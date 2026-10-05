#!/usr/bin/env python3
"""
Spatial graph statistics for the clustered Visium HD object.

Builds the neighbourhood graph per capture area (``library_key="run_id"``) so no edge
ever crosses a capture area, let alone a TMA core, then runs:

* ``sq.gr.nhood_enrichment`` -- which clusters sit next to which, as a permutation
  z-score;
* ``sq.gr.interaction_matrix`` -- the raw count of edges between cluster pairs;
* ``sq.gr.centrality_scores`` -- degree, closeness and clustering coefficient per cluster;
* ``sq.gr.co_occurrence`` -- the enrichment of one cluster near another as a function of
  radius, computed separately per group so normal and pigmented cores can be compared.
  ``src/explain_co-occurence.md`` explains how to read the score.

Visium HD bins are a square lattice, so the default ``--n-neighs 8`` is the Moore
neighbourhood (the eight touching bins), not the hexagonal 6 that Visium v2 uses.

``co_occurrence`` is bound by pairwise distances, so each group is subsampled to
``--max-bins-per-group`` before it runs. The other statistics use the full graph.

Left to itself, ``co_occurrence`` spreads its radii over the full span of the
coordinates, which for a 6.5 mm capture area means the smallest radius is already a few
hundred microns and the direct-neighbour scale is never sampled. The radii are therefore
set explicitly, in microns, and converted to pixels with ``microns_per_pixel`` from the
Space Ranger scale factors.
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq


def build_graph(adata, library_key, coord_type, n_neighs, n_rings):
    """Neighbour graph, built independently within each capture area."""
    if library_key in adata.obs:
        adata.obs[library_key] = adata.obs[library_key].astype("category")
        print(f"Building the graph per {library_key} "
              f"({adata.obs[library_key].nunique()} libraries)")
    else:
        print(f"No {library_key} column; building one graph over all bins")
        library_key = None
    sq.gr.spatial_neighbors(
        adata, library_key=library_key, coord_type=coord_type,
        n_neighs=n_neighs, n_rings=n_rings,
    )
    edges = adata.obsp["spatial_connectivities"].nnz
    print(f"Graph: {adata.n_obs} nodes, {edges} directed edges")
    return adata


def microns_per_pixel(adata):
    """The full-resolution pixel size, from any library's Space Ranger scale factors."""
    for library in adata.uns.get("spatial", {}).values():
        value = library.get("scalefactors", {}).get("microns_per_pixel")
        if value:
            return float(value)
    raise SystemExit(
        "no microns_per_pixel in uns['spatial'][*]['scalefactors']; pass "
        "--cooccurrence-max-radius-um 0 to fall back to squidpy's data-range default"
    )


def ensure_palette(adata, cluster_key):
    """squidpy's plots warn and fall back when the categorical has no stored colors."""
    if f"{cluster_key}_colors" not in adata.uns:
        sc.pl._utils._set_default_colors_for_categorical_obs(adata, cluster_key)


def write_matrix(matrix, labels, path, name):
    """Write a cluster x cluster result as a labelled CSV."""
    frame = pd.DataFrame(matrix, index=labels, columns=labels)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    frame.to_csv(path)
    print(f"Wrote {name} to {path}")
    return frame


def run_cooccurrence(adata, cluster_key, group_key, results_dir, figure_dir,
                     interval, n_splits, max_bins, seed, max_radius_um=None):
    """Co-occurrence per level of group_key, each on a subsample of that level."""
    if group_key not in adata.obs:
        print(f"No {group_key} column; skipping co-occurrence")
        return

    um_per_px = None
    if max_radius_um:
        um_per_px = microns_per_pixel(adata)
        # radii from one bin width out to max_radius_um, expressed in pixels
        intervals = np.linspace(0.0, max_radius_um / um_per_px, interval + 1)[1:]
        print(f"co-occurrence radii: {interval} steps to {max_radius_um} um "
              f"({um_per_px:.4f} um/px)")
    else:
        intervals = interval

    for level in adata.obs[group_key].unique():
        if pd.isna(level):
            continue
        subset = adata[adata.obs[group_key] == level].copy()
        if subset.n_obs > max_bins:
            sc.pp.subsample(subset, n_obs=max_bins, random_state=seed)
        # drop clusters that vanished in the subset, otherwise the result carries
        # all-NaN rows that break the plot
        subset.obs[cluster_key] = subset.obs[cluster_key].cat.remove_unused_categories()
        ensure_palette(subset, cluster_key)
        print(f"co_occurrence for {group_key}={level} on {subset.n_obs} bins, "
              f"{subset.obs[cluster_key].nunique()} clusters")

        sq.gr.co_occurrence(subset, cluster_key=cluster_key, interval=intervals,
                            n_splits=n_splits, show_progress_bar=False)

        occ = subset.uns[f"{cluster_key}_co_occurrence"]
        labels = list(subset.obs[cluster_key].cat.categories)
        # occ["occ"] is (cluster, cluster, interval); flatten to long form
        values = occ["occ"]
        intervals = occ["interval"][1:]
        radii = occ["interval"][1:]
        records = []
        for i, source in enumerate(labels):
            for j, target in enumerate(labels):
                for k, radius in enumerate(radii):
                    records.append((source, target, float(radius), float(values[i, j, k])))
        frame = pd.DataFrame(records, columns=["cluster", "conditioned_on",
                                               "radius_px", "score"])
        if um_per_px:
            frame["radius_um"] = frame["radius_px"] * um_per_px
        frame.insert(0, group_key, level)
        path = os.path.join(results_dir, f"cooccurrence_{group_key}_{level}.csv")
        os.makedirs(results_dir, exist_ok=True)
        frame.to_csv(path, index=False)
        print(f"Wrote {path}")

        os.makedirs(figure_dir, exist_ok=True)
        figure_path = os.path.join(figure_dir,
                                   f"cooccurrence_{group_key}_{level}.pdf")
        sq.pl.co_occurrence(subset, cluster_key=cluster_key, figsize=(12, 8),
                            save=None, dpi=150)
        plt.savefig(figure_path, dpi=150, bbox_inches="tight")
        plt.close("all")
        print(f"Wrote {figure_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="clustered h5ad")
    parser.add_argument("--output", required=True, help="h5ad with the graph and results")
    parser.add_argument("--results-dir", required=True, help="directory for CSV tables")
    parser.add_argument("--figure-dir", required=True, help="directory for figures")
    parser.add_argument("--cluster-key", default="leiden")
    parser.add_argument("--library-key", default="run_id")
    parser.add_argument("--coord-type", default="grid", choices=["grid", "generic"])
    parser.add_argument("--n-neighs", type=int, default=8,
                        help="8 = Moore neighbourhood on the HD square lattice")
    parser.add_argument("--n-rings", type=int, default=1)
    parser.add_argument("--n-perms", type=int, default=1000)
    parser.add_argument("--cooccurrence-interval", type=int, default=50,
                        help="number of radius steps")
    parser.add_argument("--cooccurrence-max-radius-um", type=float, default=500.0,
                        help="largest radius in microns; 0 uses squidpy's data-range "
                             "default, whose smallest radius is far above bin scale")
    parser.add_argument("--cooccurrence-splits", type=int, default=20)
    parser.add_argument("--max-bins-per-group", type=int, default=50000)
    parser.add_argument("--group-keys", nargs="*", default=["condition"])
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)
    print(f"Loaded {adata.n_obs} bins x {adata.n_vars} genes")
    if args.cluster_key not in adata.obs:
        raise SystemExit(f"{args.input} has no {args.cluster_key} column; run "
                         "src/squidpy_cluster.py first")
    adata.obs[args.cluster_key] = adata.obs[args.cluster_key].astype("category")
    labels = list(adata.obs[args.cluster_key].cat.categories)
    ensure_palette(adata, args.cluster_key)

    adata = build_graph(adata, args.library_key, args.coord_type,
                        args.n_neighs, args.n_rings)

    os.makedirs(args.results_dir, exist_ok=True)
    os.makedirs(args.figure_dir, exist_ok=True)

    print(f"nhood_enrichment ({args.n_perms} permutations)")
    sq.gr.nhood_enrichment(adata, cluster_key=args.cluster_key, n_perms=args.n_perms,
                           n_jobs=args.n_jobs, seed=args.seed, show_progress_bar=False)
    result = adata.uns[f"{args.cluster_key}_nhood_enrichment"]
    write_matrix(result["zscore"], labels,
                 os.path.join(args.results_dir, "nhood_enrichment_zscore.csv"),
                 "nhood enrichment z-scores")
    write_matrix(result["count"], labels,
                 os.path.join(args.results_dir, "nhood_enrichment_count.csv"),
                 "nhood enrichment counts")
    sq.pl.nhood_enrichment(adata, cluster_key=args.cluster_key, figsize=(8, 8))
    plt.savefig(os.path.join(args.figure_dir, "nhood_enrichment.pdf"),
                dpi=150, bbox_inches="tight")
    plt.close("all")

    print("interaction_matrix")
    sq.gr.interaction_matrix(adata, cluster_key=args.cluster_key, normalized=False)
    write_matrix(adata.uns[f"{args.cluster_key}_interactions"], labels,
                 os.path.join(args.results_dir, "interaction_matrix.csv"),
                 "interaction matrix")

    print("centrality_scores")
    sq.gr.centrality_scores(adata, cluster_key=args.cluster_key, n_jobs=args.n_jobs)
    centrality = adata.uns[f"{args.cluster_key}_centrality_scores"]
    centrality.to_csv(os.path.join(args.results_dir, "centrality_scores.csv"))
    print(f"Wrote {os.path.join(args.results_dir, 'centrality_scores.csv')}")
    sq.pl.centrality_scores(adata, cluster_key=args.cluster_key, figsize=(12, 4))
    plt.savefig(os.path.join(args.figure_dir, "centrality_scores.pdf"),
                dpi=150, bbox_inches="tight")
    plt.close("all")

    for group_key in args.group_keys:
        run_cooccurrence(adata, args.cluster_key, group_key, args.results_dir,
                         args.figure_dir, args.cooccurrence_interval,
                         args.cooccurrence_splits, args.max_bins_per_group, args.seed,
                         args.cooccurrence_max_radius_um)

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    adata.write_h5ad(args.output, compression="gzip")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
