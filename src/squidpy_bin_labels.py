#!/usr/bin/env python3
"""
Label Visium HD bins with the clusters of the segmented cells they overlap.

Bins are not integrated or clustered on their own. At 16 um the median bin holds 15
counts across the 2000 HVGs, too few for Harmony to correct: the embedding orders bins
by depth, and bins with no HVG counts collapse to one PCA point that Harmony copies
once per batch into islands around the UMAP. Cells (median ~530 counts) integrate and
cluster cleanly, so the cell level is the only one integrated, and every bin takes its
labels from the cells inside it.

The overlap comes from Space Ranger's ``outs/barcode_mappings.parquet``, which places
each 2 um square in its 8 um and 16 um bin and, when it lies in a segmented cell, that
cell. For every bin and every transferred key:

* ``<key>``: the label covering the most of the bin's area, among cells that carry a
  label in the cell object (cells dropped by cell-level QC carry none). Bins whose
  labelled-cell coverage is below ``--min-cell-fraction`` of their area are
  ``unassigned``: acellular matrix, lumen, or tissue whose cells failed QC.
* ``<key>_purity``: the winning label's share of the labelled-cell area in the bin.
* ``cell_coverage``: labelled-cell area over bin area.
* ``obsm['<cluster_key>_fractions']``: area fraction of every cluster in every bin
  (rows sum to ``cell_coverage``), for anything that wants mixtures, not a vote.

Expression is prepared exactly as ``src/squidpy_cluster.py`` does it (normalise, log,
HVGs, scale, PCA), so the downstream rules see the same ``X`` and ``.raw`` they did
when bins were clustered directly; only the labels change origin.

With ``--domains``, bins also get a ``domain``: a tissue region from their own
expression plus their neighbourhood's, which covers the bins no cell reaches (most of
the acellular dermis). It follows BANKSY: each bin's scaled HVG expression is joined
to the scaled mean of its spatial neighbours within ``--domain-n-rings`` rings,
weighted ``sqrt(1 - lambda)`` and ``sqrt(lambda)``, then PCA, Harmony on the same
batch key as the cell run, kNN and Leiden. Averaging ~25 bins at two rings gives each
bin roughly a cell's worth of counts, which is what makes it integrable when single
bins are not.
"""

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import scanpy as sc
import scipy.sparse as sp

from squidpy_cluster import integrate, preprocess, write_crosstabs

UNASSIGNED = "unassigned"
UNASSIGNED_COLOR = "#d3d3d3"


def run_id_from_mapping(path):
    """output/spaceranger/<run_id>/outs/barcode_mappings.parquet -> <run_id>"""
    return Path(path).parts[-3]


def read_overlaps(mapping_paths, bin_column):
    """One row per 2 um square that lies in a cell: bin and cell, both run-suffixed."""
    frames = []
    for path in mapping_paths:
        run_id = run_id_from_mapping(path)
        table = pq.read_table(path, columns=[bin_column, "cell_id"],
                              filters=[("in_cell", "==", True)]).to_pandas()
        frames.append(pd.DataFrame({
            "bin": table[bin_column].to_numpy(dtype=object) + f"-{run_id}",
            "cell": table["cell_id"].to_numpy(dtype=object) + f"-{run_id}",
        }))
        print(f"{run_id}: {len(table)} 2 um squares inside cells")
    return pd.concat(frames, ignore_index=True)


def vote(overlaps, bins, labels, squares_per_bin, min_cell_fraction):
    """Area-weighted label per bin from the labelled cells it overlaps.

    Returns the winning label (UNASSIGNED below the coverage floor), its purity, the
    labelled-cell coverage, and the per-label area-fraction table.
    """
    labelled = overlaps.assign(label=labels.reindex(overlaps["cell"]).to_numpy())
    labelled = labelled[labelled["label"].notna() & labelled["bin"].isin(bins)]
    area = (labelled.groupby(["bin", "label"], observed=True).size()
            .unstack(fill_value=0).reindex(index=bins, fill_value=0))
    covered = area.sum(axis=1)
    coverage = covered / squares_per_bin
    winner = area.idxmax(axis=1).where(covered > 0)
    purity = (area.max(axis=1) / covered.replace(0, np.nan))
    winner = winner.where(coverage >= min_cell_fraction, UNASSIGNED).fillna(UNASSIGNED)
    return winner, purity, coverage, area / squares_per_bin


def as_category(values, cell_categories):
    """Keep the cell object's category order, with unassigned last."""
    categories = [c for c in cell_categories if c in set(values)] + [UNASSIGNED]
    return pd.Categorical(values, categories=categories)


def neighbourhood_domains(adata, batch_key, library_key, n_rings, lam, n_pcs,
                          n_neighbors, resolution, theta, crosstab_dir, cluster_key):
    """BANKSY-style domains: own plus neighbourhood-mean expression, integrated, Leiden.

    Expects ``adata`` as ``preprocess`` leaves it: scaled HVGs in X, log-normalised
    counts for every gene in ``.raw``. Writes obs['domain'] and obsm['X_domain_pca'].
    """
    import squidpy as sq

    graph_holder = sc.AnnData(obs=adata.obs[[library_key]].copy(),
                              obsm={"spatial": adata.obsm["spatial"]})
    sq.gr.spatial_neighbors(graph_holder, library_key=library_key, coord_type="grid",
                            n_neighs=8, n_rings=n_rings)
    graph = graph_holder.obsp["spatial_connectivities"].tocsr().astype(np.float32)
    degree = np.asarray(graph.sum(axis=1)).ravel()
    print(f"Domains: {n_rings}-ring grid neighbourhood, median {np.median(degree):.0f} "
          f"neighbours per bin; {int((degree == 0).sum())} bins have none")
    degree[degree == 0] = 1
    lognorm = adata.raw[:, adata.var_names].X
    lognorm = sp.csr_matrix(lognorm, dtype=np.float32)
    neighbour_mean = np.asarray((sp.diags(1 / degree) @ graph @ lognorm).todense())
    neighbour_mean -= neighbour_mean.mean(axis=0)
    std = neighbour_mean.std(axis=0)
    std[std == 0] = 1
    neighbour_mean = np.clip(neighbour_mean / std, -10, 10)

    own = np.asarray(adata.X, dtype=np.float32)
    features = np.hstack([np.sqrt(1 - lam) * own, np.sqrt(lam) * neighbour_mean])
    del neighbour_mean, own

    holder = sc.AnnData(X=features, obs=adata.obs[[batch_key]].copy())
    sc.pp.pca(holder, n_comps=n_pcs, svd_solver="arpack")
    rep = integrate(holder, [batch_key], theta)
    sc.pp.neighbors(holder, n_neighbors=n_neighbors, use_rep=rep)
    sc.tl.leiden(holder, resolution=resolution, key_added="domain", flavor="igraph",
                 n_iterations=2, directed=False, random_state=0)
    adata.obs["domain"] = holder.obs["domain"].to_numpy()
    adata.obs["domain"] = adata.obs["domain"].astype("category")
    adata.obsm["X_domain_pca"] = holder.obsm[rep]
    sizes = adata.obs["domain"].value_counts().sort_index()
    print(f"{len(sizes)} domains at resolution {resolution} (lambda {lam}): "
          f"{sizes.to_dict()}")

    if crosstab_dir:
        os.makedirs(crosstab_dir, exist_ok=True)
        for key in (cluster_key, "broad_cell_type", "condition", batch_key, "core_id"):
            if key not in adata.obs:
                continue
            path = os.path.join(crosstab_dir, f"domain_by_{key}.csv")
            pd.crosstab(adata.obs["domain"], adata.obs[key],
                        normalize="index").round(5).to_csv(path)
            print(f"Wrote {path}")
        worst = pd.crosstab(adata.obs["domain"], adata.obs[batch_key],
                            normalize="index").max(axis=1).sort_values(ascending=False)
        print(f"Largest single-{batch_key} share per domain:")
        print(worst.round(3).head(10).to_string())


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="combined bin-level h5ad")
    parser.add_argument("--cell-clustered", required=True,
                        help="clustered cell-level h5ad whose labels are transferred")
    parser.add_argument("--barcode-mappings", nargs="+", required=True,
                        help="outs/barcode_mappings.parquet, one per capture area")
    parser.add_argument("--bin-size-um", type=int, required=True)
    parser.add_argument("--output", required=True, help="labelled bin-level h5ad")
    parser.add_argument("--summary-csv", required=True,
                        help="per-run coverage and assignment summary")
    parser.add_argument("--crosstab-dir", help="directory for cluster-composition tables")
    parser.add_argument("--cluster-key", default="leiden")
    parser.add_argument("--transfer-keys", nargs="*", default=[],
                        help="further categorical cell obs columns to transfer, if present")
    parser.add_argument("--batch-key", default="sample_id",
                        help="batch key the cell object was integrated on; crosstab only")
    parser.add_argument("--min-cell-fraction", type=float, default=0.25,
                        help="labelled-cell area below this share of the bin -> unassigned")
    parser.add_argument("--n-top-genes", type=int, default=2000)
    parser.add_argument("--hvg-batch-key", default=None)
    parser.add_argument("--library-key", default="run_id",
                        help="obs column the spatial graph is built within")
    parser.add_argument("--domains", action="store_true",
                        help="also call neighbourhood-smoothed tissue domains")
    parser.add_argument("--domain-n-rings", type=int, default=2)
    parser.add_argument("--domain-lambda", type=float, default=0.8,
                        help="weight of the neighbourhood mean against the bin itself")
    parser.add_argument("--domain-n-pcs", type=int, default=30)
    parser.add_argument("--domain-n-neighbors", type=int, default=30)
    parser.add_argument("--domain-resolution", type=float, default=0.3)
    parser.add_argument("--theta", type=float, default=2.0,
                        help="Harmony theta for the domain embedding")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)
    print(f"Loaded {adata.n_obs} bins x {adata.n_vars} genes")
    cells = sc.read_h5ad(args.cell_clustered, backed="r")
    keys = [args.cluster_key] + [k for k in args.transfer_keys
                                 if k in cells.obs and k != args.cluster_key]
    cell_obs = cells.obs[keys].astype("category")
    colors = {k: list(cells.uns[f"{k}_colors"]) for k in keys if f"{k}_colors" in cells.uns}
    cells.file.close()
    print(f"Transferring {keys} from {len(cell_obs)} cells")

    bin_column = f"square_{args.bin_size_um:03d}um"
    squares_per_bin = (args.bin_size_um // 2) ** 2
    overlaps = read_overlaps(args.barcode_mappings, bin_column)
    unmatched = ~overlaps["cell"].isin(cell_obs.index)
    print(f"{int(unmatched.sum())} of {len(overlaps)} in-cell squares belong to cells "
          f"absent from the cell object (dropped by cell QC)")

    bins = adata.obs_names
    for key in keys:
        labels = cell_obs[key].astype(object)
        winner, purity, coverage, fractions = vote(
            overlaps, bins, labels, squares_per_bin, args.min_cell_fraction)
        adata.obs[key] = as_category(winner.to_numpy(), list(cell_obs[key].cat.categories))
        adata.obs[f"{key}_purity"] = purity.to_numpy()
        if key == args.cluster_key:
            adata.obs["cell_coverage"] = coverage.to_numpy()
            adata.obsm[f"{key}_fractions"] = fractions.astype(np.float32)
            adata.obsm[f"{key}_fractions"].columns = \
                adata.obsm[f"{key}_fractions"].columns.astype(str)
        if key in colors:
            palette = dict(zip(cell_obs[key].cat.categories, colors[key]))
            adata.uns[f"{key}_colors"] = (
                [palette[c] for c in adata.obs[key].cat.categories[:-1]] + [UNASSIGNED_COLOR])

    cov = adata.obs["cell_coverage"]
    print(f"cell_coverage quantiles [0.1, 0.25, 0.5, 0.75, 0.9]: "
          f"{cov.quantile([0.1, 0.25, 0.5, 0.75, 0.9]).round(3).to_dict()}")
    for floor in (0.1, 0.25, 0.5):
        print(f"  bins below coverage {floor}: {(cov < floor).mean():.3f}")
    assigned = adata.obs[args.cluster_key] != UNASSIGNED
    print(f"{int(assigned.sum())} of {adata.n_obs} bins assigned "
          f"({assigned.mean():.3f}) at min_cell_fraction {args.min_cell_fraction}; "
          f"median purity of assigned bins "
          f"{adata.obs.loc[assigned, f'{args.cluster_key}_purity'].median():.3f}")

    summary = adata.obs.assign(assigned=assigned).groupby("run_id", observed=True).agg(
        n_bins=("assigned", "size"),
        fraction_assigned=("assigned", "mean"),
        median_coverage=("cell_coverage", "median"),
        median_purity=(f"{args.cluster_key}_purity", "median"),
    )
    os.makedirs(os.path.dirname(args.summary_csv) or ".", exist_ok=True)
    summary.round(4).to_csv(args.summary_csv)
    print(summary.round(3).to_string())
    print(f"Wrote {args.summary_csv}")

    # same X / .raw / PCA as a directly clustered bin object; every bin is kept, since
    # nothing here embeds them
    adata = preprocess(adata, args.n_top_genes, args.hvg_batch_key, min_hvg_counts=0)

    write_crosstabs(adata, args.cluster_key, [args.batch_key], args.crosstab_dir)

    if args.domains:
        neighbourhood_domains(adata, args.batch_key, args.library_key,
                              args.domain_n_rings, args.domain_lambda,
                              args.domain_n_pcs, args.domain_n_neighbors,
                              args.domain_resolution, args.theta, args.crosstab_dir,
                              args.cluster_key)

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    adata.write_h5ad(args.output, compression="gzip")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
