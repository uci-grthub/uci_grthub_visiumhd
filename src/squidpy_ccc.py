#!/usr/bin/env python3
"""
Cell-cell communication with LIANA+, scoring CellChat's formulation against CellChatDB.

``sq.gr.ligrec`` is not used. It depends on omnipath, which is broken against NumPy 2 in
this environment (``AttributeError: module 'numpy' has no attribute 'NAN'``), and it only
offers a CellPhoneDB-style cluster permutation test. LIANA+ (Dimitrov et al., Nature Cell
Biology 2024) implements CellChat's own scoring natively in Python over AnnData, ships
CellChatDB as a selectable resource, and adds a spatially-informed local score -- all
without an R dependency or a Seurat round-trip.

Two views are produced:

* **cluster level** -- ``li.mt.cellchat``: for every (source cluster, target cluster,
  ligand-receptor pair), a communication probability and a permutation p-value.
* **bin level** (``--spatial``) -- ``li.mt.bivariate``: a local interaction score per bin,
  using the spatial graph rather than cluster labels. This is the view that uses the HD
  resolution rather than averaging it away.

The CellChatDB resource is downloaded on first use. Compute nodes may have no outbound
network, so warm the cache on a login node first:

    PYTHONNOUSERSITE=1 python -c "import liana as li; li.resource.select_resource('cellchatdb')"
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc


def subset_to_group(adata, group_key, level):
    """Restrict to one level of the grouping variable, e.g. condition == 'pigmented'."""
    if group_key not in adata.obs:
        raise SystemExit(f"no {group_key} column in the input")
    mask = adata.obs[group_key].astype(str) == str(level)
    if not mask.any():
        raise SystemExit(
            f"{group_key}={level} matches no bins; available: "
            f"{sorted(adata.obs[group_key].astype(str).unique())}"
        )
    subset = adata[mask].copy()
    print(f"{group_key}={level}: {subset.n_obs} bins")
    return subset


def run_cluster_level(adata, groupby, resource, n_perms, expr_prop, n_jobs, seed):
    """li.mt.cellchat over cluster pairs; returns the flat liana_res frame."""
    import liana as li

    adata.obs[groupby] = adata.obs[groupby].astype("category").cat.remove_unused_categories()
    print(f"cellchat on {adata.obs[groupby].nunique()} {groupby} groups "
          f"against {resource} ({n_perms} permutations)")
    li.mt.cellchat(
        adata, groupby=groupby, resource_name=resource, expr_prop=expr_prop,
        n_perms=n_perms, use_raw=False, seed=seed, n_jobs=n_jobs, verbose=True,
    )
    return adata.uns["liana_res"]


def run_bin_level(adata, resource, local_name, n_perms, output, n_jobs, seed,
                  library_key="run_id", coord_type="grid", n_neighs=8):
    """li.mt.bivariate: a local LR score per bin or cell, over the spatial graph."""
    import liana as li

    if "spatial_connectivities" not in adata.obsp:
        import squidpy as sq

        print("no spatial graph on the input; building one per capture area")
        library_key = library_key if library_key in adata.obs else None
        sq.gr.spatial_neighbors(adata, library_key=library_key, coord_type=coord_type,
                                n_neighs=n_neighs)

    print(f"bivariate local scores ({local_name}, {n_perms} permutations)")
    local = li.mt.bivariate(
        adata, resource_name=resource, local_name=local_name, mask_negatives=True,
        add_categories=True, n_perms=n_perms, use_raw=False, seed=seed,
        connectivity_key="spatial_connectivities", verbose=True,
        x_name="ligand", y_name="receptor",
    )
    if local is None:
        local = adata.obsm.get("local_scores")
    if local is None:
        raise SystemExit("li.mt.bivariate returned no local scores")
    os.makedirs(os.path.dirname(output) or ".", exist_ok=True)
    local.write_h5ad(output, compression="gzip")
    print(f"Wrote {local.n_obs} bins x {local.n_vars} interactions to {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="clustered h5ad")
    parser.add_argument("--output-csv", required=True, help="liana_res table to write")
    parser.add_argument("--group-key", default="condition",
                        help="obs column naming the comparison, e.g. condition")
    parser.add_argument("--condition", required=True,
                        help="the level of --group-key to score")
    parser.add_argument("--groupby", default="leiden",
                        help="obs column giving the interacting groups")
    parser.add_argument("--resource", default="cellchatdb")
    parser.add_argument("--n-perms", type=int, default=1000)
    parser.add_argument("--expr-prop", type=float, default=0.1,
                        help="minimum fraction of bins in a group expressing the gene")
    parser.add_argument("--figure", help="dotplot of the top interactions")
    parser.add_argument("--top-n", type=int, default=25, help="interactions to plot")
    parser.add_argument("--spatial-output",
                        help="h5ad of per-bin local scores; enables --spatial")
    parser.add_argument("--local-name", default="cosine")
    parser.add_argument("--local-n-perms", type=int, default=100)
    parser.add_argument("--library-key", default="run_id",
                        help="spatial graph is built within each level of this column")
    parser.add_argument("--coord-type", default="grid", choices=["grid", "generic"],
                        help="grid for HD bins, generic (kNN) for segmented cells")
    parser.add_argument("--n-neighs", type=int, default=8)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)
    print(f"Loaded {adata.n_obs} observations x {adata.n_vars} genes")
    subset = subset_to_group(adata, args.group_key, args.condition)
    del adata

    frame = run_cluster_level(subset, args.groupby, args.resource, args.n_perms,
                              args.expr_prop, args.n_jobs, args.seed)
    frame = frame.copy()
    frame.insert(0, args.group_key, args.condition)
    os.makedirs(os.path.dirname(args.output_csv) or ".", exist_ok=True)
    frame.to_csv(args.output_csv, index=False)
    print(f"Wrote {len(frame)} interactions to {args.output_csv}")
    if frame.empty:
        print("Warning: no interactions passed the filters. --expr-prop is often too "
              "strict for sparse bins; try lowering it.")

    if args.figure and not frame.empty:
        import liana as li

        li.pl.dotplot(
            adata=subset, colour="lr_probs", size="cellchat_pvals", inverse_size=True,
            orderby="cellchat_pvals", orderby_ascending=True, top_n=args.top_n,
            figure_size=(11, 9),
        )
        os.makedirs(os.path.dirname(args.figure) or ".", exist_ok=True)
        plt.savefig(args.figure, dpi=150, bbox_inches="tight")
        plt.close("all")
        print(f"Wrote {args.figure}")

    if args.spatial_output:
        run_bin_level(subset, args.resource, args.local_name, args.local_n_perms,
                      args.spatial_output, args.n_jobs, args.seed,
                      args.library_key, args.coord_type, args.n_neighs)


if __name__ == "__main__":
    main()
