#!/usr/bin/env python3
"""
Rank spatially variable genes within one capture area.

``sq.gr.spatial_autocorr`` scores each gene by how much its expression clusters in space
rather than being scattered at random: Moran's I for spatial autocorrelation, Geary's C
for the complementary local-dissimilarity view. Unlike the neighbourhood and
co-occurrence statistics, this is linear in the number of bins, so it runs at both bin
sizes.

Each capture area is scored separately. Pooling them would let a gene that simply differs
between patients look spatially structured, since capture area and patient are the same
variable here.
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq

from squidpy_spatial_plots import plot_panels


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="clustered h5ad (all capture areas)")
    parser.add_argument("--run-id", required=True, help="capture area to score")
    parser.add_argument("--output-csv", required=True, help="ranked gene table to write")
    parser.add_argument("--figure-dir", help="directory for top-gene spatial plots")
    parser.add_argument("--library-key", default="run_id")
    parser.add_argument("--mode", default="moran", choices=["moran", "geary", "both"])
    parser.add_argument("--n-perms", type=int, default=100)
    parser.add_argument("--n-top-genes", type=int, default=20,
                        help="how many top genes to plot")
    parser.add_argument("--n-neighs", type=int, default=8)
    parser.add_argument("--coord-type", default="grid", choices=["grid", "generic"])
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)
    if args.library_key not in adata.obs:
        raise SystemExit(f"{args.input} has no {args.library_key} column")
    subset = adata[adata.obs[args.library_key] == args.run_id].copy()
    if subset.n_obs == 0:
        raise SystemExit(
            f"{args.run_id} is not present in {args.input}; available: "
            f"{sorted(adata.obs[args.library_key].unique())}"
        )
    del adata
    print(f"{args.run_id}: {subset.n_obs} bins x {subset.n_vars} genes")

    # the combined object's graph spans every capture area, so rebuild it here for this
    # one alone
    sq.gr.spatial_neighbors(subset, coord_type=args.coord_type, n_neighs=args.n_neighs)

    modes = ["moran", "geary"] if args.mode == "both" else [args.mode]
    frames = []
    for mode in modes:
        print(f"spatial_autocorr ({mode}, {args.n_perms} permutations)")
        sq.gr.spatial_autocorr(subset, mode=mode, n_perms=args.n_perms,
                               n_jobs=args.n_jobs, seed=args.seed,
                               show_progress_bar=False)
        frame = subset.uns[f"{mode}I" if mode == "moran" else f"{mode}C"].copy()
        frame["mode"] = mode
        frame.index.name = "gene"
        frames.append(frame.reset_index())

    import pandas as pd

    table = pd.concat(frames, ignore_index=True)
    table.insert(0, "run_id", args.run_id)
    os.makedirs(os.path.dirname(args.output_csv) or ".", exist_ok=True)
    table.to_csv(args.output_csv, index=False)
    print(f"Wrote {args.output_csv}")

    if args.figure_dir:
        primary = modes[0]
        ranked = table[table["mode"] == primary]
        score = "I" if primary == "moran" else "C"
        top = ranked.sort_values(score, ascending=(primary == "geary"))
        genes = top["gene"].head(args.n_top_genes).tolist()
        os.makedirs(args.figure_dir, exist_ok=True)
        path = os.path.join(args.figure_dir, f"{args.run_id}_top_{primary}.pdf")
        figure = plot_panels(subset, genes, args.library_key, size=2, ncols=4,
                             panel_inches=4.0)
        figure.savefig(path, dpi=150, bbox_inches="tight")
        plt.close("all")
        print(f"Wrote {path}")
        print("Top genes: " + ", ".join(genes[:10]))


if __name__ == "__main__":
    main()
