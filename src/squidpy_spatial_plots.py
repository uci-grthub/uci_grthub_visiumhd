#!/usr/bin/env python3
"""
Spatial scatter panels for one capture area.

Plots each bin at its position in the capture area, coloured by cluster and by the core
annotation, so the clustering can be read against the tissue. The four TMA cores are laid
out 2x2 within each capture area, which makes the core panel a direct check that the
Loupe selections and the ``metadata/sample_metadata.csv`` positions agree.
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq

#: obs columns worth a panel, when present.
DEFAULT_KEYS = ("leiden", "domain", "bin_leiden", "core_id", "condition", "timepoint_weeks",
                "broad_cell_type", "total_counts")


def plot_panels(adata, keys, library_key, size, ncols, panel_inches=5.0):
    """One spatial scatter per key, each on its own axes in a constrained layout.

    sq.pl.spatial_scatter's own grid spaces panels evenly, so a long categorical legend
    (23 Leiden clusters, Azimuth cell types, core ids) runs into the next panel. Drawing
    each key onto its own axes lets constrained layout reserve room for every legend and
    colorbar.
    """
    ncols = min(ncols, len(keys))
    nrows = -(-len(keys) // ncols)
    # room for a legend beside every panel
    figure, axes = plt.subplots(
        nrows, ncols, figsize=(ncols * panel_inches * 1.6, nrows * panel_inches),
        layout="constrained", squeeze=False,
    )
    for axis, key in zip(axes.flat, keys):
        # shape=None draws points rather than hexes: Visium HD bins are square and far
        # too numerous for patch rendering
        sq.pl.spatial_scatter(adata, color=key, library_key=library_key, shape=None,
                              size=size, ax=axis, legend_fontsize="small")
    for axis in axes.flat[len(keys):]:
        axis.set_visible(False)
    return figure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="clustered h5ad")
    parser.add_argument("--run-id", required=True, help="capture area to plot")
    parser.add_argument("--output", required=True, help="figure to write")
    parser.add_argument("--library-key", default="run_id")
    parser.add_argument("--color", nargs="*", default=None,
                        help="obs columns or genes to colour by")
    parser.add_argument("--size", type=float, default=2.0)
    parser.add_argument("--ncols", type=int, default=3)
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

    keys = args.color or [k for k in DEFAULT_KEYS
                          if k in subset.obs or k in subset.var_names]
    if not keys:
        raise SystemExit("nothing to colour by; pass --color")
    print(f"{args.run_id}: {subset.n_obs} bins, plotting {', '.join(keys)}")

    for key in keys:
        if key in subset.obs and str(subset.obs[key].dtype) == "category":
            subset.obs[key] = subset.obs[key].cat.remove_unused_categories()
            if f"{key}_colors" not in subset.uns:
                sc.pl._utils._set_default_colors_for_categorical_obs(subset, key)

    figure = plot_panels(subset, keys, args.library_key, args.size, args.ncols)
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    figure.savefig(args.output, dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
