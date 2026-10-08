#!/usr/bin/env python3
"""
Concatenate the per-capture-area Visium HD objects for one bin size and apply QC.

``spaceranger aggr`` refuses Visium HD runs, so the six capture areas are combined here
instead. Each input already carries its TMA core assignment and the per-core metadata
(patient, condition, timepoint) from ``src/squidpy_make_adata.py``.

Depth is filtered by a **low absolute floor, then relatively, per section**. At 16 um
the raw median depth runs 4 to 300 counts across the 24 sections, so a high
``min_counts`` (50) keeps 13% of TMA2_top-right against 83% of TMA4_top-left -- a 6.3x
spread in retention -- and ``condition`` is nested inside section, which makes that an
uneven filter on the contrast being tested. The floor is therefore kept low (10 counts,
0.5% of the sequencing at 16 um) and only removes bins that are uninformative in every
section: without it, bins under 10 counts fragment the kNN graph. A per-section quantile
trim of what survives then keeps 0.90-0.95 of each section. See docs/section_qc.md.

``pct_counts_mt`` is kept as a filter but is a weak axis on this probe-based FFPE assay
-- the 11 MT- probes in the panel do not report mitochondrial content the way a poly-A
assay does -- and at 25% it removes 0.55% of bins overall.

Sections are judged absolutely and bins relatively: ``section_qc`` reports per-section
interior depth and names sections below the gate, but nothing is ever excluded
automatically. Dropping a section would take both pigmented cores of a patient out of the
paired comparison, which is a decision for a person rather than a threshold.

The per-core bin counts written before and after filtering are the first real look at how
much tissue each core contributed, and the place where a bad Loupe selection shows up.

QC is reported per *section* -- one TMA core, ``core_id`` -- not per capture area. Each
capture area holds four cores from different conditions and timepoints, so a single bad
section averages away against its three neighbours when the violins are grouped by
``run_id``. Two figures come out: violins of the three QC axes with one column per
section, and a grid of per-section spatial maps of sequencing depth, which is where a
tissue fold, an edge dropout or a selection that caught background actually shows.
"""

import argparse
import os

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from squidpy_core_annotation import POSITIONS


def run_id_from_path(path):
    """output/squidpy/square_016um/per_sample/<run id>.h5ad -> <run id>"""
    return os.path.splitext(os.path.basename(path))[0]


def summarise(adata, stage):
    """Per-core bin counts and median depth at one stage of filtering."""
    grouped = adata.obs.groupby("core_id", observed=True)
    frame = pd.DataFrame(
        {
            "stage": stage,
            "n_bins": grouped.size(),
            "median_counts": grouped["total_counts"].median(),
            "median_genes": grouped["n_genes_by_counts"].median(),
            "median_pct_mt": grouped["pct_counts_mt"].median(),
        }
    )
    for column in ("run_id", "sample_id", "patient", "condition", "timepoint_weeks",
                   "tma_position", "replicate"):
        if column in adata.obs:
            frame[column] = grouped[column].first()
    return frame.reset_index()


#: the three QC axes, in the order they are plotted
QC_KEYS = ("total_counts", "n_genes_by_counts", "pct_counts_mt")

#: counts below this are treated as uninformative when summarising a section
LOW_COUNT = 10


def quantile_trim(adata, group_key, metric, lower, upper):
    """Drop the tails of ``metric`` within each group, inclusive of the bounds.

    Beyond the low absolute floor, the depth filter is relative because an absolute
    floor is a different real stringency in every section: min_counts 50 kept 13% of
    TMA2_top-right against 83% of TMA4_top-left, a 6.3x spread over a 75x spread in raw
    median depth (4 to 300 counts at 16 um). ``condition`` is nested inside section, so
    an unevenly applied filter is an uneven filter on the contrast being tested.

    The realized kept fraction is computed rather than assumed to be ``upper - lower``:
    ``metric`` is integer counts and the shallow sections pile up at low values, so a
    bound can land on a value many bins share. TMA2_top-right's 10th percentile is 1
    count, so its lower trim removes almost nothing and it keeps ~0.90. Bounds are
    inclusive, which is the conservative direction.
    """
    values = adata.obs[metric]
    groups = adata.obs[group_key].astype(str)
    bounds = values.groupby(groups, observed=True).quantile([lower, upper]).unstack()
    keep = ((values >= groups.map(bounds[lower]))
            & (values <= groups.map(bounds[upper])))

    report = pd.DataFrame({
        "n_before": groups.value_counts(),
        "n_kept": groups[keep.to_numpy()].value_counts(),
        f"{metric}_lower": bounds[lower],
        f"{metric}_upper": bounds[upper],
    })
    report["n_kept"] = report["n_kept"].fillna(0).astype(int)
    report["frac_kept"] = report["n_kept"] / report["n_before"]
    return adata[keep.to_numpy()].copy(), report.reset_index(names=group_key)


def section_qc(adata, group_key):
    """Per-section depth, including the statistics only the core interior can give.

    The overall median is contaminated by the background rim that most Loupe selections
    include, and the rim's size varies with how generously each core was drawn: the inner
    half of TMA2_top-left runs at 354 counts against 25 outside it. Reporting the
    interior separately is what distinguishes weak tissue (TMA2_bottom-right, interior
    8.5) from a wide rim on healthy tissue (TMA2_top-left, interior 354).
    """
    groups = adata.obs[group_key].astype(str).to_numpy()
    counts = adata.obs["total_counts"].to_numpy(dtype=float)
    coords = adata.obsm["spatial"]

    rows = []
    for name in sorted(set(groups)):
        mask = groups == name
        section, xy = counts[mask], coords[mask]
        radius = np.linalg.norm(xy - xy.mean(axis=0), axis=1)
        inner = radius <= np.median(radius)
        rows.append({
            group_key: name,
            "n_bins": int(mask.sum()),
            "median_counts": float(np.median(section)),
            "interior_median": float(np.median(section[inner])),
            "exterior_median": float(np.median(section[~inner])),
            "interior_frac_low": float((section[inner] < LOW_COUNT).mean()),
            "frac_low": float((section < LOW_COUNT).mean()),
        })
    frame = pd.DataFrame(rows)

    for column in ("run_id", "sample_id", "condition", "timepoint_weeks"):
        if column in adata.obs:
            first = (adata.obs[[group_key, column]].astype(str)
                     .drop_duplicates(subset=[group_key])
                     .set_index(group_key)[column])
            frame[column] = frame[group_key].map(first)
    return frame


def section_order(adata, key):
    """Sections ordered by capture area, then in TMA reading order within it.

    ``POSITIONS`` is the canonical top-left, top-right, bottom-left, bottom-right order
    of the 2x2 core layout, so with four panels to a row each row of the per-section grid
    is one capture area laid out the way the TMA is.
    """
    columns = [c for c in ("run_id", "tma_position") if c in adata.obs]
    if not columns:
        return sorted(adata.obs[key].astype(str).unique())
    frame = (adata.obs[[key] + columns].astype(str)
             .drop_duplicates(subset=[key]))
    if "tma_position" in frame:
        frame["tma_position"] = pd.Categorical(
            frame["tma_position"], categories=list(POSITIONS), ordered=True)
    return frame.sort_values(columns + [key])[key].tolist()


def plot_qc_violins(adata, key, path):
    """Violins of the three QC axes, one column per section."""
    order = section_order(adata, key)
    # tint each section by its capture area, so the run-level effect stays readable
    # behind the per-section detail
    runs = (adata.obs[[key, "run_id"]].astype(str).drop_duplicates(subset=[key])
            .set_index(key)["run_id"]) if "run_id" in adata.obs else None
    if runs is not None:
        levels = sorted(runs.unique())
        colors = plt.get_cmap("tab10")(np.linspace(0, 1, 10))
        palette = [colors[levels.index(runs[name]) % 10] for name in order]
    else:
        palette = None

    figure, axes = plt.subplots(
        3, 1, figsize=(0.45 * len(order) + 3, 9), sharex=True
    )
    for axis, qc_key in zip(axes, QC_KEYS):
        sc.pl.violin(adata, qc_key, groupby=key, order=order, palette=palette,
                     rotation=90, stripplot=False, ax=axis, show=False)
        if qc_key != "pct_counts_mt":
            # linear y cannot separate a section with a median of 7 counts from one with
            # a median of 300: tails to 12,000 compress every section into the same spike
            axis.set_yscale("log")
        axis.set_title(qc_key)
        axis.set_xlabel("")
        axis.set_ylabel("")
    figure.tight_layout()
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    figure.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(figure)
    print(f"Wrote {path}")


def plot_qc_sections(adata, key, path, ncols=4, kept=None):
    """One spatial map of depth per section, on a shared colour scale.

    Drawn on the *untrimmed* object: this figure answers "is this section good tissue",
    which a pre-trim view answers and a post-trim view mostly obscures -- after trimming,
    a weak section looks sparse because of the filter rather than because of the tissue.
    ``kept`` maps section to the fraction the quantile trim retains, shown per panel.

    Four panels per row, which is the 2x2 core layout of one capture area. Raw matplotlib
    rather than ``sq.pl.spatial_scatter``: this figure wants no H&E underlay, and the
    bins of one core are a small enough patch that the per-library image machinery buys
    nothing.
    """
    order = section_order(adata, key)
    labels = adata.obs[key].astype(str).to_numpy()
    depth = np.log1p(adata.obs["total_counts"].to_numpy(dtype=float))
    # shared scale, clipped to the 1st-99th percentile so one hot section does not
    # flatten the rest
    vmin, vmax = np.percentile(depth, [1, 99])

    nrows = int(np.ceil(len(order) / ncols))
    figure, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 3.4 * nrows))
    axes = np.atleast_1d(axes).ravel()
    scatter = None
    for axis, name in zip(axes, order):
        mask = labels == name
        coords = adata.obsm["spatial"][mask]
        scatter = axis.scatter(coords[:, 0], coords[:, 1], c=depth[mask], s=0.5,
                               linewidths=0, cmap="viridis", vmin=vmin, vmax=vmax)
        obs = adata.obs.loc[mask]
        title = str(name)
        if "sample_name" in obs:
            title += f"\n{obs['sample_name'].iloc[0]}"
        title += (f"\nn={int(mask.sum())}, "
                  f"median={int(obs['total_counts'].median())} counts")
        if kept is not None and name in kept:
            title += f", kept {kept[name]:.2f}"
        axis.set_title(title, fontsize="small")
        axis.set_aspect("equal")
        axis.invert_yaxis()  # image row order
        axis.set_xticks([])
        axis.set_yticks([])
    for axis in axes[len(order):]:
        axis.axis("off")
    if scatter is not None:
        figure.colorbar(scatter, ax=axes.tolist(), shrink=0.6,
                        label="log1p(total_counts)")
    figure.suptitle(f"Sequencing depth per section ({key})")
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    figure.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(figure)
    print(f"Wrote {path}")


def plot_section_depth(adata, key, qc_frame, path, gate=None):
    """Three views of per-section depth, drawn before any trimming.

    Panel A answers "is this section usable at all", panel B answers "is it shifted
    everywhere or does it just have a bad tail", and panel C separates the two different
    failures that an overall median cannot tell apart: weak tissue sits near the origin,
    a wide background rim on healthy tissue sits far above the diagonal.
    """
    order = section_order(adata, key)
    frame = qc_frame.set_index(key).reindex(order)
    runs = frame["run_id"] if "run_id" in frame else pd.Series("all", index=frame.index)
    levels = sorted(runs.unique())
    cmap = plt.get_cmap("tab10")
    colors = [cmap(levels.index(runs[name]) % 10) for name in frame.index]
    below = frame.index[frame["interior_median"] < gate] if gate else frame.index[:0]

    figure, axes = plt.subplots(1, 3, figsize=(19, 6))

    # A: interior depth per section, sorted, with the gate drawn
    ordered = frame.sort_values("interior_median")
    axis = axes[0]
    axis.bar(range(len(ordered)), ordered["interior_median"],
             color=[colors[list(frame.index).index(n)] for n in ordered.index])
    axis.set_yscale("log")
    axis.set_xticks(range(len(ordered)))
    axis.set_xticklabels(ordered.index, rotation=90, fontsize="x-small")
    if gate:
        axis.axhline(gate, color="crimson", linestyle="--", linewidth=1,
                     label=f"gate = {gate} counts (reported, not enforced)")
        axis.legend(fontsize="small", frameon=False)
    axis.set_ylabel("median total_counts, core interior")
    axis.set_title("A. interior depth per section")

    # B: the whole distribution, not just its median
    axis = axes[1]
    groups = adata.obs[key].astype(str).to_numpy()
    counts = adata.obs["total_counts"].to_numpy(dtype=float)
    for name, color in zip(frame.index, colors):
        values = np.sort(counts[groups == name])
        fraction = np.arange(1, len(values) + 1) / len(values)
        weak = name in below
        axis.step(np.log1p(values), fraction, where="post",
                  color="crimson" if weak else color,
                  alpha=1.0 if weak else 0.35,
                  linewidth=2.0 if weak else 0.9,
                  label=name if weak else None)
    axis.set_xlabel("log1p(total_counts)")
    axis.set_ylabel("cumulative fraction of bins")
    axis.set_title("B. depth distribution per section")
    if len(below):
        axis.legend(fontsize="small", frameon=False, loc="lower right")

    # C: interior against exterior, which tells weak tissue from a wide rim
    axis = axes[2]
    axis.scatter(frame["exterior_median"], frame["interior_median"],
                 c=colors, s=45, edgecolors="black", linewidths=0.4)
    for name, row in frame.iterrows():
        axis.annotate(name, (row["exterior_median"], row["interior_median"]),
                      fontsize="xx-small", xytext=(3, 3), textcoords="offset points")
    limits = [max(frame[["exterior_median", "interior_median"]].min().min(), 1),
              frame[["exterior_median", "interior_median"]].max().max() * 1.6]
    axis.plot(limits, limits, color="grey", linestyle=":", linewidth=1)
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("median total_counts, core exterior")
    axis.set_ylabel("median total_counts, core interior")
    axis.set_title("C. interior vs exterior\n(above the line = background rim; "
                   "near the origin = weak tissue)")

    figure.tight_layout()
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    figure.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(figure)
    print(f"Wrote {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--adatas", nargs="+", required=True,
                        help="per-capture-area h5ad files for one bin size")
    parser.add_argument("--output", required=True, help="combined h5ad to write")
    parser.add_argument("--qc-table", required=True,
                        help="per-core bin counts before and after filtering")
    parser.add_argument("--qc-figure", help="per-section QC violin figure to write")
    parser.add_argument("--qc-section-figure",
                        help="per-section spatial depth figure to write")
    parser.add_argument("--qc-group-key", default="core_id",
                        help="obs column the QC figures are grouped by; one TMA core "
                             "(a section) by default")
    parser.add_argument("--min-counts", type=int, default=0,
                        help="absolute floor on total counts; 0 disables it, which is "
                             "the default because an absolute floor is a different real "
                             "stringency in every section")
    parser.add_argument("--min-genes", type=int, default=0,
                        help="absolute floor on detected genes; 0 disables it")
    parser.add_argument("--quantile-lower", type=float, default=0.0,
                        help="per-section lower quantile of --quantile-metric to drop")
    parser.add_argument("--quantile-upper", type=float, default=1.0,
                        help="per-section upper quantile of --quantile-metric to drop")
    parser.add_argument("--quantile-metric", default="total_counts",
                        help="obs column the per-section quantile trim acts on")
    parser.add_argument("--quantile-report",
                        help="per-section bounds and realized kept fraction, as CSV")
    parser.add_argument("--section-qc-table",
                        help="per-section depth including core-interior statistics")
    parser.add_argument("--section-gate-counts", type=float,
                        help="name sections whose interior median is below this; "
                             "nothing is ever excluded automatically")
    parser.add_argument("--section-depth-figure",
                        help="three-panel per-section depth figure to write")
    parser.add_argument("--min-cells", type=int, default=10,
                        help="drop genes detected in fewer bins")
    parser.add_argument("--max-pct-mt", type=float, default=25.0,
                        help="drop bins above this percent mitochondrial counts")
    args = parser.parse_args()

    adatas = {}
    for path in sorted(args.adatas):
        run_id = run_id_from_path(path)
        adatas[run_id] = sc.read_h5ad(path)
        print(f"{run_id}: {adatas[run_id].n_obs} bins")

    combined = ad.concat(
        adatas, label="library_id", index_unique="-", merge="first", uns_merge="unique"
    )
    # ad.concat drops uns["spatial"] unless every object agrees, so rebuild it: each
    # capture area contributes its own library, keyed by run id
    combined.uns["spatial"] = {
        run_id: adata.uns["spatial"][run_id]
        for run_id, adata in adatas.items()
        if "spatial" in adata.uns and run_id in adata.uns["spatial"]
    }
    # ad.concat keeps an ordered categorical only when every capture area has the same
    # levels, and A1_TMA6 has no W16 cores, so restore it (and patient) over the union
    for column, ordered in (("timepoint_weeks", True), ("patient", False)):
        if column in combined.obs:
            values = combined.obs[column].astype(object)
            combined.obs[column] = pd.Categorical(
                values, categories=sorted(values.dropna().unique()), ordered=ordered)
    # same reason, unordered: core_id has four levels per capture area by construction,
    # so concat hands back object dtype, and both the QC grouping and the Harmony batch
    # key want a categorical
    for column in ("core_id", "sample_id", "sample_name", "condition", "tma_position",
                   "core_index", "diagnosis", "sample_type", "run_id"):
        if column in combined.obs:
            combined.obs[column] = combined.obs[column].astype(str).astype("category")
    print(f"Combined: {combined.n_obs} bins x {combined.n_vars} genes")

    group_key = args.qc_group_key
    if group_key not in combined.obs:
        raise SystemExit(
            f"--qc-group-key {group_key} is not an obs column; have "
            f"{sorted(combined.obs.columns)}"
        )
    n_sections = combined.obs[group_key].astype(str).nunique()
    print(f"QC grouped by {group_key}: {n_sections} sections")

    before = summarise(combined, "before_filter")

    # the QC figures are drawn on the untrimmed object, so a weak section looks weak
    # because of its tissue rather than because the filter already removed it. Only obs
    # and the coordinates are needed, so this view costs a frame rather than a matrix.
    qc_view = ad.AnnData(obs=combined.obs.copy(),
                         obsm={"spatial": combined.obsm["spatial"].copy()})

    # section-level QC on the untrimmed object: reported, never enforced. Dropping a
    # section automatically would take both pigmented cores of a patient and remove them
    # from the paired normal-vs-pigmented comparison; that is not a threshold's call.
    sections = section_qc(combined, group_key)
    if args.section_qc_table:
        os.makedirs(os.path.dirname(args.section_qc_table) or ".", exist_ok=True)
        sections.to_csv(args.section_qc_table, index=False)
        print(f"Wrote {args.section_qc_table}")
    if args.section_gate_counts:
        weak = sections[sections["interior_median"] < args.section_gate_counts]
        if len(weak):
            listed = ", ".join(
                f"{row[group_key]} (interior median {row['interior_median']:.0f})"
                for _, row in weak.iterrows()
            )
            print(
                f"Warning: {len(weak)} section(s) below the interior-depth gate of "
                f"{args.section_gate_counts:.0f} counts: {listed}. "
                "These were KEPT -- nothing is excluded automatically."
            )

    if args.min_counts:
        sc.pp.filter_cells(combined, min_counts=args.min_counts)
    if args.min_genes:
        sc.pp.filter_cells(combined, min_genes=args.min_genes)
    combined = combined[combined.obs["pct_counts_mt"] <= args.max_pct_mt].copy()

    report = None
    if args.quantile_lower > 0 or args.quantile_upper < 1:
        if args.quantile_metric not in combined.obs:
            raise SystemExit(
                f"--quantile-metric {args.quantile_metric} is not an obs column"
            )
        combined, report = quantile_trim(combined, group_key, args.quantile_metric,
                                         args.quantile_lower, args.quantile_upper)
        print(
            f"Per-{group_key} quantile trim of {args.quantile_metric} at "
            f"[{args.quantile_lower}, {args.quantile_upper}]: kept "
            f"{report['frac_kept'].min():.3f}-{report['frac_kept'].max():.3f} per section"
        )
        if args.quantile_report:
            os.makedirs(os.path.dirname(args.quantile_report) or ".", exist_ok=True)
            report.to_csv(args.quantile_report, index=False)
            print(f"Wrote {args.quantile_report}")

    # after the trim, so the gene set matches the bins that survive
    sc.pp.filter_genes(combined, min_cells=args.min_cells)
    print(f"After filtering: {combined.n_obs} bins x {combined.n_vars} genes")

    if combined.n_obs == 0:
        raise SystemExit(
            "Every bin was filtered out. Widen --quantile-lower/--quantile-upper, lower "
            "--min-counts/--min-genes or raise --max-pct-mt."
        )

    after = summarise(combined, "after_filter")
    table = pd.concat([before, after], ignore_index=True)
    os.makedirs(os.path.dirname(args.qc_table) or ".", exist_ok=True)
    table.to_csv(args.qc_table, index=False)
    print(f"Wrote {args.qc_table}")

    kept = None
    if report is not None:
        kept = report.set_index(group_key)["frac_kept"].to_dict()
    if args.qc_figure:
        plot_qc_violins(qc_view, group_key, args.qc_figure)
    if args.qc_section_figure:
        plot_qc_sections(qc_view, group_key, args.qc_section_figure, kept=kept)
    if args.section_depth_figure:
        plot_section_depth(qc_view, group_key, sections, args.section_depth_figure,
                           gate=args.section_gate_counts)

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    combined.write_h5ad(args.output, compression="gzip")
    print(f"Wrote {args.output}")
    summary_columns = [c for c in ("core_id", "sample_id", "condition",
                                   "timepoint_weeks", "replicate", "n_bins",
                                   "median_counts") if c in table.columns]
    print(table[table["stage"] == "after_filter"][summary_columns].to_string(index=False))


if __name__ == "__main__":
    main()
