# Section QC: how bins and cells are filtered per TMA core

The numbers below come from the `squidpy_combine` outputs of 2026-10-02
(`results/squidpy/{square_016um,square_008um,cell}/qc/`), covering six capture
areas and 24 TMA cores.

## Summary

QC runs in `src/squidpy_combine.py` (rule `squidpy_combine`) after the six
capture areas are concatenated, once per level (16 um, 8 um, cells). The unit of
QC is the **section**, meaning one TMA core (`core_id`), not the capture area.

Two rules shape the approach:

1. **Sections are reported, never dropped.** Each section's depth is measured
   in its core interior and compared against a gate. A section below the gate
   is named in the log and in `section_qc.csv`, but it stays in the object.
2. **Bins are filtered mostly relative to their own section.** A low absolute
   floor removes bins that are empty everywhere. A per-section quantile trim
   then removes the bottom 10% of each section's depth distribution.

Settings live under `squidpy.qc` in `config.yaml`.

## Why the section, not the capture area

Each capture area holds four cores, and they differ in patient, condition and
timepoint. A bad core averages away against its three neighbours when QC is
grouped by `run_id`. For example, TMA2_bottom-right (interior median of 8.5
counts at 16 um) shares capture area A1_TMA2 with TMA2_top-left (354).

Every `sample_id` is also a single condition x timepoint, so a section is the
smallest unit that a depth problem can confound with biology. QC is grouped by
`squidpy.qc.group_key: core_id`.

## Order of operations

On the concatenated object:

1. **Section report** on the unfiltered object (below). The QC figures are also
   drawn from this unfiltered view, so a weak section looks weak because of its
   tissue, not because the filter thinned it.
2. **Absolute floor**: `sc.pp.filter_cells(min_counts=...)`, with the value
   taken from `min_counts_per_level` (10 at every level).
3. **`min_genes`**: 0, so this step is off. Counts and genes are collinear, and
   one depth axis is enough (see [No gene-count filter](#no-gene-count-filter)).
4. **Mitochondrial cap**: `pct_counts_mt <= 25`, which removes 0.55% of bins.
   This is a weak axis on this probe-based FFPE panel, whose 11 MT- probes do
   not track mitochondrial content the way poly-A capture does.
5. **Per-section quantile trim** of `total_counts`, lower 0.1 and upper 1.0,
   applied to whatever survived steps 2–4.
6. **Gene filter**: `min_cells: 10`, applied once to the whole object after
   the trim, so the gene set matches the bins that survive.

## 1. Section report (never enforced)

`section_qc()` reports each section's depth twice: over the whole section and
over its **interior**. The interior is the half of the bins closest to the
section's centroid.

The interior matters because most Loupe core selections include a background
rim, and its width depends on how generously the core was drawn. The overall
median mixes tissue with rim, so it cannot tell weak tissue from a wide rim:

| Section (16 um) | Overall median | Interior | Exterior | Reading |
|---|---|---|---|---|
| TMA2_bottom-right | 4 | 8.5 | 2 | weak tissue |
| TMA2_top-left | 90 | 354 | 25 | good tissue, wide rim |

Panel C of `section_depth.pdf` plots interior against exterior depth for every
section, which makes this split visible.

### The gate

`section_gate.interior_median_counts: 20` names sections whose interior median
falls below 20. At 16 um it names two sections:

| Section | Sample | Interior median | Interior bins < 10 counts |
|---|---|---|---|
| TMA2_bottom-right | M007_pigmented_W16 | 8.5 | 53% |
| TMA2_top-right | M007_pigmented_W0 | 13 | 43% |
| next lowest: TMA3_bottom-right | M008_pigmented_W16 | 29 | 24% |

Interior depth across the 24 sections is a continuum (8.5, 13, 29, 37, 39, 49,
50, 76 … 481), so 20 is a reporting threshold, not a natural break. The column
that does separate cleanly is `interior_frac_low`, the share of interior bins
under 10 counts: 0.53 and 0.43 for these two sections, against 0.26 for the
next one.

Applied unchanged at 8 um, the gate names eight sections, because 8 um bins
hold about a quarter of the counts. Treat the 16 um list as the reference.

### Why nothing is dropped automatically

Both sections below the gate are patient M007's **pigmented** cores (W0 and
W16). Dropping them would remove M007 from the paired normal-vs-pigmented
comparison entirely. That is a study-design decision for a person, not
something a threshold should decide, so the pipeline warns and keeps them:

```
Warning: 2 section(s) below the interior-depth gate of 20 counts:
TMA2_bottom-right (interior median 8), TMA2_top-right (interior median 13).
These were KEPT -- nothing is excluded automatically.
```

## 2. Absolute floor: `min_counts` 10

A purely relative filter cannot remove a bin that is uninformative in every
section. With no floor, the 16 um object kept 32k bins under 10 counts. Those
bins formed 110 disconnected pieces of the kNN graph, which forced 273 Leiden
clusters, 213 of them under 100 bins.

The floor is set low on purpose:

- **It is cheap.** At 16 um it removes 24% of bins but only 0.5% of
  sequencing (99.51% of counts kept). At 8 um it keeps 98.15% of counts, and
  at the cell level 99.08%.
- **A high floor filters sections unevenly.** The earlier `min_counts: 50`
  removed 87% of TMA2_top-right but 17% of TMA4_top-left: a 6.3x spread in
  retention across a 75x spread in raw median depth (4 to 300 counts). Because
  condition is nested in section, an uneven filter is an uneven filter on the
  contrast being tested. At 10 the spread is 2.8x, which reflects real
  differences in tissue quality rather than hiding them.

Per level (`min_counts_per_level`):

- An 8 um bin covers a quarter of a 16 um bin's area, so the same floor is
  about 4x stricter there. A floor of 3 would match 16 um on cost (about 0.5%
  of counts) and on area scaling. It is held at 10 anyway, because the floor's
  job at 8 um is graph connectivity: at 10 the 8 um kNN graph already splits
  into 61 pieces at `n_neighbors` 30.
- Cells are segmented objects with a median of about 530 counts before QC, so
  the floor rarely applies to them.

## 3. Per-section quantile trim

After the floor, each section independently loses its bottom 10% of
`total_counts`. The bounds are inclusive, so a section whose 10th percentile
falls on a heavily shared integer keeps slightly more than 90%. Kept fractions
come out at 0.900–0.917 per section at 16 um, 0.900–0.946 at 8 um and
0.898–0.905 for cells.

The resulting depth cutoff therefore varies by section:

| Level | Lower bound range across sections (counts) |
|---|---|
| 16 um | 12–38 |
| 8 um | 11–19 |
| Cell | 28–501 |

There is **no upper trim**. These are bins, not droplets, so a high-count bin
is not a doublet. It is dense tissue such as epidermis or a cluster of dermal
cells. At 16 um the top decile holds 63.7% of all UMIs in the experiment, at a
median of 1,485 genes per bin, against 0.6% of UMIs for the bottom 30%.
Trimming the top would throw away most of the signal.

## No gene-count filter

`min_genes` is 0 at every level. A gene threshold would remove almost exactly
the bins the count filters already remove, and it would bias against the same
cores and cell types. Measured on 2026-10-08 on the post-QC objects
(`output/squidpy/<level>/clustered.h5ad`).

### Genes and counts are nearly the same number

On this probe-based assay, almost every UMI in a bin is a different gene:

| Level | Median genes | Median UMIs per gene | Spearman, genes vs counts |
|---|---|---|---|
| 8 um | 60 | 1.07 | 0.998 |
| 16 um | 108 | 1.11 | 0.999 |
| Cell | 528 | 1.19 | 0.995 |

A `min_genes` threshold is therefore a `min_counts` threshold under another name.

### What extra thresholds would remove

These fractions are on top of the current QC. "Worst core" is the lowest
retention of any core. From the second row of each level onwards, that core is
TMA2_top-right (M007_pigmented_W0); at the lowest threshold it is
TMA3_bottom-right. The cell-type columns use `broad_cell_type`, transferred
from cells for the bin levels.

| Level | `min_genes` | Kept | Worst core | Epithelial | Fibroblast | No cell overlap | Normal / pigmented | UMIs lost |
|---|---|---|---|---|---|---|---|---|
| 8 um | 10 | 99.9% | 99% | 100% | 100% | 99.9% | 100 / 99.9% | 0.0% |
| 8 um | 15 | 93% | 74% | 99.7% | 97% | 89% | 94 / 92% | 0.5% |
| 8 um | 20 | 83% | 57% | 99% | 91% | 72% | 85 / 81% | 1.3% |
| 8 um | 30 | 70% | 39% | 98% | 81% | 52% | 72 / 67% | 2.9% |
| 8 um | 50 | 55% | 23% | 95% | 64% | 32% | 57 / 52% | 5.9% |
| 16 um | 15 | 98% | 85% | 100% | 99.8% | 96% | 99 / 96% | 0.1% |
| 16 um | 25 | 86% | 60% | 99.9% | 98% | 79% | 90 / 83% | 0.5% |
| 16 um | 40 | 75% | 40% | 99.6% | 95% | 61% | 78 / 70% | 1.2% |
| 16 um | 60 | 64% | 28% | 99% | 90% | 46% | 68 / 60% | 2.1% |
| 16 um | 100 | 52% | 17% | 98% | 78% | 30% | 56 / 47% | 3.9% |
| Cell | 50 | 98% | 85% | 99% | 97% | – | 99.5 / 97% | 0.1% |
| Cell | 100 | 93% | 65% | 98% | 85% | – | 97 / 87% | 0.5% |
| Cell | 150 | 86% | 45% | 96% | 74% | – | 93 / 78% | 1.3% |
| Cell | 200 | 81% | 33% | 94% | 64% | – | 89 / 71% | 2.4% |
| Cell | 300 | 70% | 17% | 89% | 47% | – | 80 / 59% | 5.3% |

The losses are small in UMIs but uneven in what they remove. Once a threshold
starts to bite, it removes the shallow cores, fibroblasts and bins with no cell
overlap first, and pigmented cores faster than normal ones. That is the same
unevenness that ruled out the old `min_counts: 50` floor. The scRNA-seq
convention of 200 genes per cell would drop 36% of fibroblasts and two thirds
of TMA2_top-right.

### Low-complexity cells are cornified keratinocytes

A gene filter could add something new by catching observations with *fewer
genes than their counts predict*, which flags degraded or ambient-dominated
droplets in scRNA-seq. Here those outliers are real tissue. 1,025 cells fall
more than 3 SD below the log-genes vs log-counts trend:

- They are deep, not shallow: median 2,804 counts against 629 overall.
- One gene makes up a median of 16% of their counts, against 2% for other
  cells.
- That gene is usually FLG, then LOR, FLG2, KRTDAP, KRT2, DMKN, ASPRV1, CDSN and
  LCE1C, which mark the granular and cornified layers.
- 83% are `broad_cell_type` Epithelial cell.

A complexity filter would delete the cornified layer.

### If a gene filter is required anyway

If a reviewer or a convention requires a gene threshold, set it at the level
the count floor already implies, where it removes at most about 2% and does
not change which cores or cell types lose observations:

| Level | `min_genes` | Kept |
|---|---|---|
| 8 um | 10 | 99.9% |
| 16 um | 15 | 98% |
| Cell | 50 | 98% |

`squidpy.qc.min_genes` is a single value applied to every level. Unlike
`min_counts`, it has no `min_genes_per_level` form, so using these three values
together would need that option added to the Snakefile first.

## Outcome

Retention is the fraction of each section's bins kept through all steps:

| Level | Before | After | Kept overall | Kept per section | Gene count after |
|---|---|---|---|---|---|
| 16 um | 185,921 | 127,421 | 69% | 31–88% | 16,861 |
| 8 um | 660,443 | 350,892 | 53% | 20–79% | 16,828 |
| Cell | 58,660 | 52,203 | 89% | 82–90% | 16,452 |

At 16 um, the sections that lose the most are the weak ones from the report.
After QC, every section's median depth is at least 33 counts:

| Section (16 um) | Sample | Kept | Median counts before → after |
|---|---|---|---|
| TMA2_bottom-right | M007_pigmented_W16 | 31% | 4 → 42 |
| TMA2_top-right | M007_pigmented_W0 | 37% | 7 → 33 |
| TMA3_bottom-right | M008_pigmented_W16 | 46% | 10 → 43 |
| TMA2_bottom-left | M007_normal_W16 | 52% | 16 → 78 |
| TMA5_bottom-right | P011_pigmented_W16 | 56% | 17 → 79 |
| … | | | |
| TMA4_top-right | P009_pigmented_W0 | 88% | 209 → 271 |

## Caveats

- **Depth is not balanced across conditions.** At 16 um, four of the five
  lowest-retention sections are pigmented, and three of those are
  pigmented W16 cores. Median retention is 69% for pigmented sections against
  76% for normal ones. QC narrows this gap without removing it, so pigmented
  W16 contrasts should be read alongside `section_qc.csv`. A depth difference
  can look like an expression difference, especially for low-expressed genes.
- **Sections below the gate have very few cells.** Before QC, TMA2_bottom-right
  has 331 segmented cells, TMA2_top-right 361 and TMA3_bottom-right 418,
  against several thousand in most other cores. Per-cluster proportions and
  per-sample pseudobulk in those samples rest on small numbers.
- **The relative trim removes different cells in different sections.** At the
  cell level the cutoff ranges from 28 counts to 501. In a deep section the
  bottom decile can include real small or transcriptionally quiet cells, such
  as lymphocytes, whereas in a shallow one it removes mostly debris. So far
  this is only a concern; no measurement has checked it.
- **The interior is a geometric proxy.** "The half of bins nearest the
  centroid" assumes a roughly round core with tissue in the middle. A
  fragmented or eccentric core can put real tissue in the "exterior".
- **The bin-level QC is upstream of label transfer.** Bins removed here cannot
  receive cell labels or domains. Cells removed here leave their 2 um squares
  uncovered: 13% of in-cell squares belong to cells removed by cell-level QC
  (see `docs/clustering_level.md`).

## Outputs

All under `results/squidpy/<level>/qc/`:

| File | Contents |
|---|---|
| `section_qc.csv` | Per-section depth: overall, interior and exterior medians, and fractions under 10 counts. Built before any filter. |
| `section_depth.pdf` | A: interior depth per section, sorted, with the gate. B: depth distribution per section. C: interior vs exterior. |
| `qc_per_section.pdf` | Spatial map of depth per section on a shared colour scale, before filtering. Each panel shows the fraction the trim keeps. This is where folds, edge dropout and background selections show up. |
| `qc_violin.pdf` | Violins of `total_counts`, `n_genes_by_counts` and `pct_counts_mt`, one per section, on a log scale, tinted by capture area. |
| `quantile_bounds.csv` | Per-section trim bounds and the fraction actually kept. |
| `bins_per_core.csv` | Per-core bin counts and median depth, genes and %MT, before and after filtering. |

The log `output/logs/squidpy_combine_<level>.log` lists the sections below the
gate and the per-section kept range.

## Changing it

- To **exclude a section**, filter it out explicitly downstream, or remove the
  core from the Loupe region export. Do this as a deliberate decision that
  names the affected patient pairing, not by raising the gate.
- `quantile_trim.enabled: false` turns off the relative trim. `lower` sets its
  stringency.
- `min_genes` is deliberately 0; see [No gene-count filter](#no-gene-count-filter)
  for the cost of raising it.
- Changing `min_counts_per_level` at 8 um changes graph connectivity, not just
  depth. See `docs/clustering_level.md` for the 8 um floor sweep (floors of 25
  to 400).
