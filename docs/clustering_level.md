# Why clustering runs on segmented cells, not 8 or 16 um bins

Decided 2026-10-06. Measurements below are from the integrations run that day on
the six capture areas (24 TMA cores, 22 `sample_id` batches, 6 patients).

## Summary

The workflow integrates (Harmony) and clusters (Leiden) only the segmented cells
from Space Ranger's `segmented_outputs`. The 8 and 16 um bins are not embedded
or clustered on their own. Each bin takes its cluster from the cells that cover
it, and 16 um bins also get a spatially smoothed tissue `domain`.

Bins hold too few counts to integrate. A median 16 um bin has 15 counts across
the 2000 highly variable genes (HVGs) that the embedding is built from, and a
median 8 um bin has 8. At that depth, the PCA, Harmony and UMAP mostly reflect
how many molecules each bin caught, not what tissue it sits in. A segmented
cell has a median of 62 HVG counts and 635 total, which is enough to
integrate and cluster.

## Depth per observation

| Level | Observations | Median total counts | Median HVG counts | HVG counts = 0 |
|---|---|---|---|---|
| 8 um bins | 350,892 | 66 | 8 | 10,456 (3.0%) |
| 16 um bins | 127,421 | 122 | 15 | 1,853 (1.5%) |
| Segmented cells | 52,203 | 635 | 62 | 2 |

These are the objects after QC. QC already removed every bin under 10 counts
and the bottom decile per section (`squidpy.qc` in `config.yaml`), so the
remaining bins are the better ones.

## What went wrong at the bin levels

### 1. The embedding sorts bins by depth, not by tissue

The 16 um UMAP showed a dense central starburst of shallow bins with spokes
radiating out. Colouring by `total_counts` showed depth running smoothly across
the embedding; real structure (epidermis, cornified layer, deep dermis) showed
up only at the high-depth edges. With about 15 HVG counts per bin, most of the
2000 genes are zero in any given bin. After log-normalising and scaling, the
leading principal components mostly measure how many genes were detected at
all.

Sections differ several-fold in depth: TMA5_bottom-right and TMA6_bottom-right
are the shallowest. Depth is therefore partly a batch variable, and Harmony
cannot separate "shallow section" from "acellular tissue".

### 2. Empty bins form islands around the UMAP

The sparse clusters ringing both integrated UMAPs were 879 bins in 64 islands
(sample_id run) and 1,640 bins in 33 islands (patient_id run). Between 91% and
98% of those bins had **zero counts across all 2000 HVGs**: their total counts
(about 20 per bin) all fell in genes outside the HVG set.

The mechanism:

1. `sc.pp.scale` turns every all-zero row into the same vector. All 1,853
   zero-HVG bins at 16 um sat on **one identical point** in `X_pca`.
2. Harmony shifts that point differently for each batch. In `X_pca_harmony`
   it became 127 distinct points in the sample_id run and 29 in the patient_id
   run (at 8 um, 49 and 6).
3. UMAP cannot place dozens of exact duplicates inside the main body, so each
   copy became its own tight island, one per sample (or per patient).

The patient_id run had fewer but larger duplicate groups, and 13 of them broke
off the kNN graph entirely (653 bins). Leiden returns at least one cluster per
disconnected component, so those groups became clusters regardless of
resolution.

`squidpy.cluster.min_hvg_counts` (default 1) now drops such observations before
scaling. At the cell level it removes 2 cells.

### 3. The kNN graph fragments and small clusters can't be removed

At 8 um the kNN graph split into 61 disconnected pieces at `n_neighbors` 30, 33
at 50, 11 at 75 and only became connected at 100. Even then, 78 clusters held
under 100 bins at every setting. An earlier 16 um run without an absolute count
floor produced 110 components and 273 clusters, 213 of them under 100 bins and
each from a single sample. Raising `n_neighbors` connects the graph but does
not fix the underlying problem, which is that sparse bins do not resemble each
other.

### 4. Either choice about Harmony is wrong for bins

- **With Harmony**, the correction is estimated on embeddings dominated by
  depth (sections 1–3 above), so it aligns depth artefacts across batches as
  much as biology.
- **Without Harmony**, bin clusters split by section and depth. Because every
  `sample_id` batch is a single condition x timepoint, a section-specific
  cluster then reads as a condition difference in per-cluster DE, CCC and
  neighbourhood enrichment.

Neither gives bin clusters that mean the same thing in every core.

### 5. A bin is not a cell

A 16 um bin covers 256 um², about the footprint of one dermal cell, and its
grid is placed without regard to cell boundaries. A cluster of bins is
therefore a cluster of tissue patches, not a cell type. In the transfer below,
most bins that overlap labelled cells lie within a single cell cluster (median
purity 1.0), so this matters less here than the depth problem. It is still
why bins can describe *where* tissue is but not *which cell* it is.

## Would stricter QC or more HVGs rescue the 8 um level?

Not enough to change the decision. Both were measured on the 8 um object
(`output/squidpy/square_008um/combined.h5ad`, 350,892 bins, median 64 total
counts; 2026-10-08). The target is cell-like depth: a segmented cell has a
median of 635 total and 62 HVG counts.

### Stricter depth filtering

Fraction of bins kept by an absolute total-count floor, overall, across the 24
cores, and by the transferred `broad_cell_type`:

| Min total counts | Kept | Range across cores | Epithelial | Fibroblast | Unassigned | Median of kept |
|---|---|---|---|---|---|---|
| 10 (current) | 100% | 100% | 100% | 100% | 100% | 64 |
| 25 | 77% | 49–95% | 99% | 87% | 63% | 106 |
| 50 | 57% | 24–77% | 96% | 66% | 34% | 186 |
| 100 | 40% | 10–57% | 89% | 40% | 17% | 315 |
| 200 | 27% | 4.5–43% | 78% | 16% | 10% | 485 |
| 400 | 16% | 1.4–28% | 59% | 3.5% | 4.6% | 703 |

Reaching cell-like depth needs a floor of about 200, which keeps 27% of bins.
That remainder is mostly epidermis: 78% of epithelial bins survive against 16%
of fibroblast bins, and almost none of the acellular dermis. Some cores lose
95% of their bins while others keep 43%. How many bins survive then depends on
each section's sequencing depth, and since every `sample_id` batch is a single
condition x timepoint, that unevenness would show up as condition differences.
The old `min_counts: 50` floor already showed this (87% of TMA2_top-right
removed against 17% of TMA4_top-left; see `squidpy.qc` in `config.yaml`).

### More HVGs

HVG counts per bin as `n_top_genes` grows (genes ranked by the current
`seurat` dispersion):

| HVGs | Median HVG counts per bin | 90th percentile | Share of all counts | Bins with zero HVG counts | Bins with at least 30 |
|---|---|---|---|---|---|
| 500 | 0 | 4 | 1.0% | 53% | 1.3% |
| 1,000 | 2 | 11 | 2.6% | 21% | 2.9% |
| 2,000 (current) | 8 | 39 | 8.9% | 3.0% | 14% |
| 3,000 | 12 | 59 | 12.5% | 1.2% | 24% |
| 5,000 | 20 | 111 | 21% | 0.2% | 39% |
| 8,000 | 34 | 229 | 40% | 0% | 54% |
| all 16,828 | 64 | 620 | 100% | 0% | 71% |

The current 2000 HVGs cover only 9% of all counts. The `seurat` flavour ranks
genes by dispersion on log-normalised data, which favours sparse, low-expressed
genes, so most bins catch almost nothing in the HVG set. Going to 5000 would
remove nearly all the zero-HVG islands (0.2% of bins left) and raise the median
to 20. But the extra genes carry less biological signal, and even with every
gene a median bin has 64 counts.

Changing the selection method helps more than changing the number:
`flavor="seurat_v3"` (ranks on raw counts) or Pearson-residual HVG selection
picks higher-expressed, informative genes, so each bin carries more counts at
the same 2000. With Pearson residuals, 3000 HVGs and a floor of 25, the median rises
to 28 (see the comparison below). The cell level uses the same `seurat` flavour today,
and `squidpy.cluster.hvg_flavor` can switch it.

### Both together

A floor of 100 plus 5000 HVGs would embed about 40% of bins at roughly 50–60
HVG counts each. That is cell-like depth, but with the epidermis bias above,
and an 8 um bin covers 64 um², smaller than a cell, so its clusters would
still describe parts of cells rather than whole cells. The current design gets
more out of 8 um by using its resolution for area-weighted label transfer, and
uses neighbourhood smoothing at 16 um (domains) for regions defined directly
from bins.

### Four recipes compared at 8 um

To find the best option if clustering *had* to happen at 8 um, four recipes were
run on the same 8 um object (2026-10-08). All four integrate on `sample_id`, use
kNN 100 on 30 dimensions, Leiden at resolution 0.5 and UMAP `min_dist` 0.1.

| Arm | Filter | Features | Integration |
|---|---|---|---|
| baseline | floor 10, drop zero-HVG bins | 2000 `seurat` HVGs, log1p + scale, PCA 50 | Harmony (theta 2) |
| pearson | floor 25 | 3000 Pearson-residual HVGs (selected per `sample_id`); PCA 50 on analytic Pearson residuals, with depth from all genes | Harmony |
| pearson_banksy | floor 25 | as pearson, joined to the mean over two rings of neighbours (24 bins), BANKSY cell-typing lambda 0.2 | Harmony |
| scvi | floor 25 | same 3000 HVGs, raw counts | scVI (negative binomial, 30 latent dimensions, GPU) |

The floor of 25 keeps 270,642 of 350,892 bins (77%). With 3000 Pearson HVGs a
median bin has 28 HVG counts, against 8 for the baseline, and only 81 bins have
none (baseline: 10,456).

Each arm was scored against the cell-level clusters transferred onto the same
bins. That transfer uses each bin's overlap with segmented cells, independently
of its own expression. Only bins with `leiden_purity` of at least 0.75 were
scored, 138,000–147,000 per arm, so acellular areas are not tested.

| | baseline | pearson | pearson_banksy | scvi |
|---|---|---|---|---|
| Bins | 340,436 | 270,561 | 270,561 | 270,561 |
| Leiden clusters | 32 | 20 | 15 | 7 |
| Largest cluster | 9% | 15% | 25% | 56% |
| kNN graph components | 1 | 1 | 1 | 1 |
| ARI / NMI vs cell clusters | 0.17 / 0.28 | 0.33 / 0.50 | **0.44 / 0.56** | 0.28 / 0.50 |
| ARI / NMI vs `broad_cell_type` | 0.10 / 0.18 | 0.16 / 0.34 | 0.19 / 0.38 | 0.21 / 0.39 |
| Largest \|Spearman\| of log depth with dimensions 1–5 | 0.27 | 0.41 | 0.34 | **0.16** |
| Share of log-depth variance explained by cluster | 0.30 | 0.48 | 0.49 | 0.32 |
| kNN neighbours from the same sample (0.058 if fully mixed) | 0.12 | 0.13 | 0.16 | 0.14 |
| Clusters with >80% of bins from one sample | 0 | 0 | 2 | 1 |

Agreement with the cell clusters barely changes when every arm is scored on
the same 267,095 bins.

- **baseline**: reproduces the problem. The UMAP is a starburst ordered by depth,
  and the spatial maps are salt-and-pepper, with neighbouring bins in unrelated
  clusters.
- **pearson**: changing only the features roughly doubles agreement with the cell
  clusters. Follicles, epidermal layers and dermis appear as spatially coherent
  regions.
- **pearson_banksy**: agrees best with the cell clusters, and most of its
  clusters are 70–98% one cell cluster. It has three costs:
  - Two clusters come from single samples: cluster 13 (1,609 bins, 99% one
    sample) and cluster 14 (360 bins, 100% one sample). The 16 um domains
    produced the same kind of sample-specific clusters.
  - Mixing across samples is somewhat weaker.
  - One 68,000-bin cluster at a median depth of 49 collects shallow dermis.
- **scvi**: removes the most depth from the embedding, because it models library
  size directly. But it leaves 56% of bins (shallow dermis, median depth 66) in a
  single cluster. Training also stopped at the 60-epoch cap, so the model may be
  under-fitted.

Caveats on reading these scores:

- The arms found different numbers of clusters at the same resolution, which
  biases ARI, and each ran with a single seed.
- Depth still separates clusters in the Pearson arms (0.48 of log-depth
  variance), but part of that is biological: epidermis is deep and dermis
  shallow. These scores cannot separate the two.

**Conclusion.** If integration must happen at 8 um, use Pearson-residual
features with neighbour smoothing (and a floor of about 25). That recipe is
clearly better than the current one. Its clusters still describe roughly
40 um neighbourhoods rather than single cells, and it agrees with cell clusters
at ARI 0.44. Cell-level clustering with label transfer therefore remains the
source of cell identity. Smoothed 8 um clusters are best read as
finer-grained tissue domains.

Untested: scVI clustered at a higher resolution, scVI trained for more epochs,
and a sweep of lambda (0.2–0.5) and neighbour rings (1–2).

## Why segmented cells work

- **Depth.** Cells carry about four times the HVG counts of a 16 um bin, and
  only 2 of 52,203 have none. The embedding is organised by expression rather
  than detection.
- **Integration holds.** The cell kNN graph is a single component in both runs.
  In the sample_id run, no single sample makes up more than 27% of any cluster,
  except one cluster of 290 cells (50%). The condition difference survives: in
  some clusters up to 64% of cells come from one condition, rather than every
  cluster being 50/50.
- **Clusters are cell types.** The cell UMAP has 16 clusters in the sample_id
  run (15 in patient_id). Its small outlying groups are recognisable cell
  types, not artefacts: smooth muscle (MYH11, DES, TAGLN), inner root sheath
  (TCHH, KRT71, KRT25), hair shaft cortex (KRTAP genes) and cornified
  keratinocytes (FLG, LOR, LCE).
- **Biological unit.** Cell-type proportions, DE within a cell type and
  cell–cell communication are defined on cells.

## What happens to the bins instead

`src/squidpy_bin_labels.py` (rules `squidpy_bin_labels` and
`squidpy_bin_labels_alt`) writes the bin-level `clustered*.h5ad` that the
downstream rules read:

- **Transferred cell labels.** `outs/barcode_mappings.parquet` places every
  2 um square in its 8 um bin, its 16 um bin and, if it lies in one, a
  segmented cell. Each bin takes the cluster (plus `broad_cell_type` and
  `coarse_cell_type`) covering the most of its area. Bins where labelled cells
  cover less than `label_transfer.min_cell_fraction` (0.25) of the area are
  `unassigned`. The bin's sample_id or patient_id object takes labels from the
  cell run with the same batch key, so a cluster means the same thing at every
  level. Also written: `leiden_purity`, `cell_coverage`, and
  `obsm['leiden_fractions']` (each cluster's share of each bin).
- **Tissue domains (16 um).** Labels only reach bins that overlap cells: 30–42%
  of 16 um bins per capture area and 38–55% of 8 um bins. The rest are mostly
  acellular dermis. For full coverage, each 16 um bin's expression is joined to
  the mean of its 24 neighbours within two rings (BANKSY-style, lambda 0.8),
  which pools about a cell's worth of counts. That pooled signal is enough to
  run PCA, Harmony and Leiden, and gives `domain`. Settings are in
  `squidpy.domains`.
- **Expression unchanged.** Bins still go through the same normalise / HVG /
  scale step, so SVG, CCC and DE see the same `X` and `.raw` as before. Only
  the origin of the labels changed.

## Caveats

- Unassigned bins form their own group in neighbourhood enrichment, CCC and
  per-cluster DE. Read that group as acellular tissue or cells that failed
  QC, not as a cell type.
- 13% of in-cell 2 um squares belong to cells removed by cell-level QC; those
  areas count as uncovered.
- A few domains come almost entirely from one sample or one patient. In the
  sample_id run, domains 4, 10 and 14 are each under 50 bins, and domain 5
  has 512 bins, 97% from one sample. In the patient_id run, domain 7 has 583
  bins, 99.8% from one patient. Treat these as sample-specific regions, not
  shared tissue types.
- Segmentation quality now sets the ceiling on cell-type resolution. Cells
  are segmented from H&E, and dense regions such as the epidermis can merge
  neighbouring nuclei.

## Where to look

- Cell integration quality: `results/squidpy/cell/qc*/cluster_by_*.csv`, and
  the UMAPs in `results/squidpy/cell/figures*/`.
- Label transfer per capture area:
  `results/squidpy/{square_016um,square_008um}/qc*/label_transfer.csv`.
- Domain composition: `results/squidpy/square_016um/qc*/domain_by_*.csv`.
- Spatial panels with leiden, domain and cell type side by side:
  `results/squidpy/square_016um/spatial/*_panels.pdf`.
- The island investigation (scripts, tables, highlighted UMAPs) and the
  8 um QC/HVG sweep (`qc_hvg_sweep.py`): `output/umap_islands/` (not tracked
  in git).
- The 8 um recipe comparison: `output/umap_islands/compare_8um.py` and
  `score_8um.py` (with `run_compare_8um*.sh`). Results, including
  `summary.csv`, per-arm cluster tables, `umaps.png` and `spatial.png`, are in
  `output/umap_islands/compare_8um/`. The plaid grid in `spatial.png` is
  aliasing from drawing 8 um squares at low dpi, not a feature of the data.
