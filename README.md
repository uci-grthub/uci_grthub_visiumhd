# KadoK Visium HD (order 0526I-16)

Snakemake workflow for this project's Visium HD data on HPC3 (SLURM). Space
Ranger counts each capture area; the repository also carries the lab's Cell
Ranger and downstream single-cell rules, which this project does not currently
use.

Tissue: six TMAs of 2 mm human skin biopsies, two capture areas per Visium HD
6.5 mm slide, three slides. See `metadata/` for the slide layout documents.

## Layout
- `Snakefile` — the workflow. All rules live here except the includes below.
- `config.yaml` — references, paths, resources, per-sample overrides.
- `rules/` — `scenic.smk` is included; `jbrowse.smk`, `numbat.smk` and
  `pileup_and_phase.smk` are not wired into the Snakefile.
- `src/` — R, Python and shell scripts called by the rules (not `scripts/`).
- `data/` — FASTQs, one directory per sequencing run.
- `metadata/` — CytAssist run folders, H&E images, slide layout spreadsheets.
- `envs/` — conda environment YAMLs and the Space Ranger pixi manifest.
- `output/` — results; `output/logs/` for per-rule logs.
- `results/` — tables and figures meant to be read, including `results/squidpy/`.

## Running

Snakemake lives in a conda environment, not a module:

```bash
source /opt/apps/mamba/24.3.0/etc/profile.d/conda.sh
conda activate snakemake_dfs
```

Dry run first. It prints the discovered capture areas and warns about anything
it could not resolve — read those before submitting:

```bash
snakemake -n
```

Submit the Visium HD counting:

```bash
snakemake --profile slurm spaceranger_all
```

The `slurm` profile is at `~/.config/snakemake/slurm`; it submits one `sbatch`
job per rule with the partition, account and memory each rule declares. Space
Ranger takes hours per capture area, so run the Snakemake process itself inside
an allocation (or a small `sbatch` wrapper), not on a login node and not in the
background of an interactive session.

`spaceranger_all` stops after counting. The default `all` target additionally
pulls in the Cell Ranger summary and the Harmony integration, which have no
inputs in this project — use `spaceranger_all` unless you have added single-cell
samples. `src/submit_workflow.sh` is stale: it refers to a `cluster.yaml` and a
`snakemake` module that do not exist here.

To rerun one capture area:

```bash
snakemake --profile slurm output/spaceranger/A1_TMA4/outs/web_summary.html
```

## Visium HD

### How capture areas are discovered

There is no sample list to maintain. `rule spaceranger_count` builds one job per
capture area by reading, for each run folder under `spaceranger.cytassist_dir`:

- the `*-run-info.csv`, for the slide ID and the sample loaded in each capture area;
- the CytAssist `.tif` for that capture area;
- FASTQs under `spaceranger.fastq_root`, matched as `<prefix>_S*_L*_R1_*.fastq.gz`
  and searched recursively, so a sample sequenced across several runs contributes
  all of its directories;
- a microscope image from `spaceranger.microscope_image_dir`, if one matches.

A capture area with no FASTQs on disk is skipped with a warning, so the workflow
can be run while transfers are still arriving; the rest appear on the next
invocation with no edit.

Slide and capture area are read by Space Ranger from the CytAssist image
metadata, not passed on the command line. The workflow passes `--slide`, `--area`
and `--override-id` only when both are set in the config for that sample, for the
case where the values entered at the instrument were wrong.

### Run ids

The run id is the sample name entered at the CytAssist, with any character
outside `[a-zA-Z0-9_-]` replaced by `_`, since that is all Space Ranger accepts
for `--id`. If one name was used on more than one capture area, each becomes
`<name>_<slide>_<area>`. Both cases are announced at parse time. Results land in
`output/spaceranger/<run id>/`.

In this project the ids are `A1_TMA2`, `A1_TMA4`, `A1_TMA6`, `D1_TMA1`,
`D1_TMA3`, `D1_TMA5`.

### Configuration

```yaml
spaceranger:
  pixi_manifest: "envs/spaceranger/pixi.toml"
  transcriptome: "/dfs8/commondata/spaceranger/ref/refdata-gex-GRCh38-2020-A"
  probe_set: "/dfs8/commondata/spaceranger/probe/Human_V2/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv"
  fastq_root: "data"
  cytassist_dir: "metadata/H_E_images_out_0526I_16/CytAssist"
  microscope_image_dir: "metadata/H_E_images_out_0526I_16/TIFF"
  microscope_image_regex: "TMA\\s*(\\d+)$"     # captures {token} from the sample name
  microscope_image_pattern: "TMA {token}.tif"  # glob; {sample} {slide} {area} {token}
  create_bam: false
  localcores: 24
  localmem: 112       # GB
  custom_bin_size:    # extra bin size in um, on top of the standard 2/8/16
  mem_mb: 120000
  runtime: 2880       # minutes
  partition: "standard"
  account: "sbsandme_lab"
  samples: {}
```

`samples:` holds per-capture-area overrides, keyed by run id, by CytAssist sample
name, or by `<slide>_<area>`. Prefer `<slide>_<area>`: it is the only key that
stays unambiguous when a sample name repeats across slides.

```yaml
  samples:
    H1-PWCYQ3Y_A1:
      fastq_prefix: "CoreSampleName"   # when the core did not use the CytAssist name
      image: "metadata/.../section.tif"
      loupe_alignment: "metadata/.../alignment.json"
      slide: "H1-PWCYQ3Y"              # slide + area together add --override-id
      area: "A1"
```

### Outputs

Per capture area, under `output/spaceranger/<run id>/outs/`: `web_summary.html`,
`metrics_summary.csv`, and `binned_outputs/square_008um/filtered_feature_bc_matrix{,.h5}`.
The 2 µm and 16 µm bins are written alongside. Logs in `output/logs/`.

Every capture area needs an H&E microscope image (found through
`spaceranger.microscope_image_*`, or set per sample with
`spaceranger.samples.<id>.image`); the workflow stops at parse time if one is
missing. Space Ranger segments nuclei from it and expands them into cells
(default 8 µm), writing cell-level counts to
`outs/segmented_outputs/filtered_feature_cell_matrix.h5` and the polygons to
`outs/segmented_outputs/cell_segmentations.geojson`, next to the usual
`binned_outputs/`.

`rule combine_spaceranger_metrics` concatenates every capture area's
`metrics_summary.csv` into `results/spaceranger_metrics_summary.csv` (one row
per capture area, prefixed with a `run_id` column). It is part of the
`spaceranger_all` target.

`spaceranger aggr` is **not** used, and cannot be: it rejects Visium HD runs at
its first stage — "Visium HD runs ... are not currently supported by spaceranger
aggr" — in 3.1.2 and in 4.x. Aggregating HD capture areas into one object happens
downstream instead, in `rule squidpy_combine` — see [Downstream: squidpy](#downstream-squidpy)
and `snakemake --profile slurm squidpy_all`.

Space Ranger emits one matrix per capture area with no notion of which tissue
core a bin came from. Assigning the TMA cores (patient, normal/pigmented,
timepoint — see `metadata/TMA_slide_image_metadata_summary_v3.pdf`) to bins is a
separate, manual step, documented in `docs/core_annotation.md`.

## Downstream: squidpy

`snakemake --profile slurm squidpy_all` builds the combined objects and runs the spatial
analyses. It is deliberately **not** part of `spaceranger_all`: counting must stay
runnable without the manual annotation step below.

### The design, and why the batch key is `sample_id`

| run_id | slide | patient | cores | samples |
|---|---|---|---|---|
| `D1_TMA1` | H1-HBFQYB4 | M004 | normal/pigmented × W0/W16 | 4 |
| `A1_TMA2` | H1-HBFQYB4 | M007 | normal/pigmented × W0/W16 | 4 |
| `D1_TMA3` | H1-PWCYQ3Y | M008 | normal/pigmented × W0/W16 | 4 |
| `A1_TMA4` | H1-PWCYQ3Y | P009 | normal/pigmented × W0/W16 | 4 |
| `D1_TMA5` | H1-ZXYKDH7 | P011 | normal/pigmented × W0/W16 | 4 |
| `A1_TMA6` | H1-ZXYKDH7 | P012 | normal/pigmented × W0, two replicates each | 2 |
| | | | **24 cores** | **22 samples** |

`run_id` is simultaneously the capture area, the CytAssist run, the slide position and
the patient — one variable with six levels — while `condition` and `timepoint_weeks`
vary *within* each capture area.

Harmony runs with one batch per **biological sample**: `squidpy.cluster.batch_key:
sample_id`, built in `src/squidpy_make_adata.py` as
`<patient_id>_<condition>_W<timepoint_weeks>`. That is 22 levels against the
experiment's 24 cores, because the two `A1_TMA6` cores of each condition are technical
replicates of the same patient 12 W0 sample and share a `sample_id`. Harmony is
therefore not asked to remove the difference between two sections of one sample — the
only purely technical variation in this design, and the one thing worth being able to
see. It is the same unit `src/squidpy_de.py` pseudobulks on (`de.sample_key:
core_index`).

The cost is not small. Every batch is a single condition × timepoint, so `condition` and
`timepoint_weeks` are perfectly nested inside the batch key, and Harmony attenuates
genuine normal-vs-pigmented and W0-vs-W16 composition differences along with the
technical ones. Three things keep that visible:

- Harmony only rewrites the PCA embedding. `rank_genes_groups` and the pyDESeq2
  pseudobulk path read `.raw` / `layers["counts"]`, so **expression-level DE is
  unaffected**. Only clustering, UMAP and differential-abundance reading are.
- `X_pca` (uncorrected) is kept next to `X_pca_harmony`, and `use_harmony: false` gives
  the uncorrected comparison.
- `qc/cluster_by_*.csv` are written every run. A cluster that belongs to one `sample_id`
  means Harmony under-corrected: raise `squidpy.cluster.harmony_theta`. A
  `cluster_by_condition.csv` that is uniform across every cluster means it
  over-corrected: lower `harmony_theta`.

`A1_TMA6` has no W16 cores, so the timepoint contrast rests on five patients.

### QC is per section

Every QC figure is grouped by **section** — one TMA core, `core_id`, 24 of them —
not by capture area (`squidpy.qc.group_key`). A capture area holds four cores from
different conditions and timepoints, so grouping by `run_id` averages a bad section away
against its three neighbours. `qc/qc_violin.pdf` puts one violin column per section,
tinted by capture area; `qc/qc_per_section.pdf` is a grid of per-section spatial depth
maps, four to a row so each row is one capture area in TMA reading order. That grid is
where a tissue fold, an edge dropout or a Loupe selection that caught background shows
up first. Both are drawn on the **untrimmed** object, so a weak section looks weak
because of its tissue rather than because the filter already removed it.

Three figures, three questions:

| figure | answers |
|---|---|
| `qc/section_depth.pdf` | is this section usable at all |
| `qc/qc_per_section.pdf` | where inside the core does depth fall off |
| `qc/qc_violin.pdf` | the three QC axes side by side (log y on the depth axes) |

### Depth is filtered relatively, not against an absolute floor

Raw median depth runs **4 → 300 counts** across the 24 sections at 16 µm. A single
`min_counts` is therefore a different real stringency in every one of them:

| core_id | raw median | kept by `min_counts: 50` | kept by the decile trim |
|---|---|---|---|
| `TMA2_top-right` | 7 | **0.131** | 0.901 |
| `TMA2_bottom-right` | 4 | 0.144 | 0.900 |
| `TMA3_bottom-right` | 10 | 0.213 | 0.821 |
| `TMA4_top-left` | 300 | **0.828** | 0.802 |

A 6.3× spread in retention. `condition` is nested inside section, so an unevenly applied
filter is an uneven filter *on the contrast being tested*.

Depth is therefore cut in two stages: a **low absolute floor**
(`squidpy.qc.min_counts_per_level`) followed by the **per-section bottom decile**
(`squidpy.qc.quantile_trim`). The floor exists because a purely relative filter cannot
remove bins that are uninformative in *every* section — without it, 32 k bins under 10
counts broke the 16 µm kNN graph into 110 components and forced 273 Leiden clusters, 213
of them holding under 100 bins. The floor is cheap because shallow bins hold almost no
data: at 16 µm it removes 24% of bins for 0.5% of the sequencing. `min_genes` stays 0 —
counts and genes are collinear and one depth axis is enough — and `max_pct_mt: 25` is
kept only because it is a different axis (it removes 0.55% of bins).

Settings whose meaning changes with bin size are per level, falling back to the flat
value when a level is unlisted:

| level | `min_counts` | `n_neighbors` | obs | components | substantial clusters |
|---|---|---|---|---|---|
| `square_016um` | 10 | 30 | 127,421 | 1 | 30 |
| `square_008um` | 10 | 100 | 350,892 | 1 | 34 |
| `cell` | 10 | 30 | 52,203 | 1 | 16 |

An 8 µm bin covers a quarter the area of a 16 µm bin and runs at a third the depth, so
the same floor is far stricter there and the same neighbourhood is far sparser: at
`n_neighbors: 30` the 8 µm graph has 61 components, needing 100 to reach one. At 8 µm
roughly 1.7% of bins still end up in low-depth, sample-specific clusters (two of them
clear 500 bins, at median depths of 17–18 counts against the object median of 64); 16 µm
has no equivalent, where the largest single-sample share is 0.236.

Two things to know about it. The realized kept fraction is **not** exactly 80%:
`total_counts` is integer and the shallow sections pile up at low values, so
`TMA2_top-right`'s 10th percentile is 1 count and it keeps 0.90. `qc/quantile_bounds.csv`
reports the actual per-section bounds and fractions. And trimming the top decile removes
**dense tissue, not doublets** — these are bins, not cells, so a high-count bin is
epidermis or a dermal cell cluster. Set `quantile_trim.upper: 1.0` to keep that tail.

### Weak sections are reported, never excluded

`qc/section_qc.csv` carries per-section depth measured in the **core interior** — the
inner half by radius from the section centroid — because the overall median is
contaminated by the background rim that most Loupe selections include, and the rim varies
with how generously each core was drawn (`TMA2_top-left` runs 354 counts inside against
25 outside). That separates the two failures an overall median conflates: weak tissue,
and a wide rim on healthy tissue.

Interior depth is a **continuum** across the 24 sections — 8.5, 13, 29, 37, 39, 49, 50,
76 … 481 — so `squidpy.qc.section_gate` is a reporting convenience rather than a boundary
found in the data. The column that does break cleanly is `interior_frac_low`: 0.53 and
0.43 for `TMA2_bottom-right` and `TMA2_top-right` against 0.26 for the next section down.
Those two were checked against the hypothesis that `A1_TMA2` had background swept inside
the selection — they do not; their own interiors run at 8.5 and 13 counts while
`TMA2_top-left` on the same capture area runs at 354. They are low-yield tissue, and
redrawing them in Loupe recovers nothing.

**Nothing is ever excluded automatically.** Dropping those two would take patient M007's
pigmented cores at both timepoints and remove M007 from the paired normal-vs-pigmented
comparison entirely. The gate names them in the log and in the QC table; the decision is
yours, and the place to act on it is a manual edit, not a threshold.

### Manual step: core annotation

Nothing downstream runs until the TMA cores are drawn by hand in Loupe and exported to
`metadata/regions/regions_tma<n>.csv` (`squidpy.core_annotation_dir` and
`core_annotation_pattern`; `<n>` is the run's `tma` in `metadata/sample_metadata.csv`).
`rule squidpy_make_adata` declares those
files as inputs, so a missing one is a clean `MissingInputException` at dry-run time
naming the capture area. Full instructions, accepted file layouts and the checks to make
afterwards are in `docs/core_annotation.md`.

The selections double as the tissue mask, and this matters: Space Ranger's tissue
detection over-called on several capture areas. `A1_TMA2` reports all 175,561 of its
16 µm bins as tissue, including bins at negative pixel coordinates, against 49,654 for
`A1_TMA4`.

### Bin sizes

Both 16 µm and 8 µm objects are built. The permutation-heavy steps — neighbourhood
enrichment, co-occurrence, cell–cell communication and image features — run at **16 µm
only**, because 8 µm is 2.49 M bins across the six capture areas against 650 k at 16 µm.
8 µm still gets QC, PCA, Harmony, Leiden, spatially variable genes and the spatial plots.

### Segmented cells

A third level, `cell`, runs the same chain — combine, cluster, spatial statistics,
SVG, cell–cell communication, image features, plots — on Space Ranger's segmented
outputs (`segmented_outputs/filtered_feature_cell_matrix.h5`), writing to
`output/squidpy/cell/` and `results/squidpy/cell/`. It is on by default
(`squidpy.cell.enabled`).

- Cell positions are the `cell_centroid` of `cell_segmentations.geojson`.
- Cores come from the same 8 µm Loupe export: each cell takes the core that holds most of
  its 2 µm bins in `outs/barcode_mappings.parquet`. About 99% of cells land in a core.
- Cells are not a lattice, so their spatial graph is k nearest neighbours
  (`squidpy.cell.coord_type: generic`, `n_neighs: 6`) instead of the 8-neighbour grid.
- Space Ranger's Azimuth labels (`broad_cell_type`, `coarse_cell_type`, `fine_cell_type`)
  are attached to `obs` and plotted.
- There are few cells: 4.9k–17k per capture area, ~59k in all, some cores under 500.

Loupe only ships a `.cloupe` for `square_008um`, so annotation happens at 8 µm and the
workflow transfers the labels: HD barcodes encode the bin size and array coordinates and
the bins nest exactly, so the 8 µm bin at (r, c) is inside the 16 µm bin at
(r // 2, c // 2).

### Cell–cell communication uses LIANA+, not `sq.gr.ligrec`

`sq.gr.ligrec` needs omnipath, which is broken against NumPy 2 in this environment
(`AttributeError: module 'numpy' has no attribute 'NAN'`), and offers only a
CellPhoneDB-style cluster permutation test. `rule squidpy_ccc` uses
[LIANA+](https://liana-py.readthedocs.io/) instead, which implements CellChat's own
scoring natively over AnnData, ships CellChatDB as a selectable resource, and adds a
per-bin spatially-informed score (`li.mt.bivariate`) — no R and no Seurat round trip.

### Differential expression (pseudobulk, pyDESeq2)

Full notes, caveats and first-run results: `docs/differential_expression.md`.

`rule squidpy_de` (target `squidpy_de_all`, also part of `squidpy_all`) runs the eight
comparisons listed in `metadata/9-30-26 Meta data for 10X visium seq data analysis.xlsx`,
declared under `squidpy.de.comparisons` in `config.yaml`, at 16 µm and cell level. Bins
are not replicates, so raw counts are summed per core and pyDESeq2 tests across cores:
once over the whole core, once per Leiden cluster (pseudobulks under
`min_obs_per_pseudobulk` bins are dropped). The two `A1_TMA6` technical replicates share
a `core_index` and are summed into one sample. log2FC is test / reference.

| xlsx | config name | cores (No) | design |
|---|---|---|---|
| 1 | `melasma_vs_normal_W0` | 2,6,10 vs 1,5,9 | paired `~ patient_id + sample_type` |
| 2 | `melasma_vs_normal_W16` | 4,8,12 vs 3,7,11 | paired |
| 3 | `melasma_W16_vs_W0` | 4,8,12 vs 2,6,10 | paired `~ patient_id + timepoint_weeks` |
| 4 | `PIH_vs_normal_W0` | 14,18,22 vs 13,17,21 | paired |
| 5 | `PIH_vs_normal_W16` | 16,20 vs 15,19 | paired |
| 6 | `PIH_W16_vs_W0` | 16,20 vs 14,18 | paired; P012 dropped (baseline only) |
| 7 | `melasma_vs_PIH_W0` | 2,6,10 vs 14,18,22 | unpaired `~ sample_type` |
| 8 | `melasma_vs_PIH_W16` | 4,8,12 vs 16,20 | unpaired |

The patient term in the paired designs also absorbs capture area, since each capture
area holds one patient. Comparisons 7 and 8 have no such term: diagnosis is confounded
with patient and capture area, with 3 vs 3 and 3 vs 2 cores, so read them as
exploratory. No integration is involved in the test itself — Harmony only corrects the
PCs — but per-cluster results assume the Leiden clusters mean the same thing in every
core (check the UMAPs by `condition` first).

The clinical columns (`patient_id`, `age`, `sex`, `skin_type`, `diagnosis`,
`sample_type` = normal / melasma / PIH, `de_group`) come from that spreadsheet:

```bash
python src/make_sample_metadata.py \
    "metadata/Visium_Samples_Slide_042426_v4 with slide orientation.xlsx" \
    -o metadata/sample_metadata.csv \
    --cytassist-dir metadata/H_E_images_out_0526I_16/CytAssist \
    --clinical-xlsx "metadata/9-30-26 Meta data for 10X visium seq data analysis.xlsx"
```

`src/squidpy_de.py` joins them by `core_id` at run time, so they need not be in the
AnnData objects. `sample_metadata.csv` is an input of `squidpy_make_adata`, though, so
regenerating it makes Snakemake rebuild everything downstream unless its mtime is put
back.

### Prerequisites

The squidpy rules run in the pixi environment `envs/squidpy/pixi.toml`, locked in
`envs/squidpy/pixi.lock`. Every package comes from conda-forge or bioconda, so it
solves without the PyPI mapping that breaks the workspace `pixi.toml`. Install it once
(about 6 minutes; it also checks the imports and warms the CellChatDB cache, since
compute nodes may have no outbound network):

```bash
sbatch scripts/install_squidpy_pixi.sh
```

Versions worth knowing: squidpy 1.8.2, scanpy 1.11.5, anndata 0.12, liana 1.10.0,
harmonypy 2.0.2, numpy 2.0.2. liana is held at 1.x (2.0 needs scanpy >= 1.12 and
anndata >= 0.13), and its `numba < 0.61` pin is what caps numpy at 2.0.

Interactive use:

```bash
PYTHONNOUSERSITE=1 pixi run --manifest-path envs/squidpy/pixi.toml python src/...
```

The older conda env `squidpy_dfs` also has everything the rules need
(`scripts/install_squidpy_extras.sh` added harmonypy, leidenalg, python-igraph and
liana 1.10) but the workflow no longer uses it.

### Outputs

Intermediates in `output/squidpy/square_<bin>um/`, tables and figures in
`results/squidpy/square_<bin>um/`:

- `per_sample/<run_id>.h5ad`, `combined.h5ad`, `clustered.h5ad`, `spatial_stats.h5ad`
- `qc/<run_id>_core_assignment.pdf`, `qc/bins_per_core.csv`, `qc/qc_violin.pdf`,
  `qc/qc_per_section.pdf`, `qc/section_depth.pdf`, `qc/section_qc.csv`,
  `qc/quantile_bounds.csv`, `qc/cluster_by_<annotation>.csv`
- `markers.csv`, `figures/umap_*.pdf`
- `nhood_enrichment_zscore.csv`, `interaction_matrix.csv`, `centrality_scores.csv`,
  `cooccurrence_<group>_<level>.csv`
- `svg/<run_id>_spatial_autocorr.csv`, `spatial/<run_id>_panels.pdf`
- `ccc/<condition>_liana_cellchat.csv`, `ccc/<condition>_dotplot.pdf`
- `image_features/<run_id>.csv`
- `de/<comparison>/whole_core.csv`, `per_cluster.csv`, `samples.csv`, `volcano.pdf`

`src/explain_co-occurence.md` explains how to read the co-occurrence score. Co-occurrence
radii are set in microns (`squidpy.cooccurrence.max_radius_um`, default 500) rather than
left to squidpy's data-range default, which on a 6.5 mm capture area never samples below
a few hundred microns.

## Requirements

- **Space Ranger 4.1.0** — there is **no module** for it on HPC3, and the
  cluster-wide install at `/dfs8/commondata/spaceranger` stops at 3.1.2, which has
  no cell segmentation. The workflow calls it through `envs/spaceranger/pixi.toml`,
  a pixi environment that puts the install at
  `/dfs9/ucightf-lab/kstachel/TOOLS/spaceranger-4.1.0` on PATH. Space Ranger is proprietary and is not on conda-forge or bioconda, so the
  environment installs nothing. Check it with:
  ```bash
  pixi run --manifest-path envs/spaceranger/pixi.toml spaceranger --version
  ```
  This manifest is deliberately separate from the workspace `pixi.toml`, whose
  PyPI dependencies fail to solve on this cluster.
- Snakemake 9 with the `cluster-generic` executor plugin (conda env `snakemake_dfs`).
- For the downstream spatial rules: pixi (0.65 tested) and the `envs/squidpy`
  environment — see [Prerequisites](#prerequisites). Every such rule exports
  `PYTHONNOUSERSITE=1`.
- References under `/dfs8/commondata/spaceranger/`.
- CellRanger 8.0.1 (module) for the single-cell rules only.

## Other rules in this repository

Current and used: `spaceranger_count`, `combine_spaceranger_metrics`,
`fastqc`, `multiqc`, `collect_web_summaries`, and the downstream spatial rules
`squidpy_make_adata`, `squidpy_combine`, `squidpy_cluster`, `squidpy_spatial_stats`,
`squidpy_svg`, `squidpy_ccc`, `squidpy_image_features`, `squidpy_spatial_plots`
(target: `squidpy_all`).

Single-cell rules, working but idle in this project (no Cell Ranger samples
configured): `cellranger_count`, `multi_sample_summary`,
`tenx_harmony_integration`, `tenx_scvi_integration`, `tenx_harmony_notebook`,
`plot_filtering_timeline`, `cellbender`, `scrublet`.

Legacy, and broken as written — do not expect these to run without repair:

- `scvelo` and `pileup_and_phasing` call `src/compute_velocity.py` and
  `src/pileup_and_phase.R`, neither of which is in the repository.
- `velocyto` passes `config["..."]` as literal shell text rather than an
  evaluated value.
- `rules/scenic.smk` consumes `output/loom/<sample>.loom`, which was produced by
  a `loompy` rule that has been removed.
- `rules/jbrowse.smk`, `rules/numbat.smk` and `rules/pileup_and_phase.smk` are
  not included by the Snakefile.

## Troubleshooting

- **A `config.yaml` change did not rerun anything.** The slurm profile
  (`~/.config/snakemake/slurm/config.yaml`) sets `rerun-triggers: mtime`, so only input
  file timestamps are considered. `config.yaml` is not a declared input of any rule, so
  changing a threshold, a batch key or a quantile bound is invisible to it and the stale
  outputs are kept. A plain `snakemake -n` disagrees, because without the profile the
  default triggers include `params` and `code` — so the dry run promises a rerun the
  profile will not perform. Force it by name:

  ```bash
  sbatch scripts/run_snakemake.sh output/squidpy/square_016um/clustered.h5ad \
      --forcerun squidpy_combine
  ```

  Editing a `src/*.py` script does rerun its rule: every squidpy rule declares its script
  under `input:` rather than `params:` for exactly this reason.

- **A capture area is missing from the dry run.** The warning says which piece
  was not found. "No FASTQs" usually means the sequencing core named them
  differently from the CytAssist sample; set `fastq_prefix` for that
  `<slide>_<area>`.
- **"no microscope image".** The run proceeds on the CytAssist image alone, at
  much lower resolution (3200×3000 against roughly 10000×10000). Fix
  `microscope_image_pattern`/`microscope_image_regex`, or set `image` for that
  sample.
- **Tissue looks misaligned in `web_summary.html`.** Check the slide and capture
  area in the CytAssist run info against the physical layout. The slide
  spreadsheet in `metadata/` numbers the areas in the opposite order from the
  CytAssist run info, so confirm which is correct before trusting an alignment.
- **Re-running a capture area** deletes and recreates its output directory;
  Space Ranger does not resume a partial run here.
- **A package looks installed but behaves like an older version, or raises an
  `AttributeError` for something that should exist.**
  `~/.local/lib/python3.12/site-packages/` shadows the conda envs for roughly forty
  packages — `numpy`, `scanpy`, `anndata`, `matplotlib`, `h5py`, `numba`, `sklearn`,
  `umap`, `scipy`, `pandas` and more — and it wins over whatever env is active. Worse,
  `leidenalg` and `igraph` are present *only* there, so `sc.tl.leiden` looks available
  and then has no backend inside a clean env. Every squidpy rule exports
  `PYTHONNOUSERSITE=1`; set it by hand in any interactive session too:
  ```bash
  export PYTHONNOUSERSITE=1
  python -c "import numpy, scanpy; print(numpy.__file__, scanpy.__file__)"
  ```
  Both paths should be under `envs/squidpy/.pixi`, not under `~/.local`.
- **`squidpy_make_adata` says a core annotation is missing.** That file is drawn by hand;
  see `docs/core_annotation.md`. The dry run names every capture area still waiting.
- **Cores land in the wrong quadrant.** Check
  `results/squidpy/square_016um/qc/<run_id>_core_assignment.pdf` against
  `metadata/TMA_slide_image_metadata_summary_v3.pdf` before running the rest — a mirrored
  or rotated slide silently swaps normal with pigmented.
- **`squidpy_image_features` refuses an image.** It compares the TIFF against the size
  implied by `tissue_hires_scalef`; a mismatch means the capture area was registered
  against a different image. Bins whose crop window leaves the image are dropped with a
  count, which is normal — the 6.5 mm capture area is wider than these microscope images.
- **`squidpy_ccc` writes an empty table.** `expr_prop` (default 0.1) is usually too
  strict for sparse bins; lower it in `config.yaml`.
- Per-rule SLURM logs are under `logs/<rule>/`; Space Ranger's own output is in
  `output/logs/spaceranger_<run id>.log`, and the squidpy rules log to
  `output/logs/squidpy_*.log`.

## Documentation
- `rules/README.md` — downstream analysis rules.
- `docs/core_annotation.md` — drawing the TMA cores in Loupe, the required file
  layout, and the checks to make afterwards.
- `docs/prep_geo_data.md` — GEO submission staging.
