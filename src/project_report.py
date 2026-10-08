#!/usr/bin/env python3
"""
Generate the project report for the Visium HD spatial analysis: references, software
versions, what each pipeline step did and with which settings, per-capture-area and
per-level QC, the differential expression comparisons, the per-core sample sheet, and
the Snakemake rules that produced the current outputs.

Laid out after the bulk RNA-seq project report (KlasH 260923, src/generate_report.py)
so the core's reports read alike.

Sources:
  - Space Ranger: <run>/_versions, <run>/outs/web_summary.html (annotation model),
    results/spaceranger_metrics_summary.csv
  - references: <transcriptome>/reference.json and the probe set's '#key=value' header
  - squidpy env versions: `pixi list --json` of envs/squidpy, written by the rule
  - settings: config.yaml, the same file the Snakefile loads
  - QC, clusters, DE: the CSVs under results/squidpy/<level>/
  - executed steps: .snakemake/metadata, one JSON record per output file. A record
    counts only while its output still exists, and records are grouped back into jobs
    by (rule, job_hash, shellcmd), since a job with several outputs leaves one each.

Run with a Python that has reportlab, pandas and PyYAML (report.python in config.yaml).
"""

import argparse
import base64
import binascii
import html
import json
import os
import re
from collections import defaultdict
from datetime import datetime

import pandas as pd
import yaml
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import landscape, letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (BaseDocTemplate, Frame, NextPageTemplate, PageBreak,
                                PageTemplate, Paragraph, Spacer, Table, TableStyle)

BLUE = colors.HexColor("#1f4788")
LEVEL_LABELS = {"cell": "Segmented cells", "square_016um": "16 µm bins",
                "square_008um": "8 µm bins"}
SHORT_LABELS = {"cell": "cells", "square_016um": "16 µm", "square_008um": "8 µm"}

REFERENCES = [
    "10x Genomics. Space Ranger 4.1: Visium HD analysis pipelines. "
    "https://www.10xgenomics.com/support/software/space-ranger",
    "Hao et al. (2021). Integrated analysis of multimodal single-cell data (Azimuth). Cell.",
    "Mölder et al. (2021). Sustainable data analysis with Snakemake. F1000Research.",
    "Virshup et al. (2024). anndata: Access and store annotated data matrices. JOSS.",
    "Wolf et al. (2018). SCANPY: large-scale single-cell gene expression data analysis. "
    "Genome Biology.",
    "Palla et al. (2022). Squidpy: a scalable framework for spatial omics analysis. "
    "Nature Methods.",
    "Korsunsky et al. (2019). Fast, sensitive and accurate integration of single-cell data "
    "with Harmony. Nature Methods.",
    "Traag et al. (2019). From Louvain to Leiden: guaranteeing well-connected communities. "
    "Scientific Reports.",
    "McInnes et al. (2018). UMAP: Uniform Manifold Approximation and Projection for "
    "Dimension Reduction. arXiv:1802.03426.",
    "Singhal et al. (2024). BANKSY unifies cell typing and tissue domain segmentation for "
    "scalable spatial omics data analysis. Nature Genetics.",
    "Moran P.A.P. (1950). Notes on continuous stochastic phenomena. Biometrika.",
    "Dimitrov et al. (2024). LIANA+ provides an all-in-one framework for cell-cell "
    "communication inference. Nature Cell Biology.",
    "Jin et al. (2021). Inference and analysis of cell-cell communication using CellChat. "
    "Nature Communications.",
    "Muzellec et al. (2023). PyDESeq2: a python package for bulk RNA-seq differential "
    "expression analysis. Bioinformatics.",
    "Love et al. (2014). Moderated estimation of fold change and dispersion for RNA-seq "
    "data with DESeq2. Genome Biology.",
]


# --------------------------------------------------------------------------- inputs

def read_jobs(metadata_dir):
    """rule -> list of jobs, each {"outputs", "endtime"}, for outputs that still exist."""
    jobs = {}
    for name in os.listdir(metadata_dir):
        try:
            path = base64.urlsafe_b64decode(name + "=" * (-len(name) % 4)).decode()
        except (binascii.Error, UnicodeDecodeError):
            continue
        if not os.path.exists(path):
            continue
        with open(os.path.join(metadata_dir, name)) as handle:
            record = json.load(handle)
        if record.get("incomplete"):
            continue
        key = (record["rule"], record.get("job_hash"), record.get("shellcmd"))
        job = jobs.setdefault(key, {"outputs": [], "endtime": None})
        job["outputs"].append(path)
        if record.get("endtime") is not None:
            job["endtime"] = max(job["endtime"] or 0, record["endtime"])
    by_rule = defaultdict(list)
    for (rule, _, _), job in jobs.items():
        by_rule[rule].append(job)
    return by_rule


def rule_descriptions(workflow_files):
    """rule -> (position, first sentence of the '# Rule: ...' comment above it)."""
    descriptions = {}
    for path in workflow_files:
        with open(path) as handle:
            lines = handle.read().splitlines()
        for index, line in enumerate(lines):
            match = re.match(r"rule (\w+):", line)
            if not match:
                continue
            block = []
            above = index - 1
            while above >= 0 and lines[above].startswith("#"):
                block.insert(0, lines[above].lstrip("#").strip())
                above -= 1
            text = " ".join(block).split("  ")[0]
            text = re.sub(r"^Rule( to)?:?\s*", "", text).strip()
            text = re.split(r"(?<=[a-z0-9)])\. ", text)[0].split(" -- ")[0].rstrip(".")
            descriptions[match.group(1)] = (len(descriptions), text[:1].upper() + text[1:])
    return descriptions


def spaceranger_runs(run_dirs):
    """run id -> {"pipelines", "martian", "model"} from _versions and the web summary."""
    runs = {}
    for run_dir in sorted(run_dirs):
        versions_file = os.path.join(run_dir, "_versions")
        if not os.path.exists(versions_file):
            continue
        with open(versions_file) as handle:
            versions = json.load(handle)
        model = None
        summary = os.path.join(run_dir, "outs", "web_summary.html")
        if os.path.exists(summary):
            with open(summary, errors="replace") as handle:
                match = re.search(r'\["Cell annotation model","([^"]+)"\]', handle.read())
            model = match.group(1) if match else None
        runs[os.path.basename(run_dir.rstrip("/"))] = {
            "pipelines": versions.get("pipelines"), "martian": versions.get("martian"),
            "model": model}
    return runs


def probe_set_header(probe_set):
    header = {}
    with open(probe_set) as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            key, _, value = line[1:].strip().partition("=")
            header[key] = value
    return header


def read_csv(path):
    return pd.read_csv(path) if os.path.exists(path) else None


# --------------------------------------------------------------------------- report

class ReportGenerator:
    def __init__(self, args):
        self.args = args
        with open(args.config) as handle:
            self.cfg = yaml.safe_load(handle)
        self.sq = self.cfg.get("squidpy", {})
        self.results = os.path.join(self.cfg["paths"].get("results", "results"), "squidpy")
        with open(args.pixi_list) as handle:
            self.packages = {p["name"]: p["version"] for p in json.load(handle)}
        with open(os.path.join(args.transcriptome, "reference.json")) as handle:
            self.reference = json.load(handle)
        self.probes = probe_set_header(args.probe_set)
        self.runs = spaceranger_runs(args.spaceranger_runs)
        self.metrics = read_csv(args.spaceranger_metrics)
        self.samples = read_csv(self.sq.get("sample_metadata", "metadata/sample_metadata.csv"))
        self.jobs = read_jobs(args.metadata_dir)
        self.descriptions = rule_descriptions(args.workflow_files)
        self.levels = [f"square_{b}um" for b in self.sq.get("bin_sizes", [])]
        if self.sq.get("cell", {}).get("enabled", True):
            self.levels = ["cell"] + self.levels
        self._styles()

    # ---- formatting helpers

    def _styles(self):
        styles = getSampleStyleSheet()
        self.h3 = styles["Heading3"]
        self.title = ParagraphStyle("T", parent=styles["Heading1"], fontSize=22,
                                    textColor=BLUE, spaceAfter=12, alignment=TA_CENTER,
                                    fontName="Helvetica-Bold")
        self.heading = ParagraphStyle("H", parent=styles["Heading2"], fontSize=14,
                                      textColor=BLUE, spaceAfter=10, spaceBefore=12,
                                      fontName="Helvetica-Bold")
        self.body = ParagraphStyle("B", parent=styles["BodyText"], fontSize=10,
                                   alignment=TA_LEFT, spaceAfter=10, leading=14)
        self.cell = ParagraphStyle("C", parent=styles["BodyText"], fontSize=8, leading=10)
        self.path = ParagraphStyle("P", parent=styles["BodyText"], fontSize=7, leading=9,
                                   fontName="Courier", wordWrap="CJK")
        self.ref = ParagraphStyle("R", parent=styles["BodyText"], fontSize=9,
                                  leftIndent=0.2 * inch, spaceAfter=8, leading=11)

    def table(self, rows, widths, header=True, label_column=False, align="LEFT"):
        """A table in the house style: blue header row, banded body, grey grid."""
        table = Table(rows, colWidths=[w * inch for w in widths], repeatRows=1 if header else 0)
        style = [
            ("ALIGN", (0, 0), (-1, -1), align),
            ("ALIGN", (0, 0), (0, -1), "LEFT"),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("FONTSIZE", (0, 0), (-1, -1), 8.5),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ("TOPPADDING", (0, 0), (-1, -1), 4),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
        ]
        if header:
            style += [("BACKGROUND", (0, 0), (-1, 0), BLUE),
                      ("TEXTCOLOR", (0, 0), (-1, 0), colors.whitesmoke),
                      ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                      ("ROWBACKGROUNDS", (0, 1), (-1, -1),
                       [colors.white, colors.HexColor("#F5F5F5")])]
        if label_column:
            first = 1 if header else 0
            style += [("FONTNAME", (0, first), (0, -1), "Helvetica-Bold")]
            if not header:
                style += [("BACKGROUND", (0, 0), (0, -1), colors.HexColor("#E8EEF7"))]
        table.setStyle(TableStyle(style))
        return table

    def p(self, text, style=None):
        """A table cell that wraps; the text is data, so markup characters are escaped."""
        return Paragraph(html.escape(str(text)), style or self.cell)

    def version(self, package):
        return self.packages.get(package, "?")

    def level_csv(self, level, *parts):
        return read_csv(os.path.join(self.results, level, *parts))

    @staticmethod
    def n_clusters(table):
        """Rows of a cluster_by_* table, less the bins' "unassigned" label."""
        if table is None:
            return "N/A"
        return str((table.iloc[:, 0].astype(str) != "unassigned").sum())

    @staticmethod
    def number(value, digits=0):
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return "N/A"
        return f"{value:,.{digits}f}"

    # ---- sections

    def project_information(self, story):
        story.append(Paragraph("Project Information", self.heading))
        rows = [["Generation Date:", datetime.now().strftime("%B %d, %Y")],
                ["Author:", self.args.author],
                ["Project Directory:", self.p(os.getcwd(), self.path)],
                ["Workflow Commit:", self.args.git_commit]]
        story.append(self.table(rows, [2.2, 4.8], header=False, label_column=True))
        story.append(Spacer(1, 0.2 * inch))

    def summary_statistics(self, story):
        story.append(Paragraph("Summary Statistics", self.heading))
        rows = [["Metric", "Value"]]
        meta = self.samples
        if meta is not None:
            rows += [
                ["Visium HD Slides", str(meta["visium_slide_id"].nunique())],
                ["Capture Areas", str(meta["cytassist_sample"].nunique())],
                ["TMA Cores", str(meta["core_id"].nunique())],
                ["Patients", str(meta["patient_id"].nunique())],
                ["Diagnoses", ", ".join(sorted(meta["diagnosis"].dropna().unique()))],
                ["Sample Types", ", ".join(sorted(meta["sample_type"].dropna().unique()))],
                ["Timepoints (weeks)", ", ".join(str(t) for t in
                                                 sorted(meta["timepoint_weeks"].unique()))],
            ]
        rows.append(["Organism", "Human (" + ", ".join(self.reference.get("genomes", [])) + ")"])
        if self.metrics is not None:
            rows += [
                ["Total Sequenced Reads", self.number(self.metrics["Number of Reads"].sum())],
                ["Segmented Cells (Space Ranger)",
                 self.number(self.metrics["Number of Cells"].sum())],
            ]
        cells = self.level_csv("cell", "qc", "bins_per_core.csv")
        if cells is not None:
            kept = cells.loc[cells["stage"] == "after_filter", "n_bins"].sum()
            rows.append(["Cells Passing QC", self.number(kept)])
        clusters = self.level_csv("cell", "qc", "cluster_by_condition.csv")
        if clusters is not None:
            rows.append(["Leiden Clusters (cells)", self.n_clusters(clusters)])
        rows.append(["Probe Set Genes", self.number(self.metrics["Number of Genes"].iloc[0])
                     if self.metrics is not None else "N/A"])
        story.append(self.table(rows, [2.7, 3.3], align="CENTER"))
        story.append(Spacer(1, 0.2 * inch))

    def reference_genome(self, story):
        story.append(Paragraph("Reference Genome and Annotation", self.heading))
        ref = self.reference
        models = sorted({r["model"] for r in self.runs.values() if r["model"]})
        rows = [["Item", "Value"],
                ["Organism", "Homo sapiens"],
                ["Genome build", ", ".join(ref.get("genomes", []))],
                ["10x reference version", ref.get("version", "?")],
                ["Genome FASTA", self.p(", ".join(ref.get("input_fasta_files", [])), self.path)],
                ["Gene annotation", self.p(", ".join(ref.get("input_gtf_files", [])), self.path)],
                ["Reference built with", f"mkref {ref.get('mkref_version', '?')}"],
                ["Transcriptome", self.p(self.args.transcriptome, self.path)],
                ["Probe set", f"{self.probes.get('panel_name', '?')} "
                              f"({self.probes.get('reference_genome', '?')} "
                              f"{self.probes.get('reference_version', '')})"],
                ["Probe set file", self.p(self.args.probe_set, self.path)],
                ["Cell type annotation", ", ".join(models) or "N/A"],
                ["Ligand-receptor resource",
                 f"{self.sq.get('ccc', {}).get('resource', 'cellchatdb')} "
                 f"(LIANA+ {self.version('liana')}, OmniPath {self.version('omnipath')})"]]
        story.append(self.table(rows, [1.8, 5.2], label_column=True))
        story.append(Spacer(1, 0.2 * inch))

    def software_versions(self, story):
        story.append(Paragraph("Software Versions", self.heading))
        pipelines = sorted({r["pipelines"] for r in self.runs.values()})
        rows = [["Software", "Version", "Used for"],
                ["Space Ranger", ", ".join(pipelines), "Demultiplexing, probe alignment, "
                                                      "binning, cell segmentation"],
                ["Snakemake", self.args.snakemake_version, "Workflow management"],
                ["Python", self.version("python"), "Downstream analysis (pixi env "
                                                    "envs/squidpy)"]]
        uses = [("anndata", "Data structures"), ("scanpy", "Normalisation, PCA, kNN, UMAP"),
                ("squidpy", "Spatial graphs and statistics, image features"),
                ("harmonypy", "Batch integration"), ("leidenalg", "Clustering"),
                ("python-igraph", "Graph backend"), ("umap-learn", "UMAP embedding"),
                ("liana", "Cell-cell communication"), ("omnipath", "Ligand-receptor resource"),
                ("pydeseq2", "Pseudobulk differential expression"),
                ("scikit-learn", "PCA, neighbours"), ("numpy", "Numerics"),
                ("scipy", "Sparse matrices, statistics"), ("pandas", "Tables"),
                ("matplotlib", "Figures")]
        rows += [[name, self.version(name), use] for name, use in uses
                 if name in self.packages]
        story.append(self.table(rows, [1.6, 1.2, 4.2]))
        story.append(Spacer(1, 0.2 * inch))

    def pipeline_overview(self, story):
        story.append(Paragraph("Pipeline Overview", self.heading))
        sq, qc, cl = self.sq, self.sq.get("qc", {}), self.sq.get("cluster", {})
        sp, co = sq.get("spatial", {}), sq.get("cooccurrence", {})
        svg, ccc, img = sq.get("svg", {}), sq.get("ccc", {}), sq.get("image_features", {})
        dom, de = sq.get("domains", {}), sq.get("de", {})
        lt = cl.get("label_transfer", {})
        sr = self.cfg.get("spaceranger", {})
        pipelines = ", ".join(sorted({r["pipelines"] for r in self.runs.values()}))
        trim = qc.get("quantile_trim", {})
        sections = [
            (f"Space Ranger v{pipelines}",
             "Each capture area was processed with <i>spaceranger count</i> against the "
             f"{self.probes.get('panel_name', 'probe set')} and the "
             f"{self.reference.get('version', '')} GRCh38 transcriptome, using the CytAssist "
             "image for fiducial alignment and the H&amp;E microscope image for cell "
             "segmentation. Space Ranger aggregates the 2 µm barcoded squares into 8 and "
             "16 µm bins, assigns squares to segmented cells, and annotates cells with the "
             f"{', '.join(sorted({r['model'] for r in self.runs.values() if r['model']}))} "
             f"model (create-bam={str(sr.get('create_bam', False)).lower()}, "
             f"{sr.get('localcores', 24)} cores)."),
            ("Core assignment (scanpy " + self.version("scanpy") + ")",
             "Bins and cells were assigned to TMA cores from manual Loupe Browser region "
             "exports drawn on the 8 µm bins; each cell takes the core holding most of its "
             "2 µm squares. Bins and cells outside every core were dropped, and the per-core "
             "metadata (patient, diagnosis, sample type, timepoint) attached from the sample "
             "sheet. Per-bin QC metrics were computed with "
             "<i>sc.pp.calculate_qc_metrics</i>."),
            ("Quality control",
             "Capture areas were concatenated per analysis level (segmented cells, "
             f"{', '.join(LEVEL_LABELS[l] for l in self.levels if l != 'cell')}). Bins and "
             f"cells with fewer than {qc.get('min_counts', 10)} UMIs or more than "
             f"{qc.get('max_pct_mt', 25)}% mitochondrial UMIs were removed, then the lowest "
             f"{trim.get('lower', 0.1):.0%} of {trim.get('metric', 'total_counts')} was trimmed "
             f"within each TMA core (no upper trim). Genes detected in fewer than "
             f"{qc.get('min_cells', 10)} bins were removed. Cores whose interior median depth "
             f"fell below {qc.get('section_gate', {}).get('interior_median_counts', 20)} UMIs "
             "were flagged in the QC report but not excluded."),
            (f"Integration and clustering: scanpy {self.version('scanpy')}, "
             f"harmonypy {self.version('harmonypy')}, leidenalg {self.version('leidenalg')}",
             "Segmented cells were normalised to the median library size and "
             f"log-transformed; {cl.get('n_top_genes', 2000)} highly variable genes were "
             f"selected, scaled and reduced to {cl.get('n_pcs', 30)} principal components. "
             "Harmony (theta "
             f"{cl.get('harmony_theta', 2.0)}) corrected the embedding separately for each "
             f"batch key ({', '.join(cl.get('batch_keys', []))}; primary: "
             f"{cl.get('primary_batch_key', '')}). A {cl.get('n_neighbors', 30)}-nearest-"
             "neighbour graph on the corrected components was clustered with Leiden "
             f"(resolution {cl.get('resolution', 0.5)}) and embedded with UMAP (min_dist "
             f"{cl.get('umap_min_dist', 0.1)}). Cluster markers were ranked with "
             "<i>sc.tl.rank_genes_groups</i>."),
            ("Bin labels and tissue domains",
             "Bins are too shallow to integrate directly, so each bin took the Leiden "
             "cluster and Space Ranger cell types of the segmented cells covering most of "
             "its area (Space Ranger barcode mappings); bins whose cell coverage was below "
             f"{lt.get('min_cell_fraction', 0.25):.0%} were left unassigned. "
             + (f"Tissue domains were found for {', '.join(LEVEL_LABELS[l] for l in dom.get('levels', []))} "
                "with a BANKSY-style neighbourhood augmentation: each bin's scaled expression "
                f"was joined to the mean of its {dom.get('n_rings', 2)}-ring spatial "
                f"neighbourhood (lambda {dom.get('lambda', 0.8)}), followed by PCA "
                f"({dom.get('n_pcs', 30)} PCs), Harmony and Leiden (resolution "
                f"{dom.get('resolution', 0.3)})." if dom.get("enabled") else "")),
            (f"Spatial statistics: squidpy {self.version('squidpy')}",
             "Spatial neighbour graphs were built per capture area (bins: square grid, "
             f"{sp.get('n_neighs', 8)} neighbours; cells: {sq.get('cell', {}).get('n_neighs', 6)} "
             "nearest neighbours). Neighbourhood enrichment ("
             f"{sp.get('n_perms', 1000)} permutations), interaction matrices and centrality "
             "scores were computed on the Leiden clusters, and co-occurrence (up to "
             f"{co.get('max_radius_um', 500)} µm in {co.get('interval', 50)} steps) within "
             f"each {' and each '.join(sp.get('group_keys', []))}."),
            ("Spatially variable genes",
             f"Moran's I was computed per capture area with <i>sq.gr.spatial_autocorr</i> "
             f"({svg.get('n_perms', 100)} permutations); the top {svg.get('n_top_genes', 20)} "
             "genes per capture area are plotted."),
            (f"Cell-cell communication: LIANA+ {self.version('liana')}",
             "Ligand-receptor interactions between Leiden clusters were scored per condition "
             f"with LIANA+'s CellChat method against {ccc.get('resource', 'cellchatdb')} "
             f"({ccc.get('n_perms', 1000)} permutations, minimum expression proportion "
             f"{ccc.get('expr_prop', 0.1)}). "
             + ("Spatially local ligand-receptor scores were computed per bin with "
                f"<i>li.mt.bivariate</i> ({ccc.get('local_name', 'cosine')}, "
                f"{ccc.get('local_n_perms', 100)} permutations)." if ccc.get("spatial_local")
                else "")),
            ("H&amp;E image features",
             f"{' and '.join(img.get('features', [])).capitalize()} features were extracted from the "
             "microscope image under each bin and cell with "
             f"<i>sq.im.calculate_image_features</i> (image scale {img.get('scale', 0.5)}, "
             f"spot scale {img.get('spot_scale', 2.0)})."),
            (f"Differential expression: PyDESeq2 {self.version('pydeseq2')}",
             "Raw counts were summed into pseudobulk profiles per TMA core (technical "
             "replicate cores summed together), both over the whole core and per Leiden "
             f"cluster (minimum {de.get('min_obs_per_pseudobulk', 50)} bins or cells per "
             f"pseudobulk). Genes with fewer than {de.get('min_gene_counts', 10)} counts were "
             "removed. Paired comparisons used the design ~ patient_id + factor, dropping "
             "patients seen at only one level of the factor; unpaired "
             f"comparisons ~ factor. Genes were called significant at padj &lt; "
             f"{de.get('alpha', 0.05)} (Benjamini-Hochberg). Levels: "
             f"{', '.join(LEVEL_LABELS[l] for l in de.get('levels', []))}."),
        ]
        for title, text in sections:
            story.append(Paragraph(f"<b>{title}</b>", self.h3))
            story.append(Paragraph(text, self.body))

    def sequencing_statistics(self, story):
        story.append(Paragraph("Sequencing and Mapping Statistics", self.heading))
        if self.metrics is None:
            story.append(Paragraph("results/spaceranger_metrics_summary.csv was not found.",
                                   self.body))
            return
        rows = [["Capture\nArea", "Reads", "Valid\nBarcodes", "Seq.\nSaturation",
                 "Mapped to\nProbe Set", "Bins 16 µm\nUnder Tissue", "Cells",
                 "Median\nGenes/Cell", "Median\nUMIs/Cell"]]
        for _, m in self.metrics.sort_values("run_id").iterrows():
            rows.append([m["run_id"], self.number(m["Number of Reads"]),
                         f"{m['Valid Barcodes']:.1%}", f"{m['Sequencing Saturation']:.1%}",
                         f"{m['Reads Mapped Confidently to Probe Set']:.1%}",
                         self.number(m["Number of Bins Under Tissue 16 µm"]),
                         self.number(m["Number of Cells"]),
                         self.number(m["Median Genes per Cell"]),
                         self.number(m["Median UMIs per Cell"])])
        story.append(self.table(rows, [0.8, 0.9, 0.7, 0.7, 0.75, 0.85, 0.7, 0.7, 0.7],
                                align="RIGHT"))
        story.append(Paragraph("Full per-run metrics are in "
                               "results/spaceranger_metrics_summary.csv and each run's "
                               "outs/web_summary.html.", self.body))

    def qc_statistics(self, story):
        story.append(Paragraph("Filtering and Clustering Summary", self.heading))
        rows = [["Level", "Before QC", "After QC", "% Kept", "Median UMIs\n(after)",
                 "Clusters", "Bins Labelled\nfrom Cells"]]
        for level in self.levels:
            per_core = self.level_csv(level, "qc", "bins_per_core.csv")
            if per_core is None:
                continue
            before = per_core.loc[per_core["stage"] == "before_filter", "n_bins"].sum()
            after_rows = per_core[per_core["stage"] == "after_filter"]
            after = after_rows["n_bins"].sum()
            clusters = self.level_csv(level, "qc", "cluster_by_condition.csv")
            transfer = self.level_csv(level, "qc", "label_transfer.csv")
            labelled = "N/A (clustered)"
            if transfer is not None:
                labelled = (f"{(transfer['fraction_assigned'] * transfer['n_bins']).sum() / transfer['n_bins'].sum():.1%}")
            rows.append([LEVEL_LABELS.get(level, level), self.number(before),
                         self.number(after), f"{after / before:.1%}" if before else "N/A",
                         self.number(after_rows["median_counts"].median()),
                         self.n_clusters(clusters), labelled])
        story.append(self.table(rows, [1.3, 0.95, 0.95, 0.7, 0.9, 0.7, 1.1], align="RIGHT"))
        story.append(Paragraph(
            "Median UMIs is the median over cores of each core's median. Bin-level clusters "
            "are the transferred cell clusters; the remaining bins are covered by the "
            "tissue domains. Per-core tables are in results/squidpy/&lt;level&gt;/qc/.",
            self.body))

    def de_comparisons(self, story):
        story.append(Paragraph("Differential Expression Comparisons Undertaken", self.heading))
        de = self.sq.get("de", {})
        levels = de.get("levels", [])
        alpha = de.get("alpha", 0.05)
        rows = [["Comparison", "Test vs Reference", "Subset", "Paired", "Cores\n(test/ref)"]
                + [f"Sig. Genes\n{SHORT_LABELS.get(l, l)}" for l in levels]]
        for name, spec in de.get("comparisons", {}).items():
            subset = ", ".join(f"{k}={v}" for k, v in (spec.get("subset") or {}).items())
            cores, sig = "N/A", []
            for level in levels:
                result = self.level_csv(level, "de", name, "whole_core.csv")
                sig.append(self.number((result["padj"] < alpha).sum())
                           if result is not None else "N/A")
                used = self.level_csv(level, "de", name, "samples.csv")
                if used is not None and cores == "N/A":
                    used = used[(used["leiden"] == "whole_core") & used["used"]]
                    counts = used[spec["factor"]].astype(str).value_counts()
                    cores = f"{counts.get(str(spec['test']), 0)}/{counts.get(str(spec['reference']), 0)}"
            rows.append([self.p(name), self.p(f"{spec['factor']}: {spec['test']} vs "
                                              f"{spec['reference']}"),
                         self.p(subset), "yes" if spec.get("paired") else "no", cores] + sig)
        story.append(self.table(rows, [1.6, 1.4, 1.3, 0.5, 0.6] + [1.6 / max(len(levels), 1)] * len(levels),
                                align="CENTER"))
        story.append(Paragraph(
            f"Significant genes: whole-core pseudobulk, padj &lt; {alpha}. Per-cluster "
            "results and volcano plots are in results/squidpy/&lt;level&gt;/de/&lt;comparison&gt;/.",
            self.body))

    def workflow_steps(self, story):
        story.append(Paragraph("Workflow Steps Executed", self.heading))
        story.append(Paragraph(
            "Snakemake rules that produced the current outputs, from Snakemake's own job "
            "records (.snakemake/metadata), in workflow order. Exact commands are in "
            "output/logs/ and logs/&lt;rule&gt;/.", self.body))
        order = sorted(self.jobs, key=lambda r: self.descriptions.get(r, (10 ** 6, ""))[0])
        rows = [["Rule", "Description", "Jobs", "Last Finished"]]
        for rule in order:
            ends = [j["endtime"] for j in self.jobs[rule] if j["endtime"]]
            rows.append([self.p(rule), self.p(self.descriptions.get(rule, (0, ""))[1]),
                         str(len(self.jobs[rule])),
                         datetime.fromtimestamp(max(ends)).strftime("%Y-%m-%d")
                         if ends else "N/A"])
        story.append(self.table(rows, [1.7, 3.9, 0.5, 0.9]))

    def sample_details(self, story):
        story.append(NextPageTemplate("Landscape"))
        story.append(PageBreak())
        story.append(Paragraph("Sample Details", self.heading))
        if self.samples is None:
            story.append(Paragraph("No sample sheet found.", self.body))
            return
        kept = {}
        per_core = self.level_csv("cell", "qc", "bins_per_core.csv")
        if per_core is not None:
            after = per_core[per_core["stage"] == "after_filter"]
            kept = dict(zip(after["core_id"], after["n_bins"]))
        rows = [["Core", "Sample Name", "Capture\nArea", "Slide", "Patient", "Age", "Sex",
                 "Skin\nType", "Diagnosis", "Sample\nType", "Week", "Cells\n(after QC)"]]
        for _, s in self.samples.sort_values(["cytassist_sample", "core_index"]).iterrows():
            rows.append([self.p(s["core_id"]), self.p(s["sample_name"]),
                         s["cytassist_sample"], self.p(s["visium_slide_id"]),
                         s["patient_id"], str(s["age"]), s["sex"], s["skin_type"],
                         s["diagnosis"], s["sample_type"], str(s["timepoint_weeks"]),
                         self.number(kept.get(s["core_id"]))])
        story.append(self.table(rows, [1.2, 1.1, 0.75, 1.0, 0.6, 0.4, 0.35, 0.45, 0.75,
                                       0.75, 0.45, 0.75], align="CENTER"))
        story.append(Spacer(1, 0.2 * inch))
        story.append(NextPageTemplate("Portrait"))

    def references(self, story):
        story.append(PageBreak())
        story.append(Paragraph("References", self.heading))
        for number, reference in enumerate(REFERENCES, 1):
            story.append(Paragraph(f"<b>{number}.</b> {reference}", self.ref))

    def generate(self):
        margin = 0.75 * inch
        pw, ph = letter
        lw, lh = landscape(letter)
        portrait = PageTemplate(
            id="Portrait", frames=[Frame(margin, margin, pw - 2 * margin, ph - 2 * margin)],
            onPage=lambda canvas, doc: canvas.setPageSize(letter))
        wide = PageTemplate(
            id="Landscape", frames=[Frame(margin, margin, lw - 2 * margin, lh - 2 * margin)],
            onPage=lambda canvas, doc: canvas.setPageSize(landscape(letter)))
        doc = BaseDocTemplate(self.args.output, pagesize=letter,
                              pageTemplates=[portrait, wide], title=self.args.title,
                              author=self.args.author)
        story = [Paragraph(self.args.title, self.title), Spacer(1, 0.2 * inch)]
        self.project_information(story)
        self.summary_statistics(story)
        self.reference_genome(story)
        self.software_versions(story)
        self.pipeline_overview(story)
        self.sequencing_statistics(story)
        self.qc_statistics(story)
        self.de_comparisons(story)
        self.workflow_steps(story)
        self.sample_details(story)
        self.references(story)
        doc.build(story)
        print(f"Report generated: {self.args.output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--metadata-dir", default=".snakemake/metadata")
    parser.add_argument("--workflow-files", nargs="+", required=True,
                        help="Snakefile and included .smk files, in include order")
    parser.add_argument("--spaceranger-runs", nargs="+", required=True,
                        help="output/spaceranger/<run id> directories")
    parser.add_argument("--spaceranger-metrics", required=True,
                        help="combined Space Ranger metrics CSV")
    parser.add_argument("--transcriptome", required=True)
    parser.add_argument("--probe-set", required=True)
    parser.add_argument("--pixi-list", required=True,
                        help="`pixi list --json` of the squidpy env")
    parser.add_argument("--snakemake-version", required=True)
    parser.add_argument("--git-commit", required=True)
    parser.add_argument("--author", default="Kevin Stachelek")
    parser.add_argument("--title", default="Visium HD Spatial Transcriptomics Project Report")
    parser.add_argument("--output", required=True, help="PDF to write")
    ReportGenerator(parser.parse_args()).generate()


if __name__ == "__main__":
    main()
