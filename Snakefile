# Multi-sample CellRanger Snakemake workflow
# For processing multiple single-cell RNA-seq samples

EMAIL = "kstachel@uci.edu"
onstart:
   shell("mail -s 'STARTED' {EMAIL} < {log}")
onsuccess:
   shell("mail -s 'DONE' {EMAIL} < {log}")
onerror:
   shell("mail -s 'ERROR' {EMAIL} < {log}")

from collections import Counter
import csv
from datetime import datetime
import fnmatch
import glob
import json
import os
import re

# Load configuration
configfile: "config.yaml"

# Auto-detect samples from FASTQ directory
def get_samples_from_fastq_dir():
    """
    Detect sample names by inspecting the FASTQ directory structure.

    Supported layouts:
    - per-sample subfolders: <FASTQ_DIR>/<sample>/<sample>_S*_L*_R1_*.fastq.gz
    - flat files in FASTQ_DIR: <FASTQ_DIR>/<sample>_S*_L*_R1_*.fastq.gz

    Returns a sorted list of unique sample names, excluding 'Undetermined'.
    """
    # Prefer configured path; fall back to default
    fastq_path = config.get("paths", {}).get("fastqs", "data/FASTQ")
    if not os.path.exists(fastq_path):
        return []

    samples = set()

    # 1) Look inside one-level subdirectories (typical 10x demux output)
    subdirs = [d for d in glob.glob(os.path.join(fastq_path, "*")) if os.path.isdir(d)]
    for d in subdirs:
        sample_dir_name = os.path.basename(d)
        if sample_dir_name == "Undetermined":
            continue
        # R1 files inside the subdir; be lenient on prefix
        r1_in_dir = glob.glob(os.path.join(d, "*_R1_*.fastq.gz"))
        if r1_in_dir:
            samples.add(sample_dir_name)

    # 2) Also consider flat layout at root (no subdirs)
    r1_files_root = glob.glob(os.path.join(fastq_path, "*_R1_*.fastq.gz"))
    for file in r1_files_root:
        basename = os.path.basename(file)
        m = re.match(r'^(.+?)_S\d+_', basename)
        if m:
            name = m.group(1)
            if name != "Undetermined":
                samples.add(name)

    return sorted(samples)

# Extract sample list and configuration
# Use auto-detected samples if available, otherwise fall back to config
auto_samples = get_samples_from_fastq_dir()
SAMPLES = auto_samples if auto_samples else (config.get("samples") or [])
FASTQ_DIR = config.get("paths", {}).get("fastqs", "data/FASTQ")
TRANSCRIPTOME = config["references"]["transcriptome"]
OUTPUT_DIR = config["paths"]["output"]
RESULTS_DIR = config["paths"].get("results", "results")

# Visium HD / Space Ranger configuration
SPACERANGER = config.get("spaceranger", {})

# Downstream spatial analysis (squidpy / LIANA+) configuration
SQUIDPY = config.get("squidpy", {})
SQUIDPY_DIR = f"{OUTPUT_DIR}/squidpy"
SQUIDPY_RESULTS = f"{RESULTS_DIR}/squidpy"
# pixi env for every squidpy rule (envs/squidpy/pixi.toml, locked in pixi.lock)
SQUIDPY_MANIFEST = SQUIDPY.get("pixi_manifest", "envs/squidpy/pixi.toml")
# --frozen: use the installed env as locked, without re-checking or updating it, so
# parallel jobs do not race on envs/squidpy/.pixi (install it first with
# scripts/install_squidpy_pixi.sh).
# -X faulthandler: a segfault in numba/igraph/BLAS still prints the Python stack to the log
SQUIDPY_PYTHON = (f"pixi run --frozen --manifest-path {SQUIDPY_MANIFEST} "
                  "python -X faulthandler")
BIN_SIZES = [str(b) for b in SQUIDPY.get("bin_sizes", ["008", "016"])]
HEAVY_BIN = str(SQUIDPY.get("heavy_bin_size", "016"))
SQUIDPY_CELL = SQUIDPY.get("cell", {})
USE_CELLS = bool(SQUIDPY_CELL.get("enabled", True))
# analysis levels: a directory name per bin size, plus "cell" for the segmented outputs
LEVELS = [f"square_{b}um" for b in BIN_SIZES] + (["cell"] if USE_CELLS else [])
HEAVY_LEVELS = [f"square_{HEAVY_BIN}um"] + (["cell"] if USE_CELLS else [])
# graphs, SVG and plots are drawn per capture area; kept apart from the Harmony batch key
LIBRARY_KEY = SQUIDPY.get("library_key", "run_id")
CONDITIONS = SQUIDPY.get("conditions", ["normal", "pigmented"])
CORE_ANNOTATION_DIR = SQUIDPY.get("core_annotation_dir", "metadata/core_annotations")
CORE_ANNOTATION_PATTERN = SQUIDPY.get("core_annotation_pattern", "{sample}_cores.csv")
SAMPLE_METADATA = SQUIDPY.get("sample_metadata", "metadata/sample_metadata.csv")

# {level} sits next to {sample} in the squidpy paths, so it has to be constrained or
# Snakemake splits "square_016um/per_sample/A1_TMA4" in the wrong place
wildcard_constraints:
    level = r"square_\d{3}um|cell",
    condition = r"[A-Za-z0-9_-]+",
    comparison = r"[A-Za-z0-9_-]+"


def squidpy_resource(rule_name, key, default):
    """Look up one Slurm resource for a squidpy rule from config.yaml."""
    resources = SQUIDPY.get("resources", {})
    return resources.get(rule_name, {}).get(key, default)


def squidpy_partition():
    return SQUIDPY.get("resources", {}).get("partition", "standard")


def squidpy_account():
    return SQUIDPY.get("resources", {}).get("account", "sbsandme_lab")

def read_cytassist_run_info(path):
    """Parse a CytAssist *-run-info.csv into a key -> value dict."""
    info = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split(",")
            if len(fields) >= 2:
                info[fields[0].strip()] = fields[1].strip()
    return info

def sanitize_run_id(text):
    """Space Ranger --id only accepts [a-zA-Z0-9_-]+."""
    return re.sub(r"[^A-Za-z0-9_-]+", "_", text)

def get_sample_override(sample_id, name, slide, area):
    """
    Per-sample config overrides.

    Keyed by Space Ranger id, by CytAssist sample name, or by <slide>_<area>;
    the last is the only stable key when a sample name is reused across slides.
    """
    overrides = SPACERANGER.get("samples") or {}
    for key in (sample_id, name, f"{slide}_{area}"):
        if key in overrides:
            return overrides[key] or {}
    return {}

def find_fastq_dirs(prefix):
    """
    All directories holding R1 FASTQs for a FASTQ sample prefix, sorted.

    A sample sequenced across several runs (e.g. a top-up) has one directory
    per run; Space Ranger takes them as a comma-separated --fastqs list.
    """
    root = SPACERANGER.get("fastq_root", "data")
    hits = glob.glob(os.path.join(root, "**", f"{prefix}_S*_L*_R1_*.fastq.gz"),
                     recursive=True)
    return sorted({os.path.dirname(path) for path in hits})

def find_cytassist_image(run_dir, name, slide, area):
    """
    CytAssist image for a capture area, or None.

    The instrument writes <serial>_<datetime>_<slide>_<area>_<sample>.tif, but
    the sample name entered at the instrument sometimes already carries the
    area, so try the most specific patterns first. Anchoring on the slide keeps
    a sample name that is the suffix of another from matching the wrong file.
    """
    patterns = [
        f"*_{slide}_{area}_{name}.tif",
        f"*_{slide}_{name}.tif",
        f"*_{area}_{name}.tif",
        f"*_{name}.tif",
    ]
    for pattern in patterns:
        hits = sorted(glob.glob(os.path.join(run_dir, pattern)))
        if len(hits) > 1:
            print(f"Warning: {len(hits)} CytAssist images match {pattern} in "
                  f"{run_dir}; using {os.path.basename(hits[0])}")
        if hits:
            return hits[0]
    return None

def find_microscope_image(name, slide, area):
    """
    Microscope image for a capture area, or None.

    spaceranger.microscope_image_pattern is a glob relative to
    microscope_image_dir, in which {sample}, {slide} and {area} are
    substituted, along with {token}: the first group captured from the sample
    name by spaceranger.microscope_image_regex. In this project the capture
    areas are named <area>_TMA<n> while the images are one per TMA, so the
    regex pulls out <n> and the pattern is 'TMA {token}.tif'.
    """
    image_dir = SPACERANGER.get("microscope_image_dir")
    if not image_dir:
        return None

    token = ""
    token_re = SPACERANGER.get("microscope_image_regex")
    if token_re:
        m = re.search(token_re, name)
        if not m:
            return None
        token = m.group(1) if m.groups() else m.group(0)

    pattern = SPACERANGER.get("microscope_image_pattern", "*{sample}*.tif")
    pattern = pattern.format(sample=name, slide=slide, area=area, token=token)

    hits = sorted(glob.glob(os.path.join(image_dir, pattern)))
    if not hits:
        # globs are case sensitive, image extensions often are not
        hits = sorted(os.path.join(image_dir, entry)
                      for entry in os.listdir(image_dir)
                      if fnmatch.fnmatch(entry.lower(), pattern.lower()))
    return hits[0] if hits else None

def collect_cytassist_areas():
    """
    Every capture area described by the CytAssist run folders.

    Each folder holds a *-run-info.csv naming the slide and the sample loaded
    in each capture area, plus one tif per area. Returns a list of
    {slide, area, name, cytaimage} in slide order.
    """
    cytassist_dir = SPACERANGER.get("cytassist_dir")
    if not cytassist_dir or not os.path.isdir(cytassist_dir):
        return []

    areas = []
    for run_dir in sorted(glob.glob(os.path.join(cytassist_dir, "*"))):
        run_info_files = sorted(glob.glob(os.path.join(run_dir, "*-run-info.csv")))
        if not run_info_files:
            continue
        info = read_cytassist_run_info(run_info_files[0])
        slide = info.get("Visium Slide ID")

        for key, name in sorted(info.items()):
            m = re.match(r"^Sample ([A-Z]+\d+)$", key)
            if not m or not name:
                continue
            area = m.group(1)

            cytaimage = find_cytassist_image(run_dir, name, slide, area)
            if cytaimage is None:
                print(f"Warning: no CytAssist image for {name} ({slide} {area}) "
                      f"in {run_dir}; skipping")
                continue

            areas.append({"slide": slide, "area": area, "name": name,
                          "cytaimage": cytaimage})
    return areas

def get_visiumhd_samples():
    """
    Space Ranger inputs per capture area, keyed by run id.

    A capture area is identified by its slide and area, not by the sample name
    the operator typed, which is free text and may repeat across slides. The
    run id is that name when it is unique and valid, and <name>_<slide>_<area>
    otherwise. Capture areas without FASTQs on disk are skipped.
    """
    areas = collect_cytassist_areas()
    name_counts = Counter(entry["name"] for entry in areas)

    samples = {}
    missing_images = []
    for entry in areas:
        name, slide, area = entry["name"], entry["slide"], entry["area"]

        if name_counts[name] > 1:
            sample_id = sanitize_run_id(f"{name}_{slide}_{area}")
            print(f"Warning: sample name '{name}' is used by {name_counts[name]} "
                  f"capture areas; this one runs as '{sample_id}'")
        else:
            sample_id = sanitize_run_id(name)
            if sample_id != name:
                print(f"Note: sample name '{name}' runs as '{sample_id}' "
                      f"(Space Ranger --id allows [a-zA-Z0-9_-] only)")

        if sample_id in samples:
            raise ValueError(f"Duplicate Space Ranger run id '{sample_id}'; set a "
                             f"distinct id under spaceranger.samples")

        override = get_sample_override(sample_id, name, slide, area)

        fastq_prefix = override.get("fastq_prefix", name)
        fastq_dirs = find_fastq_dirs(fastq_prefix)
        if not fastq_dirs:
            print(f"Warning: no FASTQs under "
                  f"{SPACERANGER.get('fastq_root', 'data')} for '{fastq_prefix}' "
                  f"({slide} {area}); skipping. Set spaceranger.samples."
                  f"{slide}_{area}.fastq_prefix if the core named them differently")
            continue

        # the H&E image is required: Space Ranger segments nuclei from it
        image = override.get("image") or find_microscope_image(name, slide, area)
        if not image:
            missing_images.append(f"{sample_id} ({slide} {area})")

        samples[sample_id] = {
            # slide/area are read from the CytAssist image metadata; these are
            # only passed to Space Ranger when set in the config
            "slide": override.get("slide"),
            "area": override.get("area"),
            "run_slide": slide,
            "run_area": area,
            "name": name,
            "cytaimage": override.get("cytaimage", entry["cytaimage"]),
            "image": image,
            "loupe_alignment": override.get("loupe_alignment"),
            "fastq_prefix": fastq_prefix,
            "fastq_dirs": fastq_dirs,
        }
    if missing_images:
        raise ValueError(
            "No H&E microscope image for: " + ", ".join(missing_images) +
            ". Cell segmentation needs one per capture area; check "
            "spaceranger.microscope_image_* or set spaceranger.samples.<id>.image")
    return samples

VISIUMHD_SAMPLES = get_visiumhd_samples()

# Print detected samples for debugging
print(f"Visium HD samples: {sorted(VISIUMHD_SAMPLES)}")
print(f"Auto-detected samples: {auto_samples}")
print(f"Using samples: {SAMPLES}")
print(f"FASTQ directory: {FASTQ_DIR}")
print(f"Number of samples to process: {len(SAMPLES)}")

# Rule all - defines final outputs for all samples
rule all:
    input:
        # Visium HD outputs
        expand(f"{OUTPUT_DIR}/spaceranger/{{sample}}/outs/web_summary.html", sample=sorted(VISIUMHD_SAMPLES)),
        expand(f"{OUTPUT_DIR}/spaceranger/{{sample}}/outs/binned_outputs/square_008um/filtered_feature_bc_matrix", sample=sorted(VISIUMHD_SAMPLES)),
        expand(f"{OUTPUT_DIR}/spaceranger/{{sample}}/outs/segmented_outputs/filtered_feature_cell_matrix.h5", sample=sorted(VISIUMHD_SAMPLES)),
        expand(f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/web_summary.html", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/filtered_feature_bc_matrix", sample=SAMPLES),
        f"{OUTPUT_DIR}/multi_sample_summary.csv",
        f"{OUTPUT_DIR}/scanpy/combined_harmony_integrated.h5ad",
        f"{OUTPUT_DIR}/scanpy/inspect_integrated_anndata_combined.ipynb",
        # per-sample filtering timeline plots
        expand(f"{OUTPUT_DIR}/qc/filtering_timeline/{{sample}}.png", sample=SAMPLES)
        # scenic outputs
        # expand(f"{OUTPUT_DIR}/scenic/{{sample}}.loom", sample=SAMPLES),
        # # CellBender outputs
        # expand(f"{OUTPUT_DIR}/cellbender/{{sample}}_filtered.h5", sample=SAMPLES),
        # Scrublet outputs
        # expand(f"{OUTPUT_DIR}/scrublet/{{sample}}_doublets.csv", sample=SAMPLES)
        # # Velocyto outputs
        # expand(f"{OUTPUT_DIR}/velocyto/{{sample}}.loom", sample=SAMPLES)
        # # scVelo outputs
        # expand(f"{OUTPUT_DIR}/scanpy/{{sample}}_scvelo.h5ad", sample=SAMPLES),
        # # Pileup and phasing outputs
        # expand(f"{OUTPUT_DIR}/numbat/{{sample}}_allele_counts.tsv.gz", sample=SAMPLES),
        # # Numbat outputs
        # expand(f"{OUTPUT_DIR}/numbat/{{sample}}/done.txt", sample=SAMPLES)


# Target rule: Visium HD counting only, without the single-cell downstream steps
rule spaceranger_all:
    input:
        expand(f"{OUTPUT_DIR}/spaceranger/{{sample}}/outs/web_summary.html", sample=sorted(VISIUMHD_SAMPLES)),
        expand(f"{OUTPUT_DIR}/spaceranger/{{sample}}/outs/binned_outputs/square_008um/filtered_feature_bc_matrix", sample=sorted(VISIUMHD_SAMPLES)),
        expand(f"{OUTPUT_DIR}/spaceranger/{{sample}}/outs/segmented_outputs/filtered_feature_cell_matrix.h5", sample=sorted(VISIUMHD_SAMPLES)),
        f"{RESULTS_DIR}/spaceranger_metrics_summary.csv"


# Helper function to get FASTQ files for a sample
def get_fastq_files(wildcards, read):
    """Get R1 or R2 FASTQ file for a sample"""
    import glob
    # Support both flat and per-sample subfolder layouts
    pattern_flat = f"{FASTQ_DIR}/{wildcards.sample}_S*_L*_{read}_*.fastq.gz"
    pattern_dir = f"{FASTQ_DIR}/{wildcards.sample}/*_{read}_*.fastq.gz"
    files = glob.glob(pattern_flat) + glob.glob(pattern_dir)
    if not files:
        raise ValueError(f"No {read} FASTQ files found for sample {wildcards.sample}")
    return files[0]  # Return first match

# Rule: FastQC on raw FASTQ files
rule fastqc:
    input:
        r1 = lambda wildcards: get_fastq_files(wildcards, "R1"),
        r2 = lambda wildcards: get_fastq_files(wildcards, "R2")
    output:
        r1_html = f"fastqc/{{sample}}_r1_fastqc.html",
        r1_zip = f"fastqc/{{sample}}_r1_fastqc.zip",
        r2_html = f"fastqc/{{sample}}_r2_fastqc.html",
        r2_zip = f"fastqc/{{sample}}_r2_fastqc.zip"
    threads: 2
    resources:
        mem_mb = 4000,
        cpus = 2,
        partition = "standard",
        account = "sbsandme_lab"
    shell:
        """
        module load fastqc/0.11.9
        mkdir -p fastqc
        fastqc -o fastqc -t {threads} {input.r1} {input.r2}
        module unload fastqc/0.11.9
        """

# Rule: MultiQC report
rule multiqc:
    input:
        expand(f"fastqc/{{sample}}_r1_fastqc.html", sample=SAMPLES),
        expand(f"fastqc/{{sample}}_r2_fastqc.html", sample=SAMPLES)
    output:
        report = f"{OUTPUT_DIR}/multiqc_report.html"
    threads: 2
    resources:
        mem_mb = 4000,
        cpus = 2,
        partition = "standard",
        account = "sbsandme_lab"
    shell:
        """
        module load singularity/3.11.3
        singularity run /dfs9/ucightf-lab/kstachel/TOOLS/multiqc-1.20.sif multiqc {OUTPUT_DIR}/fastqc -o {OUTPUT_DIR}
        module unload singularity/3.11.3
        """

# Rule: CellRanger count for each sample
rule cellranger_count:
    input:
        transcriptome = TRANSCRIPTOME,
        fastq_dir = f"{FASTQ_DIR}/{{sample}}"
    output:
        web_summary = f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/web_summary.html",
        filtered_matrix_h5 = f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/filtered_feature_bc_matrix.h5",
        filtered_matrix_dir = directory(f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/filtered_feature_bc_matrix"),
        bam_file = f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/possorted_genome_bam.bam",
        molecule_info = f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/molecule_info.h5"
    params:
        output_dir = OUTPUT_DIR,
        localcores = config["params"]["localcores"],
        localmem = config["params"]["localmem"],
        chemistry = config["params"].get("chemistry", "auto"),
        sample_id = "{sample}"
    threads: 8
    resources:
        mem_mb = 47663,  # 48GB in MB
        mem_mb_per_cpu=6000,
        cpus = 8,
        partition = "standard",
        account = "sbsandme_lab"
    shell:
        """
        module load cellranger/8.0.1
        
        rm  -rf {params.sample_id}
        rm  -rf {params.output_dir}/cellranger/{params.sample_id}

        cellranger count --id={params.sample_id} \
        --transcriptome={input.transcriptome} \
        --fastqs={input.fastq_dir} \
        --create-bam=true \
        --sample={params.sample_id} \
        --localcores={params.localcores} \
        --localmem={params.localmem} \
        --chemistry={params.chemistry}

        mv {params.sample_id} {params.output_dir}/cellranger/{params.sample_id}
        
        module unload cellranger/8.0.1
        """


# Rule: Space Ranger count for each Visium HD capture area
def spaceranger_input(sample, key):
    """Optional per-sample input file, or [] when not configured."""
    value = VISIUMHD_SAMPLES[sample].get(key)
    return value if value else []

def spaceranger_extra_args(wildcards, input):
    """Assemble the optional Space Ranger flags for a sample."""
    meta = VISIUMHD_SAMPLES[wildcards.sample]
    args = []
    if input.loupe_alignment:
        args.append(f'--loupe-alignment="{input.loupe_alignment}"')
    if meta.get("slide") and meta.get("area"):
        # only set when the CytAssist image metadata is known to be wrong
        args.append(f'--slide={meta["slide"]} --area={meta["area"]} --override-id')
    bin_size = SPACERANGER.get("custom_bin_size")
    if bin_size:
        args.append(f"--custom-bin-size={bin_size}")
    return " ".join(args)

rule spaceranger_count:
    input:
        transcriptome = SPACERANGER.get("transcriptome", ""),
        probe_set = SPACERANGER.get("probe_set", ""),
        cytaimage = lambda wildcards: spaceranger_input(wildcards.sample, "cytaimage"),
        image = lambda wildcards: VISIUMHD_SAMPLES[wildcards.sample]["image"],
        loupe_alignment = lambda wildcards: spaceranger_input(wildcards.sample, "loupe_alignment")
    output:
        web_summary = f"{OUTPUT_DIR}/spaceranger/{{sample}}/outs/web_summary.html",
        metrics_summary = f"{OUTPUT_DIR}/spaceranger/{{sample}}/outs/metrics_summary.csv",
        filtered_matrix_h5 = f"{OUTPUT_DIR}/spaceranger/{{sample}}/outs/binned_outputs/square_008um/filtered_feature_bc_matrix.h5",
        filtered_matrix_dir = directory(f"{OUTPUT_DIR}/spaceranger/{{sample}}/outs/binned_outputs/square_008um/filtered_feature_bc_matrix"),
        cell_matrix_h5 = f"{OUTPUT_DIR}/spaceranger/{{sample}}/outs/segmented_outputs/filtered_feature_cell_matrix.h5",
        cell_segmentations = f"{OUTPUT_DIR}/spaceranger/{{sample}}/outs/segmented_outputs/cell_segmentations.geojson"
    params:
        manifest = SPACERANGER.get("pixi_manifest", "envs/spaceranger/pixi.toml"),
        outdir = f"{OUTPUT_DIR}/spaceranger/{{sample}}",
        fastq_dirs = lambda wildcards: ",".join(VISIUMHD_SAMPLES[wildcards.sample]["fastq_dirs"]),
        fastq_prefix = lambda wildcards: VISIUMHD_SAMPLES[wildcards.sample]["fastq_prefix"],
        create_bam = str(SPACERANGER.get("create_bam", False)).lower(),
        localcores = SPACERANGER.get("localcores", 24),
        localmem = SPACERANGER.get("localmem", 112),
        extra = spaceranger_extra_args
    threads: SPACERANGER.get("localcores", 24)
    resources:
        mem_mb = SPACERANGER.get("mem_mb", 120000),  # 117GB, fits any standard node
        cpus = SPACERANGER.get("localcores", 24),
        runtime = SPACERANGER.get("runtime", 2880),  # 48h
        partition = SPACERANGER.get("partition", "standard"),
        account = SPACERANGER.get("account", "sbsandme_lab")
    log:
        f"{OUTPUT_DIR}/logs/spaceranger_{{sample}}.log"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/logs
        rm -rf {params.outdir}

        pixi run --manifest-path {params.manifest} spaceranger count \
            --id={wildcards.sample} \
            --output-dir={params.outdir} \
            --transcriptome={input.transcriptome} \
            --probe-set={input.probe_set} \
            --fastqs={params.fastq_dirs} \
            --sample={params.fastq_prefix} \
            --cytaimage="{input.cytaimage}" \
            --image="{input.image}" \
            --create-bam={params.create_bam} \
            --localcores={params.localcores} \
            --localmem={params.localmem} \
            --disable-ui \
            {params.extra} > {log} 2>&1
        """


# Rule: combine the per-capture-area Space Ranger metrics into one table
#
# spaceranger aggr is not an option here: it rejects Visium HD runs
# ("Visium HD runs ... are not currently supported by spaceranger aggr").
rule combine_spaceranger_metrics:
    input:
        metrics = expand(f"{OUTPUT_DIR}/spaceranger/{{sample}}/outs/metrics_summary.csv", sample=sorted(VISIUMHD_SAMPLES))
    output:
        combined = f"{RESULTS_DIR}/spaceranger_metrics_summary.csv"
    params:
        script = "src/combine_spaceranger_metrics.py"
    threads: 1
    resources:
        mem_mb = 4000,
        cpus = 1,
        runtime = 30,
        partition = SPACERANGER.get("partition", "standard"),
        account = SPACERANGER.get("account", "sbsandme_lab")
    shell:
        """
        python {params.script} \
            --metrics {input.metrics} \
            --output {output.combined}
        """


# ---------------------------------------------------------------------------
# Downstream spatial analysis: squidpy + LIANA+
#
# spaceranger aggr rejects Visium HD, so the capture areas are combined here from the
# per-run binned matrices. Space Ranger also has no notion of which TMA core a bin came
# from, and every contrast in this experiment (normal vs pigmented, W0 vs W16) lives
# inside a capture area, so rule squidpy_make_adata depends on the manual Loupe core
# annotations named by squidpy.core_annotation_dir/_pattern -- see docs/core_annotation.md.
#
# Every shell block exports PYTHONNOUSERSITE=1: ~/.local/lib/python3.12/site-packages
# otherwise shadows numpy, scanpy, anndata, matplotlib and ~35 more packages of the
# pixi env. The rules run through `pixi run --manifest-path envs/squidpy/pixi.toml`,
# so snakemake needs no --use-conda for them.
# ---------------------------------------------------------------------------

def read_run_tma():
    """run id -> TMA number, from the cytassist_sample and tma columns of sample_metadata."""
    if not os.path.exists(SAMPLE_METADATA):
        return {}
    with open(SAMPLE_METADATA, newline="") as handle:
        return {row["cytassist_sample"]: row["tma"] for row in csv.DictReader(handle)}

RUN_TMA = read_run_tma()

def core_annotation_path(wildcards):
    """The Loupe core export for one capture area, per squidpy.core_annotation_pattern."""
    try:
        tma = RUN_TMA[wildcards.sample] if "{tma}" in CORE_ANNOTATION_PATTERN else None
    except KeyError:
        raise ValueError(f"{wildcards.sample} is not in the cytassist_sample column of "
                         f"{SAMPLE_METADATA}, so its core annotation cannot be named")
    name = CORE_ANNOTATION_PATTERN.format(sample=wildcards.sample, tma=tma)
    return os.path.join(CORE_ANNOTATION_DIR, name)

def spatial_param(key, default):
    """A spatial-graph setting per level: squidpy.cell overrides squidpy.spatial for cells."""
    def inner(wildcards):
        if wildcards.level == "cell" and key in SQUIDPY_CELL:
            return SQUIDPY_CELL[key]
        return SQUIDPY.get("spatial", {}).get(key, default)
    return inner


# Rule: one AnnData per capture area per bin size, with TMA cores assigned to bins
rule squidpy_make_adata:
    wildcard_constraints:
        level = r"square_\d{3}um"
    input:
        script = "src/squidpy_make_adata.py",
        annotation_module = "src/squidpy_core_annotation.py",
        matrix_h5 = f"{OUTPUT_DIR}/spaceranger/{{sample}}/outs/binned_outputs/{{level}}/filtered_feature_bc_matrix.h5",
        positions = f"{OUTPUT_DIR}/spaceranger/{{sample}}/outs/binned_outputs/{{level}}/spatial/tissue_positions.parquet",
        scalefactors = f"{OUTPUT_DIR}/spaceranger/{{sample}}/outs/binned_outputs/{{level}}/spatial/scalefactors_json.json",
        core_annotation = core_annotation_path,
        sample_metadata = SAMPLE_METADATA
    output:
        adata = f"{SQUIDPY_DIR}/{{level}}/per_sample/{{sample}}.h5ad",
        qc_figure = f"{SQUIDPY_RESULTS}/{{level}}/qc/{{sample}}_core_assignment.pdf"
    threads: squidpy_resource("make_adata", "cpus", 2)
    resources:
        mem_mb = squidpy_resource("make_adata", "mem_mb", 32000),
        cpus = squidpy_resource("make_adata", "cpus", 2),
        runtime = squidpy_resource("make_adata", "runtime", 60),
        partition = squidpy_partition(),
        account = squidpy_account()
    log:
        f"{OUTPUT_DIR}/logs/squidpy_make_adata_{{sample}}_{{level}}.log"
    shell:
        """
        export PYTHONNOUSERSITE=1
        export PYTHONPATH=src:${{PYTHONPATH:-}}
        {SQUIDPY_PYTHON} {input.script} \
            --matrix-h5 {input.matrix_h5} \
            --positions {input.positions} \
            --scalefactors {input.scalefactors} \
            --core-annotation {input.core_annotation} \
            --sample-metadata {input.sample_metadata} \
            --run-id {wildcards.sample} \
            --output {output.adata} \
            --qc-figure {output.qc_figure} > {log} 2>&1
        """


# Rule: one AnnData per capture area from the segmented cells (Space Ranger >= 4.0)
#
# Cores are drawn at 8 um, so the 8 um positions travel with the export; each cell takes
# the core that holds most of its 2 um bins in barcode_mappings.parquet.
rule squidpy_make_adata_cell:
    input:
        script = "src/squidpy_make_adata.py",
        annotation_module = "src/squidpy_core_annotation.py",
        matrix_h5 = f"{OUTPUT_DIR}/spaceranger/{{sample}}/outs/segmented_outputs/filtered_feature_cell_matrix.h5",
        segmentations = f"{OUTPUT_DIR}/spaceranger/{{sample}}/outs/segmented_outputs/cell_segmentations.geojson",
        scalefactors = f"{OUTPUT_DIR}/spaceranger/{{sample}}/outs/segmented_outputs/spatial/scalefactors_json.json",
        barcode_mappings = f"{OUTPUT_DIR}/spaceranger/{{sample}}/outs/barcode_mappings.parquet",
        positions = f"{OUTPUT_DIR}/spaceranger/{{sample}}/outs/binned_outputs/square_008um/spatial/tissue_positions.parquet",
        core_annotation = core_annotation_path,
        sample_metadata = SAMPLE_METADATA
    output:
        adata = f"{SQUIDPY_DIR}/cell/per_sample/{{sample}}.h5ad",
        qc_figure = f"{SQUIDPY_RESULTS}/cell/qc/{{sample}}_core_assignment.pdf"
    params:
        # optional: Space Ranger only writes it when the local Azimuth model ran
        cell_types = f"{OUTPUT_DIR}/spaceranger/{{sample}}/outs/segmented_outputs/cell_types/Azimuth/cell_types.csv"
    threads: squidpy_resource("make_adata", "cpus", 2)
    resources:
        mem_mb = squidpy_resource("make_adata", "mem_mb", 32000),
        cpus = squidpy_resource("make_adata", "cpus", 2),
        runtime = squidpy_resource("make_adata", "runtime", 60),
        partition = squidpy_partition(),
        account = squidpy_account()
    log:
        f"{OUTPUT_DIR}/logs/squidpy_make_adata_{{sample}}_cell.log"
    shell:
        """
        export PYTHONNOUSERSITE=1
        export PYTHONPATH=src:${{PYTHONPATH:-}}
        {SQUIDPY_PYTHON} {input.script} \
            --level cell \
            --matrix-h5 {input.matrix_h5} \
            --positions {input.positions} \
            --scalefactors {input.scalefactors} \
            --segmentations {input.segmentations} \
            --barcode-mappings {input.barcode_mappings} \
            --cell-types {params.cell_types} \
            --core-annotation {input.core_annotation} \
            --sample-metadata {input.sample_metadata} \
            --run-id {wildcards.sample} \
            --output {output.adata} \
            --qc-figure {output.qc_figure} > {log} 2>&1
        """


def squidpy_min_counts(wildcards):
    """Absolute count floor for one level, from squidpy.qc.min_counts_per_level.

    One number is not comparable across levels. An 8 um bin covers a quarter the area of
    a 16 um bin, so the same floor is far stricter there: at min_counts 10, 16 um kept
    127,421 of 185,921 bins (69%) against 350,892 of ~744,000 at 8 um (47%). Falls back
    to squidpy.qc.min_counts when a level is not listed.
    """
    qc = SQUIDPY.get("qc", {})
    return qc.get("min_counts_per_level", {}).get(
        wildcards.level, qc.get("min_counts", 10))


def squidpy_quantile_flags():
    """Per-section quantile trim flags, or nothing when it is switched off.

    Same shape as squidpy_subsample_flag: the script defaults to no trim, so the
    config is the only place the filter is turned on.
    """
    trim = SQUIDPY.get("qc", {}).get("quantile_trim", {})
    if not trim.get("enabled"):
        return ""
    return (f"--quantile-lower {trim.get('lower', 0.1)} "
            f"--quantile-upper {trim.get('upper', 0.9)} "
            f"--quantile-metric {trim.get('metric', 'total_counts')}")


# Rule: concatenate the capture areas for one bin size and apply QC filters
rule squidpy_combine:
    input:
        script = "src/squidpy_combine.py",
        annotation_module = "src/squidpy_core_annotation.py",
        adatas = expand(f"{SQUIDPY_DIR}/{{{{level}}}}/per_sample/{{sample}}.h5ad",
                        sample=sorted(VISIUMHD_SAMPLES))
    output:
        combined = f"{SQUIDPY_DIR}/{{level}}/combined.h5ad",
        qc_table = f"{SQUIDPY_RESULTS}/{{level}}/qc/bins_per_core.csv",
        qc_figure = f"{SQUIDPY_RESULTS}/{{level}}/qc/qc_violin.pdf",
        qc_sections = f"{SQUIDPY_RESULTS}/{{level}}/qc/qc_per_section.pdf",
        qc_quantiles = f"{SQUIDPY_RESULTS}/{{level}}/qc/quantile_bounds.csv",
        qc_section_table = f"{SQUIDPY_RESULTS}/{{level}}/qc/section_qc.csv",
        qc_section_depth = f"{SQUIDPY_RESULTS}/{{level}}/qc/section_depth.pdf"
    params:
        qc_group_key = SQUIDPY.get("qc", {}).get("group_key", "core_id"),
        quantile = squidpy_quantile_flags(),
        section_gate = SQUIDPY.get("qc", {}).get("section_gate", {}).get(
            "interior_median_counts", 50),
        min_counts = squidpy_min_counts,
        min_genes = SQUIDPY.get("qc", {}).get("min_genes", 0),
        min_cells = SQUIDPY.get("qc", {}).get("min_cells", 10),
        max_pct_mt = SQUIDPY.get("qc", {}).get("max_pct_mt", 25)
    threads: squidpy_resource("combine", "cpus", 4)
    resources:
        mem_mb = squidpy_resource("combine", "mem_mb", 128000),
        cpus = squidpy_resource("combine", "cpus", 4),
        runtime = squidpy_resource("combine", "runtime", 180),
        partition = squidpy_partition(),
        account = squidpy_account()
    log:
        f"{OUTPUT_DIR}/logs/squidpy_combine_{{level}}.log"
    shell:
        """
        export PYTHONNOUSERSITE=1
        {SQUIDPY_PYTHON} {input.script} \
            --adatas {input.adatas} \
            --output {output.combined} \
            --qc-table {output.qc_table} \
            --qc-figure {output.qc_figure} \
            --qc-section-figure {output.qc_sections} \
            --qc-group-key {params.qc_group_key} \
            --quantile-report {output.qc_quantiles} \
            --section-qc-table {output.qc_section_table} \
            --section-depth-figure {output.qc_section_depth} \
            --section-gate-counts {params.section_gate} \
            {params.quantile} \
            --min-counts {params.min_counts} \
            --min-genes {params.min_genes} \
            --min-cells {params.min_cells} \
            --max-pct-mt {params.max_pct_mt} > {log} 2>&1
        """


def squidpy_cluster_setting(key, default):
    """One squidpy.cluster setting, overridable per level.

    Bin size changes what these settings mean. An 8 um bin covers a quarter the area of a
    16 um bin and runs at a third the depth, so the neighbour count that leaves the 16 um
    graph in a single connected component leaves the 8 um graph in 61, which forces 160
    clusters holding under 100 bins each. cluster.<key>_per_level wins over cluster.<key>.
    """
    def inner(wildcards):
        cluster = SQUIDPY.get("cluster", {})
        return cluster.get(f"{key}_per_level", {}).get(
            wildcards.level, cluster.get(key, default))
    return inner


def squidpy_cluster_resource(key, default):
    """8 um carries 2.5M bins against 650k at 16 um, so it gets its own allocation."""
    def inner(wildcards):
        name = {"square_008um": "cluster_008", "cell": "cluster_cell"}.get(
            wildcards.level, "cluster_016")
        return squidpy_resource(name, key, default)
    return inner


def squidpy_subsample_flag(wildcards):
    """Cap the 8 um object, which otherwise densifies to ~20 GB for the PCA."""
    cap = SQUIDPY.get("cluster", {}).get("subsample_008")
    if wildcards.level == "square_008um" and cap:
        return f"--subsample {int(cap)}"
    return ""


# Rule: normalise, Harmony-integrate on squidpy.cluster.batch_key, and cluster
rule squidpy_cluster:
    input:
        script = "src/squidpy_cluster.py",
        combined = f"{SQUIDPY_DIR}/{{level}}/combined.h5ad"
    output:
        clustered = f"{SQUIDPY_DIR}/{{level}}/clustered.h5ad",
        markers = f"{SQUIDPY_RESULTS}/{{level}}/markers.csv"
    params:
        figure_dir = f"{SQUIDPY_RESULTS}/{{level}}/figures",
        n_top_genes = SQUIDPY.get("cluster", {}).get("n_top_genes", 2000),
        n_pcs = SQUIDPY.get("cluster", {}).get("n_pcs", 30),
        n_neighbors = squidpy_cluster_setting("n_neighbors", 15),
        resolution = squidpy_cluster_setting("resolution", 1.0),
        batch_key = SQUIDPY.get("cluster", {}).get("batch_key", "sample_id"),
        hvg_batch = (lambda k: f"--hvg-batch-key {k}" if k else "")(
            SQUIDPY.get("cluster", {}).get("hvg_batch_key")),
        theta = SQUIDPY.get("cluster", {}).get("harmony_theta", 2.0),
        crosstab_dir = f"{SQUIDPY_RESULTS}/{{level}}/qc",
        cluster_key = SQUIDPY.get("cluster", {}).get("cluster_key", "leiden"),
        harmony = "" if SQUIDPY.get("cluster", {}).get("use_harmony", True) else "--no-harmony",
        subsample = squidpy_subsample_flag
    threads: squidpy_cluster_resource("cpus", 8)
    resources:
        mem_mb = squidpy_cluster_resource("mem_mb", 128000),
        cpus = squidpy_cluster_resource("cpus", 8),
        runtime = squidpy_cluster_resource("runtime", 480),
        partition = squidpy_partition(),
        account = squidpy_account()
    log:
        f"{OUTPUT_DIR}/logs/squidpy_cluster_{{level}}.log"
    shell:
        """
        export PYTHONNOUSERSITE=1
        {SQUIDPY_PYTHON} {input.script} \
            --input {input.combined} \
            --output {output.clustered} \
            --figure-dir {params.figure_dir} \
            --markers-csv {output.markers} \
            --n-top-genes {params.n_top_genes} \
            --n-pcs {params.n_pcs} \
            --n-neighbors {params.n_neighbors} \
            --resolution {params.resolution} \
            --batch-key {params.batch_key} \
            {params.hvg_batch} \
            --theta {params.theta} \
            --crosstab-dir {params.crosstab_dir} \
            --cluster-key {params.cluster_key} \
            {params.harmony} {params.subsample} > {log} 2>&1
        """


# Rule: neighbourhood enrichment, interaction matrix, centrality and co-occurrence
rule squidpy_spatial_stats:
    input:
        script = "src/squidpy_spatial_stats.py",
        clustered = f"{SQUIDPY_DIR}/{{level}}/clustered.h5ad"
    output:
        stats = f"{SQUIDPY_DIR}/{{level}}/spatial_stats.h5ad",
        nhood = f"{SQUIDPY_RESULTS}/{{level}}/nhood_enrichment_zscore.csv",
        interactions = f"{SQUIDPY_RESULTS}/{{level}}/interaction_matrix.csv",
        centrality = f"{SQUIDPY_RESULTS}/{{level}}/centrality_scores.csv"
    params:
        results_dir = f"{SQUIDPY_RESULTS}/{{level}}",
        figure_dir = f"{SQUIDPY_RESULTS}/{{level}}/figures",
        cluster_key = SQUIDPY.get("cluster", {}).get("cluster_key", "leiden"),
        library_key = LIBRARY_KEY,
        coord_type = spatial_param("coord_type", "grid"),
        n_neighs = spatial_param("n_neighs", 8),
        n_rings = spatial_param("n_rings", 1),
        n_perms = SQUIDPY.get("spatial", {}).get("n_perms", 1000),
        group_keys = " ".join(SQUIDPY.get("spatial", {}).get("group_keys", ["condition"])),
        interval = SQUIDPY.get("cooccurrence", {}).get("interval", 50),
        max_radius_um = SQUIDPY.get("cooccurrence", {}).get("max_radius_um", 500),
        n_splits = SQUIDPY.get("cooccurrence", {}).get("n_splits", 20),
        max_bins = SQUIDPY.get("cooccurrence", {}).get("max_bins_per_group", 50000)
    threads: squidpy_resource("spatial_stats", "cpus", 16)
    resources:
        mem_mb = squidpy_resource("spatial_stats", "mem_mb", 192000),
        cpus = squidpy_resource("spatial_stats", "cpus", 16),
        runtime = squidpy_resource("spatial_stats", "runtime", 720),
        partition = squidpy_partition(),
        account = squidpy_account()
    log:
        f"{OUTPUT_DIR}/logs/squidpy_spatial_stats_{{level}}.log"
    shell:
        """
        export PYTHONNOUSERSITE=1
        {SQUIDPY_PYTHON} {input.script} \
            --input {input.clustered} \
            --output {output.stats} \
            --results-dir {params.results_dir} \
            --figure-dir {params.figure_dir} \
            --cluster-key {params.cluster_key} \
            --library-key {params.library_key} \
            --coord-type {params.coord_type} \
            --n-neighs {params.n_neighs} \
            --n-rings {params.n_rings} \
            --n-perms {params.n_perms} \
            --cooccurrence-interval {params.interval} \
            --cooccurrence-max-radius-um {params.max_radius_um} \
            --cooccurrence-splits {params.n_splits} \
            --max-bins-per-group {params.max_bins} \
            --group-keys {params.group_keys} \
            --n-jobs {threads} > {log} 2>&1
        """


# Rule: spatially variable genes, scored within each capture area
rule squidpy_svg:
    input:
        script = "src/squidpy_svg.py",
        clustered = f"{SQUIDPY_DIR}/{{level}}/clustered.h5ad"
    output:
        table = f"{SQUIDPY_RESULTS}/{{level}}/svg/{{sample}}_spatial_autocorr.csv"
    params:
        figure_dir = f"{SQUIDPY_RESULTS}/{{level}}/svg",
        library_key = LIBRARY_KEY,
        mode = SQUIDPY.get("svg", {}).get("mode", "moran"),
        n_perms = SQUIDPY.get("svg", {}).get("n_perms", 100),
        n_top_genes = SQUIDPY.get("svg", {}).get("n_top_genes", 20),
        n_neighs = spatial_param("n_neighs", 8),
        coord_type = spatial_param("coord_type", "grid")
    threads: squidpy_resource("svg", "cpus", 8)
    resources:
        mem_mb = squidpy_resource("svg", "mem_mb", 64000),
        cpus = squidpy_resource("svg", "cpus", 8),
        runtime = squidpy_resource("svg", "runtime", 240),
        partition = squidpy_partition(),
        account = squidpy_account()
    log:
        f"{OUTPUT_DIR}/logs/squidpy_svg_{{sample}}_{{level}}.log"
    shell:
        """
        export PYTHONNOUSERSITE=1
        {SQUIDPY_PYTHON} {input.script} \
            --input {input.clustered} \
            --run-id {wildcards.sample} \
            --output-csv {output.table} \
            --figure-dir {params.figure_dir} \
            --library-key {params.library_key} \
            --mode {params.mode} \
            --n-perms {params.n_perms} \
            --n-top-genes {params.n_top_genes} \
            --n-neighs {params.n_neighs} \
            --coord-type {params.coord_type} \
            --n-jobs {threads} > {log} 2>&1
        """


# Rule: cell-cell communication, LIANA+ with CellChat scoring against CellChatDB
#
# sq.gr.ligrec is not used: it needs omnipath, which is broken against NumPy 2 here, and
# it offers only a CellPhoneDB-style cluster permutation test.
rule squidpy_ccc:
    input:
        script = "src/squidpy_ccc.py",
        clustered = f"{SQUIDPY_DIR}/{{level}}/clustered.h5ad"
    output:
        table = f"{SQUIDPY_RESULTS}/{{level}}/ccc/{{condition}}_liana_cellchat.csv",
        figure = f"{SQUIDPY_RESULTS}/{{level}}/ccc/{{condition}}_dotplot.pdf"
    params:
        groupby = SQUIDPY.get("ccc", {}).get("groupby", "leiden"),
        group_key = SQUIDPY.get("ccc", {}).get("group_key", "condition"),
        resource = SQUIDPY.get("ccc", {}).get("resource", "cellchatdb"),
        n_perms = SQUIDPY.get("ccc", {}).get("n_perms", 1000),
        expr_prop = SQUIDPY.get("ccc", {}).get("expr_prop", 0.1),
        top_n = SQUIDPY.get("ccc", {}).get("top_n", 25),
        local_name = SQUIDPY.get("ccc", {}).get("local_name", "cosine"),
        local_n_perms = SQUIDPY.get("ccc", {}).get("local_n_perms", 100),
        library_key = LIBRARY_KEY,
        coord_type = spatial_param("coord_type", "grid"),
        n_neighs = spatial_param("n_neighs", 8),
        spatial_flag = (
            f"--spatial-output {SQUIDPY_DIR}/{{level}}/ccc/{{condition}}_local.h5ad"
            if SQUIDPY.get("ccc", {}).get("spatial_local", True) else ""
        )
    threads: squidpy_resource("ccc", "cpus", 16)
    resources:
        mem_mb = squidpy_resource("ccc", "mem_mb", 128000),
        cpus = squidpy_resource("ccc", "cpus", 16),
        runtime = squidpy_resource("ccc", "runtime", 480),
        partition = squidpy_partition(),
        account = squidpy_account()
    log:
        f"{OUTPUT_DIR}/logs/squidpy_ccc_{{condition}}_{{level}}.log"
    shell:
        """
        export PYTHONNOUSERSITE=1
        {SQUIDPY_PYTHON} {input.script} \
            --input {input.clustered} \
            --output-csv {output.table} \
            --figure {output.figure} \
            --group-key {params.group_key} \
            --condition {wildcards.condition} \
            --groupby {params.groupby} \
            --resource {params.resource} \
            --n-perms {params.n_perms} \
            --expr-prop {params.expr_prop} \
            --top-n {params.top_n} \
            --local-name {params.local_name} \
            --local-n-perms {params.local_n_perms} \
            --library-key {params.library_key} \
            --coord-type {params.coord_type} \
            --n-neighs {params.n_neighs} \
            --n-jobs {threads} {params.spatial_flag} > {log} 2>&1
        """


# Rule: H&E texture and summary features per bin, from the microscope TIFF
rule squidpy_image_features:
    input:
        script = "src/squidpy_image_features.py",
        clustered = f"{SQUIDPY_DIR}/{{level}}/clustered.h5ad"
    output:
        features = f"{SQUIDPY_DIR}/{{level}}/image_features/{{sample}}.csv"
    params:
        image_dir = SPACERANGER.get("microscope_image_dir", ""),
        # lambda, or Snakemake reads the literal {token} in the pattern as a wildcard
        image_pattern = lambda wildcards: SPACERANGER.get(
            "microscope_image_pattern", "TMA {token}.tif"),
        image_regex = SPACERANGER.get("microscope_image_regex", r"TMA\s*(\d+)$"),
        library_key = LIBRARY_KEY,
        features = " ".join(SQUIDPY.get("image_features", {}).get("features", ["summary", "texture"])),
        scale = SQUIDPY.get("image_features", {}).get("scale", 0.5),
        spot_scale = SQUIDPY.get("image_features", {}).get("spot_scale", 2.0)
    threads: squidpy_resource("image_features", "cpus", 8)
    resources:
        mem_mb = squidpy_resource("image_features", "mem_mb", 96000),
        cpus = squidpy_resource("image_features", "cpus", 8),
        runtime = squidpy_resource("image_features", "runtime", 240),
        partition = squidpy_partition(),
        account = squidpy_account()
    log:
        f"{OUTPUT_DIR}/logs/squidpy_image_features_{{sample}}_{{level}}.log"
    shell:
        """
        export PYTHONNOUSERSITE=1
        {SQUIDPY_PYTHON} {input.script} \
            --input {input.clustered} \
            --run-id {wildcards.sample} \
            --output-csv {output.features} \
            --image-dir "{params.image_dir}" \
            --image-pattern "{params.image_pattern}" \
            --image-regex "{params.image_regex}" \
            --library-key {params.library_key} \
            --features {params.features} \
            --scale {params.scale} \
            --spot-scale {params.spot_scale} \
            --n-jobs {threads} > {log} 2>&1
        """


# Rule: spatial scatter panels per capture area
rule squidpy_spatial_plots:
    input:
        script = "src/squidpy_spatial_plots.py",
        clustered = f"{SQUIDPY_DIR}/{{level}}/clustered.h5ad"
    output:
        figure = f"{SQUIDPY_RESULTS}/{{level}}/spatial/{{sample}}_panels.pdf"
    params:
        library_key = LIBRARY_KEY
    threads: squidpy_resource("spatial_plots", "cpus", 2)
    resources:
        mem_mb = squidpy_resource("spatial_plots", "mem_mb", 32000),
        cpus = squidpy_resource("spatial_plots", "cpus", 2),
        runtime = squidpy_resource("spatial_plots", "runtime", 60),
        partition = squidpy_partition(),
        account = squidpy_account()
    log:
        f"{OUTPUT_DIR}/logs/squidpy_spatial_plots_{{sample}}_{{level}}.log"
    shell:
        """
        export PYTHONNOUSERSITE=1
        {SQUIDPY_PYTHON} {input.script} \
            --input {input.clustered} \
            --run-id {wildcards.sample} \
            --output {output.figure} \
            --library-key {params.library_key} > {log} 2>&1
        """


# Pseudobulk differential expression: one comparison from squidpy.de.comparisons, whole
# core and per cluster, with pyDESeq2 across cores (see src/squidpy_de.py). Counts come
# from combined.h5ad because clustered.h5ad keeps only the highly variable genes.
SQUIDPY_DE = SQUIDPY.get("de", {})
DE_COMPARISONS = SQUIDPY_DE.get("comparisons", {})
DE_LEVELS = SQUIDPY_DE.get("levels", HEAVY_LEVELS)

rule squidpy_de:
    input:
        script = "src/squidpy_de.py",
        combined = f"{SQUIDPY_DIR}/{{level}}/combined.h5ad",
        clustered = f"{SQUIDPY_DIR}/{{level}}/clustered.h5ad",
        sample_metadata = SAMPLE_METADATA
    output:
        whole_core = f"{SQUIDPY_RESULTS}/{{level}}/de/{{comparison}}/whole_core.csv",
        per_cluster = f"{SQUIDPY_RESULTS}/{{level}}/de/{{comparison}}/per_cluster.csv",
        samples = f"{SQUIDPY_RESULTS}/{{level}}/de/{{comparison}}/samples.csv",
        volcano = f"{SQUIDPY_RESULTS}/{{level}}/de/{{comparison}}/volcano.pdf"
    params:
        outdir = f"{SQUIDPY_RESULTS}/{{level}}/de/{{comparison}}",
        comparison = lambda wildcards: json.dumps(DE_COMPARISONS[wildcards.comparison]),
        sample_key = SQUIDPY_DE.get("sample_key", "core_index"),
        cluster_key = SQUIDPY_DE.get("cluster_key", "leiden"),
        min_obs = SQUIDPY_DE.get("min_obs_per_pseudobulk", 50),
        min_samples = SQUIDPY_DE.get("min_samples_per_group", 2),
        min_gene_counts = SQUIDPY_DE.get("min_gene_counts", 10),
        alpha = SQUIDPY_DE.get("alpha", 0.05)
    threads: squidpy_resource("de", "cpus", 4)
    resources:
        mem_mb = squidpy_resource("de", "mem_mb", 32000),
        cpus = squidpy_resource("de", "cpus", 4),
        runtime = squidpy_resource("de", "runtime", 120),
        partition = squidpy_partition(),
        account = squidpy_account()
    log:
        f"{OUTPUT_DIR}/logs/squidpy_de_{{comparison}}_{{level}}.log"
    shell:
        """
        export PYTHONNOUSERSITE=1
        {SQUIDPY_PYTHON} {input.script} \
            --counts {input.combined} \
            --clustered {input.clustered} \
            --sample-metadata {input.sample_metadata} \
            --name {wildcards.comparison} \
            --comparison '{params.comparison}' \
            --outdir {params.outdir} \
            --sample-key {params.sample_key} \
            --cluster-key {params.cluster_key} \
            --min-obs-per-pseudobulk {params.min_obs} \
            --min-samples-per-group {params.min_samples} \
            --min-gene-counts {params.min_gene_counts} \
            --alpha {params.alpha} \
            --n-jobs {threads} > {log} 2>&1
        """


# Target rule: every DE comparison at every DE level
rule squidpy_de_all:
    input:
        expand(f"{SQUIDPY_RESULTS}/{{level}}/de/{{comparison}}/whole_core.csv",
               level=DE_LEVELS, comparison=sorted(DE_COMPARISONS))


# Target rule: the whole downstream spatial analysis
rule squidpy_all:
    input:
        expand(f"{SQUIDPY_DIR}/{{level}}/clustered.h5ad", level=LEVELS),
        expand(f"{SQUIDPY_RESULTS}/{{level}}/qc/qc_per_section.pdf", level=LEVELS),
        expand(f"{SQUIDPY_RESULTS}/{{level}}/qc/section_depth.pdf", level=LEVELS),
        expand(f"{SQUIDPY_RESULTS}/{{level}}/qc/section_qc.csv", level=LEVELS),
        expand(f"{SQUIDPY_RESULTS}/{{level}}/qc/quantile_bounds.csv", level=LEVELS),
        expand(f"{SQUIDPY_RESULTS}/{{level}}/svg/{{sample}}_spatial_autocorr.csv",
               level=LEVELS, sample=sorted(VISIUMHD_SAMPLES)),
        expand(f"{SQUIDPY_RESULTS}/{{level}}/spatial/{{sample}}_panels.pdf",
               level=LEVELS, sample=sorted(VISIUMHD_SAMPLES)),
        expand(f"{SQUIDPY_DIR}/{{level}}/spatial_stats.h5ad", level=HEAVY_LEVELS),
        expand(f"{SQUIDPY_RESULTS}/{{level}}/ccc/{{condition}}_liana_cellchat.csv",
               level=HEAVY_LEVELS, condition=CONDITIONS),
        expand(f"{SQUIDPY_DIR}/{{level}}/image_features/{{sample}}.csv",
               level=HEAVY_LEVELS, sample=sorted(VISIUMHD_SAMPLES)),
        expand(f"{SQUIDPY_RESULTS}/{{level}}/de/{{comparison}}/whole_core.csv",
               level=DE_LEVELS, comparison=sorted(DE_COMPARISONS))


# Rule: Generate multi-sample summary
rule multi_sample_summary:
    input:
        metrics_summaries = expand(f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/metrics_summary.csv", sample=SAMPLES)
    output:
        summary = f"{OUTPUT_DIR}/multi_sample_summary.csv"
    params:
        script = "src/filter_cellranger_out.py",
        mt_thresh = config.get('mt_thresh', 5),
        min_genes = config.get('min_genes', 200),
        min_cells = config.get('min_cells', 5),
        organism = config.get('organism', 'auto')
    threads: 8
    resources:
        mem_mb = 47663,  # 128GB in MB
        mem_mb_per_cpu=6000,
        cpus = 8,
        partition = "standard",
        account = "sbsandme_lab"
    run:
        import pandas as pd
        import os
        import subprocess
        import sys

        dfs = []
        for path in input.metrics_summaries:
            # Ensure the file exists and is not empty
            if not os.path.exists(path):
                print(f"Warning: metrics file not found: {path}")
                continue
            try:
                df = pd.read_csv(path)
            except pd.errors.EmptyDataError:
                print(f"Warning: metrics file empty: {path}")
                continue

            # Extract sample name from the path: .../cellranger/{sample}/outs/metrics_summary.csv
            parts = path.split(os.sep)
            sample = None
            # Try to locate the 'cellranger' dir and pick the next component as sample
            if 'cellranger' in parts:
                idx = parts.index('cellranger')
                if idx + 1 < len(parts):
                    sample = parts[idx + 1]
            if sample is None:
                # Fallback: use parent directory name two levels up
                sample = os.path.basename(os.path.dirname(os.path.dirname(path)))

            df['sample_id'] = sample

            # Call the filter script and capture the number of cells remaining
            matrix_h5 = os.path.join(OUTPUT_DIR, 'cellranger', sample, 'outs', 'filtered_feature_bc_matrix.h5')
            cmd = [sys.executable, params.script,
                   '--adata_path', matrix_h5,
                   '--min_genes', str(params.min_genes),
                   '--min_cells', str(params.min_cells),
                   '--mt_thresh', str(params.mt_thresh),
                   '--organism', str(params.organism)]
            # Run the filter script and capture stdout/stderr separately so warnings/errors
            # printed to stderr are not mixed into the stdout result.
            proc = subprocess.run(cmd, capture_output=True, text=True)
            if proc.returncode != 0:
                # Log stderr for debugging, but do not include it in the stored output
                print(f"Warning: filter script failed for sample {sample}. stderr:\n{proc.stderr}")
                out = proc.stdout.strip()
            else:
                out = proc.stdout.strip()

            df['cells_after_filtering'] = out
            df['min_cells'] = params.min_cells
            df['min_genes'] = params.min_genes
            df['mt_thresh'] = params.mt_thresh
            dfs.append(df)

        if not dfs:
            # No data to write; create an empty file with a header
            os.makedirs(os.path.dirname(output.summary), exist_ok=True)
            with open(output.summary, 'w') as fh:
                fh.write('')
            print(f"No metrics summaries found; created empty {output.summary}")
            return

        combined = pd.concat(dfs, sort=False, ignore_index=True)
        # Pop the column
        popped_column = combined.pop('sample_id')
        # Insert the column at the front (index 0)
        combined.insert(0, 'sample_id', popped_column)
        os.makedirs(os.path.dirname(output.summary), exist_ok=True)
        combined.to_csv(output.summary, index=False)
        print(f"Wrote combined metrics to {output.summary}")


# Rule: collect web_summary.html files into a single directory for easy viewing
rule collect_web_summaries:
    input:
        web_summaries = expand(f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/web_summary.html", sample=SAMPLES)
    output:
        directory(f"{OUTPUT_DIR}/web_summaries")
    params:
        outdir = f"{OUTPUT_DIR}/web_summaries"
    run:
        import os, shutil
        # Ensure output directory exists
        os.makedirs(params.outdir, exist_ok=True)
        for src in input.web_summaries:
            # Try to extract sample name from the path
            m = re.search(rf"{re.escape(OUTPUT_DIR)}/cellranger/([^/]+)/outs/web_summary.html", src)
            if m:
                sample = m.group(1)
            else:
                # Fallback: infer from parent directories
                sample = os.path.basename(os.path.dirname(os.path.dirname(src)))
            dest = os.path.join(params.outdir, f"{sample}_web_summary.html")
            shutil.copy(src, dest)
        # Touch a file to mark completion (optional; directory() is sufficient)
        open(os.path.join(params.outdir, ".done"), 'w').close()

# Rule: 10x scVI integration
rule tenx_scvi_integration:
    input:
       filtered_matrix_dirs = expand(f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/filtered_feature_bc_matrix", sample=SAMPLES)
    output:
        combined_adata = f"{OUTPUT_DIR}/scanpy/combined.h5ad",
        integration_results = f"{OUTPUT_DIR}/scanpy/combined_integrated.h5ad"
    conda: "scvi-tools"
    params:
        script = "src/tenx_scvi_integration.py",
        input_dir = f"{OUTPUT_DIR}/cellranger",
        min_genes = config.get("min_genes", 300),
        min_cells = config.get("min_cells", 5),
        n_top_genes = config.get("n_top_genes", 2000),
        batch_key = config.get("batch_key", "batch"),
        output_prefix = f"{OUTPUT_DIR}/scanpy/combined"
    threads: 4
    resources:
        mem_mb = 48000,  # 32GB in MB
        cpus = 8,
        partition = "gpu",
        account = "sbsandme_lab_gpu"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/scanpy
        python {params.script} \
            --filtered_matrix_dirs {input.filtered_matrix_dirs} \
            --output_prefix {params.output_prefix} \
            --min_genes {params.min_genes} \
            --min_cells {params.min_cells} \
            --n_top_genes {params.n_top_genes} \
            --batch_key {params.batch_key}
        """

# Rule: 10x harmony integration
rule tenx_harmony_integration:
    input:
       filtered_matrix_dirs = expand(f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/filtered_feature_bc_matrix", sample=SAMPLES)
    output:
        combined_adata = f"{OUTPUT_DIR}/scanpy/combined.h5ad",
        integration_results = f"{OUTPUT_DIR}/scanpy/combined_harmony_integrated.h5ad"
    conda: "scvi-tools"
    params:
        script = "src/tenx_harmony_integration.py",
        input_dir = f"{OUTPUT_DIR}/cellranger",
        min_genes = config.get("min_genes", 300),
        min_cells = config.get("min_cells", 5),
        n_top_genes = config.get("n_top_genes", 2000),
        batch_key = config.get("batch_key", "batch"),
        output_prefix = f"{OUTPUT_DIR}/scanpy/combined",
        metadata = config.get("metadata", None)
    threads: 8
    resources:
        mem_mb = 48000,  # 32GB in MB
        cpus = 8,
        account = "sbsandme_lab"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/scanpy
        python {params.script} \
            --filtered_matrix_dirs {input.filtered_matrix_dirs} \
            --output_prefix {params.output_prefix} \
            --min_genes {params.min_genes} \
            --min_cells {params.min_cells} \
            --n_top_genes {params.n_top_genes} \
            --batch_key {params.batch_key} \
            --metadata {params.metadata}
        """

# Rule: 10x harmony notebook
rule tenx_harmony_notebook:
    input:
       integration_results = f"{OUTPUT_DIR}/scanpy/combined_harmony_integrated.h5ad"
    output:
        integration_notebook = f"{OUTPUT_DIR}/scanpy/inspect_integrated_anndata_combined.ipynb"
    conda: "scvi-tools"
    params:
        script = "src/submit_harmony_integration.sh",
        input_dir = f"{OUTPUT_DIR}/cellranger",
        min_genes = config.get("min_genes", 300),
        min_cells = config.get("min_cells", 5),
        n_top_genes = config.get("n_top_genes", 2000),
        batch_key = config.get("batch_key", "batch"),
        output_prefix = f"{OUTPUT_DIR}/scanpy/combined"
    threads: 4
    resources:
        mem_mb = 32000,  # 32GB in MB
        cpus = 8,
        account = "sbsandme_lab"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/scanpy
        echo {input.integration_results}
        {params.script} --no-integration \
        --min-genes {params.min_genes} \
        {params.output_prefix}
        """


# Rule: per-sample filtering timeline plot
rule plot_filtering_timeline:
    input:
        matrix = f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/filtered_feature_bc_matrix.h5"
    output:
        png = f"{OUTPUT_DIR}/qc/filtering_timeline/{{sample}}.png"
    params:
        script = "src/plot_filtering_timeline.py",
        batch_key = config.get("batch_key", "batch"),
        min_genes = config.get("min_genes", 200),
        min_cells = config.get("min_cells", 3)
    threads: 4
    resources:
        mem_mb = 24000,  # 24GB in MB
        cpus = 4,
        partition = "standard",
        account = "sbsandme_lab"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/qc
        python {params.script} \
            --input {input.matrix} \
            --batch-key {params.batch_key} \
            --batch-value {wildcards.sample} \
            --min-genes {params.min_genes} \
            --min-cells {params.min_cells} \
            --out {output.png}
        """

# Rule: Velocyto
rule velocyto:
    input:
        sample_folder = f"{OUTPUT_DIR}/{{sample}}"
    output:
        loom = f"{OUTPUT_DIR}/velocyto/{{sample}}.loom"
    threads: 4
    shell:
        """
        module load velocyto/0.17.17
        velocyto run10x -m config["references"]["repeat_mask"] {input.sample_folder} config["gtf"] > {output.loom}
        module unload velocyto/0.17.17
        """

# Rule: scVelo
rule scvelo:
    input:
        loom_file = f"{OUTPUT_DIR}/velocyto/{{sample}}.loom",
        anndata_file = f"{OUTPUT_DIR}/scanpy/{{sample}}.h5ad",
        seurat_embeddings = f"{OUTPUT_DIR}/seurat/{{sample}}_embeddings.csv"
    output:
        scvelo_h5ad = f"{OUTPUT_DIR}/scanpy/{{sample}}_scvelo.h5ad"
    threads: 2
    shell:
        '''python src/compute_velocity.py {input.anndata_file} {input.loom_file} {input.seurat_embeddings}'''

# Rule: pileup and phasing
rule pileup_and_phasing:
    input:
        bam = f"{OUTPUT_DIR}/{{sample}}/outs/possorted_genome_bam.bam",
        barcodes = f"{OUTPUT_DIR}/{{sample}}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        script = "src/pileup_and_phase.R"
    output:
        allele_df = f"{OUTPUT_DIR}/numbat/{{sample}}_allele_counts.tsv.gz"
    threads: 4
    shell:
        '''Rscript {input.script} --label {wildcards.sample} --samples {wildcards.sample} --bams {input.bam} --barcodes {input.barcodes} --gmap config["gmap"] --snpvcf config["snpvcf"] --paneldir config["paneldir"] --outdir {OUTPUT_DIR}/numbat/{wildcards.sample} --ncores {threads}'''

# Rule: SCENIC analysis
include: "rules/scenic.smk"

# Rule: CellBender remove-background
rule cellbender:
    input:
        h5 = f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/raw_feature_bc_matrix.h5"
    output:
        clean_h5 = f"{OUTPUT_DIR}/cellbender/{{sample}}_filtered.h5"
    conda:
        "envs/cellbender_env.yaml"
    threads: 8
    resources:
        mem_mb = 64000,
        cpus = 8,
        partition = "standard"
    params:
        expected_cells = config.get("cellbender_expected_cells", 3000),
        total_droplets = config.get("cellbender_total_droplets", 20000)
    shell:
        """
        cellbender remove-background \
            --input {input.h5} \
            --output {output.clean_h5} \
            --expected-cells {params.expected_cells} \
            --total-droplets-included {params.total_droplets} \
            --cpu-threads {threads}
        """

# Rule: Scrublet doublet detection
rule scrublet:
    input:
        matrix = f"{OUTPUT_DIR}/cellranger/{{sample}}/outs/filtered_feature_bc_matrix.h5"
    output:
        doublets = f"{OUTPUT_DIR}/scrublet/{{sample}}_doublets.csv"
    threads: 2
    resources:
        mem_mb = 8000,
        cpus = 2,
        partition = "standard"
    params:
        script = "src/run_scrublet.py"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/scrublet
        python {params.script} --input {input.matrix} --output {output.doublets}
        """
