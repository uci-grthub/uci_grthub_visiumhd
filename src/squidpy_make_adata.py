#!/usr/bin/env python3
"""
Build one AnnData per Visium HD capture area, with TMA cores assigned to bins or cells.

``squidpy.read.visium`` cannot read Visium HD: it looks for
``spatial/tissue_positions_list.csv`` and Space Ranger 3.1.2 writes
``spatial/tissue_positions.parquet``. The object is therefore assembled here from the
binned feature matrix, the positions parquet and the scale factors.

Bins are then assigned to one of the four TMA cores using the manual Loupe annotation
(see ``src/squidpy_core_annotation.py`` and ``docs/core_annotation.md``) and joined to
``metadata/sample_metadata.csv`` to attach patient, condition and timepoint. Bins that
fall in no core are dropped: the cores occupy roughly 38% of a 6.5 mm capture area, and
on several areas Space Ranger's tissue detection called the entire capture area tissue
(A1_TMA2 reports all 175,561 of its 16 um bins as tissue, including bins at negative
pixel coordinates), so the core annotation doubles as the tissue mask.

``--level cell`` builds the object from Space Ranger's segmented outputs instead
(``segmented_outputs/filtered_feature_cell_matrix.h5``: H&E nuclei expanded into cells).
Cell positions are the ``cell_centroid`` of ``cell_segmentations.geojson``. A cell is put
in a core by majority vote over its 2 um bins in ``barcode_mappings.parquet``, each
carrying the core label of the 8 um bin that contains it, so the same 8 um Loupe export
serves both levels. The Azimuth cell types Space Ranger writes are attached when present.
"""

import argparse
import json
import os
import re

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import scanpy as sc

from squidpy_core_annotation import POSITIONS, read_core_annotation


def load_binned_matrix(matrix_h5, positions, scalefactors, run_id, image_dir=None):
    """Assemble an AnnData for one capture area at one bin size."""
    adata = sc.read_10x_h5(matrix_h5)
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()

    frame = pd.read_parquet(positions).set_index("barcode")
    missing = adata.obs_names.difference(frame.index)
    if len(missing):
        raise SystemExit(
            f"{len(missing)} barcode(s) in {matrix_h5} are absent from {positions}; "
            "the matrix and the positions file are not from the same run"
        )
    adata.obs = adata.obs.join(frame, how="left")

    # x = column, y = row, in full-resolution pixels of the microscope image
    adata.obsm["spatial"] = adata.obs[
        ["pxl_col_in_fullres", "pxl_row_in_fullres"]
    ].to_numpy(dtype=float)

    with open(scalefactors) as handle:
        factors = json.load(handle)

    images = read_library_images(image_dir) if image_dir else {}
    adata.uns["spatial"] = {run_id: {"scalefactors": factors, "images": images}}

    return adata


def read_library_images(image_dir):
    """The hires/lowres H&E PNGs Space Ranger writes next to the scale factors."""
    images = {}
    for key, name in (("hires", "tissue_hires_image.png"),
                      ("lowres", "tissue_lowres_image.png")):
        path = os.path.join(image_dir, name)
        if os.path.exists(path):
            images[key] = plt.imread(path)
        else:
            print(f"Warning: no {name} in {image_dir}; plots will have no H&E")
    return images


def load_cell_matrix(matrix_h5, segmentations, scalefactors, run_id, cell_types=None):
    """Assemble an AnnData for one capture area from the segmented cell outputs."""
    adata = sc.read_10x_h5(matrix_h5)
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()

    with open(segmentations) as handle:
        features = json.load(handle)["features"]
    rows, areas = [], []
    for feature in features:
        props = feature["properties"]
        x, y = props["cell_centroid"]
        rows.append((f"cellid_{int(props['cell_id']):09d}-1", x, y))
        ring = np.asarray(feature["geometry"]["coordinates"][0], dtype=float)
        # shoelace area of the outer ring, in full-resolution pixels
        areas.append(0.5 * abs(np.dot(ring[:-1, 0], ring[1:, 1])
                               - np.dot(ring[1:, 0], ring[:-1, 1])))
    centroids = pd.DataFrame(rows, columns=["barcode", "x", "y"]).set_index("barcode")

    missing = adata.obs_names.difference(centroids.index)
    if len(missing):
        raise SystemExit(
            f"{len(missing)} cell(s) in {matrix_h5} have no polygon in {segmentations}, "
            f"e.g. {', '.join(missing[:3])}"
        )
    # x = column, y = row, in full-resolution pixels, as for the bins
    adata.obsm["spatial"] = centroids.loc[adata.obs_names, ["x", "y"]].to_numpy(dtype=float)

    if cell_types and os.path.exists(cell_types):
        types = pd.read_csv(cell_types).set_index("barcode")
        keep = [c for c in ("broad_cell_type", "coarse_cell_type", "fine_cell_type",
                            "final_level_softmax_prob") if c in types.columns]
        adata.obs = adata.obs.join(types[keep], how="left")
        print(f"{run_id}: attached {', '.join(keep)} from {cell_types}")

    with open(scalefactors) as handle:
        factors = json.load(handle)
    # the segmented scale factors carry no spot size, which squidpy's plots and image
    # crops need: use the diameter of a circle with the median cell area
    factors.setdefault("spot_diameter_fullres", float(2 * np.sqrt(np.median(areas) / np.pi)))
    images = read_library_images(os.path.dirname(scalefactors))
    adata.uns["spatial"] = {run_id: {"scalefactors": factors, "images": images}}
    return adata


def assign_cores(adata, annotation_path, positions, sample_metadata, run_id):
    """Attach tma_position and the per-core metadata, dropping unassigned bins."""
    frame = pd.read_parquet(positions).set_index("barcode")
    assignment = read_core_annotation(annotation_path, frame)
    adata.obs["tma_position"] = assignment.reindex(adata.obs_names).to_numpy()
    return attach_core_metadata(adata, annotation_path, sample_metadata, run_id, "bins")


def assign_cores_to_cells(adata, annotation_path, bin_positions, barcode_mappings,
                          sample_metadata, run_id):
    """Put each cell in the core that holds most of its 2 um bins."""
    frame = pd.read_parquet(bin_positions).set_index("barcode")
    bin_um = int(re.match(r"s_(\d+)um_", frame.index[0]).group(1))
    labels = read_core_annotation(annotation_path, frame)

    # one row per 2 um bin inside a cell; only two columns, the full table is ~11M rows
    mapping = pq.read_table(
        barcode_mappings, columns=[f"square_{bin_um:03d}um", "cell_id"],
        filters=[("in_cell", "=", True)],
    ).to_pandas()
    mapping["tma_position"] = mapping.iloc[:, 0].map(labels)
    votes = (mapping.dropna(subset=["tma_position"])
             .groupby("cell_id")["tma_position"]
             .agg(lambda values: values.value_counts().idxmax()))
    adata.obs["tma_position"] = votes.reindex(adata.obs_names).to_numpy()
    return attach_core_metadata(adata, annotation_path, sample_metadata, run_id, "cells")


def attach_core_metadata(adata, annotation_path, sample_metadata, run_id, unit):
    """Drop observations outside every core and join the per-core sample metadata."""
    assigned = adata.obs["tma_position"].notna()
    print(
        f"{run_id}: {int(assigned.sum())} of {adata.n_obs} {unit} fall in a core "
        f"({100 * assigned.mean():.1f}%)"
    )
    if not assigned.any():
        raise SystemExit(
            f"{run_id}: no {unit[:-1]} was assigned to a core. Check that {annotation_path} "
            "covers this capture area and that the selections are named "
            f"{list(POSITIONS)}."
        )
    adata = adata[assigned].copy()

    metadata = pd.read_csv(sample_metadata)
    cores = metadata[metadata["cytassist_sample"] == run_id]
    if cores.empty:
        raise SystemExit(
            f"{run_id} does not appear in the cytassist_sample column of "
            f"{sample_metadata}"
        )

    columns = [
        "core_id", "sample_name", "tma", "tma_position", "core_index", "replicate",
        "patient", "patient_id", "condition", "timepoint_weeks", "diagnosis",
        "sample_type", "visium_slide", "visium_slide_id", "capture_area",
    ]
    cores = cores[[c for c in columns if c in cores.columns]]

    duplicated = cores["tma_position"].duplicated()
    if duplicated.any():
        # TMA6 carries two replicates of the same core, distinguished by position
        raise SystemExit(
            f"{run_id}: {sample_metadata} has more than one row for position(s) "
            f"{sorted(cores.loc[duplicated, 'tma_position'])}"
        )

    merged = adata.obs.reset_index().merge(cores, on="tma_position", how="left")
    merged.index = adata.obs_names
    unmatched = merged["core_id"].isna()
    if unmatched.any():
        positions_seen = sorted(adata.obs.loc[unmatched.to_numpy(), "tma_position"].unique())
        raise SystemExit(
            f"{run_id}: core position(s) {positions_seen} from {annotation_path} have no "
            f"row in {sample_metadata}. Expected {sorted(cores['tma_position'])}."
        )

    adata.obs = merged.drop(columns=[merged.columns[0]])
    adata.obs["run_id"] = run_id

    # One batch per biological sample, 22 across the experiment against 24 cores: the two
    # A1_TMA6 cores of each condition are technical replicates of the same sample and so
    # carry the same sample_id. This is the key Harmony integrates on.
    patient_label = (adata.obs["patient_id"] if "patient_id" in adata.obs
                     else adata.obs["patient"].astype(str).radd("P"))
    adata.obs["sample_id"] = (
        patient_label.astype(str)
        + "_" + adata.obs["condition"].astype(str)
        + "_W" + adata.obs["timepoint_weeks"].astype(str)
    )

    for column in ("core_id", "sample_name", "sample_id", "condition", "core_index",
                   "tma_position", "diagnosis", "sample_type", "run_id"):
        if column in adata.obs:
            adata.obs[column] = adata.obs[column].astype("category")
    # numeric labels, not continuous scales: categorical so plots get a legend rather
    # than a colorbar. Timepoints keep week order; patient ids are just sorted.
    for column, ordered in (("timepoint_weeks", True), ("patient", False)):
        if column in adata.obs:
            values = adata.obs[column]
            adata.obs[column] = pd.Categorical(
                values, categories=sorted(values.dropna().unique()), ordered=ordered)
    return adata


def plot_core_assignment(adata, run_id, path, unit="bins"):
    """Spatial scatter coloured by core -- the check that the selections landed right."""
    coords = adata.obsm["spatial"]
    figure, axis = plt.subplots(figsize=(6, 6))
    for core, group in adata.obs.groupby("core_id", observed=True):
        mask = adata.obs["core_id"] == core
        label = f"{core} ({group['sample_name'].iloc[0]})" if len(group) else str(core)
        axis.scatter(coords[mask, 0], coords[mask, 1], s=0.5, linewidths=0, label=label)
    axis.set_aspect("equal")
    axis.invert_yaxis()  # image row order
    axis.set_xlabel("pxl_col_in_fullres")
    axis.set_ylabel("pxl_row_in_fullres")
    axis.set_title(f"{run_id}: core assignment ({adata.n_obs} {unit})")
    axis.legend(markerscale=12, fontsize="small", loc="center left",
                bbox_to_anchor=(1.02, 0.5), frameon=False)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    figure.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(figure)
    print(f"Wrote {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level", choices=["bin", "cell"], default="bin",
                        help="binned outputs, or the segmented cell outputs")
    parser.add_argument("--matrix-h5", required=True,
                        help="binned_outputs/square_<bin>um/filtered_feature_bc_matrix.h5, "
                             "or segmented_outputs/filtered_feature_cell_matrix.h5")
    parser.add_argument("--positions", required=True,
                        help="binned_outputs/square_<bin>um/spatial/tissue_positions.parquet; "
                             "for --level cell, the one at the bin size the cores were drawn at")
    parser.add_argument("--scalefactors", required=True,
                        help="spatial/scalefactors_json.json for the level being built")
    parser.add_argument("--segmentations",
                        help="--level cell: segmented_outputs/cell_segmentations.geojson")
    parser.add_argument("--barcode-mappings",
                        help="--level cell: outs/barcode_mappings.parquet")
    parser.add_argument("--cell-types",
                        help="--level cell: segmented_outputs/cell_types/Azimuth/cell_types.csv")
    parser.add_argument("--core-annotation", required=True,
                        help="metadata/regions/regions_tma<n>.csv")
    parser.add_argument("--sample-metadata", required=True,
                        help="metadata/sample_metadata.csv")
    parser.add_argument("--run-id", required=True, help="Space Ranger run id")
    parser.add_argument("--output", required=True, help="h5ad to write")
    parser.add_argument("--qc-figure", help="core-assignment figure to write")
    args = parser.parse_args()

    if args.level == "cell":
        if not (args.segmentations and args.barcode_mappings):
            parser.error("--level cell needs --segmentations and --barcode-mappings")
        unit = "cells"
        adata = load_cell_matrix(args.matrix_h5, args.segmentations, args.scalefactors,
                                 args.run_id, args.cell_types)
        print(f"{args.run_id}: loaded {adata.n_obs} cells x {adata.n_vars} genes")
        adata = assign_cores_to_cells(adata, args.core_annotation, args.positions,
                                      args.barcode_mappings, args.sample_metadata,
                                      args.run_id)
    else:
        unit = "bins"
        image_dir = os.path.dirname(args.positions)
        adata = load_binned_matrix(
            args.matrix_h5, args.positions, args.scalefactors, args.run_id, image_dir
        )
        print(f"{args.run_id}: loaded {adata.n_obs} bins x {adata.n_vars} genes")
        adata = assign_cores(
            adata, args.core_annotation, args.positions, args.sample_metadata, args.run_id
        )

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # probe-based FFPE assay: 11 MT- probes are in the panel, so pct_counts_mt is
    # computable but is a weak QC axis next to counts and genes per bin or cell
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], inplace=True, log1p=False, percent_top=None
    )

    if args.qc_figure:
        plot_core_assignment(adata, args.run_id, args.qc_figure, unit)

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    adata.write_h5ad(args.output, compression="gzip")
    print(f"Wrote {adata.n_obs} {unit} x {adata.n_vars} genes to {args.output}")
    print(
        adata.obs.groupby(["sample_id", "core_id"], observed=True)
        .size()
        .rename(unit)
        .to_string()
    )


if __name__ == "__main__":
    main()
