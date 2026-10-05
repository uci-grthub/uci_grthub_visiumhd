#!/usr/bin/env python3
"""
Read the manual TMA core annotations and map Visium HD bins to a core position.

Space Ranger emits one matrix per capture area with no notion of which of the four
2 mm TMA cores a bin came from, and every biological contrast in this experiment
(normal vs pigmented, week 0 vs week 16) lives *inside* a capture area. The cores are
therefore drawn by hand in Loupe and exported to
``metadata/regions/regions_tma<n>.csv``.

Only ``square_008um`` carries a ``.cloupe``, so the export is made at 8 um and has to be
transferred to the other bin sizes. Two export layouts are accepted and auto-detected:

* **Barcode export** -- a ``Barcode`` column plus one category column. Space Ranger HD
  barcodes encode the bin size and the array coordinates as
  ``s_008um_<row>_<col>-1``, so bins nest exactly: the 8 um bin at (r, c) sits inside the
  16 um bin at (r // 2, c // 2). No geometry is needed and the transfer is exact.
* **Polygon export** -- a ``Selection`` column with ``X``/``Y`` vertices in full-resolution
  pixels. Assigned by point-in-polygon against ``pxl_col_in_fullres`` /
  ``pxl_row_in_fullres`` from ``spatial/tissue_positions.parquet``, which makes it
  independent of bin size. Geometry construction is delegated to
  :func:`read_polygons_with_names.read_polygons_with_names`.

Bins that fall in no core are left unassigned; the caller drops them, which is what
repairs the runaway tissue detection on some capture areas (A1_TMA2 reports every one of
its 175,561 16 um bins as tissue, including bins at negative pixel coordinates).
"""

from __future__ import annotations

import csv
import re
import warnings

import numpy as np
import pandas as pd

from read_polygons_with_names import read_polygons_with_names

#: The four core positions, matching ``tma_position`` in ``metadata/sample_metadata.csv``.
POSITIONS = ("top-left", "top-right", "bottom-left", "bottom-right")

#: ``s_008um_00269_00526-1`` -> bin size in microns, array row, array column.
BARCODE_RE = re.compile(r"^s_(\d+)um_(\d+)_(\d+)-\d+$")

_POSITION_ALIASES = {
    "tl": "top-left",
    "tr": "top-right",
    "bl": "bottom-left",
    "br": "bottom-right",
    "upperleft": "top-left",
    "upperright": "top-right",
    "lowerleft": "bottom-left",
    "lowerright": "bottom-right",
}


def normalize_position(name):
    """Normalise a Loupe selection name to one of :data:`POSITIONS`.

    Loupe names survive round trips unpredictably -- ``Top Left``, ``top_left``,
    ``topleft`` and ``TL`` all mean the same core. Anything unrecognised is returned
    lowercased and stripped so the caller can report it verbatim.
    """
    text = re.sub(r"[^a-z0-9]", "", str(name).lower())
    if text in _POSITION_ALIASES:
        return _POSITION_ALIASES[text]
    for position in POSITIONS:
        if text == position.replace("-", ""):
            return position
    return str(name).strip().lower()


def parse_hd_barcodes(barcodes):
    """Decode Space Ranger HD barcodes into bin size and array coordinates.

    Returns a frame indexed by barcode with integer ``bin_um``, ``array_row`` and
    ``array_col`` columns. Barcodes that do not match :data:`BARCODE_RE` raise, because
    silently dropping them would quietly shrink the annotation.
    """
    series = pd.Series(list(barcodes), dtype="string")
    parsed = series.str.extract(BARCODE_RE)
    bad = parsed[0].isna()
    if bad.any():
        examples = ", ".join(series[bad].head(3).tolist())
        raise ValueError(
            f"{int(bad.sum())} barcode(s) are not Space Ranger Visium HD barcodes "
            f"(expected s_<bin>um_<row>_<col>-1), e.g. {examples}"
        )
    frame = pd.DataFrame(
        {
            "bin_um": parsed[0].astype(int).to_numpy(),
            "array_row": parsed[1].astype(int).to_numpy(),
            "array_col": parsed[2].astype(int).to_numpy(),
        },
        index=pd.Index(series, name="barcode"),
    )
    return frame


def _sniff_layout(path):
    """Return ('barcode'|'polygon', skiprows, header) by looking for the header row.

    Loupe writes a short preamble before the header in some export types, so the header
    is located rather than assumed.
    """
    with open(path, newline="") as handle:
        for index, row in enumerate(csv.reader(handle)):
            if index > 10:
                break
            fields = {str(cell).strip().lower() for cell in row}
            if not fields:
                continue
            if "barcode" in fields:
                return "barcode", index, [str(cell).strip() for cell in row]
            if "selection" in fields and {"x", "y"} <= fields:
                return "polygon", index, [str(cell).strip() for cell in row]
    raise ValueError(
        f"{path}: could not find a header row with either a 'Barcode' column or "
        "'Selection' plus 'X'/'Y' columns in the first 10 lines"
    )


def _read_barcode_annotation(path, skiprows, header):
    """Read a Loupe barcode export into a barcode -> position Series."""
    frame = pd.read_csv(path, skiprows=skiprows)
    barcode_col = next(c for c in frame.columns if str(c).strip().lower() == "barcode")
    value_cols = [c for c in frame.columns if c != barcode_col]
    if not value_cols:
        raise ValueError(
            f"{path}: found a 'Barcode' column but no category column alongside it "
            f"(columns: {', '.join(map(str, header))})"
        )
    if len(value_cols) > 1:
        warnings.warn(
            f"{path}: {len(value_cols)} category columns present "
            f"({', '.join(map(str, value_cols))}); using {value_cols[0]!r}",
            stacklevel=2,
        )
    labels = frame[value_cols[0]]
    # Loupe's "include unlabeled" export tags every bin outside a selection "Unlabeled"
    stripped = labels.astype(str).str.strip()
    assigned = labels.notna() & (stripped != "") & (stripped.str.lower() != "unlabeled")
    result = pd.Series(
        labels[assigned].map(normalize_position).to_numpy(),
        index=pd.Index(frame.loc[assigned, barcode_col].astype(str), name="barcode"),
        name="tma_position",
        dtype="object",
    )
    return result[~result.index.duplicated(keep="first")]


def _transfer_between_bins(annotation, target_positions):
    """Carry a barcode annotation from the bin size it was drawn at to another.

    Visium HD bins nest exactly, so the mapping is integer arithmetic on the array
    coordinates rather than geometry. Coarsening (8 um -> 16 um) puts several source bins
    in one target bin and is resolved by majority vote; refining (8 um -> 2 um) gives every
    target bin the label of the source bin containing it.
    """
    source = parse_hd_barcodes(annotation.index)
    source["tma_position"] = annotation.to_numpy()

    source_um = sorted(source["bin_um"].unique())
    if len(source_um) != 1:
        raise ValueError(f"annotation mixes bin sizes: {source_um}")
    source_um = int(source_um[0])

    target = parse_hd_barcodes(target_positions.index)
    target_um = sorted(target["bin_um"].unique())
    if len(target_um) != 1:
        raise ValueError(f"target matrix mixes bin sizes: {target_um}")
    target_um = int(target_um[0])

    if source_um == target_um:
        return annotation.reindex(target_positions.index)

    if target_um % source_um == 0:
        # coarsening: many source bins per target bin, decided by majority vote
        factor = target_um // source_um
        keys = pd.MultiIndex.from_arrays(
            [source["array_row"] // factor, source["array_col"] // factor]
        )
        votes = source.groupby([keys.get_level_values(0), keys.get_level_values(1)])[
            "tma_position"
        ].agg(lambda values: values.value_counts().idxmax())
    elif source_um % target_um == 0:
        # refining: every target bin inherits from the source bin that contains it
        factor = source_um // target_um
        target_keys = pd.MultiIndex.from_arrays(
            [target["array_row"] // factor, target["array_col"] // factor]
        )
        votes = pd.Series(
            source["tma_position"].to_numpy(),
            index=pd.MultiIndex.from_arrays([source["array_row"], source["array_col"]]),
        )
        votes = votes[~votes.index.duplicated(keep="first")]
        return pd.Series(
            votes.reindex(target_keys).to_numpy(),
            index=target_positions.index,
            name="tma_position",
            dtype="object",
        )
    else:
        raise ValueError(
            f"cannot transfer a {source_um} um annotation to {target_um} um bins: "
            "the bin sizes do not nest"
        )

    target_keys = pd.MultiIndex.from_arrays([target["array_row"], target["array_col"]])
    return pd.Series(
        votes.reindex(target_keys).to_numpy(),
        index=target_positions.index,
        name="tma_position",
        dtype="object",
    )


def _read_polygon_annotation(path, skiprows, target_positions):
    """Assign bins by point-in-polygon against full-resolution pixel coordinates."""
    try:
        from shapely.geometry import Point
        from shapely.prepared import prep
    except ImportError as error:
        raise SystemExit(
            f"{path} is a polygon export, which needs shapely: {error}"
        ) from error

    polygons = read_polygons_with_names(path, save_csvs=False, skiprows=skiprows)

    missing = {"pxl_col_in_fullres", "pxl_row_in_fullres"} - set(target_positions.columns)
    if missing:
        raise ValueError(
            "tissue_positions is missing the full-resolution pixel columns needed for "
            f"polygon assignment: {', '.join(sorted(missing))}"
        )

    xs = target_positions["pxl_col_in_fullres"].to_numpy(dtype=float)
    ys = target_positions["pxl_row_in_fullres"].to_numpy(dtype=float)
    assignment = np.full(len(target_positions), None, dtype=object)

    for name, geometry in polygons.items():
        if geometry is None or not hasattr(geometry, "contains"):
            warnings.warn(
                f"{path}: selection {name!r} has too few vertices for a polygon; skipped",
                stacklevel=2,
            )
            continue
        minx, miny, maxx, maxy = geometry.bounds
        # cheap bounding-box prefilter, then the exact test only on survivors
        candidates = np.flatnonzero(
            (xs >= minx) & (xs <= maxx) & (ys >= miny) & (ys <= maxy)
        )
        if candidates.size == 0:
            continue
        prepared = prep(geometry)
        position = normalize_position(name)
        for index in candidates:
            if prepared.contains(Point(xs[index], ys[index])):
                assignment[index] = position

    return pd.Series(
        assignment, index=target_positions.index, name="tma_position", dtype="object"
    )


def read_core_annotation(path, target_positions=None):
    """Map the bins of one capture area to a TMA core position.

    Args:
        path: the Loupe export, ``metadata/regions/regions_tma<n>.csv``.
        target_positions: the ``tissue_positions.parquet`` frame for the bin size being
            built, indexed by barcode. Required for a polygon export and for transferring
            a barcode export between bin sizes. When omitted, the barcode export is
            returned as drawn.

    Returns:
        A ``tma_position`` Series indexed by barcode. Unassigned bins are ``None``.
    """
    layout, skiprows, header = _sniff_layout(path)

    if layout == "barcode":
        annotation = _read_barcode_annotation(path, skiprows, header)
        if target_positions is not None:
            annotation = _transfer_between_bins(annotation, target_positions)
    else:
        if target_positions is None:
            raise ValueError(
                f"{path} is a polygon export, which needs the tissue_positions frame to "
                "assign bins"
            )
        annotation = _read_polygon_annotation(path, skiprows, target_positions)

    found = sorted({p for p in annotation.dropna().unique()})
    unexpected = [p for p in found if p not in POSITIONS]
    if unexpected:
        warnings.warn(
            f"{path}: selection name(s) {unexpected} do not match the expected "
            f"tma_position values {list(POSITIONS)}; they will not join to "
            "sample_metadata.csv",
            stacklevel=2,
        )
    if len(found) < len(POSITIONS):
        warnings.warn(
            f"{path}: only {len(found)} core(s) annotated ({found}); expected "
            f"{len(POSITIONS)}",
            stacklevel=2,
        )
    return annotation


def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("annotation", help="metadata/regions/regions_tma<n>.csv")
    parser.add_argument(
        "--positions",
        help="spatial/tissue_positions.parquet for the bin size to assign",
    )
    args = parser.parse_args()

    positions = None
    if args.positions:
        positions = pd.read_parquet(args.positions).set_index("barcode")

    annotation = read_core_annotation(args.annotation, positions)
    counts = annotation.value_counts(dropna=False)
    print(f"{args.annotation}: {len(annotation)} bins")
    for position, count in counts.items():
        label = "unassigned" if position is None else position
        print(f"  {label:>14}: {count}")


if __name__ == "__main__":
    main()
