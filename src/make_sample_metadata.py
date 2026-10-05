#!/usr/bin/env python
"""Build a per-core sample metadata CSV from the Visium TMA design spreadsheet.

The spreadsheet has two halves:

* columns A-F: one row per 2 mm biopsy core, grouped in blocks of four by the
  6.5 mm capture area they were arrayed onto.
* columns J-L: for each TMA, a 2x2 map of which core index sits in which
  position of the array, read row-major (1 = top-left, 2 = top-right,
  3 = bottom-left, 4 = bottom-right).

The CytAssist run-info files map each TMA to the Visium slide serial and
capture area it was transferred to; they are joined in when available.

The clinical spreadsheet (``--clinical-xlsx``, e.g. "9-30-26 Meta data for 10X visium
seq data analysis.xlsx") has one row per core, keyed by ``No`` == core index, with
patient id (M = melasma, P = PIH), age, sex, skin type, sample type and timepoint. It
adds the patient-level diagnosis and the per-core sample_type (normal / melasma / PIH)
that the differential expression comparisons are defined on.
"""

import argparse
import csv
import re
from pathlib import Path

import openpyxl

POSITION_NAMES = {
    (1, 1): "top-left",
    (1, 2): "top-right",
    (2, 1): "bottom-left",
    (2, 2): "bottom-right",
}

FIELDS = [
    "core_id",
    "sample_name",
    "tma",
    "core_row",
    "core_col",
    "tma_position",
    "core_index",
    "replicate",
    "patient",
    "condition",
    "timepoint_weeks",
    "capture_area_label",
    "visium_slide",
    "visium_slide_id",
    "capture_area",
    "cytassist_sample",
    "tissue_slide_number",
    "slides_requested",
    "comments",
]

CLINICAL_FIELDS = [
    "patient_id",
    "age",
    "sex",
    "skin_type",
    "diagnosis",
    "sample_type",
    "de_group",
]

DIAGNOSES = {"M": "melasma", "P": "PIH"}


def read_cores(worksheet):
    """Return {core_index: {...}} from columns A-F, carrying merged cells down."""
    cores = {}
    ordered = []
    capture_label = tissue_slide = requested = None
    for row in range(2, worksheet.max_row + 1):
        name = worksheet.cell(row=row, column=2).value
        if name is None:
            continue
        # C/D/E are merged over each block of four cores; only the first row carries a value
        capture_label = worksheet.cell(row=row, column=3).value or capture_label
        tissue_slide = worksheet.cell(row=row, column=4).value or tissue_slide
        requested = worksheet.cell(row=row, column=5).value or requested
        record = {
            "sample_name": str(name).strip(),
            "capture_area_label": str(capture_label).strip(),
            "tissue_slide_number": tissue_slide,
            "slides_requested": requested,
            "comments": str(worksheet.cell(row=row, column=6).value or "").strip(),
        }
        index = worksheet.cell(row=row, column=1).value
        if index is not None:
            cores[int(index)] = record
        ordered.append(record)
    return cores


def read_tma_maps(worksheet):
    """Return {tma_number: {(row, col): core_index}} from the J/K/L layout blocks."""
    maps = {}
    for row in range(1, worksheet.max_row + 1):
        header = worksheet.cell(row=row, column=10).value
        match = re.fullmatch(r"\s*TMA\s*(\d+)\s*", str(header or ""))
        if not match:
            continue
        grid = {}
        for grid_row in (1, 2):
            for grid_col, column in ((1, 11), (2, 12)):
                value = worksheet.cell(row=row + grid_row, column=column).value
                if value is not None:
                    grid[(grid_row, grid_col)] = int(value)
        maps[int(match.group(1))] = grid
    return maps


def read_cytassist(cytassist_dir):
    """Return {tma_number: {slide_id, capture_area, cytassist_sample}}."""
    placements = {}
    if cytassist_dir is None:
        return placements
    for info in sorted(Path(cytassist_dir).glob("*/*run-info.csv")):
        slide_id = None
        samples = {}
        with open(info, newline="") as handle:
            for fields in csv.reader(handle):
                if len(fields) < 2:
                    continue
                key, value = fields[0].strip(), fields[1].strip()
                if key == "Visium Slide ID":
                    slide_id = value
                elif key.startswith("Sample "):
                    samples[key.split(None, 1)[1]] = value
        for area, sample in samples.items():
            match = re.search(r"TMA\s*(\d+)", sample)
            if not match:
                continue
            placements[int(match.group(1))] = {
                "visium_slide_id": slide_id,
                "capture_area": area,
                "cytassist_sample": sample,
            }
    return placements


def parse_comment(comment):
    """Pull patient, condition and timepoint out of e.g. 'patient #4 16 weeks pigmented'."""
    patient = re.search(r"patient\s*#?\s*(\d+)", comment, re.I)
    weeks = re.search(r"(\d+)\s*weeks?", comment, re.I)
    if re.search(r"pigment", comment, re.I):
        condition = "pigmented"
    elif re.search(r"normal", comment, re.I):
        condition = "normal"
    else:
        condition = ""
    return {
        "patient": patient.group(1) if patient else "",
        "condition": condition,
        "timepoint_weeks": weeks.group(1) if weeks else "0",
    }


def parse_capture_label(label):
    """'6.5mm x 6.5mm area 2 Visium slide 1' -> slide number."""
    slide = re.search(r"slide\s*(\d+)", label, re.I)
    return slide.group(1) if slide else ""


def build_rows(worksheet, cytassist_dir):
    cores = read_cores(worksheet)
    maps = read_tma_maps(worksheet)
    placements = read_cytassist(cytassist_dir)

    rows = []
    for tma in sorted(maps):
        seen = {}
        for (grid_row, grid_col), index in sorted(maps[tma].items()):
            core = cores.get(index)
            if core is None:
                raise KeyError(f"TMA{tma} position {grid_row},{grid_col} refers to unknown core {index}")
            seen[index] = seen.get(index, 0) + 1
            replicate = seen[index]
            row = {
                "core_id": f"TMA{tma}_{POSITION_NAMES[(grid_row, grid_col)]}",
                "sample_name": core["sample_name"],
                "tma": tma,
                "core_row": grid_row,
                "core_col": grid_col,
                "tma_position": POSITION_NAMES[(grid_row, grid_col)],
                "core_index": index,
                "replicate": replicate,
                "capture_area_label": core["capture_area_label"],
                "visium_slide": parse_capture_label(core["capture_area_label"]),
                "tissue_slide_number": core["tissue_slide_number"],
                "slides_requested": core["slides_requested"],
                "comments": core["comments"],
                "visium_slide_id": "",
                "capture_area": "",
                "cytassist_sample": "",
            }
            row.update(parse_comment(core["comments"]))
            row.update(placements.get(tma, {}))
            rows.append(row)
    return rows


def read_clinical(xlsx):
    """Return {core_index: {...}} from the clinical sheet, one row per core."""
    worksheet = openpyxl.load_workbook(xlsx, data_only=True).active
    rows = worksheet.iter_rows(values_only=True)
    header = [str(h).strip() if h is not None else "" for h in next(rows)]
    column = {name: header.index(name) for name in
              ("No", "Patient ID", "Age", "Sex", "Skin type", "Sample type",
               "Treatment timepoint")}
    clinical = {}
    for values in rows:
        index = values[column["No"]]
        if not isinstance(index, (int, float)):
            continue  # footnote rows
        patient_id = str(values[column["Patient ID"]]).strip()
        diagnosis = DIAGNOSES.get(patient_id[:1].upper())
        if diagnosis is None:
            raise SystemExit(f"{xlsx}: patient id {patient_id} is neither M (melasma) "
                             "nor P (PIH)")
        sample_type = str(values[column["Sample type"]]).strip().lower()
        timepoint = str(values[column["Treatment timepoint"]]).strip().lower()
        clinical[int(index)] = {
            "patient_id": patient_id,
            "age": values[column["Age"]],
            "sex": values[column["Sex"]],
            "skin_type": values[column["Skin type"]],
            "diagnosis": diagnosis,
            "sample_type": "normal" if sample_type.startswith("normal") else diagnosis,
            "_pigmented": sample_type.startswith("pigment"),
            "_weeks": re.search(r"(\d+)\s*weeks?", timepoint),
            "_baseline": timepoint == "baseline",
        }
    return clinical


def add_clinical(rows, clinical, xlsx):
    """Join the clinical sheet onto rows by core_index, checking the two files agree."""
    for row in rows:
        entry = clinical.get(int(row["core_index"]))
        if entry is None:
            raise SystemExit(f"core {row['core_index']} ({row['core_id']}) is not in {xlsx}")
        weeks = entry["_weeks"].group(1) if entry["_weeks"] else ("0" if entry["_baseline"] else None)
        problems = []
        if str(int(re.sub(r"\D", "", entry["patient_id"]))) != str(row["patient"]):
            problems.append(f"patient {entry['patient_id']} vs {row['patient']}")
        if ("pigmented" if entry["_pigmented"] else "normal") != row["condition"]:
            problems.append(f"sample type {entry['sample_type']} vs {row['condition']}")
        if weeks != str(row["timepoint_weeks"]):
            problems.append(f"timepoint {weeks} vs {row['timepoint_weeks']}")
        if problems:
            raise SystemExit(f"core {row['core_index']} ({row['core_id']}) disagrees with "
                             f"{xlsx}: {'; '.join(problems)}")
        row.update({key: value for key, value in entry.items() if not key.startswith("_")})
        row["de_group"] = f"{entry['sample_type']}_W{row['timepoint_weeks']}"
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("xlsx", type=Path)
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument(
        "--cytassist-dir",
        type=Path,
        default=None,
        help="directory of CytAssist run folders, each holding a *-run-info.csv",
    )
    parser.add_argument(
        "--clinical-xlsx",
        type=Path,
        default=None,
        help="per-core clinical sheet (patient id, age, sex, skin type, sample type)",
    )
    args = parser.parse_args()

    worksheet = openpyxl.load_workbook(args.xlsx, data_only=True).active
    rows = build_rows(worksheet, args.cytassist_dir)
    fields = FIELDS
    if args.clinical_xlsx is not None:
        rows = add_clinical(rows, read_clinical(args.clinical_xlsx), args.clinical_xlsx)
        fields = FIELDS + CLINICAL_FIELDS

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} cores to {args.output}")


if __name__ == "__main__":
    main()
