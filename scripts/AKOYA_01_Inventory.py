#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler inventory / schema audit
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 01 of the AKOYA analysis series. REVISION 2.

PURPOSE
    Pure read-only inventory. Reads the 8 QuPath-exported CSVs, reports what is
    actually in them, and writes summary tables. Does NOT normalize, threshold,
    re-phenotype, or modify anything.

WHAT CHANGED IN REVISION 2 (and why)
    1. SCAN IDENTITY IS NOW AN OUTPUT COLUMN.
       Script 02's between-slide versus within-slide variance decomposition keyed
       its "slide" off the treatment column of 01_file_summary.csv, because this
       script never wrote a scan identifier anywhere downstream could read it.
       Slide and treatment are perfectly confounded in this design, so that
       substitution was harmless arithmetically but it made the resulting
       quantity impossible to describe honestly: it is the fraction of variance
       explained by arm, which is batch and biology combined, not batch alone.
       This script now writes scan_image_value (the raw Image string, which is
       the join key) and scan_id (a short stable label) so script 02 can key on
       the acquisition rather than on the treatment label.

       It also warns when a section file carries more than one Image value.
       Script 02 previously took df[Image].iloc[0], which is silently wrong if a
       file spans more than one acquisition.

    2. SLIDE POSITION IS NOW AN OUTPUT COLUMN.
       The exclusion of G3_43102 and G4_43112 rests on their being the
       position-1 section of their scan, flush against the scan boundary, with
       the shortest Y span. Those quantities were reconstructed inside script 02
       every time they were needed. They are inventory facts, so they belong
       here: slide_position_rank, y_span_um, and y_offset_from_scan_min_um.
       Rank 1 is the smallest y_min within a scan, matching the convention
       script 02 already used.

    3. SCAN IS PROPAGATED INTO THE LONG TABLES.
       04_phenotype_counts_long.csv and 06_marker_stats_long.csv now carry
       scan_id and scan_image_value so they are self-describing and downstream
       scripts do not have to join back through the file summary.

    Nothing else changed. No thresholds, no flag rules, no filtering. This
    script still modifies nothing and drops nothing.

WHAT IT ANSWERS
    - How many cells per slide, per animal, per group
    - How many cells per Phenotype, per animal, per group (counts and percents)
    - Are the column sets identical across all 8 files
    - What are the exact Phenotype label strings, including invisible characters
    - Does Parent encode animal ID only, or are there sub-regions (ROIs) per slide
    - Which acquisition (scan) each section came from, and where it sat on it
    - Coordinate extents, tissue bounding box, cell density
    - Per-marker intensity distributions per slide, and any failed / near-zero markers

OUTPUTS (written to OUT_DIR/tables/)
    00_audit_report.txt            full stdout mirror
    01_file_summary.csv            one row per file, now with scan and position
    02_column_audit.csv            one row per column, presence across files
    03_phenotype_labels_repr.csv   exact label strings with codepoint inspection
    04_phenotype_counts_long.csv   tidy counts: file x phenotype
    05_phenotype_counts_wide.csv   matrix of counts and percents
    06_marker_stats_long.csv       tidy per-file per-marker distribution stats
    07_marker_flags.csv            markers flagged as failed / saturated / suspect
    08_image_parent_values.csv     unique Image and Parent values per file
    09_scan_layout.csv             one row per file: scan, position, Y geometry

USAGE
    conda activate sc_pre
    python AKOYA_01_Inventory.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================
# ALL TUNABLE PARAMETERS LIVE HERE
# =============================================================================

DATA_DIR = "/master/jlehle/WORKING/AKOYA/data"
OUT_DIR = "/master/jlehle/WORKING/AKOYA/inventory"

# Group prefix in the filename -> treatment label.
# Confirmed with Deepak: G3 = D1MT-treated, G4 = untreated control.
GROUP_MAP = {
    "G3": "D1MT",
    "G4": "Untreated",
}

# Animals that also appear in the week-15 BAL scRNA-seq libraries.
# These are the cross-modality bridge animals.
SCRNA_OVERLAP_ANIMALS = ["43102", "43109"]

# Animals suspected to come from an earlier / different cohort based on ID range.
# Flagged only, not excluded. Confirm with Annu / Bindu.
SUSPECTED_OTHER_SERIES = ["31438", "36463"]

# Column-name conventions from the QuPath export
INTENSITY_SUFFIX = ": Mean"          # marker intensity columns end with this
CENTROID_X_PREFIX = "Centroid X"     # avoids hard-coding the micron symbol
CENTROID_Y_PREFIX = "Centroid Y"
PHENOTYPE_COL = "Phenotypes"
PARENT_COL = "Parent"
IMAGE_COL = "Image"

# ---- scan identity ----------------------------------------------------------
# The Image column is treated as the acquisition identifier. Each section file is
# expected to carry exactly one Image value. If a file carries more than one, the
# modal value is used as the scan and a warning is raised, because any downstream
# per-scan analysis is then only approximately right.
SCAN_ID_PREFIX = "scan"              # short stable label: scan_01, scan_02, ...
# Position rank 1 is the section with the smallest y_min on its scan. This is the
# same convention script 02 used when it identified the position-1 sections.
POSITION_RANK_ASCENDING_Y = True

# Thresholds for flagging suspect markers (diagnostic only, nothing is dropped)
ZERO_FRAC_FAIL = 0.98      # >= this fraction of exact zeros -> likely failed cycle
ZERO_FRAC_WARN = 0.90      # >= this fraction of exact zeros -> sparse, interpret with care
LOW_MAX_FAIL = 1.0         # max intensity below this -> essentially no signal
CV_LOW_WARN = 0.10         # coefficient of variation below this -> almost no dynamic range

# Percentiles reported for every marker
PERCENTILES = [1, 5, 25, 50, 75, 95, 99, 99.9]

# Encodings tried in order when reading a CSV
ENCODINGS = ["utf-8", "utf-8-sig", "latin-1"]

# Minimum cells per phenotype per animal that we would consider workable for
# nearest-neighbour statistics downstream. Reporting only, nothing is filtered.
MIN_CELLS_FOR_NN = 50

WRITE_TABLES = True


# %% Cell 2 - imports and helpers
# =============================================================================

import os
import sys
import glob
import unicodedata
from datetime import datetime

import numpy as np
import pandas as pd


class Tee(object):
    """Mirror stdout to a log file so the console output is reproducible."""

    def __init__(self, path):
        self.terminal = sys.stdout
        self.log = open(path, "w", encoding="utf-8", newline="\n")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        self.terminal.flush()
        self.log.flush()

    def close(self):
        try:
            self.log.close()
        except Exception:
            pass


def banner(text, char="=", width=88):
    print("\n" + char * width)
    print(text)
    print(char * width)


def sub(text, char="-", width=88):
    print("\n" + text)
    print(char * len(text) if len(text) < width else char * width)


def safe_read_csv(path):
    """Try a few encodings. Warn and continue rather than dying on one bad file."""
    last_err = None
    for enc in ENCODINGS:
        try:
            df = pd.read_csv(path, encoding=enc, low_memory=False)
            if enc != ENCODINGS[0]:
                print(f"    WARNING: {os.path.basename(path)} read with encoding "
                      f"'{enc}' rather than '{ENCODINGS[0]}'")
            return df, enc, None
        except Exception as e:
            last_err = e
            continue
    return None, None, last_err


def parse_filename(path):
    """
    Expect names like G3_43102.csv -> group prefix G3, animal 43102.
    Returns dict. Warns and falls back to the raw stem if the pattern does not hold.
    """
    stem = os.path.splitext(os.path.basename(path))[0]
    group_prefix, animal_id = None, None
    if "_" in stem:
        parts = stem.split("_")
        if len(parts) == 2:
            group_prefix, animal_id = parts[0], parts[1]
    if group_prefix is None or animal_id is None:
        print(f"    WARNING: filename '{stem}' does not match the expected "
              f"<GROUP>_<ANIMAL> pattern. Using stem as sample_id.")
        group_prefix = "UNKNOWN"
        animal_id = stem
    treatment = GROUP_MAP.get(group_prefix, "UNMAPPED")
    if treatment == "UNMAPPED":
        print(f"    WARNING: group prefix '{group_prefix}' is not in GROUP_MAP. "
              f"Treatment set to UNMAPPED.")
    return {
        "sample_id": stem,
        "group_prefix": group_prefix,
        "animal_id": animal_id,
        "treatment": treatment,
    }


def describe_string(s):
    """
    Return an inspection record for a label string.
    Catches zero-width spaces, non-breaking spaces, trailing whitespace, and any
    other non-ASCII character that will silently break string matching later.
    """
    if not isinstance(s, str):
        return {
            "label_repr": repr(s),
            "n_chars": np.nan,
            "has_nonascii": False,
            "nonascii_codepoints": "",
            "nonascii_names": "",
            "has_leading_trailing_ws": False,
            "ascii_safe_label": "",
        }
    nonascii = [c for c in s if ord(c) > 127]
    codepoints = ";".join(f"U+{ord(c):04X}" for c in nonascii)
    names = ";".join(unicodedata.name(c, "UNNAMED") for c in nonascii)
    ascii_safe = "".join(c for c in s if ord(c) <= 127).strip()
    return {
        "label_repr": repr(s),
        "n_chars": len(s),
        "has_nonascii": len(nonascii) > 0,
        "nonascii_codepoints": codepoints,
        "nonascii_names": names,
        "has_leading_trailing_ws": s != s.strip(),
        "ascii_safe_label": ascii_safe,
    }


def split_marker_column(col):
    """
    'CD68: Membrane: Mean' -> ('CD68', 'Membrane')
    'Granzyme-B: Cell: Mean' -> ('Granzyme-B', 'Cell')
    Returns (None, None) if the column does not parse.
    """
    if not col.endswith(INTENSITY_SUFFIX):
        return None, None
    body = col[: -len(INTENSITY_SUFFIX)]
    parts = [p.strip() for p in body.split(":")]
    if len(parts) < 2:
        return None, None
    marker = ":".join(parts[:-1]).strip()
    compartment = parts[-1].strip()
    return marker, compartment


def find_col(df, prefix):
    """Find a column by prefix so we never have to type the micron symbol."""
    hits = [c for c in df.columns if c.startswith(prefix)]
    if len(hits) == 0:
        return None
    if len(hits) > 1:
        print(f"    WARNING: multiple columns start with '{prefix}': {hits}. "
              f"Using the first.")
    return hits[0]


def write_csv(df, path):
    """Always LF line endings so downstream bash never chokes on CRLF."""
    if not WRITE_TABLES:
        return
    df.to_csv(path, index=False, lineterminator="\n", encoding="utf-8")
    print(f"    wrote {path}  ({df.shape[0]} rows x {df.shape[1]} cols)")


# %% Cell 3 - discover files and set up output directory
# =============================================================================

os.makedirs(OUT_DIR, exist_ok=True)
TABLE_DIR = os.path.join(OUT_DIR, "tables")
os.makedirs(TABLE_DIR, exist_ok=True)

_tee = Tee(os.path.join(TABLE_DIR, "00_audit_report.txt"))
sys.stdout = _tee

banner("AKOYA INVENTORY / SCHEMA AUDIT (revision 2)")
print(f"Run time   : {datetime.now().isoformat(timespec='seconds')}")
print(f"Data dir   : {DATA_DIR}")
print(f"Output dir : {TABLE_DIR}")
print(f"Group map  : {GROUP_MAP}")
print(f"Scan key   : the '{IMAGE_COL}' column, modal value per file")

csv_paths = sorted(glob.glob(os.path.join(DATA_DIR, "*.csv")))
if len(csv_paths) == 0:
    print(f"\nERROR: no CSV files found in {DATA_DIR}. Nothing to do.")
    sys.stdout = _tee.terminal
    _tee.close()
    sys.exit(1)

print(f"\nFound {len(csv_paths)} CSV file(s):")
file_meta = {}
for p in csv_paths:
    meta = parse_filename(p)
    meta["path"] = p
    meta["file_size_mb"] = round(os.path.getsize(p) / 1e6, 2)
    meta["in_scrnaseq"] = meta["animal_id"] in SCRNA_OVERLAP_ANIMALS
    meta["suspected_other_series"] = meta["animal_id"] in SUSPECTED_OTHER_SERIES
    file_meta[meta["sample_id"]] = meta
    tags = []
    if meta["in_scrnaseq"]:
        tags.append("also in week-15 scRNA-seq")
    if meta["suspected_other_series"]:
        tags.append("SUSPECTED different cohort")
    tag_str = f"   [{'; '.join(tags)}]" if tags else ""
    print(f"  {meta['sample_id']:<14} group={meta['group_prefix']:<4} "
          f"animal={meta['animal_id']:<8} treatment={meta['treatment']:<10} "
          f"{meta['file_size_mb']:>7.2f} MB{tag_str}")


# %% Cell 4 - load every file and audit structure
# =============================================================================

banner("PER-FILE STRUCTURE AUDIT")

frames = {}
file_rows = []
image_parent_rows = []
failed_files = []

for sid in sorted(file_meta.keys()):
    meta = file_meta[sid]
    print(f"\n[{sid}]")
    df, enc, err = safe_read_csv(meta["path"])
    if df is None:
        print(f"    ERROR: could not read file. Last error: {err}")
        print(f"    Skipping and continuing.")
        failed_files.append(sid)
        continue

    frames[sid] = df
    n_rows, n_cols = df.shape
    print(f"    rows (cells) : {n_rows:,}")
    print(f"    columns      : {n_cols}")
    print(f"    encoding     : {enc}")

    # --- duplicates and missingness -----------------------------------------
    n_dup = int(df.duplicated().sum())
    n_any_na = int(df.isna().any(axis=1).sum())
    if n_dup > 0:
        print(f"    WARNING: {n_dup:,} fully duplicated rows")
    if n_any_na > 0:
        print(f"    NOTE: {n_any_na:,} rows contain at least one NaN")

    # --- Image and Parent ---------------------------------------------------
    # The modal Image value becomes this section's scan. A file carrying more
    # than one Image value is reported loudly, because every per-scan analysis
    # downstream assumes one acquisition per section file.
    img_vals, par_vals = [], []
    scan_image_value = "UNKNOWN"
    modal_image_pct = np.nan
    if IMAGE_COL in df.columns:
        img_counts = df[IMAGE_COL].dropna().value_counts()
        img_vals = sorted(img_counts.index.tolist())
        print(f"    unique Image values  : {len(img_vals)}")
        for v in img_vals[:10]:
            print(f"        {v}")
        if len(img_vals) > 10:
            print(f"        ... and {len(img_vals) - 10} more")
        if len(img_counts) > 0:
            scan_image_value = str(img_counts.index[0])
            modal_image_pct = 100.0 * float(img_counts.iloc[0]) / n_rows
            print(f"    scan (modal Image)   : {scan_image_value}")
            if len(img_counts) > 1:
                print(f"    WARNING: {len(img_counts)} distinct Image values in one "
                      f"section file. The modal value covers "
                      f"{modal_image_pct:.1f}% of cells and is being used as the "
                      f"scan. Any per-scan analysis is approximate for this file.")
        else:
            print(f"    WARNING: '{IMAGE_COL}' column present but entirely NaN")
    else:
        print(f"    WARNING: no '{IMAGE_COL}' column. Scan set to UNKNOWN, which "
              f"will collapse every file onto one scan downstream.")

    if PARENT_COL in df.columns:
        par_counts = df[PARENT_COL].fillna("<NA>").value_counts()
        par_vals = par_counts.index.tolist()
        print(f"    unique Parent values : {len(par_vals)}")
        for v, c in par_counts.head(15).items():
            print(f"        {v:<40} {c:>10,} cells")
        if len(par_vals) > 15:
            print(f"        ... and {len(par_vals) - 15} more")
        if len(par_vals) == 1:
            print(f"    -> Parent is a single value for this slide, consistent with "
                  f"it encoding animal/slide identity only (no ROI substructure).")
        else:
            print(f"    -> Parent has MULTIPLE values. Check whether these are "
                  f"separate ROIs / annotation regions rather than animal ID.")
    else:
        print(f"    WARNING: no '{PARENT_COL}' column")

    for v in img_vals:
        image_parent_rows.append({"sample_id": sid, "field": IMAGE_COL, "value": v,
                                  "n_cells": int((df[IMAGE_COL] == v).sum())})
    for v in par_vals:
        image_parent_rows.append({"sample_id": sid, "field": PARENT_COL, "value": v,
                                  "n_cells": int((df[PARENT_COL].fillna("<NA>") == v).sum())})

    # --- coordinates --------------------------------------------------------
    xcol = find_col(df, CENTROID_X_PREFIX)
    ycol = find_col(df, CENTROID_Y_PREFIX)
    x_min = x_max = y_min = y_max = np.nan
    bbox_area_mm2 = density = np.nan
    if xcol is None or ycol is None:
        print(f"    WARNING: could not find centroid columns "
              f"(x='{xcol}', y='{ycol}'). Spatial analysis will not be possible.")
    else:
        x = pd.to_numeric(df[xcol], errors="coerce")
        y = pd.to_numeric(df[ycol], errors="coerce")
        x_min, x_max = float(np.nanmin(x)), float(np.nanmax(x))
        y_min, y_max = float(np.nanmin(y)), float(np.nanmax(y))
        span_x, span_y = x_max - x_min, y_max - y_min
        bbox_area_mm2 = (span_x * span_y) / 1e6
        density = n_rows / bbox_area_mm2 if bbox_area_mm2 > 0 else np.nan
        print(f"    X range      : {x_min:>12,.1f} to {x_max:>12,.1f} um  "
              f"(span {span_x:,.1f})")
        print(f"    Y range      : {y_min:>12,.1f} to {y_max:>12,.1f} um  "
              f"(span {span_y:,.1f})")
        print(f"    bbox area    : {bbox_area_mm2:,.2f} mm^2")
        print(f"    cell density : {density:,.0f} cells / mm^2 (bounding box, not tissue)")
        if span_y > 0 and (span_x / span_y > 50 or span_y / span_x > 50):
            print(f"    WARNING: extreme aspect ratio ({span_x / span_y:.1f}). "
                  f"Check the coordinate units and whether both axes are in microns.")
        n_bad_xy = int((x.isna() | y.isna()).sum())
        if n_bad_xy > 0:
            print(f"    WARNING: {n_bad_xy:,} rows have missing or non-numeric coordinates")

    # --- phenotypes ---------------------------------------------------------
    n_phen = np.nan
    if PHENOTYPE_COL in df.columns:
        phen_counts = df[PHENOTYPE_COL].fillna("<NA>").value_counts()
        n_phen = int(len(phen_counts))
        print(f"    unique Phenotypes : {n_phen}")
        for lab, c in phen_counts.items():
            pct = 100.0 * c / n_rows
            print(f"        {repr(lab):<46} {c:>10,}  ({pct:5.2f}%)")
    else:
        print(f"    WARNING: no '{PHENOTYPE_COL}' column")

    file_rows.append({
        "sample_id": sid,
        "group_prefix": meta["group_prefix"],
        "animal_id": meta["animal_id"],
        "treatment": meta["treatment"],
        "scan_image_value": scan_image_value,
        "modal_image_pct_of_cells": modal_image_pct,
        "in_scrnaseq": meta["in_scrnaseq"],
        "suspected_other_series": meta["suspected_other_series"],
        "file_size_mb": meta["file_size_mb"],
        "encoding": enc,
        "n_cells": n_rows,
        "n_columns": n_cols,
        "n_duplicate_rows": n_dup,
        "n_rows_with_na": n_any_na,
        "n_unique_image": len(img_vals),
        "n_unique_parent": len(par_vals),
        "n_unique_phenotypes": n_phen,
        "x_min_um": x_min, "x_max_um": x_max,
        "y_min_um": y_min, "y_max_um": y_max,
        "bbox_area_mm2": bbox_area_mm2,
        "cells_per_mm2_bbox": density,
    })

file_summary = pd.DataFrame(file_rows)
if failed_files:
    print(f"\nWARNING: {len(failed_files)} file(s) failed to load: {failed_files}")

if len(frames) == 0:
    print("\nERROR: no files loaded successfully. Stopping.")
    sys.stdout = _tee.terminal
    _tee.close()
    sys.exit(1)


# %% Cell 4b - scan identity and section position on the scan
# =============================================================================
# These are inventory facts, not analysis. They are derived here once so that no
# downstream script has to reconstruct them, and so that "slide" never again has
# to be approximated by "treatment".

banner("SCAN IDENTITY AND SECTION POSITION")

# ---- short stable scan ids --------------------------------------------------
scan_values = sorted(file_summary["scan_image_value"].unique().tolist())
scan_id_map = {v: f"{SCAN_ID_PREFIX}_{i:02d}"
               for i, v in enumerate(scan_values, start=1)}
file_summary["scan_id"] = file_summary["scan_image_value"].map(scan_id_map)

print(f"    {len(scan_values)} distinct scan(s) across {len(file_summary)} section file(s)")
for v in scan_values:
    members = file_summary.loc[file_summary["scan_image_value"] == v, "sample_id"].tolist()
    print(f"      {scan_id_map[v]}  n={len(members)}  {v}")
    print(f"                 sections: {', '.join(sorted(members))}")

if len(scan_values) == 1:
    print("\n    WARNING: every file resolved to a single scan. A between-scan "
          "versus within-scan variance decomposition is not possible.")

# ---- geometry within each scan ----------------------------------------------
file_summary["y_span_um"] = file_summary["y_max_um"] - file_summary["y_min_um"]
file_summary["x_span_um"] = file_summary["x_max_um"] - file_summary["x_min_um"]
file_summary["slide_position_rank"] = np.nan
file_summary["n_sections_on_scan"] = np.nan
file_summary["y_offset_from_scan_min_um"] = np.nan

for scan, g in file_summary.groupby("scan_id"):
    file_summary.loc[g.index, "n_sections_on_scan"] = len(g)
    have_y = g.loc[g["y_min_um"].notna()]
    if not len(have_y):
        print(f"\n    WARNING: {scan} has no usable Y coordinates. Position rank "
              f"left as NaN for its sections.")
        continue
    ordered = have_y.sort_values("y_min_um", ascending=POSITION_RANK_ASCENDING_Y)
    for rank, idx in enumerate(ordered.index, start=1):
        file_summary.loc[idx, "slide_position_rank"] = rank
    file_summary.loc[have_y.index, "y_offset_from_scan_min_um"] = (
        have_y["y_min_um"] - float(have_y["y_min_um"].min()))

# ---- report -----------------------------------------------------------------
sub("Section layout on each scan (rank 1 = smallest y_min)")
print(f"    {'sample_id':<14}{'scan':<10}{'treatment':<12}{'pos':>5}"
      f"{'y_min':>12}{'y_max':>12}{'y_span':>12}{'y_offset':>12}")
print("    " + "-" * 89)
layout = file_summary.sort_values(["scan_id", "slide_position_rank"])
for _, r in layout.iterrows():
    pos = "na" if pd.isna(r["slide_position_rank"]) else f"{int(r['slide_position_rank'])}"
    print(f"    {r['sample_id']:<14}{str(r['scan_id']):<10}{r['treatment']:<12}"
          f"{pos:>5}{r['y_min_um']:>12,.1f}{r['y_max_um']:>12,.1f}"
          f"{r['y_span_um']:>12,.1f}{r['y_offset_from_scan_min_um']:>12,.1f}")

print("\n    y_offset_from_scan_min_um near zero means the section sits flush")
print("    against the bottom edge of its scan. That, together with a short")
print("    y_span, is the geometric part of the position-1 argument. Whether a")
print("    section is also globally dim is a marker question and is answered in")
print("    script 02, not here.")

# ---- is scan confounded with treatment? -------------------------------------
sub("Is scan confounded with treatment?")
ct = pd.crosstab(file_summary["scan_id"], file_summary["treatment"])
print(ct.to_string())
scans_per_treatment = file_summary.groupby("treatment")["scan_id"].nunique()
treatments_per_scan = file_summary.groupby("scan_id")["treatment"].nunique()
fully_confounded = bool((treatments_per_scan == 1).all()
                        and (scans_per_treatment == 1).all())
if fully_confounded:
    print("\n    CONFOUNDED. Each scan carries exactly one treatment arm and each")
    print("    arm sits on exactly one scan. Between-scan variance and between-arm")
    print("    variance are the same quantity in this design, so any decomposition")
    print("    that follows measures batch and biology together and must be")
    print("    described that way. It cannot isolate batch.")
elif bool((treatments_per_scan == 1).all()):
    print("\n    PARTIALLY CONFOUNDED. Every scan is single-arm, but at least one")
    print("    arm spans more than one scan, so within-arm between-scan variance")
    print("    is estimable.")
else:
    print("\n    NOT CONFOUNDED. At least one scan carries both arms, so batch and")
    print("    biology are separable.")

scan_layout = file_summary[[
    "sample_id", "animal_id", "treatment", "scan_id", "scan_image_value",
    "n_sections_on_scan", "slide_position_rank",
    "x_min_um", "x_max_um", "x_span_um",
    "y_min_um", "y_max_um", "y_span_um", "y_offset_from_scan_min_um",
    "n_cells", "bbox_area_mm2", "cells_per_mm2_bbox",
]].sort_values(["scan_id", "slide_position_rank"])

# lookups reused by the long tables below
SCAN_ID_OF = dict(zip(file_summary["sample_id"], file_summary["scan_id"]))
SCAN_IMAGE_OF = dict(zip(file_summary["sample_id"], file_summary["scan_image_value"]))
POSITION_OF = dict(zip(file_summary["sample_id"], file_summary["slide_position_rank"]))


# %% Cell 5 - cross-file column comparison
# =============================================================================

banner("CROSS-FILE COLUMN COMPARISON")

all_cols = []
for sid, df in frames.items():
    all_cols.extend(df.columns.tolist())
col_universe = sorted(set(all_cols))

col_rows = []
for col in col_universe:
    present = [sid for sid in sorted(frames) if col in frames[sid].columns]
    marker, compartment = split_marker_column(col)
    rec = {
        "column": col,
        "is_intensity": marker is not None,
        "marker": marker if marker is not None else "",
        "compartment": compartment if compartment is not None else "",
        "n_files_present": len(present),
        "present_in_all": len(present) == len(frames),
        "missing_from": ";".join(sorted(set(frames) - set(present))),
    }
    for sid in sorted(frames):
        rec[f"in_{sid}"] = col in frames[sid].columns
    col_rows.append(rec)

col_audit = pd.DataFrame(col_rows)

n_shared = int(col_audit["present_in_all"].sum())
print(f"Union of all column names across files : {len(col_universe)}")
print(f"Columns present in ALL {len(frames)} files     : {n_shared}")

if n_shared != len(col_universe):
    print("\nWARNING: column sets are NOT identical across files.")
    bad = col_audit.loc[~col_audit["present_in_all"]]
    for _, r in bad.iterrows():
        print(f"    {r['column']:<52} present in {r['n_files_present']}/{len(frames)}"
              f"  missing from: {r['missing_from']}")
else:
    print("All files share an identical column set.")

intensity_cols = col_audit.loc[col_audit["is_intensity"], "column"].tolist()
meta_cols = col_audit.loc[~col_audit["is_intensity"], "column"].tolist()

print(f"\nIntensity columns (end with '{INTENSITY_SUFFIX}') : {len(intensity_cols)}")
print(f"Non-intensity / metadata columns                : {len(meta_cols)}")
for c in meta_cols:
    print(f"    {c}")

# markers measured in more than one compartment
comp_tbl = (col_audit.loc[col_audit["is_intensity"]]
            .groupby("marker")["compartment"]
            .apply(lambda s: sorted(set(s)))
            .reset_index(name="compartments"))
multi = comp_tbl.loc[comp_tbl["compartments"].apply(len) > 1]
if len(multi) > 0:
    print(f"\nMarkers quantified in more than one compartment ({len(multi)}):")
    for _, r in multi.iterrows():
        print(f"    {r['marker']:<28} {r['compartments']}")

print(f"\nCompartment usage across the panel:")
for comp, n in col_audit.loc[col_audit["is_intensity"], "compartment"].value_counts().items():
    print(f"    {comp:<16} {n:>3} markers")


# %% Cell 6 - phenotype label inspection and counts
# =============================================================================

banner("PHENOTYPE LABELS - EXACT STRING INSPECTION")

label_universe = set()
for sid, df in frames.items():
    if PHENOTYPE_COL in df.columns:
        label_universe.update(df[PHENOTYPE_COL].dropna().unique().tolist())

label_rows = []
for lab in sorted(label_universe, key=lambda s: str(s)):
    rec = describe_string(lab)
    rec["label"] = lab
    label_rows.append(rec)
label_repr_tbl = pd.DataFrame(label_rows)

print(f"Distinct Phenotype labels across all files: {len(label_repr_tbl)}\n")
for _, r in label_repr_tbl.iterrows():
    flag = ""
    if r["has_nonascii"]:
        flag += f"  <-- NON-ASCII: {r['nonascii_codepoints']} ({r['nonascii_names']})"
    if r["has_leading_trailing_ws"]:
        flag += "  <-- LEADING/TRAILING WHITESPACE"
    print(f"    {r['label_repr']:<50} n_chars={int(r['n_chars']):>3}{flag}")

if label_repr_tbl["has_nonascii"].any() or label_repr_tbl["has_leading_trailing_ws"].any():
    print("\n    IMPORTANT: at least one label contains characters that will break "
          "naive string matching.")
    print("    Downstream scripts must match on the exact string, or map through the "
          "'ascii_safe_label' column in 03_phenotype_labels_repr.csv.")

# ---- counts -----------------------------------------------------------------
banner("PHENOTYPE COUNTS BY ANIMAL AND TREATMENT")

count_rows = []
for sid in sorted(frames):
    df = frames[sid]
    meta = file_meta[sid]
    if PHENOTYPE_COL not in df.columns:
        continue
    total = len(df)
    vc = df[PHENOTYPE_COL].fillna("<NA>").value_counts()
    for lab, c in vc.items():
        count_rows.append({
            "sample_id": sid,
            "group_prefix": meta["group_prefix"],
            "animal_id": meta["animal_id"],
            "treatment": meta["treatment"],
            "scan_id": SCAN_ID_OF.get(sid, "UNKNOWN"),
            "scan_image_value": SCAN_IMAGE_OF.get(sid, "UNKNOWN"),
            "slide_position_rank": POSITION_OF.get(sid, np.nan),
            "phenotype": lab,
            "n_cells": int(c),
            "pct_of_slide": 100.0 * c / total,
            "slide_total_cells": total,
            "meets_min_for_nn": bool(c >= MIN_CELLS_FOR_NN),
        })

pheno_long = pd.DataFrame(count_rows)

if len(pheno_long) > 0:
    wide_n = pheno_long.pivot_table(index="phenotype", columns="sample_id",
                                    values="n_cells", fill_value=0)
    wide_pct = pheno_long.pivot_table(index="phenotype", columns="sample_id",
                                      values="pct_of_slide", fill_value=0.0)

    print("\nCell counts (rows = phenotype, columns = slide):\n")
    with pd.option_context("display.width", 200, "display.max_columns", 50):
        print(wide_n.to_string())

    print("\nPercent of slide:\n")
    with pd.option_context("display.width", 200, "display.max_columns", 50,
                           "display.float_format", lambda v: f"{v:6.2f}"):
        print(wide_pct.to_string())

    # group-level rollup
    grp = (pheno_long.groupby(["treatment", "phenotype"])["n_cells"]
           .sum().reset_index())
    grp_tot = grp.groupby("treatment")["n_cells"].transform("sum")
    grp["pct_of_group"] = 100.0 * grp["n_cells"] / grp_tot
    print("\nPooled by treatment group (pooling is for inventory only, NOT for stats):\n")
    with pd.option_context("display.float_format", lambda v: f"{v:8.2f}"):
        print(grp.pivot_table(index="phenotype", columns="treatment",
                              values=["n_cells", "pct_of_group"],
                              fill_value=0).to_string())

    # power check for nearest-neighbour work
    sub("Per-animal power check for nearest-neighbour statistics")
    print(f"Phenotypes with fewer than {MIN_CELLS_FOR_NN} cells in at least one animal:")
    thin = pheno_long.loc[~pheno_long["meets_min_for_nn"]]
    if len(thin) == 0:
        print(f"    none - every phenotype clears {MIN_CELLS_FOR_NN} cells in every animal")
    else:
        for _, r in thin.sort_values(["phenotype", "sample_id"]).iterrows():
            print(f"    {r['phenotype']:<44} {r['sample_id']:<12} "
                  f"{r['n_cells']:>7,} cells")

    # phenotypes not present in every slide at all
    seen = pheno_long.groupby("phenotype")["sample_id"].nunique()
    absent = seen.loc[seen < len(frames)]
    if len(absent) > 0:
        print(f"\nPhenotypes absent from at least one slide entirely:")
        for lab, n in absent.items():
            have = sorted(pheno_long.loc[pheno_long["phenotype"] == lab, "sample_id"])
            print(f"    {lab:<44} present in {n}/{len(frames)}: {have}")
else:
    pheno_long = pd.DataFrame()
    wide_n = pd.DataFrame()
    wide_pct = pd.DataFrame()
    print("WARNING: no phenotype counts could be computed.")


# %% Cell 7 - per-marker intensity distributions
# =============================================================================

banner("PER-MARKER INTENSITY DISTRIBUTIONS")

marker_rows = []
for sid in sorted(frames):
    df = frames[sid]
    meta = file_meta[sid]
    cols = [c for c in intensity_cols if c in df.columns]
    for col in cols:
        marker, compartment = split_marker_column(col)
        v = pd.to_numeric(df[col], errors="coerce")
        n_valid = int(v.notna().sum())
        if n_valid == 0:
            marker_rows.append({
                "sample_id": sid, "treatment": meta["treatment"],
                "animal_id": meta["animal_id"],
                "scan_id": SCAN_ID_OF.get(sid, "UNKNOWN"),
                "scan_image_value": SCAN_IMAGE_OF.get(sid, "UNKNOWN"),
                "slide_position_rank": POSITION_OF.get(sid, np.nan),
                "column": col,
                "marker": marker, "compartment": compartment,
                "n_valid": 0, "n_na": int(v.isna().sum()),
                "frac_zero": np.nan, "mean": np.nan, "sd": np.nan, "cv": np.nan,
                "min": np.nan, "max": np.nan,
                **{f"p{p}": np.nan for p in PERCENTILES},
            })
            continue
        vv = v.dropna().to_numpy()
        mean = float(np.mean(vv))
        sd = float(np.std(vv, ddof=1)) if len(vv) > 1 else 0.0
        rec = {
            "sample_id": sid,
            "treatment": meta["treatment"],
            "animal_id": meta["animal_id"],
            "scan_id": SCAN_ID_OF.get(sid, "UNKNOWN"),
            "scan_image_value": SCAN_IMAGE_OF.get(sid, "UNKNOWN"),
            "slide_position_rank": POSITION_OF.get(sid, np.nan),
            "column": col,
            "marker": marker,
            "compartment": compartment,
            "n_valid": n_valid,
            "n_na": int(v.isna().sum()),
            "frac_zero": float(np.mean(vv == 0)),
            "mean": mean,
            "sd": sd,
            "cv": (sd / mean) if mean > 0 else np.nan,
            "min": float(np.min(vv)),
            "max": float(np.max(vv)),
        }
        qs = np.percentile(vv, PERCENTILES)
        for p, q in zip(PERCENTILES, qs):
            rec[f"p{p}"] = float(q)
        marker_rows.append(rec)

marker_stats = pd.DataFrame(marker_rows)
print(f"Computed distribution stats for {marker_stats['column'].nunique()} "
      f"intensity columns across {marker_stats['sample_id'].nunique()} slides "
      f"({len(marker_stats):,} rows).")


# ---- flags ------------------------------------------------------------------
def flag_row(r):
    flags = []
    if pd.isna(r["max"]):
        flags.append("ALL_NA")
        return ";".join(flags)
    if r["max"] < LOW_MAX_FAIL:
        flags.append("NO_SIGNAL")
    if r["frac_zero"] >= ZERO_FRAC_FAIL:
        flags.append("LIKELY_FAILED_CYCLE")
    elif r["frac_zero"] >= ZERO_FRAC_WARN:
        flags.append("VERY_SPARSE")
    if not pd.isna(r["cv"]) and r["cv"] < CV_LOW_WARN:
        flags.append("LOW_DYNAMIC_RANGE")
    return ";".join(flags)


marker_stats["flags"] = marker_stats.apply(flag_row, axis=1)
flagged = marker_stats.loc[marker_stats["flags"] != ""].copy()

sub("Flagged marker/slide combinations")
if len(flagged) == 0:
    print("    none - every marker has signal on every slide")
else:
    print(f"    {len(flagged)} marker x slide combination(s) flagged\n")
    for _, r in flagged.sort_values(["marker", "sample_id"]).iterrows():
        print(f"    {r['marker']:<26} {r['compartment']:<12} {r['sample_id']:<12} "
              f"max={r['max']:>10.2f}  zero_frac={r['frac_zero']:>5.3f}  "
              f"[{r['flags']}]")

# markers flagged on some slides but not others -> batch problem, not a dead antibody
if len(flagged) > 0:
    per_marker = (marker_stats.assign(is_flagged=marker_stats["flags"] != "")
                  .groupby("column")["is_flagged"].agg(["sum", "count"]))
    partial = per_marker.loc[(per_marker["sum"] > 0) &
                             (per_marker["sum"] < per_marker["count"])]
    if len(partial) > 0:
        sub("Markers flagged on SOME slides but not others (batch / cycle problem)")
        for col, r in partial.iterrows():
            bad = sorted(marker_stats.loc[(marker_stats["column"] == col) &
                                          (marker_stats["flags"] != ""), "sample_id"])
            print(f"    {col:<52} flagged in {int(r['sum'])}/{int(r['count'])}: {bad}")

# ---- dynamic range spread across slides (staining / batch drift) -------------
sub("Cross-slide spread of the 99th percentile (a proxy for staining intensity drift)")
spread = (marker_stats.pivot_table(index="column", columns="sample_id", values="p99")
          .assign(min_p99=lambda d: d.min(axis=1),
                  max_p99=lambda d: d.max(axis=1)))
spread["fold_range"] = spread["max_p99"] / spread["min_p99"].replace(0, np.nan)
worst = spread.sort_values("fold_range", ascending=False).head(25)
print(f"Top 25 markers by max/min ratio of p99 across slides:\n")
for col, r in worst.iterrows():
    print(f"    {col:<52} p99 {r['min_p99']:>9.2f} to {r['max_p99']:>10.2f}   "
          f"fold={r['fold_range']:>8.1f}")
print("\n    Large fold ranges mean a shared intensity threshold across slides would be "
      "unsafe. Per-slide normalization will need to be settled before any gating.")
print("    This spread mixes scan and section. Script 02 splits it into the")
print("    between-scan and within-scan parts, keyed on scan_id from this run.")


# %% Cell 8 - write outputs and summarize
# =============================================================================

banner("WRITING TABLES")

write_csv(file_summary, os.path.join(TABLE_DIR, "01_file_summary.csv"))
write_csv(col_audit, os.path.join(TABLE_DIR, "02_column_audit.csv"))
write_csv(label_repr_tbl, os.path.join(TABLE_DIR, "03_phenotype_labels_repr.csv"))
if len(pheno_long) > 0:
    write_csv(pheno_long, os.path.join(TABLE_DIR, "04_phenotype_counts_long.csv"))
    wide_out = wide_n.reset_index()
    pct_out = wide_pct.reset_index()
    pct_out.columns = ["phenotype"] + [f"{c}_pct" for c in pct_out.columns[1:]]
    write_csv(wide_out.merge(pct_out, on="phenotype"),
              os.path.join(TABLE_DIR, "05_phenotype_counts_wide.csv"))
write_csv(marker_stats, os.path.join(TABLE_DIR, "06_marker_stats_long.csv"))
write_csv(flagged if len(flagged) > 0 else marker_stats.head(0),
          os.path.join(TABLE_DIR, "07_marker_flags.csv"))
write_csv(pd.DataFrame(image_parent_rows),
          os.path.join(TABLE_DIR, "08_image_parent_values.csv"))
write_csv(scan_layout, os.path.join(TABLE_DIR, "09_scan_layout.csv"))

banner("SUMMARY")
print(f"Files loaded              : {len(frames)} / {len(csv_paths)}")
if failed_files:
    print(f"Files failed              : {failed_files}")
print(f"Total cells               : {file_summary['n_cells'].sum():,}")
for trt in sorted(file_summary["treatment"].unique()):
    sl = file_summary.loc[file_summary["treatment"] == trt]
    print(f"  {trt:<12} {len(sl)} slide(s), {sl['n_cells'].sum():>10,} cells "
          f"(range {sl['n_cells'].min():,} to {sl['n_cells'].max():,})")
print(f"Distinct scans            : {len(scan_values)}")
for scan, g in file_summary.groupby("scan_id"):
    arms = sorted(g["treatment"].unique())
    print(f"  {scan:<12} {len(g)} section(s), arms: {', '.join(arms)}")
print(f"Scan/treatment confounded : {fully_confounded}")
print(f"Column set identical      : {n_shared == len(col_universe)}")
print(f"Intensity columns         : {len(intensity_cols)}")
print(f"Distinct phenotype labels : {len(label_repr_tbl)}")
print(f"Flagged marker x slide    : {len(flagged)}")

sub("Open questions this run should have answered")
print("  1. Does Parent hold one value per file (animal ID) or several (ROIs)?")
print("  2. Are the Phenotype labels the final 13-phenotype call, or something coarser?")
print("  3. Which phenotypes clear the per-animal cell count needed for")
print("     nearest-neighbour statistics in every animal?")
print("  4. Which markers need per-slide normalization before any gating?")
print("  5. Does every section file resolve to exactly one scan, and does the")
print("     section layout on each scan match what was physically mounted?")

sub("What downstream now reads from here")
print("  01_file_summary.csv  scan_id, scan_image_value, slide_position_rank")
print("                       -> script 02 Q2 keys its variance decomposition on")
print("                          scan_id instead of treatment")
print("  09_scan_layout.csv   the geometric half of the position-1 exclusion")
print("                       argument, in one table")

print(f"\nAll tables written to: {TABLE_DIR}")
banner("DONE")

sys.stdout = _tee.terminal
_tee.close()
