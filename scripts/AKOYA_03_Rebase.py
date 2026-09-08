#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - RE-BASELINE ON SIX SECTIONS
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 03 of the AKOYA analysis series. REVISION 2.

WHY THIS SCRIPT EXISTS
    Scripts 01 and 02 established that the position-1 section on each scan is
    globally dim. G3_43102 (D1MT) and G4_43112 (Untreated) are excluded here,
    one from each arm, so the design stays balanced at 3 versus 3. The exclusion
    is declared once, in EXCLUDE_SECTIONS below, and every downstream script
    should read it from this file rather than re-deciding.

WHAT CHANGED IN REVISION 2 (and why)

    1. THE PLATEAU SEARCH IS CONSTRAINED BY ARM SEPARATION.
       The old search took the widest cutoff window over which every section's
       called fraction moved by less than PLATEAU_TOL. Near cutoff zero every
       section calls essentially 100 percent positive, so the curves are flat
       and identical there, and that degenerate window can win. It has an equal
       or greater width than the informative one and the loop keeps whichever it
       reaches first, which is always the low one.

       The search now runs in two steps. First it finds the SEPARATION window:
       the contiguous range of cutoffs where the arms remain completely
       separated by at least MIN_ARM_GAP in called fraction. That range is the
       actual claim, that the group difference survives any gate in it. Second,
       inside that range only, it finds the widest STABILITY window where every
       section moves less than PLATEAU_TOL. A window where every section calls
       100 percent positive can never be admissible, because the gap there is
       zero.

       Script 09 must use the identical function and the identical parameters.
       They are grouped under one header in Cell 1 for that reason, and the
       script prints every input to the search so the two runs can be diffed
       directly. Revision 1 of scripts 03 and 09 implemented the same algorithm
       on what should have been the same data and reported different windows
       (6.10 to 9.67 against 0.00 to 0.50). That divergence has never been
       explained, so this script now also prints the per-section cell counts and
       the called fraction at cutoff zero, which is where any input difference
       will show.

    2. THE VARIANCE DECOMPOSITION IS RECONCILED WITH SCRIPT 02.
       Cell 4 previously grouped by `condition`, used SECTION-WIDE MEDIANS, ran
       on six sections, and reported the naive between / (between + within)
       ratio. Script 02 groups by scan, uses p99, runs on eight, and now reports
       a corrected ICC. Four differences at once, which is why the numbers
       carried in the project notes cannot be reproduced from either script.

       This cell now computes the decomposition on BOTH statistics (p99 and
       section median), keyed on scan_id, reporting the naive ratio and the
       corrected ICC for each, all on the six retained sections. Four columns,
       one table, and the difference between the script 02 numbers and these
       becomes attributable rather than mysterious.

    3. THE EXCLUSION RATIONALE IS UPDATED TO THE COMPOSITION-CONTROLLED TEST.
       Script 02 revision 2 showed that section composition was inflating the
       original dimness numbers: 43102 falls from 68.7 percent of comparisons
       ranked dimmest on all cells to 48.8 percent when the comparison is made
       within phenotype, and 43112 from 73.1 to 56.3. Both remain roughly twice
       the 25 percent null and far outside every other section, so the decision
       is unchanged, but the reason strings now quote the defensible statistic.

       The claim that both sections have the shortest Y span on their slide is
       also removed. It is true on scan_01 (43102 at 2,797 um) and false on
       scan_02, where 36463 spans 2,575 um against 43112's 3,336 um. The flush
       Y offset, both at 0.7 um, is intact and is the geometric point.

    4. SCAN AND POSITION COME FROM THE INVENTORY.
       scan was parsed from df[Image].iloc[0]. It now reads scan_id and
       slide_position_rank from script 01, with the old parse as a fallback.
       An optional consistency check compares the EXCLUDE_SECTIONS list against
       script 02's table 12 and warns if the evidence has moved.

    5. SMALL THINGS.
       Boundary hole filling is now an explicit parameter with a comment, since
       filling holes treats airway and vessel lumens as tissue and this is lung.
       Palettes defined once. Figure and table ranges in the docstring corrected
       to what the script actually writes.

WHAT THIS SCRIPT DOES
    A  Re-baseline composition, IDO1 fraction and the batch decomposition on the
       six retained sections.
    B  Back out the IDO1 intensity threshold implied by the vendor's binary call.
    C  Threshold sensitivity with the separation-constrained window search.
    D  Rebuild the IDO1 distribution figures on six sections.
    E  Distance to the tissue boundary on CD68-lineage cells, DAPI as control.
    F  Spatial maps showing ALL phenotypes.

    Still diagnostic. Nothing is re-phenotyped or written back to source CSVs.

OUTPUTS
    figures/  F17 .. F23   (PDF + PNG, 300 DPI)
    tables/   20 .. 26     (CSV, LF line endings)

USAGE
    conda activate sc_pre
    python AKOYA_03_Rebase.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================
# ALL TUNABLE PARAMETERS LIVE HERE
# =============================================================================

DATA_DIR = "/master/jlehle/WORKING/AKOYA/data"
INVENTORY_TABLE_DIR = "/master/jlehle/WORKING/AKOYA/inventory/tables"
DIAGNOSTIC_TABLE_DIR = "/master/jlehle/WORKING/AKOYA/diagnostics/tables"
OUT_DIR = "/master/jlehle/WORKING/AKOYA/rebaseline"

# ---- THE EXCLUSION, DECLARED ONCE -------------------------------------------
# Downstream scripts should import or copy this block verbatim.
# Evidence: script 02 revision 2, table 12_slide_position_summary.csv.
# The statistic quoted is the composition-controlled rank statistic, which is
# the fraction of phenotype-by-marker comparisons in which this section is the
# dimmest of the four on its scan. Under no position effect that is 25 percent
# exactly, with no distributional assumption.
EXCLUDE_SECTIONS = {
    "G3_43102": ("Position-1 section on scan_01 (D1MT). Dimmest of its four "
                 "sections in 48.8% of composition-controlled comparisons "
                 "against a 25% null, and 68.7% on the all-cell version. No "
                 "other section on that scan exceeds 27.8%. Sits flush at "
                 "Y = 0.7 um against the scan boundary. IDO1 p99 = 8.3 on 317 "
                 "CD68-lineage cells, below the p25 of intact slide-mates, and "
                 "zero IDO1+ macrophages called. Cannot support intensity-based "
                 "claims."),
    "G4_43112": ("Position-1 section on scan_02 (Untreated). Dimmest of its "
                 "four sections in 56.3% of composition-controlled comparisons "
                 "against a 25% null, and 73.1% on the all-cell version. No "
                 "other section on that scan exceeds 23.3%. Sits flush at "
                 "Y = 0.7 um against the scan boundary. IDO1 max = 8.3 across "
                 "2,449 macrophages, and zero IDO1+ called. Cannot support "
                 "intensity-based claims."),
}
# Set False to skip the cross-check against script 02's table 12.
VERIFY_EXCLUSION_AGAINST_DIAGNOSTICS = True
EXCLUSION_RANK_NULL = 0.25       # dimmest-of-four under no position effect
EXCLUSION_RANK_FLAG = 0.40       # warn if a RETAINED section exceeds this

CONDITION_ORDER = ["D1MT", "Untreated"]
CONDITION_COLORS = {"D1MT": "#2C7FB8", "Untreated": "#D95F02"}
GROUP_MAP = {"G3": "D1MT", "G4": "Untreated"}
SCAN_COL = "scan_id"
POSITION_COL = "slide_position_rank"

# ---- column names -----------------------------------------------------------
INTENSITY_SUFFIX = ": Mean"
CENTROID_X_PREFIX = "Centroid X"
CENTROID_Y_PREFIX = "Centroid Y"
PHENOTYPE_COL = "Phenotypes"
IMAGE_COL = "Image"

IDO1_COL = "IDO1: Cytoplasm: Mean"
DAPI_COL = "DAPI: Nucleus: Mean"
HK3_COL = "3-Hydroxykynurenine: Cytoplasm: Mean"

# ---- phenotypes -------------------------------------------------------------
PHENOTYPE_ORDER = [
    "CD68+IDO1+ Macrophages", "CD68+IDO1- Macrophages", "CD163+ Macrophages",
    "Neutrophils", "Helper T cells", "CD4- T cells", "Tregs",
    "B cells", "Plasma cells", "Endothelial cells",
    "Epithelial/Tumor cells", "Other",
]

# Structural classes get real colours so tissue architecture is readable.
PHENOTYPE_COLORS = {
    "CD68+IDO1+ Macrophages": "#B2182B",
    "CD68+IDO1- Macrophages": "#EF8A62",
    "CD163+ Macrophages": "#FDBE85",
    "Neutrophils": "#7B3294",
    "Helper T cells": "#1B7837",
    "CD4- T cells": "#7FBC41",
    "Tregs": "#00441B",
    "B cells": "#2166AC",
    "Plasma cells": "#67A9CF",
    "Endothelial cells": "#E7298A",
    "Epithelial/Tumor cells": "#66C2A5",
    "Other": "#A6761D",
}

IMMUNE_PHENOTYPES = [
    "CD68+IDO1+ Macrophages", "CD68+IDO1- Macrophages", "CD163+ Macrophages",
    "Neutrophils", "Helper T cells", "CD4- T cells", "Tregs",
    "B cells", "Plasma cells",
]
STRUCTURAL_PHENOTYPES = ["Endothelial cells", "Epithelial/Tumor cells", "Other"]

IDO1_POS_PHENO = "CD68+IDO1+ Macrophages"
IDO1_NEG_PHENO = "CD68+IDO1- Macrophages"
CD68_LINEAGE = [IDO1_POS_PHENO, IDO1_NEG_PHENO]

# =============================================================================
# IDO1 CUTOFF WINDOW SEARCH
# SCRIPT 09 MUST USE THE IDENTICAL BLOCK. If either is edited, edit both.
# =============================================================================
N_CUTOFF_STEPS = 200          # grid resolution on the log1p scale
CUTOFF_GRID_HI_PERCENTILE = 99.9   # top of the grid, percentile of pooled IDO1
PLATEAU_TOL = 0.05            # stability: every section moves less than this
MIN_ARM_GAP = 0.20            # separation: min(untreated frac) - max(D1MT frac)
                              # must stay at or above this across the window.
                              # A window where every section calls 100 percent
                              # positive has a gap of zero and is inadmissible,
                              # which is what the old search kept picking.

# ---- tissue boundary distance ----------------------------------------------
BOUNDARY_GRID_UM = 50         # occupancy grid pitch for the distance transform
BOUNDARY_CLOSE_ITER = 2       # binary closing passes to fill small interior gaps
BOUNDARY_FILL_HOLES = True    # NOTE: this is lung. Filling holes treats airway
                              # and vessel lumens as tissue, so "depth" is depth
                              # into the section envelope, not into parenchyma.
                              # Set False to treat lumens as boundary.
N_BOUNDARY_BINS = 12
MIN_CELLS_PER_BOUNDARY_BIN = 30

# ---- key markers ------------------------------------------------------------
KEY_MARKERS = [
    "IDO1", "3-Hydroxykynurenine", "CD68", "CD163", "CD206", "CD4", "CD8",
    "CD3e", "CD20", "CD79a", "CD21", "MPO", "CD11b", "FoxP3", "iNOS",
    "Arginase-1", "IFNG", "Granzyme-B", "HLA-DR", "PD-1", "PD-L1", "DAPI",
]

# ---- plotting ---------------------------------------------------------------
MAX_POINTS_PER_SECTION = 80000
POINT_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]
RANDOM_SEED = 0

FONT_SIZE_BASE = 28
FONT_SIZE_TITLE = 34
FONT_SIZE_TICK = 28
FONT_SIZE_LEGEND = 28
FONT_SIZE_ANNOT = 26
DPI = 300
SAVE_PDF = True
SAVE_PNG = True

GRID_COLOR = "#DDDDDD"
AXIS_COLOR = "#333333"
TEXT_COLOR = "#000000"
FLAG_COLOR = "#B2182B"
OK_COLOR = "#1B7837"
BACKGROUND_COLOR = "#ECECEC"

USE_FLOAT32 = True


# %% Cell 2 - imports, style, helpers
# =============================================================================

import os
import sys
import glob
import gc
from datetime import datetime

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

try:
    from scipy import ndimage as ndi
    HAVE_SCIPY = True
except Exception as _e:
    HAVE_SCIPY = False
    _scipy_err = _e

plt.rcParams.update({
    "font.size": FONT_SIZE_BASE,
    "axes.titlesize": FONT_SIZE_TITLE,
    "axes.labelsize": FONT_SIZE_BASE,
    "xtick.labelsize": FONT_SIZE_TICK,
    "ytick.labelsize": FONT_SIZE_TICK,
    "legend.fontsize": FONT_SIZE_LEGEND,
    "figure.titlesize": FONT_SIZE_TITLE,
    "axes.edgecolor": AXIS_COLOR,
    "axes.labelcolor": TEXT_COLOR,
    "text.color": TEXT_COLOR,
    "xtick.color": AXIS_COLOR,
    "ytick.color": AXIS_COLOR,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

FIG_DIR = os.path.join(OUT_DIR, "figures")
TAB_DIR = os.path.join(OUT_DIR, "tables")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TAB_DIR, exist_ok=True)

rng = np.random.default_rng(RANDOM_SEED)
DIVERGING = LinearSegmentedColormap.from_list("div", ["#2166AC", "#F7F7F7", "#B2182B"])


class Tee(object):
    def __init__(self, path):
        self.terminal = sys.stdout
        self.log = open(path, "w", encoding="utf-8", newline="\n")

    def write(self, m):
        self.terminal.write(m); self.log.write(m)

    def flush(self):
        self.terminal.flush(); self.log.flush()

    def close(self):
        try:
            self.log.close()
        except Exception:
            pass


def banner(t, c="=", w=88):
    print("\n" + c * w); print(t); print(c * w)


def sub(t):
    print("\n" + t); print("-" * min(len(t), 88))


def ascii_safe(s):
    return "".join(ch for ch in str(s) if ord(ch) <= 127).strip()


def short_label(sid):
    return sid.split("_")[-1] if "_" in sid else sid


def marker_base(col):
    return col.split(":")[0].strip()


def save_fig(fig, stem):
    out = []
    if SAVE_PDF:
        p = os.path.join(FIG_DIR, f"{stem}.pdf")
        fig.savefig(p, dpi=DPI, bbox_inches="tight"); out.append(p)
    if SAVE_PNG:
        p = os.path.join(FIG_DIR, f"{stem}.png")
        fig.savefig(p, dpi=DPI, bbox_inches="tight"); out.append(p)
    plt.close(fig); gc.collect()
    for p in out:
        print(f"    wrote {p}")


def write_csv(df, name):
    p = os.path.join(TAB_DIR, name)
    df.to_csv(p, index=False, lineterminator="\n", encoding="utf-8")
    print(f"    wrote {p}  ({df.shape[0]} rows x {df.shape[1]} cols)")


def style_axes(ax, ygrid=True):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(2)
    ax.spines["bottom"].set_linewidth(2)
    if ygrid:
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, color=GRID_COLOR, linewidth=1.5)


def find_col(cols, prefix):
    hits = [c for c in cols if c.startswith(prefix)]
    return hits[0] if hits else None


def variance_components(values, groups):
    """
    One-way decomposition of `values` (already on a log scale) by `groups`.

    Returns naive_fraction and icc. naive_fraction is
    raw_between / (raw_between + within), which is what revision 1 reported and
    which is biased upward because the variance of group means already contains
    within-group variance over n. icc uses the corrected component
    (MSB - MSW) / n. Same construction as script 02 revision 2, so the two are
    directly comparable.
    """
    s = pd.Series(np.asarray(values, dtype=float))
    g = pd.Series(np.asarray(groups))
    means = s.groupby(g).mean()
    sizes = s.groupby(g).size()
    if len(means) < 2:
        return np.nan, np.nan, np.nan, np.nan
    n_bar = float(sizes.mean())
    msw = float(s.groupby(g).var(ddof=1).mean())
    raw_between = float(np.var(means.to_numpy(), ddof=1))
    naive_total = msw + raw_between
    naive = raw_between / naive_total if naive_total > 0 else np.nan
    msb = n_bar * raw_between
    sigma2_b = max(0.0, (msb - msw) / n_bar) if n_bar > 0 else np.nan
    tot = sigma2_b + msw
    icc = sigma2_b / tot if tot > 0 else np.nan
    return naive, icc, msw, raw_between


# =============================================================================
# THE CUTOFF WINDOW SEARCH. SCRIPT 09 MUST USE AN IDENTICAL COPY.
# =============================================================================

def arm_separation_gap(sw, arm_of):
    """
    sw: DataFrame, index = cutoff, columns = sample_id, values = fraction of
        CD68-lineage cells called positive at that cutoff.
    Returns a Series of the completely-separated gap at each cutoff:
        max(0, min(other arm) - max(this arm)) taken in whichever direction
        actually separates. Zero when the arms overlap.
    """
    a_cols = [c for c in sw.columns if arm_of.get(c) == CONDITION_ORDER[0]]
    b_cols = [c for c in sw.columns if arm_of.get(c) == CONDITION_ORDER[1]]
    if not a_cols or not b_cols:
        return pd.Series(np.nan, index=sw.index)
    gap_ba = sw[b_cols].min(axis=1) - sw[a_cols].max(axis=1)
    gap_ab = sw[a_cols].min(axis=1) - sw[b_cols].max(axis=1)
    return pd.concat([gap_ba, gap_ab], axis=1).max(axis=1).clip(lower=0.0)


def widest_stable_window(sw, tol, mask=None):
    """
    Widest contiguous run of cutoffs over which EVERY column of sw moves by less
    than tol. If mask is given, only cutoffs where mask is True are eligible and
    a window may not cross an ineligible cutoff.

    Returns (lo, hi, width). NaN when nothing qualifies.
    """
    idx = sw.index.to_numpy(dtype=float)
    ok = (np.ones(len(sw), dtype=bool) if mask is None
          else np.asarray(mask, dtype=bool))
    best = (np.nan, np.nan, -1.0)
    i = 0
    n = len(sw)
    while i < n:
        if not ok[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and ok[j + 1]:
            block = sw.iloc[i:j + 2]
            if float((block.max() - block.min()).max()) > tol:
                break
            j += 1
        width = idx[j] - idx[i]
        if width > best[2]:
            best = (float(idx[i]), float(idx[j]), float(width))
        i = max(j, i) + 1
    return best


def contiguous_runs(mask, index):
    """Return [(lo, hi, n_steps), ...] for each True run of mask."""
    m = np.asarray(mask, dtype=bool)
    idx = np.asarray(index, dtype=float)
    runs, i, n = [], 0, len(m)
    while i < n:
        if not m[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and m[j + 1]:
            j += 1
        runs.append((float(idx[i]), float(idx[j]), j - i + 1))
        i = j + 1
    return runs


_tee = Tee(os.path.join(TAB_DIR, "00_rebaseline_report.txt"))
sys.stdout = _tee

banner("AKOYA RE-BASELINE ON SIX SECTIONS (revision 2)")
print(f"Run time : {datetime.now().isoformat(timespec='seconds')}")
print(f"Data     : {DATA_DIR}")
print(f"Output   : {OUT_DIR}")
if not HAVE_SCIPY:
    print(f"\n    WARNING: scipy unavailable ({_scipy_err}). The tissue boundary "
          f"distance analysis will be skipped.")

sub("EXCLUDED SECTIONS")
for sid, reason in EXCLUDE_SECTIONS.items():
    print(f"    {sid}")
    print(f"        {reason}")


# %% Cell 3 - inventory metadata and exclusion cross-check
# =============================================================================

banner("INVENTORY METADATA AND EXCLUSION CROSS-CHECK")

SCAN_OF, POS_OF = {}, {}
fs_path = os.path.join(INVENTORY_TABLE_DIR, "01_file_summary.csv")
if os.path.exists(fs_path):
    INV = pd.read_csv(fs_path)
    print(f"    loaded 01_file_summary.csv  ({INV.shape[0]} x {INV.shape[1]})")
    if SCAN_COL in INV.columns:
        SCAN_OF = dict(zip(INV["sample_id"], INV[SCAN_COL].astype(str)))
    else:
        print(f"    WARNING: '{SCAN_COL}' absent, falling back to the Image column")
    if POSITION_COL in INV.columns:
        POS_OF = {r["sample_id"]: (int(r[POSITION_COL])
                                   if pd.notna(r[POSITION_COL]) else None)
                  for _, r in INV.iterrows()}
else:
    print(f"    WARNING: {fs_path} not found. Run script 01 first.")

if VERIFY_EXCLUSION_AGAINST_DIAGNOSTICS:
    p12 = os.path.join(DIAGNOSTIC_TABLE_DIR, "12_slide_position_summary.csv")
    if not os.path.exists(p12):
        print(f"    WARNING: {p12} not found, exclusion not cross-checked.")
    else:
        d12 = pd.read_csv(p12)
        col = ("frac_dimmest_by_phenotype"
               if "frac_dimmest_by_phenotype" in d12.columns
               else "frac_dimmest_of_scan")
        print(f"\n    Cross-check against script 02, column '{col}' "
              f"(null = {EXCLUSION_RANK_NULL:.2f})")
        print(f"      {'section':<14}{'pos':>4}{'dimmest':>10}   status")
        for _, r in d12.sort_values([c for c in ["scan_id", POSITION_COL]
                                     if c in d12.columns]).iterrows():
            sid = r["sample_id"]
            v = float(r[col]) if pd.notna(r[col]) else np.nan
            excluded = sid in EXCLUDE_SECTIONS
            pos = r.get(POSITION_COL, np.nan)
            pos_txt = f"{int(pos)}" if pd.notna(pos) else "na"
            if excluded:
                status = ("excluded, evidence holds"
                          if v > EXCLUSION_RANK_FLAG else
                          "EXCLUDED BUT EVIDENCE NO LONGER SUPPORTS IT")
            else:
                status = ("retained"
                          if v <= EXCLUSION_RANK_FLAG else
                          "RETAINED BUT LOOKS DIM, REVIEW")
            flag = "   <--" if "REVIEW" in status or "NO LONGER" in status else ""
            print(f"      {sid:<14}{pos_txt:>4}{100*v:>9.1f}%   {status}{flag}")


# %% Cell 4 - single pass over the retained sections
# =============================================================================

banner("LOADING RETAINED SECTIONS")

csv_paths = sorted(glob.glob(os.path.join(DATA_DIR, "*.csv")))
retained = [p for p in csv_paths
            if os.path.splitext(os.path.basename(p))[0] not in EXCLUDE_SECTIONS]
print(f"    {len(csv_paths)} files found, {len(EXCLUDE_SECTIONS)} excluded, "
      f"{len(retained)} retained")

if not retained:
    print("ERROR: nothing retained.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

hdr = pd.read_csv(retained[0], nrows=0)
ALL_COLS = hdr.columns.tolist()
XCOL = find_col(ALL_COLS, CENTROID_X_PREFIX)
YCOL = find_col(ALL_COLS, CENTROID_Y_PREFIX)
INTENSITY_COLS = [c for c in ALL_COLS if c.endswith(INTENSITY_SUFFIX)]
dtype_map = ({c: "float32" for c in INTENSITY_COLS + [XCOL, YCOL]}
             if USE_FLOAT32 else None)

sections = {}
pheno_counts = []
section_marker_stats = []   # per-section median AND p99 for every marker
pheno_profiles = []
cd68_cells = {}     # sid -> DataFrame of CD68-lineage cells

for path in retained:
    sid = os.path.splitext(os.path.basename(path))[0]
    cond = GROUP_MAP.get(sid.split("_")[0], "UNKNOWN")
    print(f"\n[{sid}]  {cond}")
    try:
        df = pd.read_csv(path, dtype=dtype_map, low_memory=False)
    except Exception as e:
        print(f"    ERROR reading: {e}. Skipping and continuing.")
        continue

    n = len(df)
    if sid in SCAN_OF:
        scan, scan_src = SCAN_OF[sid], "inventory"
    elif IMAGE_COL in df.columns and n:
        scan, scan_src = str(df[IMAGE_COL].iloc[0]), "Image column"
    else:
        scan, scan_src = f"scan_{sid.split('_')[0]}", "filename prefix"

    df["_pheno"] = df[PHENOTYPE_COL].map(ascii_safe) if PHENOTYPE_COL in df.columns \
        else "Unknown"

    x = pd.to_numeric(df[XCOL], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(df[YCOL], errors="coerce").to_numpy(dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    print(f"    n={n:,}   scan={scan} ({scan_src})   pos={POS_OF.get(sid)}")

    idx = np.flatnonzero(ok)
    idx_plot = (rng.choice(idx, size=MAX_POINTS_PER_SECTION, replace=False)
                if len(idx) > MAX_POINTS_PER_SECTION else idx)

    sections[sid] = {
        "condition": cond, "scan": scan, "position": POS_OF.get(sid),
        "n_cells": n,
        "x_plot": x[idx_plot], "y_plot": y[idx_plot],
        "pheno_plot": df["_pheno"].to_numpy()[idx_plot],
        "x_lo": float(np.nanmin(x)), "x_hi": float(np.nanmax(x)),
        "y_lo": float(np.nanmin(y)), "y_hi": float(np.nanmax(y)),
    }

    vc = df["_pheno"].value_counts()
    for ph, c in vc.items():
        pheno_counts.append({"sample_id": sid, "condition": cond,
                             "scan_id": scan, "phenotype": ph, "n_cells": int(c),
                             "pct_of_section": 100.0 * c / n})

    med = df[INTENSITY_COLS].median(axis=0)
    p99 = df[INTENSITY_COLS].quantile(0.99)
    for c in INTENSITY_COLS:
        section_marker_stats.append({
            "sample_id": sid, "condition": cond, "scan_id": scan,
            "column": c, "median": float(med[c]), "p99": float(p99[c])})

    prof = df.groupby("_pheno")[INTENSITY_COLS].mean()
    for ph, row in prof.iterrows():
        for c, v in row.items():
            pheno_profiles.append({"sample_id": sid, "condition": cond,
                                   "phenotype": ph, "column": c, "mean": float(v)})

    # ---- CD68-lineage cell table -------------------------------------------
    m = df["_pheno"].isin(CD68_LINEAGE).to_numpy() & ok
    keep = pd.DataFrame({
        "x": x[m], "y": y[m],
        "label": df.loc[m, "_pheno"].to_numpy(),
        "ido1": pd.to_numeric(df.loc[m, IDO1_COL], errors="coerce").to_numpy()
        if IDO1_COL in df.columns else np.nan,
        "dapi": pd.to_numeric(df.loc[m, DAPI_COL], errors="coerce").to_numpy()
        if DAPI_COL in df.columns else np.nan,
        "hk3": pd.to_numeric(df.loc[m, HK3_COL], errors="coerce").to_numpy()
        if HK3_COL in df.columns else np.nan,
    })
    cd68_cells[sid] = keep
    print(f"    CD68-lineage cells: {len(keep):,}  "
          f"({int((keep['label'] == IDO1_POS_PHENO).sum()):,} IDO1+)")

    sections[sid]["x_all"] = x[ok]
    sections[sid]["y_all"] = y[ok]

    del df
    gc.collect()

cond_rank = {c: i for i, c in enumerate(CONDITION_ORDER)}
SAMPLE_ORDER = sorted(sections, key=lambda s: (cond_rank.get(sections[s]["condition"], 9), s))
ARM_OF = {s: sections[s]["condition"] for s in SAMPLE_ORDER}
MARKER_OF = {s: POINT_MARKERS[i % len(POINT_MARKERS)]
             for i, s in enumerate(SAMPLE_ORDER)}

SHADES = {"D1MT": ["#08519C", "#3182BD", "#6BAED6", "#BDD7E7"],
          "Untreated": ["#A63603", "#E6550D", "#FD8D3C", "#FDBE85"]}
COLOR_OF = {}
_seen = {c: 0 for c in CONDITION_ORDER}
for s in SAMPLE_ORDER:
    c = sections[s]["condition"]
    pal = SHADES.get(c, ["#999999"])
    COLOR_OF[s] = pal[_seen.get(c, 0) % len(pal)]
    _seen[c] = _seen.get(c, 0) + 1


def section_legend(ax, **kw):
    handles = [Line2D([0], [0], color=COLOR_OF[s], marker=MARKER_OF[s],
                      markersize=16, linewidth=4,
                      label=f"{short_label(s)} ({sections[s]['condition']})")
               for s in SAMPLE_ORDER]
    ax.legend(handles=handles, frameon=False,
              fontsize=kw.pop("fontsize", FONT_SIZE_LEGEND - 10), **kw)


def sect_tick(s):
    p = sections[s]["position"]
    return (f"{short_label(s)}\n{sections[s]['condition']}\n"
            f"{sections[s]['scan']} pos{p if p else 'na'}")


# %% Cell 5 - A: re-baselined composition and IDO1 fraction
# =============================================================================

banner("A - RE-BASELINED COMPOSITION AND IDO1 FRACTION")

pc = pd.DataFrame(pheno_counts)
counts = (pc.pivot_table(index="phenotype", columns="sample_id", values="n_cells",
                         aggfunc="sum")
          .reindex(index=[p for p in PHENOTYPE_ORDER],
                   columns=SAMPLE_ORDER).fillna(0.0))
PHENO_USE = [p for p in PHENOTYPE_ORDER if p in counts.index]
counts = counts.loc[PHENO_USE]
pcts = 100.0 * counts / counts.sum(axis=0)

print(f"    Total cells retained: {int(counts.sum().sum()):,}")
for c in CONDITION_ORDER:
    mem = [s for s in SAMPLE_ORDER if sections[s]["condition"] == c]
    print(f"      {c:<12} {len(mem)} sections, "
          f"{int(counts[mem].sum().sum()):,} cells")

num = counts.loc[IDO1_POS_PHENO, SAMPLE_ORDER].to_numpy(float)
den = counts.loc[CD68_LINEAGE, SAMPLE_ORDER].sum(axis=0).to_numpy(float)
frac = np.where(den > 0, 100.0 * num / den, np.nan)

ido_frac = pd.DataFrame({
    "sample_id": SAMPLE_ORDER,
    "animal_id": [short_label(s) for s in SAMPLE_ORDER],
    "condition": [sections[s]["condition"] for s in SAMPLE_ORDER],
    "scan_id": [sections[s]["scan"] for s in SAMPLE_ORDER],
    "slide_position_rank": [sections[s]["position"] for s in SAMPLE_ORDER],
    "n_ido1_pos": num.astype(int), "n_cd68_lineage": den.astype(int),
    "pct_ido1_pos": frac,
})
sub("IDO1+ as a fraction of CD68-lineage macrophages (six sections)")
print(ido_frac.to_string(index=False))

xpos = np.arange(len(SAMPLE_ORDER))
cond_labels = [sections[s]["condition"] for s in SAMPLE_ORDER]
x_labels = [short_label(s) for s in SAMPLE_ORDER]

fig, axes = plt.subplots(1, 2, figsize=(30, 14),
                         gridspec_kw={"width_ratios": [1.35, 1.0]})
ax = axes[0]
bottom = np.zeros(len(SAMPLE_ORDER))
for ph in PHENO_USE:
    v = pcts.loc[ph, SAMPLE_ORDER].to_numpy(float)
    ax.bar(xpos, v, bottom=bottom, color=PHENOTYPE_COLORS.get(ph, "#999999"),
           edgecolor="#FFFFFF", linewidth=1.5, label=ph, zorder=3)
    bottom += v
ax.set_xticks(xpos)
ax.set_xticklabels([sect_tick(s) for s in SAMPLE_ORDER],
                   fontsize=FONT_SIZE_TICK - 12)
ax.set_ylabel("Percent of section"); ax.set_ylim(0, 100)
ax.set_title("Composition, six retained sections", fontsize=FONT_SIZE_TITLE - 6)
style_axes(ax, ygrid=False)
ax.legend(bbox_to_anchor=(1.01, 1.0), loc="upper left", frameon=False,
          fontsize=FONT_SIZE_LEGEND - 12)

ax = axes[1]
ax.bar(xpos, np.nan_to_num(frac),
       color=[CONDITION_COLORS.get(c, "#999") for c in cond_labels],
       edgecolor="#FFFFFF", linewidth=2, zorder=3)
for i in range(len(SAMPLE_ORDER)):
    ax.text(i, frac[i] + 2, f"{frac[i]:.1f}%\n{int(num[i]):,}/{int(den[i]):,}",
            ha="center", va="bottom", fontsize=FONT_SIZE_ANNOT - 8)
ax.set_xticks(xpos)
ax.set_xticklabels([sect_tick(s) for s in SAMPLE_ORDER],
                   fontsize=FONT_SIZE_TICK - 12)
ax.set_ylabel("IDO1+ % of CD68-lineage"); ax.set_ylim(0, 110)
ax.set_title("IDO1+ fraction", fontsize=FONT_SIZE_TITLE - 6)
style_axes(ax)
fig.suptitle("Re-baselined on six sections "
             "(position-1 sections of each scan excluded)", y=1.04,
             fontsize=FONT_SIZE_TITLE - 4)
save_fig(fig, "F17_rebaselined_composition")

write_csv(pc, "20_rebaselined_phenotype_counts.csv")
write_csv(ido_frac, "21_ido1_fraction_six_sections.csv")

# ---- batch decomposition on six sections, reconciled with script 02 ---------
sub("Between-scan vs within-scan variance on six sections")
print("    Computed FOUR ways so the numbers carried in the project notes can")
print("    be attributed. Revision 1 of this cell grouped by condition, used")
print("    section-wide MEDIANS, and reported the naive ratio; script 02 groups")
print("    by scan, uses p99, and reports the corrected ICC. Both statistics on")
print("    the scan grouping are given here, on the six retained sections.\n")
print("    Scan and arm are still the same factor. A high score is consistent")
print("    with batch AND with biology. Only the low end is interpretable.\n")

sms = pd.DataFrame(section_marker_stats)
rows = []
for col, g in sms.groupby("column"):
    rec = {"column": col, "marker": marker_base(col),
           "is_key_marker": marker_base(col) in KEY_MARKERS}
    for stat in ("p99", "median"):
        gg = g.loc[g[stat] > 0]
        if gg["scan_id"].nunique() < 2 or len(gg) < 4:
            rec[f"naive_{stat}"] = np.nan
            rec[f"icc_{stat}"] = np.nan
            rec[f"log2_scan_diff_{stat}"] = np.nan
            continue
        lv = np.log2(gg[stat].to_numpy(dtype=float))
        naive, icc, msw, raw_b = variance_components(lv, gg["scan_id"].to_numpy())
        rec[f"naive_{stat}"] = naive
        rec[f"icc_{stat}"] = icc
        rec[f"within_scan_var_{stat}"] = msw
        sm = pd.Series(lv).groupby(gg["condition"].to_numpy()).mean()
        rec[f"log2_scan_diff_{stat}"] = float(sm.get(CONDITION_ORDER[1], np.nan)
                                              - sm.get(CONDITION_ORDER[0], np.nan))
    rows.append(rec)
var6 = pd.DataFrame(rows).sort_values("icc_p99", ascending=False)

print(f"    {'column':<44}{'ICC p99':>9}{'naive p99':>11}{'ICC med':>9}"
      f"{'naive med':>11}{'log2 diff':>11}")
print("    " + "-" * 91)
for _, r in var6.loc[var6["is_key_marker"]].sort_values("icc_p99",
                                                        ascending=False).iterrows():
    print(f"    {r['column'][:43]:<44}{r['icc_p99']:>9.3f}{r['naive_p99']:>11.3f}"
          f"{r['icc_median']:>9.3f}{r['naive_median']:>11.3f}"
          f"{r['log2_scan_diff_p99']:>+11.2f}")

sub("Where the two statistics disagree most (key markers)")
km = var6.loc[var6["is_key_marker"]].copy()
km["abs_gap"] = (km["icc_p99"] - km["icc_median"]).abs()
for _, r in km.sort_values("abs_gap", ascending=False).head(8).iterrows():
    print(f"    {r['column'][:43]:<44} p99 {r['icc_p99']:.3f}   "
          f"median {r['icc_median']:.3f}   gap {r['abs_gap']:.3f}")
print("\n    A marker with a high p99 ICC and a low median ICC differs between")
print("    scans in its bright tail but not in its bulk. That is what a")
print("    threshold applied to the top of the distribution would see, and it")
print("    is the relevant number whenever a gate is involved.")
write_csv(var6, "22_variance_decomposition_six_sections.csv")


# %% Cell 6 - B: back out the IDO1 threshold behind the vendor's binary call
# =============================================================================

banner("B - IDO1 THRESHOLD BACK-OUT")

print("    For each section: the lowest IDO1 among cells labelled IDO1+, and the")
print("    highest IDO1 among cells labelled IDO1-. If these meet at the same")
print("    value in every section, the call was a fixed global cut. If they meet")
print("    at different values, it was per-section adaptive. If the ranges")
print("    OVERLAP, the call used more than IDO1 alone.\n")

back_rows = []
for s in SAMPLE_ORDER:
    d = cd68_cells.get(s)
    if d is None or not len(d):
        continue
    pos = d.loc[d["label"] == IDO1_POS_PHENO, "ido1"].dropna().to_numpy()
    neg = d.loc[d["label"] == IDO1_NEG_PHENO, "ido1"].dropna().to_numpy()
    rec = {"sample_id": s, "animal_id": short_label(s),
           "condition": sections[s]["condition"],
           "scan_id": sections[s]["scan"],
           "n_pos": len(pos), "n_neg": len(neg),
           "min_pos": float(np.min(pos)) if len(pos) else np.nan,
           "max_neg": float(np.max(neg)) if len(neg) else np.nan,
           "p01_pos": float(np.percentile(pos, 1)) if len(pos) else np.nan,
           "p99_neg": float(np.percentile(neg, 99)) if len(neg) else np.nan}
    if len(pos) and len(neg):
        rec["separable"] = bool(rec["min_pos"] > rec["max_neg"])
        rec["frac_neg_above_min_pos"] = float(np.mean(neg > rec["min_pos"]))
        rec["frac_pos_below_max_neg"] = float(np.mean(pos < rec["max_neg"]))
        rec["implied_cut"] = (0.5 * (rec["min_pos"] + rec["max_neg"])
                              if rec["separable"] else np.nan)
    back_rows.append(rec)

back = pd.DataFrame(back_rows)
with pd.option_context("display.width", 240, "display.max_columns", 30):
    print(back.to_string(index=False))

sub("Verdict")
if len(back) and "separable" in back.columns:
    sep = back["separable"].fillna(False)
    if sep.all():
        cuts = back["implied_cut"].dropna()
        spread = cuts.max() / cuts.min() if cuts.min() > 0 else np.inf
        print(f"    Every section is cleanly separable on IDO1 alone.")
        print(f"    Implied cut per section: "
              f"{', '.join(f'{v:.2f}' for v in cuts)}")
        if spread < 1.5:
            print(f"    Spread is {spread:.2f}x -> consistent with a FIXED GLOBAL "
                  f"cut applied across sections.")
        else:
            print(f"    Spread is {spread:.2f}x -> the cut was PER-SECTION "
                  f"ADAPTIVE. Any cross-section comparison of the binary call")
            print(f"    carries this, and it must be stated in methods.")
    else:
        bad = back.loc[~sep, "sample_id"].tolist()
        print(f"    Ranges OVERLAP in: {bad}")
        print(f"    The binary call therefore used more than IDO1 intensity")
        print(f"    alone (co-markers, spatial smoothing, or a classifier).")
        print(f"    It cannot be reproduced as a simple IDO1 gate, and we should")
        print(f"    ask Akoya what the classifier actually used before relying")
        print(f"    on the call for anything quantitative.")
write_csv(back, "23_ido1_threshold_backout.csv")


# %% Cell 7 - C: threshold sensitivity with the separation constraint
# =============================================================================

banner("C / D - THRESHOLD SENSITIVITY AND IDO1 DISTRIBUTIONS")

all_ido = np.concatenate([cd68_cells[s]["ido1"].dropna().to_numpy()
                          for s in SAMPLE_ORDER if len(cd68_cells.get(s, []))])
lo = 0.0
hi = float(np.log1p(np.percentile(all_ido, CUTOFF_GRID_HI_PERCENTILE)))
grid = np.linspace(lo, hi, N_CUTOFF_STEPS)

sub("Inputs to the cutoff search (script 09 must match these exactly)")
print(f"    grid            : {N_CUTOFF_STEPS} steps, log1p {lo:.4f} to {hi:.4f}")
print(f"    grid top        : p{CUTOFF_GRID_HI_PERCENTILE} of pooled IDO1 = "
      f"{np.percentile(all_ido, CUTOFF_GRID_HI_PERCENTILE):.4f} raw")
print(f"    pooled cells    : {len(all_ido):,}")
print(f"    PLATEAU_TOL     : {PLATEAU_TOL}")
print(f"    MIN_ARM_GAP     : {MIN_ARM_GAP}")
print(f"\n    {'section':<14}{'n cells':>10}{'frac>0 at cutoff 0':>22}")
sens_rows = []
for s in SAMPLE_ORDER:
    v = cd68_cells[s]["ido1"].dropna().to_numpy()
    if not len(v):
        continue
    lv = np.log1p(v)
    print(f"    {s:<14}{len(v):>10,}{float(np.mean(lv > 0.0)):>22.6f}")
    for g in grid:
        sens_rows.append({"sample_id": s, "condition": sections[s]["condition"],
                          "log1p_cutoff": float(g),
                          "frac_above": float(np.mean(lv > g))})
sens = pd.DataFrame(sens_rows)

sw = sens.pivot_table(index="log1p_cutoff", columns="sample_id",
                      values="frac_above")[SAMPLE_ORDER]

# ---- step 1: where do the arms stay completely separated? -------------------
gap = arm_separation_gap(sw, ARM_OF)
admissible = (gap >= MIN_ARM_GAP).to_numpy()
runs = contiguous_runs(admissible, sw.index)

sub("Step 1: arm separation")
print(f"    A cutoff is admissible when the arms are completely separated by at")
print(f"    least {MIN_ARM_GAP:.2f} in called fraction.\n")
if not runs:
    print("    NO admissible cutoff. The arms never separate by that margin.")
    print("    Do not present a cutoff-insensitivity panel. Report the")
    print("    distributions instead.")
    sep_lo = sep_hi = np.nan
else:
    for r_lo, r_hi, r_n in runs:
        print(f"    log1p {r_lo:.3f} to {r_hi:.3f}   raw IDO1 "
              f"{np.expm1(r_lo):.2f} to {np.expm1(r_hi):.2f}   "
              f"({r_n} grid steps)")
    widest = max(runs, key=lambda t: t[1] - t[0])
    sep_lo, sep_hi = widest[0], widest[1]
    peak_gap = float(gap.max())
    peak_at = float(gap.idxmax())
    print(f"\n    widest admissible run: log1p {sep_lo:.3f} to {sep_hi:.3f}, "
          f"raw IDO1 {np.expm1(sep_lo):.2f} to {np.expm1(sep_hi):.2f}")
    print(f"    maximum separation   : {peak_gap:.3f} at log1p {peak_at:.3f} "
          f"(raw {np.expm1(peak_at):.2f})")

# ---- step 2: stability inside the admissible region ------------------------
best_lo, best_hi, best_w = widest_stable_window(sw, PLATEAU_TOL, mask=admissible)
sub("Step 2: stability inside the admissible region")
if np.isfinite(best_lo):
    print(f"    Widest window where every section's called fraction moves by")
    print(f"    less than {PLATEAU_TOL:.2f} AND the arms stay separated:")
    print(f"        log1p scale : {best_lo:.3f} to {best_hi:.3f}")
    print(f"        raw IDO1    : {np.expm1(best_lo):.2f} to {np.expm1(best_hi):.2f}")
    at_lo = sw.loc[sw.index[np.argmin(np.abs(sw.index.to_numpy() - best_lo))]]
    at_hi = sw.loc[sw.index[np.argmin(np.abs(sw.index.to_numpy() - best_hi))]]
    print(f"\n    Called fraction at the two ends of that window:")
    print(f"      {'section':<14}{'at lo':>10}{'at hi':>10}")
    for s in SAMPLE_ORDER:
        print(f"      {s:<14}{100*at_lo[s]:>9.1f}%{100*at_hi[s]:>9.1f}%")
else:
    print("    No stable window inside the admissible region at this tolerance.")

# ---- what the unconstrained search would have returned ---------------------
unc_lo, unc_hi, unc_w = widest_stable_window(sw, PLATEAU_TOL, mask=None)
sub("For comparison: the unconstrained search (revision 1 behaviour)")
print(f"    log1p {unc_lo:.3f} to {unc_hi:.3f}, raw IDO1 "
      f"{np.expm1(unc_lo):.2f} to {np.expm1(unc_hi):.2f}")
if np.isfinite(unc_lo):
    g_at = float(gap.loc[gap.index[np.argmin(np.abs(gap.index.to_numpy()
                                                    - 0.5 * (unc_lo + unc_hi)))]])
    print(f"    arm separation at the centre of that window: {g_at:.3f}")
    if g_at < MIN_ARM_GAP:
        print(f"    That is below MIN_ARM_GAP, so the unconstrained window is")
        print(f"    the degenerate one: flat because every section calls the")
        print(f"    same fraction, not because the result is robust. This is")
        print(f"    the bug the constraint fixes.")

window_tbl = pd.DataFrame([{
    "n_cutoff_steps": N_CUTOFF_STEPS, "grid_hi_percentile": CUTOFF_GRID_HI_PERCENTILE,
    "grid_lo_log1p": lo, "grid_hi_log1p": hi,
    "plateau_tol": PLATEAU_TOL, "min_arm_gap": MIN_ARM_GAP,
    "separation_lo_log1p": sep_lo, "separation_hi_log1p": sep_hi,
    "separation_lo_raw": np.expm1(sep_lo), "separation_hi_raw": np.expm1(sep_hi),
    "stable_lo_log1p": best_lo, "stable_hi_log1p": best_hi,
    "stable_lo_raw": np.expm1(best_lo), "stable_hi_raw": np.expm1(best_hi),
    "unconstrained_lo_log1p": unc_lo, "unconstrained_hi_log1p": unc_hi,
    "n_pooled_cells": int(len(all_ido)),
}])
write_csv(window_tbl, "26_ido1_cutoff_windows.csv")

sens = sens.merge(gap.rename("arm_separation_gap").reset_index(),
                  on="log1p_cutoff", how="left")
write_csv(sens, "24_threshold_sensitivity.csv")

# ---- F18 ---------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(42, 14))

# --- ECDF ---
ax = axes[0]
for s in SAMPLE_ORDER:
    v = np.sort(cd68_cells[s]["ido1"].dropna().to_numpy())
    if not len(v):
        continue
    yv = np.arange(1, len(v) + 1) / len(v)
    step = max(1, len(v) // 25)
    ax.plot(np.log1p(v), yv, linewidth=4.5, color=COLOR_OF[s], alpha=0.9)
    ax.plot(np.log1p(v)[::step], yv[::step], linestyle="none",
            marker=MARKER_OF[s], markersize=15, color=COLOR_OF[s])
if np.isfinite(best_lo):
    ax.axvspan(best_lo, best_hi, color=OK_COLOR, alpha=0.12, zorder=0)
ax.set_xlabel("log(1 + IDO1 cytoplasm)")
ax.set_ylabel("Cumulative fraction of CD68+ cells")
ax.set_ylim(0, 1.02)
ax.set_title("Cumulative distribution", fontsize=FONT_SIZE_TITLE - 8)
style_axes(ax)
section_legend(ax, loc="lower right")

# --- threshold sensitivity, with the separation curve underneath -------------
ax = axes[1]
if np.isfinite(sep_lo):
    ax.axvspan(sep_lo, sep_hi, color=OK_COLOR, alpha=0.08, zorder=0)
if np.isfinite(best_lo):
    ax.axvspan(best_lo, best_hi, color=OK_COLOR, alpha=0.18, zorder=0)
    ax.text(0.5 * (best_lo + best_hi), 104,
            f"separated and stable\nraw IDO1 {np.expm1(best_lo):.1f} to "
            f"{np.expm1(best_hi):.1f}",
            ha="center", fontsize=FONT_SIZE_ANNOT - 10)
for s in SAMPLE_ORDER:
    d = sens.loc[sens["sample_id"] == s].sort_values("log1p_cutoff")
    ax.plot(d["log1p_cutoff"], 100 * d["frac_above"], linewidth=4.5,
            color=COLOR_OF[s], alpha=0.9)
    step = max(1, len(d) // 18)
    ax.plot(d["log1p_cutoff"].to_numpy()[::step],
            100 * d["frac_above"].to_numpy()[::step],
            linestyle="none", marker=MARKER_OF[s], markersize=15,
            color=COLOR_OF[s])
ax.plot(gap.index.to_numpy(), 100 * gap.to_numpy(), linewidth=5,
        color="#000000", linestyle=":", zorder=5)
ax.axhline(100 * MIN_ARM_GAP, color="#000000", linestyle="--", linewidth=2.5)
ax.set_xlabel("log(1 + IDO1) cutoff")
ax.set_ylabel("% of CD68-lineage cells above cutoff")
ax.set_ylim(0, 112)
ax.set_title("Threshold sensitivity\n"
             "dotted black = arm separation gap",
             fontsize=FONT_SIZE_TITLE - 14)
style_axes(ax)

# --- ridgeline ---
ax = axes[2]
bins = np.linspace(0, hi, 70)
for i, s in enumerate(SAMPLE_ORDER):
    v = cd68_cells[s]["ido1"].dropna().to_numpy()
    if not len(v):
        continue
    h, e = np.histogram(np.log1p(v), bins=bins, density=True)
    h = h / (h.max() if h.max() > 0 else 1.0)
    ctr = 0.5 * (e[:-1] + e[1:])
    base = len(SAMPLE_ORDER) - i
    ax.fill_between(ctr, base, base + 0.9 * h, color=COLOR_OF[s], alpha=0.75,
                    linewidth=0)
    ax.plot(ctr, base + 0.9 * h, color="#333333", linewidth=2.5)
    ax.text(hi * 1.01, base + 0.35, f"{short_label(s)}  n={len(v):,}",
            fontsize=FONT_SIZE_ANNOT - 8, va="center")
if np.isfinite(best_lo):
    ax.axvspan(best_lo, best_hi, color=OK_COLOR, alpha=0.12, zorder=0)
ax.set_yticks([]); ax.set_xlabel("log(1 + IDO1 cytoplasm)")
ax.set_title("Per-section distribution", fontsize=FONT_SIZE_TITLE - 8)
for sp in ["top", "right", "left"]:
    ax.spines[sp].set_visible(False)

fig.suptitle("IDO1 on CD68-lineage macrophages, six sections\n"
             "Shaded band = cutoffs where the arms stay completely separated "
             "AND every section's called fraction is stable", y=1.05,
             fontsize=FONT_SIZE_TITLE - 6)
save_fig(fig, "F18_ido1_threshold_and_distribution")


# %% Cell 8 - E: tissue boundary distance gradient
# =============================================================================

banner("E - TISSUE BOUNDARY DISTANCE GRADIENT")

bgrad = pd.DataFrame()
if not HAVE_SCIPY:
    print("    SKIPPED: scipy not available.")
else:
    print(f"    Occupancy grid at {BOUNDARY_GRID_UM} um, distance transform from "
          f"the tissue exterior.")
    print(f"    BOUNDARY_FILL_HOLES = {BOUNDARY_FILL_HOLES}. With holes filled, "
          f"airway and vessel")
    print(f"    lumens count as tissue, so this measures depth into the section "
          f"envelope\n    rather than depth into parenchyma. This is lung, so "
          f"that distinction matters.\n")
    rows = []
    for s in SAMPLE_ORDER:
        xa, ya = sections[s]["x_all"], sections[s]["y_all"]
        if not len(xa):
            continue
        x0, y0 = xa.min(), ya.min()
        gx = ((xa - x0) / BOUNDARY_GRID_UM).astype(int)
        gy = ((ya - y0) / BOUNDARY_GRID_UM).astype(int)
        occ = np.zeros((gy.max() + 3, gx.max() + 3), dtype=bool)
        occ[gy + 1, gx + 1] = True
        occ = ndi.binary_closing(occ, iterations=BOUNDARY_CLOSE_ITER)
        if BOUNDARY_FILL_HOLES:
            occ = ndi.binary_fill_holes(occ)
        dist = ndi.distance_transform_edt(occ) * BOUNDARY_GRID_UM

        d = cd68_cells[s]
        cgx = ((d["x"].to_numpy() - x0) / BOUNDARY_GRID_UM).astype(int) + 1
        cgy = ((d["y"].to_numpy() - y0) / BOUNDARY_GRID_UM).astype(int) + 1
        cgx = np.clip(cgx, 0, occ.shape[1] - 1)
        cgy = np.clip(cgy, 0, occ.shape[0] - 1)
        dd = dist[cgy, cgx]

        edges = np.quantile(dd, np.linspace(0, 1, N_BOUNDARY_BINS + 1))
        edges = np.unique(edges)
        if len(edges) < 3:
            print(f"    WARNING: {s} has too little depth range, skipped")
            continue
        b = np.clip(np.digitize(dd, edges[1:-1]), 0, len(edges) - 2)
        tmp = d.assign(_bin=b, _dist=dd)
        for bb, g in tmp.groupby("_bin"):
            if len(g) < MIN_CELLS_PER_BOUNDARY_BIN:
                continue
            rows.append({
                "sample_id": s, "condition": sections[s]["condition"],
                "depth_bin": int(bb), "n_cells": len(g),
                "median_depth_um": float(g["_dist"].median()),
                "median_ido1": float(g["ido1"].median()),
                "median_dapi": float(g["dapi"].median()),
                "median_hk3": float(g["hk3"].median()),
                "pct_ido1_pos": 100.0 * float((g["label"] == IDO1_POS_PHENO).mean()),
            })
        print(f"    {s:<12} depth range {dd.min():.0f} to {dd.max():.0f} um, "
              f"{len(np.unique(b))} bins")
    bgrad = pd.DataFrame(rows)

    if len(bgrad):
        fig, axes = plt.subplots(1, 3, figsize=(42, 13))
        panels = [("median_ido1", "Median IDO1 (cytoplasm)", "IDO1"),
                  ("median_dapi", "Median DAPI (nucleus)", "DAPI control"),
                  ("pct_ido1_pos", "% called IDO1+", "Binary call")]
        for ax, (colname, ylab, ttl) in zip(axes, panels):
            for s in SAMPLE_ORDER:
                d = bgrad.loc[bgrad["sample_id"] == s].sort_values("median_depth_um")
                if not len(d):
                    continue
                ax.plot(d["median_depth_um"], d[colname], linewidth=4.5,
                        marker=MARKER_OF[s], markersize=18, color=COLOR_OF[s],
                        alpha=0.9)
            ax.set_xlabel("Distance from tissue boundary (µm)")
            ax.set_ylabel(ylab)
            ax.set_title(ttl, fontsize=FONT_SIZE_TITLE - 8)
            style_axes(ax)
        section_legend(axes[0], loc="best")
        fig.suptitle("CD68-lineage macrophages by depth into the tissue\n"
                     "If IDO1 moves while DAPI stays flat, the gradient is real "
                     "signal rather than tissue quality", y=1.04,
                     fontsize=FONT_SIZE_TITLE - 4)
        save_fig(fig, "F19_boundary_distance_gradient")
        write_csv(bgrad, "25_boundary_distance_gradient.csv")
    else:
        print("    WARNING: no usable depth bins.")


# %% Cell 9 - F: spatial maps, all phenotypes
# =============================================================================

banner("F - SPATIAL MAPS")


def spatial_panel(ax, sid, show, sizes=None, bg=True):
    d = sections[sid]
    ph = d["pheno_plot"]
    if bg:
        m = ~np.isin(ph, show)
        ax.scatter(d["x_plot"][m], d["y_plot"][m], s=0.6,
                   color=BACKGROUND_COLOR, linewidths=0, rasterized=True)
    for p in show:
        m = ph == p
        if m.sum() == 0:
            continue
        ax.scatter(d["x_plot"][m], d["y_plot"][m],
                   s=(sizes or {}).get(p, 5),
                   color=PHENOTYPE_COLORS.get(p, "#000000"),
                   linewidths=0, rasterized=True)
    ax.set_title(f"{short_label(sid)}  ({d['condition']}, {d['scan']} "
                 f"pos{d['position']})\n{d['n_cells']:,} cells",
                 fontsize=FONT_SIZE_TITLE - 14)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_edgecolor(CONDITION_COLORS.get(d["condition"], "#999999"))
        sp.set_linewidth(6)


def spatial_figure(show, stem, title, sizes=None, bg=True):
    n = len(SAMPLE_ORDER)
    ncol = 3
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(15 * ncol, 13 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for k, s in enumerate(SAMPLE_ORDER):
        spatial_panel(axes[k], s, show, sizes=sizes, bg=bg)
    for k in range(n, len(axes)):
        axes[k].axis("off")
    handles = [Patch(facecolor=PHENOTYPE_COLORS.get(p, "#000"), label=p)
               for p in show]
    if bg:
        handles.append(Patch(facecolor=BACKGROUND_COLOR, label="all other cells"))
    fig.legend(handles=handles, loc="lower center",
               ncol=min(5, len(handles)), frameon=False,
               fontsize=FONT_SIZE_LEGEND - 10, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(title, y=1.02, fontsize=FONT_SIZE_TITLE - 4)
    save_fig(fig, stem)


spatial_figure(PHENO_USE, "F20_spatial_all_phenotypes",
               "All phenotypes, six retained sections\n"
               f"(downsampled to at most {MAX_POINTS_PER_SECTION:,} cells "
               "per section)",
               sizes={p: 4 for p in PHENO_USE}, bg=False)

spatial_figure([p for p in STRUCTURAL_PHENOTYPES if p in PHENO_USE],
               "F21_spatial_structural",
               "Tissue architecture: endothelium, epithelium and SMA-high stroma\n"
               'The "Other" class is SMA-positive and dim on every other marker',
               sizes={"Endothelial cells": 4, "Epithelial/Tumor cells": 4,
                      "Other": 4})

spatial_figure([p for p in IMMUNE_PHENOTYPES if p in PHENO_USE],
               "F22_spatial_immune",
               "Immune populations against structural background\n"
               "Grey = endothelium, epithelium and stroma",
               sizes={"CD68+IDO1+ Macrophages": 10, "CD68+IDO1- Macrophages": 7,
                      "CD163+ Macrophages": 10, "Neutrophils": 4,
                      "Helper T cells": 12, "CD4- T cells": 5, "Tregs": 22,
                      "B cells": 8, "Plasma cells": 5})

spatial_figure([p for p in ["B cells", "Plasma cells", "Helper T cells", "Tregs"]
                if p in PHENO_USE],
               "F23_spatial_lymphoid",
               "Lymphoid aggregate candidates\n"
               "Discrete B cell clusters with adjacent T cells would indicate "
               "BALT-like structures",
               sizes={"B cells": 12, "Plasma cells": 6, "Helper T cells": 14,
                      "Tregs": 26})


# %% Cell 10 - wrap up
# =============================================================================

banner("SUMMARY")
print(f"Sections retained : {len(SAMPLE_ORDER)}  ({', '.join(SAMPLE_ORDER)})")
print(f"Sections excluded : {', '.join(EXCLUDE_SECTIONS)}")
print(f"Cells retained    : {int(counts.sum().sum()):,}")
if np.isfinite(best_lo):
    print(f"IDO1 cutoff window: raw {np.expm1(best_lo):.2f} to "
          f"{np.expm1(best_hi):.2f}  (separated and stable)")

sub("Where to read each answer")
print("  A composition + IDO1 fraction -> F17, tables 20 / 21 / 22")
print("  B threshold back-out          -> table 23  (READ THE VERDICT ABOVE)")
print("  C threshold sensitivity       -> F18 centre panel, tables 24 / 26")
print("  D IDO1 distributions          -> F18 left and right panels")
print("  E boundary gradient           -> F19, table 25")
print("  F spatial maps                -> F20 all, F21 structural, F22 immune,")
print("                                   F23 lymphoid candidates")

sub("Hand-off to script 09")
print("  Table 26 records every input to the cutoff search and both windows.")
print("  Script 09 recomputes the same quantity from the cell assignment files")
print("  rather than the raw CSVs. Those two populations should be identical.")
print("  Diff table 26 against script 09's equivalent. If the windows differ,")
print("  the cause is the input population, not the algorithm, and the printed")
print("  per-section cell counts and cutoff-zero fractions will show where.")

sub("Decisions this run should let us close")
print("  1. Is the vendor IDO1 call reproducible as a simple intensity gate?")
print("  2. Over what cutoff range do the arms stay completely separated?")
print("  3. Is there a depth-into-tissue gradient in the retained sections?")
print("  4. Are there discrete lymphoid aggregates worth defining as structures?")
print("  5. Which statistic, p99 or section median, generated the comparability")
print("     numbers carried in the project notes? Table 22 answers this.")

banner("DONE")
sys.stdout = _tee.terminal
_tee.close()
