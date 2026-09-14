#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - FOCI DETECTION AUDIT
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 04b of the AKOYA analysis series. REVISION 2. READ-ONLY DIAGNOSTIC.

WHAT CHANGED IN REVISION 2, AND WHY

  1. THE DENSITY MAP IS NOW SCRIPT 04'S, NOT A RECONSTRUCTION OF IT.
     Revision 1 guessed at the grid and bandwidth and got both wrong: 20 um and
     50 um against script 04's actual 25 um and 75 um, on a tissue mask that
     did not fill holes where script 04 does. The result was peaks about 45
     percent too high and roughly twice as many foci per section, which is what
     under-smoothing on a finer grid produces. Every function below is now the
     same logic as script 04's make_grid, tissue_mask, density_map,
     find_local_maxima and detect_foci, with the parameters read off script 04
     rather than assumed.

     Revision 1 also invented a greedy separation step that SCRIPT 04 DOES NOT
     HAVE. Script 04 gets its peak separation entirely from the maximum_filter
     footprint of 2*sep_px+1, and every peak clearing the fold gate becomes a
     watershed marker. The "merged" rejection category is therefore removed. It
     was an artefact of the audit, not a property of the pipeline.

  2. THE CROSS-CHECK IS PAIRED AND EXACT.
     Table 35 carries background_density, peak_row and peak_col, so the
     reconstruction can be verified against script 04 focus by focus instead of
     by comparing set medians of differently sized sets. That comparison was the
     weakest part of revision 1: it could not separate "the map is wrong" from
     "the detector is more permissive", and it reported a reassuring median
     ratio of 1.05 while the one clean 1:1 section was off by 45 percent.
     Background density is a single number per section computed over the whole
     tissue mask, so it is checked FIRST and to four significant figures. If
     backgrounds match, the map matches.

  3. THE BACKGROUND STATISTIC IS TREATED AS A FINDING, NOT A SETTING.
     Section backgrounds under script 04's median-over-tissue-pixels run 326 to
     1,679 cells/mm2, a five-fold spread. At a 6x gate that means "a focus"
     is 1,957 cells/mm2 in 43106 and 10,073 in 36463, the latter more than three
     times the burden threshold of 3,000. A self-referencing gate asks lesions
     in a heavily involved section to clear a bar that the disease itself
     raised, so the most diseased sections are the ones where detection fails.
     36463 is that section, and it is the section where diffuse pathology was
     identified as missed. No fold reduction fixes this: at 4x the 36463 gate is
     still 6,715 and at 3x it is 5,037, both above burden density. Cell 6
     compares candidate background statistics and reports what each does to the
     cross-section spread.

  4. THE GATE IS A HYBRID, SWEPT IN THREE DIMENSIONS.
     A peak is kept when
         peak >= ABS_FLOOR  AND  ( peak >= fold * background  OR  peak >= ABS_RESCUE )
     ABS_FLOOR stops a low relative gate on a quiet treated background from
     admitting tissue at densities the burden definition calls normal lung.
     ABS_RESCUE admits genuinely dense tissue in a section whose own background
     is already pathological, which is the 36463 case. Both are applied
     identically to both arms. Setting ABS_FLOOR to 0 and ABS_RESCUE to infinity
     reproduces script 04 exactly, and that combination is in the grid as the
     control.

  5. COMPOSITION IS REPORTED IN THREE WAYS, AND NO COMPOSITION GATE IS IMPOSED.
     Revision 1 scored every candidate against the accepted UNTREATED foci. That
     is the right guard against circularity but it makes treated candidates
     score low by construction if residual treated foci genuinely differ, and
     the run showed exactly that: the eleven currently accepted treated foci sit
     at 0.46 cosine and 0.23 myeloid fraction against 0.78 and 0.73 untreated.
     So the cosine is now reported against BOTH an untreated reference and a
     within-arm reference, alongside the raw myeloid fraction, which depends on
     no reference at all and is the number to trust.

     Endothelial fraction is reported explicitly. The revision 1 run showed the
     endothelial population rising as the gate falls, which is the signature of
     admitting perivascular tissue rather than lesions, and it is the cleanest
     evidence that a low relative gate on a quiet background is wrong.

     Composition thresholds are SWEPT AND REPORTED, never applied. A myeloid
     floor that removes the 4x to 6x band also removes most currently accepted
     treated foci, so imposing one silently would rewrite Finding 1.

  6. THE EXPERT COUNTS ARE USED, AND ONLY WHERE THEY ARE REAL.
     Script 04 carries expert counts of 1, 2 and 2 for the treated sections and
     algorithmic targets for the untreated ones. Scoring against algorithmic
     targets is circular, so only the three treated sections contribute to the
     target error. This is the fix for the known weakness in script 04's sweep
     metric, where every top-ranked combination was fold 15 because 43109's
     algorithmic target dominated.

     Worth reading before anything else: at the CURRENT settings 43106 detects 3
     foci where the expert counted 1, and 43111 detects 1 where the expert
     counted 2. The gate is already wrong in both directions on the only
     sections with genuine ground truth. This is not a question of whether to
     loosen a correct gate.

WHAT THIS SCRIPT STILL DOES NOT DO

  It writes nothing into structures_rev4. It does not adopt a setting. The
  burden definition is untouched and appears only as a reference line. Any
  change goes into a script 04 revision 5 after this output has been read.

INPUTS
    structures_rev4/cell_assignments/*_cell_structures.csv
    structures_rev4/tables/35_foci_structures_relative.csv

OUTPUTS
    foci_audit/tables/   90_candidate_ledger, 91_background_statistics,
                         92_gate_grid, 93_reconstruction_crosscheck,
                         94_composition_references, 95_composition_gates,
                         96_expert_target_error, 00_foci_audit_report.txt
    foci_audit/figures/  F70_candidate_maps, F71_fold_vs_absolute_density,
                         F72_composition, F73_gate_grid, F74_background_statistic

USAGE
    conda activate sc_pre
    python AKOYA_04b_Foci_Audit.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================

IN_DIR = "/master/jlehle/WORKING/AKOYA/structures_rev4"
CELL_DIR = f"{IN_DIR}/cell_assignments"
FOCI_TABLE = f"{IN_DIR}/tables/35_foci_structures_relative.csv"
OUT_DIR = "/master/jlehle/WORKING/AKOYA/foci_audit"

USE_AGG = True

# ---- READ OFF SCRIPT 04. Do not change without changing script 04. ----------
GRID_UM = 25.0                  # script 04 Cell 1
DENSITY_BANDWIDTH_UM = 75.0     # script 04 Cell 1, used by both definitions
EDGE_MIN_WEIGHT = 0.25          # normalized convolution denominator floor
EDGE_CORRECTION = True
TISSUE_CLOSE_ITERATIONS = 2     # script 04 tissue_mask
TISSUE_FILL_HOLES = True        # script 04 tissue_mask fills holes
BACKGROUND_STAT = "median"      # script 04 detect_foci, median over tissue px

PEAK_SEPARATION_UM = 200.0      # -> sep_px, maximum_filter size 2*sep_px+1
FOCUS_FOLD_OVER_BACKGROUND = 6.0
FOCUS_HALF_MAX_FRACTION = 0.50
FOCUS_MIN_AREA_UM2 = 10000.0
FOCUS_MIN_CELLS = 40
FOCUS_FILL_HOLES = True
DEDUPLICATE_PLATEAU_PEAKS = True

BURDEN_THRESHOLD = 3000.0       # reference line only, never modified here

# ---- expert ground truth, from script 04 Cell 1 -----------------------------
# Only the treated sections are expert counts. The untreated targets are
# algorithmic, so scoring against them is circular and they are EXCLUDED from
# the target error. This is the fix for script 04's known sweep weakness.
EXPERT_COUNTS = {"G3_43106": 1, "G3_43111": 2, "G3_43118": 2,
                 "G4_31438": 7, "G4_36463": 5, "G4_43109": 19}
EXPERT_IS_REAL = {"G3_43106": True, "G3_43111": True, "G3_43118": True,
                  "G4_31438": False, "G4_36463": False, "G4_43109": False}

# ---- the audit --------------------------------------------------------------
PERMISSIVE_FOLD = 2.0           # candidate ledger detection floor

# The hybrid gate. keep = (peak >= floor) & ((peak >= fold*bg) | (peak >= rescue))
# fold 6 / floor 0 / rescue inf reproduces script 04 exactly and is the control.
FOLD_GRID = [6.0, 5.0, 4.0, 3.0]
ABS_FLOOR_GRID = [0.0, 1500.0, 2000.0, 2500.0]
ABS_RESCUE_GRID = [float("inf"), 3000.0]   # inf disables the rescue

# composition thresholds are REPORTED, not applied
MYELOID_FLOOR_GRID = [0.0, 0.10, 0.20, 0.30, 0.40]
ENDOTHELIAL_CEILING_GRID = [1.0, 0.30, 0.20, 0.15, 0.10]

MIN_CELLS_FOR_COMPOSITION = 20

# ---- phenotypes -------------------------------------------------------------
IDO1_POS = "CD68+IDO1+ Macrophages"
IDO1_NEG = "CD68+IDO1- Macrophages"
MYELOID_FOR_DETECTION = [IDO1_POS, IDO1_NEG, "CD163+ Macrophages", "Neutrophils"]
ENDOTHELIAL = "Endothelial cells"

PHENOTYPE_ORDER = [
    IDO1_POS, IDO1_NEG, "CD163+ Macrophages", "Neutrophils",
    "Helper T cells", "CD4- T cells", "Tregs", "B cells", "Plasma cells",
    "Endothelial cells", "Epithelial/Tumor cells", "Other",
]
PHENOTYPE_COLORS = {
    IDO1_POS: "#B2182B", IDO1_NEG: "#EF8A62", "CD163+ Macrophages": "#FDBE85",
    "Neutrophils": "#7B3294", "Helper T cells": "#1B7837",
    "CD4- T cells": "#7FBC41", "Tregs": "#00441B", "B cells": "#2166AC",
    "Plasma cells": "#67A9CF", "Endothelial cells": "#E7298A",
    "Epithelial/Tumor cells": "#66C2A5", "Other": "#A6761D",
}

CONDITION_ORDER = ["D1MT", "Untreated"]
CONDITION_COLORS = {"D1MT": "#2C7FB8", "Untreated": "#D95F02"}
CONDITION_LABELS = {"D1MT": "D1MT-treated", "Untreated": "Untreated"}
REFERENCE_ARM = "Untreated"

REASON_COLORS = {"accepted": "#1A9850", "fold": "#D73027",
                 "abs_floor": "#4575B4", "area": "#FDAE61",
                 "n_cells": "#984EA3"}

# ---- figures ----------------------------------------------------------------
FONT_SIZE_BASE = 28
FONT_SIZE_TITLE = 32
FONT_SIZE_TICK = 28
FONT_SIZE_LEGEND = 28
FONT_SIZE_PANEL = 34
AXIS_COLOR = "#333333"
TEXT_COLOR = "#000000"
GRID_COLOR = "#DDDDDD"
BACKGROUND_COLOR = "#CCCCCC"
DPI = 300
SAVE_PDF = True
SAVE_PNG = True
MAX_POINTS_MAP = 250000
RANDOM_SEED = 0


# %% Cell 2 - imports and helpers
# =============================================================================

import os
import sys
import gc
import glob
from datetime import datetime

import numpy as np
import pandas as pd
import scipy.ndimage as ndi

import matplotlib
if USE_AGG:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

try:
    from skimage.segmentation import watershed as _sk_watershed
    HAVE_SKIMAGE = True
except Exception:
    HAVE_SKIMAGE = False

pd.set_option("future.no_silent_downcasting", True)

plt.rcParams.update({
    "font.size": FONT_SIZE_BASE, "axes.titlesize": FONT_SIZE_TITLE,
    "axes.labelsize": FONT_SIZE_BASE, "xtick.labelsize": FONT_SIZE_TICK,
    "ytick.labelsize": FONT_SIZE_TICK, "legend.fontsize": FONT_SIZE_LEGEND,
    "figure.titlesize": FONT_SIZE_TITLE, "axes.edgecolor": AXIS_COLOR,
    "axes.labelcolor": TEXT_COLOR, "text.color": TEXT_COLOR,
    "xtick.color": AXIS_COLOR, "ytick.color": AXIS_COLOR,
    "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42,
})

FIG_DIR = os.path.join(OUT_DIR, "figures")
TAB_DIR = os.path.join(OUT_DIR, "tables")
for d in (FIG_DIR, TAB_DIR):
    os.makedirs(d, exist_ok=True)

rng = np.random.default_rng(RANDOM_SEED)

# Fold bands used by both the ledger summary and F72. Defined once here so the
# cells can be run out of order in Spyder.
FOLD_BANDS = [(6.0, np.inf, "6x and above"), (5.0, 6.0, "5x to 6x"),
              (4.0, 5.0, "4x to 5x"), (3.0, 4.0, "3x to 4x"),
              (PERMISSIVE_FOLD, 3.0, "below 3x")]
SEP_PX = max(1, int(round(PEAK_SEPARATION_UM / GRID_UM)))
PIXEL_AREA_MM2 = (GRID_UM * GRID_UM) / 1e6


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


def short_label(sid):
    return sid.split("_")[-1] if "_" in sid else sid


def save_fig(fig, stem):
    for flag, ext in ((SAVE_PDF, "pdf"), (SAVE_PNG, "png")):
        if not flag:
            continue
        p = os.path.join(FIG_DIR, f"{stem}.{ext}")
        fig.savefig(p, dpi=DPI, bbox_inches="tight")
        print(f"    wrote {p}")
    plt.close(fig); gc.collect()


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


def panel_letter(ax, letter, dx=-0.16, dy=1.07):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=FONT_SIZE_PANEL,
            fontweight="bold", va="top", ha="left")


_tee = Tee(os.path.join(TAB_DIR, "00_foci_audit_report.txt"))
sys.stdout = _tee

banner("AKOYA FOCI DETECTION AUDIT (script 04b revision 2, read-only)")
print(f"Run time : {datetime.now().isoformat(timespec='seconds')}")
print(f"Output   : {OUT_DIR}")
print("\nWrites nothing into structures_rev4. Adopts no setting. The burden")
print(f"definition stays at {BURDEN_THRESHOLD:,.0f} cells/mm2 applied identically")
print("to every section and appears here only as a reference line.")

sub("Reconstruction parameters, taken from script 04 Cell 1")
print(f"    grid                    {GRID_UM} um   ({PIXEL_AREA_MM2:.6f} mm2/px)")
print(f"    density bandwidth       {DENSITY_BANDWIDTH_UM} um "
      f"({DENSITY_BANDWIDTH_UM / GRID_UM:.2f} px)")
print(f"    edge min weight         {EDGE_MIN_WEIGHT}")
print(f"    tissue mask             binary_closing(iterations="
      f"{TISSUE_CLOSE_ITERATIONS}) then fill_holes={TISSUE_FILL_HOLES}")
print(f"    peak separation         {PEAK_SEPARATION_UM} um -> sep_px {SEP_PX}, "
      f"maximum_filter size {2 * SEP_PX + 1}")
print(f"    background statistic    {BACKGROUND_STAT} over tissue pixels")
print("\n    Revision 1 used 20 um and 50 um on an unfilled mask and invented a")
print("    greedy separation step script 04 does not have. That produced peaks")
print("    45 percent high and twice the foci. All of it is corrected here.")

if not HAVE_SKIMAGE:
    print("\n    WARNING: skimage unavailable. Script 04 uses watershed when it")
    print("    is present, so the EDT Voronoi fallback WILL differ from script")
    print("    04 and the cross-check will fail for that reason alone.")


# %% Cell 3 - load
# =============================================================================

banner("LOADING CELL ASSIGNMENTS")

paths = sorted(glob.glob(os.path.join(CELL_DIR, "*_cell_structures.csv")))
if not paths:
    print(f"ERROR: no cell assignment files in {CELL_DIR}")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

cells, COND_OF = {}, {}
for p in paths:
    sid = os.path.basename(p).replace("_cell_structures.csv", "")
    try:
        d = pd.read_csv(p, low_memory=False)
    except Exception as e:
        print(f"    ERROR reading {sid}: {e}")
        continue
    miss = [c for c in ("x", "y", "pheno") if c not in d.columns]
    if miss:
        print(f"    ERROR: {sid} missing {miss}, skipped")
        continue
    d = d.loc[np.isfinite(d["x"]) & np.isfinite(d["y"])].reset_index(drop=True)
    cells[sid] = d
    COND_OF[sid] = d["condition"].iloc[0] if "condition" in d.columns else "?"
    n_acc = int(d["focus_id"].gt(0).sum()) if "focus_id" in d.columns else -1
    print(f"    {sid:<14} {COND_OF[sid]:<10} {len(d):>9,} cells   "
          f"{d['focus_id'].max() if 'focus_id' in d.columns else 0:>3} foci in rev4")

SAMPLE_ORDER = sorted(cells, key=lambda s: (
    CONDITION_ORDER.index(COND_OF[s]) if COND_OF[s] in CONDITION_ORDER else 9, s))

PH_USE = [p for p in PHENOTYPE_ORDER]
myeloid_use = [p for p in MYELOID_FOR_DETECTION]

foci_ref = None
if os.path.exists(FOCI_TABLE):
    foci_ref = pd.read_csv(FOCI_TABLE)
    print(f"\n    loaded table 35  ({foci_ref.shape[0]} x {foci_ref.shape[1]})")
else:
    print(f"\n    WARNING: table 35 missing. The cross-check cannot run and the")
    print("    reconstruction is UNVERIFIED.")

sub("Ground truth available")
print(f"    {'section':<14}{'expert':>8}{'real?':>8}{'in table 35':>13}")
for s in SAMPLE_ORDER:
    n35 = (int((foci_ref["sample_id"].astype(str) == s).sum())
           if foci_ref is not None and "sample_id" in foci_ref.columns else -1)
    print(f"    {s:<14}{EXPERT_COUNTS.get(s, -1):>8}"
          f"{str(EXPERT_IS_REAL.get(s, False)):>8}{n35:>13}")
print("\n    Only the treated rows are expert counts. Note that the CURRENT")
print("    settings already miss in both directions on those three sections,")
print("    so this is not a question of whether to loosen a correct gate.")


# %% Cell 4 - script 04's density and detection, replicated
# =============================================================================
# Every function here is the same logic as the corresponding function in
# AKOYA_04_Structures.py revision 4. Differences would invalidate the audit, so
# they are kept deliberately literal rather than tidied.

def make_grid(x, y, pitch):
    x0 = float(np.floor(x.min() / pitch) * pitch)
    y0 = float(np.floor(y.min() / pitch) * pitch)
    gx = ((x - x0) / pitch).astype(int) + 1
    gy = ((y - y0) / pitch).astype(int) + 1
    return x0, y0, gx, gy, (int(gy.max()) + 2, int(gx.max()) + 2)


def tissue_mask(gx, gy, shape):
    occ = np.zeros(shape, dtype=bool)
    occ[gy, gx] = True
    occ = ndi.binary_closing(occ, iterations=TISSUE_CLOSE_ITERATIONS)
    if TISSUE_FILL_HOLES:
        occ = ndi.binary_fill_holes(occ)
    return occ


def density_map(gx, gy, shape, mask, bandwidth_um, pitch, occ=None,
                edge_correct=True, min_weight=EDGE_MIN_WEIGHT):
    counts = np.zeros(shape, dtype=float)
    np.add.at(counts, (gy[mask], gx[mask]), 1.0)
    sigma = bandwidth_um / pitch
    sm = ndi.gaussian_filter(counts, sigma=sigma, mode="constant")
    if edge_correct and occ is not None:
        w = ndi.gaussian_filter(occ.astype(float), sigma=sigma, mode="constant")
        sm = sm / np.maximum(w, min_weight)
    return sm / ((pitch * pitch) / 1e6)


def partition_by_markers(dens, markers, territory):
    if HAVE_SKIMAGE:
        return _sk_watershed(-dens, markers=markers, mask=territory), "watershed"
    _, inds = ndi.distance_transform_edt(markers == 0, return_indices=True)
    part = markers[inds[0], inds[1]]
    return np.where(territory, part, 0), "edt_voronoi"


def find_local_maxima(dens, occ, sep_px, deduplicate=True):
    mx = ndi.maximum_filter(dens, size=2 * sep_px + 1, mode="constant")
    ispeak = (dens == mx) & occ & (dens > 0)
    n_raw = int(ispeak.sum())
    if not deduplicate or n_raw == 0:
        py, px = np.nonzero(ispeak)
        return py, px, dens[py, px], n_raw, 0
    lab, n = ndi.label(ispeak, structure=np.ones((3, 3), dtype=int))
    coms = ndi.center_of_mass(ispeak, lab, list(range(1, n + 1)))
    objs = ndi.find_objects(lab)
    ys, xs = [], []
    for k, (com, sl) in enumerate(zip(coms, objs), start=1):
        sub_lab = lab[sl] == k
        yy, xx = np.nonzero(sub_lab)
        yy = yy + sl[0].start
        xx = xx + sl[1].start
        d2 = (yy - com[0]) ** 2 + (xx - com[1]) ** 2
        j = int(np.argmin(d2))
        ys.append(int(yy[j])); xs.append(int(xx[j]))
    py = np.asarray(ys, dtype=int)
    px = np.asarray(xs, dtype=int)
    return py, px, dens[py, px], n_raw, n_raw - len(py)


def background_of(dens, occ, stat=BACKGROUND_STAT):
    """
    Script 04 offers median and p75. Revision 2 adds low percentiles because
    the median of a heavily involved section IS pathology, which is the whole
    36463 problem. Adding an option changes nothing unless it is selected.
    """
    v = dens[occ]
    v = v[np.isfinite(v)]
    if not v.size:
        return np.nan
    if stat == "median":
        bg = float(np.median(v))
    elif stat == "p75":
        bg = float(np.percentile(v, 75))
    elif stat.startswith("p") and stat[1:].replace(".", "").isdigit():
        bg = float(np.percentile(v, float(stat[1:])))
    elif stat == "mode_kde":
        # densest value of the density distribution, which tracks uninvolved
        # parenchyma rather than the middle of the whole section. The bottom
        # percentile is dropped first: binary_fill_holes puts pad and lumen
        # pixels inside the tissue mask at near-zero density, and on a smoke
        # test those alone captured the modal bin and returned a background of
        # 1.0 cells/mm2.
        vv = v[v > 0]
        if vv.size:
            vv = vv[vv >= np.percentile(vv, 1.0)]
        if vv.size < 10:
            bg = float(np.median(v))
        else:
            lv = np.log10(vv)
            h, e = np.histogram(lv, bins=200)
            k = int(np.argmax(h))
            bg = float(10 ** (0.5 * (e[k] + e[k + 1])))
    else:
        bg = float(np.median(v))
    if not np.isfinite(bg) or bg <= 0:
        bg = float(np.mean(v)) or 1.0
    return bg


def detect_foci(dens, occ, gx, gy, bg, fold_min, half_frac, min_area, min_cells,
                sep_px, abs_floor=0.0, abs_rescue=float("inf"),
                fill_holes=FOCUS_FILL_HOLES, deduplicate=True):
    """
    Script 04's detect_foci with the hybrid keep-rule substituted for the plain
    fold gate. With abs_floor = 0 and abs_rescue = inf the keep-rule collapses
    to pv >= fold_min * bg, which is script 04 exactly.

    Returns (labels, records, method, diagnostics).
    """
    diag = {"n_peaks_raw": 0, "n_peaks_collapsed": 0, "n_peaks_kept": 0,
            "n_contested_pixels": 0, "n_rejected_area": 0,
            "n_rejected_cells": 0}
    if not occ.any():
        return np.zeros_like(dens, dtype=int), [], "none", diag

    py, px, pv, n_raw, n_coll = find_local_maxima(dens, occ, sep_px, deduplicate)
    diag["n_peaks_raw"] = n_raw
    diag["n_peaks_collapsed"] = n_coll

    keep = (pv >= abs_floor) & ((pv >= fold_min * bg) | (pv >= abs_rescue))
    py, px, pv = py[keep], px[keep], pv[keep]
    diag["n_peaks_kept"] = int(len(pv))
    if not len(pv):
        return np.zeros_like(dens, dtype=int), [], "none", diag

    order = np.argsort(pv)[::-1]
    py, px, pv = py[order], px[order], pv[order]

    markers = np.zeros(dens.shape, dtype=int)
    markers[py, px] = np.arange(1, len(pv) + 1)
    territory = (dens >= half_frac * pv.min()) & occ
    part, method = partition_by_markers(dens, markers, territory)

    labels = np.zeros(dens.shape, dtype=int)
    recs = []
    next_id = 1
    for i, (yy, xx, vv) in enumerate(zip(py, px, pv), start=1):
        region = (part == i) & (dens >= half_frac * vv)
        if not region[yy, xx]:
            continue
        cl, cn = ndi.label(region)
        if cn == 0:
            continue
        region = cl == cl[yy, xx]
        if fill_holes:
            region = ndi.binary_fill_holes(region)
        assigned = region & (labels == 0)
        diag["n_contested_pixels"] += int(region.sum() - assigned.sum())
        if not assigned[yy, xx]:
            continue
        cl2, cn2 = ndi.label(assigned)
        if cn2 == 0:
            continue
        assigned = cl2 == cl2[yy, xx]

        area = assigned.sum() * GRID_UM * GRID_UM
        ncell = int(assigned[gy, gx].sum())
        if area < min_area:
            diag["n_rejected_area"] += 1
            continue
        if ncell < min_cells:
            diag["n_rejected_cells"] += 1
            continue

        dt = ndi.distance_transform_edt(assigned)
        r_in = float(dt.max()) * GRID_UM
        r_eq = float(np.sqrt(area / np.pi))
        labels[assigned] = next_id
        recs.append({
            "focus_id": next_id, "peak_density": float(vv),
            "background_density": float(bg),
            "fold_over_background": float(vv / bg) if bg > 0 else np.nan,
            "peak_row": int(yy), "peak_col": int(xx),
            "area_um2": float(area), "equiv_radius_um": r_eq,
            "max_inscribed_radius_um": r_in,
            "shape_ratio": (r_eq / r_in) if r_in > 0 else np.nan,
            "n_cells": ncell})
        next_id += 1
    return labels, recs, method, diag


def composition_of(labels, gx, gy, pheno, k):
    sel = labels[gy, gx] == k
    n = int(sel.sum())
    out = {"n_cells_in_region": n}
    if not n:
        for p in PH_USE:
            out[f"frac_{p}"] = np.nan
        out["frac_myeloid"] = np.nan
        out["frac_endothelial"] = np.nan
        return out
    vc = pd.Series(pheno[sel]).value_counts()
    for p in PH_USE:
        out[f"frac_{p}"] = float(vc.get(p, 0)) / n
    out["frac_myeloid"] = float(sum(vc.get(p, 0) for p in myeloid_use)) / n
    out["frac_endothelial"] = float(vc.get(ENDOTHELIAL, 0)) / n
    return out


banner("DENSITY MAPS")

grids, dens_of, bg_of = {}, {}, {}
for sid in SAMPLE_ORDER:
    d = cells[sid]
    x = d["x"].to_numpy(dtype=float); y = d["y"].to_numpy(dtype=float)
    x0, y0, gx, gy, shape = make_grid(x, y, GRID_UM)
    occ = tissue_mask(gx, gy, shape)
    mye = d["pheno"].isin(myeloid_use).to_numpy()
    dens = density_map(gx, gy, shape, mye, DENSITY_BANDWIDTH_UM, GRID_UM,
                       occ=occ, edge_correct=EDGE_CORRECTION)
    bg = background_of(dens, occ, BACKGROUND_STAT)
    grids[sid] = dict(x0=x0, y0=y0, gx=gx, gy=gy, shape=shape, occ=occ,
                      pheno=d["pheno"].to_numpy())
    dens_of[sid] = dens
    bg_of[sid] = bg
    print(f"    {sid:<14} tissue {occ.sum() * GRID_UM * GRID_UM / 1e6:>6.1f} mm2   "
          f"bg {bg:>8.2f} /mm2   peak {float(dens[occ].max()):>9.1f}   "
          f"({float(dens[occ].max()) / bg:.1f}x)")


# %% Cell 5 - paired cross-check against table 35
# =============================================================================

banner("CROSS-CHECK AGAINST SCRIPT 04 (read this before anything else)")

print("Two levels. Background density is a single number per section computed")
print("over the whole tissue mask, so it tests the MAP alone. If backgrounds")
print("match to four significant figures the map is right and any remaining")
print("difference is in detection. Foci are then matched by peak_row and")
print("peak_col, which is exact rather than a comparison of set medians.\n")

verdict = "unverified"
cc_rows = []
if foci_ref is not None and "sample_id" in foci_ref.columns:
    labels_ctrl, recs_ctrl = {}, {}
    for sid in SAMPLE_ORDER:
        g = grids[sid]
        lab, recs, method, diag = detect_foci(
            dens_of[sid], g["occ"], g["gx"], g["gy"], bg_of[sid],
            FOCUS_FOLD_OVER_BACKGROUND, FOCUS_HALF_MAX_FRACTION,
            FOCUS_MIN_AREA_UM2, FOCUS_MIN_CELLS, SEP_PX,
            abs_floor=0.0, abs_rescue=float("inf"),
            deduplicate=DEDUPLICATE_PLATEAU_PEAKS)
        labels_ctrl[sid] = lab
        recs_ctrl[sid] = pd.DataFrame(recs)
        ref = foci_ref.loc[foci_ref["sample_id"].astype(str) == sid]
        bg35 = (float(ref["background_density"].iloc[0])
                if len(ref) and "background_density" in ref.columns else np.nan)
        rec = dict(sample_id=sid, condition=COND_OF[sid],
                   n_foci_table35=int(len(ref)), n_foci_here=int(len(recs)),
                   background_table35=bg35, background_here=bg_of[sid],
                   background_ratio=(bg_of[sid] / bg35
                                     if np.isfinite(bg35) and bg35 else np.nan),
                   partition_method=method,
                   n_peaks_raw=diag["n_peaks_raw"],
                   n_peaks_collapsed=diag["n_peaks_collapsed"],
                   n_peaks_kept=diag["n_peaks_kept"])
        # paired match on the peak pixel
        matched = 0
        ratios = []
        if len(ref) and len(recs) and {"peak_row", "peak_col"} <= set(ref.columns):
            here = pd.DataFrame(recs)
            for _, r in ref.iterrows():
                hit = here.loc[(here["peak_row"] == r["peak_row"]) &
                               (here["peak_col"] == r["peak_col"])]
                if len(hit):
                    matched += 1
                    if "peak_density" in r and float(r["peak_density"]):
                        ratios.append(float(hit["peak_density"].iloc[0])
                                      / float(r["peak_density"]))
        rec["n_foci_matched_on_peak_pixel"] = matched
        rec["peak_density_ratio_matched"] = (float(np.median(ratios))
                                             if ratios else np.nan)
        cc_rows.append(rec)

    cross = pd.DataFrame(cc_rows)
    write_csv(cross, "93_reconstruction_crosscheck.csv")
    with pd.option_context("display.width", 260):
        print(cross[["sample_id", "condition", "background_table35",
                     "background_here", "background_ratio", "n_foci_table35",
                     "n_foci_here", "n_foci_matched_on_peak_pixel",
                     "peak_density_ratio_matched"]].to_string(index=False))

    br = cross["background_ratio"].to_numpy(dtype=float)
    br = br[np.isfinite(br)]
    bg_ok = bool(len(br) and np.all(np.abs(br - 1.0) < 0.001))
    cnt_ok = bool((cross["n_foci_table35"] == cross["n_foci_here"]).all())
    match_ok = bool((cross["n_foci_matched_on_peak_pixel"]
                     == cross["n_foci_table35"]).all())

    sub("VERDICT")
    if bg_ok and cnt_ok and match_ok:
        verdict = "pass"
        print("    PASS. Backgrounds agree to better than 0.1 percent, focus")
        print("    counts match, and every table 35 focus is recovered at the")
        print("    same peak pixel. The reconstruction IS script 04.")
    elif bg_ok:
        verdict = "map_ok_detection_differs"
        print("    MAP OK, DETECTION DIFFERS. Backgrounds match, so the density")
        print("    map is right. The difference is downstream of it. Check")
        print("    partition_method above: script 04 uses watershed, and an")
        print("    edt_voronoi fallback will not reproduce it. Then check")
        print("    FOCUS_FILL_HOLES and DEDUPLICATE_PLATEAU_PEAKS.")
    else:
        verdict = "fail"
        print("    FAIL. Backgrounds do not match, so the density map itself is")
        print("    wrong and nothing below is usable.")
        worst = cross.loc[cross["background_ratio"].sub(1).abs().idxmax()]
        print(f"    Worst: {worst['sample_id']} ratio "
              f"{worst['background_ratio']:.4f}")
        print("    Ratio above 1 means under-smoothing: raise")
        print("    DENSITY_BANDWIDTH_UM or lower GRID_UM. Ratio below 1 means")
        print("    over-smoothing. If the ratio is near 1 but not within")
        print("    tolerance, EDGE_MIN_WEIGHT or the tissue mask settings are")
        print("    the remaining candidates.")
else:
    print("    Cross-check not run.")


# %% Cell 6 - background statistic diagnostic
# =============================================================================

banner("THE BACKGROUND STATISTIC IS THE 36463 PROBLEM")

print("The fold gate is relative to each section's OWN background, so what a")
print("focus MEANS in absolute terms changes section to section. In a heavily")
print("involved section the median over tissue pixels is itself pathology, and")
print("the gate then asks lesions to clear a bar the disease raised. That is")
print("why detection fails hardest in the most diseased sections, which is the")
print("opposite of what anyone wants.\n")

bg_stats = ["median", "p75", "p50", "p25", "p10", "mode_kde"]
rows = []
for sid in SAMPLE_ORDER:
    occ = grids[sid]["occ"]
    r = dict(sample_id=sid, condition=COND_OF[sid])
    for st in bg_stats:
        b = background_of(dens_of[sid], occ, st)
        r[f"bg_{st}"] = b
        r[f"gate6x_{st}"] = b * FOCUS_FOLD_OVER_BACKGROUND
    r["peak_density"] = float(dens_of[sid][occ].max())
    r["burden_threshold"] = BURDEN_THRESHOLD
    rows.append(r)
bgt = pd.DataFrame(rows)
write_csv(bgt, "91_background_statistics.csv")

sub("Background under each statistic (cells/mm2)")
with pd.option_context("display.width", 260):
    print(bgt[["sample_id", "condition"] + [f"bg_{s}" for s in bg_stats]]
          .to_string(index=False))

sub("Cross-section spread, max over min. Lower is a more comparable gate.")
for st in bg_stats:
    v = bgt[f"bg_{st}"].to_numpy(dtype=float)
    v = v[np.isfinite(v) & (v > 0)]
    print(f"    {st:<10} spread {v.max() / v.min():>6.2f}x   "
          f"range {v.min():>8.1f} to {v.max():>8.1f}")

sub(f"What a {FOCUS_FOLD_OVER_BACKGROUND:g}x gate means in absolute cells/mm2")
show = ["sample_id", "condition"] + [f"gate6x_{s}" for s in bg_stats]
with pd.option_context("display.width", 260):
    print(bgt[show].to_string(index=False))
print(f"\n    Burden threshold is {BURDEN_THRESHOLD:,.0f}. Any gate value far")
print("    ABOVE it is asking a lesion to be denser than the burden definition")
print("    requires, which is how genuine pathology gets missed. Any gate value")
print("    far BELOW it is admitting tissue the burden definition calls normal.")
print("    A statistic whose gate values straddle the burden threshold tightly")
print("    across all six sections is the one to prefer.")


# %% Cell 7 - the candidate ledger
# =============================================================================

banner("CANDIDATE LEDGER (permissive detection)")

print(f"Detection at fold {PERMISSIVE_FOLD:g} so nothing the pipeline could ever")
print("have seen is invisible. Each candidate carries absolute peak density,")
print("fold, composition, and the gate outcome at the CURRENT settings.\n")

ledger_frames = []
for sid in SAMPLE_ORDER:
    g = grids[sid]
    lab, recs, method, diag = detect_foci(
        dens_of[sid], g["occ"], g["gx"], g["gy"], bg_of[sid],
        PERMISSIVE_FOLD, FOCUS_HALF_MAX_FRACTION, FOCUS_MIN_AREA_UM2,
        FOCUS_MIN_CELLS, SEP_PX, deduplicate=DEDUPLICATE_PLATEAU_PEAKS)
    if not recs:
        continue
    df = pd.DataFrame(recs)
    comp = pd.DataFrame([composition_of(lab, g["gx"], g["gy"], g["pheno"],
                                        int(k)) for k in df["focus_id"]])
    df = pd.concat([df.reset_index(drop=True), comp], axis=1)
    df.insert(0, "sample_id", sid)
    df.insert(1, "condition", COND_OF[sid])
    ledger_frames.append(df)
    print(f"    {sid:<14} peaks raw {diag['n_peaks_raw']:>6}  collapsed "
          f"{diag['n_peaks_collapsed']:>5}  kept {diag['n_peaks_kept']:>4}  "
          f"regions {len(recs):>4}  dropped area/cells "
          f"{diag['n_rejected_area']}/{diag['n_rejected_cells']}")

ledger = (pd.concat(ledger_frames, ignore_index=True) if ledger_frames
          else pd.DataFrame())
if not len(ledger):
    print("\nERROR: no candidates. Check the cross-check verdict.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

ledger["accepted_current"] = ledger["fold_over_background"] >= FOCUS_FOLD_OVER_BACKGROUND
ledger["rejection_reason"] = np.where(ledger["accepted_current"], "accepted", "fold")

# ---- composition references, both of them ----------------------------------
frac_cols = [f"frac_{p}" for p in PH_USE]


def cosine_to(v, c):
    v = np.asarray(v, float); c = np.asarray(c, float)
    ok = np.isfinite(v) & np.isfinite(c)
    if not ok.any():
        return np.nan
    v, c = v[ok], c[ok]
    nv, nc = np.linalg.norm(v), np.linalg.norm(c)
    return float(v @ c / (nv * nc)) if nv > 0 and nc > 0 else np.nan


ref_rows = []
centroids = {}
for arm in CONDITION_ORDER + ["all"]:
    sel = ledger.loc[ledger["accepted_current"]]
    if arm != "all":
        sel = sel.loc[sel["condition"] == arm]
    if not len(sel):
        continue
    c = sel[frac_cols].mean(axis=0)
    centroids[arm] = c.to_numpy(dtype=float)
    for p in PH_USE:
        ref_rows.append(dict(reference=arm, phenotype=p,
                             mean_fraction=float(c[f"frac_{p}"]),
                             n_foci=int(len(sel))))
refs = pd.DataFrame(ref_rows)
write_csv(refs, "94_composition_references.csv")

F = ledger[frac_cols].to_numpy(dtype=float)
if REFERENCE_ARM in centroids:
    ledger["cosine_to_untreated_foci"] = [cosine_to(r, centroids[REFERENCE_ARM])
                                          for r in F]
ledger["cosine_to_own_arm_foci"] = [
    cosine_to(r, centroids.get(c, centroids.get("all")))
    for r, c in zip(F, ledger["condition"])]
for col in ("cosine_to_untreated_foci", "cosine_to_own_arm_foci"):
    if col in ledger.columns:
        ledger.loc[ledger["n_cells"] < MIN_CELLS_FOR_COMPOSITION, col] = np.nan

write_csv(ledger, "90_candidate_ledger.csv")

sub("Composition by arm and gate outcome. frac_myeloid needs no reference.")
grp = (ledger.groupby(["condition", "accepted_current"])
       .agg(n=("focus_id", "size"),
            median_peak=("peak_density", "median"),
            median_fold=("fold_over_background", "median"),
            median_myeloid=("frac_myeloid", "median"),
            median_endothelial=("frac_endothelial", "median"),
            median_cos_untreated=("cosine_to_untreated_foci", "median"),
            median_cos_own_arm=("cosine_to_own_arm_foci", "median")))
with pd.option_context("display.width", 260):
    print(grp.to_string())
print("\n    Read frac_myeloid first. It depends on no reference, so if the")
print("    treated accepted foci are already far below the untreated ones then")
print("    the untreated cosine is scoring treated tissue against a yardstick")
print("    it was never going to meet, and the within-arm column is the fair")
print("    comparison for asking whether a REJECTED treated region looks like")
print("    an ACCEPTED treated one.")

sub("Endothelial fraction by fold band, the vascular admission check")
for lo, hi, name in FOLD_BANDS:
    m = (ledger["fold_over_background"] >= lo) & (ledger["fold_over_background"] < hi)
    for c in CONDITION_ORDER:
        s = ledger.loc[m & (ledger["condition"] == c)]
        if not len(s):
            continue
        print(f"    {name:<14} {c:<10} n={len(s):>4}  myeloid "
              f"{s['frac_myeloid'].median():.3f}  endothelial "
              f"{s['frac_endothelial'].median():.3f}  peak "
              f"{s['peak_density'].median():>8.0f}")
print("\n    A rising endothelial fraction as the band falls is the signature")
print("    of admitting perivascular tissue rather than lesions.")


# %% Cell 8 - the hybrid gate grid
# =============================================================================

banner("HYBRID GATE GRID")

print("keep = (peak >= ABS_FLOOR) and ((peak >= fold * background) or")
print("       (peak >= ABS_RESCUE))")
print()
print("ABS_FLOOR blocks a low relative gate on a quiet background from")
print("admitting normal lung. ABS_RESCUE admits dense tissue in a section whose")
print("own background is already pathological, which no fold value can reach.")
print("Both apply identically to both arms. fold 6 / floor 0 / rescue inf is")
print("script 04 exactly and is the control row.")
print()
n_combo = len(FOLD_GRID) * len(ABS_FLOOR_GRID) * len(ABS_RESCUE_GRID)
print(f"    {n_combo} combinations x {len(SAMPLE_ORDER)} sections. Detection is")
print("    re-run for every one, because the surviving peak set changes the")
print("    watershed territory and therefore the partition.\n")

grid_rows = []
for fold in FOLD_GRID:
    for floor in ABS_FLOOR_GRID:
        for rescue in ABS_RESCUE_GRID:
            for sid in SAMPLE_ORDER:
                g = grids[sid]
                lab, recs, _, diag = detect_foci(
                    dens_of[sid], g["occ"], g["gx"], g["gy"], bg_of[sid],
                    fold, FOCUS_HALF_MAX_FRACTION, FOCUS_MIN_AREA_UM2,
                    FOCUS_MIN_CELLS, SEP_PX, abs_floor=floor,
                    abs_rescue=rescue,
                    deduplicate=DEDUPLICATE_PLATEAU_PEAKS)
                if recs:
                    df = pd.DataFrame(recs)
                    comp = pd.DataFrame([
                        composition_of(lab, g["gx"], g["gy"], g["pheno"], int(k))
                        for k in df["focus_id"]])
                    df = pd.concat([df.reset_index(drop=True), comp], axis=1)
                else:
                    df = pd.DataFrame()
                grid_rows.append(dict(
                    fold=fold, abs_floor=floor, abs_rescue=rescue,
                    sample_id=sid, condition=COND_OF[sid],
                    background=bg_of[sid], n_foci=len(df),
                    expert=EXPERT_COUNTS.get(sid, np.nan),
                    expert_is_real=EXPERT_IS_REAL.get(sid, False),
                    min_peak=(float(df["peak_density"].min()) if len(df) else np.nan),
                    median_peak=(float(df["peak_density"].median()) if len(df) else np.nan),
                    median_myeloid=(float(df["frac_myeloid"].median()) if len(df) else np.nan),
                    median_endothelial=(float(df["frac_endothelial"].median()) if len(df) else np.nan),
                    total_area_mm2=(float(df["area_um2"].sum()) / 1e6 if len(df) else 0.0)))
            print(f"    fold {fold:<4g} floor {floor:>6.0f} rescue "
                  f"{rescue if np.isfinite(rescue) else 'none':>6}  done")

gridf = pd.DataFrame(grid_rows)
write_csv(gridf, "92_gate_grid.csv")

# ---- expert target error, treated sections only -----------------------------
real = gridf.loc[gridf["expert_is_real"]].copy()
real["abs_error"] = (real["n_foci"] - real["expert"]).abs()
err = (real.groupby(["fold", "abs_floor", "abs_rescue"])
       .agg(total_abs_error=("abs_error", "sum"),
            max_abs_error=("abs_error", "max"))
       .reset_index())
# guards that are independent of the expert counts
guard = (gridf.groupby(["fold", "abs_floor", "abs_rescue"])
         .agg(min_peak_any=("min_peak", "min"),
              median_myeloid=("median_myeloid", "median"),
              median_endothelial=("median_endothelial", "median"),
              n_foci_total=("n_foci", "sum")).reset_index())
score = err.merge(guard, on=["fold", "abs_floor", "abs_rescue"])
score["admits_below_burden"] = score["min_peak_any"] < BURDEN_THRESHOLD
write_csv(score, "96_expert_target_error.csv")

sub("Ranked by agreement with the three EXPERT-counted treated sections")
print("    The untreated targets are algorithmic and are excluded, which is the")
print("    fix for script 04's known sweep weakness.\n")
show = score.sort_values(["total_abs_error", "max_abs_error"]).head(15)
with pd.option_context("display.width", 260):
    print(show.to_string(index=False))

sub("Control row, script 04 as it stands")
ctrl = score.loc[(score["fold"] == FOCUS_FOLD_OVER_BACKGROUND)
                 & (score["abs_floor"] == 0.0)
                 & (~np.isfinite(score["abs_rescue"]))]
with pd.option_context("display.width", 260):
    print(ctrl.to_string(index=False))

sub("Per-section foci counts, best few combinations against the control")
best = show.head(5)[["fold", "abs_floor", "abs_rescue"]].values.tolist()
best = [tuple(b) for b in best]
ctrl_key = (FOCUS_FOLD_OVER_BACKGROUND, 0.0, float("inf"))
for key in [ctrl_key] + [b for b in best if b != ctrl_key]:
    f, fl, rs = key
    sel = gridf.loc[(gridf["fold"] == f) & (gridf["abs_floor"] == fl)
                    & ((gridf["abs_rescue"] == rs)
                       if np.isfinite(rs) else ~np.isfinite(gridf["abs_rescue"]))]
    if not len(sel):
        continue
    tag = "CONTROL" if key == ctrl_key else "       "
    line = "  ".join(f"{short_label(r.sample_id)}={r.n_foci}"
                     f"({EXPERT_COUNTS.get(r.sample_id, '?')})"
                     for r in sel.itertuples())
    print(f"    {tag} fold {f:g} floor {fl:.0f} rescue "
          f"{rs if np.isfinite(rs) else 'none'}   {line}")
print("\n    Counts are shown as detected(expert). Only the three G3 sections")
print("    have genuine expert counts.")


# %% Cell 9 - composition gates, reported not applied
# =============================================================================

banner("COMPOSITION GATES: WHAT THEY WOULD COST")

print("A myeloid floor or an endothelial ceiling is an obvious way to exclude")
print("perivascular tissue. It is NOT applied here, because the currently")
print("accepted treated foci sit at a low myeloid fraction themselves, so a")
print("floor tuned to exclude the rejected band would also delete most of")
print("Finding 1's residual foci. The cost of each threshold is reported so")
print("that trade is made explicitly rather than by accident.\n")

acc = ledger.loc[ledger["accepted_current"]]
gate_rows = []
for mf in MYELOID_FLOOR_GRID:
    for ec in ENDOTHELIAL_CEILING_GRID:
        keep_acc = acc.loc[(acc["frac_myeloid"] >= mf)
                           & (acc["frac_endothelial"] <= ec)]
        band = ledger.loc[(~ledger["accepted_current"])
                          & (ledger["fold_over_background"] >= 4.0)]
        keep_band = band.loc[(band["frac_myeloid"] >= mf)
                             & (band["frac_endothelial"] <= ec)]
        r = dict(myeloid_floor=mf, endothelial_ceiling=ec)
        for c in CONDITION_ORDER:
            r[f"accepted_kept_{c}"] = int((keep_acc["condition"] == c).sum())
            r[f"accepted_total_{c}"] = int((acc["condition"] == c).sum())
            r[f"band4to6_kept_{c}"] = int((keep_band["condition"] == c).sum())
            r[f"band4to6_total_{c}"] = int((band["condition"] == c).sum())
        gate_rows.append(r)
gates = pd.DataFrame(gate_rows)
write_csv(gates, "95_composition_gates.csv")

sub("Endothelial ceiling alone, myeloid floor at 0")
sel = gates.loc[gates["myeloid_floor"] == 0.0]
with pd.option_context("display.width", 260):
    print(sel[["endothelial_ceiling"]
              + [f"accepted_kept_{c}" for c in CONDITION_ORDER]
              + [f"band4to6_kept_{c}" for c in CONDITION_ORDER]].to_string(index=False))

sub("Myeloid floor alone, endothelial ceiling at 1")
sel = gates.loc[gates["endothelial_ceiling"] == 1.0]
with pd.option_context("display.width", 260):
    print(sel[["myeloid_floor"]
              + [f"accepted_kept_{c}" for c in CONDITION_ORDER]
              + [f"band4to6_kept_{c}" for c in CONDITION_ORDER]].to_string(index=False))
print("\n    A good gate keeps the accepted columns intact while cutting the")
print("    band columns. A gate that cuts both is removing the finding along")
print("    with the noise, and an endothelial ceiling is more likely to manage")
print("    that than a myeloid floor because it targets the specific tissue")
print("    the revision 1 run implicated.")


# %% Cell 10 - figures
# =============================================================================

banner("FIGURES")

# ---- F74 background statistic ----------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(30, 13))
ax = axes[0]
xs = np.arange(len(bg_stats))
for sid in SAMPLE_ORDER:
    r = bgt.loc[bgt["sample_id"] == sid].iloc[0]
    ax.plot(xs, [r[f"bg_{s}"] for s in bg_stats], marker="o", markersize=18,
            linewidth=4, color=CONDITION_COLORS.get(COND_OF[sid], "#000000"),
            label=f"{short_label(sid)} ({CONDITION_LABELS.get(COND_OF[sid],'')})")
ax.set_xticks(xs); ax.set_xticklabels(bg_stats, rotation=30, ha="right")
ax.set_yscale("log"); ax.set_ylabel("Background density (cells / mm$^2$)")
ax.set_title("Background under each statistic", fontsize=FONT_SIZE_TITLE)
ax.legend(frameon=False, fontsize=FONT_SIZE_LEGEND - 12)
style_axes(ax); panel_letter(ax, "A")

ax = axes[1]
for sid in SAMPLE_ORDER:
    r = bgt.loc[bgt["sample_id"] == sid].iloc[0]
    ax.plot(xs, [r[f"gate6x_{s}"] for s in bg_stats], marker="s", markersize=18,
            linewidth=4, color=CONDITION_COLORS.get(COND_OF[sid], "#000000"))
ax.axhline(BURDEN_THRESHOLD, color="#000000", linestyle=":", linewidth=4)
ax.text(0, BURDEN_THRESHOLD, " burden threshold", va="bottom", ha="left",
        fontsize=FONT_SIZE_TICK - 6)
ax.set_xticks(xs); ax.set_xticklabels(bg_stats, rotation=30, ha="right")
ax.set_yscale("log")
ax.set_ylabel(f"What a {FOCUS_FOLD_OVER_BACKGROUND:g}x gate means\n(cells / mm$^2$)")
ax.set_title("The same gate, six different meanings", fontsize=FONT_SIZE_TITLE)
style_axes(ax); panel_letter(ax, "B")
fig.suptitle("A self-referencing gate is hardest on the most diseased section",
             fontsize=FONT_SIZE_TITLE, y=1.02)
fig.tight_layout()
save_fig(fig, "F74_background_statistic")

# ---- F70 candidate maps -----------------------------------------------------
ncol = 3
nrow = int(np.ceil(len(SAMPLE_ORDER) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(11 * ncol, 11 * nrow))
axes = np.atleast_1d(axes).ravel()
for ax, sid in zip(axes, SAMPLE_ORDER):
    d = cells[sid]
    n = len(d)
    idx = (rng.choice(n, MAX_POINTS_MAP, replace=False)
           if n > MAX_POINTS_MAP else np.arange(n))
    dd = d.iloc[idx]
    ax.scatter(dd["x"], dd["y"], s=0.4, color=BACKGROUND_COLOR, linewidths=0,
               rasterized=True, zorder=1)
    mm = dd["pheno"].isin(myeloid_use)
    ax.scatter(dd.loc[mm, "x"], dd.loc[mm, "y"], s=2.0, color="#8C8C8C",
               linewidths=0, rasterized=True, zorder=2)
    sl = ledger.loc[ledger["sample_id"] == sid]
    g = grids[sid]
    for acc_flag, col, mk in [(True, REASON_COLORS["accepted"], "o"),
                              (False, REASON_COLORS["fold"], "^")]:
        s = sl.loc[sl["accepted_current"] == acc_flag]
        if not len(s):
            continue
        px = g["x0"] + s["peak_col"].to_numpy() * GRID_UM
        py = g["y0"] + s["peak_row"].to_numpy() * GRID_UM
        sz = np.clip(np.sqrt(s["area_um2"].to_numpy()) * 3.0, 140, 1400)
        ax.scatter(px, py, s=sz, facecolors="none", edgecolors=col,
                   linewidths=3.5, marker=mk, zorder=4)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_edgecolor(CONDITION_COLORS.get(COND_OF[sid], "#000000"))
        sp.set_linewidth(6)
    na = int(sl["accepted_current"].sum())
    ax.set_title(f"{short_label(sid)}  {CONDITION_LABELS.get(COND_OF[sid], '')}\n"
                 f"{na} accepted (expert {EXPERT_COUNTS.get(sid, '?')}), "
                 f"{len(sl) - na} rejected   bg {bg_of[sid]:.0f}/mm$^2$",
                 fontsize=FONT_SIZE_BASE - 4)
for ax in axes[len(SAMPLE_ORDER):]:
    ax.axis("off")
handles = [Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="none",
                  markeredgecolor=REASON_COLORS["accepted"], markeredgewidth=3.5,
                  markersize=20, label="accepted at current settings"),
           Line2D([0], [0], marker="^", linestyle="none", markerfacecolor="none",
                  markeredgecolor=REASON_COLORS["fold"], markeredgewidth=3.5,
                  markersize=20, label="rejected on fold")]
axes[0].legend(handles=handles, frameon=False, loc="upper left",
               fontsize=FONT_SIZE_LEGEND - 8)
fig.suptitle("Every candidate the detector considered\nmarker area scales with "
             "the grown region", fontsize=FONT_SIZE_TITLE, y=0.995)
save_fig(fig, "F70_candidate_maps")

# ---- F71 fold vs absolute ---------------------------------------------------
fig, ax = plt.subplots(figsize=(22, 13))
for cond in CONDITION_ORDER:
    for flag, mk, lab in [(True, "o", "accepted"), (False, "^", "rejected")]:
        s = ledger.loc[(ledger["condition"] == cond)
                       & (ledger["accepted_current"] == flag)]
        if not len(s):
            continue
        ax.scatter(s["fold_over_background"], s["peak_density"], s=280,
                   marker=mk, alpha=0.75,
                   facecolors=(CONDITION_COLORS[cond] if flag else "none"),
                   edgecolors=CONDITION_COLORS[cond], linewidths=3,
                   label=f"{CONDITION_LABELS[cond]}, {lab}")
for f in FOLD_GRID:
    ax.axvline(f, color="#999999", linestyle="--", linewidth=2)
for fl in [v for v in ABS_FLOOR_GRID if v > 0]:
    ax.axhline(fl, color="#4575B4", linestyle="-.", linewidth=2)
ax.axhline(BURDEN_THRESHOLD, color="#000000", linestyle=":", linewidth=4)
ax.text(ax.get_xlim()[1], BURDEN_THRESHOLD, f"  burden {BURDEN_THRESHOLD:,.0f}",
        va="center", ha="left", fontsize=FONT_SIZE_TICK - 6)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("Fold over that section's own background")
ax.set_ylabel("Absolute peak density (cells / mm$^2$)")
ax.set_title("The admission boundary in both coordinates\n"
             "dashed vertical = fold gates, dash-dot horizontal = absolute floors",
             fontsize=FONT_SIZE_TITLE)
ax.legend(frameon=False, fontsize=FONT_SIZE_LEGEND - 8, loc="lower right")
style_axes(ax)
save_fig(fig, "F71_fold_vs_absolute_density")

# ---- F72 composition --------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(34, 12))
ax = axes[0]
for cond in CONDITION_ORDER:
    s = ledger.loc[ledger["condition"] == cond]
    ax.scatter(s["fold_over_background"], s["frac_myeloid"], s=240, alpha=0.7,
               color=CONDITION_COLORS[cond], edgecolors="#000000",
               linewidths=1.5, label=CONDITION_LABELS[cond])
ax.axvline(FOCUS_FOLD_OVER_BACKGROUND, color="#000000", linewidth=3)
ax.set_xscale("log"); ax.set_xlabel("Fold over own background")
ax.set_ylabel("Myeloid fraction of region")
ax.set_title("Myeloid fraction, no reference needed", fontsize=FONT_SIZE_TITLE - 4)
ax.legend(frameon=False, fontsize=FONT_SIZE_LEGEND - 10, loc="upper left")
style_axes(ax); panel_letter(ax, "A")

ax = axes[1]
for cond in CONDITION_ORDER:
    s = ledger.loc[ledger["condition"] == cond]
    ax.scatter(s["fold_over_background"], s["frac_endothelial"], s=240,
               alpha=0.7, color=CONDITION_COLORS[cond], edgecolors="#000000",
               linewidths=1.5)
ax.axvline(FOCUS_FOLD_OVER_BACKGROUND, color="#000000", linewidth=3)
ax.set_xscale("log"); ax.set_xlabel("Fold over own background")
ax.set_ylabel("Endothelial fraction of region")
ax.set_title("Vascular admission check", fontsize=FONT_SIZE_TITLE - 4)
style_axes(ax); panel_letter(ax, "B")

ax = axes[2]
labels_b, xs_b = [], []
bottom = None
groups = []
for lo, hi, name in FOLD_BANDS:
    m = (ledger["fold_over_background"] >= lo) & (ledger["fold_over_background"] < hi)
    if m.sum():
        groups.append((name, m))
xs_b = np.arange(len(groups))
bottom = np.zeros(len(groups))
for p in PH_USE:
    vals = np.array([float(np.nanmean(ledger.loc[m, f"frac_{p}"]))
                     if m.sum() else 0.0 for _, m in groups])
    vals = np.nan_to_num(vals)
    ax.bar(xs_b, vals, bottom=bottom, width=0.7,
           color=PHENOTYPE_COLORS.get(p, "#BBBBBB"), edgecolor="#FFFFFF",
           linewidth=1.5, label=p)
    bottom += vals
ax.set_xticks(xs_b)
ax.set_xticklabels([f"{n}\nn={int(m.sum())}" for n, m in groups],
                   fontsize=FONT_SIZE_TICK - 10)
ax.set_ylabel("Mean phenotype fraction")
ax.set_title("Composition by fold band", fontsize=FONT_SIZE_TITLE - 4)
ax.legend(frameon=False, fontsize=FONT_SIZE_LEGEND - 14, ncol=1,
          bbox_to_anchor=(1.02, 1.0), loc="upper left")
style_axes(ax); panel_letter(ax, "C")
fig.suptitle("Do the rejected regions look like the accepted foci?",
             fontsize=FONT_SIZE_TITLE, y=1.03)
fig.tight_layout()
save_fig(fig, "F72_composition")

# ---- F73 gate grid ----------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(30, 13))
ax = axes[0]
for rescue in ABS_RESCUE_GRID:
    ls = "-" if not np.isfinite(rescue) else "--"
    for floor in ABS_FLOOR_GRID:
        sel = score.loc[(score["abs_floor"] == floor)
                        & ((score["abs_rescue"] == rescue)
                           if np.isfinite(rescue)
                           else ~np.isfinite(score["abs_rescue"]))]
        if not len(sel):
            continue
        sel = sel.sort_values("fold")
        ax.plot(sel["fold"], sel["total_abs_error"], marker="o", markersize=14,
                linewidth=3, linestyle=ls,
                label=f"floor {floor:.0f}, rescue "
                      f"{'none' if not np.isfinite(rescue) else f'{rescue:.0f}'}")
ax.invert_xaxis()
ax.set_xlabel("Fold gate")
ax.set_ylabel("Total |detected - expert|\nacross the 3 expert sections")
ax.set_title("Agreement with expert counts", fontsize=FONT_SIZE_TITLE - 4)
ax.legend(frameon=False, fontsize=FONT_SIZE_LEGEND - 16, ncol=2)
style_axes(ax); panel_letter(ax, "A")

ax = axes[1]
for rescue in ABS_RESCUE_GRID:
    ls = "-" if not np.isfinite(rescue) else "--"
    for floor in ABS_FLOOR_GRID:
        sel = score.loc[(score["abs_floor"] == floor)
                        & ((score["abs_rescue"] == rescue)
                           if np.isfinite(rescue)
                           else ~np.isfinite(score["abs_rescue"]))]
        if not len(sel):
            continue
        sel = sel.sort_values("fold")
        ax.plot(sel["fold"], sel["min_peak_any"], marker="s", markersize=14,
                linewidth=3, linestyle=ls)
ax.axhline(BURDEN_THRESHOLD, color="#000000", linestyle=":", linewidth=4)
ax.text(ax.get_xlim()[0], BURDEN_THRESHOLD, " burden threshold", va="bottom",
        ha="left", fontsize=FONT_SIZE_TICK - 6)
ax.invert_xaxis(); ax.set_yscale("log")
ax.set_xlabel("Fold gate")
ax.set_ylabel("Lowest absolute peak density\nadmitted anywhere (cells / mm$^2$)")
ax.set_title("What each setting lets in", fontsize=FONT_SIZE_TITLE - 4)
style_axes(ax); panel_letter(ax, "B")
fig.suptitle("The hybrid gate: agreement against admission floor",
             fontsize=FONT_SIZE_TITLE, y=1.02)
fig.tight_layout()
save_fig(fig, "F73_gate_grid")


# %% Cell 11 - wrap up
# =============================================================================

banner("SUMMARY")

print(f"Cross-check verdict   : {verdict.upper()}")
print(f"Sections audited      : {len(SAMPLE_ORDER)}")
print(f"Candidates logged     : {len(ledger)}")
print(f"Accepted at current   : {int(ledger['accepted_current'].sum())}")
print(f"Gate combinations     : {n_combo}")

sub("Read in this order")
print("  1. The cross-check verdict. Backgrounds first: they test the map")
print("     alone. Nothing below is usable until that passes.")
print("  2. F74 and table 91. What a 6x gate means in absolute terms in each")
print("     section. This is the structural problem, and no fold value fixes")
print("     it for the section whose own background is already pathological.")
print("  3. F72 panels A and B and the fold-band table. Myeloid fraction needs")
print("     no reference, endothelial fraction is the vascular check.")
print("  4. Table 96 and F73. Which hybrid setting agrees with the three")
print("     expert counts WITHOUT admitting anything below burden density.")
print("  5. Table 95. What a composition gate would cost. It is reported, not")
print("     applied, because a myeloid floor tuned to cut the rejected band")
print("     also cuts most of Finding 1's residual treated foci.")

sub("The decision this run should let us make")
print("  Whether the fix is a lower fold, an absolute floor underneath it, a")
print("  rescue rule above it, a different background statistic, or some")
print("  combination. The control row is script 04 as it stands, so every")
print("  candidate setting is read as a change against a known baseline.")

sub("NOT changed by this script")
print("  The burden definition, script 04, the cell assignment files and every")
print("  downstream script. Any adopted change becomes script 04 revision 5,")
print("  after which 06, 07, 08 and 09 all rerun.")

banner("DONE")
sys.stdout = _tee.terminal
_tee.close()
