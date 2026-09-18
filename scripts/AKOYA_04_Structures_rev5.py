#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - STRUCTURE DEFINITION (triple definition)
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 04 of the AKOYA analysis series. REVISION 5.

REVISION 5 IS REVISION 4 WITH FIVE CHANGES. Everything else, including the whole
output schema, the marker columns, the cuff geometry, BALT detection and the
figures, is carried over unchanged, so scripts 06, 07 and 08 need only their
IN_DIR repointed from structures_rev4 to structures_rev5.

WHAT CHANGED, AND WHAT MEASUREMENT DROVE EACH CHANGE

  Parameters were calibrated against 83 expert annotations across five sections,
  rasterised onto this script's own 25 um grid (scripts 04c, 04d, 04e, 04f).
  They were not tuned to any between-arm contrast and apply to both arms.

  1. BACKGROUND_STAT  "median" -> "p25"
     Section backgrounds under the median ran 211 to 1,686 cells/mm2, so a 6x
     gate meant 1,269 cells/mm2 in one section and 10,117 in another. Annotated
     lesions on G4_36463 reached only 4.9x its median background, BELOW the
     gate, while every other section reached 7.7 to 19.8x. The median of a
     heavily involved section is itself pathology, so a self-referencing gate is
     hardest to clear in the most diseased tissue. Under p25 all five annotated
     sections clear 6x, G4_36463 at 12.5x.

  2. FOCUS_HALF_MAX_FRACTION  0.50 -> 0.25
     Annotated boundaries sit at 0.224 of the enclosed peak for complexes and
     0.258 for sharp-bordered lesions. Corroborated independently: revision 4
     structures had a median area 0.31x the regions drawn over them.

  3. MIN_PEAK_DENSITY = 1750 cells/mm2, NEW
     Window fixed at 1,565 to 1,925 by the annotations themselves: above a
     structure the pathologist called vascular, below one overlapping a drawn
     lesion. Of 18 objects it excludes NONE contain IDO1+ macrophages; of 57 it
     retains, 44 do (Fisher exact p = 2.2e-09). IDO1 is used nowhere in
     detection, so that partition is a consequence of the threshold, not an
     input to it.

  4. GRANULOMA COMPLEXES, a third definition, NEW
     G4_36463 holds 13 local maxima in its entire area, and 15 of the 24 regions
     annotated on it contain no local maximum at all. No threshold on a peak
     list recovers a region that never generated a peak. Complexes are connected
     components above a section-relative level with no peak requirement, so a
     confluent mass is one structure because it is one connected component.

  5. THE PARAMETER SWEEP IS RETIRED
     Revision 4 scored settings against six expert counts. Three free parameters
     against three informative counts is fitting rather than selection, the
     metric cannot separate over-detection from under-detection, and every
     top-ranked setting scored by reducing each treated section to one
     structure. Replaced by region-level agreement in 04d and 04e. Tables 34 and
     34b and figure F32 are therefore not produced.

  PEAK_SEPARATION_UM is UNCHANGED at 200. Measurement suggested 772, which would
  set the maximum-filter footprint to 1,575 um and merge genuinely distinct foci
  in the three sections where detection was already acceptable. That measurement
  correctly diagnoses confluence; complexes address it without the cost.

THREE DEFINITIONS, NEVER POOLED
  burden_region_id  ABSOLUTE, one fixed threshold applied identically to every
                    section. UNCHANGED. The yardstick that produces the zero in
                    the treated arm. It must not become section-relative.
  focus_id          PEAK-BASED, for internal architecture. Supplies the centre
                    and inscribed radius the radial coordinate needs. Scripts
                    06, 07 and 08 read this.
  complex_id        LEVEL-SET, for counting and extent. NEW.

CIRCULARITY, AND THE SWITCH THAT ADDRESSES IT
  Structures are defined by POOLED MYELOID DENSITY, which constrains exactly one
  quantity. It does not constrain composition within the pool, where any cell
  type sits, or any marker intensity.
    - Lymphocyte results are entirely free: those phenotypes are not in the pool.
    - Radial position of a pooled member is nearly free. The only coupling is
      that the boundary sits where density falls to a fraction of peak, so a
      rim-ward population extends the boundary it is then measured against,
      compressing its own normalised position. That is CONSERVATIVE for a
      rim-ward finding.
    - ABUNDANCE of a pooled member inside structures IS exposed.
  Set LEAVE_ONE_OUT to a phenotype to re-run with it excluded from the pool,
  writing to a separate directory, so a population can be tested inside
  structures defined without it. Run the primary first. Never mix the outputs.

CRITICAL DESIGN CONSTRAINT (unchanged from revision 4)
  No definition uses IDO1. All run on pooled myeloid density. IDO1 is measured
  about each structure afterwards, and F35 validates capture independently.

OUTPUTS
    figures/  F30, F31, F33, F34, F35, F36, F37, F38 (new)
    tables/   32, 33, 33b, 35, 35b, 36, 36b, 37 (new), 38 (schema check)
    cell_assignments/<section>_cell_structures.csv

USAGE
    conda activate sc_pre
    python AKOYA_04_Structures.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================

DATA_DIR = "/master/jlehle/WORKING/AKOYA/data"
INVENTORY_TABLE_DIR = "/master/jlehle/WORKING/AKOYA/inventory/tables"
OUT_DIR = "/master/jlehle/WORKING/AKOYA/structures_rev5"

EXCLUDE_SECTIONS = {
    "G3_43102": ("Position-1 section on scan_01 (D1MT). Dimmest of its four "
                 "sections in 48.8% of composition-controlled comparisons "
                 "against a 25% null. No other section on that scan exceeds "
                 "27.8%. Zero IDO1+ macrophages called."),
    "G4_43112": ("Position-1 section on scan_02 (Untreated). Dimmest of its "
                 "four sections in 56.3% of composition-controlled comparisons "
                 "against a 25% null. No other section on that scan exceeds "
                 "23.3%. Zero IDO1+ macrophages called."),
}

CONDITION_ORDER = ["D1MT", "Untreated"]
CONDITION_COLORS = {"D1MT": "#2C7FB8", "Untreated": "#D95F02"}
GROUP_MAP = {"G3": "D1MT", "G4": "Untreated"}
SCAN_COL = "scan_id"
POSITION_COL = "slide_position_rank"

# Reporting only. The sweep that used these as an objective is retired. Only the
# three G3 counts are expert; the G4 targets came from the revision-1 algorithm
# and agreement with them is not evidence.
EXPERT_COUNTS = {"G3_43106": 1, "G3_43111": 2, "G3_43118": 2,
                 "G4_31438": 7, "G4_36463": 5, "G4_43109": 19}
EXPERT_SOURCE = {"G3_43106": "expert", "G3_43111": "expert",
                 "G3_43118": "expert", "G4_31438": "algorithm_rev1",
                 "G4_36463": "algorithm_rev1", "G4_43109": "algorithm_rev1"}

# ---- columns ----------------------------------------------------------------
CENTROID_X_PREFIX = "Centroid X"
CENTROID_Y_PREFIX = "Centroid Y"
PHENOTYPE_COL = "Phenotypes"
IMAGE_COL = "Image"

IDO1_COL = "IDO1: Cytoplasm: Mean"
HK3_COL = "3-Hydroxykynurenine: Cytoplasm: Mean"
CD21_COL = "CD21: Membrane: Mean"
CD68_COL = "CD68: Membrane: Mean"
INOS_COL = "iNOS: Membrane: Mean"
ARG1_COL = "Arginase-1: Cytoplasm: Mean"
IFNG_COL = "IFNG: Membrane: Mean"
EXTRA_MARKER_COLS = [IDO1_COL, HK3_COL, CD21_COL, CD68_COL,
                     INOS_COL, ARG1_COL, IFNG_COL]

# ---- phenotypes -------------------------------------------------------------
IDO1_POS_PHENO = "CD68+IDO1+ Macrophages"
IDO1_NEG_PHENO = "CD68+IDO1- Macrophages"
CD68_LINEAGE = [IDO1_POS_PHENO, IDO1_NEG_PHENO]
MYELOID_FOR_DETECTION = [IDO1_POS_PHENO, IDO1_NEG_PHENO,
                         "CD163+ Macrophages", "Neutrophils"]

# Set to a phenotype in MYELOID_FOR_DETECTION to exclude it from the pool and
# write to a separate directory. None runs the primary analysis.
LEAVE_ONE_OUT = None

BALT_B, BALT_T, PLASMA = "B cells", "Helper T cells", "Plasma cells"
T_LINEAGE = ["Helper T cells", "CD4- T cells", "Tregs"]

PHENOTYPE_ORDER = [
    IDO1_POS_PHENO, IDO1_NEG_PHENO, "CD163+ Macrophages", "Neutrophils",
    "Helper T cells", "CD4- T cells", "Tregs", "B cells", "Plasma cells",
    "Endothelial cells", "Epithelial/Tumor cells", "Other",
]
PHENOTYPE_COLORS = {
    IDO1_POS_PHENO: "#B2182B", IDO1_NEG_PHENO: "#EF8A62",
    "CD163+ Macrophages": "#FDBE85", "Neutrophils": "#7B3294",
    "Helper T cells": "#1B7837", "CD4- T cells": "#7FBC41", "Tregs": "#00441B",
    "B cells": "#2166AC", "Plasma cells": "#67A9CF",
    "Endothelial cells": "#E7298A", "Epithelial/Tumor cells": "#66C2A5",
    "Other": "#A6761D",
}

# ---- spatial grid -----------------------------------------------------------
GRID_UM = 25.0
DENSITY_BANDWIDTH_UM = 75.0
EDGE_CORRECTION = True
EDGE_MIN_WEIGHT = 0.25
REPORT_EDGE_COMPARISON = True

# =============================================================================
# DEFINITION 1 - ABSOLUTE, for burden. UNCHANGED IN REVISION 5.
# =============================================================================
BURDEN_THRESHOLD = 3000.0
BURDEN_MIN_AREA_UM2 = 30000.0
BURDEN_MIN_CELLS = 100
BURDEN_FILL_HOLES = True

# =============================================================================
# DEFINITION 2 - RELATIVE FOCI, peak-based, for architecture
# =============================================================================
PEAK_SEPARATION_UM = 200.0
FOCUS_FOLD_OVER_BACKGROUND = 6.0
FOCUS_HALF_MAX_FRACTION = 0.25      # revision 4: 0.50
FOCUS_MIN_AREA_UM2 = 10000.0
FOCUS_MIN_CELLS = 40
FOCUS_FILL_HOLES = True
CUFF_WIDTH_UM = 150.0
DEDUPLICATE_PLATEAU_PEAKS = True

BACKGROUND_STAT = "p25"             # revision 4: "median"
MIN_PEAK_DENSITY = 1750.0           # NEW. Applies to foci AND complexes.

# =============================================================================
# DEFINITION 3 - GRANULOMA COMPLEXES, level-set. NEW IN REVISION 5.
# =============================================================================
COMPLEX_LEVEL_FOLD = 6.5
COMPLEX_CLOSE_UM = 0.0
COMPLEX_MIN_AREA_UM2 = 10000.0
COMPLEX_MIN_CELLS = 40

# ---- BALT detection, UNCHANGED FROM REVISION 4 ------------------------------
BALT_BANDWIDTH_UM = 75.0
BALT_B_THRESHOLD = 800.0
BALT_T_THRESHOLD = 300.0
BALT_MIN_AREA_UM2 = 15000.0
BALT_MIN_B_CELLS = 40
CD21_FOLLICLE_RATIO = 1.5
RUN_BALT_ALTERNATIVES = True
BALT_ALT_GATES = [
    ("primary_B_and_helperT", [BALT_B], [BALT_T], BALT_B_THRESHOLD, BALT_T_THRESHOLD),
    ("B_and_pooled_T", [BALT_B], T_LINEAGE, BALT_B_THRESHOLD, BALT_T_THRESHOLD),
    ("B_only", [BALT_B], None, BALT_B_THRESHOLD, None),
]

# ---- schema guard -----------------------------------------------------------
# Every column revision 4 wrote that a downstream script reads. Asserted at the
# end, because a silent rename is how a downstream script gets a column of
# zeros instead of an error.
REQUIRED_CELL_COLS = [
    "x", "y", "pheno", "condition", "sample_id", "scan_id",
    "slide_position_rank", "focus_id", "region", "radial_pos",
    "dist_to_focus_um", "burden_region_id", "in_burden_region", "balt_id",
    "myeloid_density", "section_background_density",
    "IDO1", "3-Hydroxykynurenine", "CD21", "CD68", "iNOS", "Arginase-1", "IFNG",
]
REQUIRED_FOCI_COLS = [
    "focus_id", "sample_id", "condition", "scan_id", "peak_density",
    "background_density", "fold_over_background", "boundary_density",
    "area_um2", "equiv_radius_um", "max_inscribed_radius_um", "shape_ratio",
    "n_cells", "n_cells_core", "n_cells_cuff", "centroid_x_um", "centroid_y_um",
    "pct_of_tissue_area", "geometry_matches_assignment",
    "n_macrophages_core", "n_ido1_pos_core", "pct_ido1_pos_of_mac_core",
    "median_ido1_macrophages", "median_hk3_macrophages",
    "median_inos_macrophages", "median_arg1_macrophages", "median_ifng_core",
    "overlaps_burden_region", "pct_core_in_burden_region",
]

# ---- plotting ---------------------------------------------------------------
MAX_POINTS_PER_SECTION = 80000
POINT_MARKERS = ["o", "s", "^", "D", "v", "P"]
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
BACKGROUND_COLOR = "#ECECEC"
BURDEN_OUTLINE = "#7B3294"
FOCUS_OUTLINE = "#000000"
CUFF_OUTLINE = "#666666"
BALT_OUTLINE = "#2166AC"
COMPLEX_OUTLINE = "#1A9850"
FLAG_COLOR = "#B2182B"

USE_FLOAT32 = True
WRITE_CELL_ASSIGNMENTS = True


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

try:
    from scipy import ndimage as ndi
    HAVE_SCIPY = True
except Exception as _e:
    HAVE_SCIPY = False
    _scipy_err = _e

try:
    from skimage.segmentation import watershed as _sk_watershed
    HAVE_SKIMAGE = True
except Exception as _e2:
    HAVE_SKIMAGE = False
    _skimage_err = _e2

if LEAVE_ONE_OUT:
    if LEAVE_ONE_OUT not in MYELOID_FOR_DETECTION:
        raise SystemExit(f"LEAVE_ONE_OUT='{LEAVE_ONE_OUT}' is not in "
                         f"MYELOID_FOR_DETECTION: {MYELOID_FOR_DETECTION}")
    POOL = [p for p in MYELOID_FOR_DETECTION if p != LEAVE_ONE_OUT]
    OUT_DIR = f"{OUT_DIR}_no_{LEAVE_ONE_OUT.replace(' ', '_').replace('+', '')}"
else:
    POOL = list(MYELOID_FOR_DETECTION)

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
CELL_DIR = os.path.join(OUT_DIR, "cell_assignments")
for d in (FIG_DIR, TAB_DIR, CELL_DIR):
    os.makedirs(d, exist_ok=True)

rng = np.random.default_rng(RANDOM_SEED)


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
    """Revision 4's exact normalisation. Script 01 found a zero-width space in
    the epithelial label, which breaks equality matching with no visible sign."""
    return "".join(ch for ch in str(s) if ord(ch) <= 127).strip()


def short_label(sid):
    return sid.split("_")[-1] if "_" in sid else sid


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


def write_csv(df, name, directory=None):
    p = os.path.join(directory or TAB_DIR, name)
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


def make_grid(x, y, pitch):
    x0 = float(np.floor(x.min() / pitch) * pitch)
    y0 = float(np.floor(y.min() / pitch) * pitch)
    gx = ((x - x0) / pitch).astype(int) + 1
    gy = ((y - y0) / pitch).astype(int) + 1
    return x0, y0, gx, gy, (int(gy.max()) + 2, int(gx.max()) + 2)


def density_map(gx, gy, shape, mask, bandwidth_um, pitch, occ=None,
                edge_correct=True, min_weight=EDGE_MIN_WEIGHT):
    """Normalized convolution, so density is not damped within one bandwidth of
    a tissue boundary. Unchanged from revision 4."""
    counts = np.zeros(shape, dtype=float)
    np.add.at(counts, (gy[mask], gx[mask]), 1.0)
    sigma = bandwidth_um / pitch
    sm = ndi.gaussian_filter(counts, sigma=sigma, mode="constant")
    if edge_correct and occ is not None:
        w = ndi.gaussian_filter(occ.astype(float), sigma=sigma, mode="constant")
        sm = sm / np.maximum(w, min_weight)
    return sm / ((pitch * pitch) / 1e6)


def tissue_mask(gx, gy, shape):
    occ = np.zeros(shape, dtype=bool)
    occ[gy, gx] = True
    return ndi.binary_fill_holes(ndi.binary_closing(occ, iterations=2))


def background_of(dens, occ, stat=BACKGROUND_STAT, tag=""):
    """
    Section background. REVISION 5 uses p25.

    A low percentile can legitimately be zero if more than a quarter of tissue
    pixels carry no myeloid density, which makes every fold infinite. The median
    cannot do this in practice but p25 can, and p25 is now the default, so the
    statistic escalates loudly rather than failing silently.
    """
    v = dens[occ]
    v = v[np.isfinite(v)]
    if not v.size:
        return np.nan

    def _s(st):
        if st in ("median", "p50"):
            return float(np.median(v))
        if st.startswith("p") and st[1:].replace(".", "").isdigit():
            return float(np.percentile(v, float(st[1:])))
        return float(np.median(v))

    bg = _s(stat)
    if np.isfinite(bg) and bg > 0:
        return bg
    for alt in ["p50", "p75", "p90"]:
        b = _s(alt)
        if np.isfinite(b) and b > 0:
            print(f"    WARNING [{tag}]: background '{stat}' is {bg}. Escalated "
                  f"to '{alt}' = {b:.2f} /mm2. Most tissue pixels here carry no")
            print("    myeloid density, so treat every fold-based number for "
                  "this section as unreliable.")
            return b
    pos = v[v > 0]
    if pos.size:
        b = float(pos.mean())
        print(f"    WARNING [{tag}]: no positive percentile. Using the mean of "
              f"positive pixels, {b:.2f} /mm2.")
        return b
    return np.nan


def partition_by_markers(dens, markers, territory):
    if HAVE_SKIMAGE:
        return _sk_watershed(-dens, markers=markers, mask=territory), "watershed"
    _, inds = ndi.distance_transform_edt(markers == 0, return_indices=True)
    part = markers[inds[0], inds[1]]
    return np.where(territory, part, 0), "edt_voronoi"


def find_local_maxima(dens, occ, sep_px, deduplicate=True):
    """dens == maximum_filter(dens) marks EVERY pixel of a flat maximum, not one
    per peak. Connected plateau components are collapsed to one representative."""
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
        yy, xx = np.nonzero(lab[sl] == k)
        yy = yy + sl[0].start
        xx = xx + sl[1].start
        j = int(np.argmin((yy - com[0]) ** 2 + (xx - com[1]) ** 2))
        ys.append(int(yy[j])); xs.append(int(xx[j]))
    py = np.asarray(ys, dtype=int)
    px = np.asarray(xs, dtype=int)
    return py, px, dens[py, px], n_raw, n_raw - len(py)


def detect_foci(dens, occ, gx, gy, bg, fold_min, half_frac, min_area, min_cells,
                sep_px, min_peak=MIN_PEAK_DENSITY, fill_holes=True,
                deduplicate=True):
    """
    DEFINITION 2. Revision 4's function with two changes: the background is
    passed in rather than recomputed here, and MIN_PEAK_DENSITY is applied
    alongside the fold gate.

    Geometry is computed from the pixels a focus ACTUALLY receives, so table 35
    cannot disagree with the label array or the per-cell assignments.
    """
    diag = {"n_peaks_raw": 0, "n_peaks_collapsed": 0, "n_peaks_kept": 0,
            "n_contested_pixels": 0, "n_dropped_min_peak": 0,
            "n_dropped_area_cells": 0}
    if not occ.any() or not np.isfinite(bg) or bg <= 0:
        return np.zeros_like(dens, dtype=int), [], "none", diag

    py, px, pv, n_raw, n_coll = find_local_maxima(dens, occ, sep_px, deduplicate)
    diag["n_peaks_raw"] = n_raw
    diag["n_peaks_collapsed"] = n_coll
    pass_fold = pv >= fold_min * bg
    keep = pass_fold & (pv >= min_peak)
    diag["n_dropped_min_peak"] = int((pass_fold & (pv < min_peak)).sum())
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
        contested = int(region.sum() - assigned.sum())
        diag["n_contested_pixels"] += contested
        if not assigned[yy, xx]:
            continue
        cl2, cn2 = ndi.label(assigned)
        if cn2 == 0:
            continue
        assigned = cl2 == cl2[yy, xx]

        area = assigned.sum() * GRID_UM * GRID_UM
        ncell = int(assigned[gy, gx].sum())
        if area < min_area or ncell < min_cells:
            diag["n_dropped_area_cells"] += 1
            continue

        dt = ndi.distance_transform_edt(assigned)
        r_in = float(dt.max()) * GRID_UM
        r_eq = float(np.sqrt(area / np.pi))
        labels[assigned] = next_id
        recs.append({
            "focus_id": next_id,
            "peak_density": float(vv),
            "background_density": float(bg),
            "fold_over_background": float(vv / bg),
            "boundary_density": float(half_frac * vv),
            "peak_row": int(yy), "peak_col": int(xx),
            "area_um2": float(area),
            "equiv_radius_um": r_eq,
            "max_inscribed_radius_um": r_in,
            "shape_ratio": float(r_eq / r_in) if r_in > 0 else np.nan,
            "n_contested_pixels": contested,
            "n_cells": ncell,
        })
        next_id += 1
    return labels, recs, method, diag


def detect_complexes(dens, occ, gx, gy, bg):
    """
    DEFINITION 3, NEW. Connected components above a section-relative level with
    no peak requirement, so a confluent mass is one structure because it is one
    connected component. G4_36463 needs this: 15 of its 24 annotated regions
    contain no local maximum at all.
    """
    diag = {"n_components": 0, "n_dropped_area_cells": 0,
            "n_dropped_min_peak": 0}
    if not np.isfinite(bg) or bg <= 0:
        return np.zeros_like(dens, dtype=int), [], diag
    m = (dens >= COMPLEX_LEVEL_FOLD * bg) & occ
    if COMPLEX_CLOSE_UM > 0:
        r = max(1, int(round(COMPLEX_CLOSE_UM / GRID_UM)))
        m = ndi.binary_closing(m, structure=ndi.generate_binary_structure(2, 1),
                               iterations=r) & occ
    if FOCUS_FILL_HOLES:
        m = ndi.binary_fill_holes(m)
    lab, n = ndi.label(m)
    diag["n_components"] = n
    labels = np.zeros(dens.shape, dtype=int)
    recs = []
    next_id = 1
    for i in range(1, n + 1):
        comp = lab == i
        area = comp.sum() * GRID_UM * GRID_UM
        ncell = int(comp[gy, gx].sum())
        if area < COMPLEX_MIN_AREA_UM2 or ncell < COMPLEX_MIN_CELLS:
            diag["n_dropped_area_cells"] += 1
            continue
        peak = float(dens[comp].max())
        if peak < MIN_PEAK_DENSITY:
            diag["n_dropped_min_peak"] += 1
            continue
        dt = ndi.distance_transform_edt(comp)
        r_eq = float(np.sqrt(area / np.pi))
        r_in = float(dt.max()) * GRID_UM
        yy, xx = np.unravel_index(int(np.argmax(np.where(comp, dens, -np.inf))),
                                  dens.shape)
        labels[comp] = next_id
        recs.append({
            "complex_id": next_id, "peak_density": peak,
            "background_density": float(bg),
            "fold_over_background": float(peak / bg),
            "level_density": float(COMPLEX_LEVEL_FOLD * bg),
            "peak_row": int(yy), "peak_col": int(xx),
            "area_um2": float(area), "equiv_radius_um": r_eq,
            "max_inscribed_radius_um": r_in,
            "shape_ratio": float(r_eq / r_in) if r_in > 0 else np.nan,
            "n_cells": ncell,
        })
        next_id += 1
    return labels, recs, diag


def detect_burden(dens, occ, gx, gy, threshold, min_area, min_cells,
                  fill_holes=True):
    """DEFINITION 1. Fixed absolute threshold. UNCHANGED FROM REVISION 4."""
    mask = (dens >= threshold) & occ
    if fill_holes:
        mask = ndi.binary_fill_holes(mask)
    lab, n = ndi.label(mask)
    labels = np.zeros_like(lab)
    recs = []
    next_id = 1
    for i in range(1, n + 1):
        m = lab == i
        area = m.sum() * GRID_UM * GRID_UM
        ncell = int(m[gy, gx].sum())
        if area < min_area or ncell < min_cells:
            continue
        labels[m] = next_id
        recs.append({
            "region_id": next_id, "area_um2": float(area),
            "equiv_radius_um": float(np.sqrt(area / np.pi)),
            "n_cells": ncell, "peak_density": float(dens[m].max()),
            "mean_density": float(dens[m].mean()),
        })
        next_id += 1
    return labels, recs


def detect_lymphoid(gx, gy, shape, occ, pheno, b_set, t_set, b_thr, t_thr,
                    bandwidth, min_area, min_b_cells):
    """UNCHANGED FROM REVISION 4. t_set None means the B gate alone."""
    bd = density_map(gx, gy, shape, np.isin(pheno, b_set), bandwidth, GRID_UM,
                     occ=occ, edge_correct=EDGE_CORRECTION)
    gate = (bd >= b_thr) & occ
    if t_set is not None:
        td = density_map(gx, gy, shape, np.isin(pheno, t_set), bandwidth,
                         GRID_UM, occ=occ, edge_correct=EDGE_CORRECTION)
        gate = gate & (td >= t_thr)
    gate = ndi.binary_fill_holes(gate)
    lab, n = ndi.label(gate)
    is_b = np.isin(pheno, b_set)
    keep = [i for i in range(1, n + 1)
            if ((lab == i).sum() * GRID_UM * GRID_UM >= min_area
                and int((((lab == i)[gy, gx]) & is_b).sum()) >= min_b_cells)]
    relab = np.zeros_like(lab)
    for new_i, old_i in enumerate(keep, start=1):
        relab[lab == old_i] = new_i
    return relab, len(keep)


_tee = Tee(os.path.join(TAB_DIR, "00_structures_report.txt"))
sys.stdout = _tee

banner("AKOYA STRUCTURE DEFINITION (revision 5, triple definition)")
print(f"Run time : {datetime.now().isoformat(timespec='seconds')}")
print(f"Output   : {OUT_DIR}")
print(f"Detection pool (IDO1-blind): {POOL}")
print(f"Grid {GRID_UM} um, bandwidth {DENSITY_BANDWIDTH_UM} um")
print(f"Edge correction: {EDGE_CORRECTION} (min weight {EDGE_MIN_WEIGHT})")
if LEAVE_ONE_OUT:
    print(f"\n    LEAVE-ONE-OUT RUN: '{LEAVE_ONE_OUT}' excluded from the pool.")
    print("    SENSITIVITY analysis. Output is in its own directory and must")
    print("    not be mixed with the primary run.")
if not HAVE_SCIPY:
    print(f"\nERROR: scipy required ({_scipy_err}).")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)
print("    watershed  : " + ("skimage available, using true watershed"
                             if HAVE_SKIMAGE
                             else f"NOT available ({_skimage_err}), EDT fallback"))

sub("WHAT CHANGED FROM REVISION 4")
print(f"    background statistic   median -> {BACKGROUND_STAT}")
print(f"    growth fraction        0.50 -> {FOCUS_HALF_MAX_FRACTION:.2f}")
print(f"    minimum peak density   none -> {MIN_PEAK_DENSITY:,.0f} /mm2")
print(f"    granuloma complexes    NEW, level set at {COMPLEX_LEVEL_FOLD:g}x")
print("    parameter sweep        RETIRED, see the docstring")
print(f"    peak separation        {PEAK_SEPARATION_UM:.0f} um, UNCHANGED")
print(f"    burden definition      {BURDEN_THRESHOLD:,.0f} /mm2, UNCHANGED")

sub("EXCLUDED SECTIONS")
for sid, reason in EXCLUDE_SECTIONS.items():
    print(f"    {sid}: {reason}")


# %% Cell 3 - load
# =============================================================================

banner("LOADING RETAINED SECTIONS")

SCAN_OF, POS_OF = {}, {}
fs_path = os.path.join(INVENTORY_TABLE_DIR, "01_file_summary.csv")
if os.path.exists(fs_path):
    _inv = pd.read_csv(fs_path)
    if SCAN_COL in _inv.columns:
        SCAN_OF = dict(zip(_inv["sample_id"], _inv[SCAN_COL].astype(str)))
    if POSITION_COL in _inv.columns:
        POS_OF = {r["sample_id"]: (int(r[POSITION_COL])
                                   if pd.notna(r[POSITION_COL]) else None)
                  for _, r in _inv.iterrows()}
    print(f"    inventory loaded, scan and position for {len(SCAN_OF)} section(s)")
else:
    print(f"    WARNING: {fs_path} not found, scan and position unavailable")

csv_paths = sorted(glob.glob(os.path.join(DATA_DIR, "*.csv")))
retained = [p for p in csv_paths
            if os.path.splitext(os.path.basename(p))[0] not in EXCLUDE_SECTIONS]
print(f"    {len(csv_paths)} found, {len(EXCLUDE_SECTIONS)} excluded, "
      f"{len(retained)} retained")
if not retained:
    print("ERROR: nothing to process.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

hdr = pd.read_csv(retained[0], nrows=0)
ALL_COLS = hdr.columns.tolist()
XCOL = find_col(ALL_COLS, CENTROID_X_PREFIX)
YCOL = find_col(ALL_COLS, CENTROID_Y_PREFIX)
if XCOL is None or YCOL is None or PHENOTYPE_COL not in ALL_COLS:
    print("ERROR: required columns not found.")
    print(f"    x ('{CENTROID_X_PREFIX}...')  : {XCOL}")
    print(f"    y ('{CENTROID_Y_PREFIX}...')  : {YCOL}")
    print(f"    phenotype ('{PHENOTYPE_COL}') : {PHENOTYPE_COL in ALL_COLS}")
    print(f"    columns present: {', '.join(map(str, ALL_COLS[:15]))}")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

present_extra = [c for c in EXTRA_MARKER_COLS if c in ALL_COLS]
missing_extra = [c for c in EXTRA_MARKER_COLS if c not in ALL_COLS]
if missing_extra:
    print(f"    WARNING: markers not found, will be NaN: {missing_extra}")

use_cols = [XCOL, YCOL, PHENOTYPE_COL] + present_extra
if IMAGE_COL in ALL_COLS:
    use_cols.append(IMAGE_COL)
dtype_map = ({c: "float32" for c in present_extra + [XCOL, YCOL]}
             if USE_FLOAT32 else None)

cells, geom, dens_of = {}, {}, {}
edge_rows, skipped = [], []
for path in retained:
    sid = os.path.splitext(os.path.basename(path))[0]
    cond = GROUP_MAP.get(sid.split("_")[0], "UNKNOWN")
    try:
        df = pd.read_csv(path, usecols=use_cols, dtype=dtype_map, low_memory=False)
    except Exception as e:
        print(f"    ERROR reading {sid}: {e}. Skipping.")
        skipped.append((sid, "unreadable"))
        continue
    d = pd.DataFrame({
        "x": pd.to_numeric(df[XCOL], errors="coerce"),
        "y": pd.to_numeric(df[YCOL], errors="coerce"),
        "pheno": df[PHENOTYPE_COL].map(ascii_safe),
    })
    for c in EXTRA_MARKER_COLS:
        d[c.split(":")[0].strip()] = (pd.to_numeric(df[c], errors="coerce")
                                      if c in df.columns else np.nan)
    n_bad = int((~np.isfinite(d["x"]) | ~np.isfinite(d["y"])).sum())
    if n_bad:
        print(f"    WARNING: {sid} dropped {n_bad:,} bad-coordinate cells")
    d = d.loc[np.isfinite(d["x"]) & np.isfinite(d["y"])].reset_index(drop=True)
    d["condition"] = cond

    x, y = d["x"].to_numpy(), d["y"].to_numpy()
    x0, y0, gx, gy, shape = make_grid(x, y, GRID_UM)
    occ = tissue_mask(gx, gy, shape)
    mye = d["pheno"].isin(POOL).to_numpy()
    if not mye.any():
        print(f"    ERROR: {sid} has NO cells from the detection pool.")
        print(f"      pool: {POOL}")
        print(f"      phenotypes present: {sorted(set(d['pheno']))[:8]}")
        print("      SECTION SKIPPED, no cell assignment file written.")
        skipped.append((sid, "no cells from the detection pool"))
        continue

    dens = density_map(gx, gy, shape, mye, DENSITY_BANDWIDTH_UM, GRID_UM,
                       occ=occ, edge_correct=EDGE_CORRECTION)
    bg = background_of(dens, occ, BACKGROUND_STAT, sid)
    if not np.isfinite(bg) or bg <= 0:
        # Writing a file here would hand downstream scripts a section with zero
        # structures, which is indistinguishable from a real absence of disease.
        print(f"    ERROR: {sid} has no usable background ({bg}). SECTION "
              f"SKIPPED, no file written.")
        skipped.append((sid, f"no usable background ({bg})"))
        continue

    cells[sid] = d
    geom[sid] = {"x0": x0, "y0": y0, "gx": gx, "gy": gy, "shape": shape,
                 "occ": occ, "scan": SCAN_OF.get(sid, "na"),
                 "position": POS_OF.get(sid),
                 "tissue_area_um2": float(occ.sum() * GRID_UM * GRID_UM),
                 "bg": bg, "bg_median": float(np.median(dens[occ])),
                 "peak": float(dens[occ].max())}
    dens_of[sid] = dens

    if REPORT_EDGE_COMPARISON:
        dens_raw = density_map(gx, gy, shape, mye, DENSITY_BANDWIDTH_UM,
                               GRID_UM, occ=occ, edge_correct=False)
        w = ndi.gaussian_filter(occ.astype(float),
                                sigma=DENSITY_BANDWIDTH_UM / GRID_UM,
                                mode="constant")
        edge_rows.append({
            "sample_id": sid, "condition": cond,
            "tissue_area_mm2": geom[sid]["tissue_area_um2"] / 1e6,
            "pct_tissue_edge_affected": 100.0 * float((w[occ] < 0.95).mean()),
            "median_edge_weight": float(np.median(w[occ])),
            "background_stat": BACKGROUND_STAT,
            "bg_corrected": bg,
            "bg_uncorrected": background_of(dens_raw, occ, BACKGROUND_STAT,
                                            f"{sid}/uncorrected"),
            "bg_median_corrected": geom[sid]["bg_median"],
            "peak_corrected": geom[sid]["peak"],
            "peak_uncorrected": float(dens_raw[occ].max()),
            "pct_above_burden_corrected":
                100.0 * float((dens[occ] >= BURDEN_THRESHOLD).mean()),
            "pct_above_burden_uncorrected":
                100.0 * float((dens_raw[occ] >= BURDEN_THRESHOLD).mean()),
        })
        del dens_raw, w

    print(f"    {sid:<12} {cond:<10} {len(d):>9,} cells   "
          f"tissue={geom[sid]['tissue_area_um2']/1e6:5.1f} mm^2   "
          f"bg({BACKGROUND_STAT})={bg:>6.0f}   median={geom[sid]['bg_median']:>6.0f}   "
          f"peak={geom[sid]['peak']:>7.0f} ({geom[sid]['peak']/bg:.1f}x)")
    del df
    gc.collect()

if not cells:
    print("\nERROR: no sections processed.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

cond_rank = {c: i for i, c in enumerate(CONDITION_ORDER)}
SAMPLE_ORDER = sorted(cells, key=lambda s: (cond_rank.get(cells[s]["condition"].iloc[0], 9), s))
COND_OF = {s: cells[s]["condition"].iloc[0] for s in SAMPLE_ORDER}
MARKER_OF = {s: POINT_MARKERS[i % len(POINT_MARKERS)] for i, s in enumerate(SAMPLE_ORDER)}
SHADES = {"D1MT": ["#08519C", "#3182BD", "#6BAED6"],
          "Untreated": ["#A63603", "#E6550D", "#FD8D3C"]}
COLOR_OF, _seen = {}, {c: 0 for c in CONDITION_ORDER}
for s in SAMPLE_ORDER:
    c = COND_OF[s]
    pal = SHADES.get(c, ["#999999"])
    COLOR_OF[s] = pal[_seen.get(c, 0) % len(pal)]
    _seen[c] = _seen.get(c, 0) + 1

SEP_PX = max(1, int(round(PEAK_SEPARATION_UM / GRID_UM)))

if edge_rows:
    edge = pd.DataFrame(edge_rows)
    sub("Edge correction: what it changed")
    print(f"    {'section':<12}{'edge %':>9}{'bg raw':>9}{'bg corr':>9}"
          f"{'peak raw':>10}{'peak corr':>11}{'>thr raw':>10}{'>thr corr':>11}")
    print("    " + "-" * 81)
    for _, r in edge.iterrows():
        print(f"    {r['sample_id']:<12}{r['pct_tissue_edge_affected']:>8.1f}%"
              f"{r['bg_uncorrected']:>9.0f}{r['bg_corrected']:>9.0f}"
              f"{r['peak_uncorrected']:>10.0f}{r['peak_corrected']:>11.0f}"
              f"{r['pct_above_burden_uncorrected']:>9.2f}%"
              f"{r['pct_above_burden_corrected']:>10.2f}%")
    write_csv(edge, "32_density_edge_correction.csv")


# %% Cell 4 - DEFINITION 1: absolute burden (unchanged from revision 4)
# =============================================================================

banner("DEFINITION 1 - ABSOLUTE BURDEN (unchanged from revision 4)")

print(f"    Fixed threshold {BURDEN_THRESHOLD:.0f} myeloid/mm^2 applied "
      f"identically to every section.")
print("    A near-zero result in the treated arm is the finding, not a "
      "detection failure.\n")

burden_rows, burden_labels = [], {}
for s in SAMPLE_ORDER:
    g, dens = geom[s], dens_of[s]
    lab, recs = detect_burden(dens, g["occ"], g["gx"], g["gy"],
                              BURDEN_THRESHOLD, BURDEN_MIN_AREA_UM2,
                              BURDEN_MIN_CELLS, BURDEN_FILL_HOLES)
    burden_labels[s] = lab
    total_area = sum(r["area_um2"] for r in recs)
    for r in recs:
        r.update({"sample_id": s, "condition": COND_OF[s],
                  "pct_of_tissue_area": 100.0 * r["area_um2"] / g["tissue_area_um2"]})
        burden_rows.append(r)
    pct_tissue = 100.0 * total_area / g["tissue_area_um2"]
    pct_above = 100.0 * float((dens[g["occ"]] >= BURDEN_THRESHOLD).mean())
    print(f"    {s:<12} {COND_OF[s]:<10} regions={len(recs):>3}  "
          f"burden={pct_tissue:>5.1f}%  "
          f"(tissue above threshold before filters: {pct_above:>5.2f}%)")

burden = pd.DataFrame(burden_rows)
if len(burden):
    write_csv(burden, "33_burden_regions_absolute.csv")

burden_summary = pd.DataFrame([{
    "sample_id": s, "animal_id": short_label(s), "condition": COND_OF[s],
    "scan_id": geom[s]["scan"], "slide_position_rank": geom[s]["position"],
    "tissue_area_mm2": geom[s]["tissue_area_um2"] / 1e6,
    "n_regions": int((burden["sample_id"] == s).sum()) if len(burden) else 0,
    "burden_pct": (float(burden.loc[burden["sample_id"] == s,
                                    "pct_of_tissue_area"].sum())
                   if len(burden) else 0.0),
    "background_density": geom[s]["bg"],
    "background_stat": BACKGROUND_STAT,
    "background_density_median": geom[s]["bg_median"],
    "peak_density": geom[s]["peak"],
    "peak_fold_over_background": geom[s]["peak"] / geom[s]["bg"],
} for s in SAMPLE_ORDER])
sub("Burden summary")
print(burden_summary.to_string(index=False))
write_csv(burden_summary, "33b_burden_summary.csv")


# %% Cell 5 - DEFINITIONS 2 and 3, BALT, and the cell assignment
# =============================================================================

banner("FINAL DETECTION AND CELL ASSIGNMENT")

foci_rows, complex_rows, balt_rows, balt_alt_rows, masks = [], [], [], [], {}
method_used = "none"
diag_rows = []

for s in SAMPLE_ORDER:
    d = cells[s]
    g, dens = geom[s], dens_of[s]
    gx, gy, shape, occ = g["gx"], g["gy"], g["shape"], g["occ"]
    x0, y0, bg = g["x0"], g["y0"], g["bg"]
    pheno = d["pheno"].to_numpy()

    labels, recs, method, diag = detect_foci(
        dens, occ, gx, gy, bg, FOCUS_FOLD_OVER_BACKGROUND,
        FOCUS_HALF_MAX_FRACTION, FOCUS_MIN_AREA_UM2, FOCUS_MIN_CELLS, SEP_PX,
        min_peak=MIN_PEAK_DENSITY, fill_holes=FOCUS_FILL_HOLES,
        deduplicate=DEDUPLICATE_PLATEAU_PEAKS)
    if method != "none":
        method_used = method
    n_foci = len(recs)

    clabels, crecs, cdiag = detect_complexes(dens, occ, gx, gy, bg)
    n_complex = len(crecs)

    diag.update({"sample_id": s, "condition": COND_OF[s], "n_foci": n_foci,
                 "n_complexes": n_complex,
                 **{f"complex_{k}": v for k, v in cdiag.items()}})
    diag_rows.append(diag)

    core_mask = labels > 0
    dil_px = int(np.ceil(CUFF_WIDTH_UM / GRID_UM))
    cuff_mask = (ndi.binary_dilation(core_mask, iterations=dil_px) & occ
                 & ~core_mask)
    if core_mask.any():
        _, inds = ndi.distance_transform_edt(~core_mask, return_indices=True)
        cuff_owner = labels[inds[0], inds[1]]
        dout = ndi.distance_transform_edt(~core_mask) * GRID_UM
    else:
        cuff_owner = np.zeros_like(labels)
        dout = np.zeros(shape, dtype=float)

    # radial: 0 at core centre, 1 at core boundary, 1 to 2 across the cuff
    radial = np.full(shape, np.nan, dtype=float)
    for i in range(1, n_foci + 1):
        m = labels == i
        dt = ndi.distance_transform_edt(m)
        mx = dt.max()
        if mx > 0:
            radial[m] = 1.0 - (dt[m] / mx)
    if core_mask.any():
        radial[cuff_mask] = 1.0 + np.clip(dout[cuff_mask] / CUFF_WIDTH_UM, 0, 1)

    brelab, n_balt = detect_lymphoid(
        gx, gy, shape, occ, pheno, [BALT_B], [BALT_T],
        BALT_B_THRESHOLD, BALT_T_THRESHOLD, BALT_BANDWIDTH_UM,
        BALT_MIN_AREA_UM2, BALT_MIN_B_CELLS)

    masks[s] = {"x0": x0, "y0": y0, "shape": shape, "core": labels,
                "cuff": cuff_mask, "balt": brelab, "occ": occ, "dens": dens,
                "burden": burden_labels[s], "complex": clabels}

    # ---- per-cell assignment: REVISION 4 SCHEMA, plus additive rev5 columns
    cell_core = labels[gy, gx]
    cell_cuff = np.where(cuff_mask[gy, gx], cuff_owner[gy, gx], 0)
    out = d.copy()
    out["sample_id"] = s
    out["scan_id"] = g["scan"]
    out["slide_position_rank"] = g["position"]
    out["focus_id"] = np.where(cell_core > 0, cell_core, cell_cuff)
    out["region"] = np.where(cell_core > 0, "core",
                             np.where(cell_cuff > 0, "cuff", "interstitium"))
    out["radial_pos"] = radial[gy, gx]
    out["dist_to_focus_um"] = dout[gy, gx]
    out["burden_region_id"] = burden_labels[s][gy, gx]
    out["in_burden_region"] = burden_labels[s][gy, gx] > 0
    out["balt_id"] = brelab[gy, gx]
    out["myeloid_density"] = dens[gy, gx]
    out["section_background_density"] = bg
    # NEW in revision 5, additive only
    out["focus_core_id"] = cell_core
    out["complex_id"] = clabels[gy, gx]
    out["in_complex"] = clabels[gy, gx] > 0
    out["is_detection_pool"] = np.isin(pheno, POOL)
    if WRITE_CELL_ASSIGNMENTS:
        write_csv(out, f"{s}_cell_structures.csv", directory=CELL_DIR)

    # ---- per-focus metrics: REVISION 4 SCHEMA -----------------------------
    for r in recs:
        i = r["focus_id"]
        sel_core = cell_core == i
        sel_cuff = cell_cuff == i
        core_cells = out.loc[sel_core]
        all_cells = out.loc[sel_core | sel_cuff]
        mac = core_cells.loc[core_cells["pheno"].isin(CD68_LINEAGE)]
        ys, xs = np.nonzero(labels == i)
        inside_complex = clabels[ys, xs]
        rec = dict(r)
        rec.update({
            "sample_id": s, "condition": COND_OF[s], "scan_id": g["scan"],
            "centroid_x_um": float(x0 + (xs.mean() - 1) * GRID_UM),
            "centroid_y_um": float(y0 + (ys.mean() - 1) * GRID_UM),
            "pct_of_tissue_area": 100.0 * r["area_um2"] / g["tissue_area_um2"],
            "n_cells_core": int(sel_core.sum()),
            "n_cells_cuff": int(sel_cuff.sum()),
            # consistency guard: n_cells comes from the label array,
            # n_cells_core from the per-cell assignment. They must agree.
            "geometry_matches_assignment": bool(int(sel_core.sum()) == r["n_cells"]),
            "n_macrophages_core": int(len(mac)),
            "n_ido1_pos_core": int((core_cells["pheno"] == IDO1_POS_PHENO).sum()),
            "pct_ido1_pos_of_mac_core": (
                100.0 * float((core_cells["pheno"] == IDO1_POS_PHENO).sum()) / len(mac)
                if len(mac) else np.nan),
            "median_ido1_macrophages": float(mac["IDO1"].median()) if len(mac) else np.nan,
            "median_hk3_macrophages": float(mac["3-Hydroxykynurenine"].median())
            if len(mac) else np.nan,
            "median_inos_macrophages": float(mac["iNOS"].median()) if len(mac) else np.nan,
            "median_arg1_macrophages": float(mac["Arginase-1"].median())
            if len(mac) else np.nan,
            "median_ifng_core": float(core_cells["IFNG"].median()) if len(core_cells) else np.nan,
            "overlaps_burden_region": bool((core_cells["in_burden_region"]).any()),
            "pct_core_in_burden_region": 100.0 * float(core_cells["in_burden_region"].mean())
            if len(core_cells) else np.nan,
            # NEW: which complex this focus sits inside
            "complex_id": (int(np.bincount(inside_complex[inside_complex > 0]).argmax())
                           if (inside_complex > 0).any() else 0),
        })
        for ph in PHENOTYPE_ORDER:
            n_ph = int((all_cells["pheno"] == ph).sum())
            rec[f"n_{ascii_safe(ph)}"] = n_ph
            rec[f"pct_{ascii_safe(ph)}"] = (100.0 * n_ph / len(all_cells)
                                            if len(all_cells) else np.nan)
        foci_rows.append(rec)

    # ---- per-complex metrics, NEW -----------------------------------------
    for r in crecs:
        i = r["complex_id"]
        sel = (out["complex_id"] == i).to_numpy()
        cc = out.loc[sel]
        ys, xs = np.nonzero(clabels == i)
        inner = labels[clabels == i]
        n_inside = int(len(np.unique(inner[inner > 0])))
        mac = cc.loc[cc["pheno"].isin(CD68_LINEAGE)]
        rec = dict(r)
        rec.update({
            "sample_id": s, "condition": COND_OF[s], "scan_id": g["scan"],
            "centroid_x_um": float(x0 + (xs.mean() - 1) * GRID_UM),
            "centroid_y_um": float(y0 + (ys.mean() - 1) * GRID_UM),
            "pct_of_tissue_area": 100.0 * r["area_um2"] / g["tissue_area_um2"],
            "n_cells_assigned": int(sel.sum()),
            "n_foci_inside": n_inside,
            "n_macrophages": int(len(mac)),
            "pct_ido1_pos_of_mac": (
                100.0 * float((cc["pheno"] == IDO1_POS_PHENO).sum()) / len(mac)
                if len(mac) else np.nan),
            "median_ido1_macrophages": float(mac["IDO1"].median()) if len(mac) else np.nan,
            "overlaps_burden_region": bool(cc["in_burden_region"].any()),
            "pct_in_burden_region": 100.0 * float(cc["in_burden_region"].mean())
            if len(cc) else np.nan,
        })
        for ph in PHENOTYPE_ORDER:
            n_ph = int((cc["pheno"] == ph).sum())
            rec[f"n_{ascii_safe(ph)}"] = n_ph
            rec[f"pct_{ascii_safe(ph)}"] = (100.0 * n_ph / len(cc)
                                            if len(cc) else np.nan)
        complex_rows.append(rec)

    # ---- BALT metrics, UNCHANGED FROM REVISION 4 --------------------------
    b_outside = out.loc[(out["balt_id"] == 0) & (out["pheno"] == BALT_B), "CD21"]
    cd21_bg = float(b_outside.median()) if len(b_outside) else np.nan
    for i in range(1, n_balt + 1):
        sel = out["balt_id"] == i
        s_all = out.loc[sel]
        b_in = s_all.loc[s_all["pheno"] == BALT_B, "CD21"]
        area = (brelab == i).sum() * GRID_UM * GRID_UM
        ys, xs = np.nonzero(brelab == i)
        cd21_in = float(b_in.median()) if len(b_in) else np.nan
        ratio = (cd21_in / cd21_bg) if (np.isfinite(cd21_bg) and cd21_bg > 0) else np.nan
        balt_rows.append({
            "sample_id": s, "condition": COND_OF[s], "scan_id": g["scan"],
            "balt_id": i,
            "centroid_x_um": float(x0 + (xs.mean() - 1) * GRID_UM),
            "centroid_y_um": float(y0 + (ys.mean() - 1) * GRID_UM),
            "area_um2": float(area), "n_cells": int(sel.sum()),
            "n_b_cells": int((s_all["pheno"] == BALT_B).sum()),
            "n_helper_t": int((s_all["pheno"] == BALT_T).sum()),
            "n_t_lineage": int(s_all["pheno"].isin(T_LINEAGE).sum()),
            "n_plasma": int((s_all["pheno"] == PLASMA).sum()),
            "pct_b_cells": 100.0 * float((s_all["pheno"] == BALT_B).mean()) if len(s_all) else np.nan,
            "pct_helper_t": 100.0 * float((s_all["pheno"] == BALT_T).mean()) if len(s_all) else np.nan,
            "pct_t_lineage": 100.0 * float(s_all["pheno"].isin(T_LINEAGE).mean()) if len(s_all) else np.nan,
            "pct_plasma": 100.0 * float((s_all["pheno"] == PLASMA).mean()) if len(s_all) else np.nan,
            "median_cd21_b_cells": cd21_in,
            "median_cd21_b_cells_outside": cd21_bg,
            "cd21_ratio": ratio,
            "is_follicle": bool(np.isfinite(ratio) and ratio >= CD21_FOLLICLE_RATIO),
            "min_dist_to_focus_um": float(s_all["dist_to_focus_um"].min())
            if len(s_all) else np.nan,
            "pct_in_focus_core": 100.0 * float((s_all["region"] == "core").mean())
            if len(s_all) else np.nan,
            "pct_in_focus_cuff": 100.0 * float((s_all["region"] == "cuff").mean())
            if len(s_all) else np.nan,
            "pct_in_burden_region": 100.0 * float(s_all["in_burden_region"].mean())
            if len(s_all) else np.nan,
            "pct_in_complex": 100.0 * float(s_all["in_complex"].mean())
            if len(s_all) else np.nan,
        })

    if RUN_BALT_ALTERNATIVES:
        for name, b_set, t_set, b_thr, t_thr in BALT_ALT_GATES:
            alt, n_alt = detect_lymphoid(
                gx, gy, shape, occ, pheno, b_set, t_set, b_thr, t_thr,
                BALT_BANDWIDTH_UM, BALT_MIN_AREA_UM2, BALT_MIN_B_CELLS)
            n_foll = 0
            for i in range(1, n_alt + 1):
                m_alt = alt[gy, gx] == i
                b_in = out.loc[m_alt & (out["pheno"] == BALT_B).to_numpy(), "CD21"]
                if len(b_in) and np.isfinite(cd21_bg) and cd21_bg > 0:
                    if float(b_in.median()) / cd21_bg >= CD21_FOLLICLE_RATIO:
                        n_foll += 1
            balt_alt_rows.append({
                "sample_id": s, "condition": COND_OF[s], "gate": name,
                "n_candidates": n_alt, "n_cd21_follicles": n_foll,
                "area_mm2": float((alt > 0).sum() * GRID_UM * GRID_UM / 1e6),
            })

    exp = EXPERT_COUNTS.get(s, np.nan)
    flag = f"   <-- expert {int(exp)}" if np.isfinite(exp) and n_foci != exp else ""
    print(f"    {s:<12} foci={n_foci:>3}  complexes={n_complex:>3}  "
          f"BALT={n_balt:>3}  bg={bg:>6.0f}   peaks {diag['n_peaks_raw']}->"
          f"{diag['n_peaks_raw'] - diag['n_peaks_collapsed']} after dedup, "
          f"{diag['n_peaks_kept']} kept, {diag['n_dropped_min_peak']} below "
          f"min peak{flag}")
    del out
    gc.collect()

foci = pd.DataFrame(foci_rows)
comps = pd.DataFrame(complex_rows)
balts = pd.DataFrame(balt_rows)
balt_alt = pd.DataFrame(balt_alt_rows)
diagnostics = pd.DataFrame(diag_rows)
print(f"\n    Partition method used: {method_used}")
write_csv(diagnostics, "35b_detection_diagnostics.csv")

if len(foci):
    bad = foci.loc[~foci["geometry_matches_assignment"]]
    if len(bad):
        print(f"\n    WARNING: {len(bad)} focus record(s) where the label-array "
              f"cell count disagrees with the per-cell assignment.")
        print(bad[["sample_id", "focus_id", "n_cells", "n_cells_core"]].to_string(index=False))
    else:
        print("    Geometry check: every focus record matches its per-cell "
              "assignment.")
    write_csv(foci, "35_foci_structures_relative.csv")
if len(comps):
    write_csv(comps, "37_granuloma_complexes.csv")
if len(balts):
    write_csv(balts, "36_balt_structures.csv")
if len(balt_alt):
    write_csv(balt_alt, "36b_balt_alternative_gates.csv")

sub("Structures per section")
for s in SAMPLE_ORDER:
    f = foci.loc[foci["sample_id"] == s] if len(foci) else pd.DataFrame()
    c = comps.loc[comps["sample_id"] == s] if len(comps) else pd.DataFrame()
    b = balts.loc[balts["sample_id"] == s] if len(balts) else pd.DataFrame()
    nf = int(b["is_follicle"].sum()) if len(b) else 0
    if len(f):
        print(f"    {s:<12} {COND_OF[s]:<10} foci={len(f):>3} "
              f"(expert {EXPERT_COUNTS.get(s, 'na')})  complexes={len(c):>3}  "
              f"inscribed={f['max_inscribed_radius_um'].median():>5.0f} um  "
              f"shape={f['shape_ratio'].median():>4.2f}  "
              f"peak={f['peak_density'].median():>7.0f}  "
              f"fold={f['fold_over_background'].median():>5.1f}x  "
              f"BALT={len(b):>3} ({nf} CD21+)")
    else:
        print(f"    {s:<12} {COND_OF[s]:<10} foci=  0  complexes={len(c):>3}  "
              f"BALT={len(b):>3} ({nf} CD21+)")

if len(comps):
    sub("Confluence: peak-based foci inside each complex")
    gg = comps.groupby(["condition", "sample_id"]).agg(
        n_complexes=("complex_id", "size"),
        median_foci_inside=("n_foci_inside", "median"),
        max_foci_inside=("n_foci_inside", "max"),
        total_area_mm2=("area_um2", lambda v: v.sum() / 1e6))
    print(gg.to_string())
    print("\n    'foci inside' above 1 is confluence: several density peaks")
    print("    within one connected structure. That is the quantity the")
    print("    level-set definition exists to capture.")

if len(balt_alt):
    sub("BALT gate sensitivity (diagnostic, does not change balt_id)")
    piv = balt_alt.pivot_table(index="sample_id", columns="gate",
                               values="n_candidates")
    pivf = balt_alt.pivot_table(index="sample_id", columns="gate",
                                values="n_cd21_follicles")
    gates = [gname[0] for gname in BALT_ALT_GATES if gname[0] in piv.columns]
    print(f"    {'section':<14}" + "".join(f"{gn[:22]:>24}" for gn in gates))
    for s in SAMPLE_ORDER:
        if s not in piv.index:
            continue
        row = f"    {s:<14}"
        for gn in gates:
            row += f"{int(piv.loc[s, gn]):>13} ({int(pivf.loc[s, gn])} CD21+)"[:24].rjust(24)
        print(row)

sub("SCALE WARNING - carry this into every downstream comparison")
if len(foci):
    for c in CONDITION_ORDER:
        f = foci.loc[foci["condition"] == c]
        if not len(f):
            continue
        print(f"    {c:<12} n={len(f):>3}  peak density median "
              f"{f['peak_density'].median():>8.0f}  "
              f"equiv radius {f['equiv_radius_um'].median():>5.0f} um  "
              f"inscribed radius {f['max_inscribed_radius_um'].median():>5.0f} um  "
              f"fold {f['fold_over_background'].median():>5.1f}x")
    print("\n    Treated foci are defined relative to treated background. They")
    print("    are residual myeloid foci, NOT granulomas equivalent to the")
    print("    untreated lesions. Report peak density and fold alongside any")
    print("    architecture comparison so the two are never conflated.")
    print("\n    USE max_inscribed_radius_um, NOT equiv_radius_um, to convert a")
    print("    radial-unit effect into microns.")


# %% Cell 6 - figures
# =============================================================================

banner("FIGURES")


def base_panel(sid):
    d = cells[sid]
    n = len(d)
    idx = (rng.choice(n, size=MAX_POINTS_PER_SECTION, replace=False)
           if n > MAX_POINTS_PER_SECTION else np.arange(n))
    return d.iloc[idx]


def extent_of(sid):
    mk = masks[sid]
    ny, nx = mk["shape"]
    return [mk["x0"] - GRID_UM, mk["x0"] + (nx - 1) * GRID_UM,
            mk["y0"] - GRID_UM, mk["y0"] + (ny - 1) * GRID_UM]


def contour(ax, sid, arr, color, lw, ls="solid"):
    if not np.any(arr):
        return
    ax.contour(np.asarray(arr, dtype=float), levels=[0.5], colors=[color],
               linewidths=lw, extent=extent_of(sid), origin="lower",
               linestyles=ls)


def spatial_grid(stem, title, draw_fn, handles):
    fig, axes = plt.subplots(2, 3, figsize=(45, 26))
    axes = axes.ravel()
    for k, s in enumerate(SAMPLE_ORDER):
        draw_fn(axes[k], s)
        axes[k].set_aspect("equal"); axes[k].set_xticks([]); axes[k].set_yticks([])
        for sp in axes[k].spines.values():
            sp.set_edgecolor(CONDITION_COLORS.get(COND_OF[s], "#999"))
            sp.set_linewidth(6)
    for k in range(len(SAMPLE_ORDER), len(axes)):
        axes[k].axis("off")
    fig.legend(handles=handles, loc="lower center", ncol=min(6, len(handles)),
               frameon=False, fontsize=FONT_SIZE_LEGEND - 10,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(title, y=1.02, fontsize=FONT_SIZE_TITLE - 4)
    save_fig(fig, stem)


def draw_myeloid(ax, s):
    d = base_panel(s)
    other = ~d["pheno"].isin(POOL)
    ax.scatter(d.loc[other, "x"], d.loc[other, "y"], s=0.6,
               color=BACKGROUND_COLOR, linewidths=0, rasterized=True)
    for ph in POOL:
        m = d["pheno"] == ph
        if m.sum():
            ax.scatter(d.loc[m, "x"], d.loc[m, "y"], s=5,
                       color=PHENOTYPE_COLORS[ph], linewidths=0, rasterized=True)


def draw_burden(ax, s):
    draw_myeloid(ax, s)
    contour(ax, s, masks[s]["burden"] > 0, BURDEN_OUTLINE, 4)
    n = int((burden["sample_id"] == s).sum()) if len(burden) else 0
    pct = (float(burden.loc[burden["sample_id"] == s, "pct_of_tissue_area"].sum())
           if len(burden) else 0.0)
    ax.set_title(f"{short_label(s)} ({COND_OF[s]})   {n} regions, {pct:.1f}% burden",
                 fontsize=FONT_SIZE_TITLE - 12)


spatial_grid("F30_burden_absolute_overlay",
             f"DEFINITION 1: absolute threshold {BURDEN_THRESHOLD:.0f} myeloid/mm2 "
             "applied identically to every section\n"
             "The near-absence in the treated arm is the result",
             draw_burden,
             [Patch(facecolor=PHENOTYPE_COLORS[p], label=p) for p in POOL]
             + [Line2D([0], [0], color=BURDEN_OUTLINE, linewidth=4,
                       label="granuloma-density region")])


def draw_foci(ax, s):
    draw_myeloid(ax, s)
    contour(ax, s, masks[s]["cuff"], CUFF_OUTLINE, 2.5)
    contour(ax, s, masks[s]["core"] > 0, FOCUS_OUTLINE, 4)
    f = foci.loc[foci["sample_id"] == s] if len(foci) else pd.DataFrame()
    for _, r in f.iterrows():
        ax.text(r["centroid_x_um"], r["centroid_y_um"],
                f"{r['fold_over_background']:.0f}x",
                ha="center", va="center", fontsize=FONT_SIZE_ANNOT - 10,
                color=FLAG_COLOR, fontweight="bold")
    ax.set_title(f"{short_label(s)} ({COND_OF[s]})   {len(f)} foci, "
                 f"expert {EXPERT_COUNTS.get(s, 'na')}",
                 fontsize=FONT_SIZE_TITLE - 12)


spatial_grid("F33_foci_relative_overlay",
             "DEFINITION 2: peak-based foci, relative to each section's own "
             f"{BACKGROUND_STAT} background\nLabels give fold over that background",
             draw_foci,
             [Patch(facecolor=PHENOTYPE_COLORS[p], label=p) for p in POOL]
             + [Line2D([0], [0], color=FOCUS_OUTLINE, linewidth=4, label="focus core"),
                Line2D([0], [0], color=CUFF_OUTLINE, linewidth=2.5,
                       label=f"cuff (+{CUFF_WIDTH_UM:.0f} um)")])


def draw_complex(ax, s):
    draw_myeloid(ax, s)
    contour(ax, s, masks[s]["complex"] > 0, COMPLEX_OUTLINE, 4)
    contour(ax, s, masks[s]["core"] > 0, FOCUS_OUTLINE, 2.5, ls="dashed")
    c = comps.loc[comps["sample_id"] == s] if len(comps) else pd.DataFrame()
    nf = int(c["n_foci_inside"].max()) if len(c) else 0
    ax.set_title(f"{short_label(s)} ({COND_OF[s]})   {len(c)} complexes, "
                 f"up to {nf} foci in one", fontsize=FONT_SIZE_TITLE - 12)


spatial_grid("F38_complexes_overlay",
             "DEFINITION 3: granuloma complexes, connected components above "
             f"{COMPLEX_LEVEL_FOLD:g}x background\n"
             "Green = complex, dashed black = peak-based foci inside it",
             draw_complex,
             [Line2D([0], [0], color=COMPLEX_OUTLINE, linewidth=4, label="complex"),
              Line2D([0], [0], color=FOCUS_OUTLINE, linewidth=2.5,
                     linestyle="--", label="focus core")])

# ---- F31 burden summary -----------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(42, 13))
xp = np.arange(len(SAMPLE_ORDER))
cols = [CONDITION_COLORS[COND_OF[s]] for s in SAMPLE_ORDER]
xl = [short_label(s) for s in SAMPLE_ORDER]

ax = axes[0]
ax.bar(xp, burden_summary["n_regions"], color=cols, edgecolor="#FFFFFF",
       linewidth=2, zorder=3)
ax.set_xticks(xp); ax.set_xticklabels(xl, rotation=45, ha="right",
                                      fontsize=FONT_SIZE_TICK - 8)
ax.set_ylabel("Regions")
ax.set_title("Granuloma-density regions", fontsize=FONT_SIZE_TITLE - 10)
style_axes(ax)

ax = axes[1]
ax.bar(xp, burden_summary["burden_pct"], color=cols, edgecolor="#FFFFFF",
       linewidth=2, zorder=3)
for i, v in enumerate(burden_summary["burden_pct"]):
    ax.text(i, v + 0.6, f"{v:.1f}%", ha="center", fontsize=FONT_SIZE_ANNOT - 8)
ax.set_xticks(xp); ax.set_xticklabels(xl, rotation=45, ha="right",
                                      fontsize=FONT_SIZE_TICK - 8)
ax.set_ylabel("% of tissue area")
ax.set_title("Lesion burden", fontsize=FONT_SIZE_TITLE - 10)
style_axes(ax)

ax = axes[2]
for i, s in enumerate(SAMPLE_ORDER):
    ax.scatter([i], [geom[s]["bg"]], s=420, color="#999999", marker="_",
               linewidth=6, zorder=3)
    ax.scatter([i], [geom[s]["peak"]], s=420, color=COLOR_OF[s],
               marker=MARKER_OF[s], edgecolor="#FFFFFF", linewidth=2, zorder=4)
    ax.plot([i, i], [geom[s]["bg"], geom[s]["peak"]], color=COLOR_OF[s],
            linewidth=4, alpha=0.6, zorder=2)
ax.axhline(BURDEN_THRESHOLD, color=BURDEN_OUTLINE, linestyle="--", linewidth=3.5)
ax.axhline(MIN_PEAK_DENSITY, color=FLAG_COLOR, linestyle=":", linewidth=3.5)
ax.text(len(SAMPLE_ORDER) - 0.5, BURDEN_THRESHOLD * 1.15,
        f"burden {BURDEN_THRESHOLD:.0f}", ha="right",
        fontsize=FONT_SIZE_ANNOT - 10, color=BURDEN_OUTLINE)
ax.text(len(SAMPLE_ORDER) - 0.5, MIN_PEAK_DENSITY * 0.72,
        f"min peak {MIN_PEAK_DENSITY:.0f}", ha="right",
        fontsize=FONT_SIZE_ANNOT - 10, color=FLAG_COLOR)
ax.set_yscale("log")
ax.set_xticks(xp); ax.set_xticklabels(xl, rotation=45, ha="right",
                                      fontsize=FONT_SIZE_TICK - 8)
ax.set_ylabel("Myeloid density (cells / mm$^2$)")
ax.set_title(f"{BACKGROUND_STAT} background to peak",
             fontsize=FONT_SIZE_TITLE - 12)
style_axes(ax)
fig.suptitle("Absolute burden, and why one threshold cannot serve both arms",
             y=1.03, fontsize=FONT_SIZE_TITLE - 4)
save_fig(fig, "F31_burden_summary")

# ---- F34 the scale figure ---------------------------------------------------
if len(foci):
    fig, axes = plt.subplots(1, 3, figsize=(42, 13))
    for ax, col, lab, logy in [
            (axes[0], "peak_density", "Peak myeloid density (cells / mm$^2$)", True),
            (axes[1], "fold_over_background", "Fold over own background", False),
            (axes[2], "shape_ratio", "Equivalent radius / inscribed radius", False)]:
        for i, c in enumerate(CONDITION_ORDER):
            f = foci.loc[foci["condition"] == c]
            if not len(f):
                continue
            j = rng.uniform(-0.14, 0.14, size=len(f))
            ax.scatter(np.full(len(f), i) + j, f[col], s=320,
                       color=CONDITION_COLORS[c], edgecolor="#FFFFFF",
                       linewidth=2, zorder=3)
            ax.hlines(f[col].median(), i - 0.3, i + 0.3, color="#000000",
                      linewidth=4)
        if logy:
            ax.set_yscale("log")
            ax.axhline(MIN_PEAK_DENSITY, color=FLAG_COLOR, linestyle=":",
                       linewidth=3)
            ax.axhline(BURDEN_THRESHOLD, color=BURDEN_OUTLINE, linestyle="--",
                       linewidth=3)
        if col == "shape_ratio":
            ax.axhline(1.0, color="#000000", linestyle="--", linewidth=2.5)
        ax.set_xticks(range(len(CONDITION_ORDER)))
        ax.set_xticklabels(CONDITION_ORDER)
        ax.set_ylabel(lab, fontsize=FONT_SIZE_BASE - 6)
        style_axes(ax)
    fig.suptitle("Treated foci are prominent within their own tissue but operate "
                 "at a lower absolute scale\nDotted = minimum peak density, "
                 "dashed = burden threshold", y=1.05,
                 fontsize=FONT_SIZE_TITLE - 8)
    save_fig(fig, "F34_focus_scale_comparison")


# ---- F35 IDO1 capture validation -------------------------------------------
def draw_capture(ax, s):
    d = base_panel(s)
    other = d["pheno"] != IDO1_POS_PHENO
    ax.scatter(d.loc[other, "x"], d.loc[other, "y"], s=0.6,
               color=BACKGROUND_COLOR, linewidths=0, rasterized=True)
    m = d["pheno"] == IDO1_POS_PHENO
    if m.sum():
        ax.scatter(d.loc[m, "x"], d.loc[m, "y"], s=14,
                   color=PHENOTYPE_COLORS[IDO1_POS_PHENO], linewidths=0,
                   rasterized=True)
    contour(ax, s, masks[s]["core"] > 0, FOCUS_OUTLINE, 4)
    contour(ax, s, masks[s]["complex"] > 0, COMPLEX_OUTLINE, 2.5)
    contour(ax, s, masks[s]["burden"] > 0, BURDEN_OUTLINE, 3, ls="dashed")
    full = cells[s]
    ins = full.loc[full["pheno"] == IDO1_POS_PHENO]
    pct_f = pct_b = pct_c = np.nan
    if len(ins):
        mk = masks[s]
        ggx = np.clip(((ins["x"].to_numpy() - mk["x0"]) / GRID_UM).astype(int) + 1,
                      0, mk["shape"][1] - 1)
        ggy = np.clip(((ins["y"].to_numpy() - mk["y0"]) / GRID_UM).astype(int) + 1,
                      0, mk["shape"][0] - 1)
        pct_f = 100.0 * float((mk["core"][ggy, ggx] > 0).mean())
        pct_b = 100.0 * float((mk["burden"][ggy, ggx] > 0).mean())
        pct_c = 100.0 * float((mk["complex"][ggy, ggx] > 0).mean())
    ax.set_title(f"{short_label(s)} ({COND_OF[s]})\n"
                 f"IDO1+ in foci {pct_f:.0f}%, complexes {pct_c:.0f}%, "
                 f"burden {pct_b:.0f}%", fontsize=FONT_SIZE_TITLE - 14)


spatial_grid("F35_ido1_capture_validation",
             "Independent validation: IDO1+ macrophages against boundaries drawn "
             "without IDO1\nSolid black = foci, green = complexes, dashed = burden",
             draw_capture,
             [Patch(facecolor=PHENOTYPE_COLORS[IDO1_POS_PHENO], label=IDO1_POS_PHENO),
              Line2D([0], [0], color=FOCUS_OUTLINE, linewidth=4, label="focus"),
              Line2D([0], [0], color=COMPLEX_OUTLINE, linewidth=3, label="complex"),
              Line2D([0], [0], color=BURDEN_OUTLINE, linewidth=3, linestyle="--",
                     label="burden region")])

# ---- F36 BALT ---------------------------------------------------------------
_show = [BALT_B, BALT_T, PLASMA]


def draw_balt(ax, s):
    d = base_panel(s)
    other = ~d["pheno"].isin(_show)
    ax.scatter(d.loc[other, "x"], d.loc[other, "y"], s=0.6,
               color=BACKGROUND_COLOR, linewidths=0, rasterized=True)
    for ph, sz in zip(_show, [12, 14, 6]):
        m = d["pheno"] == ph
        if m.sum():
            ax.scatter(d.loc[m, "x"], d.loc[m, "y"], s=sz,
                       color=PHENOTYPE_COLORS[ph], linewidths=0, rasterized=True)
    contour(ax, s, masks[s]["balt"] > 0, BALT_OUTLINE, 4, ls="dashed")
    contour(ax, s, masks[s]["core"] > 0, FOCUS_OUTLINE, 4)
    b = balts.loc[balts["sample_id"] == s] if len(balts) else pd.DataFrame()
    nf = int(b["is_follicle"].sum()) if len(b) else 0
    ax.set_title(f"{short_label(s)} ({COND_OF[s]})   {len(b)} candidates, {nf} CD21+",
                 fontsize=FONT_SIZE_TITLE - 12)


spatial_grid("F36_balt_candidates",
             "BALT candidates from joint B cell and helper T cell density\n"
             "Dashed blue = candidate, solid black = myeloid focus",
             draw_balt,
             [Patch(facecolor=PHENOTYPE_COLORS[p], label=p) for p in _show]
             + [Line2D([0], [0], color=BALT_OUTLINE, linewidth=4, linestyle="--",
                       label="BALT candidate"),
                Line2D([0], [0], color=FOCUS_OUTLINE, linewidth=4, label="focus")])

# ---- F37 BALT classification ------------------------------------------------
if len(balts):
    fig, axes = plt.subplots(1, 3, figsize=(45, 13))
    ax = axes[0]
    for s in SAMPLE_ORDER:
        b = balts.loc[balts["sample_id"] == s]
        if len(b):
            ax.scatter(b["pct_b_cells"], b["cd21_ratio"], s=340, color=COLOR_OF[s],
                       marker=MARKER_OF[s], edgecolor="#FFFFFF", linewidth=2,
                       zorder=3)
    ax.axhline(CD21_FOLLICLE_RATIO, color="#000000", linestyle="--", linewidth=3.5)
    ax.set_xlabel("% B cells in candidate")
    ax.set_ylabel("CD21 on B cells, inside / outside")
    ax.set_title("Follicle classification", fontsize=FONT_SIZE_TITLE - 8)
    style_axes(ax)
    ax.legend(handles=[Line2D([0], [0], color=COLOR_OF[s], marker=MARKER_OF[s],
                              markersize=16, linestyle="none",
                              label=f"{short_label(s)} ({COND_OF[s]})")
                       for s in SAMPLE_ORDER],
              frameon=False, fontsize=FONT_SIZE_LEGEND - 14, loc="best")

    ax = axes[1]
    for s in SAMPLE_ORDER:
        b = balts.loc[balts["sample_id"] == s]
        if len(b):
            ax.scatter(b["min_dist_to_focus_um"], b["pct_plasma"], s=340,
                       color=COLOR_OF[s], marker=MARKER_OF[s],
                       edgecolor="#FFFFFF", linewidth=2, zorder=3)
    ax.set_xlabel("Distance to nearest myeloid focus (um)")
    ax.set_ylabel("% plasma cells in candidate")
    ax.set_title("Proximity to foci and plasma content",
                 fontsize=FONT_SIZE_TITLE - 14)
    style_axes(ax)

    ax = axes[2]
    for s in SAMPLE_ORDER:
        b = balts.loc[balts["sample_id"] == s]
        if len(b):
            ax.scatter(b["pct_in_focus_core"] + b["pct_in_focus_cuff"],
                       b["cd21_ratio"], s=340, color=COLOR_OF[s],
                       marker=MARKER_OF[s], edgecolor="#FFFFFF", linewidth=2,
                       zorder=3)
    ax.axhline(CD21_FOLLICLE_RATIO, color="#000000", linestyle="--", linewidth=3.5)
    ax.set_xlabel("% of candidate inside a focus core or cuff")
    ax.set_ylabel("CD21 ratio")
    ax.set_title("Free-standing or granuloma-associated?",
                 fontsize=FONT_SIZE_TITLE - 16)
    style_axes(ax)
    fig.suptitle("BALT candidate characterization", y=1.04,
                 fontsize=FONT_SIZE_TITLE - 4)
    save_fig(fig, "F37_balt_classification")


# %% Cell 7 - schema guard
# =============================================================================

banner("SCHEMA CHECK")

print("Revision 5 renamed nothing. This asserts it, because a silent rename is")
print("how a downstream script ends up with a column of zeros instead of an")
print("error. Revision 5 introduced exactly this bug once already.\n")

schema_rows = []
ok_all = True
probe = sorted(glob.glob(os.path.join(CELL_DIR, "*_cell_structures.csv")))
if probe:
    pc = pd.read_csv(probe[0], nrows=5)
    miss_cell = [c for c in REQUIRED_CELL_COLS if c not in pc.columns]
    ok_all &= not miss_cell
    for c in REQUIRED_CELL_COLS:
        schema_rows.append({"table": "cell_assignments", "column": c,
                            "present": c in pc.columns})
    print(f"    cell_assignments : "
          f"{len(REQUIRED_CELL_COLS) - len(miss_cell)}"
          f"/{len(REQUIRED_CELL_COLS)} required columns present")
    if miss_cell:
        print(f"      MISSING: {miss_cell}")
    print(f"      additional (rev5): "
          f"{[c for c in pc.columns if c not in REQUIRED_CELL_COLS]}")
else:
    ok_all = False
    print("    ERROR: no cell assignment files were written.")

if len(foci):
    miss_foci = [c for c in REQUIRED_FOCI_COLS if c not in foci.columns]
    ok_all &= not miss_foci
    for c in REQUIRED_FOCI_COLS:
        schema_rows.append({"table": "35_foci", "column": c,
                            "present": c in foci.columns})
    print(f"    table 35         : "
          f"{len(REQUIRED_FOCI_COLS) - len(miss_foci)}"
          f"/{len(REQUIRED_FOCI_COLS)} required columns present")
    if miss_foci:
        print(f"      MISSING: {miss_foci}")
else:
    ok_all = False
    print("    ERROR: table 35 is empty.")

write_csv(pd.DataFrame(schema_rows), "38_schema_check.csv")
if ok_all:
    print("\n    PASS. Scripts 06, 07 and 08 need only their IN_DIR repointed")
    print("    from structures_rev4 to structures_rev5.")
else:
    print("\n    FAIL. Do NOT run downstream scripts until this passes.")


# %% Cell 8 - wrap up
# =============================================================================

banner("SUMMARY")
print(f"Partition method     : {method_used}")
print(f"Detection pool       : {', '.join(POOL)}")
print(f"Sections written     : {len(SAMPLE_ORDER)}")
print(f"Burden regions (abs) : {len(burden)}")
print(f"Foci (peak-based)    : {len(foci)}")
print(f"Complexes (level set): {len(comps)}")
print(f"BALT candidates      : {len(balts)}"
      f"  ({int(balts['is_follicle'].sum()) if len(balts) else 0} CD21+)")
print(f"Schema check         : {'PASS' if ok_all else 'FAIL'}")
if LEAVE_ONE_OUT:
    print(f"LEAVE-ONE-OUT        : {LEAVE_ONE_OUT} excluded. SENSITIVITY RUN.")
if skipped:
    print(f"Sections SKIPPED     : {len(skipped)}")
    for q, why in skipped:
        print(f"    {q}: {why}")
    print("  These have NO cell assignment file. Downstream scripts will simply")
    print("  not see them, so resolve this before re-running 06.")
for c in CONDITION_ORDER:
    bsel = burden_summary.loc[burden_summary["condition"] == c]
    fsel = foci.loc[foci["condition"] == c] if len(foci) else pd.DataFrame()
    csel = comps.loc[comps["condition"] == c] if len(comps) else pd.DataFrame()
    print(f"  {c:<12} burden {bsel['burden_pct'].mean():>5.1f}% mean, "
          f"{len(fsel):>3} foci, {len(csel):>3} complexes")

sub("Read in this order")
print("  1. table 38     : the schema check. Nothing downstream runs until it")
print("     passes.")
print("  2. table 32     : what the edge correction changed.")
print("  3. F31 right    : background to peak per section, with both the burden")
print("     threshold and the new minimum peak density drawn.")
print("  4. F33 and F38  : foci against complexes. Do the boundaries match the")
print("     slides, and does the complex definition capture confluence?")
print("  5. F35          : IDO1+ capture, the independent validation.")

sub("Downstream must be rerun")
print("  Scripts 06, 07 and 08 read cell_assignments/ or table 35. Repoint")
print("  IN_DIR from structures_rev4 to structures_rev5, re-run all three, then")
print("  09. Quote nothing from the revision 4 tables.")
print()
print("  Script 05 stays FROZEN. Its radial centre-of-mass analysis is")
print("  superseded by script 08's composition-balanced centring, which is the")
print("  correction that narrowed Finding 3 from lymphocytes to B lineage. Its")
print("  nearest-neighbour work is superseded by script 06, and its iNOS and")
print("  Arginase-1 result was dropped when the CD3e x CD20 control failed.")

sub("Expect these to move")
print("  Treated focus counts. Revision 4 reported three foci on G3_43106 with")
print("  peak densities 3,243, 1,565 and 1,333 and IDO1+ macrophage fractions")
print("  0.139, 0.000 and 0.000. Only the first clears the minimum peak")
print("  density, so two of them were vascular.")
print()
print("  Treated foci per animal fall to roughly 1, 1, 2, so two animals then")
print("  contribute a single focus each. That makes the animal variance")
print("  unidentifiable in a focus-level mixed model, because the random")
print("  intercept and the residual become the same number. Improving detection")
print("  and improving statistical power were never the same goal. Lead with")
print("  the animal-level comparison and the exact randomization test.")

banner("DONE")
sys.stdout = _tee.terminal
_tee.close()
