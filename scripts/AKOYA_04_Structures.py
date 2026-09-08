#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - STRUCTURE DEFINITION (dual definition)
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 04 of the AKOYA analysis series. REVISION 4.

WHY TWO DEFINITIONS
    The density landscape settled this. The densest point in any D1MT section is
    about 4,000 myeloid cells/mm^2. The densest in any untreated section is about
    19,600. Untreated lesions run three to five fold denser at peak, so no single
    absolute threshold can serve both arms: set it for untreated and the treated
    arm returns zero, set it for treated and the untreated sections flood.

    But every treated section does contain one or two foci sitting 8 to 16 fold
    above its own tissue background, with a sharp drop to the next peak. Those
    are real, discrete, isolated foci operating at a different absolute scale. So
    we use two definitions for two different questions:

    DEFINITION 1 - ABSOLUTE, for BURDEN
        Fixed density threshold applied identically to every section. Answers
        "how much granuloma-density tissue is there". The near-zero result in the
        treated arm IS the finding, not a failure of detection. The threshold is
        hard-coded on purpose and is not swept.

    DEFINITION 2 - RELATIVE, for ARCHITECTURE
        Local maxima that rise a set fold above their OWN section's median, each
        grown out to its own half-maximum, partitioned by watershed so adjacent
        foci do not merge. Scale-free per focus. Answers "how is a focus
        organised internally", which is the question the radial coordinate was
        built for.

    The two are never mixed. Every detected focus carries its fold-over-
    background and peak density so a treated focus can never be silently
    presented as equivalent to an untreated lesion.

CRITICAL DESIGN CONSTRAINT (unchanged)
    Neither definition uses IDO1. Both run on pooled myeloid density. IDO1 is
    measured about each structure afterwards. F35 validates capture independently.

WHAT CHANGED IN REVISION 4 (and why)

    1. THE DENSITY MAP IS EDGE-CORRECTED.
        density_map smoothed a zero-padded count array with mode="constant" and
        never renormalised by the tissue mask. Every cell within roughly one
        bandwidth of a tissue boundary therefore sat in a neighbourhood that was
        partly empty non-tissue space, and its local density was pulled down
        toward zero by however much of the kernel fell outside the tissue. A
        focus at the edge of a section is systematically under-detected, and
        because the arms differ in how fragmented their tissue is, the size of
        that bias can differ by arm.

        The fix is normalized convolution: smooth the counts, smooth the tissue
        occupancy with the same kernel, divide. Inside the tissue the occupancy
        weight is 1 and nothing changes; at a boundary it is less than 1 and the
        density is restored to what it would be if the tissue continued.
        EDGE_MIN_WEIGHT floors the denominator so thin protrusions cannot
        explode.

        Both maps are computed. Detection runs on the corrected one, and table
        32 reports the uncorrected numbers alongside so the size of the change
        is a measured quantity rather than an assumption. Expect burden
        percentages to move.

    2. EVERY FOCUS NOW CARRIES ITS MAXIMUM INSCRIBED RADIUS.
        The radial coordinate normalises by the maximum inscribed radius, that
        is dt.max() inside each focus, but downstream work converted radial
        units to microns using equiv_radius_um. For any non-circular shape the
        inscribed radius is smaller than the equivalent radius, so that
        conversion overstates the effect in microns. max_inscribed_radius_um and
        shape_ratio (equivalent over inscribed, 1.0 for a disc) are now exported
        so the conversion can be done correctly and so the shape spread can be
        checked as a source of arm-dependent noise in the radial outcome.

    3. LOCAL MAXIMA ARE DEDUPLICATED.
        find_local_maxima used dens == maximum_filter(dens), which returns EVERY
        pixel of a flat maximum. On a plateau that is several adjacent seeds for
        one peak, and watershed then splits one focus into pieces. Connected
        plateau components are now collapsed to a single representative pixel.
        The number collapsed is reported per section, which is a candidate
        explanation for untreated over-segmentation rising at fold 6.

    4. RECORDED AREA MATCHES ASSIGNED AREA.
        area_um2 and n_cells were computed from `region` but the label array was
        written with `region & (labels == 0)`, so when hole filling pushed two
        regions into contact the recorded size of the later focus exceeded what
        it actually got. Geometry is now computed from the pixels the focus
        actually receives, and the number of contested pixels is reported.

    5. THE EXCLUSION BLOCK MATCHES SCRIPT 03.
        The reason strings now quote the composition-controlled rank statistic
        from script 02 revision 2 rather than the composition-sensitive z
        fractions, and the false claim about shortest Y spans is gone.

    6. OPTIONAL BALT ALTERNATIVE GATES (RUN_BALT_ALTERNATIVES).
        The primary BALT definition is unchanged. Two alternative detections are
        run alongside and written to a separate table: B cells with POOLED T
        lineage instead of helper T only, and B cells alone. Helper T cells are
        a CD4 call, and the inventory shows CD4-positive running at 4 to 15
        percent of T cells on scan_01. If the candidate counts hold across all
        three gates the arm comparison is safe; if they do not, the primary gate
        is partly reading CD4 staining. This is a diagnostic. It does not change
        balt_id and nothing downstream reads it yet. It is here so the BALT
        session does not require another full rerun of this script.

        BALT detection depends only on lymphocyte density and each section's own
        CD21 distribution, never on FOCUS_FOLD_OVER_BACKGROUND, so candidate
        counts are invariant to the fold setting. Only min_dist_to_focus_um
        changes when the foci change.

    7. NAMING.
        The script wrote to structures_rev4 while calling itself revision 3 and
        naming its log 00_structures_rev3_report.txt. All three now agree.

OUTPUTS
    figures/  F30 .. F37
    tables/   32 .. 36b
    cell_assignments/<section>_cell_structures.csv

USAGE
    conda activate sc_pre
    python AKOYA_04_Structures.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================
# ALL TUNABLE PARAMETERS LIVE HERE
# =============================================================================

DATA_DIR = "/master/jlehle/WORKING/AKOYA/data"
INVENTORY_TABLE_DIR = "/master/jlehle/WORKING/AKOYA/inventory/tables"
OUT_DIR = "/master/jlehle/WORKING/AKOYA/structures_rev4"

# ---- exclusion, kept in step with script 03 ---------------------------------
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

EXPERT_COUNTS = {
    "G3_43106": 1, "G3_43111": 2, "G3_43118": 2,
    "G4_31438": 7, "G4_36463": 5, "G4_43109": 19,
}
EXPERT_SOURCE = {
    "G3_43106": "expert", "G3_43111": "expert", "G3_43118": "expert",
    "G4_31438": "algorithm_rev1", "G4_36463": "algorithm_rev1",
    "G4_43109": "algorithm_rev1",
}

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
DENSITY_BANDWIDTH_UM = 75.0     # used by BOTH definitions, so they are comparable

# ---- edge correction --------------------------------------------------------
# Normalized convolution. Without it, density within roughly one bandwidth of a
# tissue boundary is damped by however much of the Gaussian kernel falls on
# empty non-tissue space.
EDGE_CORRECTION = True
EDGE_MIN_WEIGHT = 0.25          # floor on the occupancy weight, so a thin
                                # protrusion cannot divide by nearly zero
REPORT_EDGE_COMPARISON = True   # compute the uncorrected map too, for table 32

# =============================================================================
# DEFINITION 1 - ABSOLUTE, for burden
# HARD-CODED ON PURPOSE. Not swept, not tuned. See the docstring.
# =============================================================================
BURDEN_THRESHOLD = 3000.0           # myeloid cells / mm^2, identical everywhere
BURDEN_MIN_AREA_UM2 = 30000.0       # 98 um equivalent radius
BURDEN_MIN_CELLS = 100
BURDEN_FILL_HOLES = True

# =============================================================================
# DEFINITION 2 - RELATIVE, for architecture
# =============================================================================
PEAK_SEPARATION_UM = 200.0          # minimum separation between local maxima
FOCUS_FOLD_OVER_BACKGROUND = 6.0    # peak must exceed this x the section median
FOCUS_HALF_MAX_FRACTION = 0.5       # boundary drawn where density falls to this
                                    # fraction of that focus's own peak
FOCUS_MIN_AREA_UM2 = 10000.0        # 56 um equivalent radius
FOCUS_MIN_CELLS = 40
FOCUS_FILL_HOLES = True             # keeps necrotic acellular centres
CUFF_WIDTH_UM = 150.0
DEDUPLICATE_PLATEAU_PEAKS = True    # collapse flat maxima to one seed each

# background statistic each peak is compared against
BACKGROUND_STAT = "median"          # "median" or "p75"

# ---- sweep for definition 2 -------------------------------------------------
RUN_FOCI_SWEEP = True
SWEEP_FOLDS = [4.0, 6.0, 8.0, 10.0, 15.0]
SWEEP_HALF_FRACTIONS = [0.3, 0.4, 0.5, 0.6]
SWEEP_MIN_AREAS_UM2 = [5000, 10000, 20000, 30000]

# ---- BALT detection ---------------------------------------------------------
# Primary definition unchanged. Detection depends only on lymphocyte density and
# each section's own CD21 distribution, so it is invariant to the focus fold
# setting; only min_dist_to_focus_um moves when the foci move.
BALT_BANDWIDTH_UM = 75.0
BALT_B_THRESHOLD = 800.0
BALT_T_THRESHOLD = 300.0
BALT_MIN_AREA_UM2 = 15000.0
BALT_MIN_B_CELLS = 40
CD21_FOLLICLE_RATIO = 1.5

# Diagnostic only. Does not change balt_id.
RUN_BALT_ALTERNATIVES = True
BALT_ALT_GATES = [
    ("primary_B_and_helperT", [BALT_B], [BALT_T], BALT_B_THRESHOLD, BALT_T_THRESHOLD),
    ("B_and_pooled_T", [BALT_B], T_LINEAGE, BALT_B_THRESHOLD, BALT_T_THRESHOLD),
    ("B_only", [BALT_B], None, BALT_B_THRESHOLD, None),
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
from matplotlib.colors import LinearSegmentedColormap

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
    """
    Kernel density of the selected cells, in cells per mm^2.

    Without edge correction the Gaussian is applied to a zero-padded count
    array, so any pixel within roughly one bandwidth of the tissue boundary
    averages in empty non-tissue space and reads low. With edge correction the
    tissue occupancy is smoothed by the same kernel and used as a denominator,
    which is normalized convolution: inside the tissue the weight is 1 and
    nothing changes, at a boundary it is less than 1 and the estimate is
    restored to what it would be if the tissue continued.
    """
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


def partition_by_markers(dens, markers, territory):
    """
    Assign every pixel of `territory` to one marker.
    Uses skimage watershed on the inverted density when available, which follows
    the density landscape properly. Falls back to a Euclidean nearest-marker
    partition (Voronoi on the seeds) using scipy alone.
    """
    if HAVE_SKIMAGE:
        return _sk_watershed(-dens, markers=markers, mask=territory), "watershed"
    _, inds = ndi.distance_transform_edt(markers == 0, return_indices=True)
    part = markers[inds[0], inds[1]]
    return np.where(territory, part, 0), "edt_voronoi"


def find_local_maxima(dens, occ, sep_px, deduplicate=True):
    """
    Returns (rows, cols, values, n_raw, n_collapsed).

    dens == maximum_filter(dens) marks every pixel of a flat maximum, not one
    per peak. On a plateau that hands watershed several adjacent seeds for a
    single peak and one focus gets split. Connected plateau components are
    collapsed to a single representative, the pixel closest to the component
    centroid.
    """
    mx = ndi.maximum_filter(dens, size=2 * sep_px + 1, mode="constant")
    ispeak = (dens == mx) & occ & (dens > 0)
    n_raw = int(ispeak.sum())
    if not deduplicate or n_raw == 0:
        py, px = np.nonzero(ispeak)
        return py, px, dens[py, px], n_raw, 0
    lab, n = ndi.label(ispeak, structure=np.ones((3, 3), dtype=int))
    ys, xs = [], []
    # centre-of-mass per component, snapped to the nearest member pixel
    coms = ndi.center_of_mass(ispeak, lab, list(range(1, n + 1)))
    objs = ndi.find_objects(lab)
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


def detect_foci(dens, occ, gx, gy, fold_min, half_frac, min_area, min_cells,
                sep_px, fill_holes=True, background_stat="median",
                deduplicate=True):
    """
    DEFINITION 2. Returns (labels, peak_records, method_used, diagnostics).

    Each surviving local maximum becomes one focus whose boundary is drawn where
    density falls to `half_frac` of THAT PEAK's value. Watershed keeps adjacent
    foci separate. Everything is relative to the section's own background.

    Geometry is computed from the pixels a focus ACTUALLY receives, not from the
    candidate region, so table 35 cannot disagree with the label array or with
    the per-cell assignments.
    """
    diag = {"n_peaks_raw": 0, "n_peaks_collapsed": 0, "n_peaks_kept": 0,
            "n_contested_pixels": 0}
    vals = dens[occ]
    if not vals.size:
        return np.zeros_like(dens, dtype=int), [], "none", diag
    bg = float(np.median(vals)) if background_stat == "median" \
        else float(np.percentile(vals, 75))
    if bg <= 0:
        bg = float(np.mean(vals)) or 1.0

    py, px, pv, n_raw, n_coll = find_local_maxima(dens, occ, sep_px, deduplicate)
    diag["n_peaks_raw"] = n_raw
    diag["n_peaks_collapsed"] = n_coll
    keep = pv >= fold_min * bg
    py, px, pv = py[keep], px[keep], pv[keep]
    diag["n_peaks_kept"] = int(len(pv))
    if not len(pv):
        return np.zeros_like(dens, dtype=int), [], "none", diag

    order = np.argsort(pv)[::-1]
    py, px, pv = py[order], px[order], pv[order]

    markers = np.zeros(dens.shape, dtype=int)
    markers[py, px] = np.arange(1, len(pv) + 1)

    # territory limit: the loosest boundary any surviving peak could ask for
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

        # only the pixels not already claimed by a brighter peak
        assigned = region & (labels == 0)
        contested = int(region.sum() - assigned.sum())
        diag["n_contested_pixels"] += contested
        if not assigned[yy, xx]:
            continue
        # keep the piece still connected to the peak after the contest
        cl2, cn2 = ndi.label(assigned)
        if cn2 == 0:
            continue
        assigned = cl2 == cl2[yy, xx]

        area = assigned.sum() * GRID_UM * GRID_UM
        ncell = int(assigned[gy, gx].sum())
        if area < min_area or ncell < min_cells:
            continue

        # geometry of the region as actually assigned
        dt = ndi.distance_transform_edt(assigned)
        r_in = float(dt.max()) * GRID_UM
        r_eq = float(np.sqrt(area / np.pi))

        labels[assigned] = next_id
        recs.append({
            "focus_id": next_id,
            "peak_density": float(vv),
            "background_density": bg,
            "fold_over_background": float(vv / bg) if bg > 0 else np.nan,
            "boundary_density": float(half_frac * vv),
            "peak_row": int(yy), "peak_col": int(xx),
            "area_um2": float(area),
            "equiv_radius_um": r_eq,
            "max_inscribed_radius_um": r_in,
            # 1.0 for a disc, larger for elongated or lobed shapes. The radial
            # coordinate normalises by the inscribed radius, so this is the
            # factor by which a microns conversion using equiv_radius overstates.
            "shape_ratio": float(r_eq / r_in) if r_in > 0 else np.nan,
            "n_contested_pixels": contested,
            "n_cells": ncell,
        })
        next_id += 1
    return labels, recs, method, diag


def detect_burden(dens, occ, gx, gy, threshold, min_area, min_cells,
                  fill_holes=True):
    """DEFINITION 1. Fixed absolute threshold, identical for every section."""
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
    """
    Generic lymphoid aggregate detection. t_set None means the B gate alone.
    Returns (relabelled mask, n_kept).
    """
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

banner("AKOYA STRUCTURE DEFINITION (revision 4, dual definition)")
print(f"Run time : {datetime.now().isoformat(timespec='seconds')}")
print(f"Output   : {OUT_DIR}")
print(f"Detection set (IDO1-blind): {MYELOID_FOR_DETECTION}")
print(f"Grid {GRID_UM} um, bandwidth {DENSITY_BANDWIDTH_UM} um (shared by both "
      f"definitions)")
print(f"Edge correction: {EDGE_CORRECTION} (min weight {EDGE_MIN_WEIGHT})")
if not HAVE_SCIPY:
    print(f"\nERROR: scipy required ({_scipy_err}).")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)
if HAVE_SKIMAGE:
    print("    watershed  : skimage available, using true watershed")
else:
    print(f"    watershed  : skimage NOT available ({_skimage_err})")
    print("                 falling back to Euclidean nearest-marker partition")

sub("DEFINITION 1 - absolute, for burden")
print(f"    threshold {BURDEN_THRESHOLD:.0f} /mm^2, min area "
      f"{BURDEN_MIN_AREA_UM2:,.0f} um^2, min cells {BURDEN_MIN_CELLS}")
print("    Hard-coded on purpose. Not swept.")
sub("DEFINITION 2 - relative, for architecture")
print(f"    peaks >= {FOCUS_FOLD_OVER_BACKGROUND:.0f}x section {BACKGROUND_STAT}, "
      f"boundary at {FOCUS_HALF_MAX_FRACTION:.0%} of each peak,")
print(f"    min area {FOCUS_MIN_AREA_UM2:,.0f} um^2, min cells {FOCUS_MIN_CELLS}, "
      f"peak separation {PEAK_SEPARATION_UM:.0f} um")
print(f"    plateau peak dedup: {DEDUPLICATE_PLATEAU_PEAKS}")

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
    print(f"    inventory loaded, scan and position available for "
          f"{len(SCAN_OF)} section(s)")
else:
    print(f"    WARNING: {fs_path} not found, scan and position unavailable")

csv_paths = sorted(glob.glob(os.path.join(DATA_DIR, "*.csv")))
retained = [p for p in csv_paths
            if os.path.splitext(os.path.basename(p))[0] not in EXCLUDE_SECTIONS]
print(f"    {len(csv_paths)} found, {len(EXCLUDE_SECTIONS)} excluded, "
      f"{len(retained)} retained")

hdr = pd.read_csv(retained[0], nrows=0)
ALL_COLS = hdr.columns.tolist()
XCOL = find_col(ALL_COLS, CENTROID_X_PREFIX)
YCOL = find_col(ALL_COLS, CENTROID_Y_PREFIX)
present_extra = [c for c in EXTRA_MARKER_COLS if c in ALL_COLS]
missing_extra = [c for c in EXTRA_MARKER_COLS if c not in ALL_COLS]
if missing_extra:
    print(f"    WARNING: markers not found, will be NaN: {missing_extra}")

use_cols = [XCOL, YCOL, PHENOTYPE_COL] + present_extra
if IMAGE_COL in ALL_COLS:
    use_cols.append(IMAGE_COL)
dtype_map = ({c: "float32" for c in present_extra + [XCOL, YCOL]}
             if USE_FLOAT32 else None)

cells, geom, dens_of, dens_raw_of = {}, {}, {}, {}
edge_rows = []
for path in retained:
    sid = os.path.splitext(os.path.basename(path))[0]
    cond = GROUP_MAP.get(sid.split("_")[0], "UNKNOWN")
    try:
        df = pd.read_csv(path, usecols=use_cols, dtype=dtype_map, low_memory=False)
    except Exception as e:
        print(f"    ERROR reading {sid}: {e}. Skipping.")
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
    cells[sid] = d

    x, y = d["x"].to_numpy(), d["y"].to_numpy()
    x0, y0, gx, gy, shape = make_grid(x, y, GRID_UM)
    occ = tissue_mask(gx, gy, shape)
    mye = d["pheno"].isin(MYELOID_FOR_DETECTION).to_numpy()

    dens = density_map(gx, gy, shape, mye, DENSITY_BANDWIDTH_UM, GRID_UM,
                       occ=occ, edge_correct=EDGE_CORRECTION)
    dens_raw = (density_map(gx, gy, shape, mye, DENSITY_BANDWIDTH_UM, GRID_UM,
                            occ=occ, edge_correct=False)
                if REPORT_EDGE_COMPARISON else None)

    geom[sid] = {"x0": x0, "y0": y0, "gx": gx, "gy": gy, "shape": shape,
                 "occ": occ, "scan": SCAN_OF.get(sid, "na"),
                 "position": POS_OF.get(sid),
                 "tissue_area_um2": float(occ.sum() * GRID_UM * GRID_UM),
                 "bg_median": float(np.median(dens[occ])),
                 "bg_p75": float(np.percentile(dens[occ], 75)),
                 "peak": float(dens[occ].max())}
    dens_of[sid] = dens
    if dens_raw is not None:
        dens_raw_of[sid] = dens_raw
        # how much of the tissue is within one bandwidth of a boundary?
        w = ndi.gaussian_filter(occ.astype(float),
                                sigma=DENSITY_BANDWIDTH_UM / GRID_UM,
                                mode="constant")
        edge_rows.append({
            "sample_id": sid, "condition": cond,
            "tissue_area_mm2": geom[sid]["tissue_area_um2"] / 1e6,
            "pct_tissue_edge_affected": 100.0 * float((w[occ] < 0.95).mean()),
            "median_edge_weight": float(np.median(w[occ])),
            "bg_median_corrected": geom[sid]["bg_median"],
            "bg_median_uncorrected": float(np.median(dens_raw[occ])),
            "peak_corrected": geom[sid]["peak"],
            "peak_uncorrected": float(dens_raw[occ].max()),
            "pct_above_burden_corrected":
                100.0 * float((dens[occ] >= BURDEN_THRESHOLD).mean()),
            "pct_above_burden_uncorrected":
                100.0 * float((dens_raw[occ] >= BURDEN_THRESHOLD).mean()),
        })

    print(f"    {sid:<12} {cond:<10} {len(d):>9,} cells   "
          f"tissue={geom[sid]['tissue_area_um2']/1e6:5.1f} mm^2   "
          f"bg={geom[sid]['bg_median']:>6.0f}   peak={geom[sid]['peak']:>7.0f}   "
          f"({geom[sid]['peak']/geom[sid]['bg_median']:.1f}x)")
    del df
    gc.collect()

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

# ---- edge correction report -------------------------------------------------
if edge_rows:
    edge = pd.DataFrame(edge_rows)
    sub("Edge correction: what it changed")
    print("    Without correction, density within one bandwidth of a tissue")
    print("    boundary is damped by the fraction of the kernel falling on")
    print("    empty space. This is how much tissue that touches and what it")
    print("    did to the numbers detection depends on.\n")
    print(f"    {'section':<12}{'edge %':>9}{'bg raw':>9}{'bg corr':>9}"
          f"{'peak raw':>10}{'peak corr':>11}{'>thr raw':>10}{'>thr corr':>11}")
    print("    " + "-" * 81)
    for _, r in edge.iterrows():
        print(f"    {r['sample_id']:<12}{r['pct_tissue_edge_affected']:>8.1f}%"
              f"{r['bg_median_uncorrected']:>9.0f}{r['bg_median_corrected']:>9.0f}"
              f"{r['peak_uncorrected']:>10.0f}{r['peak_corrected']:>11.0f}"
              f"{r['pct_above_burden_uncorrected']:>9.2f}%"
              f"{r['pct_above_burden_corrected']:>10.2f}%")
    write_csv(edge, "32_density_edge_correction.csv")


# %% Cell 4 - DEFINITION 1: absolute burden
# =============================================================================

banner("DEFINITION 1 - ABSOLUTE BURDEN")

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
    "background_density": geom[s]["bg_median"],
    "peak_density": geom[s]["peak"],
    "peak_fold_over_background": geom[s]["peak"] / geom[s]["bg_median"]
    if geom[s]["bg_median"] > 0 else np.nan,
} for s in SAMPLE_ORDER])
sub("Burden summary")
print(burden_summary.to_string(index=False))
write_csv(burden_summary, "33b_burden_summary.csv")


# %% Cell 5 - DEFINITION 2: sweep
# =============================================================================

banner("DEFINITION 2 - RELATIVE FOCI, PARAMETER SWEEP")

foci_sweep = pd.DataFrame()
if not RUN_FOCI_SWEEP:
    print("    SKIPPED")
else:
    rows = []
    for s in SAMPLE_ORDER:
        g, dens = geom[s], dens_of[s]
        exp = EXPERT_COUNTS.get(s, np.nan)
        for fold in SWEEP_FOLDS:
            for hf in SWEEP_HALF_FRACTIONS:
                lab, recs, _, _ = detect_foci(
                    dens, g["occ"], g["gx"], g["gy"], fold, hf,
                    min_area=min(SWEEP_MIN_AREAS_UM2), min_cells=FOCUS_MIN_CELLS,
                    sep_px=SEP_PX, fill_holes=FOCUS_FILL_HOLES,
                    background_stat=BACKGROUND_STAT,
                    deduplicate=DEDUPLICATE_PLATEAU_PEAKS)
                rc = pd.DataFrame(recs)
                for ma in SWEEP_MIN_AREAS_UM2:
                    k = rc.loc[rc["area_um2"] >= ma] if len(rc) else rc
                    rows.append({
                        "sample_id": s, "condition": COND_OF[s],
                        "fold_over_background": fold, "half_max_fraction": hf,
                        "min_area_um2": ma, "n_foci": int(len(k)),
                        "expert_count": exp,
                        "abs_error": abs(len(k) - exp) if np.isfinite(exp) else np.nan,
                        "median_radius_um": float(k["equiv_radius_um"].median())
                        if len(k) else np.nan,
                        "pct_tissue_area": (100.0 * float(k["area_um2"].sum())
                                            / g["tissue_area_um2"]) if len(k) else 0.0,
                    })
        print(f"    {s} swept")
    foci_sweep = pd.DataFrame(rows)
    write_csv(foci_sweep, "34_foci_parameter_sweep.csv")

    cal = (foci_sweep.groupby(["fold_over_background", "half_max_fraction",
                               "min_area_um2", "condition"])["abs_error"]
           .sum().unstack("condition").reset_index())
    have_conds = [c for c in CONDITION_ORDER if c in cal.columns]
    cal["total_error"] = cal[have_conds].sum(axis=1)
    cal = cal.sort_values(["total_error"] + have_conds[:1])
    sub("Best 15 combinations by total absolute error")
    print(f"      {'fold':>7}{'half':>7}{'min_area':>11}{'err D1MT':>11}"
          f"{'err Untr':>11}{'total':>8}")
    for _, r in cal.head(15).iterrows():
        print(f"      {r['fold_over_background']:>7.1f}{r['half_max_fraction']:>7.2f}"
              f"{int(r['min_area_um2']):>11}{r.get('D1MT', np.nan):>11.0f}"
              f"{r.get('Untreated', np.nan):>11.0f}{r['total_error']:>8.0f}")
    write_csv(cal, "34b_foci_calibration.csv")

    sub("Per-section detail at the best combination")
    b = cal.iloc[0]
    bs = foci_sweep.loc[(foci_sweep["fold_over_background"] == b["fold_over_background"]) &
                        (foci_sweep["half_max_fraction"] == b["half_max_fraction"]) &
                        (foci_sweep["min_area_um2"] == b["min_area_um2"])]
    print(bs[["sample_id", "condition", "n_foci", "expert_count",
              "median_radius_um", "pct_tissue_area"]].to_string(index=False))
    print("\n    CIRCULARITY REMINDER: the untreated targets in EXPERT_COUNTS")
    print("    came from the revision-1 algorithm, not from a pathologist, so a")
    print("    combination scoring zero error in that arm may simply be")
    print("    reproducing revision 1. Only the three treated counts are expert.")
    print("    The independent check is F35.")
    for k, v in EXPERT_SOURCE.items():
        print(f"      {k:<12} {EXPERT_COUNTS.get(k, 'na'):>3}  ({v})")

    # ---- F32 sweep heatmaps -------------------------------------------------
    n_ma = len(SWEEP_MIN_AREAS_UM2)
    fig, axes = plt.subplots(2, n_ma, figsize=(9 * n_ma, 20), squeeze=False)
    cmap = LinearSegmentedColormap.from_list("sw", ["#FFFFFF", "#9ECAE1", "#08519C"])
    for r, cond in enumerate(CONDITION_ORDER):
        target = sum(v for k, v in EXPERT_COUNTS.items() if COND_OF.get(k) == cond)
        for c, ma in enumerate(SWEEP_MIN_AREAS_UM2):
            ax = axes[r, c]
            m = (foci_sweep.loc[(foci_sweep["condition"] == cond) &
                                (foci_sweep["min_area_um2"] == ma)]
                 .pivot_table(index="fold_over_background",
                              columns="half_max_fraction",
                              values="n_foci", aggfunc="sum"))
            im = ax.imshow(m.to_numpy(float), cmap=cmap, aspect="auto")
            ax.set_xticks(np.arange(m.shape[1]))
            ax.set_xticklabels([f"{v:.1f}" for v in m.columns],
                               fontsize=FONT_SIZE_TICK - 14)
            ax.set_yticks(np.arange(m.shape[0]))
            ax.set_yticklabels([f"{v:.0f}" for v in m.index],
                               fontsize=FONT_SIZE_TICK - 14)
            for i in range(m.shape[0]):
                for j in range(m.shape[1]):
                    v = m.iloc[i, j]
                    if pd.isna(v):
                        continue
                    hit = abs(v - target) <= 1
                    ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                            fontsize=FONT_SIZE_ANNOT - 10,
                            color=FLAG_COLOR if hit else "#000000",
                            fontweight="bold" if hit else "normal")
            ax.set_xlabel("half-max fraction", fontsize=FONT_SIZE_BASE - 12)
            if c == 0:
                ax.set_ylabel("fold over background", fontsize=FONT_SIZE_BASE - 12)
            ax.set_title(f"{cond}, min area {int(ma):,} µm²\ntarget {target}",
                         fontsize=FONT_SIZE_TITLE - 14)
    fig.suptitle("Relative foci sweep: total foci across sections in each arm\n"
                 "Red bold = within 1 of the target "
                 "(untreated targets are algorithmic, not expert)", y=1.02,
                 fontsize=FONT_SIZE_TITLE - 6)
    save_fig(fig, "F32_foci_parameter_sweep")


# %% Cell 6 - DEFINITION 2: final detection, BALT, cell assignment
# =============================================================================

banner("FINAL DETECTION AND CELL ASSIGNMENT")

foci_rows, balt_rows, balt_alt_rows, masks = [], [], [], {}
method_used = "none"
diag_rows = []

for s in SAMPLE_ORDER:
    d = cells[s]
    g, dens = geom[s], dens_of[s]
    gx, gy, shape, occ = g["gx"], g["gy"], g["shape"], g["occ"]
    x0, y0 = g["x0"], g["y0"]
    pheno = d["pheno"].to_numpy()

    labels, recs, method, diag = detect_foci(
        dens, occ, gx, gy, FOCUS_FOLD_OVER_BACKGROUND, FOCUS_HALF_MAX_FRACTION,
        FOCUS_MIN_AREA_UM2, FOCUS_MIN_CELLS, SEP_PX,
        fill_holes=FOCUS_FILL_HOLES, background_stat=BACKGROUND_STAT,
        deduplicate=DEDUPLICATE_PLATEAU_PEAKS)
    if method != "none":
        method_used = method
    n_foci = len(recs)
    diag.update({"sample_id": s, "condition": COND_OF[s], "n_foci": n_foci})
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

    radial = np.full(shape, np.nan, dtype=float)
    for i in range(1, n_foci + 1):
        m = labels == i
        dt = ndi.distance_transform_edt(m)
        mx = dt.max()
        if mx > 0:
            radial[m] = 1.0 - (dt[m] / mx)
    if core_mask.any():
        radial[cuff_mask] = 1.0 + np.clip(dout[cuff_mask] / CUFF_WIDTH_UM, 0, 1)

    # ---- BALT, primary definition ------------------------------------------
    brelab, n_balt = detect_lymphoid(
        gx, gy, shape, occ, pheno, [BALT_B], [BALT_T],
        BALT_B_THRESHOLD, BALT_T_THRESHOLD, BALT_BANDWIDTH_UM,
        BALT_MIN_AREA_UM2, BALT_MIN_B_CELLS)

    masks[s] = {"x0": x0, "y0": y0, "shape": shape, "core": labels,
                "cuff": cuff_mask, "balt": brelab, "occ": occ, "dens": dens,
                "burden": burden_labels[s]}

    # ---- per-cell assignment -----------------------------------------------
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
    out["section_background_density"] = g["bg_median"]
    if WRITE_CELL_ASSIGNMENTS:
        write_csv(out, f"{s}_cell_structures.csv", directory=CELL_DIR)

    # ---- per-focus metrics -------------------------------------------------
    for r in recs:
        i = r["focus_id"]
        sel_core = cell_core == i
        sel_cuff = cell_cuff == i
        core_cells = out.loc[sel_core]
        all_cells = out.loc[sel_core | sel_cuff]
        mac = core_cells.loc[core_cells["pheno"].isin(CD68_LINEAGE)]
        ys, xs = np.nonzero(labels == i)
        rec = dict(r)
        rec.update({
            "sample_id": s, "condition": COND_OF[s], "scan_id": g["scan"],
            "centroid_x_um": float(x0 + (xs.mean() - 1) * GRID_UM),
            "centroid_y_um": float(y0 + (ys.mean() - 1) * GRID_UM),
            "pct_of_tissue_area": 100.0 * r["area_um2"] / g["tissue_area_um2"],
            "n_cells_core": int(sel_core.sum()),
            "n_cells_cuff": int(sel_cuff.sum()),
            # consistency guard: n_cells comes from the label array, n_cells_core
            # from the per-cell assignment. They must agree.
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
        })
        for ph in PHENOTYPE_ORDER:
            n_ph = int((all_cells["pheno"] == ph).sum())
            rec[f"n_{ascii_safe(ph)}"] = n_ph
            rec[f"pct_{ascii_safe(ph)}"] = (100.0 * n_ph / len(all_cells)
                                            if len(all_cells) else np.nan)
        foci_rows.append(rec)

    # ---- BALT metrics ------------------------------------------------------
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
            # a candidate sitting on a granuloma's lymphocyte cuff is a
            # different object from a free-standing follicle in uninvolved
            # parenchyma. Q7 needs to tell them apart.
            "pct_in_focus_core": 100.0 * float((s_all["region"] == "core").mean())
            if len(s_all) else np.nan,
            "pct_in_focus_cuff": 100.0 * float((s_all["region"] == "cuff").mean())
            if len(s_all) else np.nan,
            "pct_in_burden_region": 100.0 * float(s_all["in_burden_region"].mean())
            if len(s_all) else np.nan,
        })

    # ---- BALT alternative gates, diagnostic only ---------------------------
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
    print(f"    {s:<12} foci={n_foci:>3}  BALT={n_balt:>3}  "
          f"bg={g['bg_median']:>6.0f}   peaks {diag['n_peaks_raw']}->"
          f"{diag['n_peaks_raw'] - diag['n_peaks_collapsed']} after dedup, "
          f"{diag['n_peaks_kept']} above fold, {diag['n_contested_pixels']} "
          f"contested px{flag}")

foci = pd.DataFrame(foci_rows)
balts = pd.DataFrame(balt_rows)
balt_alt = pd.DataFrame(balt_alt_rows)
diagnostics = pd.DataFrame(diag_rows)
print(f"\n    Partition method used: {method_used}")
write_csv(diagnostics, "35b_detection_diagnostics.csv")

if len(foci):
    bad = foci.loc[~foci["geometry_matches_assignment"]]
    if len(bad):
        print(f"\n    WARNING: {len(bad)} focus record(s) where the label-array "
              f"cell count disagrees with the per-cell assignment. Investigate "
              f"before using table 35 geometry.")
        print(bad[["sample_id", "focus_id", "n_cells", "n_cells_core"]].to_string(index=False))
    else:
        print("    Geometry check: every focus record matches its per-cell "
              "assignment.")

sub("Foci per section, with the scale they operate at")
for s in SAMPLE_ORDER:
    f = foci.loc[foci["sample_id"] == s] if len(foci) else pd.DataFrame()
    b = balts.loc[balts["sample_id"] == s] if len(balts) else pd.DataFrame()
    nf = int(b["is_follicle"].sum()) if len(b) else 0
    if len(f):
        print(f"    {s:<12} {COND_OF[s]:<10} foci={len(f):>3} "
              f"(expert {EXPERT_COUNTS.get(s, 'na')})  "
              f"median radius={f['equiv_radius_um'].median():>5.0f} um  "
              f"inscribed={f['max_inscribed_radius_um'].median():>5.0f} um  "
              f"shape={f['shape_ratio'].median():>4.2f}  "
              f"peak={f['peak_density'].median():>7.0f}  "
              f"fold={f['fold_over_background'].median():>5.1f}x  "
              f"BALT={len(b):>3} ({nf} CD21+)")
    else:
        print(f"    {s:<12} {COND_OF[s]:<10} foci=  0 "
              f"(expert {EXPERT_COUNTS.get(s, 'na')})  "
              f"BALT={len(b):>3} ({nf} CD21+)")

if len(foci):
    write_csv(foci, "35_foci_structures_relative.csv")
if len(balts):
    write_csv(balts, "36_balt_structures.csv")
if len(balt_alt):
    write_csv(balt_alt, "36b_balt_alternative_gates.csv")
    sub("BALT gate sensitivity (diagnostic, does not change balt_id)")
    print("    Helper T cells are a CD4 call and CD4-positive runs at 4 to 15")
    print("    percent of T cells on scan_01. If candidate counts hold across")
    print("    all three gates the arm comparison is safe. If they do not, the")
    print("    primary gate is partly reading CD4 staining rather than lymphoid")
    print("    organisation.\n")
    piv = balt_alt.pivot_table(index="sample_id", columns="gate",
                               values="n_candidates")
    pivf = balt_alt.pivot_table(index="sample_id", columns="gate",
                                values="n_cd21_follicles")
    gates = [g[0] for g in BALT_ALT_GATES if g[0] in piv.columns]
    print(f"    {'section':<14}" + "".join(f"{g[:22]:>24}" for g in gates))
    print("    " + "-" * (14 + 24 * len(gates)))
    for s in SAMPLE_ORDER:
        if s not in piv.index:
            continue
        row = f"    {s:<14}"
        for g in gates:
            row += f"{int(piv.loc[s, g]):>13} ({int(pivf.loc[s, g])} CD21+)"[:24].rjust(24)
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
    print("\n    Treated foci are defined relative to treated background. They are")
    print("    residual myeloid foci, NOT granulomas equivalent to the untreated")
    print("    lesions. Report peak density and fold alongside any architecture")
    print("    comparison so the two are never conflated.")
    print("\n    USE max_inscribed_radius_um, NOT equiv_radius_um, to convert a")
    print("    radial-unit effect into microns. The radial coordinate is")
    print("    normalised by the inscribed radius. shape_ratio is the factor by")
    print("    which the equivalent radius would overstate it.")


# %% Cell 7 - figures
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


# ---- F30 burden regions -----------------------------------------------------
def draw_burden(ax, s):
    d = base_panel(s)
    other = ~d["pheno"].isin(MYELOID_FOR_DETECTION)
    ax.scatter(d.loc[other, "x"], d.loc[other, "y"], s=0.6,
               color=BACKGROUND_COLOR, linewidths=0, rasterized=True)
    for ph in MYELOID_FOR_DETECTION:
        m = d["pheno"] == ph
        if m.sum():
            ax.scatter(d.loc[m, "x"], d.loc[m, "y"], s=5,
                       color=PHENOTYPE_COLORS[ph], linewidths=0, rasterized=True)
    contour(ax, s, masks[s]["burden"] > 0, BURDEN_OUTLINE, 4)
    n = int((burden["sample_id"] == s).sum()) if len(burden) else 0
    pct = (float(burden.loc[burden["sample_id"] == s, "pct_of_tissue_area"].sum())
           if len(burden) else 0.0)
    ax.set_title(f"{short_label(s)} ({COND_OF[s]})   {n} regions, {pct:.1f}% burden",
                 fontsize=FONT_SIZE_TITLE - 12)


spatial_grid("F30_burden_absolute_overlay",
             f"DEFINITION 1: absolute threshold {BURDEN_THRESHOLD:.0f} myeloid/mm² "
             "applied identically to every section\n"
             "The near-absence in the treated arm is the result",
             draw_burden,
             [Patch(facecolor=PHENOTYPE_COLORS[p], label=p) for p in MYELOID_FOR_DETECTION]
             + [Line2D([0], [0], color=BURDEN_OUTLINE, linewidth=4,
                       label="granuloma-density region")])

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
ax.set_ylabel("Regions"); ax.set_title("Granuloma-density regions",
                                       fontsize=FONT_SIZE_TITLE - 10)
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
    ax.scatter([i], [geom[s]["bg_median"]], s=420, color="#999999",
               marker="_", linewidth=6, zorder=3)
    ax.scatter([i], [geom[s]["peak"]], s=420, color=COLOR_OF[s],
               marker=MARKER_OF[s], edgecolor="#FFFFFF", linewidth=2, zorder=4)
    ax.plot([i, i], [geom[s]["bg_median"], geom[s]["peak"]], color=COLOR_OF[s],
            linewidth=4, alpha=0.6, zorder=2)
ax.axhline(BURDEN_THRESHOLD, color=BURDEN_OUTLINE, linestyle="--", linewidth=3.5)
ax.text(len(SAMPLE_ORDER) - 0.5, BURDEN_THRESHOLD * 1.1,
        f"burden threshold {BURDEN_THRESHOLD:.0f}", ha="right",
        fontsize=FONT_SIZE_ANNOT - 10, color=BURDEN_OUTLINE)
ax.set_yscale("log")
ax.set_xticks(xp); ax.set_xticklabels(xl, rotation=45, ha="right",
                                      fontsize=FONT_SIZE_TICK - 8)
ax.set_ylabel("Myeloid density (cells / mm$^2$)")
ax.set_title("Background to peak range", fontsize=FONT_SIZE_TITLE - 12)
style_axes(ax)
fig.suptitle("Absolute burden, and why one threshold cannot serve both arms",
             y=1.03, fontsize=FONT_SIZE_TITLE - 4)
save_fig(fig, "F31_burden_summary")


# ---- F33 foci overlay -------------------------------------------------------
def draw_foci(ax, s):
    d = base_panel(s)
    other = ~d["pheno"].isin(MYELOID_FOR_DETECTION)
    ax.scatter(d.loc[other, "x"], d.loc[other, "y"], s=0.6,
               color=BACKGROUND_COLOR, linewidths=0, rasterized=True)
    for ph in MYELOID_FOR_DETECTION:
        m = d["pheno"] == ph
        if m.sum():
            ax.scatter(d.loc[m, "x"], d.loc[m, "y"], s=5,
                       color=PHENOTYPE_COLORS[ph], linewidths=0, rasterized=True)
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
             "DEFINITION 2: foci detected relative to each section's own "
             "background\nLabels give fold over that section's background density",
             draw_foci,
             [Patch(facecolor=PHENOTYPE_COLORS[p], label=p) for p in MYELOID_FOR_DETECTION]
             + [Line2D([0], [0], color=FOCUS_OUTLINE, linewidth=4, label="focus core"),
                Line2D([0], [0], color=CUFF_OUTLINE, linewidth=2.5,
                       label=f"cuff (+{CUFF_WIDTH_UM:.0f} µm)")])

# ---- F34 the scale figure ---------------------------------------------------
if len(foci):
    fig, axes = plt.subplots(1, 3, figsize=(42, 13))
    ax = axes[0]
    for i, c in enumerate(CONDITION_ORDER):
        f = foci.loc[foci["condition"] == c]
        if not len(f):
            continue
        j = rng.uniform(-0.14, 0.14, size=len(f))
        ax.scatter(np.full(len(f), i) + j, f["peak_density"], s=320,
                   color=CONDITION_COLORS[c], edgecolor="#FFFFFF", linewidth=2,
                   zorder=3)
        ax.hlines(f["peak_density"].median(), i - 0.3, i + 0.3, color="#000000",
                  linewidth=4)
    ax.set_xticks(range(len(CONDITION_ORDER)))
    ax.set_xticklabels(CONDITION_ORDER)
    ax.set_yscale("log"); ax.set_ylabel("Peak myeloid density (cells / mm$^2$)")
    ax.set_title("Absolute scale of each focus", fontsize=FONT_SIZE_TITLE - 12)
    style_axes(ax)

    ax = axes[1]
    for i, c in enumerate(CONDITION_ORDER):
        f = foci.loc[foci["condition"] == c]
        if not len(f):
            continue
        j = rng.uniform(-0.14, 0.14, size=len(f))
        ax.scatter(np.full(len(f), i) + j, f["fold_over_background"], s=320,
                   color=CONDITION_COLORS[c], edgecolor="#FFFFFF", linewidth=2,
                   zorder=3)
        ax.hlines(f["fold_over_background"].median(), i - 0.3, i + 0.3,
                  color="#000000", linewidth=4)
    ax.set_xticks(range(len(CONDITION_ORDER)))
    ax.set_xticklabels(CONDITION_ORDER)
    ax.set_ylabel("Fold over own section background")
    ax.set_title("Relative prominence", fontsize=FONT_SIZE_TITLE - 12)
    style_axes(ax)

    # shape, which is what the radial coordinate is normalised by
    ax = axes[2]
    for i, c in enumerate(CONDITION_ORDER):
        f = foci.loc[foci["condition"] == c]
        if not len(f):
            continue
        j = rng.uniform(-0.14, 0.14, size=len(f))
        ax.scatter(np.full(len(f), i) + j, f["shape_ratio"], s=320,
                   color=CONDITION_COLORS[c], edgecolor="#FFFFFF", linewidth=2,
                   zorder=3)
        ax.hlines(f["shape_ratio"].median(), i - 0.3, i + 0.3, color="#000000",
                  linewidth=4)
    ax.axhline(1.0, color="#000000", linestyle="--", linewidth=2.5)
    ax.set_xticks(range(len(CONDITION_ORDER)))
    ax.set_xticklabels(CONDITION_ORDER)
    ax.set_ylabel("Equivalent radius / inscribed radius")
    ax.set_title("Shape (1.0 = disc)", fontsize=FONT_SIZE_TITLE - 12)
    style_axes(ax)
    fig.suptitle("Treated foci are prominent within their own tissue but operate "
                 "at a lower absolute scale\nRight panel: the radial coordinate "
                 "normalises by the inscribed radius, so shape spread is noise "
                 "in the radial outcome", y=1.05,
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
    contour(ax, s, masks[s]["burden"] > 0, BURDEN_OUTLINE, 3, ls="dashed")
    full = cells[s]
    ins = full.loc[full["pheno"] == IDO1_POS_PHENO]
    pct_f = pct_b = np.nan
    if len(ins):
        mk = masks[s]
        ggx = np.clip(((ins["x"].to_numpy() - mk["x0"]) / GRID_UM).astype(int) + 1,
                      0, mk["shape"][1] - 1)
        ggy = np.clip(((ins["y"].to_numpy() - mk["y0"]) / GRID_UM).astype(int) + 1,
                      0, mk["shape"][0] - 1)
        pct_f = 100.0 * float((mk["core"][ggy, ggx] > 0).mean())
        pct_b = 100.0 * float((mk["burden"][ggy, ggx] > 0).mean())
    ax.set_title(f"{short_label(s)} ({COND_OF[s]})\n"
                 f"IDO1+ inside foci {pct_f:.0f}%, inside burden {pct_b:.0f}%",
                 fontsize=FONT_SIZE_TITLE - 14)


spatial_grid("F35_ido1_capture_validation",
             "Independent validation: IDO1+ macrophages against boundaries drawn "
             "without IDO1\nSolid = relative foci, dashed = absolute burden regions",
             draw_capture,
             [Patch(facecolor=PHENOTYPE_COLORS[IDO1_POS_PHENO], label=IDO1_POS_PHENO),
              Line2D([0], [0], color=FOCUS_OUTLINE, linewidth=4, label="focus"),
              Line2D([0], [0], color=BURDEN_OUTLINE, linewidth=3, linestyle="--",
                     label="burden region")])

# ---- F36 BALT ---------------------------------------------------------------
show = [BALT_B, BALT_T, PLASMA]


def draw_balt(ax, s):
    d = base_panel(s)
    other = ~d["pheno"].isin(show)
    ax.scatter(d.loc[other, "x"], d.loc[other, "y"], s=0.6,
               color=BACKGROUND_COLOR, linewidths=0, rasterized=True)
    for ph, sz in zip(show, [12, 14, 6]):
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
             [Patch(facecolor=PHENOTYPE_COLORS[p], label=p) for p in show]
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
    ax.set_xlabel("Distance to nearest myeloid focus (µm)")
    ax.set_ylabel("% plasma cells in candidate")
    ax.set_title("Proximity to foci and plasma content",
                 fontsize=FONT_SIZE_TITLE - 14)
    style_axes(ax)

    # is the candidate free-standing, or sitting on a granuloma?
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


# %% Cell 8 - wrap up
# =============================================================================

banner("SUMMARY")
print(f"Partition method     : {method_used}")
print(f"Edge correction      : {EDGE_CORRECTION}")
print(f"Burden regions (abs) : {len(burden)}")
print(f"Foci (relative)      : {len(foci)}")
print(f"BALT candidates      : {len(balts)}"
      f"  ({int(balts['is_follicle'].sum()) if len(balts) else 0} CD21+)")
for c in CONDITION_ORDER:
    bsel = burden_summary.loc[burden_summary["condition"] == c]
    fsel = foci.loc[foci["condition"] == c] if len(foci) else pd.DataFrame()
    print(f"  {c:<12} burden {bsel['burden_pct'].mean():>5.1f}% mean, "
          f"{len(fsel):>3} foci")

sub("Read in this order")
print("  1. table 32        : what the edge correction changed. Read this")
print("     first, because every number below moved with it.")
print("  2. F31 right panel : background-to-peak range per section, with the")
print("     absolute threshold drawn. This is the whole argument for two")
print("     definitions in one picture.")
print("  3. F30 / table 33b : burden under the absolute definition.")
print("  4. F32 / table 34b : foci sweep. Remember the untreated targets are")
print("     algorithmic, not expert.")
print("  5. F33             : do the relative boundaries match the slides?")
print("  6. F34             : the scale figure, now with the shape panel.")
print("  7. F35             : IDO1+ capture, the independent validation.")
print("  8. table 36b       : BALT gate sensitivity, for the Q7 session.")

sub("Downstream must be rerun")
print("  Every script from 05 onward reads cell_assignments/ or table 35.")
print("  The edge correction moves focus boundaries, which moves radial_pos,")
print("  which moves the Finding 3 coefficient. Treat every number carried in")
print("  the project notes as provisional until 06, 07, 08 and 09 have rerun.")

banner("DONE")
sys.stdout = _tee.terminal
_tee.close()
