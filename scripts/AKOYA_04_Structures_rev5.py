#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - STRUCTURE DETECTION, REVISION 5
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 04 of the AKOYA analysis series. WRITES THE CENTRAL ARTEFACT.

This supersedes revision 4. Every downstream script reads
structures_rev5/cell_assignments/, so 06, 07, 08 and 09 must be re-run after
this and every number in the manuscript re-baselined.

WHAT CHANGED FROM REVISION 4, AND WHAT MEASUREMENT DROVE IT

  Parameters were calibrated against 83 expert annotations across five
  sections, rasterised onto this script's own 25 um grid. They were not tuned
  to any between-arm contrast, and the same values apply to both arms.

  1. BACKGROUND_STAT  median -> p25 of tissue-pixel density
     Section backgrounds under the median ran 211 to 1,686 cells/mm2, so a 6x
     gate meant 1,269 cells/mm2 in one section and 10,117 in another. Annotated
     lesions on the most involved section reached only 4.9x its median
     background, BELOW the gate, while every other section reached 7.7 to 19.8x.
     The median of a heavily involved section is itself pathology, so a
     self-referencing gate is hardest to clear in the most diseased tissue.
     Under p25 all five annotated sections clear 6x.

  2. FOCUS_HALF_MAX_FRACTION  0.50 -> 0.25
     Annotated boundaries sit at 0.224 of the enclosed peak for complexes and
     0.258 for sharp-bordered lesions. Corroborated independently: revision 4
     structures had a median area 0.31x the regions drawn over them.

  3. NEW, granuloma complexes by level set
     The most involved section holds 13 local maxima in its entire area, and 15
     of 24 annotated regions there contain no local maximum at all. No threshold
     on a peak list recovers a region that never generated a peak, so a
     peak-free definition was added rather than the existing one retuned.

  4. NEW, MIN_PEAK_DENSITY 1750 cells/mm2
     Window fixed at 1,565 to 1,925 by the annotations: above a structure the
     expert called vascular, below one overlapping a drawn lesion. Of 18 objects
     it excludes, NONE contain IDO1+ macrophages; of 57 it retains, 44 do
     (Fisher exact p = 2.2e-09). IDO1 is used nowhere in detection, so that
     partition is a consequence of the density threshold and not an input to it.

  5. PEAK_SEPARATION_UM  unchanged at 200
     Measurement suggested 772, which would set the maximum-filter footprint to
     1,575 um and merge distinct foci in sections where detection was already
     acceptable. The measurement correctly diagnoses confluence; the level-set
     definition addresses confluence without that cost.

THREE DEFINITIONS, NEVER POOLED

  burden_*    ABSOLUTE, one fixed threshold applied identically to every
              section. UNCHANGED from revision 4. This is the yardstick that
              produces the zero in the treated animals and it must not become
              section-relative.
  focus_*     PEAK-BASED, for internal architecture. Supplies the centre and
              inscribed radius the radial coordinate needs. Scripts 06, 07 and
              08 read this.
  complex_*   LEVEL-SET, for counting and extent. A confluent mass is one
              structure because it is one connected component.

  focus_id is retained with revision 4 semantics so downstream scripts continue
  to work unchanged. complex_id is new and additive.

CIRCULARITY, AND THE SWITCH THAT ADDRESSES IT

  Structures are defined by POOLED MYELOID DENSITY. That constrains one
  quantity. It does not constrain composition within the pool, where any cell
  type sits, or any marker intensity. So:

    - Lymphocyte results (B cells, plasma cells, T cells) are entirely free:
      those phenotypes are not in the pool.
    - Radial position of a pooled member is nearly free. The only coupling is
      that the boundary sits where density falls to a fraction of peak, and a
      rim-ward population extends the boundary it is then measured against.
      That compresses its own normalised radial position, so the coupling is
      CONSERVATIVE for a rim-ward finding.
    - ABUNDANCE of a pooled member inside structures is genuinely exposed: a
      structure can reach threshold by being neutrophil-rich or macrophage-rich.

  MYELOID_POOL is therefore a parameter. Setting LEAVE_ONE_OUT to a phenotype
  re-runs detection with it excluded, writing to a separate output directory, so
  a population can be tested inside structures defined without it. Run the
  primary first, then the leave-one-out as a sensitivity analysis. Do not mix
  their outputs.

INPUTS
    AKOYA/data/*.csv                     (6 retained sections)
    inventory/tables/01_file_summary.csv (scan and position, optional)

OUTPUTS
    structures_rev5/cell_assignments/<section>_cell_structures.csv
    structures_rev5/tables/   32 to 36b, 37 (new, complexes)
    structures_rev5/figures/  F30 to F38

USAGE
    conda activate sc_pre
    python AKOYA_04_Structures_rev5.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================

DATA_DIR = "/master/jlehle/WORKING/AKOYA/data"
INVENTORY = "/master/jlehle/WORKING/AKOYA/inventory/tables/01_file_summary.csv"
OUT_DIR = "/master/jlehle/WORKING/AKOYA/structures_rev5"

USE_AGG = True

EXCLUDE_SECTIONS = ["G3_43102", "G4_43112"]   # position 1, dim, excluded in 03

# ---- density estimation, unchanged from revision 4 --------------------------
GRID_UM = 25.0
DENSITY_BANDWIDTH_UM = 75.0
EDGE_CORRECTION = True
EDGE_MIN_WEIGHT = 0.25
TISSUE_CLOSE_ITERATIONS = 2
TISSUE_FILL_HOLES = True

# ---- background, CHANGED ----------------------------------------------------
BACKGROUND_STAT = "p25"            # was "median"

# ---- ABSOLUTE burden, UNCHANGED. Do not make this section-relative. ---------
BURDEN_DENSITY = 3000.0
BURDEN_MIN_AREA_UM2 = 30000.0
BURDEN_MIN_CELLS = 100

# ---- foci, peak-based, for architecture -------------------------------------
FOCUS_FOLD_OVER_BACKGROUND = 6.0
PEAK_SEPARATION_UM = 200.0
FOCUS_HALF_MAX_FRACTION = 0.25     # was 0.50
FOCUS_MIN_AREA_UM2 = 10000.0
FOCUS_MIN_CELLS = 40
FOCUS_FILL_HOLES = True
DEDUPLICATE_PLATEAU_PEAKS = True
CUFF_DILATION_UM = 150.0

# ---- complexes, level-set, for counting and extent, NEW ---------------------
COMPLEX_LEVEL_FOLD = 6.5
COMPLEX_CLOSE_UM = 0.0
COMPLEX_MIN_AREA_UM2 = 10000.0
COMPLEX_MIN_CELLS = 40
MIN_PEAK_DENSITY = 1750.0          # NEW, applies to complexes and foci alike

# ---- the myeloid pool, parameterised for the leave-one-out ------------------
IDO1_POS = "CD68+IDO1+ Macrophages"
IDO1_NEG = "CD68+IDO1- Macrophages"
MYELOID_POOL = [IDO1_POS, IDO1_NEG, "CD163+ Macrophages", "Neutrophils"]

# Set to a phenotype name to exclude it from the pool and write to a separate
# directory, e.g. "Neutrophils" to test neutrophil abundance inside structures
# defined without them. None runs the primary analysis.
LEAVE_ONE_OUT = None

# ---- phenotypes -------------------------------------------------------------
PHENOTYPE_ORDER = [
    IDO1_POS, IDO1_NEG, "CD163+ Macrophages", "Neutrophils",
    "Helper T cells", "CD4- T cells", "Tregs", "B cells", "Plasma cells",
    "Endothelial cells", "Epithelial/Tumor cells", "Other",
]
LYMPHOID = ["Helper T cells", "CD4- T cells", "Tregs", "B cells", "Plasma cells"]

# Column names in the QuPath export, matching the conventions established in
# scripts 01, 02 and 03. Centroid columns are matched by PREFIX so the micron
# symbol is never hard-coded, and the phenotype column is 'Phenotypes', plural.
CENTROID_X_PREFIX = "Centroid X"
CENTROID_Y_PREFIX = "Centroid Y"
PHENOTYPE_COL = "Phenotypes"
PARENT_COL = "Parent"
IMAGE_COL = "Image"

# Script 01 found a zero-width space embedded in the epithelial phenotype label.
# Any label carrying one silently fails string matching, so every phenotype
# string is normalised before use. Without this a mismatch would present as an
# empty myeloid pool rather than as a naming problem.
NORMALISE_PHENOTYPE_LABELS = True

CONDITION_FROM_PREFIX = {"G3": "D1MT", "G4": "Untreated"}

# ---- appearance -------------------------------------------------------------
CONDITION_ORDER = ["D1MT", "Untreated"]
CONDITION_COLORS = {"D1MT": "#2C7FB8", "Untreated": "#D95F02"}
CONDITION_LABELS = {"D1MT": "D1MT-treated", "Untreated": "Untreated"}
PHENOTYPE_COLORS = {
    IDO1_POS: "#B2182B", IDO1_NEG: "#EF8A62", "CD163+ Macrophages": "#FDBE85",
    "Neutrophils": "#7B3294", "Helper T cells": "#1B7837",
    "CD4- T cells": "#7FBC41", "Tregs": "#00441B", "B cells": "#2166AC",
    "Plasma cells": "#67A9CF", "Endothelial cells": "#E7298A",
    "Epithelial/Tumor cells": "#66C2A5", "Other": "#A6761D",
}
BURDEN_COLOR = "#000000"
FOCUS_COLOR = "#00A0C6"
COMPLEX_COLOR = "#7B3294"
CUFF_COLOR = "#999999"
GRID_COLOR = "#DDDDDD"
AXIS_COLOR = "#333333"
TEXT_COLOR = "#000000"

FONT_SIZE_BASE = 28
FONT_SIZE_TITLE = 32
FONT_SIZE_TICK = 28
FONT_SIZE_LEGEND = 28
FONT_SIZE_PANEL = 34
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
import warnings
import unicodedata
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

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    try:
        pd.set_option("future.no_silent_downcasting", True)
    except Exception:
        pass

if LEAVE_ONE_OUT:
    OUT_DIR = f"{OUT_DIR}_no_{LEAVE_ONE_OUT.replace(' ', '_').replace('+', '')}"
    POOL = [p for p in MYELOID_POOL if p != LEAVE_ONE_OUT]
    if len(POOL) == len(MYELOID_POOL):
        raise SystemExit(f"LEAVE_ONE_OUT='{LEAVE_ONE_OUT}' is not in "
                         f"MYELOID_POOL: {MYELOID_POOL}")
else:
    POOL = list(MYELOID_POOL)

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
SEP_PX = max(1, int(round(PEAK_SEPARATION_UM / GRID_UM)))
PX_AREA_UM2 = GRID_UM * GRID_UM


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


def short_label(s):
    return s.split("_")[-1] if "_" in s else s


def save_fig(fig, stem):
    for flag, ext in ((SAVE_PDF, "pdf"), (SAVE_PNG, "png")):
        if not flag:
            continue
        p = os.path.join(FIG_DIR, f"{stem}.{ext}")
        fig.savefig(p, dpi=DPI, bbox_inches="tight")
        print(f"    wrote {p}")
    plt.close(fig); gc.collect()


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


def panel_letter(ax, letter, dx=-0.16, dy=1.07):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=FONT_SIZE_PANEL,
            fontweight="bold", va="top", ha="left")


def find_col(df, prefix):
    """Match by prefix, as scripts 01 to 03 do. Warns on ambiguity."""
    hits = [c for c in df.columns if str(c).startswith(prefix)]
    if not hits:
        hits = [c for c in df.columns
                if str(c).strip().lower().startswith(prefix.strip().lower())]
    if not hits:
        return None
    if len(hits) > 1:
        print(f"    WARNING: {len(hits)} columns start with '{prefix}': "
              f"{hits[:4]}. Using the first.")
    return hits[0]


def ascii_safe(v):
    """
    Normalise a phenotype label. Script 01 identified a zero-width space in the
    epithelial label, which breaks equality matching without any visible sign.
    """
    if not isinstance(v, str):
        return v
    v = unicodedata.normalize("NFKC", v)
    return v.replace("\u200b", "").replace("\u00a0", " ").strip()


_tee = Tee(os.path.join(TAB_DIR, "00_structures_report.txt"))
sys.stdout = _tee

banner("AKOYA STRUCTURE DETECTION, REVISION 5")
print(f"Run time : {datetime.now().isoformat(timespec='seconds')}")
print(f"Output   : {OUT_DIR}")
print("\nThis WRITES the cell assignment files every downstream script reads.")
print("Scripts 06, 07, 08 and 09 must be re-run after this, and every")
print("manuscript number re-baselined from that run.")

sub("Settings")
print(f"    grid / bandwidth        {GRID_UM} um / {DENSITY_BANDWIDTH_UM} um")
print(f"    background statistic    {BACKGROUND_STAT}          (rev4: median)")
print(f"    burden threshold        {BURDEN_DENSITY:,.0f} /mm2      UNCHANGED")
print(f"    focus fold / growth     {FOCUS_FOLD_OVER_BACKGROUND:g}x / "
      f"{FOCUS_HALF_MAX_FRACTION:.2f}        (rev4 growth: 0.50)")
print(f"    peak separation         {PEAK_SEPARATION_UM:.0f} um       UNCHANGED")
print(f"    complex level           {COMPLEX_LEVEL_FOLD:g}x           NEW")
print(f"    minimum peak density    {MIN_PEAK_DENSITY:,.0f} /mm2      NEW")
print(f"    myeloid pool            {', '.join(POOL)}")
if LEAVE_ONE_OUT:
    print(f"\n    LEAVE-ONE-OUT RUN: '{LEAVE_ONE_OUT}' excluded from the pool.")
    print("    This is a SENSITIVITY analysis. Its output lives in its own")
    print("    directory and must not be mixed with the primary run. Use it to")
    print("    test a population inside structures defined without it.")
if not HAVE_SKIMAGE:
    print("\n    WARNING: skimage missing. The EDT fallback will not reproduce")
    print("    the watershed partition. Install scikit-image.")


# %% Cell 3 - detection machinery
# =============================================================================

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


def density_map(gx, gy, shape, sel, occ):
    """
    Normalized convolution. Counts and the tissue mask are smoothed with the
    same kernel and divided, with the denominator floored, so density is not
    depressed within one bandwidth of a tissue boundary.
    """
    counts = np.zeros(shape, dtype=float)
    np.add.at(counts, (gy[sel], gx[sel]), 1.0)
    sigma = DENSITY_BANDWIDTH_UM / GRID_UM
    sm = ndi.gaussian_filter(counts, sigma=sigma, mode="constant")
    if EDGE_CORRECTION:
        w = ndi.gaussian_filter(occ.astype(float), sigma=sigma, mode="constant")
        sm = sm / np.maximum(w, EDGE_MIN_WEIGHT)
    return sm / (PX_AREA_UM2 / 1e6)


def background_of(dens, occ, stat, tag=""):
    """
    Section background. A low percentile can legitimately be zero if a quarter
    of tissue pixels carry no myeloid density, which would make every fold
    infinite, so the statistic escalates and says so rather than failing
    silently.
    """
    v = dens[occ]; v = v[np.isfinite(v)]
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
    for alt in ["p50", "p75", "p90", "p95"]:
        b = _s(alt)
        if np.isfinite(b) and b > 0:
            print(f"    WARNING [{tag}]: background '{stat}' is {bg}. Escalated")
            print(f"    to '{alt}' = {b:.2f} /mm2, meaning most tissue pixels in")
            print("    this section carry no myeloid density at all. Treat every")
            print("    fold-based number for this section as unreliable.")
            return b
    pos = v[v > 0]
    if pos.size:
        b = float(pos.mean())
        print(f"    WARNING [{tag}]: no positive percentile. Falling back to the")
        print(f"    mean of positive pixels, {b:.2f} /mm2.")
        return b
    return np.nan


def find_local_maxima(dens, occ, sep_px, deduplicate=True):
    """
    Plateau deduplication matters: dens == maximum_filter(dens) returns EVERY
    pixel of a flat maximum, so a plateau would otherwise appear as many peaks.
    """
    mx = ndi.maximum_filter(dens, size=2 * sep_px + 1, mode="constant")
    ispeak = (dens == mx) & occ & (dens > 0)
    n_raw = int(ispeak.sum())
    if not deduplicate or n_raw == 0:
        py, px = np.nonzero(ispeak)
        return py, px, dens[py, px], n_raw, 0
    lab, n = ndi.label(ispeak, structure=np.ones((3, 3), dtype=int))
    objs = ndi.find_objects(lab)
    coms = ndi.center_of_mass(ispeak, lab, list(range(1, n + 1)))
    ys, xs = [], []
    for k, (com, sl) in enumerate(zip(coms, objs), start=1):
        yy, xx = np.nonzero(lab[sl] == k)
        yy = yy + sl[0].start; xx = xx + sl[1].start
        j = int(np.argmin((yy - com[0]) ** 2 + (xx - com[1]) ** 2))
        ys.append(int(yy[j])); xs.append(int(xx[j]))
    py = np.asarray(ys, int); px = np.asarray(xs, int)
    return py, px, dens[py, px], n_raw, n_raw - len(py)


def geometry(mask):
    area = float(mask.sum()) * PX_AREA_UM2
    r_eq = float(np.sqrt(area / np.pi)) if area > 0 else 0.0
    dt = ndi.distance_transform_edt(mask)
    r_in = float(dt.max()) * GRID_UM if mask.any() else 0.0
    return area, r_eq, r_in, (r_eq / r_in if r_in > 0 else np.nan), dt


def detect_burden(dens, occ, gx, gy):
    """ABSOLUTE. One fixed threshold, identical in every section. UNCHANGED."""
    m = (dens >= BURDEN_DENSITY) & occ
    lab, n = ndi.label(m)
    labels = np.zeros(dens.shape, int); recs = []; nid = 1
    for i in range(1, n + 1):
        comp = lab == i
        area = comp.sum() * PX_AREA_UM2
        ncell = int(comp[gy, gx].sum())
        if area < BURDEN_MIN_AREA_UM2 or ncell < BURDEN_MIN_CELLS:
            continue
        labels[comp] = nid
        recs.append(dict(burden_id=nid, area_um2=float(area), n_cells=ncell,
                         peak_density=float(dens[comp].max())))
        nid += 1
    return labels, recs


def detect_foci(dens, occ, gx, gy, bg):
    """PEAK-BASED. Supplies the centre and inscribed radius the radial
    coordinate needs. focus_id keeps revision 4 semantics."""
    diag = dict(n_peaks_raw=0, n_peaks_collapsed=0, n_peaks_kept=0,
                n_contested=0, n_drop_area=0, n_drop_cells=0, n_drop_peak=0)
    py, px, pv, n_raw, n_coll = find_local_maxima(dens, occ, SEP_PX,
                                                  DEDUPLICATE_PLATEAU_PEAKS)
    diag["n_peaks_raw"] = n_raw; diag["n_peaks_collapsed"] = n_coll
    keep = (pv >= FOCUS_FOLD_OVER_BACKGROUND * bg) & (pv >= MIN_PEAK_DENSITY)
    diag["n_drop_peak"] = int(((pv >= FOCUS_FOLD_OVER_BACKGROUND * bg)
                               & (pv < MIN_PEAK_DENSITY)).sum())
    py, px, pv = py[keep], px[keep], pv[keep]
    diag["n_peaks_kept"] = int(len(pv))
    if not len(pv):
        return np.zeros(dens.shape, int), [], diag
    o = np.argsort(pv)[::-1]
    py, px, pv = py[o], px[o], pv[o]
    markers = np.zeros(dens.shape, int)
    markers[py, px] = np.arange(1, len(pv) + 1)
    territory = (dens >= FOCUS_HALF_MAX_FRACTION * pv.min()) & occ
    if HAVE_SKIMAGE:
        part = _sk_watershed(-np.where(occ, dens, 0.0), markers, mask=territory)
    else:
        _, ind = ndi.distance_transform_edt(markers == 0, return_indices=True)
        part = np.where(territory, markers[ind[0], ind[1]], 0)

    labels = np.zeros(dens.shape, int); recs = []; nid = 1
    for i, (yy, xx, vv) in enumerate(zip(py, px, pv), start=1):
        region = (part == i) & (dens >= FOCUS_HALF_MAX_FRACTION * vv)
        if not region[yy, xx]:
            continue
        cl, cn = ndi.label(region)
        if cn == 0:
            continue
        region = cl == cl[yy, xx]
        if FOCUS_FILL_HOLES:
            region = ndi.binary_fill_holes(region)
        assigned = region & (labels == 0)
        diag["n_contested"] += int(region.sum() - assigned.sum())
        if not assigned[yy, xx]:
            continue
        cl2, cn2 = ndi.label(assigned)
        if cn2 == 0:
            continue
        assigned = cl2 == cl2[yy, xx]
        area, r_eq, r_in, shape_ratio, _ = geometry(assigned)
        ncell = int(assigned[gy, gx].sum())
        if area < FOCUS_MIN_AREA_UM2:
            diag["n_drop_area"] += 1; continue
        if ncell < FOCUS_MIN_CELLS:
            diag["n_drop_cells"] += 1; continue
        labels[assigned] = nid
        recs.append(dict(focus_id=nid, peak_row=int(yy), peak_col=int(xx),
                         peak_density=float(vv),
                         background_density=float(bg),
                         fold_over_background=float(vv / bg),
                         area_um2=area, equiv_radius_um=r_eq,
                         max_inscribed_radius_um=r_in,
                         shape_ratio=shape_ratio, n_cells=ncell))
        nid += 1
    return labels, recs, diag


def detect_complexes(dens, occ, gx, gy, bg):
    """
    LEVEL-SET. Connected components above a section-relative level, no peak
    required, so a confluent mass is one structure because it is one component.
    """
    diag = dict(n_components=0, n_drop_area=0, n_drop_cells=0, n_drop_peak=0)
    m = (dens >= COMPLEX_LEVEL_FOLD * bg) & occ
    if COMPLEX_CLOSE_UM > 0:
        r = max(1, int(round(COMPLEX_CLOSE_UM / GRID_UM)))
        m = ndi.binary_closing(m, structure=ndi.generate_binary_structure(2, 1),
                               iterations=r) & occ
    if FOCUS_FILL_HOLES:
        m = ndi.binary_fill_holes(m)
    lab, n = ndi.label(m)
    diag["n_components"] = n
    labels = np.zeros(dens.shape, int); recs = []; nid = 1
    for i in range(1, n + 1):
        comp = lab == i
        area, r_eq, r_in, shape_ratio, _ = geometry(comp)
        ncell = int(comp[gy, gx].sum())
        peak = float(dens[comp].max())
        if area < COMPLEX_MIN_AREA_UM2:
            diag["n_drop_area"] += 1; continue
        if ncell < COMPLEX_MIN_CELLS:
            diag["n_drop_cells"] += 1; continue
        if peak < MIN_PEAK_DENSITY:
            diag["n_drop_peak"] += 1; continue
        yy, xx = np.unravel_index(int(np.argmax(np.where(comp, dens, -np.inf))),
                                  dens.shape)
        labels[comp] = nid
        recs.append(dict(complex_id=nid, peak_row=int(yy), peak_col=int(xx),
                         peak_density=peak, background_density=float(bg),
                         fold_over_background=float(peak / bg),
                         level_density=float(COMPLEX_LEVEL_FOLD * bg),
                         area_um2=area, equiv_radius_um=r_eq,
                         max_inscribed_radius_um=r_in,
                         shape_ratio=shape_ratio, n_cells=ncell))
        nid += 1
    return labels, recs, diag


def radial_coordinate(focus_labels, focus_recs):
    """
    Distance from each focus centre normalised by that focus's own inscribed
    radius, so the coordinate is scale-free per focus. Microns conversions
    downstream must use max_inscribed_radius_um, not the equivalent radius.
    """
    rad = np.full(focus_labels.shape, np.nan)
    for r in focus_recs:
        m = focus_labels == r["focus_id"]
        if not m.any():
            continue
        dt = ndi.distance_transform_edt(m)
        r_in_px = dt.max()
        if r_in_px <= 0:
            continue
        # 0 at the centre, 1 at the boundary
        rad[m] = 1.0 - (dt[m] / r_in_px)
    return rad


# %% Cell 4 - per-section detection
# =============================================================================

banner("DETECTION")

paths = sorted(glob.glob(os.path.join(DATA_DIR, "*.csv")))
if not paths:
    print(f"ERROR: no CSVs in {DATA_DIR}")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

burden_rows, foci_rows, complex_rows, diag_rows, edge_rows = [], [], [], [], []
section_rows = []
skipped = []

for p in paths:
    sid = os.path.basename(p).replace(".csv", "")
    try:
        d = pd.read_csv(p, low_memory=False)
    except Exception as e:
        print(f"    ERROR reading {sid}: {e}")
        continue

    cx = find_col(d, CENTROID_X_PREFIX)
    cy = find_col(d, CENTROID_Y_PREFIX)
    cp = find_col(d, PHENOTYPE_COL)
    cpar = find_col(d, PARENT_COL)
    if None in (cx, cy, cp):
        print(f"    ERROR: {sid} is missing a required column.")
        print(f"      x  ('{CENTROID_X_PREFIX}...') : {cx}")
        print(f"      y  ('{CENTROID_Y_PREFIX}...') : {cy}")
        print(f"      phenotype ('{PHENOTYPE_COL}') : {cp}")
        print(f"      columns present: {', '.join(map(str, d.columns[:12]))}"
              f"{' ...' if len(d.columns) > 12 else ''}")
        skipped.append((sid, "missing a required column"))
        continue
    if cpar is not None and d[cpar].nunique() == 1:
        sid = str(d[cpar].iloc[0])
    if sid in EXCLUDE_SECTIONS:
        print(f"    {sid:<14} EXCLUDED (position 1, reduced intensity)")
        continue

    d = d.loc[np.isfinite(d[cx]) & np.isfinite(d[cy])].reset_index(drop=True)
    x = d[cx].to_numpy(float); y = d[cy].to_numpy(float)
    pheno = d[cp].astype(str)
    if NORMALISE_PHENOTYPE_LABELS:
        pheno = pheno.map(ascii_safe)
    pheno = pheno.to_numpy()
    cond = CONDITION_FROM_PREFIX.get(sid.split("_")[0], "?")

    gx0, gy0, gx, gy, shape = make_grid(x, y, GRID_UM)
    occ = tissue_mask(gx, gy, shape)
    sel = np.isin(pheno, POOL)
    if not sel.any():
        print(f"    ERROR: {sid} contains NO cells from the myeloid pool.")
        print(f"      pool     : {', '.join(POOL)}")
        print(f"      phenotypes present: {', '.join(sorted(set(pheno))[:8])}")
        print("      Phenotype names in the export do not match MYELOID_POOL.")
        print("      SECTION SKIPPED. No cell assignment file written.")
        skipped.append((sid, "no cells from the myeloid pool"))
        continue
    dens = density_map(gx, gy, shape, sel, occ)
    bg = background_of(dens, occ, BACKGROUND_STAT, sid)
    if not np.isfinite(bg) or bg <= 0:
        # Writing a file here would silently give downstream scripts a section
        # with zero structures, which is indistinguishable from a real absence
        # of disease. Fail loudly instead.
        print(f"    ERROR: {sid} has no usable background density ({bg}).")
        print("      Every fold-based structure would be undefined, and writing")
        print("      a cell assignment file would hand downstream scripts a")
        print("      section with zero structures that looks like real absence")
        print("      of disease. SECTION SKIPPED. No file written.")
        skipped.append((sid, f"no usable background ({bg})"))
        del dens, occ
        gc.collect()
        continue

    # edge-correction audit, kept from revision 4
    dens_raw = density_map(gx, gy, shape, sel, occ) if not EDGE_CORRECTION else None
    if EDGE_CORRECTION:
        counts = np.zeros(shape); np.add.at(counts, (gy[sel], gx[sel]), 1.0)
        sigma = DENSITY_BANDWIDTH_UM / GRID_UM
        raw = ndi.gaussian_filter(counts, sigma=sigma, mode="constant") \
            / (PX_AREA_UM2 / 1e6)
        diff = np.abs(dens - raw) / np.maximum(raw, 1e-9)
        edge_rows.append(dict(sample_id=sid, condition=cond,
                              pct_tissue_changed_gt5pct=100.0 * float(
                                  ((diff > 0.05) & occ).sum()) / max(1, occ.sum()),
                              background_corrected=bg,
                              background_uncorrected=background_of(
                                  raw, occ, BACKGROUND_STAT,
                                  f"{sid}/uncorrected"),
                              peak_corrected=float(dens[occ].max()),
                              peak_uncorrected=float(raw[occ].max())))
        del raw, diff, counts

    b_lab, b_recs = detect_burden(dens, occ, gx, gy)
    f_lab, f_recs, f_diag = detect_foci(dens, occ, gx, gy, bg)
    c_lab, c_recs, c_diag = detect_complexes(dens, occ, gx, gy, bg)

    # cuff: dilation around each focus, not overlapping any focus core
    cuff_px = max(1, int(round(CUFF_DILATION_UM / GRID_UM)))
    core = f_lab > 0
    grown = ndi.binary_dilation(core, iterations=cuff_px) & occ
    cuff = grown & ~core
    cuff_lab = np.zeros(f_lab.shape, int)
    if core.any():
        _, ind = ndi.distance_transform_edt(~core, return_indices=True)
        cuff_lab[cuff] = f_lab[ind[0], ind[1]][cuff]

    rad = radial_coordinate(f_lab, f_recs)

    tissue_mm2 = occ.sum() * PX_AREA_UM2 / 1e6
    burden_area = sum(r["area_um2"] for r in b_recs) / 1e6
    section_rows.append(dict(
        sample_id=sid, condition=cond, n_cells=len(d),
        tissue_area_mm2=tissue_mm2, background_density=bg,
        background_stat=BACKGROUND_STAT,
        burden_area_mm2=burden_area,
        burden_pct_of_tissue=100.0 * burden_area / tissue_mm2 if tissue_mm2 else np.nan,
        n_burden_regions=len(b_recs), n_foci=len(f_recs),
        n_complexes=len(c_recs),
        complex_level_density=COMPLEX_LEVEL_FOLD * bg,
        focus_gate_density=FOCUS_FOLD_OVER_BACKGROUND * bg,
        min_peak_density=MIN_PEAK_DENSITY,
        peak_density=float(dens[occ].max())))

    for r in b_recs:
        r.update(sample_id=sid, condition=cond); burden_rows.append(r)
    for r in f_recs:
        m = f_lab == r["focus_id"]
        vc = pd.Series(pheno[m[gy, gx]]).value_counts()
        n = int(m[gy, gx].sum())
        r.update(sample_id=sid, condition=cond,
                 frac_myeloid=float(sum(vc.get(q, 0) for q in POOL)) / n if n else np.nan,
                 pct_ido1_pos_of_mac_core=(
                     100.0 * vc.get(IDO1_POS, 0)
                     / max(1, vc.get(IDO1_POS, 0) + vc.get(IDO1_NEG, 0))),
                 **{f"n_{q}": int(vc.get(q, 0)) for q in PHENOTYPE_ORDER})
        foci_rows.append(r)
    for r in c_recs:
        m = c_lab == r["complex_id"]
        vc = pd.Series(pheno[m[gy, gx]]).value_counts()
        n = int(m[gy, gx].sum())
        r.update(sample_id=sid, condition=cond,
                 frac_myeloid=float(sum(vc.get(q, 0) for q in POOL)) / n if n else np.nan,
                 n_foci_inside=int(len(np.unique(f_lab[m])) - (1 if 0 in f_lab[m] else 0)),
                 **{f"n_{q}": int(vc.get(q, 0)) for q in PHENOTYPE_ORDER})
        complex_rows.append(r)
    diag_rows.append(dict(sample_id=sid, condition=cond, **f_diag,
                          **{f"complex_{k}": v for k, v in c_diag.items()}))

    # ---- the central artefact ----------------------------------------------
    out = pd.DataFrame({
        "sample_id": sid, "condition": cond,
        "x": x, "y": y, "pheno": pheno,
        "grid_row": gy, "grid_col": gx,
        "in_burden_region": b_lab[gy, gx] > 0,
        "burden_id": b_lab[gy, gx],
        "focus_id": f_lab[gy, gx],
        "focus_region": np.where(f_lab[gy, gx] > 0, "core",
                                 np.where(cuff_lab[gy, gx] > 0, "cuff", "none")),
        "cuff_focus_id": cuff_lab[gy, gx],
        "radial_pos": rad[gy, gx],
        "complex_id": c_lab[gy, gx],
        "myeloid_density_at_cell": dens[gy, gx],
    })
    out["is_core"] = out["focus_id"] > 0
    out["is_myeloid_pool"] = np.isin(pheno, POOL)
    write_csv(out, f"{sid}_cell_structures.csv", CELL_DIR)

    print(f"    {sid:<14} {cond:<10} {len(d):>9,} cells  bg {bg:>8.1f}  "
          f"burden {100.0 * burden_area / tissue_mm2 if tissue_mm2 else 0:>5.1f}%  "
          f"foci {len(f_recs):>3}  complexes {len(c_recs):>3}")

    del dens, occ, b_lab, f_lab, c_lab, cuff_lab, rad, out
    gc.collect()

if not section_rows:
    print("\nERROR: no sections processed.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)


# %% Cell 5 - tables
# =============================================================================

banner("TABLES")

sections = pd.DataFrame(section_rows)
sections = sections.sort_values(
    by=["condition", "sample_id"],
    key=lambda c: c.map({v: i for i, v in enumerate(CONDITION_ORDER)})
    if c.name == "condition" else c)
write_csv(sections, "33b_burden_summary.csv")
write_csv(pd.DataFrame(burden_rows), "33_burden_regions_absolute.csv")
foci = pd.DataFrame(foci_rows)
write_csv(foci, "35_foci_structures_relative.csv")
comps = pd.DataFrame(complex_rows)
write_csv(comps, "37_granuloma_complexes.csv")
write_csv(pd.DataFrame(diag_rows), "35b_detection_diagnostics.csv")
if edge_rows:
    write_csv(pd.DataFrame(edge_rows), "32_density_edge_correction.csv")

sub("Per section")
show = ["sample_id", "condition", "n_cells", "tissue_area_mm2",
        "background_density", "burden_pct_of_tissue", "n_foci", "n_complexes"]
with pd.option_context("display.width", 260):
    print(sections[show].to_string(index=False))

sub("Burden, the ABSOLUTE definition, unchanged from revision 4")
for cond in CONDITION_ORDER:
    s = sections.loc[sections["condition"] == cond, "burden_pct_of_tissue"]
    if len(s):
        print(f"    {cond:<12} " + ", ".join(f"{v:.1f}%" for v in s))
print(f"\n    Threshold {BURDEN_DENSITY:,.0f} cells/mm2, applied identically to")
print("    every section. This definition was NOT modified in revision 5.")

sub("What the new minimum peak density removed")
dg = pd.DataFrame(diag_rows)
if "n_drop_peak" in dg.columns:
    print(f"    foci dropped for peak below {MIN_PEAK_DENSITY:,.0f} /mm2: "
          f"{int(dg['n_drop_peak'].sum())}")
if "complex_n_drop_peak" in dg.columns:
    print(f"    complexes dropped for the same reason: "
          f"{int(dg['complex_n_drop_peak'].sum())}")
print("\n    In calibration this criterion excluded 18 objects, none of which")
print("    contained IDO1+ macrophages, while retaining 57 of which 44 did.")
print("    IDO1 is used nowhere in detection.")

if len(comps):
    sub("Complexes against foci")
    g = comps.groupby(["condition", "sample_id"]).agg(
        n_complexes=("complex_id", "size"),
        median_foci_inside=("n_foci_inside", "median"),
        max_foci_inside=("n_foci_inside", "max"),
        total_area_mm2=("area_um2", lambda v: v.sum() / 1e6))
    print(g.to_string())
    print("\n    'foci inside' above 1 is confluence: several density peaks")
    print("    within one connected structure. That is the quantity the")
    print("    level-set definition exists to capture.")


# %% Cell 6 - figures
# =============================================================================

banner("FIGURES")

SECT = list(sections["sample_id"])

fig, axes = plt.subplots(1, 2, figsize=(28, 13))
ax = axes[0]
xs = np.arange(len(SECT))
ax.bar(xs, sections["burden_pct_of_tissue"],
       color=[CONDITION_COLORS.get(c, "#999") for c in sections["condition"]],
       edgecolor="#000000", linewidth=2)
ax.set_xticks(xs)
ax.set_xticklabels([f"{short_label(s)}\n{c[:4]}"
                    for s, c in zip(SECT, sections["condition"])])
ax.set_ylabel("Tissue at granuloma density (%)")
ax.set_title("Burden, absolute definition", fontsize=FONT_SIZE_TITLE - 2)
style_axes(ax); panel_letter(ax, "A")

ax = axes[1]
w = 0.38
ax.bar(xs - w / 2, sections["n_foci"], width=w, color=FOCUS_COLOR,
       edgecolor="#000000", linewidth=2, label="foci (peak-based)")
ax.bar(xs + w / 2, sections["n_complexes"], width=w, color=COMPLEX_COLOR,
       edgecolor="#000000", linewidth=2, label="complexes (level set)")
ax.set_xticks(xs)
ax.set_xticklabels([f"{short_label(s)}\n{c[:4]}"
                    for s, c in zip(SECT, sections["condition"])])
ax.set_ylabel("Structures detected")
ax.set_title("Two definitions, reported separately",
             fontsize=FONT_SIZE_TITLE - 2)
ax.legend(frameon=False, fontsize=FONT_SIZE_LEGEND - 8)
style_axes(ax); panel_letter(ax, "B")
fig.suptitle("Structure detection, revision 5", fontsize=FONT_SIZE_TITLE, y=1.02)
fig.tight_layout()
save_fig(fig, "F31_burden_summary")

if len(foci):
    fig, axes = plt.subplots(1, 3, figsize=(34, 12))
    for ax, (col, lab) in zip(axes, [
            ("peak_density", "Peak myeloid density (cells / mm$^2$)"),
            ("max_inscribed_radius_um", "Maximum inscribed radius (um)"),
            ("shape_ratio", "Shape (equivalent / inscribed)")]):
        for k, cond in enumerate(CONDITION_ORDER):
            v = foci.loc[foci["condition"] == cond, col].dropna()
            if not len(v):
                continue
            xj = k + rng.normal(0, 0.06, len(v))
            ax.scatter(xj, v, s=200, alpha=0.7, color=CONDITION_COLORS[cond],
                       edgecolors="#000000", linewidths=1.5)
            ax.hlines(v.median(), k - 0.25, k + 0.25, color="#000000",
                      linewidth=5)
        ax.set_xticks(range(len(CONDITION_ORDER)))
        ax.set_xticklabels([CONDITION_LABELS[c] for c in CONDITION_ORDER],
                           fontsize=FONT_SIZE_TICK - 8)
        ax.set_ylabel(lab, fontsize=FONT_SIZE_BASE - 6)
        if col == "peak_density":
            ax.set_yscale("log")
            ax.axhline(MIN_PEAK_DENSITY, color="#000000", linestyle="--",
                       linewidth=3)
            ax.axhline(BURDEN_DENSITY, color="#000000", linestyle=":",
                       linewidth=3)
        style_axes(ax)
    fig.suptitle("Focus properties by arm\ndashed = minimum peak density, "
                 "dotted = burden threshold", fontsize=FONT_SIZE_TITLE, y=1.03)
    fig.tight_layout()
    save_fig(fig, "F34_focus_properties")


# %% Cell 7 - wrap up
# =============================================================================

banner("SUMMARY")

print(f"Sections written   : {len(sections)}")
if skipped:
    print(f"Sections SKIPPED   : {len(skipped)}")
    for q, why in skipped:
        print(f"    {q}: {why}")
    print("  These sections have NO cell assignment file. Downstream scripts")
    print("  will simply not see them, so resolve this before re-running 06.")
print(f"Cells assigned     : {int(sections['n_cells'].sum()):,}")
print(f"Foci               : {len(foci)}")
print(f"Complexes          : {len(comps)}")
print(f"Myeloid pool       : {', '.join(POOL)}")
if LEAVE_ONE_OUT:
    print(f"LEAVE-ONE-OUT      : {LEAVE_ONE_OUT} excluded. SENSITIVITY RUN.")

sub("What must happen next")
print("  1. Scripts 06, 07, 08 and 09 re-run against structures_rev5. Their")
print("     input paths must be repointed from structures_rev4.")
print("  2. Every manuscript number re-baselined from that run. Quote nothing")
print("     from the revision 4 tables.")
print("  3. Treated focus counts and peak density summaries in particular:")
print("     revision 4 reported three foci on one treated section of which two")
print("     had peak densities below the new minimum and contained no IDO1+")
print("     macrophages.")

sub("Circularity, for whoever writes the results")
print("  Structures are defined by pooled myeloid density. Lymphocyte results")
print("  are unaffected, those phenotypes are not in the pool. Radial position")
print("  of a pooled member is nearly unaffected, and the one coupling is")
print("  conservative for a rim-ward finding. ABUNDANCE of a pooled member")
print("  inside structures IS exposed, and should be supported by a")
print("  leave-one-out run with LEAVE_ONE_OUT set to that phenotype.")

sub("Not changed by revision 5")
print(f"  The absolute burden definition, {BURDEN_DENSITY:,.0f} cells/mm2 with a")
print(f"  {BURDEN_MIN_AREA_UM2:,.0f} um2 and {BURDEN_MIN_CELLS} cell filter,")
print("  applied identically to every section.")

banner("DONE")
sys.stdout = _tee.terminal
_tee.close()
