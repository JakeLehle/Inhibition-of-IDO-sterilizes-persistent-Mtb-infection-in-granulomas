#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - EXPERT ANNOTATION MEASUREMENTS
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 04d of the AKOYA analysis series. READ-ONLY.

WHAT THIS REPLACES

    Detection parameters were previously chosen by sweeping thresholds and
    scoring against six expert COUNTS. That was three degrees of freedom tuning
    three parameters, it could not tell over-detection from under-detection, and
    every top-ranked setting reached its score by collapsing all three treated
    sections to a single focus each. It was a dead end and is recorded as one.

    This measures the same parameters directly, on regions an expert drew, over
    hundreds of grid pixels per region rather than one integer per section.

WHAT IS MEASURED, AND WHICH PARAMETER EACH ONE FIXES

    1. BOUNDARY DENSITY            -> the admission threshold
       The absolute density and the fold over background at the traced
       boundary. This is where an expert says a lesion stops.

    2. BOUNDARY AS A FRACTION OF ENCLOSED PEAK  -> FOCUS_HALF_MAX_FRACTION
       Script 04 grows each focus to 50 percent of its own peak. That number was
       assumed, never measured. If traced boundaries sit at 25 or 30 percent,
       that is the number and the sweep stops.

    3. LOCAL MAXIMA INSIDE ONE EXPERT REGION    -> PEAK_SEPARATION_UM
       Script 04 finds peaks with a maximum_filter of 2*sep_px+1 at 200 um. If
       one drawn complex contains several maxima, the separation is too small
       for confluent disease, and the spacing of those maxima gives the value it
       should be. This is the 43109 over-segmentation question.

    4. DENSITY INSIDE 'normal' REGIONS          -> BACKGROUND_STAT
       An expert-derived background, against which median, p75, p25, p10 and the
       modal estimate can be scored. This is the 36463 question: the median over
       tissue pixels in a heavily involved section is itself pathology.

    5. DENSITY PROFILE ACROSS THE BOUNDARY      -> whether 'fuzzy' is measurable
       Sampled along the outward normal, split by annotation class. If a fuzzy
       edge really does have a shallower gradient than a clean one, that is a
       property a detector can use rather than a word in a discussion.

    6. OVERLAP WITH THE CURRENT DETECTOR        -> the honest error term
       Per-region IoU, plus how many drawn regions the detector splits and how
       many detected foci fall inside one drawn region.

THE BRUSH-STAMP PROBLEM, HANDLED AUTOMATICALLY

    A brush CLICK in QuPath stamps a fixed disc; a brush DRAG paints a region.
    Stamps all share one exact area, and on G4_43109 eight annotations were
    confirmed to be the identical shape merely translated.

    A stamp records WHERE a region is but not WHERE ITS BOUNDARY IS. So stamps
    are detected here by exact duplicate area within a section and are:
        EXCLUDED from measurements 1, 2, 5 and 6, which all need a real boundary
        KEPT     for measurement 3, which needs only an enclosing location
    This is correct whether the stamps were intended as boundaries or as
    markers, and the split is reported so it is never silent.

COORDINATE HANDLING

    Annotations are in canvas pixels. Each GeoJSON carries its own transform in
    the FeatureCollection properties, written from the real cell coordinates by
    script 04c. Polygons are converted to scan microns and then onto the SAME
    25 um grid script 04 builds its density map on, using the identical origin
    convention, so a measurement here is directly a statement about script 04.

    An alignment check runs first and reports any annotation falling outside its
    section's tissue. Nothing downstream is trusted until it passes.

INPUTS
    annotations/*_annotations.geojson
    annotations/*_annotation_measurements.tsv        (optional cross-check)
    structures_rev4/cell_assignments/*_cell_structures.csv
    structures_rev4/tables/35_foci_structures_relative.csv

OUTPUTS
    annotation_measurements/tables/
        110_annotation_inventory, 111_boundary_measurements,
        112_maxima_inside_regions, 113_region_composition,
        114_expert_background, 115_detector_agreement,
        116_recommended_parameters, 00_annotation_measurements_report.txt
    annotation_measurements/figures/
        F90_boundary_threshold, F91_peak_fraction, F92_maxima_and_merge,
        F93_expert_background, F94_edge_gradient, F95_overlay_maps

USAGE
    conda activate sc_pre
    python AKOYA_04d_Annotation_Measurements.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================

ANN_DIR = "/master/jlehle/WORKING/AKOYA/annotations"
IN_DIR = "/master/jlehle/WORKING/AKOYA/structures_rev4"
CELL_DIR = f"{IN_DIR}/cell_assignments"
FOCI_TABLE = f"{IN_DIR}/tables/35_foci_structures_relative.csv"
OUT_DIR = "/master/jlehle/WORKING/AKOYA/annotation_measurements"

USE_AGG = True

# ---- must match script 04 ---------------------------------------------------
GRID_UM = 25.0
DENSITY_BANDWIDTH_UM = 75.0
EDGE_MIN_WEIGHT = 0.25
TISSUE_CLOSE_ITERATIONS = 2
TISSUE_FILL_HOLES = True
BACKGROUND_STAT = "median"
PEAK_SEPARATION_UM = 200.0
FOCUS_FOLD_OVER_BACKGROUND = 6.0
FOCUS_HALF_MAX_FRACTION = 0.50
FOCUS_MIN_AREA_UM2 = 10000.0
FOCUS_MIN_CELLS = 40
FOCUS_FILL_HOLES = True
BURDEN_THRESHOLD = 3000.0

# ---- measurement settings ---------------------------------------------------
BOUNDARY_RESAMPLE_UM = 10.0     # even spacing along the ring, removes the bias
                                # from vertices clustering where drawing slowed
NORMAL_PROFILE_UM = 300.0       # how far in and out to sample the gradient
NORMAL_STEP_UM = 10.0
MIN_REGION_AREA_UM2 = 2000.0    # below this a region is too small to measure
STAMP_AREA_TOL = 0.5            # px2; exact duplicate area means a brush stamp

CLASS_ORDER = ["complex", "missed", "clean_edge", "normal", "unclassified"]
BOUNDARY_CLASSES = ["complex", "missed", "clean_edge"]   # 'normal' has no edge

# background statistics to score against the expert regions
BG_STATS = ["median", "p75", "p50", "p25", "p10", "mode_kde"]

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

# ---- appearance -------------------------------------------------------------
CONDITION_ORDER = ["D1MT", "Untreated"]
CONDITION_COLORS = {"D1MT": "#2C7FB8", "Untreated": "#D95F02"}
CONDITION_LABELS = {"D1MT": "D1MT-treated", "Untreated": "Untreated"}
CLASS_COLORS = {"complex": "#D73027", "missed": "#4575B4",
                "clean_edge": "#1A9850", "normal": "#999999",
                "unclassified": "#000000"}
PHENOTYPE_COLORS = {
    IDO1_POS: "#B2182B", IDO1_NEG: "#EF8A62", "CD163+ Macrophages": "#FDBE85",
    "Neutrophils": "#7B3294", "Helper T cells": "#1B7837",
    "CD4- T cells": "#7FBC41", "Tregs": "#00441B", "B cells": "#2166AC",
    "Plasma cells": "#67A9CF", "Endothelial cells": "#E7298A",
    "Epithelial/Tumor cells": "#66C2A5", "Other": "#A6761D",
}
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


# %% Cell 2 - imports and helpers
# =============================================================================

import os
import sys
import gc
import json
import glob
from datetime import datetime

import numpy as np
import pandas as pd
import scipy.ndimage as ndi

import matplotlib
if USE_AGG:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.lines import Line2D

try:
    from skimage.segmentation import watershed as _sk_watershed
    HAVE_SKIMAGE = True
except Exception:
    HAVE_SKIMAGE = False

import warnings as _w
with _w.catch_warnings():          # the option is deprecated in pandas 4, where
    _w.simplefilter("ignore")      # its behaviour is already the default
    try:
        pd.set_option("future.no_silent_downcasting", True)
    except Exception:
        pass
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

SEP_PX = max(1, int(round(PEAK_SEPARATION_UM / GRID_UM)))


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


def q(v, p):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    return float(np.percentile(v, p)) if v.size else np.nan


_tee = Tee(os.path.join(TAB_DIR, "00_annotation_measurements_report.txt"))
sys.stdout = _tee

banner("AKOYA EXPERT ANNOTATION MEASUREMENTS (script 04d, read-only)")
print(f"Run time : {datetime.now().isoformat(timespec='seconds')}")
print(f"Output   : {OUT_DIR}")
print("\nWrites nothing into structures_rev4. Adopts no setting. Every number")
print("here is a measurement of what the expert drew, on script 04's own grid.")
if not HAVE_SKIMAGE:
    print("\n    WARNING: skimage missing, the detector comparison will use the")
    print("    EDT fallback and will not exactly reproduce script 04.")


# %% Cell 3 - script 04's density and detection, replicated
# =============================================================================
# Identical logic to AKOYA_04_Structures.py rev4 and to script 04b, verified
# there against table 35 to nine significant figures on background density and
# an exact peak-pixel match on all 67 foci.

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
    counts = np.zeros(shape, dtype=float)
    np.add.at(counts, (gy[sel], gx[sel]), 1.0)
    sigma = DENSITY_BANDWIDTH_UM / GRID_UM
    sm = ndi.gaussian_filter(counts, sigma=sigma, mode="constant")
    w = ndi.gaussian_filter(occ.astype(float), sigma=sigma, mode="constant")
    return (sm / np.maximum(w, EDGE_MIN_WEIGHT)) / ((GRID_UM ** 2) / 1e6)


def background_of(dens, occ, stat="median"):
    v = dens[occ]; v = v[np.isfinite(v)]
    if not v.size:
        return np.nan
    if stat in ("median", "p50"):
        return float(np.median(v))
    if stat.startswith("p") and stat[1:].replace(".", "").isdigit():
        return float(np.percentile(v, float(stat[1:])))
    if stat == "mode_kde":
        vv = v[v > 0]
        if vv.size:
            vv = vv[vv >= np.percentile(vv, 1.0)]
        if vv.size < 10:
            return float(np.median(v))
        lv = np.log10(vv)
        h, e = np.histogram(lv, bins=200)
        k = int(np.argmax(h))
        return float(10 ** (0.5 * (e[k] + e[k + 1])))
    return float(np.median(v))


def find_local_maxima(dens, occ, sep_px):
    mx = ndi.maximum_filter(dens, size=2 * sep_px + 1, mode="constant")
    ispeak = (dens == mx) & occ & (dens > 0)
    lab, n = ndi.label(ispeak, structure=np.ones((3, 3), dtype=int))
    if n == 0:
        return np.array([], int), np.array([], int), np.array([], float)
    objs = ndi.find_objects(lab)
    coms = ndi.center_of_mass(ispeak, lab, list(range(1, n + 1)))
    ys, xs = [], []
    for k, (com, sl) in enumerate(zip(coms, objs), start=1):
        yy, xx = np.nonzero(lab[sl] == k)
        yy = yy + sl[0].start; xx = xx + sl[1].start
        j = int(np.argmin((yy - com[0]) ** 2 + (xx - com[1]) ** 2))
        ys.append(int(yy[j])); xs.append(int(xx[j]))
    py = np.asarray(ys, int); px = np.asarray(xs, int)
    return py, px, dens[py, px]


def detect_foci(dens, occ, gx, gy, bg):
    """Script 04's relative foci, used only as the comparison baseline."""
    py, px, pv = find_local_maxima(dens, occ, SEP_PX)
    keep = pv >= FOCUS_FOLD_OVER_BACKGROUND * bg
    py, px, pv = py[keep], px[keep], pv[keep]
    if not len(pv):
        return np.zeros_like(dens, int), []
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
        if not assigned[yy, xx]:
            continue
        cl2, _ = ndi.label(assigned)
        assigned = cl2 == cl2[yy, xx]
        area = assigned.sum() * GRID_UM * GRID_UM
        ncell = int(assigned[gy, gx].sum())
        if area < FOCUS_MIN_AREA_UM2 or ncell < FOCUS_MIN_CELLS:
            continue
        labels[assigned] = nid
        recs.append(dict(focus_id=nid, peak_row=int(yy), peak_col=int(xx),
                         peak_density=float(vv), area_um2=float(area),
                         n_cells=ncell))
        nid += 1
    return labels, recs


# %% Cell 4 - load cells, build grids, load annotations
# =============================================================================

banner("LOADING")

ann_files = sorted(glob.glob(os.path.join(ANN_DIR, "*_annotations.geojson")))
if not ann_files:
    print(f"ERROR: no *_annotations.geojson in {ANN_DIR}")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

SECTIONS = [os.path.basename(f).replace("_annotations.geojson", "")
            for f in ann_files]
print(f"    annotation files: {', '.join(SECTIONS)}")

foci_ref = pd.read_csv(FOCI_TABLE) if os.path.exists(FOCI_TABLE) else None

sec = {}
for sid in SECTIONS:
    p = os.path.join(CELL_DIR, f"{sid}_cell_structures.csv")
    if not os.path.exists(p):
        print(f"    ERROR: no cell file for {sid}, skipped")
        continue
    d = pd.read_csv(p, low_memory=False)
    d = d.loc[np.isfinite(d["x"]) & np.isfinite(d["y"])].reset_index(drop=True)
    x = d["x"].to_numpy(float); y = d["y"].to_numpy(float)
    gx0, gy0, gx, gy, shape = make_grid(x, y, GRID_UM)
    occ = tissue_mask(gx, gy, shape)
    mye = d["pheno"].isin(MYELOID_FOR_DETECTION).to_numpy()
    dens = density_map(gx, gy, shape, mye, occ)
    bg = background_of(dens, occ, BACKGROUND_STAT)
    py, px, pv = find_local_maxima(dens, occ, SEP_PX)
    labels, recs = detect_foci(dens, occ, gx, gy, bg)
    sec[sid] = dict(d=d, x=x, y=y, gx=gx, gy=gy, gx0=gx0, gy0=gy0,
                    shape=shape, occ=occ, dens=dens, bg=bg,
                    pheno=d["pheno"].to_numpy(),
                    cond=d["condition"].iloc[0] if "condition" in d.columns else "?",
                    peak_rc=(py, px), peak_v=pv,
                    foci_labels=labels, foci_recs=recs)
    print(f"    {sid:<14} {len(d):>9,} cells  grid {shape[1]}x{shape[0]}  "
          f"bg {bg:>8.1f}/mm2  maxima {len(pv):>4}  foci {len(recs):>3}")

SECTIONS = [s for s in SECTIONS if s in sec]
SECTIONS.sort(key=lambda s: (CONDITION_ORDER.index(sec[s]["cond"])
                             if sec[s]["cond"] in CONDITION_ORDER else 9, s))


def geom_to_rings(geom):
    """GeoJSON geometry -> list of (exterior, [holes]) in canvas pixels."""
    t = geom["type"]; c = geom["coordinates"]
    polys = [c] if t == "Polygon" else (c if t == "MultiPolygon" else [])
    out = []
    for poly in polys:
        if not poly:
            continue
        out.append((np.asarray(poly[0], float),
                    [np.asarray(h, float) for h in poly[1:]]))
    return out


ann_rows = []
for f in ann_files:
    sid = os.path.basename(f).replace("_annotations.geojson", "")
    if sid not in sec:
        continue
    gj = json.load(open(f))
    pr = gj.get("properties", {})
    X0, Y0, SC = pr.get("x0_um"), pr.get("y0_um"), pr.get("um_per_px")
    if X0 is None or SC is None:
        print(f"    ERROR: {sid} geojson has no transform, skipped")
        continue
    for feat in gj["features"]:
        p = feat["properties"]
        rings = geom_to_rings(feat["geometry"])
        if not rings:
            continue
        # canvas pixels -> scan microns
        rings_um = [((r[:, 0] + 0.5) * SC + X0, (r[:, 1] + 0.5) * SC + Y0,
                     [((h[:, 0] + 0.5) * SC + X0, (h[:, 1] + 0.5) * SC + Y0)
                      for h in hs]) for r, hs in rings]
        ann_rows.append(dict(sample_id=sid, condition=sec[sid]["cond"],
                             index=p.get("index"),
                             classification=p.get("classification", "unclassified"),
                             area_px2=p.get("area_px2", np.nan),
                             area_um2=p.get("area_um2", np.nan),
                             rings_um=rings_um))

ann = pd.DataFrame(ann_rows)
print(f"\n    loaded {len(ann)} annotations across {ann['sample_id'].nunique()} sections")

# ---- brush-stamp detection --------------------------------------------------
ann["is_stamp"] = False
for sid, g in ann.groupby("sample_id"):
    a = g["area_px2"].round(1)
    dup = a.value_counts()
    dup = dup[dup > 1].index
    ann.loc[g.index[a.isin(dup)], "is_stamp"] = True

sub("Annotations by class, and the brush-stamp split")
print("    A brush CLICK stamps a fixed disc, a brush DRAG paints a region.")
print("    Stamps share one exact area. They mark WHERE a region is but not")
print("    where its boundary is, so they are excluded from every boundary")
print("    measurement and kept only for counting maxima inside a region.\n")
t = (ann.groupby(["sample_id", "classification"])
     .agg(n=("index", "size"), stamps=("is_stamp", "sum")).reset_index())
print(t.to_string(index=False))
print(f"\n    traced {int((~ann['is_stamp']).sum())}   stamped "
      f"{int(ann['is_stamp'].sum())}   of {len(ann)}")


# %% Cell 5 - rasterise onto script 04's grid, and check alignment
# =============================================================================

banner("RASTERISING ONTO SCRIPT 04's GRID")

print("Annotations are converted to scan microns using the transform each")
print("GeoJSON carries, then onto the SAME 25 um grid script 04 builds its")
print("density map on, with the identical origin convention. A measurement here")
print("is therefore directly a statement about script 04.\n")


def um_to_grid(xu, yu, s):
    return ((xu - s["gx0"]) / GRID_UM) + 1.0, ((yu - s["gy0"]) / GRID_UM) + 1.0


def rasterise(rings_um, s):
    """Boolean mask on the section grid. Holes removed explicitly rather than
    relying on ring winding order, which QuPath does not guarantee."""
    H, W = s["shape"]
    yy, xx = np.mgrid[0:H, 0:W]
    pts = np.column_stack([xx.ravel(), yy.ravel()])
    mask = np.zeros(H * W, bool)
    for ex_x, ex_y, holes in rings_um:
        cx, cy = um_to_grid(ex_x, ex_y, s)
        m = MplPath(np.column_stack([cx, cy])).contains_points(pts)
        for hx, hy in holes:
            hcx, hcy = um_to_grid(hx, hy, s)
            m &= ~MplPath(np.column_stack([hcx, hcy])).contains_points(pts)
        mask |= m
    return mask.reshape(H, W)


masks = {}
align_rows = []
for i, r in ann.iterrows():
    s = sec[r["sample_id"]]
    m = rasterise(r["rings_um"], s)
    masks[i] = m
    npx = int(m.sum())
    in_tissue = int((m & s["occ"]).sum())
    align_rows.append(dict(
        sample_id=r["sample_id"], index=r["index"],
        classification=r["classification"], is_stamp=bool(r["is_stamp"]),
        area_um2_reported=r["area_um2"],
        area_um2_rasterised=npx * GRID_UM * GRID_UM,
        n_grid_px=npx, pct_in_tissue=100.0 * in_tissue / npx if npx else np.nan))

align = pd.DataFrame(align_rows)
align["area_ratio"] = (align["area_um2_rasterised"]
                       / align["area_um2_reported"].replace(0, np.nan))
write_csv(align, "110_annotation_inventory.csv")

sub("ALIGNMENT CHECK")
bad_tissue = align.loc[align["pct_in_tissue"] < 80]
bad_area = align.loc[(align["area_ratio"] < 0.7) | (align["area_ratio"] > 1.4)]
print(f"    median rasterised / reported area ratio : "
      f"{align['area_ratio'].median():.3f}   (1.0 is perfect; the 25 um grid")
print(f"    quantises small regions, so a spread here is expected)")
print(f"    median percent of each region inside tissue : "
      f"{align['pct_in_tissue'].median():.1f}%")
if len(bad_tissue):
    print(f"\n    {len(bad_tissue)} region(s) less than 80 percent inside the")
    print("    tissue mask. Listed below; check these on the canvas.")
    print(bad_tissue[["sample_id", "index", "classification",
                      "pct_in_tissue"]].to_string(index=False))
else:
    print("    every region sits inside its section's tissue mask.")
if len(bad_area) and len(bad_area) > 0.2 * len(align):
    print(f"\n    WARNING: {len(bad_area)} regions have a rasterised area far")
    print("    from the reported area. If these are not all small regions, the")
    print("    transform may be wrong.")
    print(bad_area[["sample_id", "index", "classification", "area_um2_reported",
                    "area_um2_rasterised", "area_ratio"]].head(12).to_string(index=False))


# %% Cell 6 - measurement 1 and 2, boundary density and fraction of peak
# =============================================================================

banner("BOUNDARY DENSITY AND FRACTION OF ENCLOSED PEAK")

print("Where does an expert say a lesion stops? Sampled along the traced ring")
print("at even spacing, so vertices clustering where the drawing slowed do not")
print("bias the answer. Stamps are excluded: a stamped disc has no boundary.\n")


def resample_ring(xu, yu, step_um):
    d = np.r_[0.0, np.cumsum(np.hypot(np.diff(xu), np.diff(yu)))]
    if d[-1] <= 0:
        return xu, yu
    t = np.arange(0.0, d[-1], step_um)
    return np.interp(t, d, xu), np.interp(t, d, yu)


def sample_density(xu, yu, s):
    cx, cy = um_to_grid(xu, yu, s)
    return ndi.map_coordinates(np.nan_to_num(s["dens"]), [cy, cx],
                               order=1, mode="nearest")


bnd_rows = []
prof_rows = []
for i, r in ann.iterrows():
    if r["is_stamp"] or r["classification"] not in BOUNDARY_CLASSES:
        continue
    s = sec[r["sample_id"]]
    m = masks[i]
    if m.sum() * GRID_UM * GRID_UM < MIN_REGION_AREA_UM2:
        continue
    inside = s["dens"][m & s["occ"]]
    inside = inside[np.isfinite(inside)]
    if not inside.size:
        continue
    peak = float(inside.max())

    bx, by, bd = [], [], []
    for ex_x, ex_y, _ in r["rings_um"]:
        rx, ry = resample_ring(ex_x, ex_y, BOUNDARY_RESAMPLE_UM)
        bx.append(rx); by.append(ry); bd.append(sample_density(rx, ry, s))
    bx = np.concatenate(bx); by = np.concatenate(by); bd = np.concatenate(bd)

    bnd_rows.append(dict(
        sample_id=r["sample_id"], condition=r["condition"], index=r["index"],
        classification=r["classification"], background=s["bg"],
        area_um2=m.sum() * GRID_UM * GRID_UM, n_boundary_samples=len(bd),
        boundary_density_median=float(np.median(bd)),
        boundary_density_p25=q(bd, 25), boundary_density_p75=q(bd, 75),
        boundary_fold_median=(float(np.median(bd)) / s["bg"]
                              if np.isfinite(s["bg"]) and s["bg"] > 0 else np.nan),
        enclosed_peak_density=peak,
        enclosed_peak_fold=(peak / s["bg"]
                            if np.isfinite(s["bg"]) and s["bg"] > 0 else np.nan),
        boundary_frac_of_peak=float(np.median(bd)) / peak if peak else np.nan,
        interior_median_density=float(np.median(inside)),
        boundary_below_burden=bool(np.median(bd) < BURDEN_THRESHOLD)))

    # ---- measurement 5, gradient along the outward normal ------------------
    ex_x, ex_y, _ = r["rings_um"][0]
    rx, ry = resample_ring(ex_x, ex_y, max(BOUNDARY_RESAMPLE_UM * 4, 40.0))
    if len(rx) < 8:
        continue
    tx = np.gradient(rx); ty = np.gradient(ry)
    ln = np.hypot(tx, ty); ln[ln == 0] = 1.0
    nx, ny = ty / ln, -tx / ln          # normal; sign fixed below
    cxg, cyg = np.mean(ex_x), np.mean(ex_y)
    flip = np.sign(((rx - cxg) * nx + (ry - cyg) * ny).sum()) or 1.0
    nx, ny = nx * flip, ny * flip       # now pointing outward
    offs = np.arange(-NORMAL_PROFILE_UM, NORMAL_PROFILE_UM + 1, NORMAL_STEP_UM)
    prof = np.array([np.median(sample_density(rx + nx * o, ry + ny * o, s))
                     for o in offs])
    for o, v in zip(offs, prof):
        prof_rows.append(dict(sample_id=r["sample_id"], index=r["index"],
                              classification=r["classification"],
                              offset_um=float(o), density=float(v),
                              fold=(float(v / s["bg"])
                                    if np.isfinite(s["bg"]) and s["bg"] > 0
                                    else np.nan)))

bnd = pd.DataFrame(bnd_rows)
prof = pd.DataFrame(prof_rows)
if len(bnd):
    write_csv(bnd, "111_boundary_measurements.csv")

    sub("Boundary density, by class")
    g = (bnd.groupby("classification")
         .agg(n=("index", "size"),
              density_median=("boundary_density_median", "median"),
              density_p25=("boundary_density_median", lambda v: q(v, 25)),
              density_p75=("boundary_density_median", lambda v: q(v, 75)),
              fold_median=("boundary_fold_median", "median"),
              frac_of_peak=("boundary_frac_of_peak", "median")))
    print(g.to_string())

    sub("Boundary density, by section")
    g2 = (bnd.groupby(["condition", "sample_id"])
          .agg(n=("index", "size"), background=("background", "first"),
               boundary_density=("boundary_density_median", "median"),
               boundary_fold=("boundary_fold_median", "median"),
               frac_of_peak=("boundary_frac_of_peak", "median"),
               peak_density=("enclosed_peak_density", "median")))
    print(g2.to_string())

    sub("MEASUREMENT 2: what FOCUS_HALF_MAX_FRACTION should be")
    fp = bnd["boundary_frac_of_peak"].to_numpy(float)
    fp = fp[np.isfinite(fp)]
    print(f"    script 04 assumes {FOCUS_HALF_MAX_FRACTION:.2f}")
    print(f"    measured median   {np.median(fp):.3f}   "
          f"IQR {q(fp, 25):.3f} to {q(fp, 75):.3f}   n={len(fp)}")
    if np.median(fp) < FOCUS_HALF_MAX_FRACTION - 0.05:
        print("    Traced boundaries sit LOWER on the peak than 0.50, so script")
        print("    04 stops growing too early and its regions are smaller than")
        print("    the drawn ones.")
    elif np.median(fp) > FOCUS_HALF_MAX_FRACTION + 0.05:
        print("    Traced boundaries sit HIGHER on the peak than 0.50, so script")
        print("    04 grows too far.")
    else:
        print("    Consistent with 0.50. The growth rule is not the problem.")

    sub("MEASUREMENT 1: what an absolute admission threshold should be")
    for cls in BOUNDARY_CLASSES:
        v = bnd.loc[bnd["classification"] == cls, "boundary_density_median"]
        if len(v):
            print(f"    {cls:<12} n={len(v):>3}  median {v.median():>8.0f}  "
                  f"p10 {q(v, 10):>8.0f}  p90 {q(v, 90):>8.0f}  /mm2")
    v = bnd.loc[bnd["classification"].isin(BOUNDARY_CLASSES),
                "boundary_density_median"]
    print(f"\n    A floor at the p10 of drawn boundaries would be "
          f"{q(v, 10):,.0f} cells/mm2.")
    print(f"    Burden threshold for comparison: {BURDEN_THRESHOLD:,.0f}.")
    nb = int(bnd["boundary_below_burden"].sum())
    print(f"    {nb} of {len(bnd)} drawn boundaries sit BELOW burden density,")
    print("    which is why 'never admit below burden' was not achievable.")
else:
    print("    No traced boundary regions available.")


# %% Cell 7 - measurement 3, maxima inside regions and the merge distance
# =============================================================================

banner("LOCAL MAXIMA INSIDE ONE DRAWN REGION")

print("Script 04 separates peaks with a maximum_filter of "
      f"{2 * SEP_PX + 1} px = {(2 * SEP_PX + 1) * GRID_UM:.0f} um")
print(f"(PEAK_SEPARATION_UM = {PEAK_SEPARATION_UM:.0f}). If one drawn complex")
print("contains several maxima, the separation is too small for confluent")
print("disease, and the spacing of those maxima says what it should be.")
print("Stamps ARE included here: only an enclosing location is needed.\n")

mx_rows = []
for i, r in ann.iterrows():
    if r["classification"] == "normal":
        continue
    s = sec[r["sample_id"]]
    m = masks[i]
    py, px = s["peak_rc"]
    if not len(py):
        continue
    ins = m[py, px]
    n_in = int(ins.sum())
    pv = s["peak_v"][ins]
    # foci from the CURRENT detector whose peak falls inside
    nf = 0
    for rec in s["foci_recs"]:
        if m[rec["peak_row"], rec["peak_col"]]:
            nf += 1
    rec = dict(sample_id=r["sample_id"], condition=r["condition"],
               index=r["index"], classification=r["classification"],
               is_stamp=bool(r["is_stamp"]),
               area_um2=m.sum() * GRID_UM * GRID_UM,
               n_local_maxima_inside=n_in,
               n_current_foci_inside=nf,
               max_peak_inside=float(pv.max()) if n_in else np.nan)
    if n_in >= 2:
        yy = py[ins] * GRID_UM; xx = px[ins] * GRID_UM
        dd = np.hypot(xx[:, None] - xx[None, :], yy[:, None] - yy[None, :])
        np.fill_diagonal(dd, np.inf)
        nn = dd.min(axis=1)
        rec["min_maxima_spacing_um"] = float(nn.min())
        # the separation that would collapse this region to a single peak
        rec["merge_distance_needed_um"] = float(nn.max())
    else:
        rec["min_maxima_spacing_um"] = np.nan
        rec["merge_distance_needed_um"] = np.nan
    mx_rows.append(rec)

mx = pd.DataFrame(mx_rows)
if len(mx):
    write_csv(mx, "112_maxima_inside_regions.csv")

    sub("Maxima per drawn region")
    g = (mx.groupby(["condition", "sample_id"])
         .agg(n_regions=("index", "size"),
              maxima_total=("n_local_maxima_inside", "sum"),
              maxima_median=("n_local_maxima_inside", "median"),
              maxima_max=("n_local_maxima_inside", "max"),
              current_foci_inside=("n_current_foci_inside", "sum")))
    print(g.to_string())

    over = mx.loc[mx["n_local_maxima_inside"] >= 2]
    sub("MEASUREMENT 3: what PEAK_SEPARATION_UM should be")
    print(f"    regions containing 2 or more maxima : {len(over)} of {len(mx)}")
    if len(over):
        md = over["merge_distance_needed_um"].to_numpy(float)
        md = md[np.isfinite(md)]
        print(f"    separation needed to collapse each to one peak:")
        print(f"      median {np.median(md):>7.0f} um   p75 {q(md, 75):>7.0f}   "
              f"p90 {q(md, 90):>7.0f}   max {md.max():>7.0f}")
        print(f"    script 04 currently uses {PEAK_SEPARATION_UM:.0f} um")
        rec_sep = q(md, 90)
        print(f"\n    A separation of about {rec_sep:.0f} um would collapse 90")
        print("    percent of the drawn regions to a single peak each.")
        print("    Regions still split at that value:")
        for _, rr in over.loc[over["merge_distance_needed_um"] > rec_sep].iterrows():
            print(f"      {rr['sample_id']} idx {rr['index']} "
                  f"({rr['classification']}) needs "
                  f"{rr['merge_distance_needed_um']:.0f} um, "
                  f"{rr['n_local_maxima_inside']} maxima")
    else:
        print("    No drawn region contains more than one local maximum, so")
        print("    over-segmentation is not a peak-separation problem.")


# %% Cell 8 - measurement 4, expert background
# =============================================================================

banner("EXPERT-DERIVED BACKGROUND")

print("Density inside regions the expert called 'normal'. This is the number")
print("the fold gate should be relative to. Every candidate statistic is scored")
print("against it, per section.\n")

bg_rows = []
for i, r in ann.iterrows():
    if r["classification"] != "normal":
        continue
    s = sec[r["sample_id"]]
    v = s["dens"][masks[i] & s["occ"]]
    v = v[np.isfinite(v)]
    if v.size < 10:
        continue
    bg_rows.append(dict(sample_id=r["sample_id"], condition=r["condition"],
                        index=r["index"], n_px=int(v.size),
                        area_um2=v.size * GRID_UM * GRID_UM,
                        density_median=float(np.median(v)),
                        density_p25=q(v, 25), density_p75=q(v, 75),
                        density_p90=q(v, 90)))
bgd = pd.DataFrame(bg_rows)

bg_cmp = []
for sid in SECTIONS:
    s = sec[sid]
    row = dict(sample_id=sid, condition=s["cond"])
    e = bgd.loc[bgd["sample_id"] == sid, "density_median"]
    row["expert_background"] = float(e.median()) if len(e) else np.nan
    row["n_normal_regions"] = int(len(e))
    for st in BG_STATS:
        row[f"bg_{st}"] = background_of(s["dens"], s["occ"], st)
    bg_cmp.append(row)
bgc = pd.DataFrame(bg_cmp)
# Guarded: a section with too few 'normal' pixels, or a genuinely empty one,
# gives an expert background at or near zero and every ratio blows up. Those
# sections are reported as unusable rather than allowed to dominate the score.
_eb = bgc["expert_background"].where(bgc["expert_background"] > 0)
for st in BG_STATS:
    bgc[f"err_{st}"] = (bgc[f"bg_{st}"] / _eb).replace([np.inf, -np.inf], np.nan)
bgc["expert_background_usable"] = _eb.notna()
if len(bgd):
    write_csv(bgd, "114_expert_background.csv")

sub("Expert background against each candidate statistic (cells/mm2)")
show = ["sample_id", "condition", "n_normal_regions", "expert_background"] + \
       [f"bg_{s}" for s in BG_STATS]
with pd.option_context("display.width", 260):
    print(bgc[show].to_string(index=False))

sub("MEASUREMENT 4: which statistic tracks the expert")
print("    ratio of statistic to expert background, 1.00 is perfect\n")
with pd.option_context("display.width", 260):
    print(bgc[["sample_id"] + [f"err_{s}" for s in BG_STATS]].to_string(index=False))
n_bad = int((~bgc["expert_background_usable"]).sum())
if n_bad:
    print(f"\n    {n_bad} section(s) have no usable expert background and are")
    print("    excluded from the scoring below. Draw more 'normal' regions on")
    print("    them if this statistic matters:")
    for _, rr in bgc.loc[~bgc["expert_background_usable"]].iterrows():
        print(f"      {rr['sample_id']}  ({int(rr['n_normal_regions'])} "
              f"normal regions)")

scores = []
for st in BG_STATS:
    e = bgc.loc[bgc["expert_background_usable"], f"err_{st}"].to_numpy(float)
    e = e[np.isfinite(e) & (e > 0)]        # log2 needs strictly positive
    if not e.size:
        continue
    scores.append(dict(statistic=st, n_sections=int(e.size),
                       median_ratio=float(np.median(e)),
                       max_abs_log2_error=float(np.max(np.abs(np.log2(e)))),
                       spread=float(e.max() / e.min())))
sc = (pd.DataFrame(scores).sort_values("max_abs_log2_error")
      if scores else pd.DataFrame())
print()
print(sc.to_string(index=False))
if len(sc):
    best = sc.iloc[0]
    print(f"\n    Best tracking of the expert: {best['statistic']}   "
          f"worst error {2 ** best['max_abs_log2_error']:.2f}x  "
          f"over {int(best['n_sections'])} sections")
    print(f"    Script 04 currently uses '{BACKGROUND_STAT}'.")


# %% Cell 9 - measurement 6, agreement with the current detector
# =============================================================================

banner("AGREEMENT WITH THE CURRENT DETECTOR")

print("Per-region IoU against script 04's foci regions, plus split and merge")
print("counts. Stamps are excluded: a stamped disc cannot score a meaningful")
print("overlap against anything.\n")

ag_rows = []
for i, r in ann.iterrows():
    if r["is_stamp"] or r["classification"] == "normal":
        continue
    s = sec[r["sample_id"]]
    m = masks[i]
    det = s["foci_labels"] > 0
    inter = int((m & det).sum())
    union = int((m | det).sum())
    ids = np.unique(s["foci_labels"][m])
    ids = ids[ids > 0]
    ag_rows.append(dict(
        sample_id=r["sample_id"], condition=r["condition"], index=r["index"],
        classification=r["classification"],
        area_um2=m.sum() * GRID_UM * GRID_UM,
        iou=inter / union if union else np.nan,
        pct_region_detected=100.0 * inter / m.sum() if m.sum() else np.nan,
        n_detected_foci_overlapping=int(len(ids))))
ag = pd.DataFrame(ag_rows)
if len(ag):
    write_csv(ag, "115_detector_agreement.csv")
    sub("Overlap by class")
    print(ag.groupby("classification")
          .agg(n=("index", "size"), median_iou=("iou", "median"),
               median_pct_detected=("pct_region_detected", "median"),
               regions_split=("n_detected_foci_overlapping",
                              lambda v: int((v >= 2).sum())),
               regions_missed=("pct_region_detected",
                               lambda v: int((v < 5).sum()))).to_string())
    sub("Overlap by section")
    print(ag.groupby(["condition", "sample_id"])
          .agg(n=("index", "size"), median_iou=("iou", "median"),
               median_pct_detected=("pct_region_detected", "median"),
               regions_split=("n_detected_foci_overlapping",
                              lambda v: int((v >= 2).sum())),
               regions_missed=("pct_region_detected",
                               lambda v: int((v < 5).sum()))).to_string())
    print("\n    'regions_split' counts drawn regions overlapping 2 or more")
    print("    detected foci, which is over-segmentation. 'regions_missed'")
    print("    counts drawn regions with under 5 percent detected coverage.")


# %% Cell 10 - composition
# =============================================================================

banner("REGION COMPOSITION")

comp_rows = []
for i, r in ann.iterrows():
    s = sec[r["sample_id"]]
    m = masks[i]
    inside = m[s["gy"], s["gx"]]
    n = int(inside.sum())
    rec = dict(sample_id=r["sample_id"], condition=r["condition"],
               index=r["index"], classification=r["classification"],
               is_stamp=bool(r["is_stamp"]), n_cells=n)
    if n:
        vc = pd.Series(s["pheno"][inside]).value_counts()
        for p in PHENOTYPE_ORDER:
            rec[f"frac_{p}"] = float(vc.get(p, 0)) / n
        rec["frac_myeloid"] = float(sum(vc.get(p, 0)
                                        for p in MYELOID_FOR_DETECTION)) / n
        rec["frac_endothelial"] = float(vc.get(ENDOTHELIAL, 0)) / n
    comp_rows.append(rec)
comp = pd.DataFrame(comp_rows)
if len(comp):
    write_csv(comp, "113_region_composition.csv")
    sub("Composition by class and arm")
    print(comp.groupby(["condition", "classification"])
          .agg(n=("index", "size"), cells=("n_cells", "median"),
               myeloid=("frac_myeloid", "median"),
               endothelial=("frac_endothelial", "median")).to_string())
    print("\n    Compare 'complex' against 'normal' within each arm. A gate that")
    print("    separates them on composition is arm-blind and defensible; one")
    print("    that does not is not worth adding.")


# %% Cell 11 - recommendations
# =============================================================================

banner("RECOMMENDED PARAMETERS, MEASURED NOT SWEPT")

rec_rows = []


def add(param, current, measured, basis):
    rec_rows.append(dict(parameter=param, current=current, measured=measured,
                         basis=basis))


if len(bnd):
    fp = bnd["boundary_frac_of_peak"].replace([np.inf, -np.inf], np.nan).dropna()
    add("FOCUS_HALF_MAX_FRACTION", FOCUS_HALF_MAX_FRACTION,
        round(float(fp.median()), 3),
        f"median boundary density as a fraction of enclosed peak, "
        f"n={len(fp)} traced regions")
    v = bnd.loc[bnd["classification"].isin(BOUNDARY_CLASSES),
                "boundary_density_median"]
    add("ABSOLUTE_FLOOR (cells/mm2)", 0.0, round(q(v, 10), 0),
        f"p10 of drawn boundary density, n={len(v)}")
    add("boundary fold over background", FOCUS_FOLD_OVER_BACKGROUND,
        round(float(bnd["boundary_fold_median"].median()), 2),
        "median fold at drawn boundaries")
if len(mx):
    over = mx.loc[mx["n_local_maxima_inside"] >= 2]
    if len(over):
        add("PEAK_SEPARATION_UM", PEAK_SEPARATION_UM,
            round(q(over["merge_distance_needed_um"], 90), 0),
            f"p90 of the separation needed to collapse each drawn region to "
            f"one peak, n={len(over)} multi-maxima regions")
if len(sc):
    add("BACKGROUND_STAT", BACKGROUND_STAT, sc.iloc[0]["statistic"],
        f"smallest worst-case error against expert 'normal' regions across "
        f"{int(bgc['n_normal_regions'].sum())} regions")

recs = pd.DataFrame(rec_rows)
if len(recs):
    write_csv(recs, "116_recommended_parameters.csv")
    print(recs.to_string(index=False))

sub("How to read this")
print("  These are measurements, not a fit. Each one comes from what the")
print("  expert drew, on script 04's own grid, over hundreds of pixels per")
print("  region. None of them was chosen to improve any arm contrast, and the")
print("  same value applies to both arms.")
print()
print("  Before adopting anything, the alignment check in Cell 5 has to be")
print("  clean and the split and missed counts in Cell 9 have to make sense")
print("  against what you saw on the canvases.")


# %% Cell 12 - figures
# =============================================================================

banner("FIGURES")

# ---- F90 boundary threshold -------------------------------------------------
if len(bnd):
    fig, axes = plt.subplots(1, 2, figsize=(28, 13))
    ax = axes[0]
    for k, cls in enumerate([c for c in BOUNDARY_CLASSES
                             if (bnd["classification"] == c).any()]):
        v = bnd.loc[bnd["classification"] == cls]
        xj = k + np.random.default_rng(0).normal(0, 0.06, len(v))
        ax.scatter(xj, v["boundary_density_median"], s=320,
                   color=CLASS_COLORS[cls], edgecolors="#000000",
                   linewidths=2, alpha=0.85)
        ax.hlines(v["boundary_density_median"].median(), k - 0.25, k + 0.25,
                  color="#000000", linewidth=5)
    ax.axhline(BURDEN_THRESHOLD, color="#000000", linestyle=":", linewidth=4)
    ax.text(ax.get_xlim()[1], BURDEN_THRESHOLD, "  burden", va="center",
            fontsize=FONT_SIZE_TICK - 6)
    for sid in SECTIONS:
        ax.axhline(sec[sid]["bg"] * FOCUS_FOLD_OVER_BACKGROUND,
                   color=CONDITION_COLORS.get(sec[sid]["cond"], "#999999"),
                   linestyle="--", linewidth=2, alpha=0.7)
    ax.set_yscale("log")
    ax.set_xticks(range(len([c for c in BOUNDARY_CLASSES
                             if (bnd["classification"] == c).any()])))
    ax.set_xticklabels([c for c in BOUNDARY_CLASSES
                        if (bnd["classification"] == c).any()])
    ax.set_ylabel("Density at the drawn boundary\n(cells / mm$^2$)")
    ax.set_title("Where the expert says a lesion stops")
    style_axes(ax); panel_letter(ax, "A")

    ax = axes[1]
    for cond in CONDITION_ORDER:
        v = bnd.loc[bnd["condition"] == cond]
        if not len(v):
            continue
        ax.scatter(v["enclosed_peak_density"], v["boundary_density_median"],
                   s=320, color=CONDITION_COLORS[cond], edgecolors="#000000",
                   linewidths=2, alpha=0.85, label=CONDITION_LABELS[cond])
    lo = float(bnd["boundary_density_median"].min()) * 0.7
    hi = float(bnd["enclosed_peak_density"].max()) * 1.3
    xs = np.array([lo, hi])
    for f, ls in [(FOCUS_HALF_MAX_FRACTION, "-"), (0.25, "--"), (0.10, ":")]:
        ax.plot(xs, f * xs, ls, color="#000000", linewidth=3,
                label=f"boundary = {f:.2f} x peak")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Enclosed peak density (cells / mm$^2$)")
    ax.set_ylabel("Boundary density (cells / mm$^2$)")
    ax.set_title("Testing the 50 percent growth rule")
    ax.legend(frameon=False, fontsize=FONT_SIZE_LEGEND - 10, loc="upper left")
    style_axes(ax); panel_letter(ax, "B")
    fig.suptitle("Measurement 1 and 2: the admission threshold and the growth "
                 "rule", fontsize=FONT_SIZE_TITLE, y=1.02)
    fig.tight_layout()
    save_fig(fig, "F90_boundary_threshold")

# ---- F92 maxima and merge distance -----------------------------------------
if len(mx):
    fig, axes = plt.subplots(1, 2, figsize=(28, 13))
    ax = axes[0]
    for cond in CONDITION_ORDER:
        v = mx.loc[mx["condition"] == cond, "n_local_maxima_inside"]
        if not len(v):
            continue
        bins = np.arange(-0.5, max(6, v.max() + 1.5))
        ax.hist(v, bins=bins, alpha=0.65, color=CONDITION_COLORS[cond],
                edgecolor="#000000", linewidth=2, label=CONDITION_LABELS[cond])
    ax.axvline(1, color="#000000", linestyle="--", linewidth=3)
    ax.set_xlabel("Local maxima inside one drawn region")
    ax.set_ylabel("Regions")
    ax.set_title("Anything above 1 is over-segmentation")
    ax.legend(frameon=False, fontsize=FONT_SIZE_LEGEND - 8)
    style_axes(ax); panel_letter(ax, "A")

    ax = axes[1]
    over = mx.loc[mx["n_local_maxima_inside"] >= 2]
    if len(over):
        md = np.sort(over["merge_distance_needed_um"].dropna().to_numpy())
        ax.step(md, np.arange(1, len(md) + 1) / len(md) * 100, where="post",
                linewidth=5, color="#D73027")
        ax.axvline(PEAK_SEPARATION_UM, color="#000000", linewidth=4)
        ax.text(PEAK_SEPARATION_UM, 50, f"  current {PEAK_SEPARATION_UM:.0f} um",
                fontsize=FONT_SIZE_TICK - 6)
        ax.axhline(90, color="#999999", linestyle=":", linewidth=3)
    ax.set_xlabel("Separation needed to collapse a region to one peak (um)")
    ax.set_ylabel("Percent of multi-maxima regions collapsed")
    ax.set_title("Measurement 3: PEAK_SEPARATION_UM")
    style_axes(ax); panel_letter(ax, "B")
    fig.suptitle("Over-segmentation, measured on drawn regions",
                 fontsize=FONT_SIZE_TITLE, y=1.02)
    fig.tight_layout()
    save_fig(fig, "F92_maxima_and_merge")

# ---- F93 expert background --------------------------------------------------
if len(bgc) and bgc["expert_background"].notna().any():
    fig, ax = plt.subplots(figsize=(22, 13))
    xs = np.arange(len(bgc))
    ax.scatter(xs, bgc["expert_background"], s=600, marker="*",
               color="#000000", zorder=5, label="expert 'normal' regions")
    mk = ["o", "s", "^", "v", "D", "P"]
    for j, st in enumerate(BG_STATS):
        ax.scatter(xs, bgc[f"bg_{st}"], s=260, marker=mk[j % len(mk)],
                   alpha=0.85, label=st)
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{short_label(s)}\n{bgc['condition'].iloc[i][:4]}"
                        for i, s in enumerate(bgc["sample_id"])])
    ax.set_yscale("log")
    ax.set_ylabel("Background density (cells / mm$^2$)")
    ax.set_title("Measurement 4: which background statistic tracks the expert",
                 fontsize=FONT_SIZE_TITLE)
    ax.legend(frameon=False, fontsize=FONT_SIZE_LEGEND - 10, ncol=2)
    style_axes(ax)
    save_fig(fig, "F93_expert_background")

# ---- F94 edge gradient ------------------------------------------------------
if len(prof):
    fig, ax = plt.subplots(figsize=(22, 13))
    for cls in BOUNDARY_CLASSES:
        p = prof.loc[prof["classification"] == cls]
        if not len(p):
            continue
        g = p.groupby("offset_um")["fold"].median()
        ax.plot(g.index, g.values, linewidth=6, color=CLASS_COLORS[cls],
                label=f"{cls} (n={p['index'].nunique()})")
    ax.axvline(0, color="#000000", linewidth=4)
    ax.text(0, ax.get_ylim()[1], " drawn boundary", va="top", ha="left",
            fontsize=FONT_SIZE_TICK - 6)
    ax.axhline(1.0, color="#999999", linestyle=":", linewidth=3)
    ax.set_xlabel("Distance from the drawn boundary (um), positive is outward")
    ax.set_ylabel("Density, fold over section background")
    ax.set_title("Measurement 5: is a fuzzy edge measurably shallower?",
                 fontsize=FONT_SIZE_TITLE)
    ax.legend(frameon=False, fontsize=FONT_SIZE_LEGEND - 8)
    style_axes(ax)
    save_fig(fig, "F94_edge_gradient")

# ---- F95 overlay maps -------------------------------------------------------
ncol = 3
nrow = int(np.ceil(len(SECTIONS) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(11 * ncol, 10 * nrow))
axes = np.atleast_1d(axes).ravel()
for ax, sid in zip(axes, SECTIONS):
    s = sec[sid]
    d = np.where(s["occ"], s["dens"], np.nan)
    ax.imshow(np.log10(np.clip(d, 1.0, None)), cmap="Greys",
              interpolation="nearest")
    det = s["foci_labels"] > 0
    if det.any():
        ax.contour(det.astype(float), levels=[0.5], colors=["#00A0C6"],
                   linewidths=2.5)
    for i, r in ann.loc[ann["sample_id"] == sid].iterrows():
        m = masks[i]
        if not m.any():
            continue
        ax.contour(m.astype(float), levels=[0.5],
                   colors=[CLASS_COLORS.get(r["classification"], "#000000")],
                   linewidths=3, linestyles="--" if r["is_stamp"] else "-")
    py, px = s["peak_rc"]
    ax.scatter(px, py, s=40, color="#FF00FF", marker="+", linewidths=1.5)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_edgecolor(CONDITION_COLORS.get(s["cond"], "#000000"))
        sp.set_linewidth(6)
    ax.set_title(f"{sid}  bg {s['bg']:.0f}/mm$^2$\n"
                 f"{len(s['foci_recs'])} detected, "
                 f"{int((ann['sample_id'] == sid).sum())} drawn",
                 fontsize=FONT_SIZE_BASE - 6)
for ax in axes[len(SECTIONS):]:
    ax.axis("off")
handles = [Line2D([0], [0], color="#00A0C6", lw=4, label="current detector foci"),
           Line2D([0], [0], color="#FF00FF", lw=0, marker="+", markersize=16,
                  label="local maxima"),
           Line2D([0], [0], color="#000000", lw=3, linestyle="--",
                  label="brush stamp, boundary not usable")]
handles += [Line2D([0], [0], color=CLASS_COLORS[c], lw=4, label=f"drawn: {c}")
            for c in CLASS_ORDER if (ann["classification"] == c).any()]
fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
           fontsize=FONT_SIZE_LEGEND - 8, bbox_to_anchor=(0.5, -0.03))
fig.suptitle("Drawn regions against the current detector",
             fontsize=FONT_SIZE_TITLE, y=1.0)
fig.tight_layout()
save_fig(fig, "F95_overlay_maps")


# %% Cell 13 - wrap up
# =============================================================================

banner("SUMMARY")

print(f"Sections            : {len(SECTIONS)}")
print(f"Annotations         : {len(ann)}  "
      f"({int((~ann['is_stamp']).sum())} traced, "
      f"{int(ann['is_stamp'].sum())} stamped)")
print(f"Boundary measurable : {len(bnd)}")

sub("Read in this order")
print("  1. The alignment check in Cell 5. Nothing is usable until it passes.")
print("  2. F95. The drawn regions against the detector, on every section.")
print("     Check it against what you saw on the canvases.")
print("  3. Table 116 and F90, F92, F93. The four measured parameters.")
print("  4. Cell 9. Split and missed counts, which is the honest error term")
print("     the count-based objective could never produce.")

sub("Not changed by this script")
print("  Nothing. Any adopted parameter becomes script 04 revision 5, after")
print("  which 06, 07, 08 and 09 all rerun and every manuscript number is")
print("  re-baselined.")

banner("DONE")
sys.stdout = _tee.terminal
_tee.close()
