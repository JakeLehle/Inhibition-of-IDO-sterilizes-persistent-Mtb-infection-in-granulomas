#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - DETECTED OBJECT COMPOSITION AND REVIEW CANVAS
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 04f of the AKOYA analysis series. READ-ONLY.

THE QUESTION THIS ANSWERS

    The candidate detector finds 14 objects on G3_43106 where the pathologist
    counted 1, and 4 on G3_43118 where the count was 2. The annotation score
    cannot adjudicate that, because only one 'complex' was drawn on 43106 and
    coverage of that one region says nothing about the other thirteen objects.

    Finding 1 rests on residual treated foci being few. A setting that
    multiplies the treated count changes that sentence as a side effect rather
    than as a decision, so this has to be settled before anything is adopted.

THE NEGATIVE CONTROL WE ALREADY HAD AND HAD NOT USED

    21 regions were drawn as 'normal', meaning uninvolved parenchyma. A detected
    object sitting inside one of those is a FALSE POSITIVE by the annotator's
    own labelling.

    Everywhere else in this project precision has been refused, correctly,
    because the annotations are not exhaustive and unannotated tissue is not an
    error. Inside a 'normal' region that objection does not apply: the annotator
    has positively asserted there is no lesion there. So precision IS measurable
    on that subset, and it is the sharpest number available for the treated
    over-detection question.

    This is reported as the headline check.

THE COMPOSITION TEST

    From script 04d, measured on the drawn regions:
        treated 'normal'   5.1% myeloid, 48.3% endothelial
        treated 'complex' 15.2% myeloid, 11.5% endothelial
    Those profiles are far apart, so each detected object can be scored against
    both and the answer read off directly. If the fourteen objects on 43106 look
    like treated 'normal', they are perivascular tissue and the level is too low
    in quiet lung. If they look like treated 'complex', the count of one was
    conservative and Deepak's original point extends to the treated arm.

    CAUTION, STATED IN THE OUTPUT TOO: the treated 'complex' reference rests on
    THREE drawn regions. It is thin. Both the own-arm and the pooled reference
    are reported so the dependence on that choice is visible, and neither is
    presented as a hard call.

THE REVIEW CANVAS

    Per section, the myeloid cells with every detected object outlined and
    numbered, keyed to a table giving each object's composition, size and
    density. This is for the eye, not for annotation: no transform is written
    and nothing here is meant to be re-imported into QuPath.

WHAT IS NOT DECIDED HERE

    Nothing is adopted. If the treated objects come back looking like normal
    parenchyma, the fix is an absolute floor under the relative level, and that
    is a further run. If they look like lesions, route B at 6.5x is adopted and
    Finding 1's wording changes to match.

INPUTS
    annotations/*_annotations.geojson
    structures_rev4/cell_assignments/*_cell_structures.csv

OUTPUTS
    object_review/tables/
        130_detected_object_composition, 131_normal_region_violations,
        132_reference_profiles, 133_treated_verdict,
        134_precision_on_normal_regions, 00_object_review_report.txt
    object_review/figures/
        F100_composition_space, F101_normal_violations,
        F102_review_<section> (one per section)

USAGE
    conda activate sc_pre
    python AKOYA_04f_Object_Review.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================

ANN_DIR = "/master/jlehle/WORKING/AKOYA/annotations"
IN_DIR = "/master/jlehle/WORKING/AKOYA/structures_rev4"
CELL_DIR = f"{IN_DIR}/cell_assignments"
OUT_DIR = "/master/jlehle/WORKING/AKOYA/object_review"

USE_AGG = True

# ---- shared with script 04 --------------------------------------------------
GRID_UM = 25.0
DENSITY_BANDWIDTH_UM = 75.0
EDGE_MIN_WEIGHT = 0.25
TISSUE_CLOSE_ITERATIONS = 2
TISSUE_FILL_HOLES = True
FOCUS_FILL_HOLES = True
BURDEN_THRESHOLD = 3000.0

# ---- the three settings under review ----------------------------------------
SETTINGS = {
    "baseline":   dict(route="A", bg_stat="median", fold=6.0, half_frac=0.50,
                       sep_um=200.0, min_area=10000.0, min_cells=40),
    "candidateA": dict(route="A", bg_stat="p25", fold=6.0, half_frac=0.25,
                       sep_um=200.0, min_area=10000.0, min_cells=40),
    "candidateB": dict(route="B", bg_stat="p25", level_fold=6.5, close_um=0.0,
                       min_area=10000.0, min_cells=40),
}
PRIMARY = "candidateB"          # the one under review; others for context

# ---- composition ------------------------------------------------------------
# An object is called a violation when this much of it sits inside a region the
# annotator called 'normal'. 50 percent means the majority of the object is in
# tissue positively asserted to be uninvolved.
VIOLATION_OVERLAP_PCT = 50.0
MIN_CELLS_FOR_COMPOSITION = 20

IDO1_POS = "CD68+IDO1+ Macrophages"
IDO1_NEG = "CD68+IDO1- Macrophages"
MYELOID_FOR_DETECTION = [IDO1_POS, IDO1_NEG, "CD163+ Macrophages", "Neutrophils"]
ENDOTHELIAL = "Endothelial cells"
PHENOTYPE_ORDER = [
    IDO1_POS, IDO1_NEG, "CD163+ Macrophages", "Neutrophils",
    "Helper T cells", "CD4- T cells", "Tregs", "B cells", "Plasma cells",
    "Endothelial cells", "Epithelial/Tumor cells", "Other",
]

LESION_CLASSES = ["complex", "clean_edge"]     # 'missed' excluded: faint by
                                               # definition, would blur the
                                               # lesion reference
NORMAL_CLASS = "normal"

# ---- review canvas ----------------------------------------------------------
CANVAS_UM_PER_PX = 2.0
MAX_POINTS_CANVAS = 200000
LABEL_MIN_AREA_UM2 = 0.0        # label every object; raise to declutter

# ---- appearance -------------------------------------------------------------
CONDITION_ORDER = ["D1MT", "Untreated"]
CONDITION_COLORS = {"D1MT": "#2C7FB8", "Untreated": "#D95F02"}
CONDITION_LABELS = {"D1MT": "D1MT-treated", "Untreated": "Untreated"}
CLASS_COLORS = {"complex": "#D73027", "missed": "#4575B4",
                "clean_edge": "#1A9850", "normal": "#999999"}
LESION_COLOR = "#D73027"
NORMALLIKE_COLOR = "#4575B4"
VIOLATION_COLOR = "#000000"
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
CELL_BG_COLOR = "#D9D9D9"
MYELOID_COLOR = "#404040"

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
import warnings
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

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
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
rng = np.random.default_rng(0)


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


_tee = Tee(os.path.join(TAB_DIR, "00_object_review_report.txt"))
sys.stdout = _tee

banner("AKOYA DETECTED OBJECT REVIEW (script 04f, read-only)")
print(f"Run time : {datetime.now().isoformat(timespec='seconds')}")
print(f"Output   : {OUT_DIR}")
print(f"Primary  : {PRIMARY}")
print("\nNothing is adopted. structures_rev4 is untouched.")


# %% Cell 3 - detection machinery, identical to scripts 04b, 04d and 04e
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
    counts = np.zeros(shape, dtype=float)
    np.add.at(counts, (gy[sel], gx[sel]), 1.0)
    sigma = DENSITY_BANDWIDTH_UM / GRID_UM
    sm = ndi.gaussian_filter(counts, sigma=sigma, mode="constant")
    w = ndi.gaussian_filter(occ.astype(float), sigma=sigma, mode="constant")
    return (sm / np.maximum(w, EDGE_MIN_WEIGHT)) / ((GRID_UM ** 2) / 1e6)


def background_of(dens, occ, stat):
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
    for alt in ["p25", "p50", "median", "p75"]:
        b = _s(alt)
        if np.isfinite(b) and b > 0:
            print(f"    WARNING: background '{stat}' is {bg}; escalated to "
                  f"'{alt}' = {b:.2f}")
            return b
    return np.nan


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


def route_a(dens, occ, gx, gy, bg, fold, half_frac, sep_um, min_area, min_cells):
    sep_px = max(1, int(round(sep_um / GRID_UM)))
    py, px, pv = find_local_maxima(dens, occ, sep_px)
    keep = pv >= fold * bg
    py, px, pv = py[keep], px[keep], pv[keep]
    if not len(pv):
        return np.zeros_like(dens, int)
    o = np.argsort(pv)[::-1]
    py, px, pv = py[o], px[o], pv[o]
    markers = np.zeros(dens.shape, int)
    markers[py, px] = np.arange(1, len(pv) + 1)
    territory = (dens >= half_frac * pv.min()) & occ
    if HAVE_SKIMAGE:
        part = _sk_watershed(-np.where(occ, dens, 0.0), markers, mask=territory)
    else:
        _, ind = ndi.distance_transform_edt(markers == 0, return_indices=True)
        part = np.where(territory, markers[ind[0], ind[1]], 0)
    labels = np.zeros(dens.shape, int); nid = 1
    for i, (yy, xx, vv) in enumerate(zip(py, px, pv), start=1):
        region = (part == i) & (dens >= half_frac * vv)
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
        if (assigned.sum() * GRID_UM * GRID_UM < min_area
                or int(assigned[gy, gx].sum()) < min_cells):
            continue
        labels[assigned] = nid
        nid += 1
    return labels


def route_b(dens, occ, gx, gy, bg, level_fold, close_um, min_area, min_cells):
    m = (dens >= level_fold * bg) & occ
    if close_um > 0:
        r = max(1, int(round(close_um / GRID_UM)))
        m = ndi.binary_closing(m, structure=ndi.generate_binary_structure(2, 1),
                               iterations=r)
        m &= occ
    if FOCUS_FILL_HOLES:
        m = ndi.binary_fill_holes(m)
    lab, n = ndi.label(m)
    labels = np.zeros(dens.shape, int); nid = 1
    for i in range(1, n + 1):
        comp = lab == i
        if (comp.sum() * GRID_UM * GRID_UM < min_area
                or int(comp[gy, gx].sum()) < min_cells):
            continue
        labels[comp] = nid
        nid += 1
    return labels


def detect(s, cfg):
    bg = s["bg"][cfg["bg_stat"]]
    if not np.isfinite(bg) or bg <= 0:
        return np.zeros_like(s["dens"], int), np.nan
    if cfg["route"] == "A":
        return route_a(s["dens"], s["occ"], s["gx"], s["gy"], bg, cfg["fold"],
                       cfg["half_frac"], cfg["sep_um"], cfg["min_area"],
                       cfg["min_cells"]), bg
    return route_b(s["dens"], s["occ"], s["gx"], s["gy"], bg, cfg["level_fold"],
                   cfg["close_um"], cfg["min_area"], cfg["min_cells"]), bg


# %% Cell 4 - load
# =============================================================================

banner("LOADING")

ann_files = sorted(glob.glob(os.path.join(ANN_DIR, "*_annotations.geojson")))
SECTIONS = [os.path.basename(f).replace("_annotations.geojson", "")
            for f in ann_files]

sec = {}
for sid in SECTIONS:
    p = os.path.join(CELL_DIR, f"{sid}_cell_structures.csv")
    if not os.path.exists(p):
        continue
    d = pd.read_csv(p, low_memory=False)
    d = d.loc[np.isfinite(d["x"]) & np.isfinite(d["y"])].reset_index(drop=True)
    x = d["x"].to_numpy(float); y = d["y"].to_numpy(float)
    gx0, gy0, gx, gy, shape = make_grid(x, y, GRID_UM)
    occ = tissue_mask(gx, gy, shape)
    mye = d["pheno"].isin(MYELOID_FOR_DETECTION).to_numpy()
    dens = density_map(gx, gy, shape, mye, occ)
    sec[sid] = dict(x=x, y=y, gx=gx, gy=gy, gx0=gx0, gy0=gy0, shape=shape,
                    occ=occ, dens=dens, pheno=d["pheno"].to_numpy(),
                    is_mye=mye, n_cells=len(d),
                    cond=d["condition"].iloc[0] if "condition" in d.columns else "?",
                    bg={st: background_of(dens, occ, st)
                        for st in ["median", "p25", "p10", "p75"]})
    print(f"    {sid:<14} {len(d):>9,} cells   bg median "
          f"{sec[sid]['bg']['median']:>8.1f}   p25 {sec[sid]['bg']['p25']:>8.1f}")

SECTIONS = [s for s in SECTIONS if s in sec]
SECTIONS.sort(key=lambda s: (CONDITION_ORDER.index(sec[s]["cond"])
                             if sec[s]["cond"] in CONDITION_ORDER else 9, s))


def geom_to_rings(geom):
    t = geom["type"]; c = geom["coordinates"]
    polys = [c] if t == "Polygon" else (c if t == "MultiPolygon" else [])
    return [(np.asarray(p[0], float), [np.asarray(h, float) for h in p[1:]])
            for p in polys if p]


def um_to_grid(xu, yu, s):
    return ((xu - s["gx0"]) / GRID_UM) + 1.0, ((yu - s["gy0"]) / GRID_UM) + 1.0


def rasterise(rings_um, s):
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


ann_rows = []
for f in ann_files:
    sid = os.path.basename(f).replace("_annotations.geojson", "")
    if sid not in sec:
        continue
    gj = json.load(open(f))
    pr = gj["properties"]
    X0, Y0, SC = pr["x0_um"], pr["y0_um"], pr["um_per_px"]
    for feat in gj["features"]:
        p = feat["properties"]
        rings = geom_to_rings(feat["geometry"])
        if not rings:
            continue
        rings_um = [((r[:, 0] + 0.5) * SC + X0, (r[:, 1] + 0.5) * SC + Y0,
                     [((h[:, 0] + 0.5) * SC + X0, (h[:, 1] + 0.5) * SC + Y0)
                      for h in hs]) for r, hs in rings]
        ann_rows.append(dict(sample_id=sid, condition=sec[sid]["cond"],
                             index=p.get("index"),
                             classification=p.get("classification", "unclassified"),
                             area_px2=p.get("area_px2", np.nan),
                             rings_um=rings_um))
ann = pd.DataFrame(ann_rows)
ann["is_stamp"] = False
for sid, g in ann.groupby("sample_id"):
    a = g["area_px2"].round(1)
    dup = a.value_counts()
    ann.loc[g.index[a.isin(dup[dup > 1].index)], "is_stamp"] = True
ann_masks = {i: rasterise(r["rings_um"], sec[r["sample_id"]])
             for i, r in ann.iterrows()}

# union of 'normal' masks per section: the negative control
normal_mask = {}
for sid in SECTIONS:
    idx = ann.index[(ann["sample_id"] == sid)
                    & (ann["classification"] == NORMAL_CLASS)]
    m = np.zeros(sec[sid]["shape"], bool)
    for i in idx:
        m |= ann_masks[i]
    normal_mask[sid] = m
    print(f"    {sid:<14} 'normal' negative-control area "
          f"{m.sum() * GRID_UM ** 2 / 1e6:>6.2f} mm2 "
          f"({100.0 * m.sum() / max(1, sec[sid]['occ'].sum()):>4.1f}% of tissue)")


def composition(mask, s):
    inside = mask[s["gy"], s["gx"]]
    n = int(inside.sum())
    out = {"n_cells": n}
    if n:
        vc = pd.Series(s["pheno"][inside]).value_counts()
        for p in PHENOTYPE_ORDER:
            out[f"frac_{p}"] = float(vc.get(p, 0)) / n
        out["frac_myeloid"] = float(sum(vc.get(p, 0)
                                        for p in MYELOID_FOR_DETECTION)) / n
        out["frac_endothelial"] = float(vc.get(ENDOTHELIAL, 0)) / n
    else:
        for p in PHENOTYPE_ORDER:
            out[f"frac_{p}"] = np.nan
        out["frac_myeloid"] = np.nan
        out["frac_endothelial"] = np.nan
    return out


# %% Cell 5 - reference profiles from the drawn regions
# =============================================================================

banner("REFERENCE PROFILES")

frac_cols = [f"frac_{p}" for p in PHENOTYPE_ORDER]
ref_rows = []
for i, r in ann.iterrows():
    if r["classification"] not in LESION_CLASSES + [NORMAL_CLASS]:
        continue
    c = composition(ann_masks[i], sec[r["sample_id"]])
    if c["n_cells"] < MIN_CELLS_FOR_COMPOSITION:
        continue
    kind = "lesion" if r["classification"] in LESION_CLASSES else "normal"
    ref_rows.append(dict(sample_id=r["sample_id"], condition=r["condition"],
                         index=r["index"], classification=r["classification"],
                         kind=kind, **c))
refs = pd.DataFrame(ref_rows)

profiles = {}
prof_rows = []
for cond in CONDITION_ORDER + ["pooled"]:
    for kind in ["lesion", "normal"]:
        sel = refs.loc[refs["kind"] == kind]
        if cond != "pooled":
            sel = sel.loc[sel["condition"] == cond]
        if not len(sel):
            continue
        v = sel[frac_cols].mean(axis=0).to_numpy(float)
        profiles[(cond, kind)] = np.nan_to_num(v)
        prof_rows.append(dict(reference=cond, kind=kind, n_regions=len(sel),
                              n_cells_total=int(sel["n_cells"].sum()),
                              frac_myeloid=float(sel["frac_myeloid"].mean()),
                              frac_endothelial=float(sel["frac_endothelial"].mean()),
                              **{c: float(sel[c].mean()) for c in frac_cols}))
prof = pd.DataFrame(prof_rows)
write_csv(prof, "132_reference_profiles.csv")

show = ["reference", "kind", "n_regions", "n_cells_total", "frac_myeloid",
        "frac_endothelial"]
print(prof[show].to_string(index=False))

thin = prof.loc[(prof["kind"] == "lesion") & (prof["n_regions"] < 5)]
if len(thin):
    sub("CAUTION")
    for _, r in thin.iterrows():
        print(f"    The '{r['reference']}' lesion reference rests on "
              f"{int(r['n_regions'])} drawn region(s).")
    print("    That is thin. Both the own-arm and the pooled reference are")
    print("    carried through below so the dependence on this choice is")
    print("    visible in the output rather than buried in a centroid.")


def cosine(v, c):
    v = np.asarray(v, float); c = np.asarray(c, float)
    ok = np.isfinite(v) & np.isfinite(c)
    if not ok.any():
        return np.nan
    v, c = v[ok], c[ok]
    nv, nc = np.linalg.norm(v), np.linalg.norm(c)
    return float(v @ c / (nv * nc)) if nv > 0 and nc > 0 else np.nan


# %% Cell 6 - every detected object, for every setting
# =============================================================================

banner("DETECTED OBJECTS")

obj_rows = []
labels_of = {}
for name, cfg in SETTINGS.items():
    for sid in SECTIONS:
        s = sec[sid]
        lab, bg = detect(s, cfg)
        labels_of[(name, sid)] = lab
        n = int(lab.max())
        for oid in range(1, n + 1):
            m = lab == oid
            c = composition(m, s)
            dv = s["dens"][m]
            nm = normal_mask[sid]
            pct_in_normal = 100.0 * (m & nm).sum() / max(1, m.sum())
            # which drawn classes does it touch
            touch = set()
            for i in ann.index[ann["sample_id"] == sid]:
                if (m & ann_masks[i]).any():
                    touch.add(ann.loc[i, "classification"])
            rec = dict(setting=name, sample_id=sid, condition=s["cond"],
                       object_id=oid,
                       area_um2=float(m.sum()) * GRID_UM * GRID_UM,
                       peak_density=float(dv.max()) if dv.size else np.nan,
                       median_density=float(np.median(dv)) if dv.size else np.nan,
                       fold_over_bg=(float(dv.max() / bg)
                                     if dv.size and bg > 0 else np.nan),
                       pct_in_normal_regions=pct_in_normal,
                       touches_lesion=bool(touch & set(LESION_CLASSES)),
                       touches_missed=bool("missed" in touch),
                       touches_normal=bool(NORMAL_CLASS in touch),
                       touches_nothing=len(touch) == 0, **c)
            v = np.array([rec.get(cc, np.nan) for cc in frac_cols], float)
            for refname in [s["cond"], "pooled"]:
                pl = profiles.get((refname, "lesion"))
                pn = profiles.get((refname, "normal"))
                tag = "own_arm" if refname == s["cond"] else "pooled"
                rec[f"cos_lesion_{tag}"] = cosine(v, pl) if pl is not None else np.nan
                rec[f"cos_normal_{tag}"] = cosine(v, pn) if pn is not None else np.nan
                if pl is not None and pn is not None:
                    rec[f"call_{tag}"] = ("lesion"
                                          if rec[f"cos_lesion_{tag}"]
                                          >= rec[f"cos_normal_{tag}"] else "normal")
                else:
                    rec[f"call_{tag}"] = "undetermined"
            obj_rows.append(rec)
    print(f"    {name:<12} " + "  ".join(
        f"{short_label(q)}={int(labels_of[(name, q)].max())}" for q in SECTIONS))

objs = pd.DataFrame(obj_rows)
objs.loc[objs["n_cells"] < MIN_CELLS_FOR_COMPOSITION,
         ["call_own_arm", "call_pooled"]] = "too_few_cells"
write_csv(objs, "130_detected_object_composition.csv")


# %% Cell 7 - the negative control
# =============================================================================

banner("NEGATIVE CONTROL: DETECTIONS INSIDE REGIONS DRAWN AS 'normal'")

print("Everywhere else precision is refused, because unannotated tissue is not")
print("an error. Inside a 'normal' region that objection does not apply: the")
print("annotator has positively asserted there is no lesion there. So an object")
print("sitting mostly inside one is a false positive, and precision IS")
print("measurable on this subset.\n")
print(f"An object counts as a violation when at least "
      f"{VIOLATION_OVERLAP_PCT:.0f} percent of it lies")
print("inside a drawn 'normal' region.\n")

objs["is_violation"] = objs["pct_in_normal_regions"] >= VIOLATION_OVERLAP_PCT
viol = objs.loc[objs["is_violation"]]
write_csv(viol, "131_normal_region_violations.csv")

rows = []
for name in SETTINGS:
    o = objs.loc[objs["setting"] == name]
    for cond in CONDITION_ORDER:
        oc = o.loc[o["condition"] == cond]
        if not len(oc):
            continue
        rows.append(dict(setting=name, condition=cond, n_objects=len(oc),
                         n_violations=int(oc["is_violation"].sum()),
                         pct_violations=100.0 * oc["is_violation"].mean(),
                         violation_area_mm2=float(
                             oc.loc[oc["is_violation"], "area_um2"].sum()) / 1e6))
    rows.append(dict(setting=name, condition="ALL", n_objects=len(o),
                     n_violations=int(o["is_violation"].sum()),
                     pct_violations=100.0 * o["is_violation"].mean(),
                     violation_area_mm2=float(
                         o.loc[o["is_violation"], "area_um2"].sum()) / 1e6))
prec = pd.DataFrame(rows)
write_csv(prec, "134_precision_on_normal_regions.csv")
print(prec.to_string(index=False))

sub("Fraction of drawn 'normal' AREA covered by detections")
cov_rows = []
for name in SETTINGS:
    for sid in SECTIONS:
        nm = normal_mask[sid]
        if not nm.any():
            continue
        lab = labels_of[(name, sid)]
        cov_rows.append(dict(setting=name, sample_id=sid,
                             condition=sec[sid]["cond"],
                             normal_area_mm2=nm.sum() * GRID_UM ** 2 / 1e6,
                             pct_of_normal_detected=100.0 * ((lab > 0) & nm).sum()
                             / nm.sum()))
cov = pd.DataFrame(cov_rows)
print(cov.pivot_table(index=["condition", "sample_id"], columns="setting",
                      values="pct_of_normal_detected").to_string())
print("\n    This is the cleanest single number in the whole comparison. It is")
print("    the percentage of tissue the annotator called uninvolved that the")
print("    detector calls lesion. Lower is unambiguously better and there is no")
print("    trade against it.")


# %% Cell 8 - the treated verdict
# =============================================================================

banner("TREATED-ARM VERDICT")

print("The question: are the extra objects on the treated sections lesions or")
print("perivascular parenchyma?\n")

tre = objs.loc[(objs["setting"] == PRIMARY) & (objs["condition"] == "D1MT")]
ver_rows = []
for sid in [q for q in SECTIONS if sec[q]["cond"] == "D1MT"]:
    o = tre.loc[tre["sample_id"] == sid]
    ver_rows.append(dict(
        sample_id=sid, n_objects=len(o),
        called_lesion_own_arm=int((o["call_own_arm"] == "lesion").sum()),
        called_normal_own_arm=int((o["call_own_arm"] == "normal").sum()),
        called_lesion_pooled=int((o["call_pooled"] == "lesion").sum()),
        called_normal_pooled=int((o["call_pooled"] == "normal").sum()),
        n_violations=int(o["is_violation"].sum()),
        touching_a_drawn_lesion=int(o["touches_lesion"].sum()),
        touching_nothing_drawn=int(o["touches_nothing"].sum()),
        median_myeloid=float(o["frac_myeloid"].median()),
        median_endothelial=float(o["frac_endothelial"].median())))
ver = pd.DataFrame(ver_rows)
write_csv(ver, "133_treated_verdict.csv")
with pd.option_context("display.width", 260):
    print(ver.to_string(index=False))

sub("Against the drawn reference profiles for the same arm")
for kind in ["lesion", "normal"]:
    p = prof.loc[(prof["reference"] == "D1MT") & (prof["kind"] == kind)]
    if len(p):
        print(f"    drawn treated '{kind}':   myeloid "
              f"{p['frac_myeloid'].iloc[0]:.3f}   endothelial "
              f"{p['frac_endothelial'].iloc[0]:.3f}   "
              f"(n={int(p['n_regions'].iloc[0])} regions)")
if len(tre):
    print(f"    detected treated objects: myeloid "
          f"{tre['frac_myeloid'].median():.3f}   endothelial "
          f"{tre['frac_endothelial'].median():.3f}   (n={len(tre)} objects)")

sub("How to read this")
print("  If the detected treated objects sit near the drawn 'normal' profile,")
print("  they are perivascular tissue, the level is too low in quiet lung, and")
print("  the fix is an absolute floor underneath the relative level.")
print()
print("  If they sit near the drawn 'lesion' profile, the pathologist's count")
print("  of one on 43106 was conservative, the detector is right, and Finding 1")
print("  needs rewording rather than the detector needing changing.")
print()
print("  If they are split between the two, the composition call is not")
print("  decisive on its own and F102 is where it gets settled by eye.")


# %% Cell 8b - can a composition filter fix the treated over-detection?
# =============================================================================

banner("COMPOSITION FILTER SWEEP")

print("A myeloid floor and an endothelial ceiling applied AFTER detection,")
print("identically to both arms. The negative control makes this scoreable on")
print("two quantities that both come from the annotations, so no expert-count")
print("proxy is needed:")
print()
print("    RECALL   median coverage of drawn lesion regions      higher better")
print("    HARM     percent of drawn 'normal' AREA called lesion  lower better")
print()
print("A filter that cuts the harm without costing recall is a real fix. One")
print("that moves both together is just a stricter threshold wearing a")
print("different name.")
print()
print("Arm-blind in FORM is not arm-blind in EFFECT. Treated lesions sit at a")
print("far lower myeloid fraction than untreated ones, so a floor tuned on")
print("untreated tissue preferentially deletes treated objects. The per-arm")
print("columns exist to make that visible rather than to hide it.\n")

MYELOID_FLOOR_GRID = [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
ENDOTHELIAL_CEILING_GRID = [1.00, 0.40, 0.30, 0.25, 0.20, 0.15]

_prim = objs.loc[objs["setting"] == PRIMARY]
_targets = ann.loc[ann["classification"].isin(LESION_CLASSES + ["missed"])
                   & ~ann["is_stamp"]]


def _apply_filter(mf, ec):
    keep, rec = {}, {}
    for sid in SECTIONS:
        o = _prim.loc[_prim["sample_id"] == sid]
        k = o.loc[(o["frac_myeloid"] >= mf)
                  & (o["frac_endothelial"] <= ec), "object_id"]
        keep[sid] = set(int(v) for v in k)
        rec[sid] = len(keep[sid])
    # recall on drawn lesion regions
    cov = []
    for i, r in _targets.iterrows():
        sid = r["sample_id"]
        lab = labels_of[(PRIMARY, sid)]
        det = np.isin(lab, list(keep[sid])) if keep[sid] else np.zeros_like(lab, bool)
        m = ann_masks[i]
        cov.append(100.0 * int((m & det).sum()) / m.sum() if m.sum() else np.nan)
    # harm on drawn 'normal' regions
    harm_num = harm_den = 0
    harm_arm = {}
    for sid in SECTIONS:
        nm = normal_mask[sid]
        if not nm.any():
            continue
        lab = labels_of[(PRIMARY, sid)]
        det = np.isin(lab, list(keep[sid])) if keep[sid] else np.zeros_like(lab, bool)
        a = int((det & nm).sum()); b = int(nm.sum())
        harm_num += a; harm_den += b
        arm = sec[sid]["cond"]
        u, w = harm_arm.get(arm, (0, 0))
        harm_arm[arm] = (u + a, w + b)
    return dict(
        myeloid_floor=mf, endothelial_ceiling=ec,
        n_objects=sum(rec.values()),
        n_treated=sum(rec[q] for q in SECTIONS if sec[q]["cond"] == "D1MT"),
        n_untreated=sum(rec[q] for q in SECTIONS if sec[q]["cond"] == "Untreated"),
        recall_pct=float(np.nanmedian(cov)),
        n_missed=int(np.nansum(np.asarray(cov) < 5.0)),
        harm_pct=100.0 * harm_num / harm_den if harm_den else np.nan,
        harm_pct_D1MT=(100.0 * harm_arm["D1MT"][0] / harm_arm["D1MT"][1]
                       if harm_arm.get("D1MT", (0, 0))[1] else np.nan),
        harm_pct_Untreated=(100.0 * harm_arm["Untreated"][0]
                            / harm_arm["Untreated"][1]
                            if harm_arm.get("Untreated", (0, 0))[1] else np.nan),
        **{f"n_{short_label(q)}": rec[q] for q in SECTIONS})


filt = pd.DataFrame([_apply_filter(mf, ec)
                     for mf in MYELOID_FLOOR_GRID
                     for ec in ENDOTHELIAL_CEILING_GRID])
write_csv(filt, "135_composition_filter_sweep.csv")

_cols = (["myeloid_floor", "endothelial_ceiling"]
         + [f"n_{short_label(q)}" for q in SECTIONS]
         + ["recall_pct", "n_missed", "harm_pct", "harm_pct_D1MT",
            "harm_pct_Untreated"])

sub("Endothelial ceiling alone (myeloid floor 0)")
with pd.option_context("display.width", 260):
    print(filt.loc[filt["myeloid_floor"] == 0.0]
          .sort_values("endothelial_ceiling", ascending=False)[_cols]
          .to_string(index=False))

sub("Myeloid floor alone (endothelial ceiling 1)")
with pd.option_context("display.width", 260):
    print(filt.loc[filt["endothelial_ceiling"] == 1.0]
          .sort_values("myeloid_floor")[_cols].to_string(index=False))

base = filt.loc[(filt["myeloid_floor"] == 0.0)
                & (filt["endothelial_ceiling"] == 1.0)].iloc[0]
sub("VERDICT ON THE FILTER")
print(f"    unfiltered: recall {base['recall_pct']:.1f}%   harm "
      f"{base['harm_pct']:.1f}%   treated objects {int(base['n_treated'])}")
gain = filt.loc[(filt["recall_pct"] >= base["recall_pct"] - 5.0)
                & (filt["harm_pct"] <= base["harm_pct"] * 0.5)]
if len(gain):
    print("\n    Filters that halve the harm for under 5 points of recall:")
    with pd.option_context("display.width", 260):
        print(gain.sort_values("harm_pct")[_cols].head(8).to_string(index=False))
    print("\n    If one of these also brings the treated counts near the expert")
    print("    counts, it is the fix. Check the per-arm harm columns: a filter")
    print("    that only reduces harm in the treated arm is doing something")
    print("    arm-specific in effect and should be rejected.")
else:
    print("\n    NO filter halves the harm while keeping recall within 5 points.")
    print("    That is an important negative result. Composition cannot separate")
    print("    the extra treated objects from real lesions, because treated")
    print("    lesions are themselves only weakly myeloid-enriched (0.152")
    print("    against 0.465 untreated). So composition is the wrong lever in")
    print("    the treated arm, and the review canvases become the evidence.")
    print("    Do NOT adopt a filter that produces the expected count by")
    print("    deleting real treated lesions.")


# %% Cell 9 - review canvases
# =============================================================================

banner("REVIEW CANVASES")

print(f"One per section for {PRIMARY}: myeloid cells with every detected object")
print("outlined and numbered, keyed to table 130. For the eye only, no")
print("transform is written and these are not meant for QuPath.\n")

for sid in SECTIONS:
    s = sec[sid]
    lab = labels_of[(PRIMARY, sid)]
    o = objs.loc[(objs["setting"] == PRIMARY) & (objs["sample_id"] == sid)]
    n = int(lab.max())
    xs, ys = s["x"], s["y"]
    w_um = xs.max() - xs.min(); h_um = ys.max() - ys.min()
    fig, ax = plt.subplots(figsize=(min(40, w_um / 400.0) + 2,
                                    min(40, h_um / 400.0) + 3))
    idx = (rng.choice(len(xs), MAX_POINTS_CANVAS, replace=False)
           if len(xs) > MAX_POINTS_CANVAS else np.arange(len(xs)))
    ax.scatter(xs[idx], ys[idx], s=0.5, color=CELL_BG_COLOR, linewidths=0,
               rasterized=True)
    mm = s["is_mye"][idx]
    ax.scatter(xs[idx][mm], ys[idx][mm], s=1.6, color=MYELOID_COLOR,
               linewidths=0, rasterized=True)

    ext = [s["gx0"] - GRID_UM, s["gx0"] + (s["shape"][1] - 1) * GRID_UM,
           s["gy0"] + (s["shape"][0] - 1) * GRID_UM, s["gy0"] - GRID_UM]
    # drawn regions first, underneath
    for i in ann.index[ann["sample_id"] == sid]:
        cls = ann.loc[i, "classification"]
        ax.contour(ann_masks[i].astype(float), levels=[0.5], extent=ext,
                   colors=[CLASS_COLORS.get(cls, "#000000")], linewidths=2.5,
                   linestyles="--" if ann.loc[i, "is_stamp"] else "-",
                   origin="upper")
    for oid in range(1, n + 1):
        m = lab == oid
        r = o.loc[o["object_id"] == oid]
        if not len(r):
            continue
        r = r.iloc[0]
        col = (VIOLATION_COLOR if r["is_violation"]
               else (LESION_COLOR if r["call_own_arm"] == "lesion"
                     else NORMALLIKE_COLOR))
        ax.contour(m.astype(float), levels=[0.5], extent=ext, colors=[col],
                   linewidths=4, origin="upper")
        if r["area_um2"] >= LABEL_MIN_AREA_UM2:
            yy, xx = np.nonzero(m)
            ax.text(s["gx0"] + xx.mean() * GRID_UM,
                    s["gy0"] + yy.mean() * GRID_UM, str(oid),
                    fontsize=16, fontweight="bold", color=col,
                    ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.15", fc="#FFFFFF",
                              ec=col, lw=1.5, alpha=0.85))
    ax.set_aspect("equal"); ax.invert_yaxis()
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_edgecolor(CONDITION_COLORS.get(s["cond"], "#000000"))
        sp.set_linewidth(6)
    nv = int(o["is_violation"].sum())
    nl = int((o["call_own_arm"] == "lesion").sum())
    ax.set_title(f"{sid}  {CONDITION_LABELS.get(s['cond'], '')}  [{PRIMARY}]\n"
                 f"{n} objects: {nl} lesion-like, {n - nl - nv} normal-like, "
                 f"{nv} inside a drawn 'normal' region\n"
                 f"numbers key to table 130",
                 fontsize=18)
    handles = [Line2D([0], [0], color=LESION_COLOR, lw=4, label="detected, lesion-like"),
               Line2D([0], [0], color=NORMALLIKE_COLOR, lw=4, label="detected, normal-like"),
               Line2D([0], [0], color=VIOLATION_COLOR, lw=4,
                      label="detected inside a drawn 'normal' region")]
    handles += [Line2D([0], [0], color=CLASS_COLORS[c], lw=2.5, label=f"drawn: {c}")
                for c in ["complex", "missed", "clean_edge", "normal"]
                if ((ann["sample_id"] == sid)
                    & (ann["classification"] == c)).any()]
    ax.legend(handles=handles, loc="upper right", frameon=True, fontsize=13)
    save_fig(fig, f"F102_review_{sid}")


# %% Cell 10 - figures
# =============================================================================

banner("FIGURES")

# ---- F100 composition space -------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(28, 13))
for ax, cond in zip(axes, CONDITION_ORDER):
    o = objs.loc[(objs["setting"] == PRIMARY) & (objs["condition"] == cond)]
    r = refs.loc[refs["condition"] == cond]
    ax.scatter(o["frac_myeloid"], o["frac_endothelial"], s=240, alpha=0.75,
               color="#BBBBBB", edgecolors="#000000", linewidths=1.5,
               label="detected objects")
    v = o.loc[o["is_violation"]]
    if len(v):
        ax.scatter(v["frac_myeloid"], v["frac_endothelial"], s=340,
                   facecolors="none", edgecolors=VIOLATION_COLOR, linewidths=4,
                   label="inside a drawn 'normal' region")
    for kind, col, mk in [("lesion", LESION_COLOR, "^"),
                          ("normal", NORMALLIKE_COLOR, "s")]:
        rr = r.loc[r["kind"] == kind]
        if len(rr):
            ax.scatter(rr["frac_myeloid"], rr["frac_endothelial"], s=380,
                       marker=mk, color=col, edgecolors="#000000",
                       linewidths=2, label=f"drawn {kind} (n={len(rr)})")
        p = profiles.get((cond, kind))
        if p is not None:
            mi = PHENOTYPE_ORDER.index
            my = sum(p[mi(q)] for q in MYELOID_FOR_DETECTION)
            en = p[mi(ENDOTHELIAL)]
            ax.scatter([my], [en], s=900, marker="*", color=col,
                       edgecolors="#000000", linewidths=2.5, zorder=6,
                       label=f"{kind} centroid")
    ax.set_xlabel("Myeloid fraction")
    ax.set_ylabel("Endothelial fraction")
    ax.set_title(CONDITION_LABELS[cond], color=CONDITION_COLORS[cond])
    ax.legend(frameon=False, fontsize=FONT_SIZE_LEGEND - 12, loc="upper right")
    style_axes(ax)
panel_letter(axes[0], "A"); panel_letter(axes[1], "B")
fig.suptitle(f"Detected objects against the drawn reference profiles "
             f"[{PRIMARY}]", fontsize=FONT_SIZE_TITLE, y=1.02)
fig.tight_layout()
save_fig(fig, "F100_composition_space")

# ---- F101 negative control --------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(28, 13))
ax = axes[0]
xs = np.arange(len(SETTINGS))
for k, cond in enumerate(CONDITION_ORDER):
    vals = [prec.loc[(prec["setting"] == nm) & (prec["condition"] == cond),
                     "pct_violations"].squeeze() if
            len(prec.loc[(prec["setting"] == nm) & (prec["condition"] == cond)])
            else 0.0 for nm in SETTINGS]
    vals = [float(v) if np.isscalar(v) or np.size(v) else 0.0 for v in vals]
    ax.bar(xs + (k - 0.5) * 0.35, vals, width=0.35,
           color=CONDITION_COLORS[cond], edgecolor="#000000", linewidth=2,
           label=CONDITION_LABELS[cond])
ax.set_xticks(xs); ax.set_xticklabels(list(SETTINGS), rotation=15)
ax.set_ylabel("Objects mostly inside a drawn\n'normal' region (%)")
ax.set_title("False positives, where they are measurable")
ax.legend(frameon=False, fontsize=FONT_SIZE_LEGEND - 8)
style_axes(ax); panel_letter(ax, "A")

ax = axes[1]
piv = cov.pivot_table(index="sample_id", columns="setting",
                      values="pct_of_normal_detected")
xs = np.arange(len(piv))
for k, nm in enumerate(SETTINGS):
    if nm not in piv.columns:
        continue
    ax.bar(xs + (k - 1) * 0.26, piv[nm], width=0.26, edgecolor="#000000",
           linewidth=2, label=nm)
ax.set_xticks(xs)
ax.set_xticklabels([f"{short_label(q)}\n{sec[q]['cond'][:4]}" for q in piv.index])
ax.set_ylabel("Percent of drawn 'normal' area\ncalled lesion by the detector")
ax.set_title("Lower is unambiguously better")
ax.legend(frameon=False, fontsize=FONT_SIZE_LEGEND - 8)
style_axes(ax); panel_letter(ax, "B")
fig.suptitle("The negative control", fontsize=FONT_SIZE_TITLE, y=1.02)
fig.tight_layout()
save_fig(fig, "F101_normal_violations")


# %% Cell 11 - wrap up
# =============================================================================

banner("SUMMARY")

o = objs.loc[objs["setting"] == PRIMARY]
print(f"Setting under review : {PRIMARY}")
print(f"Objects detected     : {len(o)}")
print(f"Violations           : {int(o['is_violation'].sum())} "
      f"({100.0 * o['is_violation'].mean():.1f}%)")
tv = o.loc[o["condition"] == "D1MT"]
if len(tv):
    print(f"Treated objects      : {len(tv)}  "
          f"lesion-like {int((tv['call_own_arm'] == 'lesion').sum())}, "
          f"normal-like {int((tv['call_own_arm'] == 'normal').sum())}, "
          f"violations {int(tv['is_violation'].sum())}")

sub("Read in this order")
print("  1. Table 134 and F101. The negative control. The percent of drawn")
print("     'normal' area a setting calls lesion has no trade against it, so")
print("     it is the cleanest number in the comparison.")
print("  2. Table 133 and F100. Whether the treated objects look like drawn")
print("     lesions or drawn normal parenchyma, under both references.")
print("  3. F102_review_G3_43106. The fourteen objects, numbered. This is the")
print("     one that settles it by eye. Tell me which numbers are real.")

sub("What each outcome means")
print("  Violations low and treated objects lesion-like  -> adopt route B and")
print("     reword Finding 1 to match the higher treated count.")
print("  Violations high or treated objects normal-like  -> add an absolute")
print("     floor under the relative level and re-score.")
print("  Split                                            -> settle it on F102")
print("     by eye; the composition call is not decisive on its own.")

sub("Not changed by this script")
print("  Nothing. The burden definition, script 04, the cell assignment files")
print("  and every downstream script are untouched.")

banner("DONE")
sys.stdout = _tee.terminal
_tee.close()
