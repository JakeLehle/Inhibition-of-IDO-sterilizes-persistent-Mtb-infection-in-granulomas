#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - CANDIDATE DETECTOR, SCORED AGAINST EXPERT ANNOTATIONS
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 04e of the AKOYA analysis series. READ-ONLY.

WHAT THIS IS

    A candidate script 04 revision 5, run beside the current settings and scored
    against the 83 expert annotations. Nothing is adopted here. The output is a
    before-and-after on a real objective, so the decision to change script 04 is
    made on evidence rather than on a count.

WHAT CHANGED, AND WHAT MEASUREMENT DROVE EACH CHANGE

    1. BACKGROUND_STAT: median -> p25 of tissue-pixel density
       The regions drawn on G4_36463 peak at 4.9x its median background, BELOW
       the 6x gate, while every other section sits at 7.7 to 19.8. Under p25 all
       five sections clear 6x, 36463 at 12.5, and the cross-section spread of
       peak fold tightens from 4.0 to 3.7. The median of a heavily involved
       section is itself pathology, so a self-referencing gate is hardest on the
       most diseased tissue. That is backwards and this fixes it.

    2. FOCUS_HALF_MAX_FRACTION: 0.50 -> 0.25
       Traced boundaries sit at 0.224 of the enclosed peak for 'complex' and
       0.258 for 'clean_edge'. The pooled 0.282 in table 116 is dragged up by
       the 'missed' class at 0.637, and those regions are faint throughout so
       their boundary sits near their peak almost by construction. They are
       excluded from this estimate. Script 04 stops growing about twice too
       early.

    3. PEAK_SEPARATION_UM: UNCHANGED at 200
       The measurement said 772, from 11 multi-maxima regions. That would make
       the maximum_filter footprint 1,575 um and would merge genuinely distinct
       foci in the three sections whose detection is already acceptable. It is
       the right diagnosis of G4_43109 and the wrong instrument. Route B below
       addresses confluence directly instead.

    4. NEW ROUTE B, a level-set detector for confluent tissue
       G4_36463 contains 13 local maxima in the whole section, and the 24 drawn
       regions hold 9 of them between them with a median of ZERO maxima per
       region. Roughly fifteen drawn regions contain no local maximum at all.
       No gate creates a maximum where the density surface is a plateau, so
       peak-plus-watershed cannot see that tissue at any setting.

       Route B takes connected components of the density surface above a
       section-relative level, with no peak requirement. A confluent mass is one
       component because it is one component, and a fuzzy edge is handled by
       where the level sits rather than by a fraction of some peak that may not
       exist.

THE TWO ROUTES ARE SEPARATE OUTPUTS AND ARE NEVER MERGED

    Route A, peaks: for internal architecture. The radial coordinate needs a
        centre and an inscribed radius, both of which are meaningful for a
        compact focus and much less so for a large irregular complex.
    Route B, level sets: for counting and extent.

    This mirrors the existing absolute-for-burden against relative-for-
    architecture split. Keeping them separate matters: collapsing 43109's 44
    peaks into complexes would otherwise move the Finding 3 coefficient as a
    side effect rather than as a decision.

SCORING: RECALL ONLY, NEVER PRECISION

    The annotations are NOT exhaustive. A detected object with no annotation
    over it is not a false positive, it is unannotated tissue. So this reports
    coverage, IoU, splits and misses against drawn regions, and never a
    precision or an F-score. Reporting precision here would be dishonest.

    Objects detected outside any annotation ARE counted and shown, labelled as
    unannotated rather than as errors, because a large jump in that number is
    still worth seeing.

WHAT IS NOT TOUCHED

    The ABSOLUTE burden definition. One fixed threshold of 3,000 cells/mm2
    applied identically to both arms, which is the yardstick producing the zero
    in the treated animals. It is not modified and not scored here.

INPUTS
    annotations/*_annotations.geojson
    structures_rev4/cell_assignments/*_cell_structures.csv

OUTPUTS
    candidate_rev5/tables/
        120_settings_scored, 121_per_section_counts, 122_baseline_vs_candidate,
        123_route_b_sweep, 124_treated_foci_identifiability,
        125_region_level_detail, 00_candidate_rev5_report.txt
    candidate_rev5/figures/
        F96_score_sweep, F97_before_after_maps, F98_coverage_by_class,
        F99_counts_and_identifiability

USAGE
    conda activate sc_pre
    python AKOYA_04e_Candidate_Detector.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================

ANN_DIR = "/master/jlehle/WORKING/AKOYA/annotations"
IN_DIR = "/master/jlehle/WORKING/AKOYA/structures_rev4"
CELL_DIR = f"{IN_DIR}/cell_assignments"
OUT_DIR = "/master/jlehle/WORKING/AKOYA/candidate_rev5"

USE_AGG = True

# ---- fixed, shared by every setting -----------------------------------------
GRID_UM = 25.0
DENSITY_BANDWIDTH_UM = 75.0
EDGE_MIN_WEIGHT = 0.25
TISSUE_CLOSE_ITERATIONS = 2
TISSUE_FILL_HOLES = True
FOCUS_FILL_HOLES = True
BURDEN_THRESHOLD = 3000.0          # reference only, never modified

# ---- BASELINE, script 04 revision 4 as it stands now ------------------------
BASE = dict(bg_stat="median", fold=6.0, half_frac=0.50, sep_um=200.0,
            min_area=10000.0, min_cells=40)

# ---- CANDIDATE route A ------------------------------------------------------
CAND_A = dict(bg_stat="p25", fold=6.0, half_frac=0.25, sep_um=200.0,
              min_area=10000.0, min_cells=40)

# ---- route A sweep ----------------------------------------------------------
A_BG_STATS = ["median", "p25"]
A_HALF_FRACS = [0.50, 0.35, 0.25, 0.20]
A_FOLDS = [6.0]                    # held; the background change is the lever

# ---- route B sweep ----------------------------------------------------------
# The level IS the boundary, so the range is anchored on the measured boundary
# fold over p25: 3.4, 3.8, 6.6, 8.9, 11.8 across the five sections, median 6.6.
B_BG_STATS = ["p25"]
B_LEVEL_FOLDS = [3.0, 4.0, 5.0, 6.5, 8.0, 10.0]
B_CLOSE_UM = [0.0, 50.0, 100.0]    # bridge thin gaps between adjacent lesions
B_MIN_AREA_UM2 = 10000.0
B_MIN_CELLS = 40

# ---- scoring ----------------------------------------------------------------
MISSED_COVERAGE_PCT = 5.0          # a drawn region under this is "missed"
SCORE_CLASSES = ["complex", "missed", "clean_edge"]   # 'normal' is not a target
EXCLUDE_STAMPS_FROM_BOUNDARY_SCORING = True
MIN_PCT_IN_TISSUE = 80.0           # regions below this are flagged, not dropped

# ---- phenotypes -------------------------------------------------------------
IDO1_POS = "CD68+IDO1+ Macrophages"
IDO1_NEG = "CD68+IDO1- Macrophages"
MYELOID_FOR_DETECTION = [IDO1_POS, IDO1_NEG, "CD163+ Macrophages", "Neutrophils"]

# ---- appearance -------------------------------------------------------------
CONDITION_ORDER = ["D1MT", "Untreated"]
CONDITION_COLORS = {"D1MT": "#2C7FB8", "Untreated": "#D95F02"}
CONDITION_LABELS = {"D1MT": "D1MT-treated", "Untreated": "Untreated"}
CLASS_COLORS = {"complex": "#D73027", "missed": "#4575B4",
                "clean_edge": "#1A9850", "normal": "#999999"}
BASE_COLOR = "#999999"
CAND_A_COLOR = "#00A0C6"
CAND_B_COLOR = "#7B3294"
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


_tee = Tee(os.path.join(TAB_DIR, "00_candidate_rev5_report.txt"))
sys.stdout = _tee

banner("AKOYA CANDIDATE DETECTOR, SCORED (script 04e, read-only)")
print(f"Run time : {datetime.now().isoformat(timespec='seconds')}")
print(f"Output   : {OUT_DIR}")
print("\nNothing is adopted. structures_rev4 is untouched. The absolute burden")
print(f"definition stays at {BURDEN_THRESHOLD:,.0f} cells/mm2 and is not scored here.")
print("\nRECALL ONLY. The annotations are not exhaustive, so a detected object")
print("with no annotation over it is unannotated tissue, not a false positive.")
print("No precision or F-score is reported anywhere in this output.")


# %% Cell 3 - shared machinery
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


def background_of(dens, occ, stat, verbose=True, tag=""):
    """
    Section background, with an escalation guard.

    A low percentile can legitimately come out at ZERO: if more than a quarter
    of tissue pixels carry no myeloid density at all, p25 is 0 and every fold
    becomes infinite. The median cannot do this in practice but p25 and p10 can,
    and p25 is the candidate default, so the guard matters.

    On a zero or non-finite result the statistic escalates through higher
    percentiles until it is positive, and says so loudly. A silent fallback here
    would put an invented background under every downstream number.
    """
    v = dens[occ]; v = v[np.isfinite(v)]
    if not v.size:
        return np.nan

    def _stat(st):
        if st in ("median", "p50"):
            return float(np.median(v))
        if st.startswith("p") and st[1:].replace(".", "").isdigit():
            return float(np.percentile(v, float(st[1:])))
        return float(np.median(v))

    bg = _stat(stat)
    if np.isfinite(bg) and bg > 0:
        return bg
    for alt in ["p25", "p50", "median", "p75", "p90"]:
        b = _stat(alt)
        if np.isfinite(b) and b > 0:
            if verbose:
                print(f"    WARNING [{tag}]: background statistic '{stat}' is "
                      f"{bg}, which makes every fold infinite.")
                print(f"    Escalated to '{alt}' = {b:.2f} cells/mm2. This means "
                      f"over {stat[1:] if stat.startswith('p') else '50'} percent "
                      f"of tissue pixels carry no myeloid density, so a low")
                print(f"    percentile is not usable as a background here.")
            return b
    if verbose:
        print(f"    ERROR [{tag}]: no positive background at any percentile.")
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
    """Peak plus watershed. Script 04's method, with fold, growth fraction and
    background supplied by the caller."""
    sep_px = max(1, int(round(sep_um / GRID_UM)))
    py, px, pv = find_local_maxima(dens, occ, sep_px)
    keep = pv >= fold * bg
    py, px, pv = py[keep], px[keep], pv[keep]
    if not len(pv):
        return np.zeros_like(dens, int), []
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
    labels = np.zeros(dens.shape, int); recs = []; nid = 1
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
        area = assigned.sum() * GRID_UM * GRID_UM
        ncell = int(assigned[gy, gx].sum())
        if area < min_area or ncell < min_cells:
            continue
        dt = ndi.distance_transform_edt(assigned)
        labels[assigned] = nid
        recs.append(dict(object_id=nid, peak_row=int(yy), peak_col=int(xx),
                         peak_density=float(vv),
                         fold=float(vv / bg) if bg > 0 else np.nan,
                         area_um2=float(area), n_cells=ncell,
                         max_inscribed_radius_um=float(dt.max()) * GRID_UM))
        nid += 1
    return labels, recs


def route_b(dens, occ, gx, gy, bg, level_fold, close_um, min_area, min_cells):
    """
    Connected components of the density surface above a section-relative level.

    No peak is required, so a confluent plateau is one object because it is one
    connected component. This is what G4_36463 needs: 15 of its 24 drawn regions
    contain no local maximum at all, and no gate applied to a peak list can
    recover a region that never generated a peak.

    close_um bridges thin gaps between adjacent lesions before labelling. It is
    a morphological closing, so it can only join things already nearly touching.
    """
    level = level_fold * bg
    m = (dens >= level) & occ
    if close_um > 0:
        r = max(1, int(round(close_um / GRID_UM)))
        st = ndi.generate_binary_structure(2, 1)
        m = ndi.binary_closing(m, structure=st, iterations=r)
        m &= occ
    if FOCUS_FILL_HOLES:
        m = ndi.binary_fill_holes(m)
    lab, n = ndi.label(m)
    if n == 0:
        return np.zeros_like(dens, int), []
    labels = np.zeros(dens.shape, int); recs = []; nid = 1
    for i in range(1, n + 1):
        comp = lab == i
        area = comp.sum() * GRID_UM * GRID_UM
        if area < min_area:
            continue
        ncell = int(comp[gy, gx].sum())
        if ncell < min_cells:
            continue
        vals = dens[comp]
        dt = ndi.distance_transform_edt(comp)
        yy, xx = np.unravel_index(np.argmax(np.where(comp, dens, -np.inf)),
                                  dens.shape)
        labels[comp] = nid
        recs.append(dict(object_id=nid, peak_row=int(yy), peak_col=int(xx),
                         peak_density=float(vals.max()),
                         fold=(float(vals.max() / bg) if bg > 0 else np.nan),
                         level_density=float(level),
                         area_um2=float(area), n_cells=ncell,
                         max_inscribed_radius_um=float(dt.max()) * GRID_UM))
        nid += 1
    return labels, recs


# %% Cell 4 - load sections and annotations
# =============================================================================

banner("LOADING")

ann_files = sorted(glob.glob(os.path.join(ANN_DIR, "*_annotations.geojson")))
SECTIONS = [os.path.basename(f).replace("_annotations.geojson", "")
            for f in ann_files]

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
    sec[sid] = dict(gx=gx, gy=gy, gx0=gx0, gy0=gy0, shape=shape, occ=occ,
                    dens=dens, n_cells=len(d),
                    cond=d["condition"].iloc[0] if "condition" in d.columns else "?",
                    bg={st: background_of(dens, occ, st, tag=f"{sid}/{st}")
                        for st in sorted(set(A_BG_STATS + B_BG_STATS
                                             + ["median", "p25"]))})
    b = sec[sid]["bg"]
    print(f"    {sid:<14} {len(d):>9,} cells   bg median {b['median']:>8.1f}   "
          f"p25 {b['p25']:>8.1f}   tissue "
          f"{occ.sum() * GRID_UM ** 2 / 1e6:>5.1f} mm2")

SECTIONS = [s for s in SECTIONS if s in sec]
SECTIONS.sort(key=lambda s: (CONDITION_ORDER.index(sec[s]["cond"])
                             if sec[s]["cond"] in CONDITION_ORDER else 9, s))

# A section whose background cannot be established for a statistic is dropped
# from every run using it. Reported once here rather than once per combination.
_unusable = {}
for st in sorted(set(A_BG_STATS + B_BG_STATS)):
    bad = [q for q in SECTIONS
           if not np.isfinite(sec[q]["bg"].get(st, np.nan))
           or sec[q]["bg"].get(st, 0) <= 0]
    if bad:
        _unusable[st] = bad
if _unusable:
    print("\n    Sections with no usable background, excluded from any run")
    print("    using that statistic:")
    for st, bad in _unusable.items():
        print(f"      {st:<8} {', '.join(bad)}")
    print("    Every score below is computed on the remaining sections only,")
    print("    so scores across different background statistics are NOT")
    print("    directly comparable when this list is non-empty.")


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

masks = {}
for i, r in ann.iterrows():
    masks[i] = rasterise(r["rings_um"], sec[r["sample_id"]])
ann["n_px"] = [int(masks[i].sum()) for i in ann.index]
ann["pct_in_tissue"] = [100.0 * (masks[i] & sec[r["sample_id"]]["occ"]).sum()
                        / max(1, masks[i].sum()) for i, r in ann.iterrows()]

TARGETS = ann.loc[ann["classification"].isin(SCORE_CLASSES)].copy()
if EXCLUDE_STAMPS_FROM_BOUNDARY_SCORING:
    TARGETS = TARGETS.loc[~TARGETS["is_stamp"]]

print(f"\n    {len(ann)} annotations, {len(TARGETS)} scoreable targets")
print(f"    ({int(ann['is_stamp'].sum())} brush stamps excluded from scoring: a")
print("     stamped disc has no boundary, so an overlap against it is not")
print("     meaningful. They remain in the maxima work in script 04d.)")
low = ann.loc[ann["pct_in_tissue"] < MIN_PCT_IN_TISSUE]
if len(low):
    print(f"\n    {len(low)} region(s) under {MIN_PCT_IN_TISSUE:.0f}% inside the "
          f"tissue mask, kept but flagged:")
    print(low[["sample_id", "index", "classification",
               "pct_in_tissue"]].to_string(index=False))


# %% Cell 5 - the scorer
# =============================================================================

banner("SCORING")

print("Per drawn region: what fraction of it the detector covers, the IoU, and")
print("how many detector objects it overlaps. Per detector object: how many")
print("drawn regions it spans. Splits and over-merges are both reported,")
print("because a setting can improve one by wrecking the other.\n")


def score(labels, sid):
    """
    CORRECTED in this revision. The previous version computed each drawn
    region's IoU against `labels > 0`, which is EVERY detected object in the
    section, so the union was the whole detected area and every IoU came out at
    0.015 to 0.025 regardless of setting. That made the one metric sensitive to
    object SIZE useless, and size is exactly what the growth-fraction change is
    meant to fix.

    Two IoUs are now reported:
      best_object_iou  the single best-matching detected object. Penalises a
                       detector that covers a drawn region with many fragments.
      matched_iou      the union of objects overlapping the region. Fair when a
                       region is legitimately covered by several objects, e.g. a
                       drawn complex containing two real foci.
    Read them together: best much lower than matched means fragmentation.
    """
    tg = TARGETS.loc[TARGETS["sample_id"] == sid]
    rows = []
    for i, r in tg.iterrows():
        m = masks[i]
        ids = np.unique(labels[m]); ids = ids[ids > 0]
        if len(ids):
            matched = np.isin(labels, ids)
            inter_m = int((m & matched).sum())
            union_m = int((m | matched).sum())
            best = 0.0
            best_area = np.nan
            for oid in ids:
                om = labels == oid
                inter_o = int((m & om).sum())
                union_o = int((m | om).sum())
                v = inter_o / union_o if union_o else 0.0
                if v > best:
                    best = v
                    best_area = float(om.sum()) * GRID_UM * GRID_UM
        else:
            matched = np.zeros_like(m)
            inter_m = 0; union_m = int(m.sum())
            best = 0.0; best_area = np.nan
        drawn_area = float(m.sum()) * GRID_UM * GRID_UM
        rows.append(dict(sample_id=sid, index=r["index"],
                         classification=r["classification"],
                         pct_covered=(100.0 * int((m & (labels > 0)).sum())
                                      / m.sum()) if m.sum() else np.nan,
                         best_object_iou=float(best),
                         matched_iou=(inter_m / union_m) if union_m else np.nan,
                         n_objects_overlapping=int(len(ids)),
                         drawn_area_um2=drawn_area,
                         best_object_area_um2=best_area,
                         area_ratio=(best_area / drawn_area
                                     if np.isfinite(best_area) and drawn_area
                                     else np.nan)))
    # over-merge: one detector object spanning several drawn regions
    over = 0
    n_obj = int(labels.max())
    for oid in range(1, n_obj + 1):
        om = labels == oid
        hit = sum(1 for i in tg.index if (masks[i] & om).any())
        if hit >= 2:
            over += 1
    # objects with no annotation over them: unannotated, NOT false positives
    unann = 0
    for oid in range(1, n_obj + 1):
        om = labels == oid
        if not any((masks[i] & om).any() for i in tg.index):
            unann += 1
    return pd.DataFrame(rows), dict(n_objects=n_obj, n_over_merged=over,
                                    n_unannotated=unann)


def summarise(det_rows, obj_stats, tag, params):
    if not len(det_rows):
        return None
    d = pd.concat(det_rows, ignore_index=True)
    o = pd.DataFrame(obj_stats)
    rec = dict(route=tag, **params,
               n_targets=len(d),
               median_coverage_pct=float(d["pct_covered"].median()),
               mean_coverage_pct=float(d["pct_covered"].mean()),
               median_best_iou=float(d["best_object_iou"].median()),
               median_matched_iou=float(d["matched_iou"].median()),
               median_area_ratio=float(d["area_ratio"].median()),
               n_missed=int((d["pct_covered"] < MISSED_COVERAGE_PCT).sum()),
               n_split=int((d["n_objects_overlapping"] >= 2).sum()),
               n_objects=int(o["n_objects"].sum()),
               n_over_merged=int(o["n_over_merged"].sum()),
               n_unannotated=int(o["n_unannotated"].sum()))
    for cls in SCORE_CLASSES:
        s = d.loc[d["classification"] == cls]
        rec[f"coverage_{cls}"] = float(s["pct_covered"].median()) if len(s) else np.nan
        rec[f"missed_{cls}"] = int((s["pct_covered"] < MISSED_COVERAGE_PCT).sum())
    return rec, d


def run_a(params):
    det_rows, obj_stats, labs, counts = [], [], {}, {}
    for sid in SECTIONS:
        s = sec[sid]
        bg = s["bg"][params["bg_stat"]]
        if not np.isfinite(bg) or bg <= 0:
            continue        # reported once at load, not per combination
        lab, recs = route_a(s["dens"], s["occ"], s["gx"], s["gy"], bg,
                            params["fold"], params["half_frac"],
                            params["sep_um"], params["min_area"],
                            params["min_cells"])
        r, st = score(lab, sid)
        det_rows.append(r); obj_stats.append(st)
        labs[sid] = lab; counts[sid] = len(recs)
    return det_rows, obj_stats, labs, counts


def run_b(params):
    det_rows, obj_stats, labs, counts = [], [], {}, {}
    for sid in SECTIONS:
        s = sec[sid]
        bg = s["bg"][params["bg_stat"]]
        if not np.isfinite(bg) or bg <= 0:
            continue        # reported once at load, not per combination
        lab, recs = route_b(s["dens"], s["occ"], s["gx"], s["gy"], bg,
                            params["level_fold"], params["close_um"],
                            B_MIN_AREA_UM2, B_MIN_CELLS)
        r, st = score(lab, sid)
        det_rows.append(r); obj_stats.append(st)
        labs[sid] = lab; counts[sid] = len(recs)
    return det_rows, obj_stats, labs, counts


# ---- route A sweep ----------------------------------------------------------
sub("Route A: peaks plus watershed")
a_rows, a_counts = [], {}
for bgs in A_BG_STATS:
    for hf in A_HALF_FRACS:
        for fold in A_FOLDS:
            p = dict(bg_stat=bgs, fold=fold, half_frac=hf, sep_um=200.0,
                     min_area=10000.0, min_cells=40)
            dr, os_, labs, cnt = run_a(p)
            out = summarise(dr, os_, "A", p)
            if out:
                a_rows.append(out[0])
                a_counts[(bgs, fold, hf)] = cnt
            print(f"    bg {bgs:<7} fold {fold:g}  half {hf:.2f}   "
                  f"coverage {out[0]['median_coverage_pct']:>5.1f}%  "
                  f"IoU {out[0]['median_best_iou']:.3f}/"
                  f"{out[0]['median_matched_iou']:.3f}  "
                  f"size x{out[0]['median_area_ratio']:.2f}  missed "
                  f"{out[0]['n_missed']:>2}  split {out[0]['n_split']:>2}  "
                  f"objects {out[0]['n_objects']:>3}")

# ---- route B sweep ----------------------------------------------------------
sub("Route B: level sets, no peak required")
b_rows, b_counts = [], {}
for bgs in B_BG_STATS:
    for lf in B_LEVEL_FOLDS:
        for cu in B_CLOSE_UM:
            p = dict(bg_stat=bgs, level_fold=lf, close_um=cu)
            dr, os_, labs, cnt = run_b(p)
            out = summarise(dr, os_, "B", p)
            if out:
                b_rows.append(out[0])
                b_counts[(bgs, lf, cu)] = cnt
            print(f"    bg {bgs:<5} level {lf:>5.1f}x  close {cu:>5.0f}um   "
                  f"coverage {out[0]['median_coverage_pct']:>5.1f}%  "
                  f"IoU {out[0]['median_best_iou']:.3f}/"
                  f"{out[0]['median_matched_iou']:.3f}  "
                  f"size x{out[0]['median_area_ratio']:.2f}  missed "
                  f"{out[0]['n_missed']:>2}  split {out[0]['n_split']:>2}  "
                  f"merged {out[0]['n_over_merged']:>2}  "
                  f"objects {out[0]['n_objects']:>3}")

scored = pd.DataFrame(a_rows + b_rows)
write_csv(scored, "120_settings_scored.csv")
write_csv(pd.DataFrame(b_rows), "123_route_b_sweep.csv")


# %% Cell 6 - baseline against candidate
# =============================================================================

banner("BASELINE AGAINST CANDIDATE")

base_dr, base_os, base_labs, base_cnt = run_a(BASE)
base_sum, base_det = summarise(base_dr, base_os, "A_baseline", BASE)

cand_dr, cand_os, cand_labs, cand_cnt = run_a(CAND_A)
cand_sum, cand_det = summarise(cand_dr, cand_os, "A_candidate", CAND_A)

# best route B by coverage, with over-merging held in check
bdf = pd.DataFrame(b_rows)
bdf_ok = bdf.loc[bdf["n_over_merged"] <= max(2, int(0.1 * len(TARGETS)))]
best_b = (bdf_ok.sort_values(["median_matched_iou", "median_coverage_pct"],
                             ascending=False).iloc[0]
          if len(bdf_ok) else bdf.sort_values("median_coverage_pct",
                                              ascending=False).iloc[0])
BEST_B = dict(bg_stat=best_b["bg_stat"], level_fold=float(best_b["level_fold"]),
              close_um=float(best_b["close_um"]))
b_dr, b_os, b_labs, b_cnt = run_b(BEST_B)
b_sum, b_det = summarise(b_dr, b_os, "B_best", BEST_B)

cmp_rows = []
for tag, sm, extra in [("baseline (script 04 rev4)", base_sum, BASE),
                       ("candidate route A", cand_sum, CAND_A),
                       ("candidate route B", b_sum, BEST_B)]:
    r = dict(setting=tag)
    r.update({k: v for k, v in sm.items() if k not in ("route",)})
    cmp_rows.append(r)
cmp = pd.DataFrame(cmp_rows)
write_csv(cmp, "122_baseline_vs_candidate.csv")

show = ["setting", "median_coverage_pct", "median_best_iou",
        "median_matched_iou", "median_area_ratio", "n_missed", "n_split",
        "n_over_merged", "n_objects", "n_unannotated"]
with pd.option_context("display.width", 260):
    print(cmp[show].to_string(index=False))

sub("Coverage by annotation class (median percent of each drawn region covered)")
cls_cols = [f"coverage_{c}" for c in SCORE_CLASSES] + \
           [f"missed_{c}" for c in SCORE_CLASSES]
with pd.option_context("display.width", 260):
    print(cmp[["setting"] + cls_cols].to_string(index=False))
print("\n    The 'missed' class is the one to watch. Those are regions you")
print("    marked as lesion that the detector treats as background, and the")
print("    baseline recovers essentially none of them.")

sub("Per-section object counts")
cnt_rows = []
for sid in SECTIONS:
    cnt_rows.append(dict(sample_id=sid, condition=sec[sid]["cond"],
                         drawn_complex=int(((ann["sample_id"] == sid) &
                                            (ann["classification"] == "complex")).sum()),
                         drawn_missed=int(((ann["sample_id"] == sid) &
                                           (ann["classification"] == "missed")).sum()),
                         baseline=base_cnt.get(sid, np.nan),
                         candidate_A=cand_cnt.get(sid, np.nan),
                         candidate_B=b_cnt.get(sid, np.nan)))
cnts = pd.DataFrame(cnt_rows)
write_csv(cnts, "121_per_section_counts.csv")
print(cnts.to_string(index=False))

# per-region detail for both
base_det["setting"] = "baseline"; cand_det["setting"] = "candidate_A"
b_det["setting"] = "candidate_B"
write_csv(pd.concat([base_det, cand_det, b_det], ignore_index=True),
          "125_region_level_detail.csv")


sub("TREATED-ARM SANITY CHECK")
print("The regions drawn on treated sections are few, so the annotation score")
print("cannot police over-detection there. The expert counts can. Finding 1")
print("rests on residual foci being FEW and small in treated animals, so a")
print("setting that multiplies the treated count changes that narrative as a")
print("side effect rather than as a decision.\n")
EXPERT_COUNTS = {"G3_43106": 1, "G3_43111": 2, "G3_43118": 2,
                 "G4_31438": 7, "G4_36463": 5, "G4_43109": 19}
EXPERT_IS_REAL = {"G3_43106": True, "G3_43111": True, "G3_43118": True,
                  "G4_31438": False, "G4_36463": False, "G4_43109": False}
chk = []
for sid in SECTIONS:
    chk.append(dict(sample_id=sid, condition=sec[sid]["cond"],
                    expert=EXPERT_COUNTS.get(sid, np.nan),
                    expert_is_real=EXPERT_IS_REAL.get(sid, False),
                    drawn=int(((ann["sample_id"] == sid)
                               & (ann["classification"].isin(SCORE_CLASSES))).sum()),
                    baseline=base_cnt.get(sid, np.nan),
                    candidate_A=cand_cnt.get(sid, np.nan),
                    candidate_B=b_cnt.get(sid, np.nan)))
chkdf = pd.DataFrame(chk)
print(chkdf.to_string(index=False))
write_csv(chkdf, "126_expert_count_check.csv")
for col in ("baseline", "candidate_A", "candidate_B"):
    real = chkdf.loc[chkdf["expert_is_real"]]
    r = (real[col] / real["expert"]).replace([np.inf, -np.inf], np.nan).dropna()
    if len(r):
        print(f"\n    {col:<12} objects per expert count on the sections with")
        print(f"                 GENUINE expert counts: "
              + ", ".join(f"{s}={v:.1f}x" for s, v in
                          zip(real['sample_id'], r)))
print("\n    A ratio far above 1 on a treated section means the setting is")
print("    finding many objects where the pathologist counted one. That is not")
print("    settled by the annotation score and has to be looked at on F97.")

sub("ROUTE B FRONTIER, so the choice is not made by an arbitrary cap")
print("    The 'best' route B above is chosen subject to an over-merge cap that")
print("    is a judgement, not a measurement. The whole frontier is printed so")
print("    the trade is visible.\n")
frontier = (pd.DataFrame(b_rows)
            .sort_values(["level_fold", "close_um"])
            [["level_fold", "close_um", "median_coverage_pct",
              "median_best_iou", "median_matched_iou", "median_area_ratio",
              "n_missed", "n_split", "n_over_merged", "n_objects"]])
with pd.option_context("display.width", 260):
    print(frontier.to_string(index=False))


# %% Cell 7 - the identifiability side effect
# =============================================================================

banner("EFFECT ON THE MIXED MODEL")

print("Focus-level models fail when an animal contributes a single focus,")
print("because the animal random effect and the residual are then the same")
print("number and the variance slides to the boundary at zero. More foci per")
print("treated animal is therefore a side effect worth measuring, though it is")
print("NOT a reason to prefer a setting. The detection question is settled on")
print("agreement with the annotations alone.\n")

ident = []
for sid in SECTIONS:
    ident.append(dict(sample_id=sid, condition=sec[sid]["cond"],
                      baseline=base_cnt.get(sid, np.nan),
                      candidate_A=cand_cnt.get(sid, np.nan),
                      candidate_B=b_cnt.get(sid, np.nan)))
idf = pd.DataFrame(ident)
write_csv(idf, "124_treated_foci_identifiability.csv")
for col in ("baseline", "candidate_A", "candidate_B"):
    tr = idf.loc[idf["condition"] == "D1MT", col].dropna()
    print(f"    {col:<12} treated animals with 2+ objects: "
          f"{int((tr >= 2).sum())} of {len(tr)}   "
          f"counts {[int(v) for v in tr]}")
print("\n    Only the sections annotated so far are included, so this is")
print("    partial. G3_43111 is absent.")


# %% Cell 8 - figures
# =============================================================================

banner("FIGURES")

# ---- F96 sweep --------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(28, 13))
ax = axes[0]
adf = pd.DataFrame(a_rows)
for bgs, g in adf.groupby("bg_stat"):
    g = g.sort_values("half_frac")
    ax.plot(g["half_frac"], g["median_coverage_pct"], marker="o", markersize=18,
            linewidth=4, label=f"route A, bg {bgs}")
ax.axvline(BASE["half_frac"], color=BASE_COLOR, linestyle="--", linewidth=3)
ax.text(BASE["half_frac"], ax.get_ylim()[1], " current", va="top",
        fontsize=FONT_SIZE_TICK - 8)
ax.set_xlabel("FOCUS_HALF_MAX_FRACTION")
ax.set_ylabel("Median coverage of drawn regions (%)")
ax.set_title("Route A", fontsize=FONT_SIZE_TITLE - 2)
ax.legend(frameon=False, fontsize=FONT_SIZE_LEGEND - 8)
style_axes(ax); panel_letter(ax, "A")

ax = axes[1]
for cu, g in bdf.groupby("close_um"):
    g = g.sort_values("level_fold")
    ax.plot(g["level_fold"], g["median_coverage_pct"], marker="s",
            markersize=18, linewidth=4, label=f"closing {cu:.0f} um")
ax.axhline(base_sum["median_coverage_pct"], color=BASE_COLOR,
           linestyle="--", linewidth=4)
ax.text(ax.get_xlim()[1], base_sum["median_coverage_pct"], "  baseline",
        va="center", fontsize=FONT_SIZE_TICK - 8)
ax.set_xlabel("Level, fold over p25 background")
ax.set_ylabel("Median coverage of drawn regions (%)")
ax.set_title("Route B", fontsize=FONT_SIZE_TITLE - 2)
ax.legend(frameon=False, fontsize=FONT_SIZE_LEGEND - 8)
style_axes(ax); panel_letter(ax, "B")
fig.suptitle("Agreement with expert annotations, recall only",
             fontsize=FONT_SIZE_TITLE, y=1.02)
fig.tight_layout()
save_fig(fig, "F96_score_sweep")

# ---- F98 coverage by class --------------------------------------------------
fig, ax = plt.subplots(figsize=(22, 13))
w = 0.26
xs = np.arange(len(SCORE_CLASSES))
for k, (tag, det, col) in enumerate([("baseline", base_det, BASE_COLOR),
                                     ("candidate A", cand_det, CAND_A_COLOR),
                                     ("candidate B", b_det, CAND_B_COLOR)]):
    vals = [det.loc[det["classification"] == c, "pct_covered"].median()
            if (det["classification"] == c).any() else 0.0 for c in SCORE_CLASSES]
    ax.bar(xs + (k - 1) * w, vals, width=w, color=col, edgecolor="#000000",
           linewidth=2, label=tag)
ax.set_xticks(xs); ax.set_xticklabels(SCORE_CLASSES)
ax.set_ylabel("Median coverage of drawn regions (%)")
ax.set_title("Coverage by annotation class", fontsize=FONT_SIZE_TITLE)
ax.legend(frameon=False, fontsize=FONT_SIZE_LEGEND - 6)
style_axes(ax)
save_fig(fig, "F98_coverage_by_class")

# ---- F97 before and after maps ---------------------------------------------
nrow = len(SECTIONS)
fig, axes = plt.subplots(nrow, 3, figsize=(13 * 3, 11 * nrow))
axes = np.atleast_2d(axes)
for r_i, sid in enumerate(SECTIONS):
    s = sec[sid]
    for c_i, (tag, labs) in enumerate([("baseline", base_labs),
                                       ("candidate A", cand_labs),
                                       ("candidate B", b_labs)]):
        ax = axes[r_i, c_i]
        ax.imshow(np.log10(np.clip(np.where(s["occ"], s["dens"], np.nan), 1.0,
                                   None)), cmap="Greys", interpolation="nearest")
        if sid not in labs:
            ax.axis("off"); continue
        det = labs[sid] > 0
        if det.any():
            ax.contour(det.astype(float), levels=[0.5], colors=["#00A0C6"],
                       linewidths=2.5)
        for i, rr in ann.loc[ann["sample_id"] == sid].iterrows():
            if rr["classification"] not in SCORE_CLASSES:
                continue
            ax.contour(masks[i].astype(float), levels=[0.5],
                       colors=[CLASS_COLORS.get(rr["classification"], "#000")],
                       linewidths=2.5,
                       linestyles="--" if rr["is_stamp"] else "-")
        ax.set_xticks([]); ax.set_yticks([])
        cov = {"baseline": base_det, "candidate A": cand_det,
               "candidate B": b_det}[tag]
        cv = cov.loc[cov["sample_id"] == sid, "pct_covered"]
        ax.set_title(f"{sid}  {tag}\n{int(labs[sid].max())} objects, "
                     f"coverage {cv.median() if len(cv) else np.nan:.0f}%",
                     fontsize=FONT_SIZE_BASE - 8)
handles = [Line2D([0], [0], color="#00A0C6", lw=4, label="detected")]
handles += [Line2D([0], [0], color=CLASS_COLORS[c], lw=4, label=f"drawn: {c}")
            for c in SCORE_CLASSES]
fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
           fontsize=FONT_SIZE_LEGEND - 6, bbox_to_anchor=(0.5, -0.01))
fig.suptitle("Before and after, against the drawn regions",
             fontsize=FONT_SIZE_TITLE, y=1.0)
fig.tight_layout()
save_fig(fig, "F97_before_after_maps")

# ---- F99 counts -------------------------------------------------------------
fig, ax = plt.subplots(figsize=(22, 13))
xs = np.arange(len(SECTIONS))
for k, (col, c, lab) in enumerate([("baseline", BASE_COLOR, "baseline"),
                                   ("candidate_A", CAND_A_COLOR, "candidate A"),
                                   ("candidate_B", CAND_B_COLOR, "candidate B")]):
    ax.bar(xs + (k - 1) * 0.26, cnts[col], width=0.26, color=c,
           edgecolor="#000000", linewidth=2, label=lab)
ax.plot(xs, cnts["drawn_complex"] + cnts["drawn_missed"], "k*", markersize=30,
        linestyle="none", label="drawn regions")
ax.set_xticks(xs)
ax.set_xticklabels([f"{short_label(s)}\n{sec[s]['cond'][:4]}" for s in SECTIONS])
ax.set_ylabel("Objects detected")
ax.set_title("Object counts against what was drawn", fontsize=FONT_SIZE_TITLE)
ax.legend(frameon=False, fontsize=FONT_SIZE_LEGEND - 8)
style_axes(ax)
save_fig(fig, "F99_counts_and_identifiability")


# %% Cell 9 - wrap up
# =============================================================================

banner("SUMMARY")

d_cov = cand_sum["median_coverage_pct"] - base_sum["median_coverage_pct"]
d_miss = cand_sum["n_missed"] - base_sum["n_missed"]
print(f"Route A, background median -> p25 and growth 0.50 -> 0.25:")
print(f"    coverage {base_sum['median_coverage_pct']:.1f}% -> "
      f"{cand_sum['median_coverage_pct']:.1f}%   ({d_cov:+.1f})")
print(f"    missed   {base_sum['n_missed']} -> {cand_sum['n_missed']}   "
      f"({d_miss:+d} of {len(TARGETS)} targets)")
print(f"    splits   {base_sum['n_split']} -> {cand_sum['n_split']}")
print(f"\nRoute B best: level {BEST_B['level_fold']:g}x over "
      f"{BEST_B['bg_stat']}, closing {BEST_B['close_um']:.0f} um")
print(f"    coverage {b_sum['median_coverage_pct']:.1f}%   "
      f"missed {b_sum['n_missed']}   over-merged {b_sum['n_over_merged']}")

sub("Read in this order")
print("  1. Table 122. The three settings side by side.")
print("  2. F98. Coverage by class. The 'missed' bar is the point of route B.")
print("  3. F97. The maps. Check the candidate against what you drew,")
print("     especially G4_36463 and G4_43109.")
print("  4. Table 121 and F99. Counts, including whether route B over-merges")
print("     G4_43109 into too few objects.")

sub("What would improve this measurement")
print("  Ten brush stamps are excluded from scoring because a stamped disc has")
print("  no boundary. Six of those are 'complex' on G4_43109, which is the")
print("  section where route B matters most. Redrawing those six as polygons")
print("  would add them to the scoring and make the G4_43109 numbers")
print("  substantially more solid. Nothing else needs redoing.")

sub("Not changed by this script")
print("  Nothing. The burden definition, script 04, the cell assignment files")
print("  and every downstream script are untouched.")

banner("DONE")
sys.stdout = _tee.terminal
_tee.close()
