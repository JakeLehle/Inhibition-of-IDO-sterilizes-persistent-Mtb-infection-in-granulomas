#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - ANNOTATION CANVAS BUILDER
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 04c of the AKOYA analysis series. READ-ONLY.

WHY THIS EXISTS

    The raw Akoya imagery was never delivered. We hold the processed export:
    per-cell centroids in whole-scan micron coordinates plus Akoya's phenotype
    calls. So there is no qptiff to open in QuPath and no way to annotate the
    original image.

    This script renders the cell positions themselves to a PNG at an EXACT,
    chosen microns-per-pixel. Because the scale and origin are set here rather
    than read from an image header, mapping a polygon drawn on the PNG back to
    scan micron coordinates is exact arithmetic:

        x_um = X0_UM + (col + 0.5) * UM_PER_PX
        y_um = Y0_UM + (row + 0.5) * UM_PER_PX

    The transform is written to a sidecar JSON and baked into a generated
    QuPath export script, so the number never has to be typed by hand.

THE CIRCULARITY GUARD, WHICH IS THE POINT OF THE LAYER SPLIT

    The foci detector thresholds a Gaussian-smoothed pooled-myeloid density
    surface. That surface is a derivative of the same cells rendered here.
    Annotating ON the density map would therefore be drawing boundaries on the
    detector's own intermediate output, and any threshold calibrated against
    those boundaries would be calibrated against itself.

    So the annotation surface is the CELLS. Two layers are written for drawing:

        _phenotypes  every cell, coloured by Akoya phenotype
        _myeloid     the four pooled myeloid phenotypes only, which is the
                     population the detector actually counts

    A third image, _reference, shows the smoothed density surface with the
    CURRENT foci outlined. It is deliberately named and described as
    reference-only. Look at it AFTER drawing, to see where the detector and the
    expert disagree. Do not draw on it.

    This keeps the density definition of a focus intact, which is what we want:
    the detector stays a pure cell-density method, and the annotations only fix
    the two free parameters we have been guessing at, the boundary threshold and
    the merge distance.

WHAT THE ANNOTATIONS ARE FOR

    Not for defining foci. For measuring, on regions an expert drew:
      - the absolute density and fold at the boundary
      - what fraction of the enclosed peak the boundary sits at, which tests
        the 50 percent growth rule directly
      - how many local maxima fall inside one expert region, which quantifies
        over-segmentation and yields a merge distance
      - the density gradient across fuzzy against clean edges

    That replaces the six expert COUNTS, which were three degrees of freedom
    tuning three parameters, with region-level overlap on hundreds of pixels.

OUTPUTS
    annotation_canvas/
        <section>_phenotypes.png            draw here
        <section>_myeloid.png               or here
        <section>_reference_density.png     look at this AFTER drawing
        <section>_transform.json            the exact pixel to micron mapping
        <section>_legend.png
        QuPath_annotate_and_export.groovy   paste into QuPath, generated with
                                            the transforms already filled in
        README_annotation.txt
        00_annotation_canvas_report.txt

USAGE
    conda activate sc_pre
    python AKOYA_04c_Annotation_Canvas.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================

IN_DIR = "/master/jlehle/WORKING/AKOYA/structures_rev4"
CELL_DIR = f"{IN_DIR}/cell_assignments"
FOCI_TABLE = f"{IN_DIR}/tables/35_foci_structures_relative.csv"
OUT_DIR = "/master/jlehle/WORKING/AKOYA/annotation_canvas"

# Where the exports will land when the generated Groovy script runs. If QuPath
# is on a laptop rather than titan, change this in the .groovy after it is
# written, or here before running.
ANNOTATION_OUT_DIR = "/master/jlehle/WORKING/AKOYA/annotations"

# Sections to render. None means auto-discover every section in CELL_DIR,
# which is the six retained sections. All six are rendered on purpose: the
# sections whose detection is already acceptable supply the clean_edge contrast
# cases, and without them there is no way to tell a fuzzy boundary from a sharp
# one numerically rather than by assertion.
SECTIONS = None
EXCLUDE_SECTIONS = ["G3_43102", "G4_43112"]   # position-1, excluded upstream

# Per-section context, printed into the README so it is clear what to look for
# on each canvas. Only the G3 expert counts are genuine; the G4 targets are
# algorithmic and are labelled as such.
EXPERT_COUNTS = {"G3_43106": 1, "G3_43111": 2, "G3_43118": 2,
                 "G4_31438": 7, "G4_36463": 5, "G4_43109": 19}
EXPERT_IS_REAL = {"G3_43106": True, "G3_43111": True, "G3_43118": True,
                  "G4_31438": False, "G4_36463": False, "G4_43109": False}
SECTION_NOTES = {
    "G3_43106": "detects 3 where the expert counted 1. Look at whether the two "
                "extra are real foci or split pieces of one.",
    "G3_43111": "detects 1 where the expert counted 2. One dense lesion then a "
                "cliff. Look for the second lesion the detector misses.",
    "G3_43118": "detects 2, expert 2. Agreement case, so clean_edge examples "
                "here are especially useful.",
    "G4_31438": "detection acceptable. Mainly a source of clean_edge contrast "
                "cases.",
    "G4_36463": "the confluent one. Smallest section, highest cell density, "
                "background 1,686/mm2 so a 6x gate means 10,117 which is above "
                "burden density. Regions here never became candidates at all. "
                "Mark those 'missed' and mark uninvolved parenchyma 'normal'.",
    "G4_43109": "44 foci against an algorithmic 19. Over-segmentation. Draw "
                "granuloma complexes as you would count them, so the merge "
                "distance can be measured from the maxima inside them.",
}

# ---- canvas geometry --------------------------------------------------------
# 2.0 um/px puts 43109 at about 4955 x 4920 and 36463 at 3234 x 1387, which
# QuPath opens instantly and which resolves a 10 um cell as a 5 px dot.
# 1.0 gives finer edges at four times the file size.
UM_PER_PX = 2.0
MARGIN_UM = 100.0
CELL_RADIUS_PX = 1          # dilation radius; 1 gives a 3 px dot at 2 um/px
MYELOID_RADIUS_PX = 2       # myeloid drawn larger on its own layer

# ---- density reference layer, must match script 04 --------------------------
GRID_UM = 25.0
DENSITY_BANDWIDTH_UM = 75.0
EDGE_MIN_WEIGHT = 0.25
TISSUE_CLOSE_ITERATIONS = 2
TISSUE_FILL_HOLES = True
BURDEN_THRESHOLD = 3000.0
FOCUS_FOLD_OVER_BACKGROUND = 6.0

# ---- phenotypes -------------------------------------------------------------
IDO1_POS = "CD68+IDO1+ Macrophages"
IDO1_NEG = "CD68+IDO1- Macrophages"
MYELOID_FOR_DETECTION = [IDO1_POS, IDO1_NEG, "CD163+ Macrophages", "Neutrophils"]

# Draw order, first drawn is furthest back. Structural and unassigned cells go
# down first so immune populations read on top.
DRAW_ORDER = [
    "Other", "Epithelial/Tumor cells", "Endothelial cells",
    "B cells", "Plasma cells", "Tregs", "CD4- T cells", "Helper T cells",
    "CD163+ Macrophages", "Neutrophils", IDO1_NEG, IDO1_POS,
]
PHENOTYPE_COLORS = {
    IDO1_POS: "#B2182B", IDO1_NEG: "#EF8A62", "CD163+ Macrophages": "#FDBE85",
    "Neutrophils": "#7B3294", "Helper T cells": "#1B7837",
    "CD4- T cells": "#7FBC41", "Tregs": "#00441B", "B cells": "#2166AC",
    "Plasma cells": "#67A9CF", "Endothelial cells": "#E7298A",
    "Epithelial/Tumor cells": "#66C2A5", "Other": "#A6761D",
}
CANVAS_BACKGROUND = "#FFFFFF"
MYELOID_COLOR = "#000000"
NONMYELOID_COLOR = "#C0C0C0"

# ---- annotation classes the Groovy script will expect -----------------------
ANNOTATION_CLASSES = ["complex", "missed", "clean_edge", "normal"]

DPI = 300
SAVE_PDF = False
SAVE_PNG = True
FONT_SIZE_BASE = 28
FONT_SIZE_TITLE = 32


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
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

try:
    from PIL import Image
    HAVE_PIL = True
except Exception:
    HAVE_PIL = False

os.makedirs(OUT_DIR, exist_ok=True)


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


def hex_to_rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))


_tee = Tee(os.path.join(OUT_DIR, "00_annotation_canvas_report.txt"))
sys.stdout = _tee

banner("AKOYA ANNOTATION CANVAS BUILDER (script 04c, read-only)")
print(f"Run time : {datetime.now().isoformat(timespec='seconds')}")
print(f"Output   : {OUT_DIR}")
print(f"\nScale    : {UM_PER_PX} um per pixel, chosen here rather than read from")
print("           an image header, so the mapping back to scan microns is exact.")
if not HAVE_PIL:
    print("\n    ERROR: Pillow is required to write the PNGs.")
    print("    pip install pillow  (or conda install pillow)")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)


# %% Cell 3 - load
# =============================================================================

banner("LOADING")

if SECTIONS is None:
    found = sorted(os.path.basename(q).replace("_cell_structures.csv", "")
                   for q in glob.glob(os.path.join(CELL_DIR, "*_cell_structures.csv")))
    SECTIONS = [q for q in found if q not in EXCLUDE_SECTIONS]
    print(f"    auto-discovered {len(SECTIONS)} sections: {', '.join(SECTIONS)}")
    skipped = [q for q in found if q in EXCLUDE_SECTIONS]
    if skipped:
        print(f"    excluded (position-1, dropped upstream): {', '.join(skipped)}")

cells = {}
for sid in SECTIONS:
    p = os.path.join(CELL_DIR, f"{sid}_cell_structures.csv")
    if not os.path.exists(p):
        print(f"    ERROR: {p} not found, skipping {sid}")
        continue
    d = pd.read_csv(p, low_memory=False)
    d = d.loc[np.isfinite(d["x"]) & np.isfinite(d["y"])].reset_index(drop=True)
    cells[sid] = d
    print(f"    {sid:<14} {len(d):>9,} cells   "
          f"x {d['x'].min():>9.1f} to {d['x'].max():>9.1f} um   "
          f"y {d['y'].min():>9.1f} to {d['y'].max():>9.1f} um")

if not cells:
    print("\nERROR: nothing to render.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

foci_ref = pd.read_csv(FOCI_TABLE) if os.path.exists(FOCI_TABLE) else None
if foci_ref is not None:
    print(f"\n    loaded table 35 for the reference layer "
          f"({foci_ref.shape[0]} foci)")

print("\n    Coordinates are whole-scan, not per-section. y values in the tens")
print("    of thousands are expected: the sections tile down one slide and were")
print("    never re-zeroed. The canvas origin below removes that offset, and the")
print("    transform puts it back.")


# %% Cell 4 - render
# =============================================================================

banner("RENDERING CANVASES")


def stamp(canvas, cols, rows, rgb, radius):
    """Paint cells as small discs. Dilation is done on a mask so overlapping
    cells merge into tissue rather than into a moire of single pixels."""
    h, w = canvas.shape[:2]
    ok = (cols >= 0) & (cols < w) & (rows >= 0) & (rows < h)
    if not ok.any():
        return
    m = np.zeros((h, w), dtype=bool)
    m[rows[ok], cols[ok]] = True
    if radius > 0:
        m = ndi.binary_dilation(m, iterations=int(radius))
    canvas[m] = rgb


def density_surface(d, x0, y0, w_px, h_px):
    """Script 04's density map, rendered onto the canvas grid for reference."""
    x = d["x"].to_numpy(float); y = d["y"].to_numpy(float)
    gx0 = float(np.floor(x.min() / GRID_UM) * GRID_UM)
    gy0 = float(np.floor(y.min() / GRID_UM) * GRID_UM)
    gx = ((x - gx0) / GRID_UM).astype(int) + 1
    gy = ((y - gy0) / GRID_UM).astype(int) + 1
    shape = (int(gy.max()) + 2, int(gx.max()) + 2)

    occ = np.zeros(shape, dtype=bool)
    occ[gy, gx] = True
    occ = ndi.binary_closing(occ, iterations=TISSUE_CLOSE_ITERATIONS)
    if TISSUE_FILL_HOLES:
        occ = ndi.binary_fill_holes(occ)

    mye = d["pheno"].isin(MYELOID_FOR_DETECTION).to_numpy()
    counts = np.zeros(shape, dtype=float)
    np.add.at(counts, (gy[mye], gx[mye]), 1.0)
    sigma = DENSITY_BANDWIDTH_UM / GRID_UM
    sm = ndi.gaussian_filter(counts, sigma=sigma, mode="constant")
    wgt = ndi.gaussian_filter(occ.astype(float), sigma=sigma, mode="constant")
    dens = (sm / np.maximum(wgt, EDGE_MIN_WEIGHT)) / ((GRID_UM ** 2) / 1e6)
    dens[~occ] = np.nan
    return dens, occ, gx0, gy0


manifest = []
for sid, d in cells.items():
    x = d["x"].to_numpy(float); y = d["y"].to_numpy(float)
    pheno = d["pheno"].to_numpy()

    X0 = float(np.floor(x.min() - MARGIN_UM))
    Y0 = float(np.floor(y.min() - MARGIN_UM))
    w_px = int(np.ceil((x.max() + MARGIN_UM - X0) / UM_PER_PX))
    h_px = int(np.ceil((y.max() + MARGIN_UM - Y0) / UM_PER_PX))

    cols = ((x - X0) / UM_PER_PX).astype(int)
    rows = ((y - Y0) / UM_PER_PX).astype(int)

    print(f"\n    {sid}")
    print(f"       canvas   {w_px} x {h_px} px  at {UM_PER_PX} um/px")
    print(f"       origin   X0 {X0:.1f} um, Y0 {Y0:.1f} um")

    # ---- layer 1, phenotypes, THE ANNOTATION SURFACE ------------------------
    canvas = np.full((h_px, w_px, 3), hex_to_rgb(CANVAS_BACKGROUND), dtype=np.uint8)
    for p in DRAW_ORDER:
        m = pheno == p
        if not m.any():
            continue
        stamp(canvas, cols[m], rows[m], hex_to_rgb(PHENOTYPE_COLORS.get(p, "#999999")),
              CELL_RADIUS_PX)
    p_ph = os.path.join(OUT_DIR, f"{sid}_phenotypes.png")
    Image.fromarray(canvas).save(p_ph, optimize=True)
    print(f"       wrote    {p_ph}")

    # ---- layer 2, myeloid only, what the detector counts --------------------
    canvas2 = np.full((h_px, w_px, 3), hex_to_rgb(CANVAS_BACKGROUND), dtype=np.uint8)
    nm = ~np.isin(pheno, MYELOID_FOR_DETECTION)
    stamp(canvas2, cols[nm], rows[nm], hex_to_rgb(NONMYELOID_COLOR), CELL_RADIUS_PX)
    mm = np.isin(pheno, MYELOID_FOR_DETECTION)
    stamp(canvas2, cols[mm], rows[mm], hex_to_rgb(MYELOID_COLOR), MYELOID_RADIUS_PX)
    p_my = os.path.join(OUT_DIR, f"{sid}_myeloid.png")
    Image.fromarray(canvas2).save(p_my, optimize=True)
    print(f"       wrote    {p_my}  "
          f"({int(mm.sum()):,} myeloid of {len(d):,})")
    del canvas, canvas2
    gc.collect()

    # ---- layer 3, REFERENCE ONLY, density with current foci -----------------
    dens, occ, gx0, gy0 = density_surface(d, X0, Y0, w_px, h_px)
    bg = float(np.nanmedian(dens[occ]))
    fig, ax = plt.subplots(figsize=(w_px / 400.0, h_px / 400.0))
    ext = [gx0, gx0 + dens.shape[1] * GRID_UM,
           gy0 + dens.shape[0] * GRID_UM, gy0]
    im = ax.imshow(np.log10(np.clip(dens, 1.0, None)), extent=ext,
                   cmap="magma", interpolation="nearest")
    if foci_ref is not None and "sample_id" in foci_ref.columns:
        f = foci_ref.loc[foci_ref["sample_id"].astype(str) == sid]
        if len(f) and {"peak_row", "peak_col"} <= set(f.columns):
            fx = gx0 + f["peak_col"].to_numpy(float) * GRID_UM
            fy = gy0 + f["peak_row"].to_numpy(float) * GRID_UM
            ax.scatter(fx, fy, s=200, facecolors="none", edgecolors="#00FFFF",
                       linewidths=2.5)
    ax.set_aspect("equal"); ax.invert_yaxis()
    ax.set_xlim(X0, X0 + w_px * UM_PER_PX)
    ax.set_ylim(Y0 + h_px * UM_PER_PX, Y0)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{sid} REFERENCE ONLY - do not annotate this\n"
                 f"myeloid density, log10 cells/mm2.  background {bg:.0f}, "
                 f"6x gate {bg * FOCUS_FOLD_OVER_BACKGROUND:.0f}, "
                 f"burden {BURDEN_THRESHOLD:.0f}.  cyan = current foci",
                 fontsize=11)
    fig.colorbar(im, ax=ax, fraction=0.03)
    p_ref = os.path.join(OUT_DIR, f"{sid}_reference_density.png")
    fig.savefig(p_ref, dpi=200, bbox_inches="tight")
    plt.close(fig); gc.collect()
    print(f"       wrote    {p_ref}  (reference, background {bg:.0f}/mm2)")

    # ---- the transform ------------------------------------------------------
    tr = dict(
        sample_id=sid, um_per_px=UM_PER_PX,
        x0_um=X0, y0_um=Y0, width_px=w_px, height_px=h_px,
        n_cells=int(len(d)),
        x_min_um=float(x.min()), x_max_um=float(x.max()),
        y_min_um=float(y.min()), y_max_um=float(y.max()),
        formula="x_um = x0_um + (col + 0.5) * um_per_px ; "
                "y_um = y0_um + (row + 0.5) * um_per_px",
        coordinate_space="whole-scan microns, same frame as the cell "
                         "assignment CSV x and y columns",
        annotation_surface=[os.path.basename(p_ph), os.path.basename(p_my)],
        reference_only=os.path.basename(p_ref),
        density_background_median=bg)
    p_tr = os.path.join(OUT_DIR, f"{sid}_transform.json")
    with open(p_tr, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(tr, fh, indent=2)
    print(f"       wrote    {p_tr}")
    manifest.append(tr)
    del dens, occ
    gc.collect()

# ---- legend -----------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 8))
ax.axis("off")
handles = [Patch(facecolor=PHENOTYPE_COLORS[p], edgecolor="#000000", label=p)
           for p in reversed(DRAW_ORDER) if p in PHENOTYPE_COLORS]
ax.legend(handles=handles, loc="center", frameon=False, fontsize=14,
          title="Phenotype layer", title_fontsize=16)
p_leg = os.path.join(OUT_DIR, "legend_phenotypes.png")
fig.savefig(p_leg, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"\n    wrote {p_leg}")


# %% Cell 5 - generate the QuPath script with the transforms baked in
# =============================================================================

banner("GENERATING THE QUPATH EXPORT SCRIPT")

sec_lines = ",\n".join(
    f'    "{m["sample_id"]}": [{m["x0_um"]}, {m["y0_um"]}, {m["um_per_px"]}, '
    f'{m["width_px"]}, {m["height_px"]}]' for m in manifest)

groovy = f'''// AKOYA annotation export, generated by AKOYA_04c_Annotation_Canvas.py
// Generated {datetime.now().isoformat(timespec="seconds")}
//
// HOW TO USE
//   1. In QuPath: File > Open, and open ONE of the *_phenotypes.png or
//      *_myeloid.png canvases. Image type: Brightfield (other) is fine, or
//      Fluorescence. It does not matter, nothing is being segmented.
//   2. Set the pixel size to the canvas scale (Image tab, double-click Pixel
//      width and Pixel height). This is SAFE and recommended. QuPath stores
//      ROI coordinates in image pixels regardless of calibration, so setting
//      it changes the display and the measurement units only, never the
//      exported geometry. It makes the status bar read in real microns, which
//      matters for judging whether a region is 200 um or 2 mm across.
//      The script checks the calibration and warns if it is not the canvas
//      scale, but the export is correct either way.
//   3. Draw annotations with the polygon, brush or wand tool.
//      Classify each one: complex, missed, clean_edge, normal
//      (right-click > Set class, adding the class the first time)
//   4. Run this script. It writes the GeoJSON and a TSV carrying BOTH pixel
//      and true scan-micron coordinates, and prints a summary.
//   5. Repeat for the other section. Exports are named per image, so they
//      will not overwrite each other.

import qupath.lib.io.PathIO

def OUT = "{ANNOTATION_OUT_DIR}"
new File(OUT).mkdirs()

// sample_id : [x0_um, y0_um, um_per_px, width_px, height_px]
def TRANSFORMS = [
{sec_lines}
]

def server = getCurrentServer()
def name = server.getMetadata().getName()

// work out which section this canvas is, from the filename
def sid = TRANSFORMS.keySet().find {{ name.contains(it) }}
if (sid == null) {{
    println "ERROR: could not match the open image to a known section."
    println "  image name : " + name
    println "  known      : " + TRANSFORMS.keySet()
    println "Open one of the generated *_phenotypes.png or *_myeloid.png files."
    return
}}
def t = TRANSFORMS[sid]
def X0 = t[0], Y0 = t[1], SC = t[2], W = t[3], H = t[4]

def cal = server.getPixelCalibration()
println "section      : " + sid
println "image        : " + name
println "pixel size   : " + cal.getPixelWidthMicrons() + " um (display only)"
if (Math.abs(cal.getPixelWidthMicrons() - SC) > 1e-6) {{
    println "  NOTE: pixel size is not the canvas scale of " + SC + " um."
    println "  The export below is still correct, because QuPath stores ROI"
    println "  coordinates in image pixels regardless of calibration. Setting"
    println "  it to " + SC + " only makes the on-screen measurements read in"
    println "  real microns while you draw."
}}
println "image size   : " + server.getWidth() + " x " + server.getHeight() + " px"
println "expected     : " + W + " x " + H + " px"
if (server.getWidth() != W || server.getHeight() != H) {{
    println "  WARNING: size mismatch. The image may have been resized or"
    println "  cropped, which would break the transform. Stop and re-open the"
    println "  original PNG."
}}
println "transform    : x_um = " + X0 + " + (col + 0.5) * " + SC
println "               y_um = " + Y0 + " + (row + 0.5) * " + SC

def anns = getAnnotationObjects()
println "annotations  : " + anns.size()
if (anns.isEmpty()) {{
    println "Nothing to export. Draw some annotations first."
    return
}}

// GeoJSON, coordinates stay in canvas pixels; the python side applies the
// transform using the matching *_transform.json
def gj = new File(OUT, sid + "_annotations.geojson")
PathIO.exportObjectsAsGeoJSON(gj, anns, null)
println "wrote " + gj

// TSV carrying both frames, so the alignment can be verified without trusting
// any single file
def tsv = new File(OUT, sid + "_annotation_measurements.tsv")
tsv.withPrintWriter {{ w ->
    w.println("index\\tsample_id\\tclassification\\tcentroid_x_px\\tcentroid_y_px\\t" +
              "area_px2\\tcentroid_x_um\\tcentroid_y_um\\tarea_um2\\tn_points")
    anns.eachWithIndex {{ a, i ->
        def r = a.getROI()
        def xum = X0 + (r.getCentroidX() + 0.5) * SC
        def yum = Y0 + (r.getCentroidY() + 0.5) * SC
        def aum = r.getArea() * SC * SC
        def np = r.getAllPoints() ? r.getAllPoints().size() : 0
        w.println(String.format("%d\\t%s\\t%s\\t%.3f\\t%.3f\\t%.1f\\t%.3f\\t%.3f\\t%.1f\\t%d",
            i, sid, (a.getPathClass() ?: "unclassified").toString(),
            r.getCentroidX(), r.getCentroidY(), r.getArea(), xum, yum, aum, np))
        println String.format("  %-12s  px %7.1f, %7.1f   um %9.1f, %9.1f   area %9.0f um2",
            (a.getPathClass() ?: "unclassified").toString(),
            r.getCentroidX(), r.getCentroidY(), xum, yum, aum)
    }}
}}
println "wrote " + tsv
println ""
println "Upload both files plus this console output."
'''

p_gv = os.path.join(OUT_DIR, "QuPath_annotate_and_export.groovy")
with open(p_gv, "w", encoding="utf-8", newline="\n") as fh:
    fh.write(groovy)
print(f"    wrote {p_gv}")
print("    transforms baked in for: " + ", ".join(m["sample_id"] for m in manifest))


# %% Cell 6 - README and wrap up
# =============================================================================

banner("README")

readme = f"""AKOYA ANNOTATION CANVAS
generated {datetime.now().isoformat(timespec="seconds")} by AKOYA_04c_Annotation_Canvas.py

WHY THESE FILES EXIST

The raw Akoya imagery was never delivered, only the processed per-cell export.
So there is no qptiff to annotate. These PNGs render the cell positions at an
exact {UM_PER_PX} microns per pixel, chosen here rather than read from an image
header, which makes the mapping from a drawn polygon back to scan micron
coordinates exact arithmetic rather than an inference.

WHICH IMAGE TO DRAW ON

  <section>_phenotypes.png        every cell, coloured by phenotype
  <section>_myeloid.png           the four pooled myeloid phenotypes in black,
                                  everything else in grey. This is the exact
                                  population the detector counts.

Either is fine. The myeloid layer is closer to what the detector sees; the
phenotype layer carries more context for judging what a region is.

WHICH IMAGE NOT TO DRAW ON

  <section>_reference_density.png

This is the smoothed density surface with the CURRENT foci circled in cyan. It
is a derivative of the same myeloid cells the detector thresholds, so drawing
boundaries on it would mean calibrating a threshold against its own output.
Look at it AFTER drawing, to see where the detector and your judgment differ.

STEPS

1. QuPath: File > Open, choose one canvas PNG.
   Image type: anything. Nothing is being segmented.
2. Set the pixel size to {UM_PER_PX} um. Image tab, double-click Pixel width and
   Pixel height, enter {UM_PER_PX} for both.

   This is safe. QuPath keeps ROI coordinates in image pixels whatever the
   calibration, so this changes only what the status bar and the measurement
   table display. It does not affect the exported geometry. Do it, because
   judging whether a region is 200 um or 2 mm across by eye on an uncalibrated
   image is unnecessarily hard.

   (This corrects an earlier instruction to leave the image uncalibrated. That
   advice was over-cautious and made the drawing harder for no benefit.)
3. Draw with the polygon, brush or wand tool.
4. Classify each annotation, right-click > Set class:
      complex     a granuloma complex as you would count it
      missed      a region the detector treats as background that is lesion
      clean_edge  a focus with a sharp border, as the contrast case
      normal      parenchyma you would call uninvolved
   'normal' matters more than it looks: it tells us what background is by your
   judgment rather than by a percentile, which is the whole 36463 question.
5. Five to ten annotations per section is plenty. Two parameters are being
   calibrated, nothing is being trained.
6. Automate > Script editor, paste QuPath_annotate_and_export.groovy, Run.
   It works out which section is open from the filename, so the same script
   is used unchanged for every canvas.
7. Repeat for each remaining canvas. Exports are named per section and will
   not overwrite each other.

   Priority: G4_36463 and G4_43109 are the two the analysis is blocked on. The
   other four mainly supply clean_edge contrast cases and can follow later.

UPLOAD

  <section>_annotations.geojson
  <section>_annotation_measurements.tsv
  the console output, as text

THE TRANSFORM

  x_um = x0_um + (col + 0.5) * um_per_px
  y_um = y0_um + (row + 0.5) * um_per_px

Per section, from <section>_transform.json:
"""
for m in manifest:
    readme += (f"\n  {m['sample_id']}\n"
               f"    x0_um       {m['x0_um']}\n"
               f"    y0_um       {m['y0_um']}\n"
               f"    um_per_px   {m['um_per_px']}\n"
               f"    canvas      {m['width_px']} x {m['height_px']} px\n"
               f"    cells       {m['n_cells']:,}\n")

readme += "\n\nWHAT TO LOOK FOR ON EACH CANVAS\n\n"
for m in manifest:
    sid = m["sample_id"]
    ex = EXPERT_COUNTS.get(sid)
    real = EXPERT_IS_REAL.get(sid, False)
    readme += (f"  {sid}   background {m['density_background_median']:.0f}/mm2, "
               f"6x gate {m['density_background_median'] * FOCUS_FOLD_OVER_BACKGROUND:.0f}/mm2, "
               f"burden {BURDEN_THRESHOLD:.0f}/mm2\n"
               f"      target {ex} "
               f"({'expert' if real else 'algorithmic, not expert'})\n"
               f"      {SECTION_NOTES.get(sid, '')}\n\n")

readme += """A note on the two arms. The treated sections carry far fewer myeloid
cells, so their canvases look sparse next to the untreated ones. That is the
finding, not a rendering problem. Judge each section against its own tissue,
which is what the detector does too, and use the 'normal' class to record what
uninvolved parenchyma looks like in that particular section.
"""

readme += """
WHAT HAPPENS NEXT

The annotations are not used to define foci. The detector stays a pure myeloid
cell-density method. The annotations are used to MEASURE, on regions an expert
drew, the absolute density and fold at the boundary, what fraction of the
enclosed peak the boundary sits at (which tests the 50 percent growth rule),
how many local maxima fall inside one expert region (which quantifies
over-segmentation and gives a merge distance), and the density gradient across
fuzzy against clean edges.

That replaces an objective built on six expert counts, which was three degrees
of freedom tuning three parameters, with region overlap measured over hundreds
of pixels.
"""

p_rm = os.path.join(OUT_DIR, "README_annotation.txt")
with open(p_rm, "w", encoding="utf-8", newline="\n") as fh:
    fh.write(readme)
print(f"    wrote {p_rm}")

# ---- overview, all sections at one common scale -----------------------------
if len(manifest) > 1:
    n = len(manifest)
    ncol = 3
    nrow = int(np.ceil(n / ncol))
    # a common micron extent so section sizes are directly comparable by eye
    max_w = max(m["width_px"] * m["um_per_px"] for m in manifest)
    max_h = max(m["height_px"] * m["um_per_px"] for m in manifest)
    fig, axes = plt.subplots(nrow, ncol, figsize=(7 * ncol, 6 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for ax, m in zip(axes, manifest):
        sid = m["sample_id"]
        thumb = Image.open(os.path.join(OUT_DIR, f"{sid}_myeloid.png"))
        thumb.thumbnail((1400, 1400))
        ax.imshow(np.asarray(thumb),
                  extent=[0, m["width_px"] * m["um_per_px"],
                          m["height_px"] * m["um_per_px"], 0])
        ax.set_xlim(0, max_w); ax.set_ylim(max_h, 0)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ex = EXPERT_COUNTS.get(sid)
        real = EXPERT_IS_REAL.get(sid, False)
        ax.set_title(f"{sid}\nbg {m['density_background_median']:.0f}/mm2, "
                     f"target {ex}{'' if real else ' (algorithmic)'}",
                     fontsize=13)
        thumb.close()
    for ax in axes[len(manifest):]:
        ax.axis("off")
    fig.suptitle("All sections, myeloid layer, drawn to a common micron scale\n"
                 "reference only, annotate the full-size canvases", fontsize=16)
    fig.tight_layout()
    p_ov = os.path.join(OUT_DIR, "overview_all_sections.png")
    fig.savefig(p_ov, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    wrote {p_ov}")


banner("SUMMARY")
total_mb = 0.0
for f in glob.glob(os.path.join(OUT_DIR, "*")):
    total_mb += os.path.getsize(f) / 1e6
for m in manifest:
    print(f"    {m['sample_id']:<14} canvas {m['width_px']:>5} x {m['height_px']:>5} px   "
          f"origin ({m['x0_um']:.0f}, {m['y0_um']:.0f}) um   "
          f"{m['n_cells']:>8,} cells")

print(f"\n    total output {total_mb:.0f} MB across {len(manifest)} sections")

sub("Do next")
print("  1. Copy the annotation_canvas folder to wherever QuPath is running.")
print("  2. Read README_annotation.txt.")
print(f"  3. Annotate each of the {len(manifest)} canvases, running the generated")
print("     Groovy script once per canvas. Exports are named per section, so")
print("     they will not overwrite each other.")
print("  4. Upload every geojson, every tsv, and the console output.")
print("\n  There is no need to do all six in one sitting. 36463 and 43109 are")
print("  the ones the analysis is blocked on. The other four mainly supply")
print("  clean_edge contrast cases and can follow later.")

sub("Not changed by this script")
print("  Nothing. This is a rendering step. structures_rev4, script 04 and")
print("  every downstream script are untouched.")

banner("DONE")
sys.stdout = _tee.terminal
_tee.close()
