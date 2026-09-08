#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler inventory - FIGURE SUPPLEMENT
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 01b of the AKOYA analysis series. REVISION 2.

PURPOSE
    Turn the tables written by AKOYA_01_Inventory.py into presentation-ready
    figures for review with Deepak's group. Reads ONLY the inventory tables,
    never the raw QuPath CSVs, so the figures cannot disagree with the numbers
    already generated. Runs in seconds.

WHAT CHANGED IN REVISION 2 (and why)
    1. THE SLIDE BOUNDARY IS NOW DRAWN FROM scan_id, NOT FROM TREATMENT.
       F04 and F08 drew a black line at n_first - 0.5, where n_first counted the
       sections in the first condition. F08's title called that line the slide
       boundary. It landed in the right place only because scan and treatment
       are perfectly confounded in this cohort. It now comes from the scan_id
       column written by script 01 revision 2, so the line means what the
       caption says and it stays correct if the section set ever changes.

    2. SCAN AND SECTION POSITION APPEAR ON EVERY AXIS.
       The position-1 exclusion argument for G3_43102 and G4_43112 was invisible
       in every deck figure. Axis labels now carry scan_id and position rank, so
       the two candidate sections are identifiable on sight in F01 through F09.

    3. F08 SPLITS DRIFT INTO ITS BETWEEN-SCAN AND WITHIN-SCAN PARTS.
       The heatmap is unchanged. The ranking CSV now carries three quantities
       per marker rather than one: total spread across all sections, spread
       between the two scan medians, and the worst spread within a single scan.
       This is descriptive only. It does not decompose variance and it does not
       attribute anything to batch, because scan and arm are the same factor
       here. It exists so that script 02 is read with the right expectation
       about which markers are section-level and which are scan-level.

    4. F10 IS NEW AND OPTIONAL (RUN_LINEAGE_RATIOS).
       Within-lineage ratios that remove section composition the same way F06
       does for IDO1. It exists because the inventory counts show CD4+ running
       at 2.4 to 12.4 percent of T cells on scan_01 against 8.1 to 39.5 percent
       on scan_02, with 43106 calling 75 helper T cells in 143,790. Downstream
       BALT detection gates on helper T density, so whether that call is usable
       has to be visible before anything is built on it. Set
       RUN_LINEAGE_RATIOS = False to skip.

    5. rng AND THE SCAN LOOKUPS ARE BUILT IN CELL 3.
       rng was created inside the F05 cell and reused by F06, and n_first was
       created in the F04 cell and reused by F08. Both worked when the file ran
       top to bottom and both raised NameError when cells were run out of order
       in Spyder. They are now defined once, up front.

    Nothing about the numbers changed. No figure lost a panel.

WHAT IT MAKES (each saved as both PDF and PNG at 300 DPI)
    F01_cells_per_animal              total cells per section, by condition
    F02_tissue_area_and_density       bbox area and cell density per section
    F03_composition_stacked           phenotype composition, percent of section
    F04_composition_heatmap           log10 cell counts, phenotype x animal
    F05_phenotype_by_condition        per-phenotype percent, each animal a point
    F06_ido1_fraction_of_cd68         IDO1+ as a fraction of all CD68+ macrophages
    F07_immune_vs_structural          immune / structural / unassigned split
    F08_marker_p99_heatmap            staining intensity drift across sections
    F09_nn_feasibility                cells available per phenotype per animal
    F10_lineage_internal_ratios       CD4 fraction of T cells, B fraction of B lineage

INPUTS (from AKOYA_01_Inventory.py revision 2)
    tables/01_file_summary.csv        now carries scan_id and slide_position_rank
    tables/03_phenotype_labels_repr.csv
    tables/04_phenotype_counts_long.csv
    tables/06_marker_stats_long.csv

USAGE
    conda activate sc_pre
    python AKOYA_01b_Inventory_Figures.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================
# ALL TUNABLE PARAMETERS LIVE HERE
# =============================================================================

TABLE_DIR = "/master/jlehle/WORKING/AKOYA/inventory/tables"
FIG_DIR = "/master/jlehle/WORKING/AKOYA/inventory/figures"

# Condition ordering and colors
CONDITION_ORDER = ["D1MT", "Untreated"]
CONDITION_COLORS = {
    "D1MT": "#2C7FB8",
    "Untreated": "#D95F02",
}

# ---- scan and position annotation -------------------------------------------
# Written by script 01 revision 2. If these columns are absent the script warns
# and falls back to the old treatment-based boundary so it still runs.
SCAN_COL = "scan_id"
POSITION_COL = "slide_position_rank"
SHOW_SCAN_IN_TICKS = True        # add "scan_01 pos1" as a third axis label line
FLAG_POSITION_RANK = 1           # position rank to mark on every axis
SCAN_BOUNDARY_COLOR = "#000000"
SCAN_BOUNDARY_WIDTH = 4.0

# ---- optional diagnostic ----------------------------------------------------
RUN_LINEAGE_RATIOS = True        # F10

# Phenotype display order, coarse to fine. Labels here are the ASCII-safe forms;
# the script maps them back to the exact strings in the data, including the
# zero-width space in the epithelial label.
PHENOTYPE_ORDER = [
    "CD68+IDO1+ Macrophages",
    "CD68+IDO1- Macrophages",
    "CD163+ Macrophages",
    "Neutrophils",
    "Helper T cells",
    "CD4- T cells",
    "Tregs",
    "B cells",
    "Plasma cells",
    "Endothelial cells",
    "Epithelial/Tumor cells",
    "Other",
]

PHENOTYPE_COLORS = {
    "CD68+IDO1+ Macrophages": "#B2182B",
    "CD68+IDO1- Macrophages": "#EF8A62",
    "CD163+ Macrophages": "#FDDBC7",
    "Neutrophils": "#7B3294",
    "Helper T cells": "#1B7837",
    "CD4- T cells": "#7FBC41",
    "Tregs": "#276419",
    "B cells": "#2166AC",
    "Plasma cells": "#67A9CF",
    "Endothelial cells": "#BABABA",
    "Epithelial/Tumor cells": "#878787",
    "Other": "#4D4D4D",
}

# Functional grouping for the immune / structural summary
IMMUNE_PHENOTYPES = [
    "CD68+IDO1+ Macrophages", "CD68+IDO1- Macrophages", "CD163+ Macrophages",
    "Neutrophils", "Helper T cells", "CD4- T cells", "Tregs",
    "B cells", "Plasma cells",
]
STRUCTURAL_PHENOTYPES = ["Endothelial cells", "Epithelial/Tumor cells"]
UNASSIGNED_PHENOTYPES = ["Other"]

# Macrophage phenotypes used as the denominator for the IDO1 fraction.
# CD163+ is deliberately excluded: the call is CD68-based.
CD68_NUMERATOR = "CD68+IDO1+ Macrophages"
CD68_DENOMINATOR = ["CD68+IDO1+ Macrophages", "CD68+IDO1- Macrophages"]

# ---- lineage internal ratios for F10 ----------------------------------------
# Each entry is (panel title, numerator phenotypes, denominator phenotypes,
#                y axis label, why it matters)
LINEAGE_RATIOS = [
    ("CD4+ fraction of T cells",
     ["Helper T cells", "Tregs"],
     ["Helper T cells", "Tregs", "CD4- T cells"],
     "CD4+ % of all T cells",
     "BALT detection gates on helper T density"),
    ("B cell fraction of B lineage",
     ["B cells"],
     ["B cells", "Plasma cells"],
     "B cells % of B lineage",
     "Finding 3 is a B lineage claim"),
    ("CD163+ fraction of all macrophages",
     ["CD163+ Macrophages"],
     ["CD163+ Macrophages", "CD68+IDO1+ Macrophages", "CD68+IDO1- Macrophages"],
     "CD163+ % of macrophages",
     "CD163 is in the myeloid detection pool"),
]

# Markers highlighted in the drift heatmap because they define phenotype calls
KEY_MARKERS = [
    "IDO1", "CD68", "CD163", "CD206", "CD4", "CD8", "CD3e", "CD20", "CD79a",
    "MPO", "CD11b", "FoxP3", "iNOS", "Arginase-1", "IFNG", "Granzyme-B",
    "HLA-DR", "CD21", "PD-1", "PD-L1",
]
N_MARKERS_IN_HEATMAP = 40   # top N by cross-section spread, key markers forced in

# Threshold line drawn on the feasibility plot
MIN_CELLS_FOR_NN = 50

RANDOM_SEED = 0

# ---- figure style -----------------------------------------------------------
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
ZERO_FLAG_COLOR = "#B2182B"
POSITION_FLAG_COLOR = "#B2182B"


# %% Cell 2 - imports, style, helpers
# =============================================================================

import os
import sys
import unicodedata
from datetime import datetime

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

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
    "axes.grid": False,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def banner(text, char="=", width=88):
    print("\n" + char * width)
    print(text)
    print(char * width)


def ascii_safe(s):
    """Strip zero-width and other non-ASCII characters for display."""
    if not isinstance(s, str):
        return str(s)
    return "".join(c for c in s if ord(c) <= 127).strip()


def load_table(name, required=True):
    path = os.path.join(TABLE_DIR, name)
    if not os.path.exists(path):
        msg = f"    {'ERROR' if required else 'WARNING'}: {path} not found"
        print(msg)
        if required:
            print("    Run AKOYA_01_Inventory.py first. Stopping.")
            sys.exit(1)
        return None
    df = pd.read_csv(path)
    print(f"    loaded {name}  ({df.shape[0]} rows x {df.shape[1]} cols)")
    return df


def save_fig(fig, stem):
    """Save as both PDF and PNG at 300 DPI."""
    saved = []
    if SAVE_PDF:
        p = os.path.join(FIG_DIR, f"{stem}.pdf")
        fig.savefig(p, dpi=DPI, bbox_inches="tight")
        saved.append(p)
    if SAVE_PNG:
        p = os.path.join(FIG_DIR, f"{stem}.png")
        fig.savefig(p, dpi=DPI, bbox_inches="tight")
        saved.append(p)
    plt.close(fig)
    for p in saved:
        print(f"    wrote {p}")


def style_axes(ax, ygrid=True):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(2)
    ax.spines["bottom"].set_linewidth(2)
    if ygrid:
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, color=GRID_COLOR, linewidth=1.5)


def condition_of(sample_id, meta_map):
    return meta_map.get(sample_id, {}).get("treatment", "UNKNOWN")


def short_label(sample_id):
    """G3_43102 -> 43102 for compact axis labels."""
    return sample_id.split("_")[-1] if "_" in sample_id else sample_id


# %% Cell 3 - load tables, build ordering, scan lookups
# =============================================================================

os.makedirs(FIG_DIR, exist_ok=True)

banner("AKOYA INVENTORY FIGURES (revision 2)")
print(f"Run time  : {datetime.now().isoformat(timespec='seconds')}")
print(f"Tables    : {TABLE_DIR}")
print(f"Figures   : {FIG_DIR}\n")

file_summary = load_table("01_file_summary.csv", required=True)
pheno_long = load_table("04_phenotype_counts_long.csv", required=True)
marker_stats = load_table("06_marker_stats_long.csv", required=False)
label_tbl = load_table("03_phenotype_labels_repr.csv", required=False)

rng = np.random.default_rng(RANDOM_SEED)

# ---- animal ordering: condition first, then animal ID ------------------------
file_summary["treatment"] = file_summary["treatment"].astype(str)
cond_rank = {c: i for i, c in enumerate(CONDITION_ORDER)}
file_summary["_rank"] = file_summary["treatment"].map(
    lambda t: cond_rank.get(t, len(cond_rank))
)
file_summary = file_summary.sort_values(["_rank", "animal_id"]).reset_index(drop=True)
SAMPLE_ORDER = file_summary["sample_id"].tolist()
meta_map = file_summary.set_index("sample_id").to_dict(orient="index")

# ---- scan identity and section position -------------------------------------
HAVE_SCAN = SCAN_COL in file_summary.columns
HAVE_POSITION = POSITION_COL in file_summary.columns
if not HAVE_SCAN:
    print(f"\n    WARNING: '{SCAN_COL}' not found in 01_file_summary.csv. This "
          f"script expects script 01 revision 2 or later.")
    print(f"    Falling back to the treatment-based boundary, which is only "
          f"correct while scan and arm are perfectly confounded.")
if not HAVE_POSITION:
    print(f"    WARNING: '{POSITION_COL}' not found. Position annotation is off.")

SCAN_OF = {s: (str(meta_map[s][SCAN_COL]) if HAVE_SCAN else
               f"scan_{condition_of(s, meta_map)}") for s in SAMPLE_ORDER}
POS_OF = {}
for s in SAMPLE_ORDER:
    v = meta_map[s].get(POSITION_COL, np.nan) if HAVE_POSITION else np.nan
    POS_OF[s] = int(v) if pd.notna(v) else None

print(f"\nSection order used in all figures:")
for sid in SAMPLE_ORDER:
    m = meta_map[sid]
    pos = POS_OF[sid]
    flag = "   <-- position 1" if pos == FLAG_POSITION_RANK else ""
    print(f"    {sid:<12} {m['treatment']:<10} {SCAN_OF[sid]:<10} "
          f"pos={pos if pos is not None else 'na':<4} "
          f"{m['n_cells']:>9,} cells{flag}")

# ---- boundaries between scans, in plotting order ----------------------------
SCAN_BOUNDARIES = [i - 0.5 for i in range(1, len(SAMPLE_ORDER))
                   if SCAN_OF[SAMPLE_ORDER[i]] != SCAN_OF[SAMPLE_ORDER[i - 1]]]
print(f"\n    scan boundaries in plotting order: {SCAN_BOUNDARIES}")


def draw_scan_boundaries(ax):
    """One vertical rule wherever the scan changes along the x axis."""
    for b in SCAN_BOUNDARIES:
        ax.axvline(b, color=SCAN_BOUNDARY_COLOR, linewidth=SCAN_BOUNDARY_WIDTH,
                   zorder=6)


# ---- phenotype label reconciliation -----------------------------------------
# Map exact strings in the data (which may carry a zero-width space) to the
# ASCII-safe display names used in PHENOTYPE_ORDER.
pheno_long["phenotype_display"] = pheno_long["phenotype"].map(ascii_safe)

observed = sorted(pheno_long["phenotype_display"].unique())
missing_from_order = [p for p in observed if p not in PHENOTYPE_ORDER]
missing_from_data = [p for p in PHENOTYPE_ORDER if p not in observed]
if missing_from_order:
    print(f"\n    WARNING: phenotypes in the data but not in PHENOTYPE_ORDER, "
          f"appended at the end: {missing_from_order}")
if missing_from_data:
    print(f"    WARNING: phenotypes in PHENOTYPE_ORDER but not in the data: "
          f"{missing_from_data}")

PHENO_USE = [p for p in PHENOTYPE_ORDER if p in observed] + missing_from_order
for p in missing_from_order:
    PHENOTYPE_COLORS.setdefault(p, "#999999")

# ---- dense count and percent matrices (zeros filled in explicitly) ----------
counts = (pheno_long.pivot_table(index="phenotype_display", columns="sample_id",
                                 values="n_cells", aggfunc="sum")
          .reindex(index=PHENO_USE, columns=SAMPLE_ORDER)
          .fillna(0.0))
totals = counts.sum(axis=0)
pcts = 100.0 * counts / totals

n_zero_filled = int((counts == 0).sum().sum())
print(f"\n    count matrix: {counts.shape[0]} phenotypes x {counts.shape[1]} sections, "
      f"{n_zero_filled} zero cell(s) filled explicitly")

cond_colors = [CONDITION_COLORS.get(condition_of(s, meta_map), "#999999")
               for s in SAMPLE_ORDER]
cond_labels = [condition_of(s, meta_map) for s in SAMPLE_ORDER]
xpos = np.arange(len(SAMPLE_ORDER))


def tick_labels(n_lines=2):
    """
    n_lines = 1 : 43102
    n_lines = 2 : 43102 / D1MT
    n_lines = 3 : 43102 / D1MT / scan_01 pos1
    The third line is suppressed when SHOW_SCAN_IN_TICKS is False.
    """
    out = []
    for s in SAMPLE_ORDER:
        parts = [short_label(s)]
        if n_lines >= 2:
            parts.append(condition_of(s, meta_map))
        if n_lines >= 3 and SHOW_SCAN_IN_TICKS:
            pos = POS_OF[s]
            pos_txt = f"pos{pos}" if pos is not None else "pos na"
            parts.append(f"{SCAN_OF[s]} {pos_txt}")
        out.append("\n".join(parts))
    return out


x_labels = tick_labels(1)
x_labels2 = tick_labels(2)
x_labels3 = tick_labels(3)


def color_position_ticks(ax):
    """Colour the tick label of any position-1 section so it is visible."""
    if not HAVE_POSITION:
        return
    for tick, s in zip(ax.get_xticklabels(), SAMPLE_ORDER):
        if POS_OF[s] == FLAG_POSITION_RANK:
            tick.set_color(POSITION_FLAG_COLOR)
            tick.set_fontweight("bold")


def add_condition_legend(ax, loc="upper right"):
    handles = [Patch(facecolor=CONDITION_COLORS[c], edgecolor="none", label=c)
               for c in CONDITION_ORDER if c in cond_labels]
    ax.legend(handles=handles, loc=loc, frameon=False, fontsize=FONT_SIZE_LEGEND)


def add_condition_bands(ax):
    """Shade the background behind each condition block."""
    for c in set(cond_labels):
        idx = [i for i, lab in enumerate(cond_labels) if lab == c]
        if not idx:
            continue
        ax.axvspan(min(idx) - 0.5, max(idx) + 0.5,
                   color=CONDITION_COLORS.get(c, "#999999"), alpha=0.06, zorder=0)


# %% Cell 4 - F01 cells per animal, F02 area and density
# =============================================================================

banner("F01 / F02 - SECTION SIZE AND CELL YIELD")

fig, ax = plt.subplots(figsize=(18, 12))
bars = ax.bar(xpos, [meta_map[s]["n_cells"] for s in SAMPLE_ORDER],
              color=cond_colors, edgecolor="#FFFFFF", linewidth=2, zorder=3)
for i, s in enumerate(SAMPLE_ORDER):
    ax.text(i, meta_map[s]["n_cells"] * 1.02, f"{meta_map[s]['n_cells']:,}",
            ha="center", va="bottom", fontsize=FONT_SIZE_ANNOT, rotation=0)
ax.set_xticks(xpos)
ax.set_xticklabels(x_labels3, fontsize=FONT_SIZE_TICK - 8)
color_position_ticks(ax)
ax.set_ylabel("Segmented cells")
ax.set_xlabel("Animal")
ax.set_title("Cells recovered per lung section\n"
             "(red label = position 1 on its scan)",
             fontsize=FONT_SIZE_TITLE - 4)
ax.set_ylim(0, max(meta_map[s]["n_cells"] for s in SAMPLE_ORDER) * 1.15)
ax.yaxis.set_major_formatter(
    matplotlib.ticker.FuncFormatter(lambda v, _: f"{v/1000:,.0f}k"))
style_axes(ax)
add_condition_bands(ax)
draw_scan_boundaries(ax)
add_condition_legend(ax)
save_fig(fig, "F01_cells_per_animal")

# ---- F02: two panels, tissue footprint and density --------------------------
fig, axes = plt.subplots(1, 2, figsize=(26, 12))

ax = axes[0]
vals = [meta_map[s].get("bbox_area_mm2", np.nan) for s in SAMPLE_ORDER]
ax.bar(xpos, vals, color=cond_colors, edgecolor="#FFFFFF", linewidth=2, zorder=3)
ax.set_xticks(xpos)
ax.set_xticklabels(x_labels3, rotation=45, ha="right", fontsize=FONT_SIZE_TICK - 10)
color_position_ticks(ax)
ax.set_ylabel("Bounding box area (mm$^2$)")
ax.set_title("Section footprint")
style_axes(ax); add_condition_bands(ax); draw_scan_boundaries(ax)

ax = axes[1]
vals = [meta_map[s].get("cells_per_mm2_bbox", np.nan) for s in SAMPLE_ORDER]
ax.bar(xpos, vals, color=cond_colors, edgecolor="#FFFFFF", linewidth=2, zorder=3)
ax.set_xticks(xpos)
ax.set_xticklabels(x_labels3, rotation=45, ha="right", fontsize=FONT_SIZE_TICK - 10)
color_position_ticks(ax)
ax.set_ylabel("Cells / mm$^2$")
ax.set_title("Cell density (bounding box)")
style_axes(ax); add_condition_bands(ax); draw_scan_boundaries(ax)
add_condition_legend(ax)

fig.suptitle("Tissue footprint and cell density per section", y=1.02)
save_fig(fig, "F02_tissue_area_and_density")


# %% Cell 5 - F03 stacked composition, F04 count heatmap
# =============================================================================

banner("F03 / F04 - PHENOTYPE COMPOSITION")

fig, ax = plt.subplots(figsize=(20, 14))
bottom = np.zeros(len(SAMPLE_ORDER))
for ph in PHENO_USE:
    vals = pcts.loc[ph, SAMPLE_ORDER].to_numpy(dtype=float)
    ax.bar(xpos, vals, bottom=bottom, color=PHENOTYPE_COLORS.get(ph, "#999999"),
           edgecolor="#FFFFFF", linewidth=1.5, label=ph, zorder=3)
    bottom += vals
ax.set_xticks(xpos)
ax.set_xticklabels(x_labels3, fontsize=FONT_SIZE_TICK - 10)
color_position_ticks(ax)
ax.set_ylabel("Percent of section")
ax.set_ylim(0, 100)
ax.set_title("Phenotype composition per section")
style_axes(ax, ygrid=False)
draw_scan_boundaries(ax)
ax.legend(bbox_to_anchor=(1.01, 1.0), loc="upper left", frameon=False,
          fontsize=FONT_SIZE_LEGEND - 6)
save_fig(fig, "F03_composition_stacked")

# ---- F04: log10 count heatmap ----------------------------------------------
fig, ax = plt.subplots(figsize=(20, 15))
mat = np.log10(counts.loc[PHENO_USE, SAMPLE_ORDER].to_numpy(dtype=float) + 1.0)
cmap = LinearSegmentedColormap.from_list(
    "akoya_counts", ["#FFFFFF", "#C6DBEF", "#4292C6", "#08519C", "#08306B"])
im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=np.nanmax(mat))
ax.set_xticks(xpos)
ax.set_xticklabels(x_labels3, fontsize=FONT_SIZE_TICK - 10)
color_position_ticks(ax)
ax.set_yticks(np.arange(len(PHENO_USE)))
ax.set_yticklabels(PHENO_USE, fontsize=FONT_SIZE_TICK - 4)
for i, ph in enumerate(PHENO_USE):
    for j, s in enumerate(SAMPLE_ORDER):
        v = counts.loc[ph, s]
        txt = "0" if v == 0 else f"{int(v):,}"
        color = ZERO_FLAG_COLOR if v == 0 else (
            "#FFFFFF" if mat[i, j] > 0.65 * np.nanmax(mat) else "#000000")
        weight = "bold" if v == 0 else "normal"
        ax.text(j, i, txt, ha="center", va="center",
                fontsize=FONT_SIZE_ANNOT - 8, color=color, fontweight=weight)
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
cb.set_label("log$_{10}$(cells + 1)", fontsize=FONT_SIZE_BASE - 4)
cb.ax.tick_params(labelsize=FONT_SIZE_TICK - 8)
draw_scan_boundaries(ax)
ax.set_title("Cell counts per phenotype and section\n"
             "(red zeros = phenotype absent; black rule = scan boundary)",
             fontsize=FONT_SIZE_TITLE - 4)
save_fig(fig, "F04_composition_heatmap")


# %% Cell 6 - F05 per-phenotype percent by condition
# =============================================================================

banner("F05 - PHENOTYPE FREQUENCY BY CONDITION")

n_ph = len(PHENO_USE)
ncol = 4
nrow = int(np.ceil(n_ph / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(9 * ncol, 7.5 * nrow))
axes = np.atleast_1d(axes).ravel()

for k, ph in enumerate(PHENO_USE):
    ax = axes[k]
    for ci, cond in enumerate(CONDITION_ORDER):
        samples = [s for s in SAMPLE_ORDER if condition_of(s, meta_map) == cond]
        if not samples:
            continue
        vals = pcts.loc[ph, samples].to_numpy(dtype=float)
        jitter = rng.uniform(-0.10, 0.10, size=len(vals))
        ax.scatter(np.full(len(vals), ci) + jitter, vals,
                   s=420, color=CONDITION_COLORS[cond],
                   edgecolor="#FFFFFF", linewidth=2.5, zorder=3)
        med = float(np.median(vals))
        ax.hlines(med, ci - 0.28, ci + 0.28, color="#000000",
                  linewidth=4, zorder=4)
        for v, j, s in zip(vals, jitter, samples):
            if v == 0:
                ax.scatter([ci + j], [0], s=520, marker="x",
                           color=ZERO_FLAG_COLOR, linewidth=4, zorder=5)
            if POS_OF.get(s) == FLAG_POSITION_RANK:
                ax.scatter([ci + j], [v], s=760, facecolor="none",
                           edgecolor=POSITION_FLAG_COLOR, linewidth=3, zorder=6)
    ax.set_xticks(range(len(CONDITION_ORDER)))
    ax.set_xticklabels(CONDITION_ORDER, fontsize=FONT_SIZE_TICK - 4)
    ax.set_xlim(-0.6, len(CONDITION_ORDER) - 0.4)
    ax.set_title(ph, fontsize=FONT_SIZE_BASE - 2)
    ax.set_ylabel("% of section", fontsize=FONT_SIZE_BASE - 6)
    ax.tick_params(axis="y", labelsize=FONT_SIZE_TICK - 8)
    ax.set_ylim(bottom=0)
    style_axes(ax)

for k in range(n_ph, len(axes)):
    axes[k].axis("off")

fig.suptitle("Phenotype frequency by condition (each point = one animal, "
             "bar = group median, red ring = position-1 section)",
             y=1.005, fontsize=FONT_SIZE_TITLE - 4)
fig.tight_layout()
save_fig(fig, "F05_phenotype_by_condition")


# %% Cell 7 - F06 IDO1+ fraction of CD68+ macrophages
# =============================================================================

banner("F06 - IDO1+ FRACTION OF CD68+ MACROPHAGES")

have_num = CD68_NUMERATOR in counts.index
have_den = all(p in counts.index for p in CD68_DENOMINATOR)

if not (have_num and have_den):
    print("    WARNING: CD68 macrophage phenotypes not found. Skipping F06.")
    ido_tbl = pd.DataFrame()
else:
    num = counts.loc[CD68_NUMERATOR, SAMPLE_ORDER].to_numpy(dtype=float)
    den = counts.loc[CD68_DENOMINATOR, SAMPLE_ORDER].sum(axis=0).to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.where(den > 0, 100.0 * num / den, np.nan)

    ido_tbl = pd.DataFrame({
        "sample_id": SAMPLE_ORDER,
        "animal_id": [short_label(s) for s in SAMPLE_ORDER],
        "treatment": cond_labels,
        "scan_id": [SCAN_OF[s] for s in SAMPLE_ORDER],
        "slide_position_rank": [POS_OF[s] for s in SAMPLE_ORDER],
        "n_ido1_pos": num.astype(int),
        "n_cd68_total": den.astype(int),
        "pct_ido1_pos_of_cd68": frac,
    })
    print(ido_tbl.to_string(index=False))

    fig, axes = plt.subplots(1, 2, figsize=(28, 13),
                             gridspec_kw={"width_ratios": [1.6, 1.0]})

    # left: per animal
    ax = axes[0]
    ax.bar(xpos, np.nan_to_num(frac), color=cond_colors,
           edgecolor="#FFFFFF", linewidth=2, zorder=3)
    for i in range(len(SAMPLE_ORDER)):
        if np.isnan(frac[i]):
            continue
        lab = f"{frac[i]:.1f}%\n{int(num[i]):,}/{int(den[i]):,}"
        ax.text(i, frac[i] + 2.5, lab, ha="center", va="bottom",
                fontsize=FONT_SIZE_ANNOT - 6,
                color=ZERO_FLAG_COLOR if num[i] == 0 else TEXT_COLOR,
                fontweight="bold" if num[i] == 0 else "normal")
    ax.set_xticks(xpos)
    ax.set_xticklabels(x_labels3, fontsize=FONT_SIZE_TICK - 10)
    color_position_ticks(ax)
    ax.set_ylabel("IDO1+ % of CD68+ macrophages")
    ax.set_ylim(0, 110)
    ax.set_title("Per animal")
    style_axes(ax); add_condition_bands(ax); draw_scan_boundaries(ax)

    # right: by condition
    ax = axes[1]
    for ci, cond in enumerate(CONDITION_ORDER):
        idx = [i for i, c in enumerate(cond_labels) if c == cond]
        if not idx:
            continue
        vals = frac[idx]
        jitter = rng.uniform(-0.10, 0.10, size=len(vals))
        ax.scatter(np.full(len(vals), ci) + jitter, np.nan_to_num(vals),
                   s=620, color=CONDITION_COLORS[cond],
                   edgecolor="#FFFFFF", linewidth=3, zorder=3)
        for v, j, i0 in zip(np.nan_to_num(vals), jitter, idx):
            if POS_OF.get(SAMPLE_ORDER[i0]) == FLAG_POSITION_RANK:
                ax.scatter([ci + j], [v], s=1000, facecolor="none",
                           edgecolor=POSITION_FLAG_COLOR, linewidth=4, zorder=5)
        ax.hlines(np.nanmedian(vals), ci - 0.28, ci + 0.28,
                  color="#000000", linewidth=5, zorder=4)
    ax.set_xticks(range(len(CONDITION_ORDER)))
    ax.set_xticklabels(CONDITION_ORDER)
    ax.set_xlim(-0.6, len(CONDITION_ORDER) - 0.4)
    ax.set_ylim(0, 110)
    ax.set_ylabel("IDO1+ % of CD68+ macrophages")
    ax.set_title("By condition")
    style_axes(ax)

    fig.suptitle("IDO1+ macrophages as a fraction of all CD68+ macrophages\n"
                 "(internal normalization, removes section composition effects; "
                 "red ring = position-1 section)",
                 y=1.05, fontsize=FONT_SIZE_TITLE - 4)
    save_fig(fig, "F06_ido1_fraction_of_cd68")


# %% Cell 8 - F07 immune vs structural vs unassigned
# =============================================================================

banner("F07 - IMMUNE / STRUCTURAL / UNASSIGNED")

groups = {
    "Immune": [p for p in IMMUNE_PHENOTYPES if p in counts.index],
    "Structural": [p for p in STRUCTURAL_PHENOTYPES if p in counts.index],
    "Unassigned": [p for p in UNASSIGNED_PHENOTYPES if p in counts.index],
}
group_colors = {"Immune": "#1B7837", "Structural": "#878787", "Unassigned": "#4D4D4D"}

roll = pd.DataFrame(
    {g: counts.loc[ps, SAMPLE_ORDER].sum(axis=0) for g, ps in groups.items()}
).T
roll_pct = 100.0 * roll / roll.sum(axis=0)

fig, ax = plt.subplots(figsize=(20, 13))
bottom = np.zeros(len(SAMPLE_ORDER))
for g in ["Immune", "Structural", "Unassigned"]:
    vals = roll_pct.loc[g, SAMPLE_ORDER].to_numpy(dtype=float)
    ax.bar(xpos, vals, bottom=bottom, color=group_colors[g],
           edgecolor="#FFFFFF", linewidth=2, label=g, zorder=3)
    for i, v in enumerate(vals):
        if v > 6:
            ax.text(i, bottom[i] + v / 2, f"{v:.0f}%", ha="center", va="center",
                    color="#FFFFFF", fontsize=FONT_SIZE_ANNOT - 4)
    bottom += vals
ax.set_xticks(xpos)
ax.set_xticklabels(x_labels3, fontsize=FONT_SIZE_TICK - 10)
color_position_ticks(ax)
ax.set_ylabel("Percent of section")
ax.set_ylim(0, 100)
ax.set_title("Immune, structural and unassigned fractions")
style_axes(ax, ygrid=False)
draw_scan_boundaries(ax)
ax.legend(bbox_to_anchor=(1.01, 1.0), loc="upper left", frameon=False)
save_fig(fig, "F07_immune_vs_structural")


# %% Cell 9 - F08 marker intensity drift across sections
# =============================================================================

banner("F08 - STAINING INTENSITY DRIFT")

if marker_stats is None or "p99" not in marker_stats.columns:
    print("    WARNING: marker stats unavailable. Skipping F08.")
else:
    p99 = marker_stats.pivot_table(index="column", columns="sample_id", values="p99")
    p99 = p99.reindex(columns=SAMPLE_ORDER)

    # log2 deviation from each marker's own median across sections
    med = p99.median(axis=1).replace(0, np.nan)
    dev = np.log2(p99.div(med, axis=0).replace(0, np.nan))
    spread = dev.max(axis=1) - dev.min(axis=1)

    marker_base = {c: c.split(":")[0].strip() for c in dev.index}
    forced = [c for c in dev.index if marker_base[c] in KEY_MARKERS]
    ranked = spread.sort_values(ascending=False).index.tolist()
    selected = list(dict.fromkeys(forced + ranked))[:N_MARKERS_IN_HEATMAP]
    selected = spread.loc[selected].sort_values(ascending=False).index.tolist()

    dev_sel = dev.loc[selected]
    vmax = float(np.nanmax(np.abs(dev_sel.to_numpy())))
    vmax = max(vmax, 0.5)

    fig, ax = plt.subplots(figsize=(20, max(15, 0.65 * len(selected))))
    cmap = LinearSegmentedColormap.from_list(
        "drift", ["#2166AC", "#F7F7F7", "#B2182B"])
    im = ax.imshow(dev_sel.to_numpy(), aspect="auto", cmap=cmap,
                   norm=TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax))
    ax.set_xticks(xpos)
    ax.set_xticklabels(x_labels3, fontsize=FONT_SIZE_TICK - 12)
    color_position_ticks(ax)
    ax.set_yticks(np.arange(len(selected)))
    ylabels = [f"{c}{'  *' if marker_base[c] in KEY_MARKERS else ''}"
               for c in selected]
    ax.set_yticklabels(ylabels, fontsize=FONT_SIZE_TICK - 12)
    draw_scan_boundaries(ax)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("log$_2$(p99 / marker median)", fontsize=FONT_SIZE_BASE - 6)
    cb.ax.tick_params(labelsize=FONT_SIZE_TICK - 10)
    ax.set_title("Staining intensity drift across sections\n"
                 "(* = marker used in phenotype calling; black rule = scan boundary)",
                 fontsize=FONT_SIZE_TITLE - 4)
    save_fig(fig, "F08_marker_p99_heatmap")

    # ---- split the spread into between-scan and within-scan parts -----------
    # Descriptive only. Scan and arm are the same factor in this cohort, so this
    # separates section-level drift from scan-level drift and attributes neither
    # of them to batch. Script 02 handles the modelling.
    scan_cols = {}
    for s in SAMPLE_ORDER:
        scan_cols.setdefault(SCAN_OF[s], []).append(s)
    scan_med = pd.DataFrame({sc: p99[cols].median(axis=1)
                             for sc, cols in scan_cols.items()})
    with np.errstate(divide="ignore", invalid="ignore"):
        between_scan = np.log2(scan_med.max(axis=1)
                               / scan_med.min(axis=1).replace(0, np.nan))
        within_parts = {}
        for sc, cols in scan_cols.items():
            block = p99[cols]
            within_parts[sc] = np.log2(block.max(axis=1)
                                       / block.min(axis=1).replace(0, np.nan))
        within_scan = pd.DataFrame(within_parts).max(axis=1)

    drift_out = pd.DataFrame({
        "column": spread.index,
        "marker": [marker_base[c] for c in spread.index],
        "is_key_marker": [marker_base[c] in KEY_MARKERS for c in spread.index],
        "log2_spread_across_sections": spread.values,
        "log2_between_scan_medians": between_scan.reindex(spread.index).values,
        "log2_max_within_scan_spread": within_scan.reindex(spread.index).values,
    })
    drift_out["between_scan_share"] = (
        drift_out["log2_between_scan_medians"]
        / (drift_out["log2_between_scan_medians"]
           + drift_out["log2_max_within_scan_spread"]).replace(0, np.nan))
    drift_out = drift_out.sort_values("log2_spread_across_sections",
                                      ascending=False)
    out_path = os.path.join(FIG_DIR, "F08_marker_drift_ranking.csv")
    drift_out.to_csv(out_path, index=False, lineterminator="\n")
    print(f"    wrote {out_path}")

    print("\n    Top 12 markers by total spread, split into its parts")
    print(f"      {'column':<42}{'total':>9}{'between':>10}{'within':>9}"
          f"{'btw share':>11}")
    for _, r in drift_out.head(12).iterrows():
        share = r["between_scan_share"]
        share_txt = "na" if pd.isna(share) else f"{share:.2f}"
        print(f"      {r['column'][:41]:<42}"
              f"{r['log2_spread_across_sections']:>9.2f}"
              f"{r['log2_between_scan_medians']:>10.2f}"
              f"{r['log2_max_within_scan_spread']:>9.2f}{share_txt:>11}")
    print("\n    A high between-scan share means the marker separates the two")
    print("    acquisitions more than it separates sections within one. Because")
    print("    scan and arm are the same factor here, that is consistent with")
    print("    batch AND with a real arm difference, and this table cannot tell")
    print("    them apart. It is here to set expectations for script 02.")


# %% Cell 10 - F09 nearest-neighbour feasibility
# =============================================================================

banner("F09 - NEAREST-NEIGHBOUR FEASIBILITY")

fig, ax = plt.subplots(figsize=(20, 14))
yy = np.arange(len(PHENO_USE))
for j, s in enumerate(SAMPLE_ORDER):
    vals = counts.loc[PHENO_USE, s].to_numpy(dtype=float)
    cond = cond_labels[j]
    offset = (j - (len(SAMPLE_ORDER) - 1) / 2) * 0.085
    plotted = np.where(vals < 1, 0.7, vals)
    is_pos1 = POS_OF[s] == FLAG_POSITION_RANK
    ax.scatter(plotted, yy + offset, s=320,
               color=CONDITION_COLORS.get(cond, "#999999"),
               edgecolor=POSITION_FLAG_COLOR if is_pos1 else "#FFFFFF",
               linewidth=3 if is_pos1 else 2, zorder=3)
    zero_idx = np.where(vals == 0)[0]
    if len(zero_idx):
        ax.scatter(np.full(len(zero_idx), 0.7), zero_idx + offset, s=420,
                   marker="x", color=ZERO_FLAG_COLOR, linewidth=3.5, zorder=4)

ax.axvline(MIN_CELLS_FOR_NN, color="#000000", linestyle="--", linewidth=3.5,
           zorder=2)
ax.text(MIN_CELLS_FOR_NN * 1.15, len(PHENO_USE) - 0.4,
        f"{MIN_CELLS_FOR_NN} cells", fontsize=FONT_SIZE_ANNOT, rotation=90,
        va="top")
ax.set_xscale("log")
ax.set_yticks(yy)
ax.set_yticklabels(PHENO_USE, fontsize=FONT_SIZE_TICK - 4)
ax.invert_yaxis()
ax.set_xlabel("Cells per section (log scale)")
ax.set_title("Cells available per phenotype in each animal\n"
             "(red x = absent, red outline = position-1 section; "
             "each point is one animal)",
             fontsize=FONT_SIZE_TITLE - 6)
ax.xaxis.grid(True, color=GRID_COLOR, linewidth=1.5)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
add_condition_legend(ax, loc="lower right")
save_fig(fig, "F09_nn_feasibility")


# %% Cell 11 - F10 within-lineage ratios (optional diagnostic)
# =============================================================================

banner("F10 - WITHIN-LINEAGE RATIOS")

ratio_tbl = pd.DataFrame()
if not RUN_LINEAGE_RATIOS:
    print("    SKIPPED (RUN_LINEAGE_RATIOS is False)")
else:
    print("    Each ratio is internal to one lineage, so section composition")
    print("    cancels the same way it does in F06. A ratio that tracks the scan")
    print("    boundary rather than the animal is a warning about the marker")
    print("    that defines the numerator, not a biological result.\n")

    usable = [(t, num, den, ylab, why) for (t, num, den, ylab, why) in LINEAGE_RATIOS
              if all(p in counts.index for p in den)]
    skipped = [t for (t, num, den, _, _) in LINEAGE_RATIOS
               if not all(p in counts.index for p in den)]
    for t in skipped:
        print(f"    WARNING: '{t}' skipped, denominator phenotypes missing")

    if usable:
        rows = []
        fig, axes = plt.subplots(1, len(usable), figsize=(14 * len(usable), 13),
                                 squeeze=False)
        for k, (title, num_ph, den_ph, ylab, why) in enumerate(usable):
            num_ph_use = [p for p in num_ph if p in counts.index]
            num = counts.loc[num_ph_use, SAMPLE_ORDER].sum(axis=0).to_numpy(float)
            den = counts.loc[den_ph, SAMPLE_ORDER].sum(axis=0).to_numpy(float)
            with np.errstate(divide="ignore", invalid="ignore"):
                frac = np.where(den > 0, 100.0 * num / den, np.nan)

            for i, s in enumerate(SAMPLE_ORDER):
                rows.append({
                    "ratio": title, "sample_id": s,
                    "animal_id": short_label(s),
                    "treatment": condition_of(s, meta_map),
                    "scan_id": SCAN_OF[s], "slide_position_rank": POS_OF[s],
                    "numerator": int(num[i]), "denominator": int(den[i]),
                    "pct": frac[i],
                })

            ax = axes[0, k]
            ax.bar(xpos, np.nan_to_num(frac), color=cond_colors,
                   edgecolor="#FFFFFF", linewidth=2, zorder=3)
            for i in range(len(SAMPLE_ORDER)):
                if np.isnan(frac[i]):
                    continue
                ax.text(i, frac[i] + 1.5,
                        f"{frac[i]:.1f}%\n{int(num[i]):,}/{int(den[i]):,}",
                        ha="center", va="bottom", fontsize=FONT_SIZE_ANNOT - 10)
            ax.set_xticks(xpos)
            ax.set_xticklabels(x_labels3, fontsize=FONT_SIZE_TICK - 12)
            color_position_ticks(ax)
            ax.set_ylabel(ylab, fontsize=FONT_SIZE_BASE - 6)
            ax.set_ylim(0, max(105, float(np.nanmax(frac)) * 1.25))
            ax.set_title(f"{title}\n{why}", fontsize=FONT_SIZE_TITLE - 14)
            style_axes(ax); add_condition_bands(ax); draw_scan_boundaries(ax)

        ratio_tbl = pd.DataFrame(rows)
        out_path = os.path.join(FIG_DIR, "F10_lineage_ratios.csv")
        ratio_tbl.to_csv(out_path, index=False, lineterminator="\n")
        fig.suptitle("Within-lineage ratios: composition cancels, marker calling "
                     "does not\n(black rule = scan boundary, red label = "
                     "position-1 section)",
                     y=1.04, fontsize=FONT_SIZE_TITLE - 4)
        save_fig(fig, "F10_lineage_internal_ratios")
        print(f"    wrote {out_path}")

        for title, g in ratio_tbl.groupby("ratio", sort=False):
            print(f"\n    {title}")
            for _, r in g.iterrows():
                pos = r["slide_position_rank"]
                pos_txt = f"pos{int(pos)}" if pd.notna(pos) else "pos na"
                print(f"      {r['animal_id']:<8}{r['treatment']:<11}"
                      f"{r['scan_id']:<10}{pos_txt:<8}"
                      f"{r['pct']:>7.1f}%   "
                      f"{r['numerator']:>7,} / {r['denominator']:>8,}")


# %% Cell 12 - deck summary table and wrap up
# =============================================================================

banner("DECK SUMMARY TABLE")

deck_cols = ["sample_id", "animal_id", "treatment", "n_cells",
             "bbox_area_mm2", "cells_per_mm2_bbox", "n_unique_phenotypes"]
for c in (SCAN_COL, POSITION_COL):
    if c in file_summary.columns:
        deck_cols.insert(3, c)
deck = file_summary[deck_cols].copy()
deck = deck.set_index("sample_id").loc[SAMPLE_ORDER].reset_index()

for ph in ["CD68+IDO1+ Macrophages", "CD68+IDO1- Macrophages",
           "Helper T cells", "CD4- T cells", "B cells", "Neutrophils"]:
    if ph in counts.index:
        deck[f"n_{ascii_safe(ph)}"] = counts.loc[ph, SAMPLE_ORDER].values.astype(int)

if len(ido_tbl) > 0:
    deck = deck.merge(ido_tbl[["sample_id", "pct_ido1_pos_of_cd68"]],
                      on="sample_id", how="left")

if len(ratio_tbl) > 0:
    for title, g in ratio_tbl.groupby("ratio", sort=False):
        key = ascii_safe(title).lower().replace(" ", "_")
        deck = deck.merge(
            g[["sample_id", "pct"]].rename(columns={"pct": f"pct_{key}"}),
            on="sample_id", how="left")

deck_path = os.path.join(FIG_DIR, "F00_deck_summary_table.csv")
deck.to_csv(deck_path, index=False, lineterminator="\n")
print(f"    wrote {deck_path}")
with pd.option_context("display.width", 250, "display.max_columns", 60):
    print("\n" + deck.to_string(index=False))

banner("DONE")
print(f"Figures written to: {FIG_DIR}")
print("\nSuggested slide order for the deck:")
slides = [
    ("F01_cells_per_animal", "Cells recovered per section"),
    ("F02_tissue_area_and_density", "Section footprint and density"),
    ("F03_composition_stacked", "Phenotype composition"),
    ("F04_composition_heatmap", "Cell counts, phenotype x animal"),
    ("F05_phenotype_by_condition", "Phenotype frequency by condition"),
    ("F06_ido1_fraction_of_cd68", "IDO1+ fraction of CD68+ macrophages"),
    ("F07_immune_vs_structural", "Immune vs structural composition"),
    ("F08_marker_p99_heatmap", "Staining intensity drift"),
    ("F09_nn_feasibility", "Cells available for spatial statistics"),
]
if RUN_LINEAGE_RATIOS and len(ratio_tbl) > 0:
    slides.append(("F10_lineage_internal_ratios", "Within-lineage marker ratios"))
for i, (stem, cap) in enumerate(slides, start=1):
    print(f"   {i}. {stem:<34} {cap}")

print("\nRead F10 before script 02. If the CD4 fraction of T cells tracks the")
print("scan boundary more than it tracks the animal, then helper T cells are")
print("partly a staining readout on this panel, and any downstream step that")
print("gates on helper T density inherits that. BALT detection is one such step.")
