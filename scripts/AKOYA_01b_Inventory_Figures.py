#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler inventory - FIGURE SUPPLEMENT
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 01b of the AKOYA analysis series.

PURPOSE
    Turn the tables written by AKOYA_01_Inventory.py into presentation-ready
    figures for review with Deepak's group. Reads ONLY the inventory tables,
    never the raw QuPath CSVs, so the figures cannot disagree with the numbers
    already generated. Runs in seconds.

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

INPUTS (from AKOYA_01_Inventory.py)
    tables/01_file_summary.csv
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

# Markers highlighted in the drift heatmap because they define phenotype calls
KEY_MARKERS = [
    "IDO1", "CD68", "CD163", "CD206", "CD4", "CD8", "CD3e", "CD20", "CD79a",
    "MPO", "CD11b", "FoxP3", "iNOS", "Arginase-1", "IFNG", "Granzyme-B",
    "HLA-DR", "CD21", "PD-1", "PD-L1",
]
N_MARKERS_IN_HEATMAP = 40   # top N by cross-section spread, key markers forced in

# Threshold line drawn on the feasibility plot
MIN_CELLS_FOR_NN = 50

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


# %% Cell 3 - load tables and build ordering
# =============================================================================

os.makedirs(FIG_DIR, exist_ok=True)

banner("AKOYA INVENTORY FIGURES")
print(f"Run time  : {datetime.now().isoformat(timespec='seconds')}")
print(f"Tables    : {TABLE_DIR}")
print(f"Figures   : {FIG_DIR}\n")

file_summary = load_table("01_file_summary.csv", required=True)
pheno_long = load_table("04_phenotype_counts_long.csv", required=True)
marker_stats = load_table("06_marker_stats_long.csv", required=False)
label_tbl = load_table("03_phenotype_labels_repr.csv", required=False)

# ---- animal ordering: condition first, then animal ID ------------------------
file_summary["treatment"] = file_summary["treatment"].astype(str)
cond_rank = {c: i for i, c in enumerate(CONDITION_ORDER)}
file_summary["_rank"] = file_summary["treatment"].map(
    lambda t: cond_rank.get(t, len(cond_rank))
)
file_summary = file_summary.sort_values(["_rank", "animal_id"]).reset_index(drop=True)
SAMPLE_ORDER = file_summary["sample_id"].tolist()
meta_map = file_summary.set_index("sample_id").to_dict(orient="index")

print(f"\nSection order used in all figures:")
for sid in SAMPLE_ORDER:
    m = meta_map[sid]
    print(f"    {sid:<12} {m['treatment']:<10} {m['n_cells']:>9,} cells")

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
x_labels = [short_label(s) for s in SAMPLE_ORDER]
xpos = np.arange(len(SAMPLE_ORDER))


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

fig, ax = plt.subplots(figsize=(18, 11))
bars = ax.bar(xpos, [meta_map[s]["n_cells"] for s in SAMPLE_ORDER],
              color=cond_colors, edgecolor="#FFFFFF", linewidth=2, zorder=3)
for i, s in enumerate(SAMPLE_ORDER):
    ax.text(i, meta_map[s]["n_cells"] * 1.02, f"{meta_map[s]['n_cells']:,}",
            ha="center", va="bottom", fontsize=FONT_SIZE_ANNOT, rotation=0)
ax.set_xticks(xpos)
ax.set_xticklabels(x_labels)
ax.set_ylabel("Segmented cells")
ax.set_xlabel("Animal")
ax.set_title("Cells recovered per lung section")
ax.set_ylim(0, max(meta_map[s]["n_cells"] for s in SAMPLE_ORDER) * 1.15)
ax.yaxis.set_major_formatter(
    matplotlib.ticker.FuncFormatter(lambda v, _: f"{v/1000:,.0f}k"))
style_axes(ax)
add_condition_bands(ax)
add_condition_legend(ax)
save_fig(fig, "F01_cells_per_animal")

# ---- F02: two panels, tissue footprint and density --------------------------
fig, axes = plt.subplots(1, 2, figsize=(26, 11))

ax = axes[0]
vals = [meta_map[s].get("bbox_area_mm2", np.nan) for s in SAMPLE_ORDER]
ax.bar(xpos, vals, color=cond_colors, edgecolor="#FFFFFF", linewidth=2, zorder=3)
ax.set_xticks(xpos); ax.set_xticklabels(x_labels, rotation=45, ha="right")
ax.set_ylabel("Bounding box area (mm$^2$)")
ax.set_title("Section footprint")
style_axes(ax); add_condition_bands(ax)

ax = axes[1]
vals = [meta_map[s].get("cells_per_mm2_bbox", np.nan) for s in SAMPLE_ORDER]
ax.bar(xpos, vals, color=cond_colors, edgecolor="#FFFFFF", linewidth=2, zorder=3)
ax.set_xticks(xpos); ax.set_xticklabels(x_labels, rotation=45, ha="right")
ax.set_ylabel("Cells / mm$^2$")
ax.set_title("Cell density (bounding box)")
style_axes(ax); add_condition_bands(ax); add_condition_legend(ax)

fig.suptitle("Tissue footprint and cell density per section", y=1.02)
save_fig(fig, "F02_tissue_area_and_density")


# %% Cell 5 - F03 stacked composition, F04 count heatmap
# =============================================================================

banner("F03 / F04 - PHENOTYPE COMPOSITION")

fig, ax = plt.subplots(figsize=(20, 13))
bottom = np.zeros(len(SAMPLE_ORDER))
for ph in PHENO_USE:
    vals = pcts.loc[ph, SAMPLE_ORDER].to_numpy(dtype=float)
    ax.bar(xpos, vals, bottom=bottom, color=PHENOTYPE_COLORS.get(ph, "#999999"),
           edgecolor="#FFFFFF", linewidth=1.5, label=ph, zorder=3)
    bottom += vals
ax.set_xticks(xpos)
ax.set_xticklabels([f"{x}\n{c}" for x, c in zip(x_labels, cond_labels)],
                   fontsize=FONT_SIZE_TICK - 4)
ax.set_ylabel("Percent of section")
ax.set_ylim(0, 100)
ax.set_title("Phenotype composition per section")
style_axes(ax, ygrid=False)
ax.legend(bbox_to_anchor=(1.01, 1.0), loc="upper left", frameon=False,
          fontsize=FONT_SIZE_LEGEND - 6)
save_fig(fig, "F03_composition_stacked")

# ---- F04: log10 count heatmap ----------------------------------------------
fig, ax = plt.subplots(figsize=(20, 14))
mat = np.log10(counts.loc[PHENO_USE, SAMPLE_ORDER].to_numpy(dtype=float) + 1.0)
cmap = LinearSegmentedColormap.from_list(
    "akoya_counts", ["#FFFFFF", "#C6DBEF", "#4292C6", "#08519C", "#08306B"])
im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=np.nanmax(mat))
ax.set_xticks(xpos)
ax.set_xticklabels([f"{x}\n{c}" for x, c in zip(x_labels, cond_labels)],
                   fontsize=FONT_SIZE_TICK - 4)
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
# condition divider
n_first = sum(1 for c in cond_labels if c == CONDITION_ORDER[0])
if 0 < n_first < len(SAMPLE_ORDER):
    ax.axvline(n_first - 0.5, color="#000000", linewidth=4)
ax.set_title("Cell counts per phenotype and section\n(red zeros = phenotype absent)")
save_fig(fig, "F04_composition_heatmap")


# %% Cell 6 - F05 per-phenotype percent by condition
# =============================================================================

banner("F05 - PHENOTYPE FREQUENCY BY CONDITION")

n_ph = len(PHENO_USE)
ncol = 4
nrow = int(np.ceil(n_ph / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(9 * ncol, 7.5 * nrow))
axes = np.atleast_1d(axes).ravel()

rng = np.random.default_rng(0)
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
        for v, j in zip(vals, jitter):
            if v == 0:
                ax.scatter([ci + j], [0], s=520, marker="x",
                           color=ZERO_FLAG_COLOR, linewidth=4, zorder=5)
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
             "bar = group median)", y=1.005, fontsize=FONT_SIZE_TITLE)
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
        "animal_id": x_labels,
        "treatment": cond_labels,
        "n_ido1_pos": num.astype(int),
        "n_cd68_total": den.astype(int),
        "pct_ido1_pos_of_cd68": frac,
    })
    print(ido_tbl.to_string(index=False))

    fig, axes = plt.subplots(1, 2, figsize=(28, 12),
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
    ax.set_xticklabels([f"{x}\n{c}" for x, c in zip(x_labels, cond_labels)],
                       fontsize=FONT_SIZE_TICK - 4)
    ax.set_ylabel("IDO1+ % of CD68+ macrophages")
    ax.set_ylim(0, 110)
    ax.set_title("Per animal")
    style_axes(ax); add_condition_bands(ax)

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
                 "(internal normalization, removes section composition effects)",
                 y=1.04)
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

fig, ax = plt.subplots(figsize=(20, 12))
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
ax.set_xticklabels([f"{x}\n{c}" for x, c in zip(x_labels, cond_labels)],
                   fontsize=FONT_SIZE_TICK - 4)
ax.set_ylabel("Percent of section")
ax.set_ylim(0, 100)
ax.set_title("Immune, structural and unassigned fractions")
style_axes(ax, ygrid=False)
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

    sub = dev.loc[selected]
    vmax = float(np.nanmax(np.abs(sub.to_numpy())))
    vmax = max(vmax, 0.5)

    fig, ax = plt.subplots(figsize=(20, max(14, 0.65 * len(selected))))
    cmap = LinearSegmentedColormap.from_list(
        "drift", ["#2166AC", "#F7F7F7", "#B2182B"])
    im = ax.imshow(sub.to_numpy(), aspect="auto", cmap=cmap,
                   norm=TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax))
    ax.set_xticks(xpos)
    ax.set_xticklabels([f"{x}\n{c}" for x, c in zip(x_labels, cond_labels)],
                       fontsize=FONT_SIZE_TICK - 6)
    ax.set_yticks(np.arange(len(selected)))
    ylabels = [f"{c}{'  *' if marker_base[c] in KEY_MARKERS else ''}"
               for c in selected]
    ax.set_yticklabels(ylabels, fontsize=FONT_SIZE_TICK - 12)
    if 0 < n_first < len(SAMPLE_ORDER):
        ax.axvline(n_first - 0.5, color="#000000", linewidth=4)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("log$_2$(p99 / marker median)", fontsize=FONT_SIZE_BASE - 6)
    cb.ax.tick_params(labelsize=FONT_SIZE_TICK - 10)
    ax.set_title("Staining intensity drift across sections\n"
                 "(* = marker used in phenotype calling; black line = slide boundary)",
                 fontsize=FONT_SIZE_TITLE - 4)
    save_fig(fig, "F08_marker_p99_heatmap")

    drift_out = pd.DataFrame({
        "column": spread.index,
        "log2_spread_across_sections": spread.values,
    }).sort_values("log2_spread_across_sections", ascending=False)
    drift_out.to_csv(os.path.join(FIG_DIR, "F08_marker_drift_ranking.csv"),
                     index=False, lineterminator="\n")
    print(f"    wrote {os.path.join(FIG_DIR, 'F08_marker_drift_ranking.csv')}")


# %% Cell 10 - F09 nearest-neighbour feasibility
# =============================================================================

banner("F09 - NEAREST-NEIGHBOUR FEASIBILITY")

fig, ax = plt.subplots(figsize=(20, 13))
yy = np.arange(len(PHENO_USE))
for j, s in enumerate(SAMPLE_ORDER):
    vals = counts.loc[PHENO_USE, s].to_numpy(dtype=float)
    cond = cond_labels[j]
    offset = (j - (len(SAMPLE_ORDER) - 1) / 2) * 0.085
    plotted = np.where(vals < 1, 0.7, vals)
    ax.scatter(plotted, yy + offset, s=320,
               color=CONDITION_COLORS.get(cond, "#999999"),
               edgecolor="#FFFFFF", linewidth=2, zorder=3)
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
             "(red x = absent; each point is one animal)",
             fontsize=FONT_SIZE_TITLE - 2)
ax.xaxis.grid(True, color=GRID_COLOR, linewidth=1.5)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
add_condition_legend(ax, loc="lower right")
save_fig(fig, "F09_nn_feasibility")


# %% Cell 11 - deck summary table and wrap up
# =============================================================================

banner("DECK SUMMARY TABLE")

deck = file_summary[["sample_id", "animal_id", "treatment", "n_cells",
                     "bbox_area_mm2", "cells_per_mm2_bbox",
                     "n_unique_phenotypes"]].copy()
deck = deck.set_index("sample_id").loc[SAMPLE_ORDER].reset_index()

for ph in ["CD68+IDO1+ Macrophages", "CD68+IDO1- Macrophages",
           "Helper T cells", "CD4- T cells", "B cells", "Neutrophils"]:
    if ph in counts.index:
        deck[f"n_{ascii_safe(ph)}"] = counts.loc[ph, SAMPLE_ORDER].values.astype(int)

if len(ido_tbl) > 0:
    deck = deck.merge(ido_tbl[["sample_id", "pct_ido1_pos_of_cd68"]],
                      on="sample_id", how="left")

deck_path = os.path.join(FIG_DIR, "F00_deck_summary_table.csv")
deck.to_csv(deck_path, index=False, lineterminator="\n")
print(f"    wrote {deck_path}")
with pd.option_context("display.width", 250, "display.max_columns", 50):
    print("\n" + deck.to_string(index=False))

banner("DONE")
print(f"Figures written to: {FIG_DIR}")
print("\nSuggested slide order for the deck:")
for i, (stem, cap) in enumerate([
    ("F01_cells_per_animal", "Cells recovered per section"),
    ("F02_tissue_area_and_density", "Section footprint and density"),
    ("F03_composition_stacked", "Phenotype composition"),
    ("F04_composition_heatmap", "Cell counts, phenotype x animal"),
    ("F05_phenotype_by_condition", "Phenotype frequency by condition"),
    ("F06_ido1_fraction_of_cd68", "IDO1+ fraction of CD68+ macrophages"),
    ("F07_immune_vs_structural", "Immune vs structural composition"),
    ("F08_marker_p99_heatmap", "Staining intensity drift"),
    ("F09_nn_feasibility", "Cells available for spatial statistics"),
], start=1):
    print(f"   {i}. {stem:<34} {cap}")
