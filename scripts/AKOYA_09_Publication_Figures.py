#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - PUBLICATION FIGURES FOR THE THREE CARRIED-FORWARD FINDINGS
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 09 of the AKOYA analysis series.

This script does no new inference. It assembles the three findings that survived
the full diagnostic chain into figures, and writes a headline numbers table for
the manuscript.

    FINDING 1 - DISEASE BURDEN
        Zero granuloma-density tissue in every treated animal against 7.3, 16.0
        and 39.1 percent in untreated, under an identical absolute definition
        (3,000 myeloid cells/mm^2 applied to every section). Before any size
        filter, only 0.0 to 0.1 percent of treated tissue exceeds threshold.
        Simple, large and unarguable.

    FINDING 2 - IDO1+ MACROPHAGE FRACTION
        1.0 to 11.7 percent of CD68-lineage macrophages in treated animals
        against 49.4 to 93.7 percent untreated. Completely separated. Robust to
        cutoff across raw IDO1 6.1 to 9.7, and the vendor phenotype call is
        validated as an IDO1 intensity gate at AUC 0.977 to 1.000.

    FINDING 3 - B LINEAGE RADIAL POSITION
        B cells and plasma cells sit closer to the myeloid focus core in treated
        animals. Coefficient -0.137 radial units, p = 0.010, direction unanimous
        across all six leave-one-out fits, and unchanged when helper T cells are
        removed. Not circular: lymphocytes never enter the focus definition.

DROPPED, AND WHY (recorded here so it is not relitigated)
    3-hydroxykynurenine   no usable dynamic range (96-99% of cells nonzero,
                          p99 0.53 to 6.85). Background, not signal.
    raw / normalised
    distance outcomes     raw is bounded by focus size in principle; normalising
                          by radius over-corrects (r = -0.12 to -0.68 with
                          radius, worse than raw at -0.10 to +0.39). Delta is
                          retained but nothing survives BH.
    iNOS x Arginase-1
    co-expression         the CD3e x CD20 control pair separates the arms at 4
                          of 6 thresholds, MORE than the test pair at 3 of 6. A
                          macrophage cannot be both a T and a B cell, so that is
                          segmentation spillover. Control behaviour is too
                          inconsistent to calibrate against.

INPUTS
    structures_rev4/tables/33b_burden_summary.csv
    structures_rev4/tables/35_foci_structures_relative.csv
    structures_rev4/cell_assignments/<section>_cell_structures.csv
    lymphocyte_radial/tables/71_centered_radial_per_animal.csv
    lymphocyte_radial/tables/72_models_all.csv
    lymphocyte_radial/tables/73_radial_profiles_core.csv

OUTPUTS
    figures/  F60_finding1_burden
              F61_finding2_ido1
              F62_finding3_b_lineage
              F63_summary_three_findings
    tables/   80_headline_numbers.csv
              81_ido1_fraction.csv
              82_ido1_threshold_sensitivity.csv

USAGE
    conda activate sc_pre
    python AKOYA_09_Publication_Figures.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================

STRUCT_DIR = "/master/jlehle/WORKING/AKOYA/structures_rev4"
CELL_DIR = f"{STRUCT_DIR}/cell_assignments"
BURDEN_TABLE = f"{STRUCT_DIR}/tables/33b_burden_summary.csv"
FOCI_TABLE = f"{STRUCT_DIR}/tables/35_foci_structures_relative.csv"

LYM_DIR = "/master/jlehle/WORKING/AKOYA/lymphocyte_radial/tables"
CENTERED_TABLE = f"{LYM_DIR}/71_centered_radial_per_animal.csv"
MODELS_TABLE = f"{LYM_DIR}/72_models_all.csv"
PROFILE_TABLE = f"{LYM_DIR}/73_radial_profiles_core.csv"

OUT_DIR = "/master/jlehle/WORKING/AKOYA/publication_figures"

CONDITION_ORDER = ["D1MT", "Untreated"]
CONDITION_COLORS = {"D1MT": "#2C7FB8", "Untreated": "#D95F02"}
CONDITION_LABELS = {"D1MT": "D1MT-treated", "Untreated": "Untreated"}

# ---- phenotypes -------------------------------------------------------------
IDO1_POS = "CD68+IDO1+ Macrophages"
IDO1_NEG = "CD68+IDO1- Macrophages"
CD68_LINEAGE = [IDO1_POS, IDO1_NEG]
MYELOID = [IDO1_POS, IDO1_NEG, "CD163+ Macrophages", "Neutrophils"]
B_LINEAGE = ["B cells", "Plasma cells"]
T_LINEAGE = ["Helper T cells", "CD4- T cells", "Tregs"]
LYMPHOCYTES = T_LINEAGE + B_LINEAGE

PHENOTYPE_COLORS = {
    IDO1_POS: "#B2182B", IDO1_NEG: "#EF8A62", "CD163+ Macrophages": "#FDBE85",
    "Neutrophils": "#7B3294", "Helper T cells": "#1B7837",
    "CD4- T cells": "#7FBC41", "Tregs": "#00441B", "B cells": "#2166AC",
    "Plasma cells": "#67A9CF",
}

# ---- IDO1 threshold work ----------------------------------------------------
IDO1_COL = "IDO1"
N_CUTOFF_STEPS = 200
PLATEAU_TOL = 0.05          # window where every section moves less than this

# ---- representative sections for the spatial panels ------------------------
REPRESENTATIVE = {"D1MT": "G3_43118", "Untreated": "G4_43109"}
MAX_POINTS_MAP = 60000

# ---- plotting ---------------------------------------------------------------
FONT_SIZE_BASE = 28
FONT_SIZE_TITLE = 34
FONT_SIZE_TICK = 28
FONT_SIZE_LEGEND = 28
FONT_SIZE_ANNOT = 26
FONT_SIZE_PANEL = 40        # panel letters
DPI = 300
SAVE_PDF = True
SAVE_PNG = True
POINT_MARKERS = ["o", "s", "^", "D", "v", "P"]
RANDOM_SEED = 0

GRID_COLOR = "#DDDDDD"
AXIS_COLOR = "#333333"
TEXT_COLOR = "#000000"
FLAG_COLOR = "#B2182B"
BACKGROUND_COLOR = "#ECECEC"


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
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

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


def panel_letter(ax, letter, dx=-0.14, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes,
            fontsize=FONT_SIZE_PANEL, fontweight="bold", va="top", ha="left")


def mann_whitney_auc(pos, neg):
    pos = np.asarray(pos, float); neg = np.asarray(neg, float)
    pos = pos[np.isfinite(pos)]; neg = neg[np.isfinite(neg)]
    n1, n2 = len(pos), len(neg)
    if n1 == 0 or n2 == 0:
        return np.nan
    allv = np.concatenate([pos, neg])
    order = np.argsort(allv, kind="mergesort")
    ranks = np.empty(len(allv), float)
    ranks[order] = np.arange(1, len(allv) + 1)
    sv = allv[order]
    i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = np.mean(ranks[order[i:j + 1]])
        i = j + 1
    return float((ranks[:n1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n2))


def load_optional(path, label):
    if not os.path.exists(path):
        print(f"    WARNING: {label} not found at {path}")
        return None
    df = pd.read_csv(path)
    print(f"    loaded {label}  ({df.shape[0]} x {df.shape[1]})")
    return df


_tee = Tee(os.path.join(TAB_DIR, "00_publication_figures_report.txt"))
sys.stdout = _tee

banner("AKOYA PUBLICATION FIGURES - THREE FINDINGS")
print(f"Run time : {datetime.now().isoformat(timespec='seconds')}")
print(f"Output   : {OUT_DIR}")
print("\nNo new inference. Assembly only.")


# %% Cell 3 - load
# =============================================================================

banner("LOADING")

burden = load_optional(BURDEN_TABLE, "burden summary")
foci = load_optional(FOCI_TABLE, "foci table")
centered = load_optional(CENTERED_TABLE, "centered radial per animal")
models = load_optional(MODELS_TABLE, "models")
profiles = load_optional(PROFILE_TABLE, "radial profiles")

paths = sorted(glob.glob(os.path.join(CELL_DIR, "*_cell_structures.csv")))
if not paths:
    print(f"ERROR: no cell assignment files in {CELL_DIR}.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

cells = {}
for p in paths:
    sid = os.path.basename(p).replace("_cell_structures.csv", "")
    try:
        d = pd.read_csv(p, low_memory=False)
    except Exception as e:
        print(f"    ERROR reading {sid}: {e}")
        continue
    d["sample_id"] = sid
    cells[sid] = d
    print(f"    {sid:<12} {d['condition'].iloc[0]:<10} {len(d):>9,} cells")

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

xpos = np.arange(len(SAMPLE_ORDER))
xlab = [short_label(s) for s in SAMPLE_ORDER]
xcol = [CONDITION_COLORS[COND_OF[s]] for s in SAMPLE_ORDER]


def cond_legend(ax, **kw):
    ax.legend(handles=[Patch(facecolor=CONDITION_COLORS[c],
                             label=CONDITION_LABELS[c])
                       for c in CONDITION_ORDER],
              frameon=False, fontsize=kw.pop("fontsize", FONT_SIZE_LEGEND - 8),
              **kw)


# %% Cell 4 - IDO1 quantities, recomputed so the figure is self-contained
# =============================================================================

banner("RECOMPUTING IDO1 QUANTITIES")

ido_rows, ido_vals = [], {}
for s in SAMPLE_ORDER:
    d = cells[s]
    mac = d.loc[d["pheno"].isin(CD68_LINEAGE)]
    if not len(mac):
        continue
    v = pd.to_numeric(mac[IDO1_COL], errors="coerce").dropna().to_numpy()
    ido_vals[s] = v
    pos = pd.to_numeric(mac.loc[mac["pheno"] == IDO1_POS, IDO1_COL],
                        errors="coerce").dropna().to_numpy()
    neg = pd.to_numeric(mac.loc[mac["pheno"] == IDO1_NEG, IDO1_COL],
                        errors="coerce").dropna().to_numpy()
    ido_rows.append({
        "sample_id": s, "animal_id": short_label(s), "condition": COND_OF[s],
        "n_cd68_lineage": len(mac), "n_ido1_pos": len(pos),
        "pct_ido1_pos": 100.0 * len(pos) / len(mac),
        "auc_ido1_call": mann_whitney_auc(pos, neg),
        "median_ido1_pos": float(np.median(pos)) if len(pos) else np.nan,
        "median_ido1_neg": float(np.median(neg)) if len(neg) else np.nan,
    })
ido = pd.DataFrame(ido_rows)
write_csv(ido, "81_ido1_fraction.csv")
sub("IDO1+ fraction of CD68-lineage macrophages, and call validation")
print(ido.to_string(index=False))

# ---- threshold sensitivity --------------------------------------------------
allv = np.concatenate([v for v in ido_vals.values() if len(v)])
hi = float(np.log1p(np.percentile(allv, 99.9)))
grid = np.linspace(0.0, hi, N_CUTOFF_STEPS)
sens_rows = []
for s in SAMPLE_ORDER:
    lv = np.log1p(ido_vals.get(s, np.array([])))
    if not len(lv):
        continue
    for g in grid:
        sens_rows.append({"sample_id": s, "condition": COND_OF[s],
                          "log1p_cutoff": float(g),
                          "pct_above": 100.0 * float(np.mean(lv > g))})
sens = pd.DataFrame(sens_rows)
write_csv(sens, "82_ido1_threshold_sensitivity.csv")

sw = sens.pivot_table(index="log1p_cutoff", columns="sample_id",
                      values="pct_above")[SAMPLE_ORDER] / 100.0
best_lo = best_hi = np.nan
i = 0
while i < len(sw):
    j = i
    while j + 1 < len(sw):
        block = sw.iloc[i:j + 2]
        if (block.max() - block.min()).max() > PLATEAU_TOL:
            break
        j += 1
    span = sw.index[j] - sw.index[i]
    if np.isnan(best_lo) or span > (best_hi - best_lo):
        best_lo, best_hi = float(sw.index[i]), float(sw.index[j])
    i = max(j, i + 1)
print(f"\n    cutoff-insensitive window: log1p {best_lo:.2f} to {best_hi:.2f}, "
      f"raw IDO1 {np.expm1(best_lo):.2f} to {np.expm1(best_hi):.2f}")


# %% Cell 5 - FIGURE 1: disease burden
# =============================================================================

banner("FIGURE 1 - DISEASE BURDEN")

fig = plt.figure(figsize=(40, 26))
gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.25], hspace=0.32, wspace=0.28)

# --- A: burden percentage ---
ax = fig.add_subplot(gs[0, 0])
if burden is not None:
    b = burden.set_index("sample_id").reindex(SAMPLE_ORDER)
    ax.bar(xpos, b["burden_pct"], color=xcol, edgecolor="#FFFFFF", linewidth=2,
           zorder=3)
    for i, v in enumerate(b["burden_pct"]):
        ax.text(i, v + 1.0, f"{v:.1f}%", ha="center",
                fontsize=FONT_SIZE_ANNOT - 4,
                fontweight="bold" if v == 0 else "normal",
                color=FLAG_COLOR if v == 0 else TEXT_COLOR)
ax.set_xticks(xpos); ax.set_xticklabels(xlab, rotation=45, ha="right",
                                        fontsize=FONT_SIZE_TICK - 6)
ax.set_ylabel("Lung area in granuloma-\ndensity regions (%)",
              fontsize=FONT_SIZE_BASE - 6)
ax.set_title("Lesion burden", fontsize=FONT_SIZE_TITLE - 8)
style_axes(ax); panel_letter(ax, "A"); cond_legend(ax, loc="upper left")

# --- B: lesion count ---
ax = fig.add_subplot(gs[0, 1])
if burden is not None:
    ax.bar(xpos, b["n_regions"], color=xcol, edgecolor="#FFFFFF", linewidth=2,
           zorder=3)
    for i, v in enumerate(b["n_regions"]):
        ax.text(i, v + 0.4, f"{int(v)}", ha="center",
                fontsize=FONT_SIZE_ANNOT - 4,
                fontweight="bold" if v == 0 else "normal",
                color=FLAG_COLOR if v == 0 else TEXT_COLOR)
ax.set_xticks(xpos); ax.set_xticklabels(xlab, rotation=45, ha="right",
                                        fontsize=FONT_SIZE_TICK - 6)
ax.set_ylabel("Granuloma-density regions", fontsize=FONT_SIZE_BASE - 6)
ax.set_title("Lesion count", fontsize=FONT_SIZE_TITLE - 8)
style_axes(ax); panel_letter(ax, "B")

# --- C: background to peak, with the threshold drawn ---
ax = fig.add_subplot(gs[0, 2])
if burden is not None:
    for i, s in enumerate(SAMPLE_ORDER):
        r = b.loc[s]
        ax.plot([i, i], [r["background_density"], r["peak_density"]],
                color=COLOR_OF[s], linewidth=6, alpha=0.6, zorder=2)
        ax.scatter([i], [r["background_density"]], s=380, color="#999999",
                   marker="_", linewidth=6, zorder=3)
        ax.scatter([i], [r["peak_density"]], s=420, color=COLOR_OF[s],
                   marker=MARKER_OF[s], edgecolor="#FFFFFF", linewidth=2,
                   zorder=4)
ax.axhline(3000, color="#000000", linestyle="--", linewidth=3.5)
ax.text(len(SAMPLE_ORDER) - 0.4, 3300, "granuloma-density\nthreshold",
        ha="right", fontsize=FONT_SIZE_ANNOT - 8)
ax.set_yscale("log")
ax.set_xticks(xpos); ax.set_xticklabels(xlab, rotation=45, ha="right",
                                        fontsize=FONT_SIZE_TICK - 6)
ax.set_ylabel("Myeloid density (cells / mm$^2$)", fontsize=FONT_SIZE_BASE - 6)
ax.set_title("Background to peak range", fontsize=FONT_SIZE_TITLE - 10)
style_axes(ax); panel_letter(ax, "C")

# --- D and E: representative maps ---
for k, cond in enumerate(CONDITION_ORDER):
    sid = REPRESENTATIVE.get(cond)
    ax = fig.add_subplot(gs[1, k])
    if sid not in cells:
        ax.axis("off"); continue
    d = cells[sid]
    idx = (rng.choice(len(d), size=MAX_POINTS_MAP, replace=False)
           if len(d) > MAX_POINTS_MAP else np.arange(len(d)))
    dd = d.iloc[idx]
    other = ~dd["pheno"].isin(MYELOID)
    ax.scatter(dd.loc[other, "x"], dd.loc[other, "y"], s=0.6,
               color=BACKGROUND_COLOR, linewidths=0, rasterized=True)
    for ph in MYELOID:
        m = dd["pheno"] == ph
        if m.sum():
            ax.scatter(dd.loc[m, "x"], dd.loc[m, "y"], s=4,
                       color=PHENOTYPE_COLORS[ph], linewidths=0,
                       rasterized=True)
    if "in_burden_region" in dd.columns:
        m = dd["in_burden_region"].astype(bool)
        if m.sum():
            ax.scatter(dd.loc[m, "x"], dd.loc[m, "y"], s=2,
                       facecolor="none", edgecolor="#000000", linewidth=0.3,
                       rasterized=True)
    bp = (float(b.loc[sid, "burden_pct"]) if burden is not None else np.nan)
    ax.set_title(f"{CONDITION_LABELS[cond]}  ({short_label(sid)})\n"
                 f"burden {bp:.1f}%", fontsize=FONT_SIZE_TITLE - 12)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_edgecolor(CONDITION_COLORS[cond]); sp.set_linewidth(6)
    panel_letter(ax, "D" if k == 0 else "E", dx=-0.02, dy=1.02)

ax = fig.add_subplot(gs[1, 2]); ax.axis("off")
ax.text(0.0, 0.95,
        "Granuloma-density regions are defined by\n"
        "a single absolute threshold (3,000 myeloid\n"
        "cells/mm$^2$) applied identically to every\n"
        "section, then filtered at 30,000 µm$^2$ and\n"
        "100 cells.\n\n"
        "Before any size filter, 0.0 to 0.1% of\n"
        "treated tissue exceeds the threshold at all.\n\n"
        "Detection uses pooled myeloid density only\n"
        "and never uses IDO1.",
        transform=ax.transAxes, fontsize=FONT_SIZE_ANNOT - 4, va="top")
handles = [Patch(facecolor=PHENOTYPE_COLORS[p], label=p) for p in MYELOID]
handles.append(Patch(facecolor=BACKGROUND_COLOR, label="all other cells"))
ax.legend(handles=handles, loc="lower left", frameon=False,
          fontsize=FONT_SIZE_LEGEND - 12)

fig.suptitle("Finding 1: D1MT treatment eliminates granuloma-density tissue",
             y=0.985, fontsize=FONT_SIZE_TITLE)
save_fig(fig, "F60_finding1_burden")


# %% Cell 6 - FIGURE 2: IDO1+ macrophage fraction
# =============================================================================

banner("FIGURE 2 - IDO1+ MACROPHAGE FRACTION")

fig, axes = plt.subplots(2, 2, figsize=(34, 26))

# --- A: fraction per animal ---
ax = axes[0, 0]
ax.bar(xpos, ido.set_index("sample_id").reindex(SAMPLE_ORDER)["pct_ido1_pos"],
       color=xcol, edgecolor="#FFFFFF", linewidth=2, zorder=3)
ii = ido.set_index("sample_id").reindex(SAMPLE_ORDER)
for i, s in enumerate(SAMPLE_ORDER):
    r = ii.loc[s]
    ax.text(i, r["pct_ido1_pos"] + 2.5,
            f"{r['pct_ido1_pos']:.1f}%\n{int(r['n_ido1_pos']):,}/"
            f"{int(r['n_cd68_lineage']):,}",
            ha="center", fontsize=FONT_SIZE_ANNOT - 10)
ax.set_xticks(xpos); ax.set_xticklabels(xlab, rotation=45, ha="right",
                                        fontsize=FONT_SIZE_TICK - 6)
ax.set_ylabel("IDO1+ % of CD68-lineage\nmacrophages", fontsize=FONT_SIZE_BASE - 6)
ax.set_ylim(0, 115)
ax.set_title("IDO1+ macrophage fraction", fontsize=FONT_SIZE_TITLE - 8)
style_axes(ax); panel_letter(ax, "A"); cond_legend(ax, loc="upper left")

# --- B: by arm, showing complete separation ---
ax = axes[0, 1]
for ci, c in enumerate(CONDITION_ORDER):
    v = ido.loc[ido["condition"] == c, "pct_ido1_pos"].to_numpy()
    j = rng.uniform(-0.09, 0.09, size=len(v))
    ax.scatter(np.full(len(v), ci) + j, v, s=620, color=CONDITION_COLORS[c],
               edgecolor="#FFFFFF", linewidth=3, zorder=3)
    ax.hlines(np.median(v), ci - 0.28, ci + 0.28, color="#000000", linewidth=5)
d1 = ido.loc[ido["condition"] == "D1MT", "pct_ido1_pos"]
un = ido.loc[ido["condition"] == "Untreated", "pct_ido1_pos"]
if len(d1) and len(un):
    gap_lo, gap_hi = float(d1.max()), float(un.min())
    ax.axhspan(gap_lo, gap_hi, color="#000000", alpha=0.06, zorder=0)
    ax.text(0.5, (gap_lo + gap_hi) / 2, "complete\nseparation", ha="center",
            va="center", fontsize=FONT_SIZE_ANNOT - 4)
ax.set_xticks(range(len(CONDITION_ORDER)))
ax.set_xticklabels([CONDITION_LABELS[c] for c in CONDITION_ORDER],
                   fontsize=FONT_SIZE_TICK - 4)
ax.set_xlim(-0.6, len(CONDITION_ORDER) - 0.4); ax.set_ylim(0, 105)
ax.set_ylabel("IDO1+ % of CD68-lineage", fontsize=FONT_SIZE_BASE - 6)
ax.set_title("By treatment arm", fontsize=FONT_SIZE_TITLE - 8)
style_axes(ax); panel_letter(ax, "B")

# --- C: threshold sensitivity ---
ax = axes[1, 0]
for s in SAMPLE_ORDER:
    d = sens.loc[sens["sample_id"] == s].sort_values("log1p_cutoff")
    ax.plot(d["log1p_cutoff"], d["pct_above"], linewidth=5, color=COLOR_OF[s],
            alpha=0.9)
    step = max(1, len(d) // 14)
    ax.plot(d["log1p_cutoff"].to_numpy()[::step],
            d["pct_above"].to_numpy()[::step], linestyle="none",
            marker=MARKER_OF[s], markersize=15, color=COLOR_OF[s])
if np.isfinite(best_lo):
    ax.axvspan(best_lo, best_hi, color="#000000", alpha=0.09, zorder=0)
    ax.text((best_lo + best_hi) / 2, 103,
            f"cutoff-insensitive\nraw IDO1 {np.expm1(best_lo):.1f}"
            f" to {np.expm1(best_hi):.1f}",
            ha="center", fontsize=FONT_SIZE_ANNOT - 10)
ax.set_xlabel("log(1 + IDO1) cutoff", fontsize=FONT_SIZE_BASE - 6)
ax.set_ylabel("% of CD68+ cells called positive", fontsize=FONT_SIZE_BASE - 8)
ax.set_ylim(0, 112)
ax.set_title("The result does not depend on the cutoff",
             fontsize=FONT_SIZE_TITLE - 12)
style_axes(ax); panel_letter(ax, "C")
ax.legend(handles=[Line2D([0], [0], color=COLOR_OF[s], marker=MARKER_OF[s],
                          markersize=14, linewidth=4,
                          label=f"{short_label(s)}")
                   for s in SAMPLE_ORDER],
          frameon=False, fontsize=FONT_SIZE_LEGEND - 16, loc="upper right",
          ncol=2)

# --- D: AUC validation ---
ax = axes[1, 1]
ax.bar(xpos, ii["auc_ido1_call"], color=xcol, edgecolor="#FFFFFF", linewidth=2,
       zorder=3)
for i, s in enumerate(SAMPLE_ORDER):
    ax.text(i, ii.loc[s, "auc_ido1_call"] + 0.005,
            f"{ii.loc[s, 'auc_ido1_call']:.3f}", ha="center",
            fontsize=FONT_SIZE_ANNOT - 8)
ax.axhline(0.95, color="#000000", linestyle="--", linewidth=3)
ax.text(len(SAMPLE_ORDER) - 0.4, 0.952, "0.95", ha="right",
        fontsize=FONT_SIZE_ANNOT - 8)
ax.set_ylim(0.9, 1.02)
ax.set_xticks(xpos); ax.set_xticklabels(xlab, rotation=45, ha="right",
                                        fontsize=FONT_SIZE_TICK - 6)
ax.set_ylabel("AUC, IDO1 intensity separating\nthe phenotype call",
              fontsize=FONT_SIZE_BASE - 10)
ax.set_title("The phenotype call is an IDO1 gate",
             fontsize=FONT_SIZE_TITLE - 12)
style_axes(ax); panel_letter(ax, "D")

fig.suptitle("Finding 2: IDO1+ macrophages are near-absent in treated animals",
             y=0.98, fontsize=FONT_SIZE_TITLE)
fig.tight_layout(rect=[0, 0, 1, 0.96])
save_fig(fig, "F61_finding2_ido1")


# %% Cell 7 - FIGURE 3: B lineage radial position
# =============================================================================

banner("FIGURE 3 - B LINEAGE RADIAL POSITION")

fig, axes = plt.subplots(2, 2, figsize=(34, 26))

# --- A: radial profile ---
ax = axes[0, 0]
if profiles is not None and "B lineage" in profiles.columns:
    for s in SAMPLE_ORDER:
        d = profiles.loc[profiles["sample_id"] == s].sort_values("radial_centre")
        if not len(d):
            continue
        ax.plot(d["radial_centre"], d["B lineage"], linewidth=5,
                marker=MARKER_OF[s], markersize=16, color=COLOR_OF[s], alpha=0.9)
ax.set_xlabel("radial position (0 = core centre, 1 = boundary)",
              fontsize=FONT_SIZE_BASE - 8)
ax.set_ylabel("B lineage cells (% of bin)", fontsize=FONT_SIZE_BASE - 8)
ax.set_title("B lineage abundance across the core",
             fontsize=FONT_SIZE_TITLE - 10)
ax.set_ylim(bottom=0)
style_axes(ax); panel_letter(ax, "A")
ax.legend(handles=[Line2D([0], [0], color=COLOR_OF[s], marker=MARKER_OF[s],
                          markersize=14, linewidth=4,
                          label=f"{short_label(s)} ({CONDITION_LABELS[COND_OF[s]]})")
                   for s in SAMPLE_ORDER],
          frameon=False, fontsize=FONT_SIZE_LEGEND - 16, loc="best")

# --- B: per-animal centered position ---
ax = axes[0, 1]
show = ["B cells", "Plasma cells", "Helper T cells", "CD4- T cells", "Tregs"]
yy = np.arange(len(show))
if centered is not None:
    for i, p in enumerate(show):
        for s in SAMPLE_ORDER:
            v = centered.loc[(centered["sample_id"] == s)
                             & (centered["phenotype"] == p)]
            if not len(v):
                continue
            off = (SAMPLE_ORDER.index(s) - (len(SAMPLE_ORDER) - 1) / 2) * 0.13
            thin = int(v["n_cells"].iloc[0]) < 50
            ax.scatter(v["mean_centered_radial"].iloc[0], i + off, s=420,
                       color=COLOR_OF[s], marker=MARKER_OF[s],
                       edgecolor=FLAG_COLOR if thin else "#FFFFFF",
                       linewidth=3 if thin else 2, zorder=3)
ax.axvline(0, color="#000000", linewidth=3.5)
ax.axhspan(-0.5, 1.5, color="#2166AC", alpha=0.06, zorder=0)
ax.set_yticks(yy); ax.set_yticklabels(show, fontsize=FONT_SIZE_TICK - 6)
ax.invert_yaxis()
ax.set_xlabel("mean centred radial position\n(negative = core-ward)",
              fontsize=FONT_SIZE_BASE - 8)
ax.set_title("Per animal, core cells only", fontsize=FONT_SIZE_TITLE - 10)
ax.xaxis.grid(True, color=GRID_COLOR, linewidth=1.5); ax.set_axisbelow(True)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
panel_letter(ax, "B")

# --- C: primary models ---
ax = axes[1, 0]
if models is not None:
    keep = models.loc[models["family"].isin(["primary", "lineage",
                                             "per_phenotype"])].copy()
    keep = keep.loc[~keep.get("circular", pd.Series(False, index=keep.index))
                    .fillna(False).astype(bool)]
    keep = keep.sort_values("coef_D1MT_vs_ref")
    yy2 = np.arange(len(keep))
    for i, (_, r) in enumerate(keep.iterrows()):
        sig = r["p_value"] < 0.05
        col = FLAG_COLOR if sig else "#777777"
        ax.plot([r["ci_low"], r["ci_high"]], [i, i], color=col, linewidth=6,
                zorder=2)
        ax.scatter([r["coef_D1MT_vs_ref"]], [i], s=420, color=col,
                   edgecolor="#FFFFFF", linewidth=2, zorder=3)
        ax.text(r["ci_high"], i, f"  p={r['p_value']:.3f}", va="center",
                fontsize=FONT_SIZE_ANNOT - 12)
    ax.set_yticks(yy2)
    ax.set_yticklabels([a.replace("PRIMARY: ", "").replace(" (core)", "")[:34]
                        for a in keep["analysis"]],
                       fontsize=FONT_SIZE_TICK - 12)
    ax.invert_yaxis()
ax.axvline(0, color="#000000", linewidth=3.5)
ax.set_xlabel("D1MT effect on centred radial position (95% CI)",
              fontsize=FONT_SIZE_BASE - 10)
ax.set_title("Mixed models, non-circular populations",
             fontsize=FONT_SIZE_TITLE - 12)
ax.xaxis.grid(True, color=GRID_COLOR, linewidth=1.5); ax.set_axisbelow(True)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
panel_letter(ax, "C")

# --- D: sensitivity ---
ax = axes[1, 1]
if models is not None:
    prim = models.loc[models["family"] == "primary"]
    g = models.loc[models["family"].isin(["leave_one_out",
                                          "marker_robustness"])].copy()
    yy3 = np.arange(len(g))
    for i, (_, r) in enumerate(g.iterrows()):
        col = (CONDITION_COLORS.get(r.get("dropped_arm", ""), "#777777")
               if r["family"] == "leave_one_out" else "#444444")
        ax.plot([r["ci_low"], r["ci_high"]], [i, i], color=col, linewidth=6,
                zorder=2)
        ax.scatter([r["coef_D1MT_vs_ref"]], [i], s=380, color=col,
                   edgecolor="#FFFFFF", linewidth=2, zorder=3)
    if len(prim):
        ax.axvline(float(prim["coef_D1MT_vs_ref"].iloc[0]), color="#000000",
                   linestyle="--", linewidth=3)
    ax.set_yticks(yy3)
    ax.set_yticklabels([a[:28] for a in g["analysis"]],
                       fontsize=FONT_SIZE_TICK - 14)
    ax.invert_yaxis()
ax.axvline(0, color="#000000", linewidth=3.5)
ax.set_xlabel("D1MT effect (95% CI)", fontsize=FONT_SIZE_BASE - 10)
ax.set_title("Sensitivity: every fit stays negative\n(dashed = full model)",
             fontsize=FONT_SIZE_TITLE - 14)
ax.xaxis.grid(True, color=GRID_COLOR, linewidth=1.5); ax.set_axisbelow(True)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
panel_letter(ax, "D")

fig.suptitle("Finding 3: B lineage cells sit closer to the myeloid core in "
             "treated foci", y=0.98, fontsize=FONT_SIZE_TITLE)
fig.tight_layout(rect=[0, 0, 1, 0.96])
save_fig(fig, "F62_finding3_b_lineage")


# %% Cell 8 - FIGURE 4: the three findings side by side
# =============================================================================

banner("FIGURE 4 - SUMMARY")

fig, axes = plt.subplots(1, 3, figsize=(42, 14))

ax = axes[0]
if burden is not None:
    for ci, c in enumerate(CONDITION_ORDER):
        v = burden.loc[burden["condition"] == c, "burden_pct"].to_numpy()
        j = rng.uniform(-0.09, 0.09, size=len(v))
        ax.scatter(np.full(len(v), ci) + j, v, s=620, color=CONDITION_COLORS[c],
                   edgecolor="#FFFFFF", linewidth=3, zorder=3)
        ax.hlines(np.median(v), ci - 0.28, ci + 0.28, color="#000000",
                  linewidth=5)
ax.set_xticks(range(len(CONDITION_ORDER)))
ax.set_xticklabels([CONDITION_LABELS[c] for c in CONDITION_ORDER],
                   fontsize=FONT_SIZE_TICK - 6)
ax.set_xlim(-0.6, len(CONDITION_ORDER) - 0.4)
ax.set_ylabel("Lesion burden (% lung area)", fontsize=FONT_SIZE_BASE - 6)
ax.set_title("1. Disease burden", fontsize=FONT_SIZE_TITLE - 6)
style_axes(ax); panel_letter(ax, "A")

ax = axes[1]
for ci, c in enumerate(CONDITION_ORDER):
    v = ido.loc[ido["condition"] == c, "pct_ido1_pos"].to_numpy()
    j = rng.uniform(-0.09, 0.09, size=len(v))
    ax.scatter(np.full(len(v), ci) + j, v, s=620, color=CONDITION_COLORS[c],
               edgecolor="#FFFFFF", linewidth=3, zorder=3)
    ax.hlines(np.median(v), ci - 0.28, ci + 0.28, color="#000000", linewidth=5)
ax.set_xticks(range(len(CONDITION_ORDER)))
ax.set_xticklabels([CONDITION_LABELS[c] for c in CONDITION_ORDER],
                   fontsize=FONT_SIZE_TICK - 6)
ax.set_xlim(-0.6, len(CONDITION_ORDER) - 0.4); ax.set_ylim(0, 105)
ax.set_ylabel("IDO1+ % of CD68-lineage", fontsize=FONT_SIZE_BASE - 6)
ax.set_title("2. IDO1+ macrophages", fontsize=FONT_SIZE_TITLE - 6)
style_axes(ax); panel_letter(ax, "B")

ax = axes[2]
if centered is not None:
    bl = centered.loc[centered["phenotype"].isin(B_LINEAGE)]
    agg = (bl.groupby(["sample_id", "condition"])
           .apply(lambda g: np.average(g["mean_centered_radial"],
                                       weights=g["n_cells"]))
           .reset_index(name="b_lineage_centered"))
    for ci, c in enumerate(CONDITION_ORDER):
        v = agg.loc[agg["condition"] == c, "b_lineage_centered"].to_numpy()
        j = rng.uniform(-0.09, 0.09, size=len(v))
        ax.scatter(np.full(len(v), ci) + j, v, s=620, color=CONDITION_COLORS[c],
                   edgecolor="#FFFFFF", linewidth=3, zorder=3)
        ax.hlines(np.median(v), ci - 0.28, ci + 0.28, color="#000000",
                  linewidth=5)
    write_csv(agg, "83_b_lineage_centered_per_animal.csv")
ax.axhline(0, color="#000000", linestyle="--", linewidth=3)
ax.set_xticks(range(len(CONDITION_ORDER)))
ax.set_xticklabels([CONDITION_LABELS[c] for c in CONDITION_ORDER],
                   fontsize=FONT_SIZE_TICK - 6)
ax.set_xlim(-0.6, len(CONDITION_ORDER) - 0.4)
ax.set_ylabel("B lineage centred radial position\n(negative = core-ward)",
              fontsize=FONT_SIZE_BASE - 10)
ax.set_title("3. B lineage position", fontsize=FONT_SIZE_TITLE - 6)
style_axes(ax); panel_letter(ax, "C")

fig.suptitle("Three findings, each animal one point", y=1.02,
             fontsize=FONT_SIZE_TITLE)
fig.tight_layout(rect=[0, 0, 1, 0.95])
save_fig(fig, "F63_summary_three_findings")


# %% Cell 9 - headline numbers table
# =============================================================================

banner("HEADLINE NUMBERS")

rows = []
if burden is not None:
    for c in CONDITION_ORDER:
        v = burden.loc[burden["condition"] == c, "burden_pct"]
        n = burden.loc[burden["condition"] == c, "n_regions"]
        rows.append({"finding": "1. Lesion burden", "arm": c,
                     "metric": "% lung area in granuloma-density regions",
                     "values": ", ".join(f"{x:.1f}" for x in v),
                     "median": float(v.median()),
                     "n_animals": len(v),
                     "extra": f"regions: {', '.join(str(int(x)) for x in n)}"})
for c in CONDITION_ORDER:
    v = ido.loc[ido["condition"] == c, "pct_ido1_pos"]
    a = ido.loc[ido["condition"] == c, "auc_ido1_call"]
    rows.append({"finding": "2. IDO1+ macrophage fraction", "arm": c,
                 "metric": "% of CD68-lineage macrophages",
                 "values": ", ".join(f"{x:.1f}" for x in v),
                 "median": float(v.median()), "n_animals": len(v),
                 "extra": f"AUC of call: {a.min():.3f} to {a.max():.3f}"})
if centered is not None:
    bl = centered.loc[centered["phenotype"].isin(B_LINEAGE)]
    agg = (bl.groupby(["sample_id", "condition"])
           .apply(lambda g: np.average(g["mean_centered_radial"],
                                       weights=g["n_cells"]))
           .reset_index(name="v"))
    for c in CONDITION_ORDER:
        v = agg.loc[agg["condition"] == c, "v"]
        rows.append({"finding": "3. B lineage radial position", "arm": c,
                     "metric": "mean centred radial position (core)",
                     "values": ", ".join(f"{x:+.3f}" for x in v),
                     "median": float(v.median()), "n_animals": len(v),
                     "extra": ""})
if models is not None:
    for fam, name in [("primary", "pooled lymphocytes"),
                      ("lineage", "lineage split")]:
        for _, r in models.loc[models["family"] == fam].iterrows():
            rows.append({"finding": "3. Model", "arm": "D1MT vs Untreated",
                         "metric": r["analysis"],
                         "values": f"coef {r['coef_D1MT_vs_ref']:+.4f} "
                                   f"[{r['ci_low']:+.4f}, {r['ci_high']:+.4f}]",
                         "median": float(r["coef_D1MT_vs_ref"]),
                         "n_animals": int(r["n_animals"]),
                         "extra": f"p = {r['p_value']:.4f}, "
                                  f"{int(r['n_cells']):,} cells"})
head = pd.DataFrame(rows)
write_csv(head, "80_headline_numbers.csv")
with pd.option_context("display.width", 250, "display.max_colwidth", 60):
    print(head.to_string(index=False))

banner("DONE")
print("Figures:")
print("  F60_finding1_burden            disease burden")
print("  F61_finding2_ido1              IDO1+ macrophage fraction")
print("  F62_finding3_b_lineage         B lineage radial position")
print("  F63_summary_three_findings     all three, one point per animal")
print("\nEvery figure carries only findings that survived the diagnostic chain.")
print("3-HK, raw and normalised distance, and iNOS/Arginase-1 co-expression are")
print("deliberately absent. See the header of this script for why.")

sys.stdout = _tee.terminal
_tee.close()
