#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - WITHIN-FOCUS SPATIAL ARCHITECTURE
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 05 of the AKOYA analysis series.

WHAT THIS ANSWERS
    Q1  Where do CD4+ T cells sit relative to IDO1+ CD68+ macrophages, treated
        versus untreated? Asked as radial position within a focus and as
        within-focus nearest-neighbour distance, so structure size cannot drive
        the answer.
    Q4  Are FoxP3+ Tregs enriched in the neighbourhood of IDO1+ macrophages?
        (the canonical IDO-Treg axis)
    Q6a Do iNOS-high and Arginase-1-high macrophages occupy distinct
        neighbourhoods within a focus?

    Plus a short per-focus IDO1 vs 3-hydroxykynurenine preview, which costs
    nothing because those medians are already in table 35. The proper 3-HK
    spatial gradient is script 06.

WHY THIS DESIGN
    IDO1+ macrophages occur only inside granulomas, and untreated lesions are
    several fold larger than treated residual foci. A raw cell-to-cell distance
    comparison between arms would therefore report lesion diameter while
    appearing to report cell-cell proximity. Two things prevent that here:

    1. RADIAL POSITION is normalised per focus (0 = core centre, 1 = core
       boundary, 1-2 = cuff), so a 130 um treated focus and a 480 um untreated
       lesion both span the same range.
    2. Every statistic is compared against a WITHIN-FOCUS PERMUTATION NULL:
       phenotype labels are shuffled among the cells of that focus while
       coordinates are held fixed. The null therefore already contains that
       focus's size, shape, and cell density, and the effect size is what
       remains after all of it is accounted for.

    The unit of observation is the cell, summarised to the focus, summarised to
    the animal, then compared across arms. With three animals per arm the
    animal-level test is capped near p = 0.10, so results are reported as
    per-animal effect sizes with consistency (k of 3), not as a single p-value.

    NOTHING IS FILTERED OUT UP FRONT. Every focus and every pair is computed
    regardless of cell count, and each result carries an evidence tier
    (solid / usable / provisional / not interpretable) so stringency can be
    applied after seeing the results rather than before.

INPUTS
    structures_rev3/cell_assignments/<section>_cell_structures.csv
    structures_rev3/tables/35_foci_structures_relative.csv

OUTPUTS
    figures/  F38 .. F45
    tables/   40 .. 47

RUNTIME NOTE
    The permutation work dominates. Defaults (200 permutations, anchors capped
    at 1000 per focus, k = 100 neighbours) run in roughly 10-25 minutes on one
    node. Raise N_PERMUTATIONS once the result looks stable.

USAGE
    conda activate sc_pre
    python AKOYA_05_Radial_Architecture.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================

IN_DIR = "/master/jlehle/WORKING/AKOYA/structures_rev4"
CELL_DIR = f"{IN_DIR}/cell_assignments"
FOCI_TABLE = f"{IN_DIR}/tables/35_foci_structures_relative.csv"
OUT_DIR = "/master/jlehle/WORKING/AKOYA/radial"

CONDITION_ORDER = ["D1MT", "Untreated"]
CONDITION_COLORS = {"D1MT": "#2C7FB8", "Untreated": "#D95F02"}

# "foci"   : relative myeloid foci, identical definition in both arms (PRIMARY)
# "burden" : absolute density regions, anatomically sized but untreated-only
STRUCTURE_SOURCE = "foci"

# ---- phenotypes -------------------------------------------------------------
IDO1_POS = "CD68+IDO1+ Macrophages"
IDO1_NEG = "CD68+IDO1- Macrophages"
CD68_LINEAGE = [IDO1_POS, IDO1_NEG]
MACROPHAGE_ALL = [IDO1_POS, IDO1_NEG, "CD163+ Macrophages"]

PHENOTYPE_ORDER = [
    IDO1_POS, IDO1_NEG, "CD163+ Macrophages", "Neutrophils",
    "Helper T cells", "CD4- T cells", "Tregs", "B cells", "Plasma cells",
    "Endothelial cells", "Epithelial/Tumor cells", "Other",
]
PHENOTYPE_COLORS = {
    IDO1_POS: "#B2182B", IDO1_NEG: "#EF8A62", "CD163+ Macrophages": "#FDBE85",
    "Neutrophils": "#7B3294", "Helper T cells": "#1B7837",
    "CD4- T cells": "#7FBC41", "Tregs": "#00441B", "B cells": "#2166AC",
    "Plasma cells": "#67A9CF", "Endothelial cells": "#E7298A",
    "Epithelial/Tumor cells": "#66C2A5", "Other": "#A6761D",
}

# ---- radial binning ---------------------------------------------------------
RADIAL_MIN, RADIAL_MAX = 0.0, 2.0
N_RADIAL_BINS = 20                  # 10 across core, 10 across cuff
MIN_CELLS_PER_RADIAL_BIN = 20       # bins thinner than this are reported as NaN

# ---- nearest-neighbour pairs ------------------------------------------------
NN_ANCHORS = [IDO1_POS, IDO1_NEG]
NN_TARGETS = ["Helper T cells", "CD4- T cells", "Tregs", "B cells"]

# ---- permutation ------------------------------------------------------------
N_PERMUTATIONS = 200
MAX_ANCHORS_PER_FOCUS = 1000        # anchors subsampled above this, for speed
KNN_DEPTH = 100                     # neighbours pre-indexed per cell
RANDOM_SEED = 0

# ---- evidence tiers (REPORTING ONLY - nothing is excluded) ------------------
TIER_SOLID = 200
TIER_USABLE = 50
TIER_PROVISIONAL = 20

MIN_CELLS_PER_FOCUS = 100           # foci below this are flagged, not dropped

# ---- iNOS / Arginase-1 segregation -----------------------------------------
POLARISATION_PERCENTILES = [60, 70, 75, 80, 90]   # swept, within section
POLARISATION_PRIMARY = 75
MIXING_K = 10                       # neighbours used for the mixing score

# ---- plotting ---------------------------------------------------------------
FONT_SIZE_BASE = 28
FONT_SIZE_TITLE = 34
FONT_SIZE_TICK = 28
FONT_SIZE_LEGEND = 28
FONT_SIZE_ANNOT = 26
DPI = 300
SAVE_PDF = True
SAVE_PNG = True
POINT_MARKERS = ["o", "s", "^", "D", "v", "P"]

GRID_COLOR = "#DDDDDD"
AXIS_COLOR = "#333333"
TEXT_COLOR = "#000000"
FLAG_COLOR = "#B2182B"
CORE_BAND = "#F0F0F0"


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

try:
    from scipy.spatial import cKDTree
    HAVE_SCIPY = True
except Exception as _e:
    HAVE_SCIPY = False
    _scipy_err = _e

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


def tier_of(n):
    if n >= TIER_SOLID:
        return "solid"
    if n >= TIER_USABLE:
        return "usable"
    if n >= TIER_PROVISIONAL:
        return "provisional"
    return "not_interpretable"


def empirical_p(obs, null):
    """Two-sided empirical p with the +1 correction, and a z against the null."""
    null = np.asarray(null, float)
    null = null[np.isfinite(null)]
    if not len(null) or not np.isfinite(obs):
        return np.nan, np.nan
    mu, sd = float(np.mean(null)), float(np.std(null, ddof=1))
    z = (obs - mu) / sd if sd > 0 else np.nan
    n_ext = int(np.sum(np.abs(null - mu) >= abs(obs - mu)))
    p = (n_ext + 1) / (len(null) + 1)
    return float(min(p, 1.0)), z


_tee = Tee(os.path.join(TAB_DIR, "00_radial_report.txt"))
sys.stdout = _tee

banner("AKOYA WITHIN-FOCUS SPATIAL ARCHITECTURE")
print(f"Run time         : {datetime.now().isoformat(timespec='seconds')}")
print(f"Structure source : {STRUCTURE_SOURCE}")
print(f"Permutations     : {N_PERMUTATIONS}")
print(f"Output           : {OUT_DIR}")
if not HAVE_SCIPY:
    print(f"\nERROR: scipy required ({_scipy_err}).")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

print("\nNOTE ON UNITS: radial position is normalised per focus, and every")
print("permutation null is built within the focus with coordinates held fixed.")
print("Structure size and local cell density are therefore already inside the")
print("null and cannot drive any effect reported below.")


# %% Cell 3 - load
# =============================================================================

banner("LOADING")

paths = sorted(glob.glob(os.path.join(CELL_DIR, "*_cell_structures.csv")))
if not paths:
    print(f"ERROR: no cell assignment files in {CELL_DIR}. Run script 04 first.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

STRUCT_COL = "focus_id" if STRUCTURE_SOURCE == "foci" else "burden_region_id"

cells = {}
for p in paths:
    sid = os.path.basename(p).replace("_cell_structures.csv", "")
    try:
        d = pd.read_csv(p, low_memory=False)
    except Exception as e:
        print(f"    ERROR reading {sid}: {e}. Skipping.")
        continue
    if STRUCT_COL not in d.columns:
        print(f"    ERROR: {sid} lacks '{STRUCT_COL}'. Skipping.")
        continue
    cells[sid] = d
    n_in = int((d[STRUCT_COL] > 0).sum())
    print(f"    {sid:<12} {d['condition'].iloc[0]:<10} {len(d):>9,} cells, "
          f"{n_in:>8,} in structures")

if not cells:
    print("ERROR: nothing loaded.")
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


def structure_frames():
    """Yield (sid, struct_id, sub_dataframe) for every structure with cells."""
    for s in SAMPLE_ORDER:
        d = cells[s]
        for k, g in d.loc[d[STRUCT_COL] > 0].groupby(STRUCT_COL):
            yield s, int(k), g


inventory = []
for s, k, g in structure_frames():
    inventory.append({
        "sample_id": s, "condition": COND_OF[s], "structure_id": k,
        "n_cells": len(g),
        "n_core": int((g["region"] == "core").sum()),
        "n_cuff": int((g["region"] == "cuff").sum()),
        "n_ido1_pos": int((g["pheno"] == IDO1_POS).sum()),
        "n_cd68_lineage": int(g["pheno"].isin(CD68_LINEAGE).sum()),
        "below_min_cells": bool(len(g) < MIN_CELLS_PER_FOCUS),
    })
inv = pd.DataFrame(inventory)
sub("Structure inventory (nothing excluded, small structures flagged)")
for c in CONDITION_ORDER:
    i = inv.loc[inv["condition"] == c]
    if not len(i):
        continue
    print(f"    {c:<12} {len(i):>3} structures, "
          f"{int(i['n_cells'].sum()):>8,} cells, "
          f"{int(i['below_min_cells'].sum())} below {MIN_CELLS_PER_FOCUS} cells")
write_csv(inv, "40_structure_inventory.csv")


# %% Cell 4 - A: radial composition profiles
# =============================================================================

banner("A - RADIAL COMPOSITION PROFILES")

edges = np.linspace(RADIAL_MIN, RADIAL_MAX, N_RADIAL_BINS + 1)
centres = 0.5 * (edges[:-1] + edges[1:])

prof_rows = []
for s, k, g in structure_frames():
    r = pd.to_numeric(g["radial_pos"], errors="coerce").to_numpy()
    ok = np.isfinite(r)
    if not ok.sum():
        continue
    b = np.clip(np.digitize(r[ok], edges[1:-1]), 0, N_RADIAL_BINS - 1)
    ph = g["pheno"].to_numpy()[ok]
    for bi in range(N_RADIAL_BINS):
        m = b == bi
        n = int(m.sum())
        if n == 0:
            continue
        rec = {"sample_id": s, "condition": COND_OF[s], "structure_id": k,
               "bin": bi, "radial_centre": float(centres[bi]), "n_cells": n,
               "sufficient": n >= MIN_CELLS_PER_RADIAL_BIN}
        for p in PHENOTYPE_ORDER:
            rec[p] = 100.0 * float((ph[m] == p).mean())
        prof_rows.append(rec)

prof = pd.DataFrame(prof_rows)
write_csv(prof, "41_radial_profiles_per_structure.csv")

# animal-level: mean across structures, weighted by cells in that bin
anim_rows = []
for s in SAMPLE_ORDER:
    d = prof.loc[(prof["sample_id"] == s) & prof["sufficient"]]
    for bi in range(N_RADIAL_BINS):
        b = d.loc[d["bin"] == bi]
        if not len(b):
            continue
        w = b["n_cells"].to_numpy(float)
        rec = {"sample_id": s, "condition": COND_OF[s], "bin": bi,
               "radial_centre": float(centres[bi]),
               "n_cells": int(w.sum()), "n_structures": int(len(b))}
        for p in PHENOTYPE_ORDER:
            rec[p] = float(np.average(b[p].to_numpy(float), weights=w))
        anim_rows.append(rec)
anim_prof = pd.DataFrame(anim_rows)
write_csv(anim_prof, "41b_radial_profiles_per_animal.csv")

print(f"    {len(prof)} structure-bin rows, {len(anim_prof)} animal-bin rows")

# ---- F38 radial composition -------------------------------------------------
show = [p for p in PHENOTYPE_ORDER if p != "Other"]
ncol, nrow = 4, int(np.ceil(len(show) / 4))
fig, axes = plt.subplots(nrow, ncol, figsize=(9.5 * ncol, 8 * nrow))
axes = np.atleast_1d(axes).ravel()
for i, p in enumerate(show):
    ax = axes[i]
    ax.axvspan(0, 1, color=CORE_BAND, zorder=0)
    for s in SAMPLE_ORDER:
        d = anim_prof.loc[anim_prof["sample_id"] == s].sort_values("radial_centre")
        if not len(d):
            continue
        ax.plot(d["radial_centre"], d[p], linewidth=4, marker=MARKER_OF[s],
                markersize=12, color=COLOR_OF[s], alpha=0.9)
    ax.axvline(1.0, color="#000000", linestyle="--", linewidth=3)
    ax.set_title(p, fontsize=FONT_SIZE_BASE - 4)
    ax.set_xlim(RADIAL_MIN, RADIAL_MAX)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("radial position", fontsize=FONT_SIZE_BASE - 10)
    ax.set_ylabel("% of cells in bin", fontsize=FONT_SIZE_BASE - 12)
    ax.tick_params(labelsize=FONT_SIZE_TICK - 12)
    style_axes(ax)
for i in range(len(show), len(axes)):
    axes[i].axis("off")
handles = [Line2D([0], [0], color=COLOR_OF[s], marker=MARKER_OF[s], markersize=16,
                  linewidth=4, label=f"{short_label(s)} ({COND_OF[s]})")
           for s in SAMPLE_ORDER]
fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False,
           fontsize=FONT_SIZE_LEGEND - 10, bbox_to_anchor=(0.5, -0.02))
fig.suptitle("Radial composition within structures\n"
             "0 = core centre, 1 = core boundary (dashed), 2 = outer cuff edge. "
             "Shaded region is core.", y=1.02, fontsize=FONT_SIZE_TITLE - 4)
save_fig(fig, "F38_radial_composition_profiles")


# %% Cell 5 - B: radial centre of mass with permutation null
# =============================================================================

banner("B - RADIAL CENTRE OF MASS")

print("    For each phenotype in each structure: the mean radial position of")
print("    that phenotype, against a null built by shuffling phenotype labels")
print("    within the same structure. Negative z = pulled toward the core,")
print("    positive z = displaced toward the rim.\n")

rcm_rows = []
for s, k, g in structure_frames():
    r = pd.to_numeric(g["radial_pos"], errors="coerce").to_numpy()
    ok = np.isfinite(r)
    r, ph = r[ok], g["pheno"].to_numpy()[ok]
    n = len(r)
    if n < 2:
        continue
    for p in PHENOTYPE_ORDER:
        m = ph == p
        npos = int(m.sum())
        if npos == 0:
            continue
        obs = float(r[m].mean())
        null = np.empty(N_PERMUTATIONS, float)
        for it in range(N_PERMUTATIONS):
            idx = rng.choice(n, size=npos, replace=False)
            null[it] = r[idx].mean()
        pval, z = empirical_p(obs, null)
        rcm_rows.append({
            "sample_id": s, "condition": COND_OF[s], "structure_id": k,
            "phenotype": p, "n_cells_phenotype": npos, "n_cells_structure": n,
            "observed_rcm": obs, "null_mean_rcm": float(null.mean()),
            "z": z, "p_empirical": pval,
            "delta_rcm": obs - float(null.mean()),
            "tier": tier_of(npos),
        })

rcm = pd.DataFrame(rcm_rows)
write_csv(rcm, "42_radial_centre_of_mass.csv")

sub("Animal-level summary (median z across structures, solid + usable tiers)")
rcm_use = rcm.loc[rcm["tier"].isin(["solid", "usable"])]
rcm_anim = (rcm_use.groupby(["sample_id", "condition", "phenotype"])
            .agg(median_z=("z", "median"),
                 median_delta=("delta_rcm", "median"),
                 n_structures=("z", "size")).reset_index())
write_csv(rcm_anim, "42b_radial_centre_of_mass_by_animal.csv")

key_ph = [IDO1_POS, IDO1_NEG, "Helper T cells", "CD4- T cells", "Tregs",
          "B cells", "Neutrophils"]
print(f"    {'phenotype':<26}" + "".join(f"{short_label(s):>11}" for s in SAMPLE_ORDER))
print("    " + "-" * (26 + 11 * len(SAMPLE_ORDER)))
for p in key_ph:
    row = f"    {p:<26}"
    for s in SAMPLE_ORDER:
        v = rcm_anim.loc[(rcm_anim["sample_id"] == s) &
                         (rcm_anim["phenotype"] == p), "median_z"]
        row += f"{v.iloc[0]:>11.2f}" if len(v) else f"{'na':>11}"
    print(row)

# ---- F39 centre of mass by arm ---------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(32, 14),
                         gridspec_kw={"width_ratios": [1.3, 1.0]})
ax = axes[0]
yy = np.arange(len(key_ph))
for i, p in enumerate(key_ph):
    for s in SAMPLE_ORDER:
        v = rcm_anim.loc[(rcm_anim["sample_id"] == s) &
                         (rcm_anim["phenotype"] == p), "median_z"]
        if not len(v):
            continue
        off = (SAMPLE_ORDER.index(s) - (len(SAMPLE_ORDER) - 1) / 2) * 0.13
        ax.scatter(v.iloc[0], i + off, s=420, color=COLOR_OF[s],
                   marker=MARKER_OF[s], edgecolor="#FFFFFF", linewidth=2,
                   zorder=3)
ax.axvline(0, color="#000000", linewidth=3.5)
ax.set_yticks(yy); ax.set_yticklabels(key_ph, fontsize=FONT_SIZE_TICK - 6)
ax.invert_yaxis()
ax.set_xlabel("median z of radial centre of mass\n(negative = core, positive = rim)")
ax.set_title("Radial preference", fontsize=FONT_SIZE_TITLE - 8)
ax.xaxis.grid(True, color=GRID_COLOR, linewidth=1.5); ax.set_axisbelow(True)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.legend(handles=[Line2D([0], [0], color=COLOR_OF[s], marker=MARKER_OF[s],
                          markersize=16, linestyle="none",
                          label=f"{short_label(s)} ({COND_OF[s]})")
                   for s in SAMPLE_ORDER],
          frameon=False, fontsize=FONT_SIZE_LEGEND - 12, loc="best")

ax = axes[1]
for i, p in enumerate(key_ph):
    for ci, c in enumerate(CONDITION_ORDER):
        v = rcm_anim.loc[(rcm_anim["condition"] == c) &
                         (rcm_anim["phenotype"] == p), "median_z"]
        if not len(v):
            continue
        off = (ci - 0.5) * 0.22
        ax.scatter(v, np.full(len(v), i + off), s=300,
                   color=CONDITION_COLORS[c], edgecolor="#FFFFFF", linewidth=2,
                   zorder=3, alpha=0.9)
        ax.hlines(i + off, v.min(), v.max(), color=CONDITION_COLORS[c],
                  linewidth=4, alpha=0.5, zorder=2)
ax.axvline(0, color="#000000", linewidth=3.5)
ax.set_yticks(yy); ax.set_yticklabels([]); ax.invert_yaxis()
ax.set_xlabel("median z by arm")
ax.set_title("By arm", fontsize=FONT_SIZE_TITLE - 8)
ax.xaxis.grid(True, color=GRID_COLOR, linewidth=1.5); ax.set_axisbelow(True)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.legend(handles=[Patch(facecolor=CONDITION_COLORS[c], label=c)
                   for c in CONDITION_ORDER],
          frameon=False, fontsize=FONT_SIZE_LEGEND - 8, loc="best")
fig.suptitle("Where each population sits within a structure\n"
             "z is against a within-structure label permutation null, so "
             "structure size and cell density are already accounted for",
             y=1.04, fontsize=FONT_SIZE_TITLE - 6)
save_fig(fig, "F39_radial_centre_of_mass")


# %% Cell 6 - C: within-focus nearest-neighbour distances
# =============================================================================

banner("C - WITHIN-STRUCTURE NEAREST-NEIGHBOUR DISTANCES")

print("    Per anchor cell: distance to the nearest target cell in the same")
print("    structure. Null: reassign anchor and target labels at random among")
print("    the cells of that structure, preserving both counts and all")
print("    coordinates. Negative z = closer than the tissue arrangement would")
print("    produce by chance.\n")

nn_rows = []
for s, k, g in structure_frames():
    xy = g[["x", "y"]].to_numpy(float)
    ph = g["pheno"].to_numpy()
    n = len(xy)
    if n < 10:
        continue
    kk = int(min(KNN_DEPTH, n - 1))
    if kk < 1:
        continue
    tree = cKDTree(xy)
    dists, idxs = tree.query(xy, k=kk + 1)
    dists, idxs = dists[:, 1:], idxs[:, 1:]     # drop self

    for anchor in NN_ANCHORS:
        a_mask = ph == anchor
        n_a = int(a_mask.sum())
        if n_a == 0:
            continue
        a_rows = np.flatnonzero(a_mask)
        if len(a_rows) > MAX_ANCHORS_PER_FOCUS:
            a_rows = rng.choice(a_rows, size=MAX_ANCHORS_PER_FOCUS, replace=False)
        for target in NN_TARGETS:
            if target == anchor:
                continue
            t_mask = ph == target
            n_t = int(t_mask.sum())
            if n_t == 0:
                continue

            def nn_from(tmask, rows):
                hit = tmask[idxs[rows]]
                dm = np.where(hit, dists[rows], np.inf)
                v = dm.min(axis=1)
                return v

            obs_v = nn_from(t_mask, a_rows)
            censored = float(np.mean(~np.isfinite(obs_v)))
            obs = float(np.median(obs_v[np.isfinite(obs_v)])) \
                if np.isfinite(obs_v).any() else np.nan

            null = np.full(N_PERMUTATIONS, np.nan)
            for it in range(N_PERMUTATIONS):
                perm = rng.permutation(n)
                fake_a = perm[:n_a]
                fake_t = np.zeros(n, dtype=bool)
                fake_t[perm[n_a:n_a + n_t]] = True
                rows = (rng.choice(fake_a, size=MAX_ANCHORS_PER_FOCUS,
                                   replace=False)
                        if len(fake_a) > MAX_ANCHORS_PER_FOCUS else fake_a)
                v = nn_from(fake_t, rows)
                v = v[np.isfinite(v)]
                if len(v):
                    null[it] = np.median(v)

            pval, z = empirical_p(obs, null)
            nn_rows.append({
                "sample_id": s, "condition": COND_OF[s], "structure_id": k,
                "anchor": anchor, "target": target,
                "n_anchor": n_a, "n_target": n_t, "n_structure": n,
                "observed_median_um": obs,
                "null_median_um": float(np.nanmean(null)),
                "delta_um": obs - float(np.nanmean(null)),
                "z": z, "p_empirical": pval,
                "frac_censored": censored,
                "tier": tier_of(min(n_a, n_t)),
            })
    del tree, dists, idxs
    gc.collect()

nn = pd.DataFrame(nn_rows)
if len(nn):
    write_csv(nn, "43_within_structure_nn_distances.csv")

    sub("Censoring check (anchors with no target inside the k-neighbour list)")
    bad = nn.loc[nn["frac_censored"] > 0.2]
    if len(bad):
        print(f"    {len(bad)} rows exceed 20% censoring. Raise KNN_DEPTH if these")
        print("    matter to the conclusion:")
        for _, r in bad.head(15).iterrows():
            print(f"      {r['sample_id']:<12} struct {int(r['structure_id']):>3} "
                  f"{r['anchor']:<24} -> {r['target']:<16} "
                  f"{100*r['frac_censored']:>5.1f}%")
    else:
        print("    none above 20%")

    nn_use = nn.loc[nn["tier"].isin(["solid", "usable"])]
    nn_anim = (nn_use.groupby(["sample_id", "condition", "anchor", "target"])
               .agg(median_observed_um=("observed_median_um", "median"),
                    median_z=("z", "median"),
                    median_delta_um=("delta_um", "median"),
                    n_structures=("z", "size")).reset_index())
    write_csv(nn_anim, "43b_nn_distances_by_animal.csv")

    sub("Animal-level nearest-neighbour summary (solid + usable tiers)")
    for anchor in NN_ANCHORS:
        print(f"\n    anchor: {anchor}")
        print(f"      {'target':<18}" +
              "".join(f"{short_label(s):>13}" for s in SAMPLE_ORDER))
        for target in NN_TARGETS:
            if target == anchor:
                continue
            row = f"      {target:<18}"
            for s in SAMPLE_ORDER:
                v = nn_anim.loc[(nn_anim["sample_id"] == s) &
                                (nn_anim["anchor"] == anchor) &
                                (nn_anim["target"] == target)]
                row += (f"{v['median_observed_um'].iloc[0]:>8.0f}um"
                        if len(v) else f"{'na':>13}")
            print(row)
        print(f"      {'(z vs null)':<18}" +
              "".join(f"{'':>13}" for s in SAMPLE_ORDER))
        for target in NN_TARGETS:
            if target == anchor:
                continue
            row = f"      {target:<18}"
            for s in SAMPLE_ORDER:
                v = nn_anim.loc[(nn_anim["sample_id"] == s) &
                                (nn_anim["anchor"] == anchor) &
                                (nn_anim["target"] == target)]
                row += (f"{v['median_z'].iloc[0]:>13.2f}" if len(v) else f"{'na':>13}")
            print(row)

    # ---- F40 nn results ----------------------------------------------------
    pairs = [(a, t) for a in NN_ANCHORS for t in NN_TARGETS if a != t]
    ncol = 4
    nrow = int(np.ceil(len(pairs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(10 * ncol, 8.5 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for i, (a, t) in enumerate(pairs):
        ax = axes[i]
        d = nn_use.loc[(nn_use["anchor"] == a) & (nn_use["target"] == t)]
        for ci, c in enumerate(CONDITION_ORDER):
            v = d.loc[d["condition"] == c, "z"].to_numpy()
            if not len(v):
                continue
            j = rng.uniform(-0.13, 0.13, size=len(v))
            ax.scatter(np.full(len(v), ci) + j, v, s=240,
                       color=CONDITION_COLORS[c], edgecolor="#FFFFFF",
                       linewidth=1.8, zorder=3, alpha=0.85)
            ax.hlines(np.median(v), ci - 0.3, ci + 0.3, color="#000000",
                      linewidth=4)
        ax.axhline(0, color="#000000", linestyle="--", linewidth=3)
        ax.set_xticks(range(len(CONDITION_ORDER)))
        ax.set_xticklabels(CONDITION_ORDER, fontsize=FONT_SIZE_TICK - 10)
        ax.set_xlim(-0.6, len(CONDITION_ORDER) - 0.4)
        ax.set_title(f"{a}\nto {t}", fontsize=FONT_SIZE_BASE - 10)
        ax.set_ylabel("z vs null", fontsize=FONT_SIZE_BASE - 12)
        ax.tick_params(labelsize=FONT_SIZE_TICK - 12)
        style_axes(ax)
    for i in range(len(pairs), len(axes)):
        axes[i].axis("off")
    fig.suptitle("Within-structure proximity, each point one structure\n"
                 "Negative z = anchor and target closer than the structure's own "
                 "cell arrangement predicts", y=1.02,
                 fontsize=FONT_SIZE_TITLE - 6)
    save_fig(fig, "F40_within_structure_nn_z")
else:
    nn_use = nn_anim = pd.DataFrame()
    print("    WARNING: no nearest-neighbour rows computed.")


# %% Cell 7 - D: iNOS vs Arginase-1 segregation
# =============================================================================

banner("D - iNOS vs ARGINASE-1 MACROPHAGE SEGREGATION")

print("    Polarisation is called by a WITHIN-SECTION percentile on macrophages,")
print("    never an absolute intensity, because iNOS (between-slide fraction")
print("    0.862) and Arginase-1 (0.726) are not comparable across the two")
print("    scans. Thresholds are swept.\n")

pol_rows, mix_rows = [], []
for s in SAMPLE_ORDER:
    d = cells[s]
    mac = d.loc[d["pheno"].isin(MACROPHAGE_ALL)]
    if not len(mac):
        continue
    for pct in POLARISATION_PERCENTILES:
        inos_thr = float(np.nanpercentile(mac["iNOS"], pct))
        arg_thr = float(np.nanpercentile(mac["Arginase-1"], pct))
        hi_i = (mac["iNOS"] >= inos_thr).to_numpy()
        hi_a = (mac["Arginase-1"] >= arg_thr).to_numpy()
        pol_rows.append({
            "sample_id": s, "condition": COND_OF[s], "percentile": pct,
            "inos_threshold": inos_thr, "arg1_threshold": arg_thr,
            "n_macrophages": len(mac),
            "pct_inos_only": 100.0 * float((hi_i & ~hi_a).mean()),
            "pct_arg1_only": 100.0 * float((hi_a & ~hi_i).mean()),
            "pct_double": 100.0 * float((hi_i & hi_a).mean()),
            "pct_double_negative": 100.0 * float((~hi_i & ~hi_a).mean()),
        })

    # mixing score at the primary percentile, per structure
    inos_thr = float(np.nanpercentile(mac["iNOS"], POLARISATION_PRIMARY))
    arg_thr = float(np.nanpercentile(mac["Arginase-1"], POLARISATION_PRIMARY))
    for k, g in d.loc[(d[STRUCT_COL] > 0) &
                      d["pheno"].isin(MACROPHAGE_ALL)].groupby(STRUCT_COL):
        n = len(g)
        if n < 30:
            continue
        hi_i = (g["iNOS"] >= inos_thr).to_numpy()
        hi_a = (g["Arginase-1"] >= arg_thr).to_numpy()
        ex_i = hi_i & ~hi_a
        ex_a = hi_a & ~hi_i
        if ex_i.sum() < 5 or ex_a.sum() < 5:
            continue
        xy = g[["x", "y"]].to_numpy(float)
        kk = int(min(MIXING_K, n - 1))
        tree = cKDTree(xy)
        _, idxs = tree.query(xy, k=kk + 1)
        idxs = idxs[:, 1:]

        def mixing(a_mask, b_mask):
            rows = np.flatnonzero(a_mask)
            if not len(rows):
                return np.nan
            return float(b_mask[idxs[rows]].mean())

        obs = mixing(ex_i, ex_a)
        null = np.empty(N_PERMUTATIONS, float)
        n_i, n_a = int(ex_i.sum()), int(ex_a.sum())
        for it in range(N_PERMUTATIONS):
            perm = rng.permutation(n)
            fi = np.zeros(n, bool); fi[perm[:n_i]] = True
            fa = np.zeros(n, bool); fa[perm[n_i:n_i + n_a]] = True
            null[it] = mixing(fi, fa)
        pval, z = empirical_p(obs, null)
        mix_rows.append({
            "sample_id": s, "condition": COND_OF[s], "structure_id": int(k),
            "n_macrophages": n, "n_inos_only": n_i, "n_arg1_only": n_a,
            "observed_mixing": obs, "null_mixing": float(np.nanmean(null)),
            "z": z, "p_empirical": pval,
            "tier": tier_of(min(n_i, n_a)),
        })
        del tree, idxs

pol = pd.DataFrame(pol_rows)
mix = pd.DataFrame(mix_rows)
if len(pol):
    write_csv(pol, "44_polarisation_composition.csv")
    sub(f"Macrophage polarisation composition at the {POLARISATION_PRIMARY}th "
        f"within-section percentile")
    p0 = pol.loc[pol["percentile"] == POLARISATION_PRIMARY]
    print(p0[["sample_id", "condition", "n_macrophages", "pct_inos_only",
              "pct_arg1_only", "pct_double", "pct_double_negative"]]
          .to_string(index=False))
    print("\n    Reminder: a within-section percentile fixes the marginal counts")
    print("    by construction, so the informative quantities are the DOUBLE and")
    print("    DOUBLE-NEGATIVE fractions and the spatial mixing below, not the")
    print("    single-positive percentages.")
if len(mix):
    write_csv(mix, "45_inos_arg1_mixing.csv")
    sub("Spatial mixing of iNOS-high and Arginase-1-high macrophages")
    print("    Negative z = the two states are more segregated than chance.")
    mix_anim = (mix.loc[mix["tier"] != "not_interpretable"]
                .groupby(["sample_id", "condition"])
                .agg(median_z=("z", "median"),
                     median_observed=("observed_mixing", "median"),
                     n_structures=("z", "size")).reset_index())
    print(mix_anim.to_string(index=False))
    write_csv(mix_anim, "45b_inos_arg1_mixing_by_animal.csv")

    fig, axes = plt.subplots(1, 2, figsize=(30, 13))
    ax = axes[0]
    for ci, c in enumerate(CONDITION_ORDER):
        v = mix.loc[mix["condition"] == c, "z"].dropna().to_numpy()
        if not len(v):
            continue
        j = rng.uniform(-0.14, 0.14, size=len(v))
        ax.scatter(np.full(len(v), ci) + j, v, s=300, color=CONDITION_COLORS[c],
                   edgecolor="#FFFFFF", linewidth=2, zorder=3, alpha=0.85)
        ax.hlines(np.median(v), ci - 0.3, ci + 0.3, color="#000000", linewidth=4)
    ax.axhline(0, color="#000000", linestyle="--", linewidth=3)
    ax.set_xticks(range(len(CONDITION_ORDER)))
    ax.set_xticklabels(CONDITION_ORDER)
    ax.set_ylabel("z of mixing vs null\n(negative = segregated)")
    ax.set_title("iNOS-high / Arg1-high spatial mixing",
                 fontsize=FONT_SIZE_TITLE - 12)
    style_axes(ax)

    ax = axes[1]
    for s in SAMPLE_ORDER:
        d = pol.loc[pol["sample_id"] == s].sort_values("percentile")
        if not len(d):
            continue
        ax.plot(d["percentile"], d["pct_double"], linewidth=4,
                marker=MARKER_OF[s], markersize=14, color=COLOR_OF[s], alpha=0.9)
    ax.set_xlabel("within-section percentile threshold")
    ax.set_ylabel("% macrophages double positive")
    ax.set_title("iNOS and Arg1 co-expression", fontsize=FONT_SIZE_TITLE - 12)
    style_axes(ax)
    ax.legend(handles=[Line2D([0], [0], color=COLOR_OF[s], marker=MARKER_OF[s],
                              markersize=16, linewidth=4,
                              label=f"{short_label(s)} ({COND_OF[s]})")
                       for s in SAMPLE_ORDER],
              frameon=False, fontsize=FONT_SIZE_LEGEND - 14, loc="best")
    fig.suptitle("Macrophage polarisation state and its spatial organisation",
                 y=1.03, fontsize=FONT_SIZE_TITLE - 6)
    save_fig(fig, "F41_inos_arg1_segregation")


# %% Cell 8 - E: per-focus IDO1 vs 3-HK preview
# =============================================================================

banner("E - PER-STRUCTURE IDO1 vs 3-HK (PREVIEW ONLY)")

print("    This uses per-structure medians already computed in script 04. It is")
print("    NOT the spatial gradient analysis, which is script 06. Included here")
print("    because 3-HK has the most favourable between-slide comparability")
print("    (variance fraction 0.117) of any marker carrying a core question.\n")

foci_tbl = None
if os.path.exists(FOCI_TABLE):
    foci_tbl = pd.read_csv(FOCI_TABLE)
    keep = ["sample_id", "condition", "focus_id", "peak_density",
            "fold_over_background", "equiv_radius_um",
            "pct_ido1_pos_of_mac_core", "median_ido1_macrophages",
            "median_hk3_macrophages", "median_inos_macrophages",
            "median_arg1_macrophages", "median_ifng_core"]
    keep = [c for c in keep if c in foci_tbl.columns]
    prev = foci_tbl[keep].copy()
    print(prev.groupby("condition")[
        [c for c in ["median_ido1_macrophages", "median_hk3_macrophages",
                     "pct_ido1_pos_of_mac_core"] if c in prev.columns]]
        .median().to_string())
    write_csv(prev, "46_structure_marker_preview.csv")

    if ("median_ido1_macrophages" in prev.columns
            and "median_hk3_macrophages" in prev.columns):
        fig, axes = plt.subplots(1, 2, figsize=(30, 13))
        ax = axes[0]
        for s in SAMPLE_ORDER:
            d = prev.loc[prev["sample_id"] == s]
            if not len(d):
                continue
            ax.scatter(d["median_ido1_macrophages"], d["median_hk3_macrophages"],
                       s=340, color=COLOR_OF[s], marker=MARKER_OF[s],
                       edgecolor="#FFFFFF", linewidth=2, zorder=3)
        ax.set_xscale("log")
        ax.set_xlabel("median IDO1 in core macrophages")
        ax.set_ylabel("median 3-hydroxykynurenine in core macrophages")
        ax.set_title("Enzyme against product, per structure",
                     fontsize=FONT_SIZE_TITLE - 14)
        style_axes(ax)
        ax.legend(handles=[Line2D([0], [0], color=COLOR_OF[s],
                                  marker=MARKER_OF[s], markersize=16,
                                  linestyle="none",
                                  label=f"{short_label(s)} ({COND_OF[s]})")
                           for s in SAMPLE_ORDER],
                  frameon=False, fontsize=FONT_SIZE_LEGEND - 14, loc="best")

        ax = axes[1]
        for ci, c in enumerate(CONDITION_ORDER):
            v = prev.loc[prev["condition"] == c, "median_hk3_macrophages"].dropna()
            if not len(v):
                continue
            j = rng.uniform(-0.14, 0.14, size=len(v))
            ax.scatter(np.full(len(v), ci) + j, v, s=320,
                       color=CONDITION_COLORS[c], edgecolor="#FFFFFF",
                       linewidth=2, zorder=3, alpha=0.85)
            ax.hlines(np.median(v), ci - 0.3, ci + 0.3, color="#000000",
                      linewidth=4)
        ax.set_xticks(range(len(CONDITION_ORDER)))
        ax.set_xticklabels(CONDITION_ORDER)
        ax.set_ylabel("median 3-HK in core macrophages")
        ax.set_title("3-HK by arm", fontsize=FONT_SIZE_TITLE - 12)
        style_axes(ax)
        fig.suptitle("Preview: IDO1 protein against its downstream product\n"
                     "D1MT inhibits the enzyme, not transcription, so IDO1 "
                     "present with 3-HK reduced is the mechanistic prediction",
                     y=1.04, fontsize=FONT_SIZE_TITLE - 8)
        save_fig(fig, "F42_ido1_vs_hk3_preview")
else:
    print(f"    WARNING: {FOCI_TABLE} not found, preview skipped")


# %% Cell 9 - consistency across animals
# =============================================================================

banner("CONSISTENCY ACROSS ANIMALS")

print("    With three animals per arm, a rank test cannot go below p = 0.10.")
print("    Direction and consistency carry the weight. 'k of 3' counts animals")
print("    whose median effect runs in the stated direction.\n")

cons_rows = []


def consistency(df, group_cols, value_col, label):
    for keys, g in df.groupby(group_cols):
        keys = keys if isinstance(keys, tuple) else (keys,)
        rec = {"analysis": label}
        for cname, kv in zip(group_cols, keys):
            rec[cname] = kv
        for c in CONDITION_ORDER:
            v = g.loc[g["condition"] == c, value_col].dropna()
            n_tot = len(v)
            rec[f"{c}_n"] = n_tot
            rec[f"{c}_median"] = float(v.median()) if n_tot else np.nan
            rec[f"{c}_k_negative"] = int((v < 0).sum())
            rec[f"{c}_k_positive"] = int((v > 0).sum())
        a = rec.get(f"{CONDITION_ORDER[0]}_median", np.nan)
        b = rec.get(f"{CONDITION_ORDER[1]}_median", np.nan)
        rec["arm_difference"] = a - b if np.isfinite(a) and np.isfinite(b) else np.nan
        cons_rows.append(rec)


if len(rcm_anim):
    consistency(rcm_anim, ["phenotype"], "median_z", "radial_centre_of_mass")
if len(nn_anim):
    consistency(nn_anim, ["anchor", "target"], "median_z", "nn_proximity")
if len(mix):
    m2 = (mix.loc[mix["tier"] != "not_interpretable"]
          .groupby(["sample_id", "condition"])["z"].median().reset_index())
    m2["phenotype"] = "iNOS_vs_Arg1_mixing"
    consistency(m2, ["phenotype"], "z", "polarisation_mixing")

cons = pd.DataFrame(cons_rows)
if len(cons):
    write_csv(cons, "47_consistency_summary.csv")
    sub("Effects with a consistent direction in all animals of at least one arm")
    for _, r in cons.iterrows():
        d1_n = r.get("D1MT_n", 0) or 0
        un_n = r.get("Untreated_n", 0) or 0
        d1_consistent = d1_n > 0 and (r.get("D1MT_k_negative", 0) == d1_n
                                      or r.get("D1MT_k_positive", 0) == d1_n)
        un_consistent = un_n > 0 and (r.get("Untreated_k_negative", 0) == un_n
                                      or r.get("Untreated_k_positive", 0) == un_n)
        if not (d1_consistent or un_consistent):
            continue
        name = " / ".join(str(r[c]) for c in ["phenotype", "anchor", "target"]
                          if c in r.index and pd.notna(r.get(c)))
        print(f"    [{r['analysis']:<22}] {name:<48} "
              f"D1MT {r.get('D1MT_median', np.nan):>6.2f} "
              f"({int(r.get('D1MT_k_negative', 0))}-/{int(d1_n)})   "
              f"Untr {r.get('Untreated_median', np.nan):>6.2f} "
              f"({int(r.get('Untreated_k_negative', 0))}-/{int(un_n)})")


# %% Cell 10 - wrap up
# =============================================================================

banner("SUMMARY")
print(f"Structure source     : {STRUCTURE_SOURCE}")
print(f"Structures analysed  : {len(inv)}")
print(f"Radial COM rows      : {len(rcm)}")
print(f"NN rows              : {len(nn)}")
print(f"Mixing rows          : {len(mix)}")
print(f"Permutations         : {N_PERMUTATIONS}")

sub("Read in this order")
print("  1. F38 : radial composition. Does the coordinate behave sensibly?")
print("     Myeloid cells should peak at radial 0 by construction. If T and B")
print("     cells show structure beyond that, the coordinate is informative.")
print("  2. F39 : radial centre of mass. THE Q1 RESULT. Negative z = core,")
print("     positive z = rim, already corrected for structure size and density.")
print("  3. F40 : within-structure proximity. Q1 and Q4 as distances.")
print("  4. F41 : iNOS / Arg1 segregation. Q6a.")
print("  5. F42 : IDO1 vs 3-HK preview. Sets up script 06.")
print("  6. Table 47 : which effects are consistent across all animals of an arm.")

sub("Interpretation guardrails")
print("  - A negative z means closer or more core-ward THAN THAT STRUCTURE'S OWN")
print("    cell arrangement predicts. It is not a raw distance and should not be")
print("    described as one.")
print("  - The treated arm has few structures and few IDO1+ cells. Check the")
print("    tier column before quoting any treated-arm number.")
print("  - Helper T cell counts are unstable across identically stained sections")
print("    (26 to 4,133). Read the CD4- T cell rows alongside as a check.")
print("  - Polarisation calls use within-section percentiles, so single-positive")
print("    percentages are fixed by construction. Only co-expression and spatial")
print("    mixing are informative.")

sub("Script 06")
print("  Marker gradients: 3-HK as a function of distance to the nearest")
print("  IDO1-high macrophage (Q3), and IDO1 intensity as a function of distance")
print("  to the nearest IFNG-high cell (Q2).")

banner("DONE")
sys.stdout = _tee.terminal
_tee.close()
