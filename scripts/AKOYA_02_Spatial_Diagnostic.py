#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - SPATIAL LAYOUT AND BATCH DIAGNOSTICS
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 02 of the AKOYA analysis series.

PURPOSE
    Answer the structural questions raised by the script 01 inventory, before any
    gating, normalization or spatial statistics are committed to.

    Q1  Do the four sections per scan tile cleanly in a shared coordinate space?
    Q2  How much of the marker intensity spread is between slide versus between
        section within a slide? (Slide is perfectly confounded with condition, so
        this bounds the risk rather than resolving it.)
    Q3  Are 43102 and 43112 (bottom section of each scan) globally low across many
        markers, which would indicate a slide position artifact rather than biology?
    Q4  What does the IDO1 intensity distribution on CD68+ macrophages actually
        look like per section? Is there a distinct IDO1-high mode, and does it
        vanish in the treated arm or just shift?
    Q5  What is sitting in the "Other" bucket?
    Q6  For every anchor / target phenotype pair, how many animals have enough
        cells to support a nearest-neighbour test?

    Still diagnostic only. Nothing is re-phenotyped, gated, normalized or written
    back to the source CSVs.

OUTPUTS
    figures/  F10 .. F16   (PDF + PNG, 300 DPI)
    tables/   10 .. 15     (CSV, LF line endings)

USAGE
    conda activate sc_pre
    python AKOYA_02_Spatial_Diagnostics.py

    Or via SLURM (single job, no arrays):
    #SBATCH --partition=normal
    #SBATCH --cpus-per-task=8
    #SBATCH --mem=200G
    #SBATCH --output=/master/jlehle/WORKING/LOGS/akoya02_%j.out

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================
# ALL TUNABLE PARAMETERS LIVE HERE
# =============================================================================

DATA_DIR = "/master/jlehle/WORKING/AKOYA/data"
INVENTORY_TABLE_DIR = "/master/jlehle/WORKING/AKOYA/inventory/tables"
OUT_DIR = "/master/jlehle/WORKING/AKOYA/diagnostics"

CONDITION_ORDER = ["D1MT", "Untreated"]
CONDITION_COLORS = {"D1MT": "#2C7FB8", "Untreated": "#D95F02"}

# Scan / slide identity is parsed from the Image column. Group prefix fallback
# is used if the Image column is missing.
GROUP_MAP = {"G3": "D1MT", "G4": "Untreated"}

# ---- column names -----------------------------------------------------------
INTENSITY_SUFFIX = ": Mean"
CENTROID_X_PREFIX = "Centroid X"
CENTROID_Y_PREFIX = "Centroid Y"
PHENOTYPE_COL = "Phenotypes"
PARENT_COL = "Parent"
IMAGE_COL = "Image"

IDO1_COL = "IDO1: Cytoplasm: Mean"
CD68_COL = "CD68: Membrane: Mean"

# ---- phenotype handling -----------------------------------------------------
# ASCII-safe display names. The script maps these onto the exact strings in the
# data, which carry a zero-width space on the epithelial label.
PHENOTYPE_ORDER = [
    "CD68+IDO1+ Macrophages", "CD68+IDO1- Macrophages", "CD163+ Macrophages",
    "Neutrophils", "Helper T cells", "CD4- T cells", "Tregs",
    "B cells", "Plasma cells", "Endothelial cells",
    "Epithelial/Tumor cells", "Other",
]
PHENOTYPE_COLORS = {
    "CD68+IDO1+ Macrophages": "#B2182B", "CD68+IDO1- Macrophages": "#EF8A62",
    "CD163+ Macrophages": "#FDDBC7", "Neutrophils": "#7B3294",
    "Helper T cells": "#1B7837", "CD4- T cells": "#7FBC41", "Tregs": "#276419",
    "B cells": "#2166AC", "Plasma cells": "#67A9CF",
    "Endothelial cells": "#BABABA", "Epithelial/Tumor cells": "#878787",
    "Other": "#4D4D4D",
}

# Phenotypes counted as CD68 lineage for the IDO1 distribution work
CD68_LINEAGE = ["CD68+IDO1+ Macrophages", "CD68+IDO1- Macrophages"]

# Phenotypes highlighted on the spatial maps (everything else drawn light grey)
SPATIAL_HIGHLIGHT = [
    "CD68+IDO1+ Macrophages", "CD68+IDO1- Macrophages",
    "Helper T cells", "CD4- T cells", "B cells",
]
BACKGROUND_COLOR = "#E8E8E8"

# ---- nearest-neighbour feasibility -----------------------------------------
NN_ANCHORS = [
    "CD68+IDO1+ Macrophages", "CD68+IDO1- Macrophages", "CD163+ Macrophages",
]
NN_TARGETS = [
    "Helper T cells", "CD4- T cells", "Tregs", "B cells", "Plasma cells",
    "Neutrophils", "CD68+IDO1+ Macrophages", "CD68+IDO1- Macrophages",
]
MIN_CELLS_FOR_NN = 50          # per animal, for both anchor and target

# ---- position artifact test -------------------------------------------------
N_Y_BINS_WITHIN_SECTION = 10   # within-section Y deciles
MIN_CELLS_PER_YBIN = 200       # bins thinner than this are dropped

# ---- markers highlighted in diagnostics ------------------------------------
KEY_MARKERS = [
    "IDO1", "CD68", "CD163", "CD206", "CD4", "CD8", "CD3e", "CD20", "CD79a",
    "CD21", "MPO", "CD11b", "FoxP3", "iNOS", "Arginase-1", "IFNG",
    "Granzyme-B", "HLA-DR", "PD-1", "PD-L1", "3-Hydroxykynurenine",
]

# ---- plotting ---------------------------------------------------------------
MAX_POINTS_PER_SECTION = 60000   # downsample for scatter rendering only
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
FLAG_COLOR = "#B2182B"

USE_FLOAT32 = True   # halve memory on the intensity columns


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
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

FIG_DIR = os.path.join(OUT_DIR, "figures")
TAB_DIR = os.path.join(OUT_DIR, "tables")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TAB_DIR, exist_ok=True)

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


def banner(text, char="=", width=88):
    print("\n" + char * width); print(text); print(char * width)


def sub(text):
    print("\n" + text); print("-" * min(len(text), 88))


def ascii_safe(s):
    if not isinstance(s, str):
        return str(s)
    return "".join(c for c in s if ord(c) <= 127).strip()


def short_label(sid):
    return sid.split("_")[-1] if "_" in sid else sid


def marker_base(col):
    return col.split(":")[0].strip()


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


def find_col(cols, prefix):
    hits = [c for c in cols if c.startswith(prefix)]
    return hits[0] if hits else None


_tee = Tee(os.path.join(TAB_DIR, "00_diagnostics_report.txt"))
sys.stdout = _tee

banner("AKOYA SPATIAL LAYOUT AND BATCH DIAGNOSTICS")
print(f"Run time  : {datetime.now().isoformat(timespec='seconds')}")
print(f"Data      : {DATA_DIR}")
print(f"Figures   : {FIG_DIR}")
print(f"Tables    : {TAB_DIR}")


# %% Cell 3 - Q2: between-slide vs within-slide intensity spread
# =============================================================================
# Uses the p99 values already computed in script 01. No raw reads needed.

banner("Q2 - BETWEEN-SLIDE VS WITHIN-SLIDE INTENSITY SPREAD")

ms_path = os.path.join(INVENTORY_TABLE_DIR, "06_marker_stats_long.csv")
fs_path = os.path.join(INVENTORY_TABLE_DIR, "01_file_summary.csv")

var_tbl = pd.DataFrame()
if not (os.path.exists(ms_path) and os.path.exists(fs_path)):
    print(f"    WARNING: inventory tables not found. Skipping Q2.")
else:
    ms = pd.read_csv(ms_path)
    fs = pd.read_csv(fs_path)
    slide_of = dict(zip(fs["sample_id"], fs["treatment"]))

    ms["slide"] = ms["sample_id"].map(slide_of)
    ms = ms.loc[ms["p99"].notna() & (ms["p99"] > 0)].copy()
    ms["log2_p99"] = np.log2(ms["p99"])

    rows = []
    for col, g in ms.groupby("column"):
        slide_means = g.groupby("slide")["log2_p99"].mean()
        if len(slide_means) < 2:
            continue
        # within-slide variance, pooled across slides
        within = g.groupby("slide")["log2_p99"].var(ddof=1).mean()
        between = float(np.var(slide_means.to_numpy(), ddof=1))
        total = within + between
        rows.append({
            "column": col,
            "marker": marker_base(col),
            "is_key_marker": marker_base(col) in KEY_MARKERS,
            "within_slide_var": within,
            "between_slide_var": between,
            "frac_between_slide": between / total if total > 0 else np.nan,
            "log2_slide_diff": float(slide_means.get(CONDITION_ORDER[1], np.nan)
                                     - slide_means.get(CONDITION_ORDER[0], np.nan)),
            "within_slide_sd": np.sqrt(within),
        })
    var_tbl = pd.DataFrame(rows).sort_values("frac_between_slide", ascending=False)

    print(f"    {len(var_tbl)} markers decomposed\n")
    print("    IMPORTANT: slide is perfectly confounded with condition (one scan")
    print("    per group). A high between-slide fraction does NOT prove batch")
    print("    effect, and a low one does NOT prove absence of biology. This only")
    print("    flags which markers carry risk if a shared threshold is applied.\n")

    sub("Markers most dominated by between-slide variation (top 20)")
    for _, r in var_tbl.head(20).iterrows():
        star = " *" if r["is_key_marker"] else ""
        print(f"    {r['column']:<48} frac_between={r['frac_between_slide']:.3f}  "
              f"log2 diff={r['log2_slide_diff']:+.2f}{star}")

    sub("Key phenotype-calling markers")
    for _, r in var_tbl.loc[var_tbl["is_key_marker"]].iterrows():
        print(f"    {r['column']:<48} frac_between={r['frac_between_slide']:.3f}  "
              f"log2 diff={r['log2_slide_diff']:+.2f}  "
              f"within-slide sd={r['within_slide_sd']:.2f}")

    # ---- figure F10 ---------------------------------------------------------
    plot_tbl = var_tbl.copy()
    fig, ax = plt.subplots(figsize=(18, 14))
    ax.scatter(plot_tbl["within_slide_sd"], plot_tbl["frac_between_slide"],
               s=200, color="#BBBBBB", edgecolor="#FFFFFF", linewidth=1.5,
               zorder=3, label="all markers")
    key = plot_tbl.loc[plot_tbl["is_key_marker"]]
    ax.scatter(key["within_slide_sd"], key["frac_between_slide"],
               s=520, color="#B2182B", edgecolor="#FFFFFF", linewidth=2.5,
               zorder=4, label="phenotype-calling markers")
    for _, r in key.iterrows():
        ax.annotate(r["marker"], (r["within_slide_sd"], r["frac_between_slide"]),
                    textcoords="offset points", xytext=(12, 8),
                    fontsize=FONT_SIZE_ANNOT - 8)
    ax.axhline(0.5, color="#000000", linestyle="--", linewidth=3)
    ax.set_xlabel("Within-slide spread of log$_2$(p99), SD")
    ax.set_ylabel("Fraction of variance between slides")
    ax.set_ylim(0, 1.02)
    ax.set_title("Where marker intensity variation lives\n"
                 "High on the y axis = confounded with condition, cannot be "
                 "separated from batch",
                 fontsize=FONT_SIZE_TITLE - 6)
    style_axes(ax)
    ax.legend(loc="lower right", frameon=False, fontsize=FONT_SIZE_LEGEND - 4)
    save_fig(fig, "F10_variance_between_vs_within_slide")

    write_csv(var_tbl, "10_variance_decomposition.csv")


# %% Cell 4 - single pass over the raw CSVs
# =============================================================================

banner("LOADING RAW SECTION DATA (single pass)")

csv_paths = sorted(glob.glob(os.path.join(DATA_DIR, "*.csv")))
if not csv_paths:
    print(f"ERROR: no CSVs in {DATA_DIR}")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

# header pass to build a float32 dtype map
hdr = pd.read_csv(csv_paths[0], nrows=0)
ALL_COLS = hdr.columns.tolist()
XCOL = find_col(ALL_COLS, CENTROID_X_PREFIX)
YCOL = find_col(ALL_COLS, CENTROID_Y_PREFIX)
INTENSITY_COLS = [c for c in ALL_COLS if c.endswith(INTENSITY_SUFFIX)]
if XCOL is None or YCOL is None:
    print("ERROR: centroid columns not found.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

dtype_map = {c: "float32" for c in INTENSITY_COLS + [XCOL, YCOL]} if USE_FLOAT32 else None
print(f"    {len(INTENSITY_COLS)} intensity columns, X='{XCOL}', Y='{YCOL}'")

sections = {}          # sid -> dict of light-weight arrays and summaries
pheno_profiles = []    # tidy per-section per-phenotype mean marker profile
ybin_rows = []         # within-section Y-decile marker medians
section_marker_med = []  # per-section median for every marker
ido_cd68 = {}          # sid -> np.array of IDO1 values on CD68-lineage cells

for path in csv_paths:
    sid = os.path.splitext(os.path.basename(path))[0]
    prefix = sid.split("_")[0]
    cond = GROUP_MAP.get(prefix, "UNKNOWN")
    print(f"\n[{sid}]  condition={cond}")

    try:
        df = pd.read_csv(path, dtype=dtype_map, low_memory=False)
    except Exception as e:
        print(f"    ERROR reading file: {e}. Skipping and continuing.")
        continue

    n = len(df)
    scan = (str(df[IMAGE_COL].iloc[0]) if IMAGE_COL in df.columns and n
            else f"scan_{prefix}")
    df["_pheno"] = (df[PHENOTYPE_COL].astype(str).map(ascii_safe)
                    if PHENOTYPE_COL in df.columns else "Unknown")

    x = pd.to_numeric(df[XCOL], errors="coerce").to_numpy()
    y = pd.to_numeric(df[YCOL], errors="coerce").to_numpy()
    ok = np.isfinite(x) & np.isfinite(y)
    if (~ok).sum():
        print(f"    WARNING: {(~ok).sum():,} cells with bad coordinates, excluded "
              f"from spatial views")

    # ---- geometry ----------------------------------------------------------
    y_lo, y_hi = float(np.nanmin(y)), float(np.nanmax(y))
    x_lo, x_hi = float(np.nanmin(x)), float(np.nanmax(x))
    print(f"    n={n:,}  scan={scan}")
    print(f"    X {x_lo:,.0f} to {x_hi:,.0f}   Y {y_lo:,.0f} to {y_hi:,.0f}")

    # ---- downsampled points for plotting -----------------------------------
    idx = np.flatnonzero(ok)
    if len(idx) > MAX_POINTS_PER_SECTION:
        idx_plot = rng.choice(idx, size=MAX_POINTS_PER_SECTION, replace=False)
    else:
        idx_plot = idx
    sections[sid] = {
        "condition": cond,
        "scan": scan,
        "n_cells": n,
        "x_lo": x_lo, "x_hi": x_hi, "y_lo": y_lo, "y_hi": y_hi,
        "x_plot": x[idx_plot], "y_plot": y[idx_plot],
        "pheno_plot": df["_pheno"].to_numpy()[idx_plot],
    }

    # ---- per-section median of every marker --------------------------------
    med = df[INTENSITY_COLS].median(axis=0)
    for c, v in med.items():
        section_marker_med.append({"sample_id": sid, "condition": cond,
                                   "scan": scan, "column": c,
                                   "median": float(v)})

    # ---- per-phenotype mean marker profile ---------------------------------
    prof = df.groupby("_pheno")[INTENSITY_COLS].mean()
    for ph, row in prof.iterrows():
        for c, v in row.items():
            pheno_profiles.append({"sample_id": sid, "condition": cond,
                                   "phenotype": ph, "column": c,
                                   "mean": float(v)})

    # ---- within-section Y deciles ------------------------------------------
    if ok.sum() > 0:
        yy = y[ok]
        edges = np.quantile(yy, np.linspace(0, 1, N_Y_BINS_WITHIN_SECTION + 1))
        edges = np.unique(edges)
        if len(edges) > 2:
            bins = np.clip(np.digitize(yy, edges[1:-1]), 0, len(edges) - 2)
            sub_df = df.loc[ok, INTENSITY_COLS]
            sub_df = sub_df.assign(_bin=bins)
            gb = sub_df.groupby("_bin")
            sizes = gb.size()
            meds = gb.median()
            for b in meds.index:
                if sizes.loc[b] < MIN_CELLS_PER_YBIN:
                    continue
                for c in INTENSITY_COLS:
                    ybin_rows.append({
                        "sample_id": sid, "condition": cond,
                        "y_decile": int(b), "n_cells": int(sizes.loc[b]),
                        "column": c, "median": float(meds.loc[b, c]),
                    })

    # ---- IDO1 on CD68-lineage cells ----------------------------------------
    if IDO1_COL in df.columns:
        mask = df["_pheno"].isin(CD68_LINEAGE).to_numpy()
        vals = pd.to_numeric(df.loc[mask, IDO1_COL], errors="coerce").to_numpy()
        vals = vals[np.isfinite(vals)]
        ido_cd68[sid] = vals
        print(f"    CD68-lineage cells with IDO1 measured: {len(vals):,}")
    else:
        print(f"    WARNING: '{IDO1_COL}' not found")

    del df
    gc.collect()

if not sections:
    print("\nERROR: no sections loaded.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

# ---- ordering ---------------------------------------------------------------
cond_rank = {c: i for i, c in enumerate(CONDITION_ORDER)}
SAMPLE_ORDER = sorted(
    sections.keys(),
    key=lambda s: (cond_rank.get(sections[s]["condition"], 99), s)
)

geom = pd.DataFrame([
    {"sample_id": s, "animal_id": short_label(s),
     "condition": sections[s]["condition"], "scan": sections[s]["scan"],
     "n_cells": sections[s]["n_cells"],
     "x_lo": sections[s]["x_lo"], "x_hi": sections[s]["x_hi"],
     "y_lo": sections[s]["y_lo"], "y_hi": sections[s]["y_hi"],
     "y_mid": 0.5 * (sections[s]["y_lo"] + sections[s]["y_hi"]),
     "y_span": sections[s]["y_hi"] - sections[s]["y_lo"]}
    for s in SAMPLE_ORDER
])


# %% Cell 5 - Q1: physical slide layout
# =============================================================================

banner("Q1 - PHYSICAL SLIDE LAYOUT")

scans = list(dict.fromkeys(sections[s]["scan"] for s in SAMPLE_ORDER))
print(f"    {len(scans)} distinct scan(s):")
for sc in scans:
    members = [s for s in SAMPLE_ORDER if sections[s]["scan"] == sc]
    print(f"        {sc}")
    for s in sorted(members, key=lambda z: sections[z]["y_lo"]):
        print(f"            {s:<12} Y {sections[s]['y_lo']:>10,.0f} to "
              f"{sections[s]['y_hi']:>10,.0f}")

# overlap check within a scan
sub("Y-range overlap check within each scan")
overlap_found = False
for sc in scans:
    members = sorted([s for s in SAMPLE_ORDER if sections[s]["scan"] == sc],
                     key=lambda z: sections[z]["y_lo"])
    for a, b in zip(members[:-1], members[1:]):
        if sections[a]["y_hi"] > sections[b]["y_lo"]:
            print(f"    WARNING: {a} and {b} overlap in Y on {sc}")
            overlap_found = True
if not overlap_found:
    print("    none - sections tile without overlap, consistent with one scan "
          "containing four physically separate tissue sections")

# ---- F11 layout maps --------------------------------------------------------
n_scan = len(scans)
fig, axes = plt.subplots(1, n_scan, figsize=(13 * n_scan, 20), squeeze=False)
axes = axes.ravel()
for k, sc in enumerate(scans):
    ax = axes[k]
    members = [s for s in SAMPLE_ORDER if sections[s]["scan"] == sc]
    cond = sections[members[0]]["condition"] if members else "?"
    base = CONDITION_COLORS.get(cond, "#999999")
    shades = ["#08519C", "#3182BD", "#6BAED6", "#BDD7E7"] if cond == "D1MT" \
        else ["#A63603", "#E6550D", "#FD8D3C", "#FDBE85"]
    for j, s in enumerate(sorted(members, key=lambda z: sections[z]["y_lo"])):
        d = sections[s]
        ax.scatter(d["x_plot"], d["y_plot"], s=1.2,
                   color=shades[j % len(shades)], linewidths=0,
                   rasterized=True, label=short_label(s))
        ax.text(d["x_hi"] + 400, 0.5 * (d["y_lo"] + d["y_hi"]),
                f"{short_label(s)}\n{d['n_cells']:,}",
                fontsize=FONT_SIZE_ANNOT - 6, va="center")
    ax.set_title(f"{cond}\n{sc.split(' - ')[0]}", fontsize=FONT_SIZE_TITLE - 8)
    ax.set_xlabel("X (µm)")
    ax.set_ylabel("Y (µm)")
    ax.set_aspect("equal")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
fig.suptitle("Physical section layout on each scan\n"
             "Four tissue sections per slide, tiled in Y", y=1.01)
save_fig(fig, "F11_slide_layout")

write_csv(geom, "11_section_geometry.csv")


# %% Cell 6 - Q3: slide-position artifact test
# =============================================================================

banner("Q3 - SLIDE POSITION ARTIFACT TEST")

smed = pd.DataFrame(section_marker_med)
piv = smed.pivot_table(index="column", columns="sample_id", values="median")

# z-score each marker WITHIN its own slide, so the comparison is between the four
# sections that were stained together
z_frames = []
for sc in scans:
    members = [s for s in SAMPLE_ORDER if sections[s]["scan"] == sc]
    block = piv[members].replace(0, np.nan)
    lb = np.log2(block)
    z = lb.sub(lb.mean(axis=1), axis=0).div(lb.std(axis=1, ddof=1).replace(0, np.nan),
                                            axis=0)
    z_frames.append(z)
zmat = pd.concat(z_frames, axis=1)[SAMPLE_ORDER]

print("    Per-section median z-score across all markers")
print("    (negative = this section is globally dimmer than its slide-mates)\n")
summary_rows = []
for s in SAMPLE_ORDER:
    col = zmat[s].dropna()
    frac_low = float((col < -0.8).mean())
    summary_rows.append({
        "sample_id": s, "animal_id": short_label(s),
        "condition": sections[s]["condition"],
        "y_lo": sections[s]["y_lo"],
        "slide_position_rank": np.nan,
        "median_z_all_markers": float(col.median()),
        "frac_markers_low": frac_low,
        "n_markers": int(len(col)),
    })
pos_tbl = pd.DataFrame(summary_rows)
for sc in scans:
    members = [s for s in SAMPLE_ORDER if sections[s]["scan"] == sc]
    order = sorted(members, key=lambda z: sections[z]["y_lo"])
    for rank, s in enumerate(order, start=1):
        pos_tbl.loc[pos_tbl["sample_id"] == s, "slide_position_rank"] = rank

pos_tbl = pos_tbl.sort_values(["condition", "slide_position_rank"])
for _, r in pos_tbl.iterrows():
    flag = "   <-- globally dim" if r["median_z_all_markers"] < -0.5 else ""
    print(f"    {r['sample_id']:<12} pos={int(r['slide_position_rank'])}  "
          f"median z={r['median_z_all_markers']:+.2f}  "
          f"markers below -0.8: {100*r['frac_markers_low']:5.1f}%{flag}")

print("\n    INTERPRETATION GUIDE")
print("    If 43102 and 43112 are dim across most of the 67 markers, the zero")
print("    IDO1+ call in those two sections is a slide position artifact.")
print("    If they are dim only for IDO1 and its correlates, it is more likely")
print("    biology. Read the IDO1 row of F12 against the all-marker distribution.")

# ---- F12: two panels --------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(30, 16),
                         gridspec_kw={"width_ratios": [1.15, 1.0]})

# left: heatmap of z, key markers
ax = axes[0]
key_cols = [c for c in zmat.index if marker_base(c) in KEY_MARKERS]
key_cols = sorted(key_cols, key=lambda c: KEY_MARKERS.index(marker_base(c)))
sub_z = zmat.loc[key_cols, SAMPLE_ORDER]
vmax = float(np.nanmax(np.abs(sub_z.to_numpy()))) or 1.0
cmap = LinearSegmentedColormap.from_list("z", ["#2166AC", "#F7F7F7", "#B2182B"])
im = ax.imshow(sub_z.to_numpy(), aspect="auto", cmap=cmap,
               norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax))
ax.set_xticks(np.arange(len(SAMPLE_ORDER)))
ax.set_xticklabels([f"{short_label(s)}\npos {int(pos_tbl.loc[pos_tbl['sample_id']==s,'slide_position_rank'].iloc[0])}"
                    for s in SAMPLE_ORDER], fontsize=FONT_SIZE_TICK - 10)
ax.set_yticks(np.arange(len(key_cols)))
ax.set_yticklabels([marker_base(c) for c in key_cols], fontsize=FONT_SIZE_TICK - 10)
n_first = sum(1 for s in SAMPLE_ORDER if sections[s]["condition"] == CONDITION_ORDER[0])
ax.axvline(n_first - 0.5, color="#000000", linewidth=4)
cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
cb.set_label("z within slide", fontsize=FONT_SIZE_BASE - 8)
cb.ax.tick_params(labelsize=FONT_SIZE_TICK - 12)
ax.set_title("Marker medians, z-scored within slide", fontsize=FONT_SIZE_TITLE - 8)

# right: distribution of z across ALL markers per section
ax = axes[1]
data = [zmat[s].dropna().to_numpy() for s in SAMPLE_ORDER]
bp = ax.boxplot(data, vert=True, patch_artist=True, widths=0.6,
                medianprops=dict(color="#000000", linewidth=3.5),
                flierprops=dict(marker=".", markersize=6,
                                markerfacecolor="#999999",
                                markeredgecolor="none"))
for patch, s in zip(bp["boxes"], SAMPLE_ORDER):
    patch.set_facecolor(CONDITION_COLORS.get(sections[s]["condition"], "#999999"))
    patch.set_alpha(0.75)
    patch.set_edgecolor("#333333")
ax.axhline(0, color="#000000", linestyle="--", linewidth=3)
ax.set_xticks(np.arange(1, len(SAMPLE_ORDER) + 1))
ax.set_xticklabels([short_label(s) for s in SAMPLE_ORDER],
                   rotation=45, ha="right", fontsize=FONT_SIZE_TICK - 6)
ax.set_ylabel("z of log$_2$ median, all 67 markers")
ax.set_title("Global brightness per section", fontsize=FONT_SIZE_TITLE - 8)
style_axes(ax)

fig.suptitle("Is the bottom section of each slide globally dim?\n"
             "A section shifted down across most markers indicates a position "
             "artifact, not biology", y=1.03, fontsize=FONT_SIZE_TITLE - 4)
save_fig(fig, "F12_slide_position_artifact")

write_csv(pos_tbl, "12_slide_position_summary.csv")
write_csv(zmat.reset_index().rename(columns={"index": "column"}),
          "12b_marker_z_within_slide.csv")

# ---- within-section Y gradient for IDO1 -------------------------------------
ybin = pd.DataFrame(ybin_rows)
if len(ybin) and IDO1_COL in set(ybin["column"]):
    fig, ax = plt.subplots(figsize=(18, 12))
    for s in SAMPLE_ORDER:
        d = ybin.loc[(ybin["sample_id"] == s) & (ybin["column"] == IDO1_COL)]
        if not len(d):
            continue
        d = d.sort_values("y_decile")
        ax.plot(d["y_decile"], d["median"], marker="o", markersize=14,
                linewidth=4, label=short_label(s),
                color=CONDITION_COLORS.get(sections[s]["condition"], "#999999"),
                alpha=0.85)
    ax.set_xlabel("Within-section Y decile (0 = bottom edge)")
    ax.set_ylabel("Median IDO1 (cytoplasm)")
    ax.set_title("IDO1 signal across each section's own Y axis",
                 fontsize=FONT_SIZE_TITLE - 6)
    style_axes(ax)
    ax.legend(bbox_to_anchor=(1.01, 1.0), loc="upper left", frameon=False,
              fontsize=FONT_SIZE_LEGEND - 8, ncol=1)
    save_fig(fig, "F12b_ido1_within_section_gradient")

if len(ybin):
    write_csv(ybin, "12c_within_section_y_deciles.csv")


# %% Cell 7 - Q4: IDO1 distribution on CD68+ macrophages
# =============================================================================

banner("Q4 - IDO1 DISTRIBUTION ON CD68-LINEAGE MACROPHAGES")

if not ido_cd68:
    print("    WARNING: no IDO1 values collected. Skipping Q4.")
    ido_sum = pd.DataFrame()
else:
    rows = []
    for s in SAMPLE_ORDER:
        v = ido_cd68.get(s, np.array([]))
        if len(v) == 0:
            continue
        rows.append({
            "sample_id": s, "animal_id": short_label(s),
            "condition": sections[s]["condition"],
            "n_cd68_lineage": int(len(v)),
            "min": float(np.min(v)), "p25": float(np.percentile(v, 25)),
            "median": float(np.median(v)), "p75": float(np.percentile(v, 75)),
            "p95": float(np.percentile(v, 95)), "p99": float(np.percentile(v, 99)),
            "max": float(np.max(v)), "frac_zero": float(np.mean(v == 0)),
        })
    ido_sum = pd.DataFrame(rows)
    with pd.option_context("display.width", 220, "display.max_columns", 30):
        print(ido_sum.to_string(index=False))

    fig, axes = plt.subplots(1, 2, figsize=(30, 13))

    # left: ECDF overlay
    ax = axes[0]
    for s in SAMPLE_ORDER:
        v = ido_cd68.get(s, np.array([]))
        if len(v) == 0:
            continue
        vs = np.sort(v)
        ax.plot(np.log1p(vs), np.arange(1, len(vs) + 1) / len(vs),
                linewidth=4.5, alpha=0.85, label=short_label(s),
                color=CONDITION_COLORS.get(sections[s]["condition"], "#999999"),
                linestyle="-" if sections[s]["condition"] == "D1MT" else "--")
    ax.set_xlabel("log(1 + IDO1 cytoplasm)")
    ax.set_ylabel("Cumulative fraction of CD68+ cells")
    ax.set_title("IDO1 on CD68-lineage macrophages", fontsize=FONT_SIZE_TITLE - 8)
    ax.set_ylim(0, 1.02)
    style_axes(ax)
    ax.legend(loc="lower right", frameon=False, fontsize=FONT_SIZE_LEGEND - 10,
              ncol=2)

    # right: per-section distribution, offset ridgeline
    ax = axes[1]
    all_v = np.concatenate([v for v in ido_cd68.values() if len(v)])
    hi = np.percentile(all_v, 99.5)
    bins = np.linspace(0, np.log1p(hi), 70)
    for i, s in enumerate(SAMPLE_ORDER):
        v = ido_cd68.get(s, np.array([]))
        if len(v) == 0:
            continue
        h, e = np.histogram(np.log1p(v), bins=bins, density=True)
        h = h / (h.max() if h.max() > 0 else 1.0)
        centers = 0.5 * (e[:-1] + e[1:])
        base = len(SAMPLE_ORDER) - i
        col = CONDITION_COLORS.get(sections[s]["condition"], "#999999")
        ax.fill_between(centers, base, base + h * 0.9, color=col, alpha=0.7,
                        linewidth=0)
        ax.plot(centers, base + h * 0.9, color="#333333", linewidth=2)
        ax.text(bins[-1] * 1.01, base + 0.35,
                f"{short_label(s)}  n={len(v):,}",
                fontsize=FONT_SIZE_ANNOT - 8, va="center")
    ax.set_yticks([])
    ax.set_xlabel("log(1 + IDO1 cytoplasm)")
    ax.set_title("Per-section distribution", fontsize=FONT_SIZE_TITLE - 8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    fig.suptitle("Does the IDO1-high mode disappear in treated animals, "
                 "or does the whole distribution shift?", y=1.03,
                 fontsize=FONT_SIZE_TITLE - 4)
    save_fig(fig, "F13_ido1_distribution_cd68")

    write_csv(ido_sum, "13_ido1_cd68_summary.csv")


# %% Cell 8 - Q5: what is in "Other"
# =============================================================================

banner('Q5 - MARKER PROFILE OF THE "OTHER" BUCKET')

prof = pd.DataFrame(pheno_profiles)
other_tbl = pd.DataFrame()
if not len(prof):
    print("    WARNING: no phenotype profiles collected.")
else:
    # average across sections, then z-score each marker across phenotypes
    mean_prof = prof.groupby(["phenotype", "column"])["mean"].mean().unstack()
    mean_prof = mean_prof.replace(0, np.nan)
    lp = np.log2(mean_prof)
    z = lp.sub(lp.mean(axis=0), axis=1).div(lp.std(axis=0, ddof=1).replace(0, np.nan),
                                            axis=1)

    if "Other" in z.index:
        oz = z.loc["Other"].dropna().sort_values(ascending=False)
        sub('Markers most enriched in "Other" relative to the other phenotypes')
        for c, v in oz.head(15).items():
            print(f"    {c:<48} z={v:+.2f}")
        sub('Markers most depleted in "Other"')
        for c, v in oz.tail(10).items():
            print(f"    {c:<48} z={v:+.2f}")
        other_tbl = oz.reset_index()
        other_tbl.columns = ["column", "z_vs_other_phenotypes"]
        write_csv(other_tbl, "14_other_marker_profile.csv")
    else:
        print("    WARNING: no 'Other' phenotype found.")

    # ---- F14 heatmap: phenotype x key marker ------------------------------
    ph_use = [p for p in PHENOTYPE_ORDER if p in z.index]
    key_cols = [c for c in z.columns if marker_base(c) in KEY_MARKERS]
    key_cols = sorted(key_cols, key=lambda c: KEY_MARKERS.index(marker_base(c)))
    sub_z = z.loc[ph_use, key_cols]
    vmax = float(np.nanmax(np.abs(sub_z.to_numpy()))) or 1.0

    fig, ax = plt.subplots(figsize=(24, 14))
    im = ax.imshow(sub_z.to_numpy(), aspect="auto", cmap=cmap,
                   norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax))
    ax.set_xticks(np.arange(len(key_cols)))
    ax.set_xticklabels([marker_base(c) for c in key_cols], rotation=90,
                       fontsize=FONT_SIZE_TICK - 10)
    ax.set_yticks(np.arange(len(ph_use)))
    ax.set_yticklabels(ph_use, fontsize=FONT_SIZE_TICK - 8)
    if "Other" in ph_use:
        ax.axhline(ph_use.index("Other") - 0.5, color="#000000", linewidth=4)
        ax.axhline(ph_use.index("Other") + 0.5, color="#000000", linewidth=4)
    cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cb.set_label("z across phenotypes", fontsize=FONT_SIZE_BASE - 8)
    cb.ax.tick_params(labelsize=FONT_SIZE_TICK - 12)
    ax.set_title("Mean marker profile per phenotype\n"
                 'The boxed row is "Other", 16 to 36 percent of every section',
                 fontsize=FONT_SIZE_TITLE - 6)
    save_fig(fig, "F14_phenotype_marker_profiles")

    write_csv(prof, "14b_phenotype_marker_profiles_long.csv")


# %% Cell 9 - spatial maps of key immune populations
# =============================================================================

banner("SPATIAL MAPS OF KEY IMMUNE POPULATIONS")

fig, axes = plt.subplots(2, 4, figsize=(46, 24))
axes = axes.ravel()
for k, s in enumerate(SAMPLE_ORDER):
    ax = axes[k]
    d = sections[s]
    ph = d["pheno_plot"]
    bg = ~np.isin(ph, SPATIAL_HIGHLIGHT)
    ax.scatter(d["x_plot"][bg], d["y_plot"][bg], s=0.6, color=BACKGROUND_COLOR,
               linewidths=0, rasterized=True)
    for p in SPATIAL_HIGHLIGHT:
        m = ph == p
        if m.sum() == 0:
            continue
        ax.scatter(d["x_plot"][m], d["y_plot"][m], s=6,
                   color=PHENOTYPE_COLORS.get(p, "#000000"),
                   linewidths=0, rasterized=True)
    ax.set_title(f"{short_label(s)}  ({d['condition']})\n{d['n_cells']:,} cells",
                 fontsize=FONT_SIZE_TITLE - 10)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_edgecolor(CONDITION_COLORS.get(d["condition"], "#999999"))
        sp.set_linewidth(6)

handles = [Patch(facecolor=PHENOTYPE_COLORS.get(p, "#000"), label=p)
           for p in SPATIAL_HIGHLIGHT]
handles.append(Patch(facecolor=BACKGROUND_COLOR, label="all other cells"))
fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False,
           fontsize=FONT_SIZE_LEGEND - 6, bbox_to_anchor=(0.5, -0.03))
fig.suptitle("Spatial distribution of macrophage and lymphocyte populations\n"
             f"(plotting downsampled to at most {MAX_POINTS_PER_SECTION:,} cells "
             "per section)", y=1.02, fontsize=FONT_SIZE_TITLE - 2)
save_fig(fig, "F15_spatial_maps_key_populations")


# %% Cell 10 - Q6: nearest-neighbour pair feasibility
# =============================================================================

banner("Q6 - NEAREST-NEIGHBOUR PAIR FEASIBILITY")

cnt_path = os.path.join(INVENTORY_TABLE_DIR, "04_phenotype_counts_long.csv")
pair_tbl = pd.DataFrame()
if not os.path.exists(cnt_path):
    print(f"    WARNING: {cnt_path} not found. Skipping Q6.")
else:
    pl = pd.read_csv(cnt_path)
    pl["phenotype_display"] = pl["phenotype"].map(ascii_safe)
    cmat = (pl.pivot_table(index="phenotype_display", columns="sample_id",
                           values="n_cells", aggfunc="sum")
            .reindex(columns=SAMPLE_ORDER).fillna(0.0))

    rows = []
    for a in NN_ANCHORS:
        for t in NN_TARGETS:
            if a == t or a not in cmat.index or t not in cmat.index:
                continue
            rec = {"anchor": a, "target": t}
            for cond in CONDITION_ORDER:
                members = [s for s in SAMPLE_ORDER
                           if sections[s]["condition"] == cond]
                ok_n = sum(1 for s in members
                           if cmat.loc[a, s] >= MIN_CELLS_FOR_NN
                           and cmat.loc[t, s] >= MIN_CELLS_FOR_NN)
                rec[f"n_animals_ok_{cond}"] = ok_n
                rec[f"min_anchor_{cond}"] = int(cmat.loc[a, members].min())
                rec[f"min_target_{cond}"] = int(cmat.loc[t, members].min())
            rec["usable_both_arms"] = all(
                rec[f"n_animals_ok_{c}"] == sum(
                    1 for s in SAMPLE_ORDER if sections[s]["condition"] == c)
                for c in CONDITION_ORDER)
            rows.append(rec)
    pair_tbl = pd.DataFrame(rows)

    print(f"    Anchor/target pairs evaluated: {len(pair_tbl)}")
    print(f"    Threshold: >= {MIN_CELLS_FOR_NN} cells per animal for BOTH members\n")
    print("    n animals passing (out of 4 per arm):\n")
    hdr_line = f"    {'anchor':<26}{'target':<26}{'D1MT':>6}{'Untreated':>12}"
    print(hdr_line); print("    " + "-" * (len(hdr_line) - 4))
    for _, r in pair_tbl.iterrows():
        mark = "  <-- both arms complete" if r["usable_both_arms"] else ""
        print(f"    {r['anchor']:<26}{r['target']:<26}"
              f"{r['n_animals_ok_D1MT']:>6}{r['n_animals_ok_Untreated']:>12}{mark}")

    # ---- F16 heatmap -------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(30, 13))
    n_per_arm = {c: sum(1 for s in SAMPLE_ORDER
                        if sections[s]["condition"] == c) for c in CONDITION_ORDER}
    greens = LinearSegmentedColormap.from_list(
        "feas", ["#FFFFFF", "#D9F0D3", "#7FBC41", "#1B7837"])
    for k, cond in enumerate(CONDITION_ORDER):
        ax = axes[k]
        m = pair_tbl.pivot_table(index="anchor", columns="target",
                                 values=f"n_animals_ok_{cond}")
        m = m.reindex(index=[a for a in NN_ANCHORS if a in m.index],
                      columns=[t for t in NN_TARGETS if t in m.columns])
        im = ax.imshow(m.to_numpy(dtype=float), cmap=greens, vmin=0,
                       vmax=n_per_arm[cond], aspect="auto")
        ax.set_xticks(np.arange(m.shape[1]))
        ax.set_xticklabels(m.columns, rotation=45, ha="right",
                           fontsize=FONT_SIZE_TICK - 12)
        ax.set_yticks(np.arange(m.shape[0]))
        ax.set_yticklabels(m.index, fontsize=FONT_SIZE_TICK - 12)
        for i in range(m.shape[0]):
            for j in range(m.shape[1]):
                v = m.iloc[i, j]
                if pd.isna(v):
                    continue
                ax.text(j, i, f"{int(v)}", ha="center", va="center",
                        fontsize=FONT_SIZE_ANNOT - 4,
                        color="#FFFFFF" if v > n_per_arm[cond] * 0.6 else "#000000",
                        fontweight="bold" if v == n_per_arm[cond] else "normal")
        ax.set_title(f"{cond}  (n = {n_per_arm[cond]} animals)",
                     fontsize=FONT_SIZE_TITLE - 8)
    fig.suptitle("Animals with enough cells for a nearest-neighbour test\n"
                 f"(anchor and target both >= {MIN_CELLS_FOR_NN} cells; "
                 "rows = anchor, columns = target)", y=1.05,
                 fontsize=FONT_SIZE_TITLE - 4)
    save_fig(fig, "F16_nn_pair_feasibility")

    write_csv(pair_tbl, "15_nn_pair_feasibility.csv")


# %% Cell 11 - wrap up
# =============================================================================

banner("SUMMARY")

print(f"Sections loaded : {len(sections)}")
print(f"Distinct scans  : {len(scans)}")
print(f"Total cells     : {sum(sections[s]['n_cells'] for s in SAMPLE_ORDER):,}")

sub("Questions and where to read the answer")
print("  Q1 section tiling        -> F11_slide_layout, 11_section_geometry.csv")
print("  Q2 batch confound        -> F10_variance_between_vs_within_slide,")
print("                              10_variance_decomposition.csv")
print("  Q3 position artifact     -> F12_slide_position_artifact, F12b, 12_*.csv")
print("  Q4 IDO1 distribution     -> F13_ido1_distribution_cd68, 13_*.csv")
print('  Q5 what is in "Other"    -> F14_phenotype_marker_profiles, 14_*.csv')
print("  spatial overview         -> F15_spatial_maps_key_populations")
print("  Q6 NN feasibility        -> F16_nn_pair_feasibility, 15_*.csv")

sub("Decisions still open after this run")
print("  1. Whether 43102 and 43112 stay in the analysis, and on what grounds")
print("  2. Whether IDO1 positivity is used as the supplied binary call or as a")
print("     within-slide continuous measure")
print("  3. Which anchor/target pairs go into the nearest-neighbour work")
print('  4. Whether "Other" is excluded, split, or carried as an unknown class')

banner("DONE")
sys.stdout = _tee.terminal
_tee.close()
