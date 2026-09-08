#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - SPATIAL LAYOUT AND BATCH DIAGNOSTICS
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 02 of the AKOYA analysis series. REVISION 2.

PURPOSE
    Answer the structural questions raised by the script 01 inventory, before any
    gating, normalization or spatial statistics are committed to.

    Q1  Do the four sections per scan tile cleanly in a shared coordinate space?
    Q2  How much of the marker intensity spread is between scan versus between
        section within a scan? Scan is perfectly confounded with arm here, so
        this bounds the risk rather than resolving it.
    Q3  Are 43102 and 43112 (position 1 on each scan) globally low across many
        markers, which would indicate a slide position artifact rather than
        biology?
    Q4  What does the IDO1 intensity distribution on CD68+ macrophages actually
        look like per section? Is there a distinct IDO1-high mode, and does it
        vanish in the treated arm or just shift?
    Q5  What is sitting in the "Other" bucket?
    Q6  For every anchor / target phenotype pair, how many animals have enough
        cells to support a nearest-neighbour test?

    Still diagnostic only. Nothing is re-phenotyped, gated, normalized or written
    back to the source CSVs.

WHAT CHANGED IN REVISION 2 (and why)

    1. Q2 IS KEYED ON scan_id, NOT ON treatment.
       The previous version built its slide key as
       dict(zip(fs["sample_id"], fs["treatment"])). Because each scan carries
       exactly one arm, that produced the right grouping by accident, but it
       meant the reported quantity could never be described accurately. Script
       01 revision 2 writes a real scan_id and this script now reads it.

    2. Q2 REPORTS A CORRECTED VARIANCE COMPONENT, NOT ONLY THE NAIVE RATIO.
       frac_between_slide was between / (between + within), where "between" was
       the raw variance of the two scan means. The variance of group means
       already contains within-group variance divided by n, so that ratio is
       biased upward. The corrected between-scan component is (MSB - MSW) / n
       and the corresponding ICC is now reported alongside. The naive column is
       retained under its old name so earlier numbers can be traced.

    3. THE NAME IS HONEST NOW.
       Scan and arm are the same factor in this cohort. A high between-scan
       fraction is consistent with a staining batch difference AND with a real
       biological difference between a granuloma-rich lung and a treated lung,
       and nothing in this table can separate them. CD3e, CD163, CD4 and IFNG
       are exactly the markers you would expect to differ biologically, so
       reading their high fractions as "batch" would be wrong. Every printed
       line and the figure title now say this.

    4. Q3 GAINS A COMPOSITION-CONTROLLED VERSION AND A RANK STATISTIC.
       The old dimness test z-scored each marker's SECTION-WIDE median across
       the four sections of a scan. Section-wide medians depend on which cells
       are in the section, and composition varies enormously here: 43102 is 89
       percent endothelial, epithelial and Other, while 43109 is 36 percent
       Other and 11 percent neutrophils. A section can therefore look dim
       because of what is in it rather than how it stained.

       The test is now run two ways. The original all-cell version is kept for
       continuity. A composition-controlled version z-scores the median of each
       marker WITHIN each phenotype, using only phenotypes that clear a cell
       count in every section of that scan, so every comparison is like for
       like.

       A rank statistic is added to both. For each marker (and phenotype), the
       four sections of a scan are ranked and the fraction of comparisons where
       a section is the dimmest of its four is reported. Under no position
       effect that fraction is 0.25 with no distributional assumptions. This
       matters because with four sections per scan the minimum possible z is
       -1.5, so a median z of -1.04 is closer to the floor than it looks and
       needs a reference point.

    5. PALETTES AND THE SCAN BOUNDARY ARE DEFINED ONCE.
       cmap was created in the Q3 cell and reused by the Q5 figure, which raised
       NameError when cells were run out of order in Spyder. The condition-based
       n_first boundary is replaced by a scan-derived boundary, same fix as 01b.

    6. SMALL FIXES.
       vmax = float(...) or 1.0 returns NaN when the input is all NaN, because
       NaN is truthy. Replaced with an explicit check in all three places.
       Section geometry now carries scan_id and slide_position_rank from script
       01 rather than recomputing the rank here.

OUTPUTS
    figures/  F10 .. F16   (PDF + PNG, 300 DPI)
    tables/   10 .. 15     (CSV, LF line endings)

USAGE
    conda activate sc_pre
    python AKOYA_02_Spatial_Diagnostic.py

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

# Scan / slide identity comes from script 01 revision 2. The Image column is
# used as a fallback if those columns are absent.
GROUP_MAP = {"G3": "D1MT", "G4": "Untreated"}
SCAN_COL = "scan_id"
POSITION_COL = "slide_position_rank"

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
BRIGHTNESS_Z_LOW = -0.8        # a marker is "low" in a section below this z
BRIGHTNESS_MEDIAN_Z_FLAG = -0.5  # a section is flagged globally dim below this
# Composition-controlled version: a phenotype contributes only if it clears this
# count in EVERY section of the scan, so the comparison is like for like.
COMPOSITION_MIN_CELLS = 200
COMPOSITION_MIN_PHENOTYPES = 3   # warn if a scan has fewer qualifying phenotypes

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
SCAN_BOUNDARY_COLOR = "#000000"
SCAN_BOUNDARY_WIDTH = 4.0

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

# Palettes defined once so no cell depends on another cell having run.
CMAP_DIVERGING = LinearSegmentedColormap.from_list(
    "z", ["#2166AC", "#F7F7F7", "#B2182B"])
CMAP_FEASIBILITY = LinearSegmentedColormap.from_list(
    "feas", ["#FFFFFF", "#D9F0D3", "#7FBC41", "#1B7837"])

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


def safe_vmax(arr, floor=1.0):
    """np.nanmax(...) or 1.0 returns NaN when everything is NaN, since NaN is
    truthy. This does what that line was meant to do."""
    a = np.asarray(arr, dtype=float)
    if not np.isfinite(a).any():
        return floor
    v = float(np.nanmax(np.abs(a)))
    return v if np.isfinite(v) and v > 0 else floor


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


def zscore_within(block):
    """
    block: markers (rows) x sections (columns), positive intensities.
    Returns log2 values z-scored across the columns of that block.
    """
    b = block.replace(0, np.nan)
    lb = np.log2(b)
    sd = lb.std(axis=1, ddof=1).replace(0, np.nan)
    return lb.sub(lb.mean(axis=1), axis=0).div(sd, axis=0)


def rank_within(block):
    """
    block: markers (rows) x sections (columns).
    Rank 1 = dimmest of the block for that marker. NaN rows stay NaN.
    """
    return block.replace(0, np.nan).rank(axis=1, method="average",
                                         na_option="keep")


_tee = Tee(os.path.join(TAB_DIR, "00_diagnostics_report.txt"))
sys.stdout = _tee

banner("AKOYA SPATIAL LAYOUT AND BATCH DIAGNOSTICS (revision 2)")
print(f"Run time  : {datetime.now().isoformat(timespec='seconds')}")
print(f"Data      : {DATA_DIR}")
print(f"Inventory : {INVENTORY_TABLE_DIR}")
print(f"Figures   : {FIG_DIR}")
print(f"Tables    : {TAB_DIR}")


# %% Cell 3 - inventory metadata: scan identity and section position
# =============================================================================
# Read once, up front, so every later cell keys on the same scan definition.

banner("INVENTORY METADATA")

fs_path = os.path.join(INVENTORY_TABLE_DIR, "01_file_summary.csv")
ms_path = os.path.join(INVENTORY_TABLE_DIR, "06_marker_stats_long.csv")
cnt_path = os.path.join(INVENTORY_TABLE_DIR, "04_phenotype_counts_long.csv")

INV = None
SCAN_FROM_INVENTORY = {}
POS_FROM_INVENTORY = {}
if os.path.exists(fs_path):
    INV = pd.read_csv(fs_path)
    print(f"    loaded 01_file_summary.csv  ({INV.shape[0]} x {INV.shape[1]})")
    if SCAN_COL in INV.columns:
        SCAN_FROM_INVENTORY = dict(zip(INV["sample_id"], INV[SCAN_COL].astype(str)))
        print(f"    scan key   : '{SCAN_COL}' from script 01")
    else:
        print(f"    WARNING: '{SCAN_COL}' absent. Run script 01 revision 2. "
              f"Falling back to the Image column parsed from the raw CSVs, and "
              f"Q2 will fall back to treatment, which cannot be described "
              f"accurately.")
    if POSITION_COL in INV.columns:
        POS_FROM_INVENTORY = {r["sample_id"]: (int(r[POSITION_COL])
                                               if pd.notna(r[POSITION_COL]) else None)
                              for _, r in INV.iterrows()}
        print(f"    position   : '{POSITION_COL}' from script 01")
    else:
        print(f"    WARNING: '{POSITION_COL}' absent. Position rank will be "
              f"recomputed here from Y extents.")
else:
    print(f"    WARNING: {fs_path} not found. Run script 01 first.")


# %% Cell 4 - Q2: between-scan vs within-scan intensity spread
# =============================================================================
# Uses the p99 values already computed in script 01. No raw reads needed.

banner("Q2 - BETWEEN-SCAN VS WITHIN-SCAN INTENSITY SPREAD")

var_tbl = pd.DataFrame()
if not (os.path.exists(ms_path) and INV is not None):
    print(f"    WARNING: inventory tables not found. Skipping Q2.")
else:
    ms = pd.read_csv(ms_path)

    # Prefer the scan_id already carried in the marker table, then the file
    # summary, then treatment as a last resort with an explicit warning.
    if SCAN_COL in ms.columns:
        ms["scan"] = ms[SCAN_COL].astype(str)
        scan_source = f"'{SCAN_COL}' column of 06_marker_stats_long.csv"
    elif SCAN_FROM_INVENTORY:
        ms["scan"] = ms["sample_id"].map(SCAN_FROM_INVENTORY)
        scan_source = f"'{SCAN_COL}' column of 01_file_summary.csv"
    else:
        ms["scan"] = ms["sample_id"].map(dict(zip(INV["sample_id"],
                                                  INV["treatment"])))
        scan_source = "treatment (FALLBACK, not a real scan key)"
    print(f"    grouping variable: {scan_source}")

    arm_of_sample = dict(zip(INV["sample_id"], INV["treatment"]))
    ms["arm"] = ms["sample_id"].map(arm_of_sample)
    arm_of_scan = (ms.drop_duplicates("sample_id")
                   .groupby("scan")["arm"].agg(lambda v: sorted(set(v))))
    single_arm_scans = all(len(v) == 1 for v in arm_of_scan)
    print(f"    scans: {list(arm_of_scan.index)}")
    for sc, arms in arm_of_scan.items():
        n_sec = ms.loc[ms["scan"] == sc, "sample_id"].nunique()
        print(f"      {sc}  {n_sec} section(s)  arm(s): {', '.join(arms)}")

    ms = ms.loc[ms["p99"].notna() & (ms["p99"] > 0)].copy()
    ms["log2_p99"] = np.log2(ms["p99"])

    rows = []
    for col, g in ms.groupby("column"):
        scan_means = g.groupby("scan")["log2_p99"].mean()
        scan_sizes = g.groupby("scan")["log2_p99"].size()
        if len(scan_means) < 2:
            continue
        k = len(scan_means)
        n_bar = float(scan_sizes.mean())

        # within-scan variance, pooled across scans
        msw = float(g.groupby("scan")["log2_p99"].var(ddof=1).mean())
        # naive quantity reported by revision 1, kept for traceability
        raw_between = float(np.var(scan_means.to_numpy(), ddof=1))
        naive_total = msw + raw_between
        frac_naive = raw_between / naive_total if naive_total > 0 else np.nan

        # corrected variance components for a (near) balanced one-way layout:
        # E[MSB] = MSW + n * sigma2_between
        msb = n_bar * raw_between
        sigma2_between = max(0.0, (msb - msw) / n_bar) if n_bar > 0 else np.nan
        icc_total = sigma2_between + msw
        icc = sigma2_between / icc_total if icc_total > 0 else np.nan

        rec = {
            "column": col,
            "marker": marker_base(col),
            "is_key_marker": marker_base(col) in KEY_MARKERS,
            "n_scans": k,
            "mean_sections_per_scan": n_bar,
            "within_scan_var": msw,
            "raw_between_scan_var": raw_between,
            "sigma2_between_scan": sigma2_between,
            "frac_between_scan_naive": frac_naive,
            "icc_between_scan": icc,
            "within_scan_sd": np.sqrt(msw),
        }
        # direction, only meaningful when each scan carries one arm
        if single_arm_scans:
            per_arm = {arm_of_scan[sc][0]: v for sc, v in scan_means.items()}
            rec["log2_scan_diff"] = float(
                per_arm.get(CONDITION_ORDER[1], np.nan)
                - per_arm.get(CONDITION_ORDER[0], np.nan))
        else:
            rec["log2_scan_diff"] = np.nan

        # backwards-compatible aliases for revision 1 column names
        rec["within_slide_var"] = rec["within_scan_var"]
        rec["between_slide_var"] = rec["raw_between_scan_var"]
        rec["frac_between_slide"] = rec["frac_between_scan_naive"]
        rec["within_slide_sd"] = rec["within_scan_sd"]
        rec["log2_slide_diff"] = rec["log2_scan_diff"]
        rows.append(rec)

    var_tbl = pd.DataFrame(rows).sort_values("icc_between_scan", ascending=False)

    print(f"\n    {len(var_tbl)} markers decomposed\n")
    print("    WHAT THIS QUANTITY IS, AND IS NOT")
    print("    Scan is perfectly confounded with treatment arm in this cohort:")
    print("    one acquisition per group. A high between-scan fraction is")
    print("    therefore consistent with a staining batch difference AND with a")
    print("    real biological difference between a granuloma-rich lung and a")
    print("    treated lung. Nothing in this table separates the two. CD3e,")
    print("    CD163, CD4 and IFNG rank high partly because they are the markers")
    print("    that SHOULD differ biologically between these arms, so reading")
    print("    their scores as evidence of batch would be wrong.")
    print("    What the table does support is the opposite direction: a marker")
    print("    with a LOW between-scan component cannot be carrying much of an")
    print("    acquisition difference, so it is safe for cross-arm comparison.\n")
    print("    Two columns are reported. frac_between_scan_naive is")
    print("    between / (between + within) using the raw variance of the scan")
    print("    means, which is what revision 1 reported and which is biased")
    print("    upward because the variance of group means already contains")
    print("    within-group variance over n. icc_between_scan uses the corrected")
    print("    component (MSB - MSW) / n and is the one to quote.")

    sub("Markers most dominated by between-scan variation (top 20 by ICC)")
    print(f"    {'column':<48}{'ICC':>7}{'naive':>8}{'log2 diff':>11}")
    for _, r in var_tbl.head(20).iterrows():
        star = " *" if r["is_key_marker"] else ""
        print(f"    {r['column'][:47]:<48}{r['icc_between_scan']:>7.3f}"
              f"{r['frac_between_scan_naive']:>8.3f}"
              f"{r['log2_scan_diff']:>+11.2f}{star}")

    sub("Most comparable markers (lowest ICC), safest for cross-arm work")
    for _, r in var_tbl.tail(15).iloc[::-1].iterrows():
        star = " *" if r["is_key_marker"] else ""
        print(f"    {r['column'][:47]:<48}{r['icc_between_scan']:>7.3f}"
              f"{r['frac_between_scan_naive']:>8.3f}"
              f"{r['log2_scan_diff']:>+11.2f}{star}")

    sub("Key phenotype-calling markers")
    print(f"    {'column':<48}{'ICC':>7}{'naive':>8}{'log2 diff':>11}"
          f"{'within sd':>11}")
    for _, r in var_tbl.loc[var_tbl["is_key_marker"]].iterrows():
        print(f"    {r['column'][:47]:<48}{r['icc_between_scan']:>7.3f}"
              f"{r['frac_between_scan_naive']:>8.3f}"
              f"{r['log2_scan_diff']:>+11.2f}{r['within_scan_sd']:>11.2f}")

    # ---- figure F10 ---------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(32, 14))

    for ax, ycol, ttl in zip(
            axes,
            ["icc_between_scan", "frac_between_scan_naive"],
            ["Corrected component (ICC)", "Naive ratio (revision 1)"]):
        ax.scatter(var_tbl["within_scan_sd"], var_tbl[ycol],
                   s=200, color="#BBBBBB", edgecolor="#FFFFFF", linewidth=1.5,
                   zorder=3, label="all markers")
        key = var_tbl.loc[var_tbl["is_key_marker"]]
        ax.scatter(key["within_scan_sd"], key[ycol],
                   s=520, color=FLAG_COLOR, edgecolor="#FFFFFF", linewidth=2.5,
                   zorder=4, label="phenotype-calling markers")
        for _, r in key.iterrows():
            ax.annotate(r["marker"], (r["within_scan_sd"], r[ycol]),
                        textcoords="offset points", xytext=(12, 8),
                        fontsize=FONT_SIZE_ANNOT - 8)
        ax.axhline(0.5, color="#000000", linestyle="--", linewidth=3)
        ax.set_xlabel("Within-scan spread of log$_2$(p99), SD",
                      fontsize=FONT_SIZE_BASE - 4)
        ax.set_ylabel("Between-scan share of variance",
                      fontsize=FONT_SIZE_BASE - 4)
        ax.set_ylim(-0.02, 1.02)
        ax.set_title(ttl, fontsize=FONT_SIZE_TITLE - 8)
        style_axes(ax)
    axes[0].legend(loc="lower right", frameon=False,
                   fontsize=FONT_SIZE_LEGEND - 8)
    fig.suptitle("Where marker intensity variation lives\n"
                 "Scan and arm are the same factor here, so a high score is "
                 "consistent with batch AND with biology. A LOW score is the "
                 "informative direction.",
                 y=1.05, fontsize=FONT_SIZE_TITLE - 8)
    save_fig(fig, "F10_variance_between_vs_within_slide")

    write_csv(var_tbl, "10_variance_decomposition.csv")


# %% Cell 5 - single pass over the raw CSVs
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

sections = {}            # sid -> dict of light-weight arrays and summaries
pheno_profiles = []      # tidy per-section per-phenotype marker profile
ybin_rows = []           # within-section Y-decile marker medians
section_marker_med = []  # per-section median for every marker (all cells)
pheno_marker_med = []    # per-section per-phenotype median for every marker
pheno_counts = {}        # sid -> Series of phenotype counts
ido_cd68 = {}            # sid -> np.array of IDO1 values on CD68-lineage cells

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
    # scan identity: inventory first, Image column second, prefix last
    if sid in SCAN_FROM_INVENTORY:
        scan = SCAN_FROM_INVENTORY[sid]
        scan_src = "inventory"
    elif IMAGE_COL in df.columns and n:
        scan = str(df[IMAGE_COL].iloc[0])
        scan_src = "Image column"
    else:
        scan = f"scan_{prefix}"
        scan_src = "filename prefix"

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
    print(f"    n={n:,}  scan={scan}  ({scan_src})")
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
        "position": POS_FROM_INVENTORY.get(sid),
        "n_cells": n,
        "x_lo": x_lo, "x_hi": x_hi, "y_lo": y_lo, "y_hi": y_hi,
        "x_plot": x[idx_plot], "y_plot": y[idx_plot],
        "pheno_plot": df["_pheno"].to_numpy()[idx_plot],
    }

    # ---- per-section median of every marker, all cells ---------------------
    med = df[INTENSITY_COLS].median(axis=0)
    for c, v in med.items():
        section_marker_med.append({"sample_id": sid, "condition": cond,
                                   "scan": scan, "column": c,
                                   "median": float(v)})

    # ---- per-phenotype marker profile: mean, median and count --------------
    # The median is what the composition-controlled dimness test uses. Taking it
    # within a phenotype means the comparison between sections is like for like,
    # which the section-wide median cannot promise when composition ranges from
    # 89 percent structural to 36 percent Other.
    gb = df.groupby("_pheno")
    prof_mean = gb[INTENSITY_COLS].mean()
    prof_med = gb[INTENSITY_COLS].median()
    sizes = gb.size()
    pheno_counts[sid] = sizes
    for ph in prof_mean.index:
        for c in INTENSITY_COLS:
            pheno_profiles.append({
                "sample_id": sid, "condition": cond, "scan": scan,
                "phenotype": ph, "column": c,
                "mean": float(prof_mean.loc[ph, c]),
                "median": float(prof_med.loc[ph, c]),
                "n_cells": int(sizes.loc[ph]),
            })
            pheno_marker_med.append({
                "sample_id": sid, "scan": scan, "phenotype": ph,
                "column": c, "median": float(prof_med.loc[ph, c]),
                "n_cells": int(sizes.loc[ph]),
            })

    # ---- within-section Y deciles ------------------------------------------
    if ok.sum() > 0:
        yy = y[ok]
        edges = np.quantile(yy, np.linspace(0, 1, N_Y_BINS_WITHIN_SECTION + 1))
        edges = np.unique(edges)
        if len(edges) > 2:
            bins = np.clip(np.digitize(yy, edges[1:-1]), 0, len(edges) - 2)
            sub_df = df.loc[ok, INTENSITY_COLS]
            sub_df = sub_df.assign(_bin=bins)
            gb_y = sub_df.groupby("_bin")
            sizes_y = gb_y.size()
            meds_y = gb_y.median()
            for b in meds_y.index:
                if sizes_y.loc[b] < MIN_CELLS_PER_YBIN:
                    continue
                for c in INTENSITY_COLS:
                    ybin_rows.append({
                        "sample_id": sid, "condition": cond,
                        "y_decile": int(b), "n_cells": int(sizes_y.loc[b]),
                        "column": c, "median": float(meds_y.loc[b, c]),
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
scans = list(dict.fromkeys(sections[s]["scan"] for s in SAMPLE_ORDER))

# ---- position rank: from inventory, or recomputed as a fallback -------------
if any(sections[s]["position"] is None for s in SAMPLE_ORDER):
    print("\n    Position rank missing for at least one section. Recomputing "
          "from Y extents within each scan (rank 1 = smallest y_lo).")
    for sc in scans:
        members = sorted([s for s in SAMPLE_ORDER if sections[s]["scan"] == sc],
                         key=lambda z: sections[z]["y_lo"])
        for rank, s in enumerate(members, start=1):
            sections[s]["position"] = rank

SCAN_BOUNDARIES = [i - 0.5 for i in range(1, len(SAMPLE_ORDER))
                   if sections[SAMPLE_ORDER[i]]["scan"]
                   != sections[SAMPLE_ORDER[i - 1]]["scan"]]


def draw_scan_boundaries(ax):
    for b in SCAN_BOUNDARIES:
        ax.axvline(b, color=SCAN_BOUNDARY_COLOR, linewidth=SCAN_BOUNDARY_WIDTH,
                   zorder=6)


def sect_tick(s):
    p = sections[s]["position"]
    return f"{short_label(s)}\n{sections[s]['scan']} pos{p if p else 'na'}"


geom = pd.DataFrame([
    {"sample_id": s, "animal_id": short_label(s),
     "condition": sections[s]["condition"], "scan_id": sections[s]["scan"],
     "slide_position_rank": sections[s]["position"],
     "n_cells": sections[s]["n_cells"],
     "x_lo": sections[s]["x_lo"], "x_hi": sections[s]["x_hi"],
     "y_lo": sections[s]["y_lo"], "y_hi": sections[s]["y_hi"],
     "y_mid": 0.5 * (sections[s]["y_lo"] + sections[s]["y_hi"]),
     "y_span": sections[s]["y_hi"] - sections[s]["y_lo"]}
    for s in SAMPLE_ORDER
])


# %% Cell 6 - Q1: physical slide layout
# =============================================================================

banner("Q1 - PHYSICAL SLIDE LAYOUT")

print(f"    {len(scans)} distinct scan(s):")
for sc in scans:
    members = [s for s in SAMPLE_ORDER if sections[s]["scan"] == sc]
    print(f"        {sc}")
    for s in sorted(members, key=lambda z: sections[z]["y_lo"]):
        print(f"            {s:<12} pos {sections[s]['position']}  "
              f"Y {sections[s]['y_lo']:>10,.0f} to {sections[s]['y_hi']:>10,.0f}")

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
    shades = ["#08519C", "#3182BD", "#6BAED6", "#BDD7E7"] if cond == "D1MT" \
        else ["#A63603", "#E6550D", "#FD8D3C", "#FDBE85"]
    for j, s in enumerate(sorted(members, key=lambda z: sections[z]["y_lo"])):
        d = sections[s]
        ax.scatter(d["x_plot"], d["y_plot"], s=1.2,
                   color=shades[j % len(shades)], linewidths=0,
                   rasterized=True, label=short_label(s))
        ax.text(d["x_hi"] + 400, 0.5 * (d["y_lo"] + d["y_hi"]),
                f"{short_label(s)}\npos {d['position']}\n{d['n_cells']:,}",
                fontsize=FONT_SIZE_ANNOT - 8, va="center")
    ax.set_title(f"{cond}\n{sc}", fontsize=FONT_SIZE_TITLE - 10)
    ax.set_xlabel("X (µm)")
    ax.set_ylabel("Y (µm)")
    ax.set_aspect("equal")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
fig.suptitle("Physical section layout on each scan\n"
             "Four tissue sections per slide, tiled in Y", y=1.01)
save_fig(fig, "F11_slide_layout")

write_csv(geom, "11_section_geometry.csv")


# %% Cell 7 - Q3: slide-position artifact test, two ways
# =============================================================================

banner("Q3 - SLIDE POSITION ARTIFACT TEST")

print("    Run twice.")
print("      ALL CELLS      each marker's section-wide median, z-scored within")
print("                     scan. What revision 1 did. Composition-sensitive.")
print("      BY PHENOTYPE   each marker's median WITHIN a phenotype, z-scored")
print("                     within scan, using only phenotypes with at least")
print(f"                     {COMPOSITION_MIN_CELLS} cells in every section of")
print("                     that scan. Composition cancels.")
print("\n    Two statistics per section.")
print("      median z              how far below its slide-mates it sits. With")
print("                            four sections per scan the minimum possible")
print("                            z is -1.50, so read -1.0 as near the floor.")
print("      frac dimmest of four  fraction of comparisons where this section")
print("                            ranks lowest of its four. Under no position")
print("                            effect this is 0.25 exactly, with no")
print("                            distributional assumption. This is the")
print("                            statistic to quote.")

smed = pd.DataFrame(section_marker_med)
piv_all = smed.pivot_table(index="column", columns="sample_id", values="median")

pmed = pd.DataFrame(pheno_marker_med)

# ---- which phenotypes qualify for the composition-controlled test ------------
sub("Phenotypes qualifying for the composition-controlled test")
qualifying = {}
for sc in scans:
    members = [s for s in SAMPLE_ORDER if sections[s]["scan"] == sc]
    keep = []
    for ph in sorted(set(pmed["phenotype"])):
        counts = [int(pheno_counts.get(s, pd.Series(dtype=int)).get(ph, 0))
                  for s in members]
        if len(counts) == len(members) and min(counts) >= COMPOSITION_MIN_CELLS:
            keep.append(ph)
    qualifying[sc] = keep
    print(f"    {sc}: {len(keep)} phenotype(s)")
    for ph in keep:
        counts = [int(pheno_counts[s].get(ph, 0)) for s in members]
        print(f"        {ph:<28} min {min(counts):>7,} cells")
    if len(keep) < COMPOSITION_MIN_PHENOTYPES:
        print(f"        WARNING: fewer than {COMPOSITION_MIN_PHENOTYPES} "
              f"phenotypes qualify on {sc}. The composition-controlled test is "
              f"thin here.")

# ---- build both z matrices and both rank matrices ---------------------------
z_all_frames, r_all_frames = [], []
z_ph_frames, r_ph_frames = [], []
for sc in scans:
    members = [s for s in SAMPLE_ORDER if sections[s]["scan"] == sc]

    block = piv_all[members]
    z_all_frames.append(zscore_within(block))
    r_all_frames.append(rank_within(block))

    keep = qualifying[sc]
    if keep:
        sub_p = pmed.loc[(pmed["scan"] == sc) & (pmed["phenotype"].isin(keep))]
        wide = sub_p.pivot_table(index=["phenotype", "column"],
                                 columns="sample_id", values="median")
        wide = wide.reindex(columns=members)
        z_ph_frames.append(zscore_within(wide))
        r_ph_frames.append(rank_within(wide))

zmat = pd.concat(z_all_frames, axis=1)[SAMPLE_ORDER]
rmat = pd.concat(r_all_frames, axis=1)[SAMPLE_ORDER]
zmat_ph = (pd.concat(z_ph_frames, axis=1)[SAMPLE_ORDER]
           if z_ph_frames else pd.DataFrame())
rmat_ph = (pd.concat(r_ph_frames, axis=1)[SAMPLE_ORDER]
           if r_ph_frames else pd.DataFrame())

summary_rows = []
for s in SAMPLE_ORDER:
    col = zmat[s].dropna()
    rnk = rmat[s].dropna()
    rec = {
        "sample_id": s, "animal_id": short_label(s),
        "condition": sections[s]["condition"],
        "scan_id": sections[s]["scan"],
        "slide_position_rank": sections[s]["position"],
        "y_lo": sections[s]["y_lo"],
        "median_z_all_markers": float(col.median()) if len(col) else np.nan,
        "frac_markers_low": float((col < BRIGHTNESS_Z_LOW).mean()) if len(col) else np.nan,
        "frac_dimmest_of_scan": float((rnk == 1).mean()) if len(rnk) else np.nan,
        "n_markers": int(len(col)),
    }
    if len(zmat_ph):
        colp = zmat_ph[s].dropna()
        rnkp = rmat_ph[s].dropna()
        rec.update({
            "median_z_by_phenotype": float(colp.median()) if len(colp) else np.nan,
            "frac_low_by_phenotype": float((colp < BRIGHTNESS_Z_LOW).mean())
            if len(colp) else np.nan,
            "frac_dimmest_by_phenotype": float((rnkp == 1).mean())
            if len(rnkp) else np.nan,
            "n_phenotype_marker_pairs": int(len(colp)),
        })
    summary_rows.append(rec)
pos_tbl = pd.DataFrame(summary_rows).sort_values(["scan_id", "slide_position_rank"])

sub("Per-section brightness, both versions")
head = (f"    {'sample':<12}{'pos':>4}{'med z':>9}{'low%':>8}{'dim%':>8}"
        f"{'med z ph':>10}{'low% ph':>10}{'dim% ph':>10}")
print(head)
print("    " + "-" * (len(head) - 4))
for _, r in pos_tbl.iterrows():
    flag = ("   <-- globally dim"
            if r["median_z_all_markers"] < BRIGHTNESS_MEDIAN_Z_FLAG else "")
    mz = r.get("median_z_by_phenotype", np.nan)
    lp = r.get("frac_low_by_phenotype", np.nan)
    dp = r.get("frac_dimmest_by_phenotype", np.nan)
    print(f"    {r['sample_id']:<12}{int(r['slide_position_rank']):>4}"
          f"{r['median_z_all_markers']:>+9.2f}"
          f"{100*r['frac_markers_low']:>7.1f}%"
          f"{100*r['frac_dimmest_of_scan']:>7.1f}%"
          f"{mz:>+10.2f}{100*lp:>9.1f}%{100*dp:>9.1f}%{flag}")

print(f"\n    Reference: 'dim%' is 25.0% under no position effect, in both")
print(f"    versions. 'low%' has no clean null and is kept for continuity.")

sub("Do the two versions agree?")
if len(zmat_ph):
    for _, r in pos_tbl.iterrows():
        d_all = 100 * r["frac_dimmest_of_scan"]
        d_ph = 100 * r["frac_dimmest_by_phenotype"]
        note = ""
        if abs(d_all - d_ph) > 20:
            note = "   <-- versions disagree, composition is doing work here"
        print(f"    {r['sample_id']:<12} all cells {d_all:>5.1f}%   "
              f"by phenotype {d_ph:>5.1f}%{note}")
    print("\n    A section that is dim on ALL CELLS but not BY PHENOTYPE was")
    print("    never dim, it just holds different cells. A section dim on both")
    print("    is dim for staining reasons and is a genuine exclusion candidate.")
else:
    print("    Composition-controlled version unavailable, no qualifying "
          "phenotypes.")

print("\n    INTERPRETATION GUIDE")
print("    If 43102 and 43112 are dim across most of the 67 markers in BOTH")
print("    versions, the zero IDO1+ call in those two sections is a slide")
print("    position artifact. If they are dim only for IDO1 and its correlates,")
print("    or only in the all-cell version, it is more likely biology or")
print("    composition. Read the IDO1 row of F12 against the all-marker")
print("    distribution.")

# ---- F12: three panels ------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(42, 16),
                         gridspec_kw={"width_ratios": [1.15, 1.0, 1.0]})

# left: heatmap of z, key markers, all-cell version
ax = axes[0]
key_cols = [c for c in zmat.index if marker_base(c) in KEY_MARKERS]
key_cols = sorted(key_cols, key=lambda c: KEY_MARKERS.index(marker_base(c)))
z_key = zmat.loc[key_cols, SAMPLE_ORDER]
vmax = safe_vmax(z_key.to_numpy())
im = ax.imshow(z_key.to_numpy(), aspect="auto", cmap=CMAP_DIVERGING,
               norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax))
ax.set_xticks(np.arange(len(SAMPLE_ORDER)))
ax.set_xticklabels([sect_tick(s) for s in SAMPLE_ORDER],
                   fontsize=FONT_SIZE_TICK - 14)
ax.set_yticks(np.arange(len(key_cols)))
ax.set_yticklabels([marker_base(c) for c in key_cols], fontsize=FONT_SIZE_TICK - 10)
draw_scan_boundaries(ax)
cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
cb.set_label("z within scan", fontsize=FONT_SIZE_BASE - 8)
cb.ax.tick_params(labelsize=FONT_SIZE_TICK - 12)
ax.set_title("Marker medians, all cells", fontsize=FONT_SIZE_TITLE - 10)

# middle: distribution of z across ALL markers per section, both versions
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
if len(zmat_ph):
    for i, s in enumerate(SAMPLE_ORDER, start=1):
        v = zmat_ph[s].dropna()
        if len(v):
            ax.scatter([i], [float(v.median())], s=520, marker="D",
                       color="#000000", edgecolor="#FFFFFF", linewidth=3,
                       zorder=6)
ax.axhline(0, color="#000000", linestyle="--", linewidth=3)
ax.set_xticks(np.arange(1, len(SAMPLE_ORDER) + 1))
ax.set_xticklabels([short_label(s) for s in SAMPLE_ORDER],
                   rotation=45, ha="right", fontsize=FONT_SIZE_TICK - 8)
ax.set_ylabel("z of log$_2$ median", fontsize=FONT_SIZE_BASE - 6)
ax.set_title("Global brightness per section\n"
             "(box = all cells, diamond = median by phenotype)",
             fontsize=FONT_SIZE_TITLE - 14)
style_axes(ax)

# right: the rank statistic, with the 25 percent null drawn
ax = axes[2]
xs = np.arange(len(SAMPLE_ORDER))
w = 0.38
ax.bar(xs - w / 2, [100 * pos_tbl.set_index("sample_id").loc[s, "frac_dimmest_of_scan"]
                    for s in SAMPLE_ORDER],
       width=w, color="#BBBBBB", edgecolor="#FFFFFF", linewidth=2,
       label="all cells", zorder=3)
if len(zmat_ph):
    ax.bar(xs + w / 2,
           [100 * pos_tbl.set_index("sample_id").loc[s, "frac_dimmest_by_phenotype"]
            for s in SAMPLE_ORDER],
           width=w, color=FLAG_COLOR, edgecolor="#FFFFFF", linewidth=2,
           label="by phenotype", zorder=3)
ax.axhline(25.0, color="#000000", linestyle="--", linewidth=3.5)
ax.text(len(SAMPLE_ORDER) - 0.4, 27, "no position effect = 25%", ha="right",
        fontsize=FONT_SIZE_ANNOT - 8)
ax.set_xticks(xs)
ax.set_xticklabels([sect_tick(s) for s in SAMPLE_ORDER],
                   fontsize=FONT_SIZE_TICK - 14)
ax.set_ylabel("% of comparisons where this section\nis the dimmest of its four",
              fontsize=FONT_SIZE_BASE - 10)
ax.set_ylim(0, 105)
ax.set_title("Rank statistic", fontsize=FONT_SIZE_TITLE - 10)
style_axes(ax)
draw_scan_boundaries(ax)
ax.legend(frameon=False, fontsize=FONT_SIZE_LEGEND - 10, loc="upper center")

fig.suptitle("Is the position-1 section of each slide globally dim?\n"
             "Composition-controlled and rank-based versions included, because "
             "section composition ranges from 89% structural to 36% Other",
             y=1.04, fontsize=FONT_SIZE_TITLE - 6)
save_fig(fig, "F12_slide_position_artifact")

write_csv(pos_tbl, "12_slide_position_summary.csv")
write_csv(zmat.reset_index().rename(columns={"index": "column"}),
          "12b_marker_z_within_slide.csv")
if len(zmat_ph):
    write_csv(zmat_ph.reset_index(), "12d_marker_z_by_phenotype.csv")

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


# %% Cell 8 - Q4: IDO1 distribution on CD68+ macrophages
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
            "scan_id": sections[s]["scan"],
            "slide_position_rank": sections[s]["position"],
            "n_cd68_lineage": int(len(v)),
            "min": float(np.min(v)), "p25": float(np.percentile(v, 25)),
            "median": float(np.median(v)), "p75": float(np.percentile(v, 75)),
            "p95": float(np.percentile(v, 95)), "p99": float(np.percentile(v, 99)),
            "max": float(np.max(v)), "frac_zero": float(np.mean(v == 0)),
        })
    ido_sum = pd.DataFrame(rows)
    with pd.option_context("display.width", 240, "display.max_columns", 30):
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
        pos = sections[s]["position"]
        ax.text(bins[-1] * 1.01, base + 0.35,
                f"{short_label(s)}  pos{pos}  n={len(v):,}",
                fontsize=FONT_SIZE_ANNOT - 10, va="center")
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


# %% Cell 9 - Q5: what is in "Other"
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
    z_ph = z.loc[ph_use, key_cols]
    vmax = safe_vmax(z_ph.to_numpy())

    fig, ax = plt.subplots(figsize=(24, 14))
    im = ax.imshow(z_ph.to_numpy(), aspect="auto", cmap=CMAP_DIVERGING,
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


# %% Cell 10 - spatial maps of key immune populations
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
    ax.set_title(f"{short_label(s)}  ({d['condition']}, {d['scan']} pos{d['position']})"
                 f"\n{d['n_cells']:,} cells",
                 fontsize=FONT_SIZE_TITLE - 14)
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


# %% Cell 11 - Q6: nearest-neighbour pair feasibility
# =============================================================================

banner("Q6 - NEAREST-NEIGHBOUR PAIR FEASIBILITY")

pair_tbl = pd.DataFrame()
if not os.path.exists(cnt_path):
    print(f"    WARNING: {cnt_path} not found. Skipping Q6.")
else:
    print("    NOTE: all eight sections are counted here, including the two")
    print("    position-1 candidates. Feasibility is a property of the data as")
    print("    delivered; the exclusion decision belongs to Q3 and script 03.\n")
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
    for k, cond in enumerate(CONDITION_ORDER):
        ax = axes[k]
        m = pair_tbl.pivot_table(index="anchor", columns="target",
                                 values=f"n_animals_ok_{cond}")
        m = m.reindex(index=[a for a in NN_ANCHORS if a in m.index],
                      columns=[t for t in NN_TARGETS if t in m.columns])
        im = ax.imshow(m.to_numpy(dtype=float), cmap=CMAP_FEASIBILITY, vmin=0,
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


# %% Cell 12 - wrap up
# =============================================================================

banner("SUMMARY")

print(f"Sections loaded : {len(sections)}")
print(f"Distinct scans  : {len(scans)}")
print(f"Total cells     : {sum(sections[s]['n_cells'] for s in SAMPLE_ORDER):,}")

sub("Questions and where to read the answer")
print("  Q1 section tiling        -> F11_slide_layout, 11_section_geometry.csv")
print("  Q2 scan confound         -> F10_variance_between_vs_within_slide,")
print("                              10_variance_decomposition.csv")
print("                              quote icc_between_scan, not the naive ratio")
print("  Q3 position artifact     -> F12_slide_position_artifact, F12b,")
print("                              12_*.csv, 12d_marker_z_by_phenotype.csv")
print("                              quote frac_dimmest_by_phenotype against 25%")
print("  Q4 IDO1 distribution     -> F13_ido1_distribution_cd68, 13_*.csv")
print('  Q5 what is in "Other"    -> F14_phenotype_marker_profiles, 14_*.csv')
print("  spatial overview         -> F15_spatial_maps_key_populations")
print("  Q6 NN feasibility        -> F16_nn_pair_feasibility, 15_*.csv")

sub("Decisions still open after this run")
print("  1. Whether 43102 and 43112 stay in the analysis, and on what grounds.")
print("     The composition-controlled rank statistic is the version to decide")
print("     on. If a section is dim on all cells but not by phenotype, it was")
print("     never dim.")
print("  2. Whether IDO1 positivity is used as the supplied binary call or as a")
print("     within-scan continuous measure")
print("  3. Which anchor/target pairs go into the nearest-neighbour work")
print('  4. Whether "Other" is excluded, split, or carried as an unknown class')

sub("What NOT to conclude from Q2")
print("  A high between-scan score is not evidence of batch. Scan and arm are")
print("  the same factor in this design. The markers that score highest are the")
print("  ones expected to differ biologically between a granuloma-rich lung and")
print("  a treated lung. Only the low end of that table is interpretable, and")
print("  what it says is that those markers cannot be carrying much of an")
print("  acquisition difference.")

banner("DONE")
sys.stdout = _tee.terminal
_tee.close()
