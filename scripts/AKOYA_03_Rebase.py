#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - RE-BASELINE ON SIX SECTIONS
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 03 of the AKOYA analysis series.

WHY THIS SCRIPT EXISTS
    Scripts 01 and 02 established that the position-1 (bottom) section on each
    scan is globally dim: 67 to 69 percent of all 67 markers fall below z = -0.8
    relative to their slide-mates, versus 3 to 15 percent for every other section.
    Both also have the shortest Y spans on their slide and both start at Y = 0.7,
    flush against the scan boundary. Neither can support an intensity-based claim.

    G3_43102 (D1MT) and G4_43112 (Untreated) are therefore excluded here. One from
    each arm, so the design stays balanced at 3 versus 3. The exclusion is declared
    once, in EXCLUDE_SECTIONS below, and every downstream script should read it
    from this file rather than re-deciding.

WHAT THIS SCRIPT DOES
    A  Re-baseline composition, IDO1 fraction and the batch decomposition on the
       six retained sections.
    B  Back out the IDO1 intensity threshold implied by the vendor's binary
       CD68+IDO1+ / CD68+IDO1- call, and test whether it was a fixed global cut,
       a per-section adaptive cut, or something that used more than IDO1 alone.
    C  Threshold sensitivity: fraction of CD68-lineage cells called positive as a
       function of cutoff, to show how wide the insensitive window is.
    D  Rebuild the IDO1 distribution figures on six sections, with distinct point
       markers per section so lines can be told apart.
    E  Replace the Y-decile edge test with distance to the tissue boundary,
       restricted to CD68-lineage cells, with DAPI as a control channel.
    F  Spatial maps showing ALL phenotypes, including endothelium, epithelium and
       the SMA-high "Other" stroma, so immune positions can be read against tissue
       architecture.

    Still diagnostic. Nothing is re-phenotyped or written back to source CSVs.

OUTPUTS
    figures/  F17 .. F23   (PDF + PNG, 300 DPI)
    tables/   20 .. 26     (CSV, LF line endings)

USAGE
    conda activate sc_pre
    python AKOYA_03_Rebaseline.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================
# ALL TUNABLE PARAMETERS LIVE HERE
# =============================================================================

DATA_DIR = "/master/jlehle/WORKING/AKOYA/data"
OUT_DIR = "/master/jlehle/WORKING/AKOYA/rebaseline"

# ---- THE EXCLUSION, DECLARED ONCE -------------------------------------------
# Downstream scripts should import or copy this block verbatim.
EXCLUDE_SECTIONS = {
    "G3_43102": ("Position-1 section on the G3 scan. Globally dim: 67.2% of 67 "
                 "markers below z=-0.8 within slide, median z=-1.04. IDO1 p99=8.3 "
                 "on CD68-lineage cells, below the p25 of intact slide-mates. "
                 "Cannot support intensity-based claims."),
    "G4_43112": ("Position-1 section on the G4 scan. Globally dim: 68.7% of 67 "
                 "markers below z=-0.8 within slide, median z=-1.21. IDO1 max=8.3 "
                 "across 2,449 macrophages. Cannot support intensity-based claims."),
}

CONDITION_ORDER = ["D1MT", "Untreated"]
CONDITION_COLORS = {"D1MT": "#2C7FB8", "Untreated": "#D95F02"}
GROUP_MAP = {"G3": "D1MT", "G4": "Untreated"}

# ---- column names -----------------------------------------------------------
INTENSITY_SUFFIX = ": Mean"
CENTROID_X_PREFIX = "Centroid X"
CENTROID_Y_PREFIX = "Centroid Y"
PHENOTYPE_COL = "Phenotypes"
IMAGE_COL = "Image"

IDO1_COL = "IDO1: Cytoplasm: Mean"
DAPI_COL = "DAPI: Nucleus: Mean"
HK3_COL = "3-Hydroxykynurenine: Cytoplasm: Mean"

# ---- phenotypes -------------------------------------------------------------
PHENOTYPE_ORDER = [
    "CD68+IDO1+ Macrophages", "CD68+IDO1- Macrophages", "CD163+ Macrophages",
    "Neutrophils", "Helper T cells", "CD4- T cells", "Tregs",
    "B cells", "Plasma cells", "Endothelial cells",
    "Epithelial/Tumor cells", "Other",
]

# Structural classes now get real colours instead of grey, so tissue architecture
# is readable on the spatial maps.
PHENOTYPE_COLORS = {
    "CD68+IDO1+ Macrophages": "#B2182B",
    "CD68+IDO1- Macrophages": "#EF8A62",
    "CD163+ Macrophages": "#FDBE85",
    "Neutrophils": "#7B3294",
    "Helper T cells": "#1B7837",
    "CD4- T cells": "#7FBC41",
    "Tregs": "#00441B",
    "B cells": "#2166AC",
    "Plasma cells": "#67A9CF",
    "Endothelial cells": "#E7298A",
    "Epithelial/Tumor cells": "#66C2A5",
    "Other": "#A6761D",
}

IMMUNE_PHENOTYPES = [
    "CD68+IDO1+ Macrophages", "CD68+IDO1- Macrophages", "CD163+ Macrophages",
    "Neutrophils", "Helper T cells", "CD4- T cells", "Tregs",
    "B cells", "Plasma cells",
]
STRUCTURAL_PHENOTYPES = ["Endothelial cells", "Epithelial/Tumor cells", "Other"]

IDO1_POS_PHENO = "CD68+IDO1+ Macrophages"
IDO1_NEG_PHENO = "CD68+IDO1- Macrophages"
CD68_LINEAGE = [IDO1_POS_PHENO, IDO1_NEG_PHENO]

# ---- threshold sensitivity --------------------------------------------------
N_CUTOFF_STEPS = 200          # grid resolution on the log1p scale
PLATEAU_TOL = 0.05            # a cutoff window is "insensitive" if the called
                              # fraction moves less than this across it

# ---- tissue boundary distance ----------------------------------------------
BOUNDARY_GRID_UM = 50         # occupancy grid pitch for the distance transform
BOUNDARY_CLOSE_ITER = 2       # binary closing passes to fill small interior gaps
N_BOUNDARY_BINS = 12          # distance bins for the gradient profile
MIN_CELLS_PER_BOUNDARY_BIN = 30

# ---- key markers ------------------------------------------------------------
KEY_MARKERS = [
    "IDO1", "3-Hydroxykynurenine", "CD68", "CD163", "CD206", "CD4", "CD8",
    "CD3e", "CD20", "CD79a", "CD21", "MPO", "CD11b", "FoxP3", "iNOS",
    "Arginase-1", "IFNG", "Granzyme-B", "HLA-DR", "PD-1", "PD-L1", "DAPI",
]

# ---- plotting ---------------------------------------------------------------
MAX_POINTS_PER_SECTION = 80000
POINT_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]
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
BACKGROUND_COLOR = "#ECECEC"

USE_FLOAT32 = True


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
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

try:
    from scipy import ndimage as ndi
    HAVE_SCIPY = True
except Exception as _e:
    HAVE_SCIPY = False
    _scipy_err = _e

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
DIVERGING = LinearSegmentedColormap.from_list("div", ["#2166AC", "#F7F7F7", "#B2182B"])


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


def ascii_safe(s):
    return "".join(ch for ch in str(s) if ord(ch) <= 127).strip()


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


_tee = Tee(os.path.join(TAB_DIR, "00_rebaseline_report.txt"))
sys.stdout = _tee

banner("AKOYA RE-BASELINE ON SIX SECTIONS")
print(f"Run time : {datetime.now().isoformat(timespec='seconds')}")
print(f"Data     : {DATA_DIR}")
print(f"Output   : {OUT_DIR}")
if not HAVE_SCIPY:
    print(f"\n    WARNING: scipy unavailable ({_scipy_err}). The tissue boundary "
          f"distance analysis will be skipped.")

sub("EXCLUDED SECTIONS")
for sid, reason in EXCLUDE_SECTIONS.items():
    print(f"    {sid}")
    print(f"        {reason}")


# %% Cell 3 - single pass over the retained sections
# =============================================================================

banner("LOADING RETAINED SECTIONS")

csv_paths = sorted(glob.glob(os.path.join(DATA_DIR, "*.csv")))
retained = [p for p in csv_paths
            if os.path.splitext(os.path.basename(p))[0] not in EXCLUDE_SECTIONS]
print(f"    {len(csv_paths)} files found, {len(EXCLUDE_SECTIONS)} excluded, "
      f"{len(retained)} retained")

if not retained:
    print("ERROR: nothing retained.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

hdr = pd.read_csv(retained[0], nrows=0)
ALL_COLS = hdr.columns.tolist()
XCOL = find_col(ALL_COLS, CENTROID_X_PREFIX)
YCOL = find_col(ALL_COLS, CENTROID_Y_PREFIX)
INTENSITY_COLS = [c for c in ALL_COLS if c.endswith(INTENSITY_SUFFIX)]
dtype_map = ({c: "float32" for c in INTENSITY_COLS + [XCOL, YCOL]}
             if USE_FLOAT32 else None)

sections = {}
pheno_counts = []
section_marker_med = []
pheno_profiles = []
cd68_cells = {}     # sid -> DataFrame of CD68-lineage cells (x, y, ido1, dapi, label)

for path in retained:
    sid = os.path.splitext(os.path.basename(path))[0]
    cond = GROUP_MAP.get(sid.split("_")[0], "UNKNOWN")
    print(f"\n[{sid}]  {cond}")
    try:
        df = pd.read_csv(path, dtype=dtype_map, low_memory=False)
    except Exception as e:
        print(f"    ERROR reading: {e}. Skipping and continuing.")
        continue

    n = len(df)
    scan = str(df[IMAGE_COL].iloc[0]) if IMAGE_COL in df.columns and n else "unknown"
    df["_pheno"] = df[PHENOTYPE_COL].map(ascii_safe) if PHENOTYPE_COL in df.columns \
        else "Unknown"

    x = pd.to_numeric(df[XCOL], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(df[YCOL], errors="coerce").to_numpy(dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    print(f"    n={n:,}   scan={scan.split(' - ')[0]}")

    idx = np.flatnonzero(ok)
    idx_plot = (rng.choice(idx, size=MAX_POINTS_PER_SECTION, replace=False)
                if len(idx) > MAX_POINTS_PER_SECTION else idx)

    sections[sid] = {
        "condition": cond, "scan": scan, "n_cells": n,
        "x_plot": x[idx_plot], "y_plot": y[idx_plot],
        "pheno_plot": df["_pheno"].to_numpy()[idx_plot],
        "x_lo": float(np.nanmin(x)), "x_hi": float(np.nanmax(x)),
        "y_lo": float(np.nanmin(y)), "y_hi": float(np.nanmax(y)),
    }

    vc = df["_pheno"].value_counts()
    for ph, c in vc.items():
        pheno_counts.append({"sample_id": sid, "condition": cond,
                             "phenotype": ph, "n_cells": int(c),
                             "pct_of_section": 100.0 * c / n})

    med = df[INTENSITY_COLS].median(axis=0)
    for c, v in med.items():
        section_marker_med.append({"sample_id": sid, "condition": cond,
                                   "column": c, "median": float(v)})

    prof = df.groupby("_pheno")[INTENSITY_COLS].mean()
    for ph, row in prof.iterrows():
        for c, v in row.items():
            pheno_profiles.append({"sample_id": sid, "condition": cond,
                                   "phenotype": ph, "column": c, "mean": float(v)})

    # ---- CD68-lineage cell table -------------------------------------------
    m = df["_pheno"].isin(CD68_LINEAGE).to_numpy() & ok
    keep = pd.DataFrame({
        "x": x[m], "y": y[m],
        "label": df.loc[m, "_pheno"].to_numpy(),
        "ido1": pd.to_numeric(df.loc[m, IDO1_COL], errors="coerce").to_numpy()
        if IDO1_COL in df.columns else np.nan,
        "dapi": pd.to_numeric(df.loc[m, DAPI_COL], errors="coerce").to_numpy()
        if DAPI_COL in df.columns else np.nan,
        "hk3": pd.to_numeric(df.loc[m, HK3_COL], errors="coerce").to_numpy()
        if HK3_COL in df.columns else np.nan,
    })
    cd68_cells[sid] = keep
    print(f"    CD68-lineage cells: {len(keep):,}  "
          f"({int((keep['label'] == IDO1_POS_PHENO).sum()):,} IDO1+)")

    # keep all coordinates for the boundary transform
    sections[sid]["x_all"] = x[ok]
    sections[sid]["y_all"] = y[ok]

    del df
    gc.collect()

cond_rank = {c: i for i, c in enumerate(CONDITION_ORDER)}
SAMPLE_ORDER = sorted(sections, key=lambda s: (cond_rank.get(sections[s]["condition"], 9), s))
MARKER_OF = {s: POINT_MARKERS[i % len(POINT_MARKERS)]
             for i, s in enumerate(SAMPLE_ORDER)}

# shade within condition so sections are distinguishable but still group-coded
SHADES = {"D1MT": ["#08519C", "#3182BD", "#6BAED6", "#BDD7E7"],
          "Untreated": ["#A63603", "#E6550D", "#FD8D3C", "#FDBE85"]}
COLOR_OF = {}
_seen = {c: 0 for c in CONDITION_ORDER}
for s in SAMPLE_ORDER:
    c = sections[s]["condition"]
    pal = SHADES.get(c, ["#999999"])
    COLOR_OF[s] = pal[_seen.get(c, 0) % len(pal)]
    _seen[c] = _seen.get(c, 0) + 1


def section_legend(ax, **kw):
    handles = [Line2D([0], [0], color=COLOR_OF[s], marker=MARKER_OF[s],
                      markersize=16, linewidth=4,
                      label=f"{short_label(s)} ({sections[s]['condition']})")
               for s in SAMPLE_ORDER]
    ax.legend(handles=handles, frameon=False,
              fontsize=kw.pop("fontsize", FONT_SIZE_LEGEND - 10), **kw)


# %% Cell 4 - A: re-baselined composition and IDO1 fraction
# =============================================================================

banner("A - RE-BASELINED COMPOSITION AND IDO1 FRACTION")

pc = pd.DataFrame(pheno_counts)
counts = (pc.pivot_table(index="phenotype", columns="sample_id", values="n_cells",
                         aggfunc="sum")
          .reindex(index=[p for p in PHENOTYPE_ORDER],
                   columns=SAMPLE_ORDER).fillna(0.0))
PHENO_USE = [p for p in PHENOTYPE_ORDER if p in counts.index]
counts = counts.loc[PHENO_USE]
pcts = 100.0 * counts / counts.sum(axis=0)

print(f"    Total cells retained: {int(counts.sum().sum()):,}")
for c in CONDITION_ORDER:
    mem = [s for s in SAMPLE_ORDER if sections[s]["condition"] == c]
    print(f"      {c:<12} {len(mem)} sections, "
          f"{int(counts[mem].sum().sum()):,} cells")

num = counts.loc[IDO1_POS_PHENO, SAMPLE_ORDER].to_numpy(float)
den = counts.loc[CD68_LINEAGE, SAMPLE_ORDER].sum(axis=0).to_numpy(float)
frac = np.where(den > 0, 100.0 * num / den, np.nan)

ido_frac = pd.DataFrame({
    "sample_id": SAMPLE_ORDER,
    "animal_id": [short_label(s) for s in SAMPLE_ORDER],
    "condition": [sections[s]["condition"] for s in SAMPLE_ORDER],
    "n_ido1_pos": num.astype(int), "n_cd68_lineage": den.astype(int),
    "pct_ido1_pos": frac,
})
sub("IDO1+ as a fraction of CD68-lineage macrophages (six sections)")
print(ido_frac.to_string(index=False))

xpos = np.arange(len(SAMPLE_ORDER))
cond_labels = [sections[s]["condition"] for s in SAMPLE_ORDER]
x_labels = [short_label(s) for s in SAMPLE_ORDER]

fig, axes = plt.subplots(1, 2, figsize=(30, 13),
                         gridspec_kw={"width_ratios": [1.35, 1.0]})
ax = axes[0]
bottom = np.zeros(len(SAMPLE_ORDER))
for ph in PHENO_USE:
    v = pcts.loc[ph, SAMPLE_ORDER].to_numpy(float)
    ax.bar(xpos, v, bottom=bottom, color=PHENOTYPE_COLORS.get(ph, "#999999"),
           edgecolor="#FFFFFF", linewidth=1.5, label=ph, zorder=3)
    bottom += v
ax.set_xticks(xpos)
ax.set_xticklabels([f"{a}\n{c}" for a, c in zip(x_labels, cond_labels)],
                   fontsize=FONT_SIZE_TICK - 6)
ax.set_ylabel("Percent of section"); ax.set_ylim(0, 100)
ax.set_title("Composition, six retained sections", fontsize=FONT_SIZE_TITLE - 6)
style_axes(ax, ygrid=False)
ax.legend(bbox_to_anchor=(1.01, 1.0), loc="upper left", frameon=False,
          fontsize=FONT_SIZE_LEGEND - 12)

ax = axes[1]
ax.bar(xpos, np.nan_to_num(frac),
       color=[CONDITION_COLORS.get(c, "#999") for c in cond_labels],
       edgecolor="#FFFFFF", linewidth=2, zorder=3)
for i in range(len(SAMPLE_ORDER)):
    ax.text(i, frac[i] + 2, f"{frac[i]:.1f}%\n{int(num[i]):,}/{int(den[i]):,}",
            ha="center", va="bottom", fontsize=FONT_SIZE_ANNOT - 8)
ax.set_xticks(xpos)
ax.set_xticklabels([f"{a}\n{c}" for a, c in zip(x_labels, cond_labels)],
                   fontsize=FONT_SIZE_TICK - 6)
ax.set_ylabel("IDO1+ % of CD68-lineage"); ax.set_ylim(0, 110)
ax.set_title("IDO1+ fraction", fontsize=FONT_SIZE_TITLE - 6)
style_axes(ax)
fig.suptitle("Re-baselined on six sections "
             "(position-1 sections of each scan excluded)", y=1.03,
             fontsize=FONT_SIZE_TITLE - 4)
save_fig(fig, "F17_rebaselined_composition")

write_csv(pc, "20_rebaselined_phenotype_counts.csv")
write_csv(ido_frac, "21_ido1_fraction_six_sections.csv")

# ---- batch decomposition on six sections ------------------------------------
sub("Between-slide vs within-slide variance, recomputed on six sections")
smed = pd.DataFrame(section_marker_med)
smed = smed.loc[smed["median"] > 0].copy()
smed["log2_med"] = np.log2(smed["median"])
rows = []
for col, g in smed.groupby("column"):
    sm = g.groupby("condition")["log2_med"].mean()
    if len(sm) < 2:
        continue
    within = g.groupby("condition")["log2_med"].var(ddof=1).mean()
    between = float(np.var(sm.to_numpy(), ddof=1))
    tot = within + between
    rows.append({"column": col, "marker": marker_base(col),
                 "is_key_marker": marker_base(col) in KEY_MARKERS,
                 "within_slide_var": within, "between_slide_var": between,
                 "frac_between_slide": between / tot if tot > 0 else np.nan,
                 "log2_slide_diff": float(sm.get(CONDITION_ORDER[1], np.nan)
                                          - sm.get(CONDITION_ORDER[0], np.nan))})
var6 = pd.DataFrame(rows).sort_values("frac_between_slide", ascending=False)
print("    Key markers, recomputed:\n")
for _, r in var6.loc[var6["is_key_marker"]].iterrows():
    print(f"    {r['column']:<48} frac_between={r['frac_between_slide']:.3f}  "
          f"log2 diff={r['log2_slide_diff']:+.2f}")
print("\n    NOTE: excluding the two dim sections removes variance that was")
print("    inflating the within-slide term, so these fractions should be read")
print("    as the corrected values. Slide is still confounded with condition.")
write_csv(var6, "22_variance_decomposition_six_sections.csv")


# %% Cell 5 - B: back out the IDO1 threshold behind the vendor's binary call
# =============================================================================

banner("B - IDO1 THRESHOLD BACK-OUT")

print("    For each section: the lowest IDO1 among cells labelled IDO1+, and the")
print("    highest IDO1 among cells labelled IDO1-. If these meet at the same")
print("    value in every section, the call was a fixed global cut. If they meet")
print("    at different values, it was per-section adaptive. If the ranges")
print("    OVERLAP, the call used more than IDO1 alone.\n")

back_rows = []
for s in SAMPLE_ORDER:
    d = cd68_cells.get(s)
    if d is None or not len(d):
        continue
    pos = d.loc[d["label"] == IDO1_POS_PHENO, "ido1"].dropna().to_numpy()
    neg = d.loc[d["label"] == IDO1_NEG_PHENO, "ido1"].dropna().to_numpy()
    rec = {"sample_id": s, "animal_id": short_label(s),
           "condition": sections[s]["condition"],
           "n_pos": len(pos), "n_neg": len(neg),
           "min_pos": float(np.min(pos)) if len(pos) else np.nan,
           "max_neg": float(np.max(neg)) if len(neg) else np.nan,
           "p01_pos": float(np.percentile(pos, 1)) if len(pos) else np.nan,
           "p99_neg": float(np.percentile(neg, 99)) if len(neg) else np.nan}
    if len(pos) and len(neg):
        rec["separable"] = bool(rec["min_pos"] > rec["max_neg"])
        # fraction of negatives that sit above the dimmest positive
        rec["frac_neg_above_min_pos"] = float(np.mean(neg > rec["min_pos"]))
        rec["frac_pos_below_max_neg"] = float(np.mean(pos < rec["max_neg"]))
        # midpoint of the gap, if separable
        rec["implied_cut"] = (0.5 * (rec["min_pos"] + rec["max_neg"])
                              if rec["separable"] else np.nan)
    back_rows.append(rec)

back = pd.DataFrame(back_rows)
with pd.option_context("display.width", 240, "display.max_columns", 30):
    print(back.to_string(index=False))

sub("Verdict")
if len(back) and "separable" in back.columns:
    sep = back["separable"].fillna(False)
    if sep.all():
        cuts = back["implied_cut"].dropna()
        spread = cuts.max() / cuts.min() if cuts.min() > 0 else np.inf
        print(f"    Every section is cleanly separable on IDO1 alone.")
        print(f"    Implied cut per section: "
              f"{', '.join(f'{v:.2f}' for v in cuts)}")
        if spread < 1.5:
            print(f"    Spread is {spread:.2f}x -> consistent with a FIXED GLOBAL "
                  f"cut applied across sections.")
        else:
            print(f"    Spread is {spread:.2f}x -> the cut was PER-SECTION "
                  f"ADAPTIVE. Any cross-section comparison of the binary call")
            print(f"    carries this, and it must be stated in methods.")
    else:
        bad = back.loc[~sep, "sample_id"].tolist()
        print(f"    Ranges OVERLAP in: {bad}")
        print(f"    The binary call therefore used more than IDO1 intensity")
        print(f"    alone (co-markers, spatial smoothing, or a classifier).")
        print(f"    It cannot be reproduced as a simple IDO1 gate, and we should")
        print(f"    ask Akoya what the classifier actually used before relying")
        print(f"    on the call for anything quantitative.")
write_csv(back, "23_ido1_threshold_backout.csv")


# %% Cell 6 - C: threshold sensitivity, and D: rebuilt distributions
# =============================================================================

banner("C / D - THRESHOLD SENSITIVITY AND IDO1 DISTRIBUTIONS")

all_ido = np.concatenate([cd68_cells[s]["ido1"].dropna().to_numpy()
                          for s in SAMPLE_ORDER if len(cd68_cells.get(s, []))])
lo, hi = 0.0, float(np.log1p(np.percentile(all_ido, 99.9)))
grid = np.linspace(lo, hi, N_CUTOFF_STEPS)

sens_rows = []
for s in SAMPLE_ORDER:
    v = cd68_cells[s]["ido1"].dropna().to_numpy()
    if not len(v):
        continue
    lv = np.log1p(v)
    for g in grid:
        sens_rows.append({"sample_id": s, "condition": sections[s]["condition"],
                          "log1p_cutoff": float(g),
                          "frac_above": float(np.mean(lv > g))})
sens = pd.DataFrame(sens_rows)

# widest cutoff window over which every section's called fraction moves < tol
sw = sens.pivot_table(index="log1p_cutoff", columns="sample_id", values="frac_above")
sw = sw[SAMPLE_ORDER]
best_lo = best_hi = np.nan
i = 0
while i < len(sw):
    j = i
    while j + 1 < len(sw):
        block = sw.iloc[i:j + 2]
        if (block.max() - block.min()).max() > PLATEAU_TOL:
            break
        j += 1
    if not np.isnan(best_lo) and (sw.index[j] - sw.index[i]) <= (best_hi - best_lo):
        pass
    elif np.isnan(best_lo) or (sw.index[j] - sw.index[i]) > (best_hi - best_lo):
        best_lo, best_hi = float(sw.index[i]), float(sw.index[j])
    i = max(j, i + 1)

print(f"    Widest cutoff window where every section's called fraction moves")
print(f"    by less than {PLATEAU_TOL:.2f}:")
print(f"        log1p scale : {best_lo:.2f} to {best_hi:.2f}")
print(f"        raw IDO1    : {np.expm1(best_lo):.2f} to {np.expm1(best_hi):.2f}")
print(f"    Any cutoff inside that window gives essentially the same answer,")
print(f"    which is what makes the group difference robust to the gate.")

fig, axes = plt.subplots(1, 3, figsize=(42, 13))

# --- ECDF ---
ax = axes[0]
for s in SAMPLE_ORDER:
    v = np.sort(cd68_cells[s]["ido1"].dropna().to_numpy())
    if not len(v):
        continue
    yv = np.arange(1, len(v) + 1) / len(v)
    step = max(1, len(v) // 25)
    ax.plot(np.log1p(v), yv, linewidth=4.5, color=COLOR_OF[s], alpha=0.9)
    ax.plot(np.log1p(v)[::step], yv[::step], linestyle="none",
            marker=MARKER_OF[s], markersize=15, color=COLOR_OF[s])
if np.isfinite(best_lo):
    ax.axvspan(best_lo, best_hi, color="#000000", alpha=0.08, zorder=0)
ax.set_xlabel("log(1 + IDO1 cytoplasm)")
ax.set_ylabel("Cumulative fraction of CD68+ cells")
ax.set_ylim(0, 1.02)
ax.set_title("Cumulative distribution", fontsize=FONT_SIZE_TITLE - 8)
style_axes(ax)
section_legend(ax, loc="lower right")

# --- threshold sensitivity ---
ax = axes[1]
for s in SAMPLE_ORDER:
    d = sens.loc[sens["sample_id"] == s].sort_values("log1p_cutoff")
    ax.plot(d["log1p_cutoff"], 100 * d["frac_above"], linewidth=4.5,
            color=COLOR_OF[s], alpha=0.9)
    step = max(1, len(d) // 18)
    ax.plot(d["log1p_cutoff"].to_numpy()[::step],
            100 * d["frac_above"].to_numpy()[::step],
            linestyle="none", marker=MARKER_OF[s], markersize=15,
            color=COLOR_OF[s])
if np.isfinite(best_lo):
    ax.axvspan(best_lo, best_hi, color="#000000", alpha=0.08, zorder=0)
    ax.text(0.5 * (best_lo + best_hi), 103, "insensitive window",
            ha="center", fontsize=FONT_SIZE_ANNOT - 8)
ax.set_xlabel("log(1 + IDO1) cutoff")
ax.set_ylabel("% of CD68+ cells called positive")
ax.set_ylim(0, 108)
ax.set_title("Threshold sensitivity", fontsize=FONT_SIZE_TITLE - 8)
style_axes(ax)

# --- ridgeline ---
ax = axes[2]
bins = np.linspace(0, hi, 70)
for i, s in enumerate(SAMPLE_ORDER):
    v = cd68_cells[s]["ido1"].dropna().to_numpy()
    if not len(v):
        continue
    h, e = np.histogram(np.log1p(v), bins=bins, density=True)
    h = h / (h.max() if h.max() > 0 else 1.0)
    ctr = 0.5 * (e[:-1] + e[1:])
    base = len(SAMPLE_ORDER) - i
    ax.fill_between(ctr, base, base + 0.9 * h, color=COLOR_OF[s], alpha=0.75,
                    linewidth=0)
    ax.plot(ctr, base + 0.9 * h, color="#333333", linewidth=2.5)
    ax.text(hi * 1.01, base + 0.35, f"{short_label(s)}  n={len(v):,}",
            fontsize=FONT_SIZE_ANNOT - 8, va="center")
ax.set_yticks([]); ax.set_xlabel("log(1 + IDO1 cytoplasm)")
ax.set_title("Per-section distribution", fontsize=FONT_SIZE_TITLE - 8)
for sp in ["top", "right", "left"]:
    ax.spines[sp].set_visible(False)

fig.suptitle("IDO1 on CD68-lineage macrophages, six sections\n"
             "Shaded band = cutoff range over which the called fraction barely "
             "moves", y=1.04, fontsize=FONT_SIZE_TITLE - 4)
save_fig(fig, "F18_ido1_threshold_and_distribution")

write_csv(sens, "24_threshold_sensitivity.csv")


# %% Cell 7 - E: tissue boundary distance gradient
# =============================================================================

banner("E - TISSUE BOUNDARY DISTANCE GRADIENT")

bgrad = pd.DataFrame()
if not HAVE_SCIPY:
    print("    SKIPPED: scipy not available.")
else:
    print(f"    Occupancy grid at {BOUNDARY_GRID_UM} um, distance transform from "
          f"the tissue exterior.\n")
    rows = []
    for s in SAMPLE_ORDER:
        xa, ya = sections[s]["x_all"], sections[s]["y_all"]
        if not len(xa):
            continue
        x0, y0 = xa.min(), ya.min()
        gx = ((xa - x0) / BOUNDARY_GRID_UM).astype(int)
        gy = ((ya - y0) / BOUNDARY_GRID_UM).astype(int)
        occ = np.zeros((gy.max() + 3, gx.max() + 3), dtype=bool)
        occ[gy + 1, gx + 1] = True
        occ = ndi.binary_closing(occ, iterations=BOUNDARY_CLOSE_ITER)
        occ = ndi.binary_fill_holes(occ)
        dist = ndi.distance_transform_edt(occ) * BOUNDARY_GRID_UM

        d = cd68_cells[s]
        cgx = ((d["x"].to_numpy() - x0) / BOUNDARY_GRID_UM).astype(int) + 1
        cgy = ((d["y"].to_numpy() - y0) / BOUNDARY_GRID_UM).astype(int) + 1
        cgx = np.clip(cgx, 0, occ.shape[1] - 1)
        cgy = np.clip(cgy, 0, occ.shape[0] - 1)
        dd = dist[cgy, cgx]

        edges = np.quantile(dd, np.linspace(0, 1, N_BOUNDARY_BINS + 1))
        edges = np.unique(edges)
        if len(edges) < 3:
            print(f"    WARNING: {s} has too little depth range, skipped")
            continue
        b = np.clip(np.digitize(dd, edges[1:-1]), 0, len(edges) - 2)
        tmp = d.assign(_bin=b, _dist=dd)
        for bb, g in tmp.groupby("_bin"):
            if len(g) < MIN_CELLS_PER_BOUNDARY_BIN:
                continue
            rows.append({
                "sample_id": s, "condition": sections[s]["condition"],
                "depth_bin": int(bb), "n_cells": len(g),
                "median_depth_um": float(g["_dist"].median()),
                "median_ido1": float(g["ido1"].median()),
                "median_dapi": float(g["dapi"].median()),
                "median_hk3": float(g["hk3"].median()),
                "pct_ido1_pos": 100.0 * float((g["label"] == IDO1_POS_PHENO).mean()),
            })
        print(f"    {s:<12} depth range {dd.min():.0f} to {dd.max():.0f} um, "
              f"{len(np.unique(b))} bins")
    bgrad = pd.DataFrame(rows)

    if len(bgrad):
        fig, axes = plt.subplots(1, 3, figsize=(42, 13))
        panels = [("median_ido1", "Median IDO1 (cytoplasm)", "IDO1"),
                  ("median_dapi", "Median DAPI (nucleus)", "DAPI control"),
                  ("pct_ido1_pos", "% called IDO1+", "Binary call")]
        for ax, (colname, ylab, ttl) in zip(axes, panels):
            for s in SAMPLE_ORDER:
                d = bgrad.loc[bgrad["sample_id"] == s].sort_values("median_depth_um")
                if not len(d):
                    continue
                ax.plot(d["median_depth_um"], d[colname], linewidth=4.5,
                        marker=MARKER_OF[s], markersize=18, color=COLOR_OF[s],
                        alpha=0.9)
            ax.set_xlabel("Distance from tissue boundary (µm)")
            ax.set_ylabel(ylab)
            ax.set_title(ttl, fontsize=FONT_SIZE_TITLE - 8)
            style_axes(ax)
        section_legend(axes[0], loc="best")
        fig.suptitle("CD68-lineage macrophages by depth into the tissue\n"
                     "If IDO1 moves while DAPI stays flat, the gradient is real "
                     "signal rather than tissue quality", y=1.04,
                     fontsize=FONT_SIZE_TITLE - 4)
        save_fig(fig, "F19_boundary_distance_gradient")
        write_csv(bgrad, "25_boundary_distance_gradient.csv")
    else:
        print("    WARNING: no usable depth bins.")


# %% Cell 8 - F: spatial maps, all phenotypes
# =============================================================================

banner("F - SPATIAL MAPS")


def spatial_panel(ax, sid, show, sizes=None, bg=True):
    d = sections[sid]
    ph = d["pheno_plot"]
    if bg:
        m = ~np.isin(ph, show)
        ax.scatter(d["x_plot"][m], d["y_plot"][m], s=0.6,
                   color=BACKGROUND_COLOR, linewidths=0, rasterized=True)
    for p in show:
        m = ph == p
        if m.sum() == 0:
            continue
        ax.scatter(d["x_plot"][m], d["y_plot"][m],
                   s=(sizes or {}).get(p, 5),
                   color=PHENOTYPE_COLORS.get(p, "#000000"),
                   linewidths=0, rasterized=True)
    ax.set_title(f"{short_label(sid)}  ({d['condition']})\n{d['n_cells']:,} cells",
                 fontsize=FONT_SIZE_TITLE - 12)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_edgecolor(CONDITION_COLORS.get(d["condition"], "#999999"))
        sp.set_linewidth(6)


def spatial_figure(show, stem, title, sizes=None, bg=True):
    n = len(SAMPLE_ORDER)
    ncol = 3
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(15 * ncol, 13 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for k, s in enumerate(SAMPLE_ORDER):
        spatial_panel(axes[k], s, show, sizes=sizes, bg=bg)
    for k in range(n, len(axes)):
        axes[k].axis("off")
    handles = [Patch(facecolor=PHENOTYPE_COLORS.get(p, "#000"), label=p)
               for p in show]
    if bg:
        handles.append(Patch(facecolor=BACKGROUND_COLOR, label="all other cells"))
    fig.legend(handles=handles, loc="lower center",
               ncol=min(5, len(handles)), frameon=False,
               fontsize=FONT_SIZE_LEGEND - 10, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(title, y=1.02, fontsize=FONT_SIZE_TITLE - 4)
    save_fig(fig, stem)


# all twelve classes, nothing greyed out
spatial_figure(PHENO_USE, "F20_spatial_all_phenotypes",
               "All phenotypes, six retained sections\n"
               f"(downsampled to at most {MAX_POINTS_PER_SECTION:,} cells "
               "per section)",
               sizes={p: 4 for p in PHENO_USE}, bg=False)

# tissue architecture: endothelium, epithelium, SMA-high stroma
spatial_figure([p for p in STRUCTURAL_PHENOTYPES if p in PHENO_USE],
               "F21_spatial_structural",
               "Tissue architecture: endothelium, epithelium and SMA-high stroma\n"
               'The "Other" class is SMA-positive and dim on every other marker',
               sizes={"Endothelial cells": 4, "Epithelial/Tumor cells": 4,
                      "Other": 4})

# immune populations against a grey structural background
spatial_figure([p for p in IMMUNE_PHENOTYPES if p in PHENO_USE],
               "F22_spatial_immune",
               "Immune populations against structural background\n"
               "Grey = endothelium, epithelium and stroma",
               sizes={"CD68+IDO1+ Macrophages": 10, "CD68+IDO1- Macrophages": 7,
                      "CD163+ Macrophages": 10, "Neutrophils": 4,
                      "Helper T cells": 12, "CD4- T cells": 5, "Tregs": 22,
                      "B cells": 8, "Plasma cells": 5})

# lymphoid structure candidates: B cells, plasma cells, helper T, Tregs only
spatial_figure([p for p in ["B cells", "Plasma cells", "Helper T cells", "Tregs"]
                if p in PHENO_USE],
               "F23_spatial_lymphoid",
               "Lymphoid aggregate candidates\n"
               "Discrete B cell clusters with adjacent T cells would indicate "
               "BALT-like structures",
               sizes={"B cells": 12, "Plasma cells": 6, "Helper T cells": 14,
                      "Tregs": 26})


# %% Cell 9 - wrap up
# =============================================================================

banner("SUMMARY")
print(f"Sections retained : {len(SAMPLE_ORDER)}  ({', '.join(SAMPLE_ORDER)})")
print(f"Sections excluded : {', '.join(EXCLUDE_SECTIONS)}")
print(f"Cells retained    : {int(counts.sum().sum()):,}")

sub("Where to read each answer")
print("  A composition + IDO1 fraction -> F17, tables 20 / 21 / 22")
print("  B threshold back-out          -> table 23  (READ THE VERDICT ABOVE)")
print("  C threshold sensitivity       -> F18 centre panel, table 24")
print("  D IDO1 distributions          -> F18 left and right panels")
print("  E boundary gradient           -> F19, table 25")
print("  F spatial maps                -> F20 all, F21 structural, F22 immune,")
print("                                   F23 lymphoid candidates")

sub("Decisions this run should let us close")
print("  1. Is the vendor IDO1 call reproducible as a simple intensity gate?")
print("  2. How wide is the cutoff window that leaves the result unchanged?")
print("  3. Is there a depth-into-tissue gradient in the retained sections?")
print("  4. Are there discrete lymphoid aggregates worth defining as structures?")

banner("DONE")
sys.stdout = _tee.terminal
_tee.close()
