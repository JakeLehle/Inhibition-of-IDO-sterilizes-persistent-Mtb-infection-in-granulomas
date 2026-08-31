#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - CORRECTED DISTANCE STATISTICS AND FORMAL TREATMENT TESTS
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 06 of the AKOYA analysis series.

WHAT THIS FIXES FROM SCRIPT 05
    1. EFFECT SIZE, NOT z, IS PRIMARY. Script 05 reported z against a
       within-structure permutation null. Because the null preserves cell
       counts, its width shrinks with structure size: treated structures hold
       1,328 to 2,337 cells while untreated hold 32,707 to 99,334. A z of -59
       in 36463 reflects a very tight null, not a large biological effect.
       Reporting now leads with delta (observed minus null mean) in radial units
       or microns, which is comparable across arms. z is retained only as a
       per-structure significance flag.

    2. NO CENSORING IN NEAREST-NEIGHBOUR DISTANCES. Script 05 used a k-neighbour
       index for permutation speed; 238 of 408 rows exceeded 20% censoring and
       some reached 90%, which biases the reported median downward by
       construction. This script builds a KD-tree on the TARGET cells only and
       queries anchors against it, giving the exact nearest-target distance.
       Permutations rebuild a small target tree each time, which costs about a
       millisecond. The speed problem the k-index solved did not exist.

    3. POOLED PER-ANIMAL ANALYSIS. With one focus per treated animal, a "median
       across structures" has no within-animal variance. Cells are now pooled
       across all foci within an animal (distances still computed within each
       cell's own focus, then pooled), giving one properly supported measurement
       per animal.

    4. TREGS RECOVERED. Script 05's tier gate dropped every Treg row from the
       animal summaries, so Q4 went unanswered. Tiers are now reported as a
       column and never used to filter. Pooling within animal also raises Treg
       counts into a usable range.

    5. DEGENERATE POLARISATION COLUMNS REMOVED. A within-section percentile makes
       the iNOS-high and Arg1-high sets exactly equal in size, so pct_inos_only
       and pct_arg1_only are identical by construction and carry no information.
       Reporting now uses double-positive against the independence expectation,
       plus spatial mixing, with no tier filtering so treated structures appear.

    6. FORMAL TREATMENT TEST BY MIXED MODEL. A cell-level test comparing 5,206
       treated cells against 187,056 untreated cells is pseudoreplication: cells
       within an animal are not independent, and such a test answers "are these
       two piles of cells different" rather than "does D1MT change this". A
       linear mixed model with arm as a fixed effect, animal as a random
       intercept and structure nested within animal uses every cell while
       keeping the effective sample size at the level of the six animals.
       Expect p in the 0.05 to 0.2 range even for a strong effect. That is
       correct, not a failure.

    7. 3-HK DYNAMIC RANGE GATE. The script 05 preview showed treated foci higher
       in BOTH IDO1 (115 vs 54) and 3-HK (0.243 vs 0.105). Both 3-HK values sit
       near zero on a channel where IDO1 medians are 54 to 115. Before any
       gradient analysis is built on 3-HK, this script establishes whether the
       channel has usable dynamic range at all. If it does not, Q3 is not
       answerable with this panel and we need to know now.

PREREQUISITE
    Rerun script 04 with FOCUS_FOLD_OVER_BACKGROUND = 6.0 first, which recovers
    43118's second focus (1,925, 6.8x background) and gives the treated arm more
    than one structure per animal. This script warns if it sees the old output.

OUTPUTS
    figures/  F43 .. F49
    tables/   50 .. 58

USAGE
    conda activate sc_pre
    python AKOYA_06_Distance_Statistics.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================

IN_DIR = "/master/jlehle/WORKING/AKOYA/structures_rev4"
CELL_DIR = f"{IN_DIR}/cell_assignments"
OUT_DIR = "/master/jlehle/WORKING/AKOYA/distance_stats"

CONDITION_ORDER = ["D1MT", "Untreated"]
CONDITION_COLORS = {"D1MT": "#2C7FB8", "Untreated": "#D95F02"}
REFERENCE_ARM = "Untreated"          # model coefficient is D1MT minus this

STRUCTURE_SOURCE = "foci"            # "foci" or "burden"
STRUCT_COL = "focus_id" if STRUCTURE_SOURCE == "foci" else "burden_region_id"

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

# populations carried into the formal models
KEY_PHENOTYPES = [IDO1_POS, IDO1_NEG, "Helper T cells", "CD4- T cells",
                  "Tregs", "B cells", "Plasma cells", "Neutrophils"]

# ---- nearest-neighbour pairs ------------------------------------------------
NN_ANCHORS = [IDO1_POS, IDO1_NEG]
NN_TARGETS = ["Helper T cells", "CD4- T cells", "Tregs", "B cells",
              "Plasma cells"]

# ---- permutation ------------------------------------------------------------
N_PERMUTATIONS = 300
MAX_ANCHORS_PER_STRUCTURE = 2000
RANDOM_SEED = 0

# ---- evidence tiers (REPORTED ONLY, NEVER USED TO FILTER) -------------------
TIER_SOLID = 200
TIER_USABLE = 50
TIER_PROVISIONAL = 20

# ---- radial binning ---------------------------------------------------------
RADIAL_MIN, RADIAL_MAX = 0.0, 2.0
N_RADIAL_BINS = 20

# ---- polarisation -----------------------------------------------------------
POLARISATION_PERCENTILES = [60, 70, 75, 80, 90]
POLARISATION_PRIMARY = 75
MIXING_K = 10
MIXING_MIN_PER_CLASS = 5             # minimum cells per class to compute mixing

# ---- 3-HK dynamic range gate ------------------------------------------------
HK3_COL = "3-Hydroxykynurenine"
HK3_MIN_USABLE_P99 = 5.0             # p99 below this means no usable range
HK3_MIN_NONZERO_FRAC = 0.05          # fewer nonzero cells than this is unusable

# ---- mixed models -----------------------------------------------------------
RUN_MIXED_MODELS = True
MODEL_MAX_CELLS_PER_STRUCTURE = 3000   # subsample for tractability
MODEL_MIN_CELLS_PER_ANIMAL = 20        # animals below this are dropped from a fit

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
import warnings
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

try:
    import statsmodels.formula.api as smf
    HAVE_SM = True
except Exception as _e2:
    HAVE_SM = False
    _sm_err = _e2

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
    null = np.asarray(null, float)
    null = null[np.isfinite(null)]
    if not len(null) or not np.isfinite(obs):
        return np.nan, np.nan, np.nan
    mu, sd = float(np.mean(null)), float(np.std(null, ddof=1))
    z = (obs - mu) / sd if sd > 0 else np.nan
    n_ext = int(np.sum(np.abs(null - mu) >= abs(obs - mu)))
    p = (n_ext + 1) / (len(null) + 1)
    return float(min(p, 1.0)), z, mu


def nn_distances(anchor_xy, target_xy):
    """
    Exact nearest-target distance for every anchor. No k-neighbour truncation,
    therefore no censoring. Anchors that coincide with a target return 0.
    """
    if len(anchor_xy) == 0 or len(target_xy) == 0:
        return np.array([])
    tree = cKDTree(target_xy)
    d, _ = tree.query(anchor_xy, k=1)
    return np.asarray(d, dtype=float)


def fit_mixed(df, outcome, label, extra_note=""):
    """
    Linear mixed model: outcome ~ arm, random intercept for animal, structure
    nested within animal. Falls back to animal-only random effects if the
    nested fit fails to converge.
    """
    if not HAVE_SM:
        return None
    d = df.dropna(subset=[outcome, "condition", "sample_id"]).copy()
    if not len(d):
        return None
    counts = d.groupby("sample_id").size()
    keep = counts.loc[counts >= MODEL_MIN_CELLS_PER_ANIMAL].index
    d = d.loc[d["sample_id"].isin(keep)]
    if d["condition"].nunique() < 2 or d["sample_id"].nunique() < 3:
        return None
    d["arm"] = (d["condition"] != REFERENCE_ARM).astype(float)
    d["struct_key"] = (d["sample_id"].astype(str) + "_"
                       + d[STRUCT_COL].astype(int).astype(str))
    d["_y"] = pd.to_numeric(d[outcome], errors="coerce")
    d = d.dropna(subset=["_y"])
    if not len(d):
        return None

    res, mode = None, ""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            md = smf.mixedlm("_y ~ arm", data=d, groups=d["sample_id"],
                             re_formula="1",
                             vc_formula={"struct": "0 + C(struct_key)"})
            res = md.fit(reml=True, method="lbfgs", maxiter=200)
            mode = "animal + structure"
        except Exception:
            res = None
        if res is None or not np.isfinite(res.params.get("arm", np.nan)):
            try:
                md = smf.mixedlm("_y ~ arm", data=d, groups=d["sample_id"],
                                 re_formula="1")
                res = md.fit(reml=True, method="lbfgs", maxiter=200)
                mode = "animal only"
            except Exception:
                return None
    if res is None:
        return None

    means = d.groupby("condition")["_y"].mean()
    return {
        "analysis": label, "outcome": outcome, "random_effects": mode,
        "n_animals_D1MT": int(d.loc[d["condition"] == "D1MT", "sample_id"].nunique()),
        "n_animals_ref": int(d.loc[d["condition"] == REFERENCE_ARM, "sample_id"].nunique()),
        "n_cells": int(len(d)), "n_animals": int(d["sample_id"].nunique()),
        "n_structures": int(d["struct_key"].nunique()),
        f"mean_{CONDITION_ORDER[0]}": float(means.get(CONDITION_ORDER[0], np.nan)),
        f"mean_{CONDITION_ORDER[1]}": float(means.get(CONDITION_ORDER[1], np.nan)),
        "coef_D1MT_vs_ref": float(res.params.get("arm", np.nan)),
        "std_err": float(res.bse.get("arm", np.nan)),
        "z_stat": float(res.tvalues.get("arm", np.nan)),
        "p_value": float(res.pvalues.get("arm", np.nan)),
        "note": extra_note,
    }


_tee = Tee(os.path.join(TAB_DIR, "00_distance_stats_report.txt"))
sys.stdout = _tee

banner("AKOYA CORRECTED DISTANCE STATISTICS AND FORMAL TREATMENT TESTS")
print(f"Run time      : {datetime.now().isoformat(timespec='seconds')}")
print(f"Permutations  : {N_PERMUTATIONS}")
print(f"Output        : {OUT_DIR}")
if not HAVE_SCIPY:
    print(f"\nERROR: scipy required ({_scipy_err}).")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)
if not HAVE_SM:
    print(f"\n    WARNING: statsmodels unavailable ({_sm_err}).")
    print("    Mixed models will be skipped. Permutation and per-animal effect")
    print("    sizes still run. Install with: conda install statsmodels")

print("\nREPORTING RULE")
print("    delta (observed minus null mean) is the effect size and is primary.")
print("    z depends on the width of the null, which shrinks with structure")
print("    size, so z is NOT comparable between arms and is reported only as a")
print("    per-structure significance flag.")


# %% Cell 3 - load and prerequisite check
# =============================================================================

banner("LOADING")

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
        print(f"    ERROR reading {sid}: {e}. Skipping.")
        continue
    if STRUCT_COL not in d.columns:
        print(f"    ERROR: {sid} lacks '{STRUCT_COL}'. Skipping.")
        continue
    d["sample_id"] = sid
    cells[sid] = d
    n_struct = int(d.loc[d[STRUCT_COL] > 0, STRUCT_COL].nunique())
    print(f"    {sid:<12} {d['condition'].iloc[0]:<10} {len(d):>9,} cells, "
          f"{int((d[STRUCT_COL] > 0).sum()):>8,} in {n_struct} structures")

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

sub("PREREQUISITE CHECK")
treated_structs = {s: int(cells[s].loc[cells[s][STRUCT_COL] > 0, STRUCT_COL].nunique())
                   for s in SAMPLE_ORDER if COND_OF[s] == "D1MT"}
if treated_structs and max(treated_structs.values()) <= 1:
    print("    WARNING: every treated animal has at most ONE structure.")
    print("    Rerun script 04 with FOCUS_FOLD_OVER_BACKGROUND = 6.0 to recover")
    print("    43118's second focus (1,925, 6.8x). Per-structure summaries in")
    print("    the treated arm are single values with no within-animal variance.")
    print("    Pooled per-animal results below are still valid; per-structure")
    print("    spreads are not.")
else:
    print(f"    treated structures per animal: {treated_structs}")


def structure_frames():
    for s in SAMPLE_ORDER:
        d = cells[s]
        for k, g in d.loc[d[STRUCT_COL] > 0].groupby(STRUCT_COL):
            yield s, int(k), g


# %% Cell 4 - 3-HK dynamic range gate
# =============================================================================

banner("3-HYDROXYKYNURENINE DYNAMIC RANGE GATE")

print("    The script 05 preview showed treated foci higher in BOTH IDO1 and")
print("    3-HK, with 3-HK medians of 0.24 and 0.10 on a channel where IDO1")
print("    medians are 54 to 115. Before any gradient analysis is built on")
print("    3-HK, establish whether the channel has usable range.\n")

hk_rows = []
for s in SAMPLE_ORDER:
    d = cells[s]
    if HK3_COL not in d.columns:
        print(f"    WARNING: {s} lacks '{HK3_COL}'")
        continue
    for scope, sel in [("all_cells", d),
                       ("macrophages", d.loc[d["pheno"].isin(MACROPHAGE_ALL)]),
                       ("ido1_pos_macs", d.loc[d["pheno"] == IDO1_POS])]:
        v = pd.to_numeric(sel[HK3_COL], errors="coerce").dropna().to_numpy()
        if not len(v):
            continue
        hk_rows.append({
            "sample_id": s, "condition": COND_OF[s], "scope": scope,
            "n_cells": len(v),
            "frac_zero": float(np.mean(v == 0)),
            "frac_nonzero": float(np.mean(v > 0)),
            "median": float(np.median(v)),
            "p75": float(np.percentile(v, 75)),
            "p90": float(np.percentile(v, 90)),
            "p99": float(np.percentile(v, 99)),
            "max": float(np.max(v)),
            "median_nonzero": float(np.median(v[v > 0])) if (v > 0).any() else np.nan,
        })
hk = pd.DataFrame(hk_rows)
if len(hk):
    write_csv(hk, "50_hk3_dynamic_range.csv")
    for scope in ["all_cells", "macrophages", "ido1_pos_macs"]:
        sub(f"3-HK, scope = {scope}")
        h = hk.loc[hk["scope"] == scope]
        print(h[["sample_id", "condition", "n_cells", "frac_nonzero", "median",
                 "p90", "p99", "max", "median_nonzero"]].to_string(index=False))

    sub("VERDICT ON Q3")
    hm = hk.loc[hk["scope"] == "macrophages"]
    p99_min = float(hm["p99"].min())
    nz_min = float(hm["frac_nonzero"].min())
    print(f"    minimum p99 across sections     : {p99_min:.3f}")
    print(f"    minimum nonzero fraction        : {nz_min:.3f}")
    if p99_min < HK3_MIN_USABLE_P99 or nz_min < HK3_MIN_NONZERO_FRAC:
        print("\n    3-HK HAS NO USABLE DYNAMIC RANGE in this panel.")
        print("    The apparent arm difference is a ratio of two near-zero")
        print("    numbers and must not be interpreted. Q3 (the 3-HK halo")
        print("    around IDO1+ cells) is NOT answerable with this data.")
        print("    Options: ask Akoya whether the 3-HK antibody passed QC on")
        print("    this run, or drop Q3 and rely on IDO1 protein alone.")
        HK3_USABLE = False
    else:
        print("\n    3-HK has usable range. The correct measurement is the")
        print("    FRACTION of macrophages that are 3-HK positive plus the")
        print("    intensity among those, not the median over all macrophages.")
        HK3_USABLE = True
else:
    HK3_USABLE = False
    print("    WARNING: no 3-HK data found.")


# %% Cell 5 - radial position, delta primary, per-structure and pooled
# =============================================================================

banner("RADIAL POSITION - EFFECT SIZES")

print("    delta = observed mean radial position minus the permutation null")
print("    mean, in radial units (0 = core centre, 1 = core boundary, 2 = outer")
print("    cuff). Negative = pulled toward the core.\n")

rad_struct, rad_pool = [], []

# ---- per structure ----------------------------------------------------------
for s, k, g in structure_frames():
    r = pd.to_numeric(g["radial_pos"], errors="coerce").to_numpy()
    ok = np.isfinite(r)
    r, ph = r[ok], g["pheno"].to_numpy()[ok]
    n = len(r)
    if n < 10:
        continue
    for p in PHENOTYPE_ORDER:
        npos = int((ph == p).sum())
        if npos == 0:
            continue
        obs = float(r[ph == p].mean())
        null = np.array([r[rng.choice(n, size=npos, replace=False)].mean()
                         for _ in range(N_PERMUTATIONS)])
        pval, z, mu = empirical_p(obs, null)
        rad_struct.append({
            "sample_id": s, "condition": COND_OF[s], "structure_id": k,
            "phenotype": p, "n_cells": npos, "n_structure": n,
            "observed_radial": obs, "null_radial": mu,
            "delta_radial": obs - mu, "z": z, "p_empirical": pval,
            "tier": tier_of(npos),
        })

# ---- pooled per animal ------------------------------------------------------
for s in SAMPLE_ORDER:
    d = cells[s]
    sel = d.loc[(d[STRUCT_COL] > 0) & np.isfinite(
        pd.to_numeric(d["radial_pos"], errors="coerce"))]
    if not len(sel):
        continue
    r = pd.to_numeric(sel["radial_pos"], errors="coerce").to_numpy()
    ph = sel["pheno"].to_numpy()
    struct = sel[STRUCT_COL].to_numpy()
    n = len(r)
    for p in PHENOTYPE_ORDER:
        m = ph == p
        npos = int(m.sum())
        if npos == 0:
            continue
        obs = float(r[m].mean())
        # permute labels WITHIN each structure, then pool, so structure
        # membership is preserved in the null
        null = np.empty(N_PERMUTATIONS)
        for it in range(N_PERMUTATIONS):
            fake = np.zeros(n, dtype=bool)
            for kk in np.unique(struct):
                idx = np.flatnonzero(struct == kk)
                want = int((m & (struct == kk)).sum())
                if want:
                    fake[rng.choice(idx, size=want, replace=False)] = True
            null[it] = r[fake].mean()
        pval, z, mu = empirical_p(obs, null)
        rad_pool.append({
            "sample_id": s, "animal_id": short_label(s), "condition": COND_OF[s],
            "phenotype": p, "n_cells": npos, "n_pooled": n,
            "observed_radial": obs, "null_radial": mu,
            "delta_radial": obs - mu, "z": z, "p_empirical": pval,
            "tier": tier_of(npos),
        })
    print(f"    {s} pooled")

rad_s = pd.DataFrame(rad_struct)
rad_p = pd.DataFrame(rad_pool)
if len(rad_s):
    write_csv(rad_s, "51_radial_per_structure.csv")
if len(rad_p):
    write_csv(rad_p, "52_radial_pooled_per_animal.csv")

    sub("Pooled per-animal delta radial position (NO tier filtering)")
    print(f"    {'phenotype':<26}" + "".join(f"{short_label(s):>12}" for s in SAMPLE_ORDER))
    print("    " + "-" * (26 + 12 * len(SAMPLE_ORDER)))
    for p in KEY_PHENOTYPES:
        row = f"    {p:<26}"
        for s in SAMPLE_ORDER:
            v = rad_p.loc[(rad_p["sample_id"] == s) & (rad_p["phenotype"] == p)]
            row += (f"{v['delta_radial'].iloc[0]:>+12.3f}" if len(v) else f"{'na':>12}")
        print(row)
    print(f"\n    {'(n cells)':<26}" + "".join(f"{short_label(s):>12}" for s in SAMPLE_ORDER))
    for p in KEY_PHENOTYPES:
        row = f"    {p:<26}"
        for s in SAMPLE_ORDER:
            v = rad_p.loc[(rad_p["sample_id"] == s) & (rad_p["phenotype"] == p)]
            row += (f"{int(v['n_cells'].iloc[0]):>12,}" if len(v) else f"{'na':>12}")
        print(row)


# %% Cell 6 - nearest-neighbour distances, exact, per structure and pooled
# =============================================================================

banner("NEAREST-NEIGHBOUR DISTANCES - EXACT, NO CENSORING")

print("    KD-tree built on target cells only; every anchor gets its exact")
print("    nearest-target distance. Null: reassign anchor and target labels")
print("    within the structure, preserving counts and coordinates.\n")

nn_struct, nn_pool, nn_cells_rows = [], [], []

for s in SAMPLE_ORDER:
    d = cells[s]
    sel = d.loc[d[STRUCT_COL] > 0]
    if not len(sel):
        continue
    per_pair_obs = {}
    per_pair_null = {}

    for k, g in sel.groupby(STRUCT_COL):
        xy = g[["x", "y"]].to_numpy(float)
        ph = g["pheno"].to_numpy()
        n = len(xy)
        if n < 10:
            continue
        for anchor in NN_ANCHORS:
            a_idx = np.flatnonzero(ph == anchor)
            n_a = len(a_idx)
            if n_a == 0:
                continue
            a_use = (rng.choice(a_idx, size=MAX_ANCHORS_PER_STRUCTURE,
                                replace=False)
                     if n_a > MAX_ANCHORS_PER_STRUCTURE else a_idx)
            for target in NN_TARGETS:
                if target == anchor:
                    continue
                t_idx = np.flatnonzero(ph == target)
                n_t = len(t_idx)
                if n_t == 0:
                    continue
                obs_v = nn_distances(xy[a_use], xy[t_idx])
                obs = float(np.median(obs_v)) if len(obs_v) else np.nan

                null = np.empty(N_PERMUTATIONS)
                for it in range(N_PERMUTATIONS):
                    perm = rng.permutation(n)
                    fa = perm[:n_a]
                    ft = perm[n_a:n_a + n_t]
                    fa_use = (rng.choice(fa, size=MAX_ANCHORS_PER_STRUCTURE,
                                         replace=False)
                              if len(fa) > MAX_ANCHORS_PER_STRUCTURE else fa)
                    v = nn_distances(xy[fa_use], xy[ft])
                    null[it] = np.median(v) if len(v) else np.nan
                pval, z, mu = empirical_p(obs, null)

                nn_struct.append({
                    "sample_id": s, "condition": COND_OF[s], "structure_id": int(k),
                    "anchor": anchor, "target": target,
                    "n_anchor": n_a, "n_target": n_t, "n_structure": n,
                    "observed_median_um": obs, "null_median_um": mu,
                    "delta_um": obs - mu if np.isfinite(mu) else np.nan,
                    "z": z, "p_empirical": pval,
                    "tier": tier_of(min(n_a, n_t)),
                })

                key = (anchor, target)
                per_pair_obs.setdefault(key, []).append(obs_v)
                per_pair_null.setdefault(key, []).append(null)

                # per-cell rows for the mixed model
                take = obs_v
                if len(take) > MODEL_MAX_CELLS_PER_STRUCTURE:
                    take = rng.choice(take, size=MODEL_MAX_CELLS_PER_STRUCTURE,
                                      replace=False)
                nn_cells_rows.append(pd.DataFrame({
                    "sample_id": s, "condition": COND_OF[s],
                    STRUCT_COL: int(k), "anchor": anchor, "target": target,
                    "nn_distance_um": take,
                }))

    # pooled per animal: pool per-cell distances across structures
    for key, arrs in per_pair_obs.items():
        allv = np.concatenate([a for a in arrs if len(a)])
        if not len(allv):
            continue
        nulls = np.vstack(per_pair_null[key])          # (n_struct, n_perm)
        pooled_null = np.nanmean(nulls, axis=0)        # per permutation
        obs = float(np.median(allv))
        pval, z, mu = empirical_p(obs, pooled_null)
        nn_pool.append({
            "sample_id": s, "animal_id": short_label(s), "condition": COND_OF[s],
            "anchor": key[0], "target": key[1],
            "n_anchor_cells": int(len(allv)),
            "n_structures": int(len(arrs)),
            "observed_median_um": obs, "null_median_um": mu,
            "delta_um": obs - mu if np.isfinite(mu) else np.nan,
            "z": z, "p_empirical": pval,
            "tier": tier_of(len(allv)),
        })
    print(f"    {s} done")
    gc.collect()

nn_s = pd.DataFrame(nn_struct)
nn_p = pd.DataFrame(nn_pool)
nn_cells = (pd.concat(nn_cells_rows, ignore_index=True)
            if nn_cells_rows else pd.DataFrame())

if len(nn_s):
    write_csv(nn_s, "53_nn_per_structure.csv")
if len(nn_p):
    write_csv(nn_p, "54_nn_pooled_per_animal.csv")
    sub("Pooled per-animal nearest-neighbour distances (NO tier filtering)")
    for anchor in NN_ANCHORS:
        print(f"\n    anchor: {anchor}")
        print(f"      {'target':<18}" +
              "".join(f"{short_label(s):>14}" for s in SAMPLE_ORDER))
        for target in NN_TARGETS:
            if target == anchor:
                continue
            row = f"      {target:<18}"
            for s in SAMPLE_ORDER:
                v = nn_p.loc[(nn_p["sample_id"] == s) &
                             (nn_p["anchor"] == anchor) &
                             (nn_p["target"] == target)]
                if len(v):
                    row += (f"{v['observed_median_um'].iloc[0]:>7.0f}"
                            f"/{v['delta_um'].iloc[0]:>+6.1f}")
                else:
                    row += f"{'na':>14}"
            print(row)
        print("      (observed median um / delta vs null)")


# %% Cell 7 - polarisation, corrected
# =============================================================================

banner("iNOS vs ARGINASE-1 - CORRECTED REPORTING")

print("    pct_inos_only and pct_arg1_only are omitted: a within-section")
print("    percentile makes the two sets exactly equal in size, so those")
print("    columns are identical by construction and carry no information.")
print("    Reported instead: double-positive against the independence")
print("    expectation, and spatial mixing with NO tier filtering.\n")

pol_rows, mix_rows = [], []
for s in SAMPLE_ORDER:
    d = cells[s]
    mac = d.loc[d["pheno"].isin(MACROPHAGE_ALL)]
    if not len(mac):
        continue
    for pct in POLARISATION_PERCENTILES:
        it_ = float(np.nanpercentile(mac["iNOS"], pct))
        at_ = float(np.nanpercentile(mac["Arginase-1"], pct))
        hi_i = (mac["iNOS"] >= it_).to_numpy()
        hi_a = (mac["Arginase-1"] >= at_).to_numpy()
        p_i, p_a = float(hi_i.mean()), float(hi_a.mean())
        obs_dbl = float((hi_i & hi_a).mean())
        exp_dbl = p_i * p_a
        pol_rows.append({
            "sample_id": s, "condition": COND_OF[s], "percentile": pct,
            "n_macrophages": len(mac),
            "pct_double_positive": 100.0 * obs_dbl,
            "pct_expected_if_independent": 100.0 * exp_dbl,
            "double_positive_ratio": obs_dbl / exp_dbl if exp_dbl > 0 else np.nan,
            "pct_double_negative": 100.0 * float((~hi_i & ~hi_a).mean()),
        })

    it_ = float(np.nanpercentile(mac["iNOS"], POLARISATION_PRIMARY))
    at_ = float(np.nanpercentile(mac["Arginase-1"], POLARISATION_PRIMARY))
    for k, g in d.loc[(d[STRUCT_COL] > 0) &
                      d["pheno"].isin(MACROPHAGE_ALL)].groupby(STRUCT_COL):
        n = len(g)
        if n < 20:
            continue
        hi_i = (g["iNOS"] >= it_).to_numpy()
        hi_a = (g["Arginase-1"] >= at_).to_numpy()
        ex_i, ex_a = hi_i & ~hi_a, hi_a & ~hi_i
        n_i, n_a = int(ex_i.sum()), int(ex_a.sum())
        if n_i < MIXING_MIN_PER_CLASS or n_a < MIXING_MIN_PER_CLASS:
            continue
        xy = g[["x", "y"]].to_numpy(float)
        kk = int(min(MIXING_K, n - 1))
        tree = cKDTree(xy)
        _, idxs = tree.query(xy, k=kk + 1)
        idxs = idxs[:, 1:]

        def mixing(a_mask, b_mask):
            rows = np.flatnonzero(a_mask)
            return float(b_mask[idxs[rows]].mean()) if len(rows) else np.nan

        obs = mixing(ex_i, ex_a)
        null = np.empty(N_PERMUTATIONS)
        for it2 in range(N_PERMUTATIONS):
            perm = rng.permutation(n)
            fi = np.zeros(n, bool); fi[perm[:n_i]] = True
            fa = np.zeros(n, bool); fa[perm[n_i:n_i + n_a]] = True
            null[it2] = mixing(fi, fa)
        pval, z, mu = empirical_p(obs, null)
        mix_rows.append({
            "sample_id": s, "condition": COND_OF[s], "structure_id": int(k),
            "n_macrophages": n, "n_inos_only": n_i, "n_arg1_only": n_a,
            "observed_mixing": obs, "null_mixing": mu,
            "delta_mixing": obs - mu if np.isfinite(mu) else np.nan,
            "z": z, "p_empirical": pval, "tier": tier_of(min(n_i, n_a)),
        })
        del tree, idxs

pol = pd.DataFrame(pol_rows)
mix = pd.DataFrame(mix_rows)
if len(pol):
    write_csv(pol, "55_polarisation_composition.csv")
    sub(f"Macrophage co-expression at the {POLARISATION_PRIMARY}th percentile")
    print(pol.loc[pol["percentile"] == POLARISATION_PRIMARY][
        ["sample_id", "condition", "n_macrophages", "pct_double_positive",
         "pct_expected_if_independent", "double_positive_ratio",
         "pct_double_negative"]].to_string(index=False))
    print("\n    A ratio above 1 means iNOS and Arginase-1 are CO-EXPRESSED more")
    print("    than chance, which runs against a strict M1/M2 dichotomy.")
if len(mix):
    write_csv(mix, "56_inos_arg1_mixing.csv")
    sub("Spatial mixing, all structures, no tier filtering")
    mix_anim = (mix.groupby(["sample_id", "condition"])
                .agg(median_delta=("delta_mixing", "median"),
                     median_z=("z", "median"),
                     median_observed=("observed_mixing", "median"),
                     n_structures=("z", "size")).reset_index())
    print(mix_anim.to_string(index=False))
    write_csv(mix_anim, "56b_inos_arg1_mixing_by_animal.csv")
    missing = [s for s in SAMPLE_ORDER if s not in set(mix["sample_id"])]
    if missing:
        print(f"\n    No mixing rows for: {missing}")
        print("    Likely too few exclusively-polarised macrophages inside")
        print("    structures. Check n_inos_only / n_arg1_only in table 56.")


# %% Cell 8 - formal treatment tests by mixed model
# =============================================================================

banner("FORMAL TREATMENT TESTS - LINEAR MIXED MODELS")

print("    Cells within an animal are not independent. A plain cell-level test")
print("    across 5,206 treated and 187,056 untreated cells is")
print("    pseudoreplication and would return an implausibly small p-value for")
print("    an effect of no consequence. These models use every cell but hold")
print("    the effective sample size at the level of the six animals.")
print("    With three animals per arm, p in the 0.05 to 0.2 range is the")
print("    expected result for a real effect, not a failure.\n")

model_rows = []
if RUN_MIXED_MODELS and HAVE_SM:
    # radial position per phenotype
    for p in KEY_PHENOTYPES:
        frames = []
        for s in SAMPLE_ORDER:
            d = cells[s]
            sel = d.loc[(d[STRUCT_COL] > 0) & (d["pheno"] == p)].copy()
            sel = sel.loc[np.isfinite(pd.to_numeric(sel["radial_pos"],
                                                    errors="coerce"))]
            if not len(sel):
                continue
            out = []
            for k, g in sel.groupby(STRUCT_COL):
                out.append(g.sample(min(len(g), MODEL_MAX_CELLS_PER_STRUCTURE),
                                    random_state=RANDOM_SEED))
            if out:
                frames.append(pd.concat(out, ignore_index=True))
        if not frames:
            continue
        dd = pd.concat(frames, ignore_index=True)
        r = fit_mixed(dd, "radial_pos", f"radial position: {p}")
        if r:
            r["phenotype"] = p
            model_rows.append(r)
            print(f"    radial {p:<26} coef={r['coef_D1MT_vs_ref']:>+7.3f}  "
                  f"p={r['p_value']:.4f}  cells={r['n_cells']:,}  "
                  f"[{r['random_effects']}]")

    # nearest-neighbour distance per pair
    if len(nn_cells):
        for (a, t), g in nn_cells.groupby(["anchor", "target"]):
            r = fit_mixed(g, "nn_distance_um", f"nn distance: {a} -> {t}")
            if r:
                r["anchor"], r["target"] = a, t
                model_rows.append(r)
                print(f"    nn {a[:18]:<20} -> {t:<16} "
                      f"coef={r['coef_D1MT_vs_ref']:>+7.2f} um  "
                      f"p={r['p_value']:.4f}  cells={r['n_cells']:,}")
elif not HAVE_SM:
    print("    SKIPPED: statsmodels unavailable.")
else:
    print("    SKIPPED: RUN_MIXED_MODELS = False")

models = pd.DataFrame(model_rows)
if len(models):
    write_csv(models, "57_mixed_model_results.csv")
    sub("Interpretation")
    print("    coef is the D1MT effect relative to Untreated, in the units of")
    print("    the outcome (radial units, or microns). A negative radial coef")
    print("    means that population sits closer to the core in treated foci.")
    print("    A negative distance coef means the pair sits closer together in")
    print("    treated foci.")


# %% Cell 9 - figures
# =============================================================================

banner("FIGURES")

# ---- F43 pooled radial delta ------------------------------------------------
if len(rad_p):
    fig, axes = plt.subplots(1, 2, figsize=(32, 14),
                             gridspec_kw={"width_ratios": [1.35, 1.0]})
    ax = axes[0]
    yy = np.arange(len(KEY_PHENOTYPES))
    for i, p in enumerate(KEY_PHENOTYPES):
        for s in SAMPLE_ORDER:
            v = rad_p.loc[(rad_p["sample_id"] == s) & (rad_p["phenotype"] == p)]
            if not len(v):
                continue
            off = (SAMPLE_ORDER.index(s) - (len(SAMPLE_ORDER) - 1) / 2) * 0.13
            small = v["tier"].iloc[0] in ("provisional", "not_interpretable")
            ax.scatter(v["delta_radial"].iloc[0], i + off, s=420,
                       color=COLOR_OF[s], marker=MARKER_OF[s],
                       edgecolor=FLAG_COLOR if small else "#FFFFFF",
                       linewidth=3 if small else 2, zorder=3)
    ax.axvline(0, color="#000000", linewidth=3.5)
    ax.set_yticks(yy); ax.set_yticklabels(KEY_PHENOTYPES,
                                          fontsize=FONT_SIZE_TICK - 6)
    ax.invert_yaxis()
    ax.set_xlabel("delta radial position\n(negative = core, positive = rim)")
    ax.set_title("Pooled per animal", fontsize=FONT_SIZE_TITLE - 8)
    ax.xaxis.grid(True, color=GRID_COLOR, linewidth=1.5); ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    h = [Line2D([0], [0], color=COLOR_OF[s], marker=MARKER_OF[s], markersize=16,
                linestyle="none", label=f"{short_label(s)} ({COND_OF[s]})")
         for s in SAMPLE_ORDER]
    h.append(Line2D([0], [0], color="#FFFFFF", marker="o", markersize=16,
                    markeredgecolor=FLAG_COLOR, markeredgewidth=3,
                    linestyle="none", label="thin (<50 cells)"))
    ax.legend(handles=h, frameon=False, fontsize=FONT_SIZE_LEGEND - 14, loc="best")

    ax = axes[1]
    for i, p in enumerate(KEY_PHENOTYPES):
        for ci, c in enumerate(CONDITION_ORDER):
            v = rad_p.loc[(rad_p["condition"] == c) &
                          (rad_p["phenotype"] == p), "delta_radial"].dropna()
            if not len(v):
                continue
            off = (ci - 0.5) * 0.24
            ax.scatter(v, np.full(len(v), i + off), s=300,
                       color=CONDITION_COLORS[c], edgecolor="#FFFFFF",
                       linewidth=2, zorder=3, alpha=0.9)
            ax.hlines(i + off, v.min(), v.max(), color=CONDITION_COLORS[c],
                      linewidth=4, alpha=0.5, zorder=2)
    ax.axvline(0, color="#000000", linewidth=3.5)
    ax.set_yticks(yy); ax.set_yticklabels([]); ax.invert_yaxis()
    ax.set_xlabel("delta radial position by arm")
    ax.set_title("By arm", fontsize=FONT_SIZE_TITLE - 8)
    ax.xaxis.grid(True, color=GRID_COLOR, linewidth=1.5); ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.legend(handles=[Patch(facecolor=CONDITION_COLORS[c], label=c)
                       for c in CONDITION_ORDER],
              frameon=False, fontsize=FONT_SIZE_LEGEND - 8, loc="best")
    fig.suptitle("Radial position as an effect size, not a z score\n"
                 "delta is in radial units and is comparable between arms",
                 y=1.04, fontsize=FONT_SIZE_TITLE - 6)
    save_fig(fig, "F43_radial_delta_pooled")

# ---- F44 pooled nn distances ------------------------------------------------
if len(nn_p):
    pairs = [(a, t) for a in NN_ANCHORS for t in NN_TARGETS if a != t]
    ncol = 5
    nrow = int(np.ceil(len(pairs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(9 * ncol, 8.5 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for i, (a, t) in enumerate(pairs):
        ax = axes[i]
        d = nn_p.loc[(nn_p["anchor"] == a) & (nn_p["target"] == t)]
        for ci, c in enumerate(CONDITION_ORDER):
            v = d.loc[d["condition"] == c]
            if not len(v):
                continue
            j = rng.uniform(-0.12, 0.12, size=len(v))
            ax.scatter(np.full(len(v), ci) + j, v["observed_median_um"], s=380,
                       color=CONDITION_COLORS[c], edgecolor="#FFFFFF",
                       linewidth=2, zorder=3)
            ax.scatter(np.full(len(v), ci) + j, v["null_median_um"], s=200,
                       facecolor="none", edgecolor="#666666", linewidth=3,
                       zorder=3)
        ax.set_xticks(range(len(CONDITION_ORDER)))
        ax.set_xticklabels(CONDITION_ORDER, fontsize=FONT_SIZE_TICK - 12)
        ax.set_xlim(-0.6, len(CONDITION_ORDER) - 0.4)
        ax.set_title(f"{a[:22]}\nto {t}", fontsize=FONT_SIZE_BASE - 12)
        ax.set_ylabel("median distance (µm)", fontsize=FONT_SIZE_BASE - 14)
        ax.tick_params(labelsize=FONT_SIZE_TICK - 14)
        style_axes(ax)
    for i in range(len(pairs), len(axes)):
        axes[i].axis("off")
    fig.legend(handles=[Line2D([0], [0], marker="o", markersize=18,
                               linestyle="none", color="#444444",
                               label="observed"),
                        Line2D([0], [0], marker="o", markersize=14,
                               linestyle="none", markerfacecolor="none",
                               markeredgecolor="#666666", markeredgewidth=3,
                               label="permutation null")],
               loc="lower center", ncol=2, frameon=False,
               fontsize=FONT_SIZE_LEGEND - 8, bbox_to_anchor=(0.5, -0.03))
    fig.suptitle("Nearest-neighbour distances pooled per animal, exact and "
                 "uncensored\nOpen circles are the within-structure "
                 "permutation null", y=1.02, fontsize=FONT_SIZE_TITLE - 6)
    save_fig(fig, "F44_nn_pooled_distances")

# ---- F45 mixed model forest -------------------------------------------------
if len(models):
    m = models.copy()
    m["label"] = m["analysis"]
    m = m.sort_values("p_value")
    fig, ax = plt.subplots(figsize=(22, max(12, 0.75 * len(m))))
    yy = np.arange(len(m))
    lo = m["coef_D1MT_vs_ref"] - 1.96 * m["std_err"]
    hi = m["coef_D1MT_vs_ref"] + 1.96 * m["std_err"]
    for i, (_, r) in enumerate(m.iterrows()):
        sig = r["p_value"] < 0.05
        ax.plot([lo.iloc[i], hi.iloc[i]], [i, i],
                color=FLAG_COLOR if sig else "#666666", linewidth=5, zorder=2)
        ax.scatter([r["coef_D1MT_vs_ref"]], [i], s=420,
                   color=FLAG_COLOR if sig else "#666666",
                   edgecolor="#FFFFFF", linewidth=2, zorder=3)
        ax.text(hi.iloc[i], i, f"  p={r['p_value']:.3f}", va="center",
                fontsize=FONT_SIZE_ANNOT - 10)
    ax.axvline(0, color="#000000", linewidth=3.5)
    ax.set_yticks(yy)
    ax.set_yticklabels(m["label"], fontsize=FONT_SIZE_TICK - 12)
    ax.invert_yaxis()
    ax.set_xlabel("D1MT effect vs Untreated (95% CI)")
    ax.set_title("Mixed model treatment effects\n"
                 "animal as random intercept, structure nested within animal",
                 fontsize=FONT_SIZE_TITLE - 8)
    ax.xaxis.grid(True, color=GRID_COLOR, linewidth=1.5); ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    save_fig(fig, "F45_mixed_model_forest")

# ---- F46 3-HK dynamic range -------------------------------------------------
if len(hk):
    fig, axes = plt.subplots(1, 2, figsize=(30, 13))
    ax = axes[0]
    hm = hk.loc[hk["scope"] == "macrophages"]
    xp = np.arange(len(hm))
    for i, (_, r) in enumerate(hm.iterrows()):
        c = CONDITION_COLORS[r["condition"]]
        ax.plot([i, i], [r["median"], r["max"]], color=c, linewidth=5, alpha=0.5)
        ax.scatter([i], [r["median"]], s=380, color=c, marker="_", linewidth=6)
        ax.scatter([i], [r["p99"]], s=340, color=c, edgecolor="#FFFFFF",
                   linewidth=2, zorder=3)
        ax.scatter([i], [r["max"]], s=200, facecolor="none", edgecolor=c,
                   linewidth=3, zorder=3)
    ax.axhline(HK3_MIN_USABLE_P99, color=FLAG_COLOR, linestyle="--", linewidth=3)
    ax.text(len(hm) - 0.5, HK3_MIN_USABLE_P99 * 1.1, "usable p99 floor",
            ha="right", fontsize=FONT_SIZE_ANNOT - 10, color=FLAG_COLOR)
    ax.set_yscale("symlog", linthresh=0.1)
    ax.set_xticks(xp)
    ax.set_xticklabels([short_label(s) for s in hm["sample_id"]], rotation=45,
                       ha="right", fontsize=FONT_SIZE_TICK - 10)
    ax.set_ylabel("3-HK intensity in macrophages")
    ax.set_title("Median, p99 and max", fontsize=FONT_SIZE_TITLE - 12)
    style_axes(ax)

    ax = axes[1]
    ax.bar(xp, 100 * hm["frac_nonzero"].to_numpy(),
           color=[CONDITION_COLORS[c] for c in hm["condition"]],
           edgecolor="#FFFFFF", linewidth=2, zorder=3)
    ax.axhline(100 * HK3_MIN_NONZERO_FRAC, color=FLAG_COLOR, linestyle="--",
               linewidth=3)
    ax.set_xticks(xp)
    ax.set_xticklabels([short_label(s) for s in hm["sample_id"]], rotation=45,
                       ha="right", fontsize=FONT_SIZE_TICK - 10)
    ax.set_ylabel("% macrophages with nonzero 3-HK")
    ax.set_title("Detection rate", fontsize=FONT_SIZE_TITLE - 12)
    style_axes(ax)
    fig.suptitle("Does 3-hydroxykynurenine have usable dynamic range?\n"
                 "This decides whether Q3 is answerable with this panel",
                 y=1.04, fontsize=FONT_SIZE_TITLE - 8)
    save_fig(fig, "F46_hk3_dynamic_range")


# %% Cell 10 - consistency and wrap up
# =============================================================================

banner("CONSISTENCY ACROSS ANIMALS")

cons_rows = []


def consistency(df, group_cols, value_col, label):
    for keys, g in df.groupby(group_cols):
        keys = keys if isinstance(keys, tuple) else (keys,)
        rec = {"analysis": label, "value": value_col}
        for cname, kv in zip(group_cols, keys):
            rec[cname] = kv
        for c in CONDITION_ORDER:
            v = g.loc[g["condition"] == c, value_col].dropna()
            rec[f"{c}_n"] = len(v)
            rec[f"{c}_median"] = float(v.median()) if len(v) else np.nan
            rec[f"{c}_k_negative"] = int((v < 0).sum())
            rec[f"{c}_k_positive"] = int((v > 0).sum())
        a = rec.get(f"{CONDITION_ORDER[0]}_median", np.nan)
        b = rec.get(f"{CONDITION_ORDER[1]}_median", np.nan)
        rec["arm_difference"] = a - b if np.isfinite(a) and np.isfinite(b) else np.nan
        cons_rows.append(rec)


if len(rad_p):
    consistency(rad_p, ["phenotype"], "delta_radial", "radial_position")
if len(nn_p):
    consistency(nn_p, ["anchor", "target"], "delta_um", "nn_proximity")
if len(mix):
    m2 = (mix.groupby(["sample_id", "condition"])["delta_mixing"]
          .median().reset_index())
    m2["phenotype"] = "iNOS_vs_Arg1_mixing"
    consistency(m2, ["phenotype"], "delta_mixing", "polarisation_mixing")

cons = pd.DataFrame(cons_rows)
if len(cons):
    write_csv(cons, "58_consistency_summary.csv")
    sub("Effects separating the arms, ordered by arm difference")
    cc = cons.dropna(subset=["arm_difference"]).copy()
    cc["abs_diff"] = cc["arm_difference"].abs()
    for _, r in cc.sort_values("abs_diff", ascending=False).head(20).iterrows():
        name = " / ".join(str(r[c]) for c in ["phenotype", "anchor", "target"]
                          if c in r.index and pd.notna(r.get(c)))
        print(f"    [{r['analysis']:<20}] {name:<46} "
              f"D1MT {r.get('D1MT_median', np.nan):>+8.3f} "
              f"(n={int(r.get('D1MT_n', 0))})   "
              f"Untr {r.get('Untreated_median', np.nan):>+8.3f} "
              f"(n={int(r.get('Untreated_n', 0))})   "
              f"diff {r['arm_difference']:>+8.3f}")

banner("SUMMARY")
print(f"Radial rows (structure / pooled) : {len(rad_s)} / {len(rad_p)}")
print(f"NN rows (structure / pooled)     : {len(nn_s)} / {len(nn_p)}")
print(f"Mixing rows                      : {len(mix)}")
print(f"Mixed models fitted              : {len(models)}")
print(f"3-HK usable for Q3               : {HK3_USABLE}")

sub("Read in this order")
print("  1. F46 / table 50 : is 3-HK usable? Decides whether Q3 survives.")
print("  2. F43 / table 52 : radial position as delta. Q1 architecture.")
print("  3. F44 / table 54 : uncensored distances pooled per animal. Q1, Q4.")
print("  4. F45 / table 57 : mixed model treatment effects. The formal test.")
print("  5. Table 58       : per-animal consistency and arm differences.")

sub("What to tell Deepak about the statistics")
print("  Three levels, all reported:")
print("    - Within-animal permutation: high confidence, well powered, valid.")
print("      'In 43109, IDO1+ macrophages sit closer to CD4 T cells than that")
print("      lesion's own architecture predicts, p = 0.005, on 22,052 cells.'")
print("    - Per-animal effect sizes with k of 3 consistency: carries the")
print("      treatment claim honestly.")
print("    - Mixed model: the formal treatment test. Bounded by three animals")
print("      per arm. A p of 0.08 here is a real result, and any claim of")
print("      p < 0.001 for a treatment effect from three animals per group is")
print("      misusing the cell count.")

banner("DONE")
sys.stdout = _tee.terminal
_tee.close()
