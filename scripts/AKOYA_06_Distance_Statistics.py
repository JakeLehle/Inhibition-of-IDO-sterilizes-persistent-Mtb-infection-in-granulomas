#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - CORRECTED DISTANCE STATISTICS AND FORMAL TREATMENT TESTS
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 06 of the AKOYA analysis series. REVISION 4.

WHAT THIS FIXED FROM SCRIPT 05 (revision 1, unchanged and still correct)
    1. EFFECT SIZE, NOT z, IS PRIMARY. Script 05 reported z against a
       within-structure permutation null. Because the null preserves cell
       counts, its width shrinks with structure size, so a large z can reflect a
       tight null rather than a large effect. Reporting leads with delta,
       observed minus null mean, in radial units or microns. z is retained only
       as a per-structure significance flag.

    2. NO CENSORING IN NEAREST-NEIGHBOUR DISTANCES. Script 05 used a k-neighbour
       index; 271 of 463 rows exceeded 20 percent censoring and some reached 90
       percent, which biases the reported median downward by construction. This
       script builds a KD-tree on the TARGET cells only and queries anchors
       against it, giving the exact nearest-target distance.

    3. POOLED PER-ANIMAL ANALYSIS. With few foci per treated animal, a median
       across structures has almost no within-animal variance. Cells are pooled
       across all foci within an animal, distances still computed within each
       cell's own focus.

    4. TREGS RECOVERED. Tiers are reported as a column and never used to filter.

    5. DEGENERATE POLARISATION COLUMNS REMOVED. A within-section percentile
       makes the iNOS-high and Arg1-high sets exactly equal in size, so
       pct_inos_only and pct_arg1_only are identical by construction. Script 05
       revision 1 confirmed this: 16.48 against 16.48, 18.03 against 18.03, and
       so on. Reporting uses double-positive against the independence
       expectation plus spatial mixing.

    6. FORMAL TREATMENT TEST BY MIXED MODEL, holding the effective sample size
       at the six animals rather than at the cell count.

    7. 3-HK DYNAMIC RANGE GATE, before any gradient analysis is built on it.

WHAT CHANGED IN REVISION 2 (and why)

    1. THE POOLED PERMUTATION NULL WAS BUILT WRONG. THIS WAS THE MAIN FIX.
       Revision 1 built the pooled per-animal null as
           pooled_null = np.nanmean(np.vstack(per_structure_nulls), axis=0)
       that is, the mean across structures of each structure's own null median,
       one value per permutation index.

       THE DIRECTIONAL CLAIM MADE HERE IS WRONG. See revision 4 below. The
       FUNCTIONAL MISMATCH described next is real and is why the fix was
       correct, but revision 2 asserted that the old construction made the null
       TIGHTER in the structure-rich arm, and measurement says the opposite.

       The real, and still valid, problem: the observed
       pooled value is a cell-weighted median of concatenated distances, while
       the null was an unweighted mean of per-structure medians, so delta_um was
       not a clean observed-minus-null.

       Both are fixed the same way. For each permutation index, labels are
       shuffled independently inside every structure, the resulting distances
       are CONCATENATED across structures, and one median is taken. The old
       mean-of-medians null is computed alongside and written to table 54 so the
       size of the error is a measured quantity rather than an assertion.

    2. THE RAW RADIAL MIXED MODELS ARE FLAGGED AS SUPERSEDED by script 08's
       centred outcome, kept for continuity and excluded from the headline.

    3. BH CORRECTION IS APPLIED WITHIN DECLARED FAMILIES.

    4. THE PREREQUISITE CHECK IS CURRENT.

    5. STALE NUMBERS IN THE HEADER ARE CORRECTED.

WHAT CHANGED IN REVISION 3 (and why)

    THE TRIGGER: script 04 revision 5 re-baselined every structure. The
    background statistic moved from the section median to p25, the growth
    fraction from 0.50 to 0.25, and a minimum peak density of 1,750 cells/mm2
    was added, all calibrated against 83 expert annotations. Structures are
    therefore larger, fewer in the treated arm, and vascular objects that
    previously passed as foci are gone. Nothing in this script's METHOD was
    wrong; what changed is the input and what that input now implies.

    1. INPUT REPOINTED to structures_rev5. Every number this script produced
       against structures_rev4 is superseded.

    2. THE CELL COUNTS IN THE PSEUDOREPLICATION ARGUMENT ARE CORRECT AGAIN.
       Revision 2 quoted 14,054 treated cells in 6 structures and 197,112
       untreated in 61, which was the revision-4 structure set. On revision 5 it
       is 25,881 treated cells in 4 structures and 314,510 untreated in 76.
       Revision 2 explicitly claimed to have fixed exactly this class of stale
       number, so leaving it would have been the same error twice.

    3. fit_mixed NOW CHECKS THE STANDARD ERROR, NOT JUST THE COEFFICIENT.
       This is the load-bearing fix. The old guard was

           if res is None or not np.isfinite(res.params.get("arm", np.nan)):

       which tests the coefficient alone. statsmodels returns a result object
       with a FINITE coefficient and a NaN standard error, raising nothing, when
       the random-effect variance sits on the boundary at zero. That happens
       when an animal contributes a single structure, because the animal random
       intercept and the residual then become the same quantity and the variance
       is not identifiable from that animal's contribution.

       Under revision 5 TWO treated animals contribute exactly one structure
       each (1, 1 and 2 foci). So this is no longer a hypothetical: a NaN
       standard error would have produced a NaN p-value flowing silently into
       table 57 and into the forest plot's confidence interval. The fit is now
       required to return a finite coefficient, a finite positive standard
       error and a finite p, and a failure is announced and the row dropped
       rather than reported as NaN. fit_mode, n_structures_D1MT and
       n_structures_ref are added to the output so fallbacks are visible in the
       table rather than buried.

       This is the same error corrected in the shared statistics module in
       September, reappearing here because this script carries its own fitter.

    4. THE THIN-STRUCTURE WARNING NOW STATES THE CONSEQUENCE. Revision 2 noted
       only that per-structure spreads have no within-animal variance. With two
       single-structure animals it needs to say that the nested mixed model will
       sit on the boundary, that fallbacks to animal-only random effects are
       expected, and that this is a property of the three-per-arm design rather
       than a bug.

    NOT CHANGED, AND WHY: the nearest-neighbour censoring problem does not exist
    in this script and did not need addressing. nn_distances builds a KD-tree on
    the target cells only and queries k=1, so every anchor gets its exact
    nearest-target distance with no truncation. That was revision 2's second
    listed fix. The only fixed-k code here is MIXING_K = 10 in the polarisation
    mixing, which is a local composition score rather than a search for a
    target, so it has no censoring failure mode either. Script 05's 61 percent
    censoring under revision 5 is a property of script 05's k-neighbour index,
    and script 05 is frozen.

WHAT TO EXPECT FROM THE REVISION 5 INPUT
    - Two treated animals with a single structure each. Expect mixed model
      fallbacks to "animal only" and possibly dropped rows.
    - In table 54, null_sd_ratio_old_over_new is exactly 1.000 for the
      single-structure animals, because with one structure the old
      mean-of-medians null and the corrected pooled null are identical by
      construction. That part held.
      The prediction that the untreated animals would sit near
      1/sqrt(n_structures) was WRONG. See revision 4.
    - The 3-HK verdict should be unchanged. It is a property of the channel, not
      of the structures. If it changes, something is wrong.

WHAT CHANGED IN REVISION 4 (documentation only, NO numbers change)

    THE TRIGGER: the first structures_rev5 run measured the pooled null bug and
    it went the OPPOSITE WAY to what revision 2 asserted and revision 3 repeated.
    Nothing computed here was wrong, but the explanation printed alongside it
    was, so this revision corrects the record. Tables 50 to 58 are byte-identical
    to the revision 3 run; only the report text and the F44b labels change.

    WHAT REVISION 2 CLAIMED
        The old mean-of-medians null averaged independent permutations across
        structures, so its spread shrank as 1/sqrt(n_structures), making the
        null too TIGHT in the structure-rich untreated arm and inflating its z
        and p.

    WHAT WAS MEASURED (old null SD / corrected null SD, 300 permutations)
        G3_43106   1 structure    1.000   (predicted 1.000)
        G3_43111   1 structure    1.000   (predicted 1.000)
        G3_43118   2 structures   1.418   (predicted 0.816)
        G4_31438  12 structures   3.483   (predicted 0.289)
        G4_36463   9 structures   1.484   (predicted 0.333)
        G4_43109  53 structures   2.982   (predicted 0.137)

    The single-structure prediction held exactly, which confirms the mechanism
    at K = 1. Every K above 1 went the other way: the old null was too WIDE, not
    too tight, so the old z and p were DEFLATED and the correction makes z
    LARGER. Median |z| on G4_43109 goes from 9.90 to 48.69.

    WHY THE ARGUMENT WAS WRONG
        Simulated directly. Revision 2 treated each structure's null SD as a
        fixed quantity and divided it by sqrt(K). It is not fixed: a structure's
        null median has SD proportional to 1/sqrt(n_k), so the mean of K of them
        has SD proportional to 1/sqrt(K * n_k) = 1/sqrt(N_total). The pooled
        median over the same N_total anchors also has SD proportional to
        1/sqrt(N_total). They are the SAME ORDER, and with equal structures the
        simulated ratio is 1.0 at every K, not 1/sqrt(K).

        Above 1 arises from UNEQUAL structures, which is what this data has at
        28 to 3,984 anchors per structure. The pooled median is pinned by the
        few large structures and is therefore very stable across permutations,
        while the mean-of-medians absorbs the noise of the small ones. Simulated
        on 9 to 53 unequal structures the ratio runs 1.4 to 6.1, matching the
        observed 1.4 to 3.5.

    WHAT IS STILL TRUE, AND WHY THE FIX WAS RIGHT ANYWAY
        The functional mismatch was real: the observed pooled value is a
        cell-weighted median of concatenated distances and the old null was an
        unweighted mean of per-structure medians. Those are different
        quantities, so delta_um was not a clean observed-minus-null and the two
        were not on the same scale. Matching the functional is correct
        regardless of which way the SD moves. The fix stands; only its stated
        mechanism and direction were wrong.

        Consequence for reporting: the correction is ANTI-conservative, not
        conservative. Numbers from the corrected null carry larger z and smaller
        p than the revision 1 numbers they replace. Do not describe this
        correction as having made the analysis more cautious.

OUTPUTS
    figures/  F43 .. F46
    tables/   50 .. 58

USAGE
    conda activate sc_pre
    python AKOYA_06_Distance_Statistics.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================

IN_DIR = "/master/jlehle/WORKING/AKOYA/structures_rev5"
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
# Compute the revision-1 mean-of-medians null alongside the correct one, so the
# size of that error is measured rather than asserted. Costs nothing.
REPORT_OLD_POOLED_NULL = True

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
MIXING_K = 10                        # local composition score, not a target
                                     # search, so no censoring failure mode
MIXING_MIN_PER_CLASS = 5             # minimum cells per class to compute mixing

# ---- 3-HK dynamic range gate ------------------------------------------------
HK3_COL = "3-Hydroxykynurenine"
HK3_MIN_USABLE_P99 = 5.0             # p99 below this means no usable range
HK3_MIN_NONZERO_FRAC = 0.05          # fewer nonzero cells than this is unusable

# ---- mixed models -----------------------------------------------------------
RUN_MIXED_MODELS = True
MODEL_MAX_CELLS_PER_STRUCTURE = 3000   # subsample for tractability
MODEL_MIN_CELLS_PER_ANIMAL = 20        # animals below this are dropped from a fit
BH_ALPHA = 0.10

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
OK_COLOR = "#1B7837"
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


def benjamini_hochberg(pvals):
    p = np.asarray(pvals, float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return q
    idx = np.flatnonzero(ok)
    pv = p[idx]
    order = np.argsort(pv)
    m = len(pv)
    adj = pv[order] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(adj, 1.0)
    q[idx] = out
    return q


def empirical_p(obs, null):
    null = np.asarray(null, float)
    null = null[np.isfinite(null)]
    if not len(null) or not np.isfinite(obs):
        return np.nan, np.nan, np.nan, np.nan
    mu, sd = float(np.mean(null)), float(np.std(null, ddof=1))
    z = (obs - mu) / sd if sd > 0 else np.nan
    n_ext = int(np.sum(np.abs(null - mu) >= abs(obs - mu)))
    p = (n_ext + 1) / (len(null) + 1)
    return float(min(p, 1.0)), z, mu, sd


def nn_distances(anchor_xy, target_xy):
    """
    Exact nearest-target distance for every anchor. The KD-tree is built on the
    TARGET cells only and queried with k=1, so there is no k-neighbour
    truncation and therefore no censoring. Anchors that coincide with a target
    return 0. This is why revision 3 needed no censoring fix: script 05's
    problem was its k-neighbour index, which this script does not use.
    """
    if len(anchor_xy) == 0 or len(target_xy) == 0:
        return np.array([])
    tree = cKDTree(target_xy)
    d, _ = tree.query(anchor_xy, k=1)
    return np.asarray(d, dtype=float)


def subsample(idx, cap, generator):
    return (generator.choice(idx, size=cap, replace=False)
            if len(idx) > cap else idx)


def fit_mixed(df, outcome, label, extra_note=""):
    """
    Linear mixed model: outcome ~ arm, random intercept for animal, structure
    nested within animal. Falls back to animal-only random effects if the nested
    fit fails to converge.

    REVISION 3: the acceptance test now requires a usable STANDARD ERROR, not
    just a finite coefficient. statsmodels returns a finite coefficient with a
    NaN standard error, raising nothing, when the random-effect variance sits on
    the boundary at zero. Under structures_rev5 two treated animals contribute a
    single structure each, so the animal variance is not identifiable from their
    contribution and that failure mode is now expected rather than hypothetical.
    A NaN standard error would otherwise give a NaN p-value flowing silently
    into table 57 and a nonsense confidence interval in the forest plot.
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

    # REVISION 3 GUARD. A finite coefficient is not enough.
    _coef = float(res.params.get("arm", np.nan))
    _se = float(res.bse.get("arm", np.nan))
    _p = float(res.pvalues.get("arm", np.nan))
    if not (np.isfinite(_coef) and np.isfinite(_se) and _se > 0
            and np.isfinite(_p)):
        print(f"    NO USABLE SE: {label} ({mode}). coef={_coef} se={_se} "
              f"p={_p}.")
        print("      The random-effect variance is on the boundary at zero, "
              "which happens")
        print("      when an animal contributes a single structure. Row DROPPED "
              "rather than")
        print("      reported as NaN.")
        return None

    means = d.groupby("condition")["_y"].mean()
    n_struct_by_arm = d.groupby("condition")["struct_key"].nunique()
    return {
        "analysis": label, "outcome": outcome, "random_effects": mode,
        "fit_mode": mode,
        "n_animals_D1MT": int(d.loc[d["condition"] == "D1MT", "sample_id"].nunique()),
        "n_animals_ref": int(d.loc[d["condition"] == REFERENCE_ARM, "sample_id"].nunique()),
        "n_structures_D1MT": int(n_struct_by_arm.get("D1MT", 0)),
        "n_structures_ref": int(n_struct_by_arm.get(REFERENCE_ARM, 0)),
        "n_cells": int(len(d)), "n_animals": int(d["sample_id"].nunique()),
        "n_structures": int(d["struct_key"].nunique()),
        f"mean_{CONDITION_ORDER[0]}": float(means.get(CONDITION_ORDER[0], np.nan)),
        f"mean_{CONDITION_ORDER[1]}": float(means.get(CONDITION_ORDER[1], np.nan)),
        "coef_D1MT_vs_ref": _coef,
        "std_err": _se,
        "z_stat": float(res.tvalues.get("arm", np.nan)),
        "p_value": _p,
        "note": extra_note,
    }


_tee = Tee(os.path.join(TAB_DIR, "00_distance_stats_report.txt"))
sys.stdout = _tee

banner("AKOYA CORRECTED DISTANCE STATISTICS AND FORMAL TREATMENT TESTS (rev 4)")
print(f"Run time      : {datetime.now().isoformat(timespec='seconds')}")
print(f"Input         : {IN_DIR}")
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

print("\nREVISION 3 HEADLINE, STILL IN FORCE")
print("    Input repointed to structures_rev5. Every number this script wrote")
print("    against structures_rev4 is superseded.")
print()
print("    Two treated animals now contribute a SINGLE structure each. For")
print("    those animals the animal random intercept and the residual are the")
print("    same quantity, so the animal variance is not identifiable from their")
print("    contribution and the nested mixed model will sit on the boundary at")
print("    zero. statsmodels returns a finite coefficient with a NaN standard")
print("    error in that situation and raises nothing, so fit_mixed now tests")
print("    the standard error as well and drops the row loudly instead of")
print("    writing a NaN p-value into table 57.")

print("\nREVISION 2 FIX, STILL IN FORCE, WITH ITS RATIONALE CORRECTED")
print("    The pooled per-animal null was an unweighted mean of per-structure")
print("    null medians while the observed value is a cell-weighted median of")
print("    concatenated distances. Two different functionals, so delta_um was")
print("    not a clean observed-minus-null. It is now built by concatenating")
print("    permuted distances across structures and taking one median, exactly")
print("    as the observed value is built. Matching the functional is correct")
print("    regardless of which way the null SD moves.")
print()
print("    REVISION 4 CORRECTION: revision 2 claimed the old null was too TIGHT")
print("    in the structure-rich arm, shrinking as 1/sqrt(n_structures).")
print("    Measurement says the opposite. The old null was too WIDE for every")
print("    animal with more than one structure, so the correction makes z")
print("    LARGER, not smaller. It is ANTI-conservative. See the docstring for")
print("    the simulation. With one structure the two constructions coincide,")
print("    so a ratio of exactly 1.000 is expected there and that part held.")


# %% Cell 3 - load and prerequisite check
# =============================================================================

banner("LOADING")

paths = sorted(glob.glob(os.path.join(CELL_DIR, "*_cell_structures.csv")))
if not paths:
    print(f"ERROR: no cell assignment files in {CELL_DIR}.")
    print("Run script 04 revision 5 first, and check that IN_DIR points at")
    print("structures_rev5 rather than structures_rev4.")
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
    if "region" not in d.columns:
        print(f"    ERROR: {sid} lacks 'region'. This script expects the")
        print("    revision 4 / revision 5 schema. Skipping.")
        continue
    d["sample_id"] = sid
    cells[sid] = d
    n_struct = int(d.loc[d[STRUCT_COL] > 0, STRUCT_COL].nunique())
    print(f"    {sid:<12} {d['condition'].iloc[0]:<10} {len(d):>9,} cells, "
          f"{int((d[STRUCT_COL] > 0).sum()):>8,} in {n_struct} structures")

if not cells:
    print("ERROR: nothing loaded.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

cond_rank = {c: i for i, c in enumerate(CONDITION_ORDER)}
SAMPLE_ORDER = sorted(cells, key=lambda s: (cond_rank.get(cells[s]["condition"].iloc[0], 9), s))
COND_OF = {s: cells[s]["condition"].iloc[0] for s in SAMPLE_ORDER}
SCAN_OF = {s: (cells[s]["scan_id"].iloc[0] if "scan_id" in cells[s].columns
               else "na") for s in SAMPLE_ORDER}
MARKER_OF = {s: POINT_MARKERS[i % len(POINT_MARKERS)] for i, s in enumerate(SAMPLE_ORDER)}
SHADES = {"D1MT": ["#08519C", "#3182BD", "#6BAED6"],
          "Untreated": ["#A63603", "#E6550D", "#FD8D3C"]}
COLOR_OF, _seen = {}, {c: 0 for c in CONDITION_ORDER}
for s in SAMPLE_ORDER:
    c = COND_OF[s]
    pal = SHADES.get(c, ["#999999"])
    COLOR_OF[s] = pal[_seen.get(c, 0) % len(pal)]
    _seen[c] = _seen.get(c, 0) + 1

sub("STRUCTURE COUNTS PER ANIMAL")
n_struct_of = {}
ARM_CELLS = {}
for s in SAMPLE_ORDER:
    d = cells[s]
    n_struct_of[s] = int(d.loc[d[STRUCT_COL] > 0, STRUCT_COL].nunique())
for c in CONDITION_ORDER:
    mem = [s for s in SAMPLE_ORDER if COND_OF[s] == c]
    counts = {short_label(s): n_struct_of[s] for s in mem}
    total_cells = sum(int((cells[s][STRUCT_COL] > 0).sum()) for s in mem)
    total_struct = sum(n_struct_of[s] for s in mem)
    ARM_CELLS[c] = (total_cells, total_struct)
    print(f"    {c:<12} {counts}   {total_cells:,} cells in "
          f"{total_struct} structures")

thin = [s for s in SAMPLE_ORDER if n_struct_of[s] < 2]
if thin:
    print(f"\n    {[short_label(s) for s in thin]} contribute a SINGLE structure.")
    print("    For those animals the animal random intercept and the residual")
    print("    are the same quantity, so the animal variance is NOT")
    print("    identifiable from their contribution and a nested mixed model")
    print("    will sit on the boundary at zero. Expect fallbacks to")
    print("    'animal only', and expect some fits to return no usable standard")
    print("    error and be dropped. That is a property of the three-animals-")
    print("    per-arm design meeting a re-baselined structure set, not a bug.")
    print("    Per-structure spreads for those animals have no within-animal")
    print("    variance. Pooled per-animal results are still valid.")
print("\n    The structure counts differ enormously between arms. That is")
print("    exactly why the pooled null had to be rebuilt in revision 2: the")
print("    revision-1 construction made the null tighter in whichever arm had")
print("    more structures.")


def structure_frames():
    for s in SAMPLE_ORDER:
        d = cells[s]
        for k, g in d.loc[d[STRUCT_COL] > 0].groupby(STRUCT_COL):
            yield s, int(k), g


# %% Cell 4 - 3-HK dynamic range gate
# =============================================================================

banner("3-HYDROXYKYNURENINE DYNAMIC RANGE GATE")

print("    Before any gradient analysis is built on 3-HK, establish whether the")
print("    channel has usable range at all. Note that comparability and dynamic")
print("    range are different questions: script 03 puts 3-HK at ICC 0.250 on")
print("    p99 across scans, which is a comparability statement. This gate is")
print("    about whether the channel carries signal in the first place.")
print("\n    This verdict should be IDENTICAL to the revision 4 run. It is a")
print("    property of the channel, not of the structure definition. If it has")
print("    changed, something is wrong upstream.\n")

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
            "sample_id": s, "condition": COND_OF[s], "scan_id": SCAN_OF[s],
            "scope": scope, "n_cells": len(v),
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
HK3_USABLE = False
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
    print("    WARNING: no 3-HK data found.")


# %% Cell 5 - radial position, delta primary, per-structure and pooled
# =============================================================================

banner("RADIAL POSITION - EFFECT SIZES")

print("    delta = observed mean radial position minus the permutation null")
print("    mean, in radial units (0 = core centre, 1 = core boundary, 2 = outer")
print("    cuff). Negative = pulled toward the core.")
print("    The pooled null permutes labels WITHIN each structure and then pools,")
print("    which is the same construction the nearest-neighbour pooled null now")
print("    uses. This part was already correct in revision 1.\n")

rad_struct, rad_pool = [], []

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
        pval, z, mu, sd = empirical_p(obs, null)
        rad_struct.append({
            "sample_id": s, "condition": COND_OF[s], "structure_id": k,
            "phenotype": p, "n_cells": npos, "n_structure": n,
            "observed_radial": obs, "null_radial": mu, "null_sd": sd,
            "delta_radial": obs - mu, "z": z, "p_empirical": pval,
            "tier": tier_of(npos),
        })

for s in SAMPLE_ORDER:
    d = cells[s]
    sel = d.loc[(d[STRUCT_COL] > 0) & np.isfinite(
        pd.to_numeric(d["radial_pos"], errors="coerce"))]
    if not len(sel):
        continue
    r = pd.to_numeric(sel["radial_pos"], errors="coerce").to_numpy()
    ph = sel["pheno"].to_numpy()
    struct = sel[STRUCT_COL].to_numpy()
    uniq = np.unique(struct)
    struct_idx = {kk: np.flatnonzero(struct == kk) for kk in uniq}
    n = len(r)
    for p in PHENOTYPE_ORDER:
        m = ph == p
        npos = int(m.sum())
        if npos == 0:
            continue
        obs = float(r[m].mean())
        want = {kk: int(m[struct_idx[kk]].sum()) for kk in uniq}
        null = np.empty(N_PERMUTATIONS)
        for it in range(N_PERMUTATIONS):
            fake = np.zeros(n, dtype=bool)
            for kk in uniq:
                w = want[kk]
                if w:
                    fake[rng.choice(struct_idx[kk], size=w, replace=False)] = True
            null[it] = r[fake].mean()
        pval, z, mu, sd = empirical_p(obs, null)
        rad_pool.append({
            "sample_id": s, "animal_id": short_label(s), "condition": COND_OF[s],
            "scan_id": SCAN_OF[s], "phenotype": p,
            "n_cells": npos, "n_pooled": n, "n_structures": len(uniq),
            "observed_radial": obs, "null_radial": mu, "null_sd": sd,
            "delta_radial": obs - mu, "z": z, "p_empirical": pval,
            "tier": tier_of(npos),
        })
    print(f"    {s} pooled")

rad_s = pd.DataFrame(rad_struct)
rad_p = pd.DataFrame(rad_pool)
if len(rad_s):
    write_csv(rad_s, "51_radial_per_structure.csv")
if len(rad_p):
    rad_p["q_value"] = benjamini_hochberg(rad_p["p_empirical"].to_numpy())
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


# %% Cell 6 - nearest-neighbour distances, exact, correctly pooled
# =============================================================================

banner("NEAREST-NEIGHBOUR DISTANCES - EXACT, NO CENSORING, CORRECTLY POOLED")

print("    KD-tree built on target cells only and queried with k=1, so every")
print("    anchor gets its EXACT nearest-target distance. There is no")
print("    k-neighbour truncation here and therefore no censoring. Null:")
print("    reassign anchor and target labels within the structure, preserving")
print("    counts and coordinates.")
print("\n    THE POOLED NULL. For each permutation index, labels are shuffled")
print("    independently inside EVERY structure of that animal, the resulting")
print("    distances are concatenated across structures, and ONE median is")
print("    taken. The pooled null is therefore built exactly like the pooled")
print("    observed value: same functional, same cell weighting, and its width")
print("    does not depend on how many structures the animal happens to have.")
print("\n    With a SINGLE structure the corrected pooled null and the old")
print("    mean-of-medians null are identical by construction, so expect")
print("    null_sd_ratio_old_over_new to be exactly 1.000 for those animals.\n")

nn_struct, nn_pool, nn_cells_rows = [], [], []

for s in SAMPLE_ORDER:
    d = cells[s]
    sel = d.loc[d[STRUCT_COL] > 0]
    if not len(sel):
        continue

    # per-structure arrays, built once
    structs = {}
    for k, g in sel.groupby(STRUCT_COL):
        xy = g[["x", "y"]].to_numpy(float)
        if len(xy) < 10:
            continue
        structs[int(k)] = (xy, g["pheno"].to_numpy())

    for anchor in NN_ANCHORS:
        for target in NN_TARGETS:
            if target == anchor:
                continue

            # assemble the structures that carry both members
            items = []
            for k, (xy, ph) in structs.items():
                a_idx = np.flatnonzero(ph == anchor)
                t_idx = np.flatnonzero(ph == target)
                if not len(a_idx) or not len(t_idx):
                    continue
                a_use = subsample(a_idx, MAX_ANCHORS_PER_STRUCTURE, rng)
                items.append({"k": k, "xy": xy, "n": len(xy),
                              "n_a": len(a_idx), "n_t": len(t_idx),
                              "a_use": a_use, "t_idx": t_idx})
            if not items:
                continue

            # ---- observed ------------------------------------------------
            for item in items:
                item["obs_v"] = nn_distances(item["xy"][item["a_use"]],
                                             item["xy"][item["t_idx"]])
            pooled_obs_v = np.concatenate([it_["obs_v"] for it_ in items
                                           if len(it_["obs_v"])])
            pooled_obs = float(np.median(pooled_obs_v)) if len(pooled_obs_v) else np.nan

            # ---- one permutation pass, giving both nulls -----------------
            null_struct = {it_["k"]: np.full(N_PERMUTATIONS, np.nan)
                           for it_ in items}
            pooled_null = np.full(N_PERMUTATIONS, np.nan)
            for it in range(N_PERMUTATIONS):
                bucket = []
                for item in items:
                    n = item["n"]
                    perm = rng.permutation(n)
                    fa = perm[:item["n_a"]]
                    ft = perm[item["n_a"]:item["n_a"] + item["n_t"]]
                    fa_use = subsample(fa, MAX_ANCHORS_PER_STRUCTURE, rng)
                    v = nn_distances(item["xy"][fa_use], item["xy"][ft])
                    if len(v):
                        null_struct[item["k"]][it] = float(np.median(v))
                        bucket.append(v)
                if bucket:
                    pooled_null[it] = float(np.median(np.concatenate(bucket)))

            # ---- per-structure records -----------------------------------
            for item in items:
                k = item["k"]
                obs = (float(np.median(item["obs_v"]))
                       if len(item["obs_v"]) else np.nan)
                pval, z, mu, sd = empirical_p(obs, null_struct[k])
                nn_struct.append({
                    "sample_id": s, "condition": COND_OF[s],
                    "scan_id": SCAN_OF[s], "structure_id": k,
                    "anchor": anchor, "target": target,
                    "n_anchor": item["n_a"], "n_target": item["n_t"],
                    "n_structure": item["n"],
                    "observed_median_um": obs, "null_median_um": mu,
                    "null_sd_um": sd,
                    "delta_um": obs - mu if np.isfinite(mu) else np.nan,
                    "z": z, "p_empirical": pval,
                    "tier": tier_of(min(item["n_a"], item["n_t"])),
                })

                take = item["obs_v"]
                if len(take) > MODEL_MAX_CELLS_PER_STRUCTURE:
                    take = rng.choice(take, size=MODEL_MAX_CELLS_PER_STRUCTURE,
                                      replace=False)
                if len(take):
                    nn_cells_rows.append(pd.DataFrame({
                        "sample_id": s, "condition": COND_OF[s],
                        STRUCT_COL: k, "anchor": anchor, "target": target,
                        "nn_distance_um": take,
                    }))

            # ---- pooled record -------------------------------------------
            pval, z, mu, sd = empirical_p(pooled_obs, pooled_null)
            rec = {
                "sample_id": s, "animal_id": short_label(s),
                "condition": COND_OF[s], "scan_id": SCAN_OF[s],
                "anchor": anchor, "target": target,
                "n_anchor_cells": int(len(pooled_obs_v)),
                "n_structures": int(len(items)),
                "observed_median_um": pooled_obs,
                "null_median_um": mu, "null_sd_um": sd,
                "delta_um": pooled_obs - mu if np.isfinite(mu) else np.nan,
                "z": z, "p_empirical": pval,
                "tier": tier_of(len(pooled_obs_v)),
            }
            if REPORT_OLD_POOLED_NULL:
                old_null = np.nanmean(
                    np.vstack([null_struct[it_["k"]] for it_ in items]), axis=0)
                p_o, z_o, mu_o, sd_o = empirical_p(pooled_obs, old_null)
                rec.update({
                    "null_median_um_meanofmedians": mu_o,
                    "null_sd_meanofmedians": sd_o,
                    "delta_um_meanofmedians": (pooled_obs - mu_o
                                               if np.isfinite(mu_o) else np.nan),
                    "z_meanofmedians": z_o,
                    "p_meanofmedians": p_o,
                    "null_sd_ratio_old_over_new": (sd_o / sd if (sd and sd > 0)
                                                   else np.nan),
                })
            nn_pool.append(rec)

    print(f"    {s} done")
    gc.collect()

nn_s = pd.DataFrame(nn_struct)
nn_p = pd.DataFrame(nn_pool)
nn_cells = (pd.concat(nn_cells_rows, ignore_index=True)
            if nn_cells_rows else pd.DataFrame())

if len(nn_s):
    write_csv(nn_s, "53_nn_per_structure.csv")
if len(nn_p):
    nn_p["q_value"] = benjamini_hochberg(nn_p["p_empirical"].to_numpy())
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

    if REPORT_OLD_POOLED_NULL and "null_sd_ratio_old_over_new" in nn_p.columns:
        sub("HOW BIG WAS THE POOLED NULL BUG?")
        print("    Ratio = old mean-of-medians null SD / corrected pooled null")
        print("    SD. A ratio of exactly 1.000 means that animal has ONE")
        print("    structure, where the two constructions coincide.")
        print()
        print("    ABOVE 1 means the old null was too WIDE, so the old z and p")
        print("    were DEFLATED and the correction makes z LARGER. That is what")
        print("    this data shows, and it is the OPPOSITE of what revision 2")
        print("    predicted. Revision 2 expected roughly 1/sqrt(n_structures),")
        print("    which is wrong because a structure's own null SD is not fixed:")
        print("    it scales as 1/sqrt(n_k), so both estimators end up at")
        print("    1/sqrt(total anchors). Ratios above 1 come from UNEQUAL")
        print("    structure sizes, where the pooled median is pinned by the few")
        print("    large structures while the mean-of-medians absorbs noise from")
        print("    the small ones. The 1/sqrt column below is printed as the")
        print("    FAILED prediction, not as a target.\n")
        print(f"    {'section':<12}{'structures':>12}{'old sd/new sd':>16}"
              f"{'rev2 pred':>11}{'median |z| old':>16}{'median |z| new':>16}")
        print("    " + "-" * 83)
        for s in SAMPLE_ORDER:
            g = nn_p.loc[nn_p["sample_id"] == s]
            if not len(g):
                continue
            ns = float(g["n_structures"].median())
            ratio = float(g["null_sd_ratio_old_over_new"].median())
            print(f"    {s:<12}{ns:>12.0f}{ratio:>16.3f}"
                  f"{1/np.sqrt(ns):>11.3f}"
                  f"{float(g['z_meanofmedians'].abs().median()):>16.2f}"
                  f"{float(g['z'].abs().median()):>16.2f}")
        n_flip = int(((nn_p["p_meanofmedians"] < 0.05)
                      & (nn_p["p_empirical"] >= 0.05)).sum())
        print(f"\n    {n_flip} pooled comparison(s) were significant at p < 0.05")
        print("    under the old null and are not under the corrected one.")
        print("    Any pooled number previously quoted from table 54 is")
        print("    superseded, including the IDO1-negative to lymphocytes value")
        print("    carried in the project notes.")


# %% Cell 7 - polarisation, corrected
# =============================================================================

banner("iNOS vs ARGINASE-1 - CORRECTED REPORTING")

print("    pct_inos_only and pct_arg1_only are omitted: a within-section")
print("    percentile makes the two sets exactly equal in size, so those")
print("    columns are identical by construction. Script 05 revision 1")
print("    confirmed it empirically, 16.48 against 16.48 and so on.")
print("    Reported instead: double-positive against the independence")
print("    expectation, and spatial mixing with NO tier filtering.")
print("\n    MIXING_K is a local composition score over the k nearest")
print("    macrophages, not a search for a target, so it has no censoring")
print("    failure mode and did not need the revision 3 treatment.\n")

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
            "sample_id": s, "condition": COND_OF[s], "scan_id": SCAN_OF[s],
            "percentile": pct, "n_macrophages": len(mac),
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
        pval, z, mu, sd = empirical_p(obs, null)
        mix_rows.append({
            "sample_id": s, "condition": COND_OF[s], "structure_id": int(k),
            "n_macrophages": n, "n_inos_only": n_i, "n_arg1_only": n_a,
            "observed_mixing": obs, "null_mixing": mu, "null_sd": sd,
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
    print("    SPILLOVER CAVEAT: script 07's control pairs showed that apparent")
    print("    co-expression of any two markers rises with cell density, and")
    print("    CD3e with CD20 separated the arms more strongly than iNOS with")
    print("    Arginase-1 did. That comparison retired this as a finding. Do not")
    print("    read the ratio as biology without the control alongside.")
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
        print("    structures. Under revision 5 the treated arm holds only four")
        print("    structures in total, so this analysis is effectively dead in")
        print("    that arm. It was already retired as a finding; this simply")
        print("    confirms it cannot be revived with the current structure set.")


# %% Cell 8 - formal treatment tests by mixed model
# =============================================================================

banner("FORMAL TREATMENT TESTS - LINEAR MIXED MODELS")

_tc, _ts = ARM_CELLS.get("D1MT", (0, 0))
_uc, _us = ARM_CELLS.get("Untreated", (0, 0))
print("    Cells within an animal are not independent. A plain cell-level test")
print(f"    across {_tc:,} treated and {_uc:,} untreated cells in structures is")
print("    pseudoreplication and would return an implausibly small p-value for")
print("    an effect of no consequence. These models use every cell but hold")
print("    the effective sample size at the level of the six animals.")
print("    With three animals per arm, p in the 0.05 to 0.2 range is the")
print("    expected result for a real effect, not a failure.")
print(f"\n    Structure counts on this input: {_ts} treated, {_us} untreated.")
print("    Two treated animals contribute a single structure each, so expect")
print("    fallbacks to 'animal only' random effects, and expect some fits to")
print("    be DROPPED for returning no usable standard error. Every drop is")
print("    announced below.\n")
print("    RAW RADIAL MODELS ARE SUPERSEDED. Script 08 showed the correct")
print("    outcome is radial position CENTRED on each structure's own mean,")
print("    because that is the per-cell form of the delta effect size reported")
print("    in Cell 5. The raw-radial rows below are kept for continuity, are")
print("    flagged superseded_by_script_08, and must not be quoted.\n")

model_rows = []
n_dropped_no_se = 0
if RUN_MIXED_MODELS and HAVE_SM:
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
        r = fit_mixed(dd, "radial_pos", f"radial position: {p}",
                      extra_note="raw radial_pos, superseded by script 08")
        if r:
            r["phenotype"] = p
            r["family"] = "radial_position"
            r["superseded_by_script_08"] = True
            model_rows.append(r)
        else:
            n_dropped_no_se += 1

    if len(nn_cells):
        for (a, t), g in nn_cells.groupby(["anchor", "target"]):
            r = fit_mixed(g, "nn_distance_um", f"nn distance: {a} -> {t}",
                          extra_note="raw microns; bounded by focus size, see "
                                     "script 07")
            if r:
                r["anchor"], r["target"] = a, t
                r["family"] = "nn_distance"
                r["superseded_by_script_08"] = False
                model_rows.append(r)
                print(f"    nn {a[:18]:<20} -> {t:<16} "
                      f"coef={r['coef_D1MT_vs_ref']:>+7.2f} um  "
                      f"p={r['p_value']:.4f}  cells={r['n_cells']:,}  "
                      f"[{r['fit_mode']}]")
            else:
                n_dropped_no_se += 1
elif not HAVE_SM:
    print("    SKIPPED: statsmodels unavailable.")
else:
    print("    SKIPPED: RUN_MIXED_MODELS = False")

models = pd.DataFrame(model_rows)
if n_dropped_no_se:
    sub("MODELS DROPPED")
    print(f"    {n_dropped_no_se} fit(s) returned no usable standard error or")
    print("    could not be fitted, and were dropped rather than written as")
    print("    NaN. Revision 2 would have written those rows into table 57 with")
    print("    a NaN p-value and a nonsense confidence interval in F45.")
    print("    This is the three-animals-per-arm design meeting a structure set")
    print("    in which two treated animals contribute one structure each.")

if len(models):
    models["q_value"] = np.nan
    for fam, g in models.groupby("family"):
        models.loc[g.index, "q_value"] = benjamini_hochberg(
            g["p_value"].to_numpy())
    models["sig_q"] = models["q_value"] < BH_ALPHA
    write_csv(models, "57_mixed_model_results.csv")

    sub("Fit modes actually used")
    print(models.groupby(["family", "fit_mode"]).size().to_string())
    n_fallback = int((models["fit_mode"] == "animal only").sum())
    if n_fallback:
        print(f"\n    {n_fallback} of {len(models)} fits fell back to")
        print("    animal-only random effects, meaning the nested")
        print("    structure-within-animal term could not be estimated. Report")
        print("    the fit mode alongside any coefficient quoted from this")
        print("    table; the two modes are not the same model.")

    sup = models.loc[models["superseded_by_script_08"]]
    if len(sup):
        sub("Superseded raw-radial models (printed for the record only)")
        for _, r in sup.iterrows():
            print(f"    {r['analysis'][:44]:<46} coef={r['coef_D1MT_vs_ref']:>+7.3f}  "
                  f"p={r['p_value']:.4f}  q={r['q_value']:.4f}   DO NOT QUOTE")

    sub("Multiple testing position")
    print(f"    {int(models['sig_q'].sum())} of {len(models)} models survive BH "
          f"at q < {BH_ALPHA} within their declared family.")

    sub("Interpretation")
    print("    coef is the D1MT effect relative to Untreated, in the units of")
    print("    the outcome (radial units, or microns). A negative distance coef")
    print("    means the pair sits closer together in treated foci, but raw")
    print("    microns are bounded by focus size, which is why script 07")
    print("    re-tests distance three independent ways.")
else:
    print("\n    No mixed models were fitted.")


# %% Cell 9 - figures
# =============================================================================

banner("FIGURES")

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
        ax.set_ylabel("median distance (um)", fontsize=FONT_SIZE_BASE - 14)
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
                 "uncensored\nNull built by concatenating permuted distances "
                 "across structures, then taking one median",
                 y=1.02, fontsize=FONT_SIZE_TITLE - 8)
    save_fig(fig, "F44_nn_pooled_distances")

# ---- F44b the pooled null bug, made visible ---------------------------------
if len(nn_p) and REPORT_OLD_POOLED_NULL and "null_sd_meanofmedians" in nn_p.columns:
    fig, axes = plt.subplots(1, 2, figsize=(30, 13))
    ax = axes[0]
    for s in SAMPLE_ORDER:
        g = nn_p.loc[nn_p["sample_id"] == s]
        if not len(g):
            continue
        ax.scatter(g["n_structures"], g["null_sd_ratio_old_over_new"], s=320,
                   color=COLOR_OF[s], marker=MARKER_OF[s], edgecolor="#FFFFFF",
                   linewidth=2, zorder=3)
    ns = np.linspace(1, max(2, float(nn_p["n_structures"].max())), 100)
    ax.plot(ns, 1 / np.sqrt(ns), color="#000000", linestyle="--", linewidth=3)
    ax.axhline(1.0, color="#999999", linestyle=":", linewidth=2.5)
    ax.set_xscale("log")
    ax.set_xlabel("structures pooled in that animal")
    ax.set_ylabel("old null SD / corrected null SD")
    ax.set_title("Revision 2 predicted 1/sqrt(structures). It is wrong.\n"
                 "dashed = the failed prediction, dotted = no difference",
                 fontsize=FONT_SIZE_TITLE - 14)
    style_axes(ax)
    ax.legend(handles=[Line2D([0], [0], color=COLOR_OF[s], marker=MARKER_OF[s],
                              markersize=16, linestyle="none",
                              label=f"{short_label(s)} ({COND_OF[s]})")
                       for s in SAMPLE_ORDER],
              frameon=False, fontsize=FONT_SIZE_LEGEND - 16, loc="best")

    ax = axes[1]
    for s in SAMPLE_ORDER:
        g = nn_p.loc[nn_p["sample_id"] == s]
        if not len(g):
            continue
        ax.scatter(g["z_meanofmedians"].abs(), g["z"].abs(), s=320,
                   color=COLOR_OF[s], marker=MARKER_OF[s], edgecolor="#FFFFFF",
                   linewidth=2, zorder=3)
    lim = float(np.nanmax([nn_p["z"].abs().max(),
                           nn_p["z_meanofmedians"].abs().max(), 1.0]))
    ax.plot([0, lim], [0, lim], color="#000000", linestyle="--", linewidth=3)
    ax.set_xlabel("|z| under the old mean-of-medians null")
    ax.set_ylabel("|z| under the corrected pooled null")
    ax.set_title("Points ABOVE the line were understated by the old null",
                 fontsize=FONT_SIZE_TITLE - 14)
    style_axes(ax)
    fig.suptitle("How much the pooled null construction mattered, and in which "
                 "direction\nThe correction is ANTI-conservative: it makes z "
                 "larger, not smaller", y=1.03, fontsize=FONT_SIZE_TITLE - 8)
    save_fig(fig, "F44b_pooled_null_correction")

if len(models):
    m = models.loc[~models["superseded_by_script_08"]].copy()
    if len(m):
        m["label"] = m["analysis"]
        m = m.sort_values("p_value")
        fig, ax = plt.subplots(figsize=(22, max(12, 0.75 * len(m))))
        yy = np.arange(len(m))
        lo = m["coef_D1MT_vs_ref"] - 1.96 * m["std_err"]
        hi = m["coef_D1MT_vs_ref"] + 1.96 * m["std_err"]
        for i, (_, r) in enumerate(m.iterrows()):
            sig = bool(r["sig_q"])
            ax.plot([lo.iloc[i], hi.iloc[i]], [i, i],
                    color=FLAG_COLOR if sig else "#666666", linewidth=5, zorder=2)
            ax.scatter([r["coef_D1MT_vs_ref"]], [i], s=420,
                       color=FLAG_COLOR if sig else "#666666",
                       edgecolor="#FFFFFF", linewidth=2, zorder=3)
            ax.text(hi.iloc[i], i,
                    f"  p={r['p_value']:.3f} q={r['q_value']:.3f} "
                    f"[{r['fit_mode']}]",
                    va="center", fontsize=FONT_SIZE_ANNOT - 12)
        ax.axvline(0, color="#000000", linewidth=3.5)
        ax.set_yticks(yy)
        ax.set_yticklabels(m["label"], fontsize=FONT_SIZE_TICK - 12)
        ax.invert_yaxis()
        ax.set_xlabel("D1MT effect vs Untreated (95% CI)")
        ax.set_title("Mixed model treatment effects, superseded rows excluded\n"
                     "every interval has a finite standard error by construction",
                     fontsize=FONT_SIZE_TITLE - 10)
        ax.xaxis.grid(True, color=GRID_COLOR, linewidth=1.5); ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        save_fig(fig, "F45_mixed_model_forest")

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
print(f"Input                            : {IN_DIR}")
print(f"Radial rows (structure / pooled) : {len(rad_s)} / {len(rad_p)}")
print(f"NN rows (structure / pooled)     : {len(nn_s)} / {len(nn_p)}")
print(f"Mixing rows                      : {len(mix)}")
print(f"Mixed models fitted              : {len(models)}")
print(f"Mixed models dropped (no SE)     : {n_dropped_no_se}")
print(f"3-HK usable for Q3               : {HK3_USABLE}")
print(f"Structures                       : {_ts} treated, {_us} untreated")

sub("Read in this order")
print("  1. The MODELS DROPPED block and the fit-mode table in Cell 8. With")
print("     two single-structure treated animals, which fits survived and in")
print("     which mode decides what can be quoted at all.")
print("  2. F44b / table 54 : how much the pooled null construction mattered,")
print("     and in which direction. The single-structure animals sit at exactly")
print("     1.000. Everything above one structure sits ABOVE 1, meaning the old")
print("     null was too WIDE and the correction makes z LARGER. That is the")
print("     opposite of what revision 2 predicted; see the docstring.")
print("  3. F46 / table 50  : is 3-HK usable? Should be unchanged from the")
print("     revision 4 run.")
print("  4. F43 / table 52  : radial position as delta. Q1 architecture.")
print("  5. F44 / table 54  : uncensored distances pooled per animal. Q1, Q4.")
print("  6. F45 / table 57  : mixed model treatment effects.")
print("  7. Table 58        : per-animal consistency and arm differences.")

sub("What to tell Deepak about the statistics")
print("  Three levels, all reported:")
print("    - Within-structure permutation: high confidence, well powered, valid.")
print("      The null contains that structure's own size, shape and density.")
print("    - Per-animal effect sizes with k of 3 consistency: carries the")
print("      treatment claim honestly.")
print("    - Mixed model: the formal treatment test. Bounded by three animals")
print("      per arm. A p of 0.08 here is a real result, and any claim of")
print("      p < 0.001 for a treatment effect from three animals per group is")
print("      misusing the cell count.")
print()
print("  New with revision 5: two treated animals contribute one structure")
print("  each, so the focus-level mixed model has almost nothing to estimate")
print("  the animal variance from in that arm. Improving detection and")
print("  improving statistical power were never the same goal. Where a fit")
print("  drops or falls back, lead with the animal-level effect size and the")
print("  exact randomization test rather than reaching for a different model.")

sub("Superseded")
print("  Every number this script wrote against structures_rev4.")
print("  Every pooled p and z written by revision 1 of this script.")
print("  The STATED RATIONALE for the revision 2 pooled-null fix, corrected in")
print("  revision 4. The fix itself stands; its claimed direction did not. Any")
print("  text describing that correction as having made the analysis more")
print("  cautious is wrong and needs rewording before it reaches Methods.")
print("  Every raw-radial mixed model row, superseded by script 08's centred")
print("  outcome.")

banner("DONE")
sys.stdout = _tee.terminal
_tee.close()
