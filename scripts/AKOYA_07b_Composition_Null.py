#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - IS THE DISTANCE DELTA COMPOSITION-CONTAMINATED?
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 07b of the AKOYA analysis series. REVISION 3.
A DIAGNOSTIC. It writes to its own directory and changes nothing script 07 owns.

REVISION 3, 2 October 2026. DENOMINATOR DEGREES OF FREEDOM.

    Same change as script 08 revision 7 and script 07 revision 3.6.
    statsmodels reports a Wald z on each fixed effect with no small-sample
    correction, and res.df_resid is cells minus rank, so the arm contrast was
    tested as though it carried thousands of degrees of freedom when it is a
    between-animal contrast with about n_animals - 2 = 4.

    akoya_arm_stats.small_sample_inference refers the same statistic to t(4)
    and adds the matching interval. p_wald_z keeps what revision 2 reported as
    p_value.

    NOTHING ABOUT THIS SCRIPT'S VERDICT CHANGES. The question here is whether a
    per-animal delta separates the six animals under each of three anchor
    references, and that verdict comes from complete_separation and arm_margin
    on table 83, neither of which involves a model. The r = 0.9995 validation
    against table 53 compares distances in microns, not p-values, so it is
    untouched as well. Table 84's models are the secondary column and they are
    what moves.

    One consequence worth stating. BH in this script runs within each reference
    separately, and at BH_ALPHA = 0.10 on five targets per anchor the corrected
    p-values will clear fewer rows. That is expected and it does not bear on
    section D, which is the part of this script that answers the question.

WHY THIS SCRIPT EXISTS

    Scripts 07 and 08 agree on direction and disagree on which population
    carries the spatial effect.

        script 07, confirmatory delta : IDO1+ macrophages to PLASMA cells,
                                        margin 16.70 um, the largest of four
        script 08, composition-corrected : B CELLS, the only reportable
                                        population surviving BH; plasma cells
                                        loses per-animal separation entirely

    Script 07 argues these are two views of one spatial fact, because both
    nearest-neighbour anchors sit in the myeloid pool that defines a focus and
    therefore at the density peak, so anchor-to-non-pool-target distance is
    largely radial position in different units. If that is right, the two
    scripts should name the same population. They do not.

    THE HYPOTHESIS THIS SCRIPT TESTS. Script 07's delta is the observed
    nearest-neighbour distance minus the median of a permutation null that
    shuffles phenotype labels within the structure with coordinates and counts
    held fixed. Under that null the fake anchors are drawn UNIFORMLY from all
    cells in the structure, so the reference is "the average CELL".

    That is the same construction script 08's PRIMARY centring uses, and script
    08 has now measured what it costs. Moving from the average cell to the
    average CELL TYPE took pooled lymphocytes from -0.1454 to -0.0622 and T
    lineage from -0.0976 to -0.0119. Script 07 has never been tested for it.

    So: rebuild script 06's null with the anchor side drawn composition-balanced
    instead of uniformly, change nothing else, and see whether the plasma cell
    result is composition.

WHAT CHANGED FROM THE WRITTEN PLAN, AND WHY

    The plan proposed a closed form: reference = median nearest-neighbour
    distance from all non-target cells to the OBSERVED target cells, with no
    permutation needed. That is wrong, and it would have been wrong in a way
    that invalidated the validation step.

    Script 06's null randomises BOTH sides. It draws fake anchors and a disjoint
    fake target set, both uniformly from the structure's cells, and takes the
    median distance between them. A closed form that keeps the observed targets
    is therefore a DIFFERENT quantity from null_median_um in table 53, so it
    could not have been validated against it, and it would have changed two
    things at once: the reference type and the target handling.

    This script instead keeps script 06's null structure exactly and changes
    ONE line: how the fake anchors are drawn. Three nulls are computed from the
    same permutation stream and the same fake target set, so the target side is
    identical across all three and only the anchor reference varies.

        null_uniform    fake anchors drawn uniformly from all cells.
                        This is script 06's null. It MUST reproduce table 53.
        null_bal_adapt  fake anchors drawn with each phenotype present in that
                        structure carrying equal total weight.
        null_bal_common fake anchors drawn with each phenotype in the GLOBAL
                        common set carrying equal total weight, matching the
                        common-set centring in script 08 revision 3.1.

    Anchor COUNT barely affects the statistic, because the statistic is the
    MEDIAN over anchors, which is why a balanced draw of slightly different size
    is legitimate here. Target count is preserved exactly, because
    nearest-neighbour distance depends strongly on target density and that is
    what the existing null correctly controls.

VALIDATION BEFORE ANY RESULT

    null_uniform is recomputed here and compared against null_median_um from
    script 06 table 53, structure by structure. If they do not agree, this
    script's reconstruction of the null is wrong and the balanced versions mean
    nothing. The verdict is printed BEFORE any delta, and the balanced results
    are suppressed if the agreement fails. Nothing here may be read past a
    failed validation.

WHAT THIS SCRIPT DECIDES

    If plasma cells holds under the balanced anchor reference, scripts 07 and 08
    are measuring genuinely different things and both stand.

    If plasma cells collapses and B cells strengthens, script 07's ordering of
    the four confirmatory pairs was composition, and the spatial story is a
    B cell story throughout.

    Either way it is falsifiable, and it is the last structural question in the
    spatial analysis.

WHAT THIS SCRIPT IS NOT

    Not a replacement for script 07. Script 07 remains the owner of the
    confirmatory distance family. This script reports whether that family's
    reference is contaminated, and by how much. If it finds contamination, the
    fix belongs in script 07 as a revision, not here.

INPUTS
    structures_rev5/cell_assignments/<section>_cell_structures.csv
    structures_rev5/tables/35_foci_structures_relative.csv
    distance_stats/tables/53_nn_per_structure.csv   (the null being tested)

OUTPUTS
    composition_null/tables/   80 .. 85
    composition_null/figures/  F60 .. F62

USAGE
    conda activate sc_pre
    python AKOYA_07b_Composition_Null.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================

IN_DIR = "/master/jlehle/WORKING/AKOYA/structures_rev5"
CELL_DIR = f"{IN_DIR}/cell_assignments"
FOCI_TABLE = f"{IN_DIR}/tables/35_foci_structures_relative.csv"
NN_NULL_TABLE = "/master/jlehle/WORKING/AKOYA/distance_stats/tables/53_nn_per_structure.csv"
OUT_DIR = "/master/jlehle/WORKING/AKOYA/composition_null"

CONDITION_ORDER = ["D1MT", "Untreated"]
CONDITION_COLORS = {"D1MT": "#2C7FB8", "Untreated": "#D95F02"}
REFERENCE_ARM = "Untreated"

STRUCTURE_SOURCE = "foci"
STRUCT_COL = "focus_id" if STRUCTURE_SOURCE == "foci" else "burden_region_id"

IDO1_POS = "CD68+IDO1+ Macrophages"
IDO1_NEG = "CD68+IDO1- Macrophages"
DETECTION_POOL = [IDO1_POS, IDO1_NEG, "CD163+ Macrophages", "Neutrophils"]

# The four confirmatory pairs from script 07, plus Tregs so the evidence funnel
# stays permissive. Reportability is a label applied after the fact, never a
# pre-filter.
NN_ANCHORS = [IDO1_POS, IDO1_NEG]
NN_TARGETS = ["Plasma cells", "B cells", "CD4- T cells", "Helper T cells", "Tregs"]

# ---- the balanced anchor draw ----------------------------------------------
MIN_CELLS_FOR_BALANCE = 10        # matches script 08
BALANCE_COMMON_MIN_FRAC = 0.90    # matches script 08
MIN_ANCHORS_PER_PHENO = 3         # below this a phenotype contributes nothing

# ---- permutation -----------------------------------------------------------
N_PERMUTATIONS = 300              # matches scripts 06 and 07
MAX_CELLS_PER_STRUCTURE = 2000    # cap on any drawn anchor set, as script 07
MIN_CELLS_IN_STRUCTURE = 10
RANDOM_SEED = 0

# ---- validation, the gate on everything below ------------------------------
# null_uniform recomputed here against null_median_um from table 53.
VALIDATE_AGAINST_TABLE_53 = True
VALIDATION_MIN_CORR = 0.98        # Pearson r across matched structure-pairs
VALIDATION_MAX_MEDIAN_ABS_DIFF_UM = 3.0
REFUSE_ON_VALIDATION_FAILURE = True

# ---- reporting -------------------------------------------------------------
MIN_TREATED_TARGETS_FOR_REPORT = 50   # matches script 07
BH_ALPHA = 0.10

# ---- models ----------------------------------------------------------------
RUN_MIXED_MODELS = True
MODEL_MAX_CELLS_PER_STRUCTURE = 2000
MODEL_MIN_CELLS_PER_ANIMAL = 20
RUN_EXACT_RANDOMIZATION = True
RUN_EXACT_LMM_REFIT = True
PERM_MAX_CELLS = 20000
EXACT_P_REL_TOL = 1e-6            # see script 08 revision 3.1 item A

# ---- figures ---------------------------------------------------------------
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


# %% Cell 2 - imports, style, helpers
# =============================================================================

import os
import sys
import glob
import gc
import warnings
from datetime import datetime
from itertools import combinations

import numpy as np
import pandas as pd

# REVISION 3: the shared stats module, for small_sample_inference. This script
# previously imported nothing from it and defined its own fit_mixed, which is how
# the Wald-z defect could live here and in 06, 07 and 08 at the same time.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import akoya_arm_stats as aas

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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
    "font.size": FONT_SIZE_BASE,
    "axes.titlesize": FONT_SIZE_TITLE,
    "axes.labelsize": FONT_SIZE_BASE,
    "xtick.labelsize": FONT_SIZE_TICK,
    "ytick.labelsize": FONT_SIZE_TICK,
    "legend.fontsize": FONT_SIZE_LEGEND,
    "figure.dpi": 100,
    "savefig.dpi": DPI,
    "axes.edgecolor": AXIS_COLOR,
    "text.color": TEXT_COLOR,
    "axes.labelcolor": TEXT_COLOR,
    "xtick.color": TEXT_COLOR,
    "ytick.color": TEXT_COLOR,
})

FIG_DIR = os.path.join(OUT_DIR, "figures")
TAB_DIR = os.path.join(OUT_DIR, "tables")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TAB_DIR, exist_ok=True)

rng_perm = np.random.default_rng(RANDOM_SEED + 1)


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


def benjamini_hochberg(pvals):
    p = np.asarray(pvals, float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    if not ok.sum():
        return q
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order]
    m = len(ranked)
    adj = ranked * m / (np.arange(m) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(adj, 1.0)
    q[ok] = out
    return q


def complete_separation(a, b):
    a = [x for x in a if np.isfinite(x)]
    b = [x for x in b if np.isfinite(x)]
    if not a or not b:
        return False
    return (max(a) < min(b)) or (max(b) < min(a))


def arm_margin(t, u):
    t = [x for x in t if np.isfinite(x)]
    u = [x for x in u if np.isfinite(x)]
    if not t or not u:
        return np.nan
    if min(t) > max(u):
        return float(min(t) - max(u))
    if min(u) > max(t):
        return float(min(u) - max(t))
    return np.nan


def nn_median(anchor_xy, target_xy):
    """Median exact nearest-neighbour distance, anchors to targets."""
    if len(anchor_xy) == 0 or len(target_xy) == 0:
        return np.nan
    d, _ = cKDTree(target_xy).query(anchor_xy, k=1)
    return float(np.median(d))


def nn_all(anchor_xy, target_xy):
    if len(anchor_xy) == 0 or len(target_xy) == 0:
        return np.array([])
    d, _ = cKDTree(target_xy).query(anchor_xy, k=1)
    return np.asarray(d, float)


# ---- the balanced anchor draw ----------------------------------------------
def balanced_draw(rng, pools, n_draw):
    """
    Draw up to n_draw cell indices giving EACH PHENOTYPE equal total weight,
    rather than each cell equal weight. `pools` maps phenotype -> array of
    candidate indices, already stripped of the fake target cells.

    The returned set may be smaller than n_draw when phenotypes are thin. That
    is acceptable here and it is the reason this construction is legitimate:
    the statistic is the MEDIAN distance over anchors, which is insensitive to
    how many anchors there are, unlike the target count which is preserved
    exactly.
    """
    avail = [p for p, ix in pools.items() if len(ix) >= MIN_ANCHORS_PER_PHENO]
    if not avail:
        return None
    per = max(1, int(np.ceil(n_draw / len(avail))))
    out = []
    for p in avail:
        pool = pools[p]
        take = int(min(per, len(pool)))
        out.append(rng.choice(pool, size=take, replace=False))
    out = np.concatenate(out)
    if len(out) > n_draw:
        out = rng.choice(out, size=n_draw, replace=False)
    return out


# ---- models ----------------------------------------------------------------
def _usable(res, term="arm"):
    """Ported from script 08 revision 3.2. A finite coefficient is not enough."""
    try:
        c = float(res.params.get(term, np.nan))
        e = float(res.bse.get(term, np.nan))
        pv = float(res.pvalues.get(term, np.nan))
    except Exception:
        return False
    return bool(np.isfinite(c) and np.isfinite(e) and e > 0 and np.isfinite(pv))


def _assignment_reps(n, n_t):
    """One representative per sign-symmetric pair when the design is balanced."""
    if n % 2 == 0 and n_t * 2 == n:
        for rest in combinations(range(1, n), n_t - 1):
            yield (0,) + rest, True
    else:
        for idx in combinations(range(n), n_t):
            yield idx, False


def _finish_exact(stats, obs, label):
    stats = np.asarray([v for v in stats if np.isfinite(v)], float)
    if not len(stats):
        return np.nan, 0
    tol = EXACT_P_REL_TOL * max(1.0, abs(obs))
    hits = int(np.sum(stats >= abs(obs) - tol))
    p = hits / len(stats)
    floor = 1.0 / len(stats)
    if np.isfinite(p) and p < floor - 1e-12:
        print(f"    EXACT-P FLOOR VIOLATED ({label}): {p:.4f} < {floor:.4f}. "
              f"Raising to the floor; this indicates a defect.")
        p = floor
    return float(p), int(len(stats))


def exact_p_animal_means(d, outcome):
    g = d.dropna(subset=[outcome, "sample_id", "condition"])
    if not len(g):
        return np.nan
    m = g.groupby(["sample_id", "condition"])[outcome].mean().reset_index()
    animals = list(m["sample_id"])
    vals = m[outcome].to_numpy(float)
    arms = list(m["condition"])
    n = len(animals)
    n_t = sum(1 for a in arms if a != REFERENCE_ARM)
    if n_t == 0 or n_t == n:
        return np.nan
    obs_idx = tuple(i for i, a in enumerate(arms) if a != REFERENCE_ARM)

    def _stat(idx):
        t = vals[list(idx)]
        u = vals[[i for i in range(n) if i not in idx]]
        return float(np.mean(t) - np.mean(u))

    obs = _stat(obs_idx)
    stats = [abs(_stat(idx)) for idx, _ in _assignment_reps(n, n_t)]
    p, _ = _finish_exact(stats, obs, "p_exact_means")
    return p


def fit_mixed(df, outcome, label):
    """outcome ~ arm, animal random intercept, structure nested within animal."""
    if not (HAVE_SM and RUN_MIXED_MODELS):
        return None
    need = [outcome, "condition", "sample_id", STRUCT_COL]
    if [c for c in need if c not in df.columns]:
        return None
    d = df.dropna(subset=need).copy()
    counts = d.groupby("sample_id").size()
    d = d.loc[d["sample_id"].isin(
        counts.loc[counts >= MODEL_MIN_CELLS_PER_ANIMAL].index)]
    if not len(d) or d["condition"].nunique() < 2 or d["sample_id"].nunique() < 3:
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
        if res is None or not _usable(res):
            try:
                md = smf.mixedlm("_y ~ arm", data=d, groups=d["sample_id"],
                                 re_formula="1")
                res = md.fit(reml=True, method="lbfgs", maxiter=200)
                mode = "animal only"
            except Exception:
                return None
    if res is None or not _usable(res):
        print(f"    NO USABLE SE: {label} / {outcome}. Row DROPPED.")
        return None

    n_struct_arm = d.groupby("condition")["struct_key"].nunique()

    # REVISION 3: p and CI from the number of ANIMALS, not the number of cells.
    # The coefficient and standard error are untouched; only the reference
    # distribution changes.
    _coef = float(res.params["arm"])
    _se = float(res.bse["arm"])
    _ss = aas.small_sample_inference(_coef, _se, d["sample_id"].nunique())

    return {
        "analysis": label, "outcome": outcome, "fit_mode": mode,
        "coef_D1MT_vs_ref": _coef,
        "std_err": _se,
        "p_value": _ss["p_value"],
        "p_wald_z": float(res.pvalues["arm"]),
        "ci_low": _ss["ci_low"], "ci_high": _ss["ci_high"],
        "df": _ss["df"], "t_crit": _ss["t_crit"],
        "df_method": _ss["df_method"],
        "p_exact_means": exact_p_animal_means(d, "_y")
        if RUN_EXACT_RANDOMIZATION else np.nan,
        "n_cells": int(len(d)), "n_animals": int(d["sample_id"].nunique()),
        "n_structures": int(d["struct_key"].nunique()),
        "n_structures_D1MT": int(n_struct_arm.get("D1MT", 0)),
        "n_structures_ref": int(n_struct_arm.get(REFERENCE_ARM, 0)),
    }


_tee = Tee(os.path.join(TAB_DIR, "00_composition_null_report.txt"))
sys.stdout = _tee

banner("AKOYA COMPOSITION-AWARE NULL FOR THE DISTANCE DELTA (script 07b rev 3)")
print("Inference     : model p and CI on t(n_animals - 2) df via "
      "akoya_arm_stats.")
print("                p_wald_z retains what revision 2 reported. Section D's")
print("                verdict comes from animal-level separation, not from a")
print("                model, so the verdict is unaffected by this change.")
print(f"Run time      : {datetime.now().isoformat(timespec='seconds')}")
print(f"Input         : {IN_DIR}")
print(f"Null tested   : {NN_NULL_TABLE}")
print(f"Output        : {OUT_DIR}")
if not HAVE_SCIPY:
    print(f"\nERROR: scipy required ({_scipy_err}).")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)
if not HAVE_SM:
    print(f"\n    WARNING: statsmodels unavailable ({_sm_err}). Models skipped.")

print("""
THE QUESTION
    Script 07's delta subtracts a null whose fake anchors are drawn UNIFORMLY
    from all cells in the structure, so its reference is "the average CELL".
    Script 08 measured what that reference costs on the radial outcome: pooled
    lymphocytes went from -0.1454 to -0.0622, and T lineage to nearly zero,
    when the reference became "the average CELL TYPE". Script 07's delta has
    never been tested for the same thing, and scripts 07 and 08 currently name
    different populations: plasma cells and B cells respectively.

THREE NULLS, ONE CHANGED LINE
    All three share the same permutation stream and the SAME fake target set,
    so the target side is identical and only the anchor reference varies.
      null_uniform     anchors uniform over cells        <- script 06's null
      null_bal_adapt   anchors balanced over phenotypes present in the structure
      null_bal_common  anchors balanced over the global common set

VALIDATION GATES EVERYTHING
    null_uniform is recomputed here and must reproduce null_median_um from
    table 53. If it does not, this script's reconstruction is wrong and the
    balanced versions mean nothing. The verdict is printed before any delta.
""")


# %% Cell 3 - load, provenance, and the common phenotype set
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
    miss = [c for c in [STRUCT_COL, "x", "y", "pheno", "condition"]
            if c not in d.columns]
    if miss:
        print(f"    ERROR: {sid} lacks {miss}. Skipping.")
        continue
    d = d.loc[d[STRUCT_COL] > 0].copy()
    if not len(d):
        continue
    cells[sid] = d
    print(f"    {sid:<12} {d['condition'].iloc[0]:<10} {len(d):>8,} cells in "
          f"{int(d[STRUCT_COL].nunique()):>3} structures")

if not cells:
    print("ERROR: nothing loaded.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

cond_rank = {c: i for i, c in enumerate(CONDITION_ORDER)}
SAMPLE_ORDER = sorted(cells, key=lambda s: (
    cond_rank.get(cells[s]["condition"].iloc[0], 9), s))
COND_OF = {s: cells[s]["condition"].iloc[0] for s in SAMPLE_ORDER}
MARKER_OF = {s: POINT_MARKERS[i % len(POINT_MARKERS)]
             for i, s in enumerate(SAMPLE_ORDER)}

# ---- effective independent units --------------------------------------------
print("\n    EFFECTIVE INDEPENDENT UNITS")
print("    Treatment was assigned to ANIMALS. The number of independent units")
print("    is six, three per arm, whatever the cell counts say.")
for c in CONDITION_ORDER:
    sids = [s for s in SAMPLE_ORDER if COND_OF[s] == c]
    ns = {short_label(s): int(cells[s][STRUCT_COL].nunique()) for s in sids}
    tot = sum(len(cells[s]) for s in sids)
    print(f"      {c:<12} animals {len(sids)}   structures {sum(ns.values()):>3}"
          f"   cells {tot:>9,}")
    print(f"                   per animal: {ns}")

# ---- shared null, and its provenance ----------------------------------------
sub("Shared null provenance")
SHARED_NULL = {}
if os.path.exists(NN_NULL_TABLE):
    nulls = pd.read_csv(NN_NULL_TABLE)
    key = ["sample_id", "structure_id", "anchor", "target"]
    if all(c in nulls.columns for c in key + ["null_median_um"]):
        for _, r in nulls.iterrows():
            SHARED_NULL[(r["sample_id"], int(r["structure_id"]),
                         r["anchor"], r["target"])] = float(r["null_median_um"])
        loaded_keys = {(s, int(k)) for s in SAMPLE_ORDER
                       for k in cells[s][STRUCT_COL].unique()}
        tab_keys = {(s, k) for (s, k, _a, _t) in SHARED_NULL}
        match = len(loaded_keys & tab_keys)
        frac = match / max(len(loaded_keys), 1)
        print(f"    {match} of {len(loaded_keys)} loaded structure keys present "
              f"in table 53 ({100 * frac:.1f}%)")
        print(f"    loaded {len(SHARED_NULL)} shared per-structure nulls")
        if frac < 0.95:
            print("    MISMATCH. Table 53 is from a different structure revision.")
            if REFUSE_ON_VALIDATION_FAILURE:
                print("    REFUSING TO RUN.")
                sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)
    else:
        print(f"    WARNING: {NN_NULL_TABLE} lacks expected columns.")
else:
    print(f"    WARNING: {NN_NULL_TABLE} not found. Validation impossible.")

# ---- the global common phenotype set ----------------------------------------
sub("Common phenotype set for the balanced anchor draw, over ALL structure cells")
_rows = []
for s in SAMPLE_ORDER:
    d = cells[s]
    for k, g in d.groupby(STRUCT_COL):
        vc = g["pheno"].value_counts()
        for ph_name, n_ph in vc.items():
            _rows.append({"sample_id": s, "structure_id": int(k),
                          "pheno": ph_name, "n": int(n_ph)})
ph_counts = pd.DataFrame(_rows)
n_struct_total = int(ph_counts[["sample_id", "structure_id"]]
                     .drop_duplicates().shape[0])
_ok_ph = ph_counts.loc[ph_counts["n"] >= MIN_CELLS_FOR_BALANCE]
cover = (_ok_ph.groupby("pheno").size() / max(n_struct_total, 1)).sort_values(
    ascending=False)
COMMON_PHENOS = sorted(cover.loc[cover >= BALANCE_COMMON_MIN_FRAC].index)

print(f"    {'phenotype':<30}{'structures at threshold':>26}")
print("    " + "-" * 56)
for ph_name, fr in cover.items():
    mark = "  <-- common set" if fr >= BALANCE_COMMON_MIN_FRAC else ""
    print(f"    {ph_name:<30}{100 * fr:>24.1f}%{mark}")
print(f"\n    common set: {len(COMMON_PHENOS)} phenotypes over "
      f"{n_struct_total} structures")
print("    DERIVED OVER ALL STRUCTURE CELLS, core plus cuff, because the null")
print("    this script reconstructs is script 06's, and that null shuffles")
print("    labels over every cell in the structure. Restricting to core cells")
print("    would break the validation in section B.")
print()
print("    REVISION 2 CORRECTION. Revision 1 asserted here that this set is")
print("    'small and myeloid-weighted, exactly as script 08 found'. That was")
print("    false of this script's own output. Script 08 counts coverage over")
print("    CORE cells and got 4 phenotypes; this script counts over all")
print("    structure cells and got 8, including B cells and plasma cells,")
print("    because the cuff pushes sparse populations over the 10-cell bar.")
print("    Two different sets were sharing one name. Each is correct for its")
print("    own script, neither is 'the' common set, and a number from one must")
print("    never be quoted against the other.")
print("    The ADAPTIVE balanced draw is reported beside it either way, so the")
print("    two can be compared rather than one being trusted.")
write_csv(ph_counts, "80_phenotype_counts_per_structure.csv")


# %% Cell 4 - the three nulls, and the validation that gates them
# =============================================================================

banner("A - THREE NULLS FROM ONE PERMUTATION STREAM")

print("    For each structure, anchor and target: draw a fake target set of the")
print("    OBSERVED target size, then draw three fake anchor sets against that")
print("    same fake target set. Only the anchor draw differs.\n")

rows = []
for s in SAMPLE_ORDER:
    d = cells[s]
    for k, g in d.groupby(STRUCT_COL):
        xy = g[["x", "y"]].to_numpy(float)
        ph = g["pheno"].to_numpy()
        n = len(xy)
        if n < MIN_CELLS_IN_STRUCTURE:
            continue
        idx_by_pheno = {p: np.flatnonzero(ph == p) for p in np.unique(ph)}

        # anchor sets are fixed for this structure, so build them once
        anchors_here = {}
        for anchor in NN_ANCHORS:
            ai = np.flatnonzero(ph == anchor)
            if not len(ai):
                continue
            anchors_here[anchor] = (
                rng_perm.choice(ai, MAX_CELLS_PER_STRUCTURE, replace=False)
                if len(ai) > MAX_CELLS_PER_STRUCTURE else ai)
        if not anchors_here:
            continue

        # The permutation and the phenotype pools depend on the TARGET only,
        # not on the anchor, so they are built once per target and reused
        # across both anchors. Doing it per anchor doubled the work for no gain.
        for target in NN_TARGETS:
            t_idx = np.flatnonzero(ph == target)
            n_t = len(t_idx)
            if not n_t:
                continue
            live_anchors = {a: au for a, au in anchors_here.items()
                            if a != target}
            if not live_anchors:
                continue

            observed = {a: nn_median(xy[au], xy[t_idx])
                        for a, au in live_anchors.items()}
            acc = {a: {"uniform": [], "bal_adapt": [], "bal_common": []}
                   for a in live_anchors}
            n_ph_a, n_ph_c = [], []

            for _ in range(N_PERMUTATIONS):
                perm = rng_perm.permutation(n)
                ft = perm[:n_t]                       # fake targets
                rest = perm[n_t:]
                if not len(rest):
                    continue
                tree_xy = xy[ft]

                # boolean mask rather than a Python set: the set version was
                # O(n) per phenotype per permutation and would not finish on
                # 43109's larger structures.
                in_rest = np.zeros(n, dtype=bool)
                in_rest[rest] = True
                pools_all, pools_com = {}, {}
                for p_, ix in idx_by_pheno.items():
                    keep = ix[in_rest[ix]]
                    if len(keep) >= MIN_ANCHORS_PER_PHENO:
                        pools_all[p_] = keep
                        if p_ in COMMON_PHENOS:
                            pools_com[p_] = keep
                n_ph_a.append(len(pools_all))
                n_ph_c.append(len(pools_com))

                for a, au in live_anchors.items():
                    n_a = int(min(len(au), len(rest)))
                    if n_a < 1:
                        continue
                    fa = rng_perm.choice(rest, size=n_a, replace=False)
                    acc[a]["uniform"].append(nn_median(xy[fa], tree_xy))

                    fb = balanced_draw(rng_perm, pools_all, n_a)
                    if fb is not None and len(fb):
                        acc[a]["bal_adapt"].append(nn_median(xy[fb], tree_xy))

                    fc = balanced_draw(rng_perm, pools_com, n_a)
                    if fc is not None and len(fc):
                        acc[a]["bal_common"].append(nn_median(xy[fc], tree_xy))

            for a, au in live_anchors.items():
                def _med(key, _a=a):
                    v = [x for x in acc[_a][key] if np.isfinite(x)]
                    return float(np.median(v)) if v else np.nan
                nu, na_, nc = _med("uniform"), _med("bal_adapt"), _med("bal_common")
                obs = observed[a]
                shared = SHARED_NULL.get((s, int(k), a, target), np.nan)
                rows.append({
                    "sample_id": s, "animal_id": short_label(s),
                    "condition": COND_OF[s], "structure_id": int(k),
                    "anchor": a, "target": target,
                    "n_anchor": int(len(np.flatnonzero(ph == a))),
                    "n_anchors_used": int(len(au)),
                    "n_target": int(n_t), "n_structure": int(n),
                    "observed_median_um": obs,
                    "null_uniform_um": nu,
                    "null_bal_adapt_um": na_,
                    "null_bal_common_um": nc,
                    "n_phenos_adapt": int(np.median(n_ph_a)) if n_ph_a else 0,
                    "n_phenos_common": int(np.median(n_ph_c)) if n_ph_c else 0,
                    "table53_null_um": shared,
                    "delta_uniform_um": obs - nu if np.isfinite(nu) else np.nan,
                    "delta_bal_adapt_um": obs - na_ if np.isfinite(na_) else np.nan,
                    "delta_bal_common_um": obs - nc if np.isfinite(nc) else np.nan,
                })
    gc.collect()

nn = pd.DataFrame(rows)
if not len(nn):
    print("ERROR: no structure-pairs computed.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)
write_csv(nn, "81_three_nulls_per_structure.csv")


# ---- THE VALIDATION ---------------------------------------------------------
banner("B - VALIDATION: DOES null_uniform REPRODUCE TABLE 53?")

print("    Everything below this gate assumes the reconstruction is correct.")
print("    If it is not, the balanced nulls are meaningless and nothing in")
print("    section C may be read.\n")

VALIDATION_PASSED = False
v = nn.dropna(subset=["null_uniform_um", "table53_null_um"])
if not len(v):
    print("    CANNOT VALIDATE: no structure-pairs matched table 53.")
else:
    r = float(np.corrcoef(v["null_uniform_um"], v["table53_null_um"])[0, 1])
    diff = (v["null_uniform_um"] - v["table53_null_um"]).abs()
    med_d, p95_d = float(diff.median()), float(diff.quantile(0.95))
    print(f"    matched structure-pairs   : {len(v):,} of {len(nn):,}")
    print(f"    Pearson r                 : {r:.4f}   (threshold {VALIDATION_MIN_CORR})")
    print(f"    median |difference|       : {med_d:.2f} um   "
          f"(threshold {VALIDATION_MAX_MEDIAN_ABS_DIFF_UM})")
    print(f"    95th percentile |diff|    : {p95_d:.2f} um")
    VALIDATION_PASSED = bool(r >= VALIDATION_MIN_CORR
                             and med_d <= VALIDATION_MAX_MEDIAN_ABS_DIFF_UM)
    if VALIDATION_PASSED:
        print("\n    PASS. This script reproduces script 06's null, so the only")
        print("    thing that differs in the balanced versions is the anchor")
        print("    reference, which is what we set out to vary.")
    else:
        print("\n    FAIL. The reconstruction does not match table 53.")
        print("    The balanced nulls below are NOT interpretable, because a")
        print("    difference between them and script 07 could be the anchor")
        print("    reference or could be this script being wrong.")
        print("    Investigate before reading anything further.")

    fig, ax = plt.subplots(figsize=(15, 14))
    for s in SAMPLE_ORDER:
        g = v.loc[v["sample_id"] == s]
        if not len(g):
            continue
        ax.scatter(g["table53_null_um"], g["null_uniform_um"], s=200,
                   marker=MARKER_OF[s], alpha=0.65,
                   color=CONDITION_COLORS.get(COND_OF[s], "#777777"),
                   label=f"{short_label(s)} ({COND_OF[s]})")
    lim = [0, float(max(v["table53_null_um"].max(),
                        v["null_uniform_um"].max())) * 1.05]
    ax.plot(lim, lim, color="#000000", linewidth=3, linestyle="--", zorder=1)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("table 53 null median, um")
    ax.set_ylabel("null_uniform recomputed here, um")
    ax.set_title(f"Validation: r = {r:.4f}, median |diff| = {med_d:.2f} um",
                 fontsize=FONT_SIZE_TITLE)
    ax.legend(frameon=False, fontsize=FONT_SIZE_LEGEND - 12)
    style_axes(ax)
    save_fig(fig, "F60_validation_against_table53")

write_csv(v.assign(abs_diff_um=(v["null_uniform_um"] - v["table53_null_um"]).abs())
          if len(v) else pd.DataFrame(), "82_validation_detail.csv")


# %% Cell 5 - C: does the delta change when the reference stops being cell-weighted
# =============================================================================

banner("C - THE DELTA UNDER THREE ANCHOR REFERENCES")

if not VALIDATION_PASSED and REFUSE_ON_VALIDATION_FAILURE:
    print("    SUPPRESSED. The validation in section B did not pass, and this")
    print("    section cannot be interpreted without it. Tables 81 and 82 are")
    print("    written so the failure can be diagnosed. Stopping here.")
    banner("DONE (validation failed)")
    sys.stdout = _tee.terminal
    _tee.close()
    sys.exit(1)

OUTCOMES = [("delta_uniform_um", "uniform anchors (script 07's reference)"),
            ("delta_bal_adapt_um", "balanced, phenotypes in structure"),
            ("delta_bal_common_um", "balanced, global common set")]

# per-animal summary, cell-weighted to match scripts 06, 07 and 08
sub("Per-animal delta, cell-weighted across that animal's structures")
per_rows = []
for (a, t, s), g in nn.groupby(["anchor", "target", "sample_id"]):
    w = g["n_anchors_used"].to_numpy(float)
    if not w.sum():
        continue
    rec = {"anchor": a, "target": t, "sample_id": s,
           "animal_id": short_label(s), "condition": COND_OF[s],
           "n_anchors": int(w.sum()), "n_targets": int(g["n_target"].sum()),
           "n_targets_min_structure": int(g["n_target"].min()),
           "n_structures": int(len(g))}
    for col, _ in OUTCOMES:
        vv = g[col].to_numpy(float)
        ok = np.isfinite(vv)
        rec[col] = float(np.average(vv[ok], weights=w[ok])) if ok.any() else np.nan
    per_rows.append(rec)
per = pd.DataFrame(per_rows)
write_csv(per, "83_per_animal_delta_three_references.csv")

# reportability, matching script 07
min_treated_targets = {}
for (a, t), g in per.groupby(["anchor", "target"]):
    gt = g.loc[g["condition"] == "D1MT", "n_targets_min_structure"]
    min_treated_targets[(a, t)] = int(gt.min()) if len(gt) else 0

for col, name in OUTCOMES:
    sub(f"Per-animal delta, {name}")
    print(f"    {'anchor -> target':<44}"
          + "".join(f"{short_label(s):>11}" for s in SAMPLE_ORDER)
          + f"{'separated':>12}{'margin':>10}")
    print("    " + "-" * (44 + 11 * len(SAMPLE_ORDER) + 22))
    for a in NN_ANCHORS:
        for t in NN_TARGETS:
            g = per.loc[(per["anchor"] == a) & (per["target"] == t)]
            if not len(g):
                continue
            vals, tv, uv = [], [], []
            for s in SAMPLE_ORDER:
                r_ = g.loc[g["sample_id"] == s, col]
                x = float(r_.iloc[0]) if len(r_) else np.nan
                vals.append(x)
                (tv if COND_OF[s] == "D1MT" else uv).append(x)
            sep = complete_separation(tv, uv)
            mg = arm_margin(tv, uv)
            flag = ""
            if min_treated_targets.get((a, t), 0) < MIN_TREATED_TARGETS_FOR_REPORT:
                flag = "  [THIN]"
            label = f"{a.replace('Macrophages', 'Mac')} -> {t}"
            print(f"    {label[:43]:<44}"
                  + "".join(f"{x:>+11.2f}" if np.isfinite(x) else f"{'na':>11}"
                            for x in vals)
                  + f"{'YES' if sep else '':>12}"
                  + (f"{mg:>10.2f}" if np.isfinite(mg) else f"{'':>10}")
                  + flag)

# ---- models on each reference ----------------------------------------------
sub("Mixed models on each reference, with the exact animal-level p")
print("    The reference the model uses is the ONLY thing that differs between")
print("    these three rows. p_exact_means is floored at 0.10 by the design.\n")

percell = []
for s in SAMPLE_ORDER:
    d = cells[s]
    for k, g in d.groupby(STRUCT_COL):
        xy = g[["x", "y"]].to_numpy(float)
        ph = g["pheno"].to_numpy()
        if len(xy) < MIN_CELLS_IN_STRUCTURE:
            continue
        for a in NN_ANCHORS:
            ai = np.flatnonzero(ph == a)
            if not len(ai):
                continue
            au = (rng_perm.choice(ai, MODEL_MAX_CELLS_PER_STRUCTURE, replace=False)
                  if len(ai) > MODEL_MAX_CELLS_PER_STRUCTURE else ai)
            for t in NN_TARGETS:
                if t == a:
                    continue
                ti = np.flatnonzero(ph == t)
                if not len(ti):
                    continue
                row = nn.loc[(nn["sample_id"] == s) & (nn["structure_id"] == int(k))
                             & (nn["anchor"] == a) & (nn["target"] == t)]
                if not len(row):
                    continue
                dv = nn_all(xy[au], xy[ti])
                if not len(dv):
                    continue
                frame = {"sample_id": s, "condition": COND_OF[s], STRUCT_COL: k,
                         "anchor": a, "target": t}
                # The per-cell outcome is each anchor's own distance minus that
                # structure's null median. Table 81 stores the delta rather than
                # the null, so recover the null as observed - delta and use the
                # matching null column directly where it exists.
                for col, _ in OUTCOMES:
                    null_col = col.replace("delta_", "null_")
                    if null_col in row.columns and np.isfinite(row[null_col].iloc[0]):
                        frame[col] = dv - float(row[null_col].iloc[0])
                    else:
                        frame[col] = np.nan
                percell.append(pd.DataFrame(frame))
percell = pd.concat(percell, ignore_index=True) if percell else pd.DataFrame()

model_rows = []
if len(percell) and HAVE_SM:
    for a in NN_ANCHORS:
        for t in NN_TARGETS:
            if t == a:
                continue
            d = percell.loc[(percell["anchor"] == a) & (percell["target"] == t)]
            if not len(d):
                continue
            for col, name in OUTCOMES:
                r_ = fit_mixed(d, col, f"{a} -> {t}")
                if r_:
                    r_["anchor"], r_["target"] = a, t
                    r_["reference"] = name
                    r_["reportable"] = bool(
                        min_treated_targets.get((a, t), 0)
                        >= MIN_TREATED_TARGETS_FOR_REPORT)
                    model_rows.append(r_)

models = pd.DataFrame(model_rows)
if len(models):
    # BH within each reference separately. The three references are the SAME
    # hypothesis measured three ways, not three hypotheses, so they are never
    # pooled into one family. Matches the rule agreed for script 08.
    models["q_value"] = np.nan
    for ref, g in models.groupby("reference"):
        models.loc[g.index, "q_value"] = benjamini_hochberg(
            g["p_value"].to_numpy())
    write_csv(models, "84_models_three_references.csv")

    print(f"    {'pair':<40}{'reference':<40}{'coef um':>10}{'p':>9}"
          f"{'q':>9}{'p_exact':>9}")
    print("    " + "-" * 117)
    for a in NN_ANCHORS:
        for t in NN_TARGETS:
            g = models.loc[(models["anchor"] == a) & (models["target"] == t)]
            if not len(g):
                continue
            for _, r_ in g.iterrows():
                lab = f"{a.replace('Macrophages', 'Mac')} -> {t}"
                flag = "" if r_["reportable"] else "  [THIN]"
                print(f"    {lab[:39]:<40}{r_['reference'][:39]:<40}"
                      f"{r_['coef_D1MT_vs_ref']:>+10.2f}{r_['p_value']:>9.4f}"
                      f"{r_['q_value']:>9.4f}{r_['p_exact_means']:>9.3f}{flag}")


# %% Cell 6 - the verdict
# =============================================================================

banner("D - VERDICT")

print("    The question was whether script 07's plasma cell result is carried")
print("    by its cell-weighted anchor reference. Read the shift from")
print("    'uniform anchors' to the two balanced references, per pair.\n")

verdict_rows = []
for a in NN_ANCHORS:
    for t in NN_TARGETS:
        if t == a:
            continue
        g = per.loc[(per["anchor"] == a) & (per["target"] == t)]
        if not len(g):
            continue
        rec = {"anchor": a, "target": t,
               "n_targets_min_treated": min_treated_targets.get((a, t), 0)}
        rec["reportable"] = bool(rec["n_targets_min_treated"]
                                 >= MIN_TREATED_TARGETS_FOR_REPORT)
        for col, _ in OUTCOMES:
            tv = g.loc[g["condition"] == "D1MT", col].tolist()
            uv = g.loc[g["condition"] != "D1MT", col].tolist()
            rec[f"{col}_separated"] = complete_separation(tv, uv)
            rec[f"{col}_margin"] = arm_margin(tv, uv)
        verdict_rows.append(rec)
vd = pd.DataFrame(verdict_rows)
if len(vd):
    write_csv(vd, "85_separation_under_three_references.csv")

    print(f"    {'pair':<42}{'uniform':>12}{'bal adapt':>12}{'bal common':>12}"
          f"{'verdict':>34}")
    print("    " + "-" * 112)
    for _, r_ in vd.iterrows():
        lab = f"{r_['anchor'].replace('Macrophages', 'Mac')} -> {r_['target']}"
        su = r_["delta_uniform_um_separated"]
        sa = r_["delta_bal_adapt_um_separated"]
        sc = r_["delta_bal_common_um_separated"]
        if not r_["reportable"]:
            note = "THIN, not reportable"
        elif su and not (sa or sc):
            note = "LOST under both balanced refs"
        elif su and (sa and sc):
            note = "HOLDS under both"
        elif su and (sa or sc):
            note = "holds under one, lost under one"
        elif (not su) and (sa or sc):
            note = "APPEARS only when balanced"
        else:
            note = "no separation on any reference"
        print(f"    {lab[:41]:<42}{('YES' if su else '-'):>12}"
              f"{('YES' if sa else '-'):>12}{('YES' if sc else '-'):>12}"
              f"{note:>34}")

    rep = vd.loc[vd["reportable"]]
    lost = rep.loc[rep["delta_uniform_um_separated"]
                   & ~(rep["delta_bal_adapt_um_separated"]
                       | rep["delta_bal_common_um_separated"])]
    held = rep.loc[rep["delta_uniform_um_separated"]
                   & rep["delta_bal_adapt_um_separated"]
                   & rep["delta_bal_common_um_separated"]]
    print(f"\n    Of {len(rep)} reportable pairs separating on script 07's "
          f"reference:")
    print(f"      {len(held)} hold under both balanced references")
    print(f"      {len(lost)} are lost under both")
    if len(lost):
        print("\n    LOST: " + "; ".join(
            f"{r['anchor'].replace('Macrophages', 'Mac')} -> {r['target']}"
            for _, r in lost.iterrows()))
        print("    For these, the separation script 07 reports is a property of")
        print("    the cell-weighted anchor reference rather than of the")
        print("    IDO1-positive compartment. Script 07's confirmatory family")
        print("    needs a revision before those pairs are quoted.")
    if len(held):
        print("\n    HELD: " + "; ".join(
            f"{r['anchor'].replace('Macrophages', 'Mac')} -> {r['target']}"
            for _, r in held.iterrows()))
        print("    For these, the reference is not what carries the result.")

print("\n    HOW TO READ THIS AGAINST SCRIPT 08. Script 08 found B cells to be")
print("    the only reportable population surviving multiple testing once the")
print("    reference stopped being cell-weighted, and plasma cells to lose its")
print("    per-animal separation. If the same ordering appears above, the two")
print("    scripts agree after all and the disagreement was the reference. If")
print("    plasma cells holds here, they are measuring different things and")
print("    both results stand, which is the harder outcome to write up.")

banner("SUMMARY")
print(f"Input                     : {IN_DIR}")
print(f"Structure-pairs           : {len(nn):,}")
print(f"Permutations per null     : {N_PERMUTATIONS}")
print(f"Validation against t.53   : {'PASS' if VALIDATION_PASSED else 'FAIL'}")
print(f"Common set                : {len(COMMON_PHENOS)} phenotypes")
print(f"Models fitted             : {len(models) if len(models) else 0}")
print(f"Per-animal summary        : cell-weighted mean")
print(f"Reportability threshold   : {MIN_TREATED_TARGETS_FOR_REPORT} treated targets")
print("\nRead in this order")
print("  1. Section B, the validation. Nothing below it means anything if it")
print("     failed.")
print("  2. Section D, the verdict table, which answers the question.")
print("  3. Table 83 for the per-animal values behind every separation.")
print("  4. Table 84 for the models, remembering that the three references are")
print("     one hypothesis measured three ways and are never pooled for BH.")

banner("DONE")
sys.stdout = _tee.terminal
_tee.close()
