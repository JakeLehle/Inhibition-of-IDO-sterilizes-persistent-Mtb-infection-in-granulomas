#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - MIXED MODEL DIAGNOSTICS
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 10 of the AKOYA analysis series. READ-ONLY DIAGNOSTIC.

WHY THIS SCRIPT EXISTS

    Three focus-level outcomes fall back to a two-stage t-test on the six animal
    means because the linear mixed model returns no usable standard error: fold
    over background, IDO1+ percent of core macrophages, and shape ratio. The
    fallback is calibrated but it throws away the within-animal information and
    it puts a different test on the figure from the one used everywhere else.
    The decision is to move off it. This script establishes WHY the fits fail
    and WHICH route replaces them, with measured false positive rates rather
    than an argument.

    Nothing is modified. No table in structures_rev4, lymphocyte_radial or
    publication_figures is touched, and akoya_arm_stats.py is imported but not
    changed.

THE TWO COMPETING EXPLANATIONS, AND HOW THEY ARE TOLD APART

    A. DISTRIBUTIONAL. The outcome is skewed or heavy-tailed, a linear model on
       the raw scale is the wrong functional form, and a transform fixes it.

    B. IDENTIFIABILITY. Treated animals contribute 3, 1 and 2 foci. For the
       animal with a single focus the animal random effect and the residual are
       the same number, so there is no within-animal replication to separate
       them. The animal variance slides to the boundary at zero, the Hessian
       goes singular, and statsmodels returns a coefficient with a NaN standard
       error and raises nothing.

    These make different predictions and the script tests both. If A, the
    residual and outcome distributions will be visibly wrong on the raw scale
    and a log transform will produce a usable fit. If B, the REML profile
    likelihood will be monotonically increasing as the random-effect standard
    deviation goes to zero, the maximum will sit ON the boundary, and no
    transform will help because the covariance is singular because of the data.
    Cell 5 plots that profile directly, so this is settled by looking rather
    than by assuming. My expectation before running is B, and if the profile
    does not run to the boundary then that expectation was wrong and A is worth
    pursuing.

THE FOUR ROUTES COMPARED

    1. lmm            statsmodels MixedLM as currently used. Asymptotic z on six
                      clusters, so the optimistic end.
    2. lmm_penalised  the same random-intercept model with a weak prior on the
                      random-effect standard deviation that keeps it off the
                      boundary. This is what blme does in R, implemented here on
                      an exact profile REML for the random-intercept case so it
                      has no new dependency. Always identifiable, so it removes
                      the fallback entirely.
    3. lmm_exact      the mixed model as the test STATISTIC, with the p-value
                      from enumerating all C(6,3) = 20 assignments of animals to
                      arms. Every focus contributes to the estimate while the
                      calibration comes from the randomization rather than from
                      an asymptotic z on six clusters. Exact by construction.
    4. two_stage_t    the current fallback, carried for comparison only.

    Datta-Satten clustered rank-sum is also reported as the assumption-free
    reference.

THE FLOOR, STATED ONCE

    Treatment was assigned to six animals, three per arm. There are C(6,3) = 20
    assignments, so any test whose p-value comes from counting arrangements is
    floored at 1/20 one-sided and 2/20 two-sided. Routes 3 and 4 respect that
    floor. Routes 1 and 2 can go below it because their p comes from a
    likelihood, which is an assumption and not extra information. That is the
    whole trade and it should be stated in the Methods rather than resolved.

WHAT IS NOT DECIDED HERE

    The final model choice waits until the foci definition is settled in script
    04b. If the audit recovers additional treated foci, particularly for the
    animals now contributing one and two, the variance components become
    identifiable and several of these fits will simply work. Run this before
    and after that change.

INPUTS
    structures_rev4/tables/35_foci_structures_relative.csv
    structures_rev4/cell_assignments/*_cell_structures.csv   (conditions only)
    lymphocyte_radial/tables/72_models_all.csv               (optional)
    akoya_arm_stats.py

OUTPUTS
    model_diagnostics/tables/   100_variance_components, 101_profile_likelihood,
                                102_route_comparison, 103_null_calibration,
                                104_cluster_structure, 00_model_diagnostics_report.txt
    model_diagnostics/figures/  F80_profile_likelihood, F81_distributions,
                                F82_route_comparison, F83_null_calibration

USAGE
    conda activate sc_pre
    python AKOYA_10_Model_Diagnostics.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================

STRUCT_DIR = "/master/jlehle/WORKING/AKOYA/structures_rev4"
CELL_DIR = f"{STRUCT_DIR}/cell_assignments"
FOCI_TABLE = f"{STRUCT_DIR}/tables/35_foci_structures_relative.csv"
MODELS_TABLE = "/master/jlehle/WORKING/AKOYA/lymphocyte_radial/tables/72_models_all.csv"
OUT_DIR = "/master/jlehle/WORKING/AKOYA/model_diagnostics"

# Only used when __file__ is undefined, i.e. running cells in Spyder.
SCRIPT_DIR_FALLBACK = "/master/jlehle/WORKING/AKOYA/scripts"

USE_AGG = True

# ---- outcomes ---------------------------------------------------------------
# (column in table 35, label, fit on log10)
OUTCOMES = [
    ("peak_density",             "Peak myeloid density (cells/mm2)",   True),
    ("fold_over_background",     "Fold over own background",           False),
    ("max_inscribed_radius_um",  "Maximum inscribed radius (um)",      False),
    ("n_cells_core",             "Cells in the focus core",            True),
    ("pct_ido1_pos_of_mac_core", "IDO1+ % of core macrophages",        False),
    ("shape_ratio",              "Shape (equivalent / inscribed)",     False),
]

CONDITION_ORDER = ["D1MT", "Untreated"]     # coefficient is [0] minus [1]
CONDITION_COLORS = {"D1MT": "#2C7FB8", "Untreated": "#D95F02"}
CONDITION_LABELS = {"D1MT": "D1MT-treated", "Untreated": "Untreated"}

# ---- profile likelihood -----------------------------------------------------
# lambda = var(animal) / var(residual). The grid spans six orders of magnitude
# and includes exactly zero, which is the boundary the failing fits run to.
LAMBDA_GRID = None          # built in Cell 2
LAMBDA_MIN_LOG10 = -6.0
LAMBDA_MAX_LOG10 = 3.0
LAMBDA_N = 121

# blme's default covariance prior is gamma(shape, rate=0) on the random-effect
# SD, which adds (shape - 1) * log(sd) to the objective. shape 2.5 is blme's
# default and is weak: it only bites near the boundary.
PRIOR_SHAPE = 2.5

# ---- null calibration -------------------------------------------------------
RUN_NULL_SIMULATION = True
N_SIM = 1000                # raise to 5000 once the run time is known
SIM_ICC = [0.0, 0.25, 0.50, 0.75]   # true animal share of variance, no arm effect
SIM_LAMBDA_N = 41           # coarser grid inside the simulation, for speed
SIM_SEED = 0
NOMINAL_ALPHA = 0.05

# ---- appearance -------------------------------------------------------------
ROUTE_COLORS = {
    "lmm":           "#D73027",
    "lmm_penalised": "#1A9850",
    "lmm_exact":     "#4575B4",
    "two_stage_t":   "#999999",
    "datta_satten":  "#762A83",
}
ROUTE_LABELS = {
    "lmm":           "mixed model (statsmodels)",
    "lmm_penalised": "mixed model, weak prior on SD",
    "lmm_exact":     "mixed model, exact randomization p",
    "two_stage_t":   "two-stage t on animal means (current fallback)",
    "datta_satten":  "Datta-Satten clustered rank-sum",
}

GRID_COLOR = "#DDDDDD"
AXIS_COLOR = "#333333"
TEXT_COLOR = "#000000"

FONT_SIZE_BASE = 28
FONT_SIZE_TITLE = 32
FONT_SIZE_TICK = 28
FONT_SIZE_LEGEND = 28
FONT_SIZE_PANEL = 34

DPI = 300
SAVE_PDF = True
SAVE_PNG = True


# %% Cell 2 - imports and helpers
# =============================================================================

import os
import sys
import gc
import glob
import warnings
from datetime import datetime
from itertools import combinations

import numpy as np
import pandas as pd

import matplotlib
if USE_AGG:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from scipy import stats as sst

try:
    import statsmodels.formula.api as smf
    HAVE_SM = True
except Exception:
    HAVE_SM = False

_here = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() \
    else SCRIPT_DIR_FALLBACK
if _here not in sys.path:
    sys.path.insert(0, _here)
try:
    import akoya_arm_stats as aas
    HAVE_AAS = True
except Exception as e:
    HAVE_AAS = False
    _aas_err = e

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

LAMBDA_GRID = np.concatenate([
    [0.0], np.logspace(LAMBDA_MIN_LOG10, LAMBDA_MAX_LOG10, LAMBDA_N)])
SIM_LAMBDA_GRID = np.concatenate([
    [0.0], np.logspace(LAMBDA_MIN_LOG10, LAMBDA_MAX_LOG10, SIM_LAMBDA_N)])


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
    for flag, ext in ((SAVE_PDF, "pdf"), (SAVE_PNG, "png")):
        if not flag:
            continue
        p = os.path.join(FIG_DIR, f"{stem}.{ext}")
        fig.savefig(p, dpi=DPI, bbox_inches="tight")
        print(f"    wrote {p}")
    plt.close(fig); gc.collect()


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


def panel_letter(ax, letter, dx=-0.16, dy=1.07):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=FONT_SIZE_PANEL,
            fontweight="bold", va="top", ha="left")


def find_col(df, options):
    for c in options:
        if c in df.columns:
            return c
    low = {str(c).lower(): c for c in df.columns}
    for c in options:
        if str(c).lower() in low:
            return low[str(c).lower()]
    return None


# =============================================================================
# EXACT PROFILE REML FOR THE RANDOM-INTERCEPT MODEL
# =============================================================================
# VALIDATED against statsmodels MixedLM before this script was used. On healthy
# data, 10 clusters of 8 with a real animal variance, the fixed-effect
# coefficient agrees to 3e-16 and the standard error to 1.5e-3 on an SE of 0.47,
# which is 0.3 percent and comes from optimiser tolerance rather than from the
# formulation. On the real cluster structure with the animal variance truly
# zero, this profile returns lambda_hat = 0 and flags the boundary, while
# statsmodels raises LinAlgError: Singular matrix from np.linalg.inv on
# X'V^-1X. That is the failure mode reproduced exactly, and it is caught by the
# try / except in statsmodels_lmm rather than being allowed to stop the run.
#
# y = X beta + Z b + e,  b ~ N(0, s2b I),  e ~ N(0, s2e I)
# With one random intercept per cluster, V_i = s2e (I + lam J_i) where
# lam = s2b / s2e. Then |I + lam J_i| = 1 + n_i lam and
# (I + lam J_i)^-1 = I - lam/(1 + n_i lam) J_i, both in closed form, so the
# profile REML log-likelihood can be evaluated exactly on a grid of lam with no
# optimiser and no convergence question. That is the point: it lets us SEE
# whether the maximum sits on the boundary at lam = 0 instead of inferring it
# from a failed fit.

def _reml_at_lambda(y, X, groups, lam):
    n, p = X.shape
    XtVX = np.zeros((p, p))
    XtVy = np.zeros(p)
    ytVy = 0.0
    logdetA = 0.0
    for g in np.unique(groups):
        m = groups == g
        ni = int(m.sum())
        Xi = X[m]; yi = y[m]
        c = lam / (1.0 + ni * lam) if lam > 0 else 0.0
        sX = Xi.sum(axis=0)
        sy = float(yi.sum())
        XtVX += Xi.T @ Xi - c * np.outer(sX, sX)
        XtVy += Xi.T @ yi - c * sX * sy
        ytVy += float(yi @ yi) - c * sy * sy
        logdetA += np.log1p(ni * lam)
    try:
        XtVX_inv = np.linalg.inv(XtVX)
    except np.linalg.LinAlgError:
        return None
    beta = XtVX_inv @ XtVy
    rss = ytVy - float(XtVy @ beta)
    dfree = n - p
    if not np.isfinite(rss) or rss <= 0 or dfree <= 0:
        return None
    s2e = rss / dfree
    sign, ld_xvx = np.linalg.slogdet(XtVX)
    if sign <= 0:
        return None
    ll = -0.5 * (dfree * np.log(2.0 * np.pi * s2e) + logdetA + ld_xvx + dfree)
    cov = s2e * XtVX_inv
    return dict(lam=float(lam), beta=beta, cov=cov, s2e=float(s2e),
                s2b=float(lam * s2e), reml_loglik=float(ll), df=int(dfree),
                logdetA=float(logdetA))


def profile_reml(y, X, groups, grid=None, prior_shape=1.0):
    """
    Evaluate profile REML across the lambda grid and return the maximum.

    prior_shape > 1 adds (prior_shape - 1) * log(sd_b) to the objective, which
    is blme's gamma covariance prior. It is weak everywhere except near sd_b = 0
    where it goes to minus infinity, so the maximum can no longer sit exactly on
    the boundary. This is the whole mechanism: the fit becomes identifiable
    without changing the estimand.
    """
    grid = LAMBDA_GRID if grid is None else grid
    recs, obj = [], []
    for lam in grid:
        r = _reml_at_lambda(y, X, groups, lam)
        if r is None:
            recs.append(None); obj.append(-np.inf); continue
        o = r["reml_loglik"]
        if prior_shape > 1.0:
            sd_b = np.sqrt(max(r["s2b"], 0.0))
            o = (o + (prior_shape - 1.0) * np.log(sd_b)) if sd_b > 0 else -np.inf
        recs.append(r); obj.append(o)
    obj = np.asarray(obj, dtype=float)
    if not np.isfinite(obj).any():
        return None, None, None
    k = int(np.nanargmax(obj))
    return recs[k], obj, np.asarray(grid, dtype=float)


def fit_profile(y, X, groups, prior_shape=1.0, grid=None, coef_index=1):
    """Point estimate, SE, Wald z and t on (n_clusters - 2) df."""
    best, obj, grid = profile_reml(y, X, groups, grid=grid,
                                   prior_shape=prior_shape)
    if best is None:
        return None
    se = float(np.sqrt(max(best["cov"][coef_index, coef_index], 0.0)))
    coef = float(best["beta"][coef_index])
    n_clust = int(len(np.unique(groups)))
    out = dict(coef=coef, se=se, lam=best["lam"], s2b=best["s2b"],
               s2e=best["s2e"], reml_loglik=best["reml_loglik"],
               n_obs=int(len(y)), n_clusters=n_clust,
               on_boundary=bool(best["lam"] <= 0.0),
               icc=(best["s2b"] / (best["s2b"] + best["s2e"])
                    if (best["s2b"] + best["s2e"]) > 0 else np.nan))
    if se > 0 and np.isfinite(se):
        z = coef / se
        out["z"] = z
        out["p_wald_z"] = float(2.0 * sst.norm.sf(abs(z)))
        dfree = max(1, n_clust - 2)
        out["p_wald_t"] = float(2.0 * sst.t.sf(abs(z), dfree))
        crit = float(sst.t.ppf(0.975, dfree))
        out["ci_low"] = coef - crit * se
        out["ci_high"] = coef + crit * se
    else:
        out.update(z=np.nan, p_wald_z=np.nan, p_wald_t=np.nan,
                   ci_low=np.nan, ci_high=np.nan)
    return out


def design(arm, cluster_order=None):
    """Intercept plus a 0/1 indicator for CONDITION_ORDER[0]."""
    a = np.asarray(arm)
    x1 = (a == CONDITION_ORDER[0]).astype(float)
    return np.column_stack([np.ones(len(a)), x1])


def exact_randomization_p(y, groups, arm_of_cluster, prior_shape=1.0,
                          grid=None, statistic="z"):
    """
    Enumerate every assignment of clusters to arms, refit, and compare.

    The arm label is permuted at the ANIMAL level because that is the unit
    treatment was assigned to. With 3 v 3 there are C(6,3) = 20 assignments,
    which come in 10 sign-symmetric pairs, so the two-sided p is floored at
    2/20 = 0.10. That floor is a property of the design and no test escapes it
    without an assumption.
    """
    clusters = list(arm_of_cluster.keys())
    n_c = len(clusters)
    n_a = sum(1 for c in clusters if arm_of_cluster[c] == CONDITION_ORDER[0])
    obs_arm = np.array([arm_of_cluster[g] for g in groups])
    obs = fit_profile(y, design(obs_arm), groups, prior_shape=prior_shape,
                      grid=grid)
    if obs is None:
        return None
    key = "z" if (statistic == "z" and np.isfinite(obs.get("z", np.nan))) \
        else "coef"
    obs_stat = abs(obs[key])
    stats_all = []
    for idx in combinations(range(n_c), n_a):
        lab = {clusters[i]: (CONDITION_ORDER[0] if i in idx
                             else CONDITION_ORDER[1]) for i in range(n_c)}
        arm = np.array([lab[g] for g in groups])
        r = fit_profile(y, design(arm), groups, prior_shape=prior_shape,
                        grid=grid)
        if r is None or not np.isfinite(r.get(key, np.nan)):
            continue
        stats_all.append(abs(r[key]))
    if not stats_all:
        return None
    stats_all = np.asarray(stats_all)
    tol = 1e-9 * max(1.0, obs_stat)
    p = float(np.mean(stats_all >= obs_stat - tol))
    out = dict(obs)
    out.update(p_exact=p, n_assignments=int(len(stats_all)),
               floor_two_sided=2.0 / max(1, len(stats_all)),
               statistic_used=key)
    return out


def two_stage_t(y, groups, arm_of_cluster):
    g = pd.DataFrame({"y": y, "g": groups})
    m = g.groupby("g")["y"].mean()
    arms = pd.Series({k: arm_of_cluster[k] for k in m.index})
    a = m[arms == CONDITION_ORDER[0]].to_numpy()
    b = m[arms == CONDITION_ORDER[1]].to_numpy()
    if len(a) < 2 or len(b) < 2:
        return None
    t = sst.ttest_ind(a, b, equal_var=True)
    return dict(coef=float(a.mean() - b.mean()), p_value=float(t.pvalue),
                n_clusters_a=len(a), n_clusters_b=len(b))


def statsmodels_lmm(df, value_col, log10=False):
    if not HAVE_SM:
        return None
    d = df.dropna(subset=[value_col, "arm", "animal"]).copy()
    d["_y"] = np.log10(d[value_col].clip(lower=1e-9)) if log10 else d[value_col]
    d["_arm"] = (d["arm"] == CONDITION_ORDER[0]).astype(float)
    if d["animal"].nunique() < 3 or len(d) < 5:
        return None
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fit = smf.mixedlm("_y ~ _arm", d, groups=d["animal"]).fit()
        coef = float(fit.params.get("_arm", np.nan))
        se = float(fit.bse.get("_arm", np.nan))
        p = float(fit.pvalues.get("_arm", np.nan))
        ok = bool(np.isfinite(coef) and np.isfinite(se) and se > 0
                  and np.isfinite(p))
        return dict(coef=coef, se=se, p_value=p, usable=ok,
                    converged=bool(getattr(fit, "converged", False)),
                    group_var=float(np.ravel(fit.cov_re)[0])
                    if fit.cov_re is not None else np.nan,
                    scale=float(fit.scale))
    except Exception as e:
        return dict(coef=np.nan, se=np.nan, p_value=np.nan, usable=False,
                    converged=False, group_var=np.nan, scale=np.nan,
                    error=str(e))


_tee = Tee(os.path.join(TAB_DIR, "00_model_diagnostics_report.txt"))
sys.stdout = _tee

banner("AKOYA MIXED MODEL DIAGNOSTICS (script 10, read-only)")
print(f"Run time : {datetime.now().isoformat(timespec='seconds')}")
print(f"Output   : {OUT_DIR}")
print("\nNothing is modified. akoya_arm_stats.py is imported, not changed.")
if not HAVE_SM:
    print("\n    WARNING: statsmodels unavailable. The 'lmm' route is skipped.")
    print("    The profile REML routes do not need it and still run.")
if not HAVE_AAS:
    print(f"\n    WARNING: akoya_arm_stats not importable ({_aas_err}).")
    print("    The Datta-Satten reference column is skipped. Set")
    print("    SCRIPT_DIR_FALLBACK if running cell by cell in Spyder.")


# %% Cell 3 - load and describe the cluster structure
# =============================================================================

banner("CLUSTER STRUCTURE")

print("Treatment was assigned to six animals, three per arm. The arm contrast")
print("carries six independent pieces of information no matter how many foci")
print("or cells exist. Precision of each animal's estimate and resolution of")
print("the between-arm p are different quantities.\n")

if not os.path.exists(FOCI_TABLE):
    print(f"ERROR: table 35 not found at {FOCI_TABLE}")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

foci = pd.read_csv(FOCI_TABLE)
print(f"    loaded table 35  ({foci.shape[0]} x {foci.shape[1]})")

sid_col = find_col(foci, ["sample_id", "section", "sample"])
if sid_col is None:
    print("ERROR: table 35 has no sample identifier column.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)
foci["sample_id"] = foci[sid_col].astype(str)

cond_col = find_col(foci, ["condition", "arm", "treatment", "group"])
if cond_col is not None:
    foci["arm"] = foci[cond_col].astype(str)
else:
    print("    condition not in table 35, reading it from the cell assignments")
    cmap = {}
    for p in sorted(glob.glob(os.path.join(CELL_DIR, "*_cell_structures.csv"))):
        s = os.path.basename(p).replace("_cell_structures.csv", "")
        try:
            cmap[s] = pd.read_csv(p, usecols=["condition"],
                                  nrows=1)["condition"].iloc[0]
        except Exception as e:
            print(f"    WARNING: could not read condition for {s}: {e}")
    foci["arm"] = foci["sample_id"].map(cmap)

foci["animal"] = foci["sample_id"].map(short_label)
foci = foci.loc[foci["arm"].isin(CONDITION_ORDER)].copy()
if not len(foci):
    print("ERROR: no foci with a recognised arm label.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

struct = (foci.groupby(["arm", "animal"]).size().rename("n_foci")
          .reset_index().sort_values(["arm", "animal"]))
print()
print(struct.to_string(index=False))
write_csv(struct, "104_cluster_structure.csv")

singletons = struct.loc[struct["n_foci"] < 2]
print(f"\n    animals contributing a single focus: {len(singletons)}")
if len(singletons):
    print("    " + ", ".join(f"{r.animal} ({r.arm})" for r in singletons.itertuples()))
    print("\n    For these animals the random intercept and the residual are the")
    print("    same number. There is no within-animal replication to separate")
    print("    them, so the animal variance is not identifiable from their")
    print("    contribution and the likelihood has nothing to pull it off zero.")
    print("    This is the identifiability hypothesis, tested directly in Cell 5.")

ARM_OF = dict(zip(struct["animal"], struct["arm"]))
ANIMALS = list(struct["animal"])

prev = None
if os.path.exists(MODELS_TABLE):
    prev = pd.read_csv(MODELS_TABLE)
    print(f"\n    loaded table 72 for reference ({prev.shape[0]} rows)")
    rcol = find_col(prev, ["random_effects"])
    if rcol is not None:
        nfb = int(prev[rcol].astype(str).str.contains("fallback").sum())
        print(f"    rows currently using the fallback: {nfb} of {len(prev)}")


# %% Cell 4 - distribution diagnostics
# =============================================================================

banner("DISTRIBUTIONS: IS THE FUNCTIONAL FORM WRONG?")

print("Hypothesis A says the fits fail because a linear model on the raw scale")
print("is the wrong shape for these outcomes. If so the raw distributions will")
print("be badly skewed and the log scale will look far better. Read the")
print("skewness and the Shapiro p on both scales below, then look at F81.\n")

dist_rows = []
for col, label, use_log in OUTCOMES:
    if col not in foci.columns:
        print(f"    WARNING: {col} not in table 35, skipped")
        continue
    v = foci[col].to_numpy(dtype=float)
    v = v[np.isfinite(v)]
    if len(v) < 4:
        print(f"    WARNING: {col} has {len(v)} usable values, skipped")
        continue
    rec = dict(outcome=col, label=label, n=len(v),
               min=float(v.min()), median=float(np.median(v)),
               max=float(v.max()),
               skew_raw=float(sst.skew(v)),
               kurtosis_raw=float(sst.kurtosis(v)))
    try:
        rec["shapiro_p_raw"] = float(sst.shapiro(v).pvalue)
    except Exception:
        rec["shapiro_p_raw"] = np.nan
    pos = v[v > 0]
    if len(pos) >= 4:
        lv = np.log10(pos)
        rec["skew_log10"] = float(sst.skew(lv))
        rec["kurtosis_log10"] = float(sst.kurtosis(lv))
        try:
            rec["shapiro_p_log10"] = float(sst.shapiro(lv).pvalue)
        except Exception:
            rec["shapiro_p_log10"] = np.nan
        rec["n_nonpositive"] = int(len(v) - len(pos))
    else:
        rec.update(skew_log10=np.nan, kurtosis_log10=np.nan,
                   shapiro_p_log10=np.nan, n_nonpositive=int(len(v) - len(pos)))
    rec["log10_recommended"] = use_log
    dist_rows.append(rec)

dists = pd.DataFrame(dist_rows)
if len(dists):
    show = ["outcome", "n", "skew_raw", "skew_log10", "shapiro_p_raw",
            "shapiro_p_log10", "n_nonpositive"]
    with pd.option_context("display.width", 240):
        print(dists[show].to_string(index=False))
print("\n    A large raw skew that collapses on the log scale is evidence for")
print("    hypothesis A for that outcome. Similar skew on both scales, or a")
print("    non-positive count above zero that blocks the log, is not.")


# %% Cell 5 - the profile likelihood, the decisive diagnostic
# =============================================================================

banner("PROFILE LIKELIHOOD: WHERE DOES THE ANIMAL VARIANCE GO?")

print("lambda is var(animal) / var(residual). The profile REML log-likelihood")
print("is evaluated exactly at each lambda, with no optimiser, so there is no")
print("convergence question to confound the reading.")
print()
print("If the maximum sits at lambda = 0 the animal variance estimate is ON THE")
print("BOUNDARY. That is hypothesis B. statsmodels reports a singular random")
print("effects covariance there and returns a NaN standard error without")
print("raising, which is exactly the failure we saw. No transform fixes it,")
print("because the covariance is singular because of the DATA rather than the")
print("parameterisation.\n")

prof_rows, vc_rows, prof_curves = [], [], {}

for col, label, use_log in OUTCOMES:
    if col not in foci.columns:
        continue
    d = foci.dropna(subset=[col, "arm", "animal"]).copy()
    y_raw = d[col].to_numpy(dtype=float)
    if use_log:
        if (y_raw <= 0).any():
            print(f"    NOTE [{col}]: {int((y_raw <= 0).sum())} non-positive "
                  f"values, log10 skipped for this outcome")
            y = y_raw; scale = "raw"
        else:
            y = np.log10(y_raw); scale = "log10"
    else:
        y = y_raw; scale = "raw"
    groups = d["animal"].to_numpy()
    X = design(d["arm"].to_numpy())
    if len(np.unique(groups)) < 3 or len(y) < 5:
        print(f"    WARNING [{col}]: too few observations or clusters, skipped")
        continue

    best, obj, grid = profile_reml(y, X, groups, prior_shape=1.0)
    best_p, obj_p, _ = profile_reml(y, X, groups, prior_shape=PRIOR_SHAPE)
    if best is None:
        print(f"    WARNING [{col}]: profile REML returned nothing, skipped")
        continue
    prof_curves[col] = dict(grid=grid, obj=obj, obj_prior=obj_p,
                            lam_hat=best["lam"],
                            lam_hat_prior=best_p["lam"] if best_p else np.nan,
                            label=label, scale=scale)

    fit_free = fit_profile(y, X, groups, prior_shape=1.0)
    fit_pen = fit_profile(y, X, groups, prior_shape=PRIOR_SHAPE)
    sm = statsmodels_lmm(d.assign(**{col: d[col]}), col, log10=(scale == "log10"))

    on_b = bool(best["lam"] <= 0.0)
    vc_rows.append(dict(
        outcome=col, label=label, scale=scale, n_foci=len(y),
        n_animals=int(len(np.unique(groups))),
        lambda_hat=best["lam"], var_animal=best["s2b"], var_resid=best["s2e"],
        icc=fit_free["icc"] if fit_free else np.nan,
        on_boundary=on_b,
        lambda_hat_with_prior=best_p["lam"] if best_p else np.nan,
        icc_with_prior=fit_pen["icc"] if fit_pen else np.nan,
        statsmodels_usable=(sm["usable"] if sm else np.nan),
        statsmodels_group_var=(sm["group_var"] if sm else np.nan),
        statsmodels_se=(sm["se"] if sm else np.nan)))

    print(f"    {col:<28} scale {scale:<6} lambda_hat "
          f"{best['lam']:>10.5f}  ICC "
          f"{(fit_free['icc'] if fit_free else np.nan):>6.3f}  "
          f"boundary {'YES' if on_b else 'no ':<3}  "
          f"statsmodels usable "
          f"{('yes' if sm and sm['usable'] else 'NO ') if sm else 'n/a'}")

vc = pd.DataFrame(vc_rows)
if len(vc):
    write_csv(vc, "100_variance_components.csv")

if len(prof_curves):
    rows = []
    for col, c in prof_curves.items():
        for lam, o, op in zip(c["grid"], c["obj"], c["obj_prior"]):
            rows.append(dict(outcome=col, lam=lam, reml_loglik=o,
                             reml_loglik_with_prior=op))
    write_csv(pd.DataFrame(rows), "101_profile_likelihood.csv")

sub("VERDICT ON THE TWO HYPOTHESES")
if len(vc):
    nb = int(vc["on_boundary"].sum())
    print(f"    {nb} of {len(vc)} outcomes have the animal variance ON the")
    print(f"    boundary at zero.")
    if nb:
        names = ", ".join(vc.loc[vc["on_boundary"], "outcome"])
        print(f"    On the boundary: {names}")
        print("\n    For these the failure is IDENTIFIABILITY, not functional")
        print("    form. Check F81 anyway, but a transform will not rescue them.")
    if nb < len(vc):
        names = ", ".join(vc.loc[~vc["on_boundary"], "outcome"])
        print(f"    Off the boundary: {names}")
        print("    These fit normally and are the control showing the profile")
        print("    machinery is working.")
    bad = vc.loc[(vc["statsmodels_usable"] == False) & (~vc["on_boundary"])]  # noqa: E712
    if len(bad):
        print(f"\n    NOTE: {len(bad)} outcome(s) failed in statsmodels while")
        print("    the profile finds an interior maximum. That points at the")
        print("    optimiser rather than at the data, which would be a third")
        print("    explanation neither hypothesis covers. Names: "
              + ", ".join(bad["outcome"]))


# %% Cell 6 - route comparison on the real data
# =============================================================================

banner("ROUTE COMPARISON")

print("Same outcome, same data, four routes plus the assumption-free rank-sum.")
print("Read the p-values ACROSS a row, not down a column. Where they disagree")
print("the disagreement is the assumption, and it should be reported rather")
print("than resolved by picking the smallest.\n")

route_rows = []
for col, label, use_log in OUTCOMES:
    if col not in foci.columns:
        continue
    d = foci.dropna(subset=[col, "arm", "animal"]).copy()
    y_raw = d[col].to_numpy(dtype=float)
    scale = "raw"
    if use_log and not (y_raw <= 0).any():
        y = np.log10(y_raw); scale = "log10"
    else:
        y = y_raw
    groups = d["animal"].to_numpy()
    arm = d["arm"].to_numpy()
    X = design(arm)
    if len(np.unique(groups)) < 4 or len(y) < 6:
        print(f"    WARNING [{col}]: too small for the route comparison")
        continue

    rec = dict(outcome=col, label=label, scale=scale, n_foci=len(y),
               n_animals=int(len(np.unique(groups))))

    sm = statsmodels_lmm(d, col, log10=(scale == "log10"))
    if sm:
        rec.update(lmm_coef=sm["coef"], lmm_p=(sm["p_value"] if sm["usable"]
                                               else np.nan),
                   lmm_usable=sm["usable"])
    else:
        rec.update(lmm_coef=np.nan, lmm_p=np.nan, lmm_usable=False)

    f_free = fit_profile(y, X, groups, prior_shape=1.0)
    if f_free:
        rec.update(profile_coef=f_free["coef"], profile_se=f_free["se"],
                   profile_p_z=f_free["p_wald_z"], profile_p_t=f_free["p_wald_t"],
                   profile_on_boundary=f_free["on_boundary"])

    f_pen = fit_profile(y, X, groups, prior_shape=PRIOR_SHAPE)
    if f_pen:
        rec.update(penalised_coef=f_pen["coef"], penalised_se=f_pen["se"],
                   penalised_p_z=f_pen["p_wald_z"],
                   penalised_p_t=f_pen["p_wald_t"],
                   penalised_ci_low=f_pen["ci_low"],
                   penalised_ci_high=f_pen["ci_high"],
                   penalised_icc=f_pen["icc"])

    ex = exact_randomization_p(y, groups, ARM_OF, prior_shape=PRIOR_SHAPE)
    if ex:
        rec.update(exact_p=ex["p_exact"], exact_n_assignments=ex["n_assignments"],
                   exact_floor=ex["floor_two_sided"],
                   exact_statistic=ex["statistic_used"])

    ts = two_stage_t(y, groups, ARM_OF)
    if ts:
        rec.update(two_stage_coef=ts["coef"], two_stage_p=ts["p_value"])

    if HAVE_AAS:
        try:
            ds = aas.datta_satten_ranksum(y, groups, arm, CONDITION_ORDER[0],
                                          CONDITION_ORDER[1])
            rec["datta_satten_p"] = ds.get("p_two_sided", np.nan)
        except Exception as e:
            rec["datta_satten_p"] = np.nan
            print(f"    WARNING [{col}]: rank-sum failed ({e})")

    route_rows.append(rec)
    print(f"  {col}")
    print(f"    lmm            coef {rec.get('lmm_coef', np.nan):>+9.4f}  "
          f"p {rec.get('lmm_p', np.nan):>7.4f}"
          f"{'' if rec.get('lmm_usable') else '   [NO USABLE SE]'}")
    print(f"    profile REML   coef {rec.get('profile_coef', np.nan):>+9.4f}  "
          f"p {rec.get('profile_p_t', np.nan):>7.4f}  "
          f"boundary {'YES' if rec.get('profile_on_boundary') else 'no'}")
    print(f"    penalised      coef {rec.get('penalised_coef', np.nan):>+9.4f}  "
          f"p {rec.get('penalised_p_t', np.nan):>7.4f}  "
          f"ICC {rec.get('penalised_icc', np.nan):>5.3f}")
    print(f"    exact randomi. "
          f"                p {rec.get('exact_p', np.nan):>7.4f}  "
          f"floor {rec.get('exact_floor', np.nan):.3f}")
    print(f"    two-stage t    coef {rec.get('two_stage_coef', np.nan):>+9.4f}  "
          f"p {rec.get('two_stage_p', np.nan):>7.4f}")
    print(f"    rank-sum       "
          f"                p {rec.get('datta_satten_p', np.nan):>7.4f}\n")

routes = pd.DataFrame(route_rows)
if len(routes):
    write_csv(routes, "102_route_comparison.csv")


# %% Cell 7 - null calibration
# =============================================================================

banner("NULL CALIBRATION")

print("Every route above is simulated under a TRUE NULL on the real cluster")
print("structure, with no arm effect, across a range of true animal ICC. A")
print("route is usable only if its false positive rate is at or below the")
print("nominal 5 percent. This check comes BEFORE adopting a route, not after,")
print("which was the error made previously with the two-stage fallback.\n")

null_rows = []
if RUN_NULL_SIMULATION:
    sizes = struct.set_index("animal")["n_foci"].to_dict()
    groups_sim = np.concatenate([[a] * int(sizes[a]) for a in ANIMALS])
    arm_sim = np.array([ARM_OF[a] for a in groups_sim])
    X_sim = design(arm_sim)
    rs = np.random.default_rng(SIM_SEED)

    print(f"    cluster sizes : "
          + ", ".join(f"{a}={int(sizes[a])}" for a in ANIMALS))
    print(f"    simulations   : {N_SIM} per ICC")
    print(f"    ICC values    : {SIM_ICC}")
    # Timed on this cluster structure: the exact randomization route refits the
    # model 20 times per simulation and costs about 150 ms on a 41-point grid,
    # so it dominates. Everything else is a rounding error next to it.
    est_min = N_SIM * len(SIM_ICC) * 0.16 / 60.0
    print(f"    estimated run : about {est_min:.0f} minutes, dominated by the")
    print(f"                    exact route's 20 refits per simulation. Lower")
    print(f"                    N_SIM or SIM_LAMBDA_N if that is too long.\n")

    for icc in SIM_ICC:
        hits = {k: 0 for k in ["lmm", "profile_t", "profile_z",
                               "lmm_penalised", "lmm_exact", "two_stage_t",
                               "datta_satten"]}
        usable = {k: 0 for k in hits}
        sd_b = np.sqrt(icc)
        sd_e = np.sqrt(max(1.0 - icc, 1e-9))
        for _ in range(N_SIM):
            b = {a: rs.normal(0.0, sd_b) for a in ANIMALS}
            y = np.array([b[g] for g in groups_sim]) + rs.normal(
                0.0, sd_e, size=len(groups_sim))

            if HAVE_SM:
                dd = pd.DataFrame({"v": y, "arm": arm_sim, "animal": groups_sim})
                sm = statsmodels_lmm(dd, "v", log10=False)
                if sm and sm["usable"]:
                    usable["lmm"] += 1
                    hits["lmm"] += int(sm["p_value"] < NOMINAL_ALPHA)

            f = fit_profile(y, X_sim, groups_sim, prior_shape=1.0,
                            grid=SIM_LAMBDA_GRID)
            if f and np.isfinite(f["p_wald_t"]):
                usable["profile_t"] += 1
                hits["profile_t"] += int(f["p_wald_t"] < NOMINAL_ALPHA)
                usable["profile_z"] += 1
                hits["profile_z"] += int(f["p_wald_z"] < NOMINAL_ALPHA)

            fp = fit_profile(y, X_sim, groups_sim, prior_shape=PRIOR_SHAPE,
                             grid=SIM_LAMBDA_GRID)
            if fp and np.isfinite(fp["p_wald_t"]):
                usable["lmm_penalised"] += 1
                hits["lmm_penalised"] += int(fp["p_wald_t"] < NOMINAL_ALPHA)

            ex = exact_randomization_p(y, groups_sim, ARM_OF,
                                       prior_shape=PRIOR_SHAPE,
                                       grid=SIM_LAMBDA_GRID)
            if ex and np.isfinite(ex["p_exact"]):
                usable["lmm_exact"] += 1
                hits["lmm_exact"] += int(ex["p_exact"] < NOMINAL_ALPHA)

            ts = two_stage_t(y, groups_sim, ARM_OF)
            if ts:
                usable["two_stage_t"] += 1
                hits["two_stage_t"] += int(ts["p_value"] < NOMINAL_ALPHA)

            if HAVE_AAS:
                try:
                    ds = aas.datta_satten_ranksum(y, groups_sim, arm_sim,
                                                  CONDITION_ORDER[0],
                                                  CONDITION_ORDER[1])
                    usable["datta_satten"] += 1
                    hits["datta_satten"] += int(
                        ds.get("p_two_sided", 1.0) < NOMINAL_ALPHA)
                except Exception:
                    pass

        for k in hits:
            n = usable[k]
            null_rows.append(dict(
                true_icc=icc, route=k, n_usable=n, n_sim=N_SIM,
                fpr=(hits[k] / n) if n else np.nan,
                usable_fraction=n / N_SIM))
        print(f"    ICC {icc:.2f}   " + "   ".join(
            f"{k} {(hits[k] / usable[k] if usable[k] else np.nan):.3f}"
            for k in ["lmm", "lmm_penalised", "lmm_exact", "two_stage_t"]))

    nulls = pd.DataFrame(null_rows)
    write_csv(nulls, "103_null_calibration.csv")

    sub("VERDICT")
    piv = nulls.pivot_table(index="route", columns="true_icc", values="fpr")
    print(piv.to_string())
    worst = nulls.groupby("route")["fpr"].max().sort_values()
    print("\n    Worst-case false positive rate across the ICC range:")
    for r, v in worst.items():
        flag = "OK" if v <= NOMINAL_ALPHA + 0.01 else "ANTICONSERVATIVE"
        print(f"      {r:<16} {v:.3f}   {flag}")
    print("\n    Adopt the route with the largest power among those marked OK.")
    print("    A route that is anticonservative here cannot be rescued by")
    print("    reporting an effect size next to it.")
else:
    nulls = pd.DataFrame()
    print("    Skipped, RUN_NULL_SIMULATION is False.")


# %% Cell 8 - figures
# =============================================================================

banner("FIGURES")

# ---- F80 profile likelihood -------------------------------------------------
if prof_curves:
    keys = list(prof_curves.keys())
    ncol = 3
    nrow = int(np.ceil(len(keys) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(11 * ncol, 9 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for ax, col in zip(axes, keys):
        c = prof_curves[col]
        g, o, op = c["grid"], np.asarray(c["obj"]), np.asarray(c["obj_prior"])
        gpos = g.copy()
        gpos[gpos <= 0] = np.nan
        ok = np.isfinite(o)
        ax.plot(gpos[ok], (o - np.nanmax(o))[ok], linewidth=5,
                color=ROUTE_COLORS["lmm"], label="REML profile")
        okp = np.isfinite(op)
        if okp.any():
            ax.plot(gpos[okp], (op - np.nanmax(op[okp]))[okp], linewidth=5,
                    linestyle="--", color=ROUTE_COLORS["lmm_penalised"],
                    label="with weak prior")
        ax.set_xscale("symlog", linthresh=10 ** LAMBDA_MIN_LOG10)
        ax.axvline(10 ** LAMBDA_MIN_LOG10, color="#000000", linewidth=3)
        ax.text(10 ** LAMBDA_MIN_LOG10, ax.get_ylim()[0],
                " boundary\n at zero", va="bottom", ha="left",
                fontsize=FONT_SIZE_TICK - 8)
        if np.isfinite(c["lam_hat"]) and c["lam_hat"] > 0:
            ax.axvline(c["lam_hat"], color=ROUTE_COLORS["lmm"],
                       linestyle=":", linewidth=3)
        ax.set_title(col, fontsize=FONT_SIZE_BASE)
        ax.set_xlabel(r"$\lambda$ = var(animal) / var(residual)")
        ax.set_ylabel("REML log-likelihood\n(relative to maximum)")
        ax.legend(frameon=False, fontsize=FONT_SIZE_LEGEND - 10,
                  loc="lower right")
        style_axes(ax)
    for ax in axes[len(keys):]:
        ax.axis("off")
    fig.suptitle("Where the animal variance goes\n"
                 "a curve rising monotonically toward zero means the estimate "
                 "is on the boundary and the fit is not identifiable",
                 fontsize=FONT_SIZE_TITLE, y=1.01)
    fig.tight_layout()
    save_fig(fig, "F80_profile_likelihood")

# ---- F81 distributions ------------------------------------------------------
present = [(c, l, u) for c, l, u in OUTCOMES if c in foci.columns]
if present:
    ncol = 3
    nrow = int(np.ceil(len(present) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(11 * ncol, 9 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for ax, (col, label, use_log) in zip(axes, present):
        for arm in CONDITION_ORDER:
            v = foci.loc[foci["arm"] == arm, col].to_numpy(dtype=float)
            v = v[np.isfinite(v)]
            if not len(v):
                continue
            ax.scatter(v, np.full(len(v), CONDITION_ORDER.index(arm))
                       + np.random.default_rng(0).normal(0, 0.05, len(v)),
                       s=200, alpha=0.6, color=CONDITION_COLORS[arm],
                       edgecolors="#000000", linewidths=1.5)
        if use_log and (foci[col].dropna() > 0).all():
            ax.set_xscale("log")
        ax.set_yticks(range(len(CONDITION_ORDER)))
        ax.set_yticklabels([CONDITION_LABELS[c] for c in CONDITION_ORDER],
                           fontsize=FONT_SIZE_TICK - 6)
        ax.set_xlabel(label, fontsize=FONT_SIZE_BASE - 4)
        ax.set_title(col, fontsize=FONT_SIZE_BASE)
        style_axes(ax, ygrid=False)
    for ax in axes[len(present):]:
        ax.axis("off")
    fig.suptitle("Focus-level outcomes, one point per focus\n"
                 "the shape question, answered by looking",
                 fontsize=FONT_SIZE_TITLE, y=1.01)
    fig.tight_layout()
    save_fig(fig, "F81_distributions")

# ---- F82 route comparison ---------------------------------------------------
if len(routes):
    fig, ax = plt.subplots(figsize=(22, 13))
    cols = [("lmm_p", "lmm"), ("penalised_p_t", "lmm_penalised"),
            ("exact_p", "lmm_exact"), ("two_stage_p", "two_stage_t"),
            ("datta_satten_p", "datta_satten")]
    ys = np.arange(len(routes))
    off = np.linspace(-0.3, 0.3, len(cols))
    for (c, key), o in zip(cols, off):
        if c not in routes.columns:
            continue
        ax.scatter(routes[c], ys + o, s=420, color=ROUTE_COLORS[key],
                   edgecolors="#000000", linewidths=2,
                   label=ROUTE_LABELS[key], zorder=3)
    ax.axvline(0.05, color="#000000", linestyle="--", linewidth=3)
    ax.axvline(0.10, color="#000000", linestyle=":", linewidth=3)
    ax.text(0.10, len(routes) - 0.4, "  exact two-sided floor at 3 v 3",
            fontsize=FONT_SIZE_TICK - 8, va="center")
    ax.set_yticks(ys); ax.set_yticklabels(routes["outcome"],
                                          fontsize=FONT_SIZE_TICK - 6)
    ax.set_xscale("log")
    ax.set_xlabel("p-value")
    ax.set_title("Same data, five routes\n"
                 "where they disagree, the disagreement is the assumption",
                 fontsize=FONT_SIZE_TITLE)
    ax.legend(frameon=False, fontsize=FONT_SIZE_LEGEND - 10,
              loc="center left", bbox_to_anchor=(1.02, 0.5))
    style_axes(ax, ygrid=False)
    ax.xaxis.grid(True, color=GRID_COLOR, linewidth=1.5)
    save_fig(fig, "F82_route_comparison")

# ---- F83 null calibration ---------------------------------------------------
if len(nulls):
    fig, ax = plt.subplots(figsize=(20, 13))
    for route, g in nulls.groupby("route"):
        g = g.sort_values("true_icc")
        ax.plot(g["true_icc"], g["fpr"], marker="o", markersize=18,
                linewidth=4, color=ROUTE_COLORS.get(route, "#000000"),
                label=ROUTE_LABELS.get(route, route))
    ax.axhline(NOMINAL_ALPHA, color="#000000", linestyle="--", linewidth=3)
    ax.text(ax.get_xlim()[1], NOMINAL_ALPHA, "  nominal 5%", va="center",
            fontsize=FONT_SIZE_TICK - 6)
    ax.set_xlabel("True animal ICC (no arm effect)")
    ax.set_ylabel("False positive rate at p < 0.05")
    ax.set_title("Calibration on the real cluster structure\n"
                 "anything above the dashed line is anticonservative",
                 fontsize=FONT_SIZE_TITLE)
    ax.legend(frameon=False, fontsize=FONT_SIZE_LEGEND - 10,
              loc="center left", bbox_to_anchor=(1.02, 0.5))
    style_axes(ax)
    save_fig(fig, "F83_null_calibration")


# %% Cell 9 - wrap up
# =============================================================================

banner("SUMMARY")

if len(vc):
    print(f"Outcomes examined              : {len(vc)}")
    print(f"Animal variance on the boundary: {int(vc['on_boundary'].sum())}")
if len(nulls):
    ok = nulls.groupby("route")["fpr"].max()
    good = [r for r, v in ok.items() if v <= NOMINAL_ALPHA + 0.01]
    print(f"Routes calibrated at 5 percent : {', '.join(good) if good else 'none'}")

sub("Read in this order")
print("  1. F80 and table 100. Whether the animal variance sits on the")
print("     boundary. This decides identifiability against functional form and")
print("     it is the question that has been open.")
print("  2. F81 and the skewness columns. The shape question, for completeness")
print("     and because it is the thing a reviewer will ask about.")
print("  3. F83 and table 103. Which routes are calibrated at all. Nothing")
print("     that fails here is adoptable regardless of how well it fits.")
print("  4. F82 and table 102. Among the calibrated routes, which has the most")
print("     power on the real data.")

sub("What this run is for")
print("  Replacing the two-stage animal-means t-test with a single model used")
print("  everywhere, so that one test appears on every figure. The candidate is")
print("  the mixed model with a weak prior on the random-effect SD, which is")
print("  always identifiable, reported alongside an exact randomization p as")
print("  the assumption-free companion.")

sub("What this run does NOT settle")
print("  The design floor. Three animals per arm gives 20 assignments, so an")
print("  exact p cannot fall below 0.10 two-sided. A model-based p below that")
print("  is buying resolution with normality and asymptotics on six clusters,")
print("  which is a legitimate choice but is an assumption and not extra")
print("  information. Report the effect size and the separation first.")
print()
print("  The final choice also waits on script 04b. If the foci audit recovers")
print("  more treated foci, the animals now contributing one and two will gain")
print("  within-animal replication and several of these fits become")
print("  identifiable without any change of method. Rerun this after that.")

banner("DONE")
sys.stdout = _tee.terminal
_tee.close()
