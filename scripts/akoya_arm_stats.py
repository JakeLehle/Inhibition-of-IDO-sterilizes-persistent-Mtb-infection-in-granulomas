#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
akoya_arm_stats.py

Shared statistics and annotation helpers for the AKOYA analysis series.

WHY THIS IS A MODULE AND NOT COPIED INTO EACH SCRIPT
    The IDO1 cutoff plateau search was implemented twice, in scripts 03 and 09,
    and the two reported different windows for a year without anyone noticing.
    Anything used by more than one script lives here so that cannot happen
    again. Scripts 04, 08 and 09 import it.

WHAT IS IN HERE

    1. CLUSTERED RANK-SUM TEST (Datta and Satten, JASA 2005)
       A non-parametric test comparing two arms when the observations are cells
       or foci clustered within animals. Every observation is ranked jointly,
       ranks are averaged within cluster, and a rank-sum test is run on the
       cluster mean ranks. Unequal cluster sizes are handled, which is the point
       of the method.

       WHY NOT A PLAIN WILCOXON ON CELLS. Rank-sum handles unequal group sizes
       but it does not relax independence: its null is that all observations are
       exchangeable. Cells within one animal are not. Simulated on the real B
       lineage cell counts (122, 83, 869 treated; 1,491, 83, 778 untreated) with
       animals differing but NO arm effect, a plain cell-level Wilcoxon returns
       p < 0.05 in 85.6 percent of runs and can reach p = 4e-210 on pure noise.
       The clustered version, the animal-label permutation, and a rank-sum on
       animal means all returned 0.0 percent false positives.

       THE FLOOR. With three animals per arm there are only C(6,3) = 20 arm
       assignments, so the smallest achievable p is 1/20 = 0.05 one-sided and
       2/20 = 0.10 two-sided. That is a property of the design, not of the test.
       No amount of cells buys past it, because the randomisation happened at
       the animal.

    2. ONE-SIDED TESTING WITH PRE-SPECIFIED DIRECTION
       Every call must declare the direction expected from the biology BEFORE
       seeing the data, using D1MT_LOWER or D1MT_HIGHER. The result carries the
       declared direction, the one-sided p, the two-sided p, and a flag for
       whether the observed effect actually ran the declared way. If it ran the
       other way the one-sided p is reported as it comes out, near 1, and is not
       silently flipped. Both p-values are always in the returned record so
       nothing is hidden.

    3. CENTRED RADIAL POSITION, both centrings, identical to script 08.

    4. ANNOTATION HELPER that renders the same three-line bracket on any arm
       comparison: the model line, the clustered test line, and the sample size
       line stating animals rather than cells.

USAGE
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import akoya_arm_stats as aas

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

import numpy as np
import pandas as pd

try:
    from scipy import stats as _st
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

from itertools import combinations

# ---- pre-specified directions ----------------------------------------------
D1MT_LOWER = "D1MT_lower"     # treated expected below untreated
D1MT_HIGHER = "D1MT_higher"   # treated expected above untreated
NO_DIRECTION = "two_sided"    # no pre-specification, two-sided only


# =============================================================================
# clustered rank-sum
# =============================================================================

def per_cluster_summary(values, clusters, arms):
    """Mean value and arm for each cluster, in first-appearance order."""
    values = np.asarray(values, dtype=float)
    clusters = np.asarray(clusters)
    arms = np.asarray(arms)
    ok = np.isfinite(values)
    values, clusters, arms = values[ok], clusters[ok], arms[ok]
    out = []
    for c in pd.unique(clusters):
        m = clusters == c
        out.append({"cluster": c, "arm": arms[m][0], "n": int(m.sum()),
                    "mean": float(values[m].mean()),
                    "median": float(np.median(values[m]))})
    return pd.DataFrame(out)


def datta_satten_ranksum(values, clusters, arms, arm_a, arm_b,
                         direction=NO_DIRECTION):
    """
    Clustered rank-sum comparing arm_a against arm_b.

    values   : per-observation values (cells, foci, whatever the unit is)
    clusters : animal identifier for each observation
    arms     : arm label for each observation
    direction: D1MT_LOWER, D1MT_HIGHER or NO_DIRECTION. arm_a is treated as the
               D1MT side for the purpose of interpreting the direction.

    Returns a dict. p_value is the one-sided p when a direction was declared and
    the two-sided p otherwise; p_one_sided and p_two_sided are always present.
    """
    values = np.asarray(values, dtype=float)
    clusters = np.asarray(clusters)
    arms = np.asarray(arms)
    ok = np.isfinite(values)
    values, clusters, arms = values[ok], clusters[ok], arms[ok]

    rec = {"test": "Datta-Satten clustered rank-sum",
           "n_obs": int(len(values)), "direction_declared": direction,
           "p_one_sided": np.nan, "p_two_sided": np.nan, "p_value": np.nan,
           "n_clusters_a": 0, "n_clusters_b": 0,
           "mean_rank_a": np.nan, "mean_rank_b": np.nan,
           "ran_as_predicted": None, "floor_note": ""}
    if not len(values) or not HAVE_SCIPY:
        return rec

    ranks = _st.rankdata(values)
    rows = []
    for c in pd.unique(clusters):
        m = clusters == c
        rows.append((c, arms[m][0], float(ranks[m].mean())))
    cl = pd.DataFrame(rows, columns=["cluster", "arm", "mean_rank"])
    a = cl.loc[cl["arm"] == arm_a, "mean_rank"].to_numpy()
    b = cl.loc[cl["arm"] == arm_b, "mean_rank"].to_numpy()
    rec["n_clusters_a"], rec["n_clusters_b"] = len(a), len(b)
    if len(a) < 1 or len(b) < 1:
        return rec
    rec["mean_rank_a"], rec["mean_rank_b"] = float(a.mean()), float(b.mean())

    # EXACT permutation over cluster-label assignments.
    #
    # scipy's mannwhitneyu with method="auto" falls back to the asymptotic
    # approximation whenever the cluster mean ranks contain a tie, and the
    # asymptotic p can fall BELOW the exact design floor. Ties are not exotic
    # here: an animal contributing one focus has an integer mean rank and an
    # animal contributing two has an integer or half-integer one, so exact ties
    # do occur. Enumerating the assignments removes the problem entirely, costs
    # nothing at C(6,3) = 20, and guarantees p is never below 1/n_assignments.
    n_a, n_b = len(a), len(b)
    pooled = np.concatenate([a, b])
    n_tot = n_a + n_b
    assignments = list(combinations(range(n_tot), n_a))
    obs = float(a.sum())
    stat = np.array([pooled[list(idx)].sum() for idx in assignments])
    # small tolerance so floating-point equality counts as a tie, not a miss
    tol = 1e-12 * max(1.0, abs(obs))
    p_less = float(np.mean(stat <= obs + tol))        # arm_a ranks low
    p_greater = float(np.mean(stat >= obs - tol))     # arm_a ranks high
    two = float(min(1.0, 2.0 * min(p_less, p_greater)))
    rec["p_two_sided"] = two
    rec["n_assignments"] = len(assignments)
    if direction == D1MT_LOWER:
        one = p_less
        rec["ran_as_predicted"] = bool(a.mean() < b.mean())
    elif direction == D1MT_HIGHER:
        one = p_greater
        rec["ran_as_predicted"] = bool(a.mean() > b.mean())
    else:
        one = np.nan
    rec["p_one_sided"] = one
    rec["p_value"] = one if direction != NO_DIRECTION else two

    # the design floor, from the number of clusters actually contributing
    n_a, n_b = len(a), len(b)
    n_assign = len(list(combinations(range(n_a + n_b), n_a)))
    floor_one = 1.0 / n_assign
    floor_two = 2.0 / n_assign
    rec["floor_one_sided"] = floor_one
    rec["floor_two_sided"] = floor_two
    rec["floor_note"] = (f"floor {floor_one:.3f} one-sided, {floor_two:.3f} "
                         f"two-sided at {n_a} v {n_b} animals")
    return rec


def complete_separation(a, b):
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if not len(a) or not len(b):
        return False
    return bool(a.max() < b.min() or b.max() < a.min())


# =============================================================================
# centred radial position, identical construction to script 08
# =============================================================================

def add_centred_radial(df, struct_col="focus_id", radial_col="radial_pos",
                       pheno_col="pheno", region_col="region",
                       min_cells_for_balance=10):
    """
    Adds, for core cells only:
        radial_centered_core           reference = mean of ALL core cells
        radial_centered_core_balanced  reference = unweighted mean of the
                                       per-phenotype means, so composition
                                       cannot move the reference
        balance_reference              the balanced reference per structure
    Operates on one section at a time. Returns the frame with columns added.
    """
    d = df.copy()
    d[radial_col] = pd.to_numeric(d[radial_col], errors="coerce")
    d = d.loc[np.isfinite(d[radial_col])]
    d["is_core"] = d[region_col] == "core"
    core = d.loc[d["is_core"]]
    if not len(core):
        d["radial_centered_core"] = np.nan
        d["radial_centered_core_balanced"] = np.nan
        d["balance_reference"] = np.nan
        return d

    means = core.groupby(struct_col)[radial_col].transform("mean")
    d.loc[d["is_core"], "radial_centered_core"] = core[radial_col] - means

    bal = {}
    for k, g in core.groupby(struct_col):
        per_ph = g.groupby(pheno_col)[radial_col].agg(["mean", "size"])
        per_ph = per_ph.loc[per_ph["size"] >= min_cells_for_balance]
        bal[k] = (float(per_ph["mean"].mean()) if len(per_ph)
                  else float(g[radial_col].mean()))
    d.loc[d["is_core"], "radial_centered_core_balanced"] = (
        core[radial_col] - core[struct_col].map(bal))
    d["balance_reference"] = d[struct_col].map(bal)
    return d


def weighted_mean_by(df, group_cols, value_col, weight_col):
    """groupby.apply-free weighted mean, so no pandas deprecation warnings."""
    d = df.dropna(subset=[value_col, weight_col]).copy()
    d["_wv"] = d[value_col] * d[weight_col]
    g = d.groupby(group_cols, as_index=False)[["_wv", weight_col]].sum()
    g[value_col] = g["_wv"] / g[weight_col]
    return g.drop(columns=["_wv"])


# =============================================================================
# annotation
# =============================================================================

def format_test(res, model_row=None, unit_label="cells", extra=None,
                include_rank_sum=False):
    """include_rank_sum=None means show it only when there is no model line."""
    """
    Annotation text. The mixed model leads.

    include_rank_sum defaults to False. The exact rank-sum is bounded below by
    1/20 one-sided and 2/20 two-sided at three animals per arm, so every
    completely separated comparison returns exactly 0.05, which on a figure
    reads as six identical results rather than as one arithmetic bound. It is
    still computed and still written to the arm-tests table as the
    assumption-free check; it just does not compete for attention on the plot.
    """
    lines = []
    have_model = (model_row is not None
                  and np.isfinite(model_row.get("p_value", np.nan)))
    if include_rank_sum is None:
        # a panel whose model did not converge must still report a test
        include_rank_sum = not have_model
    if have_model:
        c = model_row.get("coef_D1MT_vs_ref", np.nan)
        lo = model_row.get("ci_low", np.nan)
        hi = model_row.get("ci_high", np.nan)
        p = model_row.get("p_value", np.nan)
        route = str(model_row.get("random_effects", ""))
        name = ("animal means t-test" if "fallback" in route else "mixed model")
        lines.append(f"{name}  {c:+.3f} [{lo:+.3f}, {hi:+.3f}]  p = {p:.3f}")
    if include_rank_sum and res is not None and np.isfinite(
            res.get("p_value", np.nan)):
        side = ("one-sided" if res.get("direction_declared") != NO_DIRECTION
                else "two-sided")
        tag = ("  (ran opposite to prediction)"
               if res.get("ran_as_predicted") is False else "")
        lines.append(f"clustered rank-sum  p = {res['p_value']:.3f} ({side})"
                     f"{tag}")
        if res.get("floor_note"):
            lines.append(res["floor_note"])
    n_obs = (model_row or {}).get("n_obs") or (res or {}).get("n_obs", 0)
    n_a = ((model_row or {}).get("n_clusters_a")
           or (res or {}).get("n_clusters_a", 0))
    n_b = ((model_row or {}).get("n_clusters_b")
           or (res or {}).get("n_clusters_b", 0))
    if n_obs:
        lines.append(f"n = {int(n_obs):,} {unit_label}, {n_a} v {n_b} animals")
    if extra:
        lines.append(extra)
    return "\n".join(lines)


def annotate_arm_test(ax, res, model_row=None, x0=0, x1=1, y=None,
                      unit_label="cells", extra=None, fontsize=18,
                      color="#000000", bracket=True, pad_frac=0.06,
                      va="bottom", include_rank_sum=False):
    """
    Draw a bracket between two x positions and write the standard three-line
    result underneath it. Returns the y used, so callers can extend the axis.
    """
    lo, hi = ax.get_ylim()
    span = hi - lo
    if y is None:
        y = hi - pad_frac * span
    if bracket:
        tick = 0.02 * span
        ax.plot([x0, x0, x1, x1], [y - tick, y, y, y - tick],
                color=color, linewidth=2.5, clip_on=False, zorder=6)
    ax.text(0.5 * (x0 + x1), y + 0.012 * span,
            format_test(res, model_row=model_row, unit_label=unit_label,
                        extra=extra, include_rank_sum=include_rank_sum),
            ha="center", va=va, fontsize=fontsize, color=color, zorder=6,
            linespacing=1.35)
    return y


def separation_band(ax, a, b, color="#000000", alpha=0.05, label=True,
                    fontsize=18):
    """Shade the gap between two completely separated groups."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if not len(a) or not len(b):
        return False
    if a.max() < b.min():
        lo, hi = a.max(), b.min()
    elif b.max() < a.min():
        lo, hi = b.max(), a.min()
    else:
        return False
    ax.axhspan(lo, hi, color=color, alpha=alpha, zorder=0)
    if label:
        xm = np.mean(ax.get_xlim())
        ax.text(xm, 0.5 * (lo + hi), "complete separation", ha="center",
                va="center", fontsize=fontsize)
    return True


# =============================================================================
# mixed linear model, the primary test from revision 3 onward
# =============================================================================

try:
    import statsmodels.formula.api as _smf
    HAVE_SM = True
except Exception:
    HAVE_SM = False



# ON. See the docstring of _cluster_mean_fallback for the calibration evidence.
USE_TWO_STAGE_FALLBACK = True


def _cluster_mean_fallback(rec, d, value_col, arm_col, cluster_col, arm_a,
                           arm_b, label, log10, direction=NO_DIRECTION):
    """
    Two-stage summary: average each animal, then a two-sample t-test on those
    six numbers.

    ENABLED, AND HERE IS THE EVIDENCE. This route was briefly disabled on the
    argument that it is assumption-laden and less sensitive than the clustered
    rank-sum. Less sensitive it is. Assumption-laden it is. But it is the only
    option in the following table that is correctly calibrated, and that was
    checked after the fact rather than before, which was the error.

    Simulated under a true null on the observed cluster structure, 3, 1 and 2
    foci against 13, 4 and 44, false positive rate at a nominal 5 percent:

        binomial GLMM, animal random effect            46.7 percent
        binomial GLMM, animal and focus random effects 46.7 percent
        binomial GEE, cluster-robust, normal reference 26.2 percent
        binomial GEE, cluster-robust, t on 4 df        13.4 percent
        two-stage t-test on the six animal means        4.4 percent

    Every method that treats cells or foci as independent observations is badly
    anticonservative at six clusters, including the mixed model families that
    are supposed to handle clustering. Adding a focus random effect changes
    nothing, because the problem is not unmodelled overdispersion. Treatment is
    assigned to six animals, so the arm contrast carries six independent pieces
    of information however the likelihood is written, and a model that computes
    its standard error as though it had tens of thousands will be wrong by the
    corresponding factor.

    This route is therefore used wherever the linear mixed model returns no
    usable standard error, which happens when the animal-level variance is
    essentially zero relative to within-animal variance. Log-transforming does
    not rescue those fits; the random effects covariance is singular because of
    the data rather than the parameterisation.

    Results from this route are labelled "cluster means t-test (fallback)" and
    rendered on figures as "animal means t-test", never as a mixed model.
    """
    if not USE_TWO_STAGE_FALLBACK:
        return rec
    try:
        from scipy import stats as _sst
    except Exception:
        return rec
    g = d.groupby([cluster_col, arm_col])["_y"].mean().reset_index()
    ga = g.loc[g[arm_col] == arm_a, "_y"].to_numpy()
    gb = g.loc[g[arm_col] == arm_b, "_y"].to_numpy()
    if len(ga) < 2 or len(gb) < 2:
        return rec
    t = _sst.ttest_ind(ga, gb, equal_var=True)
    c2 = float(ga.mean() - gb.mean())
    # scipy returns a TWO-SIDED p. The clustered rank-sum is one-sided on a
    # pre-specified direction, so leaving this two-sided applies two different
    # conventions to the same hypothesis. One-sided p is half the two-sided p
    # when the effect runs as predicted, and one minus that otherwise.
    p_two = float(t.pvalue)
    if direction == D1MT_LOWER:
        p_one = p_two / 2.0 if c2 < 0 else 1.0 - p_two / 2.0
        ran_ok = bool(c2 < 0)
    elif direction == D1MT_HIGHER:
        p_one = p_two / 2.0 if c2 > 0 else 1.0 - p_two / 2.0
        ran_ok = bool(c2 > 0)
    else:
        p_one, ran_ok = np.nan, None
    sp = np.sqrt(((len(ga) - 1) * ga.var(ddof=1) + (len(gb) - 1) * gb.var(ddof=1))
                 / (len(ga) + len(gb) - 2))
    se2 = float(sp * np.sqrt(1 / len(ga) + 1 / len(gb)))
    dfree = len(ga) + len(gb) - 2
    crit = float(_sst.t.ppf(0.975, dfree))
    ok = bool(np.isfinite(se2) and se2 > 0 and np.isfinite(t.pvalue))
    rec.update({"coef_D1MT_vs_ref": c2, "std_err": se2,
                "ci_low": c2 - crit * se2, "ci_high": c2 + crit * se2,
                "p_value": (p_one if direction != NO_DIRECTION and np.isfinite(p_one)
                            else p_two),
                "p_two_sided": p_two, "p_one_sided": p_one,
                "df": dfree, "ran_as_predicted": ran_ok,
                "n_clusters_a": int(len(ga)), "n_clusters_b": int(len(gb)),
                "random_effects": "cluster means t-test (fallback)",
                "converged": ok})
    if ok:
        print(f"    NOTE [fit_arm_lmm {label}]: mixed model gave no usable "
              f"standard error, so this quantity uses the two-stage fallback, "
              f"a t-test on the {len(ga)} v {len(gb)} animal means on "
              f"{dfree} df. Effective n is the animals.")
    return rec


def fit_arm_lmm(df, value_col, arm_col, cluster_col, arm_a, arm_b,
                struct_col=None, log10=False, label="", min_per_cluster=1,
                direction=NO_DIRECTION, covariates=None):
    """
    outcome ~ arm, random intercept for cluster (animal), structure nested
    within cluster when struct_col is given.

    This is the CELL-LEVEL (or focus-level) test. Every observation enters the
    likelihood; the random effects absorb the fact that observations from one
    animal share a baseline. Unlike the exact rank-sum it is not floored at
    1/20, because the p-value comes from a likelihood rather than from counting
    the 20 possible arm assignments.

    THE ASSUMPTION IT BUYS THAT WITH. Residuals and animal effects are taken to
    be roughly normal, and statsmodels reports asymptotic z-tests with no
    small-sample degrees-of-freedom correction. With six clusters that is
    somewhat liberal, so the p-value is the optimistic end of the range and the
    exact rank-sum is the conservative end. Report both.

    log10=True fits on log10(value) and adds fold_change = 10**coef, which is
    the right scale for densities spanning an order of magnitude.
    """
    rec = {"analysis": label, "outcome": value_col, "log10": bool(log10),
           "coef_D1MT_vs_ref": np.nan, "std_err": np.nan, "ci_low": np.nan,
           "ci_high": np.nan, "p_value": np.nan, "n_obs": 0, "n_clusters": 0,
           "n_clusters_a": 0, "n_clusters_b": 0,
           "random_effects": "", "fold_change": np.nan, "converged": False,
           "clusters_dropped": ""}
    if not HAVE_SM:
        return rec
    covariates = list(covariates or [])
    need = ([value_col, arm_col, cluster_col]
            + ([struct_col] if struct_col else []) + covariates)
    if any(c not in df.columns for c in need):
        return rec
    d = df.dropna(subset=need).copy()
    d["_y"] = pd.to_numeric(d[value_col], errors="coerce")
    if log10:
        d = d.loc[d["_y"] > 0]
        d["_y"] = np.log10(d["_y"])
    d = d.dropna(subset=["_y"])
    counts = d.groupby(cluster_col).size()
    dropped = sorted(counts.loc[counts < min_per_cluster].index.tolist())
    d = d.loc[d[cluster_col].isin(counts.loc[counts >= min_per_cluster].index)]
    rec["clusters_dropped"] = ",".join(str(x) for x in dropped)
    if dropped:
        # loud, because silently dropping an animal changes the design
        print(f"    WARNING [fit_arm_lmm {label}]: dropped cluster(s) "
              f"{dropped} for having fewer than {min_per_cluster} observations")
    if not len(d) or d[arm_col].nunique() < 2 or d[cluster_col].nunique() < 3:
        return rec
    d["_arm"] = (d[arm_col] == arm_a).astype(float)   # 1 = arm_a (D1MT)
    rec["n_obs"] = int(len(d))
    rec["n_clusters"] = int(d[cluster_col].nunique())
    rec["n_clusters_a"] = int(d.loc[d[arm_col] == arm_a, cluster_col].nunique())
    rec["n_clusters_b"] = int(d.loc[d[arm_col] == arm_b, cluster_col].nunique())

    # covariates enter as fixed effects. A pooled-lineage model adjusts for cell
    # type so that composition within the lineage cannot drive the arm effect.
    terms = ["_arm"]
    for i, c in enumerate(covariates):
        if d[c].dtype == object or str(d[c].dtype).startswith("category"):
            d[f"_cv{i}"] = d[c].astype(str)
            terms.append(f"C(_cv{i})")
        else:
            v = pd.to_numeric(d[c], errors="coerce")
            sd = v.std()
            d[f"_cv{i}"] = (v - v.mean()) / sd if sd and sd > 0 else 0.0
            terms.append(f"_cv{i}")
    formula = "_y ~ " + " + ".join(terms)
    rec["covariates"] = ",".join(covariates)

    import warnings
    res, mode = None, ""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if struct_col is not None:
            try:
                d["_sk"] = (d[cluster_col].astype(str) + "_"
                            + d[struct_col].astype(int).astype(str))
                md = _smf.mixedlm(formula, data=d, groups=d[cluster_col],
                                  re_formula="1",
                                  vc_formula={"struct": "0 + C(_sk)"})
                res = md.fit(reml=True, method="lbfgs", maxiter=300)
                mode = "animal + structure"
            except Exception:
                res = None
        if res is None or not np.isfinite(res.params.get("_arm", np.nan)):
            try:
                md = _smf.mixedlm(formula, data=d, groups=d[cluster_col],
                                  re_formula="1")
                res = md.fit(reml=True, method="lbfgs", maxiter=300)
                mode = "animal only"
            except Exception:
                res = None
    if res is None:
        return _finish_fold(
            _cluster_mean_fallback(rec, d, value_col, arm_col, cluster_col,
                                   arm_a, arm_b, label, log10,
                                   direction=direction), log10)
    c = float(res.params.get("_arm", np.nan))
    se = float(res.bse.get("_arm", np.nan))
    pv = float(res.pvalues.get("_arm", np.nan))
    rec.update({"coef_D1MT_vs_ref": c, "std_err": se,
                "ci_low": c - 1.96 * se, "ci_high": c + 1.96 * se,
                "p_value": pv, "random_effects": mode})
    # A returned result object is NOT convergence. Arm is constant within
    # animal, so with a random intercept alone the animal effect and the arm
    # effect compete and the Hessian can go singular, which shows up as a NaN
    # standard error while statsmodels reports no error at all. Silently
    # dropping the model line from a figure is the worst possible outcome, so
    # this is checked and said out loud.
    rec["converged"] = bool(np.isfinite(c) and np.isfinite(se) and se > 0
                            and np.isfinite(pv))
    if not rec["converged"]:
        rec = _cluster_mean_fallback(rec, d, value_col, arm_col, cluster_col,
                                     arm_a, arm_b, label, log10,
                                     direction=direction)
    if not rec["converged"]:
        print(f"    WARNING [fit_arm_lmm {label}]: no usable fit by any route "
              f"(coef={c:.4f}, se={se}). Report the clustered rank-sum for "
              f"this quantity instead.")
    if log10 and rec["converged"]:
        rec["fold_change"] = float(10 ** rec["coef_D1MT_vs_ref"])
    return rec


def _finish_fold(rec, log10):
    if log10 and rec.get("converged") and np.isfinite(rec["coef_D1MT_vs_ref"]):
        rec["fold_change"] = float(10 ** rec["coef_D1MT_vs_ref"])
    return rec


# =============================================================================
# violin with visible per-observation and per-animal detail
# =============================================================================

def violin_with_points(ax, frame, value_col, arm_col, cluster_col,
                       arm_order, arm_colors, cluster_colors, cluster_markers,
                       centre=0.0, spacing=1.0, violin_width=0.62,
                       max_points=1200, point_size=32, point_alpha=0.45,
                       point_offset=0.0, point_jitter=0.19,
                       median_offset=0.0, median_spread=0.11, median_size=620,
                       show_box=True, rng=None, log=False):
    """
    Per arm, at x = centre + i*spacing, drawn back to front:
        the violin, filled in the arm colour                       zorder 1
        a narrow unfilled box for the quartiles                    zorder 4
        the raw observations, jittered ON THE CENTRE LINE,
        coloured by animal                                         zorder 3
        one large black-edged marker per animal at that animal's
        MEDIAN, spread slightly so three animals do not stack      zorder 10

    Points sit on the centre because that is how a violin of this kind is
    normally read. The box is unfilled so the points show through it, and the
    animal medians are the top layer so they are never obscured.

    Returns the x positions used.
    """
    rng = rng or np.random.default_rng(0)
    positions = [centre + i * spacing for i in range(len(arm_order))]

    data, pos, cols = [], [], []
    for x, c in zip(positions, arm_order):
        v = frame.loc[frame[arm_col] == c, value_col].dropna().to_numpy()
        if len(v) >= 5:
            data.append(v); pos.append(x); cols.append(arm_colors[c])
    if data:
        parts = ax.violinplot(data, positions=pos, widths=violin_width,
                              showextrema=False, showmedians=False)
        for body, c in zip(parts["bodies"], cols):
            body.set_facecolor(c); body.set_alpha(0.32)
            body.set_edgecolor("#333333"); body.set_linewidth(1.8)
            body.set_zorder(1)

    for x, c in zip(positions, arm_order):
        sub = frame.loc[frame[arm_col] == c]
        if not len(sub):
            continue
        take = sub if len(sub) <= max_points else sub.sample(
            max_points, random_state=0)
        jx = x + point_offset + rng.uniform(-point_jitter, point_jitter,
                                            size=len(take))
        ax.scatter(jx, take[value_col],
                   s=point_size, alpha=point_alpha, linewidths=0,
                   c=[cluster_colors.get(k, "#888888") for k in take[cluster_col]],
                   zorder=3, rasterized=True)

    if data and show_box:
        ax.boxplot(data, positions=pos, widths=0.20, patch_artist=True,
                   showfliers=False, zorder=4,
                   medianprops=dict(color="#000000", linewidth=4),
                   boxprops=dict(facecolor="none", edgecolor="#000000",
                                 linewidth=2.5),
                   whiskerprops=dict(color="#000000", linewidth=2.5),
                   capprops=dict(color="#000000", linewidth=2.5))

    for x, c in zip(positions, arm_order):
        sub = frame.loc[frame[arm_col] == c]
        if not len(sub):
            continue
        meds = sub.groupby(cluster_col)[value_col].median()
        for j, (k, v) in enumerate(meds.items()):
            xo = x + median_offset + (j - (len(meds) - 1) / 2.0) * median_spread
            ax.scatter([xo], [v], s=median_size,
                       color=cluster_colors.get(k, "#888888"),
                       marker=cluster_markers.get(k, "o"),
                       edgecolor="#000000", linewidth=4.5, zorder=10)
    return positions


def cluster_median_legend(cluster_colors, cluster_markers, labels):
    from matplotlib.lines import Line2D
    return [Line2D([0], [0], color=cluster_colors[k], marker=cluster_markers[k],
                   markersize=16, markeredgecolor="#000000",
                   markeredgewidth=2.5, linestyle="none", label=labels[k])
            for k in labels]
