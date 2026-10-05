#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
akoya_arm_stats.py

Shared statistics and annotation helpers for the AKOYA analysis series.

WHY THIS IS A MODULE AND NOT COPIED INTO EACH SCRIPT
    The IDO1 cutoff plateau search was implemented twice, in scripts 03 and 09,
    and the two reported different windows for a year without anyone noticing.
    Anything used by more than one script lives here so that cannot happen
    again. Scripts 09 and 10 import it.

    CORRECTED 2 October 2026. The previous header said "Scripts 04, 08 and 09
    import it", which was wrong in both directions. 04 and 08 do not import it
    at all, each defining its own local fit_mixed, and 10 does. The four
    scripts with a local fit_mixed are 06, 07, 07b and 08, and after the
    revision 7 change below they import small_sample_inference from here.

REVISION 7, 2 October 2026
    Six changes. The first is a bug fix, the second is the degrees-of-freedom
    correction, the rest are transparency.

    1. THE REVISION 3 FALLBACK FIX, BACKPORTED. fit_arm_lmm triggered its
       fallback on a finite coefficient alone, so a nested fit returning a NaN
       standard error never reached the animal-only refit that might have
       succeeded. Script 08's fit_mixed fixed this in its own revision 3 and
       the fix never came back here, which is precisely the failure mode this
       module exists to prevent. Both the fallback trigger and the acceptance
       test now go through _usable(), matching script 08.

    2. DENOMINATOR DEGREES OF FREEDOM. See small_sample_inference below. The
       arm contrast is between-animal and is now referred to t(n_animals - 2)
       rather than to the standard normal, and the interval multiplier is the
       matching t quantile rather than 1.96. The Wald z p-value is retained as
       p_wald_z so every historical table can be reconciled against the new
       one rather than taken on trust.

    3. A SILENT FALL TO ANIMAL-ONLY IS NOW LOUD. The structure-nesting attempt
       swallowed every exception, so a failure to build the structure key, for
       instance when focus_id is not castable to int, dropped the fit to
       animal-only random effects with no message. That matters more than it
       looks: simulated at this design the animal-only route converges in 8 to
       35 percent of null replicates against 99.6 percent for the nested fit,
       so an unannounced demotion to animal-only is an unannounced move to a
       route that mostly does not work.

    4. add_centred_radial now records which structures fell back to the
       cell-weighted reference. When no phenotype in a structure clears
       min_cells_for_balance the balanced reference silently became the plain
       mean of all core cells, so a structure could sit inside an analysis
       labelled composition-balanced while carrying a cell-weighted centre.
       The behaviour is unchanged. A balance_fallback column and a warning now
       say where it happened.

    5. format_test had two consecutive string literals, so the long
       explanation was a no-op expression and never reached help(). Merged.
       The model line now also prints the t degrees of freedom, so a figure
       states its own reference.

    6. The one-sided section of this docstring described a convention the
       pipeline does not implement. See the note under section 2.

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

       THE FLOOR IS CALIBRATED, AND THE IMBALANCE DOES NOT BREAK IT. Checked at
       20,000 replicates under a true null: the exact permutation on the six
       animal means rejects at 9.6 to 10.4 percent against its nominal 0.10
       across every animal ICC tried and every cells-per-animal imbalance
       tried, including a deliberate treated 50/50/2000 against untreated
       2000/2000/50. The worry that wildly unequal cell counts make the animal
       means non-exchangeable and inflate the test was tested and is wrong.

    2. ONE-SIDED TESTING WITH PRE-SPECIFIED DIRECTION
       **UNUSED. Retained so imports do not break, and documented so the code
       and the Methods do not disagree.** D1MT_LOWER and D1MT_HIGHER are
       supported by datta_satten_ranksum and by the two-stage fallback, and no
       production script passes direction= to either. Everything reported in
       this project is TWO-SIDED with the 0.100 floor, which is the right
       convention here because effects can run in opposite directions across
       cell types, and it is the convention the context notes describe.

       Worth knowing what it would buy, because the question recurs. A
       genuinely pre-specified one-sided test has a floor of 1/20 = 0.05 rather
       than 2/20 = 0.100, so S1 would reach 0.05. Choosing that after the
       results are in is not available. It is a design decision for the next
       study, not a reporting decision for this one.

    3. SMALL-SAMPLE INFERENCE on an already-fitted fixed effect. New in
       revision 7, and the one thing in this module a statistical reviewer
       will look for.

    4. CENTRED RADIAL POSITION, both centrings, identical to script 08.

    5. ANNOTATION HELPER that renders the same three-line bracket on any arm
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
# UNUSED by every production script. See section 2 of the module docstring.
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
               D1MT side for the purpose of interpreting the direction. NOTE
               that no production script passes this; the project convention is
               two-sided. See section 2 of the module docstring.

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
# small-sample inference on an already-fitted fixed effect
# =============================================================================
#
# Which small-sample reference to use. "between_animal" is the declared
# default: denominator df = n_animals - 2, the between-animal contrast df.
# "satterthwaite" is reserved and not implemented; it requires lmerTest through
# pymer4. "none" restores the pre-revision-7 Wald z behaviour and exists only so
# a historical table can be reproduced exactly.
DF_METHOD = "between_animal"


def small_sample_inference(coef, se, n_animals, n_arms=2, alpha=0.05,
                           df_method=None):
    """
    Re-refer an already-fitted fixed effect to a small-sample t distribution.

    Does NOT refit anything. Takes the coefficient and standard error the mixed
    model produced and returns the p-value and interval that the number of
    ANIMALS supports rather than the number of cells.

    WHY THIS EXISTS
        statsmodels mixedlm reports a Wald z on each fixed effect, coefficient
        over standard error referred to the standard normal, and applies no
        denominator degrees of freedom. res.df_resid is cells minus rank, so a
        six-animal three-thousand-cell fit reports df_resid near 2998 and the
        arm effect is tested as though it carried that much information.

        The arm contrast is between-animal. Treatment was applied to six
        animals, so it carries six independent pieces of information however
        the likelihood is written, and its proper denominator is about
        n_animals - 2 = 4. Referring a t(4) statistic to a normal is
        anti-conservative, and so is building a 95 percent interval as
        coef +/- 1.96 * se when the multiplier should be t(0.975, 4) = 2.776.
        The old interval was too narrow by 42 percent.

    CALIBRATION, simulated at this design
        Six animals, 3 v 3, structures 1/1/2 against 13/9/54, cells per animal
        122/83/869 against 1491/83/778, outcome built as animal effect plus
        structure effect plus cell noise with NO arm effect, 250 replicates per
        cell, nominal alpha 0.05. The nested fit converged in 99.6 percent of
        replicates, so these are not a selected subset.

            ICC_animal     Wald z (old)      same statistic vs t(4)
            0.00                 4.4 %               0.4 %
            ~0.09                9.2 %               4.4 %
            ~0.38                8.0 %               4.8 %

        The Wald z route runs at 8 to 9 percent where it matters, reproducing
        Kahan et al., Trials 2016, 17:438, which reports 8.4 and 8.6 percent at
        six clusters. The t(4) route is calibrated or conservative throughout.
        Satterthwaite lands in the same place, because with treatment constant
        within animal the arm contrast has close to n_animals - 2 df whatever
        the cell count, so this closes the gap without an R dependency.

    NOTE ON COVARIATES
        A CELL-level covariate (local density, phenotype) does not consume
        between-animal degrees of freedom, so n_arms stays 2. An ANIMAL-level
        covariate does. If a model ever adds one, raise n_arms to match or the
        df will be too generous.

    NOTE ON THE TWO-STAGE FALLBACK
        _cluster_mean_fallback was always correct on this point. It builds its
        interval with t.ppf(0.975, n_animals - 2) and takes its p from a
        t-test on the animal means. Only the LMM route was wrong, which is why
        the fallback and the exact test agreed with each other while the model
        p was the outlier.

    Returns a dict with p_value, ci_low, ci_high, df, t_stat, t_crit and
    df_method. Every field is NaN when the inputs are unusable. This function
    does not invent a result from a NaN standard error, so the caller keeps its
    own convergence guard.
    """
    out = {"p_value": np.nan, "ci_low": np.nan, "ci_high": np.nan,
           "df": np.nan, "t_stat": np.nan, "t_crit": np.nan,
           "df_method": df_method or DF_METHOD}
    method = out["df_method"]

    if not np.isfinite(coef) or not np.isfinite(se) or se <= 0:
        return out

    if method == "none":
        if not HAVE_SCIPY:
            return out
        z = float(coef / se)
        out.update({"p_value": float(2.0 * (1.0 - _st.norm.cdf(abs(z)))),
                    "ci_low": float(coef - 1.96 * se),
                    "ci_high": float(coef + 1.96 * se),
                    "df": np.inf, "t_stat": z, "t_crit": 1.96})
        return out

    if method == "satterthwaite":
        raise NotImplementedError(
            "Satterthwaite df require lmerTest via pymer4. Use "
            "DF_METHOD='between_animal' until that route exists.")

    try:
        n_animals = int(n_animals)
    except Exception:
        return out
    dfree = n_animals - int(n_arms)
    if dfree < 1 or not HAVE_SCIPY:
        # a fit on no more animals than parameters has no inference to report
        out["df"] = float(dfree)
        return out

    t_stat = float(coef / se)
    t_crit = float(_st.t.ppf(1.0 - alpha / 2.0, dfree))
    out.update({"p_value": float(2.0 * (1.0 - _st.t.cdf(abs(t_stat), dfree))),
                "ci_low": float(coef - t_crit * se),
                "ci_high": float(coef + t_crit * se),
                "df": float(dfree), "t_stat": t_stat, "t_crit": t_crit})
    return out


def df_note(n_animals, n_arms=2):
    """One line for a figure annotation or a table footer."""
    d = int(n_animals) - int(n_arms)
    return (f"p and CI on t({d}) df, the between-animal contrast df at "
            f"{int(n_animals)} animals")


# =============================================================================
# centred radial position, identical construction to script 08
# =============================================================================

def add_centred_radial(df, struct_col="focus_id", radial_col="radial_pos",
                       pheno_col="pheno", region_col="region",
                       min_cells_for_balance=10, label=""):
    """
    Adds, for core cells only:
        radial_centered_core           reference = mean of ALL core cells
        radial_centered_core_balanced  reference = unweighted mean of the
                                       per-phenotype means, so composition
                                       cannot move the reference
        balance_reference              the balanced reference per structure
        balance_fallback               True where the balanced reference COULD
                                       NOT be built and the cell-weighted mean
                                       was used instead

    Operates on one section at a time. Returns the frame with columns added.

    REVISION 7 ON THE FALLBACK. When no phenotype in a structure has at least
    min_cells_for_balance cells, the balanced reference falls back to the plain
    mean of all core cells in that structure. That is the cell-weighted
    reference, so the affected structure sits inside an analysis labelled
    composition-balanced while carrying exactly the centring the balanced
    reference was built to avoid. The behaviour is UNCHANGED, because changing
    it would move numbers. What is new is that it is visible: the
    balance_fallback column marks the rows and a warning names the structures.
    Thin treated structures are the ones at risk, so check this column before
    quoting a balanced-reference coefficient.
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
        d["balance_fallback"] = False
        return d

    means = core.groupby(struct_col)[radial_col].transform("mean")
    d.loc[d["is_core"], "radial_centered_core"] = core[radial_col] - means

    bal, fell_back = {}, []
    for k, g in core.groupby(struct_col):
        per_ph = g.groupby(pheno_col)[radial_col].agg(["mean", "size"])
        per_ph = per_ph.loc[per_ph["size"] >= min_cells_for_balance]
        if len(per_ph):
            bal[k] = float(per_ph["mean"].mean())
        else:
            bal[k] = float(g[radial_col].mean())
            fell_back.append((k, int(len(g))))
    d.loc[d["is_core"], "radial_centered_core_balanced"] = (
        core[radial_col] - core[struct_col].map(bal))
    d["balance_reference"] = d[struct_col].map(bal)
    d["balance_fallback"] = d[struct_col].isin([k for k, _ in fell_back])

    if fell_back:
        tag = f" [{label}]" if label else ""
        detail = ", ".join(f"{k} ({n} core cells)" for k, n in fell_back)
        print(f"    WARNING [add_centred_radial{tag}]: no phenotype reached "
              f"{min_cells_for_balance} core cells in {len(fell_back)} "
              f"structure(s), so the BALANCED reference fell back to the "
              f"cell-weighted mean there: {detail}")
        print("      Those rows are marked balance_fallback=True. A balanced "
              "coefficient resting on them is not composition-corrected.")
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
    """
    Annotation text. The mixed model leads.

    include_rank_sum defaults to False. The exact rank-sum is bounded below by
    1/20 one-sided and 2/20 two-sided at three animals per arm, so every
    completely separated comparison returns exactly the same number, which on a
    figure reads as six identical results rather than as one arithmetic bound.
    It is still computed and still written to the arm-tests table as the
    assumption-free check; it just does not compete for attention on the plot.

    include_rank_sum=None means show it only when there is no model line, so a
    panel whose model did not converge still reports a test.

    REVISION 7. The model line now carries its t degrees of freedom, so a
    figure states the reference its p-value was taken against rather than
    leaving a reader to assume a normal. Rows fitted before revision 7 have no
    df field and render as before.
    """
    lines = []
    have_model = (model_row is not None
                  and np.isfinite(model_row.get("p_value", np.nan)))
    if include_rank_sum is None:
        include_rank_sum = not have_model
    if have_model:
        c = model_row.get("coef_D1MT_vs_ref", np.nan)
        lo = model_row.get("ci_low", np.nan)
        hi = model_row.get("ci_high", np.nan)
        p = model_row.get("p_value", np.nan)
        route = str(model_row.get("random_effects", ""))
        name = ("animal means t-test" if "fallback" in route else "mixed model")
        dfv = model_row.get("df", np.nan)
        # np.isfinite is False for inf, which is what DF_METHOD="none" records,
        # so the normal-reference case correctly renders no t tag.
        dftag = f" t({int(dfv)})" if np.isfinite(dfv) else ""
        lines.append(f"{name}{dftag}  {c:+.3f} [{lo:+.3f}, {hi:+.3f}]  "
                     f"p = {p:.3f}")
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


def _usable(res, term="_arm"):
    """
    A returned result object is NOT a usable fit.

    Arm is constant within animal, so the animal effect and the arm effect
    compete and the Hessian can go singular, which surfaces as a NaN standard
    error while statsmodels reports no error at all. A finite coefficient alone
    is not enough.

    REVISION 7. This is script 08's revision 3 test, backported. Before this,
    fit_arm_lmm triggered its fallback on np.isfinite(res.params) alone, so a
    nested fit with a NaN standard error never reached the animal-only refit
    that might have succeeded, and nothing caught it on the way out either.
    """
    if res is None:
        return False
    try:
        c = float(res.params.get(term, np.nan))
        e = float(res.bse.get(term, np.nan))
        p = float(res.pvalues.get(term, np.nan))
    except Exception:
        return False
    return bool(np.isfinite(c) and np.isfinite(e) and e > 0 and np.isfinite(p))


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

    NOTE, revision 7. This route was already correct on degrees of freedom. It
    builds its interval with t.ppf(0.975, n_animals - 2) and takes its p from a
    t-test on the animal means, which is why it and the exact test agreed with
    each other all along while the LMM p-value was the outlier.
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
                "df": dfree, "t_crit": crit, "df_method": "animal_means_t",
                "ran_as_predicted": ran_ok,
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

    WHAT IT ASSUMES, AND WHAT REVISION 7 CHANGED. Residuals and animal effects
    are taken to be roughly normal. statsmodels then reports an asymptotic Wald
    z with NO small-sample degrees-of-freedom correction, which at six clusters
    runs at 8 to 9 percent false positive against a nominal 5. Revision 7 refers
    the same statistic to t(n_animals - 2) through small_sample_inference and
    builds the interval with the matching t quantile. The Wald z p is retained
    as p_wald_z so historical tables can be reconciled.

    The exact animal-level test remains the conservative end of the range and
    should still be reported beside this. After the df correction the two tend
    to agree, which is the reassuring outcome.

    log10=True fits on log10(value) and adds fold_change = 10**coef, which is
    the right scale for densities spanning an order of magnitude.
    """
    rec = {"analysis": label, "outcome": value_col, "log10": bool(log10),
           "coef_D1MT_vs_ref": np.nan, "std_err": np.nan, "ci_low": np.nan,
           "ci_high": np.nan, "p_value": np.nan, "p_wald_z": np.nan,
           "df": np.nan, "t_crit": np.nan, "df_method": DF_METHOD,
           "n_obs": 0, "n_clusters": 0,
           "n_clusters_a": 0, "n_clusters_b": 0,
           "random_effects": "", "fold_change": np.nan, "converged": False,
           "clusters_dropped": "", "nesting_note": ""}
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
    # NOTE these are CELL-level covariates and do not consume between-animal
    # degrees of freedom, so n_arms stays at its default of 2 below.
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
    nest_err = ""
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
            except Exception as e:
                # REVISION 7: this was a bare `res = None`, so a failure here,
                # for instance struct_col not castable to int, demoted the fit
                # to animal-only with no message at all. The animal-only route
                # converges in 8 to 35 percent of null replicates against 99.6
                # percent for the nested fit, so an unannounced demotion is an
                # unannounced move to a route that mostly does not work.
                nest_err = f"{type(e).__name__}: {e}"
                res = None
        # REVISION 7: the trigger is _usable(), not a finite coefficient. A
        # nested fit with a NaN standard error now reaches the animal-only
        # refit, which is script 08's revision 3 behaviour.
        if not _usable(res):
            if struct_col is not None:
                print(f"    WARNING [fit_arm_lmm {label}]: structure-nested "
                      f"fit unusable, refitting with animal random intercept "
                      f"only.")
                if nest_err:
                    print(f"      nesting failed with {nest_err}")
                print("      An 'animal only' fit is a SELECTED fit at this "
                      "design. Treat the row as a diagnostic, not a result.")
            try:
                md = _smf.mixedlm(formula, data=d, groups=d[cluster_col],
                                  re_formula="1")
                res = md.fit(reml=True, method="lbfgs", maxiter=300)
                mode = "animal only"
            except Exception:
                res = None
    rec["nesting_note"] = nest_err

    if not _usable(res):
        rec = _cluster_mean_fallback(rec, d, value_col, arm_col, cluster_col,
                                     arm_a, arm_b, label, log10,
                                     direction=direction)
        if not rec["converged"]:
            print(f"    WARNING [fit_arm_lmm {label}]: no usable fit by any "
                  f"route. Report the clustered rank-sum for this quantity "
                  f"instead.")
        return _finish_fold(rec, log10)

    c = float(res.params.get("_arm", np.nan))
    se = float(res.bse.get("_arm", np.nan))

    # REVISION 7: the p-value and the interval now come from the number of
    # ANIMALS, not the number of cells. p_wald_z keeps what statsmodels said so
    # the change is auditable rather than taken on trust.
    ss = small_sample_inference(c, se, rec["n_clusters"])
    rec.update({"coef_D1MT_vs_ref": c, "std_err": se,
                "ci_low": ss["ci_low"], "ci_high": ss["ci_high"],
                "p_value": ss["p_value"],
                "p_wald_z": float(res.pvalues.get("_arm", np.nan)),
                "df": ss["df"], "t_crit": ss["t_crit"],
                "df_method": ss["df_method"],
                "random_effects": mode})
    rec["converged"] = bool(np.isfinite(c) and np.isfinite(se) and se > 0
                            and np.isfinite(ss["p_value"]))
    if not rec["converged"]:
        rec = _cluster_mean_fallback(rec, d, value_col, arm_col, cluster_col,
                                     arm_a, arm_b, label, log10,
                                     direction=direction)
    if not rec["converged"]:
        print(f"    WARNING [fit_arm_lmm {label}]: no usable fit by any route "
              f"(coef={c:.4f}, se={se}). Report the clustered rank-sum for "
              f"this quantity instead.")
    return _finish_fold(rec, log10)


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
