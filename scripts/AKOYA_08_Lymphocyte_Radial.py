#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - LYMPHOCYTE RADIAL POSITION, THE PRIMARY TEST
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 08 of the AKOYA analysis series.

THE HYPOTHESIS
    Lymphocytes sit closer to the core of a myeloid focus in D1MT-treated
    animals than in untreated animals. Because the outcome is position rather
    than intensity, it is immune to the slide confound that limits most of this
    panel: coordinates carry no batch signal. Cell IDENTITY still does, which is
    why the lineage split below matters.

    Lymphocytes are also not used to define foci (detection is myeloid-only), so
    unlike every macrophage and neutrophil result these positions are not
    circular.

WHAT THIS FIXES FROM SCRIPTS 06 AND 07
    1. CENTERED RADIAL POSITION IS THE OUTCOME. Scripts 06 and 07 modelled raw
       radial_pos while reporting delta as the effect size. Those are different
       quantities, and CD4- T cells showed it: per-animal deltas of -0.21, -0.12,
       -0.09 against -0.05, +0.01, +0.17, but a model coefficient of +0.010 with
       p = 0.94.

       delta = mean radial of population P, minus the mean of a random subset of
       the same size. A random subset's expected mean is the structure's overall
       mean, so delta is exactly (radial_pos - structure mean radial). Centering
       each cell on its structure mean therefore makes the model test the same
       quantity the effect sizes report. That is the outcome used here.

    2. CORE CELLS ONLY IN THE PRIMARY MODEL. Core cells are normalised 0 to 1 by
       each focus's own size, but cuff cells run 1 to 2 across a FIXED 150 um
       band that is not size-normalised. Mixing them reintroduces structure size
       through the back door. Core-only is primary; core-plus-cuff is secondary.

    3. NORMALISED DISTANCE FAMILY DROPPED. Dividing distance by focus radius
       correlates with radius at -0.12 to -0.68, a stronger confound than the raw
       values had (-0.10 to +0.39). It over-corrects. Delta distance
       (-0.21 to +0.24) is retained.

    4. ONE POOLED TEST INSTEAD OF TEN. Five lymphocyte populations are pooled
       into a single model with cell type as a covariate, so composition
       differences between arms cannot drive the result, and power is pooled
       rather than split across ten underpowered comparisons.

    5. SPILLOVER CONTROL FOR CO-EXPRESSION. The iNOS / Arginase-1 co-expression
       ratio separates the arms completely at the 75th percentile and above. But
       segmentation spillover inflates apparent co-expression of ANY two markers,
       and gets worse as cells pack more densely; untreated foci are three to
       five fold denser. Control pairs that should never co-occur in one
       macrophage (CD3e with Pan-Cytokeratin, CD3e with CD20, CD45 with
       Pan-Cytokeratin) are computed the same way. If those separate by arm too,
       the iNOS / Arginase-1 result is spillover and must be dropped.

WHAT IS ALREADY ESTABLISHED (fold-6 run)
    Plasma cells (q = 0.088) and B cells (q = 0.088) survive BH correction in the
    radial family and are NOT circular. All three treated animals are negative on
    five of six populations. Delta distance confirms the same axis, with IDO1- to
    plasma cells at p = 0.046 and the only two non-normalised complete
    separations both on plasma cells.

INPUTS
    structures_rev4/cell_assignments/<section>_cell_structures.csv
    structures_rev4/tables/35_foci_structures_relative.csv
    distance_stats/tables/53_nn_per_structure.csv   (for the distance null)
    AKOYA/data/<section>.csv                        (spillover control markers)

OUTPUTS
    figures/  F51 .. F56
    tables/   70 .. 78

USAGE
    conda activate sc_pre
    python AKOYA_08_Lymphocyte_Radial.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================

DATA_DIR = "/master/jlehle/WORKING/AKOYA/data"
IN_DIR = "/master/jlehle/WORKING/AKOYA/structures_rev4"
CELL_DIR = f"{IN_DIR}/cell_assignments"
FOCI_TABLE = f"{IN_DIR}/tables/35_foci_structures_relative.csv"
NN_NULL_TABLE = "/master/jlehle/WORKING/AKOYA/distance_stats/tables/53_nn_per_structure.csv"
OUT_DIR = "/master/jlehle/WORKING/AKOYA/lymphocyte_radial"

CONDITION_ORDER = ["D1MT", "Untreated"]
CONDITION_COLORS = {"D1MT": "#2C7FB8", "Untreated": "#D95F02"}
REFERENCE_ARM = "Untreated"
STRUCT_COL = "focus_id"

# ---- phenotypes -------------------------------------------------------------
IDO1_POS = "CD68+IDO1+ Macrophages"
IDO1_NEG = "CD68+IDO1- Macrophages"
MACROPHAGE_ALL = [IDO1_POS, IDO1_NEG, "CD163+ Macrophages"]

# the pooled class. None of these is used to define foci.
LYMPHOCYTES = ["Helper T cells", "CD4- T cells", "Tregs", "B cells",
               "Plasma cells"]
T_LINEAGE = ["Helper T cells", "CD4- T cells", "Tregs"]
B_LINEAGE = ["B cells", "Plasma cells"]

# populations that DEFINE the foci; their radial results are circular
DETECTION_POOL = [IDO1_POS, IDO1_NEG, "CD163+ Macrophages", "Neutrophils"]

PHENOTYPE_COLORS = {
    IDO1_POS: "#B2182B", IDO1_NEG: "#EF8A62", "CD163+ Macrophages": "#FDBE85",
    "Neutrophils": "#7B3294", "Helper T cells": "#1B7837",
    "CD4- T cells": "#7FBC41", "Tregs": "#00441B", "B cells": "#2166AC",
    "Plasma cells": "#67A9CF", "Endothelial cells": "#E7298A",
    "Epithelial/Tumor cells": "#66C2A5", "Other": "#A6761D",
}

# ---- distance work ----------------------------------------------------------
RUN_DISTANCE = True
NN_ANCHORS = [IDO1_POS, IDO1_NEG]
NN_TARGETS = LYMPHOCYTES

# ---- models -----------------------------------------------------------------
MODEL_MAX_CELLS_PER_STRUCTURE = 3000
MODEL_MIN_CELLS_PER_ANIMAL = 20
BH_ALPHA = 0.10
RANDOM_SEED = 0

# ---- radial binning for profiles -------------------------------------------
N_RADIAL_BINS = 10          # core only, 0 to 1
MIN_CELLS_PER_BIN = 15

# ---- spillover control ------------------------------------------------------
RUN_SPILLOVER_CONTROL = True
INOS_RAW = "iNOS: Membrane: Mean"
ARG1_RAW = "Arginase-1: Cytoplasm: Mean"
# pairs that should NOT co-occur within a single macrophage
CONTROL_PAIRS = [
    ("CD3e: Membrane: Mean", "Pan-Cytokeratin: Membrane: Mean"),
    ("CD3e: Membrane: Mean", "CD20: Membrane: Mean"),
    ("CD45: Membrane: Mean", "Pan-Cytokeratin: Membrane: Mean"),
]
PHENOTYPE_COL_RAW = "Phenotypes"
COEXPRESSION_PERCENTILES = [50, 60, 70, 75, 80, 90]
COEXPRESSION_PRIMARY = 75

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


def ascii_safe(s):
    return "".join(ch for ch in str(s) if ord(ch) <= 127).strip()


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


def fit_mixed(df, outcome, label, covariates=None, note=""):
    """
    outcome ~ arm [+ covariates], random intercept for animal,
    structure nested within animal. Falls back to animal-only on failure.
    """
    if not HAVE_SM:
        return None
    covariates = covariates or []
    need = [outcome, "condition", "sample_id", STRUCT_COL] + covariates
    if any(c not in df.columns for c in need):
        return None
    d = df.dropna(subset=need).copy()
    counts = d.groupby("sample_id").size()
    d = d.loc[d["sample_id"].isin(counts.loc[counts >= MODEL_MIN_CELLS_PER_ANIMAL].index)]
    if not len(d) or d["condition"].nunique() < 2 or d["sample_id"].nunique() < 3:
        return None
    d["arm"] = (d["condition"] != REFERENCE_ARM).astype(float)
    d["struct_key"] = (d["sample_id"].astype(str) + "_"
                       + d[STRUCT_COL].astype(int).astype(str))
    d["_y"] = pd.to_numeric(d[outcome], errors="coerce")
    d = d.dropna(subset=["_y"])
    if not len(d):
        return None

    terms = ["arm"]
    for i, c in enumerate(covariates):
        if d[c].dtype == object or str(d[c].dtype).startswith("category"):
            d[f"_c{i}"] = d[c].astype(str)
            terms.append(f"C(_c{i})")
        else:
            v = pd.to_numeric(d[c], errors="coerce")
            sd = v.std()
            d[f"_c{i}"] = (v - v.mean()) / sd if sd and sd > 0 else 0.0
            terms.append(f"_c{i}")
    formula = "_y ~ " + " + ".join(terms)

    res, mode = None, ""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            md = smf.mixedlm(formula, data=d, groups=d["sample_id"],
                             re_formula="1",
                             vc_formula={"struct": "0 + C(struct_key)"})
            res = md.fit(reml=True, method="lbfgs", maxiter=300)
            mode = "animal + structure"
        except Exception:
            res = None
        if res is None or not np.isfinite(res.params.get("arm", np.nan)):
            try:
                md = smf.mixedlm(formula, data=d, groups=d["sample_id"],
                                 re_formula="1")
                res = md.fit(reml=True, method="lbfgs", maxiter=300)
                mode = "animal only"
            except Exception:
                return None
    if res is None:
        return None

    means = d.groupby("condition")["_y"].mean()
    coef = float(res.params.get("arm", np.nan))
    se = float(res.bse.get("arm", np.nan))
    return {
        "analysis": label, "outcome": outcome,
        "covariates": ",".join(covariates), "random_effects": mode,
        "n_cells": int(len(d)), "n_animals": int(d["sample_id"].nunique()),
        "n_animals_D1MT": int(d.loc[d["condition"] == "D1MT", "sample_id"].nunique()),
        "n_animals_ref": int(d.loc[d["condition"] == REFERENCE_ARM, "sample_id"].nunique()),
        "n_structures": int(d["struct_key"].nunique()),
        "mean_D1MT": float(means.get("D1MT", np.nan)),
        "mean_Untreated": float(means.get("Untreated", np.nan)),
        "coef_D1MT_vs_ref": coef, "std_err": se,
        "ci_low": coef - 1.96 * se, "ci_high": coef + 1.96 * se,
        "p_value": float(res.pvalues.get("arm", np.nan)),
        "note": note,
    }


_tee = Tee(os.path.join(TAB_DIR, "00_lymphocyte_radial_report.txt"))
sys.stdout = _tee

banner("AKOYA LYMPHOCYTE RADIAL POSITION - PRIMARY TEST")
print(f"Run time : {datetime.now().isoformat(timespec='seconds')}")
print(f"Input    : {IN_DIR}")
print(f"Output   : {OUT_DIR}")
if not HAVE_SCIPY:
    print(f"\nERROR: scipy required ({_scipy_err}).")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)
if not HAVE_SM:
    print(f"\n    WARNING: statsmodels unavailable ({_sm_err}). Models skipped.")

print("\nOUTCOME")
print("    centered radial position = radial_pos minus the mean radial position")
print("    of ALL cells in that structure. This is exactly the per-cell form of")
print("    the delta effect size, because a random subset's expected mean is the")
print("    structure mean. Negative = closer to the core than that structure's")
print("    own cells generally are.")


# %% Cell 3 - load and center
# =============================================================================

banner("LOADING AND CENTERING")

paths = sorted(glob.glob(os.path.join(CELL_DIR, "*_cell_structures.csv")))
if not paths:
    print(f"ERROR: no cell assignment files in {CELL_DIR}.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

frames = []
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
    d = d.loc[d[STRUCT_COL] > 0].copy()
    d["radial_pos"] = pd.to_numeric(d["radial_pos"], errors="coerce")
    d = d.loc[np.isfinite(d["radial_pos"])]
    if not len(d):
        continue

    # centering uses ALL cells in the structure, core and cuff separately so
    # the two coordinate systems are never mixed
    d["is_core"] = d["region"] == "core"
    for scope, mask in [("core", d["is_core"]), ("all", pd.Series(True, index=d.index))]:
        sub_d = d.loc[mask]
        if not len(sub_d):
            continue
        means = sub_d.groupby(STRUCT_COL)["radial_pos"].transform("mean")
        col = "radial_centered_core" if scope == "core" else "radial_centered_all"
        d.loc[mask, col] = sub_d["radial_pos"] - means
    frames.append(d)
    print(f"    {sid:<12} {d['condition'].iloc[0]:<10} {len(d):>8,} cells in "
          f"{int(d[STRUCT_COL].nunique()):>3} structures "
          f"({int(d['is_core'].sum()):>7,} core)")

cells = pd.concat(frames, ignore_index=True)
del frames
gc.collect()

cond_rank = {c: i for i, c in enumerate(CONDITION_ORDER)}
SAMPLE_ORDER = sorted(cells["sample_id"].unique(),
                      key=lambda s: (cond_rank.get(
                          cells.loc[cells["sample_id"] == s, "condition"].iloc[0], 9), s))
COND_OF = {s: cells.loc[cells["sample_id"] == s, "condition"].iloc[0]
           for s in SAMPLE_ORDER}
MARKER_OF = {s: POINT_MARKERS[i % len(POINT_MARKERS)] for i, s in enumerate(SAMPLE_ORDER)}
SHADES = {"D1MT": ["#08519C", "#3182BD", "#6BAED6"],
          "Untreated": ["#A63603", "#E6550D", "#FD8D3C"]}
COLOR_OF, _seen = {}, {c: 0 for c in CONDITION_ORDER}
for s in SAMPLE_ORDER:
    c = COND_OF[s]
    pal = SHADES.get(c, ["#999999"])
    COLOR_OF[s] = pal[_seen.get(c, 0) % len(pal)]
    _seen[c] = _seen.get(c, 0) + 1

core = cells.loc[cells["is_core"]].copy()
core["lineage"] = np.where(core["pheno"].isin(T_LINEAGE), "T lineage",
                           np.where(core["pheno"].isin(B_LINEAGE), "B lineage",
                                    "other"))
lym = core.loc[core["pheno"].isin(LYMPHOCYTES)].copy()

sub("Lymphocyte cell counts in focus cores")
tab = (lym.groupby(["sample_id", "pheno"]).size().unstack(fill_value=0)
       .reindex(index=SAMPLE_ORDER, columns=LYMPHOCYTES, fill_value=0))
print(tab.to_string())
print(f"\n    pooled lymphocytes: "
      f"{int(lym.loc[lym['condition'] == 'D1MT'].shape[0]):,} treated, "
      f"{int(lym.loc[lym['condition'] == 'Untreated'].shape[0]):,} untreated")
write_csv(tab.reset_index(), "70_lymphocyte_counts_core.csv")


# %% Cell 4 - sanity check: does centering reproduce the delta effect sizes?
# =============================================================================

banner("SANITY CHECK - CENTERING REPRODUCES DELTA")

print("    Per-animal mean of the centered radial position should match the")
print("    permutation delta reported by script 06, since delta is exactly this")
print("    quantity. Any large disagreement means something is wrong.\n")

check_rows = []
for s in SAMPLE_ORDER:
    d = core.loc[core["sample_id"] == s]
    for p in LYMPHOCYTES + [IDO1_POS, IDO1_NEG, "Neutrophils"]:
        v = d.loc[d["pheno"] == p, "radial_centered_core"].dropna()
        if not len(v):
            continue
        check_rows.append({
            "sample_id": s, "animal_id": short_label(s),
            "condition": COND_OF[s], "phenotype": p, "n_cells": len(v),
            "mean_centered_radial": float(v.mean()),
            "median_centered_radial": float(v.median()),
        })
chk = pd.DataFrame(check_rows)
write_csv(chk, "71_centered_radial_per_animal.csv")

sub("Mean centered radial position, core cells only")
print(f"    {'phenotype':<26}" + "".join(f"{short_label(s):>12}" for s in SAMPLE_ORDER))
print("    " + "-" * (26 + 12 * len(SAMPLE_ORDER)))
for p in LYMPHOCYTES + [IDO1_POS, IDO1_NEG, "Neutrophils"]:
    row = f"    {p:<26}"
    for s in SAMPLE_ORDER:
        v = chk.loc[(chk["sample_id"] == s) & (chk["phenotype"] == p),
                    "mean_centered_radial"]
        row += f"{v.iloc[0]:>+12.3f}" if len(v) else f"{'na':>12}"
    print(row)


# %% Cell 5 - PRIMARY MODEL and the lineage split
# =============================================================================

banner("PRIMARY MODEL - POOLED LYMPHOCYTES")

model_rows = []
if HAVE_SM:
    def capped(df):
        return pd.concat(
            [g.sample(min(len(g), MODEL_MAX_CELLS_PER_STRUCTURE),
                      random_state=RANDOM_SEED)
             for _, g in df.groupby(["sample_id", STRUCT_COL])],
            ignore_index=True) if len(df) else df

    # ---- PRIMARY -----------------------------------------------------------
    r = fit_mixed(capped(lym), "radial_centered_core",
                  "PRIMARY: pooled lymphocytes (core)",
                  covariates=["pheno"],
                  note="cell type as covariate so composition cannot drive it")
    if r:
        r["family"] = "primary"
        model_rows.append(r)
        print(f"    coefficient  : {r['coef_D1MT_vs_ref']:+.4f} radial units")
        print(f"    95% CI       : [{r['ci_low']:+.4f}, {r['ci_high']:+.4f}]")
        print(f"    p            : {r['p_value']:.4f}")
        print(f"    cells        : {r['n_cells']:,} in {r['n_structures']} structures")
        print(f"    animals      : {r['n_animals_D1MT']} treated, "
              f"{r['n_animals_ref']} untreated")
        print(f"    random effect: {r['random_effects']}")
        print("\n    Negative means lymphocytes sit closer to the focus core in")
        print("    treated animals than untreated, after accounting for where")
        print("    each structure's cells sit in general.")

    # ---- lineage split -----------------------------------------------------
    sub("Lineage split")
    print("    The fold-6 run showed B lineage surviving BH while T cells did")
    print("    not. This tests whether the pooled effect is class-wide or")
    print("    carried by B lineage alone.\n")
    for name, members in [("T lineage", T_LINEAGE), ("B lineage", B_LINEAGE)]:
        d = lym.loc[lym["pheno"].isin(members)]
        r = fit_mixed(capped(d), "radial_centered_core", f"{name} (core)",
                      covariates=["pheno"])
        if r:
            r["family"] = "lineage"
            model_rows.append(r)
            print(f"    {name:<12} coef={r['coef_D1MT_vs_ref']:>+8.4f}  "
                  f"p={r['p_value']:.4f}  cells={r['n_cells']:>7,}")

    # ---- per cell type -----------------------------------------------------
    sub("Per cell type (core only, centered outcome)")
    for p in LYMPHOCYTES + [IDO1_POS, IDO1_NEG, "Neutrophils"]:
        d = core.loc[core["pheno"] == p]
        circ = p in DETECTION_POOL
        r = fit_mixed(capped(d), "radial_centered_core", f"radial: {p}",
                      note="CIRCULAR: defines the foci" if circ else "")
        if r:
            r["family"] = "per_phenotype"
            r["phenotype"] = p
            r["circular"] = circ
            model_rows.append(r)
            print(f"    {p:<26} coef={r['coef_D1MT_vs_ref']:>+8.4f}  "
                  f"p={r['p_value']:.4f}  cells={r['n_cells']:>7,}"
                  f"{'  [circular]' if circ else ''}")

    # ---- secondary: core plus cuff ----------------------------------------
    sub("Secondary: core plus cuff")
    print("    Cuff cells span a FIXED 150 um band that is not size-normalised,")
    print("    so this reintroduces structure size. Reported for completeness.")
    lym_all = cells.loc[cells["pheno"].isin(LYMPHOCYTES)]
    r = fit_mixed(capped(lym_all), "radial_centered_all",
                  "SECONDARY: pooled lymphocytes (core + cuff)",
                  covariates=["pheno"],
                  note="cuff band is not size-normalised")
    if r:
        r["family"] = "secondary"
        model_rows.append(r)
        print(f"    coef={r['coef_D1MT_vs_ref']:>+8.4f}  p={r['p_value']:.4f}  "
              f"cells={r['n_cells']:,}")


# %% Cell 6 - sensitivity analyses
# =============================================================================

banner("SENSITIVITY ANALYSES")

sens_rows = []
if HAVE_SM:
    def capped(df):
        return pd.concat(
            [g.sample(min(len(g), MODEL_MAX_CELLS_PER_STRUCTURE),
                      random_state=RANDOM_SEED)
             for _, g in df.groupby(["sample_id", STRUCT_COL])],
            ignore_index=True) if len(df) else df

    sub("Leave one animal out")
    print("    31438 behaves like a treated animal on most architectural")
    print("    measures and has the lowest untreated burden (7.3%). If the")
    print("    result depends on excluding or including one animal, that must")
    print("    be visible.\n")
    for drop in SAMPLE_ORDER:
        d = lym.loc[lym["sample_id"] != drop]
        r = fit_mixed(capped(d), "radial_centered_core",
                      f"drop {short_label(drop)}", covariates=["pheno"])
        if r:
            r["family"] = "leave_one_out"
            r["dropped"] = short_label(drop)
            r["dropped_arm"] = COND_OF[drop]
            sens_rows.append(r)
            print(f"    without {short_label(drop):<8} ({COND_OF[drop]:<10}) "
                  f"coef={r['coef_D1MT_vs_ref']:>+8.4f}  p={r['p_value']:.4f}  "
                  f"animals={r['n_animals_D1MT']}v{r['n_animals_ref']}")

    sub("Marker robustness")
    print("    CD4 has a between-slide variance fraction of 0.879, the least")
    print("    comparable marker carrying any of our populations. The effect")
    print("    should survive without helper T cells.\n")
    for name, members in [
        ("without Helper T", [p for p in LYMPHOCYTES if p != "Helper T cells"]),
        ("without Tregs", [p for p in LYMPHOCYTES if p != "Tregs"]),
        ("B + plasma only", B_LINEAGE),
        ("T cells only", T_LINEAGE),
    ]:
        d = lym.loc[lym["pheno"].isin(members)]
        r = fit_mixed(capped(d), "radial_centered_core", name,
                      covariates=["pheno"])
        if r:
            r["family"] = "marker_robustness"
            sens_rows.append(r)
            print(f"    {name:<20} coef={r['coef_D1MT_vs_ref']:>+8.4f}  "
                  f"p={r['p_value']:.4f}  cells={r['n_cells']:>7,}")

models = pd.DataFrame(model_rows + sens_rows)
if len(models):
    models["q_value"] = np.nan
    for fam, g in models.groupby("family"):
        models.loc[g.index, "q_value"] = benjamini_hochberg(g["p_value"].to_numpy())
    write_csv(models, "72_models_all.csv")


# %% Cell 7 - radial profiles and per-animal effect sizes
# =============================================================================

banner("RADIAL PROFILES, CORE ONLY")

edges = np.linspace(0, 1, N_RADIAL_BINS + 1)
centres = 0.5 * (edges[:-1] + edges[1:])
prof_rows = []
for s in SAMPLE_ORDER:
    d = core.loc[core["sample_id"] == s]
    b = np.clip(np.digitize(d["radial_pos"], edges[1:-1]), 0, N_RADIAL_BINS - 1)
    d = d.assign(_bin=b)
    for bi in range(N_RADIAL_BINS):
        g = d.loc[d["_bin"] == bi]
        if len(g) < MIN_CELLS_PER_BIN:
            continue
        rec = {"sample_id": s, "condition": COND_OF[s], "bin": bi,
               "radial_centre": float(centres[bi]), "n_cells": len(g)}
        for p in LYMPHOCYTES:
            rec[p] = 100.0 * float((g["pheno"] == p).mean())
        rec["pooled_lymphocytes"] = 100.0 * float(g["pheno"].isin(LYMPHOCYTES).mean())
        rec["B lineage"] = 100.0 * float(g["pheno"].isin(B_LINEAGE).mean())
        rec["T lineage"] = 100.0 * float(g["pheno"].isin(T_LINEAGE).mean())
        prof_rows.append(rec)
prof = pd.DataFrame(prof_rows)
if len(prof):
    write_csv(prof, "73_radial_profiles_core.csv")

fig, axes = plt.subplots(1, 3, figsize=(42, 13))
for ax, col, ttl in zip(axes,
                        ["pooled_lymphocytes", "B lineage", "T lineage"],
                        ["All lymphocytes", "B lineage", "T lineage"]):
    for s in SAMPLE_ORDER:
        d = prof.loc[prof["sample_id"] == s].sort_values("radial_centre")
        if not len(d):
            continue
        ax.plot(d["radial_centre"], d[col], linewidth=5, marker=MARKER_OF[s],
                markersize=16, color=COLOR_OF[s], alpha=0.9)
    ax.set_xlabel("radial position (0 = core centre, 1 = boundary)",
                  fontsize=FONT_SIZE_BASE - 10)
    ax.set_ylabel("% of cells in bin", fontsize=FONT_SIZE_BASE - 10)
    ax.set_title(ttl, fontsize=FONT_SIZE_TITLE - 10)
    ax.set_ylim(bottom=0)
    style_axes(ax)
axes[0].legend(handles=[Line2D([0], [0], color=COLOR_OF[s], marker=MARKER_OF[s],
                               markersize=16, linewidth=4,
                               label=f"{short_label(s)} ({COND_OF[s]})")
                        for s in SAMPLE_ORDER],
               frameon=False, fontsize=FONT_SIZE_LEGEND - 16, loc="best")
fig.suptitle("Lymphocyte abundance across the focus core\n"
             "Rising toward the left means lymphocytes concentrate at the core",
             y=1.04, fontsize=FONT_SIZE_TITLE - 6)
save_fig(fig, "F51_radial_profiles_core")

# ---- F52 per-animal centered position --------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(32, 14),
                         gridspec_kw={"width_ratios": [1.3, 1.0]})
show_ph = LYMPHOCYTES + [IDO1_POS, "Neutrophils"]
ax = axes[0]
yy = np.arange(len(show_ph))
for i, p in enumerate(show_ph):
    for s in SAMPLE_ORDER:
        v = chk.loc[(chk["sample_id"] == s) & (chk["phenotype"] == p)]
        if not len(v):
            continue
        off = (SAMPLE_ORDER.index(s) - (len(SAMPLE_ORDER) - 1) / 2) * 0.13
        thin = int(v["n_cells"].iloc[0]) < 50
        ax.scatter(v["mean_centered_radial"].iloc[0], i + off, s=420,
                   color=COLOR_OF[s], marker=MARKER_OF[s],
                   edgecolor=FLAG_COLOR if thin else "#FFFFFF",
                   linewidth=3 if thin else 2, zorder=3)
ax.axvline(0, color="#000000", linewidth=3.5)
ax.set_yticks(yy)
ax.set_yticklabels([p + (" [circ]" if p in DETECTION_POOL else "")
                    for p in show_ph], fontsize=FONT_SIZE_TICK - 8)
ax.invert_yaxis()
ax.set_xlabel("mean centered radial position\n(negative = core-ward)")
ax.set_title("Per animal", fontsize=FONT_SIZE_TITLE - 10)
ax.xaxis.grid(True, color=GRID_COLOR, linewidth=1.5); ax.set_axisbelow(True)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
h = [Line2D([0], [0], color=COLOR_OF[s], marker=MARKER_OF[s], markersize=16,
            linestyle="none", label=f"{short_label(s)} ({COND_OF[s]})")
     for s in SAMPLE_ORDER]
h.append(Line2D([0], [0], color="#FFFFFF", marker="o", markersize=16,
                markeredgecolor=FLAG_COLOR, markeredgewidth=3, linestyle="none",
                label="< 50 cells"))
ax.legend(handles=h, frameon=False, fontsize=FONT_SIZE_LEGEND - 16, loc="best")

ax = axes[1]
if len(models):
    g = models.loc[models["family"].isin(["primary", "lineage", "per_phenotype"])]
    g = g.sort_values("coef_D1MT_vs_ref")
    yy2 = np.arange(len(g))
    for i, (_, r) in enumerate(g.iterrows()):
        sig = r["p_value"] < 0.05
        col = FLAG_COLOR if sig else "#777777"
        ax.plot([r["ci_low"], r["ci_high"]], [i, i], color=col, linewidth=5,
                zorder=2)
        ax.scatter([r["coef_D1MT_vs_ref"]], [i], s=380, color=col,
                   edgecolor="#FFFFFF", linewidth=2, zorder=3)
        ax.text(r["ci_high"], i, f"  p={r['p_value']:.3f}", va="center",
                fontsize=FONT_SIZE_ANNOT - 12)
    ax.axvline(0, color="#000000", linewidth=3.5)
    ax.set_yticks(yy2)
    ax.set_yticklabels([a[:36] for a in g["analysis"]],
                       fontsize=FONT_SIZE_TICK - 14)
    ax.invert_yaxis()
    ax.set_xlabel("D1MT effect (95% CI)")
    ax.set_title("Models", fontsize=FONT_SIZE_TITLE - 10)
    ax.xaxis.grid(True, color=GRID_COLOR, linewidth=1.5); ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
fig.suptitle("Centered radial position: effect sizes and models\n"
             "Outcome is per-cell delta, so structure size and shape are "
             "already removed", y=1.04, fontsize=FONT_SIZE_TITLE - 6)
save_fig(fig, "F52_centered_radial_effects")

# ---- F53 sensitivity --------------------------------------------------------
if len(models) and (models["family"] == "leave_one_out").any():
    fig, axes = plt.subplots(1, 2, figsize=(30, 13))
    prim = models.loc[models["family"] == "primary"]
    ax = axes[0]
    g = models.loc[models["family"] == "leave_one_out"]
    yy3 = np.arange(len(g))
    for i, (_, r) in enumerate(g.iterrows()):
        col = CONDITION_COLORS.get(r.get("dropped_arm", ""), "#777777")
        ax.plot([r["ci_low"], r["ci_high"]], [i, i], color=col, linewidth=5)
        ax.scatter([r["coef_D1MT_vs_ref"]], [i], s=380, color=col,
                   edgecolor="#FFFFFF", linewidth=2, zorder=3)
        ax.text(r["ci_high"], i, f"  p={r['p_value']:.3f}", va="center",
                fontsize=FONT_SIZE_ANNOT - 12)
    if len(prim):
        ax.axvline(float(prim["coef_D1MT_vs_ref"].iloc[0]), color="#000000",
                   linestyle="--", linewidth=3)
    ax.axvline(0, color="#000000", linewidth=3.5)
    ax.set_yticks(yy3)
    ax.set_yticklabels([f"without {r['dropped']}" for _, r in g.iterrows()],
                       fontsize=FONT_SIZE_TICK - 10)
    ax.invert_yaxis()
    ax.set_xlabel("D1MT effect (95% CI)")
    ax.set_title("Leave one animal out\n(dashed = full model)",
                 fontsize=FONT_SIZE_TITLE - 12)
    ax.xaxis.grid(True, color=GRID_COLOR, linewidth=1.5); ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    ax = axes[1]
    g = models.loc[models["family"] == "marker_robustness"]
    yy4 = np.arange(len(g))
    for i, (_, r) in enumerate(g.iterrows()):
        col = FLAG_COLOR if r["p_value"] < 0.05 else "#777777"
        ax.plot([r["ci_low"], r["ci_high"]], [i, i], color=col, linewidth=5)
        ax.scatter([r["coef_D1MT_vs_ref"]], [i], s=380, color=col,
                   edgecolor="#FFFFFF", linewidth=2, zorder=3)
        ax.text(r["ci_high"], i, f"  p={r['p_value']:.3f}", va="center",
                fontsize=FONT_SIZE_ANNOT - 12)
    if len(prim):
        ax.axvline(float(prim["coef_D1MT_vs_ref"].iloc[0]), color="#000000",
                   linestyle="--", linewidth=3)
    ax.axvline(0, color="#000000", linewidth=3.5)
    ax.set_yticks(yy4)
    ax.set_yticklabels([a[:28] for a in g["analysis"]],
                       fontsize=FONT_SIZE_TICK - 10)
    ax.invert_yaxis()
    ax.set_xlabel("D1MT effect (95% CI)")
    ax.set_title("Marker robustness", fontsize=FONT_SIZE_TITLE - 12)
    ax.xaxis.grid(True, color=GRID_COLOR, linewidth=1.5); ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.suptitle("Does the result depend on one animal or one marker?", y=1.04,
                 fontsize=FONT_SIZE_TITLE - 6)
    save_fig(fig, "F53_sensitivity")


# %% Cell 8 - delta distance, normalised family dropped
# =============================================================================

banner("DELTA DISTANCE (normalised family dropped)")

dist_models = pd.DataFrame()
if RUN_DISTANCE and os.path.exists(NN_NULL_TABLE):
    nulls = pd.read_csv(NN_NULL_TABLE)
    key = ["sample_id", "structure_id", "anchor", "target"]
    if all(c in nulls.columns for c in key + ["null_median_um"]):
        rows = []
        for (s, k), g in core.groupby(["sample_id", STRUCT_COL]):
            xy = g[["x", "y"]].to_numpy(float)
            ph = g["pheno"].to_numpy()
            for a in NN_ANCHORS:
                ai = np.flatnonzero(ph == a)
                if not len(ai):
                    continue
                for t in NN_TARGETS:
                    ti = np.flatnonzero(ph == t)
                    if not len(ti):
                        continue
                    d, _ = cKDTree(xy[ti]).query(xy[ai], k=1)
                    nm = nulls.loc[(nulls["sample_id"] == s)
                                   & (nulls["structure_id"] == k)
                                   & (nulls["anchor"] == a)
                                   & (nulls["target"] == t), "null_median_um"]
                    if not len(nm):
                        continue
                    rows.append(pd.DataFrame({
                        "sample_id": s, "condition": COND_OF[s], STRUCT_COL: k,
                        "anchor": a, "target": t,
                        "distance_delta_um": np.asarray(d, float) - float(nm.iloc[0]),
                    }))
        if rows:
            dd = pd.concat(rows, ignore_index=True)
            out = []
            for (a, t), g in dd.groupby(["anchor", "target"]):
                r = fit_mixed(g, "distance_delta_um", f"{a} -> {t}",
                              note="delta outcome; normalised family dropped")
                if r:
                    r["anchor"], r["target"] = a, t
                    out.append(r)
            # pooled across lymphocyte targets
            for a in NN_ANCHORS:
                g = dd.loc[dd["anchor"] == a]
                r = fit_mixed(g, "distance_delta_um",
                              f"POOLED: {a} -> all lymphocytes",
                              covariates=["target"])
                if r:
                    r["anchor"], r["target"] = a, "pooled"
                    out.append(r)
            dist_models = pd.DataFrame(out)
            if len(dist_models):
                dist_models["q_value"] = benjamini_hochberg(
                    dist_models["p_value"].to_numpy())
                write_csv(dist_models, "74_delta_distance_models.csv")
                sub("Delta distance models, core cells only")
                print(f"    {'test':<48}{'coef um':>10}{'p':>9}{'q':>9}")
                for _, r in dist_models.sort_values("p_value").iterrows():
                    print(f"    {r['analysis'][:47]:<48}"
                          f"{r['coef_D1MT_vs_ref']:>+10.2f}"
                          f"{r['p_value']:>9.4f}{r['q_value']:>9.4f}")
    else:
        print(f"    WARNING: {NN_NULL_TABLE} lacks expected columns. Skipped.")
else:
    print(f"    SKIPPED: {NN_NULL_TABLE} not found, or RUN_DISTANCE is False.")


# %% Cell 9 - spillover control for co-expression
# =============================================================================

banner("SPILLOVER CONTROL FOR iNOS / ARGINASE-1 CO-EXPRESSION")

coex = pd.DataFrame()
if RUN_SPILLOVER_CONTROL:
    print("    Segmentation spillover inflates apparent co-expression of ANY two")
    print("    markers, and worsens with cell density. Untreated foci are three")
    print("    to five fold denser. Control pairs cannot co-occur in a single")
    print("    macrophage, so if they separate by arm the same way iNOS and")
    print("    Arginase-1 do, the result is spillover.\n")

    rows = []
    for s in SAMPLE_ORDER:
        raw_path = os.path.join(DATA_DIR, f"{s}.csv")
        if not os.path.exists(raw_path):
            print(f"    WARNING: {raw_path} not found, skipping {s}")
            continue
        need = [PHENOTYPE_COL_RAW, INOS_RAW, ARG1_RAW]
        for a, b in CONTROL_PAIRS:
            need += [a, b]
        need = list(dict.fromkeys(need))
        try:
            hdr = pd.read_csv(raw_path, nrows=0).columns.tolist()
            use = [c for c in need if c in hdr]
            miss = [c for c in need if c not in hdr]
            if miss:
                print(f"    WARNING: {s} missing columns {miss}")
            d = pd.read_csv(raw_path, usecols=use, low_memory=False)
        except Exception as e:
            print(f"    ERROR reading {raw_path}: {e}")
            continue
        d["_pheno"] = d[PHENOTYPE_COL_RAW].map(ascii_safe)
        mac = d.loc[d["_pheno"].isin(MACROPHAGE_ALL)]
        if not len(mac):
            continue
        pairs = [("iNOS x Arginase-1", INOS_RAW, ARG1_RAW, "test")]
        pairs += [(f"{a.split(':')[0]} x {b.split(':')[0]}", a, b, "control")
                  for a, b in CONTROL_PAIRS]
        for name, ca, cb, kind in pairs:
            if ca not in mac.columns or cb not in mac.columns:
                continue
            va = pd.to_numeric(mac[ca], errors="coerce")
            vb = pd.to_numeric(mac[cb], errors="coerce")
            for pct in COEXPRESSION_PERCENTILES:
                ta, tb = np.nanpercentile(va, pct), np.nanpercentile(vb, pct)
                hi_a, hi_b = (va >= ta).to_numpy(), (vb >= tb).to_numpy()
                pa, pb = float(hi_a.mean()), float(hi_b.mean())
                obs = float((hi_a & hi_b).mean())
                exp = pa * pb
                rows.append({
                    "sample_id": s, "animal_id": short_label(s),
                    "condition": COND_OF[s], "pair": name, "kind": kind,
                    "percentile": pct, "n_macrophages": len(mac),
                    "observed_pct": 100.0 * obs,
                    "expected_pct": 100.0 * exp,
                    "coexpression_ratio": obs / exp if exp > 0 else np.nan,
                })
        print(f"    {s} done")
        del d, mac
        gc.collect()

    coex = pd.DataFrame(rows)
    if len(coex):
        write_csv(coex, "75_coexpression_with_controls.csv")

        sub(f"Co-expression ratio at the {COEXPRESSION_PRIMARY}th percentile")
        p0 = coex.loc[coex["percentile"] == COEXPRESSION_PRIMARY]
        piv = p0.pivot_table(index="pair", columns="animal_id",
                             values="coexpression_ratio")
        cols = [short_label(s) for s in SAMPLE_ORDER if short_label(s) in piv.columns]
        print(piv[cols].to_string())

        sub("Arm separation by pair and threshold")
        sep_rows = []
        for (pair, kind), g in coex.groupby(["pair", "kind"]):
            for pct in COEXPRESSION_PERCENTILES:
                gg = g.loc[g["percentile"] == pct]
                a = gg.loc[gg["condition"] == "D1MT", "coexpression_ratio"].dropna()
                b = gg.loc[gg["condition"] == "Untreated", "coexpression_ratio"].dropna()
                if not len(a) or not len(b):
                    continue
                sep = (a.max() < b.min()) or (b.max() < a.min())
                sep_rows.append({"pair": pair, "kind": kind, "percentile": pct,
                                 "d1mt_min": a.min(), "d1mt_max": a.max(),
                                 "untr_min": b.min(), "untr_max": b.max(),
                                 "complete_separation": sep})
        sepdf = pd.DataFrame(sep_rows)
        write_csv(sepdf, "76_coexpression_separation_with_controls.csv")
        summ = (sepdf.groupby(["pair", "kind"])["complete_separation"]
                .sum().reset_index(name="n_thresholds_separated"))
        summ["n_thresholds"] = len(COEXPRESSION_PERCENTILES)
        print(summ.to_string(index=False))

        sub("VERDICT ON THE CO-EXPRESSION RESULT")
        test_n = int(summ.loc[summ["kind"] == "test",
                              "n_thresholds_separated"].sum())
        ctrl_n = int(summ.loc[summ["kind"] == "control",
                              "n_thresholds_separated"].sum())
        n_ctrl_pairs = int((summ["kind"] == "control").sum())
        print(f"    iNOS x Arginase-1 separates at {test_n} of "
              f"{len(COEXPRESSION_PERCENTILES)} thresholds")
        print(f"    control pairs separate at {ctrl_n} of "
              f"{n_ctrl_pairs * len(COEXPRESSION_PERCENTILES)} pair-thresholds")
        if ctrl_n == 0 and test_n > 0:
            print("\n    Controls are flat and the test pair separates. The")
            print("    co-expression result is NOT explained by segmentation")
            print("    spillover and can be reported.")
        elif ctrl_n > 0:
            print("\n    Control pairs ALSO separate by arm. Apparent")
            print("    co-expression tracks something global, most likely")
            print("    segmentation spillover driven by the density difference")
            print("    between arms. The iNOS / Arginase-1 result must be")
            print("    dropped or reported with this caveat stated plainly.")

        # ---- F54 -----------------------------------------------------------
        fig, axes = plt.subplots(1, 2, figsize=(32, 13),
                                 gridspec_kw={"width_ratios": [1.3, 1.0]})
        ax = axes[0]
        pair_names = list(dict.fromkeys(coex["pair"]))
        styles = {"test": "-", "control": "--"}
        for pair in pair_names:
            kind = coex.loc[coex["pair"] == pair, "kind"].iloc[0]
            for s in SAMPLE_ORDER:
                d = coex.loc[(coex["pair"] == pair) & (coex["sample_id"] == s)]
                if not len(d):
                    continue
                d = d.sort_values("percentile")
                ax.plot(d["percentile"], d["coexpression_ratio"],
                        linestyle=styles[kind],
                        linewidth=5 if kind == "test" else 2.5,
                        marker=MARKER_OF[s] if kind == "test" else None,
                        markersize=15, color=COLOR_OF[s],
                        alpha=0.95 if kind == "test" else 0.45)
        ax.axhline(1.0, color="#000000", linestyle=":", linewidth=3)
        ax.set_xlabel("within-section percentile")
        ax.set_ylabel("co-expression ratio")
        ax.set_title("Solid = iNOS x Arginase-1, dashed = controls",
                     fontsize=FONT_SIZE_TITLE - 14)
        style_axes(ax)

        ax = axes[1]
        p0 = coex.loc[coex["percentile"] == COEXPRESSION_PRIMARY]
        xs = np.arange(len(pair_names))
        for i, pair in enumerate(pair_names):
            for ci, c in enumerate(CONDITION_ORDER):
                v = p0.loc[(p0["pair"] == pair) & (p0["condition"] == c),
                           "coexpression_ratio"].dropna()
                if not len(v):
                    continue
                off = (ci - 0.5) * 0.3
                ax.scatter(np.full(len(v), i + off), v, s=340,
                           color=CONDITION_COLORS[c], edgecolor="#FFFFFF",
                           linewidth=2, zorder=3)
        ax.axhline(1.0, color="#000000", linestyle=":", linewidth=3)
        ax.set_xticks(xs)
        ax.set_xticklabels([p.replace(" x ", "\nx ") for p in pair_names],
                           fontsize=FONT_SIZE_TICK - 14)
        ax.set_ylabel("co-expression ratio")
        ax.set_title(f"At the {COEXPRESSION_PRIMARY}th percentile",
                     fontsize=FONT_SIZE_TITLE - 14)
        style_axes(ax)
        ax.legend(handles=[Patch(facecolor=CONDITION_COLORS[c], label=c)
                           for c in CONDITION_ORDER],
                  frameon=False, fontsize=FONT_SIZE_LEGEND - 12, loc="best")
        fig.suptitle("Is the co-expression result real, or segmentation "
                     "spillover?\nControl pairs cannot co-occur in one "
                     "macrophage", y=1.04, fontsize=FONT_SIZE_TITLE - 8)
        save_fig(fig, "F54_spillover_control")


# %% Cell 10 - wrap up
# =============================================================================

banner("SUMMARY")

if len(models):
    prim = models.loc[models["family"] == "primary"]
    if len(prim):
        r = prim.iloc[0]
        print(f"PRIMARY RESULT")
        print(f"  Pooled lymphocytes, core cells, centered radial position")
        print(f"  coef {r['coef_D1MT_vs_ref']:+.4f} "
              f"[{r['ci_low']:+.4f}, {r['ci_high']:+.4f}]  p = {r['p_value']:.4f}")
        print(f"  {r['n_cells']:,} cells, {r['n_structures']} structures, "
              f"{r['n_animals_D1MT']} vs {r['n_animals_ref']} animals")

sub("How to report this")
print("  - The outcome is position, not intensity, so the slide confound does")
print("    not apply to the measurement. It DOES apply to cell identity, which")
print("    is why the lineage split and the without-Helper-T check matter.")
print("  - Lymphocytes do not define the foci, so unlike the macrophage and")
print("    neutrophil results these positions are not circular.")
print("  - With three animals per arm the p-value is bounded regardless of cell")
print("    count. The strength of the claim comes from the effect size, the")
print("    per-animal consistency, and the sensitivity analyses, not from p.")
print("  - Read the leave-one-out panel before quoting anything. If dropping")
print("    31438 changes the conclusion, say so.")

sub("Still open")
print("  1. 43111 has one focus against an expert count of two.")
print("  2. Untreated foci went from 38 to 44 in 43109 at fold 6, so untreated")
print("     over-segmentation increased. Structure is a random effect, so this")
print("     adds noise rather than bias, but it is worth a look at F33.")
print("  3. BALT organisation (Q7) is untouched. The detection is fixed and")
print("     table 36 now has 24 candidates, 18 CD21+.")

banner("DONE")
sys.stdout = _tee.terminal
_tee.close()
