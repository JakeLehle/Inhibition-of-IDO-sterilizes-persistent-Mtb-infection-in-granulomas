#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - LYMPHOCYTE RADIAL POSITION, THE PRIMARY TEST
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 08 of the AKOYA analysis series. REVISION 2.

THE HYPOTHESIS
    Lymphocytes sit closer to the core of a myeloid focus in D1MT-treated
    animals than in untreated animals. Because the outcome is position rather
    than intensity, it is immune to the slide confound that limits most of this
    panel: coordinates carry no batch signal. Cell IDENTITY still does, which is
    why the lineage split below matters.

    Lymphocytes are also not used to define foci (detection is myeloid-only), so
    unlike every macrophage and neutrophil result these positions are not
    circular.

WHAT THIS FIXES FROM SCRIPTS 06 AND 07 (unchanged, and now better supported)
    1. CENTERED RADIAL POSITION IS THE OUTCOME. Scripts 06 and 07 modelled raw
       radial_pos while reporting delta as the effect size. Those are different
       quantities, and CD4- T cells show it on the current rev4 run: per-animal
       pooled deltas of -0.226, -0.121, -0.096 against -0.052, +0.032, +0.166,
       a clean split, but a raw-radial model coefficient of -0.034 with
       p = 0.808.

       delta = mean radial of population P, minus the mean of a random subset of
       the same size. A random subset's expected mean is the structure's overall
       mean, so delta is exactly (radial_pos - structure mean radial). Centering
       each cell on its structure mean therefore makes the model test the same
       quantity the effect sizes report.

    2. CORE CELLS ONLY IN THE PRIMARY MODEL. Core cells are normalised 0 to 1 by
       each focus's own inscribed radius, but cuff cells run 1 to 2 across a
       FIXED 150 um band that is not size-normalised. Mixing them reintroduces
       structure size. Core-only is primary; core-plus-cuff is secondary.

    3. NORMALISED DISTANCE FAMILY DROPPED. Script 07 revision 2 measured this
       properly on both radius definitions. Correlation with focus radius:
       raw median |r| = 0.10, delta 0.11, normalised by equivalent radius 0.36,
       normalised by inscribed radius 0.30. All twenty normalised correlations
       are NEGATIVE, across both definitions and all ten pairs, which is
       systematic over-correction rather than residual confounding. Delta is
       retained and is the only distance outcome used here.

    4. ONE POOLED TEST INSTEAD OF TEN, with cell type as a covariate.

    5. SPILLOVER CONTROL FOR CO-EXPRESSION, using pairs that cannot co-occur in
       a single macrophage.

WHAT CHANGED IN REVISION 2 (and why)

    1. A COMPOSITION-BALANCED CENTRING IS RUN ALONGSIDE THE PRIMARY ONE.
       radial_centered_core subtracts the mean radial position of ALL core cells
       in the structure. That mean is cell-weighted, so it is dominated by
       whichever population is most abundant, and composition differs sharply
       between arms: untreated cores are heavily myeloid, treated cores are not.
       The reference point therefore moves with the arm, and the outcome reads
       "closer to the core than the average cell in this focus" rather than
       "closer to the core".

       The balanced version subtracts the UNWEIGHTED MEAN OF THE PER-PHENOTYPE
       MEANS within that structure, over phenotypes clearing
       MIN_CELLS_FOR_BALANCE cells there. Every population contributes equally
       to the reference, so a shift in composition cannot move it. If the
       coefficient survives that, the claim is about position rather than about
       what else is in the focus.

       The exact area-weighted reference, the mean of the radial map over all
       core PIXELS, is composition-free by construction and would be better
       still. It cannot be computed from the per-cell tables and needs one
       extra column exported from script 04. If mean_radial_over_area appears in
       table 35 this script will read and use it; otherwise it uses the balanced
       version and says so.

    2. BH FAMILIES ARE SPLIT PROPERLY.
       The per_phenotype family mixed the circular myeloid tests with the
       non-circular lymphocyte tests, inflating m and blending a descriptive
       family with a confirmatory one. They are now separate families.

    3. SENSITIVITY ANALYSES NO LONGER RECEIVE q-VALUES.
       leave_one_out and marker_robustness are refits of one hypothesis on
       subsets, not independent tests. BH across them is meaningless and invites
       misreading. Their q-values are NaN.

    4. THE SPILLOVER VERDICT COMPARES TEST AGAINST CONTROL, NOT AGAINST ZERO.
       The old rule dropped the co-expression result if ANY control pair
       separated at ANY threshold. With three control pairs at six thresholds
       and three animals per arm, complete separation arises by chance at
       p = 0.10 per test, so about 1.8 spurious control separations are expected
       even with a perfectly clean panel. The verdict now states the chance
       expectation and compares the test pair's separation count against the
       worst single control pair.

    5. MICRONS CONVERSION USES THE INSCRIBED RADIUS.
       The radial coordinate normalises by each focus's maximum inscribed
       radius, so converting a radial-unit coefficient with equiv_radius_um
       overstates it. Script 04 revision 4 exports max_inscribed_radius_um. On
       the current run the medians are 125 um treated and 168 um untreated
       against 138 and 201 equivalent, so a coefficient of -0.090 is about 11
       and 15 um rather than 12 and 18.

    6. PER-ANIMAL COMPLETE SEPARATION IS CHECKED AND REPORTED.
       Script 06's pooled delta shows complete separation between arms on every
       lymphocyte population. The same check is run here on the centred outcome,
       which is the quantity this script models, so the two can be compared
       rather than conflated.

    7. STALE NUMBERS CORRECTED. The prior-run results quoted in the header
       (q = 0.088 for plasma and B cells, IDO1- to plasma at p = 0.046) were
       from an earlier structure definition and are regenerated by this run.
       CD4's comparability is ICC 0.403 on p99 and 0.874 on section medians, not
       the single 0.879 previously quoted; see script 03 table 22.

INPUTS
    structures_rev4/cell_assignments/<section>_cell_structures.csv
    structures_rev4/tables/35_foci_structures_relative.csv
    distance_stats/tables/53_nn_per_structure.csv   (shared null, from 06)
    AKOYA/data/<section>.csv                        (spillover control markers)

OUTPUTS
    figures/  F51 .. F55
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

# ---- centring ---------------------------------------------------------------
# Primary keeps the cell-weighted reference for continuity with the delta effect
# sizes. Balanced weights every phenotype equally so composition cannot move the
# reference. If script 04 ever exports mean_radial_over_area, that is preferred
# over both and will be used automatically.
MIN_CELLS_FOR_BALANCE = 10
AREA_REFERENCE_COL = "mean_radial_over_area"

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
# families that are refits of one hypothesis, not independent tests
NO_BH_FAMILIES = ["leave_one_out", "marker_robustness", "centring_check"]

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
# with 3 animals per arm, complete separation occurs by chance at this rate
CHANCE_SEPARATION_RATE = 0.10

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


def complete_separation(a, b):
    """True when the two groups' ranges do not overlap at all."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if not len(a) or not len(b):
        return False
    return bool(a.max() < b.min() or b.max() < a.min())


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

banner("AKOYA LYMPHOCYTE RADIAL POSITION - PRIMARY TEST (revision 2)")
print(f"Run time : {datetime.now().isoformat(timespec='seconds')}")
print(f"Input    : {IN_DIR}")
print(f"Output   : {OUT_DIR}")
if not HAVE_SCIPY:
    print(f"\nERROR: scipy required ({_scipy_err}).")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)
if not HAVE_SM:
    print(f"\n    WARNING: statsmodels unavailable ({_sm_err}). Models skipped.")

print("\nOUTCOME")
print("    PRIMARY   centered radial position = radial_pos minus the mean radial")
print("              position of ALL core cells in that structure. This is the")
print("              per-cell form of the delta effect size, because a random")
print("              subset's expected mean is the structure mean.")
print("    BALANCED  the same, but the reference is the UNWEIGHTED MEAN OF THE")
print("              PER-PHENOTYPE MEANS. The primary reference is cell-weighted")
print("              and therefore moves with composition, which differs sharply")
print("              between arms. If the coefficient survives the balanced")
print("              version, the claim is about position rather than about what")
print("              else is in the focus.")
print("    Negative = closer to the core than that structure's own cells are.")


# %% Cell 3 - load, center two ways
# =============================================================================

banner("LOADING AND CENTERING")

# focus geometry, for the microns conversion
inscribed_of, equiv_of, area_ref_of = {}, {}, {}
HAVE_AREA_REF = False
if os.path.exists(FOCI_TABLE):
    ft = pd.read_csv(FOCI_TABLE)
    idcol = "focus_id" if "focus_id" in ft.columns else STRUCT_COL
    HAVE_AREA_REF = AREA_REFERENCE_COL in ft.columns
    for _, r in ft.iterrows():
        key = (r["sample_id"], int(r[idcol]))
        equiv_of[key] = float(r.get("equiv_radius_um", np.nan))
        inscribed_of[key] = float(r.get("max_inscribed_radius_um", np.nan))
        if HAVE_AREA_REF:
            area_ref_of[key] = float(r[AREA_REFERENCE_COL])
    print(f"    geometry for {len(equiv_of)} structures"
          f"{'' if np.isfinite(list(inscribed_of.values())).any() else ' (no inscribed radius)'}")
    if HAVE_AREA_REF:
        print(f"    '{AREA_REFERENCE_COL}' found: the exact area-weighted "
              f"reference will be used as a third centring.")
    else:
        print(f"    '{AREA_REFERENCE_COL}' not in table 35. The balanced "
              f"centring is the composition-free check.")

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

    # ---- composition-balanced reference, core only -------------------------
    core_d = d.loc[d["is_core"]]
    bal_ref = {}
    for k, g in core_d.groupby(STRUCT_COL):
        per_ph = g.groupby("pheno")["radial_pos"].agg(["mean", "size"])
        per_ph = per_ph.loc[per_ph["size"] >= MIN_CELLS_FOR_BALANCE]
        bal_ref[k] = (float(per_ph["mean"].mean()) if len(per_ph)
                      else float(g["radial_pos"].mean()))
    d.loc[d["is_core"], "radial_centered_core_balanced"] = (
        core_d["radial_pos"] - core_d[STRUCT_COL].map(bal_ref))
    d["balance_reference"] = d[STRUCT_COL].map(bal_ref)
    d["n_phenotypes_in_balance"] = d[STRUCT_COL].map(
        {k: int((g.groupby("pheno").size() >= MIN_CELLS_FOR_BALANCE).sum())
         for k, g in core_d.groupby(STRUCT_COL)})

    # ---- exact area-weighted reference, if script 04 exported it -----------
    if HAVE_AREA_REF:
        d.loc[d["is_core"], "radial_centered_core_area"] = (
            core_d["radial_pos"]
            - core_d[STRUCT_COL].map(
                lambda k: area_ref_of.get((sid, int(k)), np.nan)))

    d["inscribed_radius_um"] = d[STRUCT_COL].map(
        lambda k: inscribed_of.get((sid, int(k)), np.nan))
    d["equiv_radius_um"] = d[STRUCT_COL].map(
        lambda k: equiv_of.get((sid, int(k)), np.nan))

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

CENTRING_OUTCOMES = [("radial_centered_core", "primary (cell-weighted)")]
if "radial_centered_core_balanced" in core.columns:
    CENTRING_OUTCOMES.append(("radial_centered_core_balanced",
                              "balanced (phenotype-weighted)"))
if HAVE_AREA_REF and "radial_centered_core_area" in core.columns:
    CENTRING_OUTCOMES.append(("radial_centered_core_area",
                              "area-weighted (exact)"))

sub("How different are the two references?")
print("    If composition were the same in both arms these would agree. They")
print("    do not have to, and the size of the gap is the size of the concern.\n")
print(f"    {'section':<12}{'structs':>9}{'cell-wtd ref':>15}{'balanced ref':>15}"
      f"{'gap':>9}{'phenos':>8}")
print("    " + "-" * 68)
ref_rows = []
for s in SAMPLE_ORDER:
    g = core.loc[core["sample_id"] == s]
    for k, gg in g.groupby(STRUCT_COL):
        ref_rows.append({
            "sample_id": s, "condition": COND_OF[s], "structure_id": int(k),
            "n_core_cells": len(gg),
            "cell_weighted_reference": float(gg["radial_pos"].mean()),
            "balanced_reference": float(gg["balance_reference"].iloc[0]),
            "n_phenotypes_in_balance": int(gg["n_phenotypes_in_balance"].iloc[0]),
        })
refs = pd.DataFrame(ref_rows)
refs["reference_gap"] = (refs["cell_weighted_reference"]
                         - refs["balanced_reference"])
for s in SAMPLE_ORDER:
    g = refs.loc[refs["sample_id"] == s]
    print(f"    {s:<12}{len(g):>9}{g['cell_weighted_reference'].median():>15.3f}"
          f"{g['balanced_reference'].median():>15.3f}"
          f"{g['reference_gap'].median():>+9.3f}"
          f"{g['n_phenotypes_in_balance'].median():>8.0f}")
for c in CONDITION_ORDER:
    g = refs.loc[refs["condition"] == c]
    print(f"    {c:<12} median gap {g['reference_gap'].median():+.3f}")
write_csv(refs, "77_centring_references.csv")

sub("Lymphocyte cell counts in focus cores")
tab = (lym.groupby(["sample_id", "pheno"]).size().unstack(fill_value=0)
       .reindex(index=SAMPLE_ORDER, columns=LYMPHOCYTES, fill_value=0))
print(tab.to_string())
print(f"\n    pooled lymphocytes: "
      f"{int(lym.loc[lym['condition'] == 'D1MT'].shape[0]):,} treated, "
      f"{int(lym.loc[lym['condition'] == 'Untreated'].shape[0]):,} untreated")
write_csv(tab.reset_index(), "70_lymphocyte_counts_core.csv")


# %% Cell 4 - per-animal effect sizes and separation
# =============================================================================

banner("PER-ANIMAL CENTERED RADIAL POSITION")

print("    Per-animal mean of the centered radial position is the same quantity")
print("    script 06 reports as a permutation delta. Both centrings are shown.\n")

check_rows = []
for s in SAMPLE_ORDER:
    d = core.loc[core["sample_id"] == s]
    for p in LYMPHOCYTES + [IDO1_POS, IDO1_NEG, "Neutrophils"]:
        v = d.loc[d["pheno"] == p]
        if not len(v):
            continue
        rec = {
            "sample_id": s, "animal_id": short_label(s),
            "condition": COND_OF[s], "phenotype": p, "n_cells": len(v),
            "mean_centered_radial": float(v["radial_centered_core"].mean()),
            "median_centered_radial": float(v["radial_centered_core"].median()),
        }
        for col, _ in CENTRING_OUTCOMES[1:]:
            rec[f"mean_{col}"] = float(v[col].mean()) if col in v.columns else np.nan
        check_rows.append(rec)
chk = pd.DataFrame(check_rows)
write_csv(chk, "71_centered_radial_per_animal.csv")

for col, name in CENTRING_OUTCOMES:
    key = "mean_centered_radial" if col == "radial_centered_core" else f"mean_{col}"
    if key not in chk.columns:
        continue
    sub(f"Mean centered radial position, {name}")
    print(f"    {'phenotype':<26}" + "".join(f"{short_label(s):>12}" for s in SAMPLE_ORDER)
          + f"{'separated':>12}")
    print("    " + "-" * (26 + 12 * len(SAMPLE_ORDER) + 12))
    for p in LYMPHOCYTES + [IDO1_POS, IDO1_NEG, "Neutrophils"]:
        row = f"    {p:<26}"
        vals = {c: [] for c in CONDITION_ORDER}
        for s in SAMPLE_ORDER:
            v = chk.loc[(chk["sample_id"] == s) & (chk["phenotype"] == p), key]
            if len(v):
                row += f"{v.iloc[0]:>+12.3f}"
                vals[COND_OF[s]].append(float(v.iloc[0]))
            else:
                row += f"{'na':>12}"
        sep = complete_separation(vals[CONDITION_ORDER[0]], vals[CONDITION_ORDER[1]])
        row += f"{'YES' if sep else '':>12}"
        if p in DETECTION_POOL:
            row += "  [circular]"
        print(row)
    print(f"\n    Complete separation with three animals per arm arises by")
    print(f"    chance with probability {CHANCE_SEPARATION_RATE:.2f} per test, so")
    print(f"    read the column as a pattern across populations, not one by one.")


# %% Cell 5 - PRIMARY MODEL, centring check, lineage split
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

    lym_capped = capped(lym)

    # ---- PRIMARY -----------------------------------------------------------
    r = fit_mixed(lym_capped, "radial_centered_core",
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

        # microns conversion on the correct scale
        ins_t = np.nanmedian([v for (sid, k), v in inscribed_of.items()
                              if COND_OF.get(sid) == "D1MT"])
        ins_u = np.nanmedian([v for (sid, k), v in inscribed_of.items()
                              if COND_OF.get(sid) == "Untreated"])
        eq_t = np.nanmedian([v for (sid, k), v in equiv_of.items()
                             if COND_OF.get(sid) == "D1MT"])
        eq_u = np.nanmedian([v for (sid, k), v in equiv_of.items()
                             if COND_OF.get(sid) == "Untreated"])
        c = abs(r["coef_D1MT_vs_ref"])
        print(f"\n    IN MICRONS. The radial coordinate normalises by the maximum")
        print(f"    INSCRIBED radius, so that is the scale to convert with.")
        print(f"      inscribed  : {c * ins_t:.0f} um treated ({ins_t:.0f} um median), "
              f"{c * ins_u:.0f} um untreated ({ins_u:.0f} um)")
        print(f"      equivalent : {c * eq_t:.0f} um / {c * eq_u:.0f} um  "
              f"(overstates, do not quote)")

    # ---- centring check ----------------------------------------------------
    sub("Centring check: does the result survive a composition-free reference?")
    print("    The primary reference is cell-weighted and therefore moves with")
    print("    composition, which differs sharply between arms. If the balanced")
    print("    coefficient collapses, the effect was partly about what else is")
    print("    in the focus rather than about lymphocyte position.\n")
    for col, name in CENTRING_OUTCOMES:
        r = fit_mixed(lym_capped, col, f"pooled lymphocytes: {name}",
                      covariates=["pheno"], note=f"centring = {name}")
        if r:
            r["family"] = "centring_check"
            r["centring"] = name
            model_rows.append(r)
            print(f"    {name:<32} coef={r['coef_D1MT_vs_ref']:>+8.4f}  "
                  f"[{r['ci_low']:>+7.4f}, {r['ci_high']:>+7.4f}]  "
                  f"p={r['p_value']:.4f}")

    # ---- lineage split -----------------------------------------------------
    sub("Lineage split")
    print("    Tests whether the pooled effect is class-wide or carried by one")
    print("    lineage. Run on both centrings.\n")
    for name, members in [("T lineage", T_LINEAGE), ("B lineage", B_LINEAGE)]:
        d = lym.loc[lym["pheno"].isin(members)]
        for col, cname in CENTRING_OUTCOMES:
            r = fit_mixed(capped(d), col, f"{name} ({cname})",
                          covariates=["pheno"])
            if r:
                r["family"] = ("lineage" if col == "radial_centered_core"
                               else "lineage_centring_check")
                r["lineage"] = name
                r["centring"] = cname
                model_rows.append(r)
                print(f"    {name:<12} {cname:<32} "
                      f"coef={r['coef_D1MT_vs_ref']:>+8.4f}  "
                      f"p={r['p_value']:.4f}  cells={r['n_cells']:>7,}")

    # ---- per cell type, circular and non-circular kept apart ---------------
    sub("Per cell type (core only, primary centring)")
    for p in LYMPHOCYTES + [IDO1_POS, IDO1_NEG, "Neutrophils"]:
        d = core.loc[core["pheno"] == p]
        circ = p in DETECTION_POOL
        r = fit_mixed(capped(d), "radial_centered_core", f"radial: {p}",
                      note="CIRCULAR: defines the foci" if circ else "")
        if r:
            r["family"] = ("per_phenotype_circular" if circ
                           else "per_phenotype_lymphocyte")
            r["phenotype"] = p
            r["circular"] = circ
            model_rows.append(r)
            print(f"    {p:<26} coef={r['coef_D1MT_vs_ref']:>+8.4f}  "
                  f"p={r['p_value']:.4f}  cells={r['n_cells']:>7,}"
                  f"{'  [circular]' if circ else ''}")
    print("\n    Circular and non-circular populations are now in SEPARATE BH")
    print("    families. Mixing them inflated m and blended a descriptive")
    print("    family with a confirmatory one.")

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

print("    These are refits of ONE hypothesis on subsets, not independent")
print("    tests, so they receive no q-values. Read them as a spread.\n")

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
    print("    measures and has the lowest untreated burden. If the result")
    print("    depends on one animal, that must be visible. Note that dropping")
    print("    an animal takes the design to 2 versus 3, so the p-value moves")
    print("    for reasons of degrees of freedom alone. Read the coefficients.\n")
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
    print("    CD4 is the marker defining helper T cells. Script 03 table 22")
    print("    puts it at ICC 0.403 on p99 and 0.874 on section medians across")
    print("    the two scans, and script 01b showed CD4+ running at 4 to 15")
    print("    percent of T cells on scan_01. The effect should survive without")
    print("    helper T cells.\n")
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
        if fam in NO_BH_FAMILIES:
            continue
        models.loc[g.index, "q_value"] = benjamini_hochberg(g["p_value"].to_numpy())
    models["gets_q_value"] = ~models["family"].isin(NO_BH_FAMILIES)
    write_csv(models, "72_models_all.csv")

    sub("BH families")
    for fam, g in models.groupby("family"):
        got = "no q (refit of one hypothesis)" if fam in NO_BH_FAMILIES else \
            f"{int((g['q_value'] < BH_ALPHA).sum())} of {len(g)} at q < {BH_ALPHA}"
        print(f"    {fam:<32} n={len(g):>3}   {got}")


# %% Cell 7 - radial profiles and figures
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
    g = models.loc[models["family"].isin(
        ["primary", "lineage", "per_phenotype_lymphocyte",
         "per_phenotype_circular"])]
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
    ax.set_title("Leave one animal out\n(dashed = full model, no q-values)",
                 fontsize=FONT_SIZE_TITLE - 14)
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

# ---- F55 the centring check -------------------------------------------------
if len(models) and (models["family"] == "centring_check").any():
    fig, axes = plt.subplots(1, 2, figsize=(30, 13),
                             gridspec_kw={"width_ratios": [1.0, 1.2]})
    ax = axes[0]
    g = models.loc[models["family"] == "centring_check"]
    yy5 = np.arange(len(g))
    for i, (_, r) in enumerate(g.iterrows()):
        col = FLAG_COLOR if r["p_value"] < 0.05 else "#777777"
        ax.plot([r["ci_low"], r["ci_high"]], [i, i], color=col, linewidth=6)
        ax.scatter([r["coef_D1MT_vs_ref"]], [i], s=420, color=col,
                   edgecolor="#FFFFFF", linewidth=2, zorder=3)
        ax.text(r["ci_high"], i, f"  p={r['p_value']:.3f}", va="center",
                fontsize=FONT_SIZE_ANNOT - 10)
    ax.axvline(0, color="#000000", linewidth=3.5)
    ax.set_yticks(yy5)
    ax.set_yticklabels([r.get("centring", "") for _, r in g.iterrows()],
                       fontsize=FONT_SIZE_TICK - 10)
    ax.invert_yaxis()
    ax.set_xlabel("D1MT effect (95% CI)")
    ax.set_title("Does the reference point matter?",
                 fontsize=FONT_SIZE_TITLE - 12)
    ax.xaxis.grid(True, color=GRID_COLOR, linewidth=1.5); ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    ax = axes[1]
    for s in SAMPLE_ORDER:
        d = refs.loc[refs["sample_id"] == s]
        if not len(d):
            continue
        ax.scatter(d["cell_weighted_reference"], d["balanced_reference"],
                   s=280, color=COLOR_OF[s], marker=MARKER_OF[s],
                   edgecolor="#FFFFFF", linewidth=2, zorder=3)
    lo = float(np.nanmin(refs[["cell_weighted_reference", "balanced_reference"]].to_numpy()))
    hi = float(np.nanmax(refs[["cell_weighted_reference", "balanced_reference"]].to_numpy()))
    ax.plot([lo, hi], [lo, hi], color="#000000", linestyle="--", linewidth=3)
    ax.set_xlabel("cell-weighted reference (primary)")
    ax.set_ylabel("phenotype-balanced reference")
    ax.set_title("Where the two references sit, per structure",
                 fontsize=FONT_SIZE_TITLE - 14)
    style_axes(ax)
    ax.legend(handles=[Line2D([0], [0], color=COLOR_OF[s], marker=MARKER_OF[s],
                              markersize=16, linestyle="none",
                              label=f"{short_label(s)} ({COND_OF[s]})")
                       for s in SAMPLE_ORDER],
              frameon=False, fontsize=FONT_SIZE_LEGEND - 16, loc="best")
    fig.suptitle("Composition-free centring check\n"
                 "The primary reference is cell-weighted and moves with "
                 "composition; the balanced one cannot", y=1.04,
                 fontsize=FONT_SIZE_TITLE - 8)
    save_fig(fig, "F55_centring_check")


# %% Cell 8 - delta distance, normalised family dropped
# =============================================================================

banner("DELTA DISTANCE (normalised family dropped)")

print("    Script 07 revision 2 measured why. Correlation of each outcome with")
print("    focus radius: raw median |r| = 0.10, delta 0.11, normalised by")
print("    equivalent radius 0.36, normalised by inscribed radius 0.30. All")
print("    twenty normalised correlations are negative, which is systematic")
print("    over-correction rather than residual confounding.")
print("    The null comes from script 06 table 53, so it is identical in")
print("    scripts 06, 07 and 08.\n")

dist_models = pd.DataFrame()
if RUN_DISTANCE and os.path.exists(NN_NULL_TABLE):
    nulls = pd.read_csv(NN_NULL_TABLE)
    key = ["sample_id", "structure_id", "anchor", "target"]
    if all(c in nulls.columns for c in key + ["null_median_um"]):
        null_lookup = {(r["sample_id"], int(r["structure_id"]),
                        r["anchor"], r["target"]): float(r["null_median_um"])
                       for _, r in nulls.iterrows()}
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
                    mu = null_lookup.get((s, int(k), a, t), np.nan)
                    if not np.isfinite(mu):
                        continue
                    d, _ = cKDTree(xy[ti]).query(xy[ai], k=1)
                    rows.append(pd.DataFrame({
                        "sample_id": s, "condition": COND_OF[s], STRUCT_COL: k,
                        "anchor": a, "target": t,
                        "distance_delta_um": np.asarray(d, float) - mu,
                    }))
        if rows:
            dd = pd.concat(rows, ignore_index=True)
            out = []
            for (a, t), g in dd.groupby(["anchor", "target"]):
                r = fit_mixed(g, "distance_delta_um", f"{a} -> {t}",
                              note="delta outcome; normalised family dropped")
                if r:
                    r["anchor"], r["target"] = a, t
                    r["family"] = "delta_distance_pair"
                    out.append(r)
            for a in NN_ANCHORS:
                g = dd.loc[dd["anchor"] == a]
                r = fit_mixed(g, "distance_delta_um",
                              f"POOLED: {a} -> all lymphocytes",
                              covariates=["target"])
                if r:
                    r["anchor"], r["target"] = a, "pooled"
                    r["family"] = "delta_distance_pooled"
                    out.append(r)
            dist_models = pd.DataFrame(out)
            if len(dist_models):
                dist_models["q_value"] = np.nan
                for fam, g in dist_models.groupby("family"):
                    dist_models.loc[g.index, "q_value"] = benjamini_hochberg(
                        g["p_value"].to_numpy())
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
        n_thr = len(COEXPRESSION_PERCENTILES)
        test_rows = summ.loc[summ["kind"] == "test"]
        ctrl_rows = summ.loc[summ["kind"] == "control"]
        test_n = int(test_rows["n_thresholds_separated"].max()) if len(test_rows) else 0
        ctrl_max = int(ctrl_rows["n_thresholds_separated"].max()) if len(ctrl_rows) else 0
        ctrl_total = int(ctrl_rows["n_thresholds_separated"].sum())
        n_ctrl_pairs = int(len(ctrl_rows))
        expected_chance = CHANCE_SEPARATION_RATE * n_ctrl_pairs * n_thr
        print(f"    iNOS x Arginase-1 separates at {test_n} of {n_thr} thresholds")
        print(f"    worst control pair separates at {ctrl_max} of {n_thr}")
        print(f"    all controls together: {ctrl_total} of "
              f"{n_ctrl_pairs * n_thr} pair-thresholds")
        print(f"    expected by chance    : {expected_chance:.1f}, because with")
        print(f"    three animals per arm complete separation arises with")
        print(f"    probability {CHANCE_SEPARATION_RATE:.2f} per test even under a")
        print(f"    clean panel. A single control separation is therefore NOT")
        print(f"    evidence of spillover on its own.\n")
        if ctrl_max >= test_n and test_n > 0:
            print("    VERDICT: at least one control pair separates the arms as")
            print("    well as or better than the test pair. A macrophage cannot")
            print("    be both a T and a B cell, so the test pair's separation is")
            print("    not established as biology. DROP IT, or report it only")
            print("    with this control alongside.")
        elif ctrl_total > expected_chance * 2 and test_n > 0:
            print("    VERDICT: controls separate well above chance. Treat the")
            print("    test pair as unresolved.")
        elif test_n > ctrl_max:
            print("    VERDICT: the test pair separates more than any control")
            print("    pair does. That is consistent with a real effect, but")
            print("    control behaviour should still be reported alongside so a")
            print("    reader can judge the margin.")
        else:
            print("    VERDICT: neither the test pair nor the controls separate.")
            print("    Nothing to report either way.")

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
    cc = models.loc[models["family"] == "centring_check"]
    if len(cc):
        print("\nCENTRING CHECK")
        for _, r in cc.iterrows():
            print(f"  {r.get('centring', ''):<34} coef {r['coef_D1MT_vs_ref']:+.4f}  "
                  f"p = {r['p_value']:.4f}")
        spread = float(cc["coef_D1MT_vs_ref"].max() - cc["coef_D1MT_vs_ref"].min())
        print(f"  coefficient spread across centrings: {spread:.4f}")
        print("  If that spread is small relative to the coefficient, the result")
        print("  is about position rather than about composition.")

sub("How to report this")
print("  - The outcome is position, not intensity, so the slide confound does")
print("    not apply to the measurement. It DOES apply to cell identity, which")
print("    is why the lineage split and the without-Helper-T check matter.")
print("  - Lymphocytes do not define the foci, so unlike the macrophage and")
print("    neutrophil results these positions are not circular.")
print("  - With three animals per arm the p-value is bounded regardless of cell")
print("    count. The strength of the claim comes from the effect size, the")
print("    per-animal consistency, and the sensitivity analyses, not from p.")
print("  - Convert radial units to microns with the INSCRIBED radius. The")
print("    equivalent radius overstates by the shape ratio, about 1.13 treated")
print("    and 1.22 untreated on this run.")
print("  - Leave-one-out p-values move because the design goes to 2 versus 3.")
print("    Quote the coefficient spread, never a single leave-one-out p.")

sub("Still open")
print("  1. 43111 has one focus against an expert count of two. Its second")
print("     structure has no myeloid density peak above the prominence gate.")
print("  2. The exact area-weighted centring needs one extra column from")
print("     script 04. The balanced centring is the stand-in until then.")
print("  3. BALT organisation (Q7) is untouched. Table 36 has the candidates and")
print("     table 36b the gate sensitivity across three definitions.")

banner("DONE")
sys.stdout = _tee.terminal
_tee.close()
