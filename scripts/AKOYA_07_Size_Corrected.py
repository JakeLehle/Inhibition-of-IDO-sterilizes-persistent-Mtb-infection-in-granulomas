#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - SIZE-CORRECTED DISTANCE MODELS AND POLARISATION
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 07 of the AKOYA analysis series. REVISION 2.

WHAT THIS FIXES FROM SCRIPT 06
    Script 06's nearest-neighbour mixed models use RAW DISTANCE IN MICRONS as
    the outcome. Distances inside a focus are bounded by that focus's size, and
    on the rev4 structures treated foci have a median equivalent radius of
    138 um against 201 um untreated (125 against 168 on the inscribed radius).
    A small focus plus its 150 um cuff cannot geometrically produce a 250 um
    nearest-neighbour distance. Coefficients such as -123 um (IDO1+ to Tregs)
    and -66 um (IDO1- to Tregs) therefore partly restate that treated foci are
    smaller. That is the exact confound the radial coordinate was built to
    avoid, reintroduced through the model outcome. Those rows of table 57 must
    not be used as a treatment claim.

    The per-animal delta values in script 06 table 54 are valid as of revision 2
    of that script, which rebuilt the pooled null. They were not valid before.

    NOTE ON DIRECTION: treated animals have far fewer Tregs inside structures
    (9, 26, 54) than untreated (59, 23, 455), and scarcity should push treated
    distances UP. They went down instead, which is the opposite of the expected
    bias. Some real effect may sit underneath the confound, which is why this
    script re-tests it several independent ways rather than discarding it.

WHAT CHANGED IN REVISION 2 (and why)

    1. ONE NULL, SHARED ACROSS SCRIPTS 06, 07 AND 08.
       Revision 1 recomputed its own per-structure permutation nulls at 200
       permutations while script 06 computed the same quantity at 300, and
       script 08 read script 06's table 53. Three scripts, two different nulls,
       and the delta outcome differed slightly between them for no reason. This
       script now READS null_median_um from script 06 table 53 and only falls
       back to recomputing if that table is absent. The delta outcome is
       therefore identical in 06, 07 and 08.

    2. TWO NORMALISATIONS, NOT ONE.
       Distance was divided by equiv_radius_um. Script 04 revision 4 also
       exports max_inscribed_radius_um, which is the scale the radial coordinate
       actually normalises by, and shape_ratio between them. Both normalisations
       are now computed and both correlations with radius are reported, so the
       choice of denominator is a measured decision rather than an assumption.

    3. THE CO-EXPRESSION RESULT IS NO LONGER PRESENTED AS THE CLEANEST FINDING.
       Revision 1's header argued that because the ratio is a within-section
       correlation between two markers with fixed marginals, global intensity
       differences between scans largely cancel, making it more robust to the
       slide confound than most of the panel. That argument is about BATCH and
       it is still correct as far as it goes. It is not the argument that
       matters. Script 08's spillover control showed that apparent co-expression
       of ANY two markers rises with cell density, and that CD3e with CD20, a
       pair that cannot co-occur in one macrophage, separated the arms at 4 of 6
       thresholds against 3 of 6 for iNOS with Arginase-1. Untreated foci are
       three to five fold denser. The separation is therefore not established as
       biology, and this script now says so at the point of reporting rather
       than leaving it to script 08 to retract.

    4. RAW RADIAL MODELS ARE FLAGGED SUPERSEDED, as in script 06 revision 2.
       Script 08 established that the correct outcome is radial position centred
       on each structure's own mean. These rows are kept and carry
       superseded_by_script_08 = True.

    5. BH IS NOT APPLIED TO THE DIAGNOSTIC FAMILY.
       distance_raw_diagnostic exists to show that an effect appearing only in
       raw microns is size. Giving it q-values invites it to be read as a test.
       Its q-values are now NaN and it is excluded from the significance counts.

    6. STALE NUMBERS AND THE PREREQUISITE CHECK ARE CURRENT.

WHAT THIS SCRIPT DOES
    A  Size confound diagnostic. Observed nearest-neighbour distance against
       focus radius, per structure, on both radius definitions.
    B  Distance outcomes, each modelled separately:
         1. NORMALISED by equivalent radius
         2. NORMALISED by maximum inscribed radius
         3. DELTA, distance minus that focus's own permutation null median
            (read from script 06 table 53)
         4. RAW + COVARIATE, diagnostic only, radius nearly collinear with arm
       An effect that survives the normalised and delta families is real; one
       that appears only in raw distance is size.
    C  Benjamini-Hochberg within declared families, diagnostic family excluded.
    D  Polarisation at the 60th percentile, which widens the exclusively
       polarised classes and recovers mixing rows for the thin treated foci.
    E  The iNOS / Arginase-1 co-expression sweep, reported WITH the spillover
       caveat attached.

CIRCULARITY FLAG
    Foci are defined by pooled myeloid density, so myeloid populations sit
    core-ward by construction. Radial results for CD68+IDO1+, CD68+IDO1-,
    CD163+ and neutrophils are flagged accordingly. What is NOT circular is the
    CONTRAST between IDO1+ and IDO1- macrophages, since both are in the same
    detection pool.

OUTPUTS
    figures/  F47 .. F50
    tables/   60 .. 66

USAGE
    conda activate sc_pre
    python AKOYA_07_Size_Corrected.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================

IN_DIR = "/master/jlehle/WORKING/AKOYA/structures_rev4"
CELL_DIR = f"{IN_DIR}/cell_assignments"
FOCI_TABLE = f"{IN_DIR}/tables/35_foci_structures_relative.csv"
NN_NULL_TABLE = "/master/jlehle/WORKING/AKOYA/distance_stats/tables/53_nn_per_structure.csv"
OUT_DIR = "/master/jlehle/WORKING/AKOYA/size_corrected"

CONDITION_ORDER = ["D1MT", "Untreated"]
CONDITION_COLORS = {"D1MT": "#2C7FB8", "Untreated": "#D95F02"}
REFERENCE_ARM = "Untreated"

STRUCTURE_SOURCE = "foci"
STRUCT_COL = "focus_id" if STRUCTURE_SOURCE == "foci" else "burden_region_id"

# Prefer the per-structure null already computed by script 06, so 06, 07 and 08
# all use one null. Set False to force recomputation here.
USE_SHARED_NULL = True

# ---- phenotypes -------------------------------------------------------------
IDO1_POS = "CD68+IDO1+ Macrophages"
IDO1_NEG = "CD68+IDO1- Macrophages"
CD68_LINEAGE = [IDO1_POS, IDO1_NEG]
MACROPHAGE_ALL = [IDO1_POS, IDO1_NEG, "CD163+ Macrophages"]

# populations used to DEFINE foci; radial results for these are circular
DETECTION_POOL = [IDO1_POS, IDO1_NEG, "CD163+ Macrophages", "Neutrophils"]

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
KEY_PHENOTYPES = [IDO1_POS, IDO1_NEG, "Helper T cells", "CD4- T cells",
                  "Tregs", "B cells", "Plasma cells", "Neutrophils"]

NN_ANCHORS = [IDO1_POS, IDO1_NEG]
NN_TARGETS = ["Helper T cells", "CD4- T cells", "Tregs", "B cells",
              "Plasma cells"]

# ---- permutation (fallback only, when the shared null is unavailable) -------
N_PERMUTATIONS = 300
MAX_ANCHORS_PER_STRUCTURE = 2000
RANDOM_SEED = 0

# ---- tiers (reported, never used to filter) --------------------------------
TIER_SOLID, TIER_USABLE, TIER_PROVISIONAL = 200, 50, 20

# ---- polarisation -----------------------------------------------------------
POLARISATION_PERCENTILES = [50, 60, 70, 75, 80, 90]
POLARISATION_PRIMARY = 60          # lowered from 75 to widen exclusive classes
MIXING_K = 10
MIXING_MIN_PER_CLASS = 5
MIXING_MIN_MACROPHAGES = 20

# ---- mixed models -----------------------------------------------------------
RUN_MIXED_MODELS = True
MODEL_MAX_CELLS_PER_STRUCTURE = 2000
MODEL_MIN_CELLS_PER_ANIMAL = 20
BH_ALPHA = 0.10
# families that are confirmatory and therefore get q-values
CONFIRMATORY_FAMILIES = ["radial_position", "distance_normalised_equiv",
                         "distance_normalised_inscribed", "distance_delta"]
DIAGNOSTIC_FAMILIES = ["distance_raw_diagnostic"]

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
    return float(min((n_ext + 1) / (len(null) + 1), 1.0)), z, mu


def benjamini_hochberg(pvals):
    """Return BH q-values, NaN-safe, preserving input order."""
    p = np.asarray(pvals, float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return q
    idx = np.flatnonzero(ok)
    pv = p[idx]
    order = np.argsort(pv)
    ranked = pv[order]
    m = len(ranked)
    adj = ranked * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(adj, 1.0)
    q[idx] = out
    return q


def nn_distances(anchor_xy, target_xy):
    if len(anchor_xy) == 0 or len(target_xy) == 0:
        return np.array([])
    d, _ = cKDTree(target_xy).query(anchor_xy, k=1)
    return np.asarray(d, float)


def safe_corr(a, b, min_n=4):
    a = np.asarray(a, float); b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < min_n or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return np.nan
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def fit_mixed(df, outcome, label, covariate=None, note=""):
    """
    outcome ~ arm [+ covariate], random intercept for animal,
    structure nested within animal. Falls back to animal-only on failure.
    """
    if not HAVE_SM:
        return None
    need = [outcome, "condition", "sample_id", STRUCT_COL]
    if covariate:
        need.append(covariate)
    missing = [c for c in need if c not in df.columns]
    if missing:
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
    if covariate:
        cv = pd.to_numeric(d[covariate], errors="coerce")
        sd = cv.std()
        d["_cov"] = (cv - cv.mean()) / sd if sd and np.isfinite(sd) and sd > 0 else 0.0
    d = d.dropna(subset=["_y"])
    if not len(d):
        return None

    formula = "_y ~ arm" + (" + _cov" if covariate else "")
    res, mode = None, ""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            md = smf.mixedlm(formula, data=d, groups=d["sample_id"],
                             re_formula="1",
                             vc_formula={"struct": "0 + C(struct_key)"})
            res = md.fit(reml=True, method="lbfgs", maxiter=200)
            mode = "animal + structure"
        except Exception:
            res = None
        if res is None or not np.isfinite(res.params.get("arm", np.nan)):
            try:
                md = smf.mixedlm(formula, data=d, groups=d["sample_id"],
                                 re_formula="1")
                res = md.fit(reml=True, method="lbfgs", maxiter=200)
                mode = "animal only"
            except Exception:
                return None
    if res is None:
        return None

    means = d.groupby("condition")["_y"].mean()
    return {
        "analysis": label, "outcome": outcome, "covariate": covariate or "",
        "random_effects": mode,
        "n_cells": int(len(d)), "n_animals": int(d["sample_id"].nunique()),
        "n_animals_D1MT": int(d.loc[d["condition"] == "D1MT", "sample_id"].nunique()),
        "n_animals_ref": int(d.loc[d["condition"] == REFERENCE_ARM, "sample_id"].nunique()),
        "n_structures": int(d["struct_key"].nunique()),
        "mean_D1MT": float(means.get("D1MT", np.nan)),
        "mean_Untreated": float(means.get("Untreated", np.nan)),
        "coef_D1MT_vs_ref": float(res.params.get("arm", np.nan)),
        "std_err": float(res.bse.get("arm", np.nan)),
        "p_value": float(res.pvalues.get("arm", np.nan)),
        "coef_covariate": float(res.params.get("_cov", np.nan)) if covariate else np.nan,
        "p_covariate": float(res.pvalues.get("_cov", np.nan)) if covariate else np.nan,
        "note": note,
    }


_tee = Tee(os.path.join(TAB_DIR, "00_size_corrected_report.txt"))
sys.stdout = _tee

banner("AKOYA SIZE-CORRECTED DISTANCE MODELS AND POLARISATION (revision 2)")
print(f"Run time      : {datetime.now().isoformat(timespec='seconds')}")
print(f"Input         : {IN_DIR}")
print(f"Polarisation  : {POLARISATION_PRIMARY}th percentile (was 75)")
print(f"Output        : {OUT_DIR}")
if not HAVE_SCIPY:
    print(f"\nERROR: scipy required ({_scipy_err}).")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)
if not HAVE_SM:
    print(f"\n    WARNING: statsmodels unavailable ({_sm_err}). Models skipped.")

print("\nOUTCOME VALIDITY")
print("    RAW distance in microns is bounded by focus size and must NOT be")
print("    used as a treatment claim. Valid outcomes here: distance divided by")
print("    focus radius (two definitions), and distance minus that focus's own")
print("    permutation null. Raw distance with radius as a covariate is a")
print("    diagnostic only, because radius is nearly collinear with arm, and it")
print("    receives no q-values.")


# %% Cell 3 - load
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

# ---- focus geometry, both radius definitions --------------------------------
radius_equiv, radius_inscribed, shape_ratio = {}, {}, {}
if os.path.exists(FOCI_TABLE):
    ft = pd.read_csv(FOCI_TABLE)
    idcol = "focus_id" if "focus_id" in ft.columns else STRUCT_COL
    have_inscribed = "max_inscribed_radius_um" in ft.columns
    for _, r in ft.iterrows():
        key = (r["sample_id"], int(r[idcol]))
        radius_equiv[key] = float(r["equiv_radius_um"])
        if have_inscribed:
            radius_inscribed[key] = float(r["max_inscribed_radius_um"])
            shape_ratio[key] = float(r.get("shape_ratio", np.nan))
    print(f"\n    loaded geometry for {len(radius_equiv)} structures"
          f"{'' if have_inscribed else '  (no inscribed radius: rerun script 04)'}")
    if have_inscribed:
        for c in CONDITION_ORDER:
            g = ft.loc[ft["condition"] == c]
            if len(g):
                print(f"      {c:<12} equiv {g['equiv_radius_um'].median():>6.0f} um  "
                      f"inscribed {g['max_inscribed_radius_um'].median():>6.0f} um  "
                      f"shape {g['shape_ratio'].median():>4.2f}")
else:
    print(f"\n    WARNING: {FOCI_TABLE} not found. Radii will be derived from")
    print("    cell coordinates, which is a rougher estimate.")
    have_inscribed = False

# ---- shared per-structure null from script 06 -------------------------------
SHARED_NULL = {}
if USE_SHARED_NULL and os.path.exists(NN_NULL_TABLE):
    nt = pd.read_csv(NN_NULL_TABLE)
    need = ["sample_id", "structure_id", "anchor", "target", "null_median_um"]
    if all(c in nt.columns for c in need):
        for _, r in nt.iterrows():
            SHARED_NULL[(r["sample_id"], int(r["structure_id"]),
                         r["anchor"], r["target"])] = float(r["null_median_um"])
        print(f"    loaded {len(SHARED_NULL)} shared per-structure nulls from")
        print(f"    script 06 table 53. The delta outcome is therefore identical")
        print(f"    in scripts 06, 07 and 08.")
    else:
        print(f"    WARNING: {NN_NULL_TABLE} lacks expected columns, will "
              f"recompute nulls here.")
elif USE_SHARED_NULL:
    print(f"    WARNING: {NN_NULL_TABLE} not found. Run script 06 first, or "
          f"nulls will be recomputed here and will not match 06 and 08.")

sub("STRUCTURE COUNTS PER ANIMAL")
for c in CONDITION_ORDER:
    mem = [s for s in SAMPLE_ORDER if COND_OF[s] == c]
    print(f"    {c:<12} "
          f"{ {short_label(s): int(cells[s].loc[cells[s][STRUCT_COL] > 0, STRUCT_COL].nunique()) for s in mem} }")


def structure_frames():
    for s in SAMPLE_ORDER:
        d = cells[s]
        for k, g in d.loc[d[STRUCT_COL] > 0].groupby(STRUCT_COL):
            yield s, int(k), g


def focus_radii(sid, k, g):
    """Return (equivalent radius, inscribed radius). NaN where unavailable."""
    key = (sid, int(k))
    r_eq = radius_equiv.get(key, np.nan)
    r_in = radius_inscribed.get(key, np.nan)
    if not (np.isfinite(r_eq) and r_eq > 0):
        core = g.loc[g["region"] == "core"]
        if len(core) >= 3:
            x, y = core["x"].to_numpy(), core["y"].to_numpy()
            r_eq = float(0.5 * np.hypot(x.max() - x.min(), y.max() - y.min()))
        else:
            r_eq = np.nan
    return r_eq, r_in


# %% Cell 4 - A: the size confound, made visible
# =============================================================================

banner("A - SIZE CONFOUND DIAGNOSTIC")

print("    Observed nearest-neighbour distance against focus radius. If distance")
print("    tracks radius, raw distance cannot be compared between arms whose")
print("    foci differ in size.\n")

percell_frames, struct_rows = [], []
n_shared, n_recomputed = 0, 0

for s, k, g in structure_frames():
    xy = g[["x", "y"]].to_numpy(float)
    ph = g["pheno"].to_numpy()
    n = len(xy)
    if n < 10:
        continue
    r_eq, r_in = focus_radii(s, k, g)
    for anchor in NN_ANCHORS:
        a_idx = np.flatnonzero(ph == anchor)
        n_a = len(a_idx)
        if n_a == 0:
            continue
        a_use = (rng.choice(a_idx, size=MAX_ANCHORS_PER_STRUCTURE, replace=False)
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

            # the null: shared from script 06 when available
            mu = SHARED_NULL.get((s, int(k), anchor, target), np.nan)
            if np.isfinite(mu):
                n_shared += 1
                pval, z = np.nan, np.nan
            else:
                n_recomputed += 1
                null = np.empty(N_PERMUTATIONS)
                for it in range(N_PERMUTATIONS):
                    perm = rng.permutation(n)
                    fa, ft_ = perm[:n_a], perm[n_a:n_a + n_t]
                    fa_use = (rng.choice(fa, size=MAX_ANCHORS_PER_STRUCTURE,
                                         replace=False)
                              if len(fa) > MAX_ANCHORS_PER_STRUCTURE else fa)
                    v = nn_distances(xy[fa_use], xy[ft_])
                    null[it] = np.median(v) if len(v) else np.nan
                pval, z, mu = empirical_p(obs, null)

            struct_rows.append({
                "sample_id": s, "condition": COND_OF[s], "structure_id": k,
                "anchor": anchor, "target": target,
                "n_anchor": n_a, "n_target": n_t, "n_structure": n,
                "focus_radius_um": r_eq,
                "focus_inscribed_radius_um": r_in,
                "shape_ratio": shape_ratio.get((s, int(k)), np.nan),
                "observed_median_um": obs, "null_median_um": mu,
                "null_source": "script_06" if np.isfinite(
                    SHARED_NULL.get((s, int(k), anchor, target), np.nan)) else "recomputed",
                "delta_um": obs - mu if np.isfinite(mu) else np.nan,
                "normalised_median_equiv": obs / r_eq if (r_eq and np.isfinite(r_eq)) else np.nan,
                "normalised_median_inscribed": obs / r_in if (r_in and np.isfinite(r_in)) else np.nan,
                "z": z, "p_empirical": pval,
                "tier": tier_of(min(n_a, n_t)),
            })

            take = obs_v
            if len(take) > MODEL_MAX_CELLS_PER_STRUCTURE:
                take = rng.choice(take, size=MODEL_MAX_CELLS_PER_STRUCTURE,
                                  replace=False)
            if len(take):
                percell_frames.append(pd.DataFrame({
                    "sample_id": s, "condition": COND_OF[s], STRUCT_COL: k,
                    "anchor": anchor, "target": target,
                    "distance_um": take,
                    "focus_radius_um": r_eq,
                    "focus_inscribed_radius_um": r_in,
                    "distance_norm_equiv": take / r_eq if (r_eq and np.isfinite(r_eq)) else np.nan,
                    "distance_norm_inscribed": take / r_in if (r_in and np.isfinite(r_in)) else np.nan,
                    "distance_delta_um": take - mu if np.isfinite(mu) else np.nan,
                }))
    gc.collect()

print(f"    nulls: {n_shared} taken from script 06, {n_recomputed} recomputed here")

nn_struct = pd.DataFrame(struct_rows)
percell = (pd.concat(percell_frames, ignore_index=True)
           if percell_frames else pd.DataFrame())
corrs = pd.DataFrame()
if len(nn_struct):
    write_csv(nn_struct, "60_nn_per_structure_with_geometry.csv")

    sub("Correlation of each distance outcome with focus radius, per pair")
    print(f"    {'pair':<42}{'raw':>8}{'norm eq':>9}{'norm in':>9}{'delta':>8}")
    print("    " + "-" * 76)
    corr_rows = []
    for (a, t), g in nn_struct.groupby(["anchor", "target"]):
        rec = {"anchor": a, "target": t, "n_structures": int(len(g))}
        rec["corr_raw_vs_radius"] = safe_corr(g["observed_median_um"],
                                              g["focus_radius_um"])
        rec["corr_norm_equiv_vs_radius"] = safe_corr(g["normalised_median_equiv"],
                                                     g["focus_radius_um"])
        rec["corr_norm_inscribed_vs_radius"] = safe_corr(
            g["normalised_median_inscribed"], g["focus_radius_um"])
        rec["corr_delta_vs_radius"] = safe_corr(g["delta_um"],
                                                g["focus_radius_um"])
        corr_rows.append(rec)
        print(f"    {(a[:20] + ' -> ' + t)[:41]:<42}"
              f"{rec['corr_raw_vs_radius']:>+8.2f}"
              f"{rec['corr_norm_equiv_vs_radius']:>+9.2f}"
              f"{rec['corr_norm_inscribed_vs_radius']:>+9.2f}"
              f"{rec['corr_delta_vs_radius']:>+8.2f}")
    corrs = pd.DataFrame(corr_rows)
    write_csv(corrs, "61_distance_radius_correlation.csv")

    sub("Which correction actually removes size?")
    for col, name in [("corr_raw_vs_radius", "raw microns"),
                      ("corr_norm_equiv_vs_radius", "normalised, equivalent r"),
                      ("corr_norm_inscribed_vs_radius", "normalised, inscribed r"),
                      ("corr_delta_vs_radius", "delta vs focus null")]:
        v = corrs[col].dropna()
        if not len(v):
            continue
        print(f"    {name:<28} median |r| = {v.abs().median():.2f}   "
              f"range {v.min():+.2f} to {v.max():+.2f}")
    print("\n    The outcome to trust is the one whose correlation with radius")
    print("    is nearest zero. Over-correction shows as a correlation that")
    print("    flips sign and grows, which is worse than the raw confound.")

    # ---- F47 the confound --------------------------------------------------
    panels = [("observed_median_um", "Raw distance (INVALID outcome)"),
              ("normalised_median_equiv", "Distance / equivalent radius"),
              ("normalised_median_inscribed", "Distance / inscribed radius"),
              ("delta_um", "Distance minus focus null")]
    fig, axes = plt.subplots(1, len(panels), figsize=(13 * len(panels), 13))
    for ax, (col, ttl) in zip(axes, panels):
        for s in SAMPLE_ORDER:
            d = nn_struct.loc[nn_struct["sample_id"] == s]
            if not len(d):
                continue
            ax.scatter(d["focus_radius_um"], d[col], s=200, color=COLOR_OF[s],
                       marker=MARKER_OF[s], edgecolor="#FFFFFF", linewidth=1.5,
                       alpha=0.8, zorder=3)
        ax.set_xscale("log")
        ax.set_xlabel("focus equivalent radius (µm)",
                      fontsize=FONT_SIZE_BASE - 12)
        ax.set_ylabel(ttl, fontsize=FONT_SIZE_BASE - 14)
        ax.set_title(ttl, fontsize=FONT_SIZE_TITLE - 16,
                     color=FLAG_COLOR if "INVALID" in ttl else TEXT_COLOR)
        style_axes(ax)
    axes[0].legend(handles=[Line2D([0], [0], color=COLOR_OF[s],
                                   marker=MARKER_OF[s], markersize=16,
                                   linestyle="none",
                                   label=f"{short_label(s)} ({COND_OF[s]})")
                            for s in SAMPLE_ORDER],
                   frameon=False, fontsize=FONT_SIZE_LEGEND - 16, loc="best")
    fig.suptitle("Why raw distance cannot be compared between arms\n"
                 "Treated foci are smaller, so their distances are bounded "
                 "smaller regardless of biology", y=1.05,
                 fontsize=FONT_SIZE_TITLE - 8)
    save_fig(fig, "F47_size_confound_diagnostic")


# %% Cell 5 - B/C: distance outcomes with BH correction
# =============================================================================

banner("B/C - DISTANCE MODELS ON FOUR OUTCOMES, WITH BH CORRECTION")

model_rows = []
if RUN_MIXED_MODELS and HAVE_SM and len(percell):
    for (a, t), g in percell.groupby(["anchor", "target"]):
        for outcome, fam, cov, note in [
            ("distance_norm_equiv", "distance_normalised_equiv", None,
             "distance / equivalent radius"),
            ("distance_norm_inscribed", "distance_normalised_inscribed", None,
             "distance / maximum inscribed radius"),
            ("distance_delta_um", "distance_delta", None,
             "distance minus focus permutation null (shared with 06 and 08)"),
            ("distance_um", "distance_raw_diagnostic", "focus_radius_um",
             "DIAGNOSTIC ONLY: radius nearly collinear with arm, no q-value"),
        ]:
            if outcome not in g.columns or not np.isfinite(g[outcome]).any():
                continue
            r = fit_mixed(g, outcome, f"{a} -> {t}", covariate=cov, note=note)
            if r:
                r["family"] = fam
                r["anchor"], r["target"] = a, t
                r["superseded_by_script_08"] = False
                model_rows.append(r)

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
            frames += [g.sample(min(len(g), MODEL_MAX_CELLS_PER_STRUCTURE),
                                random_state=RANDOM_SEED)
                       for _, g in sel.groupby(STRUCT_COL)]
        if not frames:
            continue
        dd = pd.concat(frames, ignore_index=True)
        circ = p in DETECTION_POOL
        note = "raw radial_pos, superseded by script 08's centred outcome"
        if circ:
            note += "; CIRCULAR: population defines the foci"
        r = fit_mixed(dd, "radial_pos", f"radial: {p}", note=note)
        if r:
            r["family"] = "radial_position"
            r["phenotype"] = p
            r["circular"] = circ
            r["superseded_by_script_08"] = True
            model_rows.append(r)

models = pd.DataFrame(model_rows)
if len(models):
    models["q_value"] = np.nan
    for fam, g in models.groupby("family"):
        if fam in DIAGNOSTIC_FAMILIES:
            continue
        models.loc[g.index, "q_value"] = benjamini_hochberg(g["p_value"].to_numpy())
    models["sig_p05"] = models["p_value"] < 0.05
    models["sig_q10"] = models["q_value"] < BH_ALPHA
    models["is_confirmatory"] = models["family"].isin(CONFIRMATORY_FAMILIES)
    write_csv(models, "62_mixed_models_all_outcomes.csv")

    for fam in CONFIRMATORY_FAMILIES + DIAGNOSTIC_FAMILIES:
        g = models.loc[models["family"] == fam].sort_values("p_value")
        if not len(g):
            continue
        flag = ("   [DIAGNOSTIC ONLY, no q-values]" if fam in DIAGNOSTIC_FAMILIES
                else "   [SUPERSEDED by script 08]" if fam == "radial_position"
                else "")
        sub(f"Family: {fam}  (n = {len(g)} tests){flag}")
        print(f"    {'test':<44}{'coef':>10}{'p':>9}{'q(BH)':>9}  cells")
        for _, r in g.iterrows():
            mark = " *" if r["sig_q10"] else ""
            circ = " [circular]" if r.get("circular", False) is True else ""
            qtxt = "      na" if pd.isna(r["q_value"]) else f"{r['q_value']:>9.4f}"
            print(f"    {r['analysis'][:43]:<44}{r['coef_D1MT_vs_ref']:>+10.3f}"
                  f"{r['p_value']:>9.4f}{qtxt}  "
                  f"{r['n_cells']:>7,}{mark}{circ}")

    sub("Which distance effects survive size correction")
    print("    An effect must run the same direction and reach p < 0.05 in the")
    print("    delta family AND at least one normalisation to count.\n")
    fams_needed = ["distance_normalised_equiv", "distance_normalised_inscribed",
                   "distance_delta"]
    for (a, t), g in models.loc[models["family"].isin(fams_needed)].groupby(
            ["anchor", "target"]):
        got = {f: g.loc[g["family"] == f] for f in fams_needed}
        if not len(got["distance_delta"]):
            continue
        cd = float(got["distance_delta"]["coef_D1MT_vs_ref"].iloc[0])
        pd_ = float(got["distance_delta"]["p_value"].iloc[0])
        parts, ok_norm = [], False
        for f in fams_needed[:2]:
            if not len(got[f]):
                parts.append("na")
                continue
            cn = float(got[f]["coef_D1MT_vs_ref"].iloc[0])
            pn = float(got[f]["p_value"].iloc[0])
            same = np.sign(cn) == np.sign(cd)
            parts.append(f"p={pn:.3f}{'' if same else ' (opp dir)'}")
            if pn < 0.05 and same:
                ok_norm = True
        verdict = ("SURVIVES" if (pd_ < 0.05 and ok_norm)
                   else "partial" if (pd_ < 0.05 or ok_norm)
                   else "does not survive")
        print(f"    {a[:20]:<22} -> {t:<16} delta p={pd_:.3f}  "
              f"norm-eq {parts[0]:<20} norm-in {parts[1]:<20} -> {verdict}")

    sub("Multiple testing position")
    conf = models.loc[models["is_confirmatory"]]
    n_p05 = int(conf["sig_p05"].sum())
    n_q10 = int(conf["sig_q10"].sum())
    print(f"    Confirmatory families only ({len(conf)} tests):")
    print(f"      {n_p05} have p < 0.05")
    print(f"      {n_q10} survive BH at q < {BH_ALPHA}")
    if n_q10 == 0:
        print("\n    Nothing survives BH. That does not invalidate the work, but")
        print("    it does mean no single test here is a standalone claim. The")
        print("    argument rests on coherence: consistent direction across")
        print("    populations and animals, with the internal control")
        print("    (neutrophils) inverting exactly where necrosis predicts.")

    # ---- F48 forest by family ----------------------------------------------
    fams = [f for f in ["distance_normalised_equiv",
                        "distance_normalised_inscribed", "distance_delta"]
            if (models["family"] == f).any()]
    if fams:
        fig, axes = plt.subplots(1, len(fams), figsize=(16 * len(fams), 15),
                                 squeeze=False)
        for ax, fam in zip(axes[0], fams):
            g = models.loc[models["family"] == fam].sort_values("coef_D1MT_vs_ref")
            yy = np.arange(len(g))
            lo = g["coef_D1MT_vs_ref"] - 1.96 * g["std_err"]
            hi = g["coef_D1MT_vs_ref"] + 1.96 * g["std_err"]
            for i, (_, r) in enumerate(g.iterrows()):
                col = (FLAG_COLOR if r["sig_q10"] else
                       OK_COLOR if r["sig_p05"] else "#777777")
                ax.plot([lo.iloc[i], hi.iloc[i]], [i, i], color=col, linewidth=5,
                        zorder=2)
                ax.scatter([r["coef_D1MT_vs_ref"]], [i], s=380, color=col,
                           edgecolor="#FFFFFF", linewidth=2, zorder=3)
                ax.text(hi.iloc[i], i, f"  p={r['p_value']:.3f}", va="center",
                        fontsize=FONT_SIZE_ANNOT - 12)
            ax.axvline(0, color="#000000", linewidth=3.5)
            ax.set_yticks(yy)
            ax.set_yticklabels([r["analysis"][:38] for _, r in g.iterrows()],
                               fontsize=FONT_SIZE_TICK - 14)
            ax.set_xlabel("D1MT effect vs Untreated (95% CI)",
                          fontsize=FONT_SIZE_BASE - 10)
            ax.set_title(fam.replace("_", " "), fontsize=FONT_SIZE_TITLE - 12)
            ax.xaxis.grid(True, color=GRID_COLOR, linewidth=1.5)
            ax.set_axisbelow(True)
            ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        fig.legend(handles=[Line2D([0], [0], color=FLAG_COLOR, linewidth=6,
                                   label=f"q < {BH_ALPHA}"),
                            Line2D([0], [0], color=OK_COLOR, linewidth=6,
                                   label="p < 0.05 only"),
                            Line2D([0], [0], color="#777777", linewidth=6,
                                   label="not significant")],
                   loc="lower center", ncol=3, frameon=False,
                   fontsize=FONT_SIZE_LEGEND - 8, bbox_to_anchor=(0.5, -0.03))
        fig.suptitle("Distance treatment effects by outcome family, with BH\n"
                     "Raw-distance and raw-radial models are excluded: one is "
                     "bounded by focus size, the other superseded by script 08",
                     y=1.03, fontsize=FONT_SIZE_TITLE - 8)
        save_fig(fig, "F48_model_forest_by_family")
else:
    models = pd.DataFrame()
    print("    SKIPPED (statsmodels unavailable or no per-cell data)")


# %% Cell 6 - D/E: polarisation at the 60th percentile
# =============================================================================

banner("D/E - POLARISATION AND CO-EXPRESSION")

print("    READ THE CAVEAT BEFORE THE NUMBERS. Segmentation spillover inflates")
print("    apparent co-expression of ANY two markers, and worsens as cells pack")
print("    more densely. Untreated foci run three to five fold denser at peak.")
print("    Script 08's control pairs make this concrete: CD3e with CD20, a pair")
print("    that cannot co-occur inside one macrophage, separated the arms at 4")
print("    of 6 thresholds, MORE than iNOS with Arginase-1 did at 3 of 6.")
print("    Whatever separation appears below is therefore NOT established as")
print("    biology. It is reported for completeness and to keep the sweep on")
print("    record, not as a finding.\n")

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
        obs = float((hi_i & hi_a).mean())
        exp = p_i * p_a
        pol_rows.append({
            "sample_id": s, "animal_id": short_label(s), "condition": COND_OF[s],
            "percentile": pct, "n_macrophages": len(mac),
            "pct_double_positive": 100.0 * obs,
            "pct_expected_if_independent": 100.0 * exp,
            "coexpression_ratio": obs / exp if exp > 0 else np.nan,
            "pct_double_negative": 100.0 * float((~hi_i & ~hi_a).mean()),
        })

    it_ = float(np.nanpercentile(mac["iNOS"], POLARISATION_PRIMARY))
    at_ = float(np.nanpercentile(mac["Arginase-1"], POLARISATION_PRIMARY))
    for k, g in d.loc[(d[STRUCT_COL] > 0)
                      & d["pheno"].isin(MACROPHAGE_ALL)].groupby(STRUCT_COL):
        n = len(g)
        if n < MIXING_MIN_MACROPHAGES:
            continue
        hi_i = (g["iNOS"] >= it_).to_numpy()
        hi_a = (g["Arginase-1"] >= at_).to_numpy()
        ex_i, ex_a = hi_i & ~hi_a, hi_a & ~hi_i
        n_i, n_a = int(ex_i.sum()), int(ex_a.sum())
        if n_i < MIXING_MIN_PER_CLASS or n_a < MIXING_MIN_PER_CLASS:
            continue
        xy = g[["x", "y"]].to_numpy(float)
        kk = int(min(MIXING_K, n - 1))
        _, idxs = cKDTree(xy).query(xy, k=kk + 1)
        idxs = idxs[:, 1:]

        def mixing(am, bm, _idxs=idxs):
            rows = np.flatnonzero(am)
            return float(bm[_idxs[rows]].mean()) if len(rows) else np.nan

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
            "percentile": POLARISATION_PRIMARY,
            "n_macrophages": n, "n_inos_only": n_i, "n_arg1_only": n_a,
            "observed_mixing": obs, "null_mixing": mu,
            "delta_mixing": obs - mu if np.isfinite(mu) else np.nan,
            "z": z, "p_empirical": pval, "tier": tier_of(min(n_i, n_a)),
        })
        del idxs

pol = pd.DataFrame(pol_rows)
mix = pd.DataFrame(mix_rows)
sepdf = pd.DataFrame()

if len(pol):
    write_csv(pol, "63_polarisation_sweep.csv")
    sub(f"Co-expression ratio at the {POLARISATION_PRIMARY}th percentile")
    p0 = pol.loc[pol["percentile"] == POLARISATION_PRIMARY]
    print(p0[["animal_id", "condition", "n_macrophages", "pct_double_positive",
              "pct_expected_if_independent", "coexpression_ratio"]]
          .to_string(index=False))

    sub("Separation across the whole sweep")
    print(f"    {'pct':>5}{'D1MT range':>26}{'Untreated range':>26}{'overlap':>10}")
    sep_rows = []
    for pct in POLARISATION_PERCENTILES:
        g = pol.loc[pol["percentile"] == pct]
        a = g.loc[g["condition"] == "D1MT", "coexpression_ratio"].dropna()
        b = g.loc[g["condition"] == "Untreated", "coexpression_ratio"].dropna()
        if not len(a) or not len(b):
            continue
        overlap = not (a.max() < b.min() or b.max() < a.min())
        sep_rows.append({"percentile": pct, "d1mt_min": a.min(),
                         "d1mt_max": a.max(), "untr_min": b.min(),
                         "untr_max": b.max(), "overlap": overlap,
                         "gap": (b.min() - a.max()) if not overlap else np.nan})
        print(f"    {pct:>5}{a.min():>12.3f} - {a.max():<12.3f}"
              f"{b.min():>12.3f} - {b.max():<12.3f}"
              f"{'yes' if overlap else 'NO':>10}")
    sepdf = pd.DataFrame(sep_rows)
    if len(sepdf):
        write_csv(sepdf, "64_coexpression_separation.csv")
        clean = sepdf.loc[~sepdf["overlap"], "percentile"].tolist()
        if clean:
            print(f"\n    Arms do not overlap at percentiles: {clean}")
            print("    DO NOT REPORT THIS AS A FINDING. Script 08's spillover")
            print("    control pairs separate the arms as well or better, and")
            print("    they cannot co-occur in a single macrophage. Until that")
            print("    control is calibrated, this separation is consistent with")
            print("    segmentation spillover tracking the density difference")
            print("    between arms. Revisitable with nuclear-only compartments")
            print("    or a cell-size covariate; separate work.")
        else:
            print("\n    Arms overlap at every threshold.")

    # ---- F49 co-expression sweep -------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(30, 14))
    ax = axes[0]
    for s in SAMPLE_ORDER:
        d = pol.loc[pol["sample_id"] == s].sort_values("percentile")
        if not len(d):
            continue
        ax.plot(d["percentile"], d["coexpression_ratio"], linewidth=5,
                marker=MARKER_OF[s], markersize=18, color=COLOR_OF[s], alpha=0.9)
    ax.axhline(1.0, color="#000000", linestyle="--", linewidth=3)
    ax.text(POLARISATION_PERCENTILES[0], 1.02, "independence",
            fontsize=FONT_SIZE_ANNOT - 10)
    ax.set_xlabel("within-section percentile threshold")
    ax.set_ylabel("iNOS / Arginase-1 co-expression ratio")
    ax.set_title("Co-expression above chance", fontsize=FONT_SIZE_TITLE - 12)
    style_axes(ax)
    ax.legend(handles=[Line2D([0], [0], color=COLOR_OF[s], marker=MARKER_OF[s],
                              markersize=16, linewidth=4,
                              label=f"{short_label(s)} ({COND_OF[s]})")
                       for s in SAMPLE_ORDER],
              frameon=False, fontsize=FONT_SIZE_LEGEND - 16, loc="best")

    ax = axes[1]
    p0 = pol.loc[pol["percentile"] == POLARISATION_PRIMARY]
    for ci, c in enumerate(CONDITION_ORDER):
        v = p0.loc[p0["condition"] == c, "coexpression_ratio"].dropna()
        if not len(v):
            continue
        j = rng.uniform(-0.1, 0.1, size=len(v))
        ax.scatter(np.full(len(v), ci) + j, v, s=520, color=CONDITION_COLORS[c],
                   edgecolor="#FFFFFF", linewidth=3, zorder=3)
        ax.hlines(v.median(), ci - 0.3, ci + 0.3, color="#000000", linewidth=5)
    ax.axhline(1.0, color="#000000", linestyle="--", linewidth=3)
    ax.set_xticks(range(len(CONDITION_ORDER)))
    ax.set_xticklabels(CONDITION_ORDER)
    ax.set_xlim(-0.6, len(CONDITION_ORDER) - 0.4)
    ax.set_ylabel("co-expression ratio")
    ax.set_title(f"At the {POLARISATION_PRIMARY}th percentile",
                 fontsize=FONT_SIZE_TITLE - 12)
    style_axes(ax)
    fig.suptitle("iNOS and Arginase-1 co-expression in macrophages\n"
                 "NOT A FINDING: spillover controls in script 08 separate the "
                 "arms as well or better, and cannot co-occur in one cell",
                 y=1.05, fontsize=FONT_SIZE_TITLE - 10)
    save_fig(fig, "F49_coexpression_ratio")

if len(mix):
    write_csv(mix, "65_inos_arg1_mixing.csv")
    mix_anim = (mix.groupby(["sample_id", "condition"])
                .agg(median_delta=("delta_mixing", "median"),
                     median_observed=("observed_mixing", "median"),
                     n_structures=("z", "size")).reset_index())
    sub(f"Spatial mixing at the {POLARISATION_PRIMARY}th percentile")
    print(mix_anim.to_string(index=False))
    write_csv(mix_anim, "65b_inos_arg1_mixing_by_animal.csv")
    missing = [s for s in SAMPLE_ORDER if s not in set(mix["sample_id"])]
    if missing:
        print(f"\n    Still no mixing rows for: {missing}")
        print("    Lower POLARISATION_PRIMARY further, or accept that those")
        print("    foci hold too few exclusively-polarised macrophages.")
    else:
        print("\n    All animals now contribute mixing rows.")

    fig, ax = plt.subplots(figsize=(18, 13))
    for ci, c in enumerate(CONDITION_ORDER):
        v = mix.loc[mix["condition"] == c, "delta_mixing"].dropna()
        if not len(v):
            continue
        j = rng.uniform(-0.14, 0.14, size=len(v))
        ax.scatter(np.full(len(v), ci) + j, v, s=300, color=CONDITION_COLORS[c],
                   edgecolor="#FFFFFF", linewidth=2, zorder=3, alpha=0.85)
        ax.hlines(v.median(), ci - 0.3, ci + 0.3, color="#000000", linewidth=5)
    ax.axhline(0, color="#000000", linestyle="--", linewidth=3)
    ax.set_xticks(range(len(CONDITION_ORDER)))
    ax.set_xticklabels(CONDITION_ORDER)
    ax.set_ylabel("delta mixing vs null\n(negative = more segregated)")
    ax.set_title("iNOS-high and Arginase-1-high spatial segregation\n"
                 f"each point one structure, {POLARISATION_PRIMARY}th percentile",
                 fontsize=FONT_SIZE_TITLE - 10)
    style_axes(ax)
    save_fig(fig, "F50_polarisation_mixing")


# %% Cell 7 - consistency and wrap up
# =============================================================================

banner("CONSISTENCY ACROSS ANIMALS")

cons_rows = []


def consistency(df, group_cols, value_col, label):
    for keys, g in df.groupby(group_cols):
        keys = keys if isinstance(keys, tuple) else (keys,)
        rec = {"analysis": label, "value": value_col}
        for cn, kv in zip(group_cols, keys):
            rec[cn] = kv
        for c in CONDITION_ORDER:
            v = g.loc[g["condition"] == c, value_col].dropna()
            rec[f"{c}_n"] = len(v)
            rec[f"{c}_median"] = float(v.median()) if len(v) else np.nan
            rec[f"{c}_min"] = float(v.min()) if len(v) else np.nan
            rec[f"{c}_max"] = float(v.max()) if len(v) else np.nan
        a, b = rec.get("D1MT_median", np.nan), rec.get("Untreated_median", np.nan)
        rec["arm_difference"] = a - b if np.isfinite(a) and np.isfinite(b) else np.nan
        amin, amax = rec.get("D1MT_min", np.nan), rec.get("D1MT_max", np.nan)
        bmin, bmax = rec.get("Untreated_min", np.nan), rec.get("Untreated_max", np.nan)
        rec["complete_separation"] = bool(
            np.isfinite([amin, amax, bmin, bmax]).all()
            and (amax < bmin or bmax < amin))
        cons_rows.append(rec)


if len(nn_struct):
    for col, label in [("normalised_median_equiv", "nn_normalised_equiv"),
                       ("normalised_median_inscribed", "nn_normalised_inscribed"),
                       ("delta_um", "nn_delta")]:
        if col not in nn_struct.columns or not np.isfinite(nn_struct[col]).any():
            continue
        per_animal = (nn_struct.groupby(["sample_id", "condition",
                                         "anchor", "target"])[col]
                      .median().reset_index())
        consistency(per_animal, ["anchor", "target"], col, label)
if len(pol):
    p0 = pol.loc[pol["percentile"] == POLARISATION_PRIMARY].copy()
    p0["phenotype"] = "iNOS_Arg1_coexpression_SPILLOVER_UNRESOLVED"
    consistency(p0, ["phenotype"], "coexpression_ratio", "polarisation")

cons = pd.DataFrame(cons_rows)
if len(cons):
    write_csv(cons, "66_consistency_summary.csv")
    sub("Comparisons with COMPLETE SEPARATION between arms")
    sep = cons.loc[cons["complete_separation"]]
    if not len(sep):
        print("    none")
    else:
        for _, r in sep.iterrows():
            name = " / ".join(str(r[c]) for c in ["phenotype", "anchor", "target"]
                              if c in r.index and pd.notna(r.get(c)))
            print(f"    [{r['analysis']:<24}] {name:<50} "
                  f"D1MT [{r['D1MT_min']:.3f}, {r['D1MT_max']:.3f}]   "
                  f"Untr [{r['Untreated_min']:.3f}, {r['Untreated_max']:.3f}]")
        print("\n    With three animals per arm, complete separation is the")
        print("    strongest evidence this design can produce. It corresponds to")
        print("    the minimum p a rank test can return (0.10, two-sided), and")
        print("    it arises by chance with probability 0.10, so a single")
        print("    separated comparison among many is not remarkable on its own.")

banner("SUMMARY")
print(f"Structure pair-rows      : {len(nn_struct)}")
print(f"Nulls shared from 06     : {n_shared}   recomputed here: {n_recomputed}")
print(f"Models fitted            : {len(models)}")
if len(models):
    conf = models.loc[models["is_confirmatory"]]
    print(f"  confirmatory tests     : {len(conf)}")
    print(f"  p < 0.05               : {int(conf['sig_p05'].sum())}")
    print(f"  q < {BH_ALPHA} (BH)          : {int(conf['sig_q10'].sum())}")
print(f"Polarisation rows        : {len(pol)}")
print(f"Mixing rows              : {len(mix)}")

sub("Read in this order")
print("  1. F47 / table 61 : the size confound, and which correction removes")
print("     it. Compare the four correlation columns and pick the outcome")
print("     whose correlation with radius is nearest zero.")
print("  2. F48 / table 62 : models by family with BH q-values. The distance")
print("     conclusion is whatever survives delta AND a normalisation.")
print("  3. F50 / table 65 : polarisation mixing at the 60th percentile.")
print("  4. Table 66       : complete separation between arms.")
print("  5. F49 / table 64 : co-expression sweep. ON RECORD, NOT A FINDING.")
print("     See the spillover caveat printed above it.")

sub("What is claimable after this run")
print("  - Burden: near-zero granuloma-density tissue in treated animals under")
print("    a common definition, and unchanged by the edge correction.")
print("    Strongest and simplest result.")
print("  - Architecture: lymphocytes core-ward in treated foci, rim-ward in")
print("    untreated lesions, with neutrophils inverting as necrosis predicts.")
print("    A pattern argument backed by per-animal separation, not a single p.")
print("  - NOT claimable: anything from raw distances in microns, anything from")
print("    3-hydroxykynurenine, and the iNOS / Arginase-1 co-expression")
print("    separation until the spillover control is calibrated.")

banner("DONE")
sys.stdout = _tee.terminal
_tee.close()
