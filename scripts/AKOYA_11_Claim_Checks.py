#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - CLAIM CHECKS FOR THE FOUR MANUSCRIPT STATEMENTS
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 11 of the AKOYA analysis series. REVISION 1.
READ-ONLY. Writes to its own directory and changes nothing any other script owns.

WHY THIS SCRIPT EXISTS

    Four statements are going into the manuscript. Each is assembled from
    numbers that live in different scripts, and some components have never been
    tested at the level the design supports. This script tests every component
    of every statement, in one place, at the ANIMAL level, and prints a verdict
    per component rather than per statement. Figures come after, and only for
    components that pass.

    THE FOUR STATEMENTS AS AGREED (1 October 2026)

    S1  In D1MT-treated animals the tissue does not reach the myeloid density
        used here to mark granuloma-grade disease. The residual foci that are
        found are fewer, smaller and less prominent.
    S2  IDO1-positive macrophages are nearly absent from treated lung as a
        whole, while the macrophages remaining inside residual foci are still
        IDO1-positive at a fraction that is lower but overlapping.
    S3  In untreated foci, IDO1-positive macrophages sit further from plasma
        cells than chance allows. In treated foci they sit where chance puts
        them, so the effect is the loss of an exclusion.
    S4  B cells sit closer to the core of foci in treated animals.

WHAT EACH CHECK IS FOR, AND WHICH CRITICISM IT ANSWERS

    A  S1, the three adjectives, per animal. Only "fewer" has ever been tested
       at the animal level. "Smaller" and "less prominent" have been reported as
       ARM MEDIANS only: 176 against 219 um inscribed, 250 against 283 um
       equivalent, 2,882 against 7,012 cells/mm2 peak. Equivalent radius is a
       1.13x difference and the treated arm has four structures in total, two
       animals contributing one each, so two of the three treated "medians" are
       a single focus. An adjective that has not been tested does not go in a
       sentence.

    B  S2, the within-focus IDO1 fraction read from the RIGHT COLUMN.
       Every within-focus number quoted so far was derived from focus
       membership counts, which include cuff cells. Table 35 carries
       pct_ido1_pos_of_mac_core, which is core only. This has been an open item
       since 19 September and it is now in the figure path.

    C  S2, THE BATCH CONTROL. Scan is perfectly confounded with arm: both
       acquisitions ran 20251103, one per arm, and script 02 established that no
       between-arm intensity comparison can separate staining from biology. S2
       is exactly that comparison. The defence available to us is to show that
       markers with no reason to differ between arms do not differ, while IDO1
       differs by an order of magnitude. If the control markers also separate,
       S2 is in trouble and we need to know that before a reviewer tells us.
       This is the single highest-value check in the script.

    D  S2 against S1. If treated lungs have no granuloma-grade tissue they have
       fewer activated macrophages, so the whole-section IDO1-positive fraction
       falls as an arithmetic consequence. This quantifies how much of S2 is
       explained by the macrophage compartment shrinking, so the claim can be
       stated at the level the data support rather than defended later.

    E  S3, specificity and leave-one-out. The IDO1-negative anchor must not
       separate, or the finding is about myeloid cells rather than the IDO1
       compartment. And the separation must survive dropping any single animal.

    F  S4, per-phenotype leave-one-out. Leave-one-out has only ever been run on
       the POOLED model. B cells draws 89.6 percent of its treated cells from
       43118, which is also the animal carrying the 1.59 shape ratio, and the
       confirmatory margin is 0.0010 radial units, about 0.18 um. If the
       coefficient dies without 43118 the finding is one animal.

    G  S4, BALT. B cells cluster in follicles and plasma cells are scattered,
       which is the obvious biological reason the two behave differently in the
       radial and distance analyses. If the B cell effect disappears once BALT
       B cells are excluded, the finding is about follicle placement rather than
       B cell behaviour.

VALIDATION GATES, AND WHY THEY ARE HERE

    Checks F and G recompute the radial centring that script 08 owns. Two sets
    of numbers for one quantity is how this analysis has gone wrong before: it
    froze script 05 and it removed script 08's distance cell. So anything
    recomputed here is first validated against the script that owns it, and the
    dependent result is SUPPRESSED if the validation fails. The full-sample
    B cell coefficient must reproduce script 08 table 72 before any
    leave-one-out from this script is reported.

WHAT THIS SCRIPT DELIBERATELY DOES NOT DO

    No figures for the manuscript. Figures are script 12 and are gated on this
    output, because a figure built before the check is a figure that has to be
    withdrawn.

    No burden-threshold sensitivity sweep. That needs the density field and the
    tissue mask, which live in script 04, and rebuilding them here would create
    a second burden number. It belongs in script 04 as an additive diagnostic.

    No between-arm marker intensity comparison beyond the batch control, which
    exists to measure the confound rather than to look through it.

INPUTS
    structures_rev5/cell_assignments/<section>_cell_structures.csv
    structures_rev5/tables/35_foci_structures_relative.csv
    AKOYA/data/<section>.csv                    (raw panel, for the batch control)
    lymphocyte_radial/tables/72_models_all.csv  (validation target for F and G)
    composition_null/tables/83_per_animal_delta_three_references.csv  (check E)

OUTPUTS
    claim_checks/tables/90 .. 97
    claim_checks/figures/F70 .. F72   (diagnostic only, not manuscript figures)

USAGE
    conda activate sc_pre
    python AKOYA_11_Claim_Checks.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================

IN_DIR = "/master/jlehle/WORKING/AKOYA/structures_rev5"
CELL_DIR = f"{IN_DIR}/cell_assignments"
FOCI_TABLE = f"{IN_DIR}/tables/35_foci_structures_relative.csv"
DATA_DIR = "/master/jlehle/WORKING/AKOYA/data"
MODELS_08 = "/master/jlehle/WORKING/AKOYA/lymphocyte_radial/tables/72_models_all.csv"
DELTA_07B = ("/master/jlehle/WORKING/AKOYA/composition_null/tables/"
             "83_per_animal_delta_three_references.csv")
OUT_DIR = "/master/jlehle/WORKING/AKOYA/claim_checks"

CONDITION_ORDER = ["D1MT", "Untreated"]
CONDITION_COLORS = {"D1MT": "#2C7FB8", "Untreated": "#D95F02"}
REFERENCE_ARM = "Untreated"
STRUCT_COL = "focus_id"

IDO1_POS = "CD68+IDO1+ Macrophages"
IDO1_NEG = "CD68+IDO1- Macrophages"
CD68_LINEAGE = [IDO1_POS, IDO1_NEG, "CD163+ Macrophages"]
LYMPHOCYTES = ["Helper T cells", "CD4- T cells", "Tregs", "B cells",
               "Plasma cells"]

# ---- check A, the S1 adjectives --------------------------------------------
# Each component is (label, table 35 column, higher_in_untreated_expected).
S1_COMPONENTS = [
    ("fewer: foci per animal",          "_n_foci",                  True),
    ("smaller: inscribed radius um",    "max_inscribed_radius_um",  True),
    ("smaller: equivalent radius um",   "equiv_radius_um",          True),
    ("smaller: area um2",               "area_um2",                 True),
    ("less prominent: peak density",    "peak_density",             True),
    ("less prominent: fold over bg",    "fold_over_background",     True),
]

# ---- check C, the batch control --------------------------------------------
# Markers with no reason to differ between arms. Matched in the RAW export by
# COLUMN PREFIX, never by exact name, because the micron symbol and compartment
# suffixes encode inconsistently (see the project learnings).
BATCH_CONTROL_MARKERS = ["CD68", "CD3e", "CD20", "CD45", "Pan-Cytokeratin",
                         "CD163", "Ki67"]
BATCH_TEST_MARKERS = ["IDO1"]
INTENSITY_PERCENTILES = [50, 75, 90, 99]
# A control that separates the arms as strongly as the test marker is a
# problem. This is a LABEL, never a filter.
BATCH_FLAG_IF_CONTROL_SEPARATES = True

# ---- checks F and G, the radial recomputation ------------------------------
MIN_CELLS_FOR_BALANCE = 10
BALANCE_COMMON_MIN_FRAC = 0.90
# Validation: our full-sample coefficient must match script 08 table 72.
VALIDATE_AGAINST_08 = True
VALIDATION_MAX_ABS_DIFF = 0.002        # radial units
REFUSE_ON_VALIDATION_FAILURE = True
LOO_PHENOTYPES = ["B cells", "Plasma cells", "CD4- T cells"]

# ---- reportability, matching the rest of the series ------------------------
MIN_TREATED_CELLS_FOR_REPORT = 50
EXACT_P_REL_TOL = 1e-6

# ---- models ----------------------------------------------------------------
MODEL_MAX_CELLS_PER_STRUCTURE = 3000
MODEL_MIN_CELLS_PER_ANIMAL = 20
RANDOM_SEED = 0

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
WARN_COLOR = "#E08214"


# %% Cell 2 - imports, style, helpers
# =============================================================================

import os
import sys
import glob
import gc
import re
import warnings
from datetime import datetime
from itertools import combinations

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import statsmodels.formula.api as smf
    HAVE_SM = True
except Exception as _e:
    HAVE_SM = False
    _sm_err = _e

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


def normalise_label(s):
    """Strip zero-width and non-breaking characters. See the project learnings:
    the epithelial phenotype label carries a zero-width space, so equality
    matching fails silently without this."""
    if not isinstance(s, str):
        return s
    return re.sub(r"[​‌‍ ﻿]", "", s).strip()


# ---- the exact animal-level test, with its floor ---------------------------
# Treatment was assigned to six animals, three per arm, so there are C(6,3)=20
# assignments in ten sign-symmetric pairs and the two-sided p is floored at
# 2/20 = 0.10. Enumerating one representative per pair makes that floor hold by
# construction rather than by a tolerance surviving arithmetic noise.
def _assignment_reps(n, n_t):
    if n % 2 == 0 and n_t * 2 == n:
        for rest in combinations(range(1, n), n_t - 1):
            yield (0,) + rest
    else:
        for idx in combinations(range(n), n_t):
            yield idx


def exact_p_means(treated, untreated):
    """Two-sided exact randomization p on the per-animal values, plus the floor
    the design allows so 'as good as it gets' can be told from 'does not
    separate'."""
    t = [float(x) for x in treated if np.isfinite(x)]
    u = [float(x) for x in untreated if np.isfinite(x)]
    if len(t) < 1 or len(u) < 1:
        return np.nan, np.nan
    vals = np.array(t + u, float)
    n, k = len(vals), len(t)
    obs = abs(float(np.mean(vals[:k]) - np.mean(vals[k:])))
    stats = []
    for idx in _assignment_reps(n, k):
        a = vals[list(idx)]
        b = vals[[j for j in range(n) if j not in idx]]
        stats.append(abs(float(np.mean(a) - np.mean(b))))
    stats = np.asarray(stats)
    tol = EXACT_P_REL_TOL * max(1.0, abs(obs))
    p = float(np.mean(stats >= abs(obs) - tol))
    floor = 1.0 / len(stats)
    if p < floor - 1e-12:
        print(f"    EXACT-P FLOOR VIOLATED: {p:.4f} < {floor:.4f}. Raising; "
              f"this indicates a defect.")
        p = floor
    return p, floor


def complete_separation(t, u):
    t = [x for x in t if np.isfinite(x)]
    u = [x for x in u if np.isfinite(x)]
    if not t or not u:
        return False
    return (max(t) < min(u)) or (max(u) < min(t))


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


VERDICTS = []


def record(statement, component, supported, evidence, caveat=""):
    VERDICTS.append({"statement": statement, "component": component,
                     "supported": supported, "evidence": evidence,
                     "caveat": caveat})


def per_animal_test(label, statement, t_vals, u_vals, fmt="{:.3f}",
                    caveat="", units=""):
    """One line of output and one verdict row, for any per-animal quantity."""
    sep = complete_separation(t_vals, u_vals)
    mg = arm_margin(t_vals, u_vals)
    p, floor = exact_p_means(t_vals, u_vals)
    at_floor = bool(np.isfinite(p) and np.isfinite(floor) and p <= floor + 1e-9)
    vals = "  ".join(fmt.format(v) if np.isfinite(v) else "na"
                     for v in list(t_vals) + list(u_vals))
    tag = ("SEPARATED" if sep else "overlaps")
    pstr = f"p_exact={p:.3f}" if np.isfinite(p) else "p_exact=na"
    extra = "" if at_floor else "  <-- animals do not separate"
    print(f"    {label:<34}{vals:>46}   {tag:>10}  {pstr}{extra}")
    _t = [round(float(v), 4) for v in t_vals]
    _u = [round(float(v), 4) for v in u_vals]
    record(statement, label, bool(sep and at_floor),
           f"treated {_t} untreated {_u}; margin {mg:.4g}{units}; "
           f"p_exact {p:.3f}", caveat)
    return sep, mg, p, floor


_tee = Tee(os.path.join(TAB_DIR, "00_claim_checks_report.txt"))
sys.stdout = _tee

banner("AKOYA CLAIM CHECKS FOR THE FOUR MANUSCRIPT STATEMENTS (script 11 rev 1)")
print(f"Run time : {datetime.now().isoformat(timespec='seconds')}")
print(f"Input    : {IN_DIR}")
print(f"Output   : {OUT_DIR}")
if not HAVE_SM:
    print(f"\n    WARNING: statsmodels unavailable ({_sm_err}). Checks F and G "
          f"will be skipped.")

print("""
HOW TO READ THIS SCRIPT

    Every check is at the ANIMAL level, because treatment was assigned to
    animals. Cell counts are printed beside every number because they are
    large and the number of independent units is six.

    p_exact is the two-sided exact randomization p on the six animal-level
    values. It is FLOORED AT 0.10 by the design. A value AT the floor means
    the six animals separate completely, which is the strongest this design can
    produce. A value ABOVE the floor means they do not separate, and that is
    flagged, because a model p can look convincing while the animals disagree.

    The verdict table at the end is per COMPONENT, not per statement. An
    adjective that has not been tested does not go in a sentence.
""")


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
    need = [STRUCT_COL, "pheno", "condition", "region"]
    miss = [c for c in need if c not in d.columns]
    if miss:
        print(f"    ERROR: {sid} lacks {miss}. Skipping.")
        continue
    d["pheno"] = d["pheno"].map(normalise_label)
    cells[sid] = d
    n_str = int((d.loc[d[STRUCT_COL] > 0, STRUCT_COL]).nunique())
    print(f"    {sid:<12} {d['condition'].iloc[0]:<10} {len(d):>8,} cells total, "
          f"{int((d[STRUCT_COL] > 0).sum()):>7,} in {n_str:>2} structures")

if not cells:
    print("ERROR: nothing loaded.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

cond_rank = {c: i for i, c in enumerate(CONDITION_ORDER)}
SAMPLE_ORDER = sorted(cells, key=lambda s: (
    cond_rank.get(cells[s]["condition"].iloc[0], 9), s))
COND_OF = {s: cells[s]["condition"].iloc[0] for s in SAMPLE_ORDER}
TREATED = [s for s in SAMPLE_ORDER if COND_OF[s] == "D1MT"]
UNTREATED = [s for s in SAMPLE_ORDER if COND_OF[s] != "D1MT"]
MARKER_OF = {s: POINT_MARKERS[i % len(POINT_MARKERS)]
             for i, s in enumerate(SAMPLE_ORDER)}

print(f"\n    column order for every per-animal row below:")
print(f"      {'  '.join(short_label(s) for s in SAMPLE_ORDER)}")
print(f"      treated: {[short_label(s) for s in TREATED]}   "
      f"untreated: {[short_label(s) for s in UNTREATED]}")

foci = pd.DataFrame()
if os.path.exists(FOCI_TABLE):
    foci = pd.read_csv(FOCI_TABLE)
    print(f"\n    table 35: {len(foci)} focus records, "
          f"{foci['sample_id'].nunique()} sections")
    _has_area_ref = "mean_radial_over_area" in foci.columns
    print(f"    mean_radial_over_area present: {_has_area_ref}"
          f"{'' if _has_area_ref else '   (script 04 revision 6 exports it)'}")
else:
    print(f"\n    WARNING: {FOCI_TABLE} not found. Checks A and B will be "
          f"skipped.")


# %% Cell 4 - CHECK A: the three adjectives in statement 1, per animal
# =============================================================================

banner("CHECK A - STATEMENT 1, THE THREE ADJECTIVES, PER ANIMAL")

print("    Only 'fewer' has previously been tested at the animal level. Every")
print("    size and prominence number quoted so far has been an ARM MEDIAN")
print("    across structures: 176 against 219 um inscribed, 250 against 283 um")
print("    equivalent, 2,882 against 7,012 cells/mm2 peak. Equivalent radius is")
print("    a 1.13x difference, which is not obviously 'smaller'.\n")

if len(foci):
    print("    STRUCTURES PER ANIMAL, read this before the adjectives:")
    for s in SAMPLE_ORDER:
        n = int((foci["sample_id"] == s).sum())
        flag = "   <-- a single focus, so its 'median' is one structure" \
            if n == 1 else ""
        print(f"      {short_label(s):<8}{COND_OF[s]:<12}{n:>3} foci{flag}")
    print()

    a_rows = []
    print(f"    {'component':<34}{'per-animal values (treated | untreated)':>46}"
          f"   {'':>10}")
    print("    " + "-" * 100)
    for label, col, _ in S1_COMPONENTS:
        if col == "_n_foci":
            t = [float((foci["sample_id"] == s).sum()) for s in TREATED]
            u = [float((foci["sample_id"] == s).sum()) for s in UNTREATED]
            fmt = "{:.0f}"
        elif col in foci.columns:
            t = [float(foci.loc[foci["sample_id"] == s, col].median())
                 for s in TREATED]
            u = [float(foci.loc[foci["sample_id"] == s, col].median())
                 for s in UNTREATED]
            fmt = "{:.0f}" if foci[col].median() > 20 else "{:.2f}"
        else:
            print(f"    {label:<34}  COLUMN '{col}' ABSENT from table 35, skipped")
            record("S1", label, False, f"column {col} absent", "not testable")
            continue
        cav = ("two treated animals contribute a single focus, so their median "
               "is one structure") if col != "_n_foci" else ""
        sep, mg, pv, fl = per_animal_test(label, "S1", t, u, fmt=fmt, caveat=cav)
        a_rows.append({"component": label, "column": col,
                       **{f"t_{short_label(s)}": v for s, v in zip(TREATED, t)},
                       **{f"u_{short_label(s)}": v for s, v in zip(UNTREATED, u)},
                       "separated": sep, "margin": mg, "p_exact": pv,
                       "p_floor": fl})
    if a_rows:
        write_csv(pd.DataFrame(a_rows), "90_s1_adjectives_per_animal.csv")

    sub("What this means for the wording of statement 1")
    _ok = [r["component"] for r in a_rows if r["separated"]
           and np.isfinite(r["p_exact"]) and r["p_exact"] <= r["p_floor"] + 1e-9]
    _no = [r["component"] for r in a_rows if r["component"] not in _ok]
    print(f"    components with animal-level support ({len(_ok)}):")
    for c in _ok:
        print(f"      {c}")
    print(f"    components WITHOUT animal-level support ({len(_no)}):")
    for c in _no:
        print(f"      {c}")
    print()
    print("    Write the adjectives the checks support and drop the others.")
    print("    'Fewer' and 'less prominent' are different claims from 'smaller'")
    print("    and should not be carried by one another.")
else:
    print("    SKIPPED: table 35 not available.")


# %% Cell 5 - CHECK B: statement 2's within-focus fraction, from the right column
# =============================================================================

banner("CHECK B - STATEMENT 2, WITHIN-FOCUS IDO1 FROM pct_ido1_pos_of_mac_core")

print("    Every within-focus IDO1 number quoted so far (57.9, 65.3, 3.6 against")
print("    72.2, 95.9, 51.3) was derived from focus MEMBERSHIP counts, which")
print("    include cuff cells. Table 35 carries pct_ido1_pos_of_mac_core, which")
print("    is CORE ONLY and is the column the manuscript should quote. Open")
print("    since 19 September; it is in the figure path now.\n")

whole_section = {}
for s in SAMPLE_ORDER:
    d = cells[s]
    mac = d.loc[d["pheno"].isin(CD68_LINEAGE)]
    n_mac = int(len(mac))
    n_pos = int((d["pheno"] == IDO1_POS).sum())
    whole_section[s] = (100.0 * n_pos / n_mac if n_mac else np.nan, n_mac, n_pos)

if len(foci) and "pct_ido1_pos_of_mac_core" in foci.columns:
    print(f"    {'animal':<10}{'arm':<12}{'whole section %':>17}{'n mac':>9}"
          f"{'within-focus % (core)':>24}{'n mac core':>12}")
    print("    " + "-" * 84)
    b_rows = []
    for s in SAMPLE_ORDER:
        ws, n_mac, n_pos = whole_section[s]
        f_ = foci.loc[foci["sample_id"] == s]
        wf = float(f_["pct_ido1_pos_of_mac_core"].median()) if len(f_) else np.nan
        nmc = int(f_["n_macrophages_core"].sum()) if "n_macrophages_core" \
            in f_.columns else -1
        print(f"    {short_label(s):<10}{COND_OF[s]:<12}{ws:>16.1f}%{n_mac:>9,}"
              f"{wf:>23.1f}%{nmc:>12,}")
        b_rows.append({"sample_id": s, "animal_id": short_label(s),
                       "condition": COND_OF[s], "pct_whole_section": ws,
                       "n_macrophages_section": n_mac,
                       "n_ido1_pos_section": n_pos,
                       "pct_within_focus_core_median": wf,
                       "n_macrophages_core": nmc,
                       "n_foci": int(len(f_))})
    bdf = pd.DataFrame(b_rows)
    write_csv(bdf, "91_s2_ido1_fraction_two_ways.csv")

    print()
    per_animal_test("whole section IDO1+ of mac %", "S2",
                    bdf.loc[bdf["condition"] == "D1MT", "pct_whole_section"].tolist(),
                    bdf.loc[bdf["condition"] != "D1MT", "pct_whole_section"].tolist(),
                    fmt="{:.1f}", units=" pct pts")
    per_animal_test("within-focus IDO1+ of mac %", "S2",
                    bdf.loc[bdf["condition"] == "D1MT",
                            "pct_within_focus_core_median"].tolist(),
                    bdf.loc[bdf["condition"] != "D1MT",
                            "pct_within_focus_core_median"].tolist(),
                    fmt="{:.1f}", units=" pct pts",
                    caveat="43118 is a large outlier on this row")

    sub("The wording statement 2 can carry")
    print("    If the whole-section row separates and the within-focus row does")
    print("    not, the statement is that the IDO1-positive compartment is")
    print("    nearly absent from treated lung AS A WHOLE while the macrophages")
    print("    inside residual foci are still IDO1-positive at an OVERLAPPING")
    print("    fraction. Say overlapping, not 'not significantly different':")
    print("    the second is a claim about a test and the exact p here is")
    print("    floored at 0.10, so a non-significant result is uninformative")
    print("    about equivalence.")
else:
    print("    SKIPPED: pct_ido1_pos_of_mac_core not in table 35.")


# %% Cell 6 - CHECK C: the batch control for statement 2
# =============================================================================

banner("CHECK C - STATEMENT 2, THE BATCH CONTROL")

print("    Scan is PERFECTLY CONFOUNDED with arm. Both acquisitions ran")
print("    20251103, one per arm, and script 02 established that no between-arm")
print("    intensity comparison can separate staining from biology. Statement 2")
print("    is exactly that comparison.")
print()
print("    What can be shown: markers with no reason to differ between arms")
print("    should not differ, while IDO1 differs by an order of magnitude. That")
print("    makes a global staining shift implausible without proving its")
print("    absence. If the control markers separate as strongly as IDO1 does,")
print("    statement 2 cannot be defended on intensity and we need to know now.")
print()
print("    This does NOT look through the confound. It measures it.\n")

raw_paths = sorted(glob.glob(os.path.join(DATA_DIR, "*.csv")))
c_rows = []
if not raw_paths:
    print(f"    SKIPPED: no raw exports in {DATA_DIR}.")
else:
    def match_marker_cols(cols, marker):
        """Match by PREFIX, never exact name: the micron symbol and compartment
        suffixes encode inconsistently in these exports."""
        m = marker.lower()
        return [c for c in cols
                if normalise_label(str(c)).lower().startswith(m + ":")]

    for rp in raw_paths:
        stem = os.path.basename(rp).replace(".csv", "")
        sid = next((s for s in SAMPLE_ORDER if s in stem or stem in s), None)
        if sid is None:
            continue
        try:
            rd = pd.read_csv(rp, low_memory=False)
        except Exception as e:
            print(f"    ERROR reading {stem}: {e}. Skipping.")
            continue
        rd.columns = [normalise_label(c) for c in rd.columns]
        for marker in BATCH_CONTROL_MARKERS + BATCH_TEST_MARKERS:
            mcols = match_marker_cols(rd.columns, marker)
            if not mcols:
                continue
            col = mcols[0]
            v = pd.to_numeric(rd[col], errors="coerce").dropna().to_numpy(float)
            if not len(v):
                continue
            rec = {"sample_id": sid, "animal_id": short_label(sid),
                   "condition": COND_OF[sid], "marker": marker,
                   "column_used": col, "n_cells": int(len(v)),
                   "kind": "test" if marker in BATCH_TEST_MARKERS else "control"}
            for q in INTENSITY_PERCENTILES:
                rec[f"p{q}"] = float(np.percentile(v, q))
            c_rows.append(rec)
        del rd
        gc.collect()

    if not c_rows:
        print("    SKIPPED: no marker columns matched. Check the export column")
        print("    naming against BATCH_CONTROL_MARKERS.")
    else:
        cdf = pd.DataFrame(c_rows)
        write_csv(cdf, "92_batch_control_intensities.csv")
        found = sorted(cdf["marker"].unique())
        missing = [m for m in BATCH_CONTROL_MARKERS + BATCH_TEST_MARKERS
                   if m not in found]
        print(f"    markers matched: {found}")
        if missing:
            print(f"    NOT FOUND in the export: {missing}")
            print("    A control that is absent is not a passed control. Say so.")

        for q in INTENSITY_PERCENTILES:
            sub(f"Arm comparison at the {q}th percentile of intensity")
            print(f"    {'marker':<18}{'kind':>9}"
                  + "".join(f"{short_label(s):>11}" for s in SAMPLE_ORDER)
                  + f"{'separated':>11}{'p_exact':>9}")
            print("    " + "-" * (18 + 9 + 11 * len(SAMPLE_ORDER) + 20))
            for marker in found:
                g = cdf.loc[cdf["marker"] == marker]
                vals, t, u = [], [], []
                for s in SAMPLE_ORDER:
                    r_ = g.loc[g["sample_id"] == s, f"p{q}"]
                    x = float(r_.iloc[0]) if len(r_) else np.nan
                    vals.append(x)
                    (t if COND_OF[s] == "D1MT" else u).append(x)
                sep = complete_separation(t, u)
                pv, fl = exact_p_means(t, u)
                kind = g["kind"].iloc[0]
                row = (f"    {marker:<18}{kind:>9}"
                       + "".join(f"{v:>11.1f}" if np.isfinite(v) else f"{'na':>11}"
                                 for v in vals)
                       + f"{('YES' if sep else '-'):>11}"
                       + (f"{pv:>9.3f}" if np.isfinite(pv) else f"{'na':>9}"))
                if kind == "control" and sep and BATCH_FLAG_IF_CONTROL_SEPARATES:
                    row += "   <-- CONTROL SEPARATES"
                print(row)

        sub("VERDICT ON THE BATCH CONTROL")
        print("    MAGNITUDE, NOT SEPARATION COUNT. With three animals per arm")
        print("    any marker separates by chance with probability 0.10, so")
        print("    counting separations says little: a control can separate on a")
        print("    margin of 0.1 intensity units and mean nothing. This is the")
        print("    same error flagged in script 08's spillover verdict, which")
        print("    counted separations without checking direction or size.")
        print()
        print("    The statistic that matters is the SIZE of the arm difference:")
        print("    IDO1 should differ by a large fold while the controls differ")
        print("    barely at all. Fold is untreated median over treated median")
        print("    of the per-animal values.\n")
        prim = INTENSITY_PERCENTILES[-1]
        ctrl, test = [], []
        for marker in found:
            g = cdf.loc[cdf["marker"] == marker]
            t = [float(g.loc[g["sample_id"] == s, f"p{prim}"].iloc[0])
                 for s in SAMPLE_ORDER if COND_OF[s] == "D1MT"
                 and len(g.loc[g["sample_id"] == s])]
            u = [float(g.loc[g["sample_id"] == s, f"p{prim}"].iloc[0])
                 for s in SAMPLE_ORDER if COND_OF[s] != "D1MT"
                 and len(g.loc[g["sample_id"] == s])]
            if not (t and u):
                continue
            mt, mu = float(np.median(t)), float(np.median(u))
            fold = (mu / mt) if mt > 0 else np.nan
            (test if g["kind"].iloc[0] == "test" else ctrl).append(
                (marker, complete_separation(t, u), arm_margin(t, u), fold,
                 mt, mu))

        print(f"    {'marker':<18}{'kind':>9}{'treated med':>13}"
              f"{'untreated med':>15}{'fold U/T':>11}{'separated':>11}")
        print("    " + "-" * 77)
        for marker, sep_, mg_, fold_, mt_, mu_ in sorted(
                ctrl + test, key=lambda r: -(r[3] if np.isfinite(r[3]) else 0)):
            kind = "test" if marker in BATCH_TEST_MARKERS else "control"
            print(f"    {marker:<18}{kind:>9}{mt_:>13.2f}{mu_:>15.2f}"
                  f"{fold_:>11.2f}{('YES' if sep_ else '-'):>11}")

        ctrl_folds = [f for _, _, _, f, _, _ in ctrl if np.isfinite(f)]
        test_folds = {m: f for m, _, _, f, _, _ in test if np.isfinite(f)}
        max_ctrl = max(ctrl_folds) if ctrl_folds else np.nan
        ido1_fold = test_folds.get("IDO1", np.nan)
        n_ctrl_sep = sum(1 for _, s_, _, _, _, _ in ctrl if s_)
        test_sep = [m for m, s_, _, _, _, _ in test if s_]
        print(f"\n    largest control fold : {max_ctrl:.2f}"
              f"   ({len(ctrl_folds)} controls)")
        print(f"    IDO1 fold            : {ido1_fold:.2f}")
        if np.isfinite(max_ctrl) and np.isfinite(ido1_fold) and max_ctrl > 0:
            print(f"    IDO1 exceeds the largest control by "
                  f"{ido1_fold / max_ctrl:.2f}x")
        print(f"    (separation counts, secondary: controls {n_ctrl_sep} of "
              f"{len(ctrl)}, test {len(test_sep)} of {len(test)})")
        print()
        _ido1_dominates = bool(np.isfinite(max_ctrl) and np.isfinite(ido1_fold)
                               and ido1_fold >= 2.0 * max_ctrl)
        if _ido1_dominates:
            print("    IDO1's arm difference is at least twice the largest")
            print("    control's. A global staining or acquisition shift would be")
            print("    expected to move the controls comparably, so this is the")
            print("    strongest available argument that the IDO1 difference is")
            print("    not batch.")
            print("    It is NOT proof. The scan confound is structural and")
            print("    cannot be removed by any analysis of these data. State it")
            print("    in Methods as a limitation alongside this control.")
            record("S2", "batch control", True,
                   f"IDO1 fold {ido1_fold:.2f} against largest control fold "
                   f"{max_ctrl:.2f} at p{prim}",
                   "scan is confounded with arm; this bounds the concern rather "
                   "than removing it")
        elif np.isfinite(ido1_fold) and np.isfinite(max_ctrl):
            print("    IDO1's arm difference is NOT clearly larger than the")
            print("    controls'. Markers that should not differ between arms")
            print("    differ comparably, so the IDO1 intensity difference")
            print("    cannot be attributed to biology on these data.")
            print("    Statement 2 must then rest on the PHENOTYPE CALLS with")
            print("    the vendor-classifier caveat stated, or be softened to a")
            print("    description of what the calls show without a mechanistic")
            print("    attribution. Check which control is driving this before")
            print("    deciding: a single bright-channel control behaving oddly")
            print("    is different from the whole panel shifting.")
            record("S2", "batch control", False,
                   f"IDO1 fold {ido1_fold:.2f} against largest control fold "
                   f"{max_ctrl:.2f} at p{prim}",
                   "intensity-based attribution is not available")
        else:
            print("    Could not compute folds. Read table 92 directly.")
            record("S2", "batch control", False,
                   "folds not computable", "read table 92")


# %% Cell 7 - CHECK D: how much of statement 2 is statement 1
# =============================================================================

banner("CHECK D - IS STATEMENT 2 A CONSEQUENCE OF STATEMENT 1?")

print("    If treated lungs have no granuloma-grade tissue they have fewer")
print("    activated macrophages, so the whole-section IDO1-positive FRACTION")
print("    can fall without anything changing about IDO1 regulation. This")
print("    separates the two and is worth doing before a reviewer does it.\n")

d_rows = []
for s in SAMPLE_ORDER:
    d = cells[s]
    n_all = int(len(d))
    mac = d.loc[d["pheno"].isin(CD68_LINEAGE)]
    in_str = d[STRUCT_COL] > 0
    d_rows.append({
        "sample_id": s, "animal_id": short_label(s), "condition": COND_OF[s],
        "n_cells_section": n_all,
        "n_macrophages_section": int(len(mac)),
        "pct_macrophages_of_section": 100.0 * len(mac) / n_all if n_all else np.nan,
        "n_ido1_pos_section": int((d["pheno"] == IDO1_POS).sum()),
        "n_ido1_pos_in_structures": int(((d["pheno"] == IDO1_POS) & in_str).sum()),
        "pct_ido1_pos_inside_structures": (
            100.0 * ((d["pheno"] == IDO1_POS) & in_str).sum()
            / max(int((d["pheno"] == IDO1_POS).sum()), 1)),
        "n_cells_in_structures": int(in_str.sum()),
        "pct_section_in_structures": 100.0 * in_str.mean() if n_all else np.nan,
    })
ddf = pd.DataFrame(d_rows)
write_csv(ddf, "93_s2_against_s1_decomposition.csv")

print(f"    {'animal':<10}{'arm':<12}{'mac % of section':>18}"
      f"{'IDO1+ n':>10}{'% of IDO1+ inside foci':>24}{'% section in foci':>19}")
print("    " + "-" * 93)
for _, r in ddf.iterrows():
    print(f"    {r['animal_id']:<10}{r['condition']:<12}"
          f"{r['pct_macrophages_of_section']:>17.1f}%{r['n_ido1_pos_section']:>10,}"
          f"{r['pct_ido1_pos_inside_structures']:>23.1f}%"
          f"{r['pct_section_in_structures']:>18.1f}%")

print()
per_animal_test("macrophage % of section", "S2",
                ddf.loc[ddf["condition"] == "D1MT",
                        "pct_macrophages_of_section"].tolist(),
                ddf.loc[ddf["condition"] != "D1MT",
                        "pct_macrophages_of_section"].tolist(),
                fmt="{:.1f}", units=" pct pts",
                caveat="if this separates, part of the IDO1 fraction difference "
                       "is the macrophage compartment shrinking")
per_animal_test("% of IDO1+ cells inside foci", "S2",
                ddf.loc[ddf["condition"] == "D1MT",
                        "pct_ido1_pos_inside_structures"].tolist(),
                ddf.loc[ddf["condition"] != "D1MT",
                        "pct_ido1_pos_inside_structures"].tolist(),
                fmt="{:.1f}", units=" pct pts")

sub("How to read this")
print("    If the macrophage percentage of the section does NOT separate while")
print("    the IDO1-positive fraction does, then statement 2 is about IDO1 and")
print("    not about the macrophage compartment, and it stands on its own.")
print("    If it DOES separate, part of statement 2 is statement 1 restated,")
print("    and the honest claim becomes that the IDO1-positive compartment")
print("    scales with disease while being induced at an indistinguishable")
print("    rate inside residual disease. That is a more interesting sentence")
print("    than the one currently drafted and it pre-empts the criticism.")


# %% Cell 8 - CHECK E: statement 3, specificity and leave-one-out
# =============================================================================

banner("CHECK E - STATEMENT 3, SPECIFICITY AND LEAVE-ONE-OUT")

print("    Statement 3 is about IDO1-POSITIVE macrophages. If the")
print("    IDO1-negative anchor separates too, the finding is about myeloid")
print("    cells and the word IDO1-positive is doing no work. Read from script")
print("    07b table 83 rather than recomputed, so there is one set of numbers.\n")

if os.path.exists(DELTA_07B):
    dd = pd.read_csv(DELTA_07B)
    dd["animal_id"] = dd["animal_id"].astype(str)
    refs = [c for c in dd.columns if c.startswith("delta_")]
    for ref in refs:
        sub(f"Reference: {ref}")
        for anchor in [IDO1_POS, IDO1_NEG]:
            g = dd.loc[(dd["anchor"] == anchor) & (dd["target"] == "Plasma cells")]
            if not len(g):
                continue
            t = g.loc[g["condition"] == "D1MT", ref].tolist()
            u = g.loc[g["condition"] != "D1MT", ref].tolist()
            n_t_tgt = g.loc[g["condition"] == "D1MT",
                            "n_targets_min_structure"].min() \
                if "n_targets_min_structure" in g.columns else -1
            lab = f"{anchor.replace('Macrophages', 'Mac')} -> Plasma"
            thin = "" if n_t_tgt >= MIN_TREATED_CELLS_FOR_REPORT else "  [THIN]"
            sep, mg, pv, fl = per_animal_test(lab + thin, "S3", t, u,
                                              fmt="{:+.2f}", units=" um")
            # leave-one-out on the per-animal values
            if sep:
                worst = None
                for drop in range(len(t)):
                    tt = [v for i, v in enumerate(t) if i != drop]
                    m2 = arm_margin(tt, u)
                    worst = m2 if worst is None else min(worst, m2) \
                        if np.isfinite(m2) else np.nan
                for drop in range(len(u)):
                    uu = [v for i, v in enumerate(u) if i != drop]
                    m2 = arm_margin(t, uu)
                    worst = m2 if worst is None else min(worst, m2) \
                        if np.isfinite(m2) else np.nan
                print(f"      leave-one-out: worst remaining margin "
                      f"{worst if worst is not None else float('nan'):.2f} um"
                      + ("   (separation survives dropping any one animal)"
                         if worst is not None and np.isfinite(worst) and worst > 0
                         else "   <-- separation does NOT survive"))

    sub("Specificity verdict")
    print("    Statement 3 needs the IDO1-positive anchor to separate and the")
    print("    IDO1-negative one not to. If both separate, drop the word")
    print("    IDO1-positive and describe it as a myeloid result. If neither")
    print("    does, statement 3 is not supported on this reference.")
else:
    print(f"    SKIPPED: {DELTA_07B} not found. Run script 07b first.")


# %% Cell 9 - CHECKS F and G: statement 4, leave-one-out and BALT
# =============================================================================

banner("CHECK F/G - STATEMENT 4, PER-PHENOTYPE LEAVE-ONE-OUT AND BALT")

print("    B cells draws 89.6 percent of its treated core cells from 43118,")
print("    which is also the animal carrying the 1.59 shape ratio, and the")
print("    confirmatory margin is 0.0010 radial units, about 0.18 um. Script")
print("    08 runs leave-one-out on the POOLED model only, never per phenotype.")
print()
print("    This cell recomputes script 08's centring, so it is VALIDATED against")
print("    script 08 table 72 before any leave-one-out is reported. Two sets of")
print("    numbers for one quantity is how script 05 came to be frozen.\n")

core_frames = []
for s in SAMPLE_ORDER:
    d = cells[s]
    d = d.loc[(d[STRUCT_COL] > 0) & (d["region"] == "core")].copy()
    d["radial_pos"] = pd.to_numeric(d["radial_pos"], errors="coerce")
    d = d.loc[np.isfinite(d["radial_pos"])]
    if not len(d):
        continue
    d["sample_id"] = s
    core_frames.append(d)
core = pd.concat(core_frames, ignore_index=True) if core_frames else pd.DataFrame()
del core_frames
gc.collect()

VALIDATION_F_PASSED = False
if not len(core):
    print("    SKIPPED: no core cells loaded.")
elif not HAVE_SM:
    print("    SKIPPED: statsmodels unavailable.")
else:
    # cell-weighted centring, exactly as script 08's PRIMARY outcome
    core["_ref_cells"] = core.groupby(["sample_id", STRUCT_COL])[
        "radial_pos"].transform("mean")
    core["radial_centered_core"] = core["radial_pos"] - core["_ref_cells"]

    # common-set balanced centring, over CORE cells, matching script 08
    ph = (core.groupby(["sample_id", STRUCT_COL, "pheno"])["radial_pos"]
          .agg(["mean", "size"]).reset_index())
    n_struct = int(core.groupby(["sample_id", STRUCT_COL]).ngroups)
    ok = ph.loc[ph["size"] >= MIN_CELLS_FOR_BALANCE]
    cover = (ok.groupby("pheno").size() / max(n_struct, 1))
    COMMON = sorted(cover.loc[cover >= BALANCE_COMMON_MIN_FRAC].index)
    print(f"    common set over CORE cells: {len(COMMON)} phenotypes")
    print(f"      {COMMON}")
    sel = ok.loc[ok["pheno"].isin(COMMON)]
    ref_c = sel.groupby(["sample_id", STRUCT_COL])["mean"].mean()
    key = list(zip(core["sample_id"], core[STRUCT_COL].astype(int)))
    core["_ref_common"] = [ref_c.get(k, np.nan) for k in key]
    core["radial_centered_core_balanced_common"] = (
        core["radial_pos"] - core["_ref_common"])

    if "mean_radial_over_area" in foci.columns:
        aref = {(r["sample_id"], int(r[STRUCT_COL])): float(r["mean_radial_over_area"])
                for _, r in foci.iterrows() if STRUCT_COL in foci.columns}
        core["_ref_area"] = [aref.get(k, np.nan) for k in key]
        core["radial_centered_core_area"] = core["radial_pos"] - core["_ref_area"]
        print("    area-weighted centring available from table 35")
    else:
        print("    area-weighted centring NOT available (script 04 rev 6 exports it)")

    OUTCOMES = [c for c in ["radial_centered_core_area",
                            "radial_centered_core_balanced_common",
                            "radial_centered_core"] if c in core.columns]

    def capped(df):
        if not len(df):
            return df
        return pd.concat(
            [g.sample(min(len(g), MODEL_MAX_CELLS_PER_STRUCTURE),
                      random_state=RANDOM_SEED)
             for _, g in df.groupby(["sample_id", STRUCT_COL])],
            ignore_index=True)

    def fit(df, outcome):
        d = df.dropna(subset=[outcome, "condition", "sample_id", STRUCT_COL]).copy()
        counts = d.groupby("sample_id").size()
        d = d.loc[d["sample_id"].isin(
            counts.loc[counts >= MODEL_MIN_CELLS_PER_ANIMAL].index)]
        if (not len(d) or d["condition"].nunique() < 2
                or d["sample_id"].nunique() < 3):
            return None
        d["arm"] = (d["condition"] != REFERENCE_ARM).astype(float)
        d["struct_key"] = (d["sample_id"].astype(str) + "_"
                           + d[STRUCT_COL].astype(int).astype(str))
        d["_y"] = pd.to_numeric(d[outcome], errors="coerce")
        d = d.dropna(subset=["_y"])
        if not len(d):
            return None
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res, mode = None, ""
            try:
                res = smf.mixedlm("_y ~ arm", data=d, groups=d["sample_id"],
                                  re_formula="1",
                                  vc_formula={"struct": "0 + C(struct_key)"}
                                  ).fit(reml=True, method="lbfgs", maxiter=300)
                mode = "animal + structure"
            except Exception:
                res = None
            def _usable(r):
                try:
                    return bool(np.isfinite(float(r.params["arm"]))
                                and np.isfinite(float(r.bse["arm"]))
                                and float(r.bse["arm"]) > 0
                                and np.isfinite(float(r.pvalues["arm"])))
                except Exception:
                    return False
            if res is None or not _usable(res):
                try:
                    res = smf.mixedlm("_y ~ arm", data=d, groups=d["sample_id"],
                                      re_formula="1").fit(reml=True,
                                                          method="lbfgs",
                                                          maxiter=300)
                    mode = "animal only"
                except Exception:
                    return None
            if res is None or not _usable(res):
                return None
        return {"coef": float(res.params["arm"]), "p": float(res.pvalues["arm"]),
                "n_cells": int(len(d)), "fit_mode": mode}

    # ---- VALIDATION against script 08 ------------------------------------
    sub("Validation: does this reproduce script 08 table 72?")
    if not os.path.exists(MODELS_08):
        print(f"    CANNOT VALIDATE: {MODELS_08} not found.")
    else:
        m08 = pd.read_csv(MODELS_08)
        vrows = []
        for pheno in LOO_PHENOTYPES:
            for oc in OUTCOMES:
                mine = fit(capped(core.loc[core["pheno"] == pheno]), oc)
                if mine is None:
                    continue
                cand = m08.loc[m08["analysis"].astype(str).str.contains(
                    re.escape(f"radial: {pheno}"), na=False)]
                if "outcome" in m08.columns:
                    cand = cand.loc[cand["outcome"] == oc]
                if not len(cand):
                    vrows.append({"phenotype": pheno, "outcome": oc,
                                  "mine": mine["coef"], "script08": np.nan,
                                  "abs_diff": np.nan, "matched": False})
                    continue
                theirs = float(cand["coef_D1MT_vs_ref"].iloc[0])
                diff = abs(mine["coef"] - theirs)
                vrows.append({"phenotype": pheno, "outcome": oc,
                              "mine": mine["coef"], "script08": theirs,
                              "abs_diff": diff,
                              "matched": bool(diff <= VALIDATION_MAX_ABS_DIFF)})
        vdf = pd.DataFrame(vrows)
        if len(vdf):
            write_csv(vdf, "94_validation_against_script08.csv")
            print(f"    {'phenotype':<16}{'outcome':<42}{'here':>10}{'08':>10}"
                  f"{'|diff|':>9}{'ok':>5}")
            print("    " + "-" * 92)
            for _, r in vdf.iterrows():
                print(f"    {r['phenotype']:<16}{str(r['outcome'])[:41]:<42}"
                      f"{r['mine']:>+10.4f}"
                      + (f"{r['script08']:>+10.4f}" if np.isfinite(r['script08'])
                         else f"{'na':>10}")
                      + (f"{r['abs_diff']:>9.4f}" if np.isfinite(r['abs_diff'])
                         else f"{'na':>9}")
                      + f"{('yes' if r['matched'] else 'NO'):>5}")
            comparable = vdf.loc[np.isfinite(vdf["abs_diff"])]
            VALIDATION_F_PASSED = bool(len(comparable)
                                       and comparable["matched"].all())
            if VALIDATION_F_PASSED:
                print("\n    PASS. The centring here reproduces script 08, so a")
                print("    leave-one-out computed here is a leave-one-out of")
                print("    script 08's model.")
            else:
                print("\n    FAIL. This script's centring does not reproduce")
                print("    script 08. The leave-one-out below would be a")
                print("    leave-one-out of a DIFFERENT model, which is exactly")
                print("    the second-set-of-numbers problem. Suppressed.")

    # ---- CHECK F: per-phenotype leave-one-out ----------------------------
    if VALIDATION_F_PASSED or not REFUSE_ON_VALIDATION_FAILURE:
        sub("CHECK F - per-phenotype leave-one-out")
        print("    Drop each animal in turn and refit. Read the COEFFICIENT")
        print("    spread, not the p: dropping an animal takes the design to")
        print("    2 versus 3 and the p moves on degrees of freedom alone.\n")
        f_rows = []
        for pheno in LOO_PHENOTYPES:
            dp = core.loc[core["pheno"] == pheno]
            n_t = {short_label(s): int((dp["sample_id"] == s).sum())
                   for s in TREATED}
            tot_t = sum(n_t.values())
            dom = max(n_t, key=n_t.get) if n_t else "na"
            print(f"    {pheno}   treated cells {tot_t:,}, "
                  f"{100 * n_t.get(dom, 0) / max(tot_t, 1):.1f}% from {dom}")
            for oc in OUTCOMES:
                full = fit(capped(dp), oc)
                if full is None:
                    continue
                coefs = {"full": full["coef"]}
                for s in SAMPLE_ORDER:
                    r_ = fit(capped(dp.loc[dp["sample_id"] != s]), oc)
                    coefs[f"drop_{short_label(s)}"] = r_["coef"] if r_ else np.nan
                vals = [v for k, v in coefs.items() if k != "full"
                        and np.isfinite(v)]
                spread = (max(vals) - min(vals)) if vals else np.nan
                sign_flip = any((v > 0) != (full["coef"] > 0) for v in vals)
                print(f"      {str(oc)[:44]:<46}full={full['coef']:>+8.4f}  "
                      f"spread={spread:>7.4f}"
                      + ("   <-- A LEAVE-ONE-OUT FLIPS SIGN" if sign_flip else ""))
                for k, v in coefs.items():
                    if k == "full":
                        continue
                    mark = ""
                    if np.isfinite(v) and (v > 0) != (full["coef"] > 0):
                        mark = "  <-- sign flip"
                    print(f"        {k:<16}{v:>+9.4f}{mark}")
                f_rows.append({"phenotype": pheno, "outcome": oc,
                               "treated_cells": tot_t,
                               "dominant_animal": dom,
                               "dominant_share_pct":
                                   100 * n_t.get(dom, 0) / max(tot_t, 1),
                               **coefs, "coef_spread": spread,
                               "any_sign_flip": sign_flip})
                record("S4" if pheno == "B cells" else "other",
                       f"{pheno} leave-one-out on {oc}",
                       not sign_flip,
                       f"full {full['coef']:+.4f}, spread {spread:.4f}",
                       f"{100 * n_t.get(dom, 0) / max(tot_t, 1):.0f}% of treated "
                       f"cells from {dom}")
        if f_rows:
            write_csv(pd.DataFrame(f_rows), "95_per_phenotype_leave_one_out.csv")

        # ---- CHECK G: BALT exclusion -------------------------------------
        sub("CHECK G - does the B cell result survive excluding BALT?")
        print("    B cells cluster in follicles and plasma cells are scattered,")
        print("    which is the obvious biological reason the two behave")
        print("    differently in the radial and the distance analyses. If the")
        print("    radial effect disappears once BALT B cells are dropped, the")
        print("    finding is about follicle placement, not B cell behaviour.\n")
        if "balt_id" not in core.columns:
            print("    SKIPPED: balt_id not in the cell assignments.")
        else:
            g_rows = []
            print(f"    {'animal':<10}{'arm':<12}{'B cells core':>14}"
                  f"{'in BALT':>10}{'% in BALT':>11}")
            print("    " + "-" * 57)
            for s in SAMPLE_ORDER:
                b = core.loc[(core["sample_id"] == s) & (core["pheno"] == "B cells")]
                nb = int(len(b))
                nin = int((pd.to_numeric(b["balt_id"], errors="coerce") > 0).sum())
                print(f"    {short_label(s):<10}{COND_OF[s]:<12}{nb:>14,}"
                      f"{nin:>10,}{(100 * nin / nb if nb else np.nan):>10.1f}%")
                g_rows.append({"sample_id": s, "animal_id": short_label(s),
                               "condition": COND_OF[s], "n_b_cells_core": nb,
                               "n_in_balt": nin,
                               "pct_in_balt": 100 * nin / nb if nb else np.nan})
            write_csv(pd.DataFrame(g_rows), "96_balt_b_cell_counts.csv")
            print()
            for pheno in ["B cells", "Plasma cells"]:
                dp = core.loc[core["pheno"] == pheno]
                dx = dp.loc[pd.to_numeric(dp["balt_id"], errors="coerce").fillna(0) <= 0]
                print(f"    {pheno}: {len(dp):,} core cells, "
                      f"{len(dx):,} outside BALT")
                for oc in OUTCOMES:
                    a = fit(capped(dp), oc)
                    b_ = fit(capped(dx), oc)
                    if a is None:
                        continue
                    if b_ is None:
                        print(f"      {str(oc)[:44]:<46}all={a['coef']:>+8.4f}  "
                              f"non-BALT model did not fit")
                        continue
                    died = (b_["coef"] > 0) != (a["coef"] > 0) or \
                        abs(b_["coef"]) < 0.5 * abs(a["coef"])
                    print(f"      {str(oc)[:44]:<46}all={a['coef']:>+8.4f}  "
                          f"non-BALT={b_['coef']:>+8.4f}"
                          + ("   <-- EFFECT LARGELY BALT" if died else ""))
                    record("S4" if pheno == "B cells" else "other",
                           f"{pheno} BALT exclusion on {oc}", not died,
                           f"all {a['coef']:+.4f}, non-BALT {b_['coef']:+.4f}",
                           "")
    else:
        print("\n    CHECKS F AND G SUPPRESSED: the validation above did not")
        print("    pass, and a leave-one-out of an unvalidated recomputation")
        print("    cannot be interpreted. Table 94 is written so the mismatch")
        print("    can be diagnosed.")


# %% Cell 10 - verdict
# =============================================================================

banner("VERDICT BY COMPONENT")

print("    Per COMPONENT, not per statement. A statement goes in the manuscript")
print("    with the components that passed and without the ones that did not.\n")

if VERDICTS:
    vdf = pd.DataFrame(VERDICTS)
    write_csv(vdf, "97_verdicts_by_component.csv")
    for stmt in ["S1", "S2", "S3", "S4", "other"]:
        g = vdf.loc[vdf["statement"] == stmt]
        if not len(g):
            continue
        sub(f"{stmt}")
        for _, r in g.iterrows():
            mark = "SUPPORTED    " if r["supported"] else "NOT SUPPORTED"
            print(f"    [{mark}] {r['component']}")
            print(f"        {r['evidence']}")
            if r["caveat"]:
                print(f"        caveat: {r['caveat']}")
else:
    print("    No verdicts recorded, which means every check was skipped. Read")
    print("    the warnings above.")

banner("WHAT THIS SCRIPT CANNOT TELL YOU")
print("""
    The scan confound is structural. Both acquisitions ran on one day, one per
    arm, so no analysis of these data can separate staining and acquisition
    batch from biology in a between-arm intensity comparison. Check C bounds
    the concern. It does not remove it, and Methods has to say so.

    One section per animal. Granulomatous disease is focal, so these sections
    cannot establish that treated lungs contain no granuloma-grade tissue, only
    that these sections do not. If CFU or gross pathology exists for these
    animals, agreement between it and the AKOYA burden would do more for
    statement 1 than any additional statistic here.

    The vendor IDO1 binary call is not reproducible as an intensity threshold
    and the raw imagery was never delivered, so nothing the Akoya classifier
    did can be audited. Every IDO1-positive count in this analysis inherits
    that.

    Three animals per arm. The exact p is floored at 0.10, Benjamini-Hochberg
    on a family of exact p-values is vacuous by construction because the
    smallest achievable q is 0.10 times the family size, and the model p is not
    calibrated on six clusters. Multiplicity here is handled by
    pre-specification and coherence across measurements, not by correction.
    That belongs in Methods as a stated position rather than in a rebuttal.
""")

banner("DONE")
sys.stdout = _tee.terminal
_tee.close()
