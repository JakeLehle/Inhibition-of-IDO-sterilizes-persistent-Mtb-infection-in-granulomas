#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - LEAVE-ONE-OUT DETECTION AUDIT
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 04g of the AKOYA analysis series. REVISION 1.
READ-ONLY. It opens every structures directory and writes only to its own.

WHY THIS SCRIPT EXISTS

    Foci are detected on POOLED MYELOID DENSITY, and the pool is

        MYELOID_FOR_DETECTION = [CD68+IDO1+ Macrophages,
                                 CD68+IDO1- Macrophages,
                                 CD163+ Macrophages,
                                 Neutrophils]

    So the share of any of those four INSIDE a focus is partly a restatement of
    the rule that found the focus. "Treated foci contain fewer IDO1-positive
    macrophages" is not safe as a claim if dropping IDO1-positive macrophages
    from the pool is what stops a focus being found at all. Stating the cell
    count does not address this. Circularity is not a small-sample problem.

    Script 04's LEAVE_ONE_OUT exists to test it and redirects OUT_DIR to
    structures_rev5_no_<phenotype> so a sensitivity run cannot touch the
    primary. NOTHING HAS EVER READ THOSE DIRECTORIES. This script is the
    consumer.

    Not to be confused with script 11's 95_per_phenotype_leave_one_out.csv,
    which is leave-one-ANIMAL-out on the radial outcome. Different check,
    different unit, different question.

WHAT IT ASKS

    For each phenotype X dropped from the detection pool: are these the same
    foci? The statistic is cell-level, so it needs no masks and no re-detection.

        core-cell Jaccard   |core_primary AND core_LOO| / |core_primary OR core_LOO|
                            per section. The direct question.
        focus count         per animal. A detector finding a different number of
                            foci is a different detector.
        burden fraction     per animal. Finding 1's headline component.
        peak density        per animal. Finding 1's other surviving component.
        composition share   per animal, every REMAINING phenotype inside cores.
                            The circularity question for the composition panels.
        myeloid density     per cell, median per animal. See below.

    Everything is reported PER ANIMAL, never pooled across an arm, because the
    most likely outcome here is an arm asymmetry and pooling would hide it.

THE SECOND QUESTION, WHICH IS MORE INTERESTING THAN THE FIRST

    IDO1-positive macrophages are 51 to 82 percent of the CD68 lineage in
    untreated sections and 1 to 11 percent in treated. So removing them from the
    pool should cut untreated pooled density far more than treated. If it does,
    that is a number for how much of finding 1's density gap is finding 2, which
    the 1 October note flagged as "S2 is partly S1 restated" without
    quantifying. The myeloid_density column gives it directly.

    That is not circularity. It is two findings sharing a measurement, and it
    belongs in Methods either way.

NO GATES, ONLY LABELS

    This script prints the distributions a threshold would be set from and
    labels each dropped phenotype. It hard-codes no cutoff for what counts as
    "detection survived". That is the standing evidence-funnel rule: run
    permissive, read the numbers, then choose the gate on the evidence rather
    than inheriting one. REPORTABILITY_* below are LABEL boundaries used for
    printing, not filters, and nothing is dropped on their account.

HOW TO READ THE VERDICT

    A phenotype whose removal leaves the Jaccard high and the focus count
    unchanged is NOT driving detection, so its share inside foci is reportable
    and the composition figure may test it.

    A phenotype whose removal collapses detection IS driving it, so its share
    inside foci is circular. The figure may describe it and may not test it.

    A phenotype whose removal changes detection in one arm and not the other is
    neither, and is the case to expect for IDO1-positive macrophages. The
    detector is then resting on different cells in the two arms while applying
    one threshold. Report it rather than resolve it.

INPUTS
    structures_rev5/tables/33b_burden_summary.csv
    structures_rev5/tables/35_foci_structures_relative.csv
    structures_rev5/cell_assignments/<section>_cell_structures.csv
    structures_rev5_no_<phenotype>/...  the same three, for each of the four

OUTPUTS
    leave_one_out_audit/tables/   120 .. 125
    leave_one_out_audit/figures/  F105 .. F107

USAGE
    Run the four sensitivity detections first, one at a time, in
    AKOYA_04_Structures_rev6.py Cell 1:

        LEAVE_ONE_OUT = "CD68+IDO1+ Macrophages"
        LEAVE_ONE_OUT = "CD68+IDO1- Macrophages"
        LEAVE_ONE_OUT = "CD163+ Macrophages"
        LEAVE_ONE_OUT = "Neutrophils"

    Use REV 6, not rev 5, so the siblings carry the same columns as the
    primary. Set LEAVE_ONE_OUT back to None afterwards: leaving it set and
    rerunning is how the primary directory would be overwritten by a
    sensitivity run.

    Then:
        conda activate sc_pre
        python AKOYA_04g_Leave_One_Out_Audit.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================

BASE = "/master/jlehle/WORKING/AKOYA"
PRIMARY_DIR = f"{BASE}/structures_rev5"
OUT_DIR = f"{BASE}/leave_one_out_audit"

IDO1_POS = "CD68+IDO1+ Macrophages"
IDO1_NEG = "CD68+IDO1- Macrophages"

# Must match MYELOID_FOR_DETECTION in script 04. Checked against the sibling
# directories' own is_detection_pool column, so a mismatch is caught rather
# than assumed away.
MYELOID_FOR_DETECTION = [IDO1_POS, IDO1_NEG, "CD163+ Macrophages", "Neutrophils"]

# directory suffix rule, identical to script 04's
def loo_dirname(pheno):
    return f"{PRIMARY_DIR}_no_{pheno.replace(' ', '_').replace('+', '')}"


CONDITION_ORDER = ["D1MT", "Untreated"]
CONDITION_COLORS = {"D1MT": "#2C7FB8", "Untreated": "#D95F02"}
REFERENCE_ARM = "Untreated"

STRUCT_COL = "focus_id"
CORE_REGION = "core"

# Canonical palette, same hex codes as script 03. These belong in
# akoya_arm_stats.py so there is one copy; until they move, this is a verbatim
# duplicate and must not be edited independently.
PHENOTYPE_COLORS = {
    "CD68+IDO1+ Macrophages": "#B2182B",
    "CD68+IDO1- Macrophages": "#EF8A62",
    "CD163+ Macrophages": "#FDBE85",
    "Neutrophils": "#7B3294",
    "Helper T cells": "#1B7837",
    "CD4- T cells": "#7FBC41",
    "Tregs": "#00441B",
    "B cells": "#2166AC",
    "Plasma cells": "#67A9CF",
    "Endothelial cells": "#E7298A",
    "Epithelial/Tumor cells": "#66C2A5",
    "Other": "#A6761D",
}
PHENOTYPE_ORDER = list(PHENOTYPE_COLORS)

# ---- cell matching between runs --------------------------------------------
# Script 04 builds its per-cell frame with out = d.copy(), and detection labels
# rows rather than filtering them, so row order is stable across runs. That is
# an assumption and it is checked rather than trusted: cells are keyed on
# rounded coordinates and the match rate is printed. Below MIN_MATCH_RATE the
# comparison is refused for that section, because a Jaccard on a partial match
# is not a Jaccard.
COORD_DECIMALS = 2
MIN_MATCH_RATE = 0.995

# ---- LABEL boundaries, not filters -----------------------------------------
# Used only to print a word next to a number. Nothing is excluded on these.
LABEL_JACCARD_HIGH = 0.90      # "same foci"
LABEL_JACCARD_LOW = 0.50       # "different foci"
LABEL_FOCUS_COUNT_TOL = 0.20   # fractional change treated as unchanged
LABEL_SHARE_SHIFT_PP = 5.0     # percentage points, "composition moved"

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
MID_COLOR = "#E6AB02"


# %% Cell 2 - imports, style, helpers
# =============================================================================

import os
import sys
import glob
import gc
from datetime import datetime

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

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


def short_pheno(p):
    return (p.replace("CD68+IDO1+ Macrophages", "IDO1+ Mac")
             .replace("CD68+IDO1- Macrophages", "IDO1- Mac")
             .replace("CD163+ Macrophages", "CD163+ Mac")
             .replace(" cells", ""))


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
    # lineterminator is load-bearing: the default \r\n fails silently in bash
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


def cell_key(d):
    """Stable per-cell key across runs, from rounded coordinates."""
    return (d["x"].round(COORD_DECIMALS).astype(str) + "_"
            + d["y"].round(COORD_DECIMALS).astype(str))


def jaccard(a, b):
    a, b = set(a), set(b)
    if not a and not b:
        return np.nan
    u = len(a | b)
    return float(len(a & b) / u) if u else np.nan


def read_cells(directory):
    """Per-section cell assignment frames from one structures directory."""
    cdir = os.path.join(directory, "cell_assignments")
    paths = sorted(glob.glob(os.path.join(cdir, "*_cell_structures.csv")))
    out = {}
    for p in paths:
        sid = os.path.basename(p).replace("_cell_structures.csv", "")
        try:
            d = pd.read_csv(p, low_memory=False)
        except Exception as e:
            print(f"    ERROR reading {sid} in {os.path.basename(directory)}: "
                  f"{e}. Skipping this section.")
            continue
        need = ["x", "y", "pheno", "condition", "region", STRUCT_COL]
        miss = [c for c in need if c not in d.columns]
        if miss:
            print(f"    ERROR: {sid} in {os.path.basename(directory)} lacks "
                  f"{miss}. Skipping this section.")
            continue
        out[sid] = d
    return out


def read_table(directory, name):
    p = os.path.join(directory, "tables", name)
    if not os.path.exists(p):
        print(f"    WARNING: {p} not found.")
        return None
    try:
        return pd.read_csv(p)
    except Exception as e:
        print(f"    WARNING: could not read {p}: {e}")
        return None


_tee = Tee(os.path.join(TAB_DIR, "00_leave_one_out_audit_report.txt"))
sys.stdout = _tee

banner("AKOYA LEAVE-ONE-OUT DETECTION AUDIT (script 04g rev 1)")
print(f"Run time   : {datetime.now().isoformat(timespec='seconds')}")
print(f"Primary    : {PRIMARY_DIR}")
print(f"Output     : {OUT_DIR}")
print(f"Pool       : {MYELOID_FOR_DETECTION}")
print("""
THE QUESTION
    Foci are detected on pooled myeloid density, so the share of a
    detection-pool phenotype inside a focus is partly a restatement of the
    detection rule. For each of the four, this asks whether dropping it from
    the pool changes which cells end up inside a focus core.

    The statistic is a cell-level Jaccard between the primary core set and the
    leave-one-out core set, per section. No masks, no re-detection.

NO GATES
    The distributions a threshold would be set from are printed, and each
    dropped phenotype gets a label. Nothing is filtered. Choose the gate on the
    evidence afterwards.
""")


# %% Cell 3 - load the primary and every sibling, with a provenance guard
# =============================================================================

banner("LOADING")

sub("Primary")
primary = read_cells(PRIMARY_DIR)
if not primary:
    print(f"ERROR: no cell assignments in {PRIMARY_DIR}/cell_assignments.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

cond_rank = {c: i for i, c in enumerate(CONDITION_ORDER)}
SAMPLE_ORDER = sorted(primary, key=lambda s: (
    cond_rank.get(primary[s]["condition"].iloc[0], 9), s))
COND_OF = {s: primary[s]["condition"].iloc[0] for s in SAMPLE_ORDER}
MARKER_OF = {s: POINT_MARKERS[i % len(POINT_MARKERS)]
             for i, s in enumerate(SAMPLE_ORDER)}

for s in SAMPLE_ORDER:
    d = primary[s]
    ncore = int((d["region"] == CORE_REGION).sum())
    print(f"    {s:<12} {COND_OF[s]:<10} {len(d):>8,} cells   "
          f"{ncore:>7,} core   {int(d[STRUCT_COL].nunique()):>3} structures")

print("\n    EFFECTIVE INDEPENDENT UNITS")
print("    Treatment was assigned to ANIMALS. Six units, three per arm,")
print("    whatever the cell counts say. This script makes no arm-level test;")
print("    it reports per animal and leaves the inference to the owning script.")

# ---- the siblings, and the guard -------------------------------------------
sub("Leave-one-out runs, and whether each is the run it claims to be")
print("    Script 04 writes is_detection_pool per cell. A sibling directory")
print("    that dropped X must show is_detection_pool False for every X cell")
print("    and True for the other three. If it does not, the directory is not")
print("    the run its name claims and it is refused rather than compared.\n")

loo = {}
for X in MYELOID_FOR_DETECTION:
    dirname = loo_dirname(X)
    if not os.path.isdir(dirname):
        print(f"    MISSING   {short_pheno(X):<14} no directory at {dirname}")
        print(f"              run script 04 rev6 with LEAVE_ONE_OUT = \"{X}\"")
        continue
    cells_X = read_cells(dirname)
    if not cells_X:
        print(f"    EMPTY     {short_pheno(X):<14} {dirname} has no cell "
              f"assignments")
        continue

    # provenance guard
    ok, detail = True, []
    for s in sorted(cells_X):
        d = cells_X[s]
        if "is_detection_pool" not in d.columns:
            detail.append("no is_detection_pool column")
            ok = False
            break
        for p in MYELOID_FOR_DETECTION:
            m = d["pheno"] == p
            if not m.any():
                continue
            frac_in = float(d.loc[m, "is_detection_pool"].mean())
            should_be_in = (p != X)
            if should_be_in and frac_in < 0.99:
                detail.append(f"{s}: {short_pheno(p)} should be IN the pool "
                              f"but is {100*frac_in:.0f}% in")
                ok = False
            if (not should_be_in) and frac_in > 0.01:
                detail.append(f"{s}: {short_pheno(p)} should be OUT of the "
                              f"pool but is {100*frac_in:.0f}% in")
                ok = False
    if not ok:
        print(f"    REFUSED   {short_pheno(X):<14} provenance guard failed")
        for t in detail[:4]:
            print(f"              {t}")
        print("              This directory is not the leave-one-out run its "
              "name claims.")
        continue

    loo[X] = cells_X
    print(f"    OK        {short_pheno(X):<14} {len(cells_X)} sections, "
          f"pool membership verified")

if not loo:
    print("\nERROR: no usable leave-one-out directories. Run the four "
          "sensitivity detections first; the header has the exact settings.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

print(f"\n    {len(loo)} of {len(MYELOID_FOR_DETECTION)} leave-one-out runs "
      f"available for comparison.")
if len(loo) < len(MYELOID_FOR_DETECTION):
    print("    The audit continues on what is present. A phenotype with no")
    print("    run gets no verdict, and its composition share stays circular")
    print("    by default rather than by evidence.")


# %% Cell 4 - A: are these the same foci
# =============================================================================

banner("A - CORE-CELL JACCARD, PRIMARY AGAINST EACH LEAVE-ONE-OUT")

print("    For each section: the set of cells assigned to a focus CORE in the")
print("    primary run, against the same set in the leave-one-out run.")
print("    1.00 means identical foci. Cells are keyed on rounded coordinates")
print("    and the match rate is checked before any Jaccard is computed.\n")

jac_rows = []
for X, cells_X in loo.items():
    for s in SAMPLE_ORDER:
        if s not in cells_X:
            print(f"    WARNING: {short_pheno(X)} has no section {s}. Skipped.")
            continue
        dp, dl = primary[s], cells_X[s]
        kp, kl = cell_key(dp), cell_key(dl)
        setp, setl = set(kp), set(kl)
        inter = len(setp & setl)
        rate = inter / max(len(setp), 1)
        if rate < MIN_MATCH_RATE:
            print(f"    REFUSED  {short_pheno(X):<14} {s}: only "
                  f"{100*rate:.2f}% of primary cells found in the "
                  f"leave-one-out run. A Jaccard on a partial match is not a "
                  f"Jaccard, so this section is not compared.")
            continue
        core_p = set(kp[dp["region"] == CORE_REGION])
        core_l = set(kl[dl["region"] == CORE_REGION])
        jac_rows.append({
            "dropped": X, "sample_id": s, "animal_id": short_label(s),
            "condition": COND_OF[s],
            "cell_match_rate": rate,
            "n_core_primary": len(core_p), "n_core_loo": len(core_l),
            "n_core_shared": len(core_p & core_l),
            "core_jaccard": jaccard(core_p, core_l),
            "core_recovered_frac": (len(core_p & core_l) / len(core_p)
                                    if core_p else np.nan),
            "n_struct_primary": int(dp[STRUCT_COL].nunique()),
            "n_struct_loo": int(dl[STRUCT_COL].nunique()),
        })

jac = pd.DataFrame(jac_rows)
if not len(jac):
    print("ERROR: nothing comparable. Stopping.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)
write_csv(jac, "120_core_cell_jaccard.csv")

sub("Core-cell Jaccard per section")
print(f"    {'dropped':<16}" + "".join(f"{short_label(s):>11}"
                                       for s in SAMPLE_ORDER)
      + f"{'median':>10}{'label':>22}")
print("    " + "-" * (16 + 11 * len(SAMPLE_ORDER) + 32))
for X in MYELOID_FOR_DETECTION:
    g = jac.loc[jac["dropped"] == X]
    if not len(g):
        print(f"    {short_pheno(X):<16}" + f"{'no run':>{11*len(SAMPLE_ORDER)}}")
        continue
    vals = []
    for s in SAMPLE_ORDER:
        r = g.loc[g["sample_id"] == s, "core_jaccard"]
        vals.append(float(r.iloc[0]) if len(r) else np.nan)
    med = float(np.nanmedian(vals))
    lab = ("same foci" if med >= LABEL_JACCARD_HIGH
           else "different foci" if med < LABEL_JACCARD_LOW
           else "partly changed")
    print(f"    {short_pheno(X):<16}"
          + "".join(f"{v:>11.3f}" if np.isfinite(v) else f"{'na':>11}"
                    for v in vals)
          + f"{med:>10.3f}{lab:>22}")

sub("Jaccard by arm, because an asymmetry here is the thing to look for")
print(f"    {'dropped':<16}{'D1MT median':>14}{'Untreated median':>18}"
      f"{'gap':>8}{'reading':>34}")
print("    " + "-" * 90)
arm_rows = []
for X in MYELOID_FOR_DETECTION:
    g = jac.loc[jac["dropped"] == X]
    if not len(g):
        continue
    t = float(np.nanmedian(g.loc[g["condition"] == "D1MT", "core_jaccard"]))
    u = float(np.nanmedian(g.loc[g["condition"] != "D1MT", "core_jaccard"]))
    gap = t - u
    if min(t, u) >= LABEL_JACCARD_HIGH:
        reading = "not driving detection"
    elif max(t, u) < LABEL_JACCARD_LOW:
        reading = "driving detection, both arms"
    elif abs(gap) >= 0.20:
        reading = "ARM ASYMMETRIC, see below"
    else:
        reading = "partly driving, both arms"
    arm_rows.append({"dropped": X, "jaccard_D1MT": t,
                     "jaccard_Untreated": u, "gap": gap, "reading": reading})
    print(f"    {short_pheno(X):<16}{t:>14.3f}{u:>18.3f}{gap:>+8.3f}"
          f"{reading:>34}")
write_csv(pd.DataFrame(arm_rows), "121_jaccard_by_arm.csv")

print("\n    An ARM ASYMMETRIC row is not circularity and is not cleanliness.")
print("    It means the detector rests on different cells in the two arms")
print("    while applying one threshold. That belongs in Methods as a stated")
print("    property of the detector, whichever way the composition claim goes.")


# %% Cell 5 - B: finding 1's components under each leave-one-out
# =============================================================================

banner("B - FOCUS COUNT, BURDEN AND PEAK DENSITY UNDER EACH LEAVE-ONE-OUT")

print("    Finding 1 rests on burden and peak density, and both are computed")
print("    on the pooled density the leave-one-out changes. So this is not")
print("    only a circularity check, it is how much of finding 1 depends on")
print("    each phenotype being in the pool.\n")

bur_p = read_table(PRIMARY_DIR, "33b_burden_summary.csv")
comp_rows = []
for X, cells_X in loo.items():
    bur_l = read_table(loo_dirname(X), "33b_burden_summary.csv")
    for s in SAMPLE_ORDER:
        if s not in cells_X:
            continue
        dp, dl = primary[s], cells_X[s]
        rec = {"dropped": X, "sample_id": s, "animal_id": short_label(s),
               "condition": COND_OF[s],
               "n_struct_primary": int(dp[STRUCT_COL].nunique()),
               "n_struct_loo": int(dl[STRUCT_COL].nunique()),
               "n_core_primary": int((dp["region"] == CORE_REGION).sum()),
               "n_core_loo": int((dl["region"] == CORE_REGION).sum())}
        # per-cell myeloid density is the quantity the pool change acts on
        for tag, d in (("primary", dp), ("loo", dl)):
            if "myeloid_density" in d.columns:
                v = pd.to_numeric(d["myeloid_density"], errors="coerce")
                rec[f"median_density_{tag}"] = float(np.nanmedian(v))
                rec[f"p99_density_{tag}"] = float(np.nanpercentile(
                    v.dropna(), 99)) if v.notna().any() else np.nan
            else:
                rec[f"median_density_{tag}"] = np.nan
                rec[f"p99_density_{tag}"] = np.nan
        for tab, tag in ((bur_p, "primary"), (bur_l, "loo")):
            val = np.nan
            if tab is not None and "sample_id" in tab.columns:
                r = tab.loc[tab["sample_id"] == s]
                for cand in ("burden_pct_of_tissue", "burden_fraction_pct",
                             "pct_tissue_at_granuloma_density", "burden_pct"):
                    if len(r) and cand in tab.columns:
                        val = float(r[cand].iloc[0]); break
            rec[f"burden_pct_{tag}"] = val
        comp_rows.append(rec)

comp = pd.DataFrame(comp_rows)
if len(comp):
    comp["density_ratio"] = (comp["p99_density_loo"]
                             / comp["p99_density_primary"])
    comp["struct_change_frac"] = ((comp["n_struct_loo"]
                                   - comp["n_struct_primary"])
                                  / comp["n_struct_primary"].replace(0, np.nan))
    write_csv(comp, "122_finding1_components_under_loo.csv")

    sub("Structures detected, primary against leave-one-out")
    print(f"    {'dropped':<16}" + "".join(f"{short_label(s):>13}"
                                           for s in SAMPLE_ORDER))
    print("    " + "-" * (16 + 13 * len(SAMPLE_ORDER)))
    for X in MYELOID_FOR_DETECTION:
        g = comp.loc[comp["dropped"] == X]
        if not len(g):
            continue
        cells_txt = []
        for s in SAMPLE_ORDER:
            r = g.loc[g["sample_id"] == s]
            cells_txt.append(f"{int(r['n_struct_primary'].iloc[0])}->"
                             f"{int(r['n_struct_loo'].iloc[0])}"
                             if len(r) else "na")
        print(f"    {short_pheno(X):<16}"
              + "".join(f"{t:>13}" for t in cells_txt))

    sub("99th percentile myeloid density, leave-one-out over primary")
    print("    This is the number for how much of the pooled density each")
    print("    phenotype contributes, per animal. A value near 1.00 means the")
    print("    density the detector sees barely changes without it.\n")
    print(f"    {'dropped':<16}" + "".join(f"{short_label(s):>11}"
                                           for s in SAMPLE_ORDER)
          + f"{'D1MT med':>11}{'Untr med':>11}")
    print("    " + "-" * (16 + 11 * len(SAMPLE_ORDER) + 22))
    for X in MYELOID_FOR_DETECTION:
        g = comp.loc[comp["dropped"] == X]
        if not len(g):
            continue
        vals = []
        for s in SAMPLE_ORDER:
            r = g.loc[g["sample_id"] == s, "density_ratio"]
            vals.append(float(r.iloc[0]) if len(r) else np.nan)
        tm = float(np.nanmedian([v for v, s in zip(vals, SAMPLE_ORDER)
                                 if COND_OF[s] == "D1MT"]))
        um = float(np.nanmedian([v for v, s in zip(vals, SAMPLE_ORDER)
                                 if COND_OF[s] != "D1MT"]))
        print(f"    {short_pheno(X):<16}"
              + "".join(f"{v:>11.3f}" if np.isfinite(v) else f"{'na':>11}"
                        for v in vals)
              + f"{tm:>11.3f}{um:>11.3f}")

    print("\n    HOW MUCH OF FINDING 1 IS FINDING 2. Read the IDO1+ Mac row.")
    print("    IDO1-positive macrophages are 51 to 82 percent of the CD68")
    print("    lineage in untreated sections and 1 to 11 percent in treated,")
    print("    so if the untreated density ratio is far below the treated one,")
    print("    the density gap finding 1 reports is substantially the IDO1+")
    print("    compartment that finding 2 reports. That is two findings")
    print("    sharing a measurement rather than either being wrong, and the")
    print("    1 October note flagged it without a number. This is the number.")


# %% Cell 6 - C: does the composition of the remaining phenotypes move
# =============================================================================

banner("C - COMPOSITION SHIFT INSIDE CORES, REMAINING PHENOTYPES ONLY")

print("    For every phenotype EXCEPT the dropped one: its share of core cells")
print("    in the primary run against the leave-one-out run, per animal. A")
print("    share that barely moves is robust to the detection pool.")
print("    The dropped phenotype's own share is not reported, because it is")
print("    not a comparison: it was removed from the thing being measured.\n")

shift_rows = []
for X, cells_X in loo.items():
    for s in SAMPLE_ORDER:
        if s not in cells_X:
            continue
        dp, dl = primary[s], cells_X[s]
        cp = dp.loc[dp["region"] == CORE_REGION, "pheno"].value_counts(
            normalize=True) * 100.0
        cl = dl.loc[dl["region"] == CORE_REGION, "pheno"].value_counts(
            normalize=True) * 100.0
        for p in PHENOTYPE_ORDER:
            if p == X:
                continue
            a, b = float(cp.get(p, 0.0)), float(cl.get(p, 0.0))
            shift_rows.append({
                "dropped": X, "sample_id": s, "animal_id": short_label(s),
                "condition": COND_OF[s], "pheno": p,
                "share_primary_pct": a, "share_loo_pct": b,
                "shift_pp": b - a,
                "is_detection_pool": p in MYELOID_FOR_DETECTION})

shift = pd.DataFrame(shift_rows)
if len(shift):
    write_csv(shift, "123_composition_shift.csv")
    sub("Largest absolute share shift, percentage points, per dropped phenotype")
    print(f"    {'dropped':<16}{'worst pheno':<18}{'worst shift':>13}"
          f"{'median |shift|':>16}{'pool members moved':>21}")
    print("    " + "-" * 86)
    summ = []
    for X in MYELOID_FOR_DETECTION:
        g = shift.loc[shift["dropped"] == X]
        if not len(g):
            continue
        per_ph = g.groupby("pheno")["shift_pp"].apply(
            lambda v: float(np.nanmax(np.abs(v))))
        worst_ph = per_ph.idxmax()
        worst = float(per_ph.max())
        medabs = float(np.nanmedian(np.abs(g["shift_pp"])))
        moved = int((g.loc[g["is_detection_pool"]]
                     .groupby("pheno")["shift_pp"]
                     .apply(lambda v: float(np.nanmax(np.abs(v))))
                     >= LABEL_SHARE_SHIFT_PP).sum())
        summ.append({"dropped": X, "worst_pheno": worst_ph,
                     "worst_shift_pp": worst, "median_abs_shift_pp": medabs,
                     "pool_members_moved": moved})
        print(f"    {short_pheno(X):<16}{short_pheno(worst_ph):<18}"
              f"{worst:>13.2f}{medabs:>16.2f}{moved:>21}")
    write_csv(pd.DataFrame(summ), "124_composition_shift_summary.csv")

    print("\n    A remaining phenotype whose share moves a lot when some OTHER")
    print("    phenotype leaves the pool is not independent of the detection")
    print("    rule either, even though it was never in the pool. Watch the")
    print("    lymphoid rows: those are the ones the composition figure wants")
    print("    to make claims about.")


# %% Cell 7 - D: the verdict, and the figures
# =============================================================================

banner("D - VERDICT PER DROPPED PHENOTYPE")

print("    LABELS, NOT GATES. Each phenotype gets a reading and a licence for")
print("    what the composition figure may do with it. Nothing was filtered to")
print("    produce this table.\n")

verdict = []
for X in MYELOID_FOR_DETECTION:
    gj = jac.loc[jac["dropped"] == X]
    if not len(gj):
        verdict.append({"dropped": X, "status": "NO RUN",
                        "core_jaccard_median": np.nan,
                        "jaccard_D1MT": np.nan, "jaccard_Untreated": np.nan,
                        "licence": "circular by default, no evidence either way"})
        continue
    med = float(np.nanmedian(gj["core_jaccard"]))
    t = float(np.nanmedian(gj.loc[gj["condition"] == "D1MT", "core_jaccard"]))
    u = float(np.nanmedian(gj.loc[gj["condition"] != "D1MT", "core_jaccard"]))
    gc_ = comp.loc[comp["dropped"] == X] if len(comp) else pd.DataFrame()
    struct_stable = True
    if len(gc_):
        struct_stable = bool(
            np.nanmax(np.abs(gc_["struct_change_frac"].to_numpy()))
            <= LABEL_FOCUS_COUNT_TOL)
    if med >= LABEL_JACCARD_HIGH and struct_stable:
        status = "NOT DRIVING DETECTION"
        lic = "share inside foci is reportable and may be tested"
    elif med < LABEL_JACCARD_LOW:
        status = "DRIVING DETECTION"
        lic = "share inside foci is CIRCULAR, describe only, do not test"
    elif abs(t - u) >= 0.20:
        status = "ARM ASYMMETRIC"
        lic = "describe only, and state the asymmetry in Methods"
    else:
        status = "PARTLY DRIVING"
        lic = "describe only, testing needs a stated justification"
    verdict.append({"dropped": X, "status": status,
                    "core_jaccard_median": med, "jaccard_D1MT": t,
                    "jaccard_Untreated": u,
                    "structures_stable": struct_stable, "licence": lic})

vd = pd.DataFrame(verdict)
write_csv(vd, "125_circularity_verdict.csv")
print(f"    {'dropped':<16}{'status':<24}{'J med':>8}  {'licence':<52}")
print("    " + "-" * 100)
for _, r in vd.iterrows():
    j = r["core_jaccard_median"]
    print(f"    {short_pheno(r['dropped']):<16}{r['status']:<24}"
          + (f"{j:>8.3f}" if np.isfinite(j) else f"{'na':>8}")
          + f"  {r['licence']:<52}")

print("\n    WHAT SCRIPT 12 DOES WITH THIS. The composition figure reads table")
print("    125. A phenotype labelled NOT DRIVING DETECTION may carry a test on")
print("    its share. Everything else is drawn with the circularity bracket and")
print("    the caption says descriptive. The lymphoid populations are not in")
print("    the pool at all and are unaffected by any of this.")

# ---- F105: Jaccard ---------------------------------------------------------
if len(jac):
    fig, ax = plt.subplots(figsize=(18, 13))
    xs = {X: i for i, X in enumerate(MYELOID_FOR_DETECTION)}
    # ONE generator for the whole panel. Building np.random.default_rng(0)
    # inside the comprehension re-seeds on every point, so every point gets the
    # identical offset and the six sections stack exactly on top of each other.
    # The jitter is cosmetic but a jitter that does not jitter is worse than
    # none, because it looks like the sections agree perfectly.
    jrng = np.random.default_rng(0)
    for s in SAMPLE_ORDER:
        g = jac.loc[jac["sample_id"] == s]
        if not len(g):
            continue
        ax.scatter([xs[X] + jrng.uniform(-0.12, 0.12) for X in g["dropped"]],
                   g["core_jaccard"], s=420, marker=MARKER_OF[s],
                   color=CONDITION_COLORS.get(COND_OF[s], "#777777"),
                   edgecolor="#000000", linewidth=2, alpha=0.85, zorder=3,
                   label=f"{short_label(s)} ({COND_OF[s]})")
    ax.axhline(LABEL_JACCARD_HIGH, color=OK_COLOR, linestyle="--",
               linewidth=3, zorder=1)
    ax.axhline(LABEL_JACCARD_LOW, color=FLAG_COLOR, linestyle=":",
               linewidth=3, zorder=1)
    ax.set_xticks(list(xs.values()))
    ax.set_xticklabels([short_pheno(X) for X in xs], rotation=20, ha="right")
    ax.set_ylabel("core-cell Jaccard\nprimary vs leave-one-out")
    ax.set_ylim(0, 1.02)
    ax.set_title("Does dropping a phenotype from the detection pool\n"
                 "change which cells sit inside a focus core?",
                 fontsize=FONT_SIZE_TITLE)
    # Labels on the RIGHT edge. On the left they sit on top of the leftmost
    # phenotype's points, which is exactly where the most interesting values
    # tend to be.
    ax.text(0.995, LABEL_JACCARD_HIGH + 0.015, "same foci", color=OK_COLOR,
            fontsize=FONT_SIZE_ANNOT, ha="right",
            transform=ax.get_yaxis_transform())
    ax.text(0.995, LABEL_JACCARD_LOW + 0.015, "different foci",
            color=FLAG_COLOR, fontsize=FONT_SIZE_ANNOT, ha="right",
            transform=ax.get_yaxis_transform())
    # Legend OUTSIDE the axes. Inside, it collides with any animal whose
    # Jaccard lands low, and a low Jaccard is the result worth seeing.
    ax.legend(frameon=False, fontsize=FONT_SIZE_LEGEND - 10, ncol=3,
              loc="upper center", bbox_to_anchor=(0.5, -0.22))
    style_axes(ax)
    save_fig(fig, "F105_core_cell_jaccard")

# ---- F106: structures and density ------------------------------------------
if len(comp):
    fig, axes = plt.subplots(1, 2, figsize=(26, 12))
    ax = axes[0]
    for X in MYELOID_FOR_DETECTION:
        g = comp.loc[comp["dropped"] == X]
        if not len(g):
            continue
        for _, r in g.iterrows():
            ax.plot([0, 1], [r["n_struct_primary"], r["n_struct_loo"]],
                    color=CONDITION_COLORS.get(r["condition"], "#777777"),
                    linewidth=3, alpha=0.75, zorder=2)
            ax.scatter([0, 1], [r["n_struct_primary"], r["n_struct_loo"]],
                       s=260, color=CONDITION_COLORS.get(r["condition"],
                                                         "#777777"),
                       edgecolor="#000000", linewidth=2, zorder=3)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["primary", "leave-one-out"])
    ax.set_xlim(-0.3, 1.3)
    ax.set_yscale("symlog", linthresh=1)
    ax.set_ylabel("structures detected")
    ax.set_title("Focus count", fontsize=FONT_SIZE_TITLE)
    style_axes(ax)

    ax = axes[1]
    xs = {X: i for i, X in enumerate(MYELOID_FOR_DETECTION)}
    for s in SAMPLE_ORDER:
        g = comp.loc[comp["sample_id"] == s]
        if not len(g):
            continue
        ax.scatter([xs[X] for X in g["dropped"]], g["density_ratio"],
                   s=420, marker=MARKER_OF[s],
                   color=CONDITION_COLORS.get(COND_OF[s], "#777777"),
                   edgecolor="#000000", linewidth=2, alpha=0.85, zorder=3)
    ax.axhline(1.0, color="#000000", linewidth=3, zorder=1)
    ax.set_xticks(list(xs.values()))
    ax.set_xticklabels([short_pheno(X) for X in xs], rotation=20, ha="right")
    ax.set_ylabel("p99 myeloid density\nleave-one-out / primary")
    ax.set_title("How much pooled density each phenotype carries",
                 fontsize=FONT_SIZE_TITLE)
    style_axes(ax)
    fig.legend(handles=[Patch(facecolor=CONDITION_COLORS[c], label=c)
                        for c in CONDITION_ORDER],
               frameon=False, fontsize=FONT_SIZE_LEGEND,
               loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.04))
    save_fig(fig, "F106_structures_and_density")

# ---- F107: composition shift ----------------------------------------------
if len(shift):
    present = [X for X in MYELOID_FOR_DETECTION
               if (shift["dropped"] == X).any()]
    fig, axes = plt.subplots(1, max(len(present), 1),
                             figsize=(9 * max(len(present), 1), 13),
                             squeeze=False)
    for j, X in enumerate(present):
        ax = axes[0][j]
        g = shift.loc[shift["dropped"] == X]
        phs = [p for p in PHENOTYPE_ORDER if p != X]
        for i, p in enumerate(phs):
            gg = g.loc[g["pheno"] == p]
            if not len(gg):
                continue
            ax.scatter(gg["shift_pp"], np.full(len(gg), i), s=300,
                       color=PHENOTYPE_COLORS.get(p, "#999999"),
                       edgecolor="#000000", linewidth=1.5, alpha=0.9,
                       zorder=3)
        ax.axvline(0, color="#000000", linewidth=3, zorder=1)
        for b in (-LABEL_SHARE_SHIFT_PP, LABEL_SHARE_SHIFT_PP):
            ax.axvline(b, color=MID_COLOR, linestyle=":", linewidth=2.5,
                       zorder=1)
        ax.set_yticks(range(len(phs)))
        ax.set_yticklabels([short_pheno(p) for p in phs]
                           if j == 0 else [""] * len(phs))
        ax.set_xlabel("share shift, percentage points")
        ax.set_title(f"dropped\n{short_pheno(X)}", fontsize=FONT_SIZE_TITLE - 4)
        style_axes(ax, ygrid=False)
        ax.xaxis.grid(True, color=GRID_COLOR, linewidth=1.5)
        ax.set_axisbelow(True)
    save_fig(fig, "F107_composition_shift")


banner("SUMMARY")
print(f"Primary                   : {PRIMARY_DIR}")
print(f"Leave-one-out runs found  : {len(loo)} of {len(MYELOID_FOR_DETECTION)}")
print(f"Sections compared         : {jac['sample_id'].nunique()}")
print(f"Comparisons               : {len(jac)}")
print(f"Cell match rate, minimum  : {jac['cell_match_rate'].min():.5f}"
      f"   (refusal threshold {MIN_MATCH_RATE})")
print("Gates applied             : NONE. Labels only.")
print("\nRead in this order")
print("  1. The provenance guard in LOADING. A refused directory is not the")
print("     run its name claims and nothing below it means anything for that")
print("     phenotype.")
print("  2. Section A, the Jaccard, by arm. An arm asymmetry is the outcome to")
print("     expect for IDO1+ macrophages and it is not circularity.")
print("  3. Section B, the density ratio. The IDO1+ Mac row is how much of")
print("     finding 1 is finding 2.")
print("  4. Section D, table 125, which is what script 12 reads to decide")
print("     whether a composition panel may carry a test.")
print("\nWhat this script does NOT do")
print("  No arm-level test. It reports per animal and leaves inference to the")
print("  owning script. No re-detection. No change to any structures")
print("  directory. If a number here disagrees with script 04, script 04 owns")
print("  it and this script is wrong.")

banner("DONE")
sys.stdout = _tee.terminal
_tee.close()
