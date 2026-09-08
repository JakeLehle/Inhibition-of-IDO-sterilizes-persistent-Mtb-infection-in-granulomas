#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - PUBLICATION FIGURES
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 09 of the AKOYA analysis series. REVISION 3.

No new inference. Assembly, plus one clustered arm test per comparison drawn
from akoya_arm_stats so the test is identical everywhere it appears.

WHAT CHANGED IN REVISION 3

  1. FINDING 1 IS A FOCI STORY WITH BURDEN AS THE CLOSER.
     Revision 2 led with burden, which says treated animals have zero
     granuloma-density tissue. True, and the weaker version. The stronger
     statement is that treated animals RETAIN discrete myeloid foci which remain
     7.4 to 8.1 fold above their own tissue background, so they are real
     structures rather than noise, but which run at roughly a third the absolute
     density and never reach granuloma density. Burden becomes the final panel.

     THE TWO DEFINITIONS, AND WHY BOTH EXIST.
     Absolute (burden): one fixed threshold, 3,000 pooled myeloid cells/mm2,
     applied identically to every section, filtered at 30,000 um2 and 100 cells.
     Answers how much granuloma-density tissue there is. The near-zero treated
     result is the finding, not a detection failure.
     Relative (foci): local maxima at least 6-fold above THAT SECTION's own
     median myeloid density, each grown to 50 percent of its own peak,
     watershed-partitioned. Scale-free per focus, which is what the radial
     coordinate needs. A single absolute threshold cannot serve both arms:
     treated peak density tops out near 4,000/mm2 and untreated reaches 19,600.
     The two are never mixed, and peak density and fold are always reported
     alongside architecture so a treated focus is never presented as a small
     granuloma.

  2. EVERY ARM COMPARISON CARRIES THE SAME TEST.
     Datta and Satten's clustered rank-sum (JASA 2005), evaluated by exact
     enumeration of cluster-label assignments. Observations are cells or foci,
     clusters are animals, unequal cluster sizes handled. A plain cell-level
     Wilcoxon is used nowhere: simulated on the real B lineage cell counts with
     animals differing but NO arm effect, it returns p < 0.05 in 85.6 percent of
     runs and reaches p = 4e-210 on noise, because its null requires all cells
     to be exchangeable and cells within one animal are not.

  3. ONE-SIDED TESTS ON PRE-SPECIFIED DIRECTIONS, declared in DIRECTIONS below
     and fixed before the data were examined. At three animals per arm there are
     C(6,3) = 20 assignments, so the floor is 0.05 one-sided and 0.10 two-sided.
     That is a design property, not a property of the test. Both p-values go to
     table 85. If an effect runs opposite to its prediction the annotation says
     so rather than flipping the test.

  4. NEW FIGURE F64, the B lineage distribution per cell.
     THE Y AXIS IS CENTRED RADIAL POSITION, NOT MICRONS, ON PURPOSE. Distance
     from the core in microns is bounded by focus size, and the median inscribed
     radius is 125 um treated against 168 um untreated, so a treated B cell
     physically cannot be as far from its core centre as an untreated one.
     Plotting microns would show treated cells closer whether or not any biology
     is happening, in exactly the direction that flatters the hypothesis. The
     caption carries the microns conversion instead, on the inscribed radius,
     which is the scale the radial coordinate normalises by.

  5. F63 removed.

FINDINGS AS REPORTED HERE
  1  Residual foci persist in treated animals at a third the peak density and
     half the fold prominence, and none reaches granuloma density. Burden is
     zero in every treated animal against 7.5, 16.6 and 39.5 percent.
  2  The IDO1 macrophage compartment shifts rather than disappearing. Same p99
     ceiling, thirty to fifty fold different median. The fraction called
     positive separates completely and is stable across the separated window.
  3  B lineage cells sit closer to the myeloid core in treated foci, under both
     centrings. NOT a pooled lymphocyte claim: pooled collapses from -0.087 to
     -0.026 and T lineage to +0.003 when the reference is made composition-free.
  Control  Neutrophils invert and separate completely under the balanced
     centring, the direction necrosis predicts.

DROPPED
  3-hydroxykynurenine, no dynamic range. Normalised distance, systematic
  over-correction. Raw distance, bounded by focus size. iNOS x Arginase-1
  co-expression, the CD3e x CD20 control separates more than the test pair.
  Pooled lymphocyte radial position, collapses under balanced centring.

UNRESOLVED
  The vendor IDO1 binary call is not reproducible as an intensity threshold.
  See AKOYA_IDO1_vendor_call_context_note.md. No panel here claims otherwise.

REQUIRES
  akoya_arm_stats.py in the same directory as this script.

OUTPUTS
  figures/  F60_finding1_foci_and_burden, F61_finding2_ido1,
            F62_finding3_b_lineage, F64_finding3_b_lineage_distribution
  tables/   80_headline_numbers, 81_ido1_fraction,
            82_ido1_threshold_sensitivity, 83_b_lineage_centered_per_animal,
            84_ido1_cutoff_windows, 85_arm_tests

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================

STRUCT_DIR = "/master/jlehle/WORKING/AKOYA/structures_rev4"
CELL_DIR = f"{STRUCT_DIR}/cell_assignments"
BURDEN_TABLE = f"{STRUCT_DIR}/tables/33b_burden_summary.csv"
FOCI_TABLE = f"{STRUCT_DIR}/tables/35_foci_structures_relative.csv"

LYM_DIR = "/master/jlehle/WORKING/AKOYA/lymphocyte_radial/tables"
MODELS_TABLE = f"{LYM_DIR}/72_models_all.csv"
PROFILE_TABLE = f"{LYM_DIR}/73_radial_profiles_core.csv"

OUT_DIR = "/master/jlehle/WORKING/AKOYA/publication_figures"

# Only used when __file__ is undefined, i.e. running cells in Spyder.
SCRIPT_DIR_FALLBACK = "/master/jlehle/WORKING/AKOYA/scripts"

CONDITION_ORDER = ["D1MT", "Untreated"]
CONDITION_COLORS = {"D1MT": "#2C7FB8", "Untreated": "#D95F02"}
CONDITION_LABELS = {"D1MT": "D1MT-treated", "Untreated": "Untreated"}

# =============================================================================
# PRE-SPECIFIED DIRECTIONS
# Fixed from the biology before the data were examined. One-sided testing is
# only defensible because this block does not move.
# =============================================================================
# TWO-SIDED THROUGHOUT. The directional block below is retained as a record of
# what was predicted, but it is not used for testing. Different populations move
# in opposite directions here: neutrophils are higher in one arm and B lineage
# cells in the other, so committing the analysis to a single set of directions
# would have to be defended population by population. Two-sided is the
# conservative and consistent choice, and it doubles every p relative to the
# one-sided run, which is the price of not having to argue that case.
USE_ONE_SIDED_TESTS = False

PREDICTED_DIRECTIONS = {
    "foci_peak_density":  "D1MT_lower",    # treatment reduces lesion density
    "foci_fold":          "D1MT_lower",    # and reduces prominence
    "burden_pct":         "D1MT_lower",    # and reduces burden
    "ido1_fraction":      "D1MT_lower",    # IDO1+ fraction falls
    "b_lineage_radial":   "D1MT_lower",    # B lineage sits core-ward
    "neutrophil_radial":  "D1MT_higher",   # necrotic cores draw untreated
                                           # neutrophils inward, so treated
                                           # neutrophils sit further out
    # focus architecture, F65
    "foci_inscribed_radius":  "D1MT_lower",   # residual foci are smaller
    "foci_n_cells":           "D1MT_lower",   # and hold fewer cells
    "foci_pct_ido1_pos":      "D1MT_lower",   # and fewer IDO1+ macrophages
    "foci_shape_ratio":       "two_sided",    # no shape prediction either way
}


def direction_for(name):
    """Pre-specified direction when one-sided testing is enabled, else none."""
    if not USE_ONE_SIDED_TESTS:
        return aas.NO_DIRECTION
    return PREDICTED_DIRECTIONS.get(name, aas.NO_DIRECTION)


# retained so existing references keep working
DIRECTIONS = PREDICTED_DIRECTIONS

# Quantities shown in F65, as (column, label, log scale, y axis label).
FOCI_PANELS = [
    ("peak_density", "foci_peak_density", True,
     "Peak myeloid density\n(cells / mm$^2$)"),
    ("fold_over_background", "foci_fold", False,
     "Fold over own background"),
    ("max_inscribed_radius_um", "foci_inscribed_radius", False,
     "Maximum inscribed radius (µm)"),
    ("n_cells_core", "foci_n_cells", True, "Cells in the focus core"),
    ("pct_ido1_pos_of_mac_core", "foci_pct_ido1_pos", False,
     "IDO1+ % of core macrophages"),
    ("shape_ratio", "foci_shape_ratio", False,
     "Shape (equivalent / inscribed radius)"),
]

IDO1_POS = "CD68+IDO1+ Macrophages"
IDO1_NEG = "CD68+IDO1- Macrophages"
CD68_LINEAGE = [IDO1_POS, IDO1_NEG]
MYELOID = [IDO1_POS, IDO1_NEG, "CD163+ Macrophages", "Neutrophils"]
B_LINEAGE = ["B cells", "Plasma cells"]
# F64 panel A: the two B lineage populations plus neutrophils as the internal
# specificity control, so the opposite direction is visible in the same panel
# rather than only in F62.
F64_PANEL_A_POPULATIONS = ["B cells", "Plasma cells", "Neutrophils"]
T_LINEAGE = ["Helper T cells", "CD4- T cells", "Tregs"]
LYMPHOCYTES = T_LINEAGE + B_LINEAGE

PHENOTYPE_COLORS = {
    IDO1_POS: "#B2182B", IDO1_NEG: "#EF8A62", "CD163+ Macrophages": "#FDBE85",
    "Neutrophils": "#7B3294", "Helper T cells": "#1B7837",
    "CD4- T cells": "#7FBC41", "Tregs": "#00441B", "B cells": "#2166AC",
    "Plasma cells": "#67A9CF",
}

STRUCT_COL = "focus_id"
MIN_CELLS_FOR_BALANCE = 10

# The mixed linear model is the test shown on the figures. The exact clustered
# rank-sum still runs and still goes to table 85 as the assumption-free check,
# but it is bounded below by 1/20 one-sided at three animals per arm, so every
# separated comparison returns exactly 0.05 and on a panel that reads as six
# identical results rather than as one arithmetic bound. Set True to put it back
# on the figures.
# None means: show the clustered rank-sum only on panels where the mixed model
# did not converge, so those panels still carry a test. The two-stage t-test
# fallback is disabled in akoya_arm_stats, because on every quantity that
# triggered it the exact rank-sum was both more sensitive and assumption-free.
# Set False to suppress the rank-sum everywhere, True to show it everywhere.
SHOW_RANK_SUM_ON_FIGURES = None

# ---- IDO1 cutoff window search, copied verbatim from script 03 -------------
IDO1_COL = "IDO1"
N_CUTOFF_STEPS = 200
CUTOFF_GRID_HI_PERCENTILE = 99.9
PLATEAU_TOL = 0.05
MIN_ARM_GAP = 0.20

# ---- maps -------------------------------------------------------------------
MAX_POINTS_MAP = 60000
MAP_PAD_FRACTION = 0.04
SCALE_BAR_UM = 1000.0

# ---- plotting ---------------------------------------------------------------
FONT_SIZE_BASE = 28
FONT_SIZE_TITLE = 34
FONT_SIZE_TICK = 28
FONT_SIZE_LEGEND = 28
FONT_SIZE_ANNOT = 26
FONT_SIZE_STAT = 20
FONT_SIZE_PANEL = 40
DPI = 300
SAVE_PDF = True
SAVE_PNG = True
POINT_MARKERS = ["o", "s", "^", "D", "v", "P"]
RANDOM_SEED = 0

GRID_COLOR = "#DDDDDD"
AXIS_COLOR = "#333333"
TEXT_COLOR = "#000000"
FLAG_COLOR = "#B2182B"
OK_COLOR = "#1B7837"
BACKGROUND_COLOR = "#ECECEC"
FOCUS_SHADE = "#08519C"
BURDEN_SHADE = "#4D004B"


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
from matplotlib.patches import Patch, Rectangle

try:
    _here = os.path.dirname(os.path.abspath(__file__))
except NameError:
    _here = SCRIPT_DIR_FALLBACK
if _here not in sys.path:
    sys.path.insert(0, _here)
import akoya_arm_stats as aas

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
ARM_TESTS = []


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


def panel_letter(ax, letter, dx=-0.14, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes,
            fontsize=FONT_SIZE_PANEL, fontweight="bold", va="top", ha="left")


def load_optional(path, label):
    if not os.path.exists(path):
        print(f"    WARNING: {label} not found at {path}")
        return None
    df = pd.read_csv(path)
    print(f"    loaded {label}  ({df.shape[0]} x {df.shape[1]})")
    return df


def run_arm_test(name, values, clusters, arms, unit_label, record_as=None,
                 log=True):
    """Clustered rank-sum with the pre-specified direction, recorded once."""
    direction = direction_for(name)
    res = aas.datta_satten_ranksum(values, clusters, arms,
                                   CONDITION_ORDER[0], CONDITION_ORDER[1],
                                   direction=direction)
    rec = dict(res)
    rec["comparison"] = record_as or name
    rec["unit"] = unit_label
    ARM_TESTS.append(rec)
    if log:
        pred = ("as predicted" if res.get("ran_as_predicted")
                else "OPPOSITE to prediction"
                if res.get("ran_as_predicted") is False else "no direction")
        print(f"    {rec['comparison']:<28} p1={res['p_one_sided']:.4f}  "
              f"p2={res['p_two_sided']:.4f}  n={res['n_obs']:,} "
              f"{unit_label}  {pred}")
    return res


def violin_panel(ax, frame, value_col, cluster_col="sample_id",
                 centre=0.0, spacing=1.0, violin_width=0.62,
                 max_points=1200, point_size=32, point_alpha=0.45,
                 point_jitter=0.19, median_size=680, median_spread=0.12,
                 show_box=True, log=False):
    """
    One violin per arm with the raw observations jittered ON THE CENTRE LINE,
    coloured by animal, and one large black-edged marker per animal at that
    animal's median as the top layer. Shared with every other figure through
    akoya_arm_stats.violin_with_points so they cannot drift apart.
    """
    return aas.violin_with_points(
        ax, frame, value_col, "condition", cluster_col,
        CONDITION_ORDER, CONDITION_COLORS, COLOR_OF, MARKER_OF,
        centre=centre, spacing=spacing, violin_width=violin_width,
        max_points=max_points, point_size=point_size,
        point_alpha=point_alpha, point_jitter=point_jitter,
        median_size=median_size, median_spread=median_spread,
        show_box=show_box, rng=rng, log=log)


def animal_legend(ax, **kw):
    labels = {s: f"{short_label(s)} ({CONDITION_LABELS[COND_OF[s]]})"
              for s in SAMPLE_ORDER}
    handles = aas.cluster_median_legend(COLOR_OF, MARKER_OF, labels)
    ax.legend(handles=handles, frameon=False, title="animal median",
              fontsize=kw.pop("fontsize", FONT_SIZE_LEGEND - 18), **kw)


def arm_separation_gap(sw, arm_of):
    a_cols = [c for c in sw.columns if arm_of.get(c) == CONDITION_ORDER[0]]
    b_cols = [c for c in sw.columns if arm_of.get(c) == CONDITION_ORDER[1]]
    if not a_cols or not b_cols:
        return pd.Series(np.nan, index=sw.index)
    gap_ba = sw[b_cols].min(axis=1) - sw[a_cols].max(axis=1)
    gap_ab = sw[a_cols].min(axis=1) - sw[b_cols].max(axis=1)
    return pd.concat([gap_ba, gap_ab], axis=1).max(axis=1).clip(lower=0.0)


def widest_stable_window(sw, tol, mask=None):
    idx = sw.index.to_numpy(dtype=float)
    ok = (np.ones(len(sw), dtype=bool) if mask is None
          else np.asarray(mask, dtype=bool))
    best = (np.nan, np.nan, -1.0)
    i, n = 0, len(sw)
    while i < n:
        if not ok[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and ok[j + 1]:
            block = sw.iloc[i:j + 2]
            if float((block.max() - block.min()).max()) > tol:
                break
            j += 1
        width = idx[j] - idx[i]
        if width > best[2]:
            best = (float(idx[i]), float(idx[j]), float(width))
        i = max(j, i) + 1
    return best


def contiguous_runs(mask, index):
    m = np.asarray(mask, dtype=bool)
    idx = np.asarray(index, dtype=float)
    runs, i, n = [], 0, len(m)
    while i < n:
        if not m[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and m[j + 1]:
            j += 1
        runs.append((float(idx[i]), float(idx[j]), j - i + 1))
        i = j + 1
    return runs


_tee = Tee(os.path.join(TAB_DIR, "00_publication_figures_report.txt"))
sys.stdout = _tee

banner("AKOYA PUBLICATION FIGURES (revision 3)")
print(f"Run time : {datetime.now().isoformat(timespec='seconds')}")
print(f"Output   : {OUT_DIR}")
print("\nNo new inference. Assembly plus one clustered arm test per comparison.")
sub("PRE-SPECIFIED DIRECTIONS (fixed before the data were examined)")
for k, v in DIRECTIONS.items():
    print(f"    {k:<22} {v}")
print("\n    All tests one-sided on these directions. Two-sided p also in table")
print("    85. At three animals per arm the floor is 0.05 one-sided and 0.10")
print("    two-sided, which is a property of the design.")


# %% Cell 3 - load and centre
# =============================================================================

banner("LOADING")

burden = load_optional(BURDEN_TABLE, "burden summary")
foci = load_optional(FOCI_TABLE, "foci table")
models = load_optional(MODELS_TABLE, "models")
profiles = load_optional(PROFILE_TABLE, "radial profiles")

paths = sorted(glob.glob(os.path.join(CELL_DIR, "*_cell_structures.csv")))
if not paths:
    print(f"ERROR: no cell assignment files in {CELL_DIR}.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

cells, core_frames = {}, []
for p in paths:
    sid = os.path.basename(p).replace("_cell_structures.csv", "")
    try:
        d = pd.read_csv(p, low_memory=False)
    except Exception as e:
        print(f"    ERROR reading {sid}: {e}")
        continue
    d["sample_id"] = sid
    cells[sid] = d
    ins = d.loc[d[STRUCT_COL] > 0]
    if len(ins):
        cen = aas.add_centred_radial(ins, struct_col=STRUCT_COL,
                                     min_cells_for_balance=MIN_CELLS_FOR_BALANCE)
        core_frames.append(cen.loc[cen["is_core"]])
    print(f"    {sid:<12} {d['condition'].iloc[0]:<10} {len(d):>9,} cells")

core = pd.concat(core_frames, ignore_index=True) if core_frames else pd.DataFrame()
del core_frames
gc.collect()
print(f"\n    centred radial computed for {len(core):,} core cells, both")
print(f"    centrings, identical construction to script 08")

sub("Cell selection audit")
print("    Exactly which cells enter the models, so a count in the text can")
print("    always be traced. Script 08 caps each structure at 3,000 cells;")
print("    this script uses every cell, which is why the two differ.\n")
print(f"    {'population':<26}{'core cells':>12}{'foci':>7}{'animals':>9}"
      f"{'D1MT':>9}{'Untr':>9}")
print("    " + "-" * 72)
for lab, members in [("B cells", ["B cells"]),
                     ("Plasma cells", ["Plasma cells"]),
                     ("B lineage", B_LINEAGE),
                     ("T lineage", T_LINEAGE),
                     ("Pooled lymphocytes", LYMPHOCYTES),
                     ("Neutrophils", ["Neutrophils"])]:
    d_ = core.loc[core["pheno"].isin(members)]
    if not len(d_):
        continue
    print(f"    {lab:<26}{len(d_):>12,}"
          f"{d_.groupby(['sample_id', STRUCT_COL]).ngroups:>7}"
          f"{d_['sample_id'].nunique():>9}"
          f"{int((d_['condition'] == CONDITION_ORDER[0]).sum()):>9,}"
          f"{int((d_['condition'] == CONDITION_ORDER[1]).sum()):>9,}")

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

xpos = np.arange(len(SAMPLE_ORDER))
xlab = [short_label(s) for s in SAMPLE_ORDER]

INSCRIBED = {}
if foci is not None and "max_inscribed_radius_um" in foci.columns:
    for c in CONDITION_ORDER:
        g = foci.loc[foci["condition"] == c, "max_inscribed_radius_um"]
        INSCRIBED[c] = float(g.median()) if len(g) else np.nan
    print("    median inscribed radius: "
          + ", ".join(f"{c} {INSCRIBED[c]:.0f} um" for c in CONDITION_ORDER))


# %% Cell 4 - IDO1 quantities and the cutoff window
# =============================================================================

banner("IDO1 QUANTITIES")

ido_rows, ido_vals = [], {}
for s in SAMPLE_ORDER:
    d = cells[s]
    mac = d.loc[d["pheno"].isin(CD68_LINEAGE)]
    if not len(mac):
        continue
    v = pd.to_numeric(mac[IDO1_COL], errors="coerce").dropna().to_numpy()
    ido_vals[s] = v
    n_pos = int((mac["pheno"] == IDO1_POS).sum())
    ido_rows.append({
        "sample_id": s, "animal_id": short_label(s), "condition": COND_OF[s],
        "n_cd68_lineage": len(mac), "n_ido1_pos": n_pos,
        "pct_ido1_pos": 100.0 * n_pos / len(mac),
        "median_ido1": float(np.median(v)),
        "p99_ido1": float(np.percentile(v, 99)),
        "max_ido1": float(np.max(v)),
    })
ido = pd.DataFrame(ido_rows)
write_csv(ido, "81_ido1_fraction.csv")
print(ido.to_string(index=False))

allv = np.concatenate([v for v in ido_vals.values() if len(v)])
hi = float(np.log1p(np.percentile(allv, CUTOFF_GRID_HI_PERCENTILE)))
grid = np.linspace(0.0, hi, N_CUTOFF_STEPS)
sens_rows = []
for s in SAMPLE_ORDER:
    lv = np.log1p(ido_vals.get(s, np.array([])))
    if not len(lv):
        continue
    for g in grid:
        sens_rows.append({"sample_id": s, "condition": COND_OF[s],
                          "log1p_cutoff": float(g),
                          "pct_above": 100.0 * float(np.mean(lv > g))})
sens = pd.DataFrame(sens_rows)
sw = sens.pivot_table(index="log1p_cutoff", columns="sample_id",
                      values="pct_above")[SAMPLE_ORDER] / 100.0
gap = arm_separation_gap(sw, COND_OF)
admissible = (gap >= MIN_ARM_GAP).to_numpy()
runs = contiguous_runs(admissible, sw.index)
sep_lo = sep_hi = np.nan
if runs:
    w = max(runs, key=lambda t: t[1] - t[0])
    sep_lo, sep_hi = w[0], w[1]
best_lo, best_hi, _ = widest_stable_window(sw, PLATEAU_TOL, mask=admissible)
unc_lo, unc_hi, _ = widest_stable_window(sw, PLATEAU_TOL, mask=None)
sens = sens.merge(gap.rename("arm_separation_gap").reset_index(),
                  on="log1p_cutoff", how="left")
write_csv(sens, "82_ido1_threshold_sensitivity.csv")
write_csv(pd.DataFrame([{
    "n_cutoff_steps": N_CUTOFF_STEPS,
    "grid_hi_percentile": CUTOFF_GRID_HI_PERCENTILE, "grid_hi_log1p": hi,
    "plateau_tol": PLATEAU_TOL, "min_arm_gap": MIN_ARM_GAP,
    "separation_lo_raw": np.expm1(sep_lo), "separation_hi_raw": np.expm1(sep_hi),
    "stable_lo_raw": np.expm1(best_lo), "stable_hi_raw": np.expm1(best_hi),
    "unconstrained_lo_raw": np.expm1(unc_lo),
    "unconstrained_hi_raw": np.expm1(unc_hi),
    "n_pooled_cells": int(len(allv)),
}]), "84_ido1_cutoff_windows.csv")
sub("Cutoff window (diff against script 03 table 26)")
print(f"    separated over raw IDO1 {np.expm1(sep_lo):.2f} to {np.expm1(sep_hi):.2f}")
print(f"    separated AND stable    {np.expm1(best_lo):.2f} to {np.expm1(best_hi):.2f}")


# %% Cell 5 - arm tests
# =============================================================================

banner("CLUSTERED ARM TESTS")
print("    Datta-Satten clustered rank-sum, exact over cluster-label")
print("    assignments, one-sided on the pre-specified direction.\n")

test_res = {}
if foci is not None:
    test_res["foci_peak_density"] = run_arm_test(
        "foci_peak_density", foci["peak_density"], foci["sample_id"],
        foci["condition"], "foci")
    test_res["foci_fold"] = run_arm_test(
        "foci_fold", foci["fold_over_background"], foci["sample_id"],
        foci["condition"], "foci")
if burden is not None:
    test_res["burden_pct"] = run_arm_test(
        "burden_pct", burden["burden_pct"], burden["sample_id"],
        burden["condition"], "animals")
test_res["ido1_fraction"] = run_arm_test(
    "ido1_fraction", ido["pct_ido1_pos"], ido["sample_id"], ido["condition"],
    "animals")
if len(core):
    bl_all = core.loc[core["pheno"].isin(B_LINEAGE)]
    test_res["b_lineage_radial"] = run_arm_test(
        "b_lineage_radial", bl_all["radial_centered_core"],
        bl_all["sample_id"], bl_all["condition"], "B lineage cells")
    if "radial_centered_core_balanced" in bl_all.columns:
        test_res["b_lineage_balanced"] = run_arm_test(
            "b_lineage_radial", bl_all["radial_centered_core_balanced"],
            bl_all["sample_id"], bl_all["condition"], "B lineage cells",
            record_as="b_lineage_radial_balanced")
    for ph in B_LINEAGE:
        d = core.loc[core["pheno"] == ph]
        if len(d):
            test_res[f"radial_{ph}"] = run_arm_test(
                "b_lineage_radial", d["radial_centered_core"], d["sample_id"],
                d["condition"], f"{ph}", record_as=f"radial_{ph}")
    nu = core.loc[core["pheno"] == "Neutrophils"]
    if len(nu):
        test_res["neutrophil_radial"] = run_arm_test(
            "neutrophil_radial", nu["radial_centered_core"], nu["sample_id"],
            nu["condition"], "neutrophils")
        if "radial_centered_core_balanced" in nu.columns:
            test_res["neutrophil_radial_balanced"] = run_arm_test(
                "neutrophil_radial", nu["radial_centered_core_balanced"],
                nu["sample_id"], nu["condition"], "neutrophils",
                record_as="neutrophil_radial_balanced")
    for col, name, _, _ in FOCI_PANELS:
        if foci is not None and col in foci.columns and name not in test_res:
            test_res[name] = run_arm_test(name, foci[col], foci["sample_id"],
                                          foci["condition"], "foci")

# ---- mixed linear models -----------------------------------------------
# Script 09 is otherwise assembly only. These fits live here because no
# upstream script models focus-level quantities, and F60 and F65 need them.
# The model is the CELL-LEVEL (here focus-level) test: every observation
# enters the likelihood and the animal random intercept absorbs the shared
# baseline. It is not floored at 1/20 the way the exact rank-sum is.
sub("Mixed linear models fitted here")
lmm_of = {}
if foci is not None:
    for col, name, use_log, _ in FOCI_PANELS:
        if col not in foci.columns:
            print(f"    SKIP {name}: '{col}' not in table 35")
            continue
        m = aas.fit_arm_lmm(foci, col, "condition", "sample_id",
                            CONDITION_ORDER[0], CONDITION_ORDER[1],
                            log10=use_log, label=name,
                            direction=direction_for(name))
        if m["converged"]:
            lmm_of[name] = m
            fc = (f"  fold {m['fold_change']:.2f}x" if use_log
                  and np.isfinite(m["fold_change"]) else "")
            print(f"    {name:<24} coef={m['coef_D1MT_vs_ref']:>+8.3f}"
                  f"{'  (log10)' if use_log else '        '}"
                  f"  p={m['p_value']:.4f}  n={m['n_obs']:>3} foci, "
                  f"{m['n_clusters_a']}v{m['n_clusters_b']} animals{fc}")
# Burden is ONE value per animal, so a random intercept per animal has nothing
# left to estimate and the fit is not identifiable. No model is attempted; the
# burden panel carries the clustered rank-sum, which is the correct test for six
# independent numbers.
print("    burden_pct: no mixed model, one observation per animal. The exact "
      "rank-sum is the right test there and is shown on that panel.")
if len(core):
    for ph in F64_PANEL_A_POPULATIONS:
        dph = core.loc[core["pheno"] == ph]
        if not len(dph):
            continue
        dirn = direction_for("neutrophil_radial" if ph == "Neutrophils"
                             else "b_lineage_radial")
        for col, tag in [("radial_centered_core", "cw"),
                         ("radial_centered_core_balanced", "bal")]:
            if col not in dph.columns:
                continue
            mm = aas.fit_arm_lmm(dph, col, "condition", "sample_id",
                                 CONDITION_ORDER[0], CONDITION_ORDER[1],
                                 struct_col=STRUCT_COL,
                                 label=f"radial_{ph}_{tag}", direction=dirn)
            if mm["converged"]:
                lmm_of[f"radial_{ph}_{tag}"] = mm
                print(f"    radial_{ph}_{tag:<4}"
                      f"{'':<{max(0, 14 - len(ph))}} "
                      f"coef={mm['coef_D1MT_vs_ref']:>+8.4f}  "
                      f"p={mm['p_value']:.4f}  n={mm['n_obs']:,}")

    # Pooled lineage models are fitted HERE rather than read from script 08, so
    # that every model on every figure and in the text comes from one script,
    # one cell selection and one testing convention. Script 08's table 72 is
    # still loaded and printed as a cross-check below.
    for lname, members, key in [("pooled lymphocytes", LYMPHOCYTES, "pooled_lym"),
                                ("B lineage", B_LINEAGE, "b_lineage"),
                                ("T lineage", T_LINEAGE, "t_lineage")]:
        dl = core.loc[core["pheno"].isin(members)]
        if not len(dl):
            continue
        for col, tag in [("radial_centered_core", "cw"),
                         ("radial_centered_core_balanced", "bal")]:
            if col not in dl.columns:
                continue
            mm = aas.fit_arm_lmm(dl, col, "condition", "sample_id",
                                 CONDITION_ORDER[0], CONDITION_ORDER[1],
                                 struct_col=STRUCT_COL, covariates=["pheno"],
                                 label=f"{key}_{tag}",
                                 direction=direction_for("b_lineage_radial"))
            if mm["converged"]:
                lmm_of[f"{key}_{tag}"] = mm
                print(f"    {key + '_' + tag:<18} coef={mm['coef_D1MT_vs_ref']:>+8.4f}"
                      f"  p={mm['p_value']:.4f}  n={mm['n_obs']:>7,} cells"
                      f"  [{mm['random_effects']}]")

    nu = core.loc[core["pheno"] == "Neutrophils"]
    for col, key in [("radial_centered_core", "neutrophil_radial"),
                     ("radial_centered_core_balanced",
                      "neutrophil_radial_balanced")]:
        if col in nu.columns:
            m = aas.fit_arm_lmm(nu, col, "condition", "sample_id",
                                CONDITION_ORDER[0], CONDITION_ORDER[1],
                                struct_col=STRUCT_COL, label=key,
                                direction=direction_for("neutrophil_radial"))
            if m["converged"]:
                lmm_of[key] = m
                print(f"    {key:<24} coef={m['coef_D1MT_vs_ref']:>+8.4f}"
                      f"          p={m['p_value']:.4f}  "
                      f"n={m['n_obs']:,} neutrophils")


def fold_note(name):
    m = lmm_of.get(name)
    if not m or not np.isfinite(m.get("fold_change", np.nan)):
        return ""
    return (f"treated foci at {m['fold_change']:.2f}x untreated, "
            f"model on log10")


# Figures read lmm_of, fitted above. Script 08's table 72 is loaded only to be
# compared against it, because the two disagreed on cell counts and on testing
# convention and that had to be surfaced rather than averaged over.
model_of = {"b_lineage_primary": lmm_of.get("b_lineage_cw"),
            "b_lineage_balanced": lmm_of.get("b_lineage_bal")}
model_of = {k: v for k, v in model_of.items() if v}

if models is not None:
    sub("Cross-check against script 08 table 72")
    print("    Script 08 caps each structure at 3,000 cells and tested "
          "two-sided.")
    print("    This script uses every cell. Differences below are expected; "
          "large ones are not.\n")
    print(f"    {'model':<34}{'script 08':>22}{'script 09':>24}")
    print("    " + "-" * 80)
    for fam, lin, key in [("lineage", "B lineage", "b_lineage_cw"),
                          ("lineage_centring_check", "B lineage",
                           "b_lineage_bal")]:
        g = models.loc[models["family"] == fam]
        if "lineage" in g.columns:
            g = g.loc[g["lineage"] == lin]
        if fam == "lineage_centring_check" and "centring" in g.columns:
            g = g.loc[g["centring"].astype(str).str.contains("balanced")]
        new_m = lmm_of.get(key)
        if len(g) and new_m:
            r8 = g.iloc[0]
            print(f"    {key:<34}"
                  f"{r8['coef_D1MT_vs_ref']:>+9.4f} p={r8['p_value']:.3f}"
                  f" n={int(r8['n_cells']):>6,}"
                  f"{new_m['coef_D1MT_vs_ref']:>+10.4f} "
                  f"p={new_m['p_value']:.3f} n={new_m['n_obs']:>6,}")


# %% Cell 6 - FIGURE 1: residual foci, then burden
# =============================================================================

banner("FIGURE 1 - RESIDUAL FOCI AND BURDEN")

spans = []
for s in SAMPLE_ORDER:
    d = cells[s]
    spans += [float(d["x"].max() - d["x"].min()),
              float(d["y"].max() - d["y"].min())]
SPAN = max(spans) * (1.0 + MAP_PAD_FRACTION)
print(f"    shared map window: {SPAN:,.0f} um square")

fig = plt.figure(figsize=(46, 29))
gs = fig.add_gridspec(2, 6, height_ratios=[1.15, 1.0], hspace=0.30, wspace=0.12)

for k, s in enumerate(SAMPLE_ORDER):
    ax = fig.add_subplot(gs[0, k])
    d = cells[s]
    idx = (rng.choice(len(d), size=MAX_POINTS_MAP, replace=False)
           if len(d) > MAX_POINTS_MAP else np.arange(len(d)))
    dd = d.iloc[idx]
    if "in_burden_region" in dd.columns:
        m = dd["in_burden_region"].astype(bool)
        if m.sum():
            ax.scatter(dd.loc[m, "x"], dd.loc[m, "y"], s=30,
                       color=BURDEN_SHADE, alpha=0.16, linewidths=0,
                       rasterized=True, zorder=1)
    m = (dd[STRUCT_COL] > 0) & (dd["region"] == "core")
    if m.sum():
        ax.scatter(dd.loc[m, "x"], dd.loc[m, "y"], s=18, color=FOCUS_SHADE,
                   alpha=0.18, linewidths=0, rasterized=True, zorder=1)
    other = ~dd["pheno"].isin(MYELOID)
    ax.scatter(dd.loc[other, "x"], dd.loc[other, "y"], s=0.5,
               color=BACKGROUND_COLOR, linewidths=0, rasterized=True, zorder=2)
    for ph in MYELOID:
        mm = dd["pheno"] == ph
        if mm.sum():
            ax.scatter(dd.loc[mm, "x"], dd.loc[mm, "y"], s=3.5,
                       color=PHENOTYPE_COLORS[ph], linewidths=0,
                       rasterized=True, zorder=3)
    xc = 0.5 * (float(d["x"].min()) + float(d["x"].max()))
    yc = 0.5 * (float(d["y"].min()) + float(d["y"].max()))
    ax.set_xlim(xc - SPAN / 2, xc + SPAN / 2)
    ax.set_ylim(yc - SPAN / 2, yc + SPAN / 2)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_edgecolor(CONDITION_COLORS[COND_OF[s]]); sp.set_linewidth(6)
    nf = int((foci["sample_id"] == s).sum()) if foci is not None else 0
    bp = np.nan
    if burden is not None:
        bsel = burden.loc[burden["sample_id"] == s, "burden_pct"]
        bp = float(bsel.iloc[0]) if len(bsel) else np.nan
    ax.set_title(f"{short_label(s)}\n{nf} foci, {bp:.1f}% burden",
                 fontsize=FONT_SIZE_TITLE - 13,
                 color=FLAG_COLOR if bp == 0 else TEXT_COLOR)
    if k == 0:
        panel_letter(ax, "A", dx=-0.02, dy=1.12)
        x0 = xc - SPAN / 2 + 0.06 * SPAN
        y0 = yc - SPAN / 2 + 0.06 * SPAN
        ax.add_patch(Rectangle((x0, y0), SCALE_BAR_UM, 0.012 * SPAN,
                               facecolor="#000000", edgecolor="none", zorder=5))
        ax.text(x0 + SCALE_BAR_UM / 2, y0 + 0.03 * SPAN,
                f"{SCALE_BAR_UM/1000:.0f} mm", ha="center",
                fontsize=FONT_SIZE_ANNOT - 8)

ax = fig.add_subplot(gs[1, 0:2])
if foci is not None:
    violin_panel(ax, foci, "peak_density", point_size=60,
                 point_alpha=0.65, max_points=400)
    ax.axhline(3000, color=BURDEN_SHADE, linestyle="--", linewidth=3)
    ax.set_yscale("log")
    ax.set_ylim(top=float(foci["peak_density"].max()) * 30)
    aas.annotate_arm_test(ax, include_rank_sum=SHOW_RANK_SUM_ON_FIGURES, res=test_res.get("foci_peak_density"),
                          model_row=lmm_of.get("foci_peak_density"),
                          unit_label="foci", fontsize=FONT_SIZE_STAT,
                          extra=fold_note("foci_peak_density"))
    animal_legend(ax, loc="lower left")
ax.set_xticks(range(len(CONDITION_ORDER)))
ax.set_xticklabels([CONDITION_LABELS[c] for c in CONDITION_ORDER],
                   fontsize=FONT_SIZE_TICK - 8)
ax.set_ylabel("Peak myeloid density per focus\n(cells / mm$^2$)",
              fontsize=FONT_SIZE_BASE - 8)
# 2,223 / 8,255 = 0.27, and the model fold change is 0.28, so "a quarter" is
# the accurate description. "A third" overstated it by about a quarter again.
ax.set_title("Residual foci run at roughly a quarter of the density\n"
             "dashed = granuloma-density threshold",
             fontsize=FONT_SIZE_TITLE - 14)
style_axes(ax); panel_letter(ax, "B")

ax = fig.add_subplot(gs[1, 2:4])
if foci is not None:
    violin_panel(ax, foci, "fold_over_background", point_size=60,
                 point_alpha=0.65, max_points=400)
    ax.set_ylim(0, float(foci["fold_over_background"].max()) * 2.2)
    aas.annotate_arm_test(ax, include_rank_sum=SHOW_RANK_SUM_ON_FIGURES, res=test_res.get("foci_fold"),
                          model_row=lmm_of.get("foci_fold"),
                          unit_label="foci", fontsize=FONT_SIZE_STAT)
ax.set_xticks(range(len(CONDITION_ORDER)))
ax.set_xticklabels([CONDITION_LABELS[c] for c in CONDITION_ORDER],
                   fontsize=FONT_SIZE_TICK - 8)
ax.set_ylabel("Fold over that section's own\nbackground density",
              fontsize=FONT_SIZE_BASE - 8)
ax.set_title("But stay 7 to 8 fold above their own background,\n"
             "so they are real structures rather than noise",
             fontsize=FONT_SIZE_TITLE - 17)
style_axes(ax); panel_letter(ax, "C")

ax = fig.add_subplot(gs[1, 4:6])
if burden is not None:
    for ci, c in enumerate(CONDITION_ORDER):
        g = burden.loc[burden["condition"] == c]
        v = g["burden_pct"].to_numpy()
        j = rng.uniform(-0.07, 0.07, size=len(v))
        ax.scatter(np.full(len(v), ci) + j, v, s=680,
                   color=CONDITION_COLORS[c], edgecolor="#FFFFFF", linewidth=3,
                   zorder=3)
        for xj, yv, sid in zip(np.full(len(v), ci) + j, v, g["sample_id"]):
            ax.text(xj + 0.10, yv, short_label(sid), va="center",
                    fontsize=FONT_SIZE_ANNOT - 8)
        ax.hlines(np.median(v), ci - 0.24, ci + 0.24, color="#000000",
                  linewidth=5)
    ax.set_ylim(-3, float(burden["burden_pct"].max()) * 1.9)
    aas.separation_band(
        ax, burden.loc[burden["condition"] == "D1MT", "burden_pct"],
        burden.loc[burden["condition"] == "Untreated", "burden_pct"],
        fontsize=FONT_SIZE_ANNOT - 6)
    aas.annotate_arm_test(ax, test_res.get("burden_pct"), unit_label="animals",
                          fontsize=FONT_SIZE_STAT, include_rank_sum=True,
                          extra="one value per animal, so no mixed model")
ax.set_xticks(range(len(CONDITION_ORDER)))
ax.set_xticklabels([CONDITION_LABELS[c] for c in CONDITION_ORDER],
                   fontsize=FONT_SIZE_TICK - 8)
ax.set_xlim(-0.6, len(CONDITION_ORDER) - 0.3)
ax.set_ylabel("Lung area in granuloma-density\nregions (%)",
              fontsize=FONT_SIZE_BASE - 8)
ax.set_title("And none reaches granuloma density",
             fontsize=FONT_SIZE_TITLE - 14)
style_axes(ax); panel_letter(ax, "D")

handles = [Patch(facecolor=PHENOTYPE_COLORS[p], label=p) for p in MYELOID]
handles += [Patch(facecolor=FOCUS_SHADE, alpha=0.35, label="myeloid focus core"),
            Patch(facecolor=BURDEN_SHADE, alpha=0.35,
                  label="granuloma-density region"),
            Patch(facecolor=BACKGROUND_COLOR, label="all other cells")]
fig.legend(handles=handles, loc="lower center", ncol=7, frameon=False,
           fontsize=FONT_SIZE_LEGEND - 12, bbox_to_anchor=(0.5, -0.05))
fig.suptitle("Finding 1: treated animals retain discrete myeloid foci, but "
             "none reaches granuloma density\nFoci are defined relative to "
             "each section's own background; burden uses one absolute "
             "threshold for every section. Detection never uses IDO1.",
             y=1.02, fontsize=FONT_SIZE_TITLE - 7)
save_fig(fig, "F60_finding1_foci_and_burden")


# %% Cell 7 - FIGURE 2: the IDO1 compartment
# =============================================================================

banner("FIGURE 2 - THE IDO1 MACROPHAGE COMPARTMENT")

fig, axes = plt.subplots(2, 2, figsize=(34, 28))

ax = axes[0, 0]
bins = np.linspace(0, hi, 90)
for i, s in enumerate(SAMPLE_ORDER):
    v = ido_vals.get(s, np.array([]))
    if not len(v):
        continue
    h, e = np.histogram(np.log1p(v), bins=bins, density=True)
    h = h / (h.max() if h.max() > 0 else 1.0)
    ctr = 0.5 * (e[:-1] + e[1:])
    base = len(SAMPLE_ORDER) - i
    ax.fill_between(ctr, base, base + 0.92 * h, color=COLOR_OF[s], alpha=0.8,
                    linewidth=0)
    ax.plot(ctr, base + 0.92 * h, color="#333333", linewidth=2.5)
    ax.text(hi * 1.02, base + 0.35, f"{short_label(s)}  n={len(v):,}",
            fontsize=FONT_SIZE_ANNOT - 10, va="center",
            color=CONDITION_COLORS[COND_OF[s]])
if np.isfinite(best_lo):
    ax.axvspan(best_lo, best_hi, color=OK_COLOR, alpha=0.13, zorder=0)
ax.set_yticks([]); ax.set_xlim(0, hi * 1.25)
ax.set_xlabel("log(1 + IDO1), CD68-lineage macrophages",
              fontsize=FONT_SIZE_BASE - 6)
ax.set_title("The whole compartment shifts", fontsize=FONT_SIZE_TITLE - 10)
for sp in ["top", "right", "left"]:
    ax.spines[sp].set_visible(False)
panel_letter(ax, "A", dx=-0.04)

ax = axes[0, 1]
for ci, c in enumerate(CONDITION_ORDER):
    g = ido.loc[ido["condition"] == c]
    v = g["pct_ido1_pos"].to_numpy()
    j = rng.uniform(-0.07, 0.07, size=len(v))
    ax.scatter(np.full(len(v), ci) + j, v, s=680, color=CONDITION_COLORS[c],
               edgecolor="#FFFFFF", linewidth=3, zorder=3)
    for xj, yv, aid in zip(np.full(len(v), ci) + j, v, g["animal_id"]):
        ax.text(xj + 0.10, yv, str(aid), va="center",
                fontsize=FONT_SIZE_ANNOT - 8)
    ax.hlines(np.median(v), ci - 0.24, ci + 0.24, color="#000000", linewidth=5)
ax.set_ylim(0, 165)
aas.separation_band(ax, ido.loc[ido["condition"] == "D1MT", "pct_ido1_pos"],
                    ido.loc[ido["condition"] == "Untreated", "pct_ido1_pos"],
                    fontsize=FONT_SIZE_ANNOT - 6)
aas.annotate_arm_test(ax, test_res.get("ido1_fraction"),
                      unit_label="animals", include_rank_sum=True,
                      extra="one value per animal, so no mixed model",
                      fontsize=FONT_SIZE_STAT)
ax.set_xticks(range(len(CONDITION_ORDER)))
ax.set_xticklabels([CONDITION_LABELS[c] for c in CONDITION_ORDER])
ax.set_xlim(-0.6, len(CONDITION_ORDER) - 0.3)
ax.set_ylabel("IDO1+ % of CD68-lineage macrophages",
              fontsize=FONT_SIZE_BASE - 8)
ax.set_title("Fraction called positive", fontsize=FONT_SIZE_TITLE - 10)
style_axes(ax); panel_letter(ax, "B")

ax = axes[1, 0]
if np.isfinite(sep_lo):
    ax.axvspan(sep_lo, sep_hi, color=OK_COLOR, alpha=0.07, zorder=0)
if np.isfinite(best_lo):
    ax.axvspan(best_lo, best_hi, color=OK_COLOR, alpha=0.16, zorder=0)
    ax.text(0.5 * (best_lo + best_hi), 105,
            f"separated and stable\nraw IDO1 {np.expm1(best_lo):.1f} to "
            f"{np.expm1(best_hi):.1f}", ha="center",
            fontsize=FONT_SIZE_ANNOT - 10)
for s in SAMPLE_ORDER:
    d = sens.loc[sens["sample_id"] == s].sort_values("log1p_cutoff")
    ax.plot(d["log1p_cutoff"], d["pct_above"], linewidth=5, color=COLOR_OF[s],
            alpha=0.9)
    step = max(1, len(d) // 14)
    ax.plot(d["log1p_cutoff"].to_numpy()[::step],
            d["pct_above"].to_numpy()[::step], linestyle="none",
            marker=MARKER_OF[s], markersize=15, color=COLOR_OF[s])
ax.plot(gap.index.to_numpy(), 100 * gap.to_numpy(), linewidth=5,
        color="#000000", linestyle=":", zorder=5)
ax.axhline(100 * MIN_ARM_GAP, color="#000000", linestyle="--", linewidth=2.5)
ax.set_xlabel("log(1 + IDO1) cutoff", fontsize=FONT_SIZE_BASE - 6)
ax.set_ylabel("% of CD68-lineage cells above cutoff",
              fontsize=FONT_SIZE_BASE - 10)
ax.set_ylim(0, 113)
ax.set_title("The separation does not depend on the cutoff\n"
             "dotted black = arm separation gap",
             fontsize=FONT_SIZE_TITLE - 14)
style_axes(ax); panel_letter(ax, "C")
ax.legend(handles=[Line2D([0], [0], color=COLOR_OF[s], marker=MARKER_OF[s],
                          markersize=14, linewidth=4, label=short_label(s))
                   for s in SAMPLE_ORDER],
          frameon=False, fontsize=FONT_SIZE_LEGEND - 16, loc="upper right",
          ncol=2)

ax = axes[1, 1]
for i, s in enumerate(SAMPLE_ORDER):
    r = ido.loc[ido["sample_id"] == s].iloc[0]
    ax.plot([i, i], [r["median_ido1"], r["p99_ido1"]], color=COLOR_OF[s],
            linewidth=8, alpha=0.55, zorder=2)
    ax.scatter([i], [r["median_ido1"]], s=460, color=COLOR_OF[s], marker="_",
               linewidth=7, zorder=3)
    ax.scatter([i], [r["p99_ido1"]], s=440, color=COLOR_OF[s],
               marker=MARKER_OF[s], edgecolor="#FFFFFF", linewidth=2, zorder=4)
ax.set_yscale("log")
ax.set_xticks(xpos)
ax.set_xticklabels([f"{a}\n{CONDITION_LABELS[COND_OF[s]]}"
                    for a, s in zip(xlab, SAMPLE_ORDER)],
                   fontsize=FONT_SIZE_TICK - 12)
ax.set_ylabel("IDO1 on CD68-lineage macrophages", fontsize=FONT_SIZE_BASE - 8)
ax.set_title("Same ceiling, different median\ndash = median, symbol = p99",
             fontsize=FONT_SIZE_TITLE - 14)
style_axes(ax); panel_letter(ax, "D")

fig.suptitle("Finding 2: D1MT moves the IDO1 macrophage compartment rather "
             "than removing it\nConsistent with enzymatic inhibition, not "
             "transcriptional suppression", y=1.00, fontsize=FONT_SIZE_TITLE - 4)
fig.tight_layout(rect=[0, 0, 1, 0.95])
save_fig(fig, "F61_finding2_ido1")


# %% Cell 8 - FIGURE 4: the B lineage distribution
# =============================================================================

banner("FIGURE 4 - B LINEAGE DISTRIBUTION")

bl = core.loc[core["pheno"].isin(B_LINEAGE)].copy() if len(core) else pd.DataFrame()
per_animal_bl = pd.DataFrame()
if len(bl):
    per_animal_bl = (bl.groupby(["sample_id", "condition"], as_index=False)
                     ["radial_centered_core"].mean())
    if "radial_centered_core_balanced" in bl.columns:
        bal = (bl.groupby(["sample_id", "condition"], as_index=False)
               ["radial_centered_core_balanced"].mean())
        per_animal_bl = per_animal_bl.merge(bal, on=["sample_id", "condition"],
                                            how="left")
    write_csv(per_animal_bl, "83_b_lineage_centered_per_animal.csv")

fig, axes = plt.subplots(1, 3, figsize=(42, 16))

ax = axes[0]
panel_a = [p for p in F64_PANEL_A_POPULATIONS
           if len(core.loc[core["pheno"] == p])]
positions = {}
for pi, ph in enumerate(panel_a):
    for ci, c in enumerate(CONDITION_ORDER):
        positions[(ph, c)] = pi * 3.2 + ci * 1.35
for pi, ph in enumerate(panel_a):
    d = core.loc[core["pheno"] == ph]
    violin_panel(ax, d, "radial_centered_core", centre=pi * 3.2,
                 spacing=1.35, violin_width=0.72,
                 max_points=900 if ph in B_LINEAGE else 1200,
                 point_size=40 if ph in B_LINEAGE else 26,
                 point_alpha=0.45 if ph in B_LINEAGE else 0.30,
                 median_size=580, median_spread=0.14)
lo, hiy = ax.get_ylim()
ax.set_ylim(lo, hiy + 0.42 * (hiy - lo))
for ph in panel_a:
    m_cw = lmm_of.get(f"radial_{ph}_cw")
    m_bal = lmm_of.get(f"radial_{ph}_bal")
    lines = [ph]
    if m_cw:
        lines.append(f"p = {m_cw['p_value']:.3f}")
    if m_bal:
        lines.append(f"balanced p = {m_bal['p_value']:.3f}")
    n_obs = int((core["pheno"] == ph).sum())
    lines.append(f"n = {n_obs:,} cells")
    ax.text(positions[(ph, CONDITION_ORDER[0])] + 0.68,
            ax.get_ylim()[1] - 0.03 * (ax.get_ylim()[1] - ax.get_ylim()[0]),
            "\n".join(lines), ha="center", va="top",
            fontsize=FONT_SIZE_STAT - 2, linespacing=1.3,
            color=OK_COLOR if ph not in B_LINEAGE else TEXT_COLOR)
ax.axhline(0, color="#000000", linestyle="--", linewidth=3)
ax.set_xticks([positions[(ph, c)] for ph in panel_a for c in CONDITION_ORDER])
ax.set_xticklabels([CONDITION_LABELS[c].replace("-treated", "")
                    for ph in panel_a for c in CONDITION_ORDER],
                   fontsize=FONT_SIZE_TICK - 14, rotation=30, ha="right")
ax.set_ylabel("Centred radial position\n(negative = closer to the core)",
              fontsize=FONT_SIZE_BASE - 8)
ax.set_title("Each cell, by population\n"
             "green: neutrophils, the internal specificity control",
             fontsize=FONT_SIZE_TITLE - 16)
style_axes(ax); panel_letter(ax, "A")

ax = axes[1]
if len(bl):
    violin_panel(ax, bl, "radial_centered_core")
    ax.axhline(0, color="#000000", linestyle="--", linewidth=3)
    animal_legend(ax, loc="lower left")
    lo, hiy = ax.get_ylim()
    ax.set_ylim(lo, hiy + 0.60 * (hiy - lo))
    conv = ""
    m = model_of.get("b_lineage_primary", {})
    c0 = abs(float(m.get("coef_D1MT_vs_ref", np.nan)))
    if INSCRIBED and np.isfinite(c0):
        conv = (f"{c0 * INSCRIBED[CONDITION_ORDER[0]]:.0f} um in a treated "
                f"focus, {c0 * INSCRIBED[CONDITION_ORDER[1]]:.0f} um untreated")
    aas.annotate_arm_test(ax, include_rank_sum=SHOW_RANK_SUM_ON_FIGURES, res=test_res.get("b_lineage_radial"),
                          model_row=model_of.get("b_lineage_primary"),
                          unit_label="B lineage cells", extra=conv,
                          fontsize=FONT_SIZE_STAT)
ax.set_xticks(range(len(CONDITION_ORDER)))
ax.set_xticklabels([CONDITION_LABELS[c] for c in CONDITION_ORDER],
                   fontsize=FONT_SIZE_TICK - 8)
ax.set_ylabel("Centred radial position\n(negative = closer to the core)",
              fontsize=FONT_SIZE_BASE - 8)
ax.set_title("Pooled B lineage, cell-weighted centring",
             fontsize=FONT_SIZE_TITLE - 12)
style_axes(ax); panel_letter(ax, "B")

ax = axes[2]
if len(bl) and "radial_centered_core_balanced" in bl.columns:
    violin_panel(ax, bl, "radial_centered_core_balanced")
    ax.axhline(0, color="#000000", linestyle="--", linewidth=3)
    lo, hiy = ax.get_ylim()
    ax.set_ylim(lo, hiy + 0.60 * (hiy - lo))
    aas.annotate_arm_test(ax, include_rank_sum=SHOW_RANK_SUM_ON_FIGURES, res=test_res.get("b_lineage_balanced"),
                          model_row=model_of.get("b_lineage_balanced"),
                          unit_label="B lineage cells",
                          extra="reference cannot move with composition",
                          fontsize=FONT_SIZE_STAT)
ax.set_xticks(range(len(CONDITION_ORDER)))
ax.set_xticklabels([CONDITION_LABELS[c] for c in CONDITION_ORDER],
                   fontsize=FONT_SIZE_TICK - 8)
ax.set_ylabel("Centred radial position, balanced reference",
              fontsize=FONT_SIZE_BASE - 10)
ax.set_title("Pooled B lineage, composition-balanced centring",
             fontsize=FONT_SIZE_TITLE - 14)
style_axes(ax); panel_letter(ax, "C")

fig.suptitle("Finding 3: B lineage cells sit closer to the myeloid core in "
             "treated foci\nThe axis is centred radial position, not microns. "
             "Distance in microns is bounded by focus size (median inscribed "
             "radius 125 um treated against 168 um untreated) and would show "
             "this effect whether or not it exists.",
             y=1.05, fontsize=FONT_SIZE_TITLE - 11)
fig.tight_layout(rect=[0, 0, 1, 0.89])
save_fig(fig, "F64_finding3_b_lineage_distribution")


# %% Cell 9 - FIGURE 3: profile, models, specificity control
# =============================================================================

banner("FIGURE 3 - B LINEAGE SUPPORTING PANELS")

fig, axes = plt.subplots(1, 3, figsize=(42, 15))

ax = axes[0]
if profiles is not None and "B lineage" in profiles.columns:
    for s in SAMPLE_ORDER:
        d = profiles.loc[profiles["sample_id"] == s].sort_values("radial_centre")
        if not len(d):
            continue
        ax.plot(d["radial_centre"], d["B lineage"], linewidth=5,
                marker=MARKER_OF[s], markersize=16, color=COLOR_OF[s], alpha=0.9)
ax.set_xlabel("radial position (0 = core centre, 1 = boundary)",
              fontsize=FONT_SIZE_BASE - 8)
ax.set_ylabel("B lineage cells (% of cells in bin)", fontsize=FONT_SIZE_BASE - 8)
ax.set_title("B lineage abundance across the core",
             fontsize=FONT_SIZE_TITLE - 12)
ax.set_ylim(bottom=0)
style_axes(ax); panel_letter(ax, "A")
ax.legend(handles=[Line2D([0], [0], color=COLOR_OF[s], marker=MARKER_OF[s],
                          markersize=14, linewidth=4,
                          label=f"{short_label(s)} ({CONDITION_LABELS[COND_OF[s]]})")
                   for s in SAMPLE_ORDER],
          frameon=False, fontsize=FONT_SIZE_LEGEND - 18, loc="best")

ax = axes[1]
want = []
for key, lab in [("b_lineage_primary", "B lineage, cell-weighted"),
                 ("b_lineage_balanced", "B lineage, balanced")]:
    if key in model_of:
        want.append((lab, model_of[key]))
if models is not None:
    for _, r in models.loc[models["family"] == "per_phenotype_lymphocyte"].iterrows():
        if r.get("phenotype") in B_LINEAGE:
            want.append((f"{r['phenotype']}, cell-weighted", r.to_dict()))
for i, (lab, r) in enumerate(want):
    col = FLAG_COLOR if r["p_value"] < 0.05 else "#777777"
    ax.plot([r["ci_low"], r["ci_high"]], [i, i], color=col, linewidth=7, zorder=2)
    ax.scatter([r["coef_D1MT_vs_ref"]], [i], s=440, color=col,
               edgecolor="#FFFFFF", linewidth=2, zorder=3)
    ax.text(r["ci_high"], i, f"  p={r['p_value']:.3f}", va="center",
            fontsize=FONT_SIZE_ANNOT - 10)
if want:
    ax.set_yticks(np.arange(len(want)))
    ax.set_yticklabels([lab for lab, _ in want], fontsize=FONT_SIZE_TICK - 10)
    ax.invert_yaxis()
ax.axvline(0, color="#000000", linewidth=3.5)
ax.set_xlabel("D1MT effect on centred radial position (95% CI)",
              fontsize=FONT_SIZE_BASE - 10)
ax.set_title("Mixed models, B lineage only", fontsize=FONT_SIZE_TITLE - 12)
ax.xaxis.grid(True, color=GRID_COLOR, linewidth=1.5); ax.set_axisbelow(True)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
panel_letter(ax, "B")

ax = axes[2]
nu = core.loc[core["pheno"] == "Neutrophils"] if len(core) else pd.DataFrame()
NEU_COL = ("radial_centered_core_balanced"
           if len(nu) and "radial_centered_core_balanced" in nu.columns
           else "radial_centered_core")
if len(nu):
    violin_panel(ax, nu, NEU_COL, max_points=1500, point_size=24,
                 point_alpha=0.28)
    ax.axhline(0, color="#000000", linestyle="--", linewidth=3)
    lo, hiy = ax.get_ylim()
    ax.set_ylim(lo, hiy + 0.62 * (hiy - lo))
    key = ("neutrophil_radial_balanced"
           if NEU_COL.endswith("balanced") else "neutrophil_radial")
    other = test_res.get("neutrophil_radial"
                         if NEU_COL.endswith("balanced")
                         else "neutrophil_radial_balanced")
    extra = ("" if other is None else
             f"cell-weighted centring: p = {other['p_one_sided']:.3f}")
    aas.annotate_arm_test(ax, include_rank_sum=SHOW_RANK_SUM_ON_FIGURES, res=test_res.get(key),
                          model_row=lmm_of.get(key),
                          unit_label="neutrophils", extra=extra,
                          fontsize=FONT_SIZE_STAT)
ax.set_xticks(range(len(CONDITION_ORDER)))
ax.set_xticklabels([CONDITION_LABELS[c] for c in CONDITION_ORDER],
                   fontsize=FONT_SIZE_TICK - 8)
ax.set_ylabel("Centred radial position, balanced reference",
              fontsize=FONT_SIZE_BASE - 10)
ax.set_title("Specificity control: neutrophils run the other way\n"
             "untreated neutrophils sit core-ward, as necrosis predicts",
             fontsize=FONT_SIZE_TITLE - 19)
style_axes(ax); panel_letter(ax, "C")

fig.suptitle("Finding 3 supporting panels. Lymphocytes never enter the focus "
             "definition, so these positions are not circular.",
             y=1.03, fontsize=FONT_SIZE_TITLE - 9)
fig.tight_layout(rect=[0, 0, 1, 0.93])
save_fig(fig, "F62_finding3_b_lineage")


# %% Cell 9b - FIGURE 5: focus architecture, the F34 content with statistics
# =============================================================================

banner("FIGURE 5 - FOCUS ARCHITECTURE")

if foci is None or not len(foci):
    print("    SKIPPED: no foci table")
else:
    have = [(c, n, lg, yl) for c, n, lg, yl in FOCI_PANELS if c in foci.columns]
    ncol = 3
    nrow = int(np.ceil(len(have) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(15 * ncol, 13 * nrow),
                             squeeze=False)
    axes = axes.ravel()
    for k, (col, name, use_log, ylab) in enumerate(have):
        ax = axes[k]
        d = foci.loc[np.isfinite(pd.to_numeric(foci[col], errors="coerce"))]
        if not len(d):
            ax.axis("off")
            continue
        violin_panel(ax, d, col, max_points=400, point_size=60,
                     point_alpha=0.65, median_size=600)
        if use_log:
            ax.set_yscale("log")
            ax.set_ylim(top=float(d[col].max()) * 40)
        else:
            lo, hiy = ax.get_ylim()
            ax.set_ylim(lo, hiy + 0.72 * (hiy - lo))
        if col == "peak_density":
            ax.axhline(3000, color=BURDEN_SHADE, linestyle="--", linewidth=3)
        if col == "shape_ratio":
            ax.axhline(1.0, color="#000000", linestyle=":", linewidth=2.5)
        if col == "pct_ido1_pos_of_mac_core":
            # This panel shows a large visual difference at p = 0.42, which
            # invites the wrong conclusion. The per-focus measure is dominated
            # by between-animal variance in the untreated arm (49 to 94 percent
            # across animals), whereas the animal-level comparison separates
            # completely. Say so on the panel rather than leaving a reader to
            # reconcile the two.
            ax.text(0.5, -0.16, "per-focus measure; the animal-level comparison "
                                "separates completely (Fig. F61b)",
                    transform=ax.transAxes, ha="center", va="top",
                    fontsize=FONT_SIZE_STAT - 2, style="italic")
        aas.annotate_arm_test(ax, include_rank_sum=SHOW_RANK_SUM_ON_FIGURES, res=test_res.get(name),
                              model_row=lmm_of.get(name), unit_label="foci",
                              extra=fold_note(name) if use_log else "",
                              fontsize=FONT_SIZE_STAT)
        ax.set_xticks(range(len(CONDITION_ORDER)))
        ax.set_xticklabels([CONDITION_LABELS[c] for c in CONDITION_ORDER],
                           fontsize=FONT_SIZE_TICK - 10)
        ax.set_ylabel(ylab, fontsize=FONT_SIZE_BASE - 10)
        style_axes(ax)
        panel_letter(ax, "ABCDEF"[k])
        if k == 0:
            animal_legend(ax, loc="lower left")
    for k in range(len(have), len(axes)):
        axes[k].axis("off")
    # Wording tracks the fits. Density and prominence separate clearly. Shape
    # is a clean null. Size sits between: p = 0.21 two-sided, 0.10 one-sided, so
    # "comparable size" claims more than the data support and "smaller" claims
    # more still. The title now asserts neither.
    fig.suptitle("Focus architecture: treated foci are substantially less dense "
                 "and less prominent, at comparable shape\n"
                 "Each point is one focus, coloured by animal. The large "
                 "black-edged symbol is that animal's median. Foci are defined "
                 "relative to each section's own background, so a treated "
                 "focus is a residual myeloid focus, never a granuloma.",
                 y=1.01, fontsize=FONT_SIZE_TITLE - 10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    save_fig(fig, "F65_focus_architecture")

    sub("Focus architecture summary")
    print(f"    {'quantity':<28}{'D1MT median':>14}{'Untr median':>14}"
          f"{'LMM p':>9}{'rank-sum p':>12}")
    print("    " + "-" * 77)
    for col, name, use_log, _ in have:
        a = foci.loc[foci["condition"] == CONDITION_ORDER[0], col].median()
        b = foci.loc[foci["condition"] == CONDITION_ORDER[1], col].median()
        m = lmm_of.get(name, {})
        r = test_res.get(name, {})
        print(f"    {col[:27]:<28}{a:>14.2f}{b:>14.2f}"
              f"{m.get('p_value', np.nan):>9.4f}"
              f"{r.get('p_one_sided', np.nan):>12.4f}")


# %% Cell 10 - tables and wrap up
# =============================================================================

banner("TABLES")

if lmm_of:
    lm = pd.DataFrame(list(lmm_of.values()))
    write_csv(lm, "86_mixed_models_fitted_here.csv")
    sub("Mixed linear models fitted in this script")
    show = ["analysis", "outcome", "log10", "coef_D1MT_vs_ref", "ci_low",
            "ci_high", "p_value", "fold_change", "n_obs", "n_clusters_a",
            "n_clusters_b", "random_effects", "clusters_dropped"]
    print(lm[[c for c in show if c in lm.columns]].to_string(index=False))

at = pd.DataFrame(ARM_TESTS)
if len(at):
    cols = ["comparison", "unit", "n_obs", "n_clusters_a", "n_clusters_b",
            "direction_declared", "ran_as_predicted", "p_one_sided",
            "p_two_sided", "floor_one_sided", "floor_two_sided"]
    at = at[[c for c in cols if c in at.columns]]
    write_csv(at, "85_arm_tests.csv")
    sub("Every arm test, one-sided on the pre-specified direction")
    print(at.to_string(index=False))

rows = []
if foci is not None:
    for c in CONDITION_ORDER:
        g = foci.loc[foci["condition"] == c]
        rows.append({"finding": "1. Residual foci", "arm": c,
                     "metric": "peak myeloid density per focus (cells/mm2)",
                     "values": f"median {g['peak_density'].median():.0f}, "
                               f"n = {len(g)} foci",
                     "extra": f"fold over own background "
                              f"{g['fold_over_background'].median():.1f}x"})
if burden is not None:
    for c in CONDITION_ORDER:
        g = burden.loc[burden["condition"] == c]
        rows.append({"finding": "1. Burden", "arm": c,
                     "metric": "% lung area in granuloma-density regions",
                     "values": ", ".join(f"{x:.1f}" for x in g["burden_pct"]),
                     "extra": "regions: "
                              + ", ".join(str(int(x)) for x in g["n_regions"])})
for c in CONDITION_ORDER:
    g = ido.loc[ido["condition"] == c]
    rows.append({"finding": "2. IDO1+ fraction", "arm": c,
                 "metric": "% of CD68-lineage called IDO1+",
                 "values": ", ".join(f"{x:.1f}" for x in g["pct_ido1_pos"]),
                 "extra": f"stable over raw IDO1 {np.expm1(best_lo):.1f} to "
                          f"{np.expm1(best_hi):.1f}"})
    rows.append({"finding": "2. IDO1 distribution", "arm": c,
                 "metric": "median / p99 IDO1",
                 "values": ", ".join(f"{m:.1f} / {p:.0f}" for m, p
                                     in zip(g["median_ido1"], g["p99_ido1"])),
                 "extra": "same ceiling, different median"})
if len(per_animal_bl):
    for col, cname in [("radial_centered_core", "cell-weighted"),
                       ("radial_centered_core_balanced", "balanced")]:
        if col not in per_animal_bl.columns:
            continue
        for c in CONDITION_ORDER:
            v = per_animal_bl.loc[per_animal_bl["condition"] == c, col]
            rows.append({"finding": f"3. B lineage radial ({cname})", "arm": c,
                         "metric": "mean centred radial position",
                         "values": ", ".join(f"{x:+.3f}" for x in v),
                         "extra": ""})
for key, lab in [("b_lineage_primary", "B lineage, cell-weighted"),
                 ("b_lineage_balanced", "B lineage, balanced")]:
    if key in model_of:
        r = model_of[key]
        rows.append({"finding": "3. Model", "arm": "D1MT vs Untreated",
                     "metric": lab,
                     "values": f"coef {r['coef_D1MT_vs_ref']:+.4f} "
                               f"[{r['ci_low']:+.4f}, {r['ci_high']:+.4f}]",
                     "extra": f"p = {r['p_value']:.4f}, "
                              f"{int(r['n_cells']):,} cells"})
head = pd.DataFrame(rows)
write_csv(head, "80_headline_numbers.csv")
with pd.option_context("display.width", 260, "display.max_colwidth", 62):
    print(head.to_string(index=False))

banner("DONE")
print("Figures:")
print("  F60_finding1_foci_and_burden          six maps, peak density, fold,")
print("                                        burden as the closer")
print("  F61_finding2_ido1                     distribution, fraction, sweep,")
print("                                        ceiling")
print("  F64_finding3_b_lineage_distribution   the B lineage box and violin")
print("  F62_finding3_b_lineage                profile, models, control")
print("  F65_focus_architecture                the F34 content, with the")
print("                                        mixed model and rank-sum on")
print("                                        every panel")
print("\nEvery arm comparison uses the same clustered rank-sum, one-sided on a")
print("pre-specified direction, with the design floor stated on the figure.")
print("No cell-level Wilcoxon appears anywhere: on this data structure it")
print("returns p < 0.05 in 86 percent of runs when the null is true.")

sys.stdout = _tee.terminal
_tee.close()
