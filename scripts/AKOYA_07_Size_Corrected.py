#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - SIZE-CORRECTED DISTANCE MODELS AND POLARISATION
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 07 of the AKOYA analysis series. REVISION 3.4.

WHY THIS SCRIPT EXISTS
    Script 06 measures nearest-neighbour distances in microns. A distance
    measured inside a focus cannot exceed that focus's size, so if the arms
    differ in focus size, raw distance partly restates that difference. This
    script re-tests every distance against outcomes that remove size, and the
    rule is that an effect must survive them to count.

WHAT CHANGED IN REVISION 3 (and why)

    THE TRIGGER: script 04 revision 5 re-baselined the structures and script 06
    revision 4 re-ran against them. Both the premise of this script and two of
    its numbers moved.

    1. INPUT REPOINTED to structures_rev5. Every number this script wrote
       against structures_rev4 is superseded.

    2. THE SIZE CONFOUND IS MATERIALLY SMALLER THAN WHEN THIS SCRIPT WAS
       WRITTEN, and the header numbers were stale.
           equivalent radius   rev4 138 vs 201 um (1.46x)  ->  rev5 250 vs 283 (1.13x)
           inscribed radius    rev4 125 vs 168 um (1.34x)  ->  rev5 176 vs 219 (1.24x)
       Growing foci to 0.25 of peak enlarged both arms and closed the gap. The
       confound is real but weaker, so expect smaller corrections than revision
       2 was built to expect.

    3. THE TWO NORMALISED OUTCOMES ARE DEMOTED TO DIAGNOSTIC. This is the
       substantive change.

       On the revision 4 run the median absolute correlation with focus radius
       was: raw microns 0.10, normalised by equivalent radius 0.36, normalised
       by inscribed radius 0.30, delta versus the focus null 0.11. Both
       normalisations are WORSE at removing size than doing nothing, and all
       twenty normalised correlations were negative.

       That is over-correction, not residual confounding. Dividing distance by
       radius assumes distance scales proportionally with radius. It does not:
       nearest-neighbour distance is set mostly by local target density, which
       barely changes when a focus grows, so the division injects a spurious
       negative dependence larger than the original problem.

       DELTA IS NOW THE SINGLE CONFIRMATORY OUTCOME. The normalisations are
       still computed and still reported, in the size-diagnostic group beside
       raw microns, because the demotion has to remain checkable. Cell 4 prints
       an explicit verdict: if any normalisation comes out nearer zero than
       delta on this input, the demotion is wrong and the script says so.

    4. WHY DELTA, STATED PROPERLY. Delta is the observed distance minus the
       median of a permutation null that shuffles labels INSIDE that structure
       with coordinates and counts held fixed. That null therefore contains two
       confounds at once:
         (a) the structure's size and shape;
         (b) THE DENSITY OF THE TARGET POPULATION. Fewer targets mechanically
             increase nearest-neighbour distance. Treated foci hold 4.9, 6.5 and
             0.2 percent IDO1-positive macrophages against untreated at 13.9,
             18.9 and 11.2, so scarcity alone should push treated distances UP.
             They went DOWN. The observed effect runs against its own bias, and
             the permutation null measures that rather than leaving it as an
             argument.
       Point (b) was not stated in revision 2 and is the stronger of the two.

    5. fit_mixed NOW CHECKS THE STANDARD ERROR, not just the coefficient. Same
       defect and same fix as script 06 revision 3. statsmodels returns a finite
       coefficient with a NaN standard error, raising nothing, when a variance
       component sits on the boundary at zero. Under revision 5 two treated
       animals contribute a single structure each, so this is a live failure
       mode, and a NaN standard error would put a NaN p into table 62 and a
       nonsense interval into F48.

    6. EXACT RANDOMIZATION P-VALUES ALONGSIDE THE MODEL P-VALUES.
       The mixed model stays the ESTIMATOR: every cell contributes to the
       coefficient and interval, and the random effects absorb animal and
       structure baselines. What it cannot supply is a calibrated p. statsmodels
       returns an asymptotic z with no small-sample correction on six clusters;
       a null simulation on this cluster structure previously measured false
       positive rates of 46.7 percent for GLMM variants against 4.4 percent for
       a two-stage animal-means test.

       So the p now also comes from enumerating all C(6,3) = 20 assignments of
       animals to arms. Two versions are reported:
         p_exact_means  permutes the six animal-level means. Instant, exactly
                        calibrated, always computed.
         p_exact_lmm    refits the model under each assignment. Uses ANIMAL-ONLY
                        random effects and a capped cell subsample, because 20
                        refits per model times 10 pairs is 200 nested fits and
                        that is not tractable. The nested fit is retained for
                        the point estimate.
       BOTH ARE FLOORED AT 2/20 = 0.10 TWO-SIDED. That floor is a property of
       three animals per arm and no model escapes it. Where the model p is far
       below 0.10 and the exact p sits at the floor, the gap IS the assumption.

    7. THE DISTANCE AND RADIAL RESULTS ARE NOT INDEPENDENT, AND THE SCRIPT NOW
       SAYS SO. Both nearest-neighbour anchors, IDO1-positive and IDO1-negative
       macrophages, are in the myeloid pool that defines a focus, so they
       concentrate at the density peak, which is the core. Distance from an
       anchor to a non-pool target is therefore substantially "distance from the
       core to that cell", which is radial position in different units. They are
       two views of one spatial fact, not two independent confirmations, and BH
       is no longer applied across them as if they were separate families.

    8. DISTRIBUTIONS AND AN INTERPRETABLE EFFECT ARE REPORTED. Revision 2 gave
       medians and coefficients only, from which the data cannot be seen. Added:
       the empirical cumulative distribution of distance by arm per pair with
       per-animal curves, and the fraction of anchors with a target within fixed
       radii. A coefficient of -52 um means little on its own; "68 percent of
       IDO1-positive macrophages have a plasma cell within 30 um against 35
       percent" is checkable against the images.

    9. EFFECTIVE INDEPENDENT UNITS ARE PRINTED beside every cell count, so no
       table implies that 35,863 cells is 35,863 pieces of information about a
       treatment assigned to six animals.

   10. POLARISATION IS DEMOTED TO A LABELLED APPENDIX. The co-expression result
       was retired when the CD3e by CD20 spillover control separated the arms as
       well or better, and the treated arm now holds four structures in total.
       Kept on record, excluded from every significance count.

   11. A PROVENANCE CHECK ON THE SHARED NULL. This script reads
       null_median_um from script 06 table 53, keyed by structure id. Structure
       ids changed between revision 4 and revision 5, so a stale table would
       half-match silently. The check compares the table's structure keys
       against the loaded cell assignments and refuses to proceed on a
       mismatch.

WHAT CHANGED IN REVISION 3.1 (after the first structures_rev5 run)

    Two reporting defects surfaced in that run. Neither is in the models; both
    are in how a per-animal value is formed.

    A. TABLE 66 DISAGREED WITH THE TEST IT SITS NEXT TO. Three summaries of the
       same quantity were in play for IDO1-positive macrophage to plasma cell:
         script 06 table 54, cell-weighted pooled median : complete separation
         this script's exact test, per-cell mean         : p = 0.100, the 2/20
                                                           floor, reached only
                                                           on complete separation
         this script's table 66, median ACROSS STRUCTURES: not separated
       Two agreed and one did not. The odd one out took an unweighted median
       across structures, so a 28-anchor focus counted the same as a
       3,984-anchor one. The untreated arm has 76 structures, many of them
       small, so that summary was dominated by noisy little foci. It is the
       same functional mismatch corrected in script 06's pooled null, in a
       different place.

       FIX: every per-animal summary is now CELL-WEIGHTED and built from the
       per-cell frame, matching the statistic the exact test permutes and the
       mixed model estimates. PER_ANIMAL_SUMMARY records it in the output.

    B. THE PROXIMITY FRACTIONS HAD NO NULL. Percent-within-a-radius is a raw
       thresholded distance, so on its own it carries the same objection as raw
       microns and, more importantly, no control for target density or local
       cell arrangement. Those are exactly what the delta outcome removes.

       FIX: a permutation null for the proximity fraction. Within each
       structure, labels are shuffled with coordinates and counts held fixed and
       the fraction within each radius is recomputed. Reported as observed, null
       and observed-minus-null, the last being the proximity analogue of delta.
       Observed and null come from the SAME anchors, and both are weighted by
       anchors per structure.

       Read the null-corrected row, not the observed one, before claiming a
       separation. An observed separation that vanishes once the null is
       subtracted was local arrangement, not treatment.

WHAT CHANGED IN REVISION 3.2 (after reading the 3.1 output)

    Nothing in the models moved. Four changes, three to reporting and one new
    diagnostic.

    A. ANCHOR COUNTS ARE PRINTED BEFORE THE PERCENTAGES. Revision 3.1 printed
       proximity fractions with no denominator anywhere on the console. The
       treated arm holds roughly 360 IDO1-positive anchors across three
       animals, and one animal contributes a few dozen of them while carrying
       some of the largest values in the table. n_anchors was already a column
       in table 67; it is now printed next to every percentage, with a flag
       below ANCHOR_CONTRAST_MIN_N. It is a flag and never a filter.

    B. THE PROXIMITY NULL RUNS AT 300 SHUFFLES, up from 100, matching
       N_PERMUTATIONS. Several null-corrected separations turn on margins of
       one to two percentage points, which is inside the Monte Carlo error of
       a 100-shuffle null mean. The null now carries the same precision as the
       distance null it parallels.

    C. TWO REPORTING DEFECTS IN 3.1 ARE FIXED.
       1. The multiple-testing paragraph closed by invoking the neutrophil
          control inverting where necrosis predicts. That number sits in the
          radial_position family, which this same script marks SUPERSEDED and
          tells the reader not to quote. A script may not cite a number it has
          just retracted. The clause now points at script 08, where the
          neutrophil control actually lives, on the centred outcome.
       2. WHICH EXACT P LEADS IS NOW DECLARED. p_exact_means is the headline,
          because treatment was assigned to animals and permuting the six
          animal means is the randomisation that happened. p_exact_lmm is a
          sensitivity analysis: it refits under each assignment, so its
          statistic is scaled by structure-to-structure variation again, and
          the untreated arm carries 76 heterogeneous structures against 4
          treated. That is what inflates the model p, reimported into the test
          built to avoid it. The script also now states that the count of
          tests at the 0.10 floor and the complete-separation list in table 66
          are ONE fact stated twice, because complete separation of the six
          animal means forces p_exact_means to the floor by construction.

    D. NEW DIAGNOSTIC: THE IDO1-POSITIVE MINUS IDO1-NEGATIVE ANCHOR CONTRAST,
       table 68. The 3.1 output showed that the arm separation on proximity is
       carried by the UNTREATED arm rather than the treated one. Treated
       deltas sit near zero, which is where random labelling inside their own
       focus puts them, while untreated deltas run strongly negative. The
       claim is therefore the loss of an exclusion, not the creation of a
       proximity, and "99 percent of treated IDO1-positive macrophages have a
       plasma cell within 50 um" must never be quoted without its null beside
       it, because the matching treated null is 97 percent.

       The contrast makes that testable without leaving the focus. Both
       anchors are in the detection pool, in the same foci, on the same slide,
       in the same animal, and IDO1 enters detection nowhere. The CIRCULARITY
       FLAG below already asserted this; revision 3.2 computes it.

       THE CLAIM MADE HERE IN REVISION 3.2 WAS TOO STRONG AND IS CORRECTED IN
       REVISION 3.3 BELOW. See item C there before quoting anything from
       table 68. It is a diagnostic and receives no q-values.

WHAT CHANGED IN REVISION 3.3 (after reading the 3.2 output)

    The 3.2 output answered the question it was built to answer and, in doing
    so, undermined one of its own claims. Both are recorded here.

    A. THE ANCHOR COUNTS CAME BACK, AND THE THIN ANIMAL IS NOT THE PROBLEM.
       Treated IDO1-positive anchors are 99, 218 and 39, so 356 across the
       arm. That justified the worry. But leave-one-out shows the smallest
       animal is not holding anything up: every contrast separation survives
       dropping ANY single treated animal, and the binding animal is 43111
       throughout, at 218 and 116 anchors. The revision 3.2 flag fired on the
       smallest treated anchor count, which says nothing about whether that
       animal carries the result, so it labelled every separation THIN and
       told us nothing. It is replaced by the leave-one-out itself.

       The counts did surface something else. 43118 holds 39 IDO1-positive
       anchors against 1,041 IDO1-negative, a 27-fold imbalance inside one
       animal, while 43106 and 43111 sit near parity. The IDO1-positive share
       of the macrophage compartment runs 1.4, 1.9, 0.04 treated against 2.6,
       12.4, 1.05 untreated, which does not separate by arm. That reproduces
       the section 4 nuance from an independent direction: the whole-section
       IDO1 fraction separates completely and the within-focus fraction does
       not, and those are different quantities.

    B. TARGET COUNTS ARE NOW PRINTED TOO. Revision 3.2 printed anchors only.
       For a rare target the binding constraint is the other side: the anchors
       are macrophages and run in the thousands while a treated focus may hold
       single-digit Tregs. The permutation null preserves target counts, so the
       arithmetic is honest, but with a handful of targets the null is coarse
       and observed-minus-null is one draw from a short list. Tregs is exactly
       where the contrast looks strongest, so the omission mattered.

    C. THE REVISION 3.2 CLAIM FOR THE CONTRAST WAS TOO STRONG.
       3.2 said the contrast removes focus size, focus density, animal
       baseline and circularity in one construction. That holds only if both
       anchors respond to those the same way, and the output says they do not.
       The IDO1-negative term carries its own large animal-to-animal spread
       that IDO1-positive does not share: for plasma cells at 50 um it runs
       -48.1, -21.7, -15.1 treated and -14.4, -43.1, -58.4 untreated, a
       44-point spread with the arms fully overlapping. It is not a stable
       reference point, it is a second noisy measurement.

       Measured against the IDO1-positive term alone, the contrast helps a
       lot twice, is a wash twice and is worse once:
           CD4- T cells 15 um     1.5  ->  12.8
           plasma cells 30 um     5.2  ->  16.9
           plasma cells 50 um    18.2  ->  19.2
           Helper T cells 50 um   9.1  ->   9.2
           Tregs 50 um           26.5  ->  20.6
       So it is a useful second view, not a confound-cancelling estimator, and
       it must earn its place on every row rather than by construction. The
       IDO1-positive-alone margin is now printed beside every contrast margin.

       It also does not carry to the confirmatory outcome. On distance delta
       the contrast separates only for Helper T cells, at 4.6 um, and Helper T
       is the CD4-call population already discounted. Plasma does not
       separate there, because 43109 comes in at -42.9. Proximity at 30 to 50
       um reads the body of the distribution and distance delta reads a mean
       pulled by its tail, so they measure different things and currently
       disagree. THAT DISAGREEMENT IS UNRESOLVED AND THE CONTRAST DOES NOT GO
       IN A FIGURE UNTIL IT IS.

    D. RANDOM STREAMS ARE SEPARATED, AND THE SUBSAMPLES ARE KEYED ON IDENTITY.
       Revision 3.2 raised N_PROXIMITY_PERMUTATIONS from 100 to 300, a
       parameter that controls only the precision of one null. Numbers moved
       that it has no business touching: observed 15 um proximity on 36463 to
       plasma cells went 9.3 to 8.8, p_model for IDO1+ to CD4- T cells went
       0.345 to 0.232, and p_exact_lmm for IDO1+ to Tregs went 0.333 to 1.000,
       all on identical data. Cause: the proximity permutations, the fallback
       null and the anchor subsample all drew from one generator, so tripling
       the shuffles shifted every later draw including which anchors were
       subsampled in structures above MAX_ANCHORS_PER_STRUCTURE.

       Three independent streams now, and the two subsamples are seeded from
       the structure and anchor identity rather than taken from a running
       stream, so they do not depend on iteration order or on any permutation
       count. EXPECT SMALL MOVEMENTS FROM 3.2 ON THIS RUN because the
       subsamples are drawn differently; after this run the numbers should be
       stable against any parameter that is not supposed to touch them, and
       that is worth checking once by rerunning with a different shuffle count.

       p_exact_means did not move on any test across the two 3.2 runs while
       p_exact_lmm moved between 0.333 and 1.000. That is independent support
       for the revision 3.2 rule making p_exact_means the headline.

WHAT CHANGED IN REVISION 3.4 (the lock-down pass, before script 08)

    No new analysis. This revision closes the gap between what the script
    computes and what it is allowed to say, so results can be frozen.

    A. THE SCRIPT WAS STILL ASSERTING WHAT ITS OWN HEADER RETRACTED.
       Revision 3.3 corrected the anchor-contrast claim in this docstring, but
       the runtime text in the contrast block still printed that differencing
       the two anchors "removes focus size, focus density, animal baseline AND
       the circularity in one construction". Anyone reading the log rather than
       the source got the retracted claim. This is the same defect revision 3.2
       fixed for the neutrophil control, in a second place: a script may not
       assert what it has retracted. The runtime text now carries the
       correction and says why.

       The OUTCOME STATUS banner had the matching problem from the other side.
       It quoted correlation figures from the revision 4 run as if they
       described this input. It no longer quotes any number: Cell 4 measures
       the correlations on the data actually loaded and prints them, which is
       the only place they should come from.

    B. A REPORTABILITY GATE ON TARGET COUNTS, MIN_TARGETS_FOR_REPORT = 50.
       Revision 3.3 printed target counts and they retired a result. Tregs run
       6, 31 and 68 in the three treated animals, and Tregs is where the anchor
       contrast looked strongest. The null preserves target counts so nothing
       is arithmetically wrong, but a null built on six cells is a short list
       and observed-minus-null is one draw from it.

       A count is now a LABEL, never a filter. Everything is still computed,
       written to the tables and printed. A separation whose smallest treated
       target count falls under the threshold is marked NOT REPORTABLE in the
       contrast block and THIN in the table 66 list, so it cannot be lifted
       into a figure by someone reading only the margins. Same philosophy as
       the evidence tiers: run permissive, judge on the evidence side.

    C. THE CROSS-OUTCOME CHECK IS RUN AND PRINTED, NOT LEFT TO THE READER.
       On the 3.3 run the contrast separated for plasma cells on proximity at
       30 and 50 um but not on distance delta, and that survived the random
       stream fix, so it is a property of the two statistics rather than an
       artefact. Proximity reads the body of the distribution at a fixed
       radius; distance delta reads a mean a long tail pulls. The script now
       states per target whether the outcomes AGREE, DISAGREE or separate only
       where the counts are too thin to report. Only AGREE is supported by the
       contrast alone. Script 07 cannot resolve a DISAGREE, because doing so
       needs the distribution shape rather than another summary.

    D. THE 3.3 STABILITY FIX IS NOW VERIFIED INSIDE THE RUN, NOT ASSERTED.
       Two checks. Every anchor subsample is re-derived from its identity key
       and compared against the one used, which is the invariant 3.3
       introduced. And a fingerprint is taken over every OBSERVED quantity,
       the observed medians, the subsample sizes and the observed proximity
       fractions, with no null or delta in it. No parameter that controls only
       a null's precision may move that fingerprint.

       TO VERIFY ACROSS RUNS: change N_PROXIMITY_PERMUTATIONS, rerun, and
       compare one hex string. Identical means the streams are clean. Moved
       means a stream is still crossed and the earlier run cannot be trusted.
       Record the fingerprint beside any number quoted from this script.

CIRCULARITY FLAG
    Foci are defined by pooled myeloid density, so myeloid populations sit
    core-ward by construction. Radial results for CD68+IDO1+, CD68+IDO1-,
    CD163+ and neutrophils are flagged accordingly. What is NOT circular is the
    CONTRAST between IDO1+ and IDO1- macrophages, since both are in the same
    detection pool, nor is the position of any non-pool target.

OUTPUTS
    figures/  F47 .. F52
    tables/   60 .. 68

USAGE
    conda activate sc_pre
    python AKOYA_07_Size_Corrected.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================

IN_DIR = "/master/jlehle/WORKING/AKOYA/structures_rev5"
CELL_DIR = f"{IN_DIR}/cell_assignments"
FOCI_TABLE = f"{IN_DIR}/tables/35_foci_structures_relative.csv"
NN_NULL_TABLE = "/master/jlehle/WORKING/AKOYA/distance_stats/tables/53_nn_per_structure.csv"
OUT_DIR = "/master/jlehle/WORKING/AKOYA/size_corrected"

CONDITION_ORDER = ["D1MT", "Untreated"]
CONDITION_COLORS = {"D1MT": "#2C7FB8", "Untreated": "#D95F02"}
REFERENCE_ARM = "Untreated"

STRUCTURE_SOURCE = "foci"
STRUCT_COL = "focus_id" if STRUCTURE_SOURCE == "foci" else "burden_region_id"

USE_SHARED_NULL = True
# Refuse to run if the shared null table does not describe THESE structures.
# Structure ids changed between rev4 and rev5, so a stale table half-matches.
REQUIRE_NULL_PROVENANCE = True
NULL_PROVENANCE_MIN_MATCH = 0.95      # fraction of structure keys that must match

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

# ---- outcome status ---------------------------------------------------------
# DELTA is the single confirmatory outcome. The normalisations over-corrected on
# the revision 4 run and sit with raw microns in the size-diagnostic group.
# Cell 4 measures this again and says so loudly if the demotion is wrong.
CONFIRMATORY_FAMILIES = ["distance_delta"]
SIZE_DIAGNOSTIC_FAMILIES = ["distance_normalised_equiv",
                            "distance_normalised_inscribed",
                            "distance_raw_diagnostic"]
SUPERSEDED_FAMILIES = ["radial_position"]

# ---- exact randomization ----------------------------------------------------
# All C(6,3) = 20 assignments of animals to arms. Floored at 2/20 = 0.10
# two-sided. That floor is a property of the design.
RUN_EXACT_RANDOMIZATION = True
RUN_EXACT_LMM_REFIT = True          # 20 refits per model; see PERM_* below
PERM_MAX_CELLS = 20000              # fixed subsample reused by all 20 fits
PERM_RE_MODE = "animal only"        # nested fit kept for the point estimate

# ---- interpretable effect ---------------------------------------------------
# Fraction of anchors with a target within these radii, per animal.
PROXIMITY_RADII_UM = [15.0, 30.0, 50.0, 100.0]

# REVISION 3.1: the proximity fraction is a RAW thresholded distance, so on its
# own it inherits the same objection as raw microns. It now also gets a
# permutation null: within each structure, labels are shuffled with coordinates
# and counts held fixed, and the fraction within each radius is recomputed. The
# reported quantity is then observed minus null, which is the proximity
# analogue of the delta outcome and controls for local cell arrangement and
# target density.
RUN_PROXIMITY_NULL = True
N_PROXIMITY_PERMUTATIONS = 300      # ~757 structure-pairs x this many shuffles
# REVISION 3.2: raised from 100 to match N_PERMUTATIONS. Several of the
# null-corrected separations turn on margins of 1 to 2 percentage points,
# which is inside the Monte Carlo error of a 100-shuffle null mean. The null
# is now estimated at the same precision as the distance null it parallels.

# REVISION 3.1: every per-animal summary is CELL-WEIGHTED, matching the
# statistic the exact test permutes and the mixed model estimates. Revision 3
# used an unweighted median across structures in table 66, which let a
# 28-anchor focus count the same as a 3,984-anchor one and disagreed with both
# script 06 table 54 and this script's own exact test. Same functional mismatch
# as the script 06 pooled-null bug, in a different place.
PER_ANIMAL_SUMMARY = "cell-weighted mean"

# ---- REVISION 3.2: the IDO1+ minus IDO1- anchor contrast --------------------
# Both anchors are in the detection pool, sit in the same foci, on the same
# slide, in the same animal, and IDO1 is used nowhere in detection. The
# within-animal difference between the two anchors therefore removes focus
# size, focus density, animal baseline and the core-ward pull of detection in
# one construction, because whatever detection does to a pool member it does
# to both. The header CIRCULARITY FLAG already says this; revision 3.2
# computes it.
#
# NOT a confirmatory family. It is a diagnostic for whether the arm effect is
# specific to the IDO1-positive compartment or generic to myeloid cells, and
# it is only readable where the anchor counts support it. Set False and
# nothing downstream sees it.
RUN_ANCHOR_CONTRAST = True
ANCHOR_CONTRAST_MIN_N = 50    # flag, never filter: anchors per animal

# ---- REVISION 3.4: the reportability gate -----------------------------------
# Revision 3.3 printed target counts and they retired a result: Tregs run 6, 31
# and 68 in the three treated animals. The permutation null preserves target
# counts so the arithmetic is honest, but a null built on six targets is a short
# list and observed-minus-null is one draw from it. Tregs was also where the
# anchor contrast looked strongest, which is exactly how a thin result reaches a
# figure.
#
# A count is now a LABEL on the result, not a filter on the data: everything is
# still computed, written to the tables and printed, and a separation whose
# smallest treated target count falls under this threshold is marked NOT
# REPORTABLE wherever it appears. Same philosophy as the evidence tiers: run
# permissive, judge on the evidence side.
MIN_TARGETS_FOR_REPORT = 50

# ---- REVISION 3.4: the reproducibility fingerprint --------------------------
# 3.3 separated the random streams and keyed the subsamples on identity so that
# a parameter controlling only a null's precision cannot move anything else.
# That was asserted, not verified. Two checks now run inside a single pass:
#   1. every anchor subsample is re-derived from its identity key and compared
#      against the one actually used (the invariant 3.3 introduced);
#   2. a fingerprint over every OBSERVED quantity, which no null parameter may
#      touch. Change N_PROXIMITY_PERMUTATIONS, N_PERMUTATIONS or the proximity
#      radii and rerun: the fingerprint must not move. If it does, a stream is
#      still crossed.
RUN_REPRODUCIBILITY_CHECK = True

# ---- polarisation, APPENDIX ONLY --------------------------------------------
RUN_POLARISATION_APPENDIX = True
POLARISATION_PERCENTILES = [50, 60, 70, 75, 80, 90]
POLARISATION_PRIMARY = 60
MIXING_K = 10
MIXING_MIN_PER_CLASS = 5
MIXING_MIN_MACROPHAGES = 20

# ---- mixed models -----------------------------------------------------------
RUN_MIXED_MODELS = True
MODEL_MAX_CELLS_PER_STRUCTURE = 2000
MODEL_MIN_CELLS_PER_ANIMAL = 20
BH_ALPHA = 0.10

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
import zlib
import hashlib
import warnings
from itertools import combinations
from datetime import datetime
import time

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

# REVISION 3.3: SEPARATE RANDOM STREAMS. Revision 3.2 drew the proximity
# permutations, the fallback null and the anchor subsample from ONE generator,
# so raising N_PROXIMITY_PERMUTATIONS from 100 to 300 shifted every later draw
# and moved numbers it has no business touching. Measured between the two 3.2
# runs: observed 15 um proximity on 36463 to plasma cells went 9.3 to 8.8, and
# p_model for IDO1+ to CD4- T cells went 0.345 to 0.232, on identical data.
# A parameter that only controls the precision of a null must be inert
# everywhere else.
#
# The anchor subsample is now seeded from the STRUCTURE AND ANCHOR IDENTITY
# rather than drawn from a running stream, so it does not depend on iteration
# order, on how many targets were looped, or on any permutation count. crc32
# rather than hash(), because Python salts string hashing per process and the
# subsample would not be reproducible between runs.
rng_prox = np.random.default_rng(RANDOM_SEED + 1)    # proximity null only
rng_null = np.random.default_rng(RANDOM_SEED + 2)    # fallback distance null
rng_take = np.random.default_rng(RANDOM_SEED + 3)    # per-cell model subsample


def subsample_rng(*parts):
    """Deterministic generator keyed on identity, not on stream position."""
    key = "|".join(str(x) for x in parts).encode("utf-8")
    return np.random.default_rng([RANDOM_SEED, int(zlib.crc32(key))])


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
    """BH q-values, NaN-safe, preserving input order."""
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
    outcome ~ arm [+ covariate], random intercept for animal, structure nested
    within animal. Falls back to animal-only random effects on failure.

    REVISION 3: the acceptance test requires a usable STANDARD ERROR, not just a
    finite coefficient. statsmodels returns a finite coefficient with a NaN
    standard error, raising nothing, when a variance component sits on the
    boundary at zero. Under structures_rev5 two treated animals contribute a
    single structure each, so that is a live failure mode here.
    """
    if not HAVE_SM:
        return None
    need = [outcome, "condition", "sample_id", STRUCT_COL]
    if covariate:
        need.append(covariate)
    if [c for c in need if c not in df.columns]:
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

    # REVISION 3 GUARD
    _coef = float(res.params.get("arm", np.nan))
    _se = float(res.bse.get("arm", np.nan))
    _p = float(res.pvalues.get("arm", np.nan))
    if not (np.isfinite(_coef) and np.isfinite(_se) and _se > 0 and np.isfinite(_p)):
        print(f"    NO USABLE SE: {label} / {outcome} ({mode}). "
              f"coef={_coef} se={_se} p={_p}. Row DROPPED, not reported as NaN.")
        return None

    means = d.groupby("condition")["_y"].mean()
    n_struct_arm = d.groupby("condition")["struct_key"].nunique()
    return {
        "analysis": label, "outcome": outcome, "covariate": covariate or "",
        "random_effects": mode, "fit_mode": mode,
        "n_cells": int(len(d)), "n_animals": int(d["sample_id"].nunique()),
        "n_animals_D1MT": int(d.loc[d["condition"] == "D1MT", "sample_id"].nunique()),
        "n_animals_ref": int(d.loc[d["condition"] == REFERENCE_ARM, "sample_id"].nunique()),
        "n_structures": int(d["struct_key"].nunique()),
        "n_structures_D1MT": int(n_struct_arm.get("D1MT", 0)),
        "n_structures_ref": int(n_struct_arm.get(REFERENCE_ARM, 0)),
        "mean_D1MT": float(means.get("D1MT", np.nan)),
        "mean_Untreated": float(means.get("Untreated", np.nan)),
        "coef_D1MT_vs_ref": _coef, "std_err": _se, "p_value": _p,
        "coef_covariate": float(res.params.get("_cov", np.nan)) if covariate else np.nan,
        "p_covariate": float(res.pvalues.get("_cov", np.nan)) if covariate else np.nan,
        "note": note,
    }


# =============================================================================
# EXACT RANDOMIZATION AT THE ANIMAL LEVEL
# =============================================================================
# Treatment was assigned to six animals, three per arm. There are C(6,3) = 20
# assignments, which come in 10 sign-symmetric pairs, so the two-sided p is
# floored at 2/20 = 0.10. That floor is a property of the design and no model
# escapes it: a model p far below it is buying resolution with an assumption.
#
# The mixed model remains the ESTIMATOR. These functions supply only the
# p-value, so every cell still contributes to the coefficient and interval while
# calibration comes from the units treatment was actually assigned to.

def _arm_labels(animals, idx, treated_label="D1MT"):
    return {a: (treated_label if i in idx else REFERENCE_ARM)
            for i, a in enumerate(animals)}


def exact_p_animal_means(d, outcome):
    """
    Permutes the six animal-level means. Instant and exactly calibrated.
    Statistic: difference of arm means of the per-animal means.
    """
    g = d.dropna(subset=[outcome, "sample_id", "condition"])
    if not len(g):
        return np.nan, np.nan, 0
    m = g.groupby(["sample_id", "condition"])[outcome].mean().reset_index()
    animals = list(m["sample_id"])
    vals = m[outcome].to_numpy(float)
    arms = list(m["condition"])
    n_t = sum(1 for a in arms if a != REFERENCE_ARM)
    if n_t == 0 or n_t == len(animals):
        return np.nan, np.nan, 0
    obs_idx = tuple(i for i, a in enumerate(arms) if a != REFERENCE_ARM)
    stats = []
    for idx in combinations(range(len(animals)), n_t):
        t = vals[list(idx)]
        u = vals[[i for i in range(len(animals)) if i not in idx]]
        stats.append(float(np.mean(t) - np.mean(u)))
    obs = float(np.mean(vals[list(obs_idx)])
                - np.mean(vals[[i for i in range(len(animals))
                                if i not in obs_idx]]))
    stats = np.asarray(stats)
    tol = 1e-9 * max(1.0, abs(obs))
    p = float(np.mean(np.abs(stats) >= abs(obs) - tol))
    return p, obs, len(stats)


def exact_p_lmm(d, outcome, covariate=None):
    """
    Refits the model under every assignment. Uses ANIMAL-ONLY random effects and
    a fixed capped subsample, because 20 refits per model across 10 pairs is 200
    nested fits and that is not tractable. The nested fit is retained for the
    point estimate; this supplies only the p.
    """
    if not (HAVE_SM and RUN_EXACT_LMM_REFIT):
        return np.nan, 0
    need = [outcome, "sample_id", "condition"]
    if covariate:
        need.append(covariate)
    g = d.dropna(subset=need).copy()
    if not len(g) or g["sample_id"].nunique() < 4:
        return np.nan, 0
    if len(g) > PERM_MAX_CELLS:
        g = g.sample(PERM_MAX_CELLS, random_state=RANDOM_SEED)
    animals = sorted(g["sample_id"].unique())
    arm_of = {a: g.loc[g["sample_id"] == a, "condition"].iloc[0] for a in animals}
    n_t = sum(1 for a in animals if arm_of[a] != REFERENCE_ARM)
    if n_t == 0 or n_t == len(animals):
        return np.nan, 0
    g["_y"] = pd.to_numeric(g[outcome], errors="coerce")
    if covariate:
        cv = pd.to_numeric(g[covariate], errors="coerce")
        sd = cv.std()
        g["_cov"] = (cv - cv.mean()) / sd if sd and np.isfinite(sd) and sd > 0 else 0.0
    formula = "_y ~ arm" + (" + _cov" if covariate else "")

    def _coef(labels):
        g["arm"] = g["sample_id"].map(
            lambda a: 0.0 if labels[a] == REFERENCE_ARM else 1.0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                r = smf.mixedlm(formula, data=g, groups=g["sample_id"],
                                re_formula="1").fit(reml=True, method="lbfgs",
                                                    maxiter=200)
                v = float(r.params.get("arm", np.nan))
                return v if np.isfinite(v) else np.nan
            except Exception:
                return np.nan

    obs = _coef(arm_of)
    if not np.isfinite(obs):
        return np.nan, 0
    stats = []
    for idx in combinations(range(len(animals)), n_t):
        v = _coef(_arm_labels(animals, idx))
        if np.isfinite(v):
            stats.append(abs(v))
    if not stats:
        return np.nan, 0
    stats = np.asarray(stats)
    tol = 1e-9 * max(1.0, abs(obs))
    return float(np.mean(stats >= abs(obs) - tol)), len(stats)


_tee = Tee(os.path.join(TAB_DIR, "00_size_corrected_report.txt"))
sys.stdout = _tee

banner("AKOYA SIZE-CORRECTED DISTANCE MODELS AND POLARISATION (revision 3.4)")
print(f"Run time      : {datetime.now().isoformat(timespec='seconds')}")
print(f"Input         : {IN_DIR}")
print(f"Shared null   : {NN_NULL_TABLE}")
print(f"Output        : {OUT_DIR}")
if not HAVE_SCIPY:
    print(f"\nERROR: scipy required ({_scipy_err}).")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)
if not HAVE_SM:
    print(f"\n    WARNING: statsmodels unavailable ({_sm_err}). Models skipped.")

print("\nOUTCOME STATUS (revision 3.4)")
print("    CONFIRMATORY : delta, distance minus that focus's own permutation")
print("                   null. Removes structure size AND target density in")
print("                   one construction.")
print("    DIAGNOSTIC   : raw microns, and distance divided by either radius.")
print("                   The normalisations over-correct. Cell 4 MEASURES")
print("                   this on the input actually loaded and prints the")
print("                   numbers; no figures are quoted here, because a")
print("                   banner citing a previous run is how stale numbers")
print("                   get into a paper. Read Cell 4's verdict.")
print("    NOT REPORTABLE: any separation whose treated target count is under")
print(f"                   {MIN_TARGETS_FOR_REPORT}. Labelled, never filtered.")
print("    SUPERSEDED   : raw radial position, by script 08's centred outcome.")

print("\nHOW THE P-VALUES WORK")
print("    The mixed model is the ESTIMATOR: every cell contributes to the")
print("    coefficient and the interval, and the random effects absorb animal")
print("    and structure baselines.")
print("    The p-value ALSO comes from enumerating all C(6,3) = 20 assignments")
print("    of animals to arms, because an asymptotic z on six clusters is not")
print("    calibrated. Both exact p-values are FLOORED AT 0.10 two-sided.")
print("    Where the model p is far below 0.10 and the exact p sits at the")
print("    floor, that gap is the assumption the model is making, not a")
print("    disagreement about the data.")

print("\nNOT INDEPENDENT")
print("    Both nearest-neighbour anchors are in the myeloid pool that defines")
print("    a focus, so they sit at the density peak. Distance from an anchor to")
print("    a non-pool target is largely distance from the core to that cell,")
print("    which is radial position in different units. The distance and radial")
print("    results are two views of one spatial fact, and BH is NOT applied")
print("    across them as separate families.")


# %% Cell 3 - load, geometry, and the shared-null provenance check
# =============================================================================

banner("LOADING")

paths = sorted(glob.glob(os.path.join(CELL_DIR, "*_cell_structures.csv")))
if not paths:
    print(f"ERROR: no cell assignment files in {CELL_DIR}.")
    print("Run script 04 revision 5 first, and check IN_DIR points at rev5.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

cells = {}
for p in paths:
    sid = os.path.basename(p).replace("_cell_structures.csv", "")
    try:
        d = pd.read_csv(p, low_memory=False)
    except Exception as e:
        print(f"    ERROR reading {sid}: {e}. Skipping.")
        continue
    for req in (STRUCT_COL, "region", "radial_pos", "pheno", "condition"):
        if req not in d.columns:
            print(f"    ERROR: {sid} lacks '{req}'. Skipping.")
            d = None
            break
    if d is None:
        continue
    d["sample_id"] = sid
    cells[sid] = d
    n_struct = int(d.loc[d[STRUCT_COL] > 0, STRUCT_COL].nunique())
    print(f"    {sid:<12} {d['condition'].iloc[0]:<10} {len(d):>9,} cells, "
          f"{int((d[STRUCT_COL] > 0).sum()):>8,} in {n_struct} structures")

if not cells:
    print("ERROR: nothing loaded.")
    sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)

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

N_STRUCT_OF = {s: int(cells[s].loc[cells[s][STRUCT_COL] > 0, STRUCT_COL].nunique())
               for s in SAMPLE_ORDER}

sub("EFFECTIVE INDEPENDENT UNITS")
print("    Treatment was assigned to ANIMALS. Cell counts are large; the number")
print("    of independent units is six. Every table below carries n_cells,")
print("    n_structures and n_animals side by side for this reason.\n")
for c in CONDITION_ORDER:
    mem = [s for s in SAMPLE_ORDER if COND_OF[s] == c]
    cellsum = sum(int((cells[s][STRUCT_COL] > 0).sum()) for s in mem)
    print(f"    {c:<12} animals {len(mem)}   structures "
          f"{sum(N_STRUCT_OF[s] for s in mem):>3}   cells in structures "
          f"{cellsum:>9,}")
    print(f"                 per animal: "
          f"{ {short_label(s): N_STRUCT_OF[s] for s in mem} }")
thin = [s for s in SAMPLE_ORDER if N_STRUCT_OF[s] < 2]
if thin:
    print(f"\n    {[short_label(s) for s in thin]} contribute a single structure,")
    print("    so the structure-level variance is estimated almost entirely from")
    print("    the other arm and applied to both. Watch fit_mode in table 62.")

# ---- focus geometry, both radius definitions --------------------------------
radius_equiv, radius_inscribed, shape_ratio = {}, {}, {}
have_inscribed = False
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
        print(f"    {'arm':<12}{'equiv r':>10}{'inscribed r':>14}{'shape':>8}"
              f"{'n':>6}")
        geo = {}
        for c in CONDITION_ORDER:
            g = ft.loc[ft["condition"] == c]
            if len(g):
                geo[c] = (float(g["equiv_radius_um"].median()),
                          float(g["max_inscribed_radius_um"].median()))
                print(f"    {c:<12}{geo[c][0]:>9.0f}u{geo[c][1]:>13.0f}u"
                      f"{g['shape_ratio'].median():>8.2f}{len(g):>6}")
        if len(geo) == 2:
            re_ratio = geo[REFERENCE_ARM][0] / geo["D1MT"][0]
            ri_ratio = geo[REFERENCE_ARM][1] / geo["D1MT"][1]
            print(f"\n    size gap, untreated / treated: equivalent "
                  f"{re_ratio:.2f}x, inscribed {ri_ratio:.2f}x")
            print("    (revision 4 was 1.46x and 1.34x. The confound this script")
            print("     corrects is weaker on this input than when it was written.)")
else:
    print(f"\n    WARNING: {FOCI_TABLE} not found. Radii derived from cell")
    print("    coordinates, which is a rougher estimate.")

# ---- shared per-structure null, WITH PROVENANCE CHECK -----------------------
SHARED_NULL = {}
null_ok = True
if USE_SHARED_NULL and os.path.exists(NN_NULL_TABLE):
    nt = pd.read_csv(NN_NULL_TABLE)
    need = ["sample_id", "structure_id", "anchor", "target", "null_median_um"]
    if all(c in nt.columns for c in need):
        # PROVENANCE: structure ids changed between rev4 and rev5, so a stale
        # table half-matches silently. Compare its structure keys against the
        # cell assignments actually loaded.
        have_keys = set()
        for s in SAMPLE_ORDER:
            d = cells[s]
            for k in d.loc[d[STRUCT_COL] > 0, STRUCT_COL].unique():
                have_keys.add((s, int(k)))
        tab_keys = set(zip(nt["sample_id"], nt["structure_id"].astype(int)))
        matched = tab_keys & have_keys
        frac = len(matched) / max(1, len(tab_keys))
        print(f"\n    shared null provenance: {len(matched)} of {len(tab_keys)} "
              f"structure keys in table 53 match the loaded structures "
              f"({100 * frac:.1f}%)")
        if REQUIRE_NULL_PROVENANCE and frac < NULL_PROVENANCE_MIN_MATCH:
            null_ok = False
            print("\n    ERROR: the shared null table does not describe these")
            print("    structures. Structure ids changed between structures_rev4")
            print("    and structures_rev5, so a stale table would half-match and")
            print("    silently attach the wrong null to the wrong focus.")
            print(f"    RE-RUN SCRIPT 06 against {IN_DIR} first.")
            sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)
        for _, r in nt.iterrows():
            SHARED_NULL[(r["sample_id"], int(r["structure_id"]),
                         r["anchor"], r["target"])] = float(r["null_median_um"])
        print(f"    loaded {len(SHARED_NULL)} shared per-structure nulls. The")
        print("    delta outcome is therefore identical in scripts 06, 07 and 08.")
    else:
        print(f"\n    WARNING: {NN_NULL_TABLE} lacks expected columns; nulls will")
        print("    be recomputed here and will NOT match scripts 06 and 08.")
elif USE_SHARED_NULL:
    print(f"\n    WARNING: {NN_NULL_TABLE} not found. Run script 06 first, or")
    print("    nulls are recomputed here and will not match 06 and 08.")


def structure_frames():
    for s in SAMPLE_ORDER:
        d = cells[s]
        for k, g in d.loc[d[STRUCT_COL] > 0].groupby(STRUCT_COL):
            yield s, int(k), g


def focus_radii(sid, k, g):
    """(equivalent radius, inscribed radius). NaN where unavailable."""
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


# %% Cell 4 - A: the size confound, and whether each correction removes it
# =============================================================================

banner("A - SIZE CONFOUND DIAGNOSTIC")

print("    Observed nearest-neighbour distance against focus radius. If distance")
print("    tracks radius, raw distance cannot be compared between arms whose")
print("    foci differ in size. The outcome to trust is the one whose")
print("    correlation with radius is NEAREST ZERO. A correlation that flips")
print("    sign and grows is over-correction, which is worse than the original")
print("    confound.\n")

percell_frames, struct_rows = [], []
# REVISION 3.4: reproducibility check. _fp collects only OBSERVED quantities,
# which no null parameter may influence; _subsample_mismatches counts any case
# where an anchor subsample could not be re-derived from its identity key.
_fp = hashlib.sha256()
_subsample_mismatches = 0
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
        a_use = (subsample_rng(s, k, anchor).choice(
                     a_idx, size=MAX_ANCHORS_PER_STRUCTURE, replace=False)
                 if n_a > MAX_ANCHORS_PER_STRUCTURE else a_idx)
        if RUN_REPRODUCIBILITY_CHECK and n_a > MAX_ANCHORS_PER_STRUCTURE:
            # re-derive from the identity key alone. Equal means the subsample
            # does not depend on stream position, which is the whole point of
            # the revision 3.3 change.
            _again = subsample_rng(s, k, anchor).choice(
                a_idx, size=MAX_ANCHORS_PER_STRUCTURE, replace=False)
            if not np.array_equal(np.sort(a_use), np.sort(_again)):
                _subsample_mismatches += 1
        for target in NN_TARGETS:
            if target == anchor:
                continue
            t_idx = np.flatnonzero(ph == target)
            n_t = len(t_idx)
            if n_t == 0:
                continue
            obs_v = nn_distances(xy[a_use], xy[t_idx])
            obs = float(np.median(obs_v)) if len(obs_v) else np.nan

            # ---- proximity fractions, observed and against a null ----------
            prox_obs = {R: (float(np.mean(obs_v <= R)) if len(obs_v) else np.nan)
                        for R in PROXIMITY_RADII_UM}
            prox_null = {R: np.nan for R in PROXIMITY_RADII_UM}
            if RUN_PROXIMITY_NULL and len(obs_v):
                acc = {R: [] for R in PROXIMITY_RADII_UM}
                for _ in range(N_PROXIMITY_PERMUTATIONS):
                    perm = rng_prox.permutation(n)
                    fa, ft_ = perm[:n_a], perm[n_a:n_a + n_t]
                    fa_use = (rng_prox.choice(fa, size=MAX_ANCHORS_PER_STRUCTURE,
                                              replace=False)
                              if len(fa) > MAX_ANCHORS_PER_STRUCTURE else fa)
                    v = nn_distances(xy[fa_use], xy[ft_])
                    if len(v):
                        for R in PROXIMITY_RADII_UM:
                            acc[R].append(float(np.mean(v <= R)))
                for R in PROXIMITY_RADII_UM:
                    if acc[R]:
                        prox_null[R] = float(np.mean(acc[R]))

            shared = SHARED_NULL.get((s, int(k), anchor, target), np.nan)
            if np.isfinite(shared):
                n_shared += 1
                mu, pval, z, src = shared, np.nan, np.nan, "script_06"
            else:
                n_recomputed += 1
                null = np.empty(N_PERMUTATIONS)
                for it in range(N_PERMUTATIONS):
                    perm = rng_null.permutation(n)
                    fa, ft_ = perm[:n_a], perm[n_a:n_a + n_t]
                    fa_use = (rng_null.choice(fa, size=MAX_ANCHORS_PER_STRUCTURE,
                                              replace=False)
                              if len(fa) > MAX_ANCHORS_PER_STRUCTURE else fa)
                    v = nn_distances(xy[fa_use], xy[ft_])
                    null[it] = np.median(v) if len(v) else np.nan
                pval, z, mu = empirical_p(obs, null)
                src = "recomputed"

            if RUN_REPRODUCIBILITY_CHECK:
                # observed only: the median, the anchor subsample, and each
                # observed proximity fraction. Never a null or a delta.
                _fp.update(f"{s}|{k}|{anchor}|{target}|{n_a}|{len(a_use)}|"
                           f"{n_t}|{obs:.10g}".encode("utf-8"))
                for R in PROXIMITY_RADII_UM:
                    _fp.update(f"|{R}:{prox_obs[R]:.10g}".encode("utf-8"))

            rec_prox = {}
            for R in PROXIMITY_RADII_UM:
                rec_prox[f"obs_frac_{int(R)}um"] = prox_obs[R]
                rec_prox[f"null_frac_{int(R)}um"] = prox_null[R]
                rec_prox[f"delta_frac_{int(R)}um"] = (
                    prox_obs[R] - prox_null[R]
                    if np.isfinite(prox_obs[R]) and np.isfinite(prox_null[R])
                    else np.nan)

            struct_rows.append({
                "sample_id": s, "condition": COND_OF[s], "structure_id": k,
                "anchor": anchor, "target": target,
                "n_anchor": n_a, "n_anchors_used": int(len(a_use)),
                "n_target": n_t, "n_structure": n, **rec_prox,
                "focus_radius_um": r_eq, "focus_inscribed_radius_um": r_in,
                "shape_ratio": shape_ratio.get((s, int(k)), np.nan),
                "observed_median_um": obs, "null_median_um": mu,
                "null_source": src,
                "delta_um": obs - mu if np.isfinite(mu) else np.nan,
                "normalised_median_equiv": obs / r_eq if (r_eq and np.isfinite(r_eq)) else np.nan,
                "normalised_median_inscribed": obs / r_in if (r_in and np.isfinite(r_in)) else np.nan,
                "z": z, "p_empirical": pval,
                "tier": tier_of(min(n_a, n_t)),
            })

            take = obs_v
            if len(take) > MODEL_MAX_CELLS_PER_STRUCTURE:
                take = subsample_rng(s, k, anchor, target, "take").choice(
                    take, size=MODEL_MAX_CELLS_PER_STRUCTURE,
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

print(f"    nulls: {n_shared} taken from script 06 table 53, "
      f"{n_recomputed} recomputed here")

# ---- REVISION 3.4: reproducibility verdict ----------------------------------
FINGERPRINT = _fp.hexdigest()[:16] if RUN_REPRODUCIBILITY_CHECK else "off"
if RUN_REPRODUCIBILITY_CHECK:
    sub("Reproducibility check")
    print("    Revision 3.2 raised N_PROXIMITY_PERMUTATIONS from 100 to 300, a")
    print("    parameter controlling only the precision of one null, and moved")
    print("    observed proximity fractions and model p-values with it. The")
    print("    cause was one shared generator. Revision 3.3 split the streams")
    print("    and keyed the subsamples on structure and anchor identity.")
    print("    These two checks verify that rather than asserting it.\n")
    if _subsample_mismatches:
        print(f"    FAIL: {_subsample_mismatches} anchor subsample(s) could not be")
        print("    re-derived from their identity key. A stream is still")
        print("    crossed and NOTHING IN THIS RUN IS REPRODUCIBLE. Stop here.")
    else:
        print("    PASS  every anchor subsample re-derives from its identity")
        print("          key alone, so it does not depend on stream position.")
    print(f"\n    OBSERVED-VALUE FINGERPRINT : {FINGERPRINT}")
    print("    Covers every observed median, anchor subsample size and observed")
    print("    proximity fraction. No null parameter may touch any of them.")
    print("    TO VERIFY: change N_PROXIMITY_PERMUTATIONS (300 -> 200, say) and")
    print("    rerun. This fingerprint must be IDENTICAL. If it moves, a stream")
    print("    is crossed and the run before it cannot be trusted. Record the")
    print("    fingerprint beside any number quoted from this script.")

nn_struct = pd.DataFrame(struct_rows)
percell = (pd.concat(percell_frames, ignore_index=True)
           if percell_frames else pd.DataFrame())
corrs = pd.DataFrame()
DEMOTION_SUPPORTED = True

if len(nn_struct):
    write_csv(nn_struct, "60_nn_per_structure_with_geometry.csv")

    sub("Correlation of each distance outcome with focus radius, per pair")
    print(f"    {'pair':<42}{'raw':>8}{'norm eq':>9}{'norm in':>9}{'delta':>8}")
    print("    " + "-" * 76)
    corr_rows = []
    for (a, t), g in nn_struct.groupby(["anchor", "target"]):
        rec = {"anchor": a, "target": t, "n_structures": int(len(g)),
               "corr_raw_vs_radius": safe_corr(g["observed_median_um"],
                                               g["focus_radius_um"]),
               "corr_norm_equiv_vs_radius": safe_corr(g["normalised_median_equiv"],
                                                      g["focus_radius_um"]),
               "corr_norm_inscribed_vs_radius": safe_corr(
                   g["normalised_median_inscribed"], g["focus_radius_um"]),
               "corr_delta_vs_radius": safe_corr(g["delta_um"],
                                                 g["focus_radius_um"])}
        corr_rows.append(rec)
        print(f"    {(a[:20] + ' -> ' + t)[:41]:<42}"
              f"{rec['corr_raw_vs_radius']:>+8.2f}"
              f"{rec['corr_norm_equiv_vs_radius']:>+9.2f}"
              f"{rec['corr_norm_inscribed_vs_radius']:>+9.2f}"
              f"{rec['corr_delta_vs_radius']:>+8.2f}")
    corrs = pd.DataFrame(corr_rows)
    write_csv(corrs, "61_distance_radius_correlation.csv")

    sub("VERDICT: does the demotion of the normalisations hold on this input?")
    summ = {}
    for col, name in [("corr_raw_vs_radius", "raw microns"),
                      ("corr_norm_equiv_vs_radius", "normalised, equivalent r"),
                      ("corr_norm_inscribed_vs_radius", "normalised, inscribed r"),
                      ("corr_delta_vs_radius", "delta vs focus null")]:
        v = corrs[col].dropna()
        if not len(v):
            continue
        summ[name] = float(v.abs().median())
        n_neg = int((v < 0).sum())
        print(f"    {name:<28} median |r| = {v.abs().median():.2f}   "
              f"range {v.min():+.2f} to {v.max():+.2f}   "
              f"{n_neg} of {len(v)} negative")
    d_med = summ.get("delta vs focus null", np.nan)
    better = [k for k, v in summ.items()
              if k.startswith("normalised") and np.isfinite(d_med) and v < d_med]
    print()
    if better:
        DEMOTION_SUPPORTED = False
        print(f"    OBJECTION: {better} correlate with radius LESS than delta")
        print("    does on this input. The revision 3 demotion was based on the")
        print("    revision 4 run and does not hold here. Reconsider before")
        print("    treating delta as the sole confirmatory outcome.")
    else:
        print("    Demotion holds. Delta correlates with radius as well as or")
        print("    better than every normalisation, and the normalisations")
        print("    over-correct in the direction seen on revision 4.")

    # ---- F47 the confound --------------------------------------------------
    panels = [("observed_median_um", "Raw distance (diagnostic)"),
              ("normalised_median_equiv", "/ equivalent radius (diagnostic)"),
              ("normalised_median_inscribed", "/ inscribed radius (diagnostic)"),
              ("delta_um", "Delta vs focus null (CONFIRMATORY)")]
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
        ax.set_xlabel("focus equivalent radius (um)", fontsize=FONT_SIZE_BASE - 12)
        ax.set_ylabel(ttl, fontsize=FONT_SIZE_BASE - 14)
        ax.set_title(ttl, fontsize=FONT_SIZE_TITLE - 16,
                     color=OK_COLOR if "CONFIRMATORY" in ttl else TEXT_COLOR)
        style_axes(ax)
    axes[0].legend(handles=[Line2D([0], [0], color=COLOR_OF[s],
                                   marker=MARKER_OF[s], markersize=16,
                                   linestyle="none",
                                   label=f"{short_label(s)} ({COND_OF[s]})")
                            for s in SAMPLE_ORDER],
                   frameon=False, fontsize=FONT_SIZE_LEGEND - 16, loc="best")
    fig.suptitle("Which correction actually removes focus size?\n"
                 "Flat is good. A trend that flips sign and grows is "
                 "over-correction, which is worse than the raw confound",
                 y=1.05, fontsize=FONT_SIZE_TITLE - 8)
    save_fig(fig, "F47_size_confound_diagnostic")


# %% Cell 5 - B/C: distance models, with exact randomization p-values
# =============================================================================

banner("B/C - DISTANCE MODELS, MIXED MODEL ESTIMATE + EXACT RANDOMIZATION P")

print("    The mixed model is the ESTIMATOR. Every cell contributes to the")
print("    coefficient and interval; the random effects absorb animal and")
print("    structure baselines rather than pretending cells are independent")
print("    draws from one pool.")
print()
print("    The p-value also comes from enumerating all 20 assignments of the six")
print("    animals to arms, because an asymptotic z on six clusters is not")
print("    calibrated. p_exact_means permutes the animal-level means and is")
print("    instant. p_exact_lmm refits the model under each assignment using")
print("    animal-only random effects and a capped subsample, because 20 refits")
print("    per model is not tractable with the nested term.")
print()
print("    BOTH EXACT P-VALUES ARE FLOORED AT 0.10. Where the model p is far")
print("    below that and the exact p sits at the floor, the gap is the")
print("    assumption, not a disagreement about the data.\n")

model_rows = []
n_dropped_no_se = 0
_t0 = time.time()

# RUNTIME. The nested animal-plus-structure fit is the expensive step: 10 pairs
# times 4 outcomes is 40 fits over tens of thousands of rows with one dummy per
# structure. The exact LMM refit adds 20 animal-only fits per confirmatory
# model on a capped subsample. Both are timed below. If the total is painful,
# set RUN_EXACT_LMM_REFIT = False: p_exact_means is instant, exactly calibrated,
# and carries the same 0.10 floor, so nothing essential is lost.

OUTCOME_SPEC = [
    ("distance_delta_um", "distance_delta", None,
     "distance minus focus permutation null (shared with 06 and 08). "
     "Removes structure size AND target density."),
    ("distance_norm_equiv", "distance_normalised_equiv", None,
     "DIAGNOSTIC: over-corrected on rev4, see Cell 4"),
    ("distance_norm_inscribed", "distance_normalised_inscribed", None,
     "DIAGNOSTIC: over-corrected on rev4, see Cell 4"),
    ("distance_um", "distance_raw_diagnostic", "focus_radius_um",
     "DIAGNOSTIC ONLY: bounded by focus size, radius nearly collinear with arm"),
]

if RUN_MIXED_MODELS and HAVE_SM and len(percell):
    for (a, t), g in percell.groupby(["anchor", "target"]):
        for outcome, fam, cov, note in OUTCOME_SPEC:
            if outcome not in g.columns or not np.isfinite(g[outcome]).any():
                continue
            r = fit_mixed(g, outcome, f"{a} -> {t}", covariate=cov, note=note)
            if not r:
                n_dropped_no_se += 1
                continue
            r["family"] = fam
            r["anchor"], r["target"] = a, t
            r["superseded_by_script_08"] = False
            if RUN_EXACT_RANDOMIZATION and fam in CONFIRMATORY_FAMILIES:
                pm, obs_m, n_assign = exact_p_animal_means(g, outcome)
                r["p_exact_means"] = pm
                r["stat_exact_means"] = obs_m
                r["n_assignments"] = n_assign
                r["p_exact_floor"] = 2.0 / n_assign if n_assign else np.nan
                pl, n_l = exact_p_lmm(g, outcome, covariate=cov)
                r["p_exact_lmm"] = pl
                r["n_assignments_lmm"] = n_l
            r["fit_seconds"] = round(time.time() - _t0, 1)
            _t0 = time.time()
            model_rows.append(r)
            if fam in CONFIRMATORY_FAMILIES:
                print(f"    {a[:18]:<20} -> {t:<16} "
                      f"coef={r['coef_D1MT_vs_ref']:>+8.2f}  "
                      f"p_model={r['p_value']:.4f}  "
                      f"p_exact_means={r.get('p_exact_means', np.nan):.3f}  "
                      f"p_exact_lmm={r.get('p_exact_lmm', np.nan):.3f}  "
                      f"[{r['fit_mode']}, {r['fit_seconds']:.0f}s]")

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
            frames += [gg.sample(min(len(gg), MODEL_MAX_CELLS_PER_STRUCTURE),
                                 random_state=RANDOM_SEED)
                       for _, gg in sel.groupby(STRUCT_COL)]
        if not frames:
            continue
        dd = pd.concat(frames, ignore_index=True)
        circ = p in DETECTION_POOL
        note = "raw radial_pos, superseded by script 08's centred outcome"
        if circ:
            note += "; CIRCULAR: population defines the foci"
        r = fit_mixed(dd, "radial_pos", f"radial: {p}", note=note)
        if not r:
            n_dropped_no_se += 1
            continue
        r["family"] = "radial_position"
        r["phenotype"] = p
        r["circular"] = circ
        r["superseded_by_script_08"] = True
        model_rows.append(r)

models = pd.DataFrame(model_rows)

if n_dropped_no_se:
    sub("MODELS DROPPED")
    print(f"    {n_dropped_no_se} fit(s) returned no usable standard error or")
    print("    could not be fitted, and were dropped rather than written as NaN.")
    print("    Revision 2 would have written them into table 62 with a NaN")
    print("    p-value and a nonsense interval in F48.")

if len(models):
    # BH only within the confirmatory family. Diagnostic and superseded
    # families get no q-values, and the distance and radial results are NOT
    # pooled into one correction because they are not independent tests.
    models["q_value"] = np.nan
    for fam, g in models.groupby("family"):
        if fam not in CONFIRMATORY_FAMILIES:
            continue
        models.loc[g.index, "q_value"] = benjamini_hochberg(g["p_value"].to_numpy())
    models["sig_p05"] = models["p_value"] < 0.05
    models["sig_q10"] = models["q_value"] < BH_ALPHA
    models["is_confirmatory"] = models["family"].isin(CONFIRMATORY_FAMILIES)
    write_csv(models, "62_mixed_models_all_outcomes.csv")

    for fam in CONFIRMATORY_FAMILIES + SIZE_DIAGNOSTIC_FAMILIES + SUPERSEDED_FAMILIES:
        g = models.loc[models["family"] == fam].sort_values("p_value")
        if not len(g):
            continue
        if fam in CONFIRMATORY_FAMILIES:
            flag = "   [CONFIRMATORY]"
        elif fam in SUPERSEDED_FAMILIES:
            flag = "   [SUPERSEDED by script 08, do not quote]"
        else:
            flag = "   [DIAGNOSTIC ONLY, no q-values]"
        sub(f"Family: {fam}  (n = {len(g)} tests){flag}")
        header = (f"    {'test':<40}{'coef':>10}{'p_model':>9}{'q(BH)':>8}"
                  f"{'p_exact':>9}{'cells':>9}{'str':>5}{'anim':>5}")
        print(header)
        for _, r in g.iterrows():
            mark = " *" if r["sig_q10"] else ""
            circ = " [circular]" if r.get("circular", False) is True else ""
            qtxt = "      na" if pd.isna(r["q_value"]) else f"{r['q_value']:>8.4f}"
            ex = r.get("p_exact_means", np.nan)
            extxt = "       na" if pd.isna(ex) else f"{ex:>9.3f}"
            print(f"    {r['analysis'][:39]:<40}{r['coef_D1MT_vs_ref']:>+10.3f}"
                  f"{r['p_value']:>9.4f}{qtxt}{extxt}"
                  f"{r['n_cells']:>9,}{r['n_structures']:>5}{r['n_animals']:>5}"
                  f"{mark}{circ}")

    conf = models.loc[models["is_confirmatory"]]
    if len(conf) and "p_exact_means" in conf.columns:
        sub("MODEL P AGAINST EXACT P, the confirmatory family")
        print("    The exact p is floored at 0.10 two-sided. A model p far below")
        print("    that is not extra evidence; it is the asymptotic assumption.\n")
        print(f"    {'test':<40}{'p_model':>10}{'p_means':>10}{'p_lmm':>10}"
              f"{'ratio':>10}")
        for _, r in conf.sort_values("p_value").iterrows():
            pm = r.get("p_exact_means", np.nan)
            pl = r.get("p_exact_lmm", np.nan)
            ratio = (pm / r["p_value"]) if (np.isfinite(pm)
                                            and r["p_value"] > 0) else np.nan
            print(f"    {r['analysis'][:39]:<40}{r['p_value']:>10.4f}"
                  f"{pm:>10.3f}{pl:>10.3f}{ratio:>10.1f}")
        at_floor = int((conf["p_exact_means"] <= 0.101).sum())
        print(f"\n    {at_floor} of {len(conf)} confirmatory tests sit at the")
        print("    exact floor of 0.10, which is the strongest result three")
        print("    animals per arm can produce. Lead with the effect size and")
        print("    the per-animal separation, not with the model p.")
        print()
        print("    WHICH EXACT P LEADS, DECLARED (revision 3.2).")
        print("    p_exact_means is the headline. Treatment was assigned to")
        print("    ANIMALS, so permuting the six animal-level means is the")
        print("    randomisation that actually happened and it needs no model.")
        print("    p_exact_lmm is reported as a sensitivity analysis only. It")
        print("    refits under each assignment, so its statistic is again")
        print("    scaled by structure-to-structure variation, and the untreated")
        print("    arm carries 76 heterogeneous structures against 4 treated.")
        print("    That is the same thing that inflates the model p, reimported")
        print("    into a test built to avoid it. Where the two disagree,")
        print("    p_exact_means is quoted and the disagreement is reported.")
        print()
        print("    AND IT IS NOT A SECOND RESULT. Complete separation of the six")
        print("    animal means FORCES p_exact_means to 0.10 by construction, so")
        print("    the count above and the complete-separation list in table 66")
        print("    are one fact stated twice. Never present them as two.")

    sub("Multiple testing position")
    n_p05 = int(conf["sig_p05"].sum()) if len(conf) else 0
    n_q10 = int(conf["sig_q10"].sum()) if len(conf) else 0
    print(f"    Confirmatory family only ({len(conf)} tests, delta outcome):")
    print(f"      {n_p05} have p_model < 0.05")
    print(f"      {n_q10} survive BH at q < {BH_ALPHA}")
    print("\n    BH is applied ONLY within the delta family. The distance and")
    print("    radial results are two views of one spatial fact, not independent")
    print("    tests, so they are not pooled into a single correction.")
    if n_q10 == 0:
        print("\n    Nothing survives BH. That does not invalidate the work, but")
        print("    no single test here is a standalone claim. The argument rests")
        print("    on coherence: consistent direction across populations and")
        print("    animals, plus the IDO1-positive minus IDO1-negative anchor")
        print("    contrast below, which is internal to each focus.")
        print()
        print("    REVISION 3.2 CORRECTION. Earlier revisions closed this")
        print("    paragraph by invoking the neutrophil control inverting where")
        print("    necrosis predicts. That number lives in the radial_position")
        print("    family, which THIS SCRIPT marks SUPERSEDED and tells the")
        print("    reader not to quote. It cannot be cited here. The neutrophil")
        print("    control is script 08's, on the centred outcome, and the")
        print("    coherence argument may only lean on it once 08 has run.")

    # ---- F48 forest, confirmatory family ------------------------------------
    g = models.loc[models["family"].isin(CONFIRMATORY_FAMILIES)].sort_values(
        "coef_D1MT_vs_ref")
    if len(g):
        fig, ax = plt.subplots(figsize=(22, max(12, 1.1 * len(g))))
        yy = np.arange(len(g))
        lo = g["coef_D1MT_vs_ref"] - 1.96 * g["std_err"]
        hi = g["coef_D1MT_vs_ref"] + 1.96 * g["std_err"]
        for i, (_, r) in enumerate(g.iterrows()):
            ex = r.get("p_exact_means", np.nan)
            col = (FLAG_COLOR if (np.isfinite(ex) and ex <= 0.101)
                   else OK_COLOR if r["sig_p05"] else "#777777")
            ax.plot([lo.iloc[i], hi.iloc[i]], [i, i], color=col, linewidth=6,
                    zorder=2)
            ax.scatter([r["coef_D1MT_vs_ref"]], [i], s=420, color=col,
                       edgecolor="#FFFFFF", linewidth=2, zorder=3)
            ax.text(hi.iloc[i], i,
                    f"  p={r['p_value']:.3f}  exact={ex:.2f}  n={r['n_cells']:,}"
                    f" / {r['n_structures']} str / {r['n_animals']} animals",
                    va="center", fontsize=FONT_SIZE_ANNOT - 12)
        ax.axvline(0, color="#000000", linewidth=3.5)
        ax.set_yticks(yy)
        ax.set_yticklabels([r["analysis"][:38] for _, r in g.iterrows()],
                           fontsize=FONT_SIZE_TICK - 12)
        ax.set_xlabel("D1MT effect vs Untreated, delta um (95% CI)")
        ax.set_title("Distance treatment effects, confirmatory outcome only\n"
                     "red = exact randomization p at its 0.10 floor",
                     fontsize=FONT_SIZE_TITLE - 10)
        ax.xaxis.grid(True, color=GRID_COLOR, linewidth=1.5)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        save_fig(fig, "F48_model_forest_confirmatory")
else:
    models = pd.DataFrame()
    print("    SKIPPED (statsmodels unavailable or no per-cell data)")


# %% Cell 6 - the data itself: distributions and an interpretable effect
# =============================================================================

banner("DISTRIBUTIONS AND INTERPRETABLE EFFECT SIZE")

print("    A coefficient of -52 um means little on its own. Two things are")
print("    added here so the data can be seen rather than summarised:")
print("      F51  the empirical cumulative distribution of distance by arm, one")
print("           curve per animal, for every pair.")
print("      Table 67  the fraction of anchors with a target within fixed radii,")
print("           per animal. That is checkable against the images.\n")

# Aggregated from nn_struct so the observed and null fractions come from the
# SAME anchors, and weighted by anchors per structure so a 28-anchor focus does
# not count the same as a 3,984-anchor one.
prox_rows = []
if len(nn_struct):
    for (a, t, s), g in nn_struct.groupby(["anchor", "target", "sample_id"]):
        w = g["n_anchors_used"].to_numpy(float)
        if not w.sum():
            continue
        rec = {"anchor": a, "target": t, "sample_id": s,
               "animal_id": short_label(s), "condition": COND_OF[s],
               "n_anchors": int(w.sum()),
               # REVISION 3.3: the TARGET count is the binding constraint on
               # rare populations. Tregs are the clearest case: the anchors are
               # macrophages and run in the thousands, while a treated focus may
               # hold single-digit Tregs. The permutation null preserves target
               # counts so the arithmetic is honest, but with a handful of
               # targets the null is coarse and the observed value is one draw.
               # Revision 3.2 printed anchors only, which hid exactly that.
               "n_targets": int(g["n_target"].sum()),
               "n_targets_min_structure": int(g["n_target"].min()),
               "n_structures": int(len(g)),
               "summary": PER_ANIMAL_SUMMARY}
        pv = percell.loc[(percell["anchor"] == a) & (percell["target"] == t)
                         & (percell["sample_id"] == s), "distance_um"]
        pv = pv.dropna().to_numpy(float)
        rec["median_um"] = float(np.median(pv)) if len(pv) else np.nan
        rec["p25_um"] = float(np.percentile(pv, 25)) if len(pv) else np.nan
        rec["p75_um"] = float(np.percentile(pv, 75)) if len(pv) else np.nan

        def _wm(col):
            v = g[col].to_numpy(float)
            ok = np.isfinite(v)
            return float(np.average(v[ok], weights=w[ok])) if ok.any() else np.nan

        for R in PROXIMITY_RADII_UM:
            o = _wm(f"obs_frac_{int(R)}um")
            nl = _wm(f"null_frac_{int(R)}um")
            rec[f"pct_within_{int(R)}um"] = 100.0 * o if np.isfinite(o) else np.nan
            rec[f"null_pct_within_{int(R)}um"] = 100.0 * nl if np.isfinite(nl) else np.nan
            rec[f"delta_pct_within_{int(R)}um"] = (
                100.0 * (o - nl) if np.isfinite(o) and np.isfinite(nl) else np.nan)
        prox_rows.append(rec)

prox = pd.DataFrame(prox_rows)
if len(prox):
    write_csv(prox, "67_proximity_fractions.csv")
    def _sep_list(col):
        out = []
        for (a, t), g in prox.groupby(["anchor", "target"]):
            tv = g.loc[g["condition"] == "D1MT", col].dropna()
            uv = g.loc[g["condition"] == REFERENCE_ARM, col].dropna()
            if len(tv) and len(uv) and (tv.max() < uv.min() or uv.max() < tv.min()):
                out.append(f"{a[:18]} -> {t}")
        return out

    # REVISION 3.2: anchor counts printed BEFORE the percentages. A percentage
    # computed on 39 anchors and one computed on 21,000 are not the same
    # measurement, and the treated arm is thin. These are already columns in
    # table 67; revision 3.1 simply did not print them.
    sub("Cell counts behind every percentage below")
    print("    A proximity fraction is a percentage OF THE ANCHORS, computed")
    print("    against a null that preserves the TARGET count. Both matter and")
    print("    they fail differently. Too few anchors and one cluster of them")
    print("    carries a per-animal value. Too few targets and the null is")
    print("    coarse, so observed minus null is one draw from a short list.")
    print("    Read both tables before any percentage.\n")
    order = [short_label(s) for s in SAMPLE_ORDER]
    for col, what in [("n_anchors", "ANCHORS (the denominator)"),
                      ("n_targets", "TARGETS (what the null preserves)")]:
        npiv = prox.pivot_table(index=["anchor", "target"], columns="animal_id",
                                values=col, aggfunc="max")
        npiv = npiv[[c for c in order if c in npiv.columns]]
        print(f"  {what}")
        print(npiv.fillna(0).astype(int).to_string())
        print()
    for col, word in [("n_anchors", "ANCHORS"), ("n_targets", "TARGETS")]:
        thin = prox.loc[prox[col] < ANCHOR_CONTRAST_MIN_N,
                        ["animal_id", "condition", "anchor", "target", col]]
        if len(thin):
            print(f"    BELOW {ANCHOR_CONTRAST_MIN_N} {word}, flagged not filtered:")
            for _, r in thin.iterrows():
                print(f"      {r['animal_id']} ({r['condition']}) "
                      f"{r['anchor'][:22]:<24}-> {r['target']:<16}"
                      f"{int(r[col]):>6}")
        else:
            print(f"    Every animal clears {ANCHOR_CONTRAST_MIN_N} {word.lower()} "
                  f"on every pair.")
    print()

    for R in PROXIMITY_RADII_UM:
        raw_col = f"pct_within_{int(R)}um"
        dlt_col = f"delta_pct_within_{int(R)}um"
        sub(f"Percent of anchors with a target within {int(R)} um")
        cols = [short_label(s) for s in SAMPLE_ORDER]
        piv = prox.pivot_table(index=["anchor", "target"], columns="animal_id",
                               values=raw_col)
        piv = piv[[c for c in cols if c in piv.columns]]
        print("  OBSERVED")
        print(piv.round(1).to_string())
        seps = _sep_list(raw_col)
        if seps:
            print(f"    complete separation, observed: " + "; ".join(seps))
        if RUN_PROXIMITY_NULL and prox[dlt_col].notna().any():
            piv2 = prox.pivot_table(index=["anchor", "target"],
                                    columns="animal_id", values=dlt_col)
            piv2 = piv2[[c for c in cols if c in piv2.columns]]
            print("\n  OBSERVED MINUS PERMUTATION NULL (the proximity analogue")
            print("  of delta: controls for local arrangement and target density)")
            print(piv2.round(1).to_string())
            seps2 = _sep_list(dlt_col)
            if seps2:
                print(f"    complete separation, null-corrected: " + "; ".join(seps2))
            else:
                print("    complete separation, null-corrected: none")

    # ---- F51 ECDFs ----------------------------------------------------------
    pairs = sorted(set(zip(percell["anchor"], percell["target"])))
    ncol = 5
    nrow = int(np.ceil(len(pairs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(9 * ncol, 8 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for ax, (a, t) in zip(axes, pairs):
        for s in SAMPLE_ORDER:
            g = percell.loc[(percell["anchor"] == a) & (percell["target"] == t)
                            & (percell["sample_id"] == s), "distance_um"]
            v = np.sort(g.dropna().to_numpy(float))
            if not len(v):
                continue
            ax.step(v, np.arange(1, len(v) + 1) / len(v), where="post",
                    linewidth=4, color=CONDITION_COLORS[COND_OF[s]], alpha=0.85)
        for R in PROXIMITY_RADII_UM[:2]:
            ax.axvline(R, color="#999999", linestyle=":", linewidth=2)
        ax.set_xscale("log")
        ax.set_xlabel("nearest-target distance (um)",
                      fontsize=FONT_SIZE_BASE - 14)
        ax.set_ylabel("cumulative fraction of anchors",
                      fontsize=FONT_SIZE_BASE - 16)
        ax.set_title(f"{a[:20]}\nto {t}", fontsize=FONT_SIZE_BASE - 12)
        ax.tick_params(labelsize=FONT_SIZE_TICK - 14)
        style_axes(ax)
    for ax in axes[len(pairs):]:
        ax.axis("off")
    fig.legend(handles=[Line2D([0], [0], color=CONDITION_COLORS[c], lw=6, label=c)
                        for c in CONDITION_ORDER],
               loc="lower center", ncol=2, frameon=False,
               fontsize=FONT_SIZE_LEGEND - 6, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("The distance data itself, one curve per animal\n"
                 "A curve shifted LEFT means anchors sit closer to that target. "
                 "Dotted lines mark the proximity radii in table 67",
                 y=1.02, fontsize=FONT_SIZE_TITLE - 8)
    save_fig(fig, "F51_distance_ecdf")

    # ---- F52 interpretable effect ------------------------------------------
    R0 = PROXIMITY_RADII_UM[1] if len(PROXIMITY_RADII_UM) > 1 else PROXIMITY_RADII_UM[0]
    pairs_p = sorted(set(zip(prox["anchor"], prox["target"])))
    panels = [(f"pct_within_{int(R0)}um",
               f"% of anchors with a target within {int(R0)} um",
               "OBSERVED (raw threshold)")]
    if RUN_PROXIMITY_NULL and prox[f"delta_pct_within_{int(R0)}um"].notna().any():
        panels.append((f"delta_pct_within_{int(R0)}um",
                       f"observed minus null, percentage points",
                       "NULL-CORRECTED (controls arrangement and density)"))
    fig, axes = plt.subplots(1, len(panels),
                             figsize=(16 * len(panels), max(11, 1.0 * len(pairs_p))),
                             squeeze=False)
    for ax, (col, xlab, ttl) in zip(axes[0], panels):
        for i, (a, t) in enumerate(pairs_p):
            g = prox.loc[(prox["anchor"] == a) & (prox["target"] == t)]
            for _, r in g.iterrows():
                off = -0.16 if r["condition"] == "D1MT" else 0.16
                ax.scatter(r[col], i + off, s=420,
                           color=CONDITION_COLORS[r["condition"]],
                           edgecolor="#FFFFFF", linewidth=2, zorder=3)
            for c, off in [("D1MT", -0.16), (REFERENCE_ARM, 0.16)]:
                v = g.loc[g["condition"] == c, col].dropna()
                if len(v):
                    ax.hlines(i + off, v.min(), v.max(),
                              color=CONDITION_COLORS[c], linewidth=4, alpha=0.5,
                              zorder=2)
        if "delta" in col:
            ax.axvline(0, color="#000000", linewidth=3)
        ax.set_yticks(range(len(pairs_p)))
        ax.set_yticklabels([f"{a[:20]} -> {t}" for a, t in pairs_p],
                           fontsize=FONT_SIZE_TICK - 12)
        ax.invert_yaxis()
        ax.set_xlabel(xlab, fontsize=FONT_SIZE_BASE - 8)
        ax.set_title(ttl, fontsize=FONT_SIZE_TITLE - 14)
        ax.xaxis.grid(True, color=GRID_COLOR, linewidth=1.5)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    axes[0][0].legend(handles=[Patch(facecolor=CONDITION_COLORS[c], label=c)
                               for c in CONDITION_ORDER],
                      frameon=False, fontsize=FONT_SIZE_LEGEND - 8, loc="best")
    fig.suptitle("The effect in units anyone can check against the images\n"
                 "one point per animal, cell-weighted across that animal's foci",
                 y=1.02, fontsize=FONT_SIZE_TITLE - 10)
    fig.tight_layout()
    save_fig(fig, "F52_proximity_fraction")


# %% Cell 7 - APPENDIX: polarisation and co-expression
# =============================================================================

banner("APPENDIX - POLARISATION AND CO-EXPRESSION (NOT A FINDING)")

print("    DEMOTED TO AN APPENDIX IN REVISION 3. On record, excluded from every")
print("    significance count, and not to be quoted.")
print()
print("    Two reasons. Segmentation spillover inflates apparent co-expression")
print("    of ANY two markers and worsens as cells pack more densely, and")
print("    untreated foci run denser. Script 08's control made that concrete:")
print("    CD3e with CD20, a pair that cannot co-occur inside one macrophage,")
print("    separated the arms at 4 of 6 thresholds, MORE than iNOS with")
print("    Arginase-1 did at 3 of 6.")
print("    And the treated arm now holds four structures in total, so the")
print("    spatial mixing statistic has almost nothing to stand on there.\n")

pol = pd.DataFrame(); mix = pd.DataFrame(); sepdf = pd.DataFrame()
if RUN_POLARISATION_APPENDIX:
    pol_rows, mix_rows = [], []
    for s in SAMPLE_ORDER:
        d = cells[s]
        mac = d.loc[d["pheno"].isin(MACROPHAGE_ALL)]
        if not len(mac) or "iNOS" not in d.columns:
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
                "sample_id": s, "animal_id": short_label(s),
                "condition": COND_OF[s], "percentile": pct,
                "n_macrophages": len(mac),
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
            # REVISION 3.3: this appendix keeps the base `rng`. With the
            # proximity, fallback-null and subsample draws moved onto their own
            # streams, nothing consumes `rng` before this point, so the mixing
            # null is now reproducible independently of every parameter above.
            for it2 in range(N_PERMUTATIONS):
                perm = rng.permutation(n)
                fi = np.zeros(n, bool); fi[perm[:n_i]] = True
                fa = np.zeros(n, bool); fa[perm[n_i:n_i + n_a]] = True
                null[it2] = mixing(fi, fa)
            pval, z, mu = empirical_p(obs, null)
            mix_rows.append({
                "sample_id": s, "condition": COND_OF[s], "structure_id": int(k),
                "percentile": POLARISATION_PRIMARY, "n_macrophages": n,
                "n_inos_only": n_i, "n_arg1_only": n_a,
                "observed_mixing": obs, "null_mixing": mu,
                "delta_mixing": obs - mu if np.isfinite(mu) else np.nan,
                "z": z, "p_empirical": pval, "tier": tier_of(min(n_i, n_a)),
            })
            del idxs

    pol = pd.DataFrame(pol_rows)
    mix = pd.DataFrame(mix_rows)

if len(pol):
    write_csv(pol, "63_polarisation_sweep.csv")
    sub(f"Co-expression ratio at the {POLARISATION_PRIMARY}th percentile")
    p0 = pol.loc[pol["percentile"] == POLARISATION_PRIMARY]
    print(p0[["animal_id", "condition", "n_macrophages", "pct_double_positive",
              "pct_expected_if_independent", "coexpression_ratio"]]
          .to_string(index=False))
    sep_rows = []
    for pct in POLARISATION_PERCENTILES:
        g = pol.loc[pol["percentile"] == pct]
        a = g.loc[g["condition"] == "D1MT", "coexpression_ratio"].dropna()
        b = g.loc[g["condition"] == REFERENCE_ARM, "coexpression_ratio"].dropna()
        if not len(a) or not len(b):
            continue
        overlap = not (a.max() < b.min() or b.max() < a.min())
        sep_rows.append({"percentile": pct, "d1mt_min": a.min(),
                         "d1mt_max": a.max(), "untr_min": b.min(),
                         "untr_max": b.max(), "overlap": overlap})
    sepdf = pd.DataFrame(sep_rows)
    if len(sepdf):
        write_csv(sepdf, "64_coexpression_separation.csv")
        clean = sepdf.loc[~sepdf["overlap"], "percentile"].tolist()
        print(f"\n    arms do not overlap at percentiles: "
              f"{clean if clean else 'none'}")
        if clean:
            print("    STILL NOT A FINDING. See the spillover caveat above.")

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
        print(f"\n    no mixing rows for: {missing}")
        print("    Too few exclusively polarised macrophages inside structures.")
        print("    With four treated structures this analysis cannot be revived.")


# %% Cell 8 - consistency and wrap up
# =============================================================================

banner("CONSISTENCY ACROSS ANIMALS")

print("    With three animals per arm, COMPLETE SEPARATION is the strongest")
print("    evidence this design can produce. It corresponds to the minimum p an")
print("    exact test can return, 0.10 two-sided, and it arises by chance with")
print("    probability 0.10, so one separated comparison among many is not")
print("    remarkable on its own. Several, running the same direction, is.\n")

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


print(f"    PER-ANIMAL SUMMARY: {PER_ANIMAL_SUMMARY}. Revision 3 used an")
print("    unweighted median across structures here, which let a 28-anchor")
print("    focus count the same as a 3,984-anchor one and disagreed with both")
print("    script 06 table 54 and this script's own exact test. Every summary")
print("    below now matches the statistic the exact test permutes.\n")

if len(percell):
    for col, label in [("distance_delta_um", "nn_delta_CONFIRMATORY"),
                       ("distance_norm_equiv", "nn_normalised_equiv_diag"),
                       ("distance_norm_inscribed", "nn_normalised_inscribed_diag"),
                       ("distance_um", "nn_raw_diag")]:
        if col not in percell.columns or not np.isfinite(percell[col]).any():
            continue
        per_animal = (percell.groupby(["sample_id", "condition",
                                       "anchor", "target"])[col]
                      .mean().reset_index())
        consistency(per_animal, ["anchor", "target"], col, label)
if len(prox):
    for R in PROXIMITY_RADII_UM:
        consistency(prox, ["anchor", "target"], f"pct_within_{int(R)}um",
                    f"proximity_{int(R)}um_observed")
        dcol = f"delta_pct_within_{int(R)}um"
        if dcol in prox.columns and prox[dcol].notna().any():
            consistency(prox, ["anchor", "target"], dcol,
                        f"proximity_{int(R)}um_NULL_CORRECTED")
if len(pol):
    p0 = pol.loc[pol["percentile"] == POLARISATION_PRIMARY].copy()
    p0["phenotype"] = "iNOS_Arg1_coexpression_APPENDIX_NOT_A_FINDING"
    consistency(p0, ["phenotype"], "coexpression_ratio", "polarisation_appendix")

# ---- REVISION 3.2: the IDO1+ minus IDO1- anchor contrast --------------------
# Both anchors are in the myeloid pool that defines a focus, so both are pulled
# core-ward by detection, and both sit in the SAME foci, on the same slide, in
# the same animal. IDO1 is used nowhere in detection. The within-animal
# difference between the two anchors therefore removes focus size, focus
# density, animal baseline and the circularity at once, because whatever
# detection does to a pool member it does to both of them equally.
def arm_margin(tvals, uvals):
    """Completely-separated gap between the arms, or None if they overlap."""
    tv = pd.Series(tvals).dropna()
    uv = pd.Series(uvals).dropna()
    if not (len(tv) and len(uv)):
        return None
    if tv.min() > uv.max():
        return float(tv.min() - uv.max())
    if uv.min() > tv.max():
        return float(uv.min() - tv.max())
    return None


contrast_long = pd.DataFrame()
if RUN_ANCHOR_CONTRAST:
    sub("ANCHOR CONTRAST: IDO1-positive minus IDO1-negative, within animal")
    print("    Both anchors are in the detection pool, in the same foci, on the")
    print("    same slide, in the same animal, and IDO1 enters detection")
    print("    nowhere. It asks whether an arm effect is SPECIFIC to the")
    print("    IDO1-positive compartment or generic to myeloid cells.")
    print()
    print("    WHAT THIS CONTRAST DOES NOT DO (revision 3.4 correction).")
    print("    Revisions 3.2 and 3.3 printed here that differencing the two")
    print("    anchors removes focus size, focus density, animal baseline and")
    print("    the circularity in one construction. That holds only if both")
    print("    anchors respond to those the same way, and on this input they do")
    print("    not: the IDO1-negative term carries its own large spread that")
    print("    IDO1-positive does not share, with the arms overlapping. It is a")
    print("    second noisy measurement, not a stable reference point, so the")
    print("    contrast is a SECOND VIEW and not a confound-cancelling")
    print("    estimator. The 3.3 header retracted that claim and this text")
    print("    said otherwise; a script may not assert what it has retracted.")
    print("    The IDO1+ alone margin is printed beside every contrast margin")
    print("    so the contrast has to earn its place on every row.")
    print()
    print("    SIGNS DIFFER BETWEEN THE TWO OUTCOMES. For proximity a POSITIVE")
    print("    contrast means the IDO1-positive anchor sits nearer its targets")
    print("    than the IDO1-negative anchor does, each against its own null.")
    print("    For distance the sign REVERSES: NEGATIVE means IDO1-positive is")
    print("    nearer. Both are printed with the direction spelled out.")
    print()
    print("    THIS IS A DIAGNOSTIC, NOT A CONFIRMATORY FAMILY. It gets no")
    print("    q-values. n_min is the smaller of the two ANCHOR counts and")
    print(f"    n_tgt the smaller TARGET count; a separation whose treated")
    print(f"    target count is under {MIN_TARGETS_FOR_REPORT} is labelled NOT")
    print("    REPORTABLE and may not reach a figure whatever its margin.")

    c_rows = []
    n_of, t_of = {}, {}

    # --- proximity, null-corrected -------------------------------------------
    if len(prox):
        n_of = {(r["sample_id"], r["anchor"], r["target"]): r["n_anchors"]
                for _, r in prox.iterrows()}
        # target counts: the binding constraint on rare populations
        t_of = {}
        for _, r in prox.iterrows():
            key = (r["sample_id"], r["target"])
            t_of[key] = max(t_of.get(key, 0), int(r.get("n_targets", 0)))
        for R in PROXIMITY_RADII_UM:
            dcol = f"delta_pct_within_{int(R)}um"
            if dcol not in prox.columns or not prox[dcol].notna().any():
                continue
            for t in NN_TARGETS:
                for s in SAMPLE_ORDER:
                    gp = prox.loc[(prox["sample_id"] == s)
                                  & (prox["target"] == t)]
                    vp = gp.loc[gp["anchor"] == IDO1_POS, dcol].dropna()
                    vn = gp.loc[gp["anchor"] == IDO1_NEG, dcol].dropna()
                    if not (len(vp) and len(vn)):
                        continue
                    c_rows.append({
                        "outcome": f"proximity_{int(R)}um_null_corrected",
                        "units": "percentage points",
                        "positive_means": "IDO1+ nearer",
                        "anchor": "IDO1+ minus IDO1-", "target": t,
                        "sample_id": s, "animal_id": short_label(s),
                        "condition": COND_OF[s],
                        "ido1_pos": float(vp.iloc[0]),
                        "ido1_neg": float(vn.iloc[0]),
                        "value": float(vp.iloc[0]) - float(vn.iloc[0]),
                        "n_ido1_pos": int(n_of.get((s, IDO1_POS, t), 0)),
                        "n_ido1_neg": int(n_of.get((s, IDO1_NEG, t), 0)),
                        "n_targets": int(t_of.get((s, t), 0)),
                    })

    # --- distance delta, the confirmatory outcome ----------------------------
    if len(percell) and "distance_delta_um" in percell.columns:
        pa = (percell.groupby(["sample_id", "condition", "anchor", "target"])
              ["distance_delta_um"].mean().reset_index())
        na = (percell.groupby(["sample_id", "anchor", "target"])
              ["distance_delta_um"].size().to_dict())
        for t in NN_TARGETS:
            for s in SAMPLE_ORDER:
                gp = pa.loc[(pa["sample_id"] == s) & (pa["target"] == t)]
                vp = gp.loc[gp["anchor"] == IDO1_POS, "distance_delta_um"].dropna()
                vn = gp.loc[gp["anchor"] == IDO1_NEG, "distance_delta_um"].dropna()
                if not (len(vp) and len(vn)):
                    continue
                c_rows.append({
                    "outcome": "distance_delta_um",
                    "units": "um", "positive_means": "IDO1+ FARTHER",
                    "anchor": "IDO1+ minus IDO1-", "target": t,
                    "sample_id": s, "animal_id": short_label(s),
                    "condition": COND_OF[s],
                    "ido1_pos": float(vp.iloc[0]), "ido1_neg": float(vn.iloc[0]),
                    "value": float(vp.iloc[0]) - float(vn.iloc[0]),
                    "n_ido1_pos": int(na.get((s, IDO1_POS, t), 0)),
                    "n_ido1_neg": int(na.get((s, IDO1_NEG, t), 0)),
                    "n_targets": int(t_of.get((s, t), 0)),
                })

    contrast_long = pd.DataFrame(c_rows)

if len(contrast_long):
    contrast_long["n_min"] = contrast_long[["n_ido1_pos",
                                            "n_ido1_neg"]].min(axis=1)
    write_csv(contrast_long, "68_anchor_contrast.csv")

    cols = [short_label(s) for s in SAMPLE_ORDER]
    sep_seen, sep_reportable = {}, {}
    for outcome, go in contrast_long.groupby("outcome", sort=False):
        direction = go["positive_means"].iloc[0]
        unit = go["units"].iloc[0]
        print(f"\n    {outcome}   ({unit}; positive = {direction})")
        piv = go.pivot_table(index="target", columns="animal_id",
                             values="value")
        piv = piv[[c for c in cols if c in piv.columns]]
        print(piv.round(1).to_string())
        nmin = go.pivot_table(index="target", columns="animal_id",
                              values="n_min", aggfunc="min")
        nmin = nmin[[c for c in cols if c in nmin.columns]]
        print("    n_min:")
        print(nmin.fillna(0).astype(int).to_string())
        hdr = False
        for t, gt in go.groupby("target"):
            td = gt.loc[gt["condition"] == "D1MT"]
            ud = gt.loc[gt["condition"] == REFERENCE_ARM]
            if not (len(td) == 3 and len(ud) == 3):
                continue
            gap = arm_margin(td["value"], ud["value"])
            if gap is None:
                continue
            # REVISION 3.3: the 3.2 flag fired on the smallest treated anchor
            # count, which says nothing about whether that animal is holding
            # the separation up. Leave-one-out answers the actual question.
            loo = []
            for drop in td["animal_id"]:
                keep = td.loc[td["animal_id"] != drop, "value"]
                g2 = arm_margin(keep, ud["value"])
                loo.append((drop, g2))
            survives = all(g2 is not None for _, g2 in loo)
            worst_loo = min((g2 for _, g2 in loo if g2 is not None),
                            default=float("nan"))
            broken = [d for d, g2 in loo if g2 is None]
            # REVISION 3.3: does the contrast earn its place, or does the
            # IDO1-positive term separate the arms just as well on its own?
            solo = arm_margin(td["ido1_pos"], ud["ido1_pos"])
            solo_txt = "  none" if solo is None else f"{solo:>6.1f}"
            verdict = ("LOO ok" if survives
                       else "LOO FAILS dropping " + ",".join(broken))
            if solo is None:
                earns = "contrast only"
            elif gap > solo * 1.5:
                earns = "contrast wider"
            elif solo > gap * 1.5:
                earns = "IDO1+ alone wider"
            else:
                earns = "no gain over IDO1+ alone"
            n_tgt_min = (int(td["n_targets"].min())
                         if "n_targets" in td else 0)
            # REVISION 3.4: a count is a label on the result, never a filter on
            # the data. Everything above is still computed and written out.
            reportable = n_tgt_min >= MIN_TARGETS_FOR_REPORT
            if not reportable:
                verdict = "NOT REPORTABLE"
                earns = (f"only {n_tgt_min} targets in the thinnest treated "
                         f"animal")
            if not hdr:
                print(f"      {'target':<16}{'gap':>7}{'worst LOO':>11}"
                      f"{'IDO1+ alone':>13}{'n_anch':>8}{'n_tgt':>7}  verdict")
                hdr = True
            print(f"      {t:<16}{gap:>7.1f}{worst_loo:>11.1f}{solo_txt:>13}"
                  f"{int(td['n_ido1_pos'].min()):>8}{n_tgt_min:>7}"
                  f"  {verdict}; {earns}")
            sep_seen.setdefault(t, set()).add(
                "distance" if outcome == "distance_delta_um" else "proximity")
            if reportable:
                sep_reportable.setdefault(t, set()).add(
                    "distance" if outcome == "distance_delta_um"
                    else "proximity")

    # register into table 66 so the separation list is complete
    for outcome, go in contrast_long.groupby("outcome", sort=False):
        consistency(go, ["anchor", "target"], "value",
                    f"anchor_contrast_{outcome}_DIAGNOSTIC")

    # ---- REVISION 3.4: do the two outcomes agree? --------------------------
    sub("Does the contrast agree across outcomes?")
    print("    The proximity outcome reads the BODY of the distance")
    print("    distribution at a fixed radius. The distance delta outcome reads")
    print("    a MEAN, which a long tail pulls. They are not the same question")
    print("    and on the revision 3.3 run they disagreed for plasma cells. A")
    print("    contrast that separates on one and not the other is not a")
    print("    finding yet, whatever either margin says.\n")
    print(f"      {'target':<16}{'proximity':>12}{'distance':>11}  status")
    for t in NN_TARGETS:
        got = sep_reportable.get(t, set())
        seen = sep_seen.get(t, set())
        pr = ("yes" if "proximity" in got
              else "(thin)" if "proximity" in seen else "no")
        di = ("yes" if "distance" in got
              else "(thin)" if "distance" in seen else "no")
        if pr == "yes" and di == "yes":
            status = "AGREE, both outcomes"
        elif pr == "yes" or di == "yes":
            status = "DISAGREE, unresolved, not for a figure"
        elif seen:
            status = "separates only where counts are too thin to report"
        else:
            status = "no separation"
        print(f"      {t:<16}{pr:>12}{di:>11}  {status}")
    print()
    print("    A target reaching AGREE is the only one the contrast supports")
    print("    on its own. DISAGREE stays a diagnostic until the difference")
    print("    between the two statistics is resolved, which script 07 cannot")
    print("    do: it needs the distribution shape, not another summary.")

    print()
    print("    HOW TO READ A SEPARATION HERE. With three animals per arm")
    print("    complete separation arises by chance with probability 0.10 per")
    print("    comparison, and there are several comparisons, so one is not a")
    print("    result. What would be a result is the untreated arm sitting near")
    print("    zero on every radius while the treated arm does not, because a")
    print("    control that behaves like a control is not a coin flip.")

cons = pd.DataFrame(cons_rows)
if len(cons):
    write_csv(cons, "66_consistency_summary.csv")
    sub("Comparisons with COMPLETE SEPARATION between arms")
    # REVISION 3.4: the same reportability label as the contrast block, so a
    # thin-target separation cannot be lifted out of this list into a figure.
    min_tgt = {}
    if len(prox) and "n_targets" in prox.columns:
        for _, r in prox.loc[prox["condition"] == "D1MT"].iterrows():
            key = r["target"]
            v = int(r["n_targets"])
            min_tgt[key] = min(min_tgt.get(key, v), v)
    if min_tgt:
        thin_t = sorted(t for t, v in min_tgt.items()
                        if v < MIN_TARGETS_FOR_REPORT)
        print(f"    Targets under {MIN_TARGETS_FOR_REPORT} cells in the thinnest")
        print(f"    treated animal, marked [THIN] below: "
              f"{', '.join(f'{t} ({min_tgt[t]})' for t in thin_t) or 'none'}\n")
    sep = cons.loc[cons["complete_separation"]]
    if not len(sep):
        print("    none")
    else:
        n_thin = 0
        for _, r in sep.sort_values("analysis").iterrows():
            name = " / ".join(str(r[c]) for c in ["phenotype", "anchor", "target"]
                              if c in r.index and pd.notna(r.get(c)))
            tgt = r.get("target")
            thin = (pd.notna(tgt) and tgt in min_tgt
                    and min_tgt[tgt] < MIN_TARGETS_FOR_REPORT)
            n_thin += int(bool(thin))
            print(f"    [{r['analysis']:<30}] {name:<46} "
                  f"D1MT [{r['D1MT_min']:.2f}, {r['D1MT_max']:.2f}]   "
                  f"Untr [{r['Untreated_min']:.2f}, {r['Untreated_max']:.2f}]"
                  f"{'   [THIN, NOT REPORTABLE]' if thin else ''}")
        n_conf = int(sep["analysis"].str.contains("CONFIRMATORY").sum())
        print(f"\n    {n_conf} of these are in the confirmatory delta outcome.")
        print(f"    {n_thin} of {len(sep)} are marked THIN and may not reach a")
        print(f"    figure: their treated target count is under "
              f"{MIN_TARGETS_FOR_REPORT}, so the permutation null they are")
        print("    measured against is built on too few cells to be read.")

banner("SUMMARY")
print(f"Input                    : {IN_DIR}")
print(f"Structure pair-rows      : {len(nn_struct)}")
print(f"Nulls shared from 06     : {n_shared}   recomputed here: {n_recomputed}")
print(f"Demotion of normalisations supported on this input : {DEMOTION_SUPPORTED}")
print(f"Models fitted            : {len(models)}")
print(f"Models dropped (no SE)   : {n_dropped_no_se}")
if len(models) and "fit_seconds" in models.columns:
    tot = float(models["fit_seconds"].sum())
    print(f"Model fitting time       : {tot/60:.1f} min "
          f"(exact LMM refit {'ON' if RUN_EXACT_LMM_REFIT else 'OFF'})")
if len(models):
    conf = models.loc[models["is_confirmatory"]]
    print(f"  confirmatory tests     : {len(conf)}")
    print(f"  p_model < 0.05         : {int(conf['sig_p05'].sum())}")
    print(f"  q < {BH_ALPHA} (BH)          : {int(conf['sig_q10'].sum())}")
    if "p_exact_means" in conf.columns:
        print(f"  exact p at the 0.10 floor : "
              f"{int((conf['p_exact_means'] <= 0.101).sum())}")
print(f"Proximity rows           : {len(prox)}   "
      f"null-corrected: {RUN_PROXIMITY_NULL} "
      f"({N_PROXIMITY_PERMUTATIONS} shuffles)")
print(f"Per-animal summary       : {PER_ANIMAL_SUMMARY}")
if len(prox) and "n_anchors" in prox.columns:
    tmin = prox.loc[prox["condition"] == "D1MT", "n_anchors"]
    print(f"Smallest treated anchor set : {int(tmin.min()) if len(tmin) else 0}"
          f"   (flag threshold {ANCHOR_CONTRAST_MIN_N})")
    if "n_targets" in prox.columns:
        gmin = prox.loc[prox["condition"] == "D1MT", "n_targets"]
        print(f"Smallest treated target set : "
              f"{int(gmin.min()) if len(gmin) else 0}"
              f"   (flag threshold {ANCHOR_CONTRAST_MIN_N})")
print("Random streams           : separate for proximity / null / subsample; "
      "subsamples keyed on identity")
print(f"Anchor contrast (table 68): {len(contrast_long)} rows, "
      f"{'ON' if RUN_ANCHOR_CONTRAST else 'OFF'}   DIAGNOSTIC, no q-values")
print(f"Polarisation (appendix)  : {len(pol)} rows, mixing {len(mix)} rows")

sub("Read in this order")
print("  1. Cell 4 verdict : does the demotion of the normalisations hold on")
print("     this input? Everything downstream assumes it does.")
print("  2. The cell-count tables, anchors AND targets. A percentage on 39")
print("     cells and one on 21,000 are not the same measurement, and a null")
print("     built on six targets is not the same null as one built on six")
print("     thousand. Read both before any percentage or separation claim.")
print("  3. F51 / table 67 : the distance data itself, and the effect in")
print("     percent-within-a-radius. Read the OBSERVED MINUS NULL row, not")
print("     the observed one. Look at both BEFORE any p-value.")
print("  4. Table 68       : the IDO1+ minus IDO1- anchor contrast, read")
print("     column by column. 'worst LOO' is the margin after dropping the")
print("     least helpful treated animal; if it is blank the separation")
print("     rests on one animal. 'IDO1+ alone' is what the IDO1-positive")
print("     term does by itself; where the contrast does not beat it, the")
print("     contrast is adding variance rather than removing confounding.")
print("     Revision 3.3 header item C is required reading before quoting")
print("     anything from this table.")
print("  5. F48 / table 62 : the confirmatory delta models, with the model p")
print("     and p_exact_means side by side. p_exact_means leads.")
print("  6. Table 66       : complete separation between arms. Remember this")
print("     and the count at the 0.10 floor are one fact, not two.")
print("  7. Appendix Cell 7: polarisation. On record, not a finding.")

sub("How to describe the statistics")
print("  The mixed model estimates the effect using every cell, with animal and")
print("  structure random intercepts. The p-value comes from enumerating all 20")
print("  assignments of the six animals to arms, and is floored at 0.10")
print("  two-sided by that design. Report the effect size, the per-animal")
print("  values, and the exact p. A model p below 0.10 should be reported as")
print("  what it is: the asymptotic approximation, not additional evidence.")
print()
print("  ON THE PROXIMITY NUMBERS SPECIFICALLY. Quote the null-corrected")
print("  value or quote observed and null together, never observed alone. On")
print("  the 3.1 run the treated arm sat within a point or two of its own")
print("  null at 50 and 100 um while the untreated arm ran far below its")
print("  own. The defensible sentence is that the IDO1-positive compartment")
print("  is spatially SEGREGATED from lymphocytes in untreated lesions and")
print("  is not segregated in treated foci, which is the loss of an")
print("  exclusion. It is not that treated macrophages were drawn together.")

sub("Superseded")
print("  Every number this script wrote against structures_rev4.")
print("  Raw radial-position models, by script 08's centred outcome.")
print("  The iNOS / Arginase-1 co-expression separation, by its own spillover")
print("  control.")
print("  Any observed proximity percentage quoted without its permutation")
print("  null, including the 99 / 98 / 95 against 55 / 44 / 20 figure at")
print("  50 um that appeared in the session notes.")
print("  Anything resting on Tregs, and anything resting on Helper T cells")
print("  without CD4- T cells beside it. Treated Treg target counts are 6,")
print("  31 and 68; Helper T carries the CD4-call problem on top.")
print("  The revision 3.2 and 3.3 runtime claim that the anchor contrast")
print("  removes size, density, animal baseline and circularity in one")
print("  construction. It is a second view, not a confound-cancelling")
print("  estimator.")

banner("DONE")
sys.stdout = _tee.terminal
_tee.close()
