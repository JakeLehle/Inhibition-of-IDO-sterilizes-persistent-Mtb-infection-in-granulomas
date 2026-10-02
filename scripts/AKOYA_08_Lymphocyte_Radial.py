#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AKOYA Phenocycler - LYMPHOCYTE RADIAL POSITION, THE PRIMARY TEST
Rhesus Mtb + SIV, D1MT-treated (G3) vs untreated (G4), necropsy lung sections

Script 08 of the AKOYA analysis series. REVISION 6.

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
       quantities, and CD4- T cells showed it on the REV4 run: per-animal
       pooled deltas of -0.226, -0.121, -0.096 against -0.052, +0.032, +0.166,
       a clean split, but a raw-radial model coefficient of -0.034 with
       p = 0.808. [REVISION 3: those are rev4 numbers and revision 2 called
       them "the current run". They are kept as the worked example that
       motivated centring, and they are NOT this script's output. Every figure
       in this section is superseded by the first rev5 run.]

       delta = mean radial of population P, minus the mean of a random subset of
       the same size. A random subset's expected mean is the structure's overall
       mean, so delta is exactly (radial_pos - structure mean radial). Centering
       each cell on its structure mean therefore makes the model test the same
       quantity the effect sizes report.

    2. CORE CELLS ONLY IN THE PRIMARY MODEL. Core cells are normalised 0 to 1 by
       each focus's own inscribed radius, but cuff cells run 1 to 2 across a
       FIXED 150 um band that is not size-normalised. Mixing them reintroduces
       structure size. Core-only is primary; core-plus-cuff is secondary.

    3. NORMALISED DISTANCE FAMILY DROPPED. [REVISION 3: the figures that were
       here came from script 07 revision 2, against structures_rev4, and were
       stale. Measured on rev5 by script 07 revision 3.5: raw median |r| = 0.07,
       delta 0.07, normalised by equivalent radius 0.36, normalised by inscribed
       radius 0.32, with 10 of 10 negative for both normalisations against 4 of
       10 and 3 of 10 for raw and delta. Same conclusion, different numbers.
       This script no longer uses ANY distance outcome, so the point is now
       historical; see Cell 8.]

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

WHAT CHANGED IN REVISION 3 (aligning with scripts 04 rev5, 06 rev4, 07 rev3.5)

    A. THE DELTA DISTANCE CELL IS REMOVED, and this is the substantive change.
       It computed observed distances over CORE CELLS ONLY and subtracted a
       null that script 06 builds over EVERY cell in the structure, core and
       cuff. Two different populations, differenced. The cuff is a fixed 150 um
       band while the core is size-normalised, so the cuff-to-core ratio tracks
       focus size, which puts the size confound back into the one construction
       built to remove it. Script 07 already owns this family, on matching
       populations, so repairing rather than removing would have produced a
       second set of delta numbers for the same pairs in a second script. That
       is how script 05 came to be frozen. Full reasoning at Cell 8.

       It also removes an inference problem. Revision 2 gave the distance
       family its own Benjamini-Hochberg correction alongside the radial
       families, treating the two as independent. Script 07 revision 3
       established they are not: both nearest-neighbour anchors sit in the
       myeloid pool that defines a focus and therefore at the density peak, so
       anchor-to-non-pool-target distance is largely radial position in
       different units. This script now tests one spatial fact.

       scipy was required only by that cell and is now optional, so a missing
       scipy no longer exits the run.

    B. INPUT REPOINTED to structures_rev5. Every number this script wrote
       against structures_rev4 is superseded.

    C. fit_mixed HAD NO STANDARD-ERROR GUARD AT ALL. It wrote std_err, ci_low,
       ci_high and p_value unchecked, so a variance component on the boundary
       at zero produced a NaN interval and a NaN p as an ordinary table row.
       Two treated animals contribute a single structure each under rev5, so
       this is a live failure mode. Both the fallback trigger and the
       acceptance test now go through one _usable() helper, matching script 07
       revision 3.5; revision 2 triggered the fallback on a finite coefficient
       alone, so a NaN standard error never reached the animal-only fit that
       might have succeeded.

    D. FALLBACKS AND THIN ARMS ARE NOW VISIBLE. fit_mode, n_structures_D1MT and
       n_structures_ref are in every row, a fallback list is printed, and the
       effective-independent-units block names any animal contributing a single
       structure. Revision 2 reported none of this, so a model resting on one
       treated structure looked like any other row.

    E. EXACT RANDOMIZATION P-VALUES, with p_exact_means leading and p_exact_lmm
       as a sensitivity analysis, matching script 07. Both are floored at
       2/20 = 0.10 two-sided by the design. Revision 2 leaned entirely on the
       model p while its own summary told the reader that p was bounded.
       One honest limitation is documented at the functions: p_exact_means
       permutes the per-animal mean of the RAW outcome, so it does not carry
       the cell-type covariate. Composition is answered by the balanced
       centring, not by that p.

    F. GEOMETRY PROVENANCE CHECK on table 35. Structure ids changed between
       rev4 and rev5, so a stale table half-matches in silence: unmatched keys
       map to NaN and the microns conversion quietly takes a median over
       whatever did match. Same guard script 07 item 11 applies to the shared
       null, applied to the table this script actually reads. Revision 2 read
       table 53 with a column check only and no key match; that read is gone
       with the distance cell.

    G. A THIN-POPULATION REPORT, table 70b, and deliberately NO GATE YET.
       Script 07 carries MIN_TARGETS_FOR_REPORT = 50, but that threshold was
       set for a nearest-neighbour target count and should not be imported here
       before anyone has seen the radial counts. The block prints the minimum
       per arm and, separately, whether a phenotype is absent from an animal
       altogether, which silently drops that animal for that phenotype. Set the
       gate from the first rev5 run.

    H. STALE TEXT CORRECTED. The revision 2 header quoted script 07 revision 2
       correlations from rev4. The summary told the reader to convert with the
       inscribed radius because "the equivalent radius overstates by the shape
       ratio, about 1.13 treated and 1.22 untreated", which was wrong twice:
       rev4-era figures, and 1.13 is not a shape ratio, it is the rev5
       untreated-to-treated equivalent radius gap. Rev5 shape ratios are 1.34
       and 1.24.

    NOT PORTED, AND WHY. Script 07 revision 3.5 added a reproducibility
    fingerprint and a cross-run ledger because its anchor subsample had been
    drawn from a running generator. This script's per-structure subsample uses
    a fixed random_state per group, so each draw depends only on that group's
    contents and never on iteration order or on any permutation count. It never
    had that defect and the machinery would be noise here.

WHAT CHANGED IN REVISION 3.1 (after reading the first rev5 run, 21 Sep 2026)

    No change to the hypothesis or the primary outcome. One arithmetic defect,
    one family computed on the wrong reference, and three things the first run
    made visible that the script was not reporting.

    A. p_exact_lmm RETURNED A VALUE BELOW ITS OWN FLOOR. The run printed
       p_exact_lmm = 0.050 for the primary model, two lines under text stating
       both exact p-values are floored at 0.10. That value is not possible. The
       twenty assignments come in ten sign-symmetric pairs, because flipping
       every arm label turns the indicator into 1 - arm and negates the
       coefficient exactly, so |coef| is shared within a pair and the minimum
       two-sided p is 2/20.

       CAUSE. Each assignment was a separate lbfgs refit compared against
       |obs| - 1e-9. The observed assignment's own complement refits to roughly
       1e-7 from the exact negation, fell below that threshold, and was not
       counted. p_exact_means was unaffected because it is exact arithmetic,
       which is why it sat correctly at 0.100 in the same run.

       FIX. Enumerate ONE representative per sign-symmetric pair when the design
       is balanced and even, and take the observed representative's statistic as
       |obs| by definition rather than refitting it. The floor now holds by
       construction, not by a tolerance surviving optimizer noise, and it halves
       the refits. A floor assertion prints loudly and raises any p that still
       comes back under 2/n. Unbalanced designs, such as leave-one-out at 2
       versus 3, have no pairing and enumerate in full as before.

    B. THE PER-PHENOTYPE FAMILY RAN ONLY ON THE CELL-WEIGHTED CENTRING. It
       reported "5 of 5 at q < 0.1" for the lymphocytes, the most quotable
       numbers in the output, on the reference the centring check had just shown
       moves with composition: pooled lymphocytes fell from -0.1454 to -0.0622
       between the two centrings, and T lineage from -0.0976 to -0.0119. The
       family now runs on every centring, balanced first, with primary kept
       beside it and labelled as the comparison.

    C. THE BALANCED REFERENCE WAS NOT THE SAME QUANTITY IN EVERY STRUCTURE. It
       averaged over whichever phenotypes cleared MIN_CELLS_FOR_BALANCE in that
       structure, so the first run balanced 43109 over 8 phenotypes and
       everything else over 10 or 11. A mean over 8 populations is not a mean
       over 11, so the reference meant to be composition-free was not comparable
       between arms. A COMMON set is now derived across structures and used for
       a third centring, computed alongside the adaptive one so the difference
       is measured rather than assumed.

    D. A REPORTABILITY GATE, SET FROM THIS SCRIPT'S OWN COUNTS. Revision 3
       printed the counts and deliberately applied no gate. They came back with
       minimum treated core counts of Tregs 6, Helper T cells 28, B cells 57,
       Plasma cells 78, CD4- T cells 197, so a threshold of 50 separates the two
       populations already discounted on other grounds from the three carrying
       the result. A LABEL, never a filter, and F52 will not plot a retired row.
       Tregs is why this matters: 528 cells, p = 0.0396, separated on both
       centrings, survived BH, and rests on six cells in 43106.

    E. THE TREATED ARM IS MOSTLY ONE ANIMAL, AND THE SCRIPT NOW SAYS SO. 43118
       carried 3,273 of 4,264 treated core lymphocytes, 77 percent, and holds
       the 1.59 shape ratio and the 3.6 percent IDO1 outlier. The animal random
       intercept absorbs its baseline but not its weight on the arm term, so the
       cell-weighted coefficient is largely a statement about that animal.
       Printed with the arm composition, and it is the strongest argument for
       p_exact_means as the headline, since that weights animals equally.

    F. THE SPILLOVER VERDICT NOW CHECKS DIRECTION. Spillover rises with density
       and untreated foci are three to five fold denser, so it predicts higher
       apparent co-expression in UNTREATED. On the first rev5 run iNOS with
       Arginase-1 did that, at 1.52, 1.73, 1.62 untreated against 1.36, 1.12,
       1.49 treated. CD3e with CD20, the control that retired it, ran the other
       way: 1.99, 1.91, 1.65 treated against 1.68, 1.62, 1.62 untreated.
       Revision 3 counted separations without checking sign and dropped the
       result on a control pointing away from the mechanism it was standing in
       for. The verdict now separates those cases and, where no control
       separates in the spillover direction, returns UNRESOLVED rather than
       DROP. That is not a reinstatement: a pair that cannot co-occur is still
       separating the arms, and that needs explaining before either pair is
       quoted.

WHAT CHANGED IN REVISION 4 (the multiple-testing rule, agreed 24 Sep 2026)

    A. ONE CONFIRMATORY CENTRING. Revision 3.1 pooled every balanced centring
       into one Benjamini-Hochberg family, so per_phenotype_lymphocyte_balanced
       ran at n = 10: five phenotypes measured twice. That is one hypothesis on
       two references, not ten hypotheses, and adding the area-weighted
       reference would have made it fifteen.

       It decides a result. On the 23 September run plasma cells at p = 0.0544
       sits at rank 5 of 10 and fails BH, and at rank 3 of 5 and passes, and
       plasma cells is the population at the centre of the script 07 versus
       script 08 disagreement. The rule below was therefore fixed in writing
       BEFORE this run, not chosen from the results.

       One centring is confirmatory and its families carry q-values. Every
       other centring is a sensitivity analysis and carries none, the same
       treatment leave_one_out and marker_robustness already receive, for the
       same reason. The confirmatory centring is the area-weighted one once
       script 04 exports mean_radial_over_area, and the common set until then.
       See CONFIRMATORY_CENTRING_PREFERENCE in Cell 1.

    B. THE HEADLINE IS THE CONFIRMATORY CENTRING. Revision 3.2 printed the
       cell-weighted pooled model under "PRIMARY RESULT" on the same page as a
       centring check showing that reference loses more than half the
       coefficient. The summary now leads with the confirmatory fit and keeps
       the cell-weighted one beside it, labelled as a comparison.

    C. F52 PLOTS EACH PHENOTYPE ONCE. It named the balanced families
       explicitly, so a third centring would have drawn every phenotype three
       times. It now plots the confirmatory families.

    D. The per-phenotype blocks are ordered confirmatory first, and the pandas
       FutureWarning at the reportability filter is gone.

WHAT CHANGED IN REVISION 5 (after reading the revision 4 run, 27 Sep 2026)

    No change to the multiple-testing rule. Revision 4's rule was fixed in
    advance precisely so that it would not be revised once it moved a result,
    and it did move one: B cells, Helper T cells and plasma cells now survive
    BH on the confirmatory centring, up from B cells alone. That stands.

    A. THE PER-PHENOTYPE BLOCK PRINTED THE WRONG P. It showed the MODEL p and
       not p_exact_means, in the most quotable table in the output, after this
       script had already declared that the exact animal-level p leads. It
       mattered on the revision 4 run: plasma cells survived BH at
       p_model = 0.0544 while showing NO per-animal separation on the same
       centring, which means its exact p cannot be better than the 0.10 floor.
       Anyone working from p_model alone would have quoted it as a result. Both
       p-values now print side by side, with a flag where they disagree about
       significance.

    B. THE SUMMARY DID NOT SAY WHETHER THE POOLED CLAIM SURVIVES. Revision 4
       printed, under the heading HEADLINE, a coefficient of -0.0784 with a 95
       percent interval crossing zero and p_exact_means of 0.300, and left the
       reader to work out what that meant. It means the pooled claim does not
       survive on the confirmatory reference: the six animal-level means do not
       separate. "Lymphocytes sit more core-ward in treated foci" is therefore
       not supported as a pooled confirmatory statement, and the finding is the
       PER-PHENOTYPE one. The summary now says so, and says explicitly that the
       cell-weighted pooled result is not a fallback.

    C. TWO DIFFERENT COMMON SETS WERE SHARING ONE NAME. This script derives
       coverage over CORE cells, because its outcome is core-only radial
       position. Script 07b derives its own over ALL structure cells, because it
       must reproduce script 06's null, which shuffles labels over every cell in
       the structure. Same nominal 90 percent threshold, different sets: 4
       phenotypes here and 8 in script 07b on the 27 September run, because the
       cuff pushes sparse populations over the 10-cell bar. Each is right for
       its own script. Both are now named by the population they are derived
       over, and script 07b's stale claim that its set matched this one is
       corrected there.

WHAT CHANGED IN REVISION 6 (29 Sep 2026)

    A. THE REVISION 5 FLAG FIRED ON THE WRONG ROWS. It tested
       p_model < 0.05 <= p_exact. On a design floored at 0.10 that trips
       whenever the model p is small and the exact p is AT the floor, which is
       complete separation of the six animal means and the strongest result
       this design can produce. So on the 27 September run B cells, the
       cleanest result in the analysis at p_model 0.0053 and p_exact 0.100,
       carried a warning it did not deserve. And the row that needed one did
       not get it: plasma cells came back p_model 0.0544 with p_exact 0.200,
       meaning the animals do NOT separate, and 0.0544 is not below 0.05 so
       nothing printed. Exactly backwards.

       The comparison is now against the floor the design allows, which
       fit_mixed carries as p_exact_floor, not against 0.05. A row at the floor
       is marked [animals separate]; a row above it is marked ANIMALS DO NOT
       SEPARATE with both numbers shown. The model p no longer enters the test.

    B. p_exact_floor is threaded out of the exact test and into every model
       row, so the distinction above can be made in the tables as well as in
       the printed block.

    NOTE ON THE CONFIRMATORY CENTRING. Script 04 revision 6 exports
    mean_radial_over_area, so on the next run AREA_REFERENCE_COL resolves and
    CONFIRMATORY_CENTRING_PREFERENCE promotes the area-weighted reference to
    confirmatory. That is the preference order declared in revision 4 doing
    what it was written to do, not a new decision, and the q-values will move
    with it.

INPUTS
    structures_rev5/cell_assignments/<section>_cell_structures.csv
    structures_rev5/tables/35_foci_structures_relative.csv
    AKOYA/data/<section>.csv                        (spillover control markers)
    (table 53, the shared null, is NO LONGER read: see revision 3 item A)

OUTPUTS
    figures/  F51 .. F55
    tables/   70, 70b, 71 .. 78   (74, delta distance, is gone in revision 3)

USAGE
    conda activate sc_pre
    python AKOYA_08_Lymphocyte_Radial.py

Author: Jake Lehle, Kaushal Lab, Texas Biomed
"""

# %% Cell 1 - parameters
# =============================================================================

DATA_DIR = "/master/jlehle/WORKING/AKOYA/data"
IN_DIR = "/master/jlehle/WORKING/AKOYA/structures_rev5"
CELL_DIR = f"{IN_DIR}/cell_assignments"
FOCI_TABLE = f"{IN_DIR}/tables/35_foci_structures_relative.csv"
# NN_NULL_TABLE removed in revision 3: this script no longer reads the
# shared null, because it no longer computes distance. See Cell 8.
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
# REVISION 3.1. The balanced reference in revision 3 averaged over whichever
# phenotypes cleared MIN_CELLS_FOR_BALANCE IN THAT STRUCTURE, so different
# structures averaged over different phenotype sets: the first rev5 run
# balanced 43109 over 8 phenotypes and everything else over 10 or 11. A mean
# over 8 populations is not the same quantity as a mean over 11, so a reference
# meant to be composition-free was not comparable between arms. A COMMON set is
# now derived across all structures and used for a third centring, computed
# alongside the adaptive one rather than replacing it, so the difference is
# visible rather than assumed.
BALANCE_COMMON_MIN_FRAC = 0.90   # phenotype must clear the threshold in this
                                 # fraction of structures to join the common set
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
# REVISION 3: REMOVED. Script 08 no longer computes distance at all. See the
# docstring, "WHAT CHANGED IN REVISION 3", item A. Distance belongs to script
# 07 and only to script 07.

# ---- geometry provenance ----------------------------------------------------
# REVISION 3. Structure ids changed between rev4 and rev5, so a stale table 35
# would half-match in silence: geometry would come back NaN for the missing
# keys, the microns conversion would quietly use a median over whatever did
# match, and nothing would say so. Same failure mode script 07 item 11 guards
# against on the shared null, applied to the table this script actually reads.
REQUIRE_GEOMETRY_PROVENANCE = True
GEOMETRY_PROVENANCE_MIN_MATCH = 0.95

# ---- exact randomization ----------------------------------------------------
# REVISION 3. Treatment was assigned to six animals, three per arm, so there
# are C(6,3) = 20 assignments and any test that counts them is floored at
# 2/20 = 0.10 two-sided. The mixed model stays the ESTIMATOR; these supply a
# calibrated p. Matches script 07 revision 3.5, including which one leads.
RUN_EXACT_RANDOMIZATION = True
RUN_EXACT_LMM_REFIT = True
# REVISION 3.1: relative tolerance for the exact-p comparison. 1e-9 was too
# tight for a value produced by a separate lbfgs refit and cost the floor.
EXACT_P_REL_TOL = 1e-6
PERM_MAX_CELLS = 20000              # fixed subsample reused by all 20 refits
PERM_RE_MODE = "animal only"        # nested fit kept for the point estimate

# ---- models -----------------------------------------------------------------
# NOTE ON THE CAP. Script 07 uses 2000 and this uses 3000, deliberately: 07
# caps ANCHORS per structure for a nearest-neighbour search, this caps
# LYMPHOCYTES per structure for a position model. Different quantities, so the
# numbers need not match. Recorded here so the difference is not read as drift.
MODEL_MAX_CELLS_PER_STRUCTURE = 3000
MODEL_MIN_CELLS_PER_ANIMAL = 20
BH_ALPHA = 0.10
RANDOM_SEED = 0
# families that are refits of one hypothesis, not independent tests
NO_BH_FAMILIES = ["leave_one_out", "marker_robustness", "centring_check"]

# ---- REVISION 4: ONE CONFIRMATORY CENTRING --------------------------------
# Revision 3.1 put every balanced centring into one Benjamini-Hochberg family,
# so per_phenotype_lymphocyte_balanced ran at n = 10: FIVE phenotypes measured
# TWICE. That is one hypothesis on two references, not ten hypotheses, and it
# inflates m. Adding the area-weighted centring would have made it fifteen.
#
# It is not academic. On the 23 September run, plasma cells at p = 0.0544 sits
# at rank 5 of 10 and fails BH, and at rank 3 of 5 and passes. The family
# structure therefore decides whether plasma cells survives multiple testing,
# and plasma cells is the population at the centre of the script 07 versus
# script 08 disagreement.
#
# THE RULE, DECLARED IN ADVANCE OF SEEING WHICH WAY IT FALLS.
#   One centring is CONFIRMATORY and its families get q-values.
#   Every other centring is a SENSITIVITY analysis and gets none, the same
#   treatment leave_one_out and marker_robustness already receive, and for the
#   same reason: a refit of one hypothesis on a different reference is not an
#   independent test.
#   The confirmatory centring is the AREA-WEIGHTED one once script 04 exports
#   mean_radial_over_area, because it is the only reference that cannot move
#   with composition at all. Until then it is the common set.
CONFIRMATORY_CENTRING_PREFERENCE = [
    "radial_centered_core_area",             # exact, composition-free
    "radial_centered_core_balanced_common",  # one phenotype set everywhere
    "radial_centered_core_balanced",         # per-structure phenotype set
    "radial_centered_core",                  # cell-weighted, last resort
]

# ---- reportability gate (REVISION 3.1) --------------------------------------
# Revision 3 printed the counts and applied no gate, deliberately, so the
# threshold could be set from the radial numbers rather than inherited from a
# nearest-neighbour target count. The first rev5 run gave minimum treated core
# counts of: Tregs 6, Helper T cells 28, B cells 57, Plasma cells 78,
# CD4- T cells 197. A threshold of 50 therefore separates the two populations
# that were already independently discounted, Tregs and Helper T cells, from
# the three that carry the result. It matches script 07's number, but it is set
# here on this script's own evidence rather than imported.
#
# A LABEL, NEVER A FILTER. Everything is still computed, written and printed.
MIN_TREATED_CELLS_FOR_REPORT = 50

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
from itertools import combinations

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# REVISION 3: scipy was required only by the distance cell, which is gone.
# Kept as an optional import so the environment check still reports it, but it
# is no longer a hard dependency and no longer exits the run.
try:
    from scipy.spatial import cKDTree  # noqa: F401  (retained for diagnostics)
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


# =============================================================================
# EXACT RANDOMIZATION AT THE ANIMAL LEVEL  (REVISION 3, ported from script 07)
# =============================================================================
# Treatment was assigned to six animals, three per arm. There are C(6,3) = 20
# assignments, in 10 sign-symmetric pairs, so the two-sided p is FLOORED AT
# 2/20 = 0.10. That floor is a property of the design and no model escapes it.
# A model p far below it is buying resolution with an assumption, which on six
# clusters is not a safe purchase: a null simulation on this cluster structure
# measured false positive rates of 46.7 percent for GLMM variants against 4.4
# percent for a two-stage animal-means test.
#
# WHICH ONE LEADS, DECLARED, matching script 07 revision 3.2 onward.
# p_exact_means is the headline. Treatment was assigned to ANIMALS, so
# permuting the six animal-level means is the randomisation that actually
# happened and it needs no model. p_exact_lmm is a sensitivity analysis: it
# refits under each assignment, so its statistic is scaled by
# structure-to-structure variation again, and the untreated arm carries 76
# heterogeneous structures against 4 treated. That is the same thing that
# inflates the model p, reimported into the test built to avoid it. Where the
# two disagree, quote p_exact_means and report the disagreement.
#
# ONE HONEST LIMITATION, specific to this script. p_exact_means permutes the
# per-animal mean of the RAW outcome and therefore does not carry the cell-type
# covariate the primary model adjusts for. For the pooled lymphocyte model that
# means the covariate adjustment is absent from the exact p. The composition
# question is answered by the BALANCED CENTRING outcome, not by this p, so read
# the two together and do not treat p_exact_means on the primary outcome as
# composition-adjusted.


def _arm_labels(animals, idx, treated_label="D1MT"):
    return {a: (treated_label if i in idx else REFERENCE_ARM)
            for i, a in enumerate(animals)}


def _assignment_reps(n, n_t):
    """
    REVISION 3.1. Yields (idx_tuple, paired).

    THE BUG THIS FIXES. A two-sided exact randomization p on six animals, three
    per arm, cannot go below 2/20 = 0.10, because the twenty assignments come in
    ten SIGN-SYMMETRIC PAIRS: flipping every arm label turns the binary
    indicator arm into 1 - arm, so the coefficient is negated exactly and the
    two members of a pair have identical |coef|. Revision 3 enumerated all
    twenty and compared each against |obs| with a tolerance of 1e-9. For
    p_exact_means that is exact arithmetic and the floor held. For p_exact_lmm
    each assignment is a SEPARATE lbfgs refit, which lands perhaps 1e-7 from the
    exact negation, so the observed assignment's own complement fell below
    |obs| - tol and was not counted. The run of 21 September 2026 printed
    p_exact_lmm = 0.050 for the primary model, which is below the floor the
    script declares two lines above it, and is not a possible value.

    THE FIX. When the design is balanced and even, enumerate ONE representative
    per sign-symmetric pair, pinning animal 0 to the treated set, and count each
    representative once out of n_pairs. The floor is then 1/n_pairs = 2/C(n,k)
    BY CONSTRUCTION, and no tolerance comparison can break it. It also halves
    the number of refits.

    When the design is not balanced, as in leave-one-out at 2 versus 3, the
    complement of a 2-set in 5 animals is a 3-set, so there is no pairing and
    every assignment is enumerated as before.
    """
    if n % 2 == 0 and n_t * 2 == n:
        for rest in combinations(range(1, n), n_t - 1):
            yield (0,) + rest, True
    else:
        for idx in combinations(range(n), n_t):
            yield idx, False


def _canonical(idx, n, paired):
    """The representative of idx's sign-symmetric pair, as a sorted tuple."""
    if not paired:
        return tuple(sorted(idx))
    return (tuple(sorted(idx)) if 0 in idx
            else tuple(sorted(set(range(n)) - set(idx))))


def _finish_exact(stats, obs, paired, label):
    """
    Turns a list of per-representative |statistic| values into a two-sided p,
    and REFUSES to return one below the floor the design allows. A p below the
    floor is a bug, never a result, and revision 3 printed one.
    """
    stats = np.asarray([v for v in stats if np.isfinite(v)], float)
    if not len(stats):
        return np.nan, 0, np.nan
    tol = EXACT_P_REL_TOL * max(1.0, abs(obs))
    hits = int(np.sum(stats >= abs(obs) - tol))
    n_assign = 2 * len(stats) if paired else len(stats)
    p = hits / len(stats)

    # REVISION 3.2 FIX. The floor is ONE enumerated unit out of the enumerated
    # units, always, which is 1/len(stats) whether or not those units are
    # sign-symmetric pairs. Revision 3.1 wrote it as 2/n_assign, which is right
    # for the PAIRED case, where n_assign = 2 * len(stats) and the two reduce to
    # the same thing, and WRONG for the unpaired case, where it doubled the
    # floor. Leave-one-out is unpaired, at 2 treated versus 3 untreated on five
    # animals, so C(5,2) = 10 assignments with a true floor of 1/10. The 23
    # September run therefore raised three legitimate p_exact_means values of
    # 0.100 to 0.200 and printed an alarming violation notice for each. The
    # assertion was right to exist and wrong about the threshold.
    floor = 1.0 / len(stats)
    if np.isfinite(p) and p < floor - 1e-12:
        print(f"    EXACT-P FLOOR VIOLATED ({label}): p={p:.4f} < {floor:.4f} "
              f"on {n_assign} assignments.")
        print("      This is arithmetically impossible for a two-sided exact")
        print("      randomization test and indicates a defect, not a result.")
        print(f"      Raising to the floor. Investigate before quoting {label}.")
        p = floor
    return float(p), int(n_assign), float(floor)


def exact_p_animal_means(d, outcome):
    """
    Permutes the six animal-level means. Instant and exactly calibrated.
    Statistic: difference of arm means of the per-animal means. The per-animal
    value is a CELL-WEIGHTED mean, matching how every other per-animal summary
    in this series is formed.
    """
    g = d.dropna(subset=[outcome, "sample_id", "condition"])
    if not len(g):
        return np.nan, np.nan, 0
    m = g.groupby(["sample_id", "condition"])[outcome].mean().reset_index()
    animals = list(m["sample_id"])
    vals = m[outcome].to_numpy(float)
    arms = list(m["condition"])
    n = len(animals)
    n_t = sum(1 for a in arms if a != REFERENCE_ARM)
    if n_t == 0 or n_t == n:
        return np.nan, np.nan, 0
    obs_idx = tuple(i for i, a in enumerate(arms) if a != REFERENCE_ARM)

    def _stat(idx):
        t = vals[list(idx)]
        u = vals[[i for i in range(n) if i not in idx]]
        return float(np.mean(t) - np.mean(u))

    obs = _stat(obs_idx)
    stats, paired_flag = [], False
    for idx, paired in _assignment_reps(n, n_t):
        paired_flag = paired
        stats.append(abs(_stat(idx)))
    p, n_assign, floor = _finish_exact(stats, obs, paired_flag, "p_exact_means")
    return p, obs, n_assign, floor


def exact_p_lmm(d, outcome, terms):
    """
    Refits under every assignment representative using ANIMAL-ONLY random
    effects on a fixed capped subsample, because nested refits are not
    tractable at this count. The nested fit is retained for the point estimate;
    this supplies only a p.

    REVISION 3.1: the observed assignment's representative is NOT refitted. Its
    statistic is |obs|, which is what it is by definition, so it always counts
    as a hit and the floor cannot be lost to optimizer noise. That was the
    mechanism behind the impossible 0.050.

    The subsample is drawn with a fixed random_state, so it does not depend on
    stream position or on iteration order.
    """
    if not (HAVE_SM and RUN_EXACT_LMM_REFIT):
        return np.nan, 0
    g = d.dropna(subset=[outcome, "sample_id", "condition"]).copy()
    if not len(g) or g["sample_id"].nunique() < 4:
        return np.nan, 0
    if len(g) > PERM_MAX_CELLS:
        g = g.sample(PERM_MAX_CELLS, random_state=RANDOM_SEED)
    animals = sorted(g["sample_id"].unique())
    n = len(animals)
    arm_of = {a: g.loc[g["sample_id"] == a, "condition"].iloc[0] for a in animals}
    n_t = sum(1 for a in animals if arm_of[a] != REFERENCE_ARM)
    if n_t == 0 or n_t == n:
        return np.nan, 0
    formula = f"{outcome} ~ " + " + ".join(terms)

    def _coef(labels):
        g["arm"] = g["sample_id"].map(
            lambda a: 0.0 if labels[a] == REFERENCE_ARM else 1.0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                r = smf.mixedlm(formula, data=g, groups=g["sample_id"],
                                re_formula="1").fit(reml=True, method="lbfgs",
                                                    maxiter=300)
                v = float(r.params.get("arm", np.nan))
                return v if np.isfinite(v) else np.nan
            except Exception:
                return np.nan

    obs = _coef(arm_of)
    if not np.isfinite(obs):
        return np.nan, 0
    obs_idx = tuple(i for i, a in enumerate(animals) if arm_of[a] != REFERENCE_ARM)
    reps = list(_assignment_reps(n, n_t))
    paired_flag = reps[0][1] if reps else False
    obs_rep = _canonical(obs_idx, n, paired_flag)

    stats = []
    for idx, _paired in reps:
        if _canonical(idx, n, paired_flag) == obs_rep:
            stats.append(abs(obs))          # by definition, not by refit
            continue
        stats.append(abs(_coef(_arm_labels(animals, idx))))
    p, n_assign, _floor = _finish_exact(stats, obs, paired_flag, "p_exact_lmm")
    return p, n_assign


def capped(df):
    """
    Per-structure subsample for the models. REVISION 3: was defined twice,
    identically, in Cell 5 and Cell 6. One definition now.

    NOTE ON REPRODUCIBILITY. random_state is a fixed constant, so each group's
    draw depends only on that group's contents, never on iteration order or on
    how many groups came before. This script therefore never had script 07's
    crossed-stream defect and needs none of that machinery ported.
    """
    if not len(df):
        return df
    return pd.concat(
        [g.sample(min(len(g), MODEL_MAX_CELLS_PER_STRUCTURE),
                  random_state=RANDOM_SEED)
         for _, g in df.groupby(["sample_id", STRUCT_COL])],
        ignore_index=True)


def _usable(res, term="arm"):
    """
    REVISION 3. Ported from script 07 revision 3.5. One definition of a usable
    fit, applied in BOTH places that need it: the fallback trigger and the
    acceptance guard.

    A finite coefficient is not enough. statsmodels returns a finite
    coefficient with a NaN standard error, raising nothing, when a variance
    component sits on the boundary at zero. Under structures_rev5 two treated
    animals contribute a SINGLE STRUCTURE EACH, so this is not hypothetical
    here. Revision 2 had no guard at all and would have written the NaN
    standard error, a NaN confidence interval and a NaN p into table 72 as an
    ordinary row.
    """
    try:
        c = float(res.params.get(term, np.nan))
        e = float(res.bse.get(term, np.nan))
        pv = float(res.pvalues.get(term, np.nan))
    except Exception:
        return False
    return bool(np.isfinite(c) and np.isfinite(e) and e > 0 and np.isfinite(pv))


def fit_mixed(df, outcome, label, covariates=None, note=""):
    """
    outcome ~ arm [+ covariates], random intercept for animal,
    structure nested within animal. Falls back to animal-only on failure.

    REVISION 3: the fallback trigger and the acceptance test both go through
    _usable(). Revision 2 triggered the fallback on a finite coefficient alone,
    so a nested fit with a NaN standard error never reached the animal-only fit
    that might have succeeded, and then nothing caught it on the way out.
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
        if res is None or not _usable(res):
            try:
                md = smf.mixedlm(formula, data=d, groups=d["sample_id"],
                                 re_formula="1")
                res = md.fit(reml=True, method="lbfgs", maxiter=300)
                mode = "animal only"
            except Exception:
                return None
    if res is None:
        return None

    # REVISION 3 GUARD. A finite coefficient is not enough; see _usable().
    if not _usable(res):
        coef = float(res.params.get("arm", np.nan))
        se = float(res.bse.get("arm", np.nan))
        pv = float(res.pvalues.get("arm", np.nan))
        print(f"    NO USABLE SE: {label} / {outcome} ({mode}). "
              f"coef={coef} se={se} p={pv}.")
        print("      The random-effect variance is on the boundary at zero, "
              "which happens")
        print("      when an animal contributes a single structure. Both the "
              "nested fit and")
        print("      the animal-only fallback failed this test. Row DROPPED "
              "rather than")
        print("      reported as NaN.")
        return None

    means = d.groupby("condition")["_y"].mean()
    n_struct_by_arm = d.groupby("condition")["struct_key"].nunique()
    coef = float(res.params.get("arm", np.nan))
    se = float(res.bse.get("arm", np.nan))

    # REVISION 3: exact randomization p-values beside the model p.
    p_means, p_lmm = np.nan, np.nan
    n_assign_means, n_assign_lmm = 0, 0
    p_means_floor = np.nan
    if RUN_EXACT_RANDOMIZATION:
        p_means, _obs, n_assign_means, p_means_floor = exact_p_animal_means(d, "_y")
        p_lmm, n_assign_lmm = exact_p_lmm(d, "_y", terms)

    return {
        "analysis": label, "outcome": outcome,
        "covariates": ",".join(covariates), "random_effects": mode,
        # REVISION 3: fit_mode and the per-arm structure counts were absent in
        # revision 2, so a fallback to animal-only random effects, and a model
        # resting on a single treated structure, both looked like any other row.
        "fit_mode": mode,
        "n_cells": int(len(d)), "n_animals": int(d["sample_id"].nunique()),
        "n_animals_D1MT": int(d.loc[d["condition"] == "D1MT", "sample_id"].nunique()),
        "n_animals_ref": int(d.loc[d["condition"] == REFERENCE_ARM, "sample_id"].nunique()),
        "n_structures": int(d["struct_key"].nunique()),
        "n_structures_D1MT": int(n_struct_by_arm.get("D1MT", 0)),
        "n_structures_ref": int(n_struct_by_arm.get(REFERENCE_ARM, 0)),
        "mean_D1MT": float(means.get("D1MT", np.nan)),
        "mean_Untreated": float(means.get("Untreated", np.nan)),
        "coef_D1MT_vs_ref": coef, "std_err": se,
        "ci_low": coef - 1.96 * se, "ci_high": coef + 1.96 * se,
        "p_value": float(res.pvalues.get("arm", np.nan)),
        "p_exact_means": p_means, "p_exact_lmm": p_lmm,
        # REVISION 6: the floor the design allows, so a reader can tell
        # "as good as this design gets" from "the animals do not separate".
        "p_exact_floor": p_means_floor,
        "n_assignments_means": n_assign_means, "n_assignments_lmm": n_assign_lmm,
        "note": note,
    }


_tee = Tee(os.path.join(TAB_DIR, "00_lymphocyte_radial_report.txt"))
sys.stdout = _tee

banner("AKOYA LYMPHOCYTE RADIAL POSITION - PRIMARY TEST (revision 6)")
print(f"Run time : {datetime.now().isoformat(timespec='seconds')}")
print(f"Input    : {IN_DIR}")
print(f"Output   : {OUT_DIR}")
if not HAVE_SCIPY:
    print(f"\n    NOTE: scipy unavailable ({_scipy_err}). Not required in")
    print("    revision 3, which computes no distances. Continuing.")
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

# ---- REVISION 3.1: a COMMON-SET composition-balanced reference --------------
# The adaptive reference above is kept. This one restricts every structure to
# the same phenotype list, so the reference is the same quantity everywhere.
_core_all = cells.loc[cells["is_core"]]
_ph_stats = (_core_all.groupby(["sample_id", STRUCT_COL, "pheno"])["radial_pos"]
             .agg(["mean", "size"]).reset_index())
_n_struct_total = int(_core_all.groupby(["sample_id", STRUCT_COL]).ngroups)
_ok = _ph_stats.loc[_ph_stats["size"] >= MIN_CELLS_FOR_BALANCE]
_cover = (_ok.groupby("pheno").size() / max(_n_struct_total, 1)).sort_values(ascending=False)
COMMON_BALANCE_PHENOS = sorted(_cover.loc[_cover >= BALANCE_COMMON_MIN_FRAC].index)

sub("Common-set balanced reference, derived over CORE cells (revision 3.1)")
print("    Revision 3 let each structure average over whatever phenotypes it")
print("    happened to have at least "
      f"{MIN_CELLS_FOR_BALANCE} cells of, so the reference was a mean over a")
print("    different set in different structures. This derives one set and")
print("    applies it everywhere.")
print()
print("    REVISION 5: DERIVED OVER CORE CELLS, and the name now says so.")
print("    This script's outcome is core-only radial position, so coverage is")
print("    counted over core cells, which is the correct population for it.")
print("    Script 07b derives its own common set over ALL structure cells,")
print("    core plus cuff, because it must reproduce script 06's null. Same")
print("    nominal 90 percent threshold, DIFFERENT SETS: on the 27 September")
print("    run this gave 4 phenotypes here and 8 in script 07b, because the")
print("    cuff pushes sparse populations over the 10-cell bar. Both are")
print("    right for their own script and neither is 'the' common set. Never")
print("    quote one as though it were the other.\n")
print(f"    {'phenotype':<28}{'structures at threshold':>26}")
print("    " + "-" * 54)
for _ph, _fr in _cover.items():
    _mark = "  <-- in common set" if _fr >= BALANCE_COMMON_MIN_FRAC else ""
    print(f"    {_ph:<28}{100 * _fr:>24.1f}%{_mark}")
print(f"\n    common set: {len(COMMON_BALANCE_PHENOS)} phenotypes, "
      f"threshold {100 * BALANCE_COMMON_MIN_FRAC:.0f}% of {_n_struct_total} structures")

if COMMON_BALANCE_PHENOS:
    _sel = _ok.loc[_ok["pheno"].isin(COMMON_BALANCE_PHENOS)]
    _ref_common = _sel.groupby(["sample_id", STRUCT_COL])["mean"].mean()
    _n_common = _sel.groupby(["sample_id", STRUCT_COL])["pheno"].nunique()
    _key = list(zip(cells["sample_id"], cells[STRUCT_COL].astype(int)))
    cells["balance_reference_common"] = [
        _ref_common.get(k, np.nan) for k in _key]
    cells["n_common_phenos_used"] = [int(_n_common.get(k, 0)) for k in _key]
    cells.loc[cells["is_core"], "radial_centered_core_balanced_common"] = (
        cells.loc[cells["is_core"], "radial_pos"]
        - cells.loc[cells["is_core"], "balance_reference_common"])
    _short = _n_common.loc[_n_common < len(COMMON_BALANCE_PHENOS)]
    print(f"    structures holding the FULL common set: "
          f"{int((_n_common == len(COMMON_BALANCE_PHENOS)).sum())} of {_n_struct_total}")
    if len(_short):
        print(f"    {len(_short)} structure(s) have only part of it; their reference")
        print("    is the mean over the part they have, and n_common_phenos_used")
        print("    records how many. Read that column before quoting a structure.")
else:
    print("    NO COMMON SET at this threshold. The common-set centring is")
    print("    skipped and only the adaptive one is available; say so in Methods")
    print("    rather than presenting the adaptive reference as composition-free.")

# ---- REVISION 3: geometry provenance ----------------------------------------
# Structure ids changed between rev4 and rev5. A stale table 35 half-matches in
# silence: the unmatched keys map to NaN, the microns conversion then takes a
# median over whatever did match, and no line of output says so. This is the
# same failure mode script 07 item 11 guards against on the shared null,
# applied to the table this script actually reads.
_loaded_keys = {(r_s, int(r_k)) for r_s, r_k
                in cells[["sample_id", STRUCT_COL]].drop_duplicates().itertuples(
                    index=False, name=None)}
_geom_keys = set(equiv_of.keys())
if _loaded_keys:
    _matched = len(_loaded_keys & _geom_keys)
    _frac = _matched / len(_loaded_keys)
    print(f"\n    geometry provenance: {_matched} of {len(_loaded_keys)} loaded "
          f"structure keys found in table 35 ({100 * _frac:.1f}%)")
    if _frac < GEOMETRY_PROVENANCE_MIN_MATCH:
        print(f"    MISMATCH. Under {100 * GEOMETRY_PROVENANCE_MIN_MATCH:.0f}% of the "
              f"structures this script loaded")
        print(f"    have geometry in {FOCI_TABLE}.")
        print("    That table is from a different structure revision than the cell")
        print("    assignments. The radial coordinate and every microns conversion")
        print("    below would be silently wrong.")
        if REQUIRE_GEOMETRY_PROVENANCE:
            print("    REFUSING TO RUN. Set REQUIRE_GEOMETRY_PROVENANCE = False to")
            print("    override, having understood what that means.")
            sys.stdout = _tee.terminal; _tee.close(); sys.exit(1)
        print("    CONTINUING ANYWAY because REQUIRE_GEOMETRY_PROVENANCE is False.")

# ---- REVISION 3: effective independent units --------------------------------
# Matches script 07. No table below may imply that tens of thousands of cells
# are tens of thousands of pieces of information about a treatment that was
# assigned to six animals.
print("\n    EFFECTIVE INDEPENDENT UNITS")
print("    Treatment was assigned to ANIMALS. Cell counts are large; the number")
print("    of independent units is six, three per arm.")
for _c in CONDITION_ORDER:
    _sids = [x for x in SAMPLE_ORDER if COND_OF[x] == _c]
    _sel = cells.loc[cells["sample_id"].isin(_sids)]
    _ns = {short_label(x): int(_sel.loc[_sel["sample_id"] == x, STRUCT_COL].nunique())
           for x in _sids}
    print(f"      {_c:<12} animals {len(_sids)}   structures "
          f"{sum(_ns.values()):>3}   cells in structures {len(_sel):>9,}")
    print(f"                   per animal: {_ns}")
_single = [short_label(x) for x in SAMPLE_ORDER
           if int(cells.loc[cells['sample_id'] == x, STRUCT_COL].nunique()) == 1]
if _single:
    print(f"\n    {_single} contribute a SINGLE structure, so the")
    print("    structure-level variance is estimated almost entirely from the")
    print("    other arm and applied to both. Watch fit_mode in table 72.")

# ---- REVISION 3.1: which animal is the arm? ---------------------------------
# The mixed model's animal random intercept absorbs each animal's BASELINE, but
# the arm coefficient is still estimated from cells, so an arm dominated by one
# animal is largely a statement about that animal. On the first rev5 run 43118
# carried 3,273 of 4,264 treated lymphocytes, 77 percent, and it is also the
# animal holding the 1.59 shape ratio and the 3.6 percent IDO1 outlier. That is
# the strongest argument for p_exact_means as the headline, since permuting
# animal-level means weights the three animals equally. Printed rather than
# left for a reviewer to notice.
print("\n    ARM COMPOSITION BY ANIMAL, core lymphocytes")
_lym_core = cells.loc[cells["is_core"] & cells["pheno"].isin(LYMPHOCYTES)]
for _c in CONDITION_ORDER:
    _sids = [x for x in SAMPLE_ORDER if COND_OF[x] == _c]
    _tot = int(len(_lym_core.loc[_lym_core["sample_id"].isin(_sids)]))
    if not _tot:
        continue
    _shares = []
    for _x in _sids:
        _nx = int(len(_lym_core.loc[_lym_core["sample_id"] == _x]))
        _shares.append((short_label(_x), _nx, 100.0 * _nx / _tot))
    _shares.sort(key=lambda t: -t[2])
    _txt = ",  ".join(f"{a} {n:,} ({q:.0f}%)" for a, n, q in _shares)
    print(f"      {_c:<12} {_tot:>7,} cells   {_txt}")
    if _shares and _shares[0][2] >= 60.0:
        print(f"      ^ {_shares[0][0]} is {_shares[0][2]:.0f}% of this arm. The "
              f"cell-weighted coefficient is")
        print("        largely a statement about that one animal. Lead with "
              "p_exact_means,")
        print("        which weights animals equally, and read leave-one-out "
              "for that animal.")

core = cells.loc[cells["is_core"]].copy()
core["lineage"] = np.where(core["pheno"].isin(T_LINEAGE), "T lineage",
                           np.where(core["pheno"].isin(B_LINEAGE), "B lineage",
                                    "other"))
lym = core.loc[core["pheno"].isin(LYMPHOCYTES)].copy()

CENTRING_OUTCOMES = [("radial_centered_core", "primary (cell-weighted)")]
if "radial_centered_core_balanced" in core.columns:
    CENTRING_OUTCOMES.append(("radial_centered_core_balanced",
                              "balanced (phenotype-weighted)"))
if "radial_centered_core_balanced_common" in core.columns:
    CENTRING_OUTCOMES.append(("radial_centered_core_balanced_common",
                              "balanced, common set"))
if HAVE_AREA_REF and "radial_centered_core_area" in core.columns:
    CENTRING_OUTCOMES.append(("radial_centered_core_area",
                              "area-weighted (exact)"))

# ---- REVISION 4: which centring is confirmatory -----------------------------
_available = [c for c, _ in CENTRING_OUTCOMES]
CONFIRMATORY_CENTRING = next(
    (c for c in CONFIRMATORY_CENTRING_PREFERENCE if c in _available),
    "radial_centered_core")
CONFIRMATORY_NAME = dict(CENTRING_OUTCOMES).get(CONFIRMATORY_CENTRING, "?")

sub("Confirmatory centring")
print("    One centring carries the q-values; the rest are sensitivity")
print("    analyses on the same hypothesis and get none. The rule was fixed")
print("    before this run, not chosen from the results.\n")
for c, nm in CENTRING_OUTCOMES:
    tag = "  <-- CONFIRMATORY" if c == CONFIRMATORY_CENTRING else "      sensitivity"
    print(f"    {nm:<40}{tag}")
if CONFIRMATORY_CENTRING == "radial_centered_core":
    print("\n    WARNING: falling back to the CELL-WEIGHTED centring as")
    print("    confirmatory, because no composition-free reference is")
    print("    available. That reference moves with composition and the arms")
    print("    differ sharply in composition. Treat every q-value below as")
    print("    provisional until a balanced or area-weighted reference exists.")
elif CONFIRMATORY_CENTRING != "radial_centered_core_area":
    print("\n    NOTE: the exact area-weighted reference is not available,")
    print("    because table 35 lacks mean_radial_over_area. Script 04 exports")
    print("    it from revision 6 onward. Until then the common set is")
    print("    confirmatory, and it is small and myeloid-weighted.")

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

# ---- REVISION 3: thin-population report, for setting a gate later -----------
# Script 07 carries MIN_TARGETS_FOR_REPORT = 50 and labels any separation whose
# smallest treated target count falls under it NOT REPORTABLE. This script has
# no equivalent gate YET, deliberately: the threshold should be set from the
# radial counts rather than imported from a nearest-neighbour analysis before
# anyone has seen them. This block prints what the decision needs.
sub("Thin populations, read before quoting any per-phenotype radial result")
print("    No gate is applied here. These are the counts a gate would be set")
print("    from, printed so the threshold is chosen on evidence rather than")
print("    inherited. Two things retire a population: too few cells in the")
print("    thinnest treated animal, and absence from an animal altogether,")
print("    which silently drops that animal from the model for that phenotype.\n")
print(f"    {'phenotype':<20}{'min treated':>13}{'min untreated':>15}"
      f"{'animals present':>17}")
print("    " + "-" * 65)
_thin_rows = []
for _p in LYMPHOCYTES:
    _per = {x: int(tab.loc[x, _p]) if _p in tab.columns else 0 for x in SAMPLE_ORDER}
    _t = [v for x, v in _per.items() if COND_OF[x] == "D1MT"]
    _u = [v for x, v in _per.items() if COND_OF[x] == "Untreated"]
    _present = sum(1 for v in _per.values() if v > 0)
    _flag = ""
    if _present < len(SAMPLE_ORDER):
        _flag = f"  <-- missing from {len(SAMPLE_ORDER) - _present} animal(s)"
    print(f"    {_p:<20}{min(_t) if _t else 0:>13,}{min(_u) if _u else 0:>15,}"
          f"{_present:>13} of {len(SAMPLE_ORDER)}{_flag}")
    _thin_rows.append({"phenotype": _p, "min_treated": min(_t) if _t else 0,
                       "min_untreated": min(_u) if _u else 0,
                       "animals_present": _present,
                       "n_animals_total": len(SAMPLE_ORDER)})
write_csv(pd.DataFrame(_thin_rows), "70b_thin_population_report.csv")


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
        print(f"    p (model)    : {r['p_value']:.4f}")
        print(f"    p_exact_means: {r['p_exact_means']:.3f}   <-- HEADLINE, "
              f"floored at 0.10 by the design")
        print(f"    p_exact_lmm  : {r['p_exact_lmm']:.3f}   (sensitivity only)")
        print(f"    cells        : {r['n_cells']:,} in {r['n_structures']} structures "
              f"({r['n_structures_D1MT']} treated, {r['n_structures_ref']} untreated)")
        print(f"    animals      : {r['n_animals_D1MT']} treated, "
              f"{r['n_animals_ref']} untreated")
        print(f"    random effect: {r['random_effects']}   [fit_mode: {r['fit_mode']}]")
        if np.isfinite(r['p_value']) and r['p_value'] < 0.10 <= r['p_exact_means']:
            print("\n    NOTE: the model p sits below the exact floor of 0.10 while")
            print("    the exact p is at it. That gap is the asymptotic assumption")
            print("    the model is making, not extra evidence. Lead with the")
            print("    effect size and the per-animal values.")

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
                r["family"] = ("lineage" if col == CONFIRMATORY_CENTRING
                               else "lineage_sensitivity")
                r["lineage"] = name
                r["centring"] = cname
                model_rows.append(r)
                print(f"    {name:<12} {cname:<32} "
                      f"coef={r['coef_D1MT_vs_ref']:>+8.4f}  "
                      f"p={r['p_value']:.4f}  cells={r['n_cells']:>7,}")

    # ---- per cell type, circular and non-circular kept apart ---------------
    # REVISION 3.1: run on EVERY centring, not just the cell-weighted one.
    # Revision 3 reported this family on the primary centring alone and printed
    # "5 of 5 at q < 0.1" for the lymphocytes. Those were the most quotable
    # numbers in the output and they were computed on the reference the centring
    # check had just shown moves with composition: pooled lymphocytes fell from
    # -0.1454 to -0.0622 between the two. A per-phenotype q-value on the primary
    # centring is therefore not the composition-free result it reads as.
    # Balanced leads; primary is kept beside it as the comparison.
    _min_treated = {}
    for _p in LYMPHOCYTES + [IDO1_POS, IDO1_NEG, "Neutrophils"]:
        _d = core.loc[core["pheno"] == _p]
        _per = _d.groupby("sample_id").size()
        _t = [int(_per.get(x, 0)) for x in SAMPLE_ORDER if COND_OF[x] == "D1MT"]
        _min_treated[_p] = min(_t) if _t else 0

    _centrings_for_pheno = [c for c in CENTRING_OUTCOMES
                            if c[0] != "radial_centered_core"] + \
                           [("radial_centered_core", "primary (cell-weighted)")]
    # REVISION 4: confirmatory first, so the block that carries the q-values is
    # the one read first, and the cell-weighted comparison stays last.
    _centrings_for_pheno = (
        [c for c in _centrings_for_pheno if c[0] == CONFIRMATORY_CENTRING]
        + [c for c in _centrings_for_pheno if c[0] != CONFIRMATORY_CENTRING])
    for _col, _cname in _centrings_for_pheno:
        _is_primary = _col == "radial_centered_core"
        _is_conf = _col == CONFIRMATORY_CENTRING
        sub(f"Per cell type (core only, {_cname})"
            + ("   [CONFIRMATORY, carries the q-values]" if _is_conf else
               "   [composition-contaminated, for comparison]" if _is_primary else
               "   [sensitivity, no q-values]"))
        for p_ in LYMPHOCYTES + [IDO1_POS, IDO1_NEG, "Neutrophils"]:
            d = core.loc[core["pheno"] == p_]
            circ = p_ in DETECTION_POOL
            r = fit_mixed(capped(d), _col, f"radial: {p_} ({_cname})",
                          note="CIRCULAR: defines the foci" if circ else "")
            if not r:
                continue
            _suffix = "" if _col == CONFIRMATORY_CENTRING else "_sensitivity"
            r["family"] = ("per_phenotype_circular" if circ
                           else "per_phenotype_lymphocyte") + _suffix
            r["phenotype"] = p_
            r["circular"] = circ
            r["centring"] = _cname
            # REVISION 3.1 reportability label, never a filter
            r["min_treated_cells"] = _min_treated.get(p_, 0)
            r["reportable"] = bool(_min_treated.get(p_, 0)
                                   >= MIN_TREATED_CELLS_FOR_REPORT)
            model_rows.append(r)
            _flags = ""
            if circ:
                _flags += "  [circular]"
            if not r["reportable"]:
                _flags += (f"  [NOT REPORTABLE: {r['min_treated_cells']} cells in "
                           f"the thinnest treated animal]")
            # REVISION 5: print p_exact_means here. Revision 4 printed the
            # MODEL p only, in the most quotable table in the output, after the
            # script had already declared that the exact animal-level p leads.
            # It mattered: on the 27 September run plasma cells survived BH at
            # p_model = 0.0544 while showing no per-animal separation on the
            # same centring, so its exact p cannot be better than the 0.10
            # floor. A reader working from p_model alone would have quoted it.
            # REVISION 6: THE REVISION 5 FLAG FIRED ON THE WRONG ROWS.
            # It tested p_model < 0.05 <= p_exact, which on a design floored at
            # 0.10 trips whenever the model p is small and the exact p is AT
            # the floor. At the floor the six animal means separate completely,
            # which is the strongest result this design can produce, so the
            # cleanest row in the table carried a warning. Meanwhile the row
            # that needed one did not get it: on the 27 September run plasma
            # cells came back p_model = 0.0544 with p_exact = 0.200, so the
            # animals do NOT separate, and 0.0544 is not below 0.05 so nothing
            # printed.
            #
            # What matters is not the model p. It is whether the exact p sits
            # at the floor or above it.
            _pe = r.get("p_exact_means", np.nan)
            _fl = r.get("p_exact_floor", np.nan)
            _gap = ""
            if np.isfinite(_pe) and np.isfinite(_fl):
                if _pe <= _fl + 1e-9:
                    _gap = "  [animals separate]"
                else:
                    _gap = (f"  <-- ANIMALS DO NOT SEPARATE "
                            f"(exact {_pe:.3f} > floor {_fl:.3f})")
            print(f"    {p_:<26} coef={r['coef_D1MT_vs_ref']:>+8.4f}  "
                  f"p_model={r['p_value']:.4f}  p_exact="
                  f"{_pe:.3f}  cells={r['n_cells']:>7,}{_flags}{_gap}")
    print("\n    Circular and non-circular populations are in SEPARATE BH")
    print("    families. Mixing them inflated m and blended a descriptive")
    print("    family with a confirmatory one.")
    print(f"\n    REVISION 4: only the CONFIRMATORY centring, {CONFIRMATORY_NAME},")
    print("    carries q-values. Revision 3.1 pooled every balanced centring into")
    print("    one family, which made five phenotypes measured twice look like")
    print("    ten hypotheses and inflated m. One hypothesis on several")
    print("    references is a sensitivity series, not a set of tests.")
    print(f"\n    NOT REPORTABLE below {MIN_TREATED_CELLS_FOR_REPORT} cells in the "
          f"thinnest treated animal.")
    print("    On this input that retires Tregs and Helper T cells, which were")
    print("    already discounted on other grounds: six Tregs in 43106, and the")
    print("    CD4-call problem for Helper T cells. A population can look")
    print("    significant and survive BH on a handful of cells, and Tregs does.")

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
# REVISION 4: every non-confirmatory centring is a sensitivity family and gets
# no q-values, so BH is never applied across one hypothesis measured several
# ways. See CONFIRMATORY_CENTRING_PREFERENCE in Cell 1.
NO_BH_FAMILIES = NO_BH_FAMILIES + [
    f for f in (set(models["family"]) if len(models) else [])
    if str(f).endswith("_sensitivity")]
if len(models):
    models["q_value"] = np.nan
    for fam, g in models.groupby("family"):
        if fam in NO_BH_FAMILIES:
            continue
        models.loc[g.index, "q_value"] = benjamini_hochberg(g["p_value"].to_numpy())
    models["gets_q_value"] = ~models["family"].isin(NO_BH_FAMILIES)
    write_csv(models, "72_models_all.csv")

    _fb = models.loc[models["fit_mode"] == "animal only"] if "fit_mode" in models else []
    if len(_fb):
        sub("Fits that fell back to animal-only random effects")
        print("    REVISION 3: revision 2 did not record fit_mode, so a fallback")
        print("    looked like any other row. The nested structure term could not")
        print("    be estimated for these, so they answer a slightly different")
        print("    question from the rest of the table.\n")
        for _, _r in _fb.iterrows():
            print(f"    {_r['analysis'][:58]:<60} {_r['family']}")

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
    # REVISION 3.1: the per-phenotype families now exist on more than one
    # centring, so naming them literally would plot each phenotype several
    # times. Prefer the BALANCED rows, fall back to the primary ones, and drop
    # anything the reportability label retired so a six-cell population cannot
    # reach a figure.
    # REVISION 4: the confirmatory families are the unsuffixed ones, so this
    # plots each phenotype ONCE. Revision 3.1 named the balanced families
    # explicitly, which would have drawn every phenotype three times as soon as
    # a third centring existed.
    _fams = set(models["family"])
    _pheno_fams = [f for f in ["per_phenotype_lymphocyte",
                               "per_phenotype_circular"] if f in _fams]
    g = models.loc[models["family"].isin(["primary", "lineage"] + _pheno_fams)]
    if "reportable" in g.columns:
        # REVISION 4: cast explicitly. fillna on an object column raises a
        # pandas FutureWarning about silent downcasting.
        _rep = g["reportable"].map(lambda v: True if pd.isna(v) else bool(v))
        g = g.loc[_rep.to_numpy(dtype=bool)]
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


# %% Cell 8 - REMOVED IN REVISION 3
# =============================================================================
# This cell fitted delta-distance models and wrote table 74. It is gone, and
# the reason is worth keeping because it was not obvious.
#
# WHAT IT DID WRONG. It computed observed nearest-neighbour distances over
# CORE CELLS ONLY, at core.groupby(["sample_id", STRUCT_COL]), and then
# subtracted null_median_um straight out of script 06 table 53. Script 06
# builds those nulls from d.loc[d[STRUCT_COL] > 0], which is EVERY cell in the
# structure, core and cuff. Script 07 does the same. So the delta here was
# observed(core only) minus null(core plus cuff): two different cell
# populations, and the mismatch is not random. A null that includes cuff cells
# has more targets spread over a larger area, so it does not land where a
# core-only null would. Worse, the cuff is a FIXED 150 um band while the core
# is size-normalised, so the cuff-to-core ratio varies strongly with focus
# size. That is the size confound walking back in through the one construction
# that exists to remove it.
#
# The revision 2 header said "the null comes from script 06 table 53, so it is
# identical in scripts 06, 07 and 08". The null was indeed identical. What it
# was being subtracted from was not, and the sentence hid that.
#
# WHY REMOVED RATHER THAN REPAIRED. Two reasons.
#   1. Script 07 already owns the delta-distance family and reports it as the
#      single confirmatory outcome, on all structure cells, against the
#      matching null. Repairing this cell would produce a SECOND set of
#      numbers for the same pairs in a second script. That is precisely how
#      script 05 came to be frozen, and how a number gets quoted from the
#      wrong place.
#   2. It also fixes an inference problem. Revision 2 applied
#      Benjamini-Hochberg to delta_distance_pair as its own family while the
#      radial families got their own q-values, which treats distance and radial
#      as independent. Script 07 revision 3 established they are not: both
#      nearest-neighbour anchors sit in the myeloid pool that defines a focus
#      and therefore at the density peak, so anchor-to-non-pool-target distance
#      is largely radial position in different units. With the distance cell
#      gone, this script tests one spatial fact and the question does not arise.
#
# Script 08 is now purely radial. For any distance number, read script 07.
# If a core-restricted distance view is ever genuinely wanted, the correct
# route is to have script 06 export a second, core-only null column, not to
# mix populations here.


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
                # REVISION 3.1: WHICH ARM IS HIGHER. Spillover rises with cell
                # density and untreated foci are three to five fold denser, so
                # the spillover hypothesis predicts HIGHER apparent
                # co-expression in UNTREATED. A control that separates the
                # other way is not evidence of spillover, whatever its count.
                higher = ("untreated" if b.median() > a.median() else "treated")
                sep_rows.append({"pair": pair, "kind": kind, "percentile": pct,
                                 "d1mt_min": a.min(), "d1mt_max": a.max(),
                                 "untr_min": b.min(), "untr_max": b.max(),
                                 "complete_separation": sep,
                                 "arm_higher": higher,
                                 "direction_matches_spillover": higher == "untreated"})
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
        # ---- REVISION 3.1: direction ---------------------------------------
        sub("Direction, which revision 3 did not check")
        print("    Spillover rises with cell density and untreated foci are three")
        print("    to five fold denser, so spillover predicts HIGHER apparent")
        print("    co-expression in UNTREATED. Counting separations without")
        print("    checking their direction can retire a result on a control that")
        print("    is pointing the wrong way. On the 21 September rev5 run the")
        print("    test pair ran higher in untreated, as spillover predicts, while")
        print("    CD3e x CD20 ran higher in TREATED, which spillover does not.\n")
        _sep_only = sepdf.loc[sepdf["complete_separation"]]
        print(f"    {'pair':<26}{'kind':>9}{'separations':>13}"
              f"{'higher in':>12}{'fits spillover':>16}")
        print("    " + "-" * 76)
        _ctrl_agree = 0
        for (_pair, _kind), _g in sepdf.groupby(["pair", "kind"]):
            _gs = _g.loc[_g["complete_separation"]]
            _n = int(len(_gs))
            if _n:
                _hi = _gs["arm_higher"].mode().iat[0]
                _fits = bool((_gs["direction_matches_spillover"]).mode().iat[0])
            else:
                _hi, _fits = "-", False
            if _kind == "control" and _n and _fits:
                _ctrl_agree += 1
            print(f"    {_pair:<26}{_kind:>9}{_n:>9} of {n_thr}"
                  f"{_hi:>12}{('yes' if _n and _fits else 'no'):>16}")
        _test_g = sepdf.loc[(sepdf["kind"] == "test") & sepdf["complete_separation"]]
        _test_fits = bool(_test_g["direction_matches_spillover"].mode().iat[0]) \
            if len(_test_g) else False
        print(f"\n    controls separating IN THE SPILLOVER DIRECTION: {_ctrl_agree} "
              f"of {n_ctrl_pairs}")

        if ctrl_max >= test_n and test_n > 0 and _ctrl_agree > 0:
            print("\n    VERDICT: at least one control pair separates the arms as")
            print("    well as or better than the test pair, AND does so in the")
            print("    direction spillover predicts. A macrophage cannot be both a")
            print("    T and a B cell, so the test pair's separation is not")
            print("    established as biology. DROP IT, or report it only with")
            print("    this control alongside.")
        elif ctrl_max >= test_n and test_n > 0 and _ctrl_agree == 0:
            print("\n    VERDICT: UNRESOLVED, and the revision 3 verdict was too")
            print("    strong. A control pair does separate the arms as well as or")
            print("    better than the test pair, so the panel demonstrably")
            print("    produces arm-separating co-expression artefacts and the")
            print("    test pair cannot be called clean. But NO control separates")
            print("    in the direction spillover predicts, so the separating")
            print("    control is not an estimate of the artefact this test pair")
            print("    would suffer from. Revision 3 counted separations without")
            print("    checking sign and retired the result on that count.")
            print("    Do not reinstate the finding on this either. What it means")
            print("    is that the control does not settle the question, and")
            print("    something other than density-driven spillover is separating")
            print("    the arms on a pair that cannot co-occur. That needs")
            print("    explaining before either pair is quoted.")
            if _test_fits:
                print("    Note the test pair DOES run in the spillover direction,")
                print("    which is the one piece of evidence still against it.")
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
    # REVISION 4: the headline is the CONFIRMATORY centring, not the
    # cell-weighted one. Revision 3.2 printed the cell-weighted pooled model
    # under "PRIMARY RESULT" on the same page as a centring check showing that
    # reference loses more than half the coefficient. The confirmatory fit
    # lives in the centring_check family, which fits the pooled model on every
    # reference, so it is pulled out by name here.
    _cc_all = models.loc[models["family"] == "centring_check"]
    _conf = _cc_all.loc[_cc_all.get("centring", pd.Series(dtype=object))
                        == CONFIRMATORY_NAME] if len(_cc_all) else _cc_all
    if len(_conf):
        rr = _conf.iloc[0]
        print(f"HEADLINE, CONFIRMATORY CENTRING: {CONFIRMATORY_NAME}")
        print(f"  Pooled lymphocytes, core cells, centered radial position")
        print(f"  coef {rr['coef_D1MT_vs_ref']:+.4f} "
              f"[{rr['ci_low']:+.4f}, {rr['ci_high']:+.4f}]")
        print(f"  p_model {rr['p_value']:.4f}   p_exact_means "
              f"{rr['p_exact_means']:.3f} (LEADS)   p_exact_lmm "
              f"{rr['p_exact_lmm']:.3f}")
        print(f"  {rr['n_cells']:,} cells, {rr['n_structures']} structures, "
              f"fit_mode {rr['fit_mode']}")

        # REVISION 5: state whether the POOLED claim survives on the
        # confirmatory reference, rather than leaving the reader to compare a
        # confidence interval against a p-value. On the 27 September run this
        # block printed coef -0.0784 with an interval crossing zero and
        # p_exact_means of 0.300, under the heading HEADLINE, and said nothing
        # about what that means.
        _ci_crosses = bool(np.isfinite(rr["ci_low"]) and np.isfinite(rr["ci_high"])
                           and rr["ci_low"] < 0.0 < rr["ci_high"])
        _exact_ns = bool(np.isfinite(rr["p_exact_means"])
                         and rr["p_exact_means"] > 0.10 + 1e-9)
        if _ci_crosses or _exact_ns:
            print("\n  THE POOLED CLAIM DOES NOT SURVIVE ON THIS REFERENCE.")
            if _ci_crosses:
                print("    The 95% interval crosses zero.")
            if _exact_ns:
                print(f"    p_exact_means is {rr['p_exact_means']:.3f}, above the "
                      f"0.10 floor, so the six")
                print("    animal-level means do not separate.")
            print("    'Lymphocytes sit more core-ward in treated foci' is")
            print("    therefore NOT supported as a pooled confirmatory claim.")
            print("    What the confirmatory family supports is the PER-PHENOTYPE")
            print("    result; read that block and lead with the population, not")
            print("    with pooled lymphocytes. The pooled effect on the")
            print("    cell-weighted reference below is the composition-")
            print("    contaminated version and is not a fallback.")
        else:
            print("\n  The pooled claim survives on this reference: the interval")
            print("  excludes zero and the animal-level means separate.")
        print()

    prim = models.loc[models["family"] == "primary"]
    if len(prim):
        r = prim.iloc[0]
        print(f"CELL-WEIGHTED CENTRING, for comparison only")
        print(f"  Pooled lymphocytes, core cells, centered radial position")
        print(f"  coef {r['coef_D1MT_vs_ref']:+.4f} "
              f"[{r['ci_low']:+.4f}, {r['ci_high']:+.4f}]")
        print(f"  p_model {r['p_value']:.4f}   p_exact_means "
              f"{r['p_exact_means']:.3f} (HEADLINE)   p_exact_lmm "
              f"{r['p_exact_lmm']:.3f}")
        print(f"  {r['n_cells']:,} cells, {r['n_structures']} structures "
              f"({r['n_structures_D1MT']} treated, {r['n_structures_ref']} untreated), "
              f"{r['n_animals_D1MT']} vs {r['n_animals_ref']} animals")
        print(f"  fit_mode {r['fit_mode']}")
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
print("    equivalent radius overstates by the shape ratio, which on rev5 is")
print("    1.34 treated and 1.24 untreated. REVISION 3 CORRECTION: this line")
print("    previously read '1.13 treated and 1.22 untreated', which was wrong")
print("    twice. Those were rev4-era figures, and 1.13 is not a shape ratio at")
print("    all, it is the rev5 untreated-to-treated EQUIVALENT RADIUS gap.")
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
