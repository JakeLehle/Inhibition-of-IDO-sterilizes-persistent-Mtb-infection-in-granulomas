# Inhibition-of-IDO-sterilizes-persistent-Mtb-infection-in-granulomas

# AKOYA pipeline run order

Breadcrumbs for the spatial proteomics arm of *Inhibition of IDO sterilizes
persistent Mtb infection in granulomas*. Rhesus Mtb + SIV, D1MT-treated (G3)
against untreated (G4), necropsy lung sections, Akoya Phenocycler 67-plex.

Written 2 October 2026, against the scripts in `scripts/` at that date. If a
script has been revised since, trust the script and fix this file.

**What this file is for.** Picking the pipeline back up cold, and knowing which
script to run, in what order, reading what, writing where, and which ones are
frozen or superseded. Read Part 1 for the order and Part 4 for the hazards
that have actually cost time.

---

# Part 0. Where things live

All paths are on titan under `/master/jlehle/WORKING/AKOYA/`.

| | |
|---|---|
| raw input | `data/*.csv`, the Akoya per-cell export, one CSV per section |
| scripts | `scripts/` |
| per-stage output | one directory per stage, named in the table below |
| within a stage | `tables/` for CSV and the `00_*_report.txt` log, `figures/` for PDF and PNG |

Every script writes a `00_<stage>_report.txt` that tees stdout. That report is
the first thing to read when picking a stage back up, because it carries the
run time, the parameter block and the warnings.

**Raw imagery was never delivered by the vendor.** Only the processed per-cell
export exists, which is spatial coordinates plus Akoya's own phenotype calls.
No qptiff is held on titan or in Dropbox. Nothing Akoya's classifier did can be
independently checked, which is also why the vendor IDO1 binary call cannot be
reproduced as an intensity threshold. That belongs in Methods as a stated
constraint and it bounds what any script here can validate.

**Environment.** Conda env `sc_pre`, Python 3.11. Interactive work in Spyder
using the `# %%` cell format. Pipeline runs under SLURM use
`matplotlib.use("Agg")`. Tunable parameters are in Cell 1 of each script.

---

# Part 1. The canonical order

Cold start from the raw export. Stages 1 to 3 are run once and rarely touched.
Stage 4 is the detection layer everything else depends on. Stages 6 to 8 are
the analysis layer. Stages 9 to 11 consume the analysis layer.

| # | script | reads | writes | on critical path |
|---|---|---|---|---|
| 1 | `AKOYA_01_Inventory.py` | `data/` | `inventory/` | yes, entry point |
| 1b | `AKOYA_01b_Inventory_Figures.py` | `inventory/tables/` | `inventory/figures/` | no, figures only |
| 2 | `AKOYA_02_Spatial_Diagnostic.py` | `inventory/tables/` | `diagnostics/` | yes |
| 3 | `AKOYA_03_Rebase.py` | `inventory/`, `diagnostics/` | `rebaseline/` | yes |
| 4 | `AKOYA_04_Structures_rev6.py` | `data/`, `inventory/tables/` | `structures_rev5/` | **yes, the base** |
| 6 | `AKOYA_06_Distance_Statistics.py` | `structures_rev5/` | `distance_stats/` | yes |
| 7 | `AKOYA_07_Size_Corrected.py` | `structures_rev5/` | `size_corrected/` | yes |
| 7b | `AKOYA_07b_Composition_Null.py` | `structures_rev5/`, `distance_stats/` | `composition_null/` | yes |
| 8 | `AKOYA_08_Lymphocyte_Radial.py` | `structures_rev5/` | `lymphocyte_radial/` | yes |
| 11 | `AKOYA_11_Claim_Checks.py` | `structures_rev5/`, `composition_null/`, `lymphocyte_radial/`, `data/` | `claim_checks/` | yes |
| 10 | `AKOYA_10_Model_Diagnostics.py` | structures dir, `lymphocyte_radial/tables/72_models_all.csv` | `model_diagnostics/` | no, read-only diagnostics |
| 9 | `AKOYA_09_Publication_Figures.py` | structures dir, `lymphocyte_radial/tables/` | `publication_figures/` | yes, last |

`akoya_arm_stats.py` is a module, not a stage. It is imported by scripts 09 and
10 only. It holds the Datta-Satten clustered rank-sum, the centred radial
construction, the violin plotting, and `fit_arm_lmm`.

**Short form.** 01 → 01b → 02 → 03 → 04 rev6 → 06 → 07 → 07b → 08 → 11 → 10 → 09

6, 7 and 8 all read `structures_rev5` directly and do not depend on each other,
so they can run in any order or in parallel. 7b is the exception: it validates
its reconstruction against script 06's `53_nn_per_structure.csv` at r = 0.9995,
so 06 must have run first. 11 needs both 07b and 08. 9 and 10 need 08.

---

# Part 2. Frozen, superseded and off-path scripts

Do not run these as part of the pipeline. They are kept because they document
how the detection rule was arrived at, and deleting them would make the
detector look like an assertion rather than something that was scored.

| script | status | why it is kept |
|---|---|---|
| `AKOYA_04_Structures.py` | **superseded** by rev5 and rev6 | the original dual structure definition, writes `structures_rev4`. Still the input for the 04b to 04f audit branch |
| `AKOYA_04_Structures_rev5.py` | **superseded** by rev6 | triple definition. rev6 is rev5 plus four columns and nothing else |
| `AKOYA_04b_Foci_Audit.py` | audit branch, reads rev4 | the candidate ledger, background statistics and gate grid behind the detection thresholds |
| `AKOYA_04c_Annotation_Canvas.py` | audit branch, reads rev4 | builds the canvas the expert annotations were drawn on |
| `AKOYA_04d_Annotation_Measurements.py` | audit branch, reads rev4 | measures the expert annotations |
| `AKOYA_04e_Candidate_Detector.py` | audit branch, reads rev4 | scores candidate detectors against the expert annotations. This is where the rev5 parameters come from |
| `AKOYA_04f_Object_Review.py` | audit branch, reads rev4 | composition of each detected object, review canvas |
| `AKOYA_05_Radial_Architecture.py` | **FROZEN** | superseded by 06 for distance and 08 for radial position. Do not quote anything from `radial/` |

The 04b to 04f branch reads `structures_rev4` and was the development path for
the detector. It does not need rerunning on rev5, because its job was to choose
the detection parameters that rev5 and rev6 now implement. If the detection
rule is ever changed, this branch is what has to be rerun to re-score it.

---

# Part 3. Current state, 2 October 2026

The pipeline is mid-rerun. This is what has and has not happened.

| script | revision | state |
|---|---|---|
| 01, 01b, 02, 03 | current | done, stable |
| 04 rev6 | rev 6 | **built, not yet run.** Adds `mean_radial_over_area` and three companions to table 35 |
| 06 | rev 4 | done, on rev5 |
| 07 | rev 3.5 | done and locked, on rev5 |
| 07b | rev 2 | done, on rev5 |
| 08 rev6 | rev 6 | **built, not yet run.** Needs 04 rev6 first for the area-weighted reference to resolve |
| 09 | rev 4 | **stale, and points at `structures_rev4`.** See Part 4 item 3 |
| 10 | current | built, not run. Also points at `structures_rev4` |
| 11 | rev 1 | done, on rev5 |

**Pending changes agreed 2 October**, from
`claude/AKOYA_method_decision_2026-10-02.md`.

1. The denominator degrees-of-freedom fix, in five places. `akoya_arm_stats.py`
   `fit_arm_lmm`, plus the local `fit_mixed` in 06, 07, 07b and 08. Patch in
   `akoya_small_sample_patch.py`. Scripts 06, 07, 07b and 08 need the
   `akoya_arm_stats` import added first, because they do not currently have it.
2. Repoint `STRUCT_DIR` in 09 and 10 from `structures_rev4` to
   `structures_rev5`.
3. Foci per cm² of segmented tissue in script 04.

**Rerun order after those changes**, which is not the canonical order because
detection is untouched by any of them.

```
04 rev6                  # first, writes structures_rev5. Verify numbers reproduce.
06, 07, 07b, 08          # any order, 07b after 06. All four carry the df patch.
11                       # claim checks, re-baselines the verdict table
10                       # model diagnostics, confirms the df change landed
09                       # figures last, after STRUCT_DIR is repointed
```

Stages 1 to 3 do not need rerunning. Nothing in them fits a model or depends on
structure detection.

---

# Part 4. Hazards that have actually cost time

**1. rev6 replaces rev5, it does not append to it.**
`AKOYA_04_Structures_rev6.py` reads the raw export from `data/` and recomputes
detection from scratch, then writes the full superset into
`structures_rev5/`. It is standalone, not a patch on rev5's output. So once
rev6 has run, running `AKOYA_04_Structures_rev5.py` again silently removes
`mean_radial_over_area` and the three companion columns, and script 08's
area-weighted centring stops resolving. **After rev6, never run rev5 again.**

**2. Detection must reproduce exactly when rev6 runs.** rev6 changes no
detection parameter, so every focus count, every peak density and every burden
fraction must come out identical to the rev5 run. If any of them moves, stop
and find out what changed, because something was edited that should not have
been.

**3. "Rerun 09 on rev5" is not just a rerun.** `STRUCT_DIR` is hardcoded to
`structures_rev4` at line 103 of script 09 and line 107 of script 10. Rerunning
either as-is reproduces rev4 numbers and looks like a successful rerun. Both
need the one-line path change. Checked and safe: rev5 and rev6 both write
`33b_burden_summary.csv`, `35_foci_structures_relative.csv` and
`cell_assignments/*_cell_structures.csv`, which are the three things 09 and 10
read from the structures directory, so repointing breaks nothing.

**4. Copy `00_lymphocyte_radial_report.txt` before rerunning 08.** The
common-set confirmatory results should stay on record beside the area-weighted
ones. Otherwise only the final set survives and it looks as though the centring
reference was chosen after seeing the results, which is the exact thing fixing
the preference order in advance was meant to prevent.

**5. `LEAVE_ONE_OUT` redirects the output directory.** Setting it in 04 rev5 or
rev6 appends `_no_<phenotype>` to `OUT_DIR`, so leave-one-out detection runs
land in a sibling directory and cannot clobber the main one. That is by design.
It also means a leave-one-out run does not update `structures_rev5`, so the
downstream scripts will not see it unless they are pointed at the sibling.

**6. A returned model object is not convergence.** Arm is constant within
animal, so with a random intercept alone the animal and arm effects compete and
the Hessian can go singular, which surfaces as a NaN standard error with no
error raised. Every fit path checks for a finite positive standard error and
says so out loud. A row whose `fit_mode` reads "animal only" is a selected fit:
simulated at this design the animal-only fallback converges in 8 to 35 percent
of null replicates while the nested animal-plus-structure fit converges in 99.6
percent. Read "animal only" rows as diagnostics, not results.

**7. System `samtools` at `/usr/local/bin/samtools` is broken.** Not used by
this pipeline, but if a SLURM job here ever grows a samtools call, activate
conda first and `export PATH="${CONDA_PREFIX}/bin:${PATH}"`.

**8. `csv.DictWriter` line endings.** Always pass `lineterminator="\n"`. The
default `\r\n` fails silently in bash pipelines downstream.

---

# Part 5. Which script owns which number

When a manuscript number needs tracing, this is where it comes from. Nothing
should be quoted from two places, and where two scripts compute the same thing
the owner is named.

| quantity | owner | table |
|---|---|---|
| sections retained, IDO1 cutoff window | 03 | `26_ido1_cutoff_windows.csv` |
| foci per animal, burden, peak density, inscribed and equivalent radius | 04 rev6 | `33b_burden_summary.csv`, `35_foci_structures_relative.csv` |
| BALT structures | 04 rev6 | `36_balt_structures.csv` |
| granuloma complexes | 04 rev6 | `37_granuloma_complexes.csv` |
| exact nearest-neighbour distances and the shared null | 06 | `53_nn_per_structure.csv` |
| size-corrected distance models, the delta outcome | 07 | `62_mixed_models_all_outcomes.csv`, `60_nn_per_structure_with_geometry.csv` |
| proximity fractions, anchor contrast | 07 | `67_proximity_fractions.csv`, `68_anchor_contrast.csv` |
| the three anchor references and per-animal deltas | 07b | `83_per_animal_delta_three_references.csv` |
| radial position models, the confirmatory family | 08 | `72_models_all.csv` |
| centring references | 08 | `77_centring_references.csv` |
| lymphocyte counts in core, the reportability gate | 08 | `70_lymphocyte_counts_core.csv` |
| per-component verdicts for the four statements | 11 | `97_verdicts_by_component.csv` |
| BALT B cell counts | 11 | `96_balt_b_cell_counts.csv` |
| batch control marker intensities | 11 | `92_batch_control_intensities.csv` |
| S2 against S1 decomposition | 11 | `93_s2_against_s1_decomposition.csv` |
| every publication figure | 09 | `80_headline_numbers.csv`, `85_arm_tests.csv` |

Two known duplications to watch. Script 09 refits models in
`86_mixed_models_fitted_here.csv` rather than reading 08's `72_models_all.csv`,
so it is a second place the same coefficients are computed and it inherits any
defect in the fit path. Script 07 rebuilds confidence intervals from
`std_err` at lines 1883 to 1884 instead of carrying `ci_low` and `ci_high`
through from the fit, which is the second place an interval is computed. Both
are the pattern that put two different IDO1 cutoff windows in scripts 03 and 09
for a year. Consolidating them is worth doing before submission.

---

# Part 6. Planned scripts and the reserved numbering

Three pending scripts were each proposed as a different number in different
notes, and two of those collide with `AKOYA_11_Claim_Checks.py`. Reserved block,
to be confirmed:

| # | name | purpose | gated on |
|---|---|---|---|
| 11b | `AKOYA_11b_Marker_Variability.py` | panel-wide marker fold survey across all 67 markers, holding the cell set fixed to CD68-lineage macrophages. Replaces the five-control batch check. Read-only | nothing, run any time |
| 12 | `AKOYA_12_Verdict_Figures.py` | figure script gated on the script 11 verdict table, so no figure can be produced for a component that failed | 11 |
| 13 | `AKOYA_13_BCell_Neighborhood_Markers.py` | what a core-proximal B cell sits next to, and what markers it carries against B cells outside foci on the same section. Merges the planned Neighborhoods figure with the marker contrast, since both ask the same question | 11b, and the count diagnostic |

The earlier notes proposed `AKOYA_11_Neighborhoods.py` and a separate
`AKOYA_12` figure script, and the 2 October review proposed script 12 for the
marker contrast. The block above resolves all three collisions. Nothing is
written yet, so this is cheap to change, but it should be changed here first so
the next restart reads one answer.

---

# Part 7. The constraints that travel with every result

Short list, because these come up in every conversation about this arm and they
are easier to re-read than to re-derive.

Three animals per arm, so C(6,3) = 20 arm assignments and any test treating the
animal as the unit is floored at 2/20 = 0.100 two-sided. That is a property of
the design, not of the test. Complete separation of the six animals is the
strongest evidence available and it arises by chance with probability 0.100 per
comparison, so the argument rests on effect size and coherence across
measurements that do not share a failure mode, not on any single p-value.

Acquisition run is perfectly confounded with treatment arm, both acquisitions
on 20251103, one per arm. Marker intensities are therefore compared only within
sections, and between-arm intensity contrasts are not identifiable in this
design. A difference of within-section contrasts is identifiable, because the
per-slide gain cancels, but only after background subtraction, since the
additive background does not cancel and differs between runs.

43118 carries 3,273 of 4,264 treated core lymphocytes, 77 percent, and has the
largest focus in the study on all three size measures. Every cell-weighted
coefficient in the treated arm is substantially about that one animal, which is
why the animal-level exact test leads. 43106 and 43111 contribute one structure
each, so structure-level variance is estimated almost entirely from the
untreated arm.

Mixed-model p-values and intervals are referred to n_animals − 2 degrees of
freedom, because treatment was applied to animals. statsmodels reports a Wald z
with no small-sample correction, which runs at 8 to 9 percent false positive at
this design.
