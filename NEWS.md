# CLUES.xMSannotator 1.2.0

## New features

* `advanced_annotation()` gains `feature_groups`: a data frame of `feature_id`
  and `group`, e.g. the feature groups written by Step 1 Stage 3
  (HighResolutionExposomics_Workflow v1.7.0, `*_Feature_Groups.txt`). When
  given, WGCNA peak modules and retention-time clusters are not computed: each
  group is one module and every peak has `rt_cluster` 1, so isotope candidates
  are limited by `time_tolerance` only. Needs `feature_id_column`; every peak
  must have a group. On CLU0127 (about 10,000 features x 347 injections) a run
  took about 30 s instead of 11-17 min, and precision against the confirmed
  library differed by at most 0.01. `NULL` (default) leaves the output unchanged.

## MultiDB runscript (`inst/scripts/xMSannotator_CLUES_MultiDB_Runscript.R`)

Results are not comparable with runs made with the 1.1.0 runscript.

* Reads the Step 1 Stage 3 feature groups when the file exists
  (`3-<mode>_<study>_Feature_Grouping/<mode>-<study>_Feature_Groups.txt`); their
  feature IDs then define the feature set and the `fc_threshold` argument is not
  applied. Without the file, the fold-change filter and WGCNA run as before.
* New defaults: `filter_by` M+H / M-H (was `NULL`), weight 5 for M+H / M-H
  (was 1), `mass_tolerance` and `isotope_mass_tolerance` 1.5e-6 (were 5e-6 and
  7.5e-6), `time_tolerance` 2.5 s (was 10). Chosen on CLU0127 C18pos/C18neg
  against the CLUES confirmed library; the evidence is in the Section 5 and 6
  comments.
* Selects the `*_Untargeted_RAW_Feature_Table_*.txt` file in the Unadjusted
  folder. Step 1 v1.3.2+ also writes the injection QC flags there, so the
  previous `.txt` pattern passed two files to `fread()`, which failed.
* Feature metadata columns in the exported table updated to Step 1 v1.6.0 names
  (earlier names still used when present); adds `rtmin`, `rtmax` and
  `feature_group`.

## Example runscript and SLURM wrapper

* `xMSannotator_CLUES_Runscript_Example.R`: new `feature_groups_path` (Step 1
  Stage 3 groups file; `NULL` falls back to the `fc_threshold` filter and
  WGCNA), and the MultiDB runscript's settings: `filter_by` M+H / M-H, weight 5,
  1.5 ppm mass and isotope tolerance, 2.5 s, Step 1 v1.6.0 metadata columns.
* `run_CLUES_xMSannotator.sh`: comments only (groups file, `fc_threshold` use,
  memory per task).
* Documentation: `feature_groups` in the Workflow reference, Input Formats,
  pre-processing guide and Stage 5 reference; readme links `NEWS.md`.

# CLUES.xMSannotator 1.1.0

## Results are not comparable with 1.0.0

Several changes alter Stage 4 output on the default settings. Confidence levels
produced by 1.1.0 should not be compared with those from 1.0.0; a difference
between versions will otherwise read as a change in data quality.

* `filter_by` now defaults to `NULL` (was `c("M-H", "M+H")`). Whether
  `filter_by` is set switches the confidence rules: with a primary ion named, a
  compound lacking it caps at Confidence 1 and isotope evidence promotes only to
  2; with `NULL`, isotope evidence can reach Confidence 3. On one test dataset
  this alone moved Level 3 from 40 to 136 compounds. Set
  `filter_by = c("M-H", "M+H")` to restore the previous behaviour.
* The coherence step keeps the best-supported module (Confidence, then row
  count, then a primary ion, then summed intensity) instead of the module with
  most rows. The old rule tied in most multi-module groups and then fell back to
  the alphabetical order of module labels, which systematically discarded the
  M+H row.
* Confidence is capped per group coherent in both module and retention time,
  not per compound, so one compound's rows no longer inherit a verdict from a
  peak elsewhere in the chromatogram.
* The isotope intensity window is centred on the theoretical abundance
  (symmetric `1 +/- intensity_deviation_tolerance`); it previously favoured
  isotopes reading stronger than theory.
* The 100x isotope boost applies only when the isotopes belong to the module
  that won the scoring.
* Retention-time clusters are bounded by the time tolerance from their anchor,
  so rows at the far end of a long chain of close peaks are no longer dropped
  unscored.
* `mass_defect_tolerance` is now in Da: how far the charge-corrected isotope
  spacing may sit from a whole isotope step. The previous comparison of binned
  indices discarded the M+1 peak for about a third of compounds. Values >= 0.5
  disable the check, and `compute_isotopes()` now warns when given one, since
  that is the likely result of a script written for the old units (the shipped
  runscripts used 5; they now use 0.02). The `compute_isotopes()` default
  changes from 0, which now rejects every isotope, to 0.1, matching
  `advanced_annotation()`.

## New features

* **Permanently charged compounds.** The compound table accepts an optional
  integer `charge` column. Compounds with a non-zero charge are matched only as
  their bare ion (`[M]+`, `[M]+2`, `[M]-`, ...) at `monoisotopic_mass / |z|`,
  and isotope detection uses the ion's charge. When `charge` is supplied,
  neutral compounds are no longer matched at their neutral mass ("M").
  **Tables without a `charge` column behave as before**: "M" matches are kept,
  because databases such as HMDB list permanently charged compounds (choline,
  for example) at their ion mass under a sign-less formula, and "M" is the only
  adduct that finds them. A message notes when this applies.
* `advanced_annotation(charge_type = "pos" | "neg")` drops compounds of the
  opposite charge sign; `filter_compounds_by_polarity()` applies the same rule.
* `multi_rt_clusters` (default `FALSE`): annotate every retention-time cluster
  a compound matches, each with its own confidence, instead of choosing one.
  Validated against authentic-standard confirmations on CLU0009: recall of the
  confirmed feature rose from 30 to 55 of 72. Output roughly doubles, almost
  entirely as Confidence 0. Untested in combination with
  `redundancy_filtering = TRUE`.
* `n_rt_clusters` is reported in Stage 4a/4b: how many candidate retention times
  each compound was chosen from. Confidence 3 rates rise with this count, so it
  indicates how much a given annotation can be trusted.
* `isotope_ratio_dev` (observed vs theoretical isotope intensity) is retained in
  Stage 2.
* `clustering_seed` (default 12345) pins the peak-module clustering seed rather
  than relying on a WGCNA default.

## Bug fixes

* `multilevelannotationstep4` no longer duplicates the "M" adduct row when the
  adduct table already contains it.
* `identify_isotopologues` handles charge-suffixed formulas, and uses the
  electron-corrected ion mass for bare ions.
* `boost_confidence_of_IDs()` matches on the base compound id when ids are
  qualified by retention-time cluster.
* The example runscript reads delimited-text databases and joins feature-table
  QC columns using `feature_id_column`.

## Known issues

* The end-to-end `advanced_annotation` tests, and several `get_chemscore` and
  integration-utility tests, error on this release as they did on 1.0.0 (13
  tests in total). They have not been fixed here.
