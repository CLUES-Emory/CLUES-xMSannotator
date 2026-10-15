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
