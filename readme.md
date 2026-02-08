[![R Conda](https://github.com/CLUES-Emory/CLUES-xMSannotator/actions/workflows/r-conda.yml/badge.svg?branch=main)](https://github.com/CLUES-Emory/CLUES-xMSannotator/actions/workflows/r-conda.yml)
![GitHub R package version (branch & subdirectory of monorepo)](https://img.shields.io/github/r-package/v/CLUES-Emory/CLUES-xMSannotator/main?filename=DESCRIPTION)

# CLUES.xMSannotator

An updated fork of [RECETOX/recetox-xMSannotator](https://github.com/RECETOX/recetox-xMSannotator) for automated annotation of untargeted mass spectrometry data, maintained by [CLUES-Emory](https://github.com/CLUES-Emory).

**Lineage:** [kuppal2/xMSannotator](https://github.com/kuppal2/xMSannotator) (original) &#8594; [RECETOX/recetox-xMSannotator](https://github.com/RECETOX/recetox-xMSannotator) (refactored) &#8594; **CLUES-Emory/CLUES-xMSannotator** (this repo)

For more details please have a look at the [Docs](docs/).

## Summary of CLUES-Emory Changes

### New Features
- **Custom pathways**: unified pathway scoring supports both HMDB and user-provided pathway databases via `pathway_mode` parameter
- **Permutation testing**: permutation-based significance testing with parallel processing, precomputed isotope caching, and full/streaming methods. Note: This feature is in development and should not be used.
- **Compound ID support**: string compound IDs (e.g., "HMDB0000001") flow through the entire pipeline and appear in all output files
- **Feature ID passthrough**: `feature_id_column` parameter preserves custom feature identifiers (e.g., "C0001") through all stages
- **Boosted compounds**: `boosted_compounds` parameter boosts confidence of confirmed annotations to level 4, with flexible mz/rt proximity matching
- **Isotope mass tolerance**: separate `isotope_mass_tolerance` parameter for ppm-based filtering of isotope matches
- **Stage outputs**: all intermediate results saved as tab-delimited text files (Stage1 through Stage5) for inspection
- **Adduct/isotope summaries**: console output summarizing adduct detection and isotope detection after each step
- **Abundance checks**: configurable `multimer_abundance_check` and `MplusH_abundance_ratio_check` parameters

### Bug Fixes
- Fixed confidence assignment bugs: forward-iterating row deletion in `apply_multimer_rules()`, NULL adduct weights, unreachable guards, fragile column deletion by position, hardcoded column indices, `cbind()` type coercion, inconsistent early-return confidence types
- Fixed chemical scoring: hardcoded RT tolerance now uses parameter, NA rows from empty results, duplicate rows from per-row instead of per-compound processing, isotope rows incorrectly removed by `na.omit()`
- Fixed isotope handling: isotopes now preserved through chemical scoring with 100x score boost, Stage 2 output column headers and file creation
- Fixed permutation testing: global environment pollution from `data(adduct_table)`, missing preprocessing steps, flawed null score matching, parallel function lookup failures
- Fixed `feature_id_column` validation error with non-numeric IDs, duplicate feature_id columns in stage outputs, Stage 5 output creation, `rm()` warnings in `get_confidence_stage4()`

### Code Cleanup
- Removed ~1500 lines of dead code not used by `advanced_annotation()` workflow (`get_confidence_stage2.R`, `multilevelannotationstep2.R`, `get_chemscorev1.6.71.R`, `group_by_rt_histv2.R`, `compute_confidence_levels.R`)
- Removed all `setwd()` calls from pipeline functions, replaced with `file.path()` absolute paths
- Removed unused `ISgroup` column, redundant `time.y` column, redundant `forms_valid_adduct_pair` filter, and `remove_tmp_files()` auto-cleanup

## Installation

The package can be installed from GitHub:

```r
devtools::install_github("CLUES-Emory/CLUES-xMSannotator")
```

## Reference

When using this tool, please cite the original work:

Uppal, Karan, et al. "XMSannotator: An R Package for Network-Based Annotation of High-Resolution Metabolomics Data." Analytical Chemistry, vol. 89, no. 2, Jan. 2017, pp. 1063-67, doi:10.1021/acs.analchem.6b01214.
