![GitHub R package version (branch & subdirectory of monorepo)](https://img.shields.io/github/r-package/v/CLUES-Emory/CLUES-xMSannotator/main?filename=DESCRIPTION)

# CLUES.xMSannotator

An updated fork of [RECETOX/recetox-xMSannotator](https://github.com/RECETOX/recetox-xMSannotator) for automated annotation of untargeted mass spectrometry data, maintained by [CLUES-Emory](https://github.com/CLUES-Emory).

**Lineage:** [kuppal2/xMSannotator](https://github.com/kuppal2/xMSannotator) (original) &#8594; [RECETOX/recetox-xMSannotator](https://github.com/RECETOX/recetox-xMSannotator) (refactored) &#8594; **CLUES-Emory/CLUES-xMSannotator** (this repo)

For more details please have a look at the [Docs](docs/).

## Pipeline Overview

The annotation pipeline runs through five stages via `advanced_annotation()`:

1. **Stage 1 - Mass Matching**: Matches observed m/z values to compound databases across all specified adducts. Assigns peaks to correlation-based modules and RT clusters.
2. **Stage 2 - Isotope Detection**: Identifies isotopic peaks (M+1, M+2, etc.) based on mass differences, intensity ratios, and RT agreement.
3. **Stage 3 - Chemical Scoring**: Scores each compound annotation using adduct correlation evidence, module membership, and isotope support. Optionally integrates pathway enrichment (HMDB or custom).
4. **Stage 4 - Confidence Assignment**: Assigns confidence levels (0-3) based on adduct evidence, RT coherence, and filter adduct matching. For compounds spanning multiple peak modules, filters to the largest module group before evaluation. Identifies isotopologues (e.g., 13C, 15N substitutions). Upgrades non-filter compounds with strong evidence. Optionally boosts confirmed compounds to level 4. Outputs Stage4a (all rows) and Stage4b (coherent rows only); Stage 5 uses Stage4b.
5. **Stage 5 - Redundancy Filtering**: Curates annotations by removing redundant entries, keeping the highest-confidence annotation per feature.

All intermediate results are saved as tab-delimited text files (`Stage1_*.txt` through `Stage5_*.txt`) for inspection.

## Confidence Levels

Each compound annotation is assigned a confidence level reflecting the strength of supporting evidence. After Stage 4 and the evidence upgrade step, a **post-hoc evidence cap** enforces hard evidence requirements — confidence can only stay the same or decrease based on what evidence is actually present. A `Confidence_Level` text label column is added to all output files.

| Level | Label | Meaning |
|-------|-------|---------|
| 4 | Confirmed | User-confirmed compound (via `boosted_compounds` parameter) |
| 3 | High | Isotope evidence + adduct match + RT/module coherent |
| 2 | Medium | Multiple distinct adducts + RT/module coherent |
| 1 | Low | Single primary adduct match (configurable via `level1_primary_adducts`, default: M+H, M-H) |
| 0 | None | Insufficient evidence |

### Evidence Cap

After Stage 4 confidence assignment and the upgrade step, `cap_confidence_with_evidence()` enforces maximum justified confidence based on actual evidence. Confidence 0 and Confidence 4 (user-confirmed) are skipped.

**Module coherence**: When a compound has annotations in multiple peak modules, only the rows from the most-represented module are used for evidence evaluation.

**Cap rules:**

| Evidence | Maximum Confidence |
|----------|-------------------|
| Isotope rows + 1+ base adducts + module/RT coherent | 3 (High) |
| 2+ base adducts, no isotopes + module/RT coherent | 2 (Medium) |
| Single primary adduct (M+H or M-H by default) | 1 (Low) |
| Single non-primary adduct, or incoherent multi-row | 0 (None) |

The `level1_primary_adducts` parameter (default: `c("M+H", "M-H")`) controls which adducts qualify for Confidence 1 as a single match. This is independent of `filter_by`, allowing `filter_by = NULL` for equal Stage 4 scoring while still requiring primary ion evidence for Level 1.

### Evidence-Based Upgrades for Non-Filter Adducts

Compounds initially assigned Confidence 0 (no filter adduct match) are re-evaluated using three evidence types:

- **Isotope rows**: Detected isotopic peaks (M+1, M+2) for the compound's adduct
- **Multiple base adducts**: Two or more distinct adduct types detected (ignoring isotope suffixes)
- **RT coherence**: All adduct/isotope peaks within `time_tolerance` of each other

RT coherence is **required** for any upgrade above Confidence 0.

**Upgrade rules when `filter_by` is NULL/NA (same tier as filter-matched):**

| Evidence | Confidence |
|----------|-----------|
| Isotope rows + 2+ base adducts + RT coherent | 3 (High) |
| Isotope rows + 1 base adduct + RT coherent | 2 (Medium) |
| 2+ base adducts + RT coherent | 1 (Low) |

**Upgrade rules when `filter_by` is set (one tier lower):**

| Evidence | Confidence |
|----------|-----------|
| Isotope rows + 2+ base adducts + RT coherent | 2 (Medium) |
| Isotope rows + 1 base adduct + RT coherent | 1 (Low) |
| 2+ base adducts + RT coherent | 1 (Low) |

Note: After upgrades, the evidence cap may further reduce confidence if the compound's evidence doesn't meet the cap tier requirements.

## Summary of CLUES-Emory Changes

Bug fixes, new features, and code cleanup were completed using [Claude Code](https://claude.ai/claude-code) with Claude Opus.

### New Features
- **Custom pathways**: unified pathway scoring supports both HMDB and user-provided pathway databases via `pathway_mode` parameter
- **Compound ID support**: string compound IDs (e.g., "HMDB0000001") flow through the entire pipeline and appear in all output files
- **Feature ID passthrough**: `feature_id_column` parameter preserves custom feature identifiers (e.g., "C0001") through all stages
- **Confirmed compounds**: `boosted_compounds` parameter boosts confidence of confirmed annotations to level 4 (Confirmed), with flexible mz/rt proximity matching
- **Isotope mass tolerance**: separate `isotope_mass_tolerance` parameter for ppm-based filtering of isotope matches
- **Isotopologue identification**: post-confidence step uses `enviPat` to identify which specific isotope substitution each isotope peak corresponds to (e.g., 13C:1 vs 15N:1 for M+1 peaks), adding `isotopologue` and `isotopologue_quality` columns to output
- **Evidence-based confidence upgrades**: non-filter adducts (M+Na, M+NH4, etc.) are no longer permanently capped at Confidence 0. Compounds with isotope evidence, multiple adducts, and RT coherence are upgraded (see [Confidence Levels](#confidence-levels))
- **Module coherence filtering**: compounds spanning multiple peak modules are filtered to the largest module group before confidence evaluation, preventing stray matches from killing strong evidence. Stage 4 outputs split into Stage4a (all rows) and Stage4b (coherent only); Stage 5 uses Stage4b
- **Stage outputs**: all intermediate results saved as tab-delimited text files (Stage1 through Stage5) for inspection
- **Adduct/isotope summaries**: console output summarizing adduct detection and isotope detection after each step
- **Abundance checks**: configurable `multimer_abundance_check` and `MplusH_abundance_ratio_check` parameters

### Bug Fixes
- Replaced `rcdk` (Java/rJava dependency) with `enviPat` for isotope pattern calculations, eliminating the rJava dependency entirely (no JDK, `JAVA_HOME`, or `R CMD javareconf` required)
- Fixed confidence assignment bugs: dead code in non-filter adduct path (always-FALSE condition), score zeroing preventing confidence boosts, forward-iterating row deletion in `apply_multimer_rules()`, NULL adduct weights, unreachable guards, fragile column deletion by position, hardcoded column indices, `cbind()` type coercion, inconsistent early-return confidence types
- Fixed chemical scoring: hardcoded RT tolerance now uses parameter, NA rows from empty results, duplicate rows from per-row instead of per-compound processing, isotope rows incorrectly removed by `na.omit()`
- Fixed isotope handling: isotopes now preserved through chemical scoring with 100x score boost, Stage 2 output column headers and file creation
- Fixed crash on charged molecular formulas (e.g., `C12H14N2+2` for Paraquat) after rcdk to enviPat migration
- Fixed `feature_id_column` validation error with non-numeric IDs, duplicate feature_id columns in stage outputs, Stage 5 output creation, `rm()` warnings in `get_confidence_stage4()`

### Code Cleanup
- Removed ~1500 lines of dead code not used by `advanced_annotation()` workflow (`get_confidence_stage2.R`, `multilevelannotationstep2.R`, `get_chemscorev1.6.71.R`, `group_by_rt_histv2.R`, `compute_confidence_levels.R`)
- Removed experimental permutation-based p-value testing (~650 lines, never production-ready)
- Removed all `setwd()` calls from pipeline functions, replaced with `file.path()` absolute paths
- Removed unused `ISgroup` column, redundant `time.y` column, redundant `forms_valid_adduct_pair` filter, and `remove_tmp_files()` auto-cleanup
- Added roxygen2 documentation to 13 exported functions that lacked man pages

## Installation

The package can be installed from GitHub:

```r
devtools::install_github("CLUES-Emory/CLUES-xMSannotator")
```

## Reference

When using this tool, please cite the original work:

Uppal, Karan, et al. "XMSannotator: An R Package for Network-Based Annotation of High-Resolution Metabolomics Data." Analytical Chemistry, vol. 89, no. 2, Jan. 2017, pp. 1063-67, doi:10.1021/acs.analchem.6b01214.
