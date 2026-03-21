# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Development Commands

```r
# Install dependencies
remotes::install_deps(".")

# Build and install the package
R CMD INSTALL .
# Or from R:
devtools::install()

# Run all tests
devtools::test()

# Run a single test file
testthat::test_file("tests/testthat/test_advanced_annotation.R")

# R CMD check
R CMD check --as-cran .
# Or from R:
devtools::check()

# Regenerate documentation (roxygen2 → man/ and NAMESPACE)
devtools::document()

# Auto-format code (tidyverse style)
styler::style_pkg()
```

## Architecture

CLUES.xMSannotator is an R package for automated annotation of untargeted LC-MS data. It is a customized fork of xMSannotator (Karan Uppal) via RECETOX, maintained by CLUES-Emory.

### 5-Stage Annotation Pipeline

The core workflow is orchestrated by `advanced_annotation()` in `R/advanced_annotation.R`:

1. **Stage 1 — Mass Matching:** `simple_annotation()` matches observed m/z values to a compound database across all adducts. Uses C++ via Rcpp (`src/match_by_mass.cpp`) for performance. Applies golden rules to filter chemically implausible formulas.

2. **Stage 1.5 — Network Clustering:** `compute_peak_modules()` groups co-abundant peaks via WGCNA correlation network analysis. `compute_rt_modules()` sub-clusters by retention time using kernel density estimation.

3. **Stage 2 — Isotope Detection:** `compute_isotopes.R` detects M+1/M+2 isotope peaks by comparing observed peaks to theoretical isotopic envelopes (via enviPat). Validates by mass tolerance, intensity ratio, and RT agreement.

4. **Stage 3 — Chemical Scoring:** `get_chemscore()` scores annotations combining adduct evidence, peak correlation, and isotope detection (100x boost). `multilevelannotationstep3.R` optionally adds pathway enrichment (Fisher's exact test) using HMDB or custom pathway databases.

5. **Stage 4 — Confidence Assignment:** `multilevelannotationstep4()` assigns confidence 0–4 via a decision tree:
   - **4:** User-verified (`boosted_compounds` parameter)
   - **3:** Isotope evidence + primary adduct + module/RT coherent
   - **2:** 2+ distinct adducts + coherent
   - **1:** Single primary adduct (configurable via `level1_primary_adducts`)
   - **0:** Non-primary or incoherent

   `identify_isotopologues()` then labels specific isotope substitutions (e.g., 13C:1, 15N:1).

6. **Stage 5 — Redundancy Filtering:** `multilevelannotationstep5()` resolves multiple compounds matched to the same peak, keeping the highest-confidence/highest-score annotation.

Each stage writes intermediate output files (Stage1–Stage5 .txt files) for inspection.

### Key Entry Points
- `advanced_annotation()` — full 5-stage pipeline
- `simple_annotation()` — Stage 1 mass matching only

### Module Map
| File | Role |
|------|------|
| `R/advanced_annotation.R` | Pipeline orchestrator |
| `R/simple_annotation.R` | Stage 1 mass matching |
| `R/compute_peak_modules.R` | WGCNA peak clustering |
| `R/compute_rt_modules.R` | RT density-based clustering |
| `R/compute_isotopes.R` | Isotope detection |
| `R/get_chemscore.R` | Chemical score computation |
| `R/chemscore_helpers.R` | Score helper functions |
| `R/multilevelannotationstep3.R` | Pathway enrichment (Stage 3) |
| `R/multilevelannotationstep4.R` | Confidence assignment (Stage 4, largest module) |
| `R/multilevelannotationstep5.R` | Redundancy filtering (Stage 5) |
| `R/identify_isotopologues.R` | Isotopologue labeling |
| `R/utils.R` | Input validation, data formatting |
| `R/integration_utils.R` | Format conversion utilities |
| `src/match_by_mass.cpp` | C++11 Rcpp mass matching |

## Conventions

- **Code style:** Tidyverse formatting via `styler::style_pkg()`
- **Testing:** testthat edition 3 with `patrick::cases()` for parameterized tests. Test data lives in `tests/testthat/test-data/` (Parquet and RDA formats).
- **Documentation:** roxygen2 (v7.3.3) with markdown enabled. Run `devtools::document()` after changing roxygen comments.
- **C++ code:** Rcpp with C++11. After modifying `src/*.cpp`, run `Rcpp::compileAttributes()` then `devtools::document()`.
- **Data files:** `data/*.rda` contains pre-computed HMDB databases and adduct tables. `R/sysdata.rda` holds internal package data.
