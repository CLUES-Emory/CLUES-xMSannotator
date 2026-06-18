# CLUES.xMSannotator Workflow Reference

## Overview

CLUES.xMSannotator is a network-based annotation pipeline for LC-MS metabolomics data. Starting from a peak table (m/z, RT, intensities) and a compound database, it assigns putative identities to detected features through mass matching, co-abundance network analysis, isotope pattern verification, chemical scoring, pathway enrichment, confidence assignment, and redundancy filtering.

This pipeline is a customized fork of [RECETOX/recetox-xMSannotator](https://github.com/RECETOX/recetox-xMSannotator), 
which is itself a refactored version of the original 
[xMSannotator](https://doi.org/10.1021/acs.analchem.6b01214) R package 
by Karan Uppal (Emory University). Key modifications in this fork include:

- Evidence ceiling (Approach C) for post-hoc confidence adjustment
- WGCNA module and RT coherence filtering
- Support for `boosted_compounds` (user-verified standards → Confidence 4)
- Multi-evidence chemical scoring with isotope, adduct, and correlation weighting
- Intermediate stage output files (`Stage1_*.txt` through `Stage4b_*.txt`) for inspection. `Stage4b_confidence_levels.txt` is the canonical final output; `Stage5_curated_results.txt` is written only when `redundancy_filtering = TRUE` (opt-in).

The entry point is `advanced_annotation()` in `R/advanced_annotation.R`.

**Citation:** Uppal K, Walker DI, Jones DP. xMSannotator: an R package for network-based annotation of high-resolution metabolomics data. *Analytical Chemistry*. 2017;89(2):1063-1067.

---

## Quick Start

A minimal annotation run requires a peak table (m/z, RT, and sample intensity columns) and a compound database (compound_id, name, molecular_formula, monoisotopic_mass). The following example shows a positive-mode run with the most important parameters:

```r
library(CLUES.xMSannotator)

# Load your data
peak_table     <- data.table::fread("feature_table.txt", data.table = FALSE)
compound_table <- readxl::read_xlsx("compound_database.xlsx")

# Filter adduct table to positive-mode adducts
adduct_table <- get_default_adduct_table()
adduct_table <- adduct_table[adduct_table$adduct %in%
  c("M+H", "M+2H", "M+NH4", "M+Na", "M+ACN+H", "2M+H", "M+H-H2O"), ]

# Build adduct weights (all selected adducts, equal weight)
adduct_weights <- data.frame(adduct = adduct_table$adduct, weight = rep(1, nrow(adduct_table)))

annotation <- advanced_annotation(
  peak_table        = peak_table,            # features × samples intensity matrix
  compound_table    = compound_table,        # database of candidate compounds
  adduct_table      = adduct_table,          # adduct definitions for positive mode
  adduct_weights    = adduct_weights,        # which adducts to prioritize in scoring
  mass_tolerance    = 5e-6,                  # 5 ppm (fractional form)
  time_tolerance    = 10,                    # RT grouping window (seconds)
  filter_by         = "M+H",                # expected primary adduct
  pathway_mode      = "skip",               # "HMDB", "custom", or "skip"
  outloc            = "xMSannotator_output"  # output directory
)
```

After a successful run, the output directory will contain tab-delimited text files for each pipeline stage: `Stage1_mass_matched.txt`, `Stage1_peak_clusters.txt`, `Stage2_isotope_detection.txt`, `Stage3_chemical_scores.txt`, a Stage 3 pathway file, `Stage4a_confidence_levels.txt`, and `Stage4b_confidence_levels.txt`. `Stage4b_confidence_levels.txt` is the primary output — each row is a putative annotation with a confidence level (0–4) backed by the full evidence basis used to compute it.

`Stage5_curated_results.txt` is produced only when `redundancy_filtering = TRUE` (opt-in). Stage 5 resolves multi-match features to a single compound per `(mz, time)` by deleting losing-compound rows; this is convenient for a one-row-per-peak view but discards the secondary adduct and isotope rows that contributed to other compounds' confidence values. For confidence-faithful downstream analysis, prefer Stage 4b.

---

## Pipeline Data Flow

The pipeline runs seven processing steps in sequence. Each step consumes upstream results, adds new evidence, and writes an intermediate output file:

1. **Stage 1 — Mass Matching.** Finds all compound–adduct combinations whose expected m/z falls within tolerance of each detected feature. Consumes the peak table, compound database, and adduct table. Produces `Stage1_mass_matched.txt` and `Stage1_peak_clusters.txt`.

2. **Stage 1.5 — Peak Correlation & Module Detection.** Groups features into co-abundance modules using WGCNA, then sub-clusters each module by retention time. Consumes the peak intensity matrix. Module and RT cluster assignments are written to `Stage1_peak_clusters.txt`.

3. **Stage 2 — Isotope Detection.** Searches for isotopic peaks (M+1, M+2, ...) matching each annotation's theoretical isotope envelope. Consumes Stage 1 annotations joined with module/RT assignments. Produces `Stage2_isotope_detection.txt`.

4. **Stage 3 — Chemical Scoring.** Computes a multi-evidence score for each annotation based on adduct count, correlation strength, isotope confirmation, and RT coherence. Consumes Stage 2 output and the global correlation matrix. Produces `Stage3_chemical_scores.txt`.

5. **Stage 3b — Pathway Enrichment.** Tests whether annotated compounds are enriched in known biological pathways using Fisher's exact test, and boosts scores for pathway-supported annotations. Consumes Stage 3 scores and a pathway database. Produces `Stage3_HMDB_pathways.txt` (or `Stage3_custom_pathways.txt` or `Stage3_pathway_skipped.txt`).

6. **Stage 4 — Confidence Assignment.** Translates continuous scores into discrete confidence levels (0–4) based on evidence thresholds, then filters to coherent annotations. Consumes Stage 3b output. Produces `Stage4a_confidence_levels.txt` (all rows) and `Stage4b_confidence_levels.txt` (coherent rows only).

7. **Stage 5 — Redundancy Filtering (opt-in).** Resolves features matched to multiple compounds by selecting one row per `(mz, time)` and deleting the others. Consumes Stage 4b output. Produces `Stage5_curated_results.txt`. Off by default because the deletion is destructive on multi-match `(mz, time)` evidence — it removes secondary adduct and isotope rows belonging to losing compounds, which can leave those compounds' remaining rows displaying confidence values computed from evidence that no longer appears in the file. Enable via `redundancy_filtering = TRUE` when a one-row-per-peak summary is more useful than evidence integrity.

Each stage adds an independent layer of evidence to the annotation. Stage 1 provides the initial mass-match candidates. Stage 1.5 adds co-abundance structure (which features behave similarly across samples). Stage 2 adds isotope-pattern confirmation. Stage 3 synthesizes these signals into a single score. Stage 3b adds biological context from pathway databases. Stages 4 and 5 then use the accumulated evidence to assign confidence and resolve ambiguity. A compound that survives all stages with high confidence has been validated by mass accuracy, co-elution, isotope pattern, correlation network membership, and (optionally) pathway co-occurrence — multiple independent lines of evidence converging on the same identification.

---

## Stage 1: Mass Matching

**Files:** `R/simple_annotation.R`, `src/match_by_mass.cpp`

> **Key questions this stage answers:**
> Which compounds in the database could produce this observed m/z? Under what adduct form, charge state, and multimer configuration?

### What it does

Compares every detected m/z against a compound database to find putative adduct matches within a mass tolerance window. This is a brute-force search accelerated by a C++ inner loop.

### Adduct equation

The C++ engine (`match_by_mass.cpp`, line 39) computes the expected m/z for each compound–adduct combination:

```
expected_mz = (factor × M + adduct_mass) / charge
```

- **M** = compound monoisotopic mass
- **factor** = number of molecules in the adduct species (1 for `M+H`, 2 for `2M+H`, 3 for `3M+H`)
- **adduct_mass** = mass contribution of the adduct ion (e.g., +1.007276 for protonation)
- **charge** = charge state (1, 2, or 3)

### Mass tolerance

Applied as relative (ppm) tolerance via `almost_equal()` (`match_by_mass.cpp`):

```
|observed_mz − expected_mz| ≤ max(|observed_mz|, |expected_mz|) × tolerance
```

Default: `mass_tolerance = 5e-6` (5 ppm in fractional form).

### Golden Rules filters

Before output, every match is screened by `check_golden_rules()` (lines 13–38 of `simple_annotation.R`). These are elemental ratio heuristics that reject chemically implausible molecular formulas:

**Elemental ratio checks** (all must pass, requires C > 0):

| Ratio | Allowed range |
|-------|---------------|
| H/C   | 0.1–6        |
| F/C   | ≤ 6          |
| N/C   | ≤ 4          |
| O/C   | ≤ 3          |
| P/C   | ≤ 2          |
| S/C   | ≤ 3          |

**Combined element restrictions** (upper bounds on element counts when multiple heteroatoms co-occur):

| Elements | Must satisfy (if all > threshold) |
|----------|----------------------------------|
| N,O,P,S  | N<10 AND O<20 AND P<4 AND S<3   |
| N,O,P    | N<11 AND O<22 AND P<6            |
| O,P,S    | O<14 AND P<3 AND S<3             |
| P,S,N    | N<4 AND P<3 AND S<3              |
| N,O,S    | N<19 AND O<14 AND S<8            |

**Water-loss adduct rule** (`forms_valid_adduct_pair`): Adducts involving water loss (M+H−H₂O, M+H−2H₂O, M−H₂O−H) are only allowed if the molecular formula contains oxygen.

### Adduct table

The default `sample_adduct_table` contains 110 adducts spanning both positive and negative ionization modes, including:
- Single-charge primary ions: M+H, M−H, M+Na, M+K, M+NH₄, M+Cl, M+Br
- Multi-charge ions: M+2H (z=2), M+3H (z=3), M−2H (z=2), M−3H (z=3)
- Multimers: 2M+H, 3M+H, 2M−H, 3M−H
- Adducts with common mobile-phase modifiers: formic acid, acetic acid, acetonitrile, TFA, DMSO, isopropanol

### Output

`Stage1_mass_matched.txt` — all candidate annotations with columns: peak, mz, rt, compound, adduct, expected_mass, molecular_formula, name, plus inherited compound metadata.

---

## Stage 1.5: Peak Correlation & Module Detection

**Files:** `R/compute_peak_modules.R`, `R/compute_rt_modules.R`

> **Key questions this stage answers:**
> Which features co-vary in intensity across samples and therefore likely originate from the same compound or biological process? Within each co-abundance module, which features co-elute?

### Theory

Metabolites entering the mass spectrometer generate multiple related ions: adducts, in-source fragments, and isotopes. These co-eluting, co-abundant species will show correlated intensity patterns across samples. By building a co-abundance network from the peak intensity matrix, the pipeline groups peaks into **modules** — clusters of features likely originating from the same or related metabolites. This biological signal is later used to boost annotation confidence.

### Step A: Pearson correlation matrix

**Function:** `compute_peak_correlations(peak_intensity_matrix, correlation_method = "p")`

Computes a pairwise Pearson correlation matrix across all peaks using `WGCNA::cor()`. Input is a samples × peaks intensity matrix. Output is a symmetric peaks × peaks correlation matrix.

### Step B: WGCNA module detection

**Function:** `compute_peak_modules()`

1. **Soft-threshold power selection** — `WGCNA::pickSoftThreshold()` finds the power that best approximates a scale-free network topology. Fallback: power = 6 if estimation fails.
2. **Hierarchical clustering** — `flashClust::flashClust()` clusters peaks by correlation distance (1 − correlation).
3. **Dynamic tree cutting** — `dynamicTreeCut::cutreeDynamic()` identifies modules with adaptive branch cutting (parameters: `deep_split = 2`, `min_cluster_size = 10`).
4. **Blockwise module detection** — `WGCNA::blockwiseModules()` refines modules using the soft-threshold power, then merges modules with correlation above the threshold (`mergeCutHeight = 1 − correlation_threshold`, i.e., 0.3 at default).
5. **Fallback** — if `blockwiseModules()` fails, the pipeline computes a Topological Overlap Matrix (TOM) and uses `mergeCloseModules()`.

**Key defaults:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `correlation_threshold` | 0.7 | Module merge threshold (merge if r > 0.7) |
| `deep_split` | 2 | Tree cutting depth (0=conservative, 4=aggressive) |
| `min_cluster_size` | 10 | Minimum peaks per module |
| `network_type` | "unsigned" | Considers both positive and negative correlations |

### Step C: RT clustering within modules

**Function:** `compute_rt_modules(peak_table, peak_width = 1)`

Within each module, peaks are further sub-grouped by retention time using kernel density estimation:

1. **KDE** — Gaussian kernel density with bandwidth = `peak_width` (default: 1 second), minimum 2048 evaluation points.
2. **Peak detection** — `pastecs::turnpoints()` finds local maxima in the density curve (inflection points above machine epsilon).
3. **Nearest-neighbor assignment** — `RANN::nn2()` assigns each peak to its nearest RT cluster center.

**Output:** `Stage1_peak_clusters.txt` — each peak tagged with `module`, `rt_cluster`, and the combined identifier `Module_RTclust` (e.g., "110_3").

---

## Stage 2: Isotope Detection

**File:** `R/compute_isotopes.R`

> **Key questions this stage answers:**
> Does the observed isotopic envelope (M+1, M+2, ...) match the theoretical pattern predicted by this compound's molecular formula? Are the isotope peaks present at the expected masses and intensities?

### Theory

Each molecular formula produces a characteristic isotopic envelope due to the natural abundance of heavier isotopes (¹³C, ²H, ¹⁸O, etc.). For example, a compound with 10 carbons will show an M+1 peak at ~11% of the monoisotopic intensity (from ¹³C). By matching the observed isotopic pattern against the theoretical pattern, the pipeline gains additional evidence for or against an annotation.

### Algorithm

For each annotated peak (the putative monoisotopic ion):

1. **Compute theoretical pattern** — `enviPat::isopattern()` generates the expected isotope envelope from the molecular formula. Peaks below 0.1% relative abundance (`minAbund = 0.001`) are discarded.

2. **Filter candidates** — Nearby peaks are screened by:
   - Same `rt_cluster` (must be co-eluting)
   - RT within `rt_tolerance` (default: 1 second in `compute_isotopes`, but receives `time_tolerance` = 10 from `advanced_annotation`)
   - Same `mass_defect` bin (tolerance: 0, exact match by default)
   - Not the same peak

3. **Mass matching** — Each candidate's m/z is compared to the expected isotope m/z:
   ```
   expected_isotope_mz = monoisotopic_expected_mass + exact_mass_diff
   mass_error_ppm = |observed_mz − expected_isotope_mz| / expected_isotope_mz × 10⁶
   ```
   Filtered by `isotope_mass_tolerance_ppm` (defaults to `mass_tolerance × 10⁶`, i.e., 5 ppm).

4. **Intensity matching** — The observed intensity ratio (isotope/monoisotopic) is compared to the theoretical relative abundance:
   ```
   |observed_ratio − theoretical_abundance| ≤ observed_ratio × intensity_deviation_tolerance
   ```
   Default tolerance: 0.1 (10% relative deviation).

### Key defaults

| Parameter | Default | Description |
|-----------|---------|-------------|
| `intensity_deviation_tolerance` | 0.1 | 10% allowed deviation in intensity ratio |
| `mass_defect_tolerance` | 0.1 | Mass defect window (Da) |
| `rt_tolerance` | receives `time_tolerance` | RT window (seconds) |
| `isotope_mass_tolerance_ppm` | `mass_tolerance × 10⁶` | PPM filter for isotope mass matching |
| `minAbund` (enviPat) | 0.001 | 0.1% relative abundance cutoff |

### Output

`Stage2_isotope_detection.txt` — the annotation table expanded with isotope rows. A new column `mass_number_difference` indicates: 0 = monoisotopic, 1/2/3... = M+1/M+2/M+3 isotope.

---

## Stage 3: Chemical Scoring

**Files:** `R/get_chemscore.R`, `R/chemscore_helpers.R`

> **Key questions this stage answers:**
> How much independent evidence supports this annotation? Is it a lone mass match, or do multiple adducts, isotope patterns, and intensity correlations converge on the same compound?

### Theory

Not all mass matches are equally credible. A compound with multiple corroborating adducts in the same RT window, strong intensity correlations, and matching isotope patterns is far more likely to be a true identification than a lone hit. Chemical scoring quantifies this multi-evidence support into a single numeric score.

### Hub m/z identification

**Function:** `compute_hub_mz()` (line 346, `chemscore_helpers.R`)

The "hub" is the most informative feature for a given compound — the anchor around which other evidence is evaluated. Selection follows a hierarchical criterion:

1. **Primary:** Features with adducts in the `filter.by` list AND positive correlation with other features
2. **Secondary:** Features with weighted adducts (in `adduct_weights`) AND no isotope annotations
3. **Tertiary:** Any feature with correlation > 0 and no isotope annotations

Among candidates, the hub is the feature with the **highest mean intensity**. Ties are broken by the feature with the most time neighbors within `max_diff_rt`.

### Scoring formula

**Function:** `calc_base_score()` (line 363, `chemscore_helpers.R`)

```
chemical_score = n_unique_adducts × n_weighted_adducts × max_correlation × 10^max_adduct_weight
```

- **n_unique_adducts** — count of distinct adduct types (e.g., M+H, M+Na = 2)
- **n_weighted_adducts** — count of adducts matching the `adduct_weights` table
- **max_correlation** — highest Pearson correlation between the hub and any other feature for this compound (from the upper triangle of the correlation submatrix)
- **10^max_adduct_weight** — exponential encoding of the highest adduct weight (default weights typically 0–5; e.g., weight=2 → ×100)

### Isotope evidence boost

After monoisotopic scoring, if isotopic peaks are detected for the compound, the score is multiplied by 100:

```r
if (has_isotopes) chemical_score <- chemical_score * 100
```

This 100× multiplier strongly rewards annotations with isotope confirmation.

### RT scaling

**Function:** `apply_rt_scaling()` (line 785, `chemscore_helpers.R`)

The score is penalized based on the RT spread of the compound's features:

```
chemical_score = chemical_score × 1 / ((diff_rt × 0.1) + 1)^k
```

- `diff_rt` = max(RT) − min(RT) across all features for this compound
- `k = 1` if diff_rt ≤ `max_diff_rt` (gentle penalty)
- `k = 10` if diff_rt > `max_diff_rt` (severe penalty — effectively eliminates scattered annotations)

**Example at max_diff_rt = 10s:**
- diff_rt = 5s → scaling ≈ 0.67 (mild penalty)
- diff_rt = 15s → scaling ≈ 0.00003 (effectively zero)

### Confidence-level boost

After all other adjustments, scores above the minimum threshold receive an exponential boost based on the preliminary confidence level:

```
chemical_score = chemical_score × conf_level^conf_level
```

| Confidence level | Boost factor |
|:----------------:|:------------:|
| 0 (None)         | 0 (zeroed)  |
| 1 (Low)          | 1× (no change) |
| 2 (Medium)       | 4×          |
| 3 (High)         | 27×         |

### Minimum score threshold

**Function:** `get_min_chemscore()` (line 781, `chemscore_helpers.R`)

```
min_score = 100 × 2 × corthresh × (1 / ((max_diff_rt × 0.1) + 1)^3)
```

At defaults (corthresh=0.7, max_diff_rt=10): min_score ≈ 17.5.

### Worked example

Consider a hypothetical compound detected with the following evidence:
- 2 unique adducts: M+H and M+Na
- 1 of those adducts is weighted (max weight = 2)
- Highest Pearson correlation between hub and another feature: 0.85
- Isotope peaks detected (M+1 confirmed)
- RT spread across features: diff_rt = 5 seconds (max_diff_rt = 10)
- Preliminary confidence level: 2 (Medium)

The score is computed in five steps:

| Step | Calculation | Result |
|------|-------------|--------|
| **1. Base score** | `2 × 1 × 0.85 × 10^2` | 170 |
| **2. Isotope boost** | `170 × 100` | 17,000 |
| **3. RT scaling** (k=1, since 5 ≤ 10) | `17,000 × 1/((5 × 0.1) + 1)^1` = `17,000 × 1/1.5` | 11,333 |
| **4. Min score check** | `100 × 2 × 0.7 × 1/((10 × 0.1) + 1)^3` = `140 × 1/8` = 17.5 | 11,333 > 17.5: **passes** |
| **5. Confidence boost** | `11,333 × 2^2` | **45,333** |

The final chemical score of ~45,333 places this annotation firmly in the "strong multi-evidence" range (>10,000). Without isotope confirmation (skip step 2), the score would be ~453 — still meaningful but two orders of magnitude lower, reflecting the reduced certainty.

### Output

`Stage3_chemical_scores.txt` — annotations with `cur_chem_score` column representing the multi-evidence score.

---

## Stage 3b: Pathway Enrichment

**Files:** `R/multilevelannotationstep3.R`, `R/compute_pathways.R`

> **Key questions this stage answers:**
> Are the annotated compounds over-represented in known biological pathways? Does the co-location of pathway members within the same WGCNA module suggest a real biological signal rather than random matches?

### Theory

If multiple annotated metabolites map to the same biological pathway, this co-occurrence is unlikely to be coincidental and supports the annotations' validity. Pathway enrichment tests whether the overlap between annotated compounds and known pathway members exceeds what would be expected by chance.

### Fisher's exact test — two levels

**Level 1: Global pathway enrichment** (line 72–77)

For each pathway, constructs a 2×2 contingency table:

|                    | In pathway | Not in pathway |
|--------------------|:----------:|:--------------:|
| **Significant hit** |     a      |       b        |
| **Not significant** |     c      |       d        |

"Significant" = annotation has score ≥ `scorethresh` (default: 0.1) AND adduct is in `adduct_weights`.

Fisher's exact test is applied. **Threshold: p ≤ 0.05.**

**Level 2: Module-pathway overlap** (line 115–124)

For pathways passing Level 1, tests whether pathway compounds are co-located within the same WGCNA/RT module more than expected by chance. This is a second 2×2 contingency table comparing the dominant module against all other modules.

**Threshold: p ≤ 0.2** (more permissive, since modules are smaller subsets).

### Score boost conditions

All of the following must be met for a compound to receive a pathway boost:

1. Pathway passes Level 1 Fisher's test (p ≤ 0.05)
2. Module passes Level 2 Fisher's test (p ≤ 0.2)
3. ≥ 3 pathway compounds in the module
4. Compound's chemical score ≥ 0.1
5. Compound has a valid adduct (in `adduct_weights`)

**Boost amount:** `score += n_pathway_compounds` (additive, per qualifying pathway).

### Pathway databases

| Mode | Source | Selection |
|------|--------|-----------|
| `"HMDB"` (default) | Internal `hmdbCompMZ` dataset (HMDB v3.5) | SMPDB pathway IDs parsed from column 14 |
| `"custom"` | User-provided `pathway_data` data frame | Must have `compound` and `pathway` columns |
| `"skip"` | None | Pathway step is bypassed |

### Output

- `Stage3_HMDB_pathways.txt` (HMDB mode)
- `Stage3_custom_pathways.txt` (custom mode)
- `Stage3_pathway_skipped.txt` (skip mode)

The `cur_chem_score` column is renamed to `score` at this stage.

---

## Stage 4: Confidence Assignment

**File:** `R/multilevelannotationstep4.R`

> **Key questions this stage answers:**
> How confident should we be in this annotation? Is the evidence strong enough for the annotation to be reported, or is it a speculative mass match that should be flagged as low-confidence?

### Theory

Chemical scores are continuous values that resist intuitive interpretation. Confidence levels translate multi-evidence support into a discrete 0–4 scale that communicates the strength of annotation evidence at a glance. Higher confidence requires convergent evidence from multiple independent sources (adducts, isotopes, module co-membership, RT coherence).

### Confidence levels

| Level | Label | Interpretation |
|:-----:|-------|----------------|
| 0 | None | Insufficient evidence; mass match only |
| 1 | Low | Single primary adduct match (M+H or M−H) |
| 2 | Medium | Multiple distinct adducts or isotope evidence with moderate coherence |
| 3 | High | Multiple adducts + isotope evidence + module coherence + RT coherence |
| 4 | Confirmed | User-verified compound (via `boosted_compounds` parameter) |

### Assignment algorithm

**Function:** `get_confidence_stage4()` (lines 442–542)

The algorithm evaluates each compound's annotation rows through a decision tree. The following pseudocode traces the major decision points. To determine the confidence for any annotation, start at the top and follow the first matching branch:

```
IF score < 0.1
  → CONFIDENCE 0 (None)

IF adduct not recognized in adduct_table
  → CONFIDENCE 0 (None)

IF only one unique adduct type:
  IF score < 10 AND no filter_by match
    → CONFIDENCE 0 (None)
  IF score < 10 AND has filter_by match
    → CONFIDENCE 1 (Low)
  IF nrow == 2 (single adduct but two rows, suggesting isotope)
    → CONFIDENCE 2 (Medium)
  IF nrow > 2 AND no filter_by match
    → CONFIDENCE 0 (None)

IF RT spread > max.rt.diff:
  Re-cluster rows by RT → keep the best cluster
  (continue evaluation with filtered rows)

IF features span multiple WGCNA modules:
  IF score < 10
    → CONFIDENCE 0 (None)
  IF score > 10 AND has filter_by match
    → CONFIDENCE 2 (Medium)

IF only one row remains:
  IF adduct not in adduct_weights
    → CONFIDENCE 1 (Low)
  IF adduct in adduct_weights AND score > 10 AND has filter_by match
    → CONFIDENCE 2 (Medium)
  ELSE
    → CONFIDENCE 0 (None)

IF multimer (2M, 3M) is more abundant than monomer
  → downgrade confidence

IF charge > 1
  → CONFIDENCE 1 (Low)

IF carbon count < 1 in molecular formula
  → CONFIDENCE 0 (None)

Otherwise: assign based on adduct count and score thresholds
```

Two score thresholds drive most decisions: **0.1** (minimum to be considered at all) and **10** (threshold for medium-confidence promotion). The `filter_by` parameter acts as a gatekeeper — annotations matching the expected primary adduct (e.g., M+H in positive mode) can reach Confidence 1 even with low scores, while non-filter adducts require stronger evidence.

### Confidence upgrade (post-hoc)

**Function:** `upgrade_confidence_with_evidence()` (lines 869–932)

Re-evaluates compounds below Confidence 3 using isotope and adduct evidence. Can only **upgrade**, never downgrade.

| Evidence pattern | With filter_by active | With filter_by = NULL |
|------------------|-----------------------|-----------------------|
| Isotope rows + ≥1 base row + ≥2 base adducts | → Conf 2 | → Conf 3 |
| Isotope rows + ≥1 base row (single adduct) | → Conf 1 | → Conf 2 |
| ≥2 base adducts (no isotopes) | → Conf 1 | → Conf 2 |
| No evidence | No upgrade | No upgrade |

Requirements: actual base rows must exist (not just isotope inference), and module + RT coherence must pass.

### Confidence cap (post-hoc)

**Function:** `cap_confidence_with_evidence()` (lines 959–1007)

Enforces hard evidence ceilings. Can both **raise** (rescue Level 0) and **lower** (prevent over-confidence).

| Max justified confidence | Required evidence |
|:------------------------:|-------------------|
| 3 (High) | Isotope rows + ≥1 base row + ≥1 adduct + module + RT coherence |
| 2 (Medium) | ≥2 base adducts + module + RT coherence |
| 1 (Low) | Single primary adduct (M+H or M−H) from an actual base row |
| 0 (None) | Insufficient evidence |

**Exceptions:** Confidence 4 (user-confirmed) is never modified.

### Coherence filtering

**Function:** `enforce_compound_coherence()` (lines 114–122)

Two requirements for a compound's rows to be "coherent":

1. **Module coherence** — all rows must share the same WGCNA module
2. **RT coherence** — RT range (max − min) must be ≤ `max.rt.diff` seconds

For compounds split across modules, only the rows from the **largest module group** are kept. This prevents minority stray rows from degrading confidence.

### Stage 4a vs 4b

| | Stage 4a | Stage 4b |
|-|----------|----------|
| **Content** | All rows with confidence assigned | Coherent rows only (largest module group per compound) |
| **Purpose** | Audit trail — inspect all annotations including dropped rows | Production data — feeds Stage 5 |
| **File** | `Stage4a_confidence_levels.txt` | `Stage4b_confidence_levels.txt` |

### Confidence boosting (user-confirmed compounds)

**Function:** `boost_confidence_of_IDs()` (lines 664–703)

If `boosted_compounds` is provided, matching annotations are elevated to Confidence 4 with a 100× score multiplier. Matching modes:
- **ID-only:** exact compound ID match
- **Proximity-based:** compound ID + m/z tolerance (`boost_mass_tolerance`) + RT tolerance (`boost_time_tolerance`)

---

## Stage 5: Redundancy Filtering

**File:** `R/multilevelannotationstep5.R`

> **Key questions this stage answers:**
> When a single m/z feature matches multiple compounds, which annotation is best supported? After resolving conflicts, is each remaining annotation a unique match or a winner among competitors?

### Theory

A single m/z feature may match multiple compounds in the database. Stage 5 resolves these conflicts by selecting the best-supported annotation for each feature, using adduct priority and chemical score as the arbitration criteria.

### Multi-match resolution algorithm

**Function:** `multilevelannotationstep5()` (lines 53–85)

1. Identify m/z values matching multiple compounds (`get_features(mz, "multiple")`)
2. For each multi-match m/z:
   - **Boost weighted adducts:** annotations with adducts in `adduct_weights` receive a 100× score multiplier
   - **Select winner:** annotation with the highest (boosted) score
   - **Remove losers at this m/z:** other compounds' rows at this specific m/z are dropped
   - **Recalculate loser scores:** `score = (n_competing − 1) / n_competing × winning_score`

### Adduct weighting

**Function:** `increase_multimatches_score()` (lines 14–20)

Annotations with adducts matching `adduct_weights[, 1]` receive a 100× score multiplier during conflict resolution. Default weighted adducts: M+H, M−H (when `adduct_weights = NULL`, all adducts get weight 5).

This strongly favors primary adducts (protonation/deprotonation) over exotic adduct assignments when multiple compounds compete for the same feature.

### MatchCategory

Each feature in the final output receives a `MatchCategory` label:

| Category | Meaning |
|----------|---------|
| **Unique** | This m/z matched exactly one compound |
| **Multiple** | This m/z initially matched multiple compounds (even if only the winner remains after filtering) |

`MatchCategory` is assigned **after** redundancy removal, so it reflects the initial ambiguity, not the final state.

### Output

`Stage5_curated_results.txt` — the final curated annotation table. Contains the same columns as Stage 4b plus `MatchCategory`.

---

## Output Files Summary

| Stage | File | Description |
|-------|------|-------------|
| 1 | `Stage1_mass_matched.txt` | Raw mass matching results (all putative adduct annotations) |
| 1 | `Stage1_peak_clusters.txt` | Peak module and RT cluster assignments |
| 2 | `Stage2_isotope_detection.txt` | Annotations expanded with detected isotope rows |
| 3 | `Stage3_chemical_scores.txt` | Chemical scores from multi-evidence evaluation |
| 3b | `Stage3_HMDB_pathways.txt` | Pathway-enriched scores (HMDB mode) |
| 3b | `Stage3_custom_pathways.txt` | Pathway-enriched scores (custom mode) |
| 3b | `Stage3_pathway_skipped.txt` | Scores without pathway enrichment (skip mode) |
| 4a | `Stage4a_confidence_levels.txt` | All rows with confidence levels (audit trail) |
| 4b | `Stage4b_confidence_levels.txt` | **Canonical final output** — coherent rows with confidence values backed by full evidence |
| 5 | `Stage5_curated_results.txt` | Opt-in (`redundancy_filtering = TRUE`) — one row per `(mz, time)`, with losing-compound rows deleted |

All files are tab-delimited text. Columns vary by stage but always include: mz, time/rt, compound_id, Adduct, score (from Stage 3 onward), Confidence and Confidence_Level (from Stage 4 onward).

---

## Parameter Decision Guide

Before diving into the full parameter reference table, use the sections below to identify which parameters to adjust for your specific situation.

### Adjust for your instrument

- **`mass_tolerance`** — Decrease to 2–3 ppm for Orbitrap and other high-resolution instruments where mass accuracy is tight. Increase to 10–15 ppm for QTOF or lower-resolution instruments where wider windows are needed to capture true matches.
- **`adduct_table`** — Filter the default 110-adduct table down to adducts relevant to your ionization mode and chromatographic conditions. Running positive-mode C18 data against all 110 adducts (including negative-mode and exotic species) generates unnecessary false matches.

### Adjust for your study design

- **`n_workers`** — Set to the number of available CPU cores minus one; more cores speed up the WGCNA module detection step, which is the most computationally expensive part of the pipeline.
- **`pathway_mode`** — Use `"HMDB"` for human metabolomics studies where compound IDs are HMDB identifiers. Use `"custom"` with your own pathway table for non-HMDB databases (e.g., pesticide libraries, plant metabolite databases). Use `"skip"` for exposomics or targeted panels where pathway enrichment is not meaningful.
- **`boosted_compounds`** — Provide a data frame of verified standards (compound IDs with known m/z and RT) to anchor confidence calibration; these are elevated to Confidence 4 and serve as positive controls in your results.

### Tune if results look wrong

- **`correlation_threshold`** — Decrease to 0.5–0.6 if WGCNA produces too few modules (most features land in a single module). Increase to 0.8–0.9 if modules are too large and mix unrelated features.
- **`min_cluster_size`** — Decrease below 10 if small datasets produce too few modules. Increase above 10 to reduce noise modules in large studies with thousands of features.
- **`time_tolerance`** — Decrease below 10 seconds for narrow chromatographic peaks (e.g., UHPLC C18). Increase above 10 seconds for broad peaks or HILIC separations where co-eluting adducts span a wider RT range.
- **`intensity_deviation_tolerance`** — Increase to 0.2–0.3 for noisy or low-abundance data where isotope intensity ratios are imprecise. Decrease to 0.05 for clean data from standards or high-abundance metabolites.
- **`filter_by`** — Set to the expected primary adduct for your ionization mode (`"M+H"` for positive, `"M-H"` for negative) to gate confidence levels. Set to `NULL` for unbiased scoring where no single adduct is expected to dominate.

---

## Key Parameters Reference

### `advanced_annotation()` parameters

| Parameter | Default | Stage | Description |
|-----------|---------|:-----:|-------------|
| `mass_tolerance` | 5e-6 (5 ppm) | 1 | Relative mass tolerance for adduct matching |
| `adduct_table` | 110 built-in adducts | 1 | Adduct definitions (charge, factor, mass) |
| `feature_id_column` | NULL | all | Column name for custom feature IDs to preserve in all stage outputs |
| `correlation_threshold` | 0.7 | 1.5, 3 | Module merge threshold; minimum correlation for scoring evidence |
| `deep_split` | 2 | 1.5 | WGCNA tree cutting depth (0=conservative, 4=aggressive) |
| `min_cluster_size` | 10 | 1.5 | Minimum peaks per WGCNA module |
| `network_type` | "unsigned" | 1.5 | WGCNA network topology type |
| `peak_rt_width` | 1 | 1.5 | Gaussian KDE bandwidth for RT clustering (seconds) |
| `intensity_deviation_tolerance` | 0.1 | 2 | Allowed deviation in isotope intensity ratio (10%) |
| `mass_defect_tolerance` | 0.1 | 2 | Mass defect window (Da) |
| `mass_defect_precision` | 0.01 | 2 | Bin width for mass defect calculation |
| `isotope_mass_tolerance` | NULL (→ `mass_tolerance`) | 2 | PPM tolerance for isotope mass matching |
| `maximum_isotopes` | 10 | 4 | Maximum isotope peaks per compound |
| `time_tolerance` | 10 | 2, 3, 3b, 4 | RT window for grouping features (seconds) |
| `MplusH_abundance_ratio_check` | TRUE | 3 | Require M+H/M−H to be most abundant adduct |
| `multimer_abundance_check` | TRUE | 4 | Penalize multimers (2M, 3M) more abundant than monomer |
| `adduct_weights` | NULL (all adducts, weight=5) | 3, 3b, 4, 5 | Adduct priority weights for scoring and conflict resolution |
| `filter_by` | c("M−H", "M+H") | 3, 4 | Expected primary adducts (gates certain confidence levels) |
| `min_ions_per_chemical` | 2 | 4 | Minimum annotation rows per compound |
| `level1_primary_adducts` | c("M+H", "M−H") | 4 | Adducts that qualify for Confidence Level 1 |
| `pathway_mode` | "HMDB" | 3b | Pathway analysis mode: "HMDB", "custom", or "skip" |
| `pathway_data` | NULL | 3b | Custom pathway table (required if `pathway_mode = "custom"`) |
| `excluded_pathways` | NULL | 3b | Pathway IDs to exclude from enrichment analysis |
| `excluded_pathway_compounds` | NULL | 3b | Compound IDs to exclude from pathway analysis |
| `redundancy_filtering` | FALSE | 5 | Opt-in Stage 5 multi-match resolution (writes `Stage5_curated_results.txt`). Off by default; Stage 4b is the canonical final output |
| `boosted_compounds` | NULL | 4 | User-verified compound IDs (elevated to Confidence 4) |
| `boost_match_by` | c("mz", "rt") | 4 | How to match boosted compounds: by mz, rt, or both |
| `boost_mass_tolerance` | NULL (→ `mass_tolerance`) | 4 | Mass tolerance for boost matching |
| `boost_time_tolerance` | NULL (→ `time_tolerance`) | 4 | RT tolerance for boost matching |
| `identify_isotopologues_flag` | TRUE | 4 | Add isotopologue identification labels |
| `outloc` | `tempdir()` | all | Output directory for stage files |
| `n_workers` | `parallel::detectCores()` | 1.5 | Number of parallel workers for WGCNA |

### Score interpretation guide

| Score range | Typical evidence |
|:-----------:|------------------|
| 0 | No corroborating evidence; mass match only |
| 1–10 | Single adduct, low correlation |
| 10–100 | Multiple adducts or weighted adduct with good correlation |
| 100–1,000 | Strong multi-adduct evidence or isotope confirmation |
| 1,000–10,000 | Multiple weighted adducts + isotopes + high correlation |
| >10,000 | Full evidence: weighted adducts + isotopes + high confidence level boost |

---

## Common Issues

**Symptom:** Most annotations have Confidence 0 (None), even for features you expect to annotate well.
**Likely cause:** The `filter_by` parameter does not include adducts present in your `adduct_table`, so the Stage 4 decision tree hits the "no filter match" branch and assigns Confidence 0. Alternatively, `mass_tolerance` may be too tight, producing too few mass matches for multi-evidence scoring to work.
**Fix:** Verify that `filter_by` contains an adduct that is also in your `adduct_table` (e.g., `filter_by = "M+H"` for positive mode). If mass matches are sparse, widen `mass_tolerance` (e.g., from 5e-6 to 10e-6).

---

**Symptom:** All or most features are assigned to a single WGCNA module, and `Stage1_peak_clusters.txt` shows only one or two module numbers.
**Likely cause:** `correlation_threshold` is too low (modules that should be separate get merged) or `min_cluster_size` is too large (small modules are absorbed into larger ones).
**Fix:** Increase `correlation_threshold` to 0.8 or 0.9 so only highly correlated features merge. Decrease `min_cluster_size` below 10 if your dataset has fewer features. Also check that the peak table has enough samples (columns) for meaningful correlation — WGCNA needs at least ~15–20 samples to produce stable modules.

---

**Symptom:** No pathway enrichment is occurring — scores in `Stage3_HMDB_pathways.txt` are identical to `Stage3_chemical_scores.txt`.
**Likely cause:** `pathway_mode` is set to `"skip"`, or compound IDs in your database do not match the HMDB identifier format (e.g., "HMDB0000001") used by the internal pathway table.
**Fix:** Set `pathway_mode = "HMDB"` and ensure your compound table's `compound_id` column uses HMDB-format IDs. If using a non-HMDB database, set `pathway_mode = "custom"` and provide a `pathway_data` data frame with `compound` and `pathway` columns mapping your compound IDs to pathway identifiers.

---

**Symptom:** Chemical scores are uniformly low (all below 10), resulting in most annotations stuck at Confidence 0 or 1.
**Likely cause:** Few adducts are being detected per compound (most compounds match only one adduct), correlations are weak, or `time_tolerance` is too narrow for features to be grouped together. When compounds have only a single adduct and no isotope confirmation, the base score formula produces small values.
**Fix:** Widen `time_tolerance` (e.g., from 10 to 15–20 seconds) to allow more features to be grouped. Check that your `adduct_table` includes the adducts expected for your chromatographic conditions. If correlations are weak, verify that your peak table has sufficient biological replicates (not just technical replicates or QC samples).

---

**Symptom:** `Stage5_curated_results.txt` has far fewer rows than `Stage4b_confidence_levels.txt` (only present when `redundancy_filtering = TRUE`).
**Likely cause:** Stage 5 deletes losing-compound rows at every multi-match `(mz, time)` to leave a single row per peak. The `MatchCategory` column distinguishes features that were unique matches from those that had competition. Stage 4b retains the full set of competing annotations — use Stage 4b as the basis for any analysis that depends on confidence-level evidence integrity.
**Fix:** No fix needed if you opted in to Stage 5 deliberately. To see what was dropped, diff `Stage4b_confidence_levels.txt` against `Stage5_curated_results.txt`. If the deletion was unintended, set `redundancy_filtering = FALSE` (the default).
