# xMSannotator Input File Formats Guide

This document describes the required input file formats for the CLUES.xMSannotator package.

---

## Table of Contents

1. [Compound Table](#1-compound-table-required)
2. [Peak Table](#2-peak-table-required)
3. [Adduct Table](#3-adduct-table-optional)
4. [Adduct Weights](#4-adduct-weights-optional)
5. [Pathway Database](#5-pathway-database-optional)
6. [Boosted Compounds](#6-boosted-compounds-optional)
7. [Expected Adducts](#7-expected-adducts-optional)
8. [File Format Support](#8-file-format-support)
9. [Converting Existing Databases](#9-converting-existing-databases)
10. [Advanced Annotation Function Parameters](#10-advanced-annotation-function-parameters)
11. [Complete Working Example](#11-complete-working-example)

---

## 1. Compound Table (Required)

The compound database containing metabolites to match against your peaks.

### Required Columns

| Column | Type | Description | Validation |
|--------|------|-------------|------------|
| `monoisotopic_mass` | numeric | Exact monoisotopic mass (Da) | Must be numeric |
| `molecular_formula` | character | Molecular formula | Must be character string |
| `name` | character | Compound name | Must be character string |

Plus ONE of the following identifier columns:

| Column | Type | Description | Validation |
|--------|------|-------------|------------|
| `compound_id` | character | **Recommended.** User-defined compound identifier (e.g., "HMDB0000001", "C00001") | Must be unique |
| `compound` | numeric | Legacy integer compound identifier | Must be unique, integers only |

### Identifier Behavior

- **If `compound_id` is provided**: Your identifiers flow through the entire pipeline and appear as `compound_id` in all output files. An internal integer `compound` column is auto-generated.
- **If only `compound` is provided**: Legacy mode. The `compound_id` in outputs will be formatted as "Formula_1", "Formula_2", etc.

### Example: With compound_id (Recommended)

```r
compound_table <- data.frame(
  compound_id = c("HMDB0000122", "HMDB0000243", "HMDB0000169", "HMDB0000044", "HMDB0000094"),
  monoisotopic_mass = c(180.0634, 132.0423, 146.0579, 176.0321, 192.0270),
  molecular_formula = c("C6H12O6", "C4H8O5", "C5H10O5", "C6H8O6", "C6H8O7"),
  name = c("Glucose", "Threonic acid", "Ribose", "Ascorbic acid", "Citric acid")
)
```

Output `compound_id` column will contain: "HMDB0000122", "HMDB0000243", etc.

### Example: Legacy Mode (numeric compound)

```r
compound_table <- data.frame(
  compound = c(1, 2, 3, 4, 5),
  monoisotopic_mass = c(180.0634, 132.0423, 146.0579, 176.0321, 192.0270),
  molecular_formula = c("C6H12O6", "C4H8O5", "C5H10O5", "C6H8O6", "C6H8O7"),
  name = c("Glucose", "Threonic acid", "Ribose", "Ascorbic acid", "Citric acid")
)
```

Output `compound_id` column will contain: "Formula_1", "Formula_2", etc.

### Example: CSV File with compound_id

```csv
compound_id,monoisotopic_mass,molecular_formula,name
HMDB0000122,180.0634,C6H12O6,Glucose
HMDB0000243,132.0423,C4H8O5,Threonic acid
HMDB0000169,146.0579,C5H10O5,Ribose
HMDB0000044,176.0321,C6H8O6,Ascorbic acid
HMDB0000094,192.0270,C6H8O7,Citric acid
```

### Example: Parquet Loading

```r
compound_table <- load_compound_table_parquet("compounds.parquet")
```

### Notes

- Use `compound_id` for meaningful identifiers that will appear in your output files
- `monoisotopic_mass` is the neutral exact mass, NOT the m/z value
- `molecular_formula` must follow standard notation (e.g., "C6H12O6", not "C6 H12 O6")
- The formula is used for:
  - Golden rules validation (element ratio checks)
  - Isotope pattern calculation
  - Water-loss adduct validation

---

## 2. Peak Table (Required)

Your LC-MS peak data with m/z values, retention times, and intensities.

### Required Columns

| Column | Type | Description | Validation |
|--------|------|-------------|------------|
| `mz` | numeric | Measured m/z value | Required |
| `rt` | numeric | Retention time (seconds) | Required |
| `peak` | integer | Unique peak identifier | Auto-generated if missing |

### Optional Columns (for advanced annotation)

| Column | Type | Description |
|--------|------|-------------|
| `<sample1>` | numeric | Intensity for sample 1 |
| `<sample2>` | numeric | Intensity for sample 2 |
| ... | numeric | Additional sample intensities |

### Example: Minimal Peak Table

```r
peak_table <- data.frame(
  mz = c(181.0707, 133.0496, 147.0652, 177.0394, 193.0343),
  rt = c(120.5, 85.3, 95.2, 142.8, 156.1)
)
# 'peak' column will be auto-generated
```

### Example: Full Peak Table with Intensities

```r
peak_table <- data.frame(
  mz = c(181.0707, 133.0496, 147.0652, 177.0394, 193.0343),
  rt = c(120.5, 85.3, 95.2, 142.8, 156.1),
  sample_ctrl_1 = c(50000, 32000, 18000, 45000, 28000),
  sample_ctrl_2 = c(48000, 35000, 19500, 43000, 30000),
  sample_ctrl_3 = c(52000, 30000, 17500, 47000, 26000),
  sample_treat_1 = c(25000, 38000, 22000, 42000, 35000),
  sample_treat_2 = c(27000, 40000, 20000, 44000, 33000),
  sample_treat_3 = c(24000, 36000, 23000, 40000, 37000)
)
```

### Example: CSV File

```csv
mz,rt,sample_ctrl_1,sample_ctrl_2,sample_treat_1,sample_treat_2
181.0707,120.5,50000,48000,25000,27000
133.0496,85.3,32000,35000,38000,40000
147.0652,95.2,18000,19500,22000,20000
177.0394,142.8,45000,43000,42000,44000
193.0343,156.1,28000,30000,35000,33000
```

### Notes

- Intensities are required for `advanced_annotation()` (used for correlation analysis)
- Intensities are NOT required for `simple_annotation()`
- All columns must be numeric
- The `peak` identifier is auto-generated using lexicographic ordering of (mz, rt)

### Optional: Custom Feature ID Column

If your peak table includes a custom feature identifier column (e.g., "FeatureID" with values like "C0001", "C0005"), you can preserve it through the annotation pipeline using the `feature_id_column` parameter in `advanced_annotation()`.

| Column | Type | Description |
|--------|------|-------------|
| `<custom_id>` | character/numeric | Custom feature identifier (any column name) |

Example peak table with custom feature ID:

```r
peak_table <- data.frame(
  FeatureID = c("C0001", "C0002", "C0003", "C0004", "C0005"),
  mz = c(181.0707, 133.0496, 147.0652, 177.0394, 193.0343),
  rt = c(120.5, 85.3, 95.2, 142.8, 156.1),
  sample_1 = c(50000, 32000, 18000, 45000, 28000),
  sample_2 = c(48000, 35000, 19500, 43000, 30000)
)
```

When you specify `feature_id_column = "FeatureID"` in `advanced_annotation()`, this column will be included in all stage output files (Stage 1-5).

---

## 3. Adduct Table (Optional)

Defines the adducts to consider during annotation. A default table is provided.

### Required Columns

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `adduct` | character | Adduct name | "M+H" |
| `charge` | integer | Charge state | 1, -1, 2 |
| `factor` | integer | Number of molecules | 1 for M+H, 2 for 2M+H |
| `mass` | numeric | Mass shift from neutral | 1.007276 for +H |

### m/z Calculation Formula

```
expected_mz = (factor × monoisotopic_mass + mass) / |charge|
```

### Example: Common Adducts

```r
adduct_table <- data.frame(
  adduct = c("M+H", "M+Na", "M+K", "M+NH4", "M-H", "M+Cl", "M+FA-H", "2M+H", "M-H2O+H"),
  charge = c(1L, 1L, 1L, 1L, -1L, -1L, -1L, 1L, 1L),
  factor = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L),
  mass = c(
    1.007276,    # M+H: +proton
    22.989218,   # M+Na: +sodium
    38.963158,   # M+K: +potassium
    18.033823,   # M+NH4: +ammonium
    -1.007276,   # M-H: -proton
    34.969402,   # M+Cl: +chloride
    44.998201,   # M+FA-H: +formate -H
    1.007276,    # 2M+H: dimer + proton
    -17.026549   # M-H2O+H: water loss + proton
  )
)
```

### Example: Positive Mode Only

```r
adduct_table_pos <- data.frame(
  adduct = c("M+H", "M+Na", "M+K", "M+NH4", "2M+H", "M-H2O+H"),
  charge = c(1L, 1L, 1L, 1L, 1L, 1L),
  factor = c(1L, 1L, 1L, 1L, 2L, 1L),
  mass = c(1.007276, 22.989218, 38.963158, 18.033823, 1.007276, -17.026549)
)
```

### Example: Negative Mode Only

```r
adduct_table_neg <- data.frame(
  adduct = c("M-H", "M+Cl", "M+FA-H", "M-H2O-H", "2M-H"),
  charge = c(-1L, -1L, -1L, -1L, -1L),
  factor = c(1L, 1L, 1L, 1L, 2L),
  mass = c(-1.007276, 34.969402, 44.998201, -19.01839, -1.007276)
)
```

### Common Adduct Mass Shifts Reference

| Adduct | Mass Shift | Mode |
|--------|------------|------|
| M+H | +1.007276 | Positive |
| M+Na | +22.989218 | Positive |
| M+K | +38.963158 | Positive |
| M+NH4 | +18.033823 | Positive |
| M-H | -1.007276 | Negative |
| M+Cl | +34.969402 | Negative |
| M+FA-H | +44.998201 | Negative |
| M+ACN+H | +42.033823 | Positive |
| M-H2O+H | -17.026549 | Positive |
| M-H2O-H | -19.01839 | Negative |

---

## 4. Adduct Weights (Optional)

Assigns priority weights to adducts for confidence scoring and filtering.

### Required Columns

| Column | Type | Description |
|--------|------|-------------|
| `Adduct` or `adduct` | character | Adduct name |
| `Weight` or `weight` | numeric | Priority weight (higher = more important) |

### Example: Standard Weights

```r
adduct_weights <- data.frame(
  Adduct = c("M+H", "M-H", "M+Na", "M+K", "M+NH4", "M+Cl"),
  Weight = c(5, 5, 3, 2, 3, 2)
)
```

### Example: Prioritize Protonated Forms

```r
adduct_weights <- data.frame(
  Adduct = c("M+H", "M-H"),
  Weight = c(10, 10)
)
```

### Notes

- Higher weights give higher priority during redundancy filtering
- Compounds with higher-weighted adducts are retained preferentially
- Default if not provided: M+H and M-H with weight 1

---

## 5. Pathway Database (Optional)

For pathway enrichment scoring. The `advanced_annotation()` function supports three pathway modes:

1. **HMDB mode** (default): Uses built-in HMDB pathway data
2. **Custom mode**: Uses user-provided pathway data
3. **Skip mode**: Bypasses pathway matching entirely

### Required Columns for Custom Pathway Data

| Column | Type | Description |
|--------|------|-------------|
| `compound` | character | Compound identifier (must match `compound_table$compound`) |
| `pathway` | character | Pathway identifier |

### Example: Custom Pathway Mapping

```r
# Custom pathway data for use with pathway_mode = "custom"
pathway_data <- data.frame(
  compound = c("1", "1", "2", "2", "3"),  # Must match your compound_table IDs
  pathway = c("glycolysis", "tca_cycle", "glycolysis", "pentose_phosphate", "tca_cycle")
)
```

### Example: KEGG-style Mapping

```r
pathway_data <- data.frame(
  compound = c("1", "1", "2", "2", "3"),
  pathway = c("map00010", "map00500", "map00010", "map00020", "map00020")
)
```

### Using Custom Pathways in advanced_annotation()

```r
# Option 1: Skip pathway matching entirely (simplest for custom databases)
result <- advanced_annotation(
  peak_table = my_peaks,
  compound_table = my_compounds,
  pathway_mode = "skip"
)

# Option 2: Provide custom pathway data
my_pathways <- data.frame(
  compound = c("1", "1", "2", "3"),
  pathway = c("glycolysis", "tca_cycle", "glycolysis", "pentose")
)

result <- advanced_annotation(
  peak_table = my_peaks,
  compound_table = my_compounds,
  pathway_mode = "custom",
  pathway_data = my_pathways
)

# Option 3: Use HMDB pathways (default, requires HMDB compound IDs)
result <- advanced_annotation(
  peak_table = my_peaks,
  compound_table = hmdb_compounds  # Must have HMDB IDs
)
```

### Filtering Pathways

You can exclude specific pathways or compounds from the analysis:

```r
result <- advanced_annotation(
  peak_table = my_peaks,
  compound_table = my_compounds,
  pathway_mode = "custom",
  pathway_data = my_pathways,
  excluded_pathways = c("pathway_to_skip"),           # Pathways to exclude
  excluded_pathway_compounds = c("compound_to_skip")  # Compounds to exclude
)
```

### Notes

- Same compound can appear in multiple pathways (multiple rows)
- The `compound` column values must match your `compound_table$compound` values (converted to character)
- The enrichment algorithm boosts scores for compounds in shared pathways
- Use `pathway_mode = "skip"` when using custom compound databases without pathway information

---

## 6. Boosted Compounds (Optional)

List of known/validated compounds to boost to confidence level 4. Boosted compounds receive Confidence=4 and their scores are multiplied by 100.

### Required Columns

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `compound_id` | character | **Yes** | Must match `compound_id` in compound_table (appears as `compound_id` in outputs) |
| `mz` | numeric | Conditional | m/z value for proximity matching (required if "mz" in `boost_match_by`) |
| `rt` | numeric | Conditional | Retention time in seconds (required if "rt" in `boost_match_by`) |

### Related Parameters in advanced_annotation()

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `boosted_compounds` | data.frame | `NULL` | Table of compounds to boost |
| `boost_match_by` | character | `c("mz", "rt")` | Which columns to use for matching: `c("mz")`, `c("rt")`, or `c("mz", "rt")` |
| `boost_mass_tolerance` | numeric | same as `mass_tolerance` | Fractional tolerance for mz matching (e.g., `5e-6` = 5 ppm) |
| `boost_time_tolerance` | numeric | same as `time_tolerance` | Seconds tolerance for RT matching |

### Matching Behavior

- **mz + rt matching** (`boost_match_by = c("mz", "rt")`): Annotation must match compound_id AND be within mz tolerance AND within RT tolerance
- **rt-only matching** (`boost_match_by = c("rt")`): Annotation must match compound_id AND be within RT tolerance
- **mz-only matching** (`boost_match_by = c("mz")`): Annotation must match compound_id AND be within mz tolerance

### Example: With mz and RT Matching (Default)

```r
# Boost compounds that match by ID, mz, and RT
boosted_compounds <- data.frame(
  compound_id = c("HMDB0000122", "HMDB0000158", "HMDB0000167"),
  mz = c(181.0707, 147.0652, 133.0496),
  rt = c(120.5, 95.2, 85.3)
)

result <- advanced_annotation(
  peak_table = my_peaks,
  compound_table = my_compounds,
  boosted_compounds = boosted_compounds,
  boost_match_by = c("mz", "rt"),    # default
  boost_mass_tolerance = 5e-6,       # 5 ppm (defaults to mass_tolerance if not specified)
  boost_time_tolerance = 10          # 10 seconds (defaults to time_tolerance if not specified)
)
```

### Example: RT-Only Matching

```r
# Boost compounds that match by ID and RT only (ignore mz differences)
boosted_compounds <- data.frame(
  compound_id = c("HMDB0000122", "HMDB0000158", "HMDB0000167"),
  rt = c(120.5, 95.2, 85.3)  # mz column not needed
)

result <- advanced_annotation(
  peak_table = my_peaks,
  compound_table = my_compounds,
  boosted_compounds = boosted_compounds,
  boost_match_by = c("rt"),          # RT-only matching
  boost_time_tolerance = 15          # 15 second tolerance
)
```

### Example: mz-Only Matching

```r
# Boost compounds that match by ID and mz only (ignore RT differences)
boosted_compounds <- data.frame(
  compound_id = c("HMDB0000122", "HMDB0000158", "HMDB0000167"),
  mz = c(181.0707, 147.0652, 133.0496)  # rt column not needed
)

result <- advanced_annotation(
  peak_table = my_peaks,
  compound_table = my_compounds,
  boosted_compounds = boosted_compounds,
  boost_match_by = c("mz"),          # mz-only matching
  boost_mass_tolerance = 10e-6       # 10 ppm tolerance
)
```

### Notes

- The `compound_id` column must match the `compound_id` values in your compound_table
- Boosted annotations receive Confidence=4 (highest level) and score×100
- Tolerance parameters use the same format as main parameters: fractional for mass (e.g., `5e-6` = 5 ppm), seconds for time

---

## 7. Expected Adducts (Optional)

List of adducts that are expected/required for high confidence annotations.

### Required Columns

| Column | Type | Description |
|--------|------|-------------|
| `adduct` | character | Adduct name |

### Example

```r
expected_adducts <- data.frame(
  adduct = c("M+H", "M-H")
)
```

### Loading from CSV

```r
expected_adducts <- load_expected_adducts_csv("expected_adducts.csv")
```

---

## 8. File Format Support

### Supported Formats

| Format | Read Function | Write Function | Notes |
|--------|---------------|----------------|-------|
| Parquet | `load_*_parquet()` | `save_parquet()` | Recommended for large files |
| CSV | `read.csv()` / `readr::read_csv()` | `write.csv()` | Universal compatibility |
| RDA | `load()` | `save()` | Native R format |
| Data Frame | Direct use | - | In-memory |

### Parquet Loading Functions

```r
# Compound table
compound_table <- load_compound_table_parquet("compounds.parquet")

# Peak table
peak_table <- load_peak_table_parquet("peaks.parquet")

# Adduct table
adduct_table <- load_adduct_table_parquet("adducts.parquet")
```

### Saving to Parquet

```r
save_parquet(compound_table, "compounds.parquet")
save_parquet(peak_table, "peaks.parquet")
save_parquet(annotation_results, "results.parquet")
```

---

## 9. Converting Existing Databases

### From HMDB (hmdbAllinf.rda)

```r
# Load original HMDB data
load("xMSannotator-master/data/hmdbAllinf.rda")

# Convert to required format
compound_table <- data.frame(
  compound = seq_len(nrow(hmdbAllinf)),
  monoisotopic_mass = as.numeric(hmdbAllinf$MonoisotopicMass),
  molecular_formula = as.character(hmdbAllinf$Formula),
  name = as.character(hmdbAllinf$Name)
)

# Save as Parquet
save_parquet(compound_table, "hmdb_compounds.parquet")
```

### From KEGG

```r
# Assuming KEGG data with columns: KEGG_ID, ExactMass, Formula, Name
kegg_data <- read.csv("kegg_compounds.csv")

compound_table <- data.frame(
  compound = seq_len(nrow(kegg_data)),
  monoisotopic_mass = kegg_data$ExactMass,
  molecular_formula = kegg_data$Formula,
  name = kegg_data$Name
)
```

### From LipidMaps

```r
# Assuming LipidMaps data
lipidmaps_data <- read.csv("lipidmaps.csv")

compound_table <- data.frame(
  compound = seq_len(nrow(lipidmaps_data)),
  monoisotopic_mass = lipidmaps_data$EXACT_MASS,
  molecular_formula = lipidmaps_data$FORMULA,
  name = lipidmaps_data$COMMON_NAME
)
```

### From Custom CSV

```r
# Your CSV with any column names
my_data <- read.csv("my_database.csv")

compound_table <- data.frame(
  compound = seq_len(nrow(my_data)),
  monoisotopic_mass = my_data$mass,  # Map your column name
  molecular_formula = my_data$formula,
  name = my_data$compound_name
)
```

### From MassBank/mzCloud Export

```r
# Typical spectral library export format
massbank <- read.csv("massbank_export.csv")

compound_table <- data.frame(
  compound = seq_len(nrow(massbank)),
  monoisotopic_mass = massbank$EXACT_MASS,
  molecular_formula = massbank$MOLECULAR_FORMULA,
  name = massbank$COMPOUND_NAME
)
```

---

## 10. Advanced Annotation Function Parameters

The `advanced_annotation()` function accepts the following parameters:

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `peak_table` | data.frame | Feature table with m/z, retention time, and intensity columns |
| `compound_table` | data.frame | Database of compounds with names, formulas, and monoisotopic masses |

### Optional Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adduct_table` | data.frame | `NULL` | Table of adducts with name, factor, charge, and mass columns. Uses `sample_adduct_table` if NULL |
| `adduct_weights` | data.frame | `NULL` | Weights for prioritizing specific adducts. Auto-generated with weight=5 if NULL |
| `feature_id_column` | character | `NULL` | Name of column in peak_table containing custom feature identifiers to preserve in all stage outputs |
| `intensity_deviation_tolerance` | numeric | `0.1` | Tolerance for intensity deviation in isotope matching (10%) |
| `mass_tolerance` | numeric | `5e-6` | Mass accuracy as fractional (relative) tolerance. Use `5e-6` for 5 ppm, `10e-6` for 10 ppm. Do NOT enter direct ppm values. |
| `isotope_mass_tolerance` | numeric | `NULL` | Fractional tolerance for isotope m/z matching. Defaults to `mass_tolerance`. Use `5e-6` for 5 ppm. |
| `mass_defect_tolerance` | numeric | `0.1` | Tolerance for mass defect matching |
| `mass_defect_precision` | numeric | `0.01` | Precision for binning mass defects |
| `time_tolerance` | numeric | `10` | Retention time tolerance in seconds |
| `peak_rt_width` | numeric | `1` | Expected chromatographic peak width for RT clustering |
| `correlation_threshold` | numeric | `0.7` | Minimum correlation for peak grouping |
| `MplusH_abundance_ratio_check` | logical | `TRUE` | Requires secondary adducts to have lower intensity than M+H/M-H during chemical scoring. Set to FALSE to disable. |
| `multimer_abundance_check` | logical | `TRUE` | Checks that multimer adducts (2M, 3M) have lower intensity than the monomer during confidence assignment. Set to FALSE to disable. |
| `deep_split` | integer | `2` | WGCNA parameter controlling cluster splitting (0-4, higher = more clusters) |
| `min_cluster_size` | integer | `10` | Minimum peaks per module in WGCNA clustering |
| `maximum_isotopes` | integer | `10` | Maximum isotope peaks to consider per compound |
| `min_ions_per_chemical` | integer | `2` | Minimum ions required to annotate a chemical |
| `filter_by` | character vector | `c("M-H", "M+H")` | Primary adducts for confidence scoring |
| `network_type` | character | `"unsigned"` | WGCNA network type ("unsigned", "signed", or "signed hybrid") |
| `redundancy_filtering` | logical | `TRUE` | Whether to remove redundant annotations |
| `pathway_mode` | character | `"HMDB"` | Pathway matching mode: "HMDB" (default), "custom", or "skip" |
| `pathway_data` | data.frame | `NULL` | Custom pathway-compound mappings (required if pathway_mode = "custom") |
| `excluded_pathways` | character vector | `NULL` | Pathways to exclude from analysis |
| `excluded_pathway_compounds` | character vector | `NULL` | Compounds to exclude from pathway analysis |
| `boosted_compounds` | data.frame | `NULL` | Table of confirmed compounds to boost to Confidence=4 (see Section 6) |
| `boost_match_by` | character vector | `c("mz", "rt")` | Which columns to use for boost matching: `c("mz")`, `c("rt")`, or `c("mz", "rt")` |
| `boost_mass_tolerance` | numeric | same as `mass_tolerance` | Fractional tolerance for boost mz matching (e.g., `5e-6` = 5 ppm) |
| `boost_time_tolerance` | numeric | same as `time_tolerance` | Seconds tolerance for boost RT matching |
| `enable_permutation` | logical | `FALSE` | Enable permutation-based p-value calculation for annotation scores |
| `n_permutations` | integer | `1000` | Number of permutations for significance testing |
| `permutation_method` | character | `"full"` | Permutation method: `"full"` (all permutations in parallel, faster) or `"streaming"` (chunked processing, lower memory) |
| `permutation_seed` | integer | `42` | Random seed for reproducibility of permutation results |
| `outloc` | character | `tempdir()` | Output directory for intermediate files |
| `n_workers` | integer | `detectCores()` | Number of parallel workers for WGCNA and permutation testing |

**Note on mass tolerance format:** All mass tolerance parameters (`mass_tolerance`, `isotope_mass_tolerance`) use fractional (relative) tolerance notation:
- `5e-6` = 5 ppm (5 parts per million)
- `10e-6` = 10 ppm
- The matching algorithm uses: `|observed - expected| ≤ max(|observed|, |expected|) × tolerance`
- Do NOT enter ppm values directly (e.g., do NOT use `5` for 5 ppm)

### Example Usage

```r
library(CLUES.xMSannotator)

# Run with default parameters
result <- advanced_annotation(
  peak_table = my_peaks,
  compound_table = my_compounds
)

# Run with custom parameters
result <- advanced_annotation(
  peak_table = my_peaks,
  compound_table = my_compounds,
  adduct_table = my_adducts,
  mass_tolerance = 10e-6,           # 10 ppm (fractional: 10 × 10^-6)
  isotope_mass_tolerance = 5e-6,    # 5 ppm for isotope matching
  time_tolerance = 15,              # 15 seconds
  correlation_threshold = 0.8,      # stricter correlation
  min_ions_per_chemical = 3,        # require more ions
  filter_by = c("M+H"),             # positive mode only
  n_workers = 4                     # limit parallelization
)

# Preserve custom feature IDs through the pipeline
result <- advanced_annotation(
  peak_table = my_peaks,           # has "FeatureID" column
  compound_table = my_compounds,
  feature_id_column = "FeatureID", # preserve this column in all outputs
  outloc = "output/"
)
# All Stage 1-5 output files will include the "FeatureID" column

# Skip pathway matching (for custom compound databases)
result <- advanced_annotation(
  peak_table = my_peaks,
  compound_table = custom_compounds,
  pathway_mode = "skip"            # bypass HMDB pathway matching
)

# Use custom pathway database
my_pathways <- data.frame(
  compound = c("1", "1", "2", "3"),
  pathway = c("glycolysis", "tca_cycle", "glycolysis", "pentose")
)

result <- advanced_annotation(
  peak_table = my_peaks,
  compound_table = my_compounds,
  pathway_mode = "custom",
  pathway_data = my_pathways,
  excluded_pathways = c("unwanted_pathway")  # optional filtering
)

# Boost confidence of confirmed/validated compounds
validated_compounds <- data.frame(
  compound_id = c("HMDB0000122", "HMDB0000243"),  # Must match compound_table$compound_id
  mz = c(181.0707, 133.0496),                     # Observed m/z values
  rt = c(120.5, 85.3)                             # Observed retention times (seconds)
)

result <- advanced_annotation(
  peak_table = my_peaks,
  compound_table = my_compounds,
  boosted_compounds = validated_compounds,
  boost_match_by = c("rt"),           # Match by compound_id + RT only
  boost_time_tolerance = 15,          # 15 second RT tolerance
  pathway_mode = "skip"
)
# Annotations matching validated compounds will have Confidence=4 and score×100

# Run with permutation-based significance testing (full parallel method - default)
result <- advanced_annotation(
  peak_table = my_peaks,
  compound_table = my_compounds,
  enable_permutation = TRUE,          # Enable p-value calculation
  n_permutations = 1000,              # Number of permutations (use 100 for quick tests)
  permutation_method = "full",        # All permutations in parallel (default, faster)
  permutation_seed = 42,              # For reproducibility
  n_workers = 4,                      # Parallel processing
  outloc = "output/"
)
# Result will have perm_pvalue column
# Creates Stage4_permutation_pvalues.txt in addition to Stage4_confidence_levels.txt

# For large datasets with memory constraints, use streaming method
result <- advanced_annotation(
  peak_table = my_peaks,
  compound_table = my_compounds,
  enable_permutation = TRUE,
  n_permutations = 1000,
  permutation_method = "streaming",   # Process in chunks (lower memory usage)
  n_workers = 4,
  outloc = "output/"
)
```

---

## 11. Complete Working Example

### Simple Annotation

```r
library(CLUES.xMSannotator)

# 1. Create compound database
compound_table <- data.frame(
  compound = 1:5,
  monoisotopic_mass = c(180.0634, 132.0423, 146.0579, 176.0321, 192.0270),
  molecular_formula = c("C6H12O6", "C4H8O5", "C5H10O5", "C6H8O6", "C6H8O7"),
  name = c("Glucose", "Threonic acid", "Ribose", "Ascorbic acid", "Citric acid")
)

# 2. Create peak table
peak_table <- data.frame(
  mz = c(181.0707, 133.0496, 147.0652, 177.0394, 193.0343, 203.0526),
  rt = c(120.5, 85.3, 95.2, 142.8, 156.1, 165.4)
)

# 3. Define adducts (optional - has defaults)
adduct_table <- data.frame(
  adduct = c("M+H", "M+Na", "M-H"),
  charge = c(1L, 1L, -1L),
  factor = c(1L, 1L, 1L),
  mass = c(1.007276, 22.989218, -1.007276)
)

# 4. Run simple annotation
annotation <- simple_annotation(
  peak_table = peak_table,
  compound_table = compound_table,
  adduct_table = adduct_table,
  mass_tolerance = 5e-6  # 5 ppm
)

# 5. View results
print(annotation)
```

### Advanced Annotation

```r
library(CLUES.xMSannotator)

# 1. Load compound database from Parquet
compound_table <- load_compound_table_parquet("my_compounds.parquet")

# 2. Load peak table with intensities
peak_table <- load_peak_table_parquet("my_peaks.parquet")

# 3. Define adducts and weights
adduct_table <- data.frame(
  adduct = c("M+H", "M+Na", "M+K", "M+NH4"),
  charge = c(1L, 1L, 1L, 1L),
  factor = c(1L, 1L, 1L, 1L),
  mass = c(1.007276, 22.989218, 38.963158, 18.033823)
)

adduct_weights <- data.frame(
  adduct = c("M+H", "M+Na", "M+K", "M+NH4"),
  weight = c(5, 3, 2, 3)
)

# 4. Run advanced annotation with custom compound database
# Use pathway_mode = "skip" to bypass HMDB pathway matching
result <- advanced_annotation(
  peak_table = peak_table,
  compound_table = compound_table,
  adduct_table = adduct_table,
  adduct_weights = adduct_weights,
  mass_tolerance = 5e-6,              # 5 ppm
  time_tolerance = 10,                # 10 seconds
  correlation_threshold = 0.7,
  filter_by = c("M+H"),               # positive mode
  pathway_mode = "skip",              # skip HMDB pathway matching
  outloc = "output/"
)

# 5. Or use custom pathway data
my_pathways <- data.frame(
  compound = c("1", "1", "2", "2", "3"),
  pathway = c("glycolysis", "tca_cycle", "glycolysis", "pentose_phosphate", "tca_cycle")
)

result_with_pathways <- advanced_annotation(
  peak_table = peak_table,
  compound_table = compound_table,
  adduct_table = adduct_table,
  adduct_weights = adduct_weights,
  pathway_mode = "custom",
  pathway_data = my_pathways,
  outloc = "output_with_pathways/"
)

# 6. Save results
save_parquet(result, "annotation_results.parquet")
write.csv(result, "annotation_results.csv", row.names = FALSE)
```

---

## Validation Functions Reference

The package includes validation functions that check your data:

| Function | Purpose |
|----------|---------|
| `as_peak_table()` | Validates peak table format |
| `as_compound_table()` | Validates compound table format |
| `as_adduct_table()` | Validates adduct table format |
| `as_pathway_table()` | Validates custom pathway table format |
| `as_expected_adducts_table()` | Validates expected adducts |
| `load_boost_compounds_csv()` | Loads and validates boosted compounds from CSV |
| `as_boosted_compounds_table()` | Validates boosted compounds (internal, called by `load_boost_compounds_csv()`) |

### Manual Validation Example

```r
# Check your data before annotation
tryCatch({
  validated_compounds <- as_compound_table(compound_table)
  validated_peaks <- as_peak_table(peak_table, intensities = TRUE)
  validated_adducts <- as_adduct_table(adduct_table)
  print("All data validated successfully!")
}, error = function(e) {
  print(paste("Validation error:", e$message))
})
```

---

*Document created: 2026-01-22*
*Last updated: 2026-02-08*
*For use with CLUES.xMSannotator v0.10.0*
