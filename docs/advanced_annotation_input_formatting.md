# CLUES.xMSannotator Input Data Formats & Pre-Processing

> **Important:** This document describes the input formats for **CLUES.xMSannotator**, which differ from the original recetox-xMSannotator. The pipeline now allows input of XCMS feature tables generated from the CLUES data extraction workflow, sample mapfiles, and xlsx compound databases directly — recetox-aplcms conversion and HMDB RDA-to-Parquet workflows are no longer used. Feature tables from other workflows can be modified as needed for input using the formatting described below.

---

## Overview

CLUES.xMSannotator requires three input files:

1. **Feature table** — XCMS peak-detection output (tab-delimited .txt)
2. **Sample mapfile** — sample metadata for blank identification (tab-delimited .txt)
3. **Compound database** — target compounds for annotation (.xlsx)

Before annotation, the Example Runscript applies pre-processing steps: blank sample removal, fold-change filtering, and peak table construction. This document describes each file format and the pre-processing logic.

For the full API-level input format specification (compound table, peak table, adduct table, etc.), see [`xMSannotator_Input_Formats.md`](xMSannotator_Input_Formats.md).

---

## 1. Feature Table

**Format:** Tab-delimited .txt file produced by XCMS (e.g., via `XCMSv4.8.0` QA-filtered output).

### Required Columns

Only three column types are required for annotation:

| Column | Type | Description |
|--------|------|-------------|
| `mz` | numeric | Mass-to-charge ratio |
| `time` | numeric | Retention time. **Note:** the column is `time`, not `rt` — the runscript renames it to `rt` for `advanced_annotation()`. |
| Sample intensities | numeric | One column per injection, named by `File.Name` from the mapfile (e.g., `C251203_CLU0122_C18pos_001`, ...) |

**Optional but recommended:** `feature_id` (string, e.g., "FT0001") — preserved through the pipeline via the `feature_id_column` parameter.

### QC Metadata Columns (Workflow-Specific)

The remaining columns are specific to our XCMS QA pipeline output and are **not used by `advanced_annotation()`**. They are carried through the feature table for post-processing metadata joins (e.g., merging chemical annotations back with feature-level QC stats). These columns can be modified or omitted according to your feature table output.

| Category | Example Columns | Description |
|----------|-----------------|-------------|
| **XCMS peak statistics** | `mzmin`, `mzmax`, `rtmin`, `rtmax`, `n_peaks`, `n_detected`, `n_filled`, `rt_range`, `mz_range_ppm`, `int_cv` | Peak detection and alignment metrics |
| **QA scores** | `peak_quality`, `alignment_score`, `combined_score`, `quality_category`, `m_selectivity`, `c_selectivity`, `peak_shape`, `passed_filter` | Quality assessment from XCMS QA filtering |
| **Blank/detection stats** | `Water_fc`, `Water_status`, `Overall_detect_pct`, `Blank_detect_pct`, `Water_detect_pct`, `NIST_detect_pct`, `AMAP_detect_pct`, `PFAS_QAQC_detect_pct`, `HRE_Pool_detect_pct`, `Study_Sample_detect_pct` | Fold-change vs blanks, detection percentages by sample type |
| **CV columns** | `HRE.p1Std._CV`, `HRE.p2Std._CV`, `Study_Sample_CV` | Coefficient of variation (study-specific, columns vary) |

**Note:** If your feature table includes a fold-change column (any column ending in `_fc`, e.g., `Water_fc`), it is used by the pre-processing fold-change filter (Section 4.1) but is not passed to `advanced_annotation()`.

### Example (first 3 rows, abbreviated)

```
feature_id  mz        time   mzmin     mzmax     ...  Water_fc  ...  C251203_CLU0122_C18pos_001  C251203_CLU0122_C18pos_002
F00001      85.00655  2.48   85.00653  85.00658  ...  1.467     ...  451249.51                   449728.24
F00002      85.02731  97.23  85.02728  85.02735  ...  171241.4  ...  0                           0
F00005      85.02832  11.58  85.0283   85.02834  ...  44.460    ...  149682.33                   81397.55
```

### Key Notes

- Only `mz`, `time`, and sample intensity columns are extracted for the peak table passed to `advanced_annotation()`. All other columns remain in the feature table for post-processing use.
- The retention time column is named **`time`**, not `rt`. The runscript renames it to `rt` during peak table construction (required by `advanced_annotation()`).
- The `Water_fc` column (or any column ending in `_fc`) is used by the pre-processing fold-change filter, not by annotation itself.
- Sample intensity column names must match the `File.Name` column in the mapfile (after hyphen-to-dot replacement).

---

## 2. Sample Mapfile

**Format:** Tab-delimited .txt with 4 columns.

| Column | Type | Description |
|--------|------|-------------|
| `File.Name` | character | Injection file name — must match feature table column names |
| `Sample.ID` | character | Sample identifier |
| `Batch` | integer | Batch number |
| `Sample_type` | character | Sample classification (e.g., `Blank`, `Water`, `NIST`, `AMAP`, `PFAS_QAQC`, `HRE_Pool`, `Study_Sample`) |

### Example

```
File.Name                       Sample.ID             Batch  Sample_type
C251203_CLU0122_C18pos_001      Instrument_Blank_01   1      Blank
C251203_CLU0122_C18pos_002      Water_01_01_1         1      Water
C251203_CLU0122_C18pos_003      NIST1950_01_1         1      NIST
C251203_CLU0122_C18pos_004      AM-S-Y2506_01_01_1    1      AMAP
```

### Purpose

The mapfile serves two roles in pre-processing:

1. **Blank identification** — The `Sample_type` column (column 4) is matched against `blank_names` (e.g., `c("Water", "Blank")`) to identify blank samples for removal.
2. **Sample column selection** — After removing blanks, the remaining `File.Name` values are intersected with feature table column names to identify which columns contain sample intensities.

### Hyphen-to-Dot Replacement

R's `read.table()` converts hyphens to dots in column names. To match, the runscript applies `gsub("-", ".", mapfile$File.Name)` before comparing against feature table columns.

### Without a Mapfile

If `mapfile_path` is set to `NULL`, the runscript assumes all columns except `mz`, `rt`/`time`, and the feature ID column are sample intensities. No blank removal is performed.

---

## 3. Compound Database

**Format:** Excel workbook (.xlsx) with 4 columns.

| Source Column | Renamed To | Type | Description |
|---------------|------------|------|-------------|
| `Formula_ID` | `compound_id` | character | Compound identifier (e.g., "PEST0001") |
| `Compound_names` | `name` | character | Compound name |
| `Molecular_Formula` | `molecular_formula` | character | Molecular formula (e.g., "C2H4") |
| `Monoisotopic_Mass` | `monoisotopic_mass` | numeric | Exact monoisotopic mass in Daltons |

### Column Renaming

The runscript reads the xlsx and checks whether the expected column names (`compound_id`, `name`, `molecular_formula`, `monoisotopic_mass`) are present. If not, it renames the first 4 columns to match:

```r
compound_table <- read_xlsx(compound_db_path)
expected_cols <- c("compound_id", "name", "molecular_formula", "monoisotopic_mass")
if (!all(expected_cols %in% colnames(compound_table))) {
  colnames(compound_table)[1:4] <- expected_cols
}
```

This means the xlsx columns must be in the order: identifier, name, formula, mass — regardless of their original header names.

---

## 4. Pre-Processing Steps

The Example Runscript (`xMSannotator_CLUES_Runscript_Example.R`, Section 3) applies three pre-processing steps before calling `advanced_annotation()`.

### 4.1 Fold-Change Filtering

Removes features with low fold-change relative to blanks.

**Configuration:**
- `fc_threshold <- 5` — minimum fold-change to retain a feature
- `fc_threshold <- NULL` — skip filtering entirely

**Logic:**

```r
if (!is.null(fc_threshold)) {
  fc_col <- grep("_fc$", colnames(feature_table), ignore.case = TRUE, value = TRUE)
  if (length(fc_col) == 1) {
    feature_table <- feature_table[feature_table[[fc_col]] >= fc_threshold, ]
  } else if (length(fc_col) > 1) {
    warning("Multiple fold-change columns found: ..., Skipping filter.")
  }
}
```

**Behavior:**
- Searches for columns ending in `_fc` (e.g., `Water_fc`)
- If exactly 1 match: filters rows where fold-change >= threshold
- If multiple matches: warns and skips (ambiguous)
- If 0 matches: no filtering (passes silently)
- **To skip:** set `fc_threshold <- NULL`

### 4.2 Blank Sample Removal

Removes blank/control sample columns so they are not included in annotation.

**Configuration:**
- `blank_names <- c("Water", "Blank")` — sample types to treat as blanks
- `mapfile_path <- NULL` — skip blank removal (no mapfile)

**Logic:**

```r
mapfile <- read.table(mapfile_path, sep = "\t", header = TRUE)
mapfile$File.Name <- gsub("-", ".", mapfile$File.Name)

blank_idx <- which(mapfile[, 4] %in% blank_names)
if (length(blank_idx) > 0) {
  sample_names <- mapfile$File.Name[-blank_idx]
} else {
  sample_names <- mapfile$File.Name
}

sample_cols <- intersect(sample_names, colnames(feature_table))
```

**Steps:**
1. Read mapfile and replace hyphens with dots in `File.Name`
2. Identify blank rows by matching `Sample_type` (column 4) against `blank_names`
3. Remove blank file names from the sample list
4. Intersect remaining names with feature table columns to get `sample_cols`
5. Report any mapfile samples not found in the feature table

### 4.3 Peak Table Construction

Builds the peak table that `advanced_annotation()` expects.

**Logic:**

```r
time_col <- if ("time" %in% colnames(feature_table)) "time" else "rt"
id_cols <- c("mz", time_col, sample_cols)
if (!is.null(feature_id_column) && feature_id_column %in% colnames(feature_table)) {
  id_cols <- c(feature_id_column, id_cols)
}
peak_table <- feature_table[, id_cols]
if (time_col == "time") {
  peak_table <- rename(peak_table, rt = time)
}
peak_table <- unique(peak_table)
```

**Steps:**
1. Select columns: `feature_id` (if present) + `mz` + `time` + sample intensity columns
2. Rename `time` → `rt` (required by the API)
3. Deduplicate with `unique()`

**Result:** A data frame with columns `feature_id`, `mz`, `rt`, and one column per non-blank sample — ready for `advanced_annotation()`.

---

## See Also

- [`xMSannotator_Input_Formats.md`](xMSannotator_Input_Formats.md) — Full API-level input format specification (compound table, peak table, adduct table, adduct weights, pathway database, boosted compounds, expected adducts)
- [`xMSannotator_CLUES_Runscript_Example.R`](../inst/scripts/xMSannotator_CLUES_Runscript_Example.R) — Complete working example

---

*Document created: 2026-01-22*
*Last updated: 2026-03-14*
*For use with CLUES.xMSannotator v1.0.0*
