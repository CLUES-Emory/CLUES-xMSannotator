###############################################################################
# xMSannotator CLUES Example Run Script
# Uses: CLUES.xMSannotator::advanced_annotation()
#
# This script demonstrates a complete annotation workflow including:
#   - Loading and formatting a feature table with sample mapfile
#   - Loading a compound database
#   - Configuring adducts for positive or negative mode
#   - Running advanced_annotation() with all available parameters
#   - Post-processing and exporting results
#
# Modify the Configuration section (Section 2) to match your data.
# See docs/xMSannotator_Input_Formats.md for detailed file format specs.
###############################################################################

# =============================================================================
# Section 1: Load Packages
# =============================================================================
library(CLUES.xMSannotator)
library(readxl)
library(writexl)
library(dplyr)
library(data.table)


# =============================================================================
# Section 2: Configuration
# =============================================================================

# --- Study & mode settings ---
charge_type <- "pos"                # "pos" or "neg"
study_id    <- "MyStudy"
lc_column   <- "C18"               # e.g., "C18", "HILIC"
lc_mode     <- paste0(lc_column, charge_type)  # e.g., "C18pos"

# --- File paths ---
# Feature table (tab-delimited .txt with mz, time/rt, and intensity columns)
feature_table_path <- "path/to/feature_table.txt"

# Compound database (.xlsx with compound_id, name, molecular_formula, monoisotopic_mass)
compound_db_path <- "path/to/compound_database.xlsx"

# Sample mapfile (tab-delimited .txt; set to NULL if not available)
# Must have columns: File.Name and a sample type column (column 4)
mapfile_path <- "path/to/sample_mapfile.txt"

# Output directory
outloc <- "xMSannotator_output"
dir.create(outloc, recursive = TRUE, showWarnings = FALSE)

# Blank sample names to remove (matched against column 4 of mapfile)
blank_names <- c("Water", "Blank")

# Fold-change filter threshold (set to NULL to skip filtering)
# Expects a column ending in "_fc" in the feature table
fc_threshold <- 5

# Feature ID column name (set to NULL if feature table has no ID column)
feature_id_column <- "feature_id"

# Metadata file for post-processing join (set to NULL if not needed)
# If provided, will be joined to results by compound_id
db_metafile_path <- NULL


# =============================================================================
# Section 3: Load and Format Feature Table
# =============================================================================

# Read feature table
feature_table <- data.table::fread(feature_table_path, data.table = FALSE)

# Apply fold-change filter (column name varies but always ends with _fc)
if (!is.null(fc_threshold)) {
  fc_col <- grep("_fc$", colnames(feature_table), ignore.case = TRUE, value = TRUE)
  if (length(fc_col) == 1) {
    feature_table <- feature_table[feature_table[[fc_col]] >= fc_threshold, ]
    cat(sprintf("Filtered by %s >= %d: %d features retained\n", fc_col, fc_threshold, nrow(feature_table)))
  } else if (length(fc_col) > 1) {
    warning("Multiple fold-change columns found: ", paste(fc_col, collapse = ", "), ". Skipping filter.")
  }
}

# Identify sample intensity columns from mapfile
if (!is.null(mapfile_path)) {
  mapfile <- read.table(mapfile_path, sep = "\t", header = TRUE)
  mapfile$File.Name <- gsub("-", ".", mapfile$File.Name)

  # Remove blank samples (column 4 is sample type)
  blank_idx <- which(mapfile[, 4] %in% blank_names)
  if (length(blank_idx) > 0) {
    sample_names <- mapfile$File.Name[-blank_idx]
  } else {
    sample_names <- mapfile$File.Name
  }

  # Intersect with feature table columns (some mapfile entries may not be present)
  sample_cols <- intersect(sample_names, colnames(feature_table))
  missing_samples <- setdiff(sample_names, colnames(feature_table))
  if (length(missing_samples) > 0) {
    cat(sprintf("Note: %d mapfile samples not found in feature table:\n", length(missing_samples)))
    cat(paste(" ", missing_samples, collapse = "\n"), "\n")
  }
} else {
  # If no mapfile, assume all columns after mz/rt are sample intensities
  meta_cols <- c("mz", "rt", "time", feature_id_column)
  sample_cols <- setdiff(colnames(feature_table), meta_cols[meta_cols %in% colnames(feature_table)])
}

# Build peak table: feature_id, mz, rt, and sample intensity columns
# Rename 'time' -> 'rt' if needed (the API expects 'rt')
time_col <- if ("time" %in% colnames(feature_table)) "time" else "rt"
id_cols <- c("mz", time_col, sample_cols)
if (!is.null(feature_id_column) && feature_id_column %in% colnames(feature_table)) {
  id_cols <- c(feature_id_column, id_cols)
}
peak_table <- feature_table[, id_cols]
if (time_col == "time") {
  peak_table <- rename(peak_table, rt = time)
}

# Remove duplicate rows
peak_table <- unique(peak_table)

cat(sprintf("Peak table: %d features x %d samples\n", nrow(peak_table), length(sample_cols)))


# =============================================================================
# Section 4: Load Compound Database
# =============================================================================

compound_table <- read_xlsx(compound_db_path)

# Ensure expected column names
expected_cols <- c("compound_id", "name", "molecular_formula", "monoisotopic_mass")
if (!all(expected_cols %in% colnames(compound_table))) {
  # Try renaming from common alternatives (assumes first 4 columns match expected order)
  colnames(compound_table)[1:4] <- expected_cols
}

cat(sprintf("Compound database: %d compounds\n", nrow(compound_table)))


# =============================================================================
# Section 5: Configure Adducts
# =============================================================================

# Load the package's built-in adduct table
adduct_table <- CLUES.xMSannotator:::sample_adduct_table

# Select adducts based on ionization mode
if (charge_type == "pos") {
  selected_adducts <- c("M", "M+H", "M+2H", "M+NH4", "M+Na", "M+ACN+H", "2M+H", "M+H-H2O")
  filter_by <- "M+H"
} else {
  selected_adducts <- c("M-H", "M-H2O-H", "M+Na-2H", "M+Cl", "M+Hac-H", "2M-H")
  filter_by <- "M-H"
}

# Filter adduct table to selected adducts
adduct_table <- adduct_table[adduct_table$adduct %in% selected_adducts, ]

# Build adduct weights from filtered table
adduct_weights <- data.frame(
  adduct = adduct_table$adduct,
  weight = rep(1, nrow(adduct_table))
)

cat(sprintf("Using %d adducts for %s mode\n", nrow(adduct_table), charge_type))


# =============================================================================
# Section 6: Set Annotation Parameters
# =============================================================================

# Core matching tolerances
mass_tolerance      <- 5e-6     # 5 ppm (fractional, not integer ppm)
time_tolerance      <- 10       # seconds
correlation_threshold <- 0.7

# Isotope parameters
maximum_isotopes           <- 8
mass_defect_precision      <- 0.01
isotope_mass_tolerance     <- 7.5e-6  # ppm for isotope m/z matching (defaults to mass_tolerance)
mass_defect_tolerance      <- 5
intensity_deviation_tolerance <- 0.3
identify_isotopologues_flag   <- TRUE   # use enviPat to label isotopologues (13C, 15N, etc.)

# Peak and clustering parameters
peak_rt_width     <- 4
deep_split        <- 2
min_cluster_size  <- 10
network_type      <- "unsigned"

# Annotation quality parameters
MplusH_abundance_ratio_check  <- TRUE
multimer_abundance_check      <- FALSE
redundancy_filtering          <- TRUE
min_ions_per_chemical         <- 2

# Pathway settings
pathway_mode                <- "skip"   # "skip", "HMDB", or "custom"
pathway_data                <- NULL     # required if pathway_mode = "custom"
excluded_pathways           <- NULL     # pathway IDs to exclude from scoring
excluded_pathway_compounds  <- NULL     # compound IDs to exclude from pathway scoring

# Boosted compounds
boosted_compounds    <- NULL            # data.frame with compound_id, mz, rt to boost
boost_match_by       <- c("mz", "rt")   # match boosted compounds by mz, rt, or both
boost_mass_tolerance <- NULL            # defaults to mass_tolerance
boost_time_tolerance <- NULL            # defaults to time_tolerance

# Parallel processing
n_workers <- parallel::detectCores() - 1


# =============================================================================
# Section 7: Run Annotation
# =============================================================================

cat("\n")
cat("===========================================================\n")
cat(sprintf(" Starting annotation: %s | %s | %s\n", study_id, lc_mode, charge_type))
cat(sprintf(" Features: %d | Compounds: %d | Adducts: %d\n",
            nrow(peak_table), nrow(compound_table), nrow(adduct_table)))
cat(sprintf(" Output: %s\n", outloc))
cat("===========================================================\n")
cat(format(Sys.time(), "%a %b %d %X %Y"), "\n\n")

system.time(
  annotation <- advanced_annotation(
    peak_table                    = peak_table,
    compound_table                = compound_table,
    adduct_table                  = adduct_table,
    adduct_weights                = adduct_weights,
    mass_tolerance                = mass_tolerance,
    time_tolerance                = time_tolerance,
    correlation_threshold         = correlation_threshold,
    isotope_mass_tolerance        = isotope_mass_tolerance,
    mass_defect_tolerance         = mass_defect_tolerance,
    mass_defect_precision         = mass_defect_precision,
    intensity_deviation_tolerance = intensity_deviation_tolerance,
    peak_rt_width                 = peak_rt_width,
    deep_split                    = deep_split,
    min_cluster_size              = min_cluster_size,
    maximum_isotopes              = maximum_isotopes,
    min_ions_per_chemical         = min_ions_per_chemical,
    filter_by                     = filter_by,
    network_type                  = network_type,
    MplusH_abundance_ratio_check  = MplusH_abundance_ratio_check,
    multimer_abundance_check      = multimer_abundance_check,
    redundancy_filtering          = redundancy_filtering,
    identify_isotopologues_flag   = identify_isotopologues_flag,
    pathway_mode                  = pathway_mode,
    pathway_data                  = pathway_data,
    excluded_pathways             = excluded_pathways,
    excluded_pathway_compounds    = excluded_pathway_compounds,
    boosted_compounds             = boosted_compounds,
    boost_match_by                = boost_match_by,
    boost_mass_tolerance          = boost_mass_tolerance,
    boost_time_tolerance          = boost_time_tolerance,
    feature_id_column             = feature_id_column,
    outloc                        = outloc,
    n_workers                     = n_workers
  )
)

cat("\n", format(Sys.time(), "%a %b %d %X %Y"), "\n")


# =============================================================================
# Section 8: Post-Processing & Export
# =============================================================================

# Read the final stage output
if (redundancy_filtering) {
  results_file <- file.path(outloc, "Stage5_curated_results.txt")
} else {
  results_file <- file.path(outloc, "Stage4_confidence_levels.txt")
}

results <- read.table(results_file, sep = "\t", header = TRUE)

# Join with compound metadata if a metafile is provided
if (!is.null(db_metafile_path)) {
  db_metafile <- read_xlsx(db_metafile_path)
  results <- left_join(results, db_metafile, by = "compound_id")
}

# Export to Excel
output_name <- file.path(outloc, paste0(lc_mode, "-", study_id, "_Annotation.xlsx"))

write_xlsx(results, output_name)
cat(sprintf("\nResults exported to: %s\n", output_name))
cat(sprintf("Total annotations: %d\n", nrow(results)))
