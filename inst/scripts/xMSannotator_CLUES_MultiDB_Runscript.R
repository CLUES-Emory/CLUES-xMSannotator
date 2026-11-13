###############################################################################
# xMSannotator CLUES Multi-Database SLURM Run Script
# Uses: CLUES.xMSannotator::advanced_annotation()
# Date: 2026-11-13 (CLUES.xMSannotator 1.3.0)
# Note: Designed for SLURM array jobs with database mapfile support
#       Run via: run_CLUES_xMSannotator.sh
# Output: <output_dir>/<mode>/7-<mode>_<study>_xMSannotator/
#           <mode>-<study>_xMSannotator_README.md         (written by every task, same content)
#           <File_name_out>/<mode>-<study>_<File_name_out>_CLUES_Annotation.xlsx
#           <File_name_out>/run_info.txt                   (inputs, versions, settings, counts)
#           <File_name_out>/intermediate/Stage*.txt.gz     (advanced_annotation() Stage files)
###############################################################################

run_start <- Sys.time()

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

# --- Debug/default configuration (used when no CLI args provided) ---
debug_config <- list(
  array_index        = 1,
  lc_column          = "C18",
  output_base_loc    = "/clues/CLU0046.1_STUDY/2-Data_Extraction/C18/",
  study_id           = "CLU0046.1",
  sample_mapfile_dir = "/clues/CLU0046.1_STUDY/1-mzML/C18/",
  batch_correction   = "Unadjusted",
  fc_threshold       = 5,
  n_cores            = 10,
  db_mapfile_path    = "/clues/.../Compound_Database_Mapfile_v3.xlsx",
  blank_names        = "Water,Blank"
)

# --- Parse configuration from CLI args or debug_config ---
parse_config <- function(input) {
  if (is.list(input)) {
    cfg <- input
  } else {
    cfg <- list(
      array_index        = as.numeric(input[1]),
      lc_column          = input[2],
      output_base_loc    = input[3],
      study_id           = input[4],
      sample_mapfile_dir = input[5],
      batch_correction   = input[6],
      fc_threshold       = as.numeric(input[7]),
      n_cores            = as.numeric(input[8]),
      db_mapfile_path    = input[9],
      blank_names        = input[10]
    )
  }
  cfg
}

# --- Determine config source: CLI args or debug config ---
cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) > 0) {
  cat("Using command-line arguments for configuration\n")
  config <- parse_config(cli_args)
} else {
  cat("No CLI arguments provided - using debug configuration\n")
  config <- parse_config(debug_config)
}

print(sessionInfo())
print(config)

# --- Read database mapfile to determine number of databases ---
db_ls <- read_xlsx(config$db_mapfile_path)
n_databases <- nrow(db_ls)

# --- Mode + DB index logic ---
# Array IDs 1..N = positive mode; Array IDs (N+1)..2N = negative mode
if (config$array_index <= n_databases) {
  charge_type <- "pos"
  ii <- config$array_index
} else {
  charge_type <- "neg"
  ii <- config$array_index - n_databases
}

# --- Study & mode settings ---
study_id <- config$study_id
lc_column <- config$lc_column
lc_mode <- paste0(lc_column, charge_type)

# --- File paths ---
# Sample mapfile: find {lc_mode}_sample_mapfile.txt in the mapfile directory
ls_mapfiles <- list.files(config$sample_mapfile_dir, pattern = ".txt", full.names = TRUE)
mapfile_path <- ls_mapfiles[grep(paste0(lc_mode, "_sample_mapfile.txt"), ls_mapfiles)]
cat(sprintf("Sample mapfile: %s\n", mapfile_path))

# Feature table: navigate batch correction directory structure
outloc_base <- paste0(config$output_base_loc, "/", lc_mode, "/")
bc_outloc <- paste0(outloc_base, "2-", lc_mode, "_", study_id, "_BatchCorrection/")
# The Unadjusted folder also holds the injection QC flags (*_Injection_QC_Flags.txt, Step 1 v1.3.2+),
# so match the feature table by name
ft_loc <- list.dirs(bc_outloc, full.names = TRUE)[grep("Unadjusted", list.dirs(bc_outloc, full.names = FALSE))] %>%
  list.files(pattern = "_Untargeted_RAW_Feature_Table_.*\\.txt$", full.names = TRUE)
if (length(ft_loc) != 1) {
  stop("Expected one *_Untargeted_RAW_Feature_Table_*.txt in the Unadjusted folder of ", bc_outloc,
       ", found ", length(ft_loc), ": ", paste(basename(ft_loc), collapse = ", "))
}
feature_table_path <- ft_loc
cat(sprintf("Feature table:  %s\n", feature_table_path))

# Feature groups from Step 1 Stage 3 (v1.7.0+). When present they define the feature set and replace
# xMSannotator's own WGCNA clustering; without them (older Step 1 output) WGCNA runs as before.
feature_groups_path <- paste0(outloc_base, "3-", lc_mode, "_", study_id, "_Feature_Grouping/",
                              lc_mode, "-", study_id, "_Feature_Groups.txt")
use_feature_groups <- file.exists(feature_groups_path)
cat(sprintf("Feature groups: %s\n", if (use_feature_groups) feature_groups_path else
  paste("none (", feature_groups_path, "not found): WGCNA clustering")))

# Blank sample names (comma-separated arg)
blank_names <- unlist(strsplit(config$blank_names, ","))

# Fold-change filter threshold
fc_threshold <- config$fc_threshold

# Number of parallel workers
n_workers <- config$n_cores

# Feature ID column name
feature_id_column <- "feature_id"


# =============================================================================
# Section 3: Load and Format Feature Table
# =============================================================================

# Read feature table
feature_table <- data.table::fread(feature_table_path, data.table = FALSE)

# Feature set: the Step 1 Stage 3 features (blank fold-change filter already applied there,
# YAML feature_grouping: blank_fc_threshold), or else the fold-change filter below
feature_groups <- NULL
fc_col <- grep("_fc$", colnames(feature_table), ignore.case = TRUE, value = TRUE)
if (use_feature_groups) {
  feature_groups <- data.table::fread(feature_groups_path, data.table = FALSE,
                                      colClasses = list(character = "feature_id"))
  feature_groups <- data.frame(feature_id = feature_groups$feature_id, group = feature_groups$feature_group)
  missing_ids <- setdiff(feature_groups$feature_id, feature_table$feature_id)
  if (length(missing_ids) > 0) {
    stop(length(missing_ids), " feature_id values in ", feature_groups_path, " are not in the feature table")
  }
  feature_table <- feature_table[feature_table$feature_id %in% feature_groups$feature_id, ]
  cat(sprintf("Feature set from Step 1 Stage 3: %d features in %d groups (fc_threshold argument not applied)\n",
              nrow(feature_table), length(unique(feature_groups$group))))
} else if (length(fc_col) == 1) {
  feature_table <- feature_table[feature_table[[fc_col]] >= fc_threshold, ]
  cat(sprintf("Filtered by %s >= %d: %d features retained\n", fc_col, fc_threshold, nrow(feature_table)))
} else if (length(fc_col) > 1) {
  warning("Multiple fold-change columns found: ", paste(fc_col, collapse = ", "), ". Skipping filter.")
}

# Identify sample intensity columns from mapfile
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

# Build peak table: feature_id, mz, rt, and sample intensity columns
# Rename 'time' -> 'rt' (new API expects 'rt')
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

# Read compound database from mapfile entry
db_base_dir <- dirname(config$db_mapfile_path)
compound_db_path <- file.path(db_base_dir, db_ls[[ii, "DB_name"]])
cat(sprintf("Compound DB:    %s\n", compound_db_path))
compound_table <- read_xlsx(compound_db_path)

# Ensure expected column names
expected_cols <- c("compound_id", "name", "molecular_formula", "monoisotopic_mass")
if (!all(expected_cols %in% colnames(compound_table))) {
  # Try renaming from common alternatives
  colnames(compound_table)[1:4] <- expected_cols
}

# If the database carries a `charge` column, keep only compounds observable in
# this polarity: neutrals plus permanent ions of the matching sign. Done here so
# the adduct weights built below do not cover the opposite polarity's bare ions.
# No-op for a database without a `charge` column.
compound_table <- filter_compounds_by_polarity(compound_table, charge_type)

cat(sprintf("Compound database: %d compounds (%s mode)\n",
            nrow(compound_table), charge_type))

# Load metadata file for post-processing
db_metafile_path <- file.path(db_base_dir, db_ls[[ii, "Meta_File"]])
cat(sprintf("DB metadata:    %s\n", db_metafile_path))
db_metafile <- read_xlsx(db_metafile_path)
names(db_metafile)[names(db_metafile) == "Formula_ID"] <- "compound_id"


# =============================================================================
# Section 5: Configure Adducts
# =============================================================================

# Load the package's built-in adduct table
adduct_table <- get_default_adduct_table()

# Select adducts based on ionization mode.
# "M" is deliberately absent: it is a zero mass shift, i.e. the neutral mass,
# which is not an observable species. Permanently charged compounds are handled
# separately through the bare-ion pseudo-adducts below.
#
# *** CHANGED 2026-10-22: filter_by is "M+H" / "M-H" again (NULL since 2026-08-14), ***
# *** with weight 5 on M+H / M-H (was 1) and tighter tolerances (Section 6).           ***
#
# `filter_by` names the "primary ion(s)". Setting it ACTIVATES a stricter
# confidence regime (`filter_is_active` in multilevelannotationstep4): the
# Stage 3 correlation check is centred on the primary-ion row, a compound with
# no primary-ion row skips the Stage 4 decision tree (Confidence 0, or 1 when
# its score is >= 10), and the later evidence upgrade gives at most
# Confidence 2. Compounds with the primary ion can still reach Confidence 3.
# The adduct weight multiplies the score by 10^weight of the compound's
# highest-weight adduct, so weight 5 mainly decides which compound ranks first
# on a feature.
#
# Measured on CLU0127 (C18pos, checked on C18neg from the same samples) against
# the CLUES confirmed library as compound database (2,906 compounds), with the
# Step 1 Stage 3 feature groups (5-Test_Output/34-Stage3_Grouping_Eval):
#   - filter_by alone (5 ppm): precision of confident (Confidence >= 2)
#     annotations 0.33 -> 0.51 (C18pos), 0.36 -> 0.55 (C18neg); confident
#     annotations on another confirmed compound's feature 118 -> 5, 66 -> 2.
#   - weight 5 alone: precision 0.34 / 0.35; true compound ranked first at
#     154 / 189 of 208 / 256 confirmed features (baseline 140 / 163).
#   - both with 1 ppm mass / isotope tolerance and 2.5 s: precision 0.62 / 0.55,
#     ranked first 162 / 195. Cost: correct confident annotations 194 -> 139
#     and 223 -> 168.
# Limits: one study, the two modes come from the same samples, and the
# database was the confirmed library, not the production databases.
#
# Results from the two filter_by settings are NOT comparable: re-run anything
# annotated with filter_by = NULL before comparing confidence levels. Set
# filter_by <- NULL and primary_weight <- 1 to restore the 2026-08-14 settings.
primary_weight <- 5
if (charge_type == "pos") {
  selected_adducts <- c("M+H", "M+2H", "M+NH4", "M+Na", "M+ACN+H", "2M+H", "M+H-H2O")
  filter_by <- "M+H"
} else {
  selected_adducts <- c("M-H", "M-H2O-H", "M+Na-2H", "M+Cl", "M+Hac-H", "2M-H")
  filter_by <- "M-H"
}

# Filter adduct table to selected adducts
adduct_table <- adduct_table[adduct_table$adduct %in% selected_adducts, ]

# If the compound database carries a `charge` column, advanced_annotation() adds
# one bare-ion pseudo-adduct per distinct non-zero charge ("[M]+", "[M]+2", ...)
# so those compounds are searched at monoisotopic_mass / |z|. Add the same rows
# here so the weights built below cover them; the call is a no-op otherwise.
adduct_table <- add_bare_ion_adducts(adduct_table, compound_table)

# Build adduct weights from filtered table (primary ions primary_weight, others 1)
adduct_weights <- data.frame(
  adduct = adduct_table$adduct,
  weight = ifelse(adduct_table$adduct %in% c("M+H", "M-H"), primary_weight, 1)
)

cat(sprintf("Using %d adducts for %s mode\n", nrow(adduct_table), charge_type))


# =============================================================================
# Section 6: Set Annotation Parameters
# =============================================================================

# Core matching tolerances
# *** CHANGED 2026-10-22: mass and isotope tolerance 1.5 ppm (were 5 and 7.5), ***
# *** time_tolerance 2.5 s (was 10).                                           ***
# On CLU0009 (Orbitrap) confirmed compounds sat at a median 0.33 ppm and true
# matches were exhausted by 1 ppm, so the 5f annotations used 1e-6 for both
# mass_tolerance and isotope_mass_tolerance, with time_tolerance 2.5. On
# CLU0127 (run 34, with filter_by and weight 5) precision was 0.61 / 0.60
# (C18pos) and 0.55 / 0.54 (C18neg) at 1 / 2 ppm, against 0.53 / 0.52 at 5 ppm.
# Those runs used time_tolerance 10 and isotope tolerance 7.5 ppm. With 2.5 s and
# isotope tolerance equal to mass tolerance (run 35, Step 1 Stage 3 groups),
# 1.5 ppm against 1 ppm: precision 0.607 vs 0.615 (C18pos) and 0.537 vs 0.547
# (C18neg); correct top-1 at the confirmed feature 165 vs 162 of 208 and
# 198 vs 195 of 256. Measure your own instrument's accuracy before changing these.
mass_tolerance      <- 1.5e-6   # 1.5 ppm (fractional, not integer ppm)
time_tolerance      <- 2.5      # seconds
correlation_threshold <- 0.7    # 0.5-0.9 changed CLU0127 precision by <= 0.01 (run 34, 5 ppm, filter_by NULL)

# Isotope parameters
maximum_isotopes           <- 8
mass_defect_precision      <- 0.01
isotope_mass_tolerance     <- 1.5e-6  # 1.5 ppm (was 7.5e-6); defaults to mass_tolerance when NULL
mass_defect_tolerance      <- 0.02    # Da: how far the charge-corrected spacing may
                                      # sit from a whole isotope step. A cheap prefilter;
                                      # isotope_mass_tolerance does the real work.
                                      # Values >= 0.5 disable the check entirely.
intensity_deviation_tolerance <- 0.3
identify_isotopologues_flag   <- TRUE   # use enviPat to label isotopologues (13C, 15N, etc.)

# Peak and clustering parameters (WGCNA clustering; unused when Step 1 feature groups are read)
peak_rt_width     <- 4
deep_split        <- 2
min_cluster_size  <- 10
network_type      <- "unsigned"

# Annotation quality parameters
MplusH_abundance_ratio_check  <- TRUE
multimer_abundance_check      <- FALSE
redundancy_filtering          <- FALSE  # Stage 4b is the canonical final output; set TRUE to also write Stage 5 (one row per (mz, time), destructive on multi-match evidence)
# Report every retention-time cluster a compound matches, each with its own
# confidence, instead of keeping one. Off by default in advanced_annotation();
# on CLU0009 it raised recall of authentic-standard confirmed features from
# 30 to 55 of 72. Output roughly doubles, almost all as Confidence 0, and the
# n_rt_clusters column records how many candidates each compound had.
# Untested with redundancy_filtering = TRUE.
multi_rt_clusters             <- TRUE
clustering_seed               <- 12345  # pins peak-module clustering; NULL leaves the RNG untouched
min_ions_per_chemical         <- 2
level1_primary_adducts       <- c("M+H", "M-H")  # adducts that qualify for Conf 1 as a single match

# Pathway settings
pathway_mode                <- "skip"   # "skip", "HMDB", or "custom"
pathway_data                <- NULL     # required if pathway_mode = "custom"
excluded_pathways           <- NULL     # pathway IDs to exclude from scoring
excluded_pathway_compounds  <- NULL     # compound IDs to exclude from pathway scoring

# Boosted compounds
boosted_compounds    <- NULL            # data.frame with mz (and optionally rt) to boost
boost_match_by       <- c("mz", "rt")   # match boosted compounds by mz, rt, or both
boost_mass_tolerance <- NULL            # defaults to mass_tolerance
boost_time_tolerance <- NULL            # defaults to time_tolerance


# =============================================================================
# Section 7: Run Annotation
# =============================================================================

# Output layout (1.3.0): one step folder per mode, one subfolder per database. Prefix 7- follows the
# Step 1-3 folders (4- MS2 extraction, 5- confirmed compounds, 6- MS2 matching). Before 1.3.0 each
# database had its own 7-<mode>_<study>_<DB>_CLUES_Annotation/ (1.2.2) or 4c-... folder (before 1.2.2).
db_out <- db_ls[[ii, "File_name_out"]]
step_dir <- file.path(sub("/+$", "", config$output_base_loc), lc_mode, paste0("7-", lc_mode, "_", study_id, "_xMSannotator"))
db_dir <- file.path(step_dir, db_out)
outloc <- file.path(db_dir, "intermediate")   # advanced_annotation() writes its Stage files here
dir.create(outloc, recursive = TRUE, showWarnings = FALSE)

# Check if output already exists — skip if so
output_name <- file.path(db_dir, paste0(lc_mode, "-", study_id, "_", db_out, "_CLUES_Annotation.xlsx"))
if (file.exists(output_name)) {
  cat(sprintf("Output already exists, skipping: %s\n", output_name))
  quit(save = "no", status = 0)
}

# --- Run record helpers ---
fmt_setting <- function(x) {
  if (is.null(x)) return("NULL")
  if (is.data.frame(x)) return(paste(trimws(x[[1]]), x[[2]], sep = "=", collapse = ", "))
  if (is.character(x)) return(paste(trimws(x), collapse = ", "))
  paste(format(x, scientific = FALSE, trim = TRUE), collapse = ", ")
}
runscript_path <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
runscript_md5 <- if (length(runscript_path) == 1 && file.exists(runscript_path)) unname(tools::md5sum(runscript_path)) else NA
settings <- list(
  batch_correction = config$batch_correction, blank_names = blank_names,
  feature_groups = if (use_feature_groups) "Step 1 Stage 3 groups (WGCNA not run)" else paste("none: WGCNA; fc_threshold", fc_threshold),
  selected_adducts = selected_adducts, filter_by = filter_by, primary_weight = primary_weight,
  mass_tolerance = mass_tolerance, time_tolerance = time_tolerance, correlation_threshold = correlation_threshold,
  isotope_mass_tolerance = isotope_mass_tolerance, mass_defect_tolerance = mass_defect_tolerance,
  mass_defect_precision = mass_defect_precision, intensity_deviation_tolerance = intensity_deviation_tolerance,
  maximum_isotopes = maximum_isotopes, identify_isotopologues_flag = identify_isotopologues_flag,
  peak_rt_width = peak_rt_width, deep_split = deep_split, min_cluster_size = min_cluster_size, network_type = network_type,
  MplusH_abundance_ratio_check = MplusH_abundance_ratio_check, multimer_abundance_check = multimer_abundance_check,
  redundancy_filtering = redundancy_filtering, multi_rt_clusters = multi_rt_clusters, clustering_seed = clustering_seed,
  min_ions_per_chemical = min_ions_per_chemical, level1_primary_adducts = level1_primary_adducts,
  pathway_mode = pathway_mode, boosted_compounds = if (is.null(boosted_compounds)) NULL else paste(nrow(boosted_compounds), "rows"),
  n_workers = n_workers)
settings_lines <- sprintf("| %s | %s |", names(settings), vapply(settings, fmt_setting, ""))
# adduct_weights depends on the database (bare-ion pseudo-adducts [M]+ ... for charged compounds): run_info.txt only
db_settings <- list(adduct_weights = adduct_weights)

# --- Step-folder README: same content from every task of this mode (no task-specific values or times);
#     written to a temporary file and renamed, so tasks running at the same time never leave a partial file ---
readme_path <- file.path(step_dir, paste0(lc_mode, "-", study_id, "_xMSannotator_README.md"))
final_stage <- if (redundancy_filtering) "Stage5_curated_results" else "Stage4b_confidence_levels"
readme <- c(
  sprintf("# xMSannotator annotation: %s %s", study_id, lc_mode), "",
  "The Step 1 feature table annotated against each compound database of the mapfile, one database per SLURM array task",
  "(`run_CLUES_xMSannotator.sh`, `xMSannotator_CLUES_MultiDB_Runscript.R`). Every database is annotated on its own:",
  "scores and confidence levels compare candidates within one database, not across databases. Nothing is merged.", "",
  "## Layout", "```",
  sprintf("7-%s_%s_xMSannotator/", lc_mode, study_id),
  sprintf("├── %s", basename(readme_path)),
  "└── <File_name_out>/                     one folder per database (table below)",
  sprintf("    ├── %s-%s_<File_name_out>_CLUES_Annotation.xlsx   the result", lc_mode, study_id),
  "    ├── run_info.txt                     inputs and md5, versions, settings, counts, run time of this database",
  "    └── intermediate/                    advanced_annotation() Stage files, gzipped after the xlsx is written", "```", "",
  "A database folder with no xlsx, or with uncompressed `intermediate/*.txt`, is from a task that did not finish.", "",
  "## Files",
  sprintf("- **xlsx**: `%s` (one row per compound x adduct or isotope x feature) joined with the database METADATA file", final_stage),
  "  (by `compound_id`), the feature table QA columns (by `feature_id`) and the Step 1 `feature_group`. Column reference:",
  "  `docs/Stage5_Output_Reference.md` and `docs/xMSannotator_Workflow.md` of CLUES-xMSannotator.",
  "- **intermediate/** (tab-separated, gzipped):",
  "  `Stage1_mass_matched` (mass matches), `Stage1_peak_clusters` (feature groups or WGCNA modules; same for every database),",
  "  `Stage2_isotope_detection`, `Stage3_chemical_scores`, `Stage3_pathway_*`, `Stage4_confidence_levels`,",
  "  `Stage4a_confidence_levels` (all rows, including incoherent), `Stage4b_confidence_levels` (coherent rows: the xlsx),",
  "  `Stage5_curated_results` (only with redundancy_filtering = TRUE).",
  "  Fields are quoted; a quote inside a field is written as `\\\"`. `data.table::fread()` reads them; `read.table()` needs",
  "  `\\\"` turned into `\"\"` first (as the runscript does).", "",
  "## Databases", sprintf("Mapfile: `%s` (md5 %s)", config$db_mapfile_path, unname(tools::md5sum(config$db_mapfile_path))), "",
  "| Number | File_name_out | DB_name | Meta_File |", "|---|---|---|---|",
  sprintf("| %s | %s | %s | %s |", db_ls$Number, db_ls$File_name_out, db_ls$DB_name, db_ls$Meta_File), "",
  "## Inputs",
  sprintf("- Feature table: `%s`", feature_table_path),
  sprintf("- Feature groups: `%s`", if (use_feature_groups) feature_groups_path else "none"),
  sprintf("- Sample mapfile: `%s` (blank sample types removed: %s)", mapfile_path, paste(blank_names, collapse = ", ")), "",
  "## Versions",
  sprintf("- CLUES.xMSannotator %s; %s", as.character(packageVersion("CLUES.xMSannotator")), R.version.string),
  sprintf("- Runscript: `%s` (md5 %s)", if (length(runscript_path) == 1) runscript_path else "unknown", runscript_md5), "",
  sprintf("## Settings (%s mode)", charge_type),
  "Adduct weights: primary_weight on M+H / M-H, 1 on the other adducts; databases with charged compounds also get",
  "bare-ion adducts ([M]+, [M]-, ...) with weight 1. Each database's full list is in its run_info.txt.", "",
  "| setting | value |", "|---|---|", settings_lines)
readme_tmp <- paste0(readme_path, ".tmp", Sys.getpid())
writeLines(readme, readme_tmp)
if (!file.rename(readme_tmp, readme_path)) stop("Could not write ", readme_path)

cat("\n")
cat("===========================================================\n")
cat(sprintf(" Starting annotation: %s | %s | %s\n", study_id, lc_mode, charge_type))
cat(sprintf(" Database: %s\n", db_ls[[ii, "DB_name"]]))
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
    level1_primary_adducts        = level1_primary_adducts,
    charge_type                   = charge_type,
    multi_rt_clusters             = multi_rt_clusters,
    clustering_seed               = clustering_seed,
    pathway_mode                  = pathway_mode,
    pathway_data                  = pathway_data,
    excluded_pathways             = excluded_pathways,
    excluded_pathway_compounds    = excluded_pathway_compounds,
    boosted_compounds             = boosted_compounds,
    boost_match_by                = boost_match_by,
    boost_mass_tolerance          = boost_mass_tolerance,
    boost_time_tolerance          = boost_time_tolerance,
    feature_id_column             = feature_id_column,
    feature_groups                = feature_groups,
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
  results_file <- file.path(outloc, "Stage4b_confidence_levels.txt")
}

# The Stage files are written by write.table(), which escapes a double quote inside a field as \".
# read.table() does not read that form back: one compound name with a quote (NormanNetwork,
# 'Myricitrin 2"-O-gallate') merged 12,174 rows into one field (1.2.3). Turning \" into the doubled
# quote "" that read.table accepts reads every row; files without such a quote read as before.
results <- read.table(text = gsub('\\"', '""', readLines(results_file), fixed = TRUE), sep = "\t", header = TRUE)

# Join with compound metadata
results <- left_join(results, db_metafile, by = "compound_id")

# Merge feature table QA metadata. Names from Step 1 v1.6.0+ and from earlier
# versions; intersect() below keeps those present in this table.
metadata_cols <- c("feature_id", "rtmin", "rtmax", "n_detected", "n_filled",
                   "rt_range_10_90", "mz_range_ppm", "median_snr_asari",
                   "frac_no_baseline", "median_dppm", "median_width", "median_n_scans",
                   "Water_fc", "Water_status", "Overall_detected_pct",
                   "Overall_nonzero_pct", "Study_Sample_detected_pct",
                   # before Step 1 v1.6.0
                   "peak_quality", "alignment_score", "combined_score", "quality_category",
                   "peak_shape", "passed_filter", "Overall_detect_pct", "Study_Sample_detect_pct")
# Also include any _CV columns (names vary by study, e.g., HRE.p1Std._CV, Study_Sample_CV)
cv_cols <- grep("_CV$", colnames(feature_table), value = TRUE, ignore.case = TRUE)
metadata_cols <- unique(c(metadata_cols, cv_cols))
feature_metadata <- feature_table[, intersect(metadata_cols, colnames(feature_table))]
feature_metadata <- unique(feature_metadata)
results <- left_join(results, feature_metadata, by = "feature_id")
if (!is.null(feature_groups)) {
  results <- left_join(results, rename(feature_groups, feature_group = group), by = "feature_id")
}

# Export to Excel
write_xlsx(results, output_name)
cat(sprintf("\nResults exported to: %s\n", output_name))
cat(sprintf("Total annotations: %d\n", nrow(results)))


# =============================================================================
# Section 9: Compress intermediates, write run_info.txt
# =============================================================================

# Row count of Stage 1 before compression (header excluded)
stage1_file <- file.path(outloc, "Stage1_mass_matched.txt")
n_stage1 <- if (file.exists(stage1_file)) length(readLines(stage1_file)) - 1L else NA_integer_

# gzip each Stage file; the original is removed only after the decompressed copy is byte-identical
gzip_verified <- function(f) {
  raw_in <- readBin(f, "raw", file.size(f))
  gz <- paste0(f, ".gz")
  con <- gzfile(gz, "wb"); writeBin(raw_in, con); close(con)
  con <- gzfile(gz, "rb"); raw_back <- readBin(con, "raw", length(raw_in) + 1L); close(con)
  if (!identical(raw_in, raw_back)) stop("gzip check failed for ", f, "; uncompressed file kept")
  file.remove(f)
  gz
}
stage_files <- list.files(outloc, pattern = "^Stage.*\\.txt$", full.names = TRUE)
gz_files <- vapply(stage_files, gzip_verified, "")
cat(sprintf("Compressed %d Stage files in %s\n", length(gz_files), outloc))

conf_tab <- table(factor(results$Confidence, levels = 0:4))
run_info <- c(
  sprintf("run_info: xMSannotator, %s %s, database %s", study_id, lc_mode, db_out),
  sprintf("written: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("run_time_min: %.1f", as.numeric(difftime(Sys.time(), run_start, units = "mins"))),
  "",
  "[database]",
  sprintf("array_index: %s (mapfile Number %s, %s mode)", config$array_index, db_ls[[ii, "Number"]], charge_type),
  sprintf("File_name_out: %s", db_out),
  sprintf("DB_name: %s (md5 %s)", compound_db_path, unname(tools::md5sum(compound_db_path))),
  sprintf("Meta_File: %s (md5 %s)", db_metafile_path, unname(tools::md5sum(db_metafile_path))),
  sprintf("mapfile: %s (md5 %s)", config$db_mapfile_path, unname(tools::md5sum(config$db_mapfile_path))),
  sprintf("compounds after polarity filter: %d", nrow(compound_table)),
  "",
  "[inputs]",
  sprintf("feature_table: %s", feature_table_path),
  sprintf("feature_groups: %s", if (use_feature_groups) feature_groups_path else "none"),
  sprintf("sample_mapfile: %s", mapfile_path),
  sprintf("features: %d | samples: %d", nrow(peak_table), length(sample_cols)),
  "",
  "[versions]",
  sprintf("CLUES.xMSannotator: %s", as.character(packageVersion("CLUES.xMSannotator"))),
  sprintf("R: %s", R.version.string),
  sprintf("runscript: %s (md5 %s)", if (length(runscript_path) == 1) runscript_path else "unknown", runscript_md5),
  "",
  "[settings]",
  sprintf("%s: %s", names(settings), vapply(settings, fmt_setting, "")),
  sprintf("%s: %s", names(db_settings), vapply(db_settings, fmt_setting, "")),
  "",
  "[counts]",
  sprintf("stage1_mass_matches: %s", n_stage1),
  sprintf("xlsx_rows (%s): %d", final_stage, nrow(results)),
  sprintf("annotated_features: %d | compounds: %d", length(unique(results$feature_id)), length(unique(results$compound_id))),
  sprintf("rows_by_confidence: %s", paste(names(conf_tab), conf_tab, sep = "=", collapse = " ")),
  "",
  "[files]",
  sprintf("%s (md5 %s)", basename(output_name), unname(tools::md5sum(output_name))),
  sprintf("intermediate/%s", basename(gz_files)))
writeLines(run_info, file.path(db_dir, "run_info.txt"))
cat(sprintf("Run record: %s\n", file.path(db_dir, "run_info.txt")))
