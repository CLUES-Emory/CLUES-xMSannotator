create_row <- function(adduct, query) {
  query$Adduct <- adduct
  return(query)
}

append_adduct <- function(..., adduct_names) {
  query <- tibble(...)
  rows_with_adduct <- lapply(
    adduct_names,
    create_row,
    query = query
  )
  return(rows_with_adduct)
}

create_chemCompMZ <- function(database, adduct_names) {
  database <- purrr::pmap_dfr(
    database,
    ~ append_adduct(...,
      adduct_names = adduct_names
    )
  )
  return(database)
}

#' @import dplyr
#' @importFrom rlang .data
compute_mass_defect <- function(peaks, precision) {
  mutate(peaks, mass_defect = cut(.data$mz %% 1, seq(0, 1, precision), labels = FALSE))
}

#' @import dplyr
#' @importFrom rlang .data
remove_duplicates <- function(annotation, adduct_weights) {
  is_max <- function(x) x == max(x)
  is_unique <- function(score, adduct) is_max(score * 100^is.element(adduct, adduct_weights$adduct))
  recompute_score <- function(mask, score) score * (sum(mask) / length(mask)) # FIXME: check if decreasing the score is OK

  annotation <- group_by(annotation, .data$mz)
  annotation <- mutate(annotation, unique = is_unique(.data$score, .data$adduct))
  annotation <- group_by(annotation, .data$compound)
  annotation <- mutate(annotation, score = recompute_score(.data$unique, .data$score))
  annotation <- ungroup(annotation)
  annotation <- filter(annotation, .data$unique)
  annotation <- select(annotation, -.data$unique, -.data$multiple_match)
  annotation
}

#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang .data
print_confidence_distribution <- function(annotation) {
  confidence_distribution_across_compounds <- annotation %>%
    filter(!duplicated(.data$compound_id)) %>%
    count(.data$Confidence)

  confidence_distribution_across_formulas <- annotation %>%
    filter(!duplicated(.data$Formula)) %>%
    count(.data$Confidence)

  print("Confidence level distribution for unique compounds")
  print(confidence_distribution_across_compounds)
  print("Confidence level distribution for unique formulas")
  print(confidence_distribution_across_formulas)

  invisible(annotation)
}

#' Safely join feature_id column only if it doesn't already exist
#' @description Prevents duplicate feature_id.x and feature_id.y columns
#'   when feature_id has already been joined by internal functions.
#' @param data Data frame to join to
#' @param feature_id_map Mapping table with mz, time, and feature_id columns
#' @param feature_id_column Name of the feature_id column
#' @return Data frame with feature_id column (if not already present)
#' @import dplyr
safe_join_feature_id <- function(data, feature_id_map, feature_id_column) {
  if (is.null(feature_id_map)) return(data)
  if (is.null(feature_id_column)) return(data)
  if (feature_id_column %in% names(data)) return(data)  # Already has it
  left_join(data, feature_id_map, by = c("mz", "time"))
}

#' Skip pathway matching step
#' @description Minimal transformation for users who want to skip pathway matching.
#'   Renames cur_chem_score to score (required by downstream functions) and writes output.
#' @param chemscoremat Chemical score matrix from get_chemscore
#' @param outloc Output directory
#' @param mz_rt_feature_id_map Optional feature ID mapping
#' @return chemscoremat with score column renamed
skip_pathway_step <- function(chemscoremat, outloc, mz_rt_feature_id_map = NULL) {
  # Rename cur_chem_score to score (required by downstream functions)
  names(chemscoremat)[names(chemscoremat) == "cur_chem_score"] <- "score"

  # Join feature ID if mapping provided (consistent with multilevelannotationstep3)
  if (!is.null(mz_rt_feature_id_map)) {
    chemscoremat <- dplyr::left_join(chemscoremat, mz_rt_feature_id_map, by = c("mz", "time"))
  }

  write.table(chemscoremat, file = file.path(outloc, "Stage3_pathway_skipped.txt"),
              sep = "\t", row.names = FALSE)

  return(chemscoremat)
}

#' Wrapper for the advanced annotation steps.
#' @export
#' @import dplyr
#' @importFrom magrittr %>%
advanced_annotation <- function(peak_table,
                                compound_table,
                                adduct_table = NULL,
                                adduct_weights = NULL,
                                feature_id_column = NULL,
                                intensity_deviation_tolerance = 0.1,
                                mass_tolerance = 5e-6,
                                isotope_mass_tolerance = NULL,
                                mass_defect_tolerance = 0.1,
                                mass_defect_precision = 0.01,
                                time_tolerance = 10,
                                peak_rt_width = 1,
                                correlation_threshold = 0.7,
                                MplusH_abundance_ratio_check = TRUE,
                                multimer_abundance_check = TRUE,
                                deep_split = 2,
                                min_cluster_size = 10,
                                maximum_isotopes = 10,
                                min_ions_per_chemical = 2,
                                filter_by = c("M-H", "M+H"),
                                network_type = "unsigned",
                                redundancy_filtering = TRUE,
                                pathway_mode = "HMDB",
                                pathway_data = NULL,
                                excluded_pathways = NULL,
                                excluded_pathway_compounds = NULL,
                                boosted_compounds = NULL,
                                boost_match_by = c("mz", "rt"),
                                boost_mass_tolerance = NULL,
                                boost_time_tolerance = NULL,
                                enable_permutation = FALSE,
                                n_permutations = 1000,
                                permutation_method = "full",
                                permutation_seed = 42,
                                outloc = tempdir(),
                                n_workers = parallel::detectCores()) {
  if (is.null(adduct_table)) {
    adduct_table <- sample_adduct_table
  }

  if (is.null(adduct_weights)) {
    adduct_weights <- as.data.frame(tibble(adduct = adduct_table$adduct, weight = rep_len(5, length(adduct_table$adduct))))
  }

  # Validate pathway_mode
  if (!pathway_mode %in% c("HMDB", "custom", "skip")) {
    stop("pathway_mode must be one of: 'HMDB', 'custom', 'skip'")
  }

  if (pathway_mode == "custom" && is.null(pathway_data)) {
    stop("pathway_data is required when pathway_mode = 'custom'")
  }

  # Default boost tolerances to main tolerances if not specified
  if (is.null(boost_mass_tolerance)) boost_mass_tolerance <- mass_tolerance
  if (is.null(boost_time_tolerance)) boost_time_tolerance <- time_tolerance

  # Validate and process boosted_compounds
  if (!is.null(boosted_compounds)) {
    valid_match <- c("mz", "rt")
    if (!all(boost_match_by %in% valid_match)) {
      stop("boost_match_by must be c('mz'), c('rt'), or c('mz', 'rt')")
    }
    boosted_compounds <- as_boosted_compounds_table(boosted_compounds, boost_match_by)
  }

  if (is.numeric(n_workers) && n_workers > 1) {
    WGCNA::allowWGCNAThreads(n_workers)
  }

  # Default isotope_mass_tolerance to mass_tolerance if not specified
  # Convert to ppm (mass_tolerance is in fractional form, e.g., 5e-6 = 5 ppm)
  if (is.null(isotope_mass_tolerance)) {
    isotope_mass_tolerance_ppm <- mass_tolerance * 1e6
  } else {
    isotope_mass_tolerance_ppm <- isotope_mass_tolerance * 1e6
  }

  # Store original peak table to preserve user's feature ID column
  peak_table_orig <- peak_table

  # Remove feature_id_column before validation (it's non-numeric and would fail validation)
  if (!is.null(feature_id_column) && feature_id_column %in% colnames(peak_table)) {
    peak_table <- peak_table[, colnames(peak_table) != feature_id_column, drop = FALSE]
  }

  peak_table <- as_peak_table(peak_table, intensities = TRUE)
  adduct_table <- as_adduct_table(adduct_table)
  compound_table <- as_compound_table(compound_table)

  # Create feature ID mapping if user specified a column
  feature_id_map <- NULL
  mz_rt_feature_id_map <- NULL
  if (!is.null(feature_id_column)) {
    if (!feature_id_column %in% colnames(peak_table_orig)) {
      warning(paste("feature_id_column", feature_id_column, "not found in peak_table"))
    } else {
      # Map by peak (for Stage 1 outputs before reformat)
      feature_id_map <- tibble(
        peak = peak_table$peak,
        !!feature_id_column := peak_table_orig[[feature_id_column]]
      )
      # Map by mz + time (for downstream stages where peak column is lost)
      mz_rt_feature_id_map <- tibble(
        mz = peak_table$mz,
        time = peak_table$rt,
        !!feature_id_column := peak_table_orig[[feature_id_column]]
      )
    }
  }

  # Create output directory if it doesn't exist
  if (!dir.exists(outloc)) {
    dir.create(outloc, recursive = TRUE)
  }

  # Tool 1: Simple annotation
  # ---------------------------
  annotation <- simple_annotation(
    peak_table = peak_table,
    compound_table = compound_table,
    adduct_table = adduct_table,
    mass_tolerance = mass_tolerance
  )
  # ---------------------------

  # ----------------------------
  # Summary: Adduct detection results
  # ----------------------------
  n_annotations <- nrow(annotation)
  n_unique_peaks <- length(unique(annotation$peak))
  n_unique_compounds <- length(unique(annotation$compound))

  cat("\n=== Adduct Detection Summary ===\n")
  cat(sprintf("Total annotations: %d\n", n_annotations))
  cat(sprintf("Unique peaks matched: %d\n", n_unique_peaks))
  cat(sprintf("Unique compounds matched: %d\n", n_unique_compounds))

  adduct_breakdown <- annotation %>%
    count(adduct) %>%
    arrange(desc(n))

  cat("\nAdduct breakdown:\n")
  print(adduct_breakdown, n = 20)
  cat("=================================\n\n")
  # ----------------------------

  # Output Stage1: Simple annotation results (mass matching)
  # ----------------------------
  stage1_annotation <- annotation
  if (!is.null(feature_id_map)) {
    stage1_annotation <- left_join(stage1_annotation, feature_id_map, by = "peak")
  }
  write.table(stage1_annotation, file = file.path(outloc, "Stage1_mass_matched.txt"), sep = "\t", row.names = FALSE)
  # ----------------------------

  # Tool 2: Compute mass defect
  # ----------------------------
  annotation <- compute_mass_defect(annotation, precision = mass_defect_precision)
  # ----------------------------


  # Tool 3: Compute correlations (peak intensity matrix will be the input, so doesn't need to be cut)
  # also check the galaxy tool to compute correlations
  # ----------------------------
  peak_intensity_matrix <- get_peak_intensity_matrix(peak_table)
  peak_correlation_matrix <- compute_peak_correlations(peak_intensity_matrix, correlation_method = "p")
  # ----------------------------

  # Tool 4: Compute peak modules
  # ----------------------------
  peak_modules <- compute_peak_modules(
    peak_intensity_matrix = peak_intensity_matrix,
    peak_correlation_matrix = peak_correlation_matrix,
    correlation_threshold = correlation_threshold,
    deep_split = deep_split,
    min_cluster_size = min_cluster_size,
    network_type = network_type
  )
  # ----------------------------

  # Tool 5: Compute rt modules
  # ----------------------------
  peak_rt_clusters <- compute_rt_modules(
    peak_table = inner_join(peak_table, peak_modules, by = "peak"),
    peak_width = peak_rt_width
  )
  
  # Step can be done in Galaxy using the cut and join tools
  # ----------------------------
  annotation <- inner_join(annotation,
    select(peak_rt_clusters, "peak", "mean_intensity", "module", "rt_cluster"),
    by = "peak"
  )
  # ----------------------------

  # Re-use the mass defect tool on the output of compute rt modules
  # ----------------------------
  peak_table <- peak_table %>%
    select(peak, mz, rt) %>%
    inner_join(peak_rt_clusters, by = "peak") %>%
    compute_mass_defect(precision = mass_defect_precision)
  # ----------------------------

  # Output Stage1: Peak clustering results
  # ----------------------------
  stage1_output <- peak_table %>%
    mutate(Module_RTclust = paste(module, rt_cluster, sep = "_")) %>%
    select(peak, mz, rt, mean_intensity, module, rt_cluster, Module_RTclust, mass_defect)

  if (!is.null(feature_id_map)) {
    stage1_output <- left_join(stage1_output, feature_id_map, by = "peak")
  }
  write.table(stage1_output, file = file.path(outloc, "Stage1_peak_clusters.txt"), sep = "\t", row.names = FALSE)
  # ----------------------------

  # Tool 6: Compute isotopes
  # ----------------------------
  annotation <- compute_isotopes(
    annotation = annotation,
    adduct_weights = adduct_weights,
    intensity_deviation_tolerance = intensity_deviation_tolerance,
    mass_defect_tolerance = mass_defect_tolerance,
    peak_table = peak_table,
    rt_tolerance = time_tolerance,
    isotope_mass_tolerance_ppm = isotope_mass_tolerance_ppm
  )
  # ----------------------------

  # ----------------------------
  # Summary: Isotope detection results
  # ----------------------------
  n_monoisotopic <- sum(annotation$mass_number_difference == 0)
  n_isotopes <- sum(annotation$mass_number_difference > 0)

  cat("\n=== Isotope Detection Summary ===\n")
  cat(sprintf("Monoisotopic peaks: %d\n", n_monoisotopic))
  cat(sprintf("Isotopes detected:  %d (%.1f%%)\n", n_isotopes, 100 * n_isotopes / nrow(annotation)))

  if (n_isotopes > 0) {
    isotope_breakdown <- annotation %>%
      filter(mass_number_difference > 0) %>%
      count(adduct, mass_number_difference) %>%
      arrange(adduct, mass_number_difference)

    cat("\nIsotope breakdown by adduct:\n")
    print(isotope_breakdown, n = 30)
  }
  cat("=================================\n\n")
  # ----------------------------

  # Tool 7: Reformat annotation and correlation matrix for old xmsannotator
  # ----------------------------
  annotation <- reformat_annotation_table(annotation)
  global_cor <- reformat_correlation_matrix(peak_table, peak_correlation_matrix)
  # ----------------------------

  # Output Stage2: Isotope detection results
  # ----------------------------
  stage2_output <- safe_join_feature_id(annotation, mz_rt_feature_id_map, feature_id_column)
  write.table(stage2_output, file = file.path(outloc, "Stage2_isotope_detection.txt"),
              sep = "\t", row.names = FALSE)
  # ----------------------------

  # Tool 8: Compute chemscores
  # ----------------------------
  # Get unique compound_ids and process each once.
  # NOTE: distinct() only controls how many times get_chemscore is CALLED.
  # Inside get_chemscore, it queries the FULL annotation table for all rows
  # with that compound_id, so ALL adducts are preserved.
  unique_chemicals <- annotation %>%
    dplyr::distinct(compound_id, .keep_all = TRUE)

  annotation <- purrr::pmap_dfr(
    unique_chemicals,
    ~ get_chemscore(...,
                    annotation = annotation,
                    adduct_weights = adduct_weights,
                    corthresh = correlation_threshold,
                    MplusH.abundance.ratio.check = MplusH_abundance_ratio_check,
                    global_cor = global_cor,
                    max_diff_rt = time_tolerance,
                    filter.by = filter_by,
                    adduct_table = adduct_table
    )
  )
  # ----------------------------

  # Output Stage3: Chemical score results
  # ----------------------------
  stage3_output <- safe_join_feature_id(annotation, mz_rt_feature_id_map, feature_id_column)
  write.table(stage3_output, file = file.path(outloc, "Stage3_chemical_scores.txt"),
              sep = "\t", row.names = FALSE)
  # ----------------------------

  # Tool 9: pathway matching
  # ----------------------------
  if (pathway_mode == "skip") {
    annotation <- skip_pathway_step(
      chemscoremat = annotation,
      outloc = outloc,
      mz_rt_feature_id_map = mz_rt_feature_id_map
    )
  } else if (pathway_mode == "custom") {
    annotation <- multilevelannotationstep3(
      chemscoremat = annotation,
      adduct_weights = adduct_weights,
      db_name = "custom",
      max_diff_rt = time_tolerance,
      pathwaycheckmode = "pm",
      outloc = outloc,
      mz_rt_feature_id_map = mz_rt_feature_id_map,
      pathway_data = pathway_data,
      excluded_pathways = excluded_pathways,
      excluded_pathway_compounds = excluded_pathway_compounds
    )
  } else {
    # Default: HMDB pathway matching (existing behavior)
    data(hmdbCompMZ)
    chemCompMZ <- dplyr::rename(hmdbCompMZ, compound_id = HMDBID)

    annotation <- multilevelannotationstep3(
      chemscoremat = annotation,
      adduct_weights = adduct_weights,
      db_name = "HMDB",
      max_diff_rt = time_tolerance,
      pathwaycheckmode = "pm",
      outloc = outloc,
      mz_rt_feature_id_map = mz_rt_feature_id_map,
      chemCompMZ = chemCompMZ
    )
  }
  # ----------------------------

  # Tool 10: compute confidence levels
  # ----------------------------
  annotation <- multilevelannotationstep4(
    outloc = outloc,
    chemscoremat = annotation,
    max.mz.diff = mass_tolerance,
    max.rt.diff = time_tolerance,
    filter.by = filter_by,
    adduct_weights = adduct_weights,
    max_isp = maximum_isotopes,
    min_ions_perchem = min_ions_per_chemical,
    mz_rt_feature_id_map = mz_rt_feature_id_map,
    boostIDs = if (!is.null(boosted_compounds)) boosted_compounds else NA,
    boost.mz.diff = boost_mass_tolerance,
    boost.rt.diff = boost_time_tolerance,
    adduct_table = adduct_table,
    multimer_abundance_check = multimer_abundance_check
  )
  # ----------------------------

  # Output Stage4: Confidence level results
  # ----------------------------
  stage4_output <- safe_join_feature_id(annotation, mz_rt_feature_id_map, feature_id_column)
  write.table(stage4_output, file = file.path(outloc, "Stage4_confidence_levels.txt"),
              sep = "\t", row.names = FALSE)
  # ----------------------------

  # ----------------------------
  # Permutation-based significance testing (optional)
  # Note: This feature is in developmentand not ready for use. Do not use.
  # ----------------------------
  if (FALSE) {
  #if (enable_permutation) {
    perm_start_time <- Sys.time()
    annotation <- compute_permutation_pvalues(
      annotation = annotation,
      peak_table = peak_table,
      compound_table = compound_table,
      adduct_table = adduct_table,
      adduct_weights = adduct_weights,
      mass_tolerance = mass_tolerance,
      time_tolerance = time_tolerance,
      intensity_deviation_tolerance = intensity_deviation_tolerance,
      mass_defect_tolerance = mass_defect_tolerance,
      isotope_mass_tolerance_ppm = isotope_mass_tolerance_ppm,
      correlation_threshold = correlation_threshold,
      filter_by = filter_by,
      peak_correlation_matrix = peak_correlation_matrix,
      n_permutations = n_permutations,
      seed = permutation_seed,
      n_cores = n_workers,
      method = permutation_method
    )

    # Output Stage4 with permutation p-values
    stage4_perm_output <- safe_join_feature_id(annotation, mz_rt_feature_id_map, feature_id_column)
    write.table(stage4_perm_output,
                file = file.path(outloc, "Stage4_permutation_pvalues_multi.txt"),
                sep = "\t", row.names = FALSE)

    perm_end_time <- Sys.time()
    message(sprintf("Permutation testing: %.2f minutes", difftime(perm_end_time, perm_start_time, units = "mins")))
  }
  # ----------------------------     

  # Tool 11: print confidence distribution
  # ----------------------------
  print_confidence_distribution(annotation)
  # ----------------------------

  if (redundancy_filtering) {
    # Tool 12: redundancy filtering
    # ----------------------------
    annotation <- multilevelannotationstep5(
      outloc = outloc,
      adduct_weights = adduct_weights,
      chemscoremat = annotation
    )
    # ----------------------------

    # Re-use the confidence distrivution printing tool
    # ----------------------------
    print_confidence_distribution(annotation)
    # ----------------------------

    # Output Stage5: Curated results after redundancy filtering
    # ----------------------------
    stage5_output <- safe_join_feature_id(annotation, mz_rt_feature_id_map, feature_id_column)
    write.table(stage5_output, file = file.path(outloc, "Stage5_curated_results.txt"),
                sep = "\t", row.names = FALSE)
    # ----------------------------
  }

  # Join feature ID column to final annotation (for return value)
  annotation <- safe_join_feature_id(annotation, mz_rt_feature_id_map, feature_id_column)

  annotation
}
