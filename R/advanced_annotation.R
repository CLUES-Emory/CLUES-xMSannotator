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
    filter(!duplicated(.data$chemical_ID)) %>%
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
                                mass_defect_tolerance = 0.1,
                                mass_defect_precision = 0.01,
                                time_tolerance = 10,
                                peak_rt_width = 1,
                                correlation_threshold = 0.7,
                                deep_split = 2,
                                min_cluster_size = 10,
                                maximum_isotopes = 10,
                                min_ions_per_chemical = 2,
                                filter_by = c("M-H", "M+H"),
                                network_type = "unsigned",
                                redundancy_filtering = TRUE,
                                outloc = tempdir(),
                                n_workers = parallel::detectCores()) {
  if (is.null(adduct_table)) {
    adduct_table <- sample_adduct_table
  }

  if (is.null(adduct_weights)) {
    adduct_weights <- as.data.frame(tibble(adduct = adduct_table$adduct, weight = rep_len(5, length(adduct_table$adduct))))
  }

  if (is.numeric(n_workers) && n_workers > 1) {
    WGCNA::allowWGCNAThreads(n_workers)
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
    rt_tolerance = time_tolerance
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
  stage2_output <- annotation
  if (!is.null(mz_rt_feature_id_map)) {
    stage2_output <- left_join(stage2_output, mz_rt_feature_id_map, by = c("mz", "time"))
  }
  write.table(stage2_output, file = file.path(outloc, "Stage2_isotope_detection.txt"),
              sep = "\t", row.names = FALSE)
  # ----------------------------

  # Tool 8: Compute chemscores
  # ----------------------------
  # Get unique chemical_IDs and process each once.
  # NOTE: distinct() only controls how many times get_chemscore is CALLED.
  # Inside get_chemscore, it queries the FULL annotation table for all rows
  # with that chemical_ID, so ALL adducts are preserved.
  unique_chemicals <- annotation %>%
    dplyr::distinct(chemical_ID, .keep_all = TRUE)

  annotation <- purrr::pmap_dfr(
    unique_chemicals,
    ~ get_chemscore(...,
                    annotation = annotation,
                    adduct_weights = adduct_weights,
                    corthresh = correlation_threshold,
                    global_cor = global_cor,
                    max_diff_rt = time_tolerance,
                    filter.by = filter_by
    )
  )
  # ----------------------------

  # Tool 9: pathway matching
  # ----------------------------
  data(hmdbCompMZ)
  chemCompMZ <- dplyr::rename(hmdbCompMZ, chemical_ID = HMDBID)

  annotation <- multilevelannotationstep3(
    chemCompMZ = chemCompMZ,
    chemscoremat = annotation,
    adduct_weights = adduct_weights,
    db_name = "HMDB",
    max_diff_rt = time_tolerance,
    pathwaycheckmode = "pm",
    outloc = outloc,
    mz_rt_feature_id_map = mz_rt_feature_id_map
  )
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
    mz_rt_feature_id_map = mz_rt_feature_id_map
  )
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
  }

  # Join feature ID column to final annotation
  if (!is.null(mz_rt_feature_id_map)) {
    annotation <- left_join(annotation, mz_rt_feature_id_map, by = c("mz", "time"))
    # Re-write Stage5 output with feature ID column included
    if (redundancy_filtering) {
      write.table(annotation, file = file.path(outloc, "Stage5_curated_results.txt"), sep = "\t", row.names = FALSE)
    }
  }

  annotation
}
