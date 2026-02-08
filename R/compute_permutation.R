#' Precompute isotope patterns for all unique molecular formulas
#'
#' @description
#' Computes isotopic patterns for all unique molecular formulas in a compound table.
#' This cache can be reused across multiple permutations, eliminating redundant
#' calls to rcdk::get.isotopes.pattern().
#'
#' @param compound_table Data frame with molecular_formula column
#' @param minAbund Minimum abundance threshold for isotopes (default 0.001)
#' @return Named list of isotope patterns, keyed by molecular formula
#' @keywords internal
#' @importFrom rcdk get.formula get.isotopes.pattern
#' @import dplyr
precompute_isotope_patterns <- function(compound_table, minAbund = 0.001) {
  unique_formulas <- unique(compound_table$molecular_formula)
  unique_formulas <- unique_formulas[!is.na(unique_formulas) & unique_formulas != ""]

  message(sprintf("  Precomputing isotope patterns for %d unique formulas...", length(unique_formulas)))

  patterns <- lapply(unique_formulas, function(f) {
    tryCatch({
      formula_obj <- rcdk::get.formula(f)
      isotopes <- rcdk::get.isotopes.pattern(formula_obj, minAbund)
      isotopes <- dplyr::as_tibble(isotopes)
      isotopes <- dplyr::arrange(isotopes, dplyr::desc(abund))
      isotopes <- dplyr::mutate(isotopes,
        mass_number_difference = round(mass - mass[1]),
        exact_mass_diff = mass - mass[1]
      )
      isotopes
    }, error = function(e) {
      NULL
    })
  })
  names(patterns) <- unique_formulas

  n_successful <- sum(!sapply(patterns, is.null))
  message(sprintf("  Cached %d/%d isotope patterns", n_successful, length(unique_formulas)))

  patterns
}

#' Detect isotopic peaks using precomputed pattern cache
#'
#' @description
#' Same as detect_isotopic_peaks but uses a precomputed pattern cache
#' instead of computing patterns on-the-fly.
#'
#' @param ... A single annotated feature passed as a list of columns
#' @param intensity_deviation_tolerance Intensity deviation tolerance
#' @param peaks Peak table
#' @param mass_defect_tolerance Mass defect tolerance
#' @param rt_tolerance RT tolerance
#' @param isotope_mass_tolerance_ppm Mass tolerance in ppm
#' @param pattern_cache Named list of precomputed isotope patterns
#' @return Table of identified isotopes
#' @keywords internal
#' @import dplyr
detect_isotopic_peaks_cached <- function(...,
                                         intensity_deviation_tolerance,
                                         peaks,
                                         mass_defect_tolerance,
                                         rt_tolerance,
                                         isotope_mass_tolerance_ppm = NULL,
                                         pattern_cache) {
  query <- dplyr::tibble(...)

  # Look up precomputed pattern instead of computing

  isotopic_pattern <- pattern_cache[[query$molecular_formula]]

  if (is.null(isotopic_pattern)) {
    return(dplyr::tibble())  # No pattern available for this formula
  }

  # Filter isotopes from peaks (same logic as filter_isotopes)
  isotopes <- dplyr::mutate(
    peaks,
    mass_number_difference = round(mz - query$expected_mass),
    compound = query$compound,
    compound_id = if ("compound_id" %in% names(query)) query$compound_id else NA_character_,
    adduct = query$adduct,
    molecular_formula = query$molecular_formula,
    name = if ("name" %in% names(query)) query$name else NA_character_,
    monoisotopic_mass = if ("monoisotopic_mass" %in% names(query)) query$monoisotopic_mass else NA_real_
  )
  isotopes <- dplyr::filter(
    isotopes,
    peak != query$peak,
    rt_cluster == query$rt_cluster,
    dplyr::near(rt, query$rt, rt_tolerance),
    dplyr::near(mass_defect, query$mass_defect, mass_defect_tolerance)
  )

  isotopes <- dplyr::distinct(isotopes)

  if (nrow(isotopes) == 0) {
    return(dplyr::tibble())
  }

  # Match isotopes by intensity (same logic as match_isotopes_by_intensity)
  isotopes <- dplyr::mutate(isotopes,
    relative_intensity = mean_intensity / query$mean_intensity
  )
  isotopes <- dplyr::left_join(isotopes,
    dplyr::select(isotopic_pattern, mass_number_difference, abund, exact_mass_diff),
    by = "mass_number_difference"
  )

  # Calculate expected isotope mass
  isotopes <- dplyr::mutate(isotopes,
    expected_mass = query$expected_mass + exact_mass_diff,
    mass_error_ppm = abs(mz - expected_mass) / expected_mass * 1e6
  )

  # Apply ppm mass accuracy filter if specified
  if (!is.null(isotope_mass_tolerance_ppm)) {
    isotopes <- dplyr::filter(isotopes, mass_error_ppm <= isotope_mass_tolerance_ppm)
  }

  # Apply intensity ratio filter
  isotopes <- dplyr::filter(
    isotopes,
    dplyr::near(relative_intensity, abund, relative_intensity * intensity_deviation_tolerance)
  )

  isotopes <- dplyr::select(isotopes, -c(abund, relative_intensity, exact_mass_diff, mass_error_ppm))
  isotopes
}

#' Compute isotopes using precomputed pattern cache
#'
#' @description
#' Same as compute_isotopes but uses a precomputed pattern cache for speed.
#' This avoids redundant rcdk calls when processing many annotations.
#'
#' @param annotation Annotation table
#' @param adduct_weights Adduct weights table
#' @param intensity_deviation_tolerance Intensity deviation tolerance
#' @param mass_defect_tolerance Mass defect tolerance
#' @param peak_table Peak table
#' @param rt_tolerance RT tolerance
#' @param isotope_mass_tolerance_ppm Mass tolerance in ppm
#' @param pattern_cache Named list of precomputed isotope patterns
#' @return Annotation table with isotopes added
#' @keywords internal
#' @import dplyr
#' @importFrom purrr pmap_dfr
compute_isotopes_with_cache <- function(annotation,
                                        adduct_weights,
                                        intensity_deviation_tolerance = 0.1,
                                        mass_defect_tolerance = 0,
                                        peak_table,
                                        rt_tolerance = 1,
                                        isotope_mass_tolerance_ppm = NULL,
                                        pattern_cache) {
  annotation <- dplyr::mutate(annotation, mass_number_difference = 0)
  adducts <- dplyr::semi_join(annotation, adduct_weights, by = "adduct")

  if (nrow(adducts) == 0) {
    return(annotation)
  }

  isotopes <- purrr::pmap_dfr(
    adducts,
    ~ detect_isotopic_peaks_cached(...,
        peaks = peak_table,
        rt_tolerance = rt_tolerance,
        intensity_deviation_tolerance = intensity_deviation_tolerance,
        mass_defect_tolerance = mass_defect_tolerance,
        isotope_mass_tolerance_ppm = isotope_mass_tolerance_ppm,
        pattern_cache = pattern_cache
    )
  )

  annotation <- dplyr::bind_rows(annotation, isotopes)
  annotation
}

#' Compute permutation-based p-values for annotation scores
#'
#' @description
#' Performs permutation testing to assess the statistical significance of
#' annotation scores. The null distribution is generated by permuting the
#' m/z values across peaks, which breaks the true peak-compound relationships
#' while preserving the peak correlation structure.
#'
#' The key insight is that the peak correlation matrix reflects co-elution
#' patterns (peaks that appear/disappear together across samples). This structure
#' is independent of the actual mz values. By shuffling mz but keeping the
#' correlation structure, we preserve the correlation evidence while breaking
#' the mz-to-compound relationship.
#'
#' P-values are computed using a global null distribution: for each observed
#' score, the p-value is the proportion of ALL null scores (across all
#' permutations) that are >= the observed score.
#'
#' @param annotation Data frame with annotation results including scores
#' @param peak_table Peak table with mz, rt columns
#' @param compound_table Compound database
#' @param adduct_table Adduct table for mass matching
#' @param adduct_weights Adduct weights table
#' @param mass_tolerance Mass tolerance for matching (default 5e-6 = 5 ppm)
#' @param time_tolerance Retention time tolerance for grouping
#' @param intensity_deviation_tolerance Intensity deviation tolerance for isotopes
#' @param mass_defect_tolerance Mass defect tolerance for isotopes
#' @param isotope_mass_tolerance_ppm Isotope mass tolerance in ppm
#' @param correlation_threshold Correlation threshold for scoring
#' @param filter_by Adducts to filter by for scoring
#' @param peak_correlation_matrix Pre-computed peak correlation matrix (from original data)
#' @param n_permutations Number of permutations (default 1000)
#' @param seed Random seed for reproducibility
#' @param n_cores Number of cores for parallel processing
#' @param method Permutation method: "full" (default, all permutations in parallel),
#'   or "streaming" (memory-efficient, processes in chunks)
#' @return Data frame with perm_pvalue column added
#' @export
compute_permutation_pvalues <- function(annotation,
                                        peak_table,
                                        compound_table,
                                        adduct_table,
                                        adduct_weights,
                                        mass_tolerance = 5e-6,
                                        time_tolerance = 10,
                                        intensity_deviation_tolerance = 0.1,
                                        mass_defect_tolerance = 0.1,
                                        isotope_mass_tolerance_ppm = 5,
                                        correlation_threshold = 0.7,
                                        filter_by = c("M-H", "M+H"),
                                        peak_correlation_matrix = NULL,
                                        n_permutations = 1000,
                                        seed = 42,
                                        n_cores = 1,
                                        method = "full") {
  set.seed(seed)

  observed_scores <- annotation$score
  n_annotations <- nrow(annotation)

  message(sprintf("=== Computing Permutation P-values ==="))
  message(sprintf("Permutations: %d, Cores: %d, Method: %s",
                  n_permutations, n_cores, method))

  # Capture function references for parallel processing
  # (forked processes need explicit reference, not namespace lookup)
  simple_annotation_fn <- simple_annotation
  compute_mass_defect_fn <- compute_mass_defect
  reformat_annotation_table_fn <- reformat_annotation_table
  reformat_correlation_matrix_fn <- reformat_correlation_matrix
  get_chemscore_fn <- get_chemscore

  # Explicitly capture data frames for forked process access
  local_adduct_table <- adduct_table
  local_compound_table <- compound_table
  local_adduct_weights <- adduct_weights

  # Precompute isotope patterns for all formulas (major speedup)
  # This eliminates redundant rcdk calls across permutations
  isotope_pattern_cache <- precompute_isotope_patterns(compound_table)

  # Debug: Check that peak_table has required columns
  required_cols <- c("peak", "mean_intensity", "module", "rt_cluster")
  missing_cols <- setdiff(required_cols, names(peak_table))
  if (length(missing_cols) > 0) {
    warning(sprintf("peak_table missing columns for permutation: %s. Available: %s",
                    paste(missing_cols, collapse = ", "),
                    paste(names(peak_table), collapse = ", ")))
  }

  # Debug: Check adduct_table columns
  message(sprintf("  [Debug] adduct_table cols: %s", paste(names(adduct_table), collapse = ", ")))
  message(sprintf("  [Debug] adduct_table has %d rows", nrow(adduct_table)))

  # Function to run single permutation (wrapped in tryCatch for parallel safety)
  # Returns all null scores for global null distribution
  run_permutation <- function(perm_id) {
    tryCatch({
      if (perm_id == 1) {
        message(sprintf("  [Debug] Perm 1: Starting. peak_table cols: %s", paste(names(peak_table), collapse = ", ")))
        message(sprintf("  [Debug] Perm 1: peak_table has %d rows", nrow(peak_table)))
      }

      # Shuffle mz values while preserving peak structure
      permuted_peaks <- peak_table
      permuted_peaks$mz <- sample(permuted_peaks$mz)

      # Step 1: Mass matching with shuffled mz
      null_annotation <- simple_annotation_fn(
        peak_table = permuted_peaks,
        compound_table = local_compound_table,
        adduct_table = local_adduct_table,
        mass_tolerance = mass_tolerance
      )

      if (nrow(null_annotation) == 0) {
        if (perm_id == 1) message("  [Debug] Perm 1: No matches after simple_annotation")
        return(numeric(0))
      }

      # Step 1b: Add mass_defect to annotation
      # (required by compute_isotopes and reformat_annotation_table)
      null_annotation <- compute_mass_defect_fn(null_annotation, precision = 0.01)

      # Step 1c: Join module/rt_cluster/mean_intensity from peak_table
      # (required by compute_isotopes and reformat_annotation_table)
      # The peak_table already has these columns from the main workflow
      null_annotation <- dplyr::inner_join(
        null_annotation,
        dplyr::select(permuted_peaks, peak, mean_intensity, module, rt_cluster),
        by = "peak"
      )

      if (nrow(null_annotation) == 0) {
        if (perm_id == 1) message("  [Debug] Perm 1: No rows after join with peak_table columns")
        return(numeric(0))
      }

      # Step 2: Compute isotopes (using precomputed pattern cache for speed)
      null_annotation <- compute_isotopes_with_cache(
        annotation = null_annotation,
        adduct_weights = local_adduct_weights,
        intensity_deviation_tolerance = intensity_deviation_tolerance,
        mass_defect_tolerance = mass_defect_tolerance,
        peak_table = permuted_peaks,
        rt_tolerance = time_tolerance,
        isotope_mass_tolerance_ppm = isotope_mass_tolerance_ppm,
        pattern_cache = isotope_pattern_cache
      )

      # Step 3: Reformat annotation
      null_annotation <- reformat_annotation_table_fn(null_annotation)

      if (nrow(null_annotation) == 0) {
        if (perm_id == 1) message("  [Debug] Perm 1: No rows after reformat_annotation_table")
        return(numeric(0))
      }

      # Step 4: Create re-indexed correlation matrix for permuted mz values
      # The original peak_correlation_matrix preserves co-elution patterns
      # Reformat re-indexes it by the new (permuted) mz_rt labels
      if (!is.null(peak_correlation_matrix)) {
        global_cor <- reformat_correlation_matrix_fn(permuted_peaks, peak_correlation_matrix)
      } else {
        # Fallback: identity matrix (no correlation evidence)
        unique_mz_rt <- paste0(null_annotation$mz, "_", null_annotation$time)
        unique_mz_rt <- unique(unique_mz_rt)
        global_cor <- diag(length(unique_mz_rt))
        rownames(global_cor) <- colnames(global_cor) <- unique_mz_rt
      }

      # Step 5: Compute chemical scores
      unique_chemicals <- null_annotation %>%
        dplyr::distinct(compound_id, .keep_all = TRUE)

      if (nrow(unique_chemicals) == 0) {
        if (perm_id == 1) message("  [Debug] Perm 1: No unique chemicals")
        return(numeric(0))
      }

      scored <- purrr::pmap_dfr(
        unique_chemicals,
        ~ get_chemscore_fn(...,
                          annotation = null_annotation,
                          adduct_weights = local_adduct_weights,
                          corthresh = correlation_threshold,
                          global_cor = global_cor,
                          max_diff_rt = time_tolerance,
                          filter.by = filter_by,
                          adduct_table = local_adduct_table
        )
      )

      if (is.null(scored) || nrow(scored) == 0) {
        if (perm_id == 1) message("  [Debug] Perm 1: No scores from get_chemscore")
        return(numeric(0))
      }

      # Return scores (get_chemscore returns cur_chem_score)
      if ("cur_chem_score" %in% names(scored)) {
        if (perm_id == 1) message(sprintf("  [Debug] Perm 1: Returning %d scores", length(scored$cur_chem_score)))
        return(scored$cur_chem_score)
      } else if ("score" %in% names(scored)) {
        if (perm_id == 1) message(sprintf("  [Debug] Perm 1: Returning %d scores (from 'score' column)", length(scored$score)))
        return(scored$score)
      } else {
        if (perm_id == 1) message(sprintf("  [Debug] Perm 1: No score column found. Columns: %s", paste(names(scored), collapse = ", ")))
        return(numeric(0))
      }
    }, error = function(e) {
      if (perm_id == 1) {
        message(sprintf("  [Debug] Perm 1 ERROR: %s", conditionMessage(e)))
        message(sprintf("  [Debug] Error class: %s", paste(class(e), collapse = ", ")))
        message(sprintf("  [Debug] Error call: %s", deparse(conditionCall(e))))
      }
      return(NULL)
    })
  }


  # Validate method parameter

  method <- match.arg(method, c("full", "streaming"))

  # Compute p-values using selected method
  if (method == "full") {
    p_values <- compute_full_pvalues(
      observed_scores = observed_scores,
      run_permutation = run_permutation,
      n_permutations = n_permutations,
      n_cores = n_cores
    )
  } else {
    p_values <- compute_streaming_pvalues(
      observed_scores = observed_scores,
      run_permutation = run_permutation,
      n_permutations = n_permutations,
      n_cores = n_cores
    )
  }

  annotation$perm_pvalue <- p_values

  message(sprintf("P-values computed for %d annotations", n_annotations))
  message(sprintf("Annotations with p < 0.05: %d (%.1f%%)",
                  sum(p_values < 0.05, na.rm = TRUE),
                  100 * sum(p_values < 0.05, na.rm = TRUE) / n_annotations))

  return(annotation)
}

#' Compute p-values using full parallel method
#'
#' @description
#' Runs all permutations in parallel at once. Faster than streaming but uses
#' more memory since all null score matrices are held simultaneously.
#' On Unix systems, uses mclapply; on Windows uses PSOCK clusters.
#'
#' @param observed_scores Vector of observed annotation scores
#' @param run_permutation Function that runs a single permutation
#' @param n_permutations Number of permutations to run
#' @param n_cores Number of cores for parallel processing
#' @return Vector of p-values
#' @keywords internal
compute_full_pvalues <- function(observed_scores,
                                 run_permutation,
                                 n_permutations,
                                 n_cores) {
  n_annotations <- length(observed_scores)

  if (n_cores > 1) {
    if (.Platform$OS.type == "unix") {
      all_results <- parallel::mclapply(
        seq_len(n_permutations),
        run_permutation,
        mc.cores = n_cores,
        mc.preschedule = FALSE  # Better error handling for forked processes
      )
    } else {
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      all_results <- parallel::parLapply(cl, seq_len(n_permutations), run_permutation)
      parallel::stopCluster(cl)
      on.exit(NULL)
    }

    # Check for errors (mclapply returns NULL or try-error on failure)
    failed <- sapply(all_results, function(x) {
      is.null(x) || inherits(x, "try-error") ||
        (is.list(x) && length(x) > 0 && inherits(x[[1]], "error"))
    })
    n_failed <- sum(failed)
    if (n_failed > 0) {
      warning(sprintf("%d of %d permutations failed. Using %d successful permutations.",
                      n_failed, n_permutations, n_permutations - n_failed))
      all_results <- all_results[!failed]
    }

    if (length(all_results) == 0) {
      stop("All permutations failed. Check that simple_annotation works correctly.")
    }

    n_successful <- length(all_results)
  } else {
    all_results <- lapply(seq_len(n_permutations), function(i) {
      if (i %% 10 == 0 || i == 1) {
        message(sprintf("Completed %d/%d permutations", i, n_permutations))
      }
      run_permutation(i)
    })
    n_successful <- n_permutations
  }

  message(sprintf("Completed %d/%d permutations", n_successful, n_permutations))

  # Collect all null scores into global null distribution
  all_null_scores <- unlist(all_results)
  all_null_scores <- all_null_scores[!is.na(all_null_scores)]

  message(sprintf("Global null distribution: %d scores", length(all_null_scores)))

  if (length(all_null_scores) == 0) {
    warning("No null scores generated. Returning p-value = 1 for all annotations.")
    return(rep(1, n_annotations))
  }

  # For each observed score, calculate p-value against global null distribution
  # P-value = (number of null scores >= observed + 1) / (total null scores + 1)
  p_values <- sapply(observed_scores, function(obs) {
    (sum(all_null_scores >= obs) + 1) / (length(all_null_scores) + 1)
  })

  return(p_values)
}

#' Compute p-values using streaming method (memory efficient)
#'
#' @description
#' Processes permutations in chunks to minimize memory usage. On Unix systems,
#' uses mclapply for parallel processing; on Windows uses PSOCK clusters.
#'
#' @param observed_scores Vector of observed annotation scores
#' @param run_permutation Function that runs a single permutation
#' @param n_permutations Number of permutations to run
#' @param n_cores Number of cores for parallel processing
#' @return Vector of p-values
#' @keywords internal
compute_streaming_pvalues <- function(observed_scores,
                                      run_permutation,
                                      n_permutations,
                                      n_cores) {
  n_annotations <- length(observed_scores)
  # Track exceedance counts for each observed score against global null distribution
  exceedance_counts <- rep(0, n_annotations)
  total_null_scores <- 0
  n_successful <- 0
  n_failed_total <- 0

  if (n_cores > 1) {
    chunk_size <- min(50, n_permutations)

    for (chunk_start in seq(1, n_permutations, chunk_size)) {
      chunk_end <- min(chunk_start + chunk_size - 1, n_permutations)
      chunk_ids <- chunk_start:chunk_end

      if (.Platform$OS.type == "unix") {
        chunk_results <- parallel::mclapply(
          chunk_ids,
          run_permutation,
          mc.cores = n_cores,
          mc.preschedule = FALSE  # Better error handling for forked processes
        )
      } else {
        cl <- parallel::makeCluster(n_cores)
        on.exit(parallel::stopCluster(cl), add = TRUE)
        chunk_results <- parallel::parLapply(cl, chunk_ids, run_permutation)
        parallel::stopCluster(cl)
        on.exit(NULL)
      }

      # Check for errors in this chunk
      failed <- sapply(chunk_results, function(x) {
        is.null(x) || inherits(x, "try-error") ||
          (is.list(x) && length(x) > 0 && inherits(x[[1]], "error"))
      })
      n_failed_chunk <- sum(failed)
      n_failed_total <- n_failed_total + n_failed_chunk

      # Count exceedances against global null distribution
      for (null_scores in chunk_results[!failed]) {
        if (!is.null(null_scores) && length(null_scores) > 0) {
          # For each observed score, count how many null scores are >= it
          for (obs_i in seq_along(observed_scores)) {
            exceedance_counts[obs_i] <- exceedance_counts[obs_i] +
              sum(null_scores >= observed_scores[obs_i])
          }
          total_null_scores <- total_null_scores + length(null_scores)
          n_successful <- n_successful + 1
        }
      }

      gc()
      message(sprintf("Completed %d/%d permutations", chunk_end, n_permutations))
    }

    if (n_failed_total > 0) {
      warning(sprintf("%d of %d permutations failed. Using %d successful permutations.",
                      n_failed_total, n_permutations, n_successful))
    }

    if (n_successful == 0) {
      stop("All permutations failed. Check that simple_annotation works correctly.")
    }
  } else {
    for (i in seq_len(n_permutations)) {
      null_scores <- run_permutation(i)
      if (!is.null(null_scores) && length(null_scores) > 0) {
        # For each observed score, count how many null scores are >= it
        for (obs_i in seq_along(observed_scores)) {
          exceedance_counts[obs_i] <- exceedance_counts[obs_i] +
            sum(null_scores >= observed_scores[obs_i])
        }
        total_null_scores <- total_null_scores + length(null_scores)
        n_successful <- n_successful + 1
      } else {
        n_failed_total <- n_failed_total + 1
      }

      if (i %% 100 == 0) {
        message(sprintf("Completed %d/%d permutations", i, n_permutations))
      }
    }

    if (n_failed_total > 0) {
      warning(sprintf("%d of %d permutations failed. Using %d successful permutations.",
                      n_failed_total, n_permutations, n_successful))
    }

    if (n_successful == 0) {
      stop("All permutations failed. Check that simple_annotation works correctly.")
    }
  }

  message(sprintf("Global null distribution: %d scores from %d permutations",
                  total_null_scores, n_successful))

  if (total_null_scores == 0) {
    warning("No null scores generated. Returning p-value = 1 for all annotations.")
    return(rep(1, n_annotations))
  }

  # Calculate p-values with +1 correction (prevents p=0)
  # P-value = (exceedance count + 1) / (total null scores + 1)
  p_values <- (exceedance_counts + 1) / (total_null_scores + 1)

  return(p_values)
}
