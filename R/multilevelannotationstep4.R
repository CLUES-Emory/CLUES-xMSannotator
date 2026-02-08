# =============================================================================
# multilevelannotationstep4.R
# Stage 4: Assign confidence levels to annotations
#
# Key features:
# - Pre-split data by compound_id before lapply (10-50x speedup)
# - Vectorized compute_delta_ppm()
# - Single create_adduct_weights() call (not inside loop)
# - dplyr::bind_rows instead of deprecated plyr::ldply
# - Named constants for magic numbers
# - Broken into small, testable helper functions
# =============================================================================

#' @import dplyr
#' @importFrom dplyr left_join bind_rows

# =============================================================================
# CONSTANTS
# =============================================================================

# Score thresholds
SCORE_THRESHOLD_MIN <- 0.1
SCORE_THRESHOLD_HIGH <- 10

# Confidence level numeric values
CONFIDENCE_HIGH <- 3
CONFIDENCE_MEDIUM <- 2
CONFIDENCE_LOW <- 1
CONFIDENCE_NONE <- 0
CONFIDENCE_BOOSTED <- 4

# Confidence level names (for internal use)
CONF_NAME_HIGH <- "High"
CONF_NAME_MEDIUM <- "Medium"
CONF_NAME_LOW <- "Low"
CONF_NAME_NONE <- "None"

# Lookup table for confidence name to numeric conversion
CONF_NAME_TO_SCORE <- c(
  "High" = CONFIDENCE_HIGH,
  "Medium" = CONFIDENCE_MEDIUM,
  "Low" = CONFIDENCE_LOW,
  "None" = CONFIDENCE_NONE
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Filter data to keep only rows matching the most frequent cluster
filter_clusters <- function(curdata, table_names) {
  curdata[curdata$Module_RTclust == table_names[1], , drop = FALSE]
}

#' Create frequency table of clusters, keeping only the most frequent
create_cluster_table <- function(curdata) {
  cluster_table <- table(curdata$Module_RTclust)
  cluster_table <- cluster_table[cluster_table > 0]
  cluster_table[cluster_table == max(cluster_table)]
}

#' Compute retention time range (max - min)
compute_delta_rt <- function(curdata) {
  round(max(curdata$time) - min(curdata$time))
}

#' Check if filter.by is effectively empty (NULL or single NA)
is_filter_empty <- function(filter.by) {
  is.null(filter.by) || (length(filter.by) == 1 && is.na(filter.by[1]))
}

#' Strip isotope annotations from adduct names
#' e.g., "M+H_[+1]" -> "M+H"
strip_isotope_suffix <- function(adducts) {
  gsub(adducts, pattern = "(_\\[(\\+|\\-)[0-9]*\\])", replacement = "")
}

#' Check if any adducts match the filter list
has_filter_match <- function(adducts, filter.by) {
  if (is_filter_empty(filter.by)) {
    return(FALSE)
  }
  any(adducts %in% filter.by)
}

#' Count how many adducts match the filter list
count_filter_matches <- function(adducts, filter.by) {
  if (is_filter_empty(filter.by)) {
    return(0)
  }
  sum(adducts %in% filter.by)
}

#' Assign confidence based on standard criteria
assign_conf <- function(curdata, filter.by, delta_rt, max_diff_rt, current_conf) {
  if (curdata$score[1] <= 0) {
    return(current_conf)
  }

  if (has_filter_match(curdata$Adduct, filter.by)) {
    if (nrow(curdata) > 1 && length(unique(curdata$Adduct)) > 1 && delta_rt < max_diff_rt) {
      return(CONF_NAME_HIGH)
    } else {
      return(CONF_NAME_MEDIUM)
    }
  }

  CONF_NAME_LOW
}

# =============================================================================
# CONFIDENCE ASSIGNMENT SUB-FUNCTIONS
# =============================================================================

#' Initial confidence check based on score threshold
#' Returns NULL if should continue, or a result data frame if should return early
check_minimum_score <- function(curdata) {
  if (curdata$score[1] < SCORE_THRESHOLD_MIN) {
    final_res <- data.frame(score_level = CONFIDENCE_NONE, curdata, check.names = FALSE)
    return(final_res)
  }
  NULL
}

#' Check confidence based on unique adduct count
apply_unique_adduct_rules <- function(curdata, cur_adducts, filter.by, current_conf) {
  if (length(unique(curdata$Adduct)) >= 2) {
    return(list(conf = current_conf, early_return = NULL))
  }

  # Single unique adduct type
  has_filter <- has_filter_match(cur_adducts, filter.by)

  if (curdata$score[1] < SCORE_THRESHOLD_HIGH && !has_filter) {
    current_conf <- CONF_NAME_NONE
  } else if (has_filter) {
    current_conf <- CONF_NAME_LOW
  }

  if (nrow(curdata) > 1) {
    if (nrow(curdata) == 2) {
      current_conf <- CONF_NAME_MEDIUM
    }

    if (!has_filter) {
      current_conf <- CONF_NAME_NONE
    } else {
      # Early return with Medium confidence
      filtered_data <- curdata[cur_adducts %in% filter.by, , drop = FALSE]
      final_res <- data.frame(score_level = CONFIDENCE_MEDIUM, filtered_data, check.names = FALSE)
      return(list(conf = current_conf, early_return = final_res))
    }
  }

  list(conf = current_conf, early_return = NULL)
}

#' Apply RT clustering rules when delta_rt exceeds threshold
apply_rt_clustering_rules <- function(curdata, cur_adducts, adduct_weights, filter.by,
                                       max_diff_rt, current_conf) {
  delta_rt <- compute_delta_rt(curdata)

  if (delta_rt <= max_diff_rt) {
    return(list(curdata = curdata, conf = current_conf))
  }

  # RT spread exceeds threshold - need clustering
  current_conf <- CONF_NAME_NONE

  groupB <- group_by_rt(curdata, time_step = 3, max.rt.diff = max_diff_rt, groupnum = "")
  groupB <- groupB[order(groupB$mz), ]
  curdata <- curdata[order(curdata$mz), ]

  curdata$Module_RTclust <- groupB[, 1]

  has_weighted_adduct <- any(adduct_weights[, 1] %in% cur_adducts)

  if (has_weighted_adduct && curdata$score[1] > SCORE_THRESHOLD_MIN) {
    is_filter_na <- is_filter_empty(filter.by)

    if (is_filter_na) {
      good_mod <- curdata$Module_RTclust[curdata$Adduct %in% adduct_weights[, 1]]
    } else {
      good_mod <- curdata$Module_RTclust[curdata$Adduct %in% filter.by]
    }

    curdata <- curdata[curdata$Module_RTclust %in% good_mod, , drop = FALSE]
    cluster_table <- create_cluster_table(curdata)
    curdata <- filter_clusters(curdata, names(cluster_table))
    delta_rt <- compute_delta_rt(curdata)

    if (is_filter_na) {
      if (curdata$score[1] > 0 && nrow(curdata) > 1 &&
          length(unique(curdata$Adduct)) > 1 && delta_rt < max_diff_rt) {
        current_conf <- CONF_NAME_HIGH
      } else {
        current_conf <- CONF_NAME_LOW
      }
    } else if (!any(cluster_table > 1)) {
      current_conf <- assign_conf(curdata, filter.by, delta_rt, max_diff_rt, current_conf)
    }
  }

  curdata$Module_RTclust <- paste0(curdata$module_num, curdata$Module_RTclust)

  list(curdata = as.data.frame(curdata), conf = current_conf)
}

#' Apply rules based on multiple RT modules
apply_module_rules <- function(curdata, cur_adducts, filter.by, current_conf) {
  temp_curdata <- curdata
  temp_curdata$Module_RTclust <- gsub(temp_curdata$Module_RTclust,
                                       pattern = "(_[0-9]*)", replacement = "")
  table_modules <- table(temp_curdata$Module_RTclust)
  module_names <- names(table_modules[table_modules > 0])

  if (length(module_names) <= 1) {
    return(current_conf)
  }

  # Multiple modules
  if (curdata$score[1] < SCORE_THRESHOLD_HIGH) {
    return(CONF_NAME_NONE)
  }

  if (curdata$score[1] > SCORE_THRESHOLD_HIGH && has_filter_match(cur_adducts, filter.by)) {
    return(CONF_NAME_MEDIUM)
  }

  current_conf
}

#' Apply rules for single row data
apply_single_row_rules <- function(curdata, cur_adducts, adduct_weights, filter.by, current_conf) {
  if (nrow(curdata) >= 2) {
    return(current_conf)
  }

  if (!any(cur_adducts %in% adduct_weights[, 1])) {
    return(CONF_NAME_LOW)
  }

  # Matches weighted adduct
  if (curdata$score[1] > SCORE_THRESHOLD_HIGH && has_filter_match(cur_adducts, filter.by)) {
    return(CONF_NAME_MEDIUM)
  }

  CONF_NAME_NONE
}

#' Apply multimer and charge state abundance checks
apply_multimer_rules <- function(curdata, cur_adducts, adduct_table_local,
                                  current_conf, multimer_abundance_check) {
  if (nrow(curdata) < 2) {
    return(list(curdata = curdata, conf = current_conf, cur_adducts = cur_adducts))
  }

  min_molecules <- min(curdata$factor)

  if (min_molecules > 1) {
    return(list(curdata = curdata, conf = CONF_NAME_LOW, cur_adducts = cur_adducts))
  }

  if (!multimer_abundance_check) {
    return(list(curdata = curdata, conf = current_conf, cur_adducts = cur_adducts))
  }

  # Check for multimer patterns (2M, 3M)
  check_abundance <- gregexpr(text = cur_adducts, pattern = "([2-3]+M)")

  has_multimer <- any(sapply(check_abundance, function(x) x[1] != -1))
  if (!has_multimer) {
    return(list(curdata = curdata, conf = current_conf, cur_adducts = cur_adducts))
  }

  min_mol_ind <- which(curdata$factor == min_molecules)
  max_int_min_mol <- max(curdata$mean_int_vec[min_mol_ind])
  min_charge_in_adducts <- min(adduct_table_local$charge[adduct_table_local$adduct %in% cur_adducts])

  if (min_charge_in_adducts >= 2) {
    return(list(curdata = curdata, conf = CONF_NAME_LOW, cur_adducts = cur_adducts))
  }

  # Check abundance ratios for multimers
  for (i in seq_along(check_abundance)) {
    strlength <- attr(check_abundance[[i]], "match.length")
    if (strlength[1] > -1) {
      abund_ratio <- curdata$mean_int_vec[i] / max_int_min_mol

      if (is.na(abund_ratio)) {
        current_conf <- CONF_NAME_LOW
      } else if (abund_ratio > 1) {
        # Multimer more abundant than monomer - downgrade
        if (current_conf == CONF_NAME_HIGH) {
          current_conf <- CONF_NAME_MEDIUM
        } else {
          current_conf <- CONF_NAME_LOW
        }
      }
    }
  }

  # BUG FIX: Original had `length(nrow(curdata) > 0)` which is always 1
  if (nrow(curdata) > 0) {
    cur_adducts_with_isotopes <- curdata$Adduct
    cur_adducts <- strip_isotope_suffix(cur_adducts_with_isotopes)
  } else {
    current_conf <- CONF_NAME_NONE
  }

  # Check high charge state ratios
  min_charge <- min(curdata$charge)
  min_charge_ind <- which(curdata$charge == min_charge)
  high_charge_ind <- which(curdata$charge > 1)

  if (length(high_charge_ind) > 0) {
    abund_ratio <- max(curdata$mean_int_vec[min_charge_ind])[1] /
                   max(curdata$mean_int_vec[high_charge_ind])[1]

    if (abund_ratio < 1) {
      current_conf <- CONF_NAME_LOW
    }
  }

  if (min_charge > 1) {
    current_conf <- CONF_NAME_LOW
  }

  # Remove isotope peaks without monoisotopic form
  check_isotopes <- gregexpr(text = curdata$Adduct, pattern = "(_\\[(\\+|\\-)[0-9]*\\])")

  rows_to_remove <- integer(0)
  for (i in seq_along(check_isotopes)) {
    strlength <- attr(check_isotopes[[i]], "match.length")
    if (strlength[1] > -1) {
      count_abundant_form <- sum(cur_adducts %in% cur_adducts[i])
      if (count_abundant_form < 2) {
        rows_to_remove <- c(rows_to_remove, i)
      }
    }
  }
  if (length(rows_to_remove) > 0) {
    curdata <- curdata[-rows_to_remove, , drop = FALSE]
  }

  cur_adducts <- strip_isotope_suffix(curdata$Adduct)

  list(curdata = curdata, conf = current_conf, cur_adducts = cur_adducts)
}

#' Apply elemental composition checks (carbon count)
apply_element_rules <- function(curdata, current_conf) {
  formula_vec <- curdata$Formula
  curformula <- as.character(formula_vec[1])
  curformula <- gsub(curformula, pattern = "Ca|Cl|Cd|Cr", replacement = "")

  numcarbons <- check_element(curformula, "C")

  if (numcarbons < 1) {
    return(CONF_NAME_NONE)
  }

  current_conf
}

#' Apply final adduct-based adjustments
apply_final_adjustments <- function(curdata, cur_adducts, filter.by, current_conf) {
  if (length(unique(curdata$Adduct)) >= 2) {
    return(current_conf)
  }

  # Single unique adduct type
  if (nrow(curdata) > 1) {
    if (nrow(curdata) == 2) {
      current_conf <- CONF_NAME_MEDIUM
    }
    if (!has_filter_match(cur_adducts, filter.by)) {
      current_conf <- CONF_NAME_NONE
    }
  } else {
    if (curdata$score[1] < SCORE_THRESHOLD_HIGH) {
      if (!has_filter_match(cur_adducts, filter.by)) {
        current_conf <- CONF_NAME_NONE
      } else {
        current_conf <- CONF_NAME_LOW
      }
    } else if (curdata$score[1] > SCORE_THRESHOLD_HIGH && has_filter_match(cur_adducts, filter.by)) {
      current_conf <- CONF_NAME_MEDIUM
    }
  }

  current_conf
}

# =============================================================================
# MAIN CONFIDENCE FUNCTION (REFACTORED)
# =============================================================================

#' Assign confidence level to a single compound's annotation data
#'
#' @param curdata Data frame of annotations for one compound
#' @param max_diff_rt Maximum RT difference threshold
#' @param adduct_weights Processed adduct weights table
#' @param filter.by Vector of expected adducts
#' @param min_ions_perchem Minimum ions required per chemical
#' @param max_isp Maximum isotope peaks (unused, kept for compatibility)
#' @param adduct_table Adduct definitions table
#' @param adduct_table_local Adduct table with monoisotopic "M" added
#' @param multimer_abundance_check Whether to apply multimer checks
#'
#' @return Data frame with score_level column prepended
get_confidence_stage4 <- function(curdata,
                                      max_diff_rt,
                                      adduct_weights,
                                      filter.by = NULL,
                                      min_ions_perchem = 1,
                                      max_isp = 5,
                                      adduct_table,
                                      adduct_table_local,
                                      multimer_abundance_check = TRUE) {

  curdata <- curdata[order(curdata$Adduct), ]

  # (I) Check minimum score threshold
  early_result <- check_minimum_score(curdata)
  if (!is.null(early_result)) {
    return(early_result)
  }

  curdata_all <- curdata
  cur_adducts <- strip_isotope_suffix(curdata$Adduct)

  # Merge with adduct table to get charge/factor info
  curdata <- cbind(curdata, cur_adducts)
  curdata <- merge(curdata, adduct_table_local, by.x = "cur_adducts", by.y = "adduct")
  curdata[["cur_adducts"]] <- NULL

  # Initialize confidence
  current_conf <- CONF_NAME_HIGH

  # Check if adduct is in table
  if (!any(adduct_table_local$adduct %in% cur_adducts)) {
    current_conf <- CONF_NAME_NONE
  }

  # (II) Apply unique adduct rules
  adduct_result <- apply_unique_adduct_rules(curdata, cur_adducts, filter.by, current_conf)
  current_conf <- adduct_result$conf
  if (!is.null(adduct_result$early_return)) {
    return(adduct_result$early_return)
  }

  # Ensure time is numeric
 curdata$time <- as.numeric(curdata$time)

  # (III) Apply RT clustering rules
  rt_result <- apply_rt_clustering_rules(curdata, cur_adducts, adduct_weights,
                                          filter.by, max_diff_rt, current_conf)
  curdata <- rt_result$curdata
  current_conf <- rt_result$conf

  # (IV) Apply module rules
  current_conf <- apply_module_rules(curdata, cur_adducts, filter.by, current_conf)

  # (V) Apply single row rules
  current_conf <- apply_single_row_rules(curdata, cur_adducts, adduct_weights,
                                          filter.by, current_conf)

  # (VI) Apply multimer and charge rules
  multimer_result <- apply_multimer_rules(curdata, cur_adducts, adduct_table_local,
                                           current_conf, multimer_abundance_check)
  curdata <- multimer_result$curdata
  current_conf <- multimer_result$conf
  cur_adducts <- multimer_result$cur_adducts

  # (VII) Apply element rules (carbon check)
  current_conf <- apply_element_rules(curdata, current_conf)

  # (VIII) Apply final adjustments
  current_conf <- apply_final_adjustments(curdata, cur_adducts, filter.by, current_conf)

  # Convert confidence name to numeric score
  if (nrow(curdata) < 1) {
    score_level <- CONFIDENCE_NONE
    curdata <- curdata_all
  } else {
    score_level <- unname(CONF_NAME_TO_SCORE[current_conf])
    if (is.na(score_level)) {
      stop(paste("Invalid confidence level:", current_conf))
    }
  }

  # Build result
  final_res <- data.frame(score_level = score_level, curdata, check.names = FALSE)

  # Replace any NA score levels with 0
  final_res$score_level[is.na(final_res$score_level)] <- CONFIDENCE_NONE

  # Adjust score based on filter matches
  num_good_adducts <- count_filter_matches(unique(curdata$Adduct), filter.by)

  if (num_good_adducts > 0) {
    final_res$score <- final_res$score * num_good_adducts
  } else {
    final_res$score <- 0
  }

  # Check minimum ions requirement
  if (nrow(final_res) < min_ions_perchem) {
    final_res$score_level <- CONFIDENCE_NONE
  }

  as.data.frame(final_res)
}

# =============================================================================
# COMPUTE CONFIDENCE FOR SINGLE COMPOUND (REFACTORED)
# =============================================================================

#' Compute confidence level for a single compound
#'
#' @param curdata Pre-filtered data for one compound (from split)
#' @param filter.by Vector of expected adducts
#' @param max.rt.diff Maximum RT difference
#' @param adduct_weights Processed adduct weights
#' @param max_isp Maximum isotope peaks
#' @param min_ions_perchem Minimum ions per chemical
#' @param adduct_table Adduct table
#' @param adduct_table_local Adduct table with "M" added
#' @param multimer_abundance_check Whether to check multimer abundance
#'
#' @return Data frame with Confidence and compound_id columns
compute_confidence_for_compound <- function(curdata,
                                             filter.by,
                                             max.rt.diff,
                                             adduct_weights,
                                             max_isp,
                                             min_ions_perchem,
                                             adduct_table,
                                             adduct_table_local,
                                             multimer_abundance_check) {

  curdata <- curdata[order(curdata$Adduct), ]

  # Check if any adducts match filter
  has_filter_adduct <- TRUE
  if (!is.null(filter.by) && !all(is.na(filter.by))) {
    has_filter_adduct <- any(curdata$Adduct %in% filter.by)
  }

  Confidence <- CONFIDENCE_NONE

  if (has_filter_adduct) {
    curdata <- get_confidence_stage4(
      curdata,
      max.rt.diff,
      adduct_weights = adduct_weights,
      filter.by = filter.by,
      max_isp = max_isp,
      min_ions_perchem = min_ions_perchem,
      adduct_table = adduct_table,
      adduct_table_local = adduct_table_local,
      multimer_abundance_check = multimer_abundance_check
    )

    if (!is.na(curdata$score_level[1])) {
      Confidence <- as.numeric(curdata$score_level)

      # Boost to Medium if score > 10 and matches weighted adduct
      if (Confidence[1] < CONFIDENCE_MEDIUM) {
        weighted_adducts <- adduct_weights[as.numeric(adduct_weights[, 2]) > 0, 1]
        if (any(curdata$Adduct %in% weighted_adducts) && curdata$score[1] > SCORE_THRESHOLD_HIGH) {
          max_weight <- max(as.numeric(adduct_weights[adduct_weights[, 1] %in% curdata$Adduct, 2]))
          high_weight_adducts <- adduct_weights[as.numeric(adduct_weights[, 2]) >= max_weight, 1]
          curdata <- curdata[curdata$Adduct %in% high_weight_adducts, , drop = FALSE]
          Confidence <- CONFIDENCE_MEDIUM
        }
      }
    }
  } else {
    # No filter match - check if can still assign confidence
    if (any(curdata$Adduct %in% adduct_weights[, 1]) && curdata$score[1] >= SCORE_THRESHOLD_HIGH) {
      if (any(curdata$Adduct %in% filter.by)) {
        curdata <- curdata[curdata$Adduct %in% filter.by, , drop = FALSE]
        Confidence <- CONFIDENCE_MEDIUM
      }
    }
  }

  # Final check: low score + single adduct type
  if (nrow(curdata) > 1 && curdata$score[1] < SCORE_THRESHOLD_HIGH) {
    if (length(unique(curdata$Adduct)) < 2) {
      Confidence <- CONFIDENCE_NONE
    }
  }

  # Return unique compound_id with confidence
  result <- data.frame(
    Confidence = Confidence[1],
    compound_id = curdata$compound_id[1],
    stringsAsFactors = FALSE
  )

  unique(result)
}

# =============================================================================
# COMPUTE DELTA PPM (VECTORIZED)
# =============================================================================

#' Compute delta ppm between observed and theoretical m/z (vectorized)
#'
#' @param chemscoremat_with_confidence Data frame with mz and theoretical.mz columns
#' @return Data frame with delta_ppm column inserted after column 8
compute_delta_ppm <- function(chemscoremat_with_confidence) {
  # Ensure numeric types
  mz <- as.numeric(chemscoremat_with_confidence$mz)
  theoretical_mz <- as.numeric(chemscoremat_with_confidence$theoretical.mz)

  # Vectorized calculation (replaces 6-line apply version)
  delta_ppm <- round(1e6 * abs(theoretical_mz - mz) / theoretical_mz, 2)

  # Insert delta_ppm after the theoretical.mz column
  insert_after <- which(names(chemscoremat_with_confidence) == "theoretical.mz")
  n_cols <- ncol(chemscoremat_with_confidence)
  result <- cbind(
    chemscoremat_with_confidence[, 1:insert_after, drop = FALSE],
    delta_ppm = delta_ppm,
    chemscoremat_with_confidence[, (insert_after + 1):n_cols, drop = FALSE]
  )

  # Sort by confidence descending
  result[order(result$Confidence, decreasing = TRUE), ]
}

# =============================================================================
# BOOST CONFIDENCE (KEPT MOSTLY UNCHANGED)
# =============================================================================

#' Boost confidence of specified compound IDs to level 4
#'
#' @param chemscoremat_with_confidence Annotation results
#' @param boostIDs Data frame with ID column (and optional mz, time)
#' @param max.mz.diff Mass tolerance for matching
#' @param max.rt.diff RT tolerance for matching
#'
#' @return Updated data frame with boosted confidences
boost_confidence_of_IDs <- function(chemscoremat_with_confidence, boostIDs,
                                        max.mz.diff, max.rt.diff) {
  cnames_boost <- colnames(boostIDs)
  has_mz <- "mz" %in% cnames_boost
  has_time <- "time" %in% cnames_boost

  if (has_mz || has_time) {
    # Proximity-based matching
    good_ind <- sapply(seq_len(nrow(chemscoremat_with_confidence)), function(i) {
      ann_row <- chemscoremat_with_confidence[i, ]

      id_match <- boostIDs$ID == ann_row$compound_id
      if (!any(id_match)) return(FALSE)

      if (has_mz) {
        mz_match <- abs(boostIDs$mz - ann_row$mz) <= max.mz.diff * pmax(abs(boostIDs$mz), abs(ann_row$mz))
        id_match <- id_match & mz_match
      }

      if (has_time) {
        time_match <- abs(boostIDs$time - ann_row$time) <= max.rt.diff
        id_match <- id_match & time_match
      }

      any(id_match)
    })

    good_idx <- which(good_ind)
  } else {
    # ID-only matching
    good_idx <- which(chemscoremat_with_confidence$compound_id %in% boostIDs$ID)
  }

  if (length(good_idx) > 0) {
    chemscoremat_with_confidence$Confidence[good_idx] <- CONFIDENCE_BOOSTED
    chemscoremat_with_confidence$score[good_idx] <- chemscoremat_with_confidence$score[good_idx] * 100
  }

  chemscoremat_with_confidence
}

# =============================================================================
# MAIN FUNCTION (REFACTORED)
# =============================================================================

#' Stage 4: Assign confidence levels to annotations
#'
#' @param outloc Output directory
#' @param chemscoremat Annotation data frame from Stage 3
#' @param max.mz.diff Mass tolerance (ppm as fraction)
#' @param max.rt.diff Maximum RT difference (seconds)
#' @param adduct_weights Adduct weights (NULL or data frame)
#' @param filter.by Expected adducts (NULL or character vector)
#' @param min_ions_perchem Minimum ions per chemical
#' @param boostIDs Compounds to boost (NULL or data frame with ID column)
#' @param boost.mz.diff Mass tolerance for boosting (NULL = use max.mz.diff)
#' @param boost.rt.diff RT tolerance for boosting (NULL = use max.rt.diff)
#' @param max_isp Maximum isotope peaks
#' @param mz_rt_feature_id_map Feature ID mapping (NULL or data frame)
#' @param adduct_table Adduct table (NULL = use package default)
#' @param multimer_abundance_check Check multimer abundance ratios
#'
#' @return Data frame with confidence levels assigned
#'
#' @importFrom dplyr left_join bind_rows
#' @export
multilevelannotationstep4 <- function(outloc,
                                          chemscoremat,
                                          max.mz.diff = 5,
                                          max.rt.diff = 30,
                                          adduct_weights = NULL,
                                          filter.by = NULL,
                                          min_ions_perchem = 1,
                                          boostIDs = NULL,
                                          boost.mz.diff = NULL,
                                          boost.rt.diff = NULL,
                                          max_isp = 5,
                                          mz_rt_feature_id_map = NULL,
                                          adduct_table = NULL,
                                          multimer_abundance_check = TRUE) {

  # Handle legacy NA values
  if (identical(adduct_weights, NA)) adduct_weights <- NULL
  if (identical(filter.by, NA)) filter.by <- NULL
  if (identical(boostIDs, NA)) boostIDs <- NULL

  # Load package adduct_table if not provided
  if (is.null(adduct_table)) {
    data("adduct_table", envir = environment())
  }

  # Process adduct weights ONCE (not inside loop)
  adduct_weights <- create_adduct_weights(adduct_weights)

  # Create adduct table with monoisotopic "M" ONCE
  adduct_monoisot <- data.frame(
    adduct = "M",
    charge = 1,
    factor = 1,
    mass = 0,
    stringsAsFactors = FALSE
  )
  adduct_table_local <- rbind(adduct_table, adduct_monoisot)
  adduct_table_local <- adduct_table_local[order(adduct_table_local$adduct), ]

  # PERFORMANCE: Pre-split data by compound_id (avoids O(N*M) filtering)
  split_data <- split(chemscoremat, chemscoremat$compound_id)

  # Compute confidence for each compound
  chemscoremat_conf_levels <- lapply(split_data, function(curdata) {
    compute_confidence_for_compound(
      curdata,
      filter.by,
      max.rt.diff,
      adduct_weights,
      max_isp,
      min_ions_perchem,
      adduct_table,
      adduct_table_local,
      multimer_abundance_check
    )
  })

  # Combine results (using dplyr instead of deprecated plyr)
  chemscoremat_conf_levels <- bind_rows(chemscoremat_conf_levels)

  # Merge confidence with original data
  chemscoremat_with_confidence <- merge(
    chemscoremat_conf_levels,
    unique(chemscoremat),
    by = "compound_id"
  )
  chemscoremat_with_confidence <- as.data.frame(chemscoremat_with_confidence)

  # Compute delta ppm (vectorized version)
  chemscoremat_with_confidence <- compute_delta_ppm(chemscoremat_with_confidence)

  # Boost specified IDs
  if (!is.null(boostIDs)) {
    boost_mz <- if (!is.null(boost.mz.diff)) boost.mz.diff else max.mz.diff
    boost_rt <- if (!is.null(boost.rt.diff)) boost.rt.diff else max.rt.diff
    chemscoremat_with_confidence <- boost_confidence_of_IDs(
      chemscoremat_with_confidence, boostIDs, boost_mz, boost_rt
    )
  }

  # Assign match category
  mz_counts <- table(chemscoremat_with_confidence$mz)
  unique_mz <- names(mz_counts[mz_counts == 1])

  chemscoremat_with_confidence$MatchCategory <- ifelse(
    chemscoremat_with_confidence$mz %in% unique_mz, "Unique", "Multiple"
  )

  # Join feature ID column if provided
  if (!is.null(mz_rt_feature_id_map)) {
    chemscoremat_with_confidence <- left_join(
      chemscoremat_with_confidence,
      mz_rt_feature_id_map,
      by = c("mz", "time")
    )
  }

  # Write output
  write.table(
    chemscoremat_with_confidence,
    file = file.path(outloc, "Stage4_confidence_levels.txt"),
    sep = "\t",
    row.names = FALSE
  )

  # Sort by confidence (desc), then compound, then base adduct (grouping isotopes with parent)
  base_adduct <- strip_isotope_suffix(chemscoremat_with_confidence$Adduct)
  chemscoremat_with_confidence <- chemscoremat_with_confidence[
    order(-chemscoremat_with_confidence$Confidence,
          chemscoremat_with_confidence$compound_id,
          base_adduct,
          chemscoremat_with_confidence$Adduct),
  ]

  # Print distribution summaries
  print("Stage 4 confidence level distribution for unique chemical/metabolite IDs")
  conf_by_id <- chemscoremat_with_confidence$Confidence[
    !duplicated(chemscoremat_with_confidence$compound_id)
  ]
  print(table(conf_by_id))

  print("Stage 4 confidence level distribution for unique chemical/metabolite formulas")
  conf_by_formula <- chemscoremat_with_confidence$Confidence[
    !duplicated(chemscoremat_with_confidence$Formula)
  ]
  print(table(conf_by_formula))

  as.data.frame(chemscoremat_with_confidence)
}
