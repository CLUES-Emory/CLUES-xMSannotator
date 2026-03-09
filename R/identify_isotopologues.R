#' Generate isotopologue labels for a molecular formula using enviPat.
#'
#' Computes the fine isotopologue structure and builds human-readable labels
#' identifying which heavy isotope substitutions are present (e.g., "13C:1", "15N:1").
#'
#' @param formula Molecular formula string (e.g., "C5H11N").
#' @param isotopes_data enviPat isotope data (from \code{data(isotopes, package = "enviPat")}).
#' @return Data frame with columns: mass_number_diff, neutral_mass, relative_abundance, isotopologue_label.
#' @keywords internal
get_isotopologue_labels <- function(formula, isotopes_data) {
  ep_result <- suppressMessages(enviPat::isopattern(
    isotopes_data, formula,
    threshold = 0.01,
    charge = FALSE,
    emass = 0.00054857990924,
    verbose = FALSE
  ))
  ep_matrix <- ep_result[[1]]

  element_cols <- colnames(ep_matrix)[-(1:2)]
  mono_names <- get_monoisotopic_names(element_cols)

  labels <- apply(ep_matrix[, element_cols, drop = FALSE], 1, function(row) {
    heavy <- which(row > 0 & !names(row) %in% mono_names)
    if (length(heavy) == 0) return("monoisotopic")
    paste(names(row)[heavy], row[heavy], sep = ":", collapse = " ")
  })

  neutral_masses <- as.numeric(ep_matrix[, 1])
  abundances <- as.numeric(ep_matrix[, 2])
  mass_diffs <- as.integer(round(neutral_masses - neutral_masses[1]))

  data.frame(
    mass_number_diff = mass_diffs,
    neutral_mass = neutral_masses,
    relative_abundance = abundances / max(abundances),
    isotopologue_label = labels,
    stringsAsFactors = FALSE
  )
}

#' Identify isotopologues in annotation output.
#'
#' For each isotope row (Adduct contains "_[+N]" suffix), uses enviPat to determine
#' which specific isotopologue the observed peak corresponds to (e.g., 13C:1 vs 15N:1),
#' based on closest m/z match. Adds \code{isotopologue} and \code{isotopologue_quality} columns.
#'
#' @param annotation Annotation data frame from \code{multilevelannotationstep4()}.
#'   Must contain columns: Adduct, Formula, mz, compound_id, Module_RTclust, mean_int_vec.
#' @param adduct_table Adduct parameter table with columns: adduct, charge, factor, mass.
#' @param isotope_mass_tolerance_ppm Maximum allowed mass error in ppm between observed
#'   and theoretical isotopologue m/z. If NULL, no ppm filter is applied.
#' @param intensity_deviation_tolerance Numeric threshold for abundance ratio validation.
#'   Uses the same formula as isotope detection: the absolute difference between observed
#'   and theoretical intensity ratios must be within \code{observed_ratio * tolerance}.
#'
#' @return Annotation data frame with two additional columns:
#'   \describe{
#'     \item{isotopologue}{Character. Isotopologue identity label (e.g., "13C:1", "15N:1",
#'       "13C:1 37Cl:1"). NA for monoisotopic peaks or unmatched isotopes.}
#'     \item{isotopologue_quality}{Character. Match quality: "confirmed" (m/z and abundance
#'       match), "mz_only" (m/z match but abundance outside tolerance), or NA (no match).}
#'   }
#'
#' @keywords internal
identify_isotopologues <- function(annotation,
                                   adduct_table,
                                   isotope_mass_tolerance_ppm = NULL,
                                   intensity_deviation_tolerance = 0.3) {
  # enviPat is a required dependency (Imports) so this should never trigger
  if (!requireNamespace("enviPat", quietly = TRUE)) {
    stop("enviPat package is required but not installed. ",
         "Install with: install.packages('enviPat')")
  }

  annotation$isotopologue <- NA_character_
  annotation$isotopologue_quality <- NA_character_

  # Identify isotope rows via Adduct suffix pattern (e.g., "M+H_[+1]")
  isotope_regex <- "_\\[\\+[0-9]+\\]$"
  is_isotope <- grepl(isotope_regex, annotation$Adduct)

  if (!any(is_isotope)) return(annotation)

  # Parse mass_number_difference from suffix
  mnd_parsed <- rep(0L, nrow(annotation))
  mnd_parsed[is_isotope] <- as.integer(
    gsub(".*_\\[\\+(\\d+)\\]$", "\\1", annotation$Adduct[is_isotope])
  )

  # Load enviPat isotope data once
  isotopes_data <- NULL
  data(isotopes, package = "enviPat", envir = environment())
  isotopes_data <- isotopes

  # Cache enviPat results per formula
  formula_cache <- list()

  n_matched <- 0L
  n_confirmed <- 0L

  for (i in which(is_isotope)) {
    formula <- annotation$Formula[i]
    base_adduct <- strip_isotope_suffix(annotation$Adduct[i])
    mnd <- mnd_parsed[i]
    observed_mz <- annotation$mz[i]

    # Get or compute enviPat pattern for this formula
    if (is.null(formula_cache[[formula]])) {
      tryCatch({
        formula_cache[[formula]] <- get_isotopologue_labels(formula, isotopes_data)
      }, error = function(e) {
        warning(sprintf("enviPat failed for formula '%s': %s", formula, e$message))
        formula_cache[[formula]] <<- data.frame(
          mass_number_diff = integer(0),
          neutral_mass = numeric(0),
          relative_abundance = numeric(0),
          isotopologue_label = character(0),
          stringsAsFactors = FALSE
        )
      })
    }
    ep <- formula_cache[[formula]]

    # Filter to matching mass_number_difference
    candidates <- ep[ep$mass_number_diff == mnd, ]
    if (nrow(candidates) == 0) next

    # Get adduct parameters for m/z calculation
    adduct_row <- adduct_table[adduct_table$adduct == base_adduct, ]
    if (nrow(adduct_row) == 0) next

    # Calculate theoretical m/z for each candidate isotopologue
    theo_mz <- (adduct_row$factor * candidates$neutral_mass +
                  adduct_row$mass) / abs(adduct_row$charge)

    # Match to closest by m/z
    best_idx <- which.min(abs(theo_mz - observed_mz))
    best_mz <- theo_mz[best_idx]

    # PPM check: reject if mass error exceeds threshold
    mass_error_ppm <- abs(observed_mz - best_mz) / best_mz * 1e6
    if (!is.null(isotope_mass_tolerance_ppm) && mass_error_ppm > isotope_mass_tolerance_ppm) {
      next
    }

    annotation$isotopologue[i] <- candidates$isotopologue_label[best_idx]
    n_matched <- n_matched + 1L

    # Abundance check: compare observed vs theoretical intensity ratio
    # Find monoisotopic parent: same compound_id, base adduct, Module_RTclust
    parent_idx <- which(
      annotation$compound_id == annotation$compound_id[i] &
        annotation$Adduct == base_adduct &
        annotation$Module_RTclust == annotation$Module_RTclust[i]
    )

    if (length(parent_idx) == 1) {
      observed_ratio <- annotation$mean_int_vec[i] / annotation$mean_int_vec[parent_idx]
      theoretical_ratio <- candidates$relative_abundance[best_idx]

      if (abs(observed_ratio - theoretical_ratio) <= observed_ratio * intensity_deviation_tolerance) {
        annotation$isotopologue_quality[i] <- "confirmed"
        n_confirmed <- n_confirmed + 1L
      } else {
        annotation$isotopologue_quality[i] <- "mz_only"
      }
    } else {
      annotation$isotopologue_quality[i] <- "mz_only"
    }
  }

  n_isotopes <- sum(is_isotope)
  cat(sprintf("\n=== Isotopologue Identification ===\n"))
  cat(sprintf("Isotope rows: %d | Identified: %d | Confirmed: %d | mz_only: %d\n",
              n_isotopes, n_matched, n_confirmed, n_matched - n_confirmed))
  cat(sprintf("Unique formulas processed: %d\n", length(formula_cache)))
  cat(sprintf("===================================\n\n"))

  return(annotation)
}
