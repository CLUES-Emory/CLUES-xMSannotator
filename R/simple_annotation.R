check_element <- function(formula, element) {
  stopifnot(length(element) == 1)

  pattern <- paste0(element, '(?![a-z]+)([0-9]*)')
  matches <- stringr::str_match_all(formula, pattern)

  sapply(matches, function(match) {
    counts <- as.integer(match[, 2])
    sum(counts, na.rm = TRUE) + sum(is.na(counts))
  })
}

check_golden_rules <- function(formula) {
  c <- check_element(formula, "C")
  h <- check_element(formula, "H")
  f <- check_element(formula, "F")
  n <- check_element(formula, "N")
  o <- check_element(formula, "O")
  p <- check_element(formula, "P")
  s <- check_element(formula, "S")

  ratio_check <- (0.1 <= h / c) &
    (h / c <= 6) &
    (f / c <= 6) &
    (n / c <= 4) &
    (o / c <= 3) &
    (p / c <= 2) &
    (s / c <= 3)
  ratio_check <- ifelse(c > 0, ratio_check, FALSE)

  nops <- !(n > 1 & o > 1 & p > 1 & s > 1) | (n < 10 & o < 20 & p < 4 & s < 3)
  nop <- !(n > 3 & o > 3 & p > 3) | (n < 11 & o < 22 & p < 6)
  ops <- !(o > 1 & p > 1 & s > 1) | (o < 14 & p < 3 & s < 3)
  psn <- !(n > 1 & p > 1 & s > 1) | (n < 4 & p < 3 & s < 3)
  nos <- !(n > 6 & o > 6 & s > 6) | (n < 19 & o < 14 & s < 8)

  return(ratio_check & nops & nop & ops & psn & nos)
}

forms_valid_adduct_pair <- function(molecular_formula, adduct_formula) {
  has_oxygen <- check_element(molecular_formula, "O") > 0
  is_water_adducts <- adduct_formula %in% c("M+H-H2O", "M+H-2H2O", "M-H2O-H")
  !is_water_adducts | has_oxygen
}

#' Is this compound allowed to be observed as this adduct, given its charge?
#'
#' A permanently charged compound is already an ion: it is observed as the bare
#' ion at mass / |z| and never picks up H+, Na+ or forms multimers, so it is
#' restricted to its own bare-ion pseudo-adduct. A neutral compound is the
#' converse: it can form any ionization adduct but is never observed at its
#' neutral mass, so every zero-mass-shift adduct is excluded, including the
#' stock "M". Excluding "M" for charged compounds also avoids reporting one hit
#' twice, since "M" and "[M]+" give an identical expected mass at z = 1.
#'
#' The exclusion of "M" for neutral compounds is only valid when the compound
#' table actually records charge. Databases without a `charge` column (HMDB,
#' for example) list permanently charged compounds such as choline at their ion
#' mass under a sign-less formula, and "M" is the only adduct that finds them.
#' `allow_neutral_mass = TRUE` keeps that match for such tables.
#'
#' @param compound_charge Integer vector of compound formal charges.
#' @param adduct Character vector of adduct names, same length.
#' @param allow_neutral_mass Logical. If TRUE, neutral (charge 0) compounds may
#'   also match the zero-shift "M" adduct. Set when the compound table carried
#'   no `charge` column, so charge 0 means "unknown" rather than "neutral".
#' @return Logical vector: TRUE for pairs that may be matched.
#' @keywords internal
forms_valid_charge_pair <- function(compound_charge, adduct, allow_neutral_mass = FALSE) {
  is_bare_ion <- startsWith(adduct, "[M]") | adduct == "M"
  neutral_ok <- if (isTRUE(allow_neutral_mass)) !startsWith(adduct, "[M]") else !is_bare_ion
  ifelse(
    compound_charge == 0,
    neutral_ok,
    adduct == bare_ion_adduct(compound_charge)
  )
}

#' @import dplyr
#' @importFrom rlang .data
compute_multiple_matches <- function(annotation) {
  annotation <- group_by(annotation, .data$mz)
  annotation <- mutate(annotation, multiple_match = n() > 1)
  annotation <- ungroup(annotation)
  annotation
}

#' Annotate peaks by mass matching against a compound database
#'
#' Matches peaks to compounds based on m/z values and adduct transformations,
#' then filters results using Golden Rules for molecular formula validation
#' and valid adduct pair checks.
#'
#' @param peak_table Data frame with columns: peak, mz, rt.
#' @param compound_table Data frame with columns: monoisotopic_mass,
#'   molecular_formula, compound (or compound_id), name. May also carry an
#'   integer `charge` column (default 0). Compounds with a non-zero charge are
#'   matched only as their bare ion at `monoisotopic_mass / |charge|`, under a
#'   `[M]+`-style pseudo-adduct added automatically; neutral compounds are
#'   matched across the supplied adducts but never at their neutral mass.
#' @param adduct_table Data frame with columns: adduct, charge, factor, mass.
#'   If NULL, the package default adduct table is used.
#' @param mass_tolerance Mass tolerance as a fraction (e.g., 5e-6 for 5 ppm).
#' @param charge_supplied Logical, or NULL to detect it from `compound_table`.
#'   Whether the compound table carried a `charge` column. Without one, charge
#'   is unknown and compounds keep their "M" (neutral-mass) matches as before,
#'   which is how permanently charged compounds listed at their ion mass are
#'   found. Callers that have already run [as_compound_table()], which always
#'   adds the column, must pass the value detected on the original table.
#' @return Data frame of annotations with peak, compound, adduct, and mz columns.
#' @export
#' @import dplyr
#' @importFrom rlang .data
#' @useDynLib CLUES.xMSannotator
simple_annotation <- function(peak_table, compound_table, adduct_table = NULL, mass_tolerance = 5e-6,
                              charge_supplied = NULL) {
  if (is.null(adduct_table)) {
    adduct_table <- sample_adduct_table
  }
  if (is.null(charge_supplied)) {
    charge_supplied <- is.element("charge", colnames(compound_table))
  }
  if (!charge_supplied) {
    message(
      "compound_table has no `charge` column: compounds are also matched at ",
      "their neutral mass (\"M\"), as before, which is how permanently charged ",
      "compounds listed at their ion mass are found. Add a `charge` column to ",
      "match them as bare ions and exclude neutral-mass matches."
    )
  }

  peak_table <- as_peak_table(peak_table)
  adduct_table <- as_adduct_table(adduct_table)
  compound_table <- as_compound_table(compound_table)
  adduct_table <- add_bare_ion_adducts(adduct_table, compound_table)

  annotation <- match_by_mass(peak_table, adduct_table, compound_table, mass_tolerance)
  annotation <- left_join(annotation, peak_table, by = "peak")
  annotation <- left_join(annotation, compound_table, by = "compound")

  annotation <- filter(annotation, check_golden_rules(.data$molecular_formula))
  annotation <- filter(annotation, forms_valid_adduct_pair(.data$molecular_formula, .data$adduct))
  annotation <- filter(annotation, forms_valid_charge_pair(.data$charge, .data$adduct,
                                                           allow_neutral_mass = !charge_supplied))
  # charge is an input constraint, not a result: dropping it keeps the output
  # contract stable and avoids colliding with the adduct table's own `charge`
  # in the downstream merges in get_chemscore and multilevelannotationstep4.
  annotation <- select(annotation, -any_of("charge"))
  annotation <- compute_multiple_matches(annotation)

  return(annotation)
}
