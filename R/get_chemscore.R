
#' Compute chemical score for a single compound
#'
#' Filters annotation data for a given compound within a retention time window,
#' separates isotopic from monoisotopic peaks, scores monoisotopic peaks, and
#' applies a 100x boost when isotope evidence is present.
#'
#' @param ... Named arguments passed as a query row (must include compound_id, time).
#' @param annotation Data frame of annotations from simple_annotation.
#' @param corthresh Correlation threshold for peak module filtering.
#' @param global_cor Global correlation matrix between peaks.
#' @param max_diff_rt Maximum retention time difference for grouping (seconds).
#' @param adduct_weights Data frame with adduct names and their weights.
#' @param filter.by Character vector of expected adducts.
#' @param MplusH.abundance.ratio.check Logical; if TRUE, check that secondary
#'   adducts have lower intensity than the primary M+H or M-H adduct.
#' @param adduct_table Data frame with adduct definitions.
#' @return Data frame with chemical scores, or NULL if no valid data.
#' @import tidyr
#' @import dplyr
#' @importFrom magrittr %>%
#' @export
get_chemscore <- function(...,
                          annotation,
                          corthresh,
                          global_cor,
                          max_diff_rt = 10,
                          adduct_weights,
                          filter.by,
                          MplusH.abundance.ratio.check = TRUE,
                          adduct_table) {

  query <- tibble::tibble(...)

  # Filter to this compound_id within RT window (FIX: use max_diff_rt, was hardcoded to 10)
  curmchemdata <- dplyr::filter(
    annotation,
    compound_id == query$compound_id &
      abs(time - query$time) <= max_diff_rt
  )

  if (length(curmchemdata$mz) < 1) stop("No mz data found!")

  # ============================================
  # Separate isotopes before scoring
  # ============================================
  isotope_pattern <- "_\\[(\\+|\\-)[0-9]+\\]"
  isotope_rows <- curmchemdata %>% filter(grepl(isotope_pattern, Adduct))
  monoisotopic_rows <- curmchemdata %>% filter(!grepl(isotope_pattern, Adduct))

  has_isotopes <- nrow(isotope_rows) > 0

  # Handle edge case: only isotopes, no monoisotopic peaks
  if (nrow(monoisotopic_rows) < 1) {
    return(NULL)
  }

  # Score monoisotopic peaks only
  result <- compute_chemical_score(
    mchemicaldata = monoisotopic_rows,
    adduct_weights = adduct_weights,
    global_cor = global_cor,
    corthresh = corthresh,
    filter.by = filter.by,
    max_diff_rt = max_diff_rt,
    chemicalid = query$compound_id,
    MplusH.abundance.ratio.check = MplusH.abundance.ratio.check,
    adduct_table = adduct_table
  )

  # Return NULL if no valid monoisotopic data
  if (is.null(result$filtdata) || nrow(result$filtdata) < 1) {
    return(NULL)
  }

  # ============================================
  # Apply 100x isotope boost if isotopes detected
  # ============================================
  chemical_score <- result$chemical_score
  if (has_isotopes) {
    chemical_score <- chemical_score * 100
  }

  result$filtdata <- result$filtdata[order(result$filtdata$mz), ]
  cur_chem_score <- rep_len(chemical_score, nrow(result$filtdata))
  chemscoremat <- cbind(cur_chem_score, result$filtdata)

  # ============================================
  # Re-add isotope rows with boosted score
  # ============================================
  if (has_isotopes && nrow(isotope_rows) > 0) {
    # Match isotopes to their parent adducts that survived filtering
    surviving_adducts <- unique(result$filtdata$Adduct)

    isotope_rows_matched <- isotope_rows %>%
      mutate(parent_adduct = gsub(isotope_pattern, "", Adduct)) %>%
      filter(parent_adduct %in% surviving_adducts) %>%
      select(-parent_adduct)

    if (nrow(isotope_rows_matched) > 0) {
      # Add score column to isotope rows
      isotope_rows_matched$cur_chem_score <- chemical_score

      # Ensure column alignment and bind
      common_cols <- intersect(names(chemscoremat), names(isotope_rows_matched))
      chemscoremat <- dplyr::bind_rows(
        chemscoremat[, common_cols],
        isotope_rows_matched[, common_cols]
      )
    }
  }

  # Filter on critical columns only
  critical_cols <- c("mz", "time", "compound_id", "Adduct")
  chemscoremat <- chemscoremat[complete.cases(chemscoremat[, critical_cols]), ]

  if (nrow(chemscoremat) < 1) {
    return(NULL)
  }

  return(chemscoremat)
}
