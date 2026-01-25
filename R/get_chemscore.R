
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
                          mass_defect_window = 0.01,
                          MplusH.abundance.ratio.check = TRUE) {

  query <- tibble::tibble(...)

  # Filter to this chemical_ID within RT window (FIX: use max_diff_rt, was hardcoded to 10)
  curmchemdata <- dplyr::filter(
    annotation,
    chemical_ID == query$chemical_ID &
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
    chemicalid = query$chemical_ID,
    MplusH.abundance.ratio.check = MplusH.abundance.ratio.check
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
  critical_cols <- c("mz", "time", "chemical_ID", "Adduct")
  chemscoremat <- chemscoremat[complete.cases(chemscoremat[, critical_cols]), ]

  if (nrow(chemscoremat) < 1) {
    return(NULL)
  }

  return(chemscoremat)
}
