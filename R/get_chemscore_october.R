
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

  curmchemdata <- dplyr::filter(
    annotation,
    chemical_ID == query$chemical_ID &
    abs(time - query$time) <= 10
  )

  if (length(curmchemdata$mz) < 1) stop("No mz data found!")

  result <- compute_chemical_score(
    mchemicaldata = curmchemdata,
    adduct_weights = adduct_weights,
    global_cor = global_cor,
    corthresh = corthresh,
    filter.by = filter.by,
    max_diff_rt = max_diff_rt,
    chemicalid = query$chemical_ID,
    MplusH.abundance.ratio.check = MplusH.abundance.ratio.check
  )

  # Return NULL if no valid data (pmap_dfr will skip NULL values)
  if (is.null(result$filtdata) || nrow(result$filtdata) < 1) {
    return(NULL)
  }

  result$filtdata <- result$filtdata[order(result$filtdata$mz), ]
  cur_chem_score <- rep_len(result$chemical_score, nrow(result$filtdata))
  chemscoremat <- cbind(cur_chem_score, result$filtdata)

  # Remove rows only if critical columns are NA (not isotope-specific columns)
  # Isotopes have NA in theoretical.mz, Name, MonoisotopicMass - which is expected
  critical_cols <- c("mz", "time", "chemical_ID", "Adduct")
  chemscoremat <- chemscoremat[complete.cases(chemscoremat[, critical_cols]), ]

  # Return NULL if no rows remain after filtering
  if (nrow(chemscoremat) < 1) {
    return(NULL)
  }

  return(chemscoremat)
}
