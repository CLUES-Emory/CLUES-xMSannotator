lexicographic_rank <- function(...) {
  .o <- order(...)
  .x <- cbind(...)[.o,]
  cumsum(!duplicated(.x))[order(.o)]
}

#' Get only intensity part of peak table with peaks as columns and samples on rows.
#' @param peak_table Peak table from which to extract the intensities.
#' @return Intensities with samples on rows and peaks on columns.
#' @export
get_peak_intensity_matrix <- function(peak_table) {
  peak_intensity_matrix <- t(select(peak_table, -any_of(c("peak", "mz", "rt"))))
  peak_intensity_matrix <- magrittr::set_colnames(peak_intensity_matrix, peak_table$peak)
}

#' @import dplyr
as_peak_table <- function(data, intensities = FALSE) {
  optional <- rlang::quo(any_of("peak"))
  required <- rlang::quo(all_of(c("mz", "rt")))

  if (!intensities) {
    data <- select(data, !!optional, !!required)
  }

  if (!is.element("peak", colnames(data))) {
    data$peak <- lexicographic_rank(data$mz, data$rt)
  }

  stopifnot(anyDuplicated(data$peak) == 0)
  stopifnot(is.integer(data$peak))
  for (column in colnames(data)) {
    stopifnot(is.numeric(data[[column]]))
  }

  return(data)
}

#' @import dplyr
as_adduct_table <- function(data) {
  data <- select(data, 'adduct', 'charge', 'factor', 'mass')

  stopifnot(anyDuplicated(data$adduct) == 0)
  stopifnot(is.integer(data$charge))
  stopifnot(is.integer(data$factor))
  stopifnot(is.numeric(data$mass))

  return(data)
}

#' Select columns from compound table and ensure data types.
#' @param data Compound database. Can contain either:
#'   - `compound_id` (character): User-defined compound identifiers (e.g., "HMDB0000001")
#'   - `compound` (numeric): Legacy integer compound IDs
#'   If `compound_id` is provided, an integer `compound` column is auto-generated.
#'   If only `compound` is provided, `compound_id` is created as "Formula_<compound>".
#' @return Compound table with required columns and expected types.
#' @import dplyr
as_compound_table <- function(data) {
  has_compound_id <- "compound_id" %in% names(data)

  if (has_compound_id) {
    # User provided compound_id - validate and auto-generate integer compound
    stopifnot(is.character(data$compound_id) || is.numeric(data$compound_id))
    data$compound_id <- as.character(data$compound_id)
    stopifnot(anyDuplicated(data$compound_id) == 0)

    # Auto-generate integer compound column for internal C++ processing
    data$compound <- seq_len(nrow(data))

    required <- c("monoisotopic_mass", "molecular_formula", "compound", "compound_id", "name")
  } else {
    # Legacy mode: require numeric compound column
    required <- c("monoisotopic_mass", "molecular_formula", "compound", "name")
    stopifnot(is.numeric(data$compound))
    stopifnot(anyDuplicated(data$compound) == 0)

    # Create compound_id from compound for consistency in downstream processing
    data$compound_id <- paste0("Formula_", data$compound)
  }

  data <- select(data, all_of(required))

  stopifnot(is.numeric(data$compound))
  stopifnot(is.numeric(data$monoisotopic_mass))
  stopifnot(is.character(data$molecular_formula))
  stopifnot(is.character(data$name))

  return(data)
}

#' @import dplyr
as_expected_adducts_table <- function (data) {
  data <- select(data, 'adduct')
  return(data)
}

#' Validate boosted compounds table
#' @param data Boosted compounds with columns: compound_id (required), mz (optional), rt (optional)
#' @param match_by Which columns to use for matching: c("mz"), c("rt"), or c("mz", "rt")
#' @return Validated table with ID, mz, time columns for boost_confidence_of_IDs()
#' @import dplyr
as_boosted_compounds_table <- function(data, match_by = c("mz", "rt")) {
  # Require compound_id column
  if (!"compound_id" %in% names(data)) {
    stop("boosted_compounds must contain 'compound_id' column")
  }

  # Rename compound_id to ID for internal function compatibility
  data$ID <- data$compound_id

  # Rename rt to time if present (internal function uses 'time')
  if ("rt" %in% names(data)) {
    data$time <- data$rt
  }

  # Validate required columns based on match_by
  if ("mz" %in% match_by && !"mz" %in% names(data)) {
    stop("boosted_compounds must contain 'mz' column when match_by includes 'mz'")
  }
  if ("rt" %in% match_by && !"time" %in% names(data)) {
    stop("boosted_compounds must contain 'rt' column when match_by includes 'rt'")
  }

  # Select columns - ID required, mz/time based on match_by
  cols_to_select <- c("ID")
  if ("mz" %in% match_by) cols_to_select <- c(cols_to_select, "mz")
  if ("rt" %in% match_by) cols_to_select <- c(cols_to_select, "time")

  data <- select(data, all_of(cols_to_select))

  # Validate numeric columns
  if ("mz" %in% names(data)) stopifnot(is.numeric(data$mz))
  if ("time" %in% names(data)) stopifnot(is.numeric(data$time))

  return(data)
}

#' @import dplyr
load_parquet <- function (file, columns) {
  rlang::with_handlers(
    arrow::read_parquet(file, col_select = any_of(columns)),
    error = ~ rlang::abort(paste("The file", toString(file), "does not seem to be a valid Parquet file."), parent = .)
  )
}

load_csv <- function (file, columns) {
  data <- readr::read_csv(file)
  data <- select(data, any_of(columns))
}

#' @export
load_peak_table_parquet <- function(file) {
  data <- arrow::read_parquet(file)
  as_peak_table(data, intensities = TRUE)
}

#' @export
load_adduct_table_parquet <- function(file) {
  data <- arrow::read_parquet(file)
  as_adduct_table(data)
}

#' @export
load_compound_table_parquet <- function(file) {
  data <- arrow::read_parquet(file)
  as_compound_table(data)
}

#' @export
load_expected_adducts_csv <- function (file) {
  data <- load_csv(file, columns = "adduct")
  as_expected_adducts_table(data)
}

#' @export
load_boost_compounds_csv <- function (file) {
  data <- load_csv(file, columns = c("compound", "mz", "rt"))
  as_boosted_compounds_table(data)
}

#' @export
save_parquet <- function(data, file) {
  invisible(arrow::write_parquet(data, file))
}

#' @export
create_adduct_weights <- function(adduct_weights, weight = 1) {
  if (any(is.na(adduct_weights))) {
    adduct_weights <- data.frame(
      Adduct = c("M+H", "M-H"),
      Weight = c(weight, weight)
    )
  }
  return(adduct_weights)
}

#' Validate pathway table format
#' @param data Pathway data frame with compound-pathway mappings
#' @return Validated pathway table with required columns
#' @import dplyr
#' @export
as_pathway_table <- function(data) {
  required <- c("compound", "pathway")

  if (!all(required %in% colnames(data))) {
    stop(paste("pathway_data must contain columns:", paste(required, collapse = ", ")))
  }

  data <- select(data, all_of(required), everything())
  data$compound <- as.character(data$compound)
  data$pathway <- as.character(data$pathway)

  return(data)
}
