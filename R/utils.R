#' Get the default adduct table shipped with the package.
#'
#' Returns a copy of the built-in 110-row adduct table covering positive and
#' negative ionization modes. Use this instead of accessing the internal
#' \code{sample_adduct_table} object directly.
#'
#' @return A data.frame with columns \code{adduct}, \code{charge},
#'   \code{mass}, and \code{factor}.
#' @export
get_default_adduct_table <- function() {
  as.data.frame(sample_adduct_table)
}


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

#' Name of the bare-ion pseudo-adduct for a formal charge.
#'
#' A permanently charged compound (e.g. a quaternary ammonium herbicide) is
#' already an ion and is observed at mass / |z| with no ionization adduct.
#' Names are "[M]" plus the sign, plus the magnitude when |z| > 1, so e.g.
#' 1 -> "[M]+", 2 -> "[M]+2", -5 -> "[M]-5". Returns NA for a neutral compound.
#'
#' The "[M]" prefix keeps these names clear of every pattern the pipeline
#' matches on adduct names: the isotope suffix `_[+1]`, the `_.*` strip in
#' multilevelannotationstep3, and the `[2-3]+M` multimer test in step 4.
#'
#' @param charge Integer vector of formal charges.
#' @return Character vector of pseudo-adduct names, NA where charge is 0.
#' @keywords internal
bare_ion_adduct <- function(charge) {
  ifelse(
    charge == 0,
    NA_character_,
    paste0(
      "[M]",
      ifelse(charge > 0, "+", "-"),
      ifelse(abs(charge) > 1, abs(charge), "")
    )
  )
}

#' Add the bare-ion pseudo-adducts required by permanently charged compounds.
#'
#' Appends one row per distinct non-zero compound charge, with
#' \code{mass = 0}, \code{factor = 1} and \code{charge = |z|}, so that
#' \code{match_by_mass} computes \code{monoisotopic_mass / |z|}. No-op for a
#' neutral-only database or a compound table without a \code{charge} column,
#' and idempotent.
#'
#' \code{advanced_annotation()} and \code{simple_annotation()} call this
#' automatically. Call it directly only when building an \code{adduct_weights}
#' table by hand, so the weights cover the bare-ion rows as well.
#'
#' @param adduct_table Adduct table to extend.
#' @param compound_table Compound table, optionally carrying \code{charge}.
#' @return The adduct table, with any missing bare-ion rows appended.
#' @export
add_bare_ion_adducts <- function(adduct_table, compound_table) {
  if (!is.element("charge", colnames(compound_table))) {
    return(adduct_table)
  }

  # Coerced here because this runs before as_compound_table() in the advanced
  # pipeline, so the column may still be double or character (e.g. "+2").
  charges <- suppressWarnings(as.integer(unique(compound_table$charge)))
  charges <- charges[!is.na(charges) & charges != 0]
  if (length(charges) == 0) {
    return(adduct_table)
  }

  rows <- data.frame(
    adduct = bare_ion_adduct(charges),
    charge = as.integer(abs(charges)),
    factor = 1L,
    mass = 0,
    stringsAsFactors = FALSE
  )
  rows <- rows[!is.element(rows$adduct, adduct_table$adduct), , drop = FALSE]
  if (nrow(rows) == 0) {
    return(adduct_table)
  }

  rbind(as.data.frame(adduct_table), rows)
}

#' Restrict a compound table to one ionization polarity.
#'
#' A permanently charged compound is observed only in the polarity matching its
#' own charge sign: a quaternary ammonium cation cannot be seen in negative
#' mode. Neutral compounds ionize in both and are always kept.
#'
#' With \code{charge_type = NULL} the table is returned unchanged, which is the
#' behaviour of every caller written before this argument existed. A table
#' carrying both signs then warns, because a single-polarity run would otherwise
#' silently search m/z values its instrument cannot produce.
#'
#' No-op for a compound table without a \code{charge} column, and idempotent.
#'
#' @param compound_table Compound table, optionally carrying \code{charge}.
#' @param charge_type \code{"pos"}, \code{"neg"}, or NULL for no filtering.
#' @return The compound table, keeping only compounds observable in that polarity.
#' @export
filter_compounds_by_polarity <- function(compound_table, charge_type = NULL) {
  if (!is.element("charge", colnames(compound_table))) {
    return(compound_table)
  }

  # Coerced because this runs before as_compound_table() in the advanced
  # pipeline, so the column may still be double or character (e.g. "+2").
  charge <- suppressWarnings(as.integer(compound_table$charge))
  charge[is.na(charge)] <- 0L

  if (is.null(charge_type)) {
    if (any(charge > 0) && any(charge < 0)) {
      warning(
        "compound_table contains both positively and negatively charged ",
        "compounds and charge_type was not supplied, so both polarities will ",
        "be searched. Set charge_type = 'pos' or 'neg' for a single-polarity run."
      )
    }
    return(compound_table)
  }

  keep <- if (charge_type == "pos") charge >= 0 else charge <= 0
  compound_table[keep, , drop = FALSE]
}

#' Select columns from compound table and ensure data types.
#' @param data Compound database. Can contain either:
#'   - `compound_id` (character): User-defined compound identifiers (e.g., "HMDB0000001")
#'   - `compound` (numeric): Legacy integer compound IDs
#'   If `compound_id` is provided and `compound` is also present, numeric, and
#'   unique, the supplied `compound` values are preserved. If `compound_id` is
#'   provided but `compound` is missing or invalid, an integer `compound` column
#'   is auto-generated. If only `compound` is provided, `compound_id` is created
#'   as `"Formula_<compound>"`.
#'
#'   May also carry `charge` (integer): the compound's permanent formal charge.
#'   Defaults to 0 when absent. A non-zero charge means `monoisotopic_mass` is
#'   the electron-corrected ion mass and the compound is matched only as its
#'   bare ion (see [bare_ion_adduct()]), never across ionization adducts.
#' @return Compound table with required columns and expected types.
#' @import dplyr
as_compound_table <- function(data) {
  has_compound_id <- "compound_id" %in% names(data)

  if (has_compound_id) {
    # User provided compound_id - validate and preserve or auto-generate compound
    stopifnot(is.character(data$compound_id) || is.numeric(data$compound_id))
    data$compound_id <- as.character(data$compound_id)
    stopifnot(anyDuplicated(data$compound_id) == 0)

    # Auto-generate integer compound column only when not already valid
    has_valid_compound <- "compound" %in% names(data) &&
      is.numeric(data$compound) &&
      anyDuplicated(data$compound) == 0
    if (!has_valid_compound) {
      data$compound <- seq_len(nrow(data))
    }

    required <- c("monoisotopic_mass", "molecular_formula", "compound", "compound_id", "name")
  } else {
    # Legacy mode: require numeric compound column
    required <- c("monoisotopic_mass", "molecular_formula", "compound", "name")
    stopifnot(is.numeric(data$compound))
    stopifnot(anyDuplicated(data$compound) == 0)

    # Create compound_id from compound for consistency in downstream processing
    data$compound_id <- paste0("Formula_", data$compound)
  }

  # Optional formal charge. Non-zero means the compound is already an ion and is
  # matched only as its bare ion; absent or NA means neutral.
  if (!is.element("charge", colnames(data))) {
    data$charge <- 0L
  } else {
    stopifnot(is.numeric(data$charge))
    data$charge <- as.integer(data$charge)
    data$charge[is.na(data$charge)] <- 0L
  }
  required <- c(required, "charge")

  data <- select(data, all_of(required))

  stopifnot(is.numeric(data$compound))
  stopifnot(is.numeric(data$monoisotopic_mass))
  stopifnot(is.character(data$molecular_formula))
  stopifnot(is.character(data$name))
  stopifnot(is.integer(data$charge))

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

#' Load peak table from a Parquet file
#'
#' Reads a Parquet file and validates it as a peak table with required columns.
#'
#' @param file Path to the Parquet file.
#' @return Validated peak table data frame with peak, mz, rt, and intensity columns.
#' @export
load_peak_table_parquet <- function(file) {
  data <- arrow::read_parquet(file)
  as_peak_table(data, intensities = TRUE)
}

#' Load adduct table from a Parquet file
#'
#' Reads a Parquet file and validates it as an adduct table.
#'
#' @param file Path to the Parquet file.
#' @return Validated adduct table data frame with adduct, charge, factor, mass columns.
#' @export
load_adduct_table_parquet <- function(file) {
  data <- arrow::read_parquet(file)
  as_adduct_table(data)
}

#' Load compound table from a Parquet file
#'
#' Reads a Parquet file and validates it as a compound table.
#'
#' @param file Path to the Parquet file.
#' @return Validated compound table data frame.
#' @export
load_compound_table_parquet <- function(file) {
  data <- arrow::read_parquet(file)
  as_compound_table(data)
}

#' Load expected adducts from a CSV file
#'
#' Reads a CSV file and extracts the adduct column.
#'
#' @param file Path to the CSV file containing an adduct column.
#' @return Data frame with a single adduct column.
#' @export
load_expected_adducts_csv <- function (file) {
  data <- load_csv(file, columns = "adduct")
  as_expected_adducts_table(data)
}

#' Load boost compounds from a CSV file
#'
#' Reads a CSV file with compound identifiers for confidence boosting.
#'
#' @param file Path to the CSV file with compound_id, mz, and rt columns.
#' @return Validated boosted compounds table.
#' @export
load_boost_compounds_csv <- function (file) {
  data <- load_csv(file, columns = c("compound_id", "mz", "rt"))
  as_boosted_compounds_table(data)
}

#' Save data frame to a Parquet file
#'
#' Writes a data frame to disk in Apache Parquet format.
#'
#' @param data Data frame to save.
#' @param file Path to the output Parquet file.
#' @return Invisible NULL.
#' @export
save_parquet <- function(data, file) {
  invisible(arrow::write_parquet(data, file))
}

#' Create default adduct weights table
#'
#' Returns the provided adduct weights, or creates a default table with M+H and
#' M-H if the input is NULL or NA.
#'
#' @param adduct_weights Existing adduct weights data frame, or NULL/NA.
#' @param weight Default weight value for M+H and M-H adducts.
#' @return Data frame with Adduct and Weight columns.
#' @export
create_adduct_weights <- function(adduct_weights, weight = 1) {
  if (is.null(adduct_weights) || any(is.na(adduct_weights))) {
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

#' Choose one representative peak per compound, on evidence rather than row order
#'
#' `advanced_annotation()` scores each compound against a single query row, and
#' that row's `time` becomes the centre of the +/- `max_diff_rt` window that
#' [get_chemscore()] scores within. The choice therefore decides which
#' chromatographic peak a compound can be annotated to.
#'
#' Previously the row was taken by `dplyr::distinct(compound_id, .keep_all =
#' TRUE)`, i.e. whichever row happened to come first. Row order originates in
#' `match_by_mass.cpp`, whose loop nesting is adduct index -> compound index ->
#' peak index, so the anchor followed the order of the adduct table and
#' reordering that table changed which peaks were annotated.
#'
#' Candidates are the monoisotopic rows, since [get_chemscore()] returns NULL
#' when its window contains none. Each candidate is ranked on the evidence
#' co-located with it -- same `compound_id`, same module prefix, and within
#' `max_diff_rt` seconds:
#'
#' 1. isotope support: how many co-eluting isotope rows there are
#' 2. ratio fidelity: the smallest `isotope_ratio_dev` among them
#' 3. primary adduct: whether the base adduct is in `level1_primary_adducts`
#' 4. abundance: `mean_int_vec`
#' 5. `mz` then `time`, so the result never depends on input row order
#'
#' Criteria 1-2 are what make the choice evidence-based; criterion 5 is what
#' makes it deterministic. Output rows are ordered by `compound_id` for the same
#' reason.
#'
#' `isotope_ratio_dev` is read only from isotope-suffixed rows. Monoisotopic
#' rows can carry a non-NA value -- [reformat_annotation_table()] leaves the
#' adduct unsuffixed when `mass_number_difference` is 0 -- and such a value
#' does not describe the row it sits on.
#'
#' @param annotation Annotation table after [reformat_annotation_table()], with
#'   columns `compound_id`, `Adduct`, `time`, `mz`, `mean_int_vec`,
#'   `Module_RTclust`, and optionally `isotope_ratio_dev`.
#' @param level1_primary_adducts Adducts treated as primary evidence.
#' @param max_diff_rt Retention time window in seconds, matching the window
#'   [get_chemscore()] applies.
#' @return One row per `compound_id`, ordered by `compound_id`.
#' @export
select_query_anchor <- function(annotation,
                                level1_primary_adducts = c("M+H", "M-H"),
                                max_diff_rt = 10) {
  if (nrow(annotation) == 0) {
    return(annotation)
  }

  isotope_pattern <- "_\\[(\\+|\\-)[0-9]+\\]"
  is_isotope <- grepl(isotope_pattern, annotation$Adduct)
  base_adduct <- sub(isotope_pattern, "", annotation$Adduct)
  # same expression as check_compound_coherence(), so "co-located" means the
  # same thing here as it does downstream
  module <- gsub(as.character(annotation$Module_RTclust),
                 pattern = "_[0-9]*", replacement = "")
  rt <- as.numeric(annotation$time)
  mz <- as.numeric(annotation$mz)
  intensity <- as.numeric(annotation$mean_int_vec)
  ratio_dev <- if ("isotope_ratio_dev" %in% names(annotation)) {
    as.numeric(annotation$isotope_ratio_dev)
  } else {
    rep(NA_real_, nrow(annotation))
  }

  rows_by_compound <- split(seq_len(nrow(annotation)),
                            as.character(annotation$compound_id))

  chosen <- vapply(rows_by_compound, function(idx) {
    candidates <- idx[!is_isotope[idx]]
    # a compound with only isotope rows cannot be scored; keep one anyway so the
    # caller still sees it, as distinct() did
    if (length(candidates) == 0) candidates <- idx
    if (length(candidates) == 1) return(candidates)

    isotopes <- idx[is_isotope[idx]]
    support <- integer(length(candidates))
    fidelity <- rep(Inf, length(candidates))

    for (k in seq_along(candidates)) {
      i <- candidates[k]
      co_eluting <- isotopes[module[isotopes] == module[i] &
                               abs(rt[isotopes] - rt[i]) <= max_diff_rt]
      support[k] <- length(co_eluting)
      observed <- ratio_dev[co_eluting]
      observed <- observed[!is.na(observed)]
      if (length(observed) > 0) fidelity[k] <- min(observed)
    }

    primary <- as.integer(base_adduct[candidates] %in% level1_primary_adducts)
    ranked <- order(-support, fidelity, -primary, -intensity[candidates],
                    mz[candidates], rt[candidates])
    candidates[ranked[1]]
  }, integer(1))

  annotation[chosen, , drop = FALSE]
}
