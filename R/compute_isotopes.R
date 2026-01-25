#' Compute isotopic pattern of a given molecule.
#'
#' @param formula A string containing molecular or empirical formula of a compound.
#' @param minAbund A cutoff abundance ratio: only isotopes with abundance ratio above this value are returned.
#'
#' @return A dataframe containing isotopic pattern of a given molecule. The pattern is represented by three columns:
#' mass, mass number difference, and abundance. Abundance of the isotopes is normalized with the most abundant isotope
#' having a unit value and abundances of the other isotopes are represented as a share of the most abundant isotope.
#' Mass number difference is 0 for the most abundant isotope.
#'
#' @export
#' @import dplyr
#' @importFrom rcdk get.formula get.isotopes.pattern
compute_isotopic_pattern <- function(formula, minAbund = 0.001) {
  formula <- get.formula(formula)
  isotopes <- get.isotopes.pattern(formula, minAbund)
  isotopes <- as_tibble(isotopes)
  isotopes <- arrange(isotopes, desc(abund))
  isotopes <- mutate(isotopes,
    mass_number_difference = round(mass - mass[1]),
    exact_mass_diff = mass - mass[1]
  )
}

#' Match isotopes from peak table to annotated peaks based on rt cluster, rt difference, and mass defect.
#'
#' @param query A single annotated peak.
#' @param intensity_deviation_tolerance A numeric threshold by which an intensity ratio of two isotopic peaks may differ
#'  from their actual abundance ratio.
#' @param peaks A peak table containing a peak identifier (unique number), mean intensity, module, and rt cluster
#'  of each identified peak.
#' @param mass_defect_tolerance A number. Maximum difference in mass defect between two peaks of the same compound.
#' @param rt_tolerance A number. Maximum rt difference for two peaks of the same substance.
#'
#' @return A table of identified isotopes with an extra column: `mass_number_difference`.
#'  The extra column represents a nominal mass (mass number) difference between a given isotope and the most abundant
#'  isotope of the same molecule.
#'
#' @import dplyr
filter_isotopes <- function(query,
                            intensity_deviation_tolerance,
                            peaks,
                            mass_defect_tolerance,
                            rt_tolerance) {
  isotopes <- mutate(
    peaks,
    mass_number_difference = round(mz - query$expected_mass),
    compound = query$compound,
    compound_id = if ("compound_id" %in% names(query)) query$compound_id else NA_character_,
    adduct = query$adduct,
    molecular_formula = query$molecular_formula,
    name = if ("name" %in% names(query)) query$name else NA_character_,
    monoisotopic_mass = if ("monoisotopic_mass" %in% names(query)) query$monoisotopic_mass else NA_real_
  )
  isotopes <- filter(
    isotopes,
    peak != query$peak,
    rt_cluster == query$rt_cluster,
    near(rt, query$rt, rt_tolerance),
    near(mass_defect, query$mass_defect, mass_defect_tolerance)
  )
}

#' Match prefiltered isotopes by comparing their relative intensity to their relative natural abundance.
#' The relative values are computed as a ratio of an isotopic peak intensity to the intensity of its monoisotopic peak,
#' respectively as a ratio of a natural abundance of the isotope to the abundance of the most abundant isotope.
#'
#' @param query A single annotated peak.
#' @param isotopes Isotope candidates for the annotated molecule. The candidates are identified in [compute_isotopes].
#' @param pattern Isotopic pattern of a given molecule as produced by [compute_isotopic_pattern].
#' @param intensity_deviation_tolerance A numeric threshold by which an intensity ratio of two isotopic peaks may differ
#'  from their actual abundance ratio.
#' @param isotope_mass_tolerance_ppm Maximum allowed mass error in ppm between observed and theoretical isotope m/z.
#'  If NULL, no ppm filter is applied.
#'
#' @return A table of matched isotopes.
#'
#' @import dplyr
match_isotopes_by_intensity <- function(query,
                                        isotopes,
                                        pattern,
                                        intensity_deviation_tolerance,
                                        isotope_mass_tolerance_ppm = NULL) {
  isotopes <- mutate(isotopes,
    relative_intensity = mean_intensity / query$mean_intensity
  )
  isotopes <- left_join(isotopes,
    select(pattern, mass_number_difference, abund, exact_mass_diff),
    by = "mass_number_difference"
  )

  # Calculate expected isotope mass
  isotopes <- mutate(isotopes,
    expected_mass = query$expected_mass + exact_mass_diff,
    mass_error_ppm = abs(mz - expected_mass) / expected_mass * 1e6
  )

  # Apply ppm mass accuracy filter if specified
  if (!is.null(isotope_mass_tolerance_ppm)) {
    isotopes <- filter(isotopes, mass_error_ppm <= isotope_mass_tolerance_ppm)
  }

  # Apply intensity ratio filter
  isotopes <- filter(
    isotopes,
    near(relative_intensity, abund, relative_intensity * intensity_deviation_tolerance)
  )

  isotopes <- select(isotopes, -c(abund, relative_intensity, exact_mass_diff, mass_error_ppm))
}

#' Find all possible isotopes of a given annotated peak by filtering peaks based on several criteria
#' (see [filter_isotopes]) and matching the intensities with natural abundances of the isotopes.
#'
#' @param ... A single annotated feature passed as a list of columns. See [purrr::pmap()] for more details.
#' @param intensity_deviation_tolerance A numeric threshold by which an intensity ratio of two isotopic peaks may differ
#'  from their actual abundance ratio.
#' @param peaks A peak table containing a peak identifier (unique number), mean intensity, module, and rt cluster
#'  of each identified peak.
#' @param mass_defect_tolerance A number. Maximum difference in mass defect between two peaks of the same compound.
#' @param rt_tolerance A number. Maximum rt difference for two peaks of the same substance.
#' @param isotope_mass_tolerance_ppm Maximum allowed mass error in ppm between observed and theoretical isotope m/z.
#'  If NULL, no ppm filter is applied.
#'
#' @return A table with peaks that have been identified as isotopes of a given molecule from the annotation table.
#'
#' @import dplyr
#' @importFrom rlang .data
detect_isotopic_peaks <- function(...,
                                  intensity_deviation_tolerance,
                                  peaks,
                                  mass_defect_tolerance,
                                  rt_tolerance,
                                  isotope_mass_tolerance_ppm = NULL) {
  query <- tibble(...)
  isotopic_pattern <- compute_isotopic_pattern(query$molecular_formula)

  isotopes <- filter_isotopes(
    query,
    intensity_deviation_tolerance,
    peaks,
    mass_defect_tolerance,
    rt_tolerance
  )
  isotopes <- distinct(isotopes)
  isotopes <- match_isotopes_by_intensity(
    query,
    isotopes,
    isotopic_pattern,
    intensity_deviation_tolerance,
    isotope_mass_tolerance_ppm
  )
}

#' For each annotated feature find all possible isotopes from peak table.
#'
#' @param annotation A table with annotated peaks.
#' @param adduct_weights A weight-by-adduct table.
#' @param intensity_deviation_tolerance A numeric threshold by which an intensity ratio of two isotopic peaks may differ
#'  from their actual abundance ratio. (Cl35 1.0 and Cl37 0.95 -> absolute abundances are about 0.45 and 0.42)
#' @param mass_defect_tolerance A number. Maximum difference in mass defect between two peaks of the same compound.
#' @param peak_table A peak table containing a peak identifier (unique number), mean intensity, module, and rt cluster
#'  of each identified peak.
#' @param rt_tolerance A number. Maximum rt difference for two peaks of the same substance.
#' @param isotope_mass_tolerance_ppm Maximum allowed mass error in ppm between observed and theoretical isotope m/z.
#'  If NULL, no ppm filter is applied. Recommended: 5-10 ppm for high-resolution MS data.
#'
#' @return Annotation table expanded by annotated isotopic peaks.
#'
#' @import dplyr
#' @importFrom purrr pmap_dfr
compute_isotopes <- function(annotation,
                             adduct_weights,
                             intensity_deviation_tolerance = 0.1,
                             mass_defect_tolerance = 0,
                             peak_table,
                             rt_tolerance = 1,
                             isotope_mass_tolerance_ppm = NULL) {
  annotation <- mutate(annotation, mass_number_difference = 0)
  adducts <- semi_join(annotation, adduct_weights, by = "adduct")

  # This can be parallelized on `group_split(group_by(isotopes, molecular_formula))`
  isotopes <- purrr::pmap_dfr(
    adducts,
    ~ detect_isotopic_peaks(...,
        peaks = peak_table,
        rt_tolerance = rt_tolerance,
        intensity_deviation_tolerance = intensity_deviation_tolerance,
        mass_defect_tolerance = mass_defect_tolerance,
        isotope_mass_tolerance_ppm = isotope_mass_tolerance_ppm
    )
  )
  annotation <- bind_rows(annotation, isotopes)
}
