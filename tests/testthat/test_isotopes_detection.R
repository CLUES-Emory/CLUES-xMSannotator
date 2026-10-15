# `isotope_ratio_dev` is a deliberate addition that postdates the stored
# fixtures. Assert it on its own terms, then compare every other column against
# the fixture unchanged -- regenerating the fixtures instead would hide any
# unrelated drift they still guard against.
expect_ratio_dev <- function(actual) {
  expect_true("isotope_ratio_dev" %in% names(actual))
  expect_type(actual$isotope_ratio_dev, "double")
  observed <- actual$isotope_ratio_dev[!is.na(actual$isotope_ratio_dev)]
  expect_true(all(observed >= 0))
}

drop_ratio_dev <- function(actual) {
  actual[, setdiff(names(actual), "isotope_ratio_dev"), drop = FALSE]
}

test_that("Isotopes-computing: isotopes annotation.", {
  peak_table <- readRDS("test-data/isotopes_detection/peaks.rds")
  annotation <- readRDS("test-data/isotopes_detection/annotation.rds")
  adduct_weights <- tibble(adduct = c("M+H", "M+NH4"), weight = c(5, 5))

  actual <- compute_isotopes(
    annotation = annotation,
    adduct_weights = adduct_weights,
    intensity_deviation_tolerance = 0.25,
    peak_table = peak_table,
    rt_tolerance = 10,
    mass_defect_tolerance = 0.1
  )
  expected <- readRDS("test-data/isotopes_detection/annotation_with_isp.rds")
  expect_ratio_dev(actual)
  expect_equal(drop_ratio_dev(actual), expected)
})

test_that("isotopes near a mass-defect bin boundary are not discarded", {
  # The old prefilter compared binned mass-defect indices, so an isotope whose
  # m/z fell in the neighbouring 0.01 bin was rejected: the 13C shift moves the
  # fractional part by 0.00336 and crossed a boundary about a third of the time.
  # These counts are at the 5 ppm accuracy the pipeline actually runs with, so
  # every recovered peak also satisfies the ppm and intensity-ratio checks.
  peak_table <- readRDS("test-data/isotopes_detection/peaks.rds")
  annotation <- readRDS("test-data/isotopes_detection/annotation.rds")
  adduct_weights <- tibble(adduct = c("M+H", "M+NH4"), weight = c(5, 5))

  actual <- compute_isotopes(
    annotation = annotation,
    adduct_weights = adduct_weights,
    intensity_deviation_tolerance = 0.25,
    peak_table = peak_table,
    rt_tolerance = 10,
    mass_defect_tolerance = 0.1,
    isotope_mass_tolerance_ppm = 5
  )

  # 21 before the mass-defect fix, 37 after. Now 43: the intensity band was
  # made symmetric about the theoretical abundance, which admits isotopes
  # reading weak (ratio 0.7-0.769 of theory) that the observed-scaled band
  # rejected. Both fixes add real isotopes; neither relaxes the tolerance.
  expect_equal(sum(actual$mass_number_difference != 0), 43)
})

patrick::with_parameters_test_that(
  "Isotopes-computing: isotopes detection (single annotation).",
  {
    peaks <- readRDS("test-data/isotopes_detection/peaks.rds")
    actual <- detect_isotopic_peaks(
      peaks = peaks,
      intensity_deviation_tolerance = 0.25,
      rt_tolerance = 10,
      mass_defect_tolerance = 0.1,
      query
    )
    expect_ratio_dev(actual)
    expect_equal(drop_ratio_dev(actual), expected, ignore_attr = TRUE)
  },
  patrick::cases(
    isoproturon = list(
      query = readRDS("test-data/isotopes_detection/isoproturon.rds"),
      expected = readRDS("test-data/isotopes_detection/isoproturon_isp.rds")
    ),
    desisopropylatrazine = list(
      query = readRDS("test-data/isotopes_detection/desisopropylatrazine.rds"),
      expected = readRDS("test-data/isotopes_detection/desisopropylatrazine_isp.rds")
    ),
    acetochlor = list(
      query = readRDS("test-data/isotopes_detection/acetochlor.rds"),
      expected = readRDS("test-data/isotopes_detection/acetochlor_isp.rds")
    )
  )
)

test_that("Isotopes-computing: isotope matching by intensity.", {
  query <- readRDS("test-data/isotopes_detection/query.Rda")
  isotopes <- readRDS("test-data/isotopes_detection/isotopes.Rda")
  pattern <- readRDS("test-data/isotopes_detection/adenine_pattern.rda")

  expected <- readRDS("test-data/isotopes_detection/matched_isotopes.rda")
  actual <- match_isotopes_by_intensity(
    query = query,
    isotopes = isotopes,
    pattern = pattern,
    intensity_deviation_tolerance = 0.4
  )
  expect_ratio_dev(actual)
  expect_equal(drop_ratio_dev(actual), expected)
})

patrick::with_parameters_test_that(
  "Isotopes-computing: theoretical isotopic pattern",
  {
    actual <- compute_isotopic_pattern(formula, minAbund = 0.001)
    expect_equal(actual, expected)
  },
  patrick::cases(
    carbon = list(
      formula = "C",
      expected = readRDS("test-data/isotopes_detection/C_pattern.rda")
    ),
    chlorophorm = list(
      formula = "CHCl3",
      expected = readRDS("test-data/isotopes_detection/CHCl3_pattern.rda")
    ),
    chlorophyll = list(
      formula = "C55H72O5N4Mg",
      expected = readRDS("test-data/isotopes_detection/CHLA_pattern.rda")
    )
  )
)

# --- intensity acceptance band ---
# The tolerance is a fraction of the theoretical abundance, so the accepted
# band is symmetric about it: [1 - tol, 1 + tol] in observed/theoretical terms.
# Scaling by the observed ratio instead gave [1/(1+tol), 1/(1-tol)] =
# [0.769, 1.429] at tol = 0.3 -- centred on 1.099, admitting isotopes reading
# strong (interference) and rejecting cleanly measured weak ones.

# accepts a pair whose observed/theoretical ratio is `ratio`
accepted_at <- function(ratio, abund = 0.644, tol = 0.3) {
  query <- tibble(mean_intensity = 1000, expected_mass = 160.9567)
  pattern <- tibble(mass_number_difference = 2, abund = abund,
                    exact_mass_diff = 1.99705)
  isotopes <- tibble(
    mz = query$expected_mass + pattern$exact_mass_diff,
    mean_intensity = 1000 * abund * ratio,
    mass_number_difference = 2
  )
  nrow(match_isotopes_by_intensity(query, isotopes, pattern, tol)) == 1
}

test_that("band is symmetric about the theoretical abundance", {
  expect_true(accepted_at(1.00))
  expect_true(accepted_at(0.75))
  expect_true(accepted_at(1.25))
  # edges
  expect_true(accepted_at(0.71))
  expect_true(accepted_at(1.29))
  expect_false(accepted_at(0.69))
  expect_false(accepted_at(1.31))
})

test_that("a weak but clean M+2 is accepted (was rejected by the observed-scaled band)", {
  # 2,4-dichlorophenol, C18neg rt 123.4: observed 0.463 against 0.644 theory,
  # i.e. ratio 0.719 -- inside [0.7, 1.3] but below the old 0.769 floor.
  expect_true(accepted_at(0.719))
})

test_that("an M+2 reading far stronger than theory is rejected", {
  # ratio 1.40 sat inside the old band, which stretched to 1.429
  expect_false(accepted_at(1.40))
})

test_that("the band scales with the tolerance argument", {
  expect_true(accepted_at(0.55, tol = 0.5))
  expect_false(accepted_at(0.45, tol = 0.5))
  expect_true(accepted_at(0.92, tol = 0.1))
  expect_false(accepted_at(0.88, tol = 0.1))
})
