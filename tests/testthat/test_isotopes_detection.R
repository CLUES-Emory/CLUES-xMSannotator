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
  expect_equal(actual, expected)
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

  # 21 before the fix, 37 after
  expect_equal(sum(actual$mass_number_difference != 0), 37)
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
    expect_equal(actual, expected, ignore_attr = TRUE)
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
  expect_equal(actual, expected)
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
