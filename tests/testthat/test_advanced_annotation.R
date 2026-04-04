patrick::with_parameters_test_that("basic advanced_annotation functionality", {
  if (exists("skip_function") && is.function(skip_function)) {
    skip_function()
  }

  testthat_wd <- getwd()
  outloc <- file.path(tempdir(), "advanced_annotation")
  dir.create(outloc, recursive = TRUE)

  expected <- read.csv(expected_output)

  peaks <- arrow::read_parquet(peak_table)
  DB <- arrow::read_parquet(database)
  adduct_table <- sample_adduct_table %>% filter(adduct %in% adduct_names)

  annotation <- advanced_annotation(
    peaks,
    DB,
    adduct_table = adduct_table,
    peak_rt_width = peak_rt_width,
    time_tolerance = time_tolerance,
    intensity_deviation_tolerance = intensity_deviation_tolerance,
    outloc = outloc
  )

  actual <- read.table(file.path(outloc, "Stage5_curated_results.txt"), sep = "\t", header = TRUE)
  setwd(testthat_wd)

  expect_equal(actual, expected)
},
  patrick::cases(
    simple_case = list(test_identifier = "simple_case",
                       expected_output = "test-data/advanced_annotation/expected_small.csv",
                       peak_table = "test-data/advanced_annotation/peak_table.parquet",
                       database = "test-data/advanced_annotation/database_small.parquet",
                       adduct_names = c("M-H", "M+H"),
                       peak_rt_width = 1,
                       time_tolerance = 10,
                       intensity_deviation_tolerance = 0.1),
    large_case = list(test_identifier = "large_case",
                      expected_output = "test-data/advanced_annotation/expected.csv",
                      peak_table = "test-data/batch1_neg.parquet",
                      database = "test-data/advanced_annotation/hmdb.parquet",
                      adduct_names = c("M-H", "2M-H", "M-2H"),
                      peak_rt_width = 15,
                      time_tolerance = 15,
                      intensity_deviation_tolerance = 0.2,
                      skip_function = skip_on_ci)
  )
)

test_that("reformat_annotation_table carries peak column through", {
  annotation <- data.frame(
    peak = c(1L, 2L, 1L),
    mz = c(100.05, 100.05, 100.05),
    rt = c(60.0, 60.0, 60.0),
    module = c(1L, 1L, 1L),
    rt_cluster = c(1L, 1L, 1L),
    multiple_match = c(TRUE, TRUE, TRUE),
    expected_mass = c(100.05, 100.05, 100.05),
    compound_id = c("C1", "C1", "C2"),
    name = c("A", "A", "B"),
    molecular_formula = c("C5H10", "C5H10", "C5H10"),
    monoisotopic_mass = c(70.0, 70.0, 70.0),
    adduct = c("M+H", "M+H", "M+H"),
    mass_number_difference = c(0L, 0L, 0L),
    mean_intensity = c(1000.0, 800.0, 1000.0),
    mass_defect = c(5L, 5L, 5L),
    stringsAsFactors = FALSE
  )

  result <- reformat_annotation_table(annotation)
  expect_true("peak" %in% names(result))
  expect_equal(result$peak, annotation$peak)
})

test_that("safe_join_feature_id uses peak key, no row duplication", {
  # Two rows share (mz, time) but have different peaks
  data <- data.frame(
    mz = c(100.05, 100.05, 200.10),
    time = c(60.0, 60.0, 120.0),
    peak = c(1L, 2L, 3L),
    compound = c("A", "B", "C"),
    stringsAsFactors = FALSE
  )

  # Feature ID map keyed on peak (each peak maps to a unique feature)
  fid_map <- data.frame(
    peak = c(1L, 2L, 3L),
    mz = c(100.05, 100.05, 200.10),
    time = c(60.0, 60.0, 120.0),
    my_fid = c("F1", "F2", "F3"),
    stringsAsFactors = FALSE
  )

  result <- safe_join_feature_id(data, fid_map, "my_fid")

  # No row duplication
  expect_equal(nrow(result), 3L)
  # Correct feature IDs assigned via peak (not ambiguous mz/time)
  expect_equal(result$my_fid, c("F1", "F2", "F3"))
})

test_that("safe_join_feature_id falls back to (mz,time) without peak", {
  # Data without peak column (post-scoring stage)
  data <- data.frame(
    mz = c(100.05, 200.10),
    time = c(60.0, 120.0),
    compound = c("A", "B"),
    stringsAsFactors = FALSE
  )

  fid_map <- data.frame(
    peak = c(1L, 2L, 3L),
    mz = c(100.05, 100.05, 200.10),
    time = c(60.0, 60.0, 120.0),
    my_fid = c("F1", "F2", "F3"),
    stringsAsFactors = FALSE
  )

  result <- safe_join_feature_id(data, fid_map, "my_fid")

  # Fallback deduplication: no row multiplication
  expect_equal(nrow(result), 2L)
  expect_true("my_fid" %in% names(result))
})
