patrick::with_parameters_test_that(
  "multilevelannotationstep5",
  {
    if (exists("skip_function") && is.function(skip_function)) {
      skip_function()
    }

    testdata_dir <- file.path("test-data", subfolder)
    load(file.path(testdata_dir, "tempobjects.Rda"))

    testthat_wd <- getwd()

    outloc <- file.path(tempdir(), "multilevelannotationstep5", subfolder)
    dir.create(outloc, recursive = TRUE)
    # Copy Stage4b data and convert to new format for the test
    stage4_data <- read.csv(file.path(testdata_dir, "Stage4.csv"))
    write.table(stage4_data, file.path(outloc, "Stage4b_confidence_levels.txt"),
                sep = "\t", row.names = FALSE)
    expected <- read.csv(file.path(testdata_dir, "Stage5.csv"))

    actual <- multilevelannotationstep5(
      outloc = outloc,
      adduct_weights = adduct_weights)

    setwd(testthat_wd)

    actual <- dplyr::arrange_all(actual)
    expected <- dplyr::arrange_all(expected)

    comparison <- dataCompareR::rCompare(
      actual,
      expected,
      keys = names(actual)
    )

    dataCompareR::saveReport(
      comparison,
      reportName = subfolder,
      reportLocation = outloc,
      showInViewer = FALSE,
      mismatchCount = 1000
    )

    expect_equal(actual, expected)
  },
  patrick::cases(
    qc_solvent = list(subfolder = "qc_solvent"),
    qc_matrix = list(subfolder = "qc_matrix", skip_function = skip_on_ci),
    batch1_neg = list(subfolder = "batch1_neg", skip_function = skip_on_ci),
    sourceforge = list(subfolder = "sourceforge", skip_function = skip_on_ci)
  )
)

# --- Unit tests for init_chemscoremat and column preservation ---

test_that("init_chemscoremat passes through data frame with NA cells", {
  df <- data.frame(
    mz = c(100.0, 200.0, 300.0),
    time = c(1.0, 2.0, 3.0),
    compound_id = c("C001", "C002", "C003"),
    Adduct = c("M+H", "M+H", "M-H"),
    score = c(50, 60, 70),
    Confidence = c(2, 3, 0),
    Confidence_Level = c("Medium", "High", "None"),
    isotopologue = c(NA, "13C:1", NA),
    isotopologue_quality = c(NA, "confirmed", NA),
    stringsAsFactors = FALSE
  )

  result <- init_chemscoremat(df)

  expect_equal(nrow(result), 3)
  expect_true("Confidence_Level" %in% names(result))
  expect_true("isotopologue" %in% names(result))
  expect_equal(result$Confidence, c(2, 3, 0))
  expect_equal(result$mz, c(100.0, 200.0, 300.0))
})

test_that("init_chemscoremat reads file when passed scalar NA", {
  outloc <- file.path(tempdir(), "test_init_chemscoremat")
  dir.create(outloc, recursive = TRUE, showWarnings = FALSE)

  df <- data.frame(
    mz = c(100.0, 200.0),
    compound_id = c("C001", "C002"),
    Adduct = c("M+H", "M-H"),
    score = c(50, 60),
    Confidence = c(2, 3),
    stringsAsFactors = FALSE
  )
  write.table(df, file.path(outloc, "Stage4b_confidence_levels.txt"),
              sep = "\t", row.names = FALSE)

  result <- init_chemscoremat(NA, outloc)

  expect_equal(nrow(result), 2)
  expect_equal(result$Confidence, c(2, 3))
  expect_equal(result$mz, c(100.0, 200.0))
})

test_that("multilevelannotationstep5 preserves extra columns", {
  outloc <- file.path(tempdir(), "test_step5_columns")
  dir.create(outloc, recursive = TRUE, showWarnings = FALSE)

  df <- data.frame(
    compound_id = c("C001", "C002", "C003"),
    mz = c(100.0, 100.0, 200.0),
    time = c(1.0, 1.0, 2.0),
    Adduct = c("M+H", "M+H", "M-H"),
    score = c(50, 60, 70),
    Confidence = c(2, 2, 3),
    Module_RTclust = c("1_1", "1_1", "2_1"),
    Confidence_Level = c("Medium", "Medium", "High"),
    isotopologue = c(NA, NA, "13C:1"),
    isotopologue_quality = c(NA, NA, "confirmed"),
    stringsAsFactors = FALSE
  )

  result <- multilevelannotationstep5(
    outloc = outloc,
    chemscoremat = df
  )

  expect_true("Confidence_Level" %in% names(result))
  expect_true("isotopologue" %in% names(result))
  expect_true("isotopologue_quality" %in% names(result))
  expect_true("MatchCategory" %in% names(result))
})