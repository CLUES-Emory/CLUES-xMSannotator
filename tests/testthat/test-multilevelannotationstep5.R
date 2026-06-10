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

test_that("same mz, different time are not collapsed", {
  df <- data.frame(
    compound_id = c("C001", "C002", "C003"),
    mz = c(100.05, 100.05, 200.10),
    time = c(60.0, 120.0, 60.0),
    Adduct = c("M+H", "M+H", "M+H"),
    score = c(50, 30, 40),
    Confidence = c(2, 1, 2),
    Module_RTclust = c("1_1", "1_1", "2_1"),
    stringsAsFactors = FALSE
  )

  result <- multilevelannotationstep5(
    outloc = tempdir(),
    chemscoremat = df
  )

  expect_equal(nrow(result), 3)
  expect_true(all(c("C001", "C002", "C003") %in% result$compound_id))
  expect_true(all(result$MatchCategory == "Unique"))
})

test_that("same mz and time collapses to best annotation", {
  df <- data.frame(
    compound_id = c("C001", "C002"),
    mz = c(100.05, 100.05),
    time = c(60.0, 60.0),
    Adduct = c("M+H", "M+H"),
    score = c(50, 30),
    Confidence = c(2, 1),
    Module_RTclust = c("1_1", "1_1"),
    stringsAsFactors = FALSE
  )

  result <- multilevelannotationstep5(
    outloc = tempdir(),
    chemscoremat = df
  )

  expect_equal(nrow(result), 1)
  expect_equal(result$compound_id, "C001")
  expect_equal(result$MatchCategory, "Unique")
})