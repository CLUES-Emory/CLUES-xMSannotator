test_that("multilevelannotationstep4 works with adduct_table = NULL", {
  chemscoremat <- data.frame(
    compound_id = "C001",
    Module_RTclust = "1_1",
    mz = 100.05,
    time = 60.0,
    MatchCategory = "Unique",
    theoretical.mz = 100.05,
    Name = "test",
    Formula = "C5H10O2",
    MonoisotopicMass = 100.05,
    Adduct = "M+H",
    mean_int_vec = 1000,
    MD = 0.05,
    score = 50,
    stringsAsFactors = FALSE
  )

  result <- multilevelannotationstep4(
    outloc = tempdir(),
    chemscoremat = chemscoremat,
    adduct_table = NULL
  )

  expect_s3_class(result, "data.frame")
  expect_true("Confidence" %in% names(result))
  expect_true(nrow(result) >= 1)
})

test_that("RT reclustering preserves module prefix without module_num", {
  # Wide RT spread triggers reclustering inside apply_rt_clustering_rules.
  # Input has Module_RTclust = "5_2" but no module_num column.
  curdata <- data.frame(
    Module_RTclust = c("5_2", "5_2", "5_2"),
    mz = c(100.0, 101.0, 102.0),
    time = c(10.0, 50.0, 90.0),
    Adduct = c("M+H", "M-H", "M+Na"),
    score = c(50, 50, 50),
    mean_int_vec = c(1000, 500, 300),
    stringsAsFactors = FALSE
  )

  adduct_weights <- data.frame(
    Adduct = c("M+H", "M-H"),
    Weight = c(1, 1),
    stringsAsFactors = FALSE
  )

  result <- CLUES.xMSannotator:::apply_rt_clustering_rules(
    curdata, c("M+H", "M-H", "M+Na"), adduct_weights,
    c("M+H", "M-H"), max_diff_rt = 30, "High"
  )

  # Module prefix "5" must be preserved, not lost to "_"
  expect_true(all(grepl("^5", result$curdata$Module_RTclust)))
  expect_false(any(grepl("^_", result$curdata$Module_RTclust)))
})