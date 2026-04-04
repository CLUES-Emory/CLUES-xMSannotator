patrick::with_parameters_test_that(
  "multilevelannotation step 4 works",
  {
    if (exists("skip_function") && is.function(skip_function)) {
      skip_function()
    }

    # load data needed during step 4
    testdata_dir <- file.path(getwd(), "test-data", subfolder)
    load(file.path(testdata_dir, "tempobjects.Rda"))

    testthat_wd <- getwd()
    outloc <- file.path(
      tempdir(),
      "multilevelannotationstep4",
      subfolder
    )

    # create test folder
    dir.create(outloc, recursive = TRUE)
    
    chemscoremat <- read.csv(file.path(testdata_dir, "Stage3.csv"))

    # load expected results
    expected <- read.csv(file.path(testdata_dir, "Stage4.csv"))

    # compute annotation step 4
    result <- multilevelannotationstep4(
      outloc = outloc,
      chemscoremat = chemscoremat,
      max.mz.diff = max.mz.diff,
      max.rt.diff = max_diff_rt,
      filter.by = filter.by,
      adduct_weights = adduct_weights,
      max_isp = max_isp,
      min_ions_perchem = min_ions_perchem
    )
    actual <- read.table(file.path(outloc, "Stage4_confidence_levels.txt"), sep = "\t", header = TRUE)

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