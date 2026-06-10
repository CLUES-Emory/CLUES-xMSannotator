patrick::with_parameters_test_that("load_compound_table_parquet:", {
  filename <- file.path("test-data", "utils", paste0(db_name, ".parquet"))
  actual <- load_compound_table_parquet(filename)

  expected <- arrow::read_parquet(
    filename,
    col_select = c("monoisotopic_mass", "molecular_formula", "compound", "name")
  )
  expect_equal(actual, expected)
},
  db_name = c("hmdb", "pubchem")
)

test_that("Utils data reading: peak table processing.", {
  peak_table <- readRDS("test-data/utils/peak_table.rds")
  actual <- as_peak_table(peak_table, intensities = FALSE)
  peak_table <- select(peak_table, c("peak", "mz", "rt"))

  expect_equal(actual, peak_table)
})

test_that("load_boost_compounds_csv reads compound_id,mz,rt CSV", {
  tmpf <- tempfile(fileext = ".csv")
  on.exit(unlink(tmpf))
  writeLines("compound_id,mz,rt\nHMDB0000001,100.05,60.2\nHMDB0000002,200.10,120.5", tmpf)

  result <- load_boost_compounds_csv(tmpf)

  expect_s3_class(result, "data.frame")
  expect_true("ID" %in% names(result))
  expect_equal(result$ID, c("HMDB0000001", "HMDB0000002"))
  expect_equal(result$mz, c(100.05, 200.10))
  expect_equal(result$time, c(60.2, 120.5))
})

test_that("as_compound_table preserves caller-supplied numeric compound IDs", {
  data <- data.frame(
    compound_id = c("HMDB001", "HMDB002", "HMDB003"),
    compound = c(101L, 202L, 303L),
    monoisotopic_mass = c(100.05, 200.10, 300.15),
    molecular_formula = c("C5H10O2", "C10H20O4", "C15H30O6"),
    name = c("A", "B", "C"),
    stringsAsFactors = FALSE
  )

  result <- as_compound_table(data)

  expect_equal(result$compound, c(101L, 202L, 303L))
  expect_equal(result$compound_id, c("HMDB001", "HMDB002", "HMDB003"))
})

test_that("Utils data reading: bad peak table processing.", {
  peak_table <- readRDS("test-data/utils/peak_table.rds")

  expect_error(as_peak_table(peak_table, intensities = TRUE), "is\\.numeric\\(data\\[\\[.*\\]\\]\\).*not TRUE")
})
