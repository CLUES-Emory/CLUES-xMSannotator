# Tests for check_compound_coherence()

test_that("Single row is trivially coherent", {
  curdata <- data.frame(
    Module_RTclust = "73_3",
    time = 100,
    stringsAsFactors = FALSE
  )
  expect_true(check_compound_coherence(curdata, max.rt.diff = 10))
})

test_that("Same module, RT within threshold -> coherent", {
  curdata <- data.frame(
    Module_RTclust = c("73_3", "73_5", "73_3"),
    time = c(100, 105, 102),
    stringsAsFactors = FALSE
  )
  expect_true(check_compound_coherence(curdata, max.rt.diff = 10))
})

test_that("Different modules -> incoherent", {
  curdata <- data.frame(
    Module_RTclust = c("73_3", "85_2", "73_3"),
    time = c(100, 102, 101),
    stringsAsFactors = FALSE
  )
  expect_false(check_compound_coherence(curdata, max.rt.diff = 10))
})

test_that("Same module, RT exceeds threshold -> incoherent", {
  curdata <- data.frame(
    Module_RTclust = c("73_3", "73_5"),
    time = c(100, 200),
    stringsAsFactors = FALSE
  )
  expect_false(check_compound_coherence(curdata, max.rt.diff = 10))
})

test_that("Module check runs before RT check (different modules, close RT)", {
  curdata <- data.frame(
    Module_RTclust = c("73_3", "85_3"),
    time = c(100, 101),
    stringsAsFactors = FALSE
  )
  expect_false(check_compound_coherence(curdata, max.rt.diff = 10))
})

test_that("Module extraction handles compound Module_RTclust formats", {
  # Module "120" with RT cluster "15"
  curdata <- data.frame(
    Module_RTclust = c("120_15", "120_3"),
    time = c(100, 105),
    stringsAsFactors = FALSE
  )
  expect_true(check_compound_coherence(curdata, max.rt.diff = 10))
})

test_that("RT exactly at threshold is coherent", {
  curdata <- data.frame(
    Module_RTclust = c("73_3", "73_5"),
    time = c(100, 110),
    stringsAsFactors = FALSE
  )
  expect_true(check_compound_coherence(curdata, max.rt.diff = 10))
})

test_that("RT just over threshold is incoherent", {
  curdata <- data.frame(
    Module_RTclust = c("73_3", "73_5"),
    time = c(100, 110.01),
    stringsAsFactors = FALSE
  )
  expect_false(check_compound_coherence(curdata, max.rt.diff = 10))
})

# Tests for enforce_compound_coherence()

test_that("enforce: single row returned unchanged", {
  curdata <- data.frame(
    Module_RTclust = "73_3",
    Adduct = "M+H",
    time = 100,
    stringsAsFactors = FALSE
  )
  result <- enforce_compound_coherence(curdata)
  expect_equal(nrow(result), 1)
  expect_equal(result$Adduct, "M+H")
})

test_that("enforce: single module returned unchanged", {
  curdata <- data.frame(
    Module_RTclust = c("73_3", "73_5", "73_3"),
    Adduct = c("M+H", "M+Na", "M+H_[+1]"),
    time = c(100, 102, 101),
    stringsAsFactors = FALSE
  )
  result <- enforce_compound_coherence(curdata)
  expect_equal(nrow(result), 3)
})

test_that("enforce: multi-module keeps largest group", {
  curdata <- data.frame(
    Module_RTclust = c("110_1", "110_2", "62_1"),
    Adduct = c("M+H", "M+H_[+1]", "M+H-H2O"),
    time = c(14.51, 14.32, 11.53),
    stringsAsFactors = FALSE
  )
  result <- enforce_compound_coherence(curdata)
  expect_equal(nrow(result), 2)
  # Should keep module 110 rows (2 rows) over module 62 (1 row)
  modules <- gsub(as.character(result$Module_RTclust), pattern = "_[0-9]*", replacement = "")
  expect_true(all(modules == "110"))
})

test_that("enforce: tie-breaking keeps first max module", {
  curdata <- data.frame(
    Module_RTclust = c("73_3", "85_2"),
    Adduct = c("M+H", "M+Na"),
    time = c(100, 200),
    stringsAsFactors = FALSE
  )
  result <- enforce_compound_coherence(curdata)
  # Both modules have 1 row; which.max returns first max
  expect_equal(nrow(result), 1)
})

test_that("enforce: preserves all columns", {
  curdata <- data.frame(
    Module_RTclust = c("110_1", "110_2", "62_1"),
    Adduct = c("M+H", "M+H_[+1]", "M+H-H2O"),
    time = c(14.51, 14.32, 11.53),
    score = c(100, 5000, 50),
    compound_id = rep("PEST0681", 3),
    stringsAsFactors = FALSE
  )
  result <- enforce_compound_coherence(curdata)
  expect_equal(ncol(result), ncol(curdata))
  expect_equal(names(result), names(curdata))
})
