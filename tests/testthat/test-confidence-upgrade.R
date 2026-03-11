# Tests for upgrade_confidence_with_evidence()

# Helper: create minimal annotation row(s) for a compound
make_annotation <- function(compound_id, adducts, scores, times, confidence = 0,
                            module_rt_clust = "73_3") {
  n <- length(adducts)
  data.frame(
    compound_id = rep(compound_id, n),
    Adduct = adducts,
    score = scores,
    time = times,
    Confidence = rep(confidence, n),
    mz = seq(100, by = 0.5, length.out = n),
    Formula = rep("C6H12O6", n),
    Name = rep("test_compound", n),
    Module_RTclust = rep(module_rt_clust, n),
    stringsAsFactors = FALSE
  )
}

# Default adduct_weights for tests
test_adduct_weights <- data.frame(
  adduct = c("M+H", "M-H", "M+Na", "M+NH4", "M+ACN+H"),
  weight = c(5, 5, 3, 3, 2),
  stringsAsFactors = FALSE
)

test_that("filter_by NULL: isotope rows + 2 base adducts + RT coherent -> Conf 3", {
  ann <- make_annotation(
    "HMDB001",
    adducts = c("M+Na", "M+NH4", "M+Na_[+1]"),
    scores = c(50, 50, 5000),
    times = c(100, 102, 101)
  )

  result <- upgrade_confidence_with_evidence(
    ann, filter.by = NULL, adduct_weights = test_adduct_weights, max.rt.diff = 10
  )

  expect_equal(unique(result$Confidence), 3)
})

test_that("filter_by set: isotope rows + 2 base adducts + RT coherent -> Conf 2", {
  ann <- make_annotation(
    "HMDB001",
    adducts = c("M+Na", "M+NH4", "M+Na_[+1]"),
    scores = c(50, 50, 5000),
    times = c(100, 102, 101)
  )

  result <- upgrade_confidence_with_evidence(
    ann, filter.by = c("M+H", "M-H"), adduct_weights = test_adduct_weights, max.rt.diff = 10
  )

  expect_equal(unique(result$Confidence), 2)
})

test_that("filter_by NULL: isotope rows + 1 base adduct + RT coherent -> Conf 2", {
  ann <- make_annotation(
    "HMDB002",
    adducts = c("M+Na", "M+Na_[+1]"),
    scores = c(50, 5000),
    times = c(100, 101)
  )

  result <- upgrade_confidence_with_evidence(
    ann, filter.by = NULL, adduct_weights = test_adduct_weights, max.rt.diff = 10
  )

  expect_equal(unique(result$Confidence), 2)
})

test_that("filter_by set: isotope rows + 1 base adduct + RT coherent -> Conf 1", {
  ann <- make_annotation(
    "HMDB002",
    adducts = c("M+Na", "M+Na_[+1]"),
    scores = c(50, 5000),
    times = c(100, 101)
  )

  result <- upgrade_confidence_with_evidence(
    ann, filter.by = c("M+H", "M-H"), adduct_weights = test_adduct_weights, max.rt.diff = 10
  )

  expect_equal(unique(result$Confidence), 1)
})

test_that("2 base adducts + RT coherent, no isotopes -> Conf 1 (both filter modes)", {
  ann <- make_annotation(
    "HMDB003",
    adducts = c("M+Na", "M+NH4"),
    scores = c(50, 50),
    times = c(100, 105)
  )

  # filter_by NULL
  result_null <- upgrade_confidence_with_evidence(
    ann, filter.by = NULL, adduct_weights = test_adduct_weights, max.rt.diff = 10
  )
  expect_equal(unique(result_null$Confidence), 1)

  # filter_by set
  result_set <- upgrade_confidence_with_evidence(
    ann, filter.by = c("M+H", "M-H"), adduct_weights = test_adduct_weights, max.rt.diff = 10
  )
  expect_equal(unique(result_set$Confidence), 1)
})

test_that("2 base adducts + high score + RT coherent -> Conf 1 (both filter modes)", {
  ann <- make_annotation(
    "HMDB004",
    adducts = c("M+Na", "M+NH4"),
    scores = c(150, 50),
    times = c(100, 105)
  )

  result_null <- upgrade_confidence_with_evidence(
    ann, filter.by = NULL, adduct_weights = test_adduct_weights, max.rt.diff = 10
  )
  expect_equal(unique(result_null$Confidence), 1)

  result_set <- upgrade_confidence_with_evidence(
    ann, filter.by = c("M+H", "M-H"), adduct_weights = test_adduct_weights, max.rt.diff = 10
  )
  expect_equal(unique(result_set$Confidence), 1)
})

test_that("Single adduct, no isotopes, low score -> stays Conf 0", {
  ann <- make_annotation(
    "HMDB005",
    adducts = c("M+Na"),
    scores = c(5),
    times = c(100)
  )

  result <- upgrade_confidence_with_evidence(
    ann, filter.by = c("M+H", "M-H"), adduct_weights = test_adduct_weights, max.rt.diff = 10
  )

  expect_equal(unique(result$Confidence), 0)
})

test_that("RT incoherent -> stays Conf 0 regardless of other evidence", {
  ann <- make_annotation(
    "HMDB006",
    adducts = c("M+Na", "M+NH4", "M+Na_[+1]"),
    scores = c(50, 50, 5000),
    times = c(100, 200, 150)  # 100-second spread >> 10s tolerance

)

  result <- upgrade_confidence_with_evidence(
    ann, filter.by = NULL, adduct_weights = test_adduct_weights, max.rt.diff = 10
  )

  expect_equal(unique(result$Confidence), 0)
})

test_that("Existing Conf 3 compound is untouched", {
  ann <- make_annotation(
    "HMDB007",
    adducts = c("M+H", "M+Na"),
    scores = c(200, 150),
    times = c(100, 102),
    confidence = 3
  )

  result <- upgrade_confidence_with_evidence(
    ann, filter.by = c("M+H", "M-H"), adduct_weights = test_adduct_weights, max.rt.diff = 10
  )

  expect_equal(unique(result$Confidence), 3)
})

test_that("Never downgrades existing confidence", {
  ann <- make_annotation(
    "HMDB008",
    adducts = c("M+Na"),
    scores = c(5),
    times = c(100),
    confidence = 2
  )

  result <- upgrade_confidence_with_evidence(
    ann, filter.by = c("M+H", "M-H"), adduct_weights = test_adduct_weights, max.rt.diff = 10
  )

  # Evidence would give Conf 0, but existing Conf 2 is preserved
  expect_equal(unique(result$Confidence), 2)
})

test_that("Empty annotation handled gracefully", {
  ann <- data.frame(
    compound_id = character(0),
    Adduct = character(0),
    score = numeric(0),
    time = numeric(0),
    Confidence = numeric(0),
    stringsAsFactors = FALSE
  )

  result <- upgrade_confidence_with_evidence(
    ann, filter.by = NULL, adduct_weights = test_adduct_weights, max.rt.diff = 10
  )

  expect_equal(nrow(result), 0)
})

test_that("Single adduct + high score + no isotopes -> stays Conf 0", {
  ann <- make_annotation(
    "HMDB009",
    adducts = c("M+Na"),
    scores = c(150),
    times = c(100)
  )

  result_null <- upgrade_confidence_with_evidence(
    ann, filter.by = NULL, adduct_weights = test_adduct_weights, max.rt.diff = 10
  )
  expect_equal(unique(result_null$Confidence), 0)

  result_set <- upgrade_confidence_with_evidence(
    ann, filter.by = c("M+H", "M-H"), adduct_weights = test_adduct_weights, max.rt.diff = 10
  )
  expect_equal(unique(result_set$Confidence), 0)
})

test_that("Multiple compounds upgraded independently", {
  ann1 <- make_annotation(
    "HMDB010",
    adducts = c("M+Na", "M+NH4", "M+Na_[+1]"),
    scores = c(50, 50, 5000),
    times = c(100, 102, 101)
  )
  ann2 <- make_annotation(
    "HMDB011",
    adducts = c("M+ACN+H"),
    scores = c(5),
    times = c(200)
  )
  ann <- rbind(ann1, ann2)

  result <- upgrade_confidence_with_evidence(
    ann, filter.by = NULL, adduct_weights = test_adduct_weights, max.rt.diff = 10
  )

  conf_10 <- unique(result$Confidence[result$compound_id == "HMDB010"])
  conf_11 <- unique(result$Confidence[result$compound_id == "HMDB011"])

  expect_equal(conf_10, 3)  # isotope + 2 adducts
  expect_equal(conf_11, 0)  # single adduct, low score
})

test_that("Module incoherence: pre-filter keeps best module, upgrades on coherent subset", {
  ann <- make_annotation(
    "HMDB012",
    adducts = c("M+Na", "M+NH4", "M+Na_[+1]"),
    scores = c(50, 50, 5000),
    times = c(100, 102, 101)
  )
  # Put rows in different modules: 2 in module 73, 1 in module 85
  ann$Module_RTclust <- c("73_3", "85_2", "73_3")

  result <- upgrade_confidence_with_evidence(
    ann, filter.by = NULL, adduct_weights = test_adduct_weights, max.rt.diff = 10
  )

  # Pre-filter keeps module 73 (M+Na + M+Na_[+1]): isotope rows + 1 base adduct -> Conf 2
  expect_equal(unique(result$Confidence), 2)
})

test_that("Orphan isotope row (no base row) stays Conf 0", {
  # Single isotope-suffixed adduct with NO base adduct row present
  ann <- make_annotation(
    "PEST1802",
    adducts = c("M+ACN+H_[+1]"),
    scores = c(5000),
    times = c(100)
  )

  result <- upgrade_confidence_with_evidence(
    ann, filter.by = NULL, adduct_weights = test_adduct_weights, max.rt.diff = 10
  )

  expect_equal(unique(result$Confidence), 0)
})

test_that("Same module with different RT clusters is coherent", {
  ann <- make_annotation(
    "HMDB013",
    adducts = c("M+Na", "M+NH4", "M+Na_[+1]"),
    scores = c(50, 50, 5000),
    times = c(100, 102, 101)
  )
  # Same module (73), different RT clusters (_3, _5) — still coherent
  ann$Module_RTclust <- c("73_3", "73_5", "73_3")

  result <- upgrade_confidence_with_evidence(
    ann, filter.by = NULL, adduct_weights = test_adduct_weights, max.rt.diff = 10
  )

  expect_equal(unique(result$Confidence), 3)
})
