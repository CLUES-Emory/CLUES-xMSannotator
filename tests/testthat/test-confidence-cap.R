# Tests for cap_confidence_with_evidence() and add_confidence_labels()

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

# --- cap_confidence_with_evidence tests ---

test_that("Conf 3 with isotope rows + adduct stays 3", {
  ann <- make_annotation(
    "HMDB001",
    adducts = c("M+H", "M+H_[+1]"),
    scores = c(50, 5000),
    times = c(100, 101),
    confidence = 3
  )

  result <- cap_confidence_with_evidence(ann, max.rt.diff = 10)
  expect_equal(unique(result$Confidence), 3)
})

test_that("Conf 3 with 2 adducts, no isotopes -> capped to 2", {
  ann <- make_annotation(
    "HMDB002",
    adducts = c("M+H", "M+Na"),
    scores = c(50, 50),
    times = c(100, 102),
    confidence = 3
  )

  result <- cap_confidence_with_evidence(ann, max.rt.diff = 10)
  expect_equal(unique(result$Confidence), 2)
})

test_that("Conf 3 with single M+H, no isotopes -> capped to 1", {
  ann <- make_annotation(
    "HMDB003",
    adducts = c("M+H"),
    scores = c(50),
    times = c(100),
    confidence = 3
  )

  result <- cap_confidence_with_evidence(ann, max.rt.diff = 10)
  expect_equal(unique(result$Confidence), 1)
})

test_that("Conf 2 with 2 adducts, no isotopes stays 2", {
  ann <- make_annotation(
    "HMDB004",
    adducts = c("M+H", "M+Na"),
    scores = c(50, 50),
    times = c(100, 105),
    confidence = 2
  )

  result <- cap_confidence_with_evidence(ann, max.rt.diff = 10)
  expect_equal(unique(result$Confidence), 2)
})

test_that("Conf 2 with single M+H, no isotopes -> capped to 1", {
  ann <- make_annotation(
    "HMDB005",
    adducts = c("M+H"),
    scores = c(50),
    times = c(100),
    confidence = 2
  )

  result <- cap_confidence_with_evidence(ann, max.rt.diff = 10)
  expect_equal(unique(result$Confidence), 1)
})

test_that("Conf 2 with single M+Na (non-primary) -> capped to 0", {
  ann <- make_annotation(
    "HMDB006",
    adducts = c("M+Na"),
    scores = c(50),
    times = c(100),
    confidence = 2
  )

  result <- cap_confidence_with_evidence(ann, max.rt.diff = 10)
  expect_equal(unique(result$Confidence), 0)
})

test_that("Conf 4 (confirmed) is untouched", {
  ann <- make_annotation(
    "HMDB007",
    adducts = c("M+Na"),
    scores = c(50),
    times = c(100),
    confidence = 4
  )

  result <- cap_confidence_with_evidence(ann, max.rt.diff = 10)
  expect_equal(unique(result$Confidence), 4)
})

test_that("Conf 0 is untouched", {
  ann <- make_annotation(
    "HMDB008",
    adducts = c("M+Na"),
    scores = c(5),
    times = c(100),
    confidence = 0
  )

  result <- cap_confidence_with_evidence(ann, max.rt.diff = 10)
  expect_equal(unique(result$Confidence), 0)
})

test_that("Incoherent multi-row at Conf 2 -> capped to 0", {
  ann <- make_annotation(
    "HMDB009",
    adducts = c("M+H", "M+Na"),
    scores = c(50, 50),
    times = c(100, 200),  # 100-second spread >> 10s tolerance
    confidence = 2
  )

  result <- cap_confidence_with_evidence(ann, max.rt.diff = 10)
  expect_equal(unique(result$Confidence), 0)
})

test_that("Custom primary adducts: M+Na primary -> M+Na gets 1, M+H gets 0", {
  ann_na <- make_annotation(
    "HMDB010",
    adducts = c("M+Na"),
    scores = c(50),
    times = c(100),
    confidence = 2
  )
  ann_h <- make_annotation(
    "HMDB011",
    adducts = c("M+H"),
    scores = c(50),
    times = c(200),
    confidence = 2
  )
  ann <- rbind(ann_na, ann_h)

  result <- cap_confidence_with_evidence(
    ann,
    level1_primary_adducts = c("M+Na"),
    max.rt.diff = 10
  )

  conf_na <- unique(result$Confidence[result$compound_id == "HMDB010"])
  conf_h <- unique(result$Confidence[result$compound_id == "HMDB011"])

  expect_equal(conf_na, 1)  # M+Na is primary -> capped to 1
  expect_equal(conf_h, 0)   # M+H is NOT primary -> capped to 0
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

  result <- cap_confidence_with_evidence(ann, max.rt.diff = 10)
  expect_equal(nrow(result), 0)
})

# --- add_confidence_labels tests ---

test_that("Labels added correctly for all confidence levels", {
  ann <- rbind(
    make_annotation("C1", "M+H", 50, 100, confidence = 3),
    make_annotation("C2", "M+Na", 50, 200, confidence = 0),
    make_annotation("C3", "M+H", 50, 300, confidence = 4),
    make_annotation("C4", "M-H", 50, 400, confidence = 1),
    make_annotation("C5", "M+Na", 50, 500, confidence = 2)
  )

  result <- add_confidence_labels(ann)

  expect_true("Confidence_Level" %in% names(result))
  expect_equal(result$Confidence_Level[result$compound_id == "C1"], "High")
  expect_equal(result$Confidence_Level[result$compound_id == "C2"], "None")
  expect_equal(result$Confidence_Level[result$compound_id == "C3"], "Confirmed")
  expect_equal(result$Confidence_Level[result$compound_id == "C4"], "Low")
  expect_equal(result$Confidence_Level[result$compound_id == "C5"], "Medium")
})
