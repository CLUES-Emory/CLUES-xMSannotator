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

test_that("Conf 0 with non-primary adduct stays 0", {
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

test_that("Conf 0 single M+H raised to 1", {
  ann <- make_annotation(
    "HMDB020",
    adducts = c("M+H"),
    scores = c(10),
    times = c(100),
    confidence = 0
  )

  result <- cap_confidence_with_evidence(ann, max.rt.diff = 10)
  expect_equal(unique(result$Confidence), 1)
})

test_that("Conf 0 single M-H raised to 1", {
  ann <- make_annotation(
    "HMDB021",
    adducts = c("M-H"),
    scores = c(10),
    times = c(100),
    confidence = 0
  )

  result <- cap_confidence_with_evidence(ann, max.rt.diff = 10)
  expect_equal(unique(result$Confidence), 1)
})

test_that("Orphan isotope row at Conf 2 capped to 0", {
  # Single isotope-suffixed adduct with NO base adduct row
  ann <- make_annotation(
    "PEST1802",
    adducts = c("M+ACN+H_[+1]"),
    scores = c(5000),
    times = c(100),
    confidence = 2
  )

  result <- cap_confidence_with_evidence(ann, max.rt.diff = 10)
  expect_equal(unique(result$Confidence), 0)
})

test_that("Incoherent multi-row is capped per RT group, not as one compound", {
  ann <- make_annotation(
    "HMDB009",
    adducts = c("M+H", "M+Na"),
    scores = c(50, 50),
    times = c(100, 200),  # 100-second spread >> 10s tolerance
    confidence = 2
  )

  # These two rows are 100 s apart: they are separate chromatographic entities,
  # so each is judged on the evidence co-located with it.
  #   M+H  at rt 100 -> Level 1 on its own primary adduct
  #   M+Na at rt 200 -> Level 0; alone, non-primary, no isotope
  #
  # This previously asserted that BOTH rows were capped to 1, i.e. the isolated
  # M+Na row inherited Level 1 from an M+H row 100 s away. That was the scope
  # bug: evidence was read from one coherent subset and written to every row of
  # the compound. See coherent_groups().
  result <- cap_confidence_with_evidence(ann, max.rt.diff = 10)
  expect_equal(result$Confidence[result$Adduct == "M+H"], 1)
  expect_equal(result$Confidence[result$Adduct == "M+Na"], 0)
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

# --- coherent_groups() and per-group capping ---
# Regression tests for the scope bug: evidence was read from the module surviving
# enforce_compound_coherence() but written back to EVERY row of the compound, so a
# single RT-outlier row could demote that compound's well-supported rows.

test_that("coherent_groups splits a module on RT beyond tolerance", {
  d <- make_annotation("X", c("M+Na", "M+ACN+H", "M+Na"), c(50, 50, 50),
                       times = c(45, 46, 90), module_rt_clust = "12_1")
  g <- coherent_groups(d, max.rt.diff = 10)
  expect_equal(length(unique(g)), 2)
  expect_equal(g[1], g[2])
  expect_false(g[1] == g[3])
})

test_that("coherent_groups keeps a coherent module in one group", {
  d <- make_annotation("X", c("M+Na", "M+ACN+H"), c(50, 50),
                       times = c(45, 46), module_rt_clust = "12_1")
  expect_equal(length(unique(coherent_groups(d, max.rt.diff = 10))), 1)
})

test_that("coherent_groups never lets a group span more than the tolerance", {
  # chained small gaps must not accumulate past max.rt.diff
  d <- make_annotation("X", rep("M+Na", 5), rep(50, 5),
                       times = c(0, 6, 12, 18, 24), module_rt_clust = "12_1")
  g <- coherent_groups(d, max.rt.diff = 10)
  spans <- tapply(d$time, g, function(x) max(x) - min(x))
  expect_true(all(spans <= 10))
})

test_that("coherent_groups separates different modules", {
  d <- make_annotation("X", c("M+Na", "M+Na"), c(50, 50), times = c(45, 46))
  d$Module_RTclust <- c("12_1", "88_2")
  expect_equal(length(unique(coherent_groups(d, max.rt.diff = 10))), 2)
})

test_that("an RT-outlier row does not demote the compound's coherent rows", {
  # Two M+Na/M+ACN+H rows at rt 45-46 justify Confidence 2 on their own.
  # A third row in the SAME module at rt 90 previously made the whole compound
  # incoherent; with no M+H present the Level-1 rescue could not apply and every
  # row fell to 0.
  ann <- make_annotation("TERP", c("M+Na", "M+ACN+H", "M+Na"), c(50, 50, 50),
                         times = c(45, 46, 90), confidence = 2,
                         module_rt_clust = "12_1")
  result <- cap_confidence_with_evidence(ann, max.rt.diff = 10)

  # the two co-eluting rows keep Confidence 2 on their own evidence
  expect_equal(result$Confidence[result$time %in% c(45, 46)], c(2, 2))
  # the lone outlier is judged on its own evidence: single non-primary adduct
  expect_equal(result$Confidence[result$time == 90], 0)
})

test_that("per-group capping does not rescue a group that lacks evidence", {
  ann <- make_annotation("Y", c("M+H", "M+H_[+1]", "M+Na"), c(50, 5000, 50),
                         times = c(45, 46, 90), confidence = 3,
                         module_rt_clust = "12_1")
  result <- cap_confidence_with_evidence(ann, max.rt.diff = 10)
  # isotope + base adduct co-eluting -> 3 retained
  expect_equal(result$Confidence[result$time %in% c(45, 46)], c(3, 3))
  # isolated M+Na, non-primary, no isotope -> 0
  expect_equal(result$Confidence[result$time == 90], 0)
})
