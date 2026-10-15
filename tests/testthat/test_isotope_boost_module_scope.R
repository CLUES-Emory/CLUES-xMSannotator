# get_chemscore() multiplies a compound's score by 100 when isotope evidence is
# present. compute_best_score() scores each Module_RTclust separately and keeps
# the best, so the score belongs to one module -- but the isotope rows were
# gathered from the whole +/- max_diff_rt window, which can span several.
#
# Peaks that share an m/z but come from different molecules land in different
# WGCNA modules (within-module r ~ 0.98, between-module r ~ 0.10). Counting
# their isotopes as evidence therefore boosted a score on another molecule's
# data. Measured on C18pos: 137 module groups (141 rows) were boosted this way.
#
# These tests fix the boost to the module that actually won the scoring.

mono_rows <- function(extra = NULL) {
  # three adducts in module 2_1 -- this module wins the scoring
  d <- data.frame(
    compound_id = rep("CPD1", 3),
    Adduct = c("M+H", "M+Na", "M+NH4"),
    mz = c(200, 222, 217),
    time = rep(105, 3),
    Module_RTclust = rep("2_1", 3),
    mean_int_vec = c(5000, 4000, 3000),
    MonoisotopicMass = rep(199, 3),
    theoretical.mz = c(200, 222, 217),
    MatchCategory = rep("Unique", 3),
    Name = rep("test compound", 3),
    Formula = rep("C10H13NO", 3),
    MD = rep(0.1, 3),
    stringsAsFactors = FALSE
  )
  if (!is.null(extra)) d <- rbind(d, extra)
  d
}

isotope_row <- function(module, rt) {
  data.frame(
    compound_id = "CPD1", Adduct = "M+H_[+1]", mz = 201.0034, time = rt,
    Module_RTclust = module, mean_int_vec = 300, MonoisotopicMass = 199,
    theoretical.mz = 201.0034, MatchCategory = "Unique",
    Name = "test compound", Formula = "C10H13NO", MD = 0.1,
    stringsAsFactors = FALSE
  )
}

run_score <- function(ann) {
  ad <- c("M+H", "M+Na", "M+NH4")
  aw <- data.frame(Adduct = ad, Weight = rep(1, 3), stringsAsFactors = FALSE)
  at <- as.data.frame(get_default_adduct_table())
  at <- at[at$adduct %in% ad, ]
  ids <- paste0(ann$mz, "_", ann$time)
  gc <- matrix(0.99, nrow(ann), nrow(ann), dimnames = list(ids, ids))
  diag(gc) <- 1
  get_chemscore(
    compound_id = "CPD1", time = 105, mz = 200, annotation = ann,
    corthresh = 0.7, global_cor = gc, max_diff_rt = 10, adduct_weights = aw,
    filter.by = NULL, MplusH.abundance.ratio.check = FALSE, adduct_table = at
  )
}

score_of <- function(res) max(as.numeric(res$cur_chem_score))
iso_count <- function(res) sum(grepl("_[[]", res$Adduct))

test_that("the module that wins is the one scored", {
  res <- run_score(mono_rows())
  expect_equal(unique(res$Module_RTclust), "2_1")
  expect_equal(nrow(res), 3)
})

test_that("an isotope in a DIFFERENT module does not boost the score", {
  plain <- run_score(mono_rows())
  cross <- run_score(mono_rows(isotope_row("1_1", 100)))
  # the isotope sits in module 1_1, 5 s away and in a different correlation
  # module: it is a different molecule and must not count as evidence
  expect_equal(score_of(cross), score_of(plain))
  expect_equal(iso_count(cross), 0)
})

test_that("an isotope in the WINNING module still boosts the score 100x", {
  plain <- run_score(mono_rows())
  own <- run_score(mono_rows(isotope_row("2_1", 105)))
  expect_equal(score_of(own), score_of(plain) * 100)
  expect_equal(iso_count(own), 1)
})

test_that("only the winning module's isotope rows are carried into the output", {
  both <- run_score(mono_rows(rbind(isotope_row("1_1", 100),
                                    isotope_row("2_1", 105))))
  iso <- both[grepl("_[[]", both$Adduct), , drop = FALSE]
  expect_equal(nrow(iso), 1)
  expect_equal(iso$Module_RTclust, "2_1")
})

test_that("no isotope anywhere leaves the score unboosted", {
  plain <- run_score(mono_rows())
  # the unboosted score is whatever the adduct weights give; assert only that
  # adding no isotope changes nothing, and that it is not a multiple of 100
  expect_equal(score_of(plain) %% 100 == 0, FALSE)
  expect_equal(iso_count(plain), 0)
})
