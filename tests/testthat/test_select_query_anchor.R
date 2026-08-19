# Tests for select_query_anchor()
#
# The property being bought is that the chosen anchor depends on isotope
# evidence, not on where a row happens to sit in the table.

make_rows <- function(compound_id, Adduct, time, mz, mean_int_vec,
                      isotope_ratio_dev = NA_real_, Module_RTclust = "73_3") {
  n <- length(Adduct)
  data.frame(
    compound_id = rep(compound_id, length.out = n),
    Adduct = Adduct,
    time = time,
    mz = mz,
    mean_int_vec = mean_int_vec,
    Module_RTclust = rep(Module_RTclust, length.out = n),
    isotope_ratio_dev = rep(isotope_ratio_dev, length.out = n),
    stringsAsFactors = FALSE
  )
}

test_that("returns exactly one row per compound", {
  ann <- rbind(
    make_rows("A", c("M-H", "M-H"), c(10, 50), c(100, 100), c(5, 5)),
    make_rows("B", c("M+H"), 20, 200, 7),
    make_rows("C", c("M+H", "M+Na", "M+Cl"), c(30, 31, 32), c(300, 301, 302), c(1, 2, 3))
  )
  out <- select_query_anchor(ann, max_diff_rt = 10)
  expect_equal(nrow(out), 3)
  expect_equal(sort(out$compound_id), c("A", "B", "C"))
})

test_that("prefers the candidate with co-eluting isotope support", {
  # two M-H peaks; only the second has an isotope beside it
  ann <- make_rows("X",
    Adduct = c("M-H", "M-H", "M-H_[+2]"),
    time = c(35.1, 47.8, 47.9),
    mz = c(160.9566, 160.9566, 162.9537),
    mean_int_vec = c(999999, 100, 50),
    isotope_ratio_dev = c(NA, NA, 0.06)
  )
  out <- select_query_anchor(ann, max_diff_rt = 10)
  expect_equal(nrow(out), 1)
  # chosen despite being ~10000x less intense
  expect_equal(out$time, 47.8)
})

test_that("with equal support, better ratio fidelity wins (the dichlorophenol case)", {
  # both retention times carry one 37Cl isotope; rt 47.8 matches theory better.
  # Values are the observed HILICneg deviations for 2,4-dichlorophenol.
  ann <- make_rows("HFZWRUODUSTPEG",
    Adduct = c("M-H_[+2]", "M-H", "M-H", "M-H_[+2]"),
    time = c(35.04, 35.13, 47.81, 47.88),
    mz = c(162.9536, 160.9566, 160.9566, 162.9537),
    mean_int_vec = c(99889, 180435, 59771, 36002),
    isotope_ratio_dev = c(0.13488065, NA, NA, 0.05874152)
  )
  out <- select_query_anchor(ann, max_diff_rt = 10)
  expect_equal(nrow(out), 1)
  expect_equal(out$time, 47.81)
  expect_false(out$time == 35.13)  # what row order currently gives
})

test_that("intensity does not override ratio fidelity", {
  # the worse-fidelity peak is also the more abundant one
  ann <- make_rows("Y",
    Adduct = c("M-H", "M-H_[+2]", "M-H", "M-H_[+2]"),
    time = c(10, 10.1, 60, 60.1),
    mz = c(100, 102, 100, 102),
    mean_int_vec = c(1e6, 5e5, 10, 5),
    isotope_ratio_dev = c(NA, 0.40, NA, 0.01)
  )
  out <- select_query_anchor(ann, max_diff_rt = 10)
  expect_equal(out$time, 60)
})

test_that("falls back to primary adduct, then to intensity", {
  # no isotopes anywhere: M-H beats M+Cl
  ann <- make_rows("P", c("M+Cl", "M-H"), c(10, 60), c(100, 101), c(900, 5))
  out <- select_query_anchor(ann, level1_primary_adducts = c("M+H", "M-H"),
                             max_diff_rt = 10)
  expect_equal(out$Adduct, "M-H")

  # neither is primary: the more abundant wins
  ann2 <- make_rows("Q", c("M+Cl", "M+Na"), c(10, 60), c(100, 101), c(5, 900))
  out2 <- select_query_anchor(ann2, level1_primary_adducts = c("M+H", "M-H"),
                              max_diff_rt = 10)
  expect_equal(out2$Adduct, "M+Na")
})

test_that("level1_primary_adducts is honoured rather than hardcoded", {
  ann <- make_rows("R", c("M+Na", "M-H"), c(10, 60), c(100, 101), c(5, 5))
  out <- select_query_anchor(ann, level1_primary_adducts = c("M+Na"),
                             max_diff_rt = 10)
  expect_equal(out$Adduct, "M+Na")
})

test_that("isotope support is not counted across a module boundary", {
  ann <- make_rows("S",
    Adduct = c("M-H", "M-H", "M-H_[+2]"),
    time = c(10, 60, 60.1),
    mz = c(100, 100, 102),
    mean_int_vec = c(5, 5, 5),
    isotope_ratio_dev = c(NA, NA, 0.01)
  )
  ann$Module_RTclust <- c("12_1", "88_1", "12_1")  # isotope is in a DIFFERENT module
  out <- select_query_anchor(ann, max_diff_rt = 10)
  # the rt-60 row cannot claim the rt-60.1 isotope: different module. Neither
  # candidate has support, so the tie falls through to mz then time.
  expect_equal(out$time, 10)
})

test_that("isotope support is not counted beyond max_diff_rt", {
  ann <- make_rows("T",
    Adduct = c("M-H", "M-H", "M-H_[+2]"),
    time = c(10, 60, 90),
    mz = c(100, 100, 102),
    mean_int_vec = c(5, 5, 5),
    isotope_ratio_dev = c(NA, NA, 0.01)
  )
  out <- select_query_anchor(ann, max_diff_rt = 10)
  expect_equal(out$time, 10)  # rt-60 gains nothing from an isotope 30 s away
})

test_that("a non-NA dev on a monoisotopic row is ignored", {
  # reformat_annotation_table() leaves the adduct unsuffixed when
  # mass_number_difference is 0, so such rows can carry a value that does not
  # describe them. 41 rows in HILICneg do.
  ann <- make_rows("U",
    Adduct = c("M-H", "M-H2O-H"),
    time = c(10, 60),
    mz = c(100, 101),
    mean_int_vec = c(5, 5),
    isotope_ratio_dev = c(NA, 0.001)  # excellent-looking, but not an isotope row
  )
  out <- select_query_anchor(ann, level1_primary_adducts = c("M-H"),
                             max_diff_rt = 10)
  expect_equal(out$Adduct, "M-H")  # decided by primary adduct, not by the dev
})

test_that("result is invariant to input row order", {
  ann <- rbind(
    make_rows("HFZWRUODUSTPEG",
      Adduct = c("M-H_[+2]", "M-H", "M-H", "M-H_[+2]"),
      time = c(35.04, 35.13, 47.81, 47.88),
      mz = c(162.9536, 160.9566, 160.9566, 162.9537),
      mean_int_vec = c(99889, 180435, 59771, 36002),
      isotope_ratio_dev = c(0.13488065, NA, NA, 0.05874152)),
    make_rows("Z", c("M+H", "M+Na", "M+H_[+1]"), c(5, 6, 5.1),
              c(200, 222, 201), c(10, 20, 3), c(NA, NA, 0.02)),
    make_rows("W", c("M+Cl", "M+Na"), c(80, 90), c(300, 301), c(4, 4))
  )
  reference <- select_query_anchor(ann, max_diff_rt = 10)

  set.seed(1)
  for (i in 1:25) {
    shuffled <- ann[sample(nrow(ann)), , drop = FALSE]
    out <- select_query_anchor(shuffled, max_diff_rt = 10)
    expect_equal(out, reference, ignore_attr = TRUE)
  }
})

test_that("ties are broken deterministically by mz then time", {
  # two rows identical in every ranked criterion except mz
  ann <- make_rows("D", c("M+Na", "M+Na"), c(10, 10), c(150, 120), c(5, 5))
  out <- select_query_anchor(ann, max_diff_rt = 10)
  expect_equal(out$mz, 120)
})

test_that("a table without isotope_ratio_dev still works", {
  ann <- make_rows("E", c("M+Cl", "M-H"), c(10, 60), c(100, 101), c(900, 5))
  ann$isotope_ratio_dev <- NULL
  out <- select_query_anchor(ann, max_diff_rt = 10)
  expect_equal(nrow(out), 1)
  expect_equal(out$Adduct, "M-H")  # primary adduct still decides
  expect_false("isotope_ratio_dev" %in% names(out))
})

test_that("isotope support still ranks when dev is absent", {
  ann <- make_rows("F",
    Adduct = c("M+Na", "M+Na", "M+Na_[+1]"),
    time = c(10, 60, 60.1),
    mz = c(100, 100, 101),
    mean_int_vec = c(900, 5, 5)
  )
  ann$isotope_ratio_dev <- NULL
  out <- select_query_anchor(ann, max_diff_rt = 10)
  expect_equal(out$time, 60)  # support beats intensity
})

test_that("compound with only isotope rows still returns a row", {
  ann <- make_rows("G", c("M+ACN+H_[+1]"), 100, 500, 5000,
                   isotope_ratio_dev = 0.02)
  out <- select_query_anchor(ann, max_diff_rt = 10)
  expect_equal(nrow(out), 1)
})

test_that("empty annotation returns empty", {
  ann <- make_rows("H", character(0), numeric(0), numeric(0), numeric(0))
  out <- select_query_anchor(ann, max_diff_rt = 10)
  expect_equal(nrow(out), 0)
})

test_that("all columns are preserved for the caller", {
  ann <- make_rows("I", c("M+H", "M+Na"), c(10, 12), c(100, 122), c(5, 6))
  ann$Name <- "test"; ann$Formula <- "C6H12O6"
  out <- select_query_anchor(ann, max_diff_rt = 10)
  expect_equal(sort(names(out)), sort(names(ann)))
})
