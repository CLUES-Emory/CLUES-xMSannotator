# enforce_compound_coherence_by_evidence() replaces the row-count rule used to
# decide which module survives into Stage4b.
#
# The old rule, enforce_compound_coherence(), takes names(which.max(table(mods))).
# When the modules are the same size -- 365 of 473 multi-module groups in CLU0009
# C18pos -- table() orders its names lexically and which.max() returns the first,
# so the survivor is chosen by string-sorting the module label. Labels are
# numbers held as text, so "12" beats "52" beats "7". In all 98 C18pos groups
# where the weaker module was kept, the discarded row was an M+H.

rowset <- function(...) {
  d <- do.call(rbind, list(...))
  rownames(d) <- NULL
  d
}

annrow <- function(Adduct, module, Confidence, mz = 100, time = 20,
                   mean_int_vec = 1000, compound_id = "CPD1") {
  data.frame(
    compound_id = compound_id, Adduct = Adduct, mz = mz, time = time,
    Module_RTclust = module, Confidence = Confidence,
    mean_int_vec = mean_int_vec, stringsAsFactors = FALSE
  )
}

kept <- function(d, ...) {
  out <- enforce_compound_coherence_by_evidence(d, ...)
  unique(gsub(as.character(out$Module_RTclust), pattern = "_[0-9]*", replacement = ""))
}

# --- pass-through cases ---

test_that("a single row is returned unchanged", {
  d <- annrow("M+H", "52_1", 1)
  expect_equal(enforce_compound_coherence_by_evidence(d), d)
})

test_that("a single module is returned unchanged", {
  d <- rowset(annrow("M+H", "52_1", 2), annrow("M+Na", "52_1", 2))
  expect_equal(nrow(enforce_compound_coherence_by_evidence(d)), 2)
})

test_that("rt sub-clusters of one module are not treated as different modules", {
  # 52_1 and 52_3 share the WGCNA module; the suffix is module-local
  d <- rowset(annrow("M+H", "52_1", 1), annrow("M+Na", "52_3", 1))
  expect_equal(nrow(enforce_compound_coherence_by_evidence(d)), 2)
})

# --- the ranking, criterion by criterion ---

test_that("higher Confidence beats a module with more rows", {
  d <- rowset(
    annrow("M+NH4", "12_1", 0), annrow("M+ACN+H", "12_1", 0), annrow("2M+H", "12_1", 0),
    annrow("M+H", "52_1", 1)
  )
  expect_equal(kept(d), "52")
})

test_that("on equal Confidence the module with more rows wins", {
  d <- rowset(
    annrow("M+NH4", "12_1", 1), annrow("M+ACN+H", "12_1", 1),
    annrow("M+Na", "52_1", 1)
  )
  expect_equal(kept(d), "12")
})

test_that("on equal Confidence and size the primary ion wins", {
  d <- rowset(annrow("M+H-H2O", "12_1", 3), annrow("M+H", "52_1", 3))
  expect_equal(kept(d), "52")
})

test_that("with all else equal the more intense module wins", {
  d <- rowset(
    annrow("M+NH4", "12_1", 0, mean_int_vec = 100),
    annrow("M+Na", "52_1", 0, mean_int_vec = 900)
  )
  expect_equal(kept(d), "52")
})

test_that("level1_primary_adducts is honoured for negative mode", {
  d <- rowset(annrow("M+Cl", "12_1", 2), annrow("M-H", "52_1", 2))
  expect_equal(kept(d, level1_primary_adducts = c("M+H", "M-H")), "52")
})

test_that("an isotope row counts toward its module's base adduct", {
  # M-H_[+1] must strip to M-H and satisfy the primary-ion criterion
  d <- rowset(annrow("M+Cl", "12_1", 3), annrow("M-H_[+1]", "52_1", 3))
  expect_equal(kept(d), "52")
})

# --- the regression: no dependence on the module label ---

test_that("the choice does not depend on the module label's sort order", {
  # same evidence, labels swapped: the M+H module must win either way
  a <- rowset(annrow("M+H-H2O", "12_3", 0), annrow("M+H", "52_1", 1))
  b <- rowset(annrow("M+H-H2O", "52_1", 0), annrow("M+H", "12_3", 1))
  expect_equal(kept(a), "52")
  expect_equal(kept(b), "12")
})

test_that("the old rule and the new one disagree exactly where evidence does", {
  d <- rowset(annrow("M+H-H2O", "12_3", 0), annrow("M+H", "52_1", 1))
  old <- enforce_compound_coherence(d)
  new <- enforce_compound_coherence_by_evidence(d)
  expect_equal(old$Adduct, "M+H-H2O")   # lexical: "12" sorts first
  expect_equal(new$Adduct, "M+H")       # evidence: Conf 1 beats Conf 0
})

test_that("the tartrate case keeps the M-H", {
  # C18neg tartrate, rt cluster 1: three single-row modules, one holding the
  # M-H at -0.3 ppm with Confidence 1. The row-count rule keeps module 12.
  d <- rowset(
    annrow("M-H2O-H", "12_3",  0, mz = 130.9986, time = 14.6, mean_int_vec = 202520),
    annrow("M+Hac-H", "43_1",  0, mz = 209.0303, time = 15.6, mean_int_vec = 429707),
    annrow("M-H",     "52_1",  1, mz = 149.0092, time = 15.3, mean_int_vec = 475285)
  )
  expect_equal(enforce_compound_coherence(d)$Adduct, "M-H2O-H")
  out <- enforce_compound_coherence_by_evidence(d)
  expect_equal(out$Adduct, "M-H")
  expect_equal(out$mz, 149.0092)
})

test_that("the oxoproline case is decided by the primary-ion criterion", {
  # C18pos: both modules Confidence 3 with two rows each, so only the third
  # criterion separates them. Lexically "58" sorts before "7".
  d <- rowset(
    annrow("M+NH4",      "58_1", 3, mz = 147.0764, time = 17.2, mean_int_vec = 6840592),
    annrow("M+NH4_[+1]", "58_1", 3, mz = 148.0798, time = 17.2, mean_int_vec = 288204),
    annrow("M+H",        "7_2",  3, mz = 130.0499, time = 19.3, mean_int_vec = 14865197),
    annrow("M+H_[+1]",   "7_2",  3, mz = 131.0532, time = 19.2, mean_int_vec = 642194)
  )
  expect_equal(unique(enforce_compound_coherence(d)$Module_RTclust), "58_1")
  expect_equal(kept(d), "7")
})

# --- robustness ---

test_that("a table without mean_int_vec still ranks on the first three criteria", {
  d <- rowset(annrow("M+NH4", "12_1", 0), annrow("M+H", "52_1", 1))
  d$mean_int_vec <- NULL
  expect_equal(kept(d), "52")
})

test_that("the result does not depend on input row order", {
  d <- rowset(
    annrow("M-H2O-H", "12_3", 0), annrow("M+Hac-H", "43_1", 0), annrow("M-H", "52_1", 1)
  )
  ref <- enforce_compound_coherence_by_evidence(d)$Adduct
  set.seed(7)
  for (i in 1:20) {
    shuffled <- d[sample(nrow(d)), , drop = FALSE]
    expect_equal(enforce_compound_coherence_by_evidence(shuffled)$Adduct, ref)
  }
})

test_that("every column and row of the winning module is preserved", {
  d <- rowset(
    annrow("M+NH4", "12_1", 0),
    annrow("M+H", "52_1", 1, mz = 200), annrow("M+Na", "52_1", 1, mz = 222)
  )
  out <- enforce_compound_coherence_by_evidence(d)
  expect_equal(nrow(out), 2)
  expect_equal(names(out), names(d))
  expect_equal(sort(out$mz), c(200, 222))
})
