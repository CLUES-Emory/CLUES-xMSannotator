# Tests for the multi_rt_clusters helpers: select_cluster_anchors(),
# qualify_rt_clusters(), unqualify_rt_clusters().

make_rows <- function(compound_id, Adduct, time, mz, mean_int_vec,
                      isotope_ratio_dev = NA_real_, Module_RTclust = "73_3") {
  n <- length(Adduct)
  data.frame(
    compound_id = rep(compound_id, length.out = n),
    Adduct = Adduct, time = time, mz = mz, mean_int_vec = mean_int_vec,
    Module_RTclust = rep(Module_RTclust, length.out = n),
    isotope_ratio_dev = rep(isotope_ratio_dev, length.out = n),
    stringsAsFactors = FALSE
  )
}

# --- select_cluster_anchors ---

test_that("one anchor per retention-time cluster, not one per compound", {
  ann <- make_rows("A", rep("M-H", 3), c(20, 120, 220), c(100, 100, 100), c(5, 5, 5))
  out <- select_cluster_anchors(ann, max_diff_rt = 10)
  expect_equal(nrow(out), 3)
  expect_equal(sort(out$rt_cluster), 1:3)
  expect_equal(sort(out$time), c(20, 120, 220))
})

test_that("peaks within the window remain a single cluster", {
  ann <- make_rows("A", rep("M-H", 3), c(20, 25, 28), c(100, 100, 100), c(5, 9, 5))
  out <- select_cluster_anchors(ann, max_diff_rt = 10)
  expect_equal(nrow(out), 1)
})

test_that("within a cluster the isotope-supported candidate wins", {
  ann <- make_rows("X",
    Adduct = c("M-H", "M-H", "M-H_[+2]"),
    time = c(20, 120, 121), mz = c(100, 100, 102),
    mean_int_vec = c(999999, 50, 25),
    isotope_ratio_dev = c(NA, NA, 0.05))
  out <- select_cluster_anchors(ann, max_diff_rt = 10)
  expect_equal(nrow(out), 2)
  # the rt-120 cluster picks the monoisotopic row that has the isotope beside it
  expect_true(120 %in% out$time)
  # and the far more intense rt-20 peak is still reported, as its own cluster
  expect_true(20 %in% out$time)
})

test_that("every compound is represented at least once", {
  ann <- rbind(make_rows("A", "M-H", 10, 100, 5),
               make_rows("B", rep("M-H", 2), c(10, 300), c(200, 200), c(5, 5)))
  out <- select_cluster_anchors(ann, max_diff_rt = 10)
  expect_equal(sort(unique(out$compound_id)), c("A", "B"))
  expect_equal(nrow(out), 3)
})

test_that("result does not depend on input row order", {
  ann <- make_rows("A", c("M-H", "M-H", "M-H_[+2]", "M+Cl"),
                   c(20, 120, 121, 22), c(100, 100, 102, 135),
                   c(10, 50, 25, 7), c(NA, NA, 0.05, NA))
  ref <- select_cluster_anchors(ann, max_diff_rt = 10)
  set.seed(3)
  for (i in 1:20) {
    out <- select_cluster_anchors(ann[sample(nrow(ann)), , drop = FALSE], max_diff_rt = 10)
    expect_equal(out[order(out$time), c("time", "mz")],
                 ref[order(ref$time), c("time", "mz")], ignore_attr = TRUE)
  }
})

test_that("a compound with only isotope rows still yields an anchor", {
  ann <- make_rows("G", "M+ACN+H_[+1]", 100, 500, 5000, isotope_ratio_dev = 0.02)
  expect_equal(nrow(select_cluster_anchors(ann, max_diff_rt = 10)), 1)
})

test_that("empty input is handled", {
  ann <- make_rows("H", character(0), numeric(0), numeric(0), numeric(0))
  out <- select_cluster_anchors(ann, max_diff_rt = 10)
  expect_equal(nrow(out), 0)
  expect_true("rt_cluster" %in% names(out))
})

# --- qualify / unqualify round trip ---

test_that("qualify assigns one suffix per retention-time cluster", {
  ann <- data.frame(compound_id = rep("A", 4), time = c(10, 12, 200, 202),
                    stringsAsFactors = FALSE)
  q <- qualify_rt_clusters(ann, max_diff_rt = 10)
  expect_equal(q$compound_id, c("A__rt1", "A__rt1", "A__rt2", "A__rt2"))
})

test_that("qualify keeps compounds separate", {
  ann <- data.frame(compound_id = c("A", "B"), time = c(10, 10), stringsAsFactors = FALSE)
  q <- qualify_rt_clusters(ann, max_diff_rt = 10)
  expect_equal(q$compound_id, c("A__rt1", "B__rt1"))
})

test_that("unqualify restores the id and records the cluster", {
  ann <- data.frame(compound_id = rep("A", 4), time = c(10, 12, 200, 202),
                    stringsAsFactors = FALSE)
  u <- unqualify_rt_clusters(qualify_rt_clusters(ann, max_diff_rt = 10))
  expect_equal(u$compound_id, rep("A", 4))
  expect_equal(u$rt_cluster, c(1L, 1L, 2L, 2L))
})

test_that("unqualify is a no-op on unqualified input", {
  ann <- data.frame(compound_id = c("A", "B"), time = c(1, 2), stringsAsFactors = FALSE)
  expect_equal(unqualify_rt_clusters(ann), ann)
})

test_that("a compound id containing rt-like text is not mangled", {
  ann <- data.frame(compound_id = c("ABCrt12DEF", "ABCrt12DEF"), time = c(10, 200),
                    stringsAsFactors = FALSE)
  u <- unqualify_rt_clusters(qualify_rt_clusters(ann, max_diff_rt = 10))
  expect_equal(u$compound_id, rep("ABCrt12DEF", 2))
  expect_equal(u$rt_cluster, c(1L, 2L))
})
