# Tests for count_rt_clusters()
#
# The count answers "how many candidate retention times did the anchor choose
# between", so clusters are separated on adjacent gaps in retention time alone.
# This is deliberately a different rule from coherent_groups(), which is module
# aware and measures each row against a group anchor so a chain of small gaps
# cannot accumulate past the tolerance. Here a chain of small gaps IS one
# candidate: the peaks are contiguous in time, so there is nothing to choose.

ann <- function(compound_id, time) {
  data.frame(compound_id = compound_id, time = time, stringsAsFactors = FALSE)
}

test_that("a single peak is one cluster", {
  out <- count_rt_clusters(ann("A", 100), max_diff_rt = 10)
  expect_equal(out$n_rt_clusters, 1L)
  expect_equal(out$compound_id, "A")
})

test_that("peaks within the window are one cluster", {
  out <- count_rt_clusters(ann(rep("A", 3), c(100, 105, 108)), max_diff_rt = 10)
  expect_equal(out$n_rt_clusters, 1L)
})

test_that("peaks beyond the window are separate clusters", {
  out <- count_rt_clusters(ann(rep("A", 2), c(100, 160)), max_diff_rt = 10)
  expect_equal(out$n_rt_clusters, 2L)
})

test_that("the real dichlorophenol case counts three", {
  # C18neg: rt 123.4, 183.8, 195.9 -- the three mass matches the anchor chose
  # between, and the reason its Confidence deserves less weight than a
  # single-candidate compound's.
  out <- count_rt_clusters(ann(rep("HFZWRUODUSTPEG", 3), c(123.4, 183.8, 195.9)),
                           max_diff_rt = 10)
  expect_equal(out$n_rt_clusters, 3L)
})

test_that("duplicate times do not inflate the count", {
  # several adducts and isotopes of one peak share a retention time
  out <- count_rt_clusters(ann(rep("A", 5), c(100, 100, 100, 160, 160)),
                           max_diff_rt = 10)
  expect_equal(out$n_rt_clusters, 2L)
})

test_that("counts are independent per compound", {
  d <- rbind(ann(rep("A", 2), c(10, 200)),
             ann(rep("B", 3), c(50, 52, 54)),
             ann("C", 300))
  out <- count_rt_clusters(d, max_diff_rt = 10)
  out <- out[order(out$compound_id), ]
  expect_equal(out$compound_id, c("A", "B", "C"))
  expect_equal(out$n_rt_clusters, c(2L, 1L, 1L))
})

test_that("the window argument is honoured", {
  d <- ann(rep("A", 2), c(100, 130))
  expect_equal(count_rt_clusters(d, max_diff_rt = 10)$n_rt_clusters, 2L)
  expect_equal(count_rt_clusters(d, max_diff_rt = 50)$n_rt_clusters, 1L)
})

test_that("a chain of small gaps is one cluster", {
  # contiguous peaks: no gap exceeds the window, so there is one candidate
  out <- count_rt_clusters(ann(rep("A", 5), c(0, 6, 12, 18, 24)), max_diff_rt = 10)
  expect_equal(out$n_rt_clusters, 1L)
})

test_that("empty annotation returns an empty frame with the right columns", {
  out <- count_rt_clusters(ann(character(0), numeric(0)), max_diff_rt = 10)
  expect_equal(nrow(out), 0)
  expect_equal(names(out), c("compound_id", "n_rt_clusters"))
})

test_that("every compound appears exactly once", {
  d <- rbind(ann(rep("A", 4), c(1, 2, 300, 301)), ann(rep("B", 2), c(5, 6)))
  out <- count_rt_clusters(d, max_diff_rt = 10)
  expect_equal(nrow(out), 2)
  expect_false(any(duplicated(out$compound_id)))
})
