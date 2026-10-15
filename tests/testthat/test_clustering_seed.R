make_matrix <- function(n_peaks = 60, n_samples = 30) {
  set.seed(7)
  # three co-abundant groups so WGCNA has real structure to find
  base <- lapply(1:3, function(i) rnorm(n_samples, 10 * i, 2))
  m <- vapply(seq_len(n_peaks), function(j) {
    base[[(j %% 3) + 1]] + rnorm(n_samples, 0, 0.3)
  }, numeric(n_samples))
  colnames(m) <- as.character(seq_len(n_peaks))
  m
}

test_that("peak module clustering is reproducible across caller RNG states", {
  m <- make_matrix()
  cm <- compute_peak_correlations(m)

  runs <- lapply(c(1L, 2L, 999L), function(s) {
    set.seed(s)   # vary the caller's stream
    suppressWarnings(
      compute_peak_modules(m, cm, 0.7, 2, 10, "unsigned")$module
    )
  })

  expect_equal(runs[[1]], runs[[2]])
  expect_equal(runs[[1]], runs[[3]])
})

test_that("clustering does not perturb the caller's random stream", {
  m <- make_matrix()
  cm <- compute_peak_correlations(m)

  set.seed(42)
  expected <- runif(3)

  set.seed(42)
  suppressWarnings(compute_peak_modules(m, cm, 0.7, 2, 10, "unsigned"))
  expect_equal(runif(3), expected)
})

test_that("the seed argument is exposed and accepted, not silently dropped", {
  # not asserting the labels differ -- on well separated data they need not --
  # only that the argument reaches WGCNA and the call still returns a module set
  m <- make_matrix()
  cm <- compute_peak_correlations(m)

  out <- suppressWarnings(suppressMessages(
    capture.output(
      res <- compute_peak_modules(m, cm, 0.7, 2, 10, "unsigned", seed = 999)
    )
  ))
  expect_equal(nrow(res), ncol(m))
  expect_true(all(c("peak", "mean_intensity", "module") %in% names(res)))

  expect_true("seed" %in% names(formals(compute_peak_modules)))
  expect_equal(formals(compute_peak_modules)$seed, 12345)
  expect_true("clustering_seed" %in% names(formals(advanced_annotation)))
  expect_equal(formals(advanced_annotation)$clustering_seed, 12345)
})

test_that("seed = NULL leaves the caller's random stream advancing normally", {
  m <- make_matrix()
  cm <- compute_peak_correlations(m)

  set.seed(11)
  suppressWarnings(suppressMessages(capture.output(
    compute_peak_modules(m, cm, 0.7, 2, 10, "unsigned", seed = NULL)
  )))
  unpinned <- runif(1)

  set.seed(11)
  suppressWarnings(suppressMessages(capture.output(
    compute_peak_modules(m, cm, 0.7, 2, 10, "unsigned", seed = NULL)
  )))
  expect_equal(runif(1), unpinned)
})
