test_that("compute_pathways works on normal annotation tables", {
  annotations <- data.frame(
    compound = c("C001", "C001", "C002"),
    adduct = c("M+H", "M-H", "M+H"),
    peak = c(1L, 2L, 3L),
    module = c(1L, 1L, 2L),
    score = c(5.0, 3.0, 1.0),
    stringsAsFactors = FALSE
  )

  pathways <- data.frame(
    compound = c("C001", "C002"),
    pathway = c("P1", "P1"),
    stringsAsFactors = FALSE
  )

  adduct_weights <- data.frame(
    adduct = c("M+H", "M-H"),
    weight = c(5, 5),
    stringsAsFactors = FALSE
  )

  result <- compute_pathways(annotations, pathways, adduct_weights)

  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)
  expect_true("score" %in% names(result))
  expect_true(all(result$score >= 0.1))
})
