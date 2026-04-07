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


test_that("compute_pathway_scores assigns positive scores to well-supported pathways", {
  # 4 compounds in the same enriched pathway — exceeds the >= 3 threshold
  pathways <- data.frame(
    pathway = rep("P1", 4),
    compound = c("C001", "C002", "C003", "C004"),
    stringsAsFactors = FALSE
  )

  result <- CLUES.xMSannotator:::compute_pathway_scores(pathways)

  expect_true(all(result$pathway_score == 4))
  expect_equal(nrow(result), 4)
})


test_that("compute_pathway_scores filters out pathways with < 3 compounds", {
  pathways <- data.frame(
    pathway = rep("P_small", 2),
    compound = c("C001", "C002"),
    stringsAsFactors = FALSE
  )

  result <- CLUES.xMSannotator:::compute_pathway_scores(pathways)

  expect_equal(nrow(result), 0)
})


test_that("pathway-supported annotations receive a positive score increase", {
  # 6 significant annotations: 5 in pathway P1, 1 in background pathway
  annotations <- data.frame(
    compound = c(paste0("C0", sprintf("%02d", 1:5)), "C106"),
    adduct = rep("M+H", 6),
    peak = 1:6,
    module = rep(1L, 6),
    score = rep(1.0, 6),
    stringsAsFactors = FALSE
  )

  # P1 has 5 compounds (all significant); P_bg provides background
  # for Fisher's test to detect P1 enrichment
  pathways <- data.frame(
    pathway = c(rep("P1", 5), rep("P_bg", 15)),
    compound = c(
      paste0("C0", sprintf("%02d", 1:5)),
      paste0("C1", sprintf("%02d", 1:15))
    ),
    stringsAsFactors = FALSE
  )

  adduct_weights <- data.frame(
    adduct = "M+H", weight = 5, stringsAsFactors = FALSE
  )

  result <- compute_pathways(annotations, pathways, adduct_weights,
                             score_threshold = 0.1)

  p1_results <- result[result$compound %in%
    paste0("C0", sprintf("%02d", 1:5)), ]

  # P1 compounds should be boosted: original 1.0 + 5 pathway compounds
  expect_true(all(p1_results$score == 6.0))
  # C106 is in P_bg which has only 1 significant compound (< 3 threshold),
  # so it should keep its original score of 1.0
  expect_equal(result$score[result$compound == "C106"], 1.0)
})
