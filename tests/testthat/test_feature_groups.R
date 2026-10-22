load_feature_group_case <- function() {
  peaks <- arrow::read_parquet("test-data/advanced_annotation/peak_table.parquet")
  peaks$feature_id <- sprintf("F%04d", peaks$peak)
  db <- arrow::read_parquet("test-data/advanced_annotation/database_small.parquet")
  # The test database has no compound_id, which get_chemscore() filters on
  db$compound_id <- db$hmdbid
  list(
    peaks = peaks,
    db = db,
    adducts = sample_adduct_table %>% filter(adduct %in% c("M-H", "M+H")),
    # 30 s retention-time bins stand in for externally computed groups
    groups = data.frame(feature_id = peaks$feature_id, group = paste0("g", floor(peaks$rt / 30)))
  )
}

run_with_groups <- function(case, feature_groups, feature_id_column = "feature_id", outloc = tempfile()) {
  advanced_annotation(
    case$peaks,
    case$db,
    adduct_table = case$adducts,
    peak_rt_width = 1,
    time_tolerance = 10,
    intensity_deviation_tolerance = 0.1,
    pathway_mode = "skip",
    feature_id_column = feature_id_column,
    feature_groups = feature_groups,
    outloc = outloc
  )
}

test_that("feature_groups needs feature_id_column", {
  case <- load_feature_group_case()
  case$peaks$feature_id <- NULL  # a non-numeric column fails as_peak_table() unless named as feature_id_column
  expect_error(run_with_groups(case, case$groups, feature_id_column = NULL), "needs feature_id_column")
})

patrick::with_parameters_test_that("feature_groups is validated:", {
  case <- load_feature_group_case()
  expect_error(run_with_groups(case, modify(case$groups)), message)
},
  patrick::cases(
    missing_columns = list(modify = function(g) dplyr::rename(g, grp = group),
                           message = "needs columns feature_id and group"),
    duplicated_ids = list(modify = function(g) rbind(g, g[1, ]), message = "duplicated feature_id"),
    na_group = list(modify = function(g) { g$group[1] <- NA; g }, message = "have no group \\(NA\\)"),
    peak_without_group = list(modify = function(g) g[-1, ], message = "1 peaks have no group")
  )
)

test_that("feature_groups become modules with one RT cluster", {
  case <- load_feature_group_case()
  outloc <- tempfile()
  run_with_groups(case, case$groups, outloc = outloc)

  clusters <- read.table(file.path(outloc, "Stage1_peak_clusters.txt"), sep = "\t", header = TRUE)
  expected_module <- as.integer(factor(case$groups$group))[match(clusters$feature_id, case$groups$feature_id)]

  expect_gt(nrow(clusters), 0)
  expect_equal(clusters$module, expected_module)
  expect_true(all(clusters$rt_cluster == 1))
  expect_equal(clusters$Module_RTclust, paste(expected_module, 1, sep = "_"))
})
