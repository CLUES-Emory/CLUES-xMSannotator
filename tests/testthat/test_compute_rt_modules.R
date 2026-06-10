test_that("RT_clust: Integration test: RT based clustering works.", {
    skip_on_ci()
    
    peaks <- readRDS("test-data/rt_modules/qc_solvent.rda")
    peaks <- dplyr::rename(peaks, peak = feature)
    peaks <- dplyr::rename_with(
        peaks,
        ~ paste0("intensity_", .x),
        starts_with("Tribrid")
    )
    peaks$peak <- as.integer(peaks$peak)

    intensity_matrix <- t(select(peaks, starts_with("intensity")))
    intensity_matrix <- magrittr::set_colnames(intensity_matrix, peaks$peak)
    correlation_matrix <- WGCNA::cor(intensity_matrix, use = "p", method = "p")

    correlation_threshold <- 0.7
    deep_split <- 2
    min_cluster_size <- 10
    network_type <- "unsigned"

    peak_modules <- compute_peak_modules(
        peak_intensity_matrix = intensity_matrix,
        peak_correlation_matrix = correlation_matrix,
        correlation_threshold = correlation_threshold,
        deep_split = deep_split,
        min_cluster_size = min_cluster_size,
        network_type = network_type
    )

    peaks <- inner_join(peaks, peak_modules, by = "peak")
    rt_modules <- compute_rt_modules(peaks)
    expect_equal(nrow(peaks), nrow(rt_modules))
})

test_that("RT_clust: Cluster assignments are comparable.", {
    peaks <- readRDS("test-data/rt_modules/xmsannotator_qc_matrix_stage1.rda")
    peaks$peak <- as.integer(rownames(peaks))
    peaks <- dplyr::rename(
        peaks,
        rt = time,
        module = Module,
        RTclust_old = RTclust
    )

    actual <- compute_rt_modules(peaks, peak_width = 1)
    y2d <- gplots::hist2d(
        peaks$RTclust_old,
        actual$rt_cluster,
        show = FALSE
    )

    mutual_information <- entropy::mi.empirical(y2d$counts) # maximum similarity is 5.1
    expect_true(mutual_information > 2.5)
})

test_that("RT_clust: Cluster assignments are comparable.", {
    peaks <- readRDS("test-data/rt_modules/xmsannotator_qc_matrix_stage1.rda")
    peaks$peak <- as.integer(rownames(peaks))
    peaks <- dplyr::rename(
        peaks,
        rt = time,
        module = Module,
        RTclust_old = RTclust
    )

    actual <- compute_rt_modules(peaks, peak_width = 1) %>% select(rt_cluster)

    expected <- readRDS("test-data/rt_modules/expected.rds")
    expect_equal(actual, expected)
})

patrick::with_parameters_test_that("compute_cluster_assignments can handle various 1d data:",{
    actual <- compute_cluster_assignments(clusters, data)
    expect_equal(actual, expected)
}, patrick::cases(
    list_1d = list(
        clusters = c(1, 4),
        data = list(1.2, 4.2),
        expected = c(1, 2)
    ),
    array_1d = list(
        clusters = c(1, 4),
        data = array(c(1.2, 4.2)),
        expected = c(1, 2)
    ),
    vector_1d = list(
        clusters = c(1, 4),
        data = as.vector(c(1.2, 0.9)),
        expected = c(1, 1)
    )
))

test_that("estimate_kernel_density returns empty on empty input", {
  result <- CLUES.xMSannotator:::estimate_kernel_density(numeric(0))
  expect_equal(length(result$x), 0)
  expect_equal(length(result$y), 0)
})

test_that("compute_cluster_positions returns empty on empty input", {
  result <- CLUES.xMSannotator:::compute_cluster_positions(numeric(0))
  expect_equal(length(result), 0)
  expect_is(result, "numeric")
})

test_that("compute_cluster_assignments returns empty on empty input", {
  result <- compute_cluster_assignments(numeric(0), numeric(0))
  expect_equal(length(result), 0)
  expect_is(result, "integer")
})

test_that("compute_rt_modules returns correct empty table on empty input", {
  empty_table <- data.frame(
    peak = integer(0), rt = numeric(0),
    module = integer(0), mean_intensity = numeric(0)
  )
  result <- compute_rt_modules(empty_table)
  expect_equal(nrow(result), 0)
  expect_true(all(c("peak", "mean_intensity", "module", "rt_cluster") %in% names(result)))
})