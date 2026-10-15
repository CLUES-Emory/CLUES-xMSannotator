#' Compute correlation matrix for intensity table.
#' @import WGCNA
#' @param peak_intensity_matrix Intensity matrix.
#' @param correlation_method Correlation method, default "p" -> pearson
#' @return Peak correlation matrix.
#' @export
compute_peak_correlations <- function(peak_intensity_matrix, correlation_method = "p") {
  peak_correlation_matrix <- WGCNA::cor(
    peak_intensity_matrix,
    use = correlation_method,
    method = correlation_method
  )
  return(peak_correlation_matrix)
}

#' Performs clustering analysis on top of a peak intensity matrix
#'
#' @param peak_intensity_matrix matrix of peak intensities. Each column
#' represents intensities of one peak in different samples, each row represent
#' intensitiess of different peaks in one sample
#' @param peak_correlation_matrix correlation matrix of peak intensities
#' @param correlation_threshold see the \emph{cutHeight} parameter of
#' \link{WGCNA::mergeCloseModules}
#' @param deep_split see \link{WGCNA::mergeCloseModules}
#' @param min_cluster_size see \link{WGCNA::blockwiseModules},
#' \link{WGCNA::mergeCloseModules}, and \link{dynamicTreeCut::cutreeDynamic}
#' @param network_type see \link{WGCNA::adjacency}
#' @param seed Integer RNG seed pinning the clustering, or NULL to leave the
#'   caller's random stream untouched. Module labels feed the coherence checks
#'   that set confidence, so the value is fixed here rather than inherited.
#'
#' @return data frame with 3 columns: \emph{peak}, \emph{mean_intensity}, and
#' \emph{module}
#' @export
#' @import WGCNA
compute_peak_modules <- function(
  peak_intensity_matrix,
  peak_correlation_matrix,
  correlation_threshold,
  deep_split,
  min_cluster_size,
  network_type,
  seed = 12345
) {
  # Pin the RNG for the whole clustering step, restoring the caller's stream on
  # exit so this function does not perturb anything downstream of it.
  #
  # The primary path is already reproducible without this: WGCNA::blockwiseModules
  # saves .Random.seed, calls set.seed(randomSeed) with its own default of 54321,
  # and restores afterwards, so the caller's state cannot leak in. Verified by
  # running identical inputs from three separate R sessions under three different
  # caller seeds and getting byte-identical module vectors. What this adds is that
  # the seed becomes ours rather than a WGCNA default that could change upstream,
  # and the guarantee extends to the mergeCloseModules fallback below, which has
  # no randomSeed argument of its own.
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv())) {
      saved_seed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", saved_seed, envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }

  peak_correlation_matrix <- 1 - peak_correlation_matrix

  power_estimate <- WGCNA::pickSoftThreshold(peak_intensity_matrix)$powerEstimate
  if (is.na(power_estimate)) {
    warning("Unable to estimate soft threshold. Using fallback value.")
    power_estimate <- 6
  }

  # NOTE: When WGCNA is loaded incorectly then WGCNA::blockwiseModules may fail
  # due to conflicting functions WGCNA::cor and stats::cor.
  # if fails when debugging run cor <- WGCNA::cor before the call to the function
  modules <- WGCNA::blockwiseModules(
    datExpr = peak_intensity_matrix,
    checkMissingData = FALSE,
    blocks = dynamicTreeCut::cutreeDynamic(
      dendro = flashClust::flashClust(as.dist(peak_correlation_matrix)),
      distM = peak_correlation_matrix,
      minClusterSize = min_cluster_size,
      pamRespectsDendro = FALSE
    ),
    blockSizePenaltyPower = 100,
    power = power_estimate,
    deepSplit = deep_split,
    detectCutHeight = NULL,
    minModuleSize = min_cluster_size,
    pamRespectsDendro = FALSE,
    mergeCutHeight = 1 - correlation_threshold,
    randomSeed = seed
  )

  if (!modules$MEsOK) {
    warning("Unable to perform peak WGCNA clustering. Using fallback method.")

    distance_matrix <- WGCNA::TOMdist(
      WGCNA::adjacency(
        datExpr = peak_intensity_matrix,
        type = network_type,
        power = power_estimate,
        corOptions = list(use = "p", method = "pearson")
      )
    )

    modules <- WGCNA::mergeCloseModules(
      exprData = peak_intensity_matrix,
      colors = dynamicTreeCut::cutreeDynamic(
        dendro = flashClust::flashClust(as.dist(distance_matrix)),
        distM = distance_matrix,
        deepSplit = deep_split,
        minClusterSize = min_cluster_size,
        pamRespectsDendro = FALSE
      ),
      cutHeight = 1 - correlation_threshold,
      trapErrors = TRUE
    )
  }

  tibble(
    peak = as.integer(colnames(peak_intensity_matrix)),
    mean_intensity = colMeans(peak_intensity_matrix),
    module = as.integer(as.factor(modules$colors)),
  )
}
