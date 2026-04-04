#' Plot the clustering and density function.
#'
#' @param pdf Density estimate using stats::density(...) function. Black line.
#' @param data Raw data used to estimate kernel density.
#' @param peak_indices Computed indices where to plot the peaks.
#' @param plot_lines Whether to plot vertical lines at the cluster centers or not.
#' Red circles & lines.
#'
plot_clustering <- function(pdf, data, peak_indices, plot_lines = FALSE) {
  plot(pdf, type = "l", xaxt = "n")
  axis(1, at = data, las = 2)
  points(pdf$x[peak_indices], pdf$y[peak_indices], col = "red")

  if (plot_lines) abline(v = pdf$x[peak_indices], col = "red")
}


#' Compute the positions of peak indices for given kernel density estimate.
#'
#' @param density y attribute from kernel density estimates [stats::density()].
#'
#' @return Indices at which kernel density estimate position is a cluster.
#' @importFrom pastecs turnpoints
#' @importFrom rlist list.filter
compute_peak_indices <- function(density) {
  turning_points <- turnpoints(ts(density), calc.proba = FALSE)
  peak_indices <- turning_points$pos[turning_points$peaks]
  filtered <- list.filter(peak_indices, density[.] > .Machine$double.eps)
  return(filtered)
}


#' Compute number of equally spaced kernel points for density estimation.
#'
#' @param data Data for which to estimate the number of kernel points.
#'
#' @return Nearest power of 2 to length(data), but at least 2048.
compute_num_kernel_points <- function(data) {
  num_points <- length(data)
  nearest_power_of_two <- 2**ceiling(log2(num_points))
  max(nearest_power_of_two, 2048)
}


#' Compute kernel density estimate.
#'
#' @param data Data for which to estimate the kernel density.
#' @param width Kernel width to use. Default 1. [0.5;2]
#' @param kernel Kernel to use. See stats::density for details.
#' Recommended ["gaussian","biweight","cosine","epanechnikov"].
#'
#' @return Kernel density estimate [stats::density()].
estimate_kernel_density <- function(data, width = 1, kernel = "gaussian") {
  if (length(data) == 0) return(list(x = numeric(0), y = numeric(0)))

  start <- max(min(data) - 1, 0)
  end <- max(data) + 1

  pdf <- density(
    data,
    from = start,
    to = end,
    n = compute_num_kernel_points(data),
    width = width,
    kernel = kernel
  )
}

#' Compute cluster positions in data
#'
#' @param data Data for which to compute the clustering.
#' @param width Bandwidth of the clustering.
#' @param kernel Kernel to use. See [stats::density()] for details.
#' @param show If TRUE, plot the resulting clustering.
#' Scan rate should be consulted to optimize clustering. [0.5;2.0]
#'
#' @return Positions of dense clusters.
compute_cluster_positions <- function(data, width = 1, kernel = "gaussian", show = FALSE) {
  if (length(data) == 0) return(numeric(0))

  pdf <- estimate_kernel_density(data, width = width, kernel = kernel)
  peak_indices <- compute_peak_indices(pdf$y)
  if (show) plot_clustering(pdf, data, peak_indices)

  peak_positions <- pdf$x[peak_indices]
  return(peak_positions)
}

#' Compute the assignment of cluster ID based on nn-algorithm.
#'
#' @param clusters Position of cluster centers.
#' @param data Data to assign to clusters
#'
#' @return Indices of closest cluster for each query point computed using [RANN::nn2()].
#'
#' @importFrom RANN nn2
compute_cluster_assignments <- function(clusters, data) {
  if (length(clusters) == 0 || length(data) == 0) return(integer(0))

  query <- as.vector(data)
  model <- nn2(clusters, query, k = 1)
  return(model$nn.idx[, 1])
}

#' Compute RT modules based on intensity modules and chromatographic peak width.
#'
#' @param peak_table Feature table with columns ['rt', 'module','peak'].
#' @param peak_width Estimated chromatographic peak width.
#' @param show If TRUE, plot the resulting clustering.
#'
#' @return data frame with 4 columns: \emph{peak}, \emph{mean_intensity},
#' \emph{module} and \emph{rt_cluster}
#' @export
#' @import dplyr
compute_rt_modules <- function(peak_table, peak_width = 1, show = FALSE) {
  if (nrow(peak_table) == 0) {
    return(tibble(peak = integer(0), mean_intensity = numeric(0),
                  module = integer(0), rt_cluster = integer(0)))
  }

  rt_cluster <- {}

  for (subdata in peak_table %>% group_split(module)) {
    cluster_positions <- compute_cluster_positions(subdata$rt, width = peak_width, show = show)
    subdata$rt_cluster <- compute_cluster_assignments(cluster_positions, subdata$rt)
    rt_cluster <- bind_rows(rt_cluster, subdata %>% select(peak, rt_cluster))
  }
  peaks_without_sample_intensities <- peak_table %>%
    select(any_of(c("peak", "mean_intensity", "module")))

  output_table <- left_join(
    peaks_without_sample_intensities,
    rt_cluster,
    by = "peak"
  )
  return(output_table)
}
