#'@title Efficiency of Cluster Sampling for Crop Surveys
#' @param x Datasets
#' @param N Number of clusters
#' @import stats  dplyr
#' @return
#' \itemize{
#'   \item results: Results
#' }
#' @export
#'
#' @examples
#'N_clusters <- 105
#'orchards_per_cluster <- 4
#'data <- matrix(rnorm(N_clusters * orchards_per_cluster), nrow = orchards_per_cluster, byrow = TRUE)
#'colnames(data) <- paste0("Cluster_", 1:N_clusters)
#'demo_data <- data.frame(data)
#'result_imcluster <- ImCluster(demo_data, N_clusters)
#'
#' @references
#' \itemize{
#'\item Iqbal, J. M., Faizan, D and  Mansha,  G. (2018) . A Review on the Recent Development on the Cluster Sampling. Biostatistics and Biometrics.  5(5): 555673. DOI: 10.19080/BBOAJ.2018.05.555673
#'}
ImCluster <- function(x, N = NULL) {
  x = data.frame(x)
  n = nrow(x)
  m = ncol(x)

  # Calculate cluster mean
  yibar = apply(x, 1, mean)
  ynbar = sum(yibar) / n

  # Calculate cluster variance
  yibar2 = sum(yibar^2)
  nybar2 = ynbar^2
  if (!is.null(N)) {
    vynbar = ((1/n) - (1/N)) * (1/(n-1)) * (yibar2 - n * nybar2)
    se = sqrt(vynbar)
  } else {
    se = NULL
  }

  # Calculate mean square between clusters and mean square within clusters
  msb2 = (1/(n-1)) * (yibar2 - n * nybar2)
  si2 = apply(x, 1, var)
  sw2 = sum(si2/n)

  # Calculate total variance
  s2 = ((n-1) * m * msb2) + (n * (m-1) * sw2)
  S2 = s2 / ((n * m) - 1)

  # Calculate relative efficiency and intra-class correlation coefficient
  re = S2 / (m * msb2)
  rho = (1 - re) / ((m-1) * re)

  # Calculate Design effect
  DEF<-1+rho*(m-1)

  result <- list(
    clusterMean = ynbar,
    ClusterVariance = ifelse(!is.null(N), vynbar, NULL),
    se = se,
    DFBetweenClusters = n - 1,
    DFWithinClusters = n * (m - 1),
    DFTotal = n * m - 1,
    betweenCluSSq = msb2,
    WithinCluSSq = sw2,
    Total = S2,
    RelativeEfficiency = re,
    EstimatorRho = rho,
    DesignEffect= DEF
  )

  return(result)
}


