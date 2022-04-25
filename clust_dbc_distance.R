packages <- c("mclust", "cluster")

## Now load or install&load all
package_check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

source("find_optimal_theta.R")
source("00_subsetlist.R")
source("01_algorithm_distance.R")
source("02_algorithm_clustering.R")
source("03_algorithm_find_outlier.R")

clust_dbc_distance <- function(distance_matrix, standard_label, minprop, alpha_c, optim_index) {
  theta_change <- find_optimal_theta(distance_matrix, standard_label, minprop, optim_index, alpha_c)
  ### optimize parameters
  optim_theta <- theta_change$optimal_theta
  neighbor_mat_alg1 <- findneighbours(distance_matrix, optim_theta)  # algorithm 1
  cluster_alg2 <- clustering_warping(neighbor_mat_alg1)     ## algorithm 2 
  cluster_alg3 <- include_isolatepoints(cluster_alg2[[2]], minprop, distance_matrix, alpha_c) ## algorithm 3
  n <- nrow(distance_matrix)
  ##############
  empirical_label <- rep(0, n)
  for (clustindex in 1:length(cluster_alg3$clust)) {
    empirical_label[cluster_alg3$clust[[clustindex]]$element] <- clustindex
  }
  empirical_label[empirical_label == 0] <- length(cluster_alg3$clust) + 1:length(which(empirical_label==0))
    
  
  value <- adjustedRandIndex(standard_label, empirical_label)
    names(value) <- "adjusted rand index"
  #}
  return (list(empirical_label, optim_theta, value, length(cluster_alg3$clust)))
}

