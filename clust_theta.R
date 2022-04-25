packages <- c("mclust", "fda", "cluster", "purrr",
             "funHDDC", "mfaces", "NbClust")

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
source("00_subsetlist.R")
source("01_algorithm_distance.R")
source("02_algorithm_clustering.R")
source("03_algorithm_find_outlier.R")
source("silhouette_index.R")
source("simulation_parameter.R")

clust_theta <- function(distance_matrix, standard_label, minprop, theta, alpha_c) {
  neighb_matrix <- findneighbours(distance_matrix, theta)
  n <- nrow(distance_matrix)
  cluster_result <- clustering_warping(neighb_matrix)
  final_cluster <- include_isolatepoints(cluster_result[[2]], minprop, distance_matrix, alpha_c)
  clusterlabel <- function(cluster) {
    clust <- final_cluster$clust
    isolate <- final_cluster$isolate
    empirical_label <- rep(0, nrow(distance_matrix))
    for (cl in 1:length(clust)) {
      empirical_label[clust[[cl]]$element] <- cl
    }
    if (length(isolate) > 0) {
      empirical_label[isolate] <- (length(clust) + 1):(length(clust) + length(isolate))
    }
    return (empirical_label)
  }
  
  empirical_label <- clusterlabel(final_cluster)
  #clust_label <- as.numeric(names(table(empirical_label)[which(table(empirical_label) > minprop * n)]))
  
  if (length(standard_label) > 0) {
  ari <- adjustedRandIndex(standard_label,
                           empirical_label)
  } else {
    ari <- 0
  }
   silh_value <- silhouette_index(empirical_label,
                            distance_matrix, minprop)
  if (length(silh_value) == 1 && is.na(silh_value) ) {
    return (c(ari, 0))
  } else {
    mean_silhouette <- mean(silh_value[which(silh_value[, 1] != 0), 3])
    return (c(ari, mean_silhouette))
  }
}

