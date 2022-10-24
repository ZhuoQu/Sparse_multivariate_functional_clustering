packages = c("mclust", "fda", "cluster", "purrr",
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
source("sparsify_scenario.R")
source("01_algorithm_distance.R")
source("find_optimal_theta.R")
source("clust_dbc_distance.R")
source("clust_hierarchical_distance.R")
source("clust_kmedoids_distance.R")
source("clust_funhddc.R")
source("compare_clustering_method.R")

clust_trial <- function(listresult, minprop, alpha_c, assess_methods, optim_index) {
  n <- nrow(scenario[[1]][[1]])
  
  if (assess_methods == "multivariate") {
    pd <- listresult[[1]]
  } else if (assess_methods == "univariate 1") {
    pd <- listresult[[2]]
  } else if (assess_methods == "univariate 2") {
    pd <- listresult[[3]]
  } else if (assess_methods == "univariate 3") {
    pd <- listresult[[4]]
  } 
  
  outlier_nb <- scenario[[2]]
  #################################### calculation of ETD distance
  start_time_etd <- Sys.time()
  distance_matrix_etd <- distance_mfd(pd, distance_measure = "ETD")
  end_time_etd <- Sys.time()
  run_time_etd <- difftime(end_time_etd, start_time_etd, units = "secs")
  
  result_etd <- compare_clustering_method(distance_matrix_etd, outlier_nb, run_time_etd,
                                      pd, pc, minprop, alpha_c, optim_index)
  
  return (result_etd)
}

# multi_result<-clust_trial(scenario, pc, "multivariate", "silhouette")
# var1_result<-clust_trial(scenario, pc, "univariate 1", "silhouette")
# var2_result<-clust_trial(scenario, pc, "univariate 2", "silhouette")
# var3_result<-clust_trial(scenario, pc, "univariate 3", "silhouette")
