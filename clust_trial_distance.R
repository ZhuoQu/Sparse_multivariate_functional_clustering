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

clust_trial_distance <- function(listresult, minprop, alpha_c, assess_methods, optim_index) {
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
  #################################### calculation of ETD distance: pointwise L^2
  start_time_l2 <- Sys.time()
  distance_matrix_l2 <- distance_mfd(pd, distance_measure = "L^2")
  end_time_l2 <- Sys.time()
  run_time_l2 <- difftime(end_time_l2, start_time_l2, units = "secs")
  
  result_l2 <- compare_clustering_method(distance_matrix_l2, outlier_nb, run_time_etd,
                                          pd, pc, minprop, alpha_c, optim_index)
  
  #################################### calculation of ETD distance: pointwise L^1
  start_time_l1 <- Sys.time()
  distance_matrix_l1 <- distance_mfd(pd, distance_measure = "L^1")
  end_time_l1 <- Sys.time()
  run_time_l1 <- difftime(end_time_l1, start_time_l1, units = "secs")
  
  result_l1 <- compare_clustering_method(distance_matrix_l1, outlier_nb, run_time_l1,
                                         pd, pc, minprop, alpha_c, optim_index)
  
  #################################### calculation of ETD distance: pointwise L^inf
  start_time_linf <- Sys.time()
  distance_matrix_linf <- distance_mfd(pd, distance_measure = "L^inf")
  end_time_linf <- Sys.time()
  run_time_linf <- difftime(end_time_linf, start_time_linf, units = "secs")
  
  result_linf <- compare_clustering_method(distance_matrix_linf, outlier_nb, 
                                                run_time_linf, pd, pc, minprop, 
                                                alpha_c, optim_index)
  return (list(etd_l2 = result_l2, etd_l1 = result_l1, etd_linf = result_linf))
}

# multi_result<-clust_trial(scenario, pc, "multivariate", "silhouette")
# var1_result<-clust_trial(scenario, pc, "univariate 1", "silhouette")
# var2_result<-clust_trial(scenario, pc, "univariate 2", "silhouette")
# var3_result<-clust_trial(scenario, pc, "univariate 3", "silhouette")
