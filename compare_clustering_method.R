packages = c("dbscan", "gdata")
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
source("clust_dbc_distance.R")
source("clust_hierarchical_distance.R")
source("clust_kmedoids_distance.R")
source("clust_funhddc.R")

compare_clustering_method <- function(distance_matrix, outlier_nb, run_time_dist, 
                                      pd, pc, minprop, alpha_c, optim_index) {
  
    standard_label <- rep(1:clusternumber, each = n / clusternumber)
    if (length(outlier_nb) > 0) {
      standard_label[outlier_nb] <- clusternumber + 1:length(outlier_nb)
    }

  ### Three distance based methods: dbc, hierarchical, kmedoids, and funhddc
  ################################### dbc
  start_time_dbc <- Sys.time()
  
  clust_procedure <- clust_dbc_distance(distance_matrix, standard_label, minprop, 
                                        alpha_c, optim_index)
  empiricallabel_dbc <- clust_procedure[[1]]
  end_time_dbc <- Sys.time()
  value_dbc <- clust_procedure[[3]]
  run_time_dbc <- difftime(end_time_dbc, start_time_dbc, units = "secs") + run_time_dist
  p_c_dbc <- length(intersect(which(empiricallabel_dbc > clust_procedure[[4]]), outlier_nb)) / length(outlier_nb)
  p_f_dbc <- length(setdiff(which(empiricallabel_dbc > clust_procedure[[4]]), outlier_nb)) / (n - length(outlier_nb))
  
  #################################### hierarchical clustering
  start_time_hierarchical <- Sys.time()
  
  optimalresult <- clust_hierarchical_distance(distance_matrix) 
  
  value_hierarchical <- adjustedRandIndex(standard_label, optimalresult$Best.partition)
  
  end_time_hierarchical <- Sys.time()
  run_time_hierarchical <- difftime(end_time_hierarchical, start_time_hierarchical, 
                                    units = "secs") + run_time_dist
  #######################################
  ####################### k-medoids
  start_time_k_medoids <- Sys.time()
  
  # function to compute coefficient
  
  pam_result <- clust_kmedoids_distance(distance_matrix)
  
 
    value_kmedoids <- adjustedRandIndex(standard_label, pam_result$clustering)
  
  end_time_k_medoids <- Sys.time()
  run_time_k_medoids <- difftime(end_time_k_medoids, start_time_k_medoids, units = "secs") + 
    run_time_dist
  ########################## DBSCAN
  start_time_dbscan <- Sys.time()
  optim_theta <- clust_procedure[[2]]
  eps <- quantile(unlist(upperTriangle(distance_matrix)), optim_theta)
  x <- matrix(rnorm(n * 5), nrow = n)
  y <- dist(x)
  lower_dis <- lowerTriangle(distance_matrix, diag = FALSE)
  for (dis_ind in 1:length(lower_dis)) {
    y[dis_ind] <- lower_dis[dis_ind]
  }
  dbscan_result <- dbscan(y, eps, minPts = nrow(distance_matrix) * minprop)
  cluster_result <- dbscan_result$cluster
  sign_out <- which(cluster_result == 0)
  clust_label <- max(cluster_result)
  cluster_result[sign_out] <- (clust_label + 1): (clust_label + length(sign_out))
  value_dbscan <- adjustedRandIndex(standard_label, cluster_result)
  end_time_dbscan <- Sys.time()
  
  run_time_dbscan <- end_time_dbscan - start_time_dbscan + run_time_dbc
  
  p_c_dbscan <- length(intersect(which(dbscan_result$cluster == 0), outlier_nb)) / length(outlier_nb)
  p_f_dbscan <- length(setdiff(which(dbscan_result$cluster == 0), outlier_nb)) / (n - length(outlier_nb))
  
  ############################ funHDDC method
  start_time_funhddc <- Sys.time()
  
  funhddc_result <- clust_funhddc(pd, pc)
  funhddc_clustlabel <- funhddc_result$class
  
  if (length(funhddc_clustlabel) == 0) {
    value_funhddc <- 0
  } else {
    value_funhddc <- adjustedRandIndex(standard_label, funhddc_clustlabel)
  }
  
  end_time_funhddc <- Sys.time()
  
  run_time_funhddc <- difftime(end_time_funhddc, start_time_funhddc, unit = "secs") 
  
  ################## 
  result <- c(value_dbc, 
              value_dbscan,
              value_hierarchical,
              value_kmedoids,
              value_funhddc,
              run_time_dbc, 
              run_time_dbscan,
              run_time_hierarchical,  
              run_time_k_medoids, 
              run_time_funhddc,
              p_c_dbc, p_f_dbc,
              p_c_dbscan, p_f_dbscan)
  
  names(result) <- c(paste("ARI", "_dbc", sep = ""),
                     paste("ARI", "_dbscan", sep = ""), 
                     paste("ARI", "_hierarchical", sep = ""), 
                     paste("ARI", "_kmedoids", sep=""), 
                     paste("ARI", "_funhddc", sep=""),
                     "run_time_dbc", 
                     "run_time_dbscan",
                     "run_time_hierarchical", 
                     "run_time_kmedoids", 
                     "run_time_funhddc", 
                     "pc_dbc", "pf_dbc", 
                     "pc_dbscan", "pf_dbscan")
  return (result)
}


