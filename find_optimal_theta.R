source("clust_theta.R")

find_optimal_theta <- function(distance_matrix, standard_label, minprop, optim_index, alpha_c) {
  theta_candidate <- seq(0.01, 0.3, by = 0.01)
  clust_dataframe <- t(sapply(theta_candidate, function(theta) {
    clust_result <- clust_theta(distance_matrix, standard_label, minprop, theta, alpha_c)
  }))
  if (optim_index != "ARI") {
    value_index <- unlist(clust_dataframe[, 2])
  } else if (optim_index == "ARI") {
    value_index <- unlist(clust_dataframe[, 1])
  }
  names(value_index) <- theta_candidate
  optimal_theta <- as.numeric(names(value_index)[which(value_index == max(value_index))[1]])
  return (list(simulation = data.frame(theta = theta_candidate, value_index = value_index), 
               optimal_theta = optimal_theta) )
  }

# 
# find_optimal_theta <- function(distance_matrix, assess_index) {
#   theta_candidate <- seq(0.01, 0.25, by = 0.01)
#   clust_dataframe <- t(sapply(theta_candidate, function(k) {
#     clust_result <- clust_theta(distance_matrix, k)
#   }))
#   if (assess_index == "silhouette")
#   {
#     silhouette <- unlist(clust_dataframe[ , 2])
#     names(silhouette) <- theta_candidate
#     optimal_theta <- as.numeric(names(silhouette)[which(silhouette == max(silhouette))[1]])
#   } else if (assess_index == "ARI")
#   {
#     ari_value <- unlist(clust_dataframe[ , 1])
#     names(ari_value) <- theta_candidate
#     optimal_theta <- as.numeric(names(ari_value)[which(ari_value == max(ari_value))[1]])
#   }
# }
