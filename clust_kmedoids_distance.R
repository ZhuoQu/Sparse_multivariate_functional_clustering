packages = c("NbClust", "cluster", "purrr")
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

clust_kmedoids_distance <- function(distance_matrix) {
  ac <- function(k) {
    result<-pam(distance_matrix, k, diss=TRUE, stand = FALSE)
    return(result$silinfo$avg.width)
  }
  k_candidates <- 2:(sqrt(n / 2) + 5)
  silhouette <- map_dbl(k_candidates, ac)
  optimalk <- k_candidates[order(silhouette, decreasing = T)[1]]
  pam_result <- pam(distance_matrix, k = optimalk, diss = TRUE, stand = FALSE)
  return (pam_result)
}
