packages = c("NbClust", "cluster")
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

clust_hierarchical_distance <- function(distance_matrix) {
m <- c("average", "single", "complete", "ward", "weighted","flexible", "gaverage")
names(m) <- c("average", "single", "complete", "ward", "weighted","flexible", "gaverage")

# function to compute coefficient
ac <- function(x) {
  agnes(distance_matrix, method = x, par.method = c(0.3,0.3,0.4,0))$ac
}

distance_method <- m[which(map_dbl(m, ac) == max(map_dbl(m,ac)))]
method_Clust <- ifelse(distance_method %in% c("ward", "weighted", "flexible", "gaverage"), "ward.D", distance_method)
optimalresult <- NbClust(diss = dist(distance_matrix, diag = FALSE, upper = FALSE), distance = NULL, 
                         min.nc = 2, max.nc = sqrt(n / 2) + 5, method = method_Clust, index = "silhouette")

return (optimalresult)
}

