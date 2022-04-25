packages <- c( "gdata")

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

findneighbours <- function(distance_matrix, theta) {
  numer <- distance_matrix[lower.tri(distance_matrix, diag = FALSE)]
  cutoffvalue <- quantile(numer[numer != 0], probs = theta)
  neighbour_matrix <- distance_matrix < cutoffvalue
  return (neighbour_matrix)
}

clustering_warping <- function(neighbour_matrix) {
  ################### group partition ########################
  M <- 0
  S <- 1:nrow(neighbour_matrix)
  group <- list()
  while (length(S) > 0) {
    M <- M + 1
    if (length(S) == 1) {
      element <- S
      group[[M]] <- list(center = S, element = S)
    } else {
      center <- S[which(apply(neighbour_matrix[S, S], 1, sum) 
                      == max(apply(neighbour_matrix[S, S], 1, sum)))]
      if (length(center) > 1) {
        center <- center[1]
      } 
      element <- intersect(which(neighbour_matrix[center, ] == TRUE), S)
      
      group[[M]] <- list(center = center, element = element)
    }
    S <- setdiff(S, element)
  }
##################### group partition  ###################################  

##################### Cluster Partition #################################    
  P <- group
  I <- 0
  cluster_iteration <- list()
  while (length(P) > 0) {
    I <- I + 1
    cluster_iteration[[I]] <- P[[1]]$element
    remove_index <- 1
    neighbour_bound <- c()
    for (l in cluster_iteration[[I]]) {
      neighbour_bound <- union(neighbour_bound, which(neighbour_matrix[l, ] == TRUE))
    }
    if (length(P) > 1) {
      for (l in 2:length(P)) {
        for (cen in P[[l]]$center) {
          if (cen %in% neighbour_bound) {
            cluster_iteration[[I]] <- union(cluster_iteration[[I]], P[[l]]$element)
            remove_index <- c(remove_index, l)
            break
          }
        }
      }
    }
    P <- subsetlist(P, setdiff(1:length(P), remove_index))
  }
##################### Cluster Partition #############################  
  clust <- vector(mode = "list", length = length(cluster_iteration))
  clust <- lapply(1:length(cluster_iteration), 
                  function(l) {
                    element = cluster_iteration[[l]]
                    if (length(element) > 1) {
                      neighbour_within_clust <- apply(neighbour_matrix[element, element], 1, sum)
                      center <- element[which(neighbour_within_clust == max(neighbour_within_clust))]
                    } else {
                      center = element
                    }
                    return (list(center = center, element = element))
                  })
  return (list(group, clust))
}
