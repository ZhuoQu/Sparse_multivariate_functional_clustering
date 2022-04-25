source("00_subsetlist.R")
include_isolatepoints <- function(cluster_result, minprop, distance_matrix, alpha_c = 0.85) {
  ### cluster_result is the list of the clustering function. Each list includes center, core, and element
  ###################################################
  basis_cluster_index <- which(lapply(cluster_result, function(k) { 
    length(k$element) > minprop * nrow(distance_matrix)}) == T)
  
  isolatecenter <- subsetlist(cluster_result, setdiff(1:length(cluster_result), basis_cluster_index))
  
  maincenter <- subsetlist(cluster_result, basis_cluster_index)
  
  if (length(isolatecenter) == 0) {
    return (list(clust = maincenter))
  } else if (length(maincenter) == 0) {
    return (list(isolate = isolatecenter))
  } else {
    Distance_to_gravity <- lapply(maincenter, function(l) {
    distance_matrix[l$center[1], l$element] })
    
    distance_threshold <- sapply(Distance_to_gravity, function(l) {
      sort(l)[ceiling(alpha_c * length(l))]
    })
    isolate <- c()
    ###############################################
    if (length(maincenter) > 1) {
      D_g <- t(sapply(isolatecenter, function(nr) {
        sapply(maincenter, function(nc) {
          return (distance_matrix[nr$center[1], nc$center[1]])
        })
      }))
      D_g_bool <- t(apply(D_g, 1, function(k) {
        temp <- c()
        for (l in 1:length(k)) {
          temp <- c(temp, k[l] <= distance_threshold[l])
        }
        return (temp)
      }))
    } else {
      D_g <- sapply(isolatecenter, function(nr) {
      return (distance_matrix[nr$center[1], maincenter[[1]]$center[1]])
     }) 
      D_g <- as.matrix(D_g)
      D_g_bool <- D_g <= distance_threshold
    }
    rownames(D_g) <- lapply(isolatecenter, function(k) {
      return (k$center[1])
    })
    for (l in 1:nrow(D_g)) {
      if (length(unique(D_g_bool[l, ])) == 1 && 
          unique(D_g_bool[l, ]) == FALSE) {
        isolate <- c(isolate, isolatecenter[[l]]$element)
      } else {
        cand_index <- which(D_g_bool[l, ] == TRUE)
        main_quantile <- sapply(cand_index, function(num) {
          return (rank(c(D_g[l, num], Distance_to_gravity[[num]])[1] / 
                         (length(Distance_to_gravity[[num]]) + 1) ) )
        })
        index <- which(main_quantile == min(main_quantile))
        if (length(index) > 1) {
          index <- which(D_g[l, ] == min(D_g[l, ]))
        }
        maincenter[[index]]$element <- c(maincenter[[index]]$element, 
                                         isolatecenter[[index]]$element)
      }
    }
    s <- rbind(D_g, distance_threshold)
    return (list(clust = maincenter,
                 isolate = sort(isolate),
                 criteria_matrix = s))
  }
}



#############################################
# plot(Distance_to_gravity[[1]], rank(Distance_to_gravity[[1]]) / length(Distance_to_gravity[[1]]),
#      xlim = range(c(unlist(Distance_to_gravity), D_g)), ylim = c(0, 1), type = "n", xlab ="Distance to the Center", 
#      ylab ="Empirical cdf")
# if (length(maincenter) > 1) {
#   lapply(1:length(maincenter), function(l) {
#   candidate_curve <- c(Distance_to_gravity[[l]], D_g[, l])
#   points(candidate_curve, rank(candidate_curve) / length(1 + candidate_curve), col = l)
#   #points(D_g[,l],rank(D_g[,l])/length(1+candidate_curve),col=l, pch=4)
# })
# } else {
#   candidate_curve <- c(Distance_to_gravity[[1]], D_g)
#   points(candidate_curve, rank(candidate_curve) / length(1 + candidate_curve))
#   
# }
#############################################
