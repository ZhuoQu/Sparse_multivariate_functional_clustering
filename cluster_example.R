source("simulation_parameter.R")
source("cluster_generation.R")
source("sparsify_scenario.R")
source("00_subsetlist.R")
source("01_algorithm_distance.R")
source("02_algorithm_clustering.R")
source("03_algorithm_find_outlier.R")
source("find_optimal_theta.R")

packages <- c("randomcoloR")

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

nbtime <- 20
K <- 4
clusterN <- 30
grid <- seq(0, 1, length.out = nbtime * K)
nbvar <- 2
n <- clusterN * K
outlier_index <- sort(sample(1:n, 0.1 * n))
timegrid <- seq(0, 1, length.out = nbtime)
theta <- nbvar * grid
radius <- 10


group <- list() 
group[[1]] <- matrix(rep(radius * sinpi(2 * theta) * cos(2 * pi * grid), clusterN) + runif(length(grid) * clusterN, -1, 1), nrow = clusterN, byrow = TRUE)
group[[2]] <- matrix(rep(radius * sinpi(2 * theta) * sin(2 * pi * grid), clusterN) + runif(length(grid) * clusterN, -1, 1), nrow = clusterN, byrow = TRUE)

cluster_sample <- vector(mode = "list", length = nbvar)
for (clustnumber in 1:K) {
  cluster_sample[[1]] <- rbind(cluster_sample[[1]], 
                               group[[1]][, ((clustnumber - 1) * nbtime + 1):((clustnumber - 1) * nbtime + nbtime)])
  cluster_sample[[2]] <- rbind(cluster_sample[[2]], 
                               group[[2]][, ((clustnumber - 1) * nbtime + 1):((clustnumber - 1) * nbtime + nbtime)])
}

# plot(group[[1]][1, ], group[[2]][1, ], type = "n", 
#      ylim = range(cluster_sample[[2]]), xlim = range(cluster_sample[[1]]),
#      xlab = "Variable 1", ylab = "Variable 2", main = "Bivariate Functional Samples")
# lapply(1:(clusterN * K), function(gg) {
#   for (seg in 1:(nbtime - 1) ) {
#     segments(cluster_sample[[1]][gg, seg], cluster_sample[[2]][gg, seg], 
#              cluster_sample[[1]][gg, seg + 1], cluster_sample[[2]][gg, seg + 1], 
#              col = grey(seg / nbtime) )
#   }
# })

#### introduce outliers 
for (l in 1:nbvar) {
  cluster_sample[[l]][outlier_index[1], ] <- cluster_sample[[l]][outlier_index[1], ] + 5
  cluster_sample[[l]][outlier_index[2], ] <- cluster_sample[[l]][outlier_index[2], ] - 5
  cluster_sample[[l]][outlier_index[3], 4:7] <- cluster_sample[[l]][outlier_index[3], 4:7] + 5
  cluster_sample[[l]][outlier_index[4], 10:13] <- cluster_sample[[l]][outlier_index[4], 10:13] - 5
  cluster_sample[[l]][outlier_index[5], 14:20] <- cluster_sample[[l]][outlier_index[5], 14:20] + 5
  cluster_sample[[l]][outlier_index[6], 12:20] <- cluster_sample[[l]][outlier_index[6], 12:20] - 5
  a <- range(cluster_sample[[l]][outlier_index[7], ])[1]
  b <- range(cluster_sample[[l]][outlier_index[7], ])[2] - a
  cluster_sample[[l]][outlier_index[7], ] <- a + b * timegrid
  a <- range(cluster_sample[[l]][outlier_index[8], ])[1]
  b <- range(cluster_sample[[l]][outlier_index[8], ])[2] - a
  cluster_sample[[l]][outlier_index[8], ] <- a + b * timegrid
  h1 <- max(abs(range(cluster_sample[[l]][outlier_index[9], ]))) / 2
  h2 <- max(abs(range(cluster_sample[[l]][outlier_index[10], ]))) / 2
  h3 <- max(abs(range(cluster_sample[[l]][outlier_index[11], ]))) / 2
  h4 <- max(abs(range(cluster_sample[[l]][outlier_index[12], ]))) / 2
  if (l == 1) {
    cluster_sample[[l]][outlier_index[9], ] <- h1 * cos(pi * timegrid)
    cluster_sample[[l]][outlier_index[10], ] <- h2 * cos(pi * timegrid)
    cluster_sample[[l]][outlier_index[11], ] <- h3 * cos(30 * pi * timegrid)
    cluster_sample[[l]][outlier_index[12], ] <- h4 * cos(30 * pi * timegrid)
    } else if (l == 2) {
    cluster_sample[[l]][outlier_index[9], ] <- h1 * sin(pi * timegrid)
    cluster_sample[[l]][outlier_index[10], ] <- h2 * sin(pi * timegrid)
    cluster_sample[[l]][outlier_index[11], ] <- h3 * sin(30 * pi * timegrid)
    cluster_sample[[l]][outlier_index[12], ] <- h4 * sin(30 * pi * timegrid)
  }
}

############ sparsify data #################################
scenario <- list(cluster_sample, outlier_index)
sparse_scenario <- sparsify_scenario(0.8, 0.3, scenario)
###########################################################

##########################################
etd <- distance_mfd(sparse_scenario[[1]])
##### First of all we need to optimize parameters
simulation_result <- find_optimal_theta(etd, standard_label = NULL, minprop = 0.05, 
                                        optim_index = "silhouette", alpha_c = 0.87) 
optimal_theta <- simulation_result$optimal_theta
#################### Second, we implement Algorithm 2 under the optimal theta.
sample_neighbour <- findneighbours(etd, optimal_theta)
clust <- clustering_warping(sample_neighbour)
##################### Third, we implement Algorithm 3 under the optimal theta.
final_clust <- include_isolatepoints(clust[[2]], 0.05,
                                     etd, 0.87)
#################################################
######## return the silhouette value for the group #############################################
clusterlabel <- function(final_cluster) {
  clust <- final_cluster$clust
  isolate <- final_cluster$isolate
  empirical_label <- rep(0, nrow(sample_neighbour))
  for (cl in 1:length(clust)) {
    empirical_label[clust[[cl]]$element] <- cl
  }
  empirical_label[isolate] <- (length(clust) + 1):(length(clust) + length(isolate))
  
  return (empirical_label)
}

empirical_label <- clusterlabel(final_clust)
silh_value <- silhouette_index(empirical_label,
                               etd, minprop)

###################################################
pdf("visualize_1.pdf", width = 4.8, height = 5)
par(mfrow = c(1, 1), mai = c(0.7, 0.8, 0.4, 0.1), 
    mar = c(4, 3.5, 2, 0.5), mgp = c(2.5, 1.1, 0))
plot(group[[1]][1, ], group[[2]][1, ], type = "n",
     ylim = range(cluster_sample[[2]]), xlim = range(cluster_sample[[1]]),
     xlab = "Variable 1", ylab = "Variable 2",
     main = bquote(bold(paste("(a) Functional Samples (", .(K)," clusters, ", .(length(outlier_index))," outliers)", sep = ""))),
     cex.main = 1.1, cex.lab = 1.1)
lapply(setdiff(1:(clusterN * K), outlier_index), function(gg) {
  nbtime <- length(sparse_scenario[[2]][[gg]]$argvals)
  for (seg in 1:(nbtime - 1) ) {
    segments(sparse_scenario[[2]][[gg]]$y[seg], sparse_scenario[[3]][[gg]]$y[seg],
             sparse_scenario[[2]][[gg]]$y[seg + 1], sparse_scenario[[3]][[gg]]$y[seg + 1],
             col = grey(seg / nbtime) )
  }
})

lapply(outlier_index, function(gg) {
  nbtime <- length(sparse_scenario[[2]][[gg]]$argvals)
  for (seg in 1:(nbtime - 1) ) {
    segments(sparse_scenario[[2]][[gg]]$y[seg], sparse_scenario[[3]][[gg]]$y[seg],
             sparse_scenario[[2]][[gg]]$y[seg + 1], sparse_scenario[[3]][[gg]]$y[seg + 1],
             col = rgb(1, 140/250, 0, alpha = seg / nbtime), lty = 2)
  }
})

dev.off()

pdf("visualize_2.pdf", width = 4.8, height = 5)
par(mfrow = c(1, 1), mai = c(0.7, 0.8, 0.4, 0.1), 
    mar = c(4, 3.5, 2, 0.5), mgp = c(2.5, 1.1, 0))

plot(simulation_result$simulation$theta, simulation_result$simulation$value_index, 
     xlab = expression(theta), ylab = "Average Silhouette value", type = "l", 
     xlim = c(0, 0.3), ylim = c(0, 1), font.main = 2,
     main = expression(bold(paste("(b) Average Silhouette value versus ", theta, sep = ""))), 
     cex.main = 1.1, cex.lab = 1.1)
points(simulation_result$simulation$theta, simulation_result$simulation$value_index, cex = 0.5)
points(simulation_result$optimal_theta, max(simulation_result$simulation$value_index), col = "red", pch = 2)
text(simulation_result$optimal_theta, max(simulation_result$simulation$value_index) + 0.05,
     labels = bquote(paste(theta*'='*.(simulation_result$optimal_theta))), cex = 1.2)
dev.off()
######### plot groups from Criteria 1 in Algorithm 2 #############################################
figure_color <- c("black", "red", "green", "blue", "orange")
figure_color <- c(figure_color, distinctColorPalette(length(clust[[1]]) - 5))

pdf("visualize_3.pdf", width = 4.8, height = 5)
par(mfrow = c(1, 1), mai = c(0.7, 0.8, 0.4, 0.1), 
    mar = c(4, 3.5, 2, 0.5), mgp = c(2.5, 1.1, 0))
plot(group[[1]][1, ], group[[2]][1, ], type = "n",
     ylim = range(cluster_sample[[2]]), xlim = range(cluster_sample[[1]]),
     xlab = "Variable 1", ylab = "Variable 2", 
     main = bquote(bold(paste("(c) The first-layer Partition ", bolditalic(G), " (", .(length(clust[[1]])), " groups)", sep = ""))),
     cex.main = 1.1, cex.lab = 1.1)

# lapply(1:length(clust[[1]]), function(gg) {
#   if (length(clust[[1]][[gg]]$element) == 1) {
#     samp <- clust[[1]][[gg]]$element
#     points(sparse_scenario[[2]][[samp]]$y, sparse_scenario[[3]][[samp]]$y,
#            col = "gold", lty = gg, type = "l")
#   }
# })
# 
# lapply(1:length(clust[[1]]), function(gg) {
#   if (length(clust[[1]][[gg]]$element) == 2) {
#     for (samp in clust[[1]][[gg]]$element) {
#       points(sparse_scenario[[2]][[samp]]$y, sparse_scenario[[3]][[samp]]$y,
#              col = "purple", lty = gg, type = "l")
#     }
#   }
# })

lapply(1:length(clust[[1]]), function(gg) {
  #if (length(clust[[1]][[gg]]$element) > 2) {
    for (samp in clust[[1]][[gg]]$element) {
      points(sparse_scenario[[2]][[samp]]$y, sparse_scenario[[3]][[samp]]$y,
             col = figure_color[gg], lty = 1, type = "l")
    }
  #}
})

dev.off()

figure_color <- c("black", "red", "green", "blue", "orange")
figure_color <- c(figure_color, distinctColorPalette(length(clust[[2]]) - 5))
pdf("visualize_4.pdf", width = 4.8, height = 5)
par(mfrow = c(1, 1), mai = c(0.7, 0.8, 0.4, 0.1), 
    mar = c(4, 3.5, 2, 0.5), mgp = c(2.5, 1.1, 0))
######## plot clusters from Criteria 2 in Algorithm 2 #####################################################
plot(group[[1]][1, ], group[[2]][1, ], type = "n",
     ylim = range(cluster_sample[[2]]), xlim = range(cluster_sample[[1]]),
     xlab = "Variable 1", ylab = "Variable 2", 
     main = bquote(bold(paste("(d) The second-layer Partition ", bolditalic(C), " (", .(length(clust[[2]])), " groups)", sep = ""))),
     cex.main = 1.1, cex.lab = 1.1)


lapply(1:length(clust[[2]]), function(gg) {
  #if (length(clust[[2]][[gg]]$element) > 2) {
    for (samp in clust[[2]][[gg]]$element) {
      points(sparse_scenario[[2]][[samp]]$y, sparse_scenario[[3]][[samp]]$y,
             col = figure_color[gg], lty = 1, type = "l")
    }
  #}
})
# lapply(1:length(clust[[2]]), function(gg) {
#   if (length(clust[[2]][[gg]]$element) == 1) {
#     samp <- clust[[2]][[gg]]$element
#     points(sparse_scenario[[2]][[samp]]$y, sparse_scenario[[3]][[samp]]$y,
#            col = "gold", lty = gg, type = "l")
#   }
# })

# lapply(1:length(clust[[2]]), function(gg) {
#   if (length(clust[[2]][[gg]]$element) == 2) {
#     for (samp in clust[[2]][[gg]]$element) {
#       points(sparse_scenario[[2]][[samp]]$y, sparse_scenario[[3]][[samp]]$y,
#              col = "purple", lty = gg, type = "l")
#     }
#   }
# })
dev.off()

########## plot primary clusters and potential outliers from Algorithm 3####################################################
primary_cluster_index <- sapply(clust[[2]], function(l) {length(l$element)}) > n * minprop
primary_cluster_index <- which(primary_cluster_index == TRUE)
potential_outlier <- unlist(lapply(subsetlist(clust[[2]], setdiff(1:length(clust[[2]]), primary_cluster_index)), 
                            function(k) {k$element}))
#figure_color <- distinctColorPalette(length(primary_cluster_index) + 1)
figure_color <- c("black", "red", "green", "blue", "orange")
pdf("visualize_5.pdf", width = 4.8, height = 5)
par(mfrow = c(1, 1), mai = c(0.7, 0.8, 0.4, 0.1), 
    mar = c(4, 3.5, 2, 0.5), mgp = c(2.5, 1.1, 0))
plot(group[[1]][1, ], group[[2]][1, ], type = "n",
     ylim = range(cluster_sample[[2]]), xlim = range(cluster_sample[[1]]),
     xlab = "Variable 1", ylab = "Variable 2", 
     main = paste("(e) Primary Clusters (", length(primary_cluster_index),") and Potential Outliers (", length(potential_outlier), ")", sep = ""),
     cex.main = 1.02, cex.lab = 1.1)

lapply(1:length(clust[[2]]), function(gg) {
  if (gg %in% primary_cluster_index) {
    for (samp in clust[[2]][[gg]]$element) {
      points(sparse_scenario[[2]][[samp]]$y, sparse_scenario[[3]][[samp]]$y,
             col = figure_color[gg], lty = 1, type = "l")
    }
  }
})

lapply(1:length(clust[[2]]), function(gg) {
  if (gg %in% primary_cluster_index == FALSE) {
    for (samp in clust[[2]][[gg]]$element) {
      points(sparse_scenario[[2]][[samp]]$y, sparse_scenario[[3]][[samp]]$y,
             col = rev(figure_color)[1], lty = 2, type = "l")
    }
  }
})

dev.off()

#######################Primary clusters and outliers from Algorithm 3 ##############################################
pdf("visualize_6.pdf", width = 4.8, height = 5)
par(mfrow = c(1, 1), mai = c(0.7, 0.8, 0.4, 0.1), 
    mar = c(4, 3.5, 2, 0.5), mgp = c(2.5, 1.1, 0))
plot(group[[1]][1, ], group[[2]][1, ], type = "n",
     ylim = range(cluster_sample[[2]]), xlim = range(cluster_sample[[1]]),
     xlab = "Variable 1", ylab = "Variable 2", 
     main = paste("(f) Primary Clusters (", length(final_clust$clust), ")", " and Outliers (", length(final_clust$isolate), ")", sep = ""),
     cex.main = 1.1, cex.lab = 1.1)
lapply(1:length(final_clust$clust), function(gg) {
  for (samp in final_clust$clust[[gg]]$element) {
    points(sparse_scenario[[2]][[samp]]$y, sparse_scenario[[3]][[samp]]$y,
           col = figure_color[gg], lty = 1,
           type = "l")
  }
})
lapply(final_clust$isolate, function(gg) {
  points(sparse_scenario[[2]][[gg]]$y, sparse_scenario[[3]][[gg]]$y,
        col = rev(figure_color)[1], lty = 2, type = "l")

})
dev.off()
###################################################################
