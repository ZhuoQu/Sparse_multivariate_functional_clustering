source("cluster_generation.R")
source("simulation_parameter.R")

packages <- c("scatterplot3d")
package_check <- lapply(
  packages,
  FUN <- function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

fig_color <- c(rgb(0, 0, 0, maxColorValue = 255), 
               rgb(255, 0, 0, maxColorValue = 255), 
               rgb(0, 128, 0, maxColorValue = 255))
figure_color <- c("black", "red", "green", "orange")
plot_clust <- function(result, rotate_angle, sc) {
  ## result is the result of function cluster_generation
  var <- result[[1]]
  nbvar <- length(var)
  outlier_nb <- result[[2]]
  par(mfrow = c(1, ifelse(nbvar <= 3, nbvar + 1, 3)), 
      mar = c(4.5, 3.8, 3, 1),
      mgp = c(2.6, 1, 0))
  for (k in 1:nbvar) {
    plot(t, var[[k]][1, ], ylim = range(var[[k]]), type = "n", 
         ylab = "Values", xlab = "Time", main = paste("Variable ", k), 
         cex.main = 2, cex.lab = 2, cex.axis = 2)
    lapply(setdiff(1:n, outlier_nb), function (ll) {
      for (cc in 1:clusternumber) {
        if (cc < clusternumber && 
            round((cc-1) * n / clusternumber) < ll && 
            ll <= round((cc) * n / clusternumber)) {
          lines(t, var[[k]][ll, ], col = figure_color[cc])
        } else if (cc == clusternumber && 
                  (clusternumber - 1) * round(n / clusternumber) < ll) {
          lines(t, var[[k]][ll, ], col = figure_color[cc])
        } 
      }
    })
    if (length(outlier_nb) > 0) {
      apply(var[[k]][outlier_nb, ], 1, function (ll) {
        lines(t, ll, 
              col = "orange", 
              lty = 2)
      })
      text(t[1], var[[k]][outlier_nb, 1], 
           labels = outlier_nb, 
           col = rgb(128, 0, 128),
           #col = "purple", 
           cex = 0.7)
    }
  }
  if (nbvar == 3) {
    s3d <- scatterplot3d(x = var[[1]][1, ],
                         y = var[[2]][1, ],
                         z = var[[3]][1, ],
                         xlim = range(var[[1]]),
                         ylim = range(var[[2]]),
                         zlim = range(var[[3]]),
                         type = "n",
                         xlab = "Variable 1",
                         ylab = "Variable 2",
                         zlab = "Variable 3",
                         main = "Trivarate Variables", 
                         grid = TRUE, scale.y = sc, 
                         angle = rotate_angle, cex.main = 2, 
                         cex.lab = 1.6, cex.axis = 1.3,
                         mar = c(3, 3.5, 3, 1.7),
                         mgp = c(2.2, 0.5, 0))
    for (cc in 1:clusternumber) {
      if (cc < clusternumber) {
        normal_curve_index <- setdiff ((round((cc - 1) * n / clusternumber) + 1):round((cc) * n / clusternumber), 
                                       outlier_nb)
      } else {
        normal_curve_index <- setdiff(((clusternumber - 1) * round(n / clusternumber) + 1):n, outlier_nb)
      }
      lapply(normal_curve_index, function(rowindex) {
        s3d$points3d(var[[1]][rowindex, ], 
                     var[[2]][rowindex, ], 
                     var[[3]][rowindex, ], type = "l", 
                     col = figure_color[cc])
      })
    }
    if (length(outlier_nb) > 0) {
      lapply(outlier_nb, function (ll) {
        s3d$points3d(var[[1]][ll, ], var[[2]][ll, ],var[[3]][ll, ], 
                     #col = "orange", 
                     #col = rgb(255, 165, 0, maxColorValue = 255),
                     col = figure_color[4],
                     lty = 2, type = "l")
      })
      #text(var[[1]][outlier_nb,1],var[[2]][outlier_nb,1],var[[3]][outlier_nb,1],labels=outlier_nb,col="purple",cex=0.7)
    }
  }
  invisible(nbvar)
}


rotation_angle <- c(45, 70, 45, 45, 160, 45)
scaley <- c(1, 1, 1, 1, 1, 0.5)
lapply(1:length(scenario_type), function(k) {
  pdf(file = paste("scenario_", k,".pdf", sep = ""), width = 15.2, height = 4.01)
  cat(k, "\n")
  variat <- scenario_type[k]
  contam <- contam_types[1]
  scenario <- cluster_generation(3, 150, 50, 
                                 0, 
                                 contam, 
                                 3, variat)
  plot_clust(scenario, rotation_angle[k], scaley[k])
  dev.off()
})
