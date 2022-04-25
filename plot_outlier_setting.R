source("simulation_parameter.R")
source("cluster_generation.R")
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

figure_color <- c("black", "red", "green", "orange")
plot_outlier_setting <- function(rotate_angle, sc) {
  ## scenario is the result of function cluster_generation
  par(mfrow = c(1, 3), mai = c(4, 4, 4, 1))
  outlier_nb <- sample(1:n, n * 0.1)
  t <- seq(0, 1, length.out = nbtime)
  var <- generate_cluster(nbvar, n, nbtime, clusternumber, variat)

  lapply(1:6, function (con) {
  contam <- contam_types[con]
  var_outl <- introduce_outlier(nbvar, nbtime, t,
                                var, outlier_nb, contam)
  # if (con %in% 1:3) {
  #   xrange <- c(-10, 15)
  #   yrange <- c(-15, 10)
  #   zrange <- c(-15, 15)
  # } else {
  #   xrange <- c(-6, 8)
  #   yrange <- c(-8, 8)
  #   zrange <- c(-10, 10)
  # }
  ###########################################  
  scenario <- list(var_outl, outlier_nb)
  var <- scenario[[1]]
  nbvar <- length(var)
  outlier_nb <- scenario[[2]]
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
                         main = paste("Contamination ", con), 
                         grid = TRUE, scale.y = sc, 
                         angle = rotate_angle, 
                         cex.main = 1.6, 
                         cex.lab = 1.1, cex.axis = 0.9, 
                         mar = c(3, 3.5, 3, 1.7),
                         mgp = c(2.2, 0.5, 0))
    #text(x = 0, y = min(var[[2]]), "Variable 1", srt = 45)
   # text(x = max(var[[1]]) + 1, y = 0, "Variable 2", srt = 45, cex = 1.6)
    for (cc in 1:clusternumber) {
      if (cc < clusternumber) {
        normal_curve_index <- setdiff ((round((cc - 1) * n / clusternumber) + 1):round((cc) * n / clusternumber), 
                                       outlier_nb)
        lapply(normal_curve_index, function(rowindex) {
          s3d$points3d(var[[1]][rowindex, ], var[[2]][rowindex, ], var[[3]][rowindex, ], 
                       type = "l", col = figure_color[cc])
        })
      } else {
        normal_curve_index <- setdiff(((clusternumber - 1) * round(n / clusternumber)+1):n, outlier_nb)
        lapply(normal_curve_index, function(rowindex) {
          s3d$points3d(var[[1]][rowindex,], var[[2]][rowindex,], var[[3]][rowindex,], 
                       type = "l", col = figure_color[cc])
        }) 
      }
    }
    
    if (length(outlier_nb) > 0) {
      #outl_ind <- 1:length(outlier_nb)
      outl_ind <- 1:length(outlier_nb)
      if (con == 6) {
        outl_ind <- 1:5
      }
      lapply(outlier_nb[outl_ind], function (ll) {
        s3d$points3d(var[[1]][ll, ], var[[2]][ll, ],var[[3]][ll, ], 
                     col = "orange", type = "l", lwd = 1.3, lty = 2)
      })
      #text(var[[1]][outlier_nb,1],var[[2]][outlier_nb,1],var[[3]][outlier_nb,1],labels=outlier_nb,col="purple",cex=0.7)
    }
  }
  invisible(nbvar)
  })
}

rotate_angle <- c(45, 70, 45, 45, 160, 45)
sc <- c(1, 1, 1, 1, 1, 0.5)
for (l in 1:6) {
  variat <- scenario_type[l]
  plot_outlier_setting (rotate_angle[l], sc[l])
}
