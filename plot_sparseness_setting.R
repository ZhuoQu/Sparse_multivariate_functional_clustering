source("simulation_parameter.R")
source("cluster_generation.R")
source("sparsify_scenario.R")
packages <- c("scatterplot3d")

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
fig_color <- c("black", "red", "green", "orange")
plot_sparseness_setting <- function(con, rotate_angle, sc) {
  ## scenario is the result of function cluster_generation
  par(mfrow = c(1, 3))
  var <- generate_cluster(nbvar, n, nbtime, clusternumber, variat)
  outlier_nb <- sample(n, 0.07 * n, replace = FALSE)
  var_outl <- introduce_outlier(nbvar, nbtime, t, var, outlier_nb, con)
    ###########################################  
  scenario <- list(var_outl, outlier_nb)
  var <- scenario[[1]]
  nbvar <- length(var)
  outlier_nb <- scenario[[2]]
    
  ps <- 1
    #scenario_sparse <- sparsify_scenario(ps, pc, scenario)
  sparse_curve <- sample(1:n, ps * n)
  complete_curve <- setdiff(1:n, sparse_curve)
  nbvar <- length(scenario[[1]])
  outlier_nb <- scenario[[2]]
  if (nbvar == 3) {
    lapply(c(0, 0.3, 0.6), function (pc) {
    s3d <- scatterplot3d(x = scenario[[1]][[1]][1, ],
                         y = scenario[[1]][[2]][1, ],
                         z = scenario[[1]][[3]][1, ],
                         xlim = range(scenario[[1]][[1]]),
                         ylim = range(scenario[[1]][[2]]),
                         zlim = range(scenario[[1]][[3]]),
                         type = "n",
                         xlab = "Variable 1",
                         ylab = "Variable 2",
                         zlab = "Variable 3",
                         main = bquote(p[size] == .(100 * ps)*"% and"~p[curve] == .(100 * pc) * "%"), 
                         grid = TRUE, scale.y = sc, 
                         angle = rotate_angle, cex.main = 2, 
                         cex.lab = 1.3, cex.axis = 1,
                         mar = c(3, 3.5, 3, 1.7),
                         mgp = c(2.2, 0.5, 0))
        for (cc in 1:clusternumber) {
          if (cc < clusternumber) {
            normal_curve_index <- setdiff ((round((cc-1) * n / clusternumber) + 1):round((cc) * n / clusternumber), 
                                         outlier_nb)
          } else {
            normal_curve_index <- setdiff(((clusternumber - 1) * round(n / clusternumber) + 1):n, outlier_nb)
          }
          lapply(normal_curve_index, function(rowindex) {
            s3d$points3d(var[[1]][rowindex, ], var[[2]][rowindex, ], var[[3]][rowindex, ], 
                         type = "l", col = figure_color[cc])
            if (rowindex %in% sparse_curve) {
              sparse_index <- sample(1:nbtime, pc * nbtime)
              s3d$points3d(var[[1]][rowindex, sparse_index], 
                           var[[2]][rowindex, sparse_index], 
                           var[[3]][rowindex, sparse_index], 
                           type = "p", col = "white", pch = 16, cex = 0.8)
            }
          })
          if (length(outlier_nb) > 0) {
            lapply(outlier_nb, function (ll) {
              s3d$points3d(var[[1]][ll, ], var[[2]][ll, ], var[[3]][ll, ], 
                           col = "orange", 
                           lty = 2, type = "l", lwd = 1.5)
              if (ll %in% sparse_curve) {
                sparse_index <- sample(1:nbtime, pc*nbtime)
                s3d$points3d(var[[1]][ll, sparse_index], var[[2]][ll, sparse_index], 
                             var[[3]][ll, sparse_index], col = "white", 
                             type = "p", pch = 16, cex = 0.8)
              }
            })
            #text(var[[1]][outlier_nb,1],var[[2]][outlier_nb,1],var[[3]][outlier_nb,1],labels=outlier_nb,col="purple",cex=0.7)
          }
         } 
      })
    }
  invisible(nbvar)
}

rotate_angle <- c(45, 70, 45, 45, 160, 45)
sc <- c(1, 1, 1, 1, 1, 0.5)
for (l in 1:6) {
  variat <- scenario_type[l]
  con <- contam_types[1]
  pdf(file = paste("sparseness_sce_", l, ".pdf", sep = ""), 
      width = 12.17, height = 4)
  plot_sparseness_setting (con, rotate_angle[l], sc[l])
  dev.off()
}
