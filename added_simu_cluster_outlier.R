packages = c( "caret", "fda",
              "funHDDC", "mclust", "parallel")


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
#setwd("~/Dropbox/Mac/Desktop/KAUST/project/Clustering/finalcode")
source("cross_covariance.R")
source("cluster_generation.R")
source("simulation_parameter.R")
source("clust_trial.R")

# simr <- 30
# sce <- 4
# con_sce <- 4
max_tries <- 4
pc <- 0.6

added_sim_cluster_outlier <- lapply(1:repeattimes, function (simr) {
  sim_scenario <- lapply(4:5, function(sce) {
    variat <- scenario_type[sce]
    result <- sapply(5:6, function (con_sce) {
      cat(sce, "-", con_sce,"\n", sep = "")
      #### generate one simulation scenario
      set.seed(simr)
      scenario <- cluster_generation(nbvar, n, nbtime,
                                     contamination_level, contam_types[con_sce], 
                                     clusternumber, variat)
      ##############################==============================
      
      res <- simpleError("Fake error to start")
      counter <- 0
      # Sys.sleep(2*counter)
      while(inherits(res, "error") & counter < max_tries) { 
        res <- tryCatch({ 
          multi_result <- clust_trial(scenario, pc, minprop, alpha_c, "multivariate", "silh")
        },  error = function(e) e)
        counter <- counter + 1
      }
      
      res <- simpleError("Fake error to start")
      counter <- 0
      # Sys.sleep(2*counter)
      while(inherits(res, "error") & counter < max_tries) { 
        res <- tryCatch({ 
          var1_result <- clust_trial(scenario, pc, minprop, alpha_c, "univariate 1", "silh")
        },  error = function(e) e)
        counter <- counter + 1
      }
      
      res <- simpleError("Fake error to start")
      counter <- 0
      # Sys.sleep(2*counter)
      while (inherits(res, "error") & counter < max_tries) { 
        res <- tryCatch({ 
          var2_result <- clust_trial(scenario, pc, minprop, alpha_c, "univariate 2", "silh")
        },  error = function(e) e)
        counter <- counter + 1
      }
      
      res <- simpleError("Fake error to start")
      counter <- 0
      # Sys.sleep(2*counter)
      while (inherits(res, "error") & counter < max_tries) { 
        res <- tryCatch({ 
          var3_result <- clust_trial(scenario, pc, minprop, alpha_c, "univariate 3", "silh")
        },  error = function(e) e)
        counter <- counter + 1
      }
      
      temp <- c(multi_result, var1_result, var2_result, var3_result)
      return (temp)
    }
    )
    return (result)
  })
  return (sim_scenario)
})


#apply(sim_cluster_outlier_multi_compare, 1, mean)
pc = 0.3
load(paste("added_sim_cluster_outlier_0924_", 100 * pc, "%.RData", sep = ""))
load(paste("sim_cluster_outlier_0924_", 100 * pc, "%.RData", sep = ""))
sim_cluster_outlier
for (k in 1:repeattimes) {
  if (length(sim_cluster_outlier[[k]]) > 1) {
  sim_cluster_outlier[[k]][[4]][, c(5, 6)] <- added_sim_cluster_outlier[[k]][[1]]
  sim_cluster_outlier[[k]][[5]][, c(5, 6)] <- added_sim_cluster_outlier[[k]][[2]]
  }
}

save(sim_cluster_outlier, file = paste("sim_cluster_outlier_0924_", 100 * pc, "%.RData", sep = ""))
