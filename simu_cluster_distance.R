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
source("clust_trial_distance.R")

# simr <- 30
# sce <- 4
# con_sce <- 4
max_tries <- 4


sim_cluster_outlier_distance <- mclapply(1:repeattimes, function (simr) {
  sim_scenario <- lapply(1:length(scenario_type), function(sce) {
    variat <- scenario_type[sce]
    result <- lapply(1: length(contam_types), function (con_sce) {
      cat(sce, "-", con_sce,"\n", sep = "")
      #### generate one simulation scenario
      set.seed(simr)
      scenario <- cluster_generation(nbvar, n, nbtime,
                                     contamination_level, contam_types[con_sce], 
                                     clusternumber, variat)
      listresult <- sparsify_scenario(ps, pc, scenario)
      ##############################==============================
      
      res <- simpleError("Fake error to start")
      counter <- 0
      # Sys.sleep(2*counter)
      while(inherits(res, "error") & counter < max_tries) { 
        res <- tryCatch({ 
          multi_result <- clust_trial_distance(listresult, minprop, alpha_c,
                                               "multivariate", "silh")
        },  error = function(e) e)
        counter <- counter + 1
      }
      
      res <- simpleError("Fake error to start")
      counter <- 0
      # Sys.sleep(2*counter)
      while(inherits(res, "error") & counter < max_tries) { 
        res <- tryCatch({ 
          var1_result <- clust_trial_distance(listresult, minprop, alpha_c, 
                                              "univariate 1", "silh")
        },  error = function(e) e)
        counter <- counter + 1
      }
      
      res <- simpleError("Fake error to start")
      counter <- 0
      # Sys.sleep(2*counter)
      while (inherits(res, "error") & counter < max_tries) { 
        res <- tryCatch({ 
          var2_result <- clust_trial_distance(listresult, minprop, alpha_c,
                                              "univariate 2", "silh")
        },  error = function(e) e)
        counter <- counter + 1
      }
      
      res <- simpleError("Fake error to start")
      counter <- 0
      # Sys.sleep(2*counter)
      while (inherits(res, "error") & counter < max_tries) { 
        res <- tryCatch({ 
          var3_result <- clust_trial_distance(listresult, minprop, alpha_c, 
                                              "univariate 3", "silh")
        },  error = function(e) e)
        counter <- counter + 1
      }
      
      etd_l2_result <- c(multi_result[[1]], var1_result[[1]], 
                      var2_result[[1]], var3_result[[1]])
      etd_l1_result <- c(multi_result[[2]], var1_result[[2]], 
                        var2_result[[2]], var3_result[[2]])
      etd_linf_result <- c(multi_result[[3]], var1_result[[3]], 
                         var2_result[[3]], var3_result[[3]])
      temp <- cbind(etd_l2_result, etd_l1_result, etd_linf_result)
      return (temp)
    }
    )
    return (result)
  })
  return (sim_scenario)
}, mc.cores = 6)


#apply(sim_cluster_outlier_multi_compare, 1, mean)
save(sim_cluster_outlier_diatance, 
     file = paste("sim_cluster_outlier_distance_1009_", 100 * pc, "%.RData", sep = ""))

