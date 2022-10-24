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


sim_cluster_outlier <- mclapply(1:repeattimes, function (simr) {
  sim_scenario <- lapply(1:length(scenario_type), function(sce) {
    variat <- scenario_type[sce]
  result <- sapply(1: length(contam_types), function (con_sce) {
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
        multi_result <- clust_trial(listresult, minprop, alpha_c, "multivariate", "silh")
      },  error = function(e) e)
      counter <- counter + 1
    }
    
    res <- simpleError("Fake error to start")
    counter <- 0
    # Sys.sleep(2*counter)
    while(inherits(res, "error") & counter < max_tries) { 
      res <- tryCatch({ 
        var1_result <- clust_trial(listresult, minprop, alpha_c, "univariate 1", "silh")
      },  error = function(e) e)
      counter <- counter + 1
    }
    
    res <- simpleError("Fake error to start")
    counter <- 0
    # Sys.sleep(2*counter)
    while (inherits(res, "error") & counter < max_tries) { 
      res <- tryCatch({ 
        var2_result <- clust_trial(listresult, minprop, alpha_c, "univariate 2", "silh")
      },  error = function(e) e)
      counter <- counter + 1
    }
    
    res <- simpleError("Fake error to start")
    counter <- 0
    # Sys.sleep(2*counter)
    while (inherits(res, "error") & counter < max_tries) { 
      res <- tryCatch({ 
        var3_result <- clust_trial(listresult, minprop, alpha_c, "univariate 3", "silh")
      },  error = function(e) e)
      counter <- counter + 1
    }
    
    temp <- c(multi_result, var1_result, 
              var2_result, var3_result)
    return (temp)
  }
  )
  return (result)
})
  return (sim_scenario)
}, mc.cores = 6)


#apply(sim_cluster_outlier_multi_compare, 1, mean)
save(sim_cluster_outlier, 
     file = paste("sim_cluster_outlier_0924_", 100 * pc, "%.RData", sep = ""))

# lapply(1: length(sim_cluster_outlier), function(k) {
#   lapply(1:length(contam_types), function(c) {
#   pdf(file = paste("scenario_",k,"_",contam_types[c], ".pdf",sep=""), width = 15.2, height = 4.01)
#   par(mfrow=c(1, 3), mar= c(4.3,4.3,3,1.5))
#   plot(1:repeattimes, sim_cluster_outlier[[k]][[c]][1,], type = "n", ylim=c(0,1), xlab= "Repeat Times", ylab = "ARI", main = paste("Scenario", k,":",scenario_type_title[k], "(",contam_types_title[c],")", sep =" ") )
#     lines(1:repeattimes, sim_cluster_outlier[[k]][[c]][1,], col = 1)
#     lines(1:repeattimes, ifelse(!is.nan(sim_cluster_outlier[[k]][[c]][2,]), sim_cluster_outlier[[k]][[c]][2,], 0), col = 2)
#   
#   
#   plot(1:repeattimes, sim_cluster_outlier[[k]][[c]][3,], type = "n", ylim=range(sim_cluster_outlier[[k]][[c]][3:4, ]), xlab= "Repeat Times", ylab = "Running Time (s)", main = paste("Scenario", k,":",scenario_type_title[k], "(",contam_types_title[c],")", sep =" ") )
#   lines(1:repeattimes, sim_cluster_outlier[[k]][[c]][3,], col = 1)
#   lines(1:repeattimes, ifelse(!is.nan(sim_cluster_outlier[[k]][[c]][4,]), sim_cluster_outlier[[k]][[c]][4,], 0), col = 2)
#    
#   plot(1:repeattimes, sim_cluster_outlier[[k]][[c]][5,], type = "n", ylim=c(0,1), xlab= "Repeat Times", ylab = "Probability", main = paste("Scenario", k,":",scenario_type_title[k], "(",contam_types_title[c],")", sep =" ") )
#   lines(1:repeattimes, sim_cluster_outlier[[k]][[c]][5,], col = "purple")
#   lines(1:repeattimes, ifelse(!is.nan(sim_cluster_outlier[[k]][[c]][6,]), sim_cluster_outlier[[k]][[c]][6,], 0), col = "blue")
#   
#   dev.off()
#   })
# 
#  
# })
