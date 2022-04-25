library("gdata")
source("00_subsetlist.R")
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

ps <- 1
plotsimulation_meanmarginal <- function(sce, pc, index, includelegend, savefile = TRUE) { 
  if (savefile == TRUE && includelegend == TRUE) {
    pdf(file = paste("sce", sce, "_pc", 100 * pc, ".pdf", sep = ""), width = 11, height = 5.2)
  } else if (savefile == TRUE) {
    pdf(file = paste("sce", sce, "_pc", 100 * pc, ".pdf", sep = ""), width = 11, height = 4.95)
  }
  sequence <- unlist(sapply(1:130, function(k){
    ifelse(length(sim_cluster_outlier[[k]][[sce]]) > 1, return (k), return())}))
  percent <- list()
  for (j in 1:6) {
    available_list <- lapply(subsetlist(sim_cluster_outlier, sequence), function(k) {
      k[[sce]]
    } )
    if (index == "ARI") {
      mult_dbc <- sapply(available_list, function(mat) {
        mat[1, j]})
      mult_hierarchical <- sapply(available_list, function(mat) {
        mat[2, j]})
      mult_kmedoids <- sapply(available_list, function(mat) {
        mat[3, j]})
      mult_funhddc <- sapply(available_list, function(mat) {
        mat[4, j]})
      
      uni_dbc <- sapply(available_list, function(mat) {
        mean(mat[c(11, 21, 31), j])})
      uni_hierarchical <- sapply(available_list, function(mat) {
        mean(mat[c(12, 22, 32), j])})
      uni_kmedoids <- sapply(available_list, function(mat) {
        mean(mat[c(13, 23, 33), j])})
      uni_funhddc <- sapply(available_list, function(mat) {
        mean(mat[c(14, 24, 34), j])})
      
      multi <- unmatrix(sapply(available_list, function(mat) {
        mat[1:4, j]}), byrow = T)
      uni <- unmatrix(sapply(available_list, function(mat) {
        c(mean(mat[c(11, 21, 31), j]), mean(mat[c(12, 22, 32), j]), 
          mean(mat[c(13, 23, 33), j]), mean(mat[c(14, 24, 34), j]))
      }), byrow = T)
      
      
      ### multi_DBC, uni_DBC, multifunHDDC, unifunHDDC, multihierarchical, unihierarchical, multikmedoids,unikmedoids
    } else if (index == "time") {
      mult_dbc <- sapply(available_list, function(mat) {
        mat[5, j]})
      mult_hierarchical <- sapply(available_list, function(mat) {
        mat[6, j]})
      mult_kmedoids <- sapply(available_list, function(mat) {
        mat[7, j]})
      mult_funhddc <- sapply(available_list, function(mat) {
        mat[8, j]})
      
      uni_dbc <- sapply(available_list, function(mat) {
        mean(mat[c(15, 25, 35), j])})
      uni_hierarchical <- sapply(available_list, function(mat) {
        mean(mat[c(16, 26, 36), j])})
      uni_kmedoids <- sapply(available_list, function(mat) {
        mean(mat[c(17, 27, 37), j])})
      uni_funhddc <- sapply(available_list, function(mat) {
        mean(mat[c(18, 28, 38), j])})
      
      multi <- unmatrix(sapply(available_list, function(mat) {
        mat[5:8, j]}), byrow = T)
      uni <- unmatrix(sapply(available_list, function(mat) {
        c(mean(mat[c(15, 25, 35), j]), mean(mat[c(16, 26, 36), j]), 
          mean(mat[c(17, 27, 37), j]), mean(mat[c(18, 28, 38), j]))
      }), byrow = T)
    } else if (index == "percentage") {
      multi_pc <- sapply(available_list, function(mat) {
        mat[9, j]})
      multi_pf <- sapply(available_list, function(mat) {
        mat[10, j]})
      uni_pc <- sapply(available_list, function(mat) {
        mean(mat[c(19, 29, 39), j])
      })
      uni_pf <- sapply(available_list, function(mat) {
        mean(mat[c(20, 30, 40), j])
      })
      pc_union <- cbind(multi_pc, uni_pc)
      pf_union <- cbind(multi_pf, uni_pf)
      percent[[j]] <- cbind(pc_union, pf_union)
    }
    
   
    if (index != "percentage") {
      value <- c(multi,uni)
      method <- c(rep("multi_RTLP", length(sequence)), rep("multi_harchi", length(sequence)),
                rep("multi_kmedoids", length(sequence)), rep("multi_funHDDC", length(sequence)),
                rep("uni_RTLP", length(sequence)), rep("uni_harchi", length(sequence)),
                rep("uni_kmedoids", length(sequence)), rep("uni_funHDDC", length(sequence)))
      outliers <- rep(paste("Contamination", j, sep = " "), 8 * length(sequence))
      
      if (j == 1) {
        plotmatrix <- data.frame(value, method, outliers)
      } else {
        plotmatrix <- rbind(plotmatrix, data.frame(value, method, outliers))
      }
      } else {
      method <- c(rep("multi_RTLP", length(sequence)), rep("uni_RTLP", length(sequence)))
    }
      }

  if (index == "ARI") {
    formulas <- value ^ 2 ~ method * outliers
    titlename <- bquote("Boxplots of" ~ .(index)^"2" ~ "with Outlier Contaminations")
    yrange <- c(0, 1)
    ylabname <- expression("ARI"^"2")
  } else if (index == "time"){
    formulas <- value ~ method * outliers
    titlename <- "Boxplots of Running Time with Outlier Contaminations"
    ylabname <- "Time (seconds)"
    yrange <- range(plotmatrix$value)
  } else if (index == "percentage") {
    return (percent)
  }
  
  if (index != "percentage") {
    my_names <- c("contamination 1", "contamination 2", "contamination 3",
                "contamination 4", "contamination 5", "contamination 6")
    plotmatrix$method <- factor(plotmatrix$method, levels = 
                                  c("multi_RTLP", "multi_harchi", "multi_kmedoids", "multi_funHDDC",
                                    "uni_RTLP", "uni_harchi", "uni_kmedoids", "uni_funHDDC"))
    
    if (includelegend == T) {
      par(mar = c(5, 3.3, 2.5, 1), mgp = c(2.1, 1, 0))
      } else {
      par(mar = c(2.7, 3.3, 2.5, 1), mgp = c(2.1, 1, 0)) 
      }
    myplot <- boxplot(formulas, data = plotmatrix, col = c(2:8, "pink"),
                    main = bquote(.(titlename) ~ "(Scenario" ~ .(sce)*"," ~ p[size] == .(100*ps) * "% and"
                                  ~ p[curve] == .(100 * pc) * "%)"),
                    varwidth = TRUE, xlim = c(2, 47), ylim = yrange, xaxt = "n", 
                    xlab = NULL, ylab = ylabname, 
                    cex.name = 1, cex.lab = 1.2, cex.main = 1.5)
  
    axis(1, 
        at = c(seq(3, 44, length.out = 6)), 
        labels = my_names, 
        tick = FALSE, cex.axis = 0.9)
  
  # Add the grey vertical lines
    for(i in seq(8.5, 44, by = 8)) { 
      abline(v = i, lty = 1, col = "grey")
    }
    if (includelegend == TRUE) {
      par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
      # Add a legend
      legend("bottom", title = NULL, 
             legend = c("Methods:", "multi_RTLP", "multi_HIER",
                      "multi_KMED", "multi_HDDC",
                      "uni_RTLP", "uni_HIER", 
                      "uni_KMED", "uni_HDDC"), 
            col = c("white", 2:8, "pink"),
            pch = 15, bty = "o",  
            cex = 0.95, text.font = 1, 
            horiz = T, inset = c(0.01, 0.01))
    
   }
  }
    
  if (savefile == TRUE) {
    dev.off()
  }
  return (percent)
}
################ ARI
######################### Figure
load("sim_cluster_outlier_0308_0%.RData")
for (l in 1:6) {
  plotsimulation_meanmarginal(l, 0, "ARI", includelegend = F, savefile = TRUE)
}
for (l in 1:6) {
  plotsimulation_meanmarginal(l, 0, "time", includelegend = F, savefile = FALSE)
}

load("../sim_cluster_outlier_0114_10%.RData")
for (l in 1:6) {
  plotsimulation_meanmarginal(l, 0, "ARI", includelegend = F, savefile = TRUE)
}
for (l in 1:6) {
  plotsimulation_meanmarginal(l, 0, "time", includelegend = F, savefile = FALSE)
}

for (l in 1:6) {
  table_pc0 <- plotsimulation_meanmarginal(l, 0, "percentage", FALSE, FALSE)
  sapply(table_pc0, function(k) {round(100 * apply(k, 2, mean), digits = 1)})
  sapply(table_pc0, function(k) {round(100 * apply(k, 2, sd), digits = 1)})
}
load("../sim_cluster_outlier_0114_30%.RData")
for (l in 1:6) {
  plotsimulation_meanmarginal(l, 0.3, "ARI", includelegend = F, savefile = TRUE)
}
for (l in 1:6) {
  plotsimulation_meanmarginal(l, 0.3, "time", includelegend = F, savefile = FALSE)
}
for (l in 1:6) {
  table_pc30 <- plotsimulation_meanmarginal(l, 0, "percentage", FALSE, FALSE)
  sapply(table_pc30, function(k) {round(100 * apply(k, 2, mean), digits = 1)})
  sapply(table_pc30, function(k) {round(100 * apply(k, 2, sd), digits = 1)})
}

load("../sim_cluster_outlier_0114_60%.RData")
for (l in 1:6) {
  plotsimulation_meanmarginal(l, 0.6, "ARI", includelegend = T, savefile = TRUE)
}

for (l in 1:6) {
  plotsimulation_meanmarginal(l, 0.6, "time", includelegend = F, savefile = FALSE)
}

for (l in 1:6) {
  table_pc60 <- plotsimulation_meanmarginal(l, 0, "percentage", FALSE, FALSE)
  sapply(table_pc60, function(k) {round(100 * apply(k, 2, mean), digits = 1)})
  sapply(table_pc60, function(k) {round(100 * apply(k, 2, sd), digits = 1)})
}
###############################################

