library("gdata")
source("00_subsetlist.R")
packages <- c("randomcoloR", "RColorBrewer")

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

color_palette <- brewer.pal(n = 10, name = "Paired")
color_palette[c(6, 1, 7, 2, 8, 3, 9, 4, 10, 5)] <- color_palette

ps <- 1
plotsimulation_meanmarginal <- function(sce, pc, index, includelegend, savefile = TRUE) { 
  if (savefile == TRUE && includelegend == TRUE) {
    pdf(file = paste("sce", sce, "_pc", 100 * pc, ".pdf", sep = ""), width = 11, height = 5.2)
  } else if (savefile == TRUE) {
    pdf(file = paste("sce", sce, "_pc", 100 * pc, ".pdf", sep = ""), width = 11, height = 4.95)
  }
  
  sequence <- unlist(sapply(1:repeattimes, function(k){
    ifelse(length(sim_cluster_outlier[[k]]) > 1, return (k), return())}))
  
  percent <- list()
  for (j in 1:6) {
    available_list <- lapply(subsetlist(sim_cluster_outlier, sequence), function(k) {
      k[[sce]]
    } )
    if (index == "ARI") {
      mult_dbc <- sapply(available_list, function(mat) {
        mat[1, j]})
      mult_dbscan <- sapply(available_list, function(mat) {
        mat[2, j]})
      mult_hierarchical <- sapply(available_list, function(mat) {
        mat[3, j]})
      mult_kmedoids <- sapply(available_list, function(mat) {
        mat[4, j]})
      mult_funhddc <- sapply(available_list, function(mat) {
        mat[5, j]})
      
      uni_dbc <- sapply(available_list, function(mat) {
        mean(mat[c(14 + 1, 14 * 2 + 1, 14 * 3 + 1), j])})
      uni_dbscan <- sapply(available_list, function(mat) {
        mean(mat[c(14 + 2, 14 * 2 + 2, 14 * 3 + 2), j])})
      uni_hierarchical <- sapply(available_list, function(mat) {
        mean(mat[c(14 + 3, 14 * 2 + 3, 14 * 3 + 3), j])})
      uni_kmedoids <- sapply(available_list, function(mat) {
        mean(mat[c(14 + 4, 14 * 2 + 4, 14 * 3 + 4), j])})
      uni_funhddc <- sapply(available_list, function(mat) {
        mean(mat[c(14 + 5, 14 * 2 + 5, 14 * 3 + 5), j])})
      
      multi <- unmatrix(sapply(available_list, function(mat) {
        mat[1:5, j]}), byrow = T)
      uni <- unmatrix(sapply(available_list, function(mat) {
        c(mean(mat[c(14 + 1, 14 * 2 + 1, 14 * 3 + 1), j]), 
          mean(mat[c(14 + 2, 14 * 2 + 2, 14 * 3 + 2), j]), 
          mean(mat[c(14 + 3, 14 * 2 + 3, 14 * 3 + 3), j]), 
          mean(mat[c(14 + 4, 14 * 2 + 4, 14 * 3 + 4), j]),
          mean(mat[c(14 + 5, 14 * 2 + 5, 14 * 3 + 5), j]))
      }), byrow = T)
      
      
      ### multi_DBC, uni_DBC, multifunHDDC, unifunHDDC, multihierarchical, unihierarchical, multikmedoids,unikmedoids
    } else if (index == "time") {
      mult_dbc <- sapply(available_list, function(mat) {
        mat[6, j]})
      mult_dbscan <- sapply(available_list, function(mat) {
        mat[7, j]})
      mult_hierarchical <- sapply(available_list, function(mat) {
        mat[8, j]})
      mult_kmedoids <- sapply(available_list, function(mat) {
        mat[9, j]})
      mult_funhddc <- sapply(available_list, function(mat) {
        mat[10, j]})
      
      uni_dbc <- sapply(available_list, function(mat) {
        mean(mat[c(14 + 6, 14 * 2 + 6, 14 * 3 + 6), j])})
      uni_dbscan <- sapply(available_list, function(mat) {
        mean(mat[c(14 + 7, 14 * 2 + 7, 14 * 3 + 7), j])})
      uni_hierarchical <- sapply(available_list, function(mat) {
        mean(mat[c(14 + 8, 14 * 2 + 8, 14 * 3 + 8), j])})
      uni_kmedoids <- sapply(available_list, function(mat) {
        mean(mat[c(14 + 9, 14 * 2 + 9, 14 * 3 + 9), j])})
      uni_funhddc <- sapply(available_list, function(mat) {
        mean(mat[c(14 + 10, 14 * 2 + 10, 14 * 3 + 10), j])})
      
      multi <- unmatrix(sapply(available_list, function(mat) {
        mat[6:10, j]}), byrow = T)
      uni <- unmatrix(sapply(available_list, function(mat) {
        c(mean(mat[c(14 + 6, 14 * 2 + 6, 14 * 3 + 6), j]), 
          mean(mat[c(14 + 7, 14 * 2 + 7, 14 * 3 + 7), j]), 
          mean(mat[c(14 + 8, 14 * 2 + 8, 14 * 3 + 8), j]), 
          mean(mat[c(14 + 9, 14 * 2 + 9, 14 * 3 + 9), j]),
          mean(mat[c(14 + 10, 14 * 2 + 10, 14 * 3 + 10), j]))
      }), byrow = T)
    } else if (index == "percentage") {
      multi_dbc_pc <- sapply(available_list, function(mat) {
        mat[11, j]})
      multi_dbc_pf <- sapply(available_list, function(mat) {
        mat[12, j]})
      multi_dbscan_pc <- sapply(available_list, function(mat) {
        mat[13, j]})
      multi_dbscan_pf <- sapply(available_list, function(mat) {
        mat[14, j]})
      
      uni_dbc_pc <- sapply(available_list, function(mat) {
        mean(mat[c(14 + 11, 14 * 2 + 11, 14 * 3 + 11), j])
      })
      uni_dbc_pf <- sapply(available_list, function(mat) {
        mean(mat[c(14 + 12, 14 * 2 + 12, 14 * 3 + 12), j])
      })
      
      uni_dbscan_pc <- sapply(available_list, function(mat) {
        mean(mat[c(14 + 13, 14 * 2 + 13, 14 * 3 + 13), j])
      })
      uni_dbscan_pf <- sapply(available_list, function(mat) {
        mean(mat[c(14 + 14, 14 * 2 + 14, 14 * 3 + 14), j])
      })
      
      dbc_pc_union <- cbind(multi_dbc_pc, uni_dbc_pc)
      dbc_pf_union <- cbind(multi_dbc_pf, uni_dbc_pf)
      
      dbscan_pc_union <- cbind(multi_dbscan_pc, uni_dbscan_pc)
      dbscan_pf_union <- cbind(multi_dbscan_pf, uni_dbscan_pf)
      
      percent[[j]] <- cbind(dbc_pc_union, dbc_pf_union, dbscan_pc_union, dbscan_pf_union)
      colnames(percent[[j]]) <- c("multi_pc_dbc", "uni_pc_dbc",
                                  "multi_pf_dbc", "uni_pf_dbc",
                                  "multi_pc_dbscan", "uni_pc_dbscan",
                                  "multi_pf_dbscan", "uni_pf_dbscan")
      rownames(percent[[j]]) <- 1 : length(sequence)
    }
    
   
    if (index != "percentage") {
      value <- c(multi, uni)
      method <- c(rep("multi_RTLP", length(sequence)), rep("multi_dbscan", length(sequence)),
                  rep("multi_harchi", length(sequence)), rep("multi_kmedoids", length(sequence)), 
                  rep("multi_funHDDC", length(sequence)),
                rep("uni_RTLP", length(sequence)), rep("uni_dbscan", length(sequence)), 
                rep("uni_harchi", length(sequence)), rep("uni_kmedoids", length(sequence)), 
                rep("uni_funHDDC", length(sequence)))
      outliers <- rep(paste("Contamination", j, sep = " "), 10 * length(sequence))
      
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
                                  c("multi_RTLP", "multi_dbscan", "multi_harchi", "multi_kmedoids", "multi_funHDDC",
                                    "uni_RTLP", "uni_dbscan", "uni_harchi", "uni_kmedoids", "uni_funHDDC"))
    
    if (includelegend == T) {
      par(mar = c(5, 3.3, 2.5, 1), mgp = c(2.1, 1, 0))
      } else {
      par(mar = c(2.7, 3.3, 2.5, 1), mgp = c(2.1, 1, 0)) 
      }
    myplot <- boxplot(formulas, data = plotmatrix, col = color_palette,
                    main = bquote(.(titlename) ~ "(Scenario" ~ .(sce)*"," ~ p[size] == .(100 * ps) * "% and"
                                  ~ p[curve] == .(100 * pc) * "%)"),
                    varwidth = TRUE, xlim = c(2, 59), ylim = yrange, xaxt = "n", 
                    xlab = NULL, ylab = ylabname, 
                    cex.name = 1, cex.lab = 1.2, cex.main = 1.5)
  
    axis(1, 
        at = c(seq(3, 58, length.out = 6)), 
        labels = my_names, 
        tick = FALSE, cex.axis = 0.9)
  
  # Add the grey vertical lines
    for(i in seq(10.5, 50.5, length.out = 5)) { 
      abline(v = i, lty = 1, col = "grey")
    }
    if (includelegend == TRUE) {
      par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
      # Add a legend
      legend("bottom", title = NULL, 
             legend = c("Methods:", "multi_RTLP", 
                        "multi_DBSCAN", "multi_HIER",
                      "multi_KMED", "multi_HDDC",
                      "uni_RTLP", "uni_DBSCAN", "uni_HIER", 
                      "uni_KMED", "uni_HDDC"), 
            col = c("white", color_palette),
            pch = 15, bty = "o",  
            text.font = 1, cex = 0.65, 
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
pc <- 0
load(paste("sim_cluster_outlier_0924_", 100 * pc,"%.RData", sep = ""))
for (l in 1:6) {
  plotsimulation_meanmarginal(sce = l, pc, index = "ARI", 
                              includelegend = FALSE, savefile = TRUE)
}

pc <- 0.3
load(paste("sim_cluster_outlier_0924_", 100 * pc,"%.RData", sep = ""))
for (l in 1:6) {
  plotsimulation_meanmarginal(sce = l, pc, index = "ARI", 
                              includelegend = FALSE, savefile = TRUE)
}

pc <- 0.6
load(paste("sim_cluster_outlier_0924_", 100 * pc,"%.RData", sep = ""))
for (l in 1:6) {
  plotsimulation_meanmarginal(sce = l, pc, index = "ARI", 
                              includelegend = TRUE, savefile = TRUE)
}

for (l in 1:6) {
  plotsimulation_meanmarginal(sce = l, pc, index = "ARI", 
                              includelegend = FALSE, savefile = FALSE)
}

pc <- 0.3
load(paste("sim_cluster_outlier_0924_", 100 * pc,"%.RData", sep = ""))
for (l in 1:6) {
  table_outlier <- plotsimulation_meanmarginal(l, pc, "percentage", FALSE, FALSE)
  sapply(table_outlier, function(k) {round(100 * apply(k, 2, mean), digits = 1)})
  sapply(table_outlier, function(k) {round(100 * apply(k, 2, sd), digits = 1)})
}

###############################################

