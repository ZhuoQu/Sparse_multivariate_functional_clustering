library("gdata")
source("subsetList.R")
plotsimulation_meanmarginal <- function(sce, index, includelegend){ 
  sequence <- unlist(sapply(1:130, function(k){
    ifelse(length(sim_cluster_outlier[[k]][[sce]]) > 1, return (k), return())}))
  for (j in 1:6) {
    if (index != "time" && index != "percentage") {
      available_list<-lapply(subsetList(sim_cluster_outlier, sequence), function(k) {
        k[[sce]]
      } )
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
      
      mult <- unmatrix(sapply(available_list, function(mat) {
        mat[1:4, j]}), byrow = T)
      uni <- unmatrix(sapply(available_list, function(mat) {
        c(mean(mat[c(11, 21, 31), j]), mean(mat[c(12, 22, 32), j]), 
          mean(mat[c(13, 23, 33), j]), mean(mat[c(14, 24, 34), j]))
      }), byrow = T)
      
      
      ### multi_DBC, uni_DBC, multifunHDDC, unifunHDDC, multihierarchical, unihierarchical, multikmedoids,unikmedoids
    }else if (index == "time") {
      mult <- unmatrix(sapply(available_list, function(mat) {
        mat[c(5, 6, 7, 8), j]
      }), byrow = T)
      uni <- unmatrix(sapply(available_list, function(mat) {
        mat[c(15, 16, 17, 18), j]
      }), byrow = T)
    } 
    
    value <- c(mult,uni)  
    method <- c(rep("multi_DBC", length(sequence)), rep("multi_harchi", length(sequence)),
                rep("multi_kmedoids", length(sequence)), rep("multi_funHDDC", length(sequence)),
                rep("uni_DBC", length(sequence)), rep("uni_harchi", length(sequence)),
                rep("uni_kmedoids", length(sequence)), rep("uni_funHDDC", length(sequence)))
    outliers <- rep(paste("Contamination", j, sep = " "), 8 * length(sequence))
    if (j == 1) {
      plotmatrix <- data.frame(value, method, outliers)
    } else {
      plotmatrix <- rbind(plotmatrix, data.frame(value, method, outliers))
    }
  }
  
  if (index != "time" && index != "percentage") {
    titlename <- paste("Boxplots of ", index, " with Outlier Contaminations", sep = "")
    yrange <- c(0, 1)
    ylabname <- index
    ### multi_DBC, uni_DBC, multifunHDDC, unifunHDDC, multihierarchical, unihierarchical, multikmedoids,unikmedoids
  } else if (index == "time"){
    titlename <- "Boxplots of Running Time with Outlier Contaminations"
    ylabname <- "Time (seconds)"
    yrange <- range(plotmatrix$value)
  }
  my_names <- c("contamination 1", "contamination 2", "contamination 3",
                "contamination 4", "contamination 5", "contamination 6")
  
  plotmatrix$method <- factor(plotmatrix$method, levels = 
                                c("multi_DBC", "multi_harchi", "multi_kmedoids", "multi_funHDDC",
                                  "uni_DBC", "uni_harchi", "uni_kmedoids", "uni_funHDDC"))
  if (includelegend == T) {
    par(mar = c(5, 3.3, 2.5, 1), mgp = c(2.2,1,0))
  } else {
    par(mar = c(2.7,3.3,2.5,1), mgp = c(2.2,1,0)) 
  }
  myplot <- boxplot(value~method * outliers, data = plotmatrix, col = c(2:8, "pink"),
                    main = paste(titlename, " (Scenario ", sce, ")", sep = ""),
                    varwidth = TRUE, xlim = c(2, 47), ylim = yrange, xaxt = "n", 
                    xlab = NULL, ylab = ylabname, cex.names = 0.6, cex.lab = 0.9, cex = 0.8)
  
  axis(1, 
       at = c(seq(3, 44, length.out=6)), 
       labels = my_names , 
       tick=FALSE , cex.axis=0.8)
  
  # Add the grey vertical lines
  for(i in seq(8.5, 44, by = 8)){ 
    abline(v = i, lty = 1, col = "grey")
  }
  if (includelegend == TRUE) {
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    # Add a legend
    legend("bottom", title = "Methods",legend = c("multi_DBC","multi_hierarchical",
                                                  "multi_kmedoids", "multi_funHDDC",
                                                  "uni_DBC","uni_hierarchical", 
                                                  "uni_kmedoids", "uni_funHDDC"), 
           col = c(2:8, "pink"),
           pch = 15, bty = "o",  cex = 0.65, text.font=0.9, horiz = T, inset = c(0.01, 0.01))
    
  }
}
################ ARI
pdf(paste("simulation_outlier_pc60_meanmarginal.pdf",sep=""),width=9.25,height = 3.5)
par(mfrow=c(1,1))
for (sce in 1:6) {
  plotsimulation_meanmarginal(sce, "ARI", includelegend=T)
}
dev.off()
#### Running time
plotsimulation_meanmarginal(1, "ARI", includelegend = T)
plotsimulation_meanmarginal(2, "ARI", includelegend = T)
plotsimulation_meanmarginal(3, "ARI", includelegend = T)
plotsimulation_meanmarginal(4, "ARI", includelegend = T)
plotsimulation_meanmarginal(5, "ARI", includelegend = T)
plotsimulation_meanmarginal(6, "ARI", includelegend = T)

apply(sim_cluster_outlier[[1]][[1]][5:6, ], 1, mean)
apply(sim_cluster_outlier[[1]][[1]][17:18, ], 1 mean)
