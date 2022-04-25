packages <- c("mapdata", "ggplot2", "readr","dplyr", "maps", 
             "rnaturalearth", "rnaturalearthdata", "rgeos",
             "ggspatial", "plotrix")

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

source("00_subsetlist.R")
source("01_algorithm_distance.R")
source("02_algorithm_clustering.R")
source("03_algorithm_find_outlier.R")
source("find_optimal_theta.R")
#source("final_cyclone_00.R")

###################################
load("clust_NP.RData") ## RData should be under the application directory
outlier_index <- clust_NP$isolate

#############################################
worldmap <- map_data ("world", wrap = c(0, 360))
############################################# figures of clustering result
lapply(1:2, function(clustl){
  loc <- ifelse(clustl == 1, "(b)", "(d)")
  index <- unlist(lapply(subsetlist(clust_NP$clust, clustl), 
                         function(l) {l$element}))
  neighbour_number <- apply(neighbours_NP[index, index], 2, sum)
  median <- index[which(neighbour_number == max(neighbour_number))][1]
  firstquarter <- index[order(neighbour_number, decreasing = T)[1:round(0.25 * (length(neighbour_number)))]]
  secondquarter <- index[order(neighbour_number, decreasing = T)[round(0.25 * (length(neighbour_number))+1):round(0.5 * (length(neighbour_number)))] ]
  thirdquarter <- index[order(neighbour_number, decreasing = T)[round(0.5 * (length(neighbour_number)) + 1):round(0.75 * (length(neighbour_number)))] ]
# shift coordinates to recenter worldmap
  pdf(file = paste("clust_", clustl, ".pdf", sep = ""), width = 6, height = 4)
  ggplot(aes(x = long, y = lat), data = worldmap) + 
  geom_path(aes(group = group), 
            #fill = "#f9f9f9", 
            colour = "grey65") + 
    scale_y_continuous(limits = c(0, 70)) +
    scale_x_continuous(limits = c(235, 350)) +
  coord_equal() +  theme_bw() +
  geom_path(aes(x = cal_lon, y = lat,group = id, alpha = 1 - time),
            data = subset(whole_track_df, id %in% westpacific[thirdquarter]), color = "pink") +
  geom_path(aes(x = cal_lon, y = lat, group = id, alpha = 1 - time),
            data = subset(whole_track_df, id %in% westpacific[secondquarter]), color = "magenta") +
  geom_path(aes(x = cal_lon, y = lat,group = id, alpha = 1 - time),
            data = subset(whole_track_df, id %in% westpacific[firstquarter]), color = "purple") +
  geom_path(aes(x = cal_lon, y = lat, group = id, alpha = 1 - time), 
            data = subset(whole_track_df, id %in% westpacific[outlier_index]), color = "red") +
  geom_path(aes(x = cal_lon, y = lat,group = id, alpha = 1 - time),
            data = subset(whole_track_df, id %in% westpacific[median], col = grey(time)), color = "black") + 
  
  #guides(fill=guide_legend(title="Time"))+
  theme(legend.position = "none") +
  
  labs(title = paste(loc," West Pacific Cyclone Tracks: Cluster ", clustl, sep = ""),
       x = "Longitude", y = "Latitude") +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 13),
        plot.margin = margin(1, 1.2, 1, 1.2))
  #scale_colour_gradient(low="red",high="yellow",na.value=NA)
  dev.off()
})
############################################# storm speed polar plot
#####################################
for (clustl in 1:2) {
  loc <- ifelse(clustl == 1, "(c)", "(e)")
  pdf(paste("polar_plot_", clustl,".pdf", sep = ""), width = 4, height = 4)
  par(mfrow = c(1, 1), mgp = c(0.9, 0, 0))
  index <- unlist(lapply(subsetlist(clust_NP$clust, clustl), 
                         function(l) {l$element}))
  neighbour_number <- apply(neighbours_NP[index, index], 2, sum)
  median <- index[which(neighbour_number == max(neighbour_number))][1]
  firstquarter <- index[order(neighbour_number, decreasing = T)[1:round(0.25 * (length(neighbour_number)))]]
  secondquarter <- index[order(neighbour_number, decreasing = T)[round(0.25 * (length(neighbour_number))+1):round(0.5 * (length(neighbour_number)))] ]
  thirdquarter <- index[order(neighbour_number, decreasing = T)[round(0.5 * (length(neighbour_number)) + 1):round(0.75 * (length(neighbour_number)))] ]
  outlier_index <- clust_NP$isolate
  
  plotrix::polar.plot(whole_track_df$storm_speed[1], 
                      polar.pos = whole_track_df$storm_speed_dir[1],
                      label.prop = 1.1, xlab = "Storm Direction (degrees)",
                      ylab = "Storm Speed (knots)",
                      radial.lim = c(0, 100), clockwise = TRUE,
                      #show.grid.labels = seq(0, 100, length.out = 3), 
                      mar = c(2.5, 1, 3, 0), 
                      show.radial.grid = T, explab = FALSE,
                      radial.labels = seq(0, 100, by = 20), grid.left = TRUE,
                      start = 90,
                      rp.type = "p", line.col = 1, 
                      lwd = 2 ,
                      main = paste(loc," Storm Speed Polar Plot: Cluster ", clustl, sep = ""),
                      cex.main = 1.15, cex.lab = 1.15)
  
  # plotrix::polar.plot(whole_track_df$storm_speed[1], 
  #                     label.prop = 1.1, 
  #                     xlab = "Storm Direction (degrees)",
  #                     ylab = "Storm Speed (knots)",
  #                     radial.lim = c(0, 100), clockwise = TRUE,
  #                     radial.pos=seq(0, 20 * 90 / 11.1, length.out = 8),
  #                     label.pos=seq(0, 20 * 90 / 11.1, length.out = 8),
  #                     start = 90, 
  #                     labels = seq(0, 315, by = 45),
  #                     mar = c(2.5, 1, 3, 0), 
  #                     show.radial.grid = T, explab = FALSE,
  #                     grid.left = TRUE,
  #                     rp.type = "p", line.col = 1, 
  #                     lwd = 2,
  #                     main = paste(loc," Storm Speed Polar Plot: Cluster ", clustl, sep = ""),
  #                     cex.main = 1.15, cex.lab = 1.15)
  
  ### pink rgb(255, 192, 203)
  lapply(westpacific[thirdquarter], function(k){
    cat (k, "\n")
    index <- which(whole_track_df$id == k)
    if (length(index) > 10) {
      for (l in 1: max(1,(length(index) - 1))) {
        cat(l,"\n")
        plotrix::polar.plot(whole_track_df$storm_speed[index[l]:(index[l] + 1)], 
                            polar.pos = whole_track_df$storm_speed_dir[index[l]:(index[l] + 1)],
                            rp.type = "p", line.col = rgb(255, 192, 203, alpha = 255 - round(250 * l / (length(index))), maxColorValue = 255),  
                            label.prop = 1.1,
                            radial.lim = c(0, 120), clockwise = TRUE,
                            show.grid.labels = seq(0, 130, length.out = 3),
                            show.radial.grid = T, start = 90,
                            lwd = 1.5,
                            add = TRUE)
      }
    }
  })
  
  ### magenta rgb(255,0,255)
  lapply(westpacific[secondquarter], function(k){
    cat (k, "\n")
    index <- which(whole_track_df$id == k)
    if (length(index) > 10) {
      for (l in 1: max(1,(length(index) - 1))) {
        cat(l,"\n")
        plotrix::polar.plot(whole_track_df$storm_speed[index[l]:(index[l] + 1)], 
                            polar.pos = whole_track_df$storm_speed_dir[index[l]:(index[l] + 1)],
                            rp.type = "p", line.col = rgb(255, 0, 255, alpha = 255 - round(250 * l / (length(index))), maxColorValue = 255),  
                            label.prop = 1.1,
                            radial.lim = c(0, 120), clockwise = TRUE,
                            show.grid.labels = seq(0, 130, length.out = 3),
                            show.radial.grid = T, start = 90,
                            lwd = 1.5,
                            add = TRUE)
      }
    }
  })
  
  ### purple rgb(148,0,211)
  lapply(westpacific[firstquarter], function(k){
    cat (k, "\n")
    index <- which(whole_track_df$id == k)
    if (length(index) > 10) {
      for (l in 1: max(1,(length(index) - 1))) {
        cat(l,"\n")
        plotrix::polar.plot(whole_track_df$storm_speed[index[l]:(index[l] + 1)], 
                            polar.pos = whole_track_df$storm_speed_dir[index[l]:(index[l] + 1)],
                            rp.type = "p", line.col = rgb(148, 0, 211, alpha = 255 - round(250 * l / (length(index))), maxColorValue = 255),  
                            label.prop = 1.1,
                            radial.lim = c(0, 120), clockwise = TRUE,
                            show.grid.labels = seq(0, 130, length.out = 3),
                            show.radial.grid = T, start = 90,
                            lwd = 1.5,
                            add = TRUE)
      }
    }
  })
  
  lapply(westpacific[outlier_index], function(k){
    cat (k, "\n")
    index <- which(whole_track_df$id == k)
    if (length(index) > 10) {
      for (l in 1: max(1,(length(index) - 1))) {
        cat(l,"\n")
        plotrix::polar.plot(whole_track_df$storm_speed[index[l]:(index[l] + 1)], 
                            polar.pos = whole_track_df$storm_speed_dir[index[l]:(index[l] + 1)],
                            rp.type = "p", line.col = rgb(255, 0, 0, 
                                                          alpha = 255 - round(250 * l / (length(index))), 
                                                          maxColorValue = 255),  
                            label.prop = 1.1,
                            radial.lim = c(0, 120), clockwise = TRUE,
                            show.grid.labels = seq(0, 130, length.out = 3),
                            show.radial.grid = T, start = 90,
                            lwd = 1.5,
                            add = TRUE)
      }
    }
  })
  
  index <- which(whole_track_df$id == westpacific[median])
  for (l in 1: (length(index) - 1)) {
    plotrix::polar.plot(whole_track_df$storm_speed[index[l]:(index[l] + 1)], 
                        polar.pos = whole_track_df$storm_speed_dir[index[l]:(index[l] + 1)],
                        rp.type = "p", line.col = rgb(0, 0, 0, alpha = 255 - round(250 * l / (length(index))), maxColorValue = 255), 
                        label.prop = 1.1,
                        radial.lim = c(0, 120), clockwise = TRUE,
                        show.grid.labels = seq(0, 130, length.out = 3),
                        show.radial.grid = T, start = 90,
                        lwd = 1.5,
                        add = TRUE)
  }
  
  dev.off()
}



