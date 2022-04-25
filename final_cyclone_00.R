
packages = c("ggmap",
             "ggplot2", "readr","dplyr")

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
source("clust_theta.R")

whole_track_df <- read_csv('application/whole_track_cyclone.csv')
whole_track_df_copy <- whole_track_df
whole_track_df_copy$vis_lon[intersect(which(whole_track_df$vis_lon < 6),
                                      which(whole_track_df$vis_lon >= 0))] <- 0
# myLocation <- c(-180, -65, 179.99, 65)
# 
# myMap <- get_map(location = myLocation,
#                  source = "stamen", 
#                  maptype = "watercolor", 
#                  crop = FALSE)
# 



# cyclone_calculation<- lapply(1:length(unique(whole_track_df$id)),
#                              function(i){
#                                index <- which(whole_track_df$id == unique(whole_track_df$id)[i])
#                                len_time <- length(which(whole_track_df$id == unique(whole_track_df$id)[i]))
#                               
#                                return(list(
#                                  argvals = whole_track_df$time[index], 
#                                  subj = rep(unique(whole_track_df$id)[i], len_time),
#                                  y = cbind(whole_track_df$cal_lon[index], whole_track_df$lat[index])
#                                ))
#                              })
#distance_global <- distance_mfd(cyclone_calculation)

#save(distance_global, file = "distance_gb.Rdata")
load("application/distance_gb.Rdata")

#theta_global <- find_optimal_theta(distance_global, NULL, 0.05, 0.75, "silhouette")
  
#clust_global<-clustresult(distance_global, theta)

neighbours_global <- findneighbours(distance_global, theta = 0.05)
cluster_global <- clustering_warping(neighbours_global)  
clust_global <- include_isolatepoints(cluster_global[[2]], minprop = 0.03, 
                                     distance_global, alpha_c = 0.85)

###############
###################
#### west Pacific
clustl <- c(1, 8)
westpacific <- unlist(lapply(subsetlist(clust_global$clust, clustl), function(l){l$element}))
west_pacific_subset <- subset(whole_track_df, id %in% westpacific)
# ggmap(myMap) +
#   #geom_point(aes(x = lon, y = lat),
#   #           data = rbind(SI_track_df,NI_track_df), alpha = .5, color="black") +
#   
#   geom_path(aes(x = vis_lon, y = lat, group = name, colour = time),
#             data = west_pacific_subset, alpha = .5) +
#   
#   labs(title = "Global Cyclone Tracks (1842-2021)", x = "Longitude", y = "Latitude") +
#   theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(0.8, 1.2, 0.6, 1.2)) +
#   scale_colour_gradient(low = "red", high = "yellow", na.value = NA)
# 
# 
# 

################################ clustering in west pacific cyclone tracks
distance_westpacific <- distance_global[westpacific, westpacific]
#theta <- find_optimal_theta(distance_westpacific, NULL, 0.05,  0.75, "silh")
theta_candidate <- seq(0.01, 0.30, by = 0.01)
clust_dataframe <- t(sapply(theta_candidate, function(theta) {
  clust_result <- clust_theta(distance_westpacific, NULL, 0.05, theta, 0.87)
}))

theta <- theta_candidate[clust_dataframe[, 2] == max(clust_dataframe[, 2])]
neighbours_westpacific <- findneighbours(distance_westpacific, theta)
cluster_westpacific <- clustering_warping(neighbours_westpacific)  
clust_westpacific <- include_isolatepoints(cluster_westpacific[[2]], minprop = 0.05,
                                           distance_westpacific, alpha_c = 0.87)
directory <- getwd()
setwd(paste(directory, "/application", sep = ""))
save(clust_westpacific, file = "clust_westpacific.RData")


###############################
worldmap <- map_data ("world", wrap = c(0, 360))
whole_track_df$Time <- whole_track_df$time
pdf("wp.pdf", width = 6.5, height = 4)
ggplot(aes(x = long, y = lat), data = worldmap) + 
  geom_path(aes(group = group), 
            #fill = "#f9f9f9", 
            colour = "grey65") + 
  scale_y_continuous(limits = c(0, 65)) +
  scale_x_continuous(limits = c(78, 198)) +
  coord_equal() +  theme_bw() +
  geom_path(aes(x = cal_lon, y = lat, group = id, colour = Time),
            data = west_pacific_subset) + 
  labs(fill = "Time") +
  #guides(fill = guide_legend(title = "Time")) + 
  #theme(legend.position = "none") +
  
  labs(title = expression(bold("(a) West Pacific Cyclone Tracks")), x = "Longitude",y = "Latitude") +
  theme(plot.title = element_text(hjust = 0.5, size = 15), 
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 13),
        plot.margin = margin(1, 1.2, 0.5, 1)) + 
  scale_colour_gradient(low = "red",high = "yellow", na.value = NA,
                        breaks = seq(0, 1, by = 0.2), labels = seq(0.0, 1.0, by = 0.2),
                        limits = c(0, 1))
dev.off()


# 
# ############# South Pacific
# clustl = 5
# ggmap(myMap) +
#   #geom_point(aes(x = lon, y = lat),
#   #           data = rbind(SI_track_df,NI_track_df), alpha = .5, color="black") +
#   
#   geom_path(aes(x = vis_lon, y = lat,group=name,colour=time),data = subset(cyclone_subset, id %in% unlist(lapply(subsetlist(clust_global$clust,clustl), function(l){l$element}))), alpha = .5) +
#   
#   #guides(fill=guide_legend(title="Time"))+
#   #theme(legend.position = "none")+
#   
#   labs(title = "Global Cyclone Tracks (1842-2021)",x = "Longitude",y = "Latitude") +
#   theme(plot.title = element_text(hjust = 0.5),plot.margin = margin(0.8, 1.2, 0.6, 1.2))+
#   scale_colour_gradient(low="red",high="yellow",na.value=NA)
# 
# southpacific <- unlist(lapply(subsetlist(clust_global$clust,clustl), function(l){l$element}))
# distance_southpacific<- distance_global[southpacific,southpacific]
# neighbours_southpacific<-findneighbours(distance_southpacific,theta=0.05)
# cluster_southpacific<-clustering_warping(neighbours_southpacific)  
# clust_southpacific<-includeisolatepoints(cluster_southpacific,minprop=140/nrow(distance_southpacific),distance_southpacific,thres=1,outliers=T)
# ######
# clustl = 1
# index <- southpacific[unlist(lapply(subsetlist(clust_southpacific$clust,clustl), function(l){l$element}))]
# ggmap(myMap) +
#   #geom_point(aes(x = lon, y = lat),
#   #           data = rbind(SI_track_df,NI_track_df), alpha = .5, color="black") +
#   
#   geom_path(aes(x = vis_lon, y = lat,group=name,colour=time),data = subset(cyclone_subset, id %in% unique(cyclone_subset$id)[index]), alpha = .5) +
#   
#   #guides(fill=guide_legend(title="Time"))+
#   #theme(legend.position = "none")+
#   
#   labs(title = "Global Cyclone Tracks (1842-2021)",x = "Longitude",y = "Latitude") +
#   theme(plot.title = element_text(hjust = 0.5),plot.margin = margin(0.8, 1.2, 0.6, 1.2))+
#   scale_colour_gradient(low="red",high="yellow",na.value=NA)
# 
# ########### Indian ocean
# clustl = c(2,3,7)
# ggmap(myMap) +
#   #geom_point(aes(x = lon, y = lat),
#   #           data = rbind(SI_track_df,NI_track_df), alpha = .5, color="black") +
#   
#   geom_path(aes(x = vis_lon, y = lat,group=name,colour=time),data = subset(cyclone_subset, id %in% unique(cyclone_subset$id)[unlist(lapply(subsetlist(clust_global$clust,clustl), function(l){l$element}))]), alpha = .5) +
#   
#   #guides(fill=guide_legend(title="Time"))+
#   #theme(legend.position = "none")+
#   
#   labs(title = "Global Cyclone Tracks (1842-2021)",x = "Longitude",y = "Latitude") +
#   theme(plot.title = element_text(hjust = 0.5),plot.margin = margin(0.8, 1.2, 0.6, 1.2))+
#   scale_colour_gradient(low="red",high="yellow",na.value=NA)
# ###
# indianocean <- unlist(lapply(subsetlist(clust_global$clust,clustl), function(l){l$element}))
# distance_indianocean<- distance_global[indianocean,indianocean]
# neighbours_indianocean<-findneighbours(distance_indianocean,theta=0.05)
# cluster_indianocean<-clustering_warping(neighbours_indianocean)  
# clust_indianocean<-includeisolatepoints(cluster_indianocean,minprop=209/nrow(distance_indianocean),distance_indianocean,thres=1,outliers=T)
# ###
# clustl = 5
# index <- indianocean[unlist(lapply(subsetlist(clust_indianocean$clust,clustl), function(l){l$element}))]
# ggmap(myMap) +
#   #geom_point(aes(x = lon, y = lat),
#   #           data = rbind(SI_track_df,NI_track_df), alpha = .5, color="black") +
#   
#   geom_path(aes(x = vis_lon, y = lat,group=name,colour=time),data = subset(cyclone_subset, id %in% unique(cyclone_subset$id)[index]), alpha = .5) +
#   
#   #guides(fill=guide_legend(title="Time"))+
#   #theme(legend.position = "none")+
#   
#   labs(title = "Global Cyclone Tracks (1842-2021)",x = "Longitude",y = "Latitude") +
#   theme(plot.title = element_text(hjust = 0.5),plot.margin = margin(0.8, 1.2, 0.6, 1.2))+
#   scale_colour_gradient(low="red",high="yellow",na.value=NA)
# 
# #### North America
# clustl = c(4,6,9)
# ggmap(myMap) +
#   #geom_point(aes(x = lon, y = lat),
#   #           data = rbind(SI_track_df,NI_track_df), alpha = .5, color="black") +
#   
#   geom_path(aes(x = vis_lon, y = lat,group=name,colour=time),data = subset(cyclone_subset, id %in% unique(cyclone_subset$id)[unlist(lapply(subsetlist(clust_global$clust,clustl), function(l){l$element}))]), alpha = .5) +
#   
#   #guides(fill=guide_legend(title="Time"))+
#   #theme(legend.position = "none")+
#   
#   labs(title = "Global Cyclone Tracks (1842-2021)",x = "Longitude",y = "Latitude") +
#   theme(plot.title = element_text(hjust = 0.5),plot.margin = margin(0.8, 1.2, 0.6, 1.2))+
#   scale_colour_gradient(low="red",high="yellow",na.value=NA)
# ###
# ###
# atlanticocean <- unlist(lapply(subsetlist(clust_global$clust,clustl), function(l){l$element}))
# distance_atlanticocean<- distance_global[atlanticocean,atlanticocean]
# neighbours_atlanticocean<-findneighbours(distance_atlanticocean,theta=0.05)
# cluster_atlanticocean<-clustering_warping(neighbours_atlanticocean)  
# clust_atlanticocean<-includeisolatepoints(cluster_atlanticocean,minprop=209/nrow(distance_atlanticocean),distance_atlanticocean,thres=1,outliers=T)
# 
# clustl = 1
# index <- atlanticocean[unlist(lapply(subsetlist(cluster_atlanticocean,clustl), function(l){l$element}))]
# ggmap(myMap) +
#   #geom_point(aes(x = lon, y = lat),
#   #           data = rbind(SI_track_df,NI_track_df), alpha = .5, color="black") +
#   
#   geom_path(aes(x = vis_lon, y = lat,group=name,colour=time),data = subset(cyclone_subset, id %in% unique(cyclone_subset$id)[index]), alpha = .5) +
#   
#   #guides(fill=guide_legend(title="Time"))+
#   #theme(legend.position = "none")+
#   
#   labs(title = "Global Cyclone Tracks (1842-2021)",x = "Longitude",y = "Latitude") +
#   theme(plot.title = element_text(hjust = 0.5),plot.margin = margin(0.8, 1.2, 0.6, 1.2))+
#   scale_colour_gradient(low="red",high="yellow",na.value=NA)
# 
# #### trying to find east pacific ocean
# elem<- c(setdiff(1:14,c(3,5,12)), 15:672, 973:1038,1080:1141)
# core_elem <- setdiff(1:14,c(3,5,12))
# cluster_eastpacific <- list(list(center=atlanticocean[cluster_atlanticocean[[1]]$center],core = atlanticocean[cluster_atlanticocean[[1]]$core[core_elem]],element = atlanticocean[cluster_atlanticocean[[1]]$element[elem]]),
#                     list(center = atlanticocean[cluster_atlanticocean[[3]]$center],
#                          core = atlanticocean[cluster_atlanticocean[[3]]$core],
#                          element = atlanticocean[cluster_atlanticocean[[3]]$element]),
#                     list(center = atlanticocean[cluster_atlanticocean[[9]]$center],
#                          core = atlanticocean[cluster_atlanticocean[[9]]$core],
#                          element = atlanticocean[cluster_atlanticocean[[9]]$element]),
#                     list(center = atlanticocean[cluster_atlanticocean[[14]]$center],
#                          core = atlanticocean[cluster_atlanticocean[[14]]$core],
#                          element = atlanticocean[cluster_atlanticocean[[14]]$element])
# )
# 
# eastpacific <- unlist(lapply(cluster_eastpacific, function(k){k$element}))
# distance_eastpacific<- distance_global[eastpacific,eastpacific]
# neighbours_eastpacific<-findneighbours(distance_eastpacific,theta=0.05)
# cluster_eastpacific<-clustering_warping(neighbours_eastpacific)  
# clust_eastpacific<-includeisolatepoints(cluster_eastpacific,minprop=100/nrow(distance_eastpacific),distance_eastpacific,thres=1,outliers=T)
# 
# #################################
# clustl = 1
# index <- eastpacific[unlist(lapply(subsetlist(clust_eastpacific$clust,clustl), function(l){l$element}))]
# index <-eastpacific
# ggmap(myMap) +
#   #geom_point(aes(x = lon, y = lat),
#   #           data = rbind(SI_track_df,NI_track_df), alpha = .5, color="black") +
#   
#   geom_path(aes(x = vis_lon, y = lat,group=name,colour=time),data = subset(cyclone_subset, id %in% unique(cyclone_subset$id)[index]), alpha = .5) +
#   
#   #guides(fill=guide_legend(title="Time"))+
#   #theme(legend.position = "none")+
#   
#   labs(title = "Global Cyclone Tracks (1842-2021)",x = "Longitude",y = "Latitude") +
#   theme(plot.title = element_text(hjust = 0.5),plot.margin = margin(0.8, 1.2, 0.6, 1.2))+
#   scale_colour_gradient(low="red",high="yellow",na.value=NA)
# 
# ############## update the atlantic ocean
# 
# atlanticocean <- setdiff(atlanticocean,eastpacific)
# distance_atlanticocean<- distance_global[atlanticocean,atlanticocean]
# neighbours_atlanticocean<-findneighbours(distance_atlanticocean,theta=0.05)
# cluster_atlanticocean<-clustering_warping(neighbours_atlanticocean)  
# clust_atlanticocean<-includeisolatepoints(cluster_atlanticocean,minprop=150/nrow(distance_atlanticocean),distance_atlanticocean,thres=1,outliers=T)
# 
# clustl<-3
# index <- atlanticocean[unlist(lapply(subsetlist(clust_atlanticocean$clust,clustl), function(l){l$element}))]
# ggmap(myMap) +
#   #geom_point(aes(x = lon, y = lat),
#   #           data = rbind(SI_track_df,NI_track_df), alpha = .5, color="black") +
#   
#   geom_path(aes(x = vis_lon, y = lat,group=name,colour=time),data = subset(cyclone_subset, id %in% unique(cyclone_subset$id)[index]), alpha = .5) +
#   
#   #guides(fill=guide_legend(title="Time"))+
#   #theme(legend.position = "none")+
#   
#   labs(title = "Global Cyclone Tracks (1842-2021)",x = "Longitude",y = "Latitude") +
#   theme(plot.title = element_text(hjust = 0.5),plot.margin = margin(0.8, 1.2, 0.6, 1.2))+
#   scale_colour_gradient(low="red",high="yellow",na.value=NA)
# 
# 
# ######################################
# isolatel = 1:50
# ggmap(myMap) +
#   #geom_point(aes(x = lon, y = lat),
#   #           data = rbind(SI_track_df,NI_track_df), alpha = .5, color="black") +
#   
#   geom_path(aes(x = vis_lon, y = lat,group=name,colour=time),data = subset(cyclone_subset, id %in% unique(cyclone_subset$id)[unlist(lapply(subsetlist(clust_global$isolate, isolatel),function (k){k$element}))]), alpha = .5) +
#   
#   #guides(fill=guide_legend(title="Time"))+
#   #theme(legend.position = "none")+
#   
#   labs(title = "Global Cyclone Tracks (1842-2021)",x = "Longitude",y = "Latitude") +
#   theme(plot.title = element_text(hjust = 0.5),plot.margin = margin(0.8, 1.2, 0.6, 1.2))+
#   scale_colour_gradient(low="red",high="yellow",na.value=NA)
# 


### Try to run the code with SASFunclust method but it shows

##### Error in bsplineS(evalarg, breaks, norder, nderiv) : 
###  Knots do not span the values of X
library(sasfunclust)
timegrid <- unique(west_pacific_subset$time[order(west_pacific_subset$time, 
                                                  decreasing = F)])
timeindex <- sapply(west_pacific_subset$time, function(k) {
  which(timegrid == k)
})
mod <- sasfclust(X = west_pacific_subset$lon,
                 grid = timeindex,
                 curve = factor(west_pacific_subset$SID, labels = 1:length(westpacific)),
                 timeindex = timegrid,
                 lambda_s = 10^-6,
                 lambda_l = 5, 
                 G = 2,
                 maxit = 5, q = 10)
