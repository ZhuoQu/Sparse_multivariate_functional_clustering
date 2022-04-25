library(readr)
library(ggplot2)  # ggplot() fortify()
library(dplyr)  # %>% select() filter() bind_rows()
library(rgdal)  # readOGR() spTransform()
#library(raster)  # intersect()
#library(ggsn)  # north2() scalebar()
#library(rworldmap)  # getMap()
library(ggmap)
library(rlist)
library(sjmisc)
library(lubridate)
library(gdata)


maxColorValue = 100
palette <- colorRampPalette(c("red","yellow"))(maxColorValue)
ibtracs <- read_csv("ibtracs.ALL.list.v04r00.csv")
ibtracs<-ibtracs[-1,]
unique(ibtracs$NATURE)
unique(ibtracs$BASIN)


ibtracs$LAT<-as.numeric(ibtracs$LAT)
ibtracs$LON<-as.numeric(ibtracs$LON)

ibtracs$calculate_LON<-ifelse(ibtracs$LON<0,ibtracs$LON+360,ibtracs$LON)
ibtracs$visualize_LON<-ifelse(ibtracs$LON>=180,ibtracs$LON-360,ibtracs$LON)
range(ibtracs$LON)
range(ibtracs$calculate_LON)
range(ibtracs$visualize_LON)
plot(ibtracs$calculate_LON,ibtracs$visualize_LON)
#######################################
range(ibtracs$LAT)
range(ibtracs$LON)
whole_track<-list()
wholecyclone_SID <- unique(ibtracs$SID)
wholecyclone_length <- c()
wholecyclone_month <- list()
ID = 0
for (i in wholecyclone_SID)
{
  index<-which(ibtracs$SID==i)
    wholecyclone_length<-list.append(wholecyclone_length,length(index))
    wholecyclone_month<-list.append(wholecyclone_month,unique(month(ibtracs$ISO_TIME[index])))
    wholecyclone_year<-ibtracs$SEASON[index]
    lat=ibtracs$LAT[index]
    lon=ibtracs$LON[index]
    calculate_lon=ibtracs$calculate_LON[index]
    visualize_lon=ibtracs$visualize_LON[index]    
    sd_time=1:length(lat)/length(lat)
    name1=rep(i,length(index))
    name2=ibtracs$NAME[index]
    name3<-ifelse(visualize_lon>0,"pos","neg")
    basin=ibtracs$BASIN[index]
    nature = ibtracs$NATURE[index]
    ID = ID + 1
    #########################
    year=as.numeric(ibtracs$SEASON[index])
    month=month(ibtracs$ISO_TIME)[index]
    duration=ibtracs$ISO_TIME[index[length(index)]]-ibtracs$ISO_TIME[index[1]]
    track_type=ibtracs$TRACK_TYPE[index]
    storm_speed=as.numeric(ibtracs$STORM_SPEED[index])
    storm_speed_dir=as.numeric(ibtracs$STORM_DIR[index])
    wind_speed=as.numeric(ibtracs$WMO_WIND[index])
    wind_pressure=as.numeric(ibtracs$WMO_PRES[index])
    ###########################
    whole_track = list.append(whole_track,list(lat=lat,lon=lon,cal_lon=calculate_lon,vis_lon=visualize_lon,time=sd_time,
                                       SID=rep(i,length(index)),
                                       name= paste(name1,name2,name3),
                                       nature = nature,
                                       basin = basin,
                                       id = rep(ID,length(index)),
                                       year = year,
                                       month = month,
                                       duration = rep(duration,length(index)),
                                       track_type = track_type,
                                       storm_speed = storm_speed,
                                       storm_speed_dir = storm_speed_dir,
                                       wind_speed = wind_speed,
                                       wind_pressure = wind_pressure)
                            )
    
}


whole_track_df <- data.frame(
                        lon=unlist(lapply(whole_track, function(k){k$lon})),
                        cal_lon=unlist(lapply(whole_track, function(k){k$cal_lon})),
                        vis_lon=unlist(lapply(whole_track, function(k){k$vis_lon})),
                        lat=unlist(lapply(whole_track, function(k){k$lat})),
                        time=unlist(lapply(whole_track, function(k){k$time})),
                        name=unlist(lapply(whole_track, function(k){k$name})),
                        nature=unlist(lapply(whole_track, function(k){k$nature})),
                        basin=unlist(lapply(whole_track, function(k){k$basin})),
                        SID=unlist(lapply(whole_track, function(k){k$SID})),
                        id = unlist(lapply(whole_track, function(k){k$id})),
                        year = unlist(lapply(whole_track, function(k){k$year})),
                        month = unlist(lapply(whole_track, function(k){k$month})),
                        duration = unlist(lapply(whole_track, function(k){k$duration})),
                        track_type = unlist(lapply(whole_track, function(k){k$track_type})),
                        storm_speed = unlist(lapply(whole_track, function(k){k$storm_speed})),
                        storm_speed_dir = unlist(lapply(whole_track, function(k){k$storm_speed_dir})),
                        wind_speed = unlist(lapply(whole_track, function(k){k$wind_speed})),
                        wind_pressure = unlist(lapply(whole_track, function(k){k$wind_pressure}))
                        
                        
)



############################################
myLocation <- c(-180, -65, 179.99, 65)

myMap <- get_map(location=myLocation,
                 source="stamen", maptype="watercolor", crop=FALSE)

#### If we see all the whole data together
 ggmap(myMap) +
  #geom_point(aes(x = lon, y = lat),
  #           data = rbind(SI_track_df,NI_track_df), alpha = .5, color="black") +
  
  geom_path(aes(x = vis_lon, y = lat,group=name,colour=time),data = whole_track_df, alpha = .5) +
  
  #guides(fill=guide_legend(title="Time"))+
  #theme(legend.position = "none")+
  
  labs(title = "Track of Global Cyclone track (1842-2021)",x = "Longitude",y = "Latitude") +
  theme(plot.title = element_text(hjust = 0.5),plot.margin = margin(0.8, 1.2, 0.6, 1.2))+
  scale_colour_gradient(low="red",high="yellow",na.value=NA) 

###############################################  

time_length<-whole_track_df%>%group_by(nature,id)%>%summarise(timelength=length(time))


whole_length<-time_length$timelength
