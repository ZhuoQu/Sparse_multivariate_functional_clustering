source("../00_subsetlist.R")
wp_index1 <- westpacific[unlist(lapply(subsetlist(clust_westpacific$clust, 1), function(l) {l$element}))]
wp_index2 <- westpacific[unlist(lapply(subsetlist(clust_westpacific$clust, 2), function(l) {l$element}))]
wp_1_data <- subset(whole_track_df, id %in% wp_index1)
wp_2_data <- subset(whole_track_df, id %in% wp_index2)

############################
dur1 <- wp_1_data %>% group_by(id) %>% summarize(mean(duration))
dur2 <- wp_2_data %>% group_by(id) %>% summarize(mean(duration))

month1 <- wp_1_data %>% group_by(id) %>% summarize(unique(month))
month2 <- wp_2_data %>% group_by(id) %>% summarize(unique(month))

year1 <- wp_1_data %>% group_by(id) %>% summarize(unique(year))
year2 <- wp_2_data %>% group_by(id) %>% summarize(unique(year))

# year_c1<-hist(subset(wp_1_data, year > 1940)$year, 
#               xlim = c(1940, 2030), ylim = c(0, 0.018),
#               freq = FALSE, breaks = seq(1940, 2030, by = 10),
#               cex.lab = 1.1, cex.main = 1.1, xlab = "Year",
#               main = "Year of Cyclone (Cluster 1)")
# year_c2<-hist(subset(wp_2_data, year>1940)$year, 
#               xlim = c(1940, 2030), ylim = c(0, 0.018),
#               freq = FALSE, breaks = seq(1940, 2030, by=10),
#               cex.lab = 1.1, cex.main = 1.1, xlab = "Year",
#               main = "Year of Cyclone (Cluster 2)")

#par(mfrow = c(1, 3), mar = c(1, 3, 0.5, 0.5), mgp = c(1.8, 0.8, 0))
pdf(file = "hist.pdf", width = 10.3, height = 3.5)
par(mfrow = c(1, 3), mar = c(6.3, 3.3, 2.5, 1), mgp = c(2.2, 1, 0))
#par(mfrow = c(1, 3), mai = c(0.5, 0.55, 0.35, 0.07), mar = c(3.5, 3.5, 2.5, 1), mgp=c(2.3, 1, 0))

# dat <- data.frame(dens = c(dur1$`mean(duration)`, dur2$`mean(duration)`)
#                   , cluster = c(rep(1, length(dur1$`mean(duration)`)),
#                                 rep(2, length(dur2$`mean(duration)`)) )
#                   )
#Plot.
##ggplot(dat, aes(x = dens, fill = cluster)) + geom_density(alpha = 0.5)

# d1 <- hist(dur1$`mean(duration)`, breaks = 12,
#      plot = FALSE)
# d2 <- hist(dur2$`mean(duration)`, breaks = 12, plot = FALSE)
# plot(d1, freq = FALSE, ylim = c(0, 0.15), cex.lab = 1.3, 
#      cex.main = 1.3, xlim = c(0,25), 
#       xlab = "Time (Days)", 
#       main = "Density Estimate of Cyclone Duration", col = "red") 
# plot(d2, freq = FALSE, add = TRUE, col = "green") 
# #plot(d1, freq = FALSE, add = TRUE, col = "red") 
# lines(x = density(dur1$`mean(duration)`), lwd = 3, lty = 2, col = "red")
# lines(x = density(dur2$`mean(duration)`), lwd = 3, lty = 2, col = "green")

plot(NA, xlim = range(dur1$`mean(duration)`), ylim = c(0, 0.12), 
     main = "Density Estimate of Cyclone Duration", cex.main = 1.5, 
     xlab = "Time (Days)", ylab = "Density", cex.lab = 1.4,
     cex.axis = 1.4)
lines(x = density(dur1$`mean(duration)`), lwd = 2, col = "black")
lines(x = density(dur2$`mean(duration)`), lwd = 2, col = "grey")

y1 <- hist(year1$`unique(year)`, breaks = seq(1860, 2030, by = 10), plot = FALSE)
           
y2 <- hist(year2$`unique(year)`, 
                breaks = seq(1860, 2030, by = 10), plot = FALSE)

# plot(y1, xlim = c(1880, 2030), 
#      ylim = c(0, 0.018), freq = FALSE, 
#      cex.lab = 1.3, 
#      cex.main = 1.3, xlab = "Year", 
#      main = "Density Estimate of Cyclone Year", col = "red")
# plot(y2, freq = FALSE, add = TRUE, col = "green") 
# lines(x = density(year1$`unique(year)`), col = "red", lwd = 3, lty = 2)
# lines(x = density(year2$`unique(year)`), col ="green", lwd = 3, lty = 2)

plot(NA, xlim = c(1850, 2021), ylim = c(0, 0.017), 
     main = "Density Estimate of Cyclone Year", cex.main = 1.5, 
     xlab = "Year", ylab = "Density", cex.lab = 1.4, cex.axis = 1.4)
lines(x = density(year1$`unique(year)`), col = 1, lwd = 2)
lines(x = density(year2$`unique(year)`), col = "grey", lwd = 2)

m1 <- hist(month1$`unique(month)`, breaks = seq(0, 12, by = 1),
           plot = FALSE)
m2 <- hist(month2$`unique(month)`, breaks = seq(0, 12, by = 1),
           plot = FALSE)
# plot(m1, xlim = c(0, 12), ylim = c(0, 0.35), freq = FALSE, cex.lab = 1.3, cex.main = 1.3, 
#      xlab= "Months", main="Density Estimate of Cyclone Month", col = "red") 
# plot(m2, freq = FALSE, add = TRUE, col = "green") 
# #plot(d1, freq = FALSE, add = TRUE, col = "red") 
# lines(x = density(month1$`unique(month)`), lwd = 3, lty = 2, col = "red")
# lines(x = density(month2$`unique(month)`), lwd = 3, lty = 2, col = "green")

plot(NA, xlim = c(0, 12), ylim = c(0, 0.32), 
     main = "Density Estimate of Cyclone Month", cex.main = 1.5,
     xlab = "Months", ylab = "Density", cex.lab = 1.4, cex.axis = 1.4)
lines(x = density(month1$`unique(month)`), lwd = 2, col = 1)
lines(x = density(month2$`unique(month)`), lwd = 2, col = "grey")

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
# Add a legend
legend("bottom", 
       legend = c("Cluster 1", "Cluster 2"), 
       col = c("black", "grey"), lwd = rep(1, 2), 
       cex = 1.3, text.font = 1, 
       horiz = T, inset = c(0.01, 0.01))
dev.off()
