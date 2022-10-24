source("cluster_generation.R")
source("sparsify_scenario.R")
source("00_subsetlist.R")

nbtime <- 15
n <- 20

example <- generate_cluster(nbvar = 2, n = 20, nbtime = 15, 
                            clusternumber = 2, variat = "phase")

example[[2]][1, ] <- example[[2]][1, ] + 4
example[[2]][2, ] <- example[[2]][2, ] - 2
example[[2]][3, ] <- example[[2]][3, ] + 3
example[[2]][5, ] <- example[[2]][5, ] -5
common_grid <- seq(0, 1, length.out = nbtime)

observed_index <- sapply(1:n, function(l) {
  sort(unique(sample(nbtime, nbtime, replace = TRUE)))})
pdf(file = "etd_example.pdf", width = 8, height = 4)
par(mfrow = c(1, 2), 
    #mar = c(6, 4, 3, 3), 
    mai = c(0.7, 0.8, 0.4, 0.1), mar = c(4, 3.5, 2, 0.5),
    mgp = c(2.5, 1.1, 0))
    #mai = c(0.8, 0.8, 0.4, 0.1))  
plot(common_grid, example[[2]][1, ], type = "n",
     ylim = range(example[[2]]), xlim = c(0, 1.03),
     xlab = "Time", ylab = "Value", main = "Univariate Functional Observations",
     cex.main = 0.95)

select_index <- 1:5
lapply(select_index, function(l) {
  time_index <- observed_index[[l]]
  lines(common_grid[time_index], example[[2]][l, time_index], col = l + 1, lty = l)})

lapply(select_index, function(l) {
  time_index <- observed_index[[l]]
  points(common_grid[time_index], example[[2]][l, time_index], col = "black", cex = 0.5)
  text(1.03, y = example[[2]][l, rev(time_index)[1]], 
       labels = l, cex = 0.7)
  })


standard_grid <- seq(0, 1, length.out = max(unlist(lapply(subsetlist(observed_index, select_index), length))))
  
plot(standard_grid, example[[2]][1, 1:length(standard_grid)], type = "n",
     ylim = range(example[[2]]), xlim = c(0, 1.03),
     xlab = "Time", ylab = "Value", main = "Interpolated Univariate Functional Observations",
     cex.main = 0.87)

lapply(select_index, function(l) {
  time_index <- observed_index[[l]]
  obs_time <- common_grid[time_index]
  update_time <- sapply(1:length(standard_grid), function(l) {
    which(abs(standard_grid[l] - obs_time) == min(abs(standard_grid[l] - obs_time)))
    })
  new_time <- c(obs_time, standard_grid)
  rematch_index <- order(new_time, decreasing = FALSE)
  new_obs <- c(example[[2]][l, time_index], example[[2]][l, time_index[update_time]])
  #lines(new_time[rematch_index], new_obs[rematch_index], col = l + 1, lty = l)
  lines(standard_grid, example[[2]][l, time_index[update_time]], col = l + 1, lty = l)
  })

lapply(select_index, function(l) {
  time_index <- observed_index[[l]]
  obs_time <- common_grid[time_index]
  update_time <- sapply(1:length(standard_grid), function(l) {
    which(abs(standard_grid[l] - obs_time) == min(abs(standard_grid[l] - obs_time)))
  })
  new_time <- c(obs_time, standard_grid)
  rematch_index <- order(new_time, decreasing = FALSE)
  new_obs <- c(example[[2]][l, time_index], example[[2]][l, time_index[update_time]])
  #points(new_time[rematch_index], new_obs[rematch_index], col = "black", cex = 0.5)
  points(standard_grid, example[[2]][l, time_index[update_time]], col = "black", cex = 0.5)
  text(1.03, y = rev(example[[2]][l, observed_index[[l]][update_time]])[1], 
       labels = l, cex = 0.7)
  })
dev.off()

