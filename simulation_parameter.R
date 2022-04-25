
nbvar <- 3
n <- 150
nbtime <- 50
clusternumber <- 3
t <- seq(0, 1, length.out = nbtime)
repeattimes <- 100
scenario_type <- c("amplitude", "phase", "shift", "clover", "cyclone", "helix")
scenario_type_title <- c("Amplitude", "Phase", "Shift", "Clover", "Cyclone", "Helix")

contamination_level <- 0.1
contam_types <- c("pure magnitude outlier", "peak magnitude outlier", "partial magnitude outlier", "shape outlier I", "shape outlier II", "shape outlier III")
contam_types_title <- c("Pure Magnitude Outlier", "Peak Magnitude Outlier", "Partial Magnitude Outlier", "Shape Outlier I", "Shape Outlier II", "Shape Outlier III")

ps <- 1

minprop <- 0.05
alpha_c <- 0.85

pc <- 0
