sparsify_scenario <- function(ps, pc, scenario) {
  ## ps is a parameter specifying the proportion of sparseness among samples
  ## pc is a parameter specifying the proportion of sparseness among curves with missing values
  ## scenario is a list with the data and the outlier index
  ## the first list in scenario is with length p, and each sublist is a n * nbtime matrix.
  n <- nrow(scenario[[1]][[1]])
  nbvar <- length(scenario[[1]])
  nbtime <- ncol(scenario[[1]][[1]])
  t <- seq(0, 1, length.out = nbtime)
  
  PD <- lapply(1:n, function(i){
    sample_curve <- sample(1:n, ps * n)
    if (i %in% setdiff(1:n, sample_curve) ) {
      observed_index <- 1:nbtime
    } else {
      observed_index <- sort(sample(1:nbtime, (1 - pc) * nbtime))
    }
    values <- scenario[[1]][[1]][i, ]
    if (nbvar >= 2) {
      for (kk in 2:nbvar) {
        values <- cbind(values, scenario[[1]][[kk]][i, ])
      }
    }
    return (list(argvals = t[observed_index],
                 subj = rep(i, length(observed_index)), 
                 y = values[observed_index, ]))
    
  })
  
  marginal_PD <- function(nbvar) {
    result <- PD
    for (l in 1:length(result)) {
      result[[l]]$y <- PD[[l]]$y[, nbvar]
    }
    return (result)
  }
  
  result <- vector(mode = "list", length = nbvar + 1)
  result[[1]] <- PD
  for (l in 2:(nbvar + 1)) {
    result[[l]] <- PD
    result[[l]] <- marginal_PD(l - 1)
  }
  return (result)
}
