
distance_mfd <- function(pd) {
  ### We assume that pd is a list. Each sublist is a list with three components: argvals, subj and y ( a matrix of length|t| * variable numbers)
  DO <- lapply(pd, function(k) {k$y})
  n <- length(pd)
  if ( "numeric" %in% class(pd[[1]]$y) ) {
    nbvar <- 1
  } else {
    nbvar <- ncol(pd[[1]]$y)
  }
  ## define the standard time grid
  sample_time <- seq(0, 1, length.out = max(sapply(pd, function(k) {length(k$argvals)})))
  
  ## define the time grid for each subject
  time_warping_index <- function(index) {
    time_orig <- pd[[index]]$argvals
    time_warping <- sapply(1:length(sample_time), function(k) {
      ind <- which(
        abs(time_orig-sample_time[k]) == min(abs(time_orig-sample_time[k]))
      )
      if (length(ind) > 1) {
        ind <- ind[1]
      }
      return (ind)
    })
    return (time_warping)
  }
  ## define the distance between the rowindex-th curve and the colindex-th curve [rowindex, colindex]
  distance_warping <- function(rowindex, colindex) {
    time_r_index <- time_warping_index(rowindex)
    time_c_index <- time_warping_index(colindex)
    if (nbvar > 1) {
      difference <- DO[[rowindex]][time_r_index, ] - DO[[colindex]][time_c_index, ]
    } else {
      difference <- DO[[rowindex]][time_r_index] - DO[[colindex]][time_c_index]
    }
    if (class(difference) == "numeric" && nbvar > 1) {
      val <- sqrt(sum(difference^2))
    } else if (class(difference) == "numeric" && nbvar == 1) {
      val <- max(abs(difference))
    } else {
      val <- max(apply(difference, 1, function(k) { sqrt(sum(k^2)) } ), na.rm = TRUE)
    }
    return (val)
  }
  
  ############################# calculate the distance matrix
  upt <- lapply(1:(n-1), function (nr) {
    sapply((nr+1):n, function (nc) {
      distance_warping(nr, nc)
    })
  })
  
  numer <- unlist(upt)
  
  distance_matrix <- matrix(0, nrow = n, ncol = n)
  
  distance_matrix[lower.tri(distance_matrix, diag = FALSE)] <- numer
  
  distance_matrix <- t(distance_matrix)
  
  lowerTriangle(distance_matrix) <- upperTriangle(distance_matrix, byrow = TRUE)
  
  return (distance_matrix)
}
