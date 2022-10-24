packages = c("funHDDC", "fda", "stats")
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
clust_funhddc <- function (pd, pc) {
  n <- length(pd)
  nbvar <- ifelse(class(pd[[1]]$y)[1] == "numeric", 1, length(pd[[1]]$y[1, ]))
  if (pc < 0.2) {
    basisobj <- create.bspline.basis(c(0,1), nbasis = 15)
  } else if (pc < 0.5) {
    basisobj <- create.bspline.basis(c(0,1), nbasis = 10)  
  } else if (pc < 0.8) {
    basisobj <- create.bspline.basis(c(0,1), nbasis = 6)
  } else {
    basisobj <- create.bspline.basis(c(0,1), nbasis = 4)
  }
  
  if (length(ncol(pd[[1]]$y)) == 0) {
    y <- sapply(1:n, function(l){   
      result1 <- smooth.basis(pd[[l]]$argvals, pd[[l]]$y, basisobj)
      val <- predict(result1,t)
      return (val)
    })
    yfd <- Data2fd(t, y)
    
  } else {
    
  yfd <- lapply(1:nbvar, function(num) {
    temp <- sapply(1:n, function(l){   
    result1 <- smooth.basis(pd[[l]]$argvals, pd[[l]]$y[, num], basisobj)
    val <- predict(result1, t)
    return (val)
  })
    Data2fd(t, temp)
    })
  }
  
  funhddc_result <- funHDDC(yfd, K = 2:(sqrt(n / 2) + 5))
  
  return (funhddc_result)
}
