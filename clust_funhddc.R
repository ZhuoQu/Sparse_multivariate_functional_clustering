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
    funhddc_result <- funHDDC(yfd, K = 2:(sqrt(n / 2) + 5))
    
  } else {
    
  y1 <- sapply(1:n, function(l){   
    result1 <- smooth.basis(pd[[l]]$argvals, pd[[l]]$y[,1], basisobj)
    val <- predict(result1,t)
    return (val)
  })
  
  y2 <- sapply(1:n, function(l){   
    result1 <- smooth.basis(pd[[l]]$argvals, pd[[l]]$y[, 2], basisobj)
    val <- predict(result1,t)
    return (val)
  })
  
  y3 <- sapply(1:n, function(l){   
    result1 <- smooth.basis(pd[[l]]$argvals, pd[[l]]$y[, 3], basisobj)
    val <- predict(result1,t)
    return (val)
  })
  yfd1 <- Data2fd(t, y1)
  yfd2 <- Data2fd(t, y2)
  yfd3 <- Data2fd(t, y3)
  
  funhddc_result <- funHDDC(list(yfd1, yfd2, yfd3), K = 2:(sqrt(n / 2) + 5))
  }
  return (funhddc_result)
}
