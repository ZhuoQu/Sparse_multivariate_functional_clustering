silhouette_index <- function (x, dmatrix, minprop, ...) {
  cll <- match.call()
  if (is.list(x) && !is.null(cl <- x$clustering)) 
    x <- cl
  n <- length(x)
  if (!all(x == round(x))) 
    stop("'x' must only have integer codes")
  whole_cluster <- unique(x)
  if (0 %in% whole_cluster) {
    basis_cluster <- unique(as.numeric(x[x != 0]))
  } else {
    basis_cluster <- as.numeric(names(table(x))[as.numeric(table(x)) > minprop * n])
  }
  k <- length(unique(x))
  if (k <= 1 || k >= n) {
    return(NA)
  }
  
  if (missing(dmatrix)) 
    stop("Need a diss.matrix 'dmatrix'")
  if (is.null(dm <- dim(dmatrix)) || length(dm) != 2 || 
      !all(n == dm)) 
    stop("'dmatrix' is not a dissimilarity matrix compatible to 'x'")
  
  wds <- matrix(NA, n, 3, dimnames = list(names(x), c("cluster", 
                                                      "neighbor", "sil_width")))
  for (j in whole_cluster) {
    #cat(j, "\n")
    Nj <- sum(iC <- x == j)
    wds[iC, "cluster"] <- j
    diC <- rbind(apply(dmatrix[x %in% setdiff(basis_cluster, j), iC, drop = FALSE], 2, 
                       function(r) tapply(r, x[x %in% setdiff(basis_cluster, j)], mean)))
    if (length(diC) == 0) {
      s.i <- 0
    } else {
      minC <- apply(diC, 2, which.min)
      if (j %in% basis_cluster) {
        wds[iC, "neighbor"] <- setdiff(basis_cluster, j)[minC]
      } else {
        wds[iC, "neighbor"] <- basis_cluster[minC]
      }
      s.i <- if (Nj > 1) {
        a.i <- colSums(dmatrix[iC, iC])/(Nj - 1)
        b.i <- diC[cbind(minC, seq(along = minC))]
        ifelse(a.i != b.i, (b.i - a.i)/pmax(b.i, a.i), 0)
      }
      else 0
    }

    wds[iC, "sil_width"] <- s.i
  }
  attr(wds, "Ordered") <- FALSE
  attr(wds, "call") <- cll
  class(wds) <- "silhouette"
  wds
}
