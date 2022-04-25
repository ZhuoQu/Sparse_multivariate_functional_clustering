packages <- c("caret","sjmisc", "pcaPP", "dplyr", "roahd", "MASS")
package_check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
nu_construct <- function(lwl, upl, nbvar){
  nu_sp <- matrix(NA, nrow = nbvar, ncol = nbvar)
  diag(nu_sp) <- runif(nbvar, lwl, upl)
  for (i in 1:nbvar) {
    for (j in 1:nbvar) {
      nu_sp[i, j] <- 1 / 2 * (diag(nu_sp)[i] + diag(nu_sp)[j])
    }
  }
  return (nu_sp)
}

beta_construct <- function(nbvar){
  B <- matrix(runif(nbvar^2) / sqrt(nbvar), nrow = nbvar)
  beta <- t(B) %*% B
  diag(beta) <- 1
  return (beta)
} 

rho_construct <- function(beta, nu_sp, nbvar){
  rho <- matrix(NA, nrow = nbvar, ncol = nbvar)
  for (i in 1:nbvar) {
    for (j in 1:nbvar) {
      rho[i,j] <- beta[i, j] * (gamma(diag(nu_sp)[i] + nbvar / 2))^{1 / 2} * 
        (gamma(diag(nu_sp)[j] + nbvar / 2)) ^ {1 / 2} / 
        (gamma(diag(nu_sp)[i]) * diag(nu_sp)[j])^{1 / 2}
      rho[i,j] <- rho[i, j] * 
        gamma(1 / 2 * (diag(nu_sp)[i] + diag(nu_sp)[j])) / 
        gamma(1 / 2 * (diag(nu_sp)[i] + diag(nu_sp)[j] + nbvar))
    }
  }
  return (rho)
}

a_construct <- function(values, nbvar){
  a <- matrix(rep(values, nbvar ^ 2), nrow = nbvar,
              ncol = nbvar, byrow = TRUE) 
  return (a)
}

cross_covariance <- function(nbvar, t, var = c(0.5, 0.5),
                             rho = matrix(c(1, 0.5, 0.5, 1),
                                          nrow = 2, ncol = 2,
                                          byrow = TRUE),
                           nu = matrix(c(0.1, (0.1 + 0.2) / 2,
                                         (0.1 + 0.2) / 2, 0.2),
                                       nrow = 2, ncol = 2, byrow = TRUE),
                          
                            a = matrix(c(2.5, 2.5, 2.5, 2.5),
                                       nrow = 2, ncol = 2, byrow = TRUE)){
  p <- length(t)
  if (length(var) != nbvar || nrow(rho) != nbvar ||
      nrow(nu) != nbvar || nrow(a) != nbvar) {
      paste ("Error!")}
  else {
    H <- matrix(0, nrow = p, ncol = p)
    for (i in 1:p) {
      for (j in 1:p) {
        H[i, j] <- abs(t[i] - t[j]) 
      }
    }

  matern <- array(NA, c(nbvar, nbvar, p, p))

  for (i in 1:nbvar) {
    for (j in 1:nbvar) {
      matern[i, j, , ] <- rho[i, j] * sqrt(var[i] * var[j]) * 
        (a[i, j] * H)^nu[i, j] * besselK(a[i, j] * H, nu[i, j]) /
        (2^(nu[i, j] - 1) * gamma(nu[i, j]))
      diag(matern[i, j, , ]) <- rho[i, j] * sqrt(var[i] * var[j])
    }
  }
   
  #covariance function in time
  cov <- c()
  for (i in 1:nbvar) {
    cov_i <- c()
    for (j in 1:nbvar) {
      cov_i <- cbind(cov_i, matern[i, j, , ])
    }
    cov <- rbind(cov, cov_i)
  }
  return (cov)
 }
}
