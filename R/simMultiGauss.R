#' Simulate the centers for each Gaussian in the mutlivariate Gaussian mixture model
#' 
#' @param k       Number of different Gaussians
#' @param d       Number of dimensions for multivariate Gaussian mixture model
#' @return mus  : A (d x k) matrix of mean parameters
#' @export
createMus <- function(k, d) {
  mus = replicate(k, runif(d,1,60*(1/d)))
  return(mus)
}
#' Simulate the array of covariance parameters for the mutlivariate Gaussian mixture model. This function is adapted from https://stat.ethz.ch/pipermail/r-help/2008-February/153708
#' 
#' @param k       Number of different Gaussians
#' @param d       Number of dimensions for multivariate Gaussian mixture model
#' @return mus  : A (d x k) matrix of mean parameters
#' @export
createCovs <- function(k, d) {
  covList = array(dim=c(d, d, k))
  for (i in 1:k) {
    # https://stat.ethz.ch/pipermail/r-help/2008-February/153708
    Z <- matrix(ncol=d, rnorm(d^2))
    decomp <- qr(Z)
    Q <- qr.Q(decomp) 
    R <- qr.R(decomp)
    diag1 <- diag(R)
    ph <- diag1 / abs(diag1)
    O <- Q %*% diag(ph)
    ev = runif(d, 0, 10)
    Z <- t(O) %*% diag(ev) %*% O
    covList[,,i] = Z
  }
  return(covList)
}
#' Simulate the mixing proportions for the Gaussian mixture model
#' 
#' @param k       Number of different Gaussians
#' @return mix_prop : k-dimensional Vector of mixing proportions (to be normalized)
#' @export
createMixProp <- function(k) {
  mix_prop = sample(1:k, k, replace=T)
  return(mix_prop)
}
#' Simulate Random Samples from Multivariate Gaussian Mixture Model
#' 
#' @param n       Number of random samples to simulate
#' @param mix_prop k-dimensional Vector of mixing proportions (to be normalized)
#' @param mus     (d x k) Matrix of mean parameters
#' @param covs    (d x d x k) Multi-way array of covariance parameters
#' @param perc_out  Percent (\%) of random samples that will be outliers
#' @return A list containing the following values
#' * x : A matrix of simulated random samples from multi-dimensional Gaussian mixtures
#' * z : A vector of hidden variables indicating the source components
#' @export
simMultiGauss <- function(n, mix_prop, mus, covs, perc_out) {  # function definition // HL
  k   <- length(mix_prop)          # k is the number of components
  d   <- nrow(mus)                # d is the dimension of data
  
  stopifnot(k == ncol(mus))       # sanity checking
  stopifnot(dim(covs) == c(d,d,k)) # sanity checking
  
  cov.chols <- array(NA, c(d,d,k))
  for(i in 1:k) { cov.chols[,,i] <- chol(covs[,,i]) }
  
  z  <- apply(rmultinom(n, 1, mix_prop), 2, which.max) # simulate hidden labels // HL
  y <- matrix(rnorm(n*d), n, d)  # simulate i.i.d. from standard normal // HL
  x  <- matrix(NA, n, d)
  for(i in 1:(n-ceiling(n*perc_out/100))) {
    x[i,] <- (y[i,] %*% cov.chols[,,z[i]]) + mus[,z[i]] # convert to MVN // HL
  }
  if (perc_out != 0) {
    for(i in (n-ceiling(n*perc_out/100)+1):n) {
      x[i,] <- (y[i,] %*% (3*cov.chols[,,z[i]])) + mus[,z[i]] # convert to MVN // HL
    }
  }
  return(list(x=x,z=z))             # return x and z as a list
}
