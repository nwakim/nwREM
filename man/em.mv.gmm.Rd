% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/em.mv.gmm.R
\name{em.mv.gmm}
\alias{em.mv.gmm}
\title{E-M Algorithm for Multivariate Gaussian Mixture Model  // HL}
\usage{
em.mv.gmm(x, k, max.iter = 10000, tol = 1e-08)
}
\arguments{
\item{x}{A (n x d) matrix of observed data   // HL}

\item{k}{The number of mixing components  // HL}

\item{max.iter}{Maximum number of iterations}

\item{tol}{Relative tolerance of likelihood at convergence.}
}
\value{
A list containing the following values:
* llk     : log-likelihood value
* lambdas : size k vector of estimated mixing proportions
* mus     : (d x k) of estimated means
* covs    : (d x d x k) vectors of estimated standard deviations
* iter    : number of iterations
}
\description{
E-M Algorithm for Multivariate Gaussian Mixture Model  // HL
}
