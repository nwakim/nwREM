% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{hc_fn}
\alias{hc_fn}
\title{This is an Rcpp implementation of the hc function from the mclust package. This will be used while implementing the Mclust function.}
\usage{
hc_fn(x)
}
\arguments{
\item{x}{A (n x d) matrix of observed data}
}
\value{
A numeric two-column matrix in which the ith row gives the minimum index for observations in each of the two clusters merged at the ith stage of agglomerative hierarchical clustering.
}
\description{
This is an Rcpp implementation of the hc function from the mclust package. This will be used while implementing the Mclust function.
}
