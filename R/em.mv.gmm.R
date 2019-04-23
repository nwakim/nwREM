#' Function to evaluate likelihoods under multivariate normal  // HL
#'
#' @param x    (n x d) matrix of observed data    // HL
#' @param mus  d-dimensional vector of means      // HL
#' @param covs (d x d) matrix of covariance       // HL
#' @return a size n vector of log-likelihoods     // HL
#' @export
mv.dnorm <- function(x, mus, covs) {   # function definition // HL
    inv.chol.covs <- chol(solve(covs)) # compute cholesky decomposition of inverse covariance // HL
    det.covs <- det(covs)              # compute determinant
    n <- nrow(x)
    d <- ncol(x)
    m <- x - matrix(mus, n, d, byrow=TRUE) # computes (x - mu) // HL
    xt <- m %*% t(inv.chol.covs)       # transformation to evaluate likelihoods as a dot product // HL
    return( exp(-0.5*rowSums(xt*xt)) / sqrt(det.covs * (2*pi)^d) ) # evaluates likelihoods at once // HL
}

#' E-M Algorithm for Multivariate Gaussian Mixture Model  // HL
#'
#' @param x        A (n x d) matrix of observed data   // HL
#' @param k        The number of mixing components  // HL
#' @param max.iter Maximum number of iterations
#' @param tol      Relative tolerance of likelihood at convergence.
#' @return A list containing the following values:
#' * llk     : log-likelihood value
#' * lambdas : size k vector of estimated mixing proportions
#' * mus     : (d x k) of estimated means
#' * covs    : (d x d x k) vectors of estimated standard deviations
#' * iter    : number of iterations
#' @export
em.mv.gmm <- function(x, k, max.iter = 10000, tol = 1e-8) {
    n       <- nrow(x)                # n = sample size // HL
    d       <- ncol(x)                # d = dimension of data // HL
    
    
    
    Mclust1 = Mclust(x, verbose=F, G=k, modelNames = "VVV", 
                     initialization = list(hcPairs = hc(x)), 
                     control = emControl(tol=c(1.e-5, 1.e-6), itmax=15))
    
    #mix_prop <- rep(1/k, k)            # uniform mixing proportions // HL
    lambdas = Mclust1$parameters$pro
    
    #mus     <- t(x[sample(1:n,k),])   # random starting point as mean // HL
    mus = Mclust1$parameters$mean
    
    #covs    <- array(cov(x),c(d,d,k)) # pooled covariance across all components // HL
    covs = array(unlist(Mclust1$parameters$variance$sigma), c(d,d,k))
    
    
    prevLLK <- -1e300                 # start of E-M loop // HL
    for(i in 1:max.iter) {
        ## E-step to evaluate the fractional counts // HL
        W <- matrix(NA, n, k)
        for(j in 1:k) {               # for each component // HL
            W[,j] <- lambdas[j] * mv.dnorm(x, mus[,j], covs[,,j])  # evaluate likelihoods // HL
        }
        Wsum <- rowSums(W)            # calculate normalization factor // HL
        W    <- W/matrix(Wsum, n, k)  # normalize fractional counts // HL
        llk  <- sum(log(Wsum))        # Evaluate the observed data log-likelihood // HL
        #print(llk)
        if ( llk - prevLLK < tol ) { break }  # check convergence // HL
        prevLLK <- llk

        ## M-step to maximize the expected log-likelihood // HL
        lambdas <- colSums(W) / n                        # update mixing proportions // HL
        mus     <- crossprod(x, W) / matrix(lambdas * n, d, k, byrow=TRUE) # update means // HL
        for(j in 1:k) {                                  # for each component, update covariance // HL
            xc <- (x - matrix(mus[,j], n, d, byrow=TRUE)) * matrix(sqrt(W[,j]),n,d)         # // HL
            covs[,,j] <- crossprod(xc,xc)/(lambdas[j]*n)      # efficiently comput in O(n^2d) // HL
        }
    }                                 # end of E-M loop // HL

    classification = apply(W, 1, which.max)
    
    return ( list(llk=llk,            # log-likelihood
                  lambdas=lambdas,    # mixing proportion
                  mus=mus,            # (d x k) matrix of means
                  covs=covs,          # (d x d x k) matrices of covariance
                  probs=W,            # normalized probabilities
                  classification = classification, 
                  iter=i) )           # iterations
}
