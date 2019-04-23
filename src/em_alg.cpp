#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

//' This is an Rcpp implementation of the hc function from the mclust package. This will be used while implementing the Mclust function.
//' @param x    A (n x d) matrix of observed data 
//' @return A numeric two-column matrix in which the ith row gives the minimum index for observations in each of the two clusters merged at the ith stage of agglomerative hierarchical clustering.
// [[Rcpp::export]]
NumericMatrix hc_fn(NumericMatrix x){
  
  // Obtaining namespace of mclust package
  Environment pkg = Environment::namespace_env("mclust");
  
  // Picking up hc() function from Matrix package
  Function f = pkg["hc"];
  
  // Executing hc(x)
  return f( Named("data", x) );
}
//' This is an Rcpp implementation of the emControl function from the mclust package. This will be used while implementing the Mclust function.
//' @param itmax  Maximum number of iterations for Mclust()
//' @return A named list in which the names are the names of the arguments and the values are the values supplied to the arguments.
// [[Rcpp::export]]
List emControl_fn(int itmax){
  
  // Obtaining namespace of mclust package
  Environment pkg = Environment::namespace_env("mclust");
  
  // Picking up hc() function from Matrix package
  Function f = pkg["emControl"];
  
  // Executing hc(x)
  return f( Named("itmax", itmax) );
}

//' Implements the Mclust() function from library mclust.
//' @param x    A (n x d) matrix of observed data 
//' @param k    The number of mixing components
//' @return A list of parameter estimates that will be used for the initial parameter estimates in the EM algorithm.
// [[Rcpp::export]]
List Mclust_fn(NumericMatrix x, int k){

  // Obtaining namespace of mclust package
  Environment pkg = Environment::namespace_env("mclust");
  
  // Picking up Mclust() function from Matrix package
  Function f = pkg["Mclust"];
  
  // Executing Mclust(x, verbose=F, G=k, modelNames = "VVV", 
  //                  initialization = list(hcPairs = hc(x)), 
  //                  control = emControl(tol=c(1.e-5, 1.e-6), itmax=15)))
  int itmax = 15;
  return f( Named("data", x), Named("verbose", false), 
            Named("G", k), Named("modelNames", "VVV"), 
            Named("initialization", List::create(Named("hcPairs") = hc_fn(x))), 
            Named("control", emControl_fn( itmax ) ) );
}
//' Function to evaluate likelihood or log-likelihood of multivariate normal using armadillo
//' @param x    (n x d) matrix of observed data 
//' @param mean  Row vector (d) of means   
//' @param sigma (d x d) matrix of covariance
//' @param logd boolean indicating log-likelihood (true) or likelihood (false) is returned
//' @return a size n vector of log-likelihoods
// [[Rcpp::export]]
arma::vec dmvnrm_arma(arma::mat x,  
                      arma::rowvec mean,  
                      arma::mat sigma, 
                      bool logd = false) { 
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  
  const double log2pi = std::log(2.0 * M_PI);
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;    
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;     
  }  
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}
//' Function to calculate the Mahalanobis distance. Credit: Rcpp gallery.
//' @param x    (n x d) matrix of observed data 
//' @param center  Row vector (d) of means   
//' @param cov (d x d) matrix of covariance
// [[Rcpp::export]]
arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov) {
  int n = x.n_rows;
  arma::mat x_cen;
  x_cen.copy_size(x);
  for (int i=0; i < n; i++) {
    x_cen.row(i) = x.row(i) - center;
  }
  return sum((x_cen * cov.i()) % x_cen, 1);    
}

//' Implements the EM algorithm on mvGMM data
//' @param x        A (n x d) matrix of observed data 
//' @param k        The number of mixing components
//' @param lambda   Final value for regularization parameter - will go lambda^10 to lambda. Lambda will automatically equal infinity, if not specified. This is equivalent to a standard EM algorithm
//' @param max_it   Maximum number of iterations
//' @param tol      Relative tolerance of likelihood at convergence
//' @return A list containing the following values:
//' * llk     : Likelihood of final iteration
//' * mus     : (d x k) of estimated means
//' * mix_prop : size k vector of estimated mixing proportions
//' * covs    : (d x d x k) vectors of estimated standard deviations
//' * probs : Matrix (n x k) of the probabilties of observed values in each cluster
//' * classification Size n vector of most likely cluster for each observed value
//' * iter    : number of iterations
// [[Rcpp::export]]
List em_alg_GMM(NumericMatrix& x, int k, int lambda, 
                int max_it, double tol) {
  
  int n = x.nrow();
  int d = x.ncol();
  
  arma::mat x1 = Rcpp::as<arma::mat>(x);
  
  // ######### Initial values - random sigma and tau ###########
  // Initialize mix_prop
  arma::vec mix_prop(k);  

  // Initialize errors
  arma::mat err(n, d, fill::zeros); 

  // Initialize mus
  arma::mat mus(d, k);

  // Run Mclust to get initial parameter values
  List Mclust1 = Mclust_fn(x, k);
  List param = Mclust1["parameters"];
  
  // Create mixing proportions from Mclust
  mix_prop = Rcpp::as<arma::vec>(param["pro"]);
  // Create mus from Mclust
  mus = Rcpp::as<arma::mat>(param[1]);
  // Create covs cube from Mclust
  List var = param["variance"];
  NumericVector vecArray = var["sigma"];
  IntegerVector arrayDims = vecArray.attr("dim");
  arma::cube covs(vecArray.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
  
  double prevLLK = -1e300;
  double llk;
  int iter; int j;
  mat xc(n, d);
  mat mal(n, k);
  mat W(n, k);
  
  for(iter = 0; iter < max_it; ++iter) {
    
    // E-step to evaluate the fractional counts
    for(j = 0; j < k; ++j) {
      rowvec rowMus = conv_to<rowvec>::from(mus.col(j));

      vec mvdnorm = dmvnrm_arma(x1 - err, rowMus, covs.slice(j), false);
      
      W.col(j) = mix_prop(j) * mvdnorm; 
    }
    
    vec Wsum(n, 1);
    Wsum = sum(W, 1);
    mat Wsum_mat(n,k);
    Wsum_mat.each_col() = Wsum;
    W = W/Wsum_mat;
    
    vec err_rowSum = sum(abs(err), 1);
    int sum_err = 0; 
    for ( int p = 0; p < n; ++p) {
      if ( err_rowSum(p) > 0 ) { sum_err += 1;}
    }
    llk = sum( log(Wsum) ) - (lambda/2) * sum_err;       
    if ( llk - prevLLK < tol ) { break; }  
    prevLLK = llk;
    
    //M-step to maximize the expected log-likelihood
    mix_prop = conv_to<vec>::from( sum(W, 0) / n );    
    mat mu_denom(d,k);
    mu_denom.each_row() = conv_to<rowvec>::from(mix_prop * n);
    mus = (x1-err).t() * W / mu_denom; 
    
    mat mu_rep(n, d);
    mat sqrt_W(n, d);
    rowvec mu_rowvec(d);
    for(j = 0; j < k; ++j ) { 
      mu_rowvec = conv_to<rowvec>::from(mus(span::all, j));
      mu_rep.each_row() = mu_rowvec;
      sqrt_W.each_col() = sqrt(W(span::all, j));
      xc = ((x1-err) - mu_rep) % sqrt_W;
      covs.slice(j) = xc.t()*xc / (mix_prop(j)*n);
      mal(span::all, j) = Mahalanobis(x1, mu_rowvec, covs.slice(j));
    }
    err.zeros();

    uvec k_min(n);
    k_min = index_min(mal, 1);
    vec min_dist(n);
    min_dist = min( mal, 1 ); 
    uvec ind = find( min_dist > lambda );
    mat mu_err = mus.cols(k_min(ind));
    err.rows(ind) = x1.rows(ind) - mu_err.t();
  }
  
  uvec classification = index_max(W, 1) + 1;

  return(List::create( Named("llk") = llk,
                       Named("mus") = mus,
                       Named("mix_prop") = mix_prop,
                       Named("covs") = covs,
                       Named("probs") = W,
                       Named("classification") = classification,
                       Named("iter") = iter));
}