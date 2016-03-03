#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
//' @param X vector of test.statistics at time t
//' @param n_iter number of iterations of the gibbs sampler
//' @param n_burnin number of burn in iterations
//' @param sampling_interval interval of sampling times
//' @return 
// [[Rcpp::export]]
arma::mat gibbs_sampler(arma::vec X, int n_iter, int n_burnin, int sampling_interval) {
  int n = X.size();
  int a = 1; // gamma prior hyper param
  int b = 1; // gamma prior hyper param
  int m = 0; // normal prior hyper param
  int alpha = 1;
  int beta = 1;
  int alpha_norm = 2; // normal hyper param 
  double pi0 = R::rbeta(alpha, beta);
  double phi1 = R::rgamma(a/2, b/2);
  int mu1 = R::rnorm(m, 1 / (alpha_norm * phi1));
  arma::mat S(n_iter, 3);
  
  //for (i=1; i<n_iter; i++){
    //p_z_num = pi0 * exp(-.5*(X^2));
  //}
  return(S);
}