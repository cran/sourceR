#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector like_calc_pois(int no_J, int no_I, int no_T, int no_L, List lambda_can, double d, List dat) { 
                
  double ll = 0.0; //nest times then locations

  for (int t = 0; t < no_T; ++t) {
    
    List lambda_can_t = lambda_can[t];
    List dat_t = dat[t]; //dat is the human data
    
    for (int l = 0; l < no_L; ++l) {
      
      arma::vec data_t_l = Rcpp::as<arma::vec>(dat_t[l]);
      arma::vec lambda_can_t_l = Rcpp::as<arma::vec>(lambda_can_t[l]);
      ll += sum((data_t_l % log(lambda_can_t_l)) - lambda_can_t_l);
      
    }
  }

  return wrap(ll);
}
