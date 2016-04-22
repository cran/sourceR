#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector like_calc_nbinom(int no_J, int no_I, int no_T, int no_L, List lambda_can, double d, List dat) { 

  double ll = 0.0; //nest times then locations

  for (int t = 0; t < no_T; ++t) {
    List lambda_can_t = lambda_can[t];
    List dat_t = dat[t]; //dat is the human data
    
    for (int l = 0; l < no_L; ++l) {
      
      NumericVector data_t_l = dat_t[l];
      NumericVector lambda_can_t_l = lambda_can_t[l];
      double sum = 0.0;
      
      for (int i = 0; i < no_I; ++i) {
        
        sum += Rf_dnbinom_mu(data_t_l[i], d, lambda_can_t_l[i], true);
        
      }
      
      ll += sum;
      
    }
  }

  return wrap(ll);
}
