#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]

List calc_li(int no_J, int no_I, int no_T, int no_L, List r, List a, List prev, NumericVector q) { 

  List lambdai(no_T); //nest times then locations

  arma::vec q1 = Rcpp::as<arma::vec>(q);

  for (int t = 0; t < no_T; ++t) {

    arma::vec prev1 = Rcpp::as<arma::vec>(prev[t]);
    lambdai[t] = List(no_L);

    for (int l = 0; l < no_L; ++l) {

      List tmp = lambdai[t]; // this has a pointer back to the original List (lambdai), so updating tmp[l], updates lambdai[t][l]
      List at = a[t];

      arma::mat r1 = Rcpp::as<arma::mat>(r[t]);
      arma::vec a1 = Rcpp::as<arma::vec>(at[l]);

      tmp[l] = q1 % (r1 * (a1 % prev1)); // * does matrix multiplication; % does element wise multiplication

    }
  }

  return lambdai;
}

