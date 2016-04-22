#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]

List calc_lj(int no_J, int no_I, int no_T, int no_L, List r, List a, List prev, NumericVector q) { 

  List lambdaj(no_T); //nest times then locations
  NumericVector colSums(no_I);

  arma::vec q1 = Rcpp::as<arma::vec>(q);

  for (int t = 0; t < no_T; ++t) {

    arma::vec prev1 = Rcpp::as<arma::vec>(prev[t]);
    lambdaj[t] = List(no_L);

    for (int l = 0; l < no_L; ++l) {

      List tmp = lambdaj[t]; // this has a pointer back to the original List (lambdai), so updating tmp[l], updates lambdai[t][l]
      List at = a[t];

      arma::mat r1 = Rcpp::as<arma::mat>(r[t]);
      arma::vec a1 = Rcpp::as<arma::vec>(at[l]);

      NumericMatrix qr(no_I, no_J);
      arma::mat qr1 = Rcpp::as<arma::mat>(qr);
    
      for (int j = 0; j < no_J; ++j) {
        qr1.col(j) = r1.col(j) % q1;
      }

      colSums = arma::sum(qr1);
      arma::vec colSums1 = Rcpp::as<arma::vec>(colSums);

      tmp[l] = (a1 % prev1) % colSums1; // * does matrix multiplication; % does element wise multiplication
    }
  }

  return wrap(lambdaj);
}
