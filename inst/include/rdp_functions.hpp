#ifndef rdp_functions
#define rdp_functions

#include <RcppArmadillo.h>

//
arma::vec ritnb_cpp(const int & n, const arma::vec & mu, const arma::vec & theta, const arma::vec & p, const arma::vec & i, const arma::vec & t, const int & seed);

//
double ditnb_cpp(const int & x, const double & mu, const double & theta, const double & p, const int & i, const int & t);
arma::vec ditnb_cpp(const arma::vec & x, const arma::vec & mu, const arma::vec & theta, const arma::vec & p, const arma::vec & i, const arma::vec & t);

//
arma::vec pitnb_cpp(const arma::vec & x, const arma::vec & mu, const arma::vec & theta, const arma::vec & p, const arma::vec & i, const arma::vec & t);

#endif //rdp_functions
