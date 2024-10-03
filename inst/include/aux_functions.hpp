#ifndef aux_functions
#define aux_functions

#include <RcppArmadillo.h>
#include "link_functions.hpp"

//// Log-likelihood
double loglikelihood(const arma::mat & X, const arma::vec & y, const arma::vec & beta, const double & theta, const double & p, const int & i, const int & t, const int & N, Link & LO);

//// Derivative of beta function
double beta_derivative_rectangle(const double & a, const double & b, const double & theta, const int & t, const int & N);
double beta_derivative_boole(const double & a, const double & b, const double & theta, const int & t, const int & N);

#endif //aux_functions
