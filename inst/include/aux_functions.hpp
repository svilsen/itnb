#ifndef aux_functions
#define aux_functions

#include <RcppArmadillo.h>
#include "links.hpp"

//// Derivative of beta function
double beta_derivative_rectangle(const double & a, const double & b, const double & theta, const int & t, const int & N);
double beta_derivative_boole(const double & a, const double & b, const double & theta, const int & t, const int & N);

//// EM-steps
void update_z(arma::vec & z, const arma::mat & X, const arma::vec & y, const arma::vec & yi, const arma::vec & beta, const double & theta, const double & p, const int & i, const int & t, const int & N, Link & LO);
void update_p(double & p, const arma::vec & yi, const arma::vec & z, const int & N);

//// Complete log-likelihoods
double complete_loglikelihood_itpois(const arma::mat & X, const arma::vec & y, const arma::vec & beta, const double & p, const int & i, const int & t, Link & LO);
double complete_loglikelihood_itnb(const arma::mat & X, const arma::vec & y, const arma::vec & beta, const double & theta, const double & p, const int & i, const int & t, Link & LO);

#endif //aux_functions
