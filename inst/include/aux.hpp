#ifndef aux_functions
#define aux_functions

#include <RcppArmadillo.h>
#include "links.hpp"

//// Log-likelihood
double loglike_pois(const arma::mat & X, const arma::vec & y, const arma::vec & beta, const double & p, const int & i, const int & t, const int & N, Link & LO);
double loglike_nb(const arma::mat & X, const arma::vec & y, const arma::vec & beta, const double & theta, const double & p, const int & i, const int & t, const int & N, Link & LO);

//// Derivative of beta function
double beta_derivative_rectangle(const double & a, const double & b, const double & theta, const int & t, const int & N);
double beta_derivative_boole(const double & a, const double & b, const double & theta, const int & t, const int & N);

//// EM-steps
void update_z(arma::vec & z, const arma::mat & X, const arma::vec & y, const arma::vec & yi, const arma::vec & beta, const double & theta, const double & p, const int & i, const int & t, const int & N, Link & LO, const bool & is_poisson);
void update_p(double & p, const arma::vec & yi, const arma::vec & z, const int & N);

#endif //aux_functions
