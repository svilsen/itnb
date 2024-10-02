#ifndef mle_functions
#define mle_functions

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

void optimise_tnb(
        arma::vec & beta_j, double & theta_j, const double & p_0, double & loglike_j,
        int & j, bool & not_converged, std::string & convergence_flag,
        const arma::mat & X, const arma::vec & y, const arma::vec yi,
        const int & i, const int & t, const int & N, const int & M,
        const double & tolerance, const arma::vec & lambda,
        const int & steps, const bool & exact
);

#endif //mle_functions
