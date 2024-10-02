#ifndef em_functions
#define em_functions

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

void optimise_itnb(
        arma::vec & beta_j, double & theta_j, double & p_j, double & loglike_j,
        int & j, bool & not_converged, std::string & convergence_flag,
        arma::mat & beta_trace, arma::vec & theta_trace, arma::vec & p_trace, arma::vec & loglike_trace,
        const arma::mat & X, const arma::vec & y, const arma::vec yi,
        const int & i, const int & t, const int & N, const int & M,
        const int & iteration_min, const int & iteration_max,
        const double & tolerance, const arma::vec & lambda,
        const int & steps, const bool & exact,
        const int & trace, const bool & save_trace
);

#endif //em_functions
