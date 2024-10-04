#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "rdp_functions.hpp"
#include "link_functions.hpp"

//// Log-likelihood
double loglikelihood(const arma::mat & X, const arma::vec & y, const arma::vec & beta, const double & theta, const double & p, const int & i, const int & t, const int & N, Link & LO) {
    //
    double d = 0.0;
    for (int n = 0; n < N; n++) {
        const arma::rowvec & x_n = X.row(n);
        const arma::vec eta_n = x_n * beta;
        const double & mu_n = LO.link_inv(eta_n[0]);

        d += ditnb_cpp(y[n], mu_n, theta, p, i, t);
    }

    return d;
}

//// Derivative of beta function
double beta_derivative_integrand(const double & s, const double & t, const double & theta) {
    //
    const double s_rev = 1.0 - s;
    double r = 0.0;
    if ((s > 1e-8) & (s < (1.0 - 1e-8))) {
        r = std::pow(s, t) * std::pow(s_rev, theta - 1) * std::log(s_rev);
    }

    return r;
}

double beta_derivative_rectangle(const double & a, const double & b, const double & theta, const int & t, const int & N) {
    //
    const double delta_x = (b - a) / N;
    double i = 0.0;
    for (int n = 0; n < N; n++) {
        i += beta_derivative_integrand(a + n * delta_x, t, theta) * delta_x;
    }

    return i;
}

double beta_derivative_boole(const double & a, const double & b, const double & theta, const int & t, const int & N) {
    //
    const double delta_x = (b - a) / N;

    //
    double i = 0.0;
    for (int n = 0; n < (N - 3); n += 4) {
        i += 7.0 * beta_derivative_integrand(a + n * delta_x, t, theta);
        i += 32.0 * beta_derivative_integrand(a + (n + 1) * delta_x, t, theta);
        i += 12.0 * beta_derivative_integrand(a + (n + 2) * delta_x, t, theta);
        i += 32.0 * beta_derivative_integrand(a + (n + 3) * delta_x, t, theta);
        i += 7.0 * beta_derivative_integrand(a + (n + 4) * delta_x, t, theta);
    }

    //
    const double scale_x = (2.0 / 45.0) * delta_x;
    return scale_x * i;
}


////
// [[Rcpp::export]]
arma::vec lm_fast_cpp(const arma::mat & X, const arma::vec & y) {
    arma::vec beta = arma::solve(X, y, arma::solve_opts::fast);
    return beta;
}

