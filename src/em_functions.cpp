#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <roptim.h>
// [[Rcpp::depends(roptim)]]

#include "rdp_functions.hpp"

using namespace roptim;

////
// [[Rcpp::export]]
double complete_loglikelihood(const arma::vec & x, const arma::vec & xi, const arma::vec & z, const arma::vec & mu, const arma::vec & theta, const arma::vec & p, const arma::vec & i, const arma::vec & t) {
    //
    const int & N_x = x.size();

    const int & N_mu = mu.size();
    double mu_n = mu[0];

    const int & N_theta = theta.size();
    double theta_n = theta[0];

    const int & N_p = p.size();
    double p_n = p[0];

    const int & N_i = i.size();
    int i_n = i[0];

    const int & N_t = t.size();
    int t_n = t[0];

    //
    double d = 0.0;
    for (int n = 0; n < N_x; n++) {
        if (N_mu > 1) {
            mu_n = mu[n];
        }
        if (N_theta > 1) {
            theta_n = theta[n];
        }
        if (N_p > 1) {
            p_n = p[n];
        }
        if (N_i > 1) {
            i_n = i[n];
        }
        if (N_t > 1) {
            t_n = t[n];
        }

        double log_d = ditnb_cpp(x[n], mu_n, theta_n, 0.0, i_n, t_n);
        if (p_n == 0) {
            log_d = (1 - z[n]) * log_d;
        }
        else {
            log_d = z[n] * std::log(p_n) * xi[n] + (1.0 - z[n]) * (std::log(1.0 - p_n) + log_d);
        }

        d += log_d;
    }

    return d;
}

// [[Rcpp::export]]
double restricted_loglikelihood(const arma::vec & x, const arma::vec & mu, const arma::vec & theta, const arma::vec & p, const arma::vec & i, const arma::vec & t) {
    //
    const int & N_x = x.size();

    const int & N_mu = mu.size();
    double mu_n = mu[0];

    const int & N_theta = theta.size();
    double theta_n = theta[0];

    const int & N_p = p.size();
    double p_n = p[0];

    const int & N_i = i.size();
    int i_n = i[0];

    const int & N_t = t.size();
    int t_n = t[0];

    //
    double d = 0.0;
    for (int n = 0; n < N_x; n++) {
        if (N_mu > 1) {
            mu_n = mu[n];
        }
        if (N_theta > 1) {
            theta_n = theta[n];
        }
        if (N_p > 1) {
            p_n = p[n];
        }
        if (N_i > 1) {
            i_n = i[n];
        }
        if (N_t > 1) {
            t_n = t[n];
        }

        d += ditnb_cpp(x[n], mu_n, theta_n, p_n, i_n, t_n);
    }

    return d;
}

////
arma::vec update_z(const arma::vec & x, const arma::vec & xi, const arma::vec & mu, const arma::vec & theta, const arma::vec & p, const arma::vec & i, const arma::vec & t) {
    //
    const int & N_x = x.size();

    const int & N_mu = mu.size();
    double mu_n = mu[0];

    const int & N_theta = theta.size();
    double theta_n = theta[0];

    const int & N_p = p.size();
    double p_n = p[0];

    const int & N_i = i.size();
    int i_n = i[0];

    const int & N_t = t.size();
    int t_n = t[0];

    //
    arma::vec z(N_x);
    for (int n = 0; n < N_x; n++) {
        if (N_mu > 1) {
            mu_n = mu[n];
        }
        if (N_theta > 1) {
            theta_n = theta[n];
        }
        if (N_p > 1) {
            p_n = p[n];
        }
        if (N_i > 1) {
            i_n = i[n];
        }
        if (N_t > 1) {
            t_n = t[n];
        }

        double d_n = std::exp(ditnb_cpp(x[n], mu_n, theta_n, p_n, i_n, t_n));
        z[n] = p_n * xi[n] / (p_n * xi[n] + (1.0 - p_n) * d_n);
    }

    return z;
}

class MSTEP : public Functor {
public:
    MSTEP(const arma::vec & x_, const arma::vec & xi_, const arma::vec & p_, const arma::vec & i_, const arma::vec & t_, const int & N_mu_, const int & N_theta_) : x(x_), xi(xi_), p(p_), i(i_), t(t_) {

    }

    double operator()(const arma::vec &x) override {
        return 0.0;
    }

    void Gradient(const arma::vec & x, arma::vec & gr) override {

    }

private:
    arma::vec x;
    arma::vec xi;
    arma::vec p;
    arma::vec i;
    arma::vec t;

    int N_mu;
    int N_theta;
};


////
// [[Rcpp::export]]
Rcpp::List em_itnb(const arma::vec & x, const arma::vec & xi, const arma::vec & mu, const arma::vec & theta, const arma::vec & p, const arma::vec & i, const arma::vec & t) {

    //
    arma::vec z_star = update_z(x, xi, mu, theta, p, i, t);
    double c_log_likelihood = complete_loglikelihood(x, xi, z_star, mu, theta, p, i, t);

    //
    MSTEP restricted(x, xi, p, i, t, mu.size(), theta.size());
    Roptim<MSTEP> opt("L-BFGS-B");

    arma::vec pars = arma::join_cols(mu, theta);
    arma::vec lb = 1e-8 * arma::ones(pars.size());

    opt.set_lower(lb);
    opt.set_hessian(true);
    opt.minimize(restricted, pars);

    return Rcpp::List::create(
        Rcpp::Named("LogLikelihood") = 0.0,
        Rcpp::Named("mu") = 0.0,
        Rcpp::Named("theta") = 0.0,
        Rcpp::Named("p") = 0.0,
        Rcpp::Named("Converged") = 0.0
    );
}
