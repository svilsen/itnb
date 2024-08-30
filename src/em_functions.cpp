#include <Rcpp.h>

#include "rdp_functions.hpp"

////
// [[Rcpp::export]]
double complete_loglikelihood(const std::vector<int> & x, const std::vector<int> & xi, const std::vector<double> & z, const std::vector<double> & mu, const std::vector<double> & theta, const std::vector<double> & p, std::vector<int> & i, std::vector<int> & t) {
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
double restricted_loglikelihood(const std::vector<int> & x, const std::vector<double> & mu, const std::vector<double> & theta, const std::vector<double> & p, std::vector<int> & i, std::vector<int> & t) {
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
std::vector<double> update_z(const std::vector<int> & x, const std::vector<int> & xi, const std::vector<double> & mu, const std::vector<double> & theta, const std::vector<double> & p, std::vector<int> & i, std::vector<int> & t) {
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
    std::vector<double> z(N_x);
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

////
// [[Rcpp::export]]
Rcpp::List em_itnb(const std::vector<int> & x, const std::vector<int> & xi, const std::vector<double> & mu, const std::vector<double> & theta, const std::vector<double> & p, std::vector<int> & i, std::vector<int> & t) {

    std::vector<double> z_star = update_z(x, xi, mu, theta, p, i, t);
    double c_log_likelihood = complete_loglikelihood(x, xi, z_star, mu, theta, p, i, t);

    return Rcpp::List::create(
        Rcpp::Named("LogLikelihood") = 0.0,
        Rcpp::Named("mu") = 0.0,
        Rcpp::Named("theta") = 0.0,
        Rcpp::Named("p") = 0.0,
        Rcpp::Named("Converged") = 0.0
    );
}
