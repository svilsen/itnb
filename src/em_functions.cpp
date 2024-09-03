#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <roptim.h>
// [[Rcpp::depends(roptim)]]

#include "rdp_functions.hpp"

using namespace roptim;

////
// [[Rcpp::export]]
double complete_loglikelihood(const arma::vec & x, const arma::vec & xi, const arma::vec & z, const double & mu, const double & theta, const double & p, const int & i, const int & t) {
    //
    const int & N = x.size();

    //
    double d = 0.0;
    for (int n = 0; n < N; n++) {
        double log_d = ditnb_cpp(x[n], mu, theta, 0.0, i, t);
        if (p == 0) {
            log_d = (1 - z[n]) * log_d;
        }
        else {
            log_d = z[n] * std::log(p) * xi[n] + (1.0 - z[n]) * (std::log(1.0 - p) + log_d);
        }

        d += log_d;
    }

    return d;
}

// [[Rcpp::export]]
double restricted_loglikelihood(const arma::vec & x, const double & mu, const double & theta, const double & p, const int & i, const int & t) {
    //
    const int & N = x.size();

    //
    double d = 0.0;
    for (int n = 0; n < N; n++) {
        d += ditnb_cpp(x[n], mu, theta, p, i, t);
    }

    return d;
}

////
arma::vec update_z(const arma::vec & x, const arma::vec & xi, const double & mu, const double & theta, const double & p, const int & i, const int & t) {
    //
    const int & N = x.size();

    //
    arma::vec z(N);
    for (int n = 0; n < N; n++) {
        double d_n = std::exp(ditnb_cpp(x[n], mu, theta, p, i, t));
        z[n] = p * xi[n] / (p * xi[n] + (1.0 - p) * d_n);
    }

    return z;
}

double update_p(const arma::vec & xi, const arma::vec & z_star) {
    return 0.0;
}

class RL : public Functor {
public:
    // Parameters
    double p_;

    // Constructor
    RL(const arma::vec & x, const arma::vec & xi, const double & p, const int & i, const int & t) : p_(p), x_(x), xi_(xi), i_(i), t_(t) {

    }

    // Functor / objective function
    double operator()(const arma::vec & x) override {
        double mu_ = x[0];
        double theta_ = x[0];

        double rl = restricted_loglikelihood(x_, mu_, theta_, p_, i_, t_);
        return -rl;
    }

    // Gradient function
    void Gradient(const arma::vec & x, arma::vec & gr) override {

    }

private:
    // Data
    arma::vec x_;
    arma::vec xi_;

    // Parameters
    int i_;
    int t_;
};


////
// [[Rcpp::export]]
Rcpp::List em_itnb(
        const arma::vec & x, const arma::vec & xi,
        const double & mu_0, const double & theta_0, const double & p_0,
        const int & i, const int & t,
        const int & iteration_min, const int & iteration_max, const double & tolerance,
        const int & trace, const bool & save_trace
) {

    //
    arma::vec z_star = update_z(x, xi, mu_0, theta_0, p_0, i, t);
    double c_log_likelihood_old = HUGE_VAL;
    double c_log_likelihood = complete_loglikelihood(x, xi, z_star, mu_0, theta_0, p_0, i, t);

    //
    Roptim<RL> opt("L-BFGS-B");
    RL r_log_likelihood(x, xi, p_0, i, t);

    arma::vec pars_j = {mu_0, theta_0};
    arma::vec lb = 1e-8 * arma::ones(pars_j.size());

    opt.set_lower(lb);
    opt.set_hessian(true);

    //
    arma::vec mu_trace, theta_trace, p_trace, loglike_trace;
    if (save_trace) {
        mu_trace = arma::vec(iteration_max, mu_0);
        theta_trace = arma::vec(iteration_max, theta_0);
        p_trace = arma::vec(iteration_max, p_0);
        loglike_trace = arma::vec(iteration_max, p_0);
    }

    //
    bool not_converged = true;
    std::string convergence_flag = "";

    double mu_j = mu_0;
    double theta_j = theta_0;
    double p_j = p_0;

    //
    int j = 1;
    while (not_converged) {
        //// E-step
        z_star = update_z(x, xi, mu_j, theta_j, p_j, i, t);

        //// M-step
        // Inflation proportion
        p_j = update_p(xi, z_star);
        r_log_likelihood.p_ = p_j;

        // NB parameters
        pars_j = {mu_j, theta_j};
        opt.minimize(r_log_likelihood, pars_j);
        arma::vec pars_opt = opt.par();

        mu_j = pars_opt[0];
        theta_j = pars_opt[1];

        //// Convergence
        c_log_likelihood_old = c_log_likelihood;
        c_log_likelihood = complete_loglikelihood(x, xi, z_star, mu_j, theta_j, p_j, i, t);
        double abs_delta_c_log_likelihood = std::abs(c_log_likelihood - c_log_likelihood_old);

        if (j > iteration_min) {
            if (j < iteration_max) {

                if (abs_delta_c_log_likelihood < tolerance) {
                    not_converged = false;
                }
            }
            else {
                not_converged = false;
                convergence_flag = "MAXIMUM NUMBER OF ITERATIONS EXCEEDED.";
            }
        }

        //// Trace
        if ((trace > 0) & (((j % trace) == 0) | (!not_converged))) {
            Rcpp::Rcout << "Iteration: " <<  j << "\t Current log-likelihood: " << c_log_likelihood << "\t Absolute change in log-likelihood: " << abs_delta_c_log_likelihood << "\n"
                        << "\t Parameters: " << "\t mu = " << mu_j << "\t theta = " << theta_j << "\t p =" << p_j << "\n";
        }

        if (save_trace) {
            mu_trace[j] = mu_j;
            theta_trace[j] = theta_j;
            p_trace[j] = theta_j;
            loglike_trace[j] = c_log_likelihood;
        }

        ////
        j++;
    }

    //
    Rcpp::List res;
    if (save_trace) {
        res = Rcpp::List::create(
            Rcpp::Named("LogLikelihood") = loglike_trace,
            Rcpp::Named("mu") = mu_trace,
            Rcpp::Named("theta") = theta_trace,
            Rcpp::Named("p") = p_trace,
            Rcpp::Named("Converged") = !not_converged,
            Rcpp::Named("Flag") = convergence_flag
        );
    }
    else {
        res = Rcpp::List::create(
            Rcpp::Named("LogLikelihood") = c_log_likelihood,
            Rcpp::Named("mu") = mu_j,
            Rcpp::Named("theta") = theta_j,
            Rcpp::Named("p") = p_j,
            Rcpp::Named("Converged") = !not_converged,
            Rcpp::Named("Flag") = convergence_flag
        );
    }

    return res;
}
