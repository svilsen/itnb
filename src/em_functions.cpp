#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <roptim.h>
// [[Rcpp::depends(roptim)]]

#include "boost/math/special_functions/beta.hpp"
#include "boost/math/special_functions/digamma.hpp"
// [[Rcpp::depends(BH)]]

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
double restricted_loglikelihood(const arma::vec & x, const arma::vec & z, const double & mu, const double & theta, const double & p, const int & i, const int & t) {
    //
    const int & N = x.size();

    //
    double d = 0.0;
    for (int n = 0; n < N; n++) {
        d += (1.0 - z[n]) * ditnb_cpp(x[n], mu, theta, p, i, t);
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
        double d_n = std::exp(ditnb_cpp(x[n], mu, theta, 0.0, i, t));
        z[n] = (p * xi[n]) / (p * xi[n] + (1.0 - p) * d_n);
    }

    return z;
}

double update_p(const arma::vec & xi, const arma::vec & z_star) {
    //
    const int & N = xi.size();

    double num = 0.0;
    double denom = 0.0;
    for (int n = 0; n < N; n++) {
        num += (z_star[n] * xi[n]);
        denom += (1.0 - z_star[n] * (1.0 - xi[n]));
    }

    return num / denom;
}

//[[Rcpp::export]]
double theta_trapz(const double & a, const double & b, const double & theta, const int & t, const int & N) {
    //
    const double b_rev = 1.0 - b;
    const double theta_1 = theta - 1.0;

    //
    const double delta_x = (b - a) / N;
    double i = 0.5 * (std::pow(b, t) * std::pow(b_rev, theta_1) * std::log(b_rev));
    for (int n = 1; n < N; n++) {
        //
        const double x_n = a + n * delta_x;
        const double x_n_rev = 1.0 - x_n;

        //
        i += (std::pow(x_n, t) * std::pow(x_n_rev, theta_1) * std::log(x_n_rev));
    }

    //
    i *= delta_x;
    return i;
}

class RL : public Functor {
public:
    // Parameters
    double p_;
    arma::vec z_;

    // Constructor
    RL(const arma::vec & x, const arma::vec & xi, const arma::vec & z, const double & p, const int & i, const int & t, const double & e) : p_(p), z_(z), x_(x), xi_(xi), N_(x_.size()), i_(i), t_(t), e_(e) { }

    // Functor / objective function
    double operator()(const arma::vec & x) override {
        double mu_ = x[0];
        double theta_ = x[1];

        //
        double d = -restricted_loglikelihood(x_, z_, mu_, theta_, p_, i_, t_);
        return d;
    }

    // Gradient function
    // void Gradient(const arma::vec & x, arma::vec & gr) override {
    //     const double & mu_ = x[0];
    //     const double & theta_ = x[1];
    //
    //     gr[0] = (restricted_loglikelihood(x_, z_, mu_ - e_, theta_, p_, i_, t_) - restricted_loglikelihood(x_, z_, mu_ + e_, theta_, p_, i_, t_)) / (2.0 * e_);
    //     gr[1] = (restricted_loglikelihood(x_, z_, mu_, theta_ - e_, p_, i_, t_) - restricted_loglikelihood(x_, z_, mu_, theta_ + e_, p_, i_, t_)) / (2.0 * e_);
    // }

    void Gradient(const arma::vec & x, arma::vec & gr) override {
        //
        const double & mu_ = x[0];
        const double & theta_ = x[1];

        //
        const double & t_p_t_ = std::pow(theta_, theta_);

        //
        gr[0] = 0.0;
        gr[1] = 0.0;
        for (int n = 0; n < N_; n++) {
            ////
            const double & mt_ = mu_ + theta_;
            const double & mti_ = 1.0 / mt_;
            const double & xt_ = x_[n] + theta_;

            //// Mean parameters
            // Derivatives of the regularised beta function
            const double & b_mu_n = std::pow(mu_, t_) * t_p_t_ / std::pow(mt_, t_ + 1 + theta_);
            const double & b_mu_d = boost::math::beta(t_ + 1, theta_, mu_ * mti_);
            const double & I_mu_ = b_mu_n / b_mu_d;

            // Gradient update
            gr[0] += (1.0 - z_[n]) * (x_[n] / mu_ - xt_ * mti_ - I_mu_);

            //// Overdispersion parameter
            // Derivatives of the regularised beta function
            const double & b_t_f = b_mu_n * mu_ / theta_;

            const double & b_t_i = theta_trapz(0.0, mu_ * mti_, theta_, t_, 100);

            const double & b_t_d = boost::math::digamma(t_ + 1.0 + theta_) - boost::math::digamma(theta_);
            const double & I_t_ = -b_t_f + b_t_i + b_t_d;

            // Gradient update
            gr[1] += (1.0 - z_[n]) * (1.0 + std::log(theta_) - xt_ * mti_ + boost::math::digamma(xt_) - boost::math::digamma(theta_) - I_t_);
        }
    }

private:
    // Data
    arma::vec x_;
    arma::vec xi_;

    int N_;

    // Parameters
    int i_;
    int t_;
    double e_;
};


////
// [[Rcpp::export]]
Rcpp::List em_itnb_cpp(
        const arma::vec & x, const arma::vec & xi,
        const double & mu_0, const double & theta_0, const double & p_0,
        const int & i, const int & t,
        const int & iteration_min, const int & iteration_max, const double & tolerance,
        const int & trace, const bool & save_trace
) {
    // Cube root of 10^(-16)
    const double e = 4.6416 * std::pow(10, -6);

    //
    arma::vec z_star = update_z(x, xi, mu_0, theta_0, p_0, i, t);
    double c_log_likelihood_old = HUGE_VAL;
    double c_log_likelihood = complete_loglikelihood(x, xi, z_star, mu_0, theta_0, p_0, i, t);

    //
    Roptim<RL> opt("L-BFGS-B");
    RL r_log_likelihood(x, xi, z_star, p_0, i, t, e);

    arma::vec pars_j = {mu_0, theta_0};

    arma::vec lb = arma::ones(pars_j.size());
    lb[0] = t + e;
    lb[1] = 1e-8;

    opt.set_lower(lb);
    opt.set_hessian(true);

    //
    arma::vec mu_trace, theta_trace, p_trace, loglike_trace;
    if (save_trace) {
        mu_trace.set_size(iteration_max);
        mu_trace[0] = mu_0;

        theta_trace.set_size(iteration_max);
        theta_trace[0] = theta_0;

        p_trace.set_size(iteration_max);
        p_trace[0] = p_0;

        loglike_trace.set_size(iteration_max);
        loglike_trace[0] = c_log_likelihood;
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
        r_log_likelihood.z_ = z_star;

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
        if (trace > 0) {
            if (((j % trace) == 0) | (!not_converged)) {
                Rcpp::Rcout << "Iteration: " <<  j << "\t Current log-likelihood: " << c_log_likelihood << "\t Absolute change in log-likelihood: " << abs_delta_c_log_likelihood << "\n"
                            << "\t Parameters: " << "\t mu = " << mu_j << "\t theta = " << theta_j << "\t p =" << p_j << "\n";
            }
        }

        if (save_trace) {
            mu_trace[j] = mu_j;
            theta_trace[j] = theta_j;
            p_trace[j] = p_j;
            loglike_trace[j] = c_log_likelihood;
        }

        ////
        j++;
    }

    //
    Rcpp::List trace_list;
    if (save_trace) {
        trace_list = Rcpp::List::create(
            Rcpp::Named("LogLikelihood") = loglike_trace.head(j),
            Rcpp::Named("mu") = mu_trace.head(j),
            Rcpp::Named("theta") = theta_trace.head(j),
            Rcpp::Named("p") = p_trace.head(j)
        );
    }

    return Rcpp::List::create(
        Rcpp::Named("LogLikelihood") = c_log_likelihood,
        Rcpp::Named("mu") = mu_j,
        Rcpp::Named("theta") = theta_j,
        Rcpp::Named("p") = p_j,
        Rcpp::Named("trace") = trace_list,
        Rcpp::Named("Converged") = !not_converged,
        Rcpp::Named("Flag") = convergence_flag
    );
}
