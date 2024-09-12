#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <roptim.h>
// [[Rcpp::depends(roptim)]]

#include "rdp_functions.hpp"

using namespace roptim;

////
double complete_loglikelihood(const arma::vec & x, const arma::vec & xi, const arma::vec & z, const double & mu, const double & theta, const double & p, const int & i, const int & t) {
    //
    const int & N = x.size();

    //
    double d = 0.0;
    for (int n = 0; n < N; n++) {
        double log_d = ditnb_cpp(x[n], mu, theta, 0.0, i, t);
        d += z[n] * std::log(p) * xi[n] + (1.0 - z[n]) * (std::log(1.0 - p) + log_d);
    }

    return d;
}

//// E-step
void update_z(arma::vec & z, const arma::vec & x, const arma::vec & xi, const double & mu, const double & theta, const double & p, const int & i, const int & t, const int & N) {
    //
    for (int n = 0; n < N; n++) {
        double d_n = std::exp(ditnb_cpp(x[n], mu, theta, 0.0, i, t));
        z[n] = (p * xi[n]) / (p * xi[n] + (1.0 - p) * d_n);
    }
}

//// M-step
// Inflation proportion
void update_p(double & p, const arma::vec & xi, const arma::vec & z, const int & N) {
    double num = 0.0;
    double denom = 0.0;
    for (int n = 0; n < N; n++) {
        num += (z[n] * xi[n]);
        denom += (1.0 - z[n] * (1.0 - xi[n]));
    }

    p = num / denom;

    if (p < 1e-16) {
        p = 4e-16;
    }
}

// NB-parameters
double beta_derivative_integrand(const double & s, const double & t, const double & theta) {
    const double s_rev = 1.0 - s;
    double r = 0.0;
    if ((s > 1e-8) & (s < (1.0 - 1e-8))) {
        r = std::pow(s, t) * std::pow(s_rev, theta - 1) * std::log(s_rev);
    }

    return r;
}

double beta_derivative_trapz(const double & a, const double & b, const double & theta, const int & t, const int & N) {
    //
    const double delta_x = (b - a) / N;

    double i = 0.5 * (beta_derivative_integrand(a, t, theta) + beta_derivative_integrand(b, t, theta));

    for (int n = 1; n < N; n++) {
        i += beta_derivative_integrand(a + n * delta_x, t, theta);
    }

    return i * delta_x;
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

class RL : public Functor {
public:
    // Parameters
    double p_;
    arma::vec z_;

    // Constructor
    RL(const arma::vec & x, const arma::vec & xi, const arma::vec & z, const double & p, const int & i, const int & t, const double & e) : p_(p), z_(z), x_(x), xi_(xi), N_(x.size()), i_(i), t_(t), e_(e) { }

    // Functor / objective function
    double operator()(const arma::vec & x) override {
        const double & mu_ = x[0];
        const double & theta_ = x[1];

        //
        double d = 0.0;
        for (int n = 0; n < N_; n++) {
            d += (1.0 - z_[n]) * ditnb_cpp(x_[n], mu_, theta_, p_, i_, t_);
        }

        return -d;
    }

    // Gradient function
    void Gradient(const arma::vec & x, arma::vec & gr) override {
        //
        const double & mu_ = x[0];
        const double & theta_ = x[1];

        //
        double l = 0.0;
        double l_mu_p = 0.0;
        double l_theta_p = 0.0;
        for (int n = 0; n < N_; n++) {
            l += (1.0 - z_[n]) * ditnb_cpp(x_[n], mu_, theta_, p_, i_, t_);

            //// Location
            l_mu_p += (1.0 - z_[n]) * ditnb_cpp(x_[n], mu_ + e_, theta_, p_, i_, t_);

            //// Overdispersion
            l_theta_p += (1.0 - z_[n]) * ditnb_cpp(x_[n], mu_, theta_ + e_, p_, i_, t_);
        }

        //
        const double ei = 1.0 / e_;
        gr[0] = (l - l_mu_p) * ei;
        gr[1] = (l - l_theta_p) * ei;

        //     gr[0] = 0.0;
        //     gr[1] = 0.0;
        //     for (int n = 0; n < N_; n++) {
        //         ////
        //         const double & mu_theta_ = mu_ + theta_;
        //         const double & mu_theta_i_ = 1.0 / mu_theta_;
        //
        //         const double & x_theta_ = x_[n] + theta_;
        //         const double b_ = std::exp(R::pbeta(mu_ * mu_theta_i_, t_ + 1, theta_, true, true) + R::lbeta(t_ + 1, theta_));
        //
        //         //// Location
        //         const double & b_mu_n_ = std::pow(mu_, t_) * std::pow(theta_, theta_) / std::pow(mu_theta_, t_ + 1 + theta_);
        //         double I_mu_ = b_mu_n_ / b_;
        //
        //         gr[0] += (z_[n] - 1.0) * (x_[n] / mu_ - x_theta_ * mu_theta_i_ - I_mu_);
        //
        //         //// Overdispersion
        //         // Derivatives of the regularised beta function
        //         const double & b_t_f_ = std::pow(mu_, t_ + 1.0) * std::pow(theta_, theta_ - 1.0) / std::pow(mu_theta_, t_ + 1 + theta_);
        //         const double & b_t_i_ = beta_derivative_boole(0.0, mu_ * mu_theta_i_, theta_, t_, M_);
        //         const double & b_t_d_ = R::digamma(theta_) - R::digamma(t_ + 1.0 + theta_);
        //
        //         double I_theta_ = (b_t_i_ - b_t_f_) / b_ - b_t_d_;
        //         gr[1] += (z_[n] - 1.0) * (1.0 + std::log(theta_) - std::log(mu_theta_) - x_theta_ * mu_theta_i_ + R::digamma(x_theta_) - R::digamma(theta_) - I_theta_);
        //     }
        // }
    }

private:
    // Data
    arma::vec x_;
    arma::vec xi_;

    // Number of observations
    int N_;

    // Parameters
    int i_;
    int t_;

    // Precision
    double e_;
};

////
// [[Rcpp::export]]
Rcpp::List em_itnb_cpp(
        const arma::vec & x, const arma::vec & xi,
        const double & mu_0, const double & theta_0, const double & p_0,
        const int & i, const int & t,
        const int & iteration_min, const int & iteration_max,
        const double & tolerance,
        const int & trace, const bool & save_trace
) {
    //
    const int & N = x.size();
    const double e = 4.6416 * 1e-6; // std::pow(2e-16, 0.33);

    //
    double p_j = p_0;
    if (p_j < 1e-16) {
        p_j = 4e-16;
    }

    arma::vec z(N);
    update_z(z, x, xi, mu_0, theta_0, p_0, i, t, N);

    double c_log_likelihood_old = HUGE_VAL;
    double c_log_likelihood = complete_loglikelihood(x, xi, z, mu_0, theta_0, p_j, i, t);
    double abs_delta_c_log_likelihood = std::abs(c_log_likelihood - c_log_likelihood_old);

    //
    Roptim<RL> opt("L-BFGS-B");
    RL r_log_likelihood(x, xi, z, p_0, i, t, e);

    //
    arma::vec lb = arma::ones(2);
    lb[0] = t + 1;
    lb[1] = 1e-8;

    const double x_max = x.max();
    arma::vec ub = arma::ones(2);
    ub[0] = x_max;
    ub[1] = x_max * x_max;

    opt.set_lower(lb);
    opt.set_upper(ub);

    //
    arma::vec pars_j = {mu_0, theta_0};
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

    if (trace > 0) {
        Rcpp::Rcout << "Iteration: " <<  0 << "\t Current log-likelihood: " << c_log_likelihood << "\t Absolute change in log-likelihood: " << abs_delta_c_log_likelihood << "\n"
                    << "\t Parameters: " << "\t mu = " << mu_j << "\t theta = " << theta_j << "\t p =" << p_j << "\n";
    }

    //
    int j = 1;
    while (not_converged) {
        //// E-step
        update_z(z, x, xi, mu_j, theta_j, p_j, i, t, N);

        //// M-step
        // Inflation proportion
        update_p(p_j, xi, z, N);

        // NB parameters
        r_log_likelihood.p_ = p_j;
        r_log_likelihood.z_ = z;

        pars_j = {mu_j, theta_j};

        opt.minimize(r_log_likelihood, pars_j);
        arma::vec pars_opt = opt.par();

        mu_j = pars_opt[0];
        theta_j = pars_opt[1];

        //// Convergence
        c_log_likelihood_old = c_log_likelihood;
        c_log_likelihood = complete_loglikelihood(x, xi, z, mu_j, theta_j, p_j, i, t);
        abs_delta_c_log_likelihood = std::abs(c_log_likelihood - c_log_likelihood_old);

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
