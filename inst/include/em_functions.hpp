#ifndef em_functions
#define em_functions

#include <RcppArmadillo.h>

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <roptim.h>
// [[Rcpp::depends(roptim)]]

#include "rdp_functions.hpp"

using namespace roptim;

////
double loglikelihood(const arma::vec & x, const double & mu, const double & theta, const double & p, const int & i, const int & t, const int & N);

//// E-step
void update_z(arma::vec & z, const arma::vec & x, const arma::vec & xi, const double & mu, const double & theta, const double & p, const int & i, const int & t, const int & N);

//// M-step
// Inflation proportion
void update_p(double & p, const arma::vec & xi, const arma::vec & z, const int & N);

// NB-parameters
double beta_derivative_integrand(const double & s, const double & t, const double & theta);

double beta_derivative_boole(const double & a, const double & b, const double & theta, const int & t, const int & N);

class RL : public Functor {
public:
    // Parameters
    arma::vec z_;

    // Constructor
    RL(const arma::vec & x, const arma::vec & xi, const arma::vec & z, const int & i, const int & t, const int & steps, const bool & fd, const double & steps_fd) : z_(z), x_(x), xi_(xi), N_(x.size()), i_(i), t_(t), steps_(steps), fd_(fd), steps_fd_(steps_fd)  { }

    // Functor / objective function
    double operator()(const arma::vec & x) override {
        const double & mu_ = x[0];
        const double & theta_ = x[1];

        //
        double d = 0.0;
        for (int n = 0; n < N_; n++) {
            d += (1.0 - z_[n]) * ditnb_cpp(x_[n], mu_, theta_, 0.0, i_, t_);
        }

        return -d;
    }

    // Gradient function
    void Gradient(const arma::vec & x, arma::vec & gr) override {
        //
        const double & mu_ = x[0];
        const double & theta_ = x[1];

        //
        if (fd_) {
            double l = 0.0;
            double l_mu_p = 0.0;
            double l_theta_p = 0.0;
            for (int n = 0; n < N_; n++) {
                l += (1.0 - z_[n]) * ditnb_cpp(x_[n], mu_, theta_, 0.0, i_, t_);

                //// Location
                l_mu_p += (1.0 - z_[n]) * ditnb_cpp(x_[n], mu_ + steps_fd_, theta_, 0.0, i_, t_);

                //// Overdispersion
                l_theta_p += (1.0 - z_[n]) * ditnb_cpp(x_[n], mu_, theta_ + steps_fd_, 0.0, i_, t_);
            }

            //
            const double ei = 1.0 / steps_fd_;
            gr[0] = (l  - l_mu_p) * ei;
            gr[1] = (l - l_theta_p) * ei;
        }
        else {
            gr[0] = 0.0;
            gr[1] = 0.0;
            for (int n = 0; n < N_; n++) {
                ////
                const double & mu_theta_ = mu_ + theta_;
                const double & mu_theta_i_ = 1.0 / mu_theta_;

                const double & x_theta_ = x_[n] + theta_;
                const double b_ = std::exp(R::pbeta(mu_ * mu_theta_i_, t_ + 1, theta_, true, true) + R::lbeta(t_ + 1, theta_));

                //// Location
                const double & b_mu_n_ = std::pow(mu_, t_) * std::pow(theta_, theta_) / std::pow(mu_theta_, t_ + 1 + theta_);
                double I_mu_ = b_mu_n_ / b_;

                gr[0] += (z_[n] - 1.0) * (x_[n] / mu_ - x_theta_ * mu_theta_i_ - I_mu_);

                //// Overdispersion
                // Derivatives of the regularised beta function
                const double & b_t_f_ = std::pow(mu_, t_ + 1.0) * std::pow(theta_, theta_ - 1.0) / std::pow(mu_theta_, t_ + 1 + theta_);
                const double & b_t_i_ = beta_derivative_boole(0.0, mu_ * mu_theta_i_, theta_, t_, steps_);
                const double & b_t_d_ = R::digamma(theta_) - R::digamma(t_ + 1.0 + theta_);

                double I_theta_ = (b_t_i_ - b_t_f_) / b_ - b_t_d_;
                gr[1] += (z_[n] - 1.0) * (1.0 + std::log(theta_) - std::log(mu_theta_) - x_theta_ * mu_theta_i_ + R::digamma(x_theta_) - R::digamma(theta_) - I_theta_);
            }
        }
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

    // Integral precision
    int steps_;

    // Finite difference
    bool fd_;
    double steps_fd_;
};

////
void optimise_itnb(
        double mu_j, double theta_j, double p_j, double loglike_j,
        int j, bool not_converged, std::string convergence_flag,
        arma::vec mu_trace, arma::vec theta_trace, arma::vec p_trace, arma::vec loglike_trace,
        const arma::vec & x, const arma::vec xi,
        const int & i, const int & t, const int & N,
        const int & iteration_min, const int & iteration_max,
        const double & tolerance,
        const int & steps,
        const bool & fd, const double & steps_fd,
        const int & trace, const bool & save_trace
);

//
Rcpp::List em_itnb_cpp(
        const arma::vec & x, const arma::vec & xi,
        const double & mu_0, const double & theta_0, const double & p_0,
        const int & i, const int & t,
        const int & iteration_min, const int & iteration_max,
        const double & tolerance,
        const int & steps,
        const bool & fd, const double & steps_fd,
        const int & trace, const bool & save_trace
);


#endif //em_functions
