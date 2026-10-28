#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <roptim.h>
// [[Rcpp::depends(roptim)]]

#include "rdp.hpp"
#include "links.hpp"
#include "aux_functions.hpp"

using namespace roptim;

// Truncated negative binomial parameters
class EM_ITPOIS : public Functor {
public:
    // Parameters
    arma::vec z;

    // Constructor
    EM_ITPOIS(const arma::mat & X_, const arma::vec & y_, const arma::vec & yi_, const arma::vec & z_, const int & i_, const int & t_, const int & steps_, const bool & exact_, Link & LO_) :
        z(z_), X(X_), y(y_), yi(yi_), N(y.size()), M(X.n_cols), i(i_), t(t_), steps(steps_), exact(exact_), LO(LO_) { }

    // Functor / objective function
    double operator()(const arma::vec & par) override {
        const arma::vec & beta = par;

        //
        double d = 0.0;
        for (int n = 0; n < N; n++) {
            //
            const arma::rowvec & x_n = X.row(n);
            const double & eta_n = arma::dot(x_n, beta);
            const double & mu_n = LO.link_inv(eta_n);

            d -= (1.0 - z[n]) * ditpois_cpp(y[n], mu_n, 0.0, i, t);
        }

        if (std::isinf(d) || std::isnan(d)) {
            d = 1e16;
        }

        return d;
    }

    // Gradient function
    void Gradient(const arma::vec & par, arma::vec & gr) override {
        //
        gr = arma::zeros(M);
        if (exact) {
            const arma::vec & beta = par;
            for (int n = 0; n < N; n++) {
                //
                const arma::rowvec & x_n = X.row(n);
                const double & eta_n = arma::dot(x_n, beta);
                const double & mu_n = LO.link_inv(eta_n);

                //
                double gr_beta = (z[n] - 1.0) * (y[n] / mu_n - 1.0);
                if (t > -1) {
                    //
                    const double & b = R::ppois(t, mu_n, 0, 0);
                    const double & b_mu_n = R::dpois(t, mu_n, 0);
                    const double & I_mu = b_mu_n / b;

                    //
                    gr_beta += (z[n] - 1.0) * (-I_mu);
                }

                //
                for (int m = 0; m < M; m++) {
                    const double & gr_mu = LO.link_inv_dev(mu_n, x_n, m);
                    gr[m] += gr_beta * gr_mu;
                }
            }
        }
        else{
            ApproximateGradient(par, gr);
        }
    }

private:
    // Data
    arma::mat X;
    arma::vec y;
    arma::vec yi;

    // Number of observations
    int N;
    int M;

    // Parameters
    int i;
    int t;

    // Integral precision
    int steps;

    //
    bool exact;

    // Link
    Link LO;
};

//// Optimiser
void optimise_itpois(
        arma::vec & beta_j, double & p_j, double & loglike_j,
        int & j, bool & not_converged, std::string & convergence_flag,
        arma::mat & beta_trace, arma::vec & alpha_trace, arma::vec & p_trace,
        arma::vec & loglike_trace, arma::mat & approx_hessian,
        const arma::mat & X, const arma::vec & y, const arma::vec yi,
        const int & i, const int & t, const std::string & link,
        const int & N, const int & M,
        const int & iteration_min, const int & iteration_max,
        const double & tolerance,
        const int & steps, const bool & exact,
        const int & trace, const bool & save_trace
) {
    //
    not_converged = true;
    convergence_flag = "";

    //
    Link LO(link);

    //
    arma::vec z = arma::zeros(N);
    update_z(z, X, y, yi, beta_j, HUGE_VAL, p_j, i, t, N, LO);

    //
    arma::vec pars_j = arma::vec(M);
    arma::vec pars_j_old = pars_j;

    pars_j = beta_j;

    //
    Roptim<EM_ITPOIS> opt("L-BFGS-B");
    EM_ITPOIS log_likelihood(X, y, yi, z, i, t, steps, exact, LO);

    //
    loglike_j = opt.value();
    double loglike_j_old = HUGE_VAL;
    double delta_loglike_j = loglike_j - loglike_j_old;

    if (trace > 0) {
        Rcpp::Rcout << "Iteration: " << j << "\t Current log-likelihood: " << loglike_j << "\t Change in log-likelihood: " << delta_loglike_j << "\n"
                    << "\t Parameters: " << "\t beta = " << beta_j.t() << "\t alpha = " << 0.0 << "\t p = " << p_j << "\n";
    }

    //
    j = 1;
    while (not_converged) {
        ////
        //
        const double p_j_old = p_j;
        pars_j_old = beta_j;

        //// E-step
        update_z(z, X, y, yi, beta_j, HUGE_VAL, p_j, i, t, N, LO);

        //// M-step
        // Inflation proportion
        update_p(p_j, yi, z, N);

        // NB parameters
        log_likelihood.z = z;

        opt.minimize(log_likelihood, pars_j_old);
        pars_j = opt.par();
        beta_j = pars_j;

        //// Convergence
        loglike_j_old = loglike_j;
        loglike_j = opt.value();
        delta_loglike_j = loglike_j - loglike_j_old;

        if (delta_loglike_j < 0.0) {
            not_converged = false;
            convergence_flag = "DECREASING LIKELIHOOD DETECTED.";

            p_j = p_j_old;
            beta_j = pars_j_old;

            loglike_j = loglike_j_old;
        }
        else if (j > iteration_min) {
            if (j < iteration_max) {
                if (std::abs(delta_loglike_j) < tolerance) {
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
                Rcpp::Rcout << "Iteration: " <<  j << "\t Current log-likelihood: " << loglike_j << "\t Absolute change in log-likelihood: " << delta_loglike_j << "\n"
                            << "\t Parameters: " << "\t beta = " << beta_j.t() << "\t alpha = " << 0.0 << "\t p = " << p_j << "\n";
            }
        }

        if (save_trace) {
            beta_trace.row(j) = beta_j.t();
            alpha_trace[j] = 0.0;
            p_trace[j] = p_j;
            loglike_trace[j] = loglike_j;
        }

        ////
        j++;
    }

    log_likelihood.ApproximateHessian(opt.par(), approx_hessian);
}

//// R interface
// [[Rcpp::export]]
Rcpp::List em_itpois_cpp(
        const arma::mat & X, const arma::vec & y, const arma::vec & yi,
        const arma::vec & beta_0, const double & p_0,
        const int & i, const int & t, const std::string & link,
        const int & iteration_min, const int & iteration_max,
        const double & tolerance,
        const int & steps, const bool & exact,
        const int & trace, const bool & save_trace
) {
    //
    const int & N = y.size();
    const int & M = X.n_cols;

    //
    int j = 0;
    bool not_converged;
    std::string convergence_flag = "";

    //
    double p_j = p_0;

    //
    arma::vec beta_j = beta_0;

    //
    double loglike_j = 0.0;

    //
    arma::mat beta_trace;
    arma::vec alpha_trace, p_trace, loglike_trace;
    if (save_trace) {
        beta_trace.set_size(iteration_max + 1, M);
        beta_trace.row(0) = beta_0.t();

        alpha_trace.set_size(iteration_max + 1);
        alpha_trace[0] = 0.0;

        p_trace.set_size(iteration_max + 1);
        p_trace[0] = p_0;

        loglike_trace.set_size(iteration_max + 1);
        loglike_trace[0] = HUGE_VAL;
    }

    //
    arma::mat approx_hessian = arma::mat(beta_j.size(), beta_j.size());

    //
    optimise_itpois(
        beta_j, p_j, loglike_j,
        j, not_converged, convergence_flag,
        beta_trace, alpha_trace, p_trace,
        loglike_trace, approx_hessian,
        X, y, yi, i, t, link, N, M,
        iteration_min, iteration_max,
        tolerance,
        steps, exact,
        trace, save_trace
    );

    //
    Rcpp::List trace_list;
    if (save_trace) {
        trace_list = Rcpp::List::create(
            Rcpp::Named("LogLikelihood") = loglike_trace.head(j),
            Rcpp::Named("beta") = beta_trace.rows(0, j - 1),
            Rcpp::Named("alpha") = alpha_trace.head(j),
            Rcpp::Named("p") = p_trace.head(j)
        );
    }

    // HOT FIX
    approx_hessian += 1e-8 * arma::eye(beta_j.size(), beta_j.size());

    //
    arma::mat vcov = arma::inv(approx_hessian);

    //
    return Rcpp::List::create(
        Rcpp::Named("formula") = 0,
        Rcpp::Named("data") = 0,
        Rcpp::Named("i") = i,
        Rcpp::Named("t") = t,
        Rcpp::Named("link") = link,
        Rcpp::Named("loglikelihood") = -loglike_j,
        Rcpp::Named("beta") = beta_j,
        Rcpp::Named("alpha") = 0.0,
        Rcpp::Named("p") = p_j,
        Rcpp::Named("vcov") = vcov,
        Rcpp::Named("logtheta") = 0,
        Rcpp::Named("trace") = trace_list,
        Rcpp::Named("converged") = !not_converged,
        Rcpp::Named("iterations") = j,
        Rcpp::Named("flag") = convergence_flag
    );
}
