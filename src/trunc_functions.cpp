#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <roptim.h>
// [[Rcpp::depends(roptim)]]

#include "rdp_functions.hpp"
#include "aux_functions.hpp"

using namespace roptim;

// Negative binomial
class MLE: public Functor {
public:
    double lambda;

    // Constructor
    MLE(const arma::mat & X_, const arma::vec & y_, const int & i_, const int & t_, const int & steps_, const bool & exact_, const double & lambda_) : lambda(lambda_), X(X_), y(y_), N(y.size()), M(X.n_cols), i(i_), t(t_), steps(steps_), exact(exact_) { }

    // Functor / objective function
    double operator()(const arma::vec & par) override {
        const arma::vec & beta = par.head(M);
        const double & theta = par[M];

        //
        double d = 0.0;
        double penalisation = 0.0;
        for (int n = 0; n < N; n++) {
            //
            const arma::rowvec & x_n = X.row(n);
            const arma::vec mu_n = x_n * beta;

            //
            double mu = mu_n[0];
            if (mu < 1e-6) {
                mu = 1e-6;
            }

            const double e = (mu_n[0] - (t + 1));
            double g = 0.0;
            if (e < 0) {
                g = e * e;
            }

            const double d_n = ditnb_cpp(y[n], mu, theta, 0.0, i, t);

            d += d_n;
            penalisation += lambda * g;
        }

        return -d + penalisation;
    }

    // Gradient function
    void Gradient(const arma::vec & par, arma::vec & gr) override {
        //
        gr = arma::zeros(M + 1);
        if (exact) {
            const arma::vec & beta = par.head(M);
            const double & theta = par[M];

            for (int n = 0; n < N; n++) {
                //
                const arma::rowvec & x_n = X.row(n);
                const arma::vec mu_n = x_n * beta;

                double mu = mu_n[0];
                if (mu < 1e-6) {
                    mu = 1e-6;
                }

                //
                const double & mu_theta = mu + theta;
                const double & mu_theta_i = 1.0 / mu_theta;

                const double & y_theta = y[n] + theta;

                //
                double gr_beta = (y[n] / mu - y_theta * mu_theta_i);
                double gr_theta = (1.0 + std::log(theta) - std::log(mu_theta) - y_theta * mu_theta_i + R::digamma(y_theta) - R::digamma(theta));
                if (t > -1) {
                    const double & b = std::exp(R::pbeta(mu * mu_theta_i, t + 1, theta, true, true) + R::lbeta(t + 1, theta));

                    // Location
                    const double & b_mu_n = std::pow(mu, t) * std::pow(theta, theta) / std::pow(mu_theta, t + 1 + theta);
                    const double & I_mu = b_mu_n / b;

                    gr_beta += (-I_mu);

                    // Overdispersion
                    const double & b_t_f = std::pow(mu, t + 1.0) * std::pow(theta, theta - 1.0) / std::pow(mu_theta, t + 1 + theta);
                    const double & b_t_i = beta_derivative_rectangle(0.0, mu * mu_theta_i, theta, t, steps);
                    const double & b_t_d = R::digamma(theta) - R::digamma(t + 1.0 + theta);

                    const double & I_theta = (b_t_i - b_t_f) / b - b_t_d;
                    gr_theta += (-I_theta);
                }

                //
                for (int m = 0; m < M; m++) {
                    double gr_g = 0.0;
                    if (mu_n[0] < (t + 1)) {
                        gr_g = 2.0 * x_n[m] * (mu_n[0] - (t + 1));
                    }

                    gr[m] += (-gr_beta * x_n[m]) + lambda * gr_g;
                }

                //
                gr[M] += (-gr_theta);
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
};

//// Optimiser
void optimise_tnb(
        arma::vec & beta_j, double & theta_j, const double & p_0, double & loglike_j,
        int & j, bool & not_converged, std::string & convergence_flag,
        const arma::mat & X, const arma::vec & y,
        const int & i, const int & t, const int & N, const int & M,
        const double & tolerance, const arma::vec & lambda,
        const int & steps, const bool & exact,
        const int & trace
) {
    not_converged = true;
    convergence_flag = "";

    //
    arma::vec pars_j = arma::vec(M + 1);
    pars_j.head(M) = beta_j;
    pars_j[M] = theta_j;

    //
    Roptim<MLE> opt("L-BFGS-B");
    MLE r_log_likelihood(X, y, i, t, steps, exact, lambda[0]);

    //
    arma::vec lb = (-HUGE_VAL) * arma::ones(M + 1);
    lb[M] = 1e-8;

    opt.set_lower(lb);
    opt.control.trace = trace;

    //
    opt.minimize(r_log_likelihood, pars_j);
    pars_j = opt.par();

    //
    beta_j = pars_j.head(M);
    theta_j = pars_j[M];

    //
    loglike_j = loglikelihood(X, y, beta_j, theta_j, p_0, i, t, N);
}


//// R interface
// [[Rcpp::export]]
Rcpp::List mle_itnb_cpp(
        const arma::mat & X, const arma::vec & y,
        const arma::vec & beta_0, const double & theta_0, const double & p_0,
        const int & i, const int & t,
        const double & tolerance, const arma::vec & lambda,
        const int & steps, const bool & exact,
        const int & trace
) {
    //
    const int & N = y.size();
    const int & M = X.n_cols;

    //
    int j = 0;
    bool not_converged;
    std::string convergence_flag;

    //
    arma::vec beta_j = beta_0;
    double theta_j = theta_0;

    //
    double loglike_j = 0.0;

    //
    optimise_tnb(
        beta_j, theta_j, p_0, loglike_j,
        j, not_converged, convergence_flag,
        X, y, i, t, N, M,
        tolerance, lambda,
        steps, exact, trace
    );

    //
    return Rcpp::List::create(
        Rcpp::Named("formula") = 0,
        Rcpp::Named("data") = 0,
        Rcpp::Named("i") = i,
        Rcpp::Named("t") = t,
        Rcpp::Named("loglikelihood") = loglike_j,
        Rcpp::Named("beta") = beta_j,
        Rcpp::Named("theta") = theta_j,
        Rcpp::Named("p") = p_0,
        Rcpp::Named("trace") = 0,
        Rcpp::Named("converged") = !not_converged,
        Rcpp::Named("flag") = convergence_flag
    );
}
