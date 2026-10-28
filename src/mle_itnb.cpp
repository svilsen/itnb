#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <roptim.h>
// [[Rcpp::depends(roptim)]]

#include "rdp.hpp"
#include "links.hpp"
#include "aux_functions.hpp"

using namespace roptim;

//// Truncated negative binomial parameters
class MLE_ITNB : public Functor {
public:
    // Constructor
    MLE_ITNB(const arma::mat & X_, const arma::vec & y_, const int & i_, const int & t_, const int & steps_, const bool & exact_, Link & LO_) :
        X(X_), y(y_), N(y.size()), M(X.n_cols), i(i_), t(t_), steps(steps_), exact(exact_), LO(LO_) { }

    // Functor / objective function
    double operator()(const arma::vec & par) override {
        const arma::vec & beta = par.head(M);
        const double & theta = std::exp(par[M]);
        const double & p = 1.0 / (1.0 + std::exp(-par[M + 1]));

        //
        double d = 0.0;
        for (int n = 0; n < N; n++) {
            //
            const arma::rowvec & x_n = X.row(n);
            const double & eta_n = arma::dot(x_n, beta);
            const double & mu_n = LO.link_inv(eta_n);

            //
            d -= ditnb_cpp(y[n], mu_n, theta, p, i, t);
        }

        if (std::isinf(d) || std::isnan(d)) {
            d = 1e16;
        }

        return d;
    }

    // Gradient function
    void Gradient(const arma::vec & par, arma::vec & gr) override {
        //
        gr = arma::zeros(M + 2);
        if (exact) {
            const arma::vec & beta = par.head(M);
            const double & theta = std::exp(par[M]);
            const double & p = 1.0 / (1.0 + std::exp(-par[M + 1]));

            for (int n = 0; n < N; n++) {
                //
                const arma::rowvec & x_n = X.row(n);
                const arma::vec eta_n = x_n * beta;
                const double mu_n = LO.link_inv(eta_n[0]);

                //
                const double y_n = y[n];
                const double f_n = ditnb_cpp(y[n], mu_n, theta, 0.0, i, t);

                //
                const double & mu_theta = mu_n + theta;
                const double & mu_theta_i = 1.0 / mu_theta;

                const double & y_theta = y_n + theta;

                //
                double gr_beta = -(y_n / mu_n - y_theta * mu_theta_i);
                double gr_theta = -(1.0 + std::log(theta) - std::log(mu_theta) - y_theta * mu_theta_i + R::digamma(y_theta) - R::digamma(theta));

                if (t > -1) {
                    //
                    const double & b = std::exp(R::pbeta(mu_n * mu_theta_i, t + 1, theta, 1, 1) + R::lbeta(t + 1, theta));

                    // Location
                    const double & b_mu_n = std::pow(mu_n, t) * std::pow(theta, theta) / std::pow(mu_theta, t + 1 + theta);
                    const double & I_mu = b_mu_n / b;

                    gr_beta += (I_mu);

                    // Overdispersion
                    const double & b_t_f = std::pow(mu_n, t + 1.0) * std::pow(theta, theta - 1.0) / std::pow(mu_theta, t + 1 + theta);
                    const double & b_t_i = beta_derivative_boole(0.0, mu_n * mu_theta_i, theta, t, steps);
                    const double & b_t_d = R::digamma(theta) - R::digamma(t + 1.0 + theta);

                    const double & I_theta = (b_t_i - b_t_f) / b - b_t_d;
                    gr_theta += (I_theta);
                }


                double gr_p;
                if (y_n == i) {
                    const double numerator = std::exp(f_n);
                    const double denominator = p + (1.0 - p) * std::exp(f_n);

                    gr_beta *= ((1.0 - p) * numerator / denominator);
                    gr_theta *= ((1.0 - p) * numerator / denominator);

                    gr_p = (numerator - 1.0) / denominator;
                }
                else {
                    gr_p = 1.0 / (1.0 - p);
                }

                //
                for (int m = 0; m < M; m++) {
                    const double & gr_mu = LO.link_inv_dev(mu_n, x_n, m);
                    gr[m] += gr_beta * gr_mu;
                }

                //
                gr[M] += gr_theta * theta;

                //
                gr[M + 1] += gr_p * p * (1.0 - p);
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

    // Link
    Link LO;
};

//// Optimiser
void optimise_itnb(
        arma::vec & beta_j, double & theta_j, double & p_j, double & loglike_j,
        int & j, bool & not_converged, std::string & convergence_flag, arma::mat & approx_hessian,
        const arma::mat & X, const arma::vec & y,
        const int & i, const int & t, const std::string & link,
        const int & N, const int & M,
        const double & tolerance,
        const int & steps, const bool & exact,
        const int & trace
) {
    //
    not_converged = true;
    convergence_flag = "";

    //
    Link LO(link);

    //
    arma::vec pars_j = arma::vec(M + 2);
    pars_j.head(M) = beta_j;
    pars_j[M] = std::log(theta_j);
    pars_j[M + 1] = std::log(p_j / (1 - p_j));

    //
    Roptim<MLE_ITNB> opt("L-BFGS-B");
    MLE_ITNB log_likelihood(X, y, i, t, steps, exact, LO);

    opt.control.trace = trace;

    //
    opt.minimize(log_likelihood, pars_j);
    pars_j = opt.par();

    //
    beta_j = pars_j.head(M);
    theta_j = std::exp(pars_j[M]);
    p_j = 1.0 / (1.0 + std::exp(-pars_j[M + 1]));

    //
    loglike_j = opt.value();
    log_likelihood.ApproximateHessian(opt.par(), approx_hessian);
}


//// R interface
// [[Rcpp::export]]
Rcpp::List mle_itnb_cpp(
        const arma::mat & X, const arma::vec & y,
        const arma::vec & beta_0, const double & theta_0, const double & p_0,
        const int & i, const int & t, const std::string & link,
        const double & tolerance,
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
    double p_j = p_0;

    //
    double loglike_j = 0.0;

    //
    arma::mat approx_hessian = arma::mat(beta_j.size() + 2, beta_j.size() + 2);

    //
    optimise_itnb(
        beta_j, theta_j, p_j, loglike_j,
        j, not_converged, convergence_flag,
        approx_hessian,
        X, y, i, t, link, N, M,
        tolerance,
        steps, exact, trace
    );

    // HOT FIX
    approx_hessian += 1e-8 * arma::eye(beta_j.size() + 2, beta_j.size() + 2);

    //
    arma::mat vcov = arma::inv(approx_hessian);
    double se_logtheta = std::sqrt(vcov(beta_j.size(), beta_j.size()));

    vcov.shed_col(beta_j.size());
    vcov.shed_row(beta_j.size());

    //
    return Rcpp::List::create(
        Rcpp::Named("formula") = 0,
        Rcpp::Named("data") = 0,
        Rcpp::Named("i") = i,
        Rcpp::Named("t") = t,
        Rcpp::Named("link") = link,
        Rcpp::Named("loglikelihood") = -loglike_j,
        Rcpp::Named("beta") = beta_j,
        Rcpp::Named("alpha") = 1.0 / theta_j,
        Rcpp::Named("p") = p_j,
        Rcpp::Named("vcov") = vcov,
        Rcpp::Named("logtheta") = se_logtheta,
        Rcpp::Named("trace") = 0,
        Rcpp::Named("converged") = !not_converged,
        Rcpp::Named("iterations") = 0,
        Rcpp::Named("flag") = convergence_flag
    );
}
