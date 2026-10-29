#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <roptim.h>
// [[Rcpp::depends(roptim)]]

#include "rdp.hpp"
#include "links.hpp"
#include "aux_functions.hpp"

using namespace roptim;

////
//
class MLE_ITPOIS : public Functor {
public:
    // Constructor
    MLE_ITPOIS(const arma::mat & X_, const arma::vec & y_, const int & i_, const int & t_, const int & steps_, const bool & exact_, Link & LO_) :
        X(X_), y(y_), N(y.size()), M(X.n_cols), i(i_), t(t_), steps(steps_), exact(exact_), LO(LO_) { }

    // Functor / objective function
    double operator()(const arma::vec & par) override {
        const arma::vec & beta = par.head(M);
        const double & p = 1.0 / (1.0 + std::exp(-par[M]));

        //
        double d = -complete_loglikelihood_itpois(X, y, beta, p, i, t, LO);
        return d;
    }

    // Gradient function
    void Gradient(const arma::vec & par, arma::vec & gr) override {
        //
        gr = arma::zeros(M + 1);
        if (exact) {
            const arma::vec & beta = par.head(M);
            const double & p = 1.0 / (1.0 + std::exp(-par[M]));

            for (int n = 0; n < N; n++) {
                //
                const arma::rowvec & x_n = X.row(n);
                const double eta_n = arma::dot(x_n, beta);
                const double mu_n = LO.link_inv(eta_n);

                //
                const double y_n = y[n];
                const double f_n = ditpois_cpp(y_n, mu_n, 0.0, i, t);

                //
                double gr_beta = -(y_n / mu_n - 1.0);
                if (t > -1) {
                    //
                    const double & b = R::ppois(t, mu_n, 0, 0);
                    const double & b_mu_n = R::dpois(t, mu_n, 0);
                    const double & I_mu = b_mu_n / b;

                    //
                    gr_beta += (I_mu);
                }

                //
                double gr_p;
                if (y_n == i) {
                    const double numerator = std::exp(f_n);
                    const double denominator = p + (1.0 - p) * std::exp(f_n);

                    gr_beta *= ((1.0 - p) * numerator / denominator);
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
                gr[M] += gr_p * p * (1.0 - p);
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

//
void optimise_itpois(
        arma::vec & beta_j, double & p_j, double & loglike_j,
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
    arma::vec pars_j = arma::vec(M + 1);
    pars_j.head(M) = beta_j;
    pars_j[M] = std::log(p_j / (1 - p_j));

    //
    Roptim<MLE_ITPOIS> opt("L-BFGS-B");
    MLE_ITPOIS log_likelihood(X, y, i, t, steps, exact, LO);

    opt.control.trace = trace;

    //
    opt.minimize(log_likelihood, pars_j);
    pars_j = opt.par();

    //
    beta_j = pars_j.head(M);
    p_j = 1.0 / (1.0 + std::exp(-pars_j[M]));

    //
    loglike_j = opt.value();
    log_likelihood.ApproximateHessian(opt.par(), approx_hessian);
}


//// R interface
// [[Rcpp::export]]
Rcpp::List mle_itpois_cpp(
        const arma::mat & X, const arma::vec & y,
        const arma::vec & beta_0, const double & p_0,
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
    double p_j = p_0;

    //
    double loglike_j = 0.0;

    //
    arma::mat approx_hessian = arma::mat(beta_j.size() + 1, beta_j.size() + 1);

    //
    optimise_itpois(
        beta_j, p_j, loglike_j,
        j, not_converged, convergence_flag,
        approx_hessian,
        X, y, i, t, link, N, M,
        tolerance,
        steps, exact, trace
    );

    // HOT FIX
    approx_hessian += 1e-8 * arma::eye(beta_j.size() + 1, beta_j.size() + 1);

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
        Rcpp::Named("logtheta") = 0.0,
        Rcpp::Named("trace") = 0,
        Rcpp::Named("converged") = !not_converged,
        Rcpp::Named("iterations") = 0,
        Rcpp::Named("flag") = convergence_flag
    );
}
