#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <roptim.h>
// [[Rcpp::depends(roptim)]]

#include "rdp_functions.hpp"

using namespace roptim;

////
double loglikelihood(const arma::mat & X, const arma::vec & y, const arma::vec & beta, const double & theta, const double & p, const int & i, const int & t, const int & N) {
    //
    double d = 0.0;
    for (int n = 0; n < N; n++) {
        //
        const arma::rowvec & x_n = X.row(n);
        const arma::vec mu_n = x_n * beta;
        const double & mu = mu_n[0];

        //
        d += ditnb_cpp(y[n], mu, theta, p, i, t);
    }

    return d;
}

//// E-step
void update_z(arma::vec & z, const arma::mat & X, const arma::vec & y, const arma::vec & yi, const arma::vec & beta, const double & theta, const double & p, const int & i, const int & t, const int & N) {
    //
    for (int n = 0; n < N; n++) {
        //
        const arma::rowvec & x_n = X.row(n);
        const arma::vec mu_n = x_n * beta;
        const double & mu = mu_n[0];

        //
        double d_n = std::exp(ditnb_cpp(y[n], mu, theta, 0.0, i, t));
        z[n] = (p * yi[n]) / (p * yi[n] + (1.0 - p) * d_n);
    }
}

//// M-step
// Inflation proportion
void update_p(double & p, const arma::vec & yi, const arma::vec & z, const int & N) {
    double num = 0.0;
    double denom = 0.0;
    for (int n = 0; n < N; n++) {
        num += (z[n] * yi[n]);
        denom += (1.0 + z[n] * (yi[n] - 1.0));
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

double beta_derivative_rectangle(const double & a, const double & b, const double & theta, const int & t, const int & N) {
    //
    const double delta_x = (b - a) / N;
    double i = 0.0;
    for (int n = 0; n < N; n++) {
        i += beta_derivative_integrand(a + n * delta_x, t, theta) * delta_x;
    }

    return i;
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

class EM : public Functor {
public:
    // Parameters
    arma::vec z;

    // Constructor
    EM(const arma::mat & X_, const arma::vec & y_, const arma::vec & yi_, const arma::vec & z_, const int & i_, const int & t_, const int & steps_, const bool & exact_) : z(z_), X(X_), y(y_), yi(yi_), N(y.size()), M(X.n_cols), i(i_), t(t_), steps(steps_), exact(exact_) { }

    // Functor / objective function
    double operator()(const arma::vec & par) override {
        const arma::vec & beta = par.head(M);
        const double & theta = par[M];

        //
        double d = 0.0;
        for (int n = 0; n < N; n++) {
            //
            const arma::rowvec & x_n = X.row(n);
            const arma::vec mu_n = x_n * beta;
            double mu = mu_n[0];

            if (mu < 1e-8) {
                mu = 1e-8; // t + 1;
            }

            d += (z[n] - 1.0) * ditnb_cpp(y[n], mu, theta, 0.0, i, t);
        }

        return d;
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

                if (mu < 1e-8) {
                    mu = 1e-8; // t + 1;
                }

                ////
                const double & mu_theta = mu + theta;
                const double & mu_theta_i = 1.0 / mu_theta;

                const double & y_theta = y[n] + theta;

                ////
                double gr_beta = (z[n] - 1.0) * (y[n] / mu - y_theta * mu_theta_i);
                double gr_theta = (z[n] - 1.0) * (1.0 + std::log(theta) - std::log(mu_theta) - y_theta * mu_theta_i + R::digamma(y_theta) - R::digamma(theta));
                if (t > -1) {
                    const double & b = std::exp(R::pbeta(mu * mu_theta_i, t + 1, theta, true, true) + R::lbeta(t + 1, theta));

                    //// Location
                    const double & b_mu_n = std::pow(mu, t) * std::pow(theta, theta) / std::pow(mu_theta, t + 1 + theta);
                    const double & I_mu = b_mu_n / b;

                    gr_beta += (z[n] - 1.0) * (-I_mu);

                    //// Overdispersion
                    const double & b_t_f = std::pow(mu, t + 1.0) * std::pow(theta, theta - 1.0) / std::pow(mu_theta, t + 1 + theta);
                    const double & b_t_i = beta_derivative_boole(0.0, mu * mu_theta_i, theta, t, steps);
                    const double & b_t_d = R::digamma(theta) - R::digamma(t + 1.0 + theta);

                    const double & I_theta = (b_t_i - b_t_f) / b - b_t_d;
                    gr_theta += (z[n] - 1.0) * (-I_theta);
                }

                for (int m = 0; m < M; m++) {
                    gr[m] += gr_beta * x_n[m];
                }

                gr[M] += gr_theta;
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
};

////
void optimise_itnb(
        arma::vec & beta_j, double & theta_j, double & p_j, double & loglike_j,
        int & j, bool & not_converged, std::string & convergence_flag,
        arma::mat & beta_trace, arma::vec & theta_trace, arma::vec & p_trace, arma::vec & loglike_trace,
        const arma::mat & X, const arma::vec & y, const arma::vec yi,
        const int & i, const int & t, const int & N, const int & M,
        const int & iteration_min, const int & iteration_max,
        const double & tolerance,
        const int & steps, const bool & exact,
        const int & trace, const bool & save_trace
) {
    not_converged = true;
    convergence_flag = "";

    loglike_j = loglikelihood(X, y, beta_j, theta_j, p_j, i, t, N);
    double loglike_j_old = HUGE_VAL;
    double delta_loglike_j = loglike_j - loglike_j_old;

    arma::vec z(N);
    update_z(z, X, y, yi, beta_j, theta_j, p_j, i, t, N);

    //
    arma::vec pars_j = arma::vec(M + 1);
    arma::vec pars_j_old = pars_j;

    pars_j.head(M) = beta_j;
    pars_j[M] = theta_j;

    //
    Roptim<EM> opt("L-BFGS-B");
    EM r_log_likelihood(X, y, yi, z, i, t, steps, exact);

    //
    arma::vec lb = (-HUGE_VAL) * arma::ones(M + 1);
    lb[M] = 1e-8;

    opt.set_lower(lb);

    if (trace > 0) {
        Rcpp::Rcout << "Iteration: " << j << "\t Current log-likelihood: " << loglike_j << "\t Change in log-likelihood: " << delta_loglike_j << "\n"
                    << "\t Parameters: " << "\t beta = " << beta_j.t() << "\t theta = " << theta_j << "\t p = " << p_j << "\n";
    }

    //
    j = 1;
    while (not_converged) {
        ////
        //
        const double p_j_old = p_j;

        //
        pars_j_old.head(M) = beta_j;
        pars_j_old[M] = theta_j;

        //// E-step
        update_z(z, X, y, yi, beta_j, theta_j, p_j, i, t, N);

        //// M-step
        // Inflation proportion
        update_p(p_j, yi, z, N);

        // NB parameters
        r_log_likelihood.z = z;

        opt.minimize(r_log_likelihood, pars_j_old);
        arma::vec pars_j = opt.par();

        beta_j = pars_j.head(M);
        theta_j = pars_j[M];

        //// Convergence
        loglike_j_old = loglike_j;
        loglike_j = loglikelihood(X, y, beta_j, theta_j, p_j, i, t, N);
        delta_loglike_j = loglike_j - loglike_j_old;

        if (delta_loglike_j < 0.0) {
            not_converged = false;
            convergence_flag = "DECREASING LIKELIHOOD DETECTED.";

            p_j = p_j_old;
            beta_j = pars_j_old.head(M);
            theta_j = pars_j_old[M];

            loglike_j = loglike_j_old;
            break;
        }

        if (j > iteration_min) {
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
                            << "\t Parameters: " << "\t beta = " << beta_j.t() << "\t theta = " << theta_j << "\t p = " << p_j << "\n";
            }
        }

        if (save_trace) {
            beta_trace.row(j) = beta_j.t();
            theta_trace[j] = theta_j;
            p_trace[j] = p_j;
            loglike_trace[j] = loglike_j;
        }

        ////
        j++;
    }
}


////
// [[Rcpp::export]]
Rcpp::List em_itnb_cpp(
        const arma::mat & X, const arma::vec & y, const arma::vec & yi,
        const arma::vec & beta_0, const double & theta_0, const double & p_0,
        const int & i, const int & t,
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
    std::string convergence_flag;

    //
    double p_j = p_0;
    if (p_j < 1e-16) {
        p_j = 4e-16;
    }

    //
    arma::vec beta_j = beta_0;
    double theta_j = theta_0;

    //
    double loglike_j = 0.0;

    //
    arma::mat beta_trace;
    arma::vec theta_trace, p_trace, loglike_trace;
    if (save_trace) {
        beta_trace.set_size(iteration_max, M);
        beta_trace.row(0) = beta_0.t();

        theta_trace.set_size(iteration_max);
        theta_trace[0] = theta_0;

        p_trace.set_size(iteration_max);
        p_trace[0] = p_0;

        loglike_trace.set_size(iteration_max);
        loglike_trace[0] = HUGE_VAL;
    }

    //
    optimise_itnb(
        beta_j, theta_j, p_j, loglike_j,
        j, not_converged, convergence_flag,
        beta_trace, theta_trace, p_trace, loglike_trace,
        X, y, yi, i, t, N, M,
        iteration_min, iteration_max, tolerance,
        steps, exact,
        trace, save_trace
    );

    //
    Rcpp::List trace_list;
    if (save_trace) {
        trace_list = Rcpp::List::create(
            Rcpp::Named("LogLikelihood") = loglike_trace.head(j),
            Rcpp::Named("beta") = beta_trace.rows(0, j - 1),
            Rcpp::Named("theta") = theta_trace.head(j),
            Rcpp::Named("p") = p_trace.head(j)
        );
    }

    return Rcpp::List::create(
        Rcpp::Named("formula") = 0,
        Rcpp::Named("data") = 0,
        Rcpp::Named("i") = i,
        Rcpp::Named("t") = t,
        Rcpp::Named("loglikelihood") = loglike_j,
        Rcpp::Named("beta") = beta_j,
        Rcpp::Named("theta") = theta_j,
        Rcpp::Named("p") = p_j,
        Rcpp::Named("trace") = trace_list,
        Rcpp::Named("converged") = !not_converged,
        Rcpp::Named("flag") = convergence_flag
    );
}
