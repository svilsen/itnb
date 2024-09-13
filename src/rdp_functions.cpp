#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//// ritnb function
// Generate random variates from a truncated negative binomial distribution
int rtnbinom_cpp(const double & mu, const double & theta, const int & t) {
    //
    int x = t;
    while (x <= t){
        x = R::rnbinom(theta, theta / (mu + theta));
    }

    //
    return x;
}

// [[Rcpp::export]]
arma::vec ritnb_cpp(const int & n, const arma::vec & mu, const arma::vec & theta, const arma::vec & p, const arma::vec & i, const arma::vec & t) {
    //
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
    arma::vec x(n);
    for (int j = 0; j < n; j++) {
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

        int val;
        const double u = R::runif(0, 1);
        if (u < p_n) {
            val = i_n;
        }
        else {
            if (t_n < 0) {
                val = R::rnbinom(theta_n, theta_n / (mu_n + theta_n));
            }
            else {
                val = rtnbinom_cpp(mu_n, theta_n, t_n);
            }
        }

        x[j] = val;
    }

    //
    return x;
}

//// ditnb function
double ditnb_cpp(const int & x, const double & mu, const double & theta, const double & p, const int & i, const int & t) {
    //
    const double tm = theta + mu;
    double log_d = theta * std::log(theta) - theta * std::log(theta + mu) + std::lgamma(theta + x) - std::lgamma(theta) - std::lgamma(x + 1) + x * std::log(mu) - x * std::log(theta + mu);

    //
    if (t > -1) {
        if (x <= t) {
            log_d = -HUGE_VAL;
        }
        else {
            //
            const double pb = R::pbeta(mu / tm, t + 1, theta, true, true);
            log_d = log_d - pb;

            //
            if (p > 1e-16) {
                log_d = std::log(1.0 - p) + log_d;
                if (x == i) {
                    log_d = std::log1p(p + std::expm1(log_d));
                }
            }
        }
    }

    //
    return log_d;
}

// [[Rcpp::export]]
arma::vec ditnb_cpp(const arma::vec & x, const arma::vec & mu, const arma::vec & theta, const arma::vec & p, const arma::vec & i, const arma::vec & t) {
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
    arma::vec d(N_x);
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

        d[n] = ditnb_cpp(x[n], mu_n, theta_n, p_n, i_n, t_n);
    }

    return d;
}


//// pitnb function
// [[Rcpp::export]]
arma::vec pitnb_cpp(const arma::vec & x, const arma::vec & mu, const arma::vec & theta, const arma::vec & p, const arma::vec & i, const arma::vec & t) {
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
    int m_n = std::max(t_n, 0);

    //
    arma::vec d(N_x);
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
            m_n = std::max(t_n, 0);
        }

        d[n] = 0.0;
        for (int m = m_n; m < (x[n] + 1); m++) {
            const double log_d_nm = ditnb_cpp(m, mu_n, theta_n, p_n, i_n, t_n);
            d[n] += std::exp(log_d_nm);
        }

        d[n] = std::log(d[n]);
    }

    return d;
}

