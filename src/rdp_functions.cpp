#include "Rcpp.h"

#include "boost/random/mersenne_twister.hpp"
#include "boost/random/gamma_distribution.hpp"
#include "boost/random/poisson_distribution.hpp"
#include "boost/random/uniform_01.hpp"

//// ritnb function
// Generate random variates from a negative binomial distribution
int rnbinom_mu(const double & mu, const double & theta, const int & seed) {
    //
    boost::random::mt19937 rng(seed);

    //
    boost::random::gamma_distribution<> gamma(theta, mu / theta);

    //
    double g = gamma(rng);
    boost::random::poisson_distribution poisson(g);

    //
    int x = poisson(rng);
    return x;
}

// Generate random variates from a truncated negative binomial distribution
int rtnbinom_mu(const double & mu, const double & theta, const int & t, const int & seed) {
    //
    boost::random::mt19937 rng(seed);

    //
    boost::random::gamma_distribution<> gamma(theta, mu / theta);

    //
    int x = t;

    while (x <= t){
        double g = gamma(rng);
        boost::random::poisson_distribution poisson(g);

        x = poisson(rng);
    }

    return x;
}

//[[Rcpp::export]]
std::vector<int> ritnb_mu(const int & n, const double & mu, const double & theta, const double & p, const int & i, const int & t, const int & seed) {
    //
    boost::random::mt19937 rng(seed);

    //
    boost::random::uniform_01<> uniform_01;
    boost::random::uniform_int_distribution uniform(1, 100000);

    //
    std::vector<int> x(n);
    for (int j = 0; j < n; j++) {
        int val;
        const double u = uniform_01(rng);
        if (u < p) {
            val = i;
        }
        else {
            if (t < 0) {
                val = rnbinom_mu(mu, theta, seed + uniform(rng));
            }
            else {
                val = rtnbinom_mu(mu, theta, t, seed + uniform(rng));
            }
        }

        x[j] = val;
    }

    return x;
}

//// ditnb function


//// pitnb function


