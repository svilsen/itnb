#include <RcppArmadillo.h>
#include "links.hpp"

//
Link::Link(const std::string & s) {
    if (s == "identity") {
        type = 0;
    }
    else if (s == "sqrt") {
        type = 1;
    }
    else if (s == "log") {
        type = 2;
    }
    else { // reciprocal
        type = 3;
    }
}

double Link::link(const double & mu) {
    switch (type) {
        case 0 : // identity
            return mu;
        case 1 : // sqrt
            return std::sqrt(mu);
        case 2 : // log
            return std::log(mu);
        case 3 : // reciprocal
            return 1.0 / mu;
    }

    return mu;
}

double Link::link_inv(const double & eta) {
    switch (type) {
        case 0 : // identity
            return std::max(eta, 1e-6);
        case 1 : // sqrt
            return eta * eta;
        case 2 : // log
            return std::exp(eta);
        case 3 : // reciprocal
            return 1.0 / eta;
    }

    return eta;
}

double Link::link_inv_dev(const double & mu, const arma::rowvec & x, const int & m) {
    switch (type) {
        case 0 : // identity
            return x[m];
        case 1 : // sqrt
            return 2.0 * std::sqrt(mu) * x[m];
        case 2 : // log
            return mu * x[m];
        case 3 : // reciprocal
            return -mu * mu;
    }

    return x[m];
}

