#include <RcppArmadillo.h>
#include "link_functions.hpp"

//
Link::Link(const std::string & s_) : s(s_) { }

double Link::link(const double & mu) {
    if (s == "identity") {
        return mu;
    }
    else if (s == "sqrt") {
        return std::sqrt(mu);
    }
    else {
        return std::log(mu);
    }
}

double Link::link_inv(const double & eta) {
    if (s == "identity") {
        double mu = eta;
        if (mu < 1e-6) {
            mu = 1e-6;
        }

        return mu;
    }
    else if (s == "sqrt") {
        return eta * eta;
    }
    else {
        return std::exp(eta);
    }
}

double Link::link_inv_dev(const double & mu, const arma::rowvec & x, const int & m) {
    if (s == "identity") {
        return x[m];
    }
    else if (s == "sqrt") {
        double eta = link_inv(mu);
        return 2.0 * eta * x[m];
    }
    else {
        return mu * x[m];
    }
}

// // Identity
// double link_identity::link(const double & mu) {
//     return mu;
// }
//
// double link_identity::link_inv(const double & eta) {
//     double mu = eta;
//     if (mu < 1e-6) {
//         mu = 1e-6;
//     }
//
//     return mu;
// }
//
// double link_identity::link_inv_dev(const double & mu, const arma::vec & x, const int & m) {
//     return x[m];
// }
//
// // Sqrt.
// double link_sqrt::link(const double & mu) {
//     return std::sqrt(mu);
// }
//
// double link_sqrt::link_inv(const double & eta) {
//     return eta * eta;
// }
//
// double link_sqrt::link_inv_dev(const double & mu, const arma::vec & x, const int & m) {
//     double eta = link_inv(mu);
//     return 2.0 * eta * x[m];
// }
//
// // Log.
// double link_log::link(const double & mu) {
//     return std::log(mu);
// }
//
// double link_log::link_inv(const double & eta) {
//     return std::exp(eta);
// }
//
// double link_log::link_inv_dev(const double & mu, const arma::vec & x, const int & m) {
//     return mu * x[m];
// }

