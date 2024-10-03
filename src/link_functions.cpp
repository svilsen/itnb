#include <RcppArmadillo.h>
#include "link_functions.hpp"

//
Link::Link() {

}

double Link::link(const double & mu) {
    return 0.0;
}

double Link::link_inv(const double & eta) {
    return 0.0;
}

double Link::link_inv_dev(const double & mu, const arma::vec & x, const int & m) {
    return 0.0;
}


// Identity
double link_identity::link(const double & mu) {
    return mu;
}

double link_identity::link_inv(const double & eta) {
    double mu = eta;
    if (mu < 1e-6) {
        mu = 1e-6;
    }

    return mu;
}

double link_identity::link_inv_dev(const double & mu, const arma::vec & x, const int & m) {
    return x[m];
}

// Log.
double link_log::link(const double & mu) {
    return std::log(mu);
}

double link_log::link_inv(const double & eta) {
    return std::exp(eta);
}

double link_log::link_inv_dev(const double & mu, const arma::vec & x, const int & m) {
    return mu * x[m];
}

// Sqrt.
double link_sqrt::link(const double & mu) {
    return std::sqrt(mu);
}

double link_sqrt::link_inv(const double & eta) {
    return eta * eta;
}

double link_sqrt::link_inv_dev(const double & mu, const arma::vec & x, const int & m) {
    double eta = link_inv(mu);
    return 2.0 * eta * x[m];
}
