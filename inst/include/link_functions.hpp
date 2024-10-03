#ifndef link_functions
#define link_functions

#include <RcppArmadillo.h>

//
class Link {
public:
    // Constructor
    Link();

    // Functions
    double link(const double & mu);
    double link_inv(const double & eta);
    double link_inv_dev(const double & mu, const arma::vec & x, const int & m);
};

class link_identity : public Link {
    double link(const double & mu);
    double link_inv(const double & eta);
    double link_inv_dev(const double & mu, const arma::vec & x, const int & m);
};

class link_log : public Link {
    double link(const double & mu);
    double link_inv(const double & eta);
    double link_inv_dev(const double & mu, const arma::vec & x, const int & m);
};

class link_sqrt : public Link {
    double link(const double & mu);
    double link_inv(const double & eta);
    double link_inv_dev(const double & mu, const arma::vec & x, const int & m);
};

#endif //link_functions
