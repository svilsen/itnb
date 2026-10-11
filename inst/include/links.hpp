#ifndef links
#define links

#include <RcppArmadillo.h>

//
class Link {
public:
    // Constructor
    Link(const std::string & s);

    // Functions
    double link(const double & mu);
    double link_inv(const double & eta);
    double link_inv_dev(const double & mu, const arma::rowvec & x, const int & m);

private:
    int type;
};

#endif //links
