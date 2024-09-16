#include <RcppArmadillo.h>

#include "rdp_functions.hpp"
#include "em_functions.hpp"

Rcpp::List bs_itnb_cpp() {
    return Rcpp::List::create(
        Rcpp::Named("TEST") = 0.0
    );
}
