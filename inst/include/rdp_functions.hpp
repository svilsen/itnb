#ifndef rdp_functions
#define rdp_functions

#include <Rcpp.h>

//
std::vector<int> ritnb_cpp(const int & n, const std::vector<double> & mu, const std::vector<double> & theta, const std::vector<double> & p, const std::vector<int> & i, const std::vector<int> & t, const int & seed);

//
double ditnb_cpp(const int & x, const double & mu, const double & theta, const double & p, const int & i, const int & t);
std::vector<double> ditnb_cpp(const std::vector<int> & x, const std::vector<double> & mu, const std::vector<double> & theta, const std::vector<double> & p, const std::vector<int> & i, const std::vector<int> & t);

//
std::vector<double> pitnb_cpp(const std::vector<int> & x, const std::vector<double> & mu, const std::vector<double> & theta, const std::vector<double> & p, const std::vector<int> & i, const std::vector<int> & t);

#endif //rdp_functions
