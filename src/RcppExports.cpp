// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// em_itnb_cpp
Rcpp::List em_itnb_cpp(const arma::mat& X, const arma::vec& y, const arma::vec& yi, const arma::vec& beta_0, const double& theta_0, const double& p_0, const int& i, const int& t, const std::string& link, const int& iteration_min, const int& iteration_max, const double& tolerance, const arma::vec& lambda, const int& steps, const bool& exact, const int& trace, const bool& save_trace);
RcppExport SEXP _itnb_em_itnb_cpp(SEXP XSEXP, SEXP ySEXP, SEXP yiSEXP, SEXP beta_0SEXP, SEXP theta_0SEXP, SEXP p_0SEXP, SEXP iSEXP, SEXP tSEXP, SEXP linkSEXP, SEXP iteration_minSEXP, SEXP iteration_maxSEXP, SEXP toleranceSEXP, SEXP lambdaSEXP, SEXP stepsSEXP, SEXP exactSEXP, SEXP traceSEXP, SEXP save_traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type yi(yiSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta_0(beta_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type theta_0(theta_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type p_0(p_0SEXP);
    Rcpp::traits::input_parameter< const int& >::type i(iSEXP);
    Rcpp::traits::input_parameter< const int& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type link(linkSEXP);
    Rcpp::traits::input_parameter< const int& >::type iteration_min(iteration_minSEXP);
    Rcpp::traits::input_parameter< const int& >::type iteration_max(iteration_maxSEXP);
    Rcpp::traits::input_parameter< const double& >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const int& >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type exact(exactSEXP);
    Rcpp::traits::input_parameter< const int& >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< const bool& >::type save_trace(save_traceSEXP);
    rcpp_result_gen = Rcpp::wrap(em_itnb_cpp(X, y, yi, beta_0, theta_0, p_0, i, t, link, iteration_min, iteration_max, tolerance, lambda, steps, exact, trace, save_trace));
    return rcpp_result_gen;
END_RCPP
}
// ritnb_cpp
arma::vec ritnb_cpp(const int& n, const arma::vec& mu, const arma::vec& theta, const arma::vec& p, const arma::vec& i, const arma::vec& t);
RcppExport SEXP _itnb_ritnb_cpp(SEXP nSEXP, SEXP muSEXP, SEXP thetaSEXP, SEXP pSEXP, SEXP iSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type i(iSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(ritnb_cpp(n, mu, theta, p, i, t));
    return rcpp_result_gen;
END_RCPP
}
// ditnb_cpp
arma::vec ditnb_cpp(const arma::vec& x, const arma::vec& mu, const arma::vec& theta, const arma::vec& p, const arma::vec& i, const arma::vec& t);
RcppExport SEXP _itnb_ditnb_cpp(SEXP xSEXP, SEXP muSEXP, SEXP thetaSEXP, SEXP pSEXP, SEXP iSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type i(iSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(ditnb_cpp(x, mu, theta, p, i, t));
    return rcpp_result_gen;
END_RCPP
}
// pitnb_cpp
arma::vec pitnb_cpp(const arma::vec& x, const arma::vec& mu, const arma::vec& theta, const arma::vec& p, const arma::vec& i, const arma::vec& t);
RcppExport SEXP _itnb_pitnb_cpp(SEXP xSEXP, SEXP muSEXP, SEXP thetaSEXP, SEXP pSEXP, SEXP iSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type i(iSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(pitnb_cpp(x, mu, theta, p, i, t));
    return rcpp_result_gen;
END_RCPP
}
// mle_itnb_cpp
Rcpp::List mle_itnb_cpp(const arma::mat& X, const arma::vec& y, const arma::vec& beta_0, const double& theta_0, const double& p_0, const int& i, const int& t, const std::string& link, const double& tolerance, const arma::vec& lambda, const int& steps, const bool& exact, const int& trace);
RcppExport SEXP _itnb_mle_itnb_cpp(SEXP XSEXP, SEXP ySEXP, SEXP beta_0SEXP, SEXP theta_0SEXP, SEXP p_0SEXP, SEXP iSEXP, SEXP tSEXP, SEXP linkSEXP, SEXP toleranceSEXP, SEXP lambdaSEXP, SEXP stepsSEXP, SEXP exactSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta_0(beta_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type theta_0(theta_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type p_0(p_0SEXP);
    Rcpp::traits::input_parameter< const int& >::type i(iSEXP);
    Rcpp::traits::input_parameter< const int& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type link(linkSEXP);
    Rcpp::traits::input_parameter< const double& >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const int& >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type exact(exactSEXP);
    Rcpp::traits::input_parameter< const int& >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(mle_itnb_cpp(X, y, beta_0, theta_0, p_0, i, t, link, tolerance, lambda, steps, exact, trace));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_itnb_em_itnb_cpp", (DL_FUNC) &_itnb_em_itnb_cpp, 17},
    {"_itnb_ritnb_cpp", (DL_FUNC) &_itnb_ritnb_cpp, 6},
    {"_itnb_ditnb_cpp", (DL_FUNC) &_itnb_ditnb_cpp, 6},
    {"_itnb_pitnb_cpp", (DL_FUNC) &_itnb_pitnb_cpp, 6},
    {"_itnb_mle_itnb_cpp", (DL_FUNC) &_itnb_mle_itnb_cpp, 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_itnb(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
