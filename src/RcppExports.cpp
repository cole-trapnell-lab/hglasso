// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// timesTwo
NumericVector timesTwo(NumericVector x);
RcppExport SEXP hglasso_timesTwo(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    __result = Rcpp::wrap(timesTwo(x));
    return __result;
END_RCPP
}
// BB_logistic
SEXP BB_logistic(SEXP R_X, SEXP R_A, SEXP R_rho);
RcppExport SEXP hglasso_BB_logistic(SEXP R_XSEXP, SEXP R_ASEXP, SEXP R_rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type R_X(R_XSEXP);
    Rcpp::traits::input_parameter< SEXP >::type R_A(R_ASEXP);
    Rcpp::traits::input_parameter< SEXP >::type R_rho(R_rhoSEXP);
    __result = Rcpp::wrap(BB_logistic(R_X, R_A, R_rho));
    return __result;
END_RCPP
}