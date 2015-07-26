
#ifndef _hglasso_HGLASSO_H
#define _hglasso_HGLASSO_H

#include <Rcpp.h>
#include <RcppEigen.h>

//#define EIGEN_NO_DEBUG

using namespace Eigen;

//RcppExport SEXP BB_logistic(SEXP R_X, SEXP R_A, SEXP R_rho);

void gradient_eval_symmetric(const MatrixXd& Theta, 
                             const MatrixXd& Xhat, 
                             const MatrixXd& X, 
                             const MatrixXd& A,
                             const double& rho);
// functions to be used within C++
void BB_logistic_cpp(const MatrixXd& X, const MatrixXd& A, const double& rho,
                     MatrixXd& Theta);

#endif
