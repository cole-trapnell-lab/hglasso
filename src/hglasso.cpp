#include "hglasso.h"

using namespace Rcpp;
using namespace Eigen;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

void gradient_eval_symmetric(const MatrixXd& Theta, 
                             const MatrixXd& Xhat, 
                             const MatrixXd& A, 
                             const MatrixXd& X,
                             const double& rho,
                             MatrixXd& grad)
{
  int p = Theta.rows();
  int n = X.rows();
  
  grad = MatrixXd(p, p);
  
  MatrixXd mat_temp = ArrayXXd(p, n);

  MatrixXd theta_diag = MatrixXd::Zero(p, p);
  theta_diag.diagonal() = Theta.diagonal();
  
  MatrixXd theta_X_t = Theta * X.transpose();
  
  mat_temp =  (Theta.diagonal() * MatrixXd::Constant(1, n, 1.0) + 
          theta_X_t - 
          theta_diag * X.transpose()).array().exp();
  ArrayXXd temp = mat_temp;

  mat_temp = ((Theta.diagonal() * MatrixXd::Constant(1, n, 1.0) + 
          theta_X_t -
          theta_diag * X.transpose()).array().exp() + 1);
  temp /= mat_temp.array();

  for (size_t k = 0; k < p; ++k)
  {
    VectorXd rep_elm = ( (  Theta.row(k) * X.transpose() - X.col(k).transpose() * Theta(k,k)).array() + Theta(k,k)).exp() /
                       ( ( (Theta.row(k) * X.transpose() - X.col(k).transpose() * Theta(k,k)).array() + Theta(k,k)).exp() + 1);

    ArrayXXd temp_mat(p, n);
    for (size_t j = 0; j < p; j++){
      temp_mat.row(j) = rep_elm;
    }
    temp_mat *= X.transpose().array();

    ArrayXXd temp_mat2(p, n);
    for (size_t j = 0; j < p; j++){
      temp_mat2.row(j) = X.col(k).transpose();
    }
    
    VectorXd grad_row_tmp = temp_mat.transpose().colwise().sum() + 
                            (temp*temp_mat2).transpose().colwise().sum();

    grad_row_tmp +=  -Xhat.row(k) - Xhat.col(k).transpose(); 

    grad_row_tmp += (Theta.row(k) - A.row(k))*rho + (Theta.col(k).transpose() - A.col(k).transpose())*rho;

    grad.row(k) = grad_row_tmp;

    ArrayXd vec_row_tmp = (  (Theta.row(k) * X.transpose() - X.col(k).transpose() * Theta(k,k)).array() + Theta(k,k)).array().exp();
    vec_row_tmp /= (  (  (Theta.row(k) * X.transpose() - X.col(k).transpose() * Theta(k,k)).array() + Theta(k,k)).array().exp() + 1);

    double update_k = vec_row_tmp.sum() + 
                      rho*(Theta(k,k) - A(k,k)) - Xhat(k,k);

    grad(k,k) = update_k;

  }
  
  
  
}

void BB_logistic_cpp(const MatrixXd& X, const MatrixXd& A, const double& rho,
                     MatrixXd& Theta)
{
  MatrixXd Xhat = X.transpose() * X;
  //Rcpp::Rcout << Xhat << std::endl;
  int n = X.rows();
  int p = X.cols();
  Theta = MatrixXd::Identity(p, p);
  
  MatrixXd Theta_old = 2 * Theta;
  
  MatrixXd D = MatrixXd(p, p);
  MatrixXd D_old = MatrixXd(p, p);
  
  gradient_eval_symmetric(Theta_old, Xhat, A, X, rho, D_old);
  
  MatrixXd S = MatrixXd(p, p);
  MatrixXd Y = MatrixXd(p, p);
  double alpha_bb = 0.0;
  
  // Algorithm parameters
  size_t max_iter = 400;
  double tol = 1e-5;
  
  size_t iter = 0;
  // Main algorithm
  for (size_t ind = 0; ind < max_iter; ++ind)
  {
    S = Theta - Theta_old;
    Theta_old = Theta;
    
    gradient_eval_symmetric(Theta_old, Xhat, A, X, rho, D);
    
    Y = D - D_old;
    D_old = D;
    
    
    alpha_bb = std::min( (S.array()*S.array()).sum()/(S.array()*Y.array()).sum(), (S.array()*Y.array()).sum()/(Y.array()*Y.array()).sum() ); // BB step size
    
    if (isnan(alpha_bb)){
      break;
    }
    
    //Rcpp::Rcout << alpha_bb << std::endl;
    Theta = Theta - alpha_bb*D; // Gradient descent step
      
    double rel_err = (Theta - Theta_old).norm()/Theta_old.norm();
    if (rel_err < tol)
    {
      iter = ind;
      break;		
    }	
  }
}

// [[Rcpp::export]]
SEXP BB_logistic(SEXP R_X, SEXP R_A, SEXP R_rho) {
  NumericMatrix Rcpp_X(R_X);   
  

  const int n = Rcpp_X.nrow(), p = Rcpp_X.ncol();
  Map<MatrixXd> X(Rcpp_X.begin(), n, p); 
  
  //Rcpp::Rcout << X << std::endl;
  
  NumericMatrix Rcpp_A(R_A);              
  Map<MatrixXd> A(Rcpp_A.begin(), p, p);	 
  
  //Rcpp::Rcout << A << std::endl;
  
  double rho = as<double>(R_rho);
  
  MatrixXd Theta;
  
  BB_logistic_cpp(X, A, rho, Theta);
  //return wrap(1);
  return wrap(Theta);
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
