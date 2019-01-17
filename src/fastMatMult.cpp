#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;

//fast matrix multiplication

// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;

  return Rcpp::wrap(C);
}

