#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// // [[Rcpp::export]]
// NumericVector timesTwo(NumericVector x) {
//   return x * 2;
// }


// [[Rcpp::export]]
double cppCondNumber(Eigen::MatrixXd X) {
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(X);
    const Eigen::VectorXd D = svd.singularValues();
    return D[0] / D[D.size() - 1];
}

// [[Rcpp::export]]
double cppCdev(const Eigen::Map<Eigen::MatrixXd>& X, const Eigen::Map<Eigen::MatrixXd>& Y){
    // fast.condNumber((X %*% MASS::ginv(Y))))
    Eigen::MatrixXd Xinv;
    Eigen::MatrixXd Z;
    if (X.cols() > X.rows()) {
        Xinv = X.transpose().completeOrthogonalDecomposition().pseudoInverse();
        Z = Xinv * Y.transpose();
    } else {
        Xinv = X.completeOrthogonalDecomposition().pseudoInverse();
        Z = Xinv * Y;
    }
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(Z);
    const auto D = svd.singularValues();
    return D[0] / D[D.size() - 1];
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

// /*** R
// timesTwo(42)
// */
