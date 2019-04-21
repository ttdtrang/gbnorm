// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
Rcpp::List eigenJacobiSVD(const Eigen::Map<Eigen::MatrixXd> X){
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeThinU|Eigen::ComputeThinV);
    const Eigen::MatrixXd U = svd.matrixU();
    const Eigen::MatrixXd V = svd.matrixV();
    const Eigen::VectorXd D = svd.singularValues();
    return List::create(
            _["u"] = U,
            _["v"] = V,
            _["d"] = D);
}

// [[Rcpp::export]]
Rcpp::List eigenBDCSVD(const Eigen::Map<Eigen::MatrixXd> X){
    const Eigen::BDCSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeThinU|Eigen::ComputeThinV);
    const Eigen::MatrixXd U = svd.matrixU();
    const Eigen::MatrixXd V = svd.matrixV();
    const Eigen::VectorXd D = svd.singularValues();
    return List::create(
            _["u"] = U,
            _["v"] = V,
            _["d"] = D);
}

// [[Rcpp::export]]
SEXP eigenJacobiSV(const Eigen::Map<Eigen::MatrixXd> X){
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(X);
    const Eigen::VectorXd D = svd.singularValues();
    return Rcpp::wrap(D);
}

// [[Rcpp::export]]
SEXP eigenBDCSV(const Eigen::Map<Eigen::MatrixXd> X){
    const Eigen::BDCSVD<Eigen::MatrixXd> svd(X);
    const Eigen::VectorXd D = svd.singularValues();
    return Rcpp::wrap(D);
}

// [[Rcpp::export]]
SEXP eigenPseudoInverse(const Eigen::Map<Eigen::MatrixXd> X){
    const Eigen::MatrixXd Xinv = X.completeOrthogonalDecomposition().pseudoInverse();
    return Rcpp::wrap(Xinv); 
}
