// source: https://stackoverflow.com/questions/35923787/fast-large-matrix-multiplication-in-r

// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>


// [[Rcpp::export]]
SEXP armaMatMult(arma::mat A, arma::mat B){
    arma::mat C = A * B;

    return Rcpp::wrap(C);
}

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
SEXP transClose_W(Rcpp::NumericMatrix x){
    int nr = INTEGER(Rf_getAttrib(x, R_DimSymbol))[0];
    int nc = INTEGER(Rf_getAttrib(x, R_DimSymbol))[1];
    double *px = REAL(x);
    int i, j, k;
    for (k = 0; k < nr; k++) {
        for (i = 0; i < nc; i++) {
            for (j = 0; j < nr; j++) {
                px[i*nc + j] = (round(px[i*nc + j]) ||
                              (round(px[i*nc + k]) && round(px[k*nc + j])));
            }
        }
    }
    return Rcpp::wrap(x);
}

// [[Rcpp::export]]
SEXP transClose_Del(Rcpp::NumericMatrix x, Rcpp::IntegerVector u, Rcpp::IntegerVector v){
    int nr = INTEGER(Rf_getAttrib(x, R_DimSymbol))[0];
    int nc = INTEGER(Rf_getAttrib(x, R_DimSymbol))[1];
    int i = v[0]-1;
    int j = u[0]-1;
    double *px = REAL(x);
    int k;
    for (k = 0; k < nr; k++) {
        px[i*nc + j] = (round(px[i*nc + j]) ||
			(round(px[i*nc + k]) && round(px[k*nc + j])));
	if (px[i*nc + j]==1) {
	    break;
	}
    }
    return Rcpp::wrap(x);
}

// [[Rcpp::export]]
SEXP transClose_Ins(Rcpp::NumericMatrix x, Rcpp::IntegerVector u, Rcpp::IntegerVector v){
    int nr = INTEGER(Rf_getAttrib(x, R_DimSymbol))[0];
    int nc = INTEGER(Rf_getAttrib(x, R_DimSymbol))[1];
    int k = u[0]-1;
    int l = v[0]-1;
    double *px = REAL(x);
    int i, j;
    for (i = 0; i < nc; i++) {
	for (j = 0; j < nr; j++) {
	    px[j*nc + i] = (round(px[j*nc + i]) ||
			    (round(px[j*nc + l]) && round(px[k*nc + i])));
	}
    }
    return Rcpp::wrap(x);
}
// [[Rcpp::export]]
SEXP maxCol_row(Rcpp::NumericMatrix x){
    int nr = INTEGER(Rf_getAttrib(x, R_DimSymbol))[0];
    int nc = INTEGER(Rf_getAttrib(x, R_DimSymbol))[1];
    double *px = REAL(x), *buf = (double *) R_alloc(nr, sizeof(double));
    for(int i = 0; i < nr; i++) buf[i] = R_NegInf;
    SEXP ans = PROTECT(Rf_allocVector(INTSXP, nr));
    int *pans = INTEGER(ans);
    for(int i = 0; i < nr; i++) {
        for(int j = 0; j < nc; j++) {
            if(px[i + j*nr] > buf[i]) {
                buf[i] = px[i + j*nr];
                pans[i] = j + 1;
            }
        }
    }
    UNPROTECT(1);
    return Rcpp::wrap(ans);
}
