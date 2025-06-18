// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
void generate_snp_matrix(int n_population, int n_loci, Rcpp::NumericVector snp_mafs) {
    arma::umat snp_matrix(n_population, n_loci);
    for (int i = 0; i < n_population; i++) {
        for (int j = 0; j < n_loci; j++) {
            snp_matrix(i, j) = R::rbinom(2, )
        }
    }
}