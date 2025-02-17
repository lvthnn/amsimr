// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat swap_sol(const arma::mat &sol, const arma::uvec &swap,
                   const arma::uvec &cf_idx) {
  arma::mat new_sol = sol;
  new_sol.submat(swap, cf_idx) = sol.submat(arma::reverse(swap), cf_idx);
  return new_sol;
}

// [[Rcpp::export]]
double psi(const arma::mat &sol, const arma::uvec &phi_m,
           const arma::uvec &phi_f, const arma::vec &w) {
  arma::mat delta = sol.cols(phi_m) - sol.cols(phi_f);
  arma::rowvec score = arma::sum(arma::square(delta), 0) % w.t();
  return arma::sum(score);
}

// [[Rcpp::export]]
double dpsi(const arma::mat &sol_prop, const arma::mat &sol,
            const arma::uvec &phi_m, const arma::uvec &phi_f,
            const arma::vec w) {
  return psi(sol_prop, phi_m, phi_f, w) - psi(sol, phi_m, phi_f, w);
}

// [[Rcpp::export]]
Rcpp::List optim_matching(arma::mat &sol, const arma::mat &psi_vec,
                          const arma::uvec &cf_idx, int n_iter = 10000,
                          double alpha = 0.9995, double temp0 = 5,
                          bool eval = false, bool progress = false) {
  // Set up some variables used in the annealing algorithm
  int n_pairs = sol.n_rows;
  arma::uvec phi_m = arma::conv_to<arma::uvec>::from(psi_vec.col(0));
  arma::uvec phi_f = arma::conv_to<arma::uvec>::from(psi_vec.col(1));
  arma::vec w = psi_vec.col(2);
  double temp = temp0;

  // Simulated annealing algorithm
  for (int i = 0; i < n_iter; i++) {
    arma::uvec swap = arma::randperm(n_pairs, 2);
    arma::mat sol_prop = swap_sol(sol, swap, cf_idx);
    double dpsi_prop = dpsi(sol_prop, sol, phi_m, phi_f, w);
    double rho_prop = std::min(1.0, std::exp(-dpsi_prop / temp));
    double u = Rcpp::as<double>(Rcpp::runif(1));

    if (u < rho_prop) {
      sol = sol_prop;
    }

    temp *= alpha;
  }

  return Rcpp::List::create(Rcpp::Named("sol") = sol);
}
