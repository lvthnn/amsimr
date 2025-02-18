// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
double psi(const arma::mat &sol, const arma::uvec &phi_m,
           const arma::uvec &phi_f, const arma::vec &w) {
  arma::mat delta = sol.cols(phi_m) - sol.cols(phi_f);
  arma::rowvec score = arma::sum(arma::square(delta), 0) % w.t();
  return arma::sum(score);
}

// [[Rcpp::export]]
double dpsi(const arma::mat &sol, const arma::uvec &swap,
            const arma::uvec &phi_m, const arma::uvec &phi_f,
            const arma::vec w) {
  arma::mat sol_eff = sol.rows(swap);
  arma::uvec swap_eff = {0, 1};
  double psi_sol = psi(sol_eff, phi_m, phi_f, w);
  sol_eff.submat(swap_eff, phi_f) =
      sol_eff.submat(arma::reverse(swap_eff), phi_f);
  double psi_prop = psi(sol_eff, phi_m, phi_f, w);
  return psi_prop - psi_sol;
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

  arma::vec psi_eval, dpsi_eval, rho_eval, temp_eval;

  if (eval) {
    psi_eval.set_size(n_iter);
    dpsi_eval.set_size(n_iter);
    rho_eval.set_size(n_iter);
    temp_eval.set_size(n_iter);

    psi_eval[0] = psi(sol, phi_m, phi_f, w);
    dpsi_eval[0] = 0.0;
    rho_eval[0] = 1.0;
    temp_eval[0] = temp0;
  }

  // Simulated annealing algorithm
  for (int i = 1; i < n_iter; i++) {
    arma::uvec swap =
        arma::randi<arma::uvec>(2, arma::distr_param(0, n_pairs - 1));

    double dpsi_prop = dpsi(sol, swap, phi_m, phi_f, w);
    double rho_prop = std::min(1.0, std::exp(-dpsi_prop / temp));
    double u = Rcpp::as<double>(Rcpp::runif(1));

    if (u < rho_prop) {
      sol.submat(swap, phi_f) = sol.submat(arma::reverse(swap), phi_f);
    }

    if (eval) {
      psi_eval[i] = psi(sol, phi_m, phi_f, w);
      dpsi_eval[i] = dpsi_prop;
      rho_eval[i] = rho_prop;
      temp_eval[i] = temp;
    }

    temp *= alpha;
  }

  if (eval) {
    return Rcpp::List::create(Rcpp::Named("sol") = sol,
                              Rcpp::Named("psi") = Rcpp::wrap(psi_eval),
                              Rcpp::Named("dpsi") = Rcpp::wrap(dpsi_eval),
                              Rcpp::Named("rho") = Rcpp::wrap(rho_eval),
                              Rcpp::Named("temp") = Rcpp::wrap(temp_eval));
  } else {
    return Rcpp::List::create(Rcpp::Named("sol") = sol);
  }
}
