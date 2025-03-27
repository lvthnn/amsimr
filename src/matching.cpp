// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Computes correlation update when swapping pairs
//'
//' @name compute_delta
//' @param sol_mat Solution matrix
//' @param swap_idx Indices to swap
//' @param male_snp_idx Male SNP indices
//' @param female_snp_idx Female SNP indices
//' @return Correlation vector update
//' @keywords internal
arma::vec compute_delta(
  const arma::mat &sol_mat,
  const arma::uvec &swap_idx,
  const arma::uvec &male_snp_idx,
  const arma::uvec &female_snp_idx
) {
  const int pop_size = sol_mat.n_rows;
  arma::vec cor_update = arma::vec(male_snp_idx.n_elem);

  for (size_t i = 0; i < male_snp_idx.n_elem; i++) {
    double m0 = sol_mat(swap_idx(0), male_snp_idx[i]);
    double m1 = sol_mat(swap_idx(1), male_snp_idx[i]);
    double f0 = sol_mat(swap_idx(0), female_snp_idx[i]);
    double f1 = sol_mat(swap_idx(1), female_snp_idx[i]);

    cor_update[i] = (1.0 / pop_size) * (m0 * (f1 - f0) + m1 * (f0 - f1));
  }

  return cor_update;
}

//' Computes energy differential between states
//'
//' @name compute_dpsi
//' @param curr_cor Current correlation
//' @param target_cor Target correlation
//' @param delta_cor Correlation update
//' @return Energy differential
//' @keywords internal
double compute_dpsi(
  const arma::vec &curr_cor,
  const arma::vec &target_cor,
  const arma::vec &delta_cor
) {
  const arma::vec diff_cor = curr_cor - target_cor;
  return arma::dot(delta_cor, delta_cor) + 2.0 * arma::dot(diff_cor, delta_cor);
}

//' Evaluates energy of current state
//'
//' @name compute_psi
//' @param curr_cor Current correlation
//' @param target_cor Target correlation
//' @return Energy value
//' @keywords internal
double compute_psi(
  const arma::vec &curr_cor,
  const arma::vec &target_cor
) {
  return std::pow(arma::norm(curr_cor - target_cor, 2), 2);
}

//' Automatically determine initial temperature for simulated annealing
//'
//' @name auto_init_temp
//' @param sol_mat Current solution matrix
//' @param snp_pairs SNP pair data (male_idx, female_idx, target_cor)
//' @param female_swap_idx Columns to swap in solution
//' @param num_samples Number of random moves to sample
//' @param accept_ratio Target initial acceptance ratio (default: 0.8)
//' @return Recommended initial temperature
//' @keywords internal
double auto_init_temp(
  const arma::mat &sol_mat,
  const arma::mat &snp_pairs,
  const arma::uvec &female_swap_idx,
  const int num_samples,
  const double accept_ratio
) {
  const int pop_size = sol_mat.n_rows;
  
  // Extract SNP indices and target correlation
  const arma::uvec male_snp_idx =
    arma::conv_to<arma::uvec>::from(snp_pairs.col(0));
  const arma::uvec female_snp_idx =
    arma::conv_to<arma::uvec>::from(snp_pairs.col(1));
  const arma::vec target_correlation = snp_pairs.col(2);
  
  // Distribution parameters for swap indices
  const arma::distr_param swap_distr = arma::distr_param(0, pop_size - 1);
  
  // Compute current correlation vector
  arma::vec curr_cor = arma::vec(male_snp_idx.n_elem);
  for (size_t i = 0; i < male_snp_idx.n_elem; i++) {
    curr_cor[i] = (1.0 / pop_size) * dot(
      sol_mat.col(male_snp_idx[i]),
      sol_mat.col(female_snp_idx[i])
    );
  }
  
  // Sample random moves and collect energy deltas
  arma::vec energy_deltas = arma::vec(num_samples);
  
  for (int i = 0; i < num_samples; i++) {
    arma::uvec swap_idx = arma::randi<arma::uvec>(2, swap_distr);
    arma::vec delta_cor = compute_delta(
      sol_mat,
      swap_idx,
      male_snp_idx,
      female_snp_idx
    );
    energy_deltas[i] = compute_dpsi(
      curr_cor,
      target_correlation,
      delta_cor
    );
  }
  
  // Keep only positive energy deltas (uphill moves)
  arma::vec pos_deltas = energy_deltas.elem(arma::find(energy_deltas > 0));
  
  // If there are no positive deltas, return a small default value
  if (pos_deltas.n_elem == 0) return 1e-9;

  double median_delta = arma::median(pos_deltas);
  return -median_delta / std::log(accept_ratio);
}

//' Performs simulated annealing to optimize matching
//'
//' @param sol_mat Solution matrix to optimize
//' @param snp_pairs SNP pair data [male_idx, female_idx, target_cor]
//' @param female_swap_idx Columns to swap in solution
//' @param max_iterations Maximum iterations to run
//' @param temp_decay Temperature decay rate
//' @param init_temp Initial temperature
//' @param auto_temp_samples Number of samples to draw for determining initial
//'   temperature of simulated annealing algorithm
//' @param auto_accept_ratio Desired initial acceptance ratio for calibration of
//'   initial temperature value in simulated annealing algorithm
//' @param collect_metrics Whether to collect optimization metrics
//' @param quietly Print diagnostics while running annealing
//' @return Optimized solution matrix and optional metrics
// [[Rcpp::export]]
Rcpp::List optim_matching(
  arma::mat &sol_mat,
  const arma::mat &snp_pairs,
  const arma::uvec &female_swap_idx,
  const int num_iterations = 10000,
  const double temp_decay = 0.995,
  const double init_temp = 1e-9,
  const int auto_temp_samples = 100000,
  const double auto_accept_ratio = 0.995,
  const bool collect_metrics = false,
  const bool quietly = true 
) {
  const int pop_size = sol_mat.n_rows;

  // Extract SNP genotype indices and target correlation
  const arma::uvec male_snp_idx =
    arma::conv_to<arma::uvec>::from(snp_pairs.col(0));
  const arma::uvec female_snp_idx =
    arma::conv_to<arma::uvec>::from(snp_pairs.col(1));
  const arma::vec target_correlation = snp_pairs.col(2);

  // Distribution parameters for generation of swap indices
  const arma::distr_param swap_distr = arma::distr_param(0, pop_size - 1);
  
  double curr_temp;
  if (init_temp <= 0.0) {
    curr_temp = auto_init_temp(
      sol_mat, 
      snp_pairs, 
      female_swap_idx, 
      auto_temp_samples, 
      auto_accept_ratio
    );
    if (!quietly) {
      Rcpp::Rcout << "Auto-determined initial temperature: "
                  << curr_temp << "\n";
    }
  } else {
    curr_temp = init_temp;
  }

  // Vectors for storing evaluation metrics if requested
  arma::vec energy_vals, energy_deltas, acc_probs, temp_vals;

  // Compute the assortative SNP correlation vector for the initial state
  arma::vec curr_cor = arma::vec(male_snp_idx.n_elem);
  for (size_t i = 0; i < male_snp_idx.n_elem; i++) {
    curr_cor[i] = (1.0 / pop_size) * dot(
      sol_mat.col(male_snp_idx[i]),
      sol_mat.col(female_snp_idx[i])
    );
  }

  // If metric collection is requested, initialize metric vectors
  if (collect_metrics) {
    energy_vals.set_size(num_iterations);
    energy_deltas.set_size(num_iterations);
    acc_probs.set_size(num_iterations);
    temp_vals.set_size(num_iterations);

    energy_vals[0] = compute_psi(curr_cor, target_correlation);
    energy_deltas[0] = 0.0;
    acc_probs[0] = 1.0;
    temp_vals[0] = curr_temp;
  }

  // Run the simulated annealing algorithm
  for (int i = 0; i < num_iterations; i++) {
    arma::uvec swap_idx = arma::randi<arma::uvec>(2, swap_distr);
    arma::vec delta_cor = compute_delta(
      sol_mat,
      swap_idx,
      male_snp_idx,
      female_snp_idx
    );
    double energy_delta = compute_dpsi(
      curr_cor,
      target_correlation,
      delta_cor
    );
    arma::vec prop_cor = curr_cor + delta_cor;
    double acc_prob = std::min(1.0, std::exp(-energy_delta / curr_temp));
    double unf = Rcpp::as<double>(Rcpp::runif(1));

    if (unf < acc_prob) {
      sol_mat.submat(swap_idx, female_swap_idx) = sol_mat.submat(
        reverse(swap_idx), female_swap_idx
      );
      curr_cor = prop_cor;
    }

    curr_temp *= temp_decay;

    if (collect_metrics) {
      energy_vals[i] = compute_psi(curr_cor, target_correlation);
      energy_deltas[i] = energy_delta;
      acc_probs[i] = acc_prob;
      temp_vals[i] = curr_temp;
    }
  }

  if (collect_metrics) {
    return Rcpp::List::create(
      Rcpp::Named("sol_mat") = sol_mat,
      Rcpp::Named("energy_vals") = Rcpp::wrap(energy_vals),
      Rcpp::Named("energy_deltas") = Rcpp::wrap(energy_deltas),
      Rcpp::Named("acceptance_probs") = Rcpp::wrap(acc_probs),
      Rcpp::Named("temp_vals") = Rcpp::wrap(temp_vals)
    );
  } else {
    return Rcpp::List::create(Rcpp::Named("sol_mat") = sol_mat);
  }
}
