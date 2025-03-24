// [[depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Computes correlation update when swapping pairs
//'
//' @param sol_mat Solution matrix
//' @param swap_idx Indices to swap
//' @param male_snp_idx Male SNP indices
//' @param female_snp_idx Female SNP indices
//' @return Correlation vector update
//[[Rcpp::export]]
vec compute_delta(
  const mat &sol_mat,
  const uvec &swap_idx,
  const uvec &male_snp_idx,
  const uvec &female_snp_idx
) {
  const int pop_size = sol_mat.n_rows;
  vec cor_update = vec(male_snp_idx.n_elem);

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
//' @param curr_cor Current correlation
//' @param target_cor Target correlation
//' @param delta_cor Correlation update
//' @return Energy differential
//[[Rcpp::export]]
double compute_dpsi(
  const vec &curr_cor,
  const vec &target_cor,
  const vec &delta_cor
) {
  const vec diff_cor = curr_cor - target_cor;
  return dot(delta_cor, delta_cor) + 2.0 * dot(diff_cor, delta_cor);
}

//' Evaluates energy of current state
//'
//' @param curr_cor Current correlation
//' @param target_cor Target correlation
//' @return Energy value
//[[Rcpp::export]]
double compute_psi(
  const vec &curr_cor,
  const vec &target_cor
) {
  return std::pow(norm(curr_cor - target_cor, 2), 2);
}

//' Performs simulated annealing to optimize matching
//'
//' @param sol_mat Solution matrix to optimize
//' @param snp_pairs SNP pair data [male_idx, female_idx, target_cor]
//' @param female_swap_idx Columns to swap in solution
//' @param max_iterations Maximum iterations to run
//' @param temp_decay Temperature decay rate
//' @param init_temp Initial temperature
//' @param collect_metrics Whether to collect optimization metrics
//' @return Optimized solution matrix and optional metrics
// [[Rcpp::export]]
List optim_matching(
  mat &sol_mat,
  const mat &snp_pairs,
  const uvec &female_swap_idx,
  const int num_iterations = 10000,
  const double temp_decay = 1.0,
  const double init_temp = 1.0,
  const bool collect_metrics = false
) {
  const int pop_size = sol_mat.n_rows;

  // Extract SNP genotype indices and target correlation
  const uvec male_snp_idx = conv_to<uvec>::from(snp_pairs.col(0));
  const uvec female_snp_idx = conv_to<uvec>::from(snp_pairs.col(1));
  const vec target_correlation = snp_pairs.col(2);

  // Distribution parameters for generation of swap indices
  const distr_param swap_distr = distr_param(0, pop_size - 1);

  double curr_temp = init_temp;

  // Vectors for storing evaluation metrics if requested
  vec energy_vals, energy_deltas, acc_probs, temp_vals;

  // Compute the assortative SNP correlation vector for the initial state
  vec curr_cor = vec(male_snp_idx.n_elem);
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
    uvec swap_idx = randi<uvec>(2, swap_distr);
    vec delta_cor = compute_delta(
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
    vec prop_cor = curr_cor + delta_cor;
    double acc_prob = std::min(1.0, std::exp(-energy_delta / curr_temp));
    double unf = as<double>(runif(1));

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
    return List::create(
      Named("sol_mat") = sol_mat,
      Named("energy_vals") = wrap(energy_vals),
      Named("energy_deltas") = wrap(energy_deltas),
      Named("acceptance_probs") = wrap(acc_probs),
      Named("temp_vals") = wrap(temp_vals)
    );
  } else {
    return List::create(Named("sol_mat") = sol_mat);
  }
}
