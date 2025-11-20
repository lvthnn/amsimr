// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <amsim/simulation.h>
#include <amsim/simulation_config.h>

#include "bindings_utils.h"

//' Create a new SimulationConfig object reference
//'
//' @return A new SimulationConfig instance
//'
// [[Rcpp::export(".SimulationConfig__new")]]
SEXP SimulationConfig__new() {
  auto config = std::make_unique<amsim::SimulationConfig>();
  Rcpp::XPtr<amsim::SimulationConfig> ptr(config.release(), true);
  return ptr;
}

//' Configure simulation parameters.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//' @param n_generations Integer. Number of generations to simulate.
//' @param n_individuals Integer. Number of individuals in the population.
//' @param output_dir String. Output directory to write data out to.
//' @param random_seed Integer. Seed used by random number generator.
//'
//' @return Self, for method chaining.
//'
// [[Rcpp::export(".SimulationConfig__simulation")]]
SEXP SimulationConfig__simulation(
    SEXP ptr,
    std::size_t n_generations,
    std::size_t n_individuals,
    std::string output_dir,
    SEXP random_seed) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  config->simulation(
      n_generations,
      n_individuals,
      output_dir,
      _optional<std::size_t>(random_seed));
  return config;
}

//' Configure genome parameters.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//' @param n_loci Number of genetic loci simulated in the population.
//' @param locus_maf Vector of locus MAFs
//' @param locus_recombination Vector of locus recombination probabilities
//' @param locus_mutation Vector of locus mutation probabilities
//'
//' @return Self, for method chaining
//'
// [[Rcpp::export(".SimulationConfig__genome")]]
SEXP SimulationConfig__genome(
    SEXP ptr,
    std::size_t n_loci,
    std::vector<double> locus_maf,
    std::vector<double> locus_recombination,
    std::vector<double> locus_mutation) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  config->genome(n_loci, locus_maf, locus_recombination, locus_mutation);
  return config;
}

//' Configure phenome parameters.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//' @param n_phenotypes Number of phenotypes modelled in the population.
//' @param names Vector of phenotype names.
//' @param n_causal_loci Number of causal loci for each phenotype.
//' @param h2_genetic Narrow-sense phenotype heritabilities.
//' @param h2_environmental Phenotype environmental component variances.
//' @param h2_vertical Phenotype vertical component variances.
//'
//' @return Self, for method chaining
//'
// [[Rcpp::export(".SimulationConfig__phenome")]]
SEXP SimulationConfig__phenome(
    SEXP ptr,
    std::size_t n_phenotypes,
    std::vector<std::string> names,
    std::vector<std::size_t> n_causal_loci,
    std::vector<double> h2_genetic,
    std::vector<double> h2_environmental,
    std::vector<double> h2_vertical,
    std::vector<double> genetic_cor,
    std::vector<double> environmental_cor) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  config->phenome(
      n_phenotypes,
      names,
      n_causal_loci,
      h2_genetic,
      h2_environmental,
      h2_vertical,
      genetic_cor,
      environmental_cor);
  return config;
}

//' Specify a random mating regime.
//'
//' @return Self, for method chaining.
//'
// [[Rcpp::export(".SimulationConfig__random_mating")]]
SEXP SimulationConfig__random_mating(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  config->random_mating();
  return config;
}

//' Specify an assortative mating regime.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//' @param mate_cor Column-major double trait cross-correlation matrix
//' @param tol_inf Termination threshold for annealing routine.
//' @param n_iterations Maximum number of steps performed in annealing routine.
//' @param temp_init Initial temperature of annealing routine
//' @param temp_decay Temperature decay parameter of annealing routine
//'
//' @return Self, for method chaining
//'
// [[Rcpp::export(".SimulationConfig__assortative_mating")]]
SEXP SimulationConfig__assortative_mating(
    SEXP ptr,
    std::vector<double> mate_cor,
    SEXP tol_inf = R_NilValue,
    SEXP n_iterations = R_NilValue,
    SEXP temp_init = R_NilValue,
    SEXP temp_decay = R_NilValue) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  config->assortative_mating(
      mate_cor,
      _optional<double>(tol_inf),
      _optional<std::size_t>(n_iterations),
      _optional<double>(temp_init),
      _optional<double>(temp_decay));
  return config;
}

//' Configure simulation metrics
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//' @param metrics List of metric objects to attach to the simulation.
//'
//' @return Self, for method chaining.
//'
// [[Rcpp::export(".SimulationConfig__metrics")]]
SEXP SimulationConfig__metrics(SEXP ptr, Rcpp::List metrics) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);

  std::size_t n_specs = metrics.size();
  std::vector<amsim::MetricSpec> metric_specs;
  metric_specs.reserve(n_specs);

  for (std::size_t sp = 0; sp < n_specs; ++sp) {
    Rcpp::XPtr<amsim::MetricSpec> spec(metrics[sp]);
    metric_specs.push_back(*spec);
  }

  config->metrics(metric_specs);
  return config;
}

//' Get the number of generations.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
// [[Rcpp::export(".SimulationConfig__n_generations")]]
std::size_t SimulationConfig__n_generations(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->n_gen;
}

//' Get the number of individuals.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
// [[Rcpp::export(".SimulationConfig__n_individuals")]]
std::size_t SimulationConfig__n_individuals(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->n_ind;
}

//' Get the number of loci.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
// [[Rcpp::export(".SimulationConfig__n_loci")]]
std::size_t SimulationConfig__n_loci(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->n_loc;
}

//' Get the output directory.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
// [[Rcpp::export(".SimulationConfig__output_dir")]]
std::string SimulationConfig__output_dir(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->out_dir;
}

//' Get the RNG seed.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
// [[Rcpp::export(".SimulationConfig__random_seed")]]
std::size_t SimulationConfig__random_seed(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->rng_seed;
}

//' Get the locus MAF vector.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
// [[Rcpp::export(".SimulationConfig__locus_maf")]]
std::vector<double> SimulationConfig__locus_maf(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->v_maf;
}

//' Get the locus recombination probabilities.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
// [[Rcpp::export(".SimulationConfig__locus_recombination")]]
std::vector<double> SimulationConfig__locus_recombination(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->v_rec;
}

//' Get the locus mutation probabilities.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
// [[Rcpp::export(".SimulationConfig__locus_mutation")]]
std::vector<double> SimulationConfig__locus_mutation(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->v_mut;
}

//' Get the number of phenotypes.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
// [[Rcpp::export(".SimulationConfig__n_phenotypes")]]
std::size_t SimulationConfig__n_phenotypes(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->n_pheno;
}

//' Get the phenotype names.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
// [[Rcpp::export(".SimulationConfig__phenotype_names")]]
std::vector<std::string> SimulationConfig__phenotype_names(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->v_name;
}

//' Get the phenotype names.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
// [[Rcpp::export(".SimulationConfig__n_causal_loci")]]
std::vector<std::size_t> SimulationConfig__n_causal_loci(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->v_n_loc;
}

//' Get the narrow-sense phenotype heritabilities.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
// [[Rcpp::export(".SimulationConfig__h2_genetic")]]
std::vector<double> SimulationConfig__h2_genetic(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->v_h2_gen;
}

//' Get the phenotype environmental component variances.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
// [[Rcpp::export(".SimulationConfig__h2_environmental")]]
std::vector<double> SimulationConfig__h2_environmental(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->v_h2_env;
}

//' Get the phenotype vertical component variances.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
// [[Rcpp::export(".SimulationConfig__h2_vertical")]]
std::vector<double> SimulationConfig__h2_vertical(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->v_h2_vert;
}

//' Get the phenotype genetic component correlation matrix.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
// [[Rcpp::export(".SimulationConfig__genetic_cor")]]
std::vector<double> SimulationConfig__genetic_cor(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->gen_cor;
}

//' Get the phenotype environmental component correlation matrix.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
// [[Rcpp::export(".SimulationConfig__environmental_cor")]]
std::vector<double> SimulationConfig__environmental_cor(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->env_cor;
}

//' Get the mating correlation matrix.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
// [[Rcpp::export(".SimulationConfig__mate_cor")]]
std::vector<double> SimulationConfig__mate_cor(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->mate_cor;
}

//' Get the annealing tolerance threshold.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
// [[Rcpp::export(".SimulationConfig__tol_inf")]]
double SimulationConfig__tol_inf(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->tol_inf;
}

//' Get the maximum number of annealing iterations.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
// [[Rcpp::export(".SimulationConfig__n_iterations")]]
std::size_t SimulationConfig__n_iterations(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->n_itr;
}

//' Get the initial annealing temperature.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
// [[Rcpp::export(".SimulationConfig__temp_init")]]
double SimulationConfig__temp_init(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->temp_init;
}

//' Get the annealing temperature decay parameter.
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
// [[Rcpp::export(".SimulationConfig__temp_decay")]]
double SimulationConfig__temp_decay(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->temp_decay;
}

//' Run (multithreaded) simulations
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//' @param n_replicates Integer. Number of independent simulation replicates to
//'   run.
//' @param n_threads Integer. Number of parallel threads to use for execution.
//' @param summarise Logical. Whether to summarize results across replicates
//'   (default: TRUE).
//' @param log_file Logical. Whether to write log output to a file (default:
//'   FALSE).
//' @param log_level String. Logging verbosity level.
//'
// [[Rcpp::export(".SimulationConfig__run_simulations")]]
void run_simulations(
    SEXP ptr,
    std::size_t n_replicates,
    std::size_t n_threads,
    bool summarise = true,
    bool log_file = false,
    std::string log_level = "info") {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  amsim::LogLevel level = _s_LogLevel(log_level);
  amsim::run_simulations(
      *config, n_replicates, n_threads, summarise, log_file, level);
}
