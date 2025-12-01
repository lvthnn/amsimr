// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <amsim/simulation.h>
#include <amsim/simulation_config.h>

#include <utility>
#include <cstdint>
#include <filesystem>

#include "bindings_utils.h"

//' Create a new SimulationConfig object reference
//'
//' @return A new SimulationConfig instance.
//'
//' @noRd
//'
// [[Rcpp::export]]
SEXP SimulationConfig_new() {
  auto config = std::make_unique<amsim::SimulationConfig>();
  Rcpp::XPtr<amsim::SimulationConfig> ptr(config.release(), true);
  return ptr;
}

//' Configure simulation parameters
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//' @param n_generations Integer. Number of generations to simulate.
//' @param n_individuals Integer. Number of individuals in the population.
//' @param output_dir String. Output directory to write data out to.
//' @param random_seed Integer. Seed used by random number generator.
//'
//' @return Self, for method chaining.
//'
//' @noRd
//'
// [[Rcpp::export]]
SEXP SimulationConfig_simulation(
    SEXP ptr,
    std::size_t n_generations,
    std::size_t n_individuals,
    const std::string& output_dir,
    SEXP random_seed) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  config->simulation(
      n_generations,
      n_individuals,
      output_dir,
      rOptional<std::size_t>(random_seed));
  return config;
}

//' Configure genome parameters
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//' @param n_loci Number of genetic loci simulated in the population.
//' @param locus_maf Vector of locus MAFs.
//' @param locus_recombination Vector of locus recombination probabilities.
//' @param locus_mutation Vector of locus mutation probabilities.
//'
//' @return Self, for method chaining.
//'
//' @noRd
//'
// [[Rcpp::export]]
SEXP SimulationConfig_genome(
    SEXP ptr,
    std::size_t n_loci,
    std::vector<double> locus_maf,
    std::vector<double> locus_recombination,
    std::vector<double> locus_mutation) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  config->genome(
      n_loci,
      std::move(locus_maf),
      std::move(locus_recombination),
      std::move(locus_mutation));
  return config;
}

//' Configure phenome parameters
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//' @param n_phenotypes Number of phenotypes modelled in the population.
//' @param names Vector of phenotype names.
//' @param n_causal_loci Number of causal loci for each phenotype.
//' @param h2_genetic Narrow-sense phenotype heritabilities.
//' @param h2_environmental Phenotype environmental component variances.
//' @param h2_vertical Phenotype vertical component variances.
//' @param genetic_cor Genetic correlation matrix (column-major).
//' @param environmental_cor Environmental correlation matrix (column-major).
//'
//' @return Self, for method chaining.
//'
//' @noRd
//'
// [[Rcpp::export]]
SEXP SimulationConfig_phenome(
    SEXP ptr,
    std::size_t n_phenotypes,
    const std::vector<std::string>& names,
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
      std::move(n_causal_loci),
      std::move(h2_genetic),
      std::move(h2_environmental),
      std::move(h2_vertical),
      std::move(genetic_cor),
      std::move(environmental_cor));
  return config;
}

//' Specify a random mating regime
//'
//' @return Self, for method chaining.
//'
//' @noRd
//'
// [[Rcpp::export]]
SEXP SimulationConfig_random_mating(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  config->random_mating();
  return config;
}

//' Specify an assortative mating regime
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//' @param mate_cor Column-major double trait cross-correlation matrix.
//' @param tol_inf Termination threshold for annealing routine.
//' @param n_iterations Maximum number of steps performed in annealing routine.
//' @param temp_init Initial temperature of annealing routine.
//' @param temp_decay Temperature decay parameter of annealing routine.
//'
//' @return Self, for method chaining
//'
//' @noRd
//'
// [[Rcpp::export]]
SEXP SimulationConfig_assortative_mating(
    SEXP ptr,
    std::vector<double> mate_cor,
    SEXP tol_inf = R_NilValue,
    SEXP n_iterations = R_NilValue,
    SEXP temp_init = R_NilValue,
    SEXP temp_decay = R_NilValue) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  config->assortative_mating(
      std::move(mate_cor),
      rOptional<double>(tol_inf),
      rOptional<std::size_t>(n_iterations),
      rOptional<double>(temp_init),
      rOptional<double>(temp_decay));
  return config;
}

//' Configure simulation metrics
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//' @param metrics List of metric objects to attach to the simulation.
//'
//' @return Self, for method chaining.
//'
//' @noRd
//'
// [[Rcpp::export]]
SEXP SimulationConfig_metrics(SEXP ptr, Rcpp::List metrics) {
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

//' Get the number of generations
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
//' @return Integer specifying the number of generations to simulate.
//'
//' @noRd
//'
// [[Rcpp::export]]
std::size_t SimulationConfig_n_generations(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->n_gen;
}

//' Get the number of individuals
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
//' @return Integer specifying the population size.
//'
//' @noRd
//'
// [[Rcpp::export]]
std::size_t SimulationConfig_n_individuals(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->n_ind;
}

//' Get the number of loci
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
//' @return Integer specifying the number of genetic loci.
//'
//' @noRd
//'
// [[Rcpp::export]]
std::size_t SimulationConfig_n_loci(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->n_loc;
}

//' Get the output directory
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
//' @return String specifying the output directory path.
//'
//' @noRd
//'
// [[Rcpp::export]]
std::string SimulationConfig_output_dir(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->out_dir;
}

//' Get the RNG seed
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
//' @return Integer specifying the random number generator seed.
//'
//' @noRd
//'
// [[Rcpp::export]]
std::size_t SimulationConfig_random_seed(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->rng_seed;
}

//' Shuffle RNG seed for replicate simulations
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//' @param rep_id New seed 
//'
//' @noRd
//'
// [[Rcpp::export]]
void SimulationConfig_shuffle_random_seed(SEXP ptr, std::size_t rep_id) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  config->rng_seed = amsim::shuffle_seed(config->rng_seed, rep_id);
}


//' Get the locus MAF vector
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
//' @return Numeric vector of minor allele frequencies for each locus.
//'
//' @noRd
//'
// [[Rcpp::export]]
std::vector<double> SimulationConfig_locus_maf(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->v_maf;
}

//' Get the locus recombination probabilities
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
//' @return Numeric vector of recombination probabilities for each locus.
//'
//' @noRd
//'
// [[Rcpp::export]]
std::vector<double> SimulationConfig_locus_recombination(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->v_rec;
}

//' Get the locus mutation probabilities
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
//' @return Numeric vector of mutation probabilities for each locus.
//'
//' @noRd
//'
// [[Rcpp::export]]
std::vector<double> SimulationConfig_locus_mutation(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->v_mut;
}

//' Get the number of phenotypes
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
//' @return Integer specifying the number of phenotypic traits.
//'
//' @noRd
//'
// [[Rcpp::export]]
std::size_t SimulationConfig_n_phenotypes(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->n_pheno;
}

//' Get the phenotype names
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
//' @return Character vector of phenotype names.
//'
//' @noRd
//'
// [[Rcpp::export]]
std::vector<std::string> SimulationConfig_phenotype_names(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->v_name;
}

//' Get the number of causal loci for each phenotype
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
//' @return Integer vector specifying the number of causal loci per phenotype.
//'
//' @noRd
//'
// [[Rcpp::export]]
std::vector<std::size_t> SimulationConfig_n_causal_loci(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->v_n_loc;
}

//' Get the narrow-sense phenotype heritabilities
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
//' @return Numeric vector of narrow-sense heritabilities for each phenotype.
//'
//' @noRd
//'
// [[Rcpp::export]]
std::vector<double> SimulationConfig_h2_genetic(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->v_h2_gen;
}

//' Get the phenotype environmental component variances
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
//' @return Numeric vector of environmental variance components for each
//'   phenotype.
//'
//' @noRd
//'
// [[Rcpp::export]]
std::vector<double> SimulationConfig_h2_environmental(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->v_h2_env;
}

//' Get the phenotype vertical component variances
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
//' @return Numeric vector of vertical transmission components for each
//'   phenotype.
//'
//' @noRd
//'
// [[Rcpp::export]]
std::vector<double> SimulationConfig_h2_vertical(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->v_h2_vert;
}

//' Get the phenotype genetic component correlation matrix
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
//' @return Numeric vector representing the genetic correlation matrix in
//'   column-major order.
//'
//' @noRd
//'
// [[Rcpp::export]]
std::vector<double> SimulationConfig_genetic_cor(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->gen_cor;
}

//' Get the phenotype environmental component correlation matrix
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
//' @return Numeric vector representing the environmental correlation matrix in
//'   column-major order.
//'
//' @noRd
//'
// [[Rcpp::export]]
std::vector<double> SimulationConfig_environmental_cor(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->env_cor;
}

//' Get the mating correlation matrix
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
//' @return Numeric vector representing the mating correlation matrix in
//'   column-major order.
//'
//' @noRd
//'
// [[Rcpp::export]]
std::vector<double> SimulationConfig_mate_cor(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->mate_cor;
}

//' Get the annealing tolerance threshold
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
//' @return Numeric value specifying the termination threshold for the
//'   annealing routine.
//'
//' @noRd
//'
// [[Rcpp::export]]
double SimulationConfig_tol_inf(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->tol_inf;
}

//' Get the maximum number of annealing iterations
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
//' @return Integer specifying the maximum number of annealing iterations.
//'
//' @noRd
//'
// [[Rcpp::export]]
std::size_t SimulationConfig_n_iterations(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->n_itr;
}

//' Get the initial annealing temperature
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
//' @return Numeric value specifying the initial annealing temperature.
//'
//' @noRd
//'
// [[Rcpp::export]]
double SimulationConfig_temp_init(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->temp_init;
}

//' Get the annealing temperature decay parameter
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//'
//' @return Numeric value specifying the temperature decay rate.
//'
//' @noRd
//'
// [[Rcpp::export]]
double SimulationConfig_temp_decay(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  return config->temp_decay;
}


//' Run a single simulation
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//' @param output_dir Output directory to overwrite configuration.
//' @param ransom_seed Random seed to overwrite configuration
//' @param log_file Whether to write log output to file.
//' @param log_level Logging verbosity level.
//'
//' @noRd
//'
// [[Rcpp::export]]
void run_simulation(
    SEXP ptr,
    SEXP output_dir,
    bool log_file,
    const std::string& log_level) {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  amsim::LogLevel level = strLogLevel(log_level);
  amsim::run_simulation(
      *config,
      rOptional<std::string>(output_dir),
      log_file,
      level);
}


//' Run (multithreaded) simulations
//'
//' @param ptr An external pointer to a SimulationConfig instance.
//' @param n_replicates Number of independent simulation replicates to
//'   run.
//' @param n_threads Number of parallel threads to use for execution.
//' @param summarise Whether to summarize results across replicates.
//' @param log_file Whether to write log output to a file.
//' @param log_level String. Logging verbosity level.
//'
//' @noRd
//'
// [[Rcpp::export]]
void run_simulations(
    SEXP ptr,
    std::size_t n_replicates,
    std::size_t n_threads,
    bool summarise = true,
    bool log_file = false,
    const std::string& log_level = "info") {
  Rcpp::XPtr<amsim::SimulationConfig> config(ptr);
  amsim::LogLevel level = strLogLevel(log_level);
  amsim::run_simulations(
      *config, n_replicates, n_threads, summarise, log_file, level);
}
