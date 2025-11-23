// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <amsim/metric.h>
#include <amsim/metricspec.h>

//' Copy a C-style array into a std::vector
//'
//' @param ptr Pointer to the first element of the source array.
//' @param n_elem Number of elements to copy.
//'
//' @return A std::vector containing a copy of the array elements.
//'
//' @noRd
template <typename T>
std::vector<T> wrap_array(const T* ptr, const std::size_t n_elem) {
  std::vector<T> result(n_elem);
  std::copy_n(ptr, n_elem, result.data());
  return result;
}

//' Wrap a SimulationContext into an R-compatible list structure
//'
//' @details
//' Converts the internal C++ simulation context into a nested Rcpp::List
//' that can be passed to user-defined R metric functions. The resulting
//' list contains:
//' - n_individuals: Population size
//' - n_loci: Number of genetic loci
//' - n_phenotypes: Number of phenotypic traits
//' - n_sex: Number of individuals per sex
//' - genome: List with locus MAFs, means, and variances
//' - phenotypes: Named list of phenotype data (components, means, variances)
//' - mating: List containing mate matching state
//'
//' @param ctx Reference to the SimulationContext to wrap.
//'
//' @return An Rcpp::List containing the wrapped context data.
//'
//' @noRd
Rcpp::List wrap_context(const amsim::SimulationContext& ctx) {
  Rcpp::List genome = Rcpp::List::create(
      Rcpp::Named("locus_mafs") = ctx.genome.v_lmaf(),
      Rcpp::Named("locus_means") = ctx.genome.v_lmean(),
      Rcpp::Named("locus_vars") = ctx.genome.v_lvar());

  Rcpp::List phenotypes;
  for (const auto& pheno : ctx.phenotypes)
    phenotypes[pheno.name()] = Rcpp::List::create(
        Rcpp::Named("name") = pheno.name(),
        Rcpp::Named("genetic") =
            wrap_array<double>(pheno(amsim::ComponentType::GENETIC), ctx.n_ind),
        Rcpp::Named("environmental") = wrap_array<double>(
            pheno(amsim::ComponentType::ENVIRONMENTAL), ctx.n_ind),
        Rcpp::Named("vertical") = wrap_array<double>(
            pheno(amsim::ComponentType::VERTICAL), ctx.n_ind),
        Rcpp::Named("total") =
            wrap_array<double>(pheno(amsim::ComponentType::TOTAL), ctx.n_ind),
        Rcpp::Named("mean") = Rcpp::List::create(
            Rcpp::Named("genetic") =
                pheno.comp_mean(amsim::ComponentType::GENETIC),
            Rcpp::Named("environmental") =
                pheno.comp_mean(amsim::ComponentType::ENVIRONMENTAL),
            Rcpp::Named("vertical") =
                pheno.comp_mean(amsim::ComponentType::VERTICAL),
            Rcpp::Named("total") =
                pheno.comp_mean(amsim::ComponentType::TOTAL)),
        Rcpp::Named("var") = Rcpp::List::create(
            Rcpp::Named("genetic") =
                pheno.comp_var(amsim::ComponentType::GENETIC),
            Rcpp::Named("environmental") =
                pheno.comp_var(amsim::ComponentType::ENVIRONMENTAL),
            Rcpp::Named("vertical") =
                pheno.comp_var(amsim::ComponentType::VERTICAL),
            Rcpp::Named("total") =
                pheno.comp_var(amsim::ComponentType::TOTAL)));

  Rcpp::List mating =
      Rcpp::List::create(Rcpp::Named("mate_matching") = ctx.model.state);

  Rcpp::List r_ctx = Rcpp::List::create(
      Rcpp::Named("n_individuals") = ctx.n_ind,
      Rcpp::Named("n_loci") = ctx.n_loc,
      Rcpp::Named("n_phenotypes") = ctx.n_pheno,
      Rcpp::Named("n_sex") = ctx.n_sex,
      Rcpp::Named("genome") = genome,
      Rcpp::Named("phenotypes") = phenotypes,
      Rcpp::Named("mating") = mating);

  return r_ctx;
}

//' Wrap an R function into a thread-safe C++ MetricFunc
//'
//' @details
//' Creates a C++ lambda that wraps the provided R function, converting the
//' SimulationContext to an R list before invocation. A static mutex ensures
//' thread-safe access to the R runtime when running parallel simulations.
//'
//' @param r_func The R function to wrap. Must accept a list (simulation
//'   context) and return a numeric vector.
//'
//' @return A MetricFunc callable that invokes the R function with the
//'   wrapped context.
//'
//' @noRd
amsim::MetricFunc wrap_metric(const Rcpp::Function& r_func) {
  static std::mutex r_mutex;

  return [r_func](const amsim::SimulationContext& ctx) -> std::vector<double> {
    std::lock_guard<std::mutex> lock(r_mutex);
    Rcpp::List r_ctx = wrap_context(ctx);
    Rcpp::NumericVector result = r_func(r_ctx);
    return Rcpp::as<std::vector<double>>(result);
  };
}

//' Create a custom metric specification
//'
//' @description
//' Defines a custom metric using an R function that computes values from the
//' simulation context. This allows users to extend the built-in metrics with
//' arbitrary calculations.
//'
//' @details
//' The provided R function receives a list containing the current simulation
//' state, including genome statistics, phenotype values and summary statistics,
//' and mating information. The function must return a numeric vector of length
//' `n_rows * n_cols`.
//'
//' The context list passed to the metric function contains:
//' \describe{
//'   \item{n_individuals}{Population size}
//'   \item{n_loci}{Number of genetic loci}
//'   \item{n_phenotypes}{Number of phenotypic traits}
//'   \item{n_sex}{Number of individuals per sex}
//'   \item{genome}{List with locus_mafs, locus_means, locus_vars}
//'   \item{phenotypes}{Named list of phenotype data with component vectors
//'     (genetic, environmental, vertical, total) and summary statistics}
//'   \item{mating}{List containing mate_matching state}
//' }
//'
//' @param name String. Unique identifier for the metric.
//' @param r_metric_func Function. An R function that accepts the simulation
//'   context list and returns a numeric vector.
//' @param r_n_rows Integer. Number of rows in the metric output.
//' @param r_n_cols Integer or NULL. Number of columns in the metric output
//'   (default: 1).
//' @param r_labels Character vector or NULL. Labels for metric components
//'   (default: empty).
//' @param require_lat Logical. Whether the metric requires linkage
//'   disequilibrium tracking (default: FALSE).
//'
//' @return An external pointer to a MetricSpec object for use with
//'   Simulation$metrics().
//'
// [[Rcpp::export]]
SEXP custom_metric(
    const std::string& name,
    const Rcpp::Function& r_metric_func,
    SEXP r_n_rows,
    SEXP r_n_cols = R_NilValue,
    SEXP r_labels = R_NilValue,
    const bool require_lat = false) {
  std::size_t n_cols;
  std::vector<std::string> labels;
  if (r_n_cols == R_NilValue) n_cols = 1;
  else n_cols = Rcpp::as<std::size_t>(r_n_cols);

  if (r_labels == R_NilValue) labels = {};
  else labels = Rcpp::as<std::vector<std::string>>(r_labels);

  auto n_rows = Rcpp::as<std::size_t>(r_n_rows);

  amsim::MetricFunc metric_func = wrap_metric(r_metric_func);

  amsim::MetricSetup metric_setup =
      [name, metric_func, n_rows, n_cols, labels, require_lat](
          const amsim::SimulationContext& /*ctx*/) -> amsim::Metric {
    return amsim::Metric(
        metric_func, name, n_rows, n_cols, labels, require_lat);
  };

  amsim::MetricSpec spec(name, metric_func, metric_setup, require_lat);
  auto ptr_unique = std::make_unique<amsim::MetricSpec>(spec);

  Rcpp::XPtr<amsim::MetricSpec> ptr(ptr_unique.release(), true);

  return ptr;
}
