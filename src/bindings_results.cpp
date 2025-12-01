// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <amsim/simulation_results.h>

#include "bindings_utils.h"

//' Create a new SimulationResults object reference
//'
//' @return A new SimulationResults instance.
//'
//' @noRd
//'
// [[Rcpp::export]]
SEXP SimulationResults_new(const std::string& output_dir) {
  auto results = std::make_unique<amsim::SimulationResults>(output_dir);
  Rcpp::XPtr<amsim::SimulationResults> ptr(results.release(), true);
  return ptr;
}

//' Perform summarisation on replicate metric data
//'
//' @param ptr An external pointer to a SimulationResults instance.
//'
//' @noRd
//'
// [[Rcpp::export]]
void SimulationResults_summarise(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationResults> results(ptr);
  results->summarise();
}

//' Save summarised results to disk
//'
//' @param ptr An external pointer to a SimulationResults instance.
//' @param metrics A list of metric names to be saved to file.
//' @param out_dir Output directory to write data out to.
//' @param overwrite Specifies whether directory should be overwritten.
//'
//' @noRd
//'
// [[Rcpp::export]]
void SimulationResults_save(
    SEXP ptr,
    SEXP metrics = R_NilValue,
    SEXP out_dir = R_NilValue,
    bool overwrite = false) {
  Rcpp::XPtr<amsim::SimulationResults> results(ptr);
  results->save(
      rOptional<std::vector<std::string>>(metrics),
      rOptional<std::string>(out_dir),
      overwrite);
}

//' Load a summarised metric into the binding layer.
//'
//' @param ptr An external pointer to a SimulationResults instance.
//' @param metric_name Name of the metric to be loaded into the environment.
//'
//' @return A DataFrame comprising the summarised results data for the specified
//'   metric.
//'
//' @noRd
//'
// [[Rcpp::export]]
Rcpp::DataFrame SimulationResults_load(
    SEXP ptr, const std::string& metric_name) {
  Rcpp::XPtr<amsim::SimulationResults> results(ptr);
  Rcpp::DataFrame metric_results;

  amsim::ResultsTable table = (*results)(metric_name);
  std::size_t n_rows = table.data[0].size() / table.label_map.size();
  std::size_t n_cols = table.label_map.size();
  std::vector<std::size_t> index(n_rows * n_cols);
  std::vector<std::string> names(n_rows * n_cols);

  for (std::size_t col = 0; col < n_cols; ++col) {
    for (std::size_t row = 0; row < n_rows; ++row) {
      std::size_t idx = (row * n_cols) + col;
      index[idx] = row + 1;
      names[idx] = table.inv_label_map[col];
    }
  }

  metric_results["gen"] = index;
  metric_results["name"] = names;
  metric_results["mean"] = table.data[0];
  metric_results["median"] = table.data[1];
  metric_results["stddev"] = table.data[2];
  metric_results["stderr"] = table.data[3];
  metric_results["lower_ci95"] = table.data[4];
  metric_results["upper_ci95"] = table.data[5];
  metric_results["quant_025"] = table.data[6];
  metric_results["quant_975"] = table.data[7];

  return metric_results;
}

//' Get metric names
//'
//' @param ptr An external pointer to a SimulationResults instance.
//'
//' @return A vector of the metric names recognised by the results instance.
//'
//' @noRd
//'
// [[Rcpp::export]]
std::vector<std::string> SimulationResults_metric_names(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationResults> results(ptr);
  return results->metric_names();
}
