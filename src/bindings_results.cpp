// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <amsim/simulation_results.h>

#include "bindings_utils.h"

SEXP SimulationResults__new(const std::string& output_dir) {
  auto results = std::make_unique<amsim::SimulationResults>(output_dir);
  Rcpp::XPtr<amsim::SimulationResults> ptr(results.release(), true);
  return ptr;
}

void SimulationResults__summarise(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationResults> results(ptr);
  results->summarise();
}

void SimulationResults__save(
    SEXP ptr,
    SEXP metrics = R_NilValue,
    SEXP out_dir = R_NilValue,
    bool overwrite = false) {
  Rcpp::XPtr<amsim::SimulationResults> results(ptr);
  results->save(
      _optional<std::vector<std::string>>(metrics),
      _optional<std::string>(out_dir),
      overwrite);
}

Rcpp::DataFrame SimulationResults__load(
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
      std::size_t idx = row * n_cols + col;
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

std::vector<std::string> SimulationResults__metric_names(SEXP ptr) {
  Rcpp::XPtr<amsim::SimulationResults> results(ptr);
  return results->metric_names();
}
