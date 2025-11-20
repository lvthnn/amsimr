// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <amsim/metricspec.h>

#include "bindings_utils.h"

// [[Rcpp::export]]
SEXP pheno_h2() {
  amsim::MetricSpec pheno_h2 = amsim::pheno_h2();
  auto metric = std::make_unique<amsim::MetricSpec>(pheno_h2);
  Rcpp::XPtr<amsim::MetricSpec> ptr(metric.release(), true);
  return ptr;
}

// [[Rcpp::export]]
SEXP pheno_comp_mean(const std::string& component_type) {
  amsim::ComponentType type = _s_ComponentType(component_type);
  amsim::MetricSpec pheno_comp_mean = amsim::pheno_comp_mean(type);
  auto metric = std::make_unique<amsim::MetricSpec>(pheno_comp_mean);
  Rcpp::XPtr<amsim::MetricSpec> ptr(metric.release(), true);
  return ptr;
}

// [[Rcpp::export]]
SEXP pheno_comp_var(const std::string& component_type) {
  amsim::ComponentType type = _s_ComponentType(component_type);
  amsim::MetricSpec pheno_comp_var = amsim::pheno_comp_var(type);
  auto metric = std::make_unique<amsim::MetricSpec>(pheno_comp_var);
  Rcpp::XPtr<amsim::MetricSpec> ptr(metric.release(), true);
  return ptr;
}

// [[Rcpp::export]]
SEXP pheno_comp_cor(const std::string& component_type) {
  amsim::ComponentType type = _s_ComponentType(component_type);
  amsim::MetricSpec pheno_comp_cor = amsim::pheno_comp_cor(type);
  auto metric = std::make_unique<amsim::MetricSpec>(pheno_comp_cor);
  Rcpp::XPtr<amsim::MetricSpec> ptr(metric.release(), true);
  return ptr;
}

// [[Rcpp::export]]
SEXP pheno_comp_xcor(const std::string& component_type) {
  amsim::ComponentType type = _s_ComponentType(component_type);
  amsim::MetricSpec pheno_comp_xcor = amsim::pheno_comp_xcor(type);
  auto metric = std::make_unique<amsim::MetricSpec>(pheno_comp_xcor);
  Rcpp::XPtr<amsim::MetricSpec> ptr(metric.release(), true);
  return ptr;
}

// [[Rcpp::export]]
SEXP pheno_latent_h2() {
  amsim::MetricSpec pheno_latent_h2 = amsim::pheno_latent_h2();
  auto metric = std::make_unique<amsim::MetricSpec>(pheno_latent_h2);
  Rcpp::XPtr<amsim::MetricSpec> ptr(metric.release(), true);
  return ptr;
}

// [[Rcpp::export]]
SEXP pheno_latent_comp_mean(const std::string& component_type) {
  amsim::ComponentType type = _s_ComponentType(component_type);
  amsim::MetricSpec pheno_comp_mean = amsim::pheno_latent_comp_mean(type);
  auto metric = std::make_unique<amsim::MetricSpec>(pheno_comp_mean);
  Rcpp::XPtr<amsim::MetricSpec> ptr(metric.release(), true);
  return ptr;
}

// [[Rcpp::export]]
SEXP pheno_latent_comp_var(const std::string& component_type) {
  amsim::ComponentType type = _s_ComponentType(component_type);
  amsim::MetricSpec pheno_comp_var = amsim::pheno_latent_comp_var(type);
  auto metric = std::make_unique<amsim::MetricSpec>(pheno_comp_var);
  Rcpp::XPtr<amsim::MetricSpec> ptr(metric.release(), true);
  return ptr;
}

// [[Rcpp::export]]
SEXP pheno_latent_comp_cor(const std::string& component_type) {
  amsim::ComponentType type = _s_ComponentType(component_type);
  amsim::MetricSpec pheno_comp_cor = amsim::pheno_latent_comp_cor(type);
  auto metric = std::make_unique<amsim::MetricSpec>(pheno_comp_cor);
  Rcpp::XPtr<amsim::MetricSpec> ptr(metric.release(), true);
  return ptr;
}

// [[Rcpp::export]]
SEXP pheno_latent_comp_xcor(const std::string& component_type) {
  amsim::ComponentType type = _s_ComponentType(component_type);
  amsim::MetricSpec pheno_comp_xcor = amsim::pheno_latent_comp_xcor(type);
  auto metric = std::make_unique<amsim::MetricSpec>(pheno_comp_xcor);
  Rcpp::XPtr<amsim::MetricSpec> ptr(metric.release(), true);
  return ptr;
}
