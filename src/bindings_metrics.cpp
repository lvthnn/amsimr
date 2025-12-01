// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <amsim/metricspec.h>

#include "bindings_utils.h"

//' Phenotype heritability metric
//'
//' @description
//' Creates a metric specification for tracking narrow-sense heritability
//' (h²) of each phenotype across generations.
//'
//' @return An external pointer to a MetricSpec object for use with
//'   Simulation$metrics().
//'
//' @noRd
//'
// [[Rcpp::export]]
SEXP MetricSpec_pheno_h2() {
  amsim::MetricSpec pheno_h2 = amsim::pheno_h2();
  auto metric = std::make_unique<amsim::MetricSpec>(pheno_h2);
  Rcpp::XPtr<amsim::MetricSpec> ptr(metric.release(), true);
  return ptr;
}

//' Phenotype component mean metric
//'
//' @description
//' Creates a metric specification for tracking the mean of a specific
//' phenotype component across generations.
//'
//' @param component_type String. The phenotype component to track. One of
//'   "genetic", "environmental", "vertical", or "total".
//'
//' @return An external pointer to a MetricSpec object for use with
//'   Simulation$metrics().
//'
//' @noRd
//'
// [[Rcpp::export]]
SEXP MetricSpec_pheno_comp_mean(const std::string& component_type) {
  amsim::ComponentType type = strComponentType(component_type);
  amsim::MetricSpec pheno_comp_mean = amsim::pheno_comp_mean(type);
  auto metric = std::make_unique<amsim::MetricSpec>(pheno_comp_mean);
  Rcpp::XPtr<amsim::MetricSpec> ptr(metric.release(), true);
  return ptr;
}

//' Phenotype component variance metric
//'
//' @description
//' Creates a metric specification for tracking the variance of a specific
//' phenotype component across generations.
//'
//' @param component_type String. The phenotype component to track. One of
//'   "genetic", "environmental", "vertical", or "total".
//'
//' @return An external pointer to a MetricSpec object for use with
//'   Simulation$metrics().
//'
//' @noRd
//'
// [[Rcpp::export]]
SEXP MetricSpec_pheno_comp_var(const std::string& component_type) {
  amsim::ComponentType type = strComponentType(component_type);
  amsim::MetricSpec pheno_comp_var = amsim::pheno_comp_var(type);
  auto metric = std::make_unique<amsim::MetricSpec>(pheno_comp_var);
  Rcpp::XPtr<amsim::MetricSpec> ptr(metric.release(), true);
  return ptr;
}

//' Phenotype component correlation metric
//'
//' @description
//' Creates a metric specification for tracking within-phenotype correlations
//' of a specific component across generations. This measures the correlation
//' matrix of the specified component across all phenotypes.
//'
//' @param component_type String. The phenotype component to track. One of
//'   "genetic", "environmental", "vertical", or "total".
//'
//' @return An external pointer to a MetricSpec object for use with
//'   Simulation$metrics().
//'
//' @noRd
//'
// [[Rcpp::export]]
SEXP MetricSpec_pheno_comp_cor(const std::string& component_type) {
  amsim::ComponentType type = strComponentType(component_type);
  amsim::MetricSpec pheno_comp_cor = amsim::pheno_comp_cor(type);
  auto metric = std::make_unique<amsim::MetricSpec>(pheno_comp_cor);
  Rcpp::XPtr<amsim::MetricSpec> ptr(metric.release(), true);
  return ptr;
}

//' Phenotype component cross-correlation metric
//'
//' @description
//' Creates a metric specification for tracking cross-mate correlations of a
//' specific phenotype component across generations. This measures the
//' correlation between mated pairs for the specified component.
//'
//' @param component_type String. The phenotype component to track. One of
//'   "genetic", "environmental", "vertical", or "total".
//'
//' @return An external pointer to a MetricSpec object for use with
//'   Simulation$metrics().
//'
//' @noRd
//'
// [[Rcpp::export]]
SEXP MetricSpec_pheno_comp_xcor(const std::string& component_type) {
  amsim::ComponentType type = strComponentType(component_type);
  amsim::MetricSpec pheno_comp_xcor = amsim::pheno_comp_xcor(type);
  auto metric = std::make_unique<amsim::MetricSpec>(pheno_comp_xcor);
  Rcpp::XPtr<amsim::MetricSpec> ptr(metric.release(), true);
  return ptr;
}

//' Latent phenotype heritability metric
//'
//' @description
//' Creates a metric specification for tracking narrow-sense heritability
//' (h²) of latent phenotype values across generations. Latent values
//' represent the underlying trait values before observation noise.
//'
//' @return An external pointer to a MetricSpec object for use with
//'   Simulation$metrics().
//'
//' @noRd
//'
// [[Rcpp::export]]
SEXP MetricSpec_pheno_latent_h2() {
  amsim::MetricSpec pheno_latent_h2 = amsim::pheno_latent_h2();
  auto metric = std::make_unique<amsim::MetricSpec>(pheno_latent_h2);
  Rcpp::XPtr<amsim::MetricSpec> ptr(metric.release(), true);
  return ptr;
}

//' Latent phenotype component mean metric
//'
//' @description
//' Creates a metric specification for tracking the mean of a specific
//' latent phenotype component across generations.
//'
//' @param component_type String. The phenotype component to track. One of
//'   "genetic", "environmental", "vertical", or "total".
//'
//' @return An external pointer to a MetricSpec object for use with
//'   Simulation$metrics().
//'
//' @noRd
//'
// [[Rcpp::export]]
SEXP MetricSpec_pheno_latent_comp_mean(const std::string& component_type) {
  amsim::ComponentType type = strComponentType(component_type);
  amsim::MetricSpec pheno_comp_mean = amsim::pheno_latent_comp_mean(type);
  auto metric = std::make_unique<amsim::MetricSpec>(pheno_comp_mean);
  Rcpp::XPtr<amsim::MetricSpec> ptr(metric.release(), true);
  return ptr;
}

//' Latent phenotype component variance metric
//'
//' @description
//' Creates a metric specification for tracking the variance of a specific
//' latent phenotype component across generations.
//'
//' @param component_type String. The phenotype component to track. One of
//'   "genetic", "environmental", "vertical", or "total".
//'
//' @return An external pointer to a MetricSpec object for use with
//'   Simulation$metrics().
//'
//' @noRd
//'
// [[Rcpp::export]]
SEXP MetricSpec_pheno_latent_comp_var(const std::string& component_type) {
  amsim::ComponentType type = strComponentType(component_type);
  amsim::MetricSpec pheno_comp_var = amsim::pheno_latent_comp_var(type);
  auto metric = std::make_unique<amsim::MetricSpec>(pheno_comp_var);
  Rcpp::XPtr<amsim::MetricSpec> ptr(metric.release(), true);
  return ptr;
}

//' Latent phenotype component correlation metric
//'
//' @description
//' Creates a metric specification for tracking within-phenotype correlations
//' of a specific latent component across generations.
//'
//' @param component_type String. The phenotype component to track. One of
//'   "genetic", "environmental", "vertical", or "total".
//'
//' @return An external pointer to a MetricSpec object for use with
//'   Simulation$metrics().
//'
//' @noRd
//'
// [[Rcpp::export]]
SEXP MetricSpec_pheno_latent_comp_cor(const std::string& component_type) {
  amsim::ComponentType type = strComponentType(component_type);
  amsim::MetricSpec pheno_comp_cor = amsim::pheno_latent_comp_cor(type);
  auto metric = std::make_unique<amsim::MetricSpec>(pheno_comp_cor);
  Rcpp::XPtr<amsim::MetricSpec> ptr(metric.release(), true);
  return ptr;
}

//' Latent phenotype component cross-correlation metric
//'
//' @description
//' Creates a metric specification for tracking cross-mate correlations of a
//' specific latent phenotype component across generations.
//'
//' @param component_type String. The phenotype component to track. One of
//'   "genetic", "environmental", "vertical", or "total".
//'
//' @return An external pointer to a MetricSpec object for use with
//'   Simulation$metrics().
//'
//' @noRd
//'
// [[Rcpp::export]]
SEXP MetricSpec_pheno_latent_comp_xcor(const std::string& component_type) {
  amsim::ComponentType type = strComponentType(component_type);
  amsim::MetricSpec pheno_comp_xcor = amsim::pheno_latent_comp_xcor(type);
  auto metric = std::make_unique<amsim::MetricSpec>(pheno_comp_xcor);
  Rcpp::XPtr<amsim::MetricSpec> ptr(metric.release(), true);
  return ptr;
}
