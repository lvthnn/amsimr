#' Metric Specification
#'
#' @description
#' An R6 class for specifying metrics to collect during genetic simulations.
#' Metrics are statistics calculated and recorded at each generation, such as
#' heritabilities, phenotypic component means and variances, correlations, and
#' mate cross-correlations.
#'
#' @details
#' Metric objects are created using convenience functions:
#' * `pheno_h2()` - Phenotype heritabilities
#' * `pheno_comp_mean()` - Phenotype component means
#' * `pheno_comp_var()` - Phenotype component variances
#' * `pheno_comp_cor()` - Phenotype component correlations
#' * `pheno_comp_xcor()` - Between-mate phenotype component cross-correlations
#' * `pheno_latent_h2()` - Latent phenotype heritabilities
#' * `pheno_latent_comp_mean()` - Latent phenotype component means
#' * `pheno_latent_comp_var()` - Latent phenotype component variances
#' * `pheno_latent_comp_cor()` - Latent phenotype component correlations
#' * `pheno_latent_comp_xcor()` - Between-mate latent phenotype cross-correlations
#' * `custom_metric()` - User-defined metrics
#'
#' @examples
#' \dontrun{
#' sim <- Simulation$new()
#' sim$simulation(
#'   n_generations = 100,
#'   n_individuals = 1000,
#'   output_dir = "output",
#'   random_seed = 42
#' )$metrics(
#'   pheno_h2(),
#'   pheno_comp_mean("genetic"),
#'   pheno_comp_xcor("total")
#' )
#' }
#'
#' @export
Metric <- R6::R6Class(
  classname = "Metric",

  public = list(
    #' @description
    #' Create a new Metric object
    #'
    #' @param name String. Name of the metric.
    #' @param func Function. Metric specification function.
    #' @param params List. Parameters for the metric function.
    initialize = function(name = NULL, func = NULL, params = list()) {
      checkmate::assert_string(name)
      checkmate::assert_function(func)

      private$.name <- name
      private$.func <- func
      private$.params <- params
    },

    #' @description
    #' Build the metric specification
    #'
    #' @return Metric specification object for use in simulations.
    build = function() {
      do.call(private$.func, private$.params)
    }
  ),

  private = list(
    .name = NULL,
    .func = NULL,
    .params = list()
  )
)

#' Phenotype heritability metric
#'
#' @export
pheno_h2 <- function() {
  Metric$new(name = "pheno_h2", func = MetricSpec_pheno_h2)
}

#' Phenotype component means metric
#'
#' @param component_type Target component type.
#'
#' @export
pheno_comp_mean <- function(component_type) {
  Metric$new(
    name = "pheno_comp_mean",
    func = MetricSpec_pheno_comp_mean,
    params = list(component_type = component_type)
  )
}

#' Phenotype component variances metric
#'
#' @param component_type Target component type.
#'
#' @export
pheno_comp_var <- function(component_type) {
  Metric$new(
    name = "pheno_comp_var",
    func = MetricSpec_pheno_comp_var,
    params = list(component_type = component_type)
  )
}

#' Phenotype component correlations metric
#'
#' @param component_type Target component type.
#'
#' @export
pheno_comp_cor <- function(component_type) {
  Metric$new(
    name = "pheno_comp_cor",
    func = MetricSpec_pheno_comp_cor,
    params = list(component_type = component_type)
  )
}

#' Between-mate phenotype component cross-correlation metric
#'
#' @param component_type Target component type.
#'
#' @export
pheno_comp_xcor <- function(component_type) {
  Metric$new(
    name = "pheno_comp_xcor",
    func = MetricSpec_pheno_comp_xcor,
    params = list(component_type = component_type)
  )
}

#' Latent phenotype heritabilities metric
#'
#' @param component_type Target component type.
#'
#' @export
pheno_latent_h2 <- function() {
  Metric$new(name = "pheno_latent_h2", func = MetricSpec_pheno_latent_h2)
}

#' Latent phenotype component means metric
#'
#' @param component_type Target component type.
#'
#' @export
pheno_latent_comp_mean <- function(component_type) {
  Metric$new(
    name = "pheno_latent_comp_mean",
    func = MetricSpec_pheno_latent_comp_mean,
    params = list(component_type = component_type)
  )
}

#' Latent phenotype component variances metric
#'
#' @param component_type Target component type.
#'
#' @export
pheno_latent_comp_var <- function(component_type) {
  Metric$new(
    name = "pheno_latent_comp_var",
    func = MetricSpec_pheno_latent_comp_var,
    params = list(component_type = component_type)
  )
}

#' Latent phenotype component correlation metric
#'
#' @param component_type Target component type.
#'
#' @export
pheno_latent_comp_cor <- function(component_type) {
  Metric$new(
    name = "pheno_latent_comp_cor",
    func = MetricSpec_pheno_latent_comp_cor,
    params = list(component_type = component_type)
  )
}

#' Between-mate latent phenotype cross-correlation metric
#'
#' @param component_type Target component type.
#'
#' @export
pheno_latent_comp_xcor <- function(component_type) {
  Metric$new(
    name = "pheno_latent_comp_xcor",
    func = MetricSpec_pheno_latent_comp_xcor,
    params = list(component_type = component_type)
  )
}

#' Create a custom Metric
#'
#' @param name Name of the custom metric function
#' @param metric_func The R function defining the metric
#' @param n_rows Number of rows in output from single `metric_func` call
#' @param n_cols Number of columns in output from single `metric_func` call
#' @param labels Labels of output data produced by single `metric_func` call
#' @param require_lat Does this metric require a latent phenotype buffer?
#'
#' @export
custom_metric <- function(
  name,
  metric_func,
  n_rows = NULL,
  n_cols = NULL,
  labels = NULL,
  require_lat = FALSE
) {
  checkmate::assert_string(name)
  checkmate::assert_function(metric_func)
  checkmate::assert_count(n_rows, null.ok = TRUE)
  checkmate::assert_count(n_cols, null.ok = TRUE)
  checkmate::assert_character(labels)
  checkmate::assert_logical(require_lat)

  Metric$new(
    name = name,
    func = MetricSpec_custom_metric,
    params = list(
      name = name,
      metric_func = metric_func,
      n_rows = n_rows,
      n_cols = n_cols,
      require_lat = require_lat
    )
  )
}
