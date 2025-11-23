#' Simulation Results
#' @description
#' An R6 class providing an interface for accessing and analyzing summarised
#' simulation output data. This class automatically aggregates replicate data
#' into summary statistics and loads metrics into R-native data structures for
#' further analysis.
#'
#' @details
#' The SimulationResults class provides functionality for:
#' * Loading simulation output from a specified directory
#' * Aggregating replicate data into summary statistics
#' * Accessing metric data as data frames
#'
#' Upon initialization, the class automatically performs summarisation across
#' all replicates and loads all available metrics. Each metric is stored as a
#' data frame containing the following summary statistics per generation:
#' * `gen`: Generation number
#' * `name`: Metric component name
#' * `mean`: Arithmetic mean across replicates
#' * `median`: Median value across replicates
#' * `stddev`: Standard deviation across replicates
#' * `stderr`: Standard error of the mean
#' * `lower_ci95`: Lower bound of 95% confidence interval
#' * `upper_ci95`: Upper bound of 95% confidence interval
#' * `quant_025`: 2.5th percentile
#' * `quant_975`: 97.5th percentile
#'
#' @examples
#' \dontrun{
#' # Load results from a completed simulation
#' results <- SimulationResults$new("output/my_simulation")
#'
#' # Get available metric names
#' results$metric_names
#'
#' # Retrieve heritability estimates
#' h2_data <- results$get("pheno_h2")
#'
#' # Plot mean heritability over generations
#' library(ggplot2)
#' ggplot(h2_data, aes(x = gen, y = mean, color = name)) +
#'   geom_line() +
#'   geom_ribbon(aes(ymin = lower_ci95, ymax = upper_ci95), alpha = 0.2)
#' }
#'
#' @export
SimulationResults <- R6::R6Class(
  classname = "SimulationResults",

  public = list(
    #' @description
    #' Create a new SimulationResults object
    #'
    #' @details
    #' Initializes the results object by reading simulation output from the
    #' specified directory, performing summarisation across all replicates,
    #' and loading all available metrics into memory as data frames.
    #'
    #' @param output_dir String. Path to the simulation output directory
    #'   containing replicate data files.
    #'
    #' @return A new `SimulationResults` instance with all metrics loaded.
    initialize = function(output_dir) {
      checkmate::assert_string(output_dir)
      private$.results_obj <- .SimulationResults_new(output_dir)
      .SimulationResults_summarise(private$.results_obj)

      metric_names <- .SimulationResults_metric_names(private$.results_obj)
      for (metric_name in metric_names) {
        private$.results_dat[[metric_name]] <- .SimulationResults_load(
          private$.results_obj,
          metric_name
        )
      }
    },

    #' @description
    #' Retrieve summarised data for a specific metric
    #'
    #' @param metric_name String. Name of the metric to retrieve. Must be one
    #'   of the available metrics (see `metric_names` field).
    #'
    #' @return A data frame containing summarised metric data with columns for
    #'   generation, name, and summary statistics (mean, median, stddev,
    #'   stderr, confidence intervals, and quantiles).
    get = function(metric_name) {
      checkmate::assert_choice(
        metric_name,
        choices = names(private$.results_dat)
      )
      private$.results_dat[[metric_name]]
    }
  ),

  private = list(
    .results_obj = NULL,
    .results_dat = list()
  ),

  active = list(
    #' @field metric_names Character vector of available metric names
    metric_names = function() {
      .SimulationResults_metric_names(private$.results_obj)
    }
  )
)
