#' Simulation Configuration
#' @description
#' An R6 class for configuring and running genetic simulations with assortative
#' mating. This class provides a builder pattern interface for setting update
#' simulation parameters including genome architecture, phenotypic traits,
#' mating strategies, as well as the metrics and statistics used to produce
#' output data for analysis.
#'
#' @details
#' The Simulation class allows configuration of:
#' * Basic simulation parameters (generations, population size, output)
#' * Genome architecture (loci, MAF, recombination, mutation rates)
#' * Phenome specification (traits, heritability, correlations)
#' * Mating strategies (random or assortative)
#' * Metrics collection
#'
#' All configuration methods return `self` for method chaining.
#'
#' @examples
#' \dontrun{
#' sim <- Simulation$new()
#' sim$simulation(
#'   n_generations = 100,
#'   n_individuals = 1000,
#'   output_dir = "output",
#'   random_seed = 42
#' )$genome(
#'   n_loci = 1000,
#'   locus_maf = 0.3,
#'   locus_recombination = 0.5,
#'   locus_mutation = 1e-8
#' )$phenome(
#'   n_phenotypes = 2,
#'   names = c("height", "weight"),
#'   n_causal_loci = c(100, 100),
#'   h2_genetic = c(0.5, 0.6),
#'   h2_environmental = c(0.3, 0.2),
#'   h2_vertical = c(0.2, 0.2),
#'   genetic_cor = diag(2),
#'   environmental_cor = diag(2)
#' )$random_mating()$metrics(
#'   pheno_h2(),
#'   pheno_comp_xcor("total")
#' )
#' }
#'
#' @export
Simulation <- R6::R6Class(
  classname = "Simulation",

  public = list(
    #' @description
    #' Create a new Simulation object
    initialize = function() {
      private$.config <- .SimulationConfig_new()
    },

    #' @description
    #' Configure basic simulation parameters
    #'
    #' @param n_generations Integer. Number of generations to simulate.
    #' @param n_individuals Integer. Population size.
    #' @param output_dir String. Directory path for output files.
    #' @param random_seed Integer or NULL. Random seed for reproducibility.
    #' @return Self, invisibly, for method chaining.
    simulation = function(
      n_generations,
      n_individuals,
      output_dir,
      random_seed = NULL
    ) {
      checkmate::assert_count(n_generations, positive = TRUE)
      checkmate::assert_count(n_individuals, positive = TRUE)
      checkmate::assert_string(output_dir)
      checkmate::assert_count(random_seed, positive = TRUE, null.ok = TRUE)

      .SimulationConfig_simulation(
        ptr = private$.config,
        n_generations = n_generations,
        n_individuals = n_individuals,
        output_dir = output_dir,
        random_seed = random_seed
      )

      invisible(self)
    },

    #' @description
    #' Configure genome architecture
    #'
    #' @param n_loci Integer. Number of genetic loci.
    #' @param locus_maf Numeric vector. Minor allele frequencies for each locus.
    #'   Can be a single value (recycled) or a vector of length `n_loci`.
    #' @param locus_recombination Numeric vector. Recombination probabilities.
    #'   Can be a single value (recycled) or a vector of length `n_loci`.
    #' @param locus_mutation Numeric vector. Mutation probabilities.
    #'   Can be a single value (recycled) or a vector of length `n_loci`.
    #' @return Self, invisibly, for method chaining.
    genome = function(
      n_loci,
      locus_maf,
      locus_recombination,
      locus_mutation
    ) {
      checkmate::assert_count(n_loci, positive = TRUE)

      if (length(locus_maf) == 1) {
        locus_maf <- rep(locus_maf, n_loci)
      }
      if (length(locus_recombination) == 1) {
        locus_recombination <- rep(locus_recombination, n_loci)
      }
      if (length(locus_mutation) == 1) {
        locus_mutation <- rep(locus_mutation, n_loci)
      }

      assert_probs(locus_maf, len = n_loci)
      assert_probs(locus_recombination, len = n_loci)
      assert_probs(locus_mutation, len = n_loci)

      .SimulationConfig_genome(
        ptr = private$.config,
        n_loci = n_loci,
        locus_maf = locus_maf,
        locus_recombination = locus_recombination,
        locus_mutation = locus_mutation
      )

      invisible(self)
    },

    #' @description
    #' Configure phenome (trait) architecture
    #'
    #' @param n_phenotypes Integer. Number of phenotypic traits.
    #' @param names Character vector. Names for each phenotype.
    #' @param n_causal_loci Integer vector. Number of causal loci per phenotype.
    #' @param h2_genetic Numeric vector. Narrow-sense heritabilities.
    #' @param h2_environmental Numeric vector. Environmental variance
    #'   components.
    #' @param h2_vertical Numeric vector. Vertical transmission components.
    #' @param genetic_cor Numeric matrix. Genetic correlation matrix.
    #' @param environmental_cor Numeric matrix. Environmental correlation
    #'   matrix.
    #' @return Self, invisibly, for method chaining.
    phenome = function(
      n_phenotypes,
      names,
      n_causal_loci,
      h2_genetic,
      h2_environmental,
      h2_vertical,
      genetic_cor,
      environmental_cor
    ) {
      checkmate::assert_count(n_phenotypes)
      checkmate::assert_character(names, len = n_phenotypes, unique = TRUE)
      checkmate::assert_integerish(
        n_causal_loci,
        len = n_phenotypes,
        lower = 1
      )
      assert_probs(h2_genetic)
      assert_probs(h2_environmental)
      assert_probs(h2_vertical)
      checkmate::assert_matrix(
        genetic_cor,
        nrows = n_phenotypes,
        ncols = n_phenotypes,
        any.missing = FALSE,
        all.missing = FALSE
      )
      checkmate::assert_matrix(
        environmental_cor,
        nrows = n_phenotypes,
        ncols = n_phenotypes,
        any.missing = FALSE,
        all.missing = FALSE
      )

      .SimulationConfig_phenome(
        private$.config,
        n_phenotypes = n_phenotypes,
        names = names,
        n_causal_loci = n_causal_loci,
        h2_genetic = h2_genetic,
        h2_environmental = h2_environmental,
        h2_vertical = h2_vertical,
        genetic_cor = genetic_cor,
        environmental_cor = environmental_cor
      )

      invisible(self)
    },

    #' @description
    #' Configure random mating strategy
    #'
    #' @return Self, invisibly, for method chaining.
    random_mating = function() {
      .SimulationConfig_random_mating(private$.config)

      invisible(self)
    },

    #' @description
    #' Configure assortative mating strategy
    #'
    #' @param mate_cor Numeric matrix. Mating correlation matrix.
    #' @param tol_inf Numeric or NULL. Tolerance threshold for annealing.
    #' @param n_iterations Integer or NULL. Maximum annealing iterations.
    #' @param temp_init Numeric or NULL. Initial annealing temperature.
    #' @param temp_decay Numeric or NULL. Temperature decay rate.
    #' @return Self, invisibly, for method chaining.
    assortative_mating = function(
      mate_cor = mate_cor,
      tol_inf = NULL,
      n_iterations = NULL,
      temp_init = NULL,
      temp_decay = NULL
    ) {
      checkmate::assert_matrix(
        mate_cor,
        nrows = self$n_phenotypes,
        ncols = self$n_phenotypes
      )
      checkmate::assert_numeric(tol_inf, lower = 0, null.ok = TRUE)
      checkmate::assert_count(n_iterations, null.ok = TRUE)
      checkmate::assert_numeric(temp_init, lower = 0, null.ok = TRUE)
      checkmate::assert_numeric(
        temp_decay,
        lower = 0,
        upper = 1,
        null.ok = TRUE
      )

      .SimulationConfig_assortative_mating(
        private$.config,
        mate_cor = mate_cor,
        tol_inf = tol_inf,
        n_iterations = n_iterations,
        temp_init = temp_init,
        temp_decay = temp_decay
      )

      invisible(self)
    },

    #' @description
    #' Configure simulation metrics
    #'
    #' @param metrics List. List of metric specification objects.
    #' @return Self, invisibly, for method chaining.
    metrics = function(metrics = list()) {
      checkmate::assert_list(metrics)
      .SimulationConfig_metrics(private$.config, metrics = metrics)

      invisible(self)
    },

    #' @description
    #' Run the configured genetic simulation
    #'
    #' @details
    #' Executes the simulation pipeline with the configured parameters. The
    #' simulation runs for the specified number of replicates, optionally using
    #' multiple threads for parallel execution. Results are written to the
    #' configured output directory.
    #'
    #' The simulation proceeds through the following steps for each generation:
    #' 1. Evaluate phenotypes for all individuals
    #' 2. Perform mating according to the configured strategy
    #' 3. Generate offspring through reproduction (recombination and mutation)
    #' 4. Collect configured metrics
    #' 5. Write output data
    #'
    #' @param n_replicates Integer. Number of independent simulation replicates.
    #' @param n_threads Integer. Number of parallel threads for execution.
    #' @param summarise Logical. If TRUE, aggregate results across replicates
    #'   (default: TRUE).
    #' @param log_file Logical. If TRUE, write log messages to file
    #'   (default: TRUE).
    #' @param log_level String. Logging verbosity: "trace", "debug", "info",
    #'   "warn", "error", or "critical" (default: "info").
    #'
    #' @return Self (invisibly). Results are written to the output directory. If
    #'   the summarise argument is TRUE, the data are aggregated into summary
    #'   statistic tables and are made accesssible under the results field.
    #'
    #' @examples
    #' \dontrun{
    #' sim <- Simulation$new()
    #' # ... configure simulation ...
    #' sim$run(
    #'   n_replicates = 10,
    #'   n_threads = 4,
    #'   summarise = TRUE,
    #'   log_level = "info"
    #' )
    #' }
    run = function(
      n_replicates,
      n_threads,
      summarise = TRUE,
      log_file = TRUE,
      log_level = "info"
    ) {
      checkmate::assert_count(n_replicates)
      checkmate::assert_count(n_threads)
      checkmate::assert_logical(summarise)
      checkmate::assert_logical(log_file)
      checkmate::assert_string(log_level)

      .run_simulations(
        private$.config,
        n_replicates = n_replicates,
        n_threads = n_threads,
        summarise = summarise,
        log_file = log_file,
        log_level = log_level
      )

      if (summarise) {
        private$.results <- SimulationResults$new(self$output_dir)
      }

      invisible(self)
    }
  ),

  private = list(
    .config = NULL,
    .results = NULL
  ),

  active = list(
    #' @field n_generations Number of generations
    n_generations = function() {
      .SimulationConfig_n_generations(private$.config)
    },
    #' @field n_individuals Population size
    n_individuals = function() {
      .SimulationConfig_n_individuals(private$.config)
    },
    #' @field output_dir Output directory path
    output_dir = function() {
      .SimulationConfig_output_dir(private$.config)
    },
    #' @field random_seed Random seed value
    random_seed = function() {
      .SimulationConfig_random_seed(private$.config)
    },
    #' @field n_loci Number of genetic loci
    n_loci = function() {
      .SimulationConfig_n_loci(private$.config)
    },
    #' @field locus_maf Vector of locus minor allele frequencies
    locus_maf = function() {
      .SimulationConfig_locus_maf(private$.config)
    },
    #' @field locus_recombination Vector of locus recombination probabilities
    locus_recombination = function() {
      .SimulationConfig_locus_recombination(private$.config)
    },
    #' @field locus_mutation Vector of locus mutation probabilities
    locus_mutation = function() {
      .SimulationConfig_locus_mutation(private$.config)
    },
    #' @field n_phenotypes Number of phenotypes
    n_phenotypes = function() {
      .SimulationConfig_n_phenotypes(private$.config)
    },
    #' @field phenotype_names Vector of phenotype names
    phenotype_names = function() {
      .SimulationConfig_phenotype_names(private$.config)
    },
    #' @field n_causal_loci Vector of causal loci counts per phenotype
    n_causal_loci = function() {
      .SimulationConfig_n_causal_loci(private$.config)
    },
    #' @field h2_genetic Vector of narrow-sense heritabilities
    h2_genetic = function() {
      .SimulationConfig_h2_genetic(private$.config)
    },
    #' @field h2_environmental Vector of environmental variance components
    h2_environmental = function() {
      .SimulationConfig_h2_environmental(private$.config)
    },
    #' @field h2_vertical Vector of vertical transmission components
    h2_vertical = function() {
      .SimulationConfig_h2_vertical(private$.config)
    },
    #' @field genetic_cor Genetic correlation matrix
    genetic_cor = function() {
      .SimulationConfig_genetic_cor(private$.config)
    },
    #' @field environmental_cor Environmental correlation matrix
    environmental_cor = function() {
      .SimulationConfig_environmental_cor(private$.config)
    },
    #' @field mate_cor Mating correlation matrix
    mate_cor = function() {
      .SimulationConfig_mate_cor(private$.config)
    },
    #' @field tol_inf Annealing tolerance threshold
    tol_inf = function() {
      .SimulationConfig_tol_inf(private$.config)
    },
    #' @field n_iterations Maximum number of annealing iterations
    n_iterations = function() {
      .SimulationConfig_n_iterations(private$.config)
    },
    #' @field temp_init Initial annealing temperature
    temp_init = function() {
      .SimulationConfig_temp_init(private$.config)
    },
    #' @field temp_decay Annealing temperature decay rate
    temp_decay = function() {
      .SimulationConfig_temp_decay(private$.config)
    },
    #' @field results Summary statistic of metrics from simulation
    results = function() {
      private$.results
    }
  )
)

#' @include utils.R
