Population <- R6::R6Class(
  "Population",
  public = list(
    initialize = function(n_population, n_loci, phenotypes = NULL,
                          mating_model = NULL) {
      checkmate::assert_count(n_population)
      checkmate::assert_count(n_loci)
      checkmate::assert_class(phenotypes, "Phenotype", null.ok = TRUE)
      checkmate::assert_class(mating_model, "MatingModel", null.ok = TRUE)

      private$n_population <- n_population
      private$n_loci <- n_loci
      private$phenotypes <- if (!is.null(phenotypes)) phenotypes else list()

      # @TODO — Validate whether `mating_model` pairs matrix matches number of
      #         phenotypes and whether names are in accordance with names of
      #         phenotypes
      private$mating_model <- if (!is.null(mating_model)) mating_model else MatingModel$new()
    },

    get_n_population = function() return(private$n_population),
    get_n_loci = function() return(private$n_loci),

    get_mating_model = function() return(private$mating_model),
    set_mating_model = function(mating_model) { ... },

    get_phenotypes = function() return(private$phenotypes),
    add_phenotype = function(phenotype) {
      checkmate::assert_class(phenotype, "Phenotype")
      phenotypes <- append(phenotypes, phenotype)
    },
    remove_phenotype = function(phenotype) { ... },

  ),
  private = list(
    n_population = NULL,
    n_loci       = NULL,
    data         = NULL,
    phenotypes   = NULL,
    mating_model = NULL,

    initialise_data = function() {

    }
  )
)