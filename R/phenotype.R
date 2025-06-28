#' Phenotype Class for Genetic Simulation
#'
#' @description
#' The `Phenotype` R6 class defines and manages phenotypic traits in genetic
#' simulations. It allows specification of trait heritability and causal
#' genetic architecture through either a defined number of random causal loci
#' or explicit causal SNPs with their associated effect sizes.
#'
#' @details
#' The Phenotype class supports two modes of operation: (1) specifying the
#' number of causal loci for random selection, or (2) providing explicit
#' causal SNPs with their effect sizes. The class computes phenotypic values
#' by scaling SNP data and calculating the genetic component based on the
#' specified causal architecture. The heritability parameter controls the
#' proportion of phenotypic variance explained by genetic factors.
#'
#' @section Methods:
#' - `initialize()`: Create a new phenotype with specified parameters.
#' - `compute_values()`: Calculate phenotypic values from SNP matrix data.
#'
#' @seealso
#' Other related classes and functions: `Population`, `MatingModel`
#'
#' @export
Phenotype <- R6::R6Class(
  "Phenotype",

  public = list(
    #' @description Initialise a new Phenotype object
    #' @param name String. Name of the phenotype.
    #' @param heritability Numeric. The narrow-sense heritability of the
    #'   phenotype.
    #' @param n_loci Integer. Number of causal loci affecting the phenotype.
    #' @param causal_snps Integer vector of length equal to `n_loci` specifying
    #'   the indices of causal SNPs.
    #' @param causal_effects Numeric vector of length equal to `n_loci`
    #'   specifying the effect sizes of causal SNPs.
    initialize = function(name, heritability, n_loci = NULL, causal_snps = NULL,
                          causal_effects = NULL) {
      checkmate::assert_string(name, null.ok = FALSE)
      checkmate::assert_numeric(heritability, len = 1, lower = 0, upper = 1,
                                null.ok = FALSE)
      if (missing(n_loci)) {
        checkmate::assert_integer(causal_snps, any.missing = FALSE,
                                  all.missing = FALSE)
        checkmate::assert_numeric(causal_effects, len = length(causal_snps))
        private$.causal_snps    <- gsub("rs", "", causal_snps)
        private$.causal_effects <- causal_effects
        private$.n_loci         <- length(causal_snps)
      } else if (missing(causal_snps) && missing(causal_effects)) {
        checkmate::assert_count(n_loci)
        private$.n_loci <- n_loci
      } else {
        cli::cli_abort(c(
          "x" = "Invalid parameter combination: provide either `n_loci` or both
                 `causal_snps` and `causal_effects`."
        ))
      }
      private$.name         <- name
      private$.heritability <- heritability
    },

    #' @description Compute the data vector of a Phenotype.
    #' @param snp_matrix Matrix of SNPs from a Population instance.
    compute_values = function(snp_matrix) {
      scale_fbm   <- bigstatsr::big_scale()
      scale_data  <- scale_fbm(snp_matrix)[self$causal_snps, ]
      self$values <- bigstatsr::big_prodVec(
        X       = snp_matrix,
        y.col   = self$causal_effects,
        ind.col = self$causal_snps,
        center  = scale_data$center,
        scale   = scale_data$scale
      )
      # + rnorm(nrow(snp_matrix), sd = sqrt(1 - self$heritability))
    }
  ),

  private = list(
    .name           = NULL,
    .values         = NULL,
    .n_loci         = NULL,
    .causal_snps    = NULL,
    .causal_effects = NULL,
    .heritability   = NULL
  ),

  active = list(
    #' @field name Get the phenotype name.
    name = function() return(private$.name),

    #' @field values Get the values for the phenotype.
    values = function(value) {
      if (missing(value)) {
        return(private$.values)
      } else {
        private$.values <- value
      }
    },

    #' @field n_loci Get the number of causal loci.
    n_loci = function(value) {
      if (missing(value)) {
        return(private$.n_loci)
      } else {
        private$.n_loci <- value
      }
    },

    #' @field causal_snps Get the indices of causal SNPs.
    causal_snps = function(value) {
      if (missing(value)) {
        return(private$.causal_snps)
      } else {
        private$.causal_snps <- value
      }
    },

    #' @field causal_effects Get the causal effects.
    causal_effects = function(value) {
      if (missing(value)) {
        return(private$.causal_effects)
      } else {
        private$.causal_effects <- value
      }
    },

    #' @field heritability Get the narrow-sense heritability.
    heritability = function() return(private$.heritability)
  )
)
