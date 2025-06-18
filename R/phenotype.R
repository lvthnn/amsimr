#' Phenotype Class
#'
#' Represents a phenotype with its associated genetic properties.
#'
#' @importFrom R6 R6Class
#' @importFrom checkmate assert_string assert_count assert_numeric assert_number
#'
#' @export
Phenotype <- R6Class(
  "Phenotype",
  public = list(
    #' @field name Name of phenotype
    name = NULL,

    #' @field loci Indices of causal loci
    loci = NULL,

    #' @field n_loci Number of causal loci
    n_loci = NULL,

    #' @field loci_effects Genic effects of causal loci
    loci_effects = NULL,

    #' @field heritability Narrow-sense heritability of phenotype
    heritability = NULL,

    #' Initialise a Phenotype object
    #'
    #' @param name Name of the phenotype
    #' @param n_loci Number of causal loci
    #' @param heritability Narrow-sense heritability of phenotype
    #' @param loci_effects Vector of genic effects of causal loci
    initialize = function(name, heritability, loci = NULL, n_loci = NULL,
                          loci_effects = NULL) {
      assert_string(name)
      assert_count()
      assert_
      assert_numeric(heritability, lower = 0, upper = 1, len = 1)

      if (is.null(loci_effects)) {
        loci_effects <- rep(heritability / n_loci, n_loci)
      } else {
        assert_numeric(loci_effects, len = n_loci, any.missing = FALSE)
      }

      if (!is.null(loci) && !is.null(n_loci))

      self$name <- name
      self$n_loci <- n_loci
      self$loci_effects <- loci_effects
      self$heritability <- heritability
    }
  )
)
