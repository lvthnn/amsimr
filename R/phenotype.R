Phenotype <- R6::R6Class(
  "Phenotype",
  public = list(

    initialize = function(name, n_loci, heritability) {
      checkmate::assert_string(name, null.ok = FALSE)
      checkmate::assert_count(n_loci, null.ok = FALSE)
      checkmate::assert_numeric(heritability, len = 1, lower = 0,
                                upper = 1, null.ok = FALSE)

      private$name <- name
      private$n_loci <- n_loci
      private$heritability <- heritability
    },

    get_name = function() return(private$name),
    get_n_loci = function() return(private$n_loci),
    get_causal_snps = function() return(private$causal_effects)
  ),
  private = list(
    name = NULL,
    n_loci = NULL,
    causal_snps = NULL,
    causal_effects = NULL,
    heritability = NULL
  )
)