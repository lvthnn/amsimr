Phenotype <- R6::R6Class(
  "Phenotype",

  public = list(
    initialize = function(name, heritability, n_loci = NULL, causal_snps = NULL,
                          causal_effects = NULL) {
      checkmate::assert_string(name, null.ok = FALSE)
      checkmate::assert_numeric(heritability, len = 1, lower = 0, upper = 1,
                                null.ok = FALSE)
      if (missing(n_loci)) {
        checkmate::assert_character(causal_snps, pattern = "rs[0-9]+$")
        checkmate::assert_numeric(causal_effects, len = length(causal_snps))
        private$.causal_snps    <- causal_snps
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

    # Getters
    name           = function() return(private$.name),
    data           = function() return(private$.data),
    n_loci         = function() return(private$.n_loci),
    causal_snps    = function() return(private$.causal_snps),
    causal_effects = function() return(private$.causal_effects),
    heritability   = function() return(private$.heritability),

    # Setters
    set_n_loci = function(n_loci) private$.n_loci <- n_loci,
    set_causal_snps = function(causal_snps) {
      private$.causal_snps <- causal_snps
    },
    set_causal_effects = function(causal_effects) {
      private$.causal_effects <- causal_effects
    },

    # Other
    compute_data = function(snp_matrix) {
      checkmate::assert_matrix(snp_matrix, mode = "integer",
                               any.missing = FALSE, all.missing = FALSE,
                               ncols = self$n_loci())
      private$.data <- snp_matrix %*% self$causal_effects() +
        rnorm(nrow(snp_matrix), 0, 1 - self$heritability())
    }
  ),

  private = list(
    .name           = NULL,
    .data           = NULL,
    .n_loci         = NULL,
    .causal_snps    = NULL,
    .causal_effects = NULL,
    .heritability   = NULL
  )
)