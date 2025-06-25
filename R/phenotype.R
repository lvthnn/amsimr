Phenotype <- R6::R6Class(
  "Phenotype",

  public = list(
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
      scale_fbm  <- bigstatsr::big_scale()
      scale_data <- scale_fbm(snp_matrix)[self$causal_snps(), ]

      private$.data  <- bigstatsr::big_prodVec(X       = snp_matrix,
                                               y.col   = self$causal_effects(),
                                               ind.col = self$causal_snps(),
                                               center  = scale_data$center,
                                               scale   = scale_data$scale) # +
        # rnorm(nrow(snp_matrix), sd = sqrt(1 - self$heritability()))
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