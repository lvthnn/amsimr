Population <- R6::R6Class(
  "Population",

  public = list(
    initialize = function(n_population, n_loci, snp_dist = "beta",
                          snp_params = c(1.5, 2.5), phenotypes = NULL,
                          mating_model = NULL) {
      checkmate::assert_count(n_population)
      checkmate::assert_count(n_loci)
      checkmate::assert_choice(snp_dist, choices = c("beta", "unif"))
      checkmate::assert_numeric(snp_params, len = 2)
      checkmate::assert_class(phenotypes, "Phenotype", null.ok = TRUE)
      checkmate::assert_class(mating_model, "MatingModel", null.ok = TRUE)

      rdist  <- get(paste0("r", snp_dist), mode = "function")
      param1 <- snp_params[[1]]
      param2 <- snp_params[[2]]

      snp_mafs <- rdist(n_loci, param1, param2)
      names(snp_mafs) <- paste0("rs", 1:n_loci)

      checkmate::assert_numeric(snp_mafs, lower = 0, upper = 1, finite = TRUE,
                                any.missing = FALSE, all.missing = FALSE)

      private$.n_population <- n_population
      private$.n_loci       <- n_loci
      private$.snp_mafs     <- snp_mafs
      private$.snp_dist     <- snp_dist
      private$.snp_params   <- snp_params
      private$.phenotypes   <- if (!is.null(phenotypes)) phenotypes else list()
      private$.mating_model <- if (!is.null(mating_model)) mating_model
    },

    # Getter methods
    data         = function() return(private$.data),
    n_population = function() return(private$.n_population),
    n_loci       = function() return(private$.n_loci),
    mating_model = function() return(private$.mating_model),
    snp_mafs     = function() return(private$.snp_mafs),
    phenotypes   = function() return(private$.phenotypes),

    # Setter methods
    set_mating_model = function(mating_model) {
      private$.check_compiled()
      checkmate::assert_class(mating_model, "MatingModel")
      private$.mating_model <- mating_model
    },

    # Other
    add_phenotype = function(phenotype) {
      private$.check_compiled()
      checkmate::assert_class(phenotype, "Phenotype")
      private$.phenotypes[[phenotype$name()]] <- phenotype
    },
    remove_phenotype = function(phenotype_name) {
      private$.check_compiled()
      phenotype_names <- names(self$phenotypes())
      checkmate::assert_choice(phenotype_name, choices = phenotype_names)
      private$.phenotypes[[phenotype_name]] <- NULL
    },
    compile = function() {
      # Compiles the Population object, generating data and locking further
      # modifications to configuration
      private$.setup_phenotypes()
      private$.initialise_data()
    }
  ),

  private = list(
    # Fields
    .genotype_data  = NULL,
    .phenotype_data = NULL,
    .n_population   = NULL,
    .n_loci         = NULL,
    .data           = NULL,
    .snp_mafs       = NULL,
    .snp_dist       = NULL,
    .snp_params     = NULL,
    .phenotypes     = NULL,
    .mating_model   = NULL,
    .compiled       = FALSE,

    # Utility methods
    .setup_phenotypes = function() {
      for (phenotype in self$phenotypes()) {
        if (is.null(phenotype$causal_snps())) {
          if (phenotype$n_loci() > self$n_loci()) {
            cli::cli_abort(c("x" = "Number of causal SNPs specified for
                                    phenotype `{phenotype$name()}` exceeds
                                    number of SNPs modelled in the
                                    population."))
          }
          causal_idx <- sample(1:self$n_loci(), size = phenotype$n_loci())
          causal_eff <- rep(phenotype$heritability() / phenotype$n_loci(),
                            phenotype$n_loci())
          phenotype$set_causal_snps(paste0("rs", causal_idx))
          phenotype$set_causal_effects(causal_eff)
        }
      }
    },

    .initialise_data = function() {
      # - Generate ancestral haplotype matrix based on MAFs (Rcpp),
      #   and use to create genetic matrix. Store as an arma matrix
      #   in the .genotype_data field.

      # - Compute trait vectors using genetic matrix and store in
      #   .phenotype_data field.
    },

    .check_compiled = function() {
      if (private$.compiled) cli::cli_abort("Cannot modify object after
        compilation. Create a new instance to make changes.")
    }
  )
)