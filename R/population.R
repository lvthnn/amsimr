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
      private$.n_pairs      <- n_population / 2
      private$.n_loci       <- n_loci
      private$.snp_mafs     <- snp_mafs
      private$.snp_dist     <- snp_dist
      private$.snp_params   <- snp_params
      private$.phenotypes   <- if (!is.null(phenotypes)) phenotypes else list()
      private$.mating_model <- if (!is.null(mating_model)) mating_model
    },

    # Getter methods
    genotype_data  = function() return(private$.genotype_data),
    phenotype_data = function() return(private$.phenotype_data),
    n_population   = function() return(private$.n_population),
    n_pairs        = function() return(private$.n_population),
    n_loci         = function() return(private$.n_loci),
    mating_model   = function() return(private$.mating_model),
    snp_mafs       = function() return(private$.snp_mafs),
    phenotypes     = function() return(private$.phenotypes),

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
      private$.setup_phenotypes()
      private$.initialise_data()
      private$.compiled <- TRUE
    }
  ),

  private = list(
    # Fields
    .genotype_data  = NULL,
    .phenotype_data = NULL,
    .n_population   = NULL,
    .n_pairs        = NULL,
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

    .initial_ancestral_snps = function() {
      mat_dim    <- self$n_pairs() * self$n_loci()
      snp_probs  <- rep(self$snp_mafs(), each = self$n_pairs())
      snp_counts <- rbinom(mat_dim, size = 2, prob = snp_probs)
      return(matrix(snp_counts, self$n_pairs() * self$n_loci(),
                    nrow = self$n_pairs(), ncol = self$n_loci()))
    },

    .generate_haplotype = function(snp_parental) {
      prob_parental      <- snp_parental / 2
      prob_haplotype     <- rep(prob_parental, each = 2)
      mat_dim            <- self$n_pairs() * self$n_loci()
      haplotype_parental <- rbinom(mat_dim, size = 1, prob = prob_haplotype)
      return(matrix(haplotype_parental, self$n_pairs() * self$n_loci(),
                    nrow = self$n_population(), ncol = self$n_loci()))
    },

    .generate_snp_matrix = function(snps_paternal, snps_maternal) {
      haplotype_paternal <- private$.generate_haplotype(snps_paternal)
      haplotype_maternal <- private$.generate_haplotype(snps_maternal)
      snpmat <- haplotype_paternal + haplotype_maternal
      return(snpmat)
    },

    .initialise_data = function() {
      sex <- rep(0:1, times = self$n_pairs())
      sibling_id <- rep(1:self$n_pairs(), each = 2)
      snps_paternal <- private$.initial_ancestral_snps()
      snps_maternal <- private$.initial_ancestral_snps()
      snps_init <- private$.generate_snp_matrix(snps_paternal, snps_maternal)
      private$.genotype_data <- snps_init
    },

    .check_compiled = function() {
      if (private$.compiled) cli::cli_abort("Cannot modify object after
        compilation. Create a new instance to make changes.")
    }
  )
)