Population <- R6::R6Class(
  "Population",

  public = list(
    initialize = function(n_population, n_loci, snp_dist = "beta",
                          snp_params = c(5, 5), phenotypes = NULL,
                          mating_model = NULL) {
      checkmate::assert_count(n_population)
      checkmate::assert_count(n_loci)
      checkmate::assert_choice(snp_dist, choices = c("beta", "unif"))
      checkmate::assert_numeric(snp_params, len = 2)
      checkmate::assert_class(phenotypes, "Phenotype", null.ok = TRUE)
      checkmate::assert_class(mating_model, "MatingModel", null.ok = TRUE)

      rdist           <- get(paste0("r", snp_dist), mode = "function")
      param1          <- snp_params[[1]]
      param2          <- snp_params[[2]]
      snp_mafs        <- rdist(n_loci, param1, param2)
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
    generation     = function() return(private$.generation),

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
    compile = function(verbose = TRUE) {
      cli::cli_h1("Compiling Population object:")
      cli::cli_h2("Initialising population data")
      private$.initialise_data()
      cli::cli_h2("Generating sibling pairs")
      private$.tick_generation()
      cli::cli_h2("Setting up phenotypes")
      private$.setup_phenotypes()
      private$.compiled <- TRUE
      cli::cli_alert_info("Done!")
    }
  ),

  private = list(
    # Fields
    .genotype_data  = NULL,
    .phenotype_data = NULL,
    .n_population   = NULL,
    .n_pairs        = NULL,
    .n_loci         = NULL,
    .snp_mafs       = NULL,
    .snp_dist       = NULL,
    .snp_params     = NULL,
    .summariser     = NULL,
    .phenotypes     = NULL,
    .mating_model   = NULL,
    .compiled       = FALSE,
    .generation     = 0,

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
          sd_eff     <- sqrt(phenotype$heritability() / phenotype$n_loci())
          causal_idx <- sample(1:self$n_loci(), size = phenotype$n_loci())
          causal_eff <- rnorm(phenotype$n_loci(), mean = 0, sd = sd_eff)
          phenotype$set_causal_snps(paste0("rs", causal_idx))
          phenotype$set_causal_effects(causal_eff)
        }
        causal_idx  <- sort(as.integer(gsub("rs", "", phenotype$causal_snps())))
        causal_snps <- self$genotype_data()[, causal_idx]
        phenotype$set_causal_snps(paste0("rs", causal_idx))
        phenotype$compute_data(causal_snps)
      }
    },

    .initialise_data = function() {
      private$.genotype_data <- .init_ancestral_snps(
        self$n_population(),
        self$n_loci(),
        self$snp_mafs()
      )
    },

    .tick_generation = function(ncores = parallel::detectCores() - 1) {
      private$.generation <- private$.generation + 1
      genotypes_new <- bigstatsr::FBM(
        nrow = self$n_population(),
        ncol = self$n_loci(),
        type = "unsigned short",
        backingfile = paste0("genotypes_gen", self$generation())
      )

      n_sex <- population$n_population() / 2
      chunk_size <- 10000

      # Create a list of chunks to fill the matrix
      start_rows <- seq(1, n_sex, by = chunk_size)
      end_rows   <- pmin(start_rows + chunk_size - 1, n_sex)
      chunk_data <- list(start = start_rows, end = end_rows)

      future::plan(future::multisession, workers = ncores)

      furrr::future_walk(seq_along(chunk_data$start), ~ {
        start_row      <- chunk_data$start[.x]
        end_row        <- chunk_data$end[.x]
        n_row          <- end_row - start_row + 1

        row_idx_male   <- start_row:end_row
        row_idx_female <- (n_sex + start_row):(n_sex + end_row)

        chunk_male     <- population$genotype_data()[row_idx_male, ]
        chunk_female   <- population$genotype_data()[row_idx_female, ]
        prob_male      <- rep(as.vector(chunk_male / 2), times = 2)
        prob_female    <- rep(as.vector(chunk_female / 2), times = 2)

        male_gametes <- matrix(
          rbinom(
            n_row * population$n_loci() * 2,
            size = 1,
            prob = prob_male
          ),
          nrow = n_row * 2,
          ncol = population$n_loci(),
        )

        female_gametes <- matrix(
          rbinom(
            n_row * population$n_loci() * 2,
            size = 1,
            prob = prob_female
          ),
          nrow = n_row * 2,
          ncol = population$n_loci()
        )

        male_gametes_m   <- male_gametes[1:n_row, ]
        male_gametes_f   <- female_gametes[1:n_row, ]
        female_gametes_m <- male_gametes[(n_row + 1):(n_row * 2), ]
        female_gametes_f <- female_gametes[(n_row + 1):(n_row * 2), ]

        genotypes_new[row_idx_male, ] <- male_gametes_m + male_gametes_f
        genotypes_new[row_idx_female, ] <- female_gametes_m + female_gametes_f
      }, .options = furrr::furrr_options(seed = TRUE))
      future::plan(future::sequential)
      private$.genotype_data <- genotypes_new
    },

    .check_compiled = function() {
      if (private$.compiled) cli::cli_abort("Cannot modify object after
        compilation. Create a new instance to make changes.")
    }
  )
)