#' Population Class for Genomic Simulation
#'
#' @description
#' The `Population` R6 class encapsulates the state and behavior of a simulated
#' genomic population. It manages genotype and phenotype data, supports
#' initialization, updating, and querying of individuals, and integrates with
#' phenotype and mating model classes for flexible simulation workflows.
#'
#' @details
#' The Population class simulates genomic populations with realistic genetic
#' structure. Generates correlated genotype data, supports multiple phenotypes,
#' and can advance through generations via mating simulations. The class uses
#' parallel processing and file-backed matrices for scalability. Once compiled,
#' the object becomes immutable.
#'
#' @section Methods:
#' - `initialize()`: Set up the initial population state.
#' - `genotype_data()`: Get the genotype data matrix.
#' - `phenotype_data()`: Get the phenotype data.
#' - `haplotype_data()`: Get the haplotype data.
#' - `recombination_rate()`: Get the recombination rate.
#' - `n_population()`: Get the number of individuals in the population.
#' - `n_loci()`: Get the number of loci.
#' - `mating_model()`: Get the mating model.
#' - `snp_mafs()`: Get the SNP minor allele frequencies.
#' - `phenotypes()`: Get the list of phenotypes.
#' - `generation()`: Get the current generation.
#' - `set_mating_model()`: Set the mating model.
#' - `add_phenotype()`: Add a phenotype to the population.
#' - `remove_phenotype()`: Remove a phenotype from the population.
#' - `compile()`: Compile and initialise the population data.
#'
#' @seealso
#' Other related classes and functions: `Phenotype`, `MatingModel`

Population <- R6::R6Class(
  "Population",

  public = list(

    #' @description Initialise a new Population object
    #' @param n_population Integer. The number of individuals in the population.
    #'   Must be even.
    #' @param n_loci Integer. The number of genetic loci to simulate.
    #' @param recombination_rate Numeric. The recombination rate between loci
    #'   (default: 1e-8).
    #' @param snp_dist Character. Distribution for SNP minor allele frequencies
    #'   ("beta" or "unif").
    #' @param snp_params Numeric vector of length 2. Parameters for the SNP
    #'   distribution.
    #' @param phenotypes Phenotype object or NULL. Optional phenotype to
    #'   associate with the population.
    #' @param mating_model MatingModel object or NULL. Optional mating model for
    #'   the population.
    #' @return A new Population object.
    initialize = function(n_population, n_loci, recombination_rate = 1e-8,
                          snp_dist = "beta", snp_params = c(5, 5),
                          phenotypes = NULL, mating_model = NULL) {
      checkmate::assert_count(n_population)
      checkmate::assert_count(n_loci)
      checkmate::assert_numeric(recombination_rate, lower = 0, upper = 1)
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
      private$.n_loci       <- n_loci
      private$.snp_mafs     <- snp_mafs
      private$.snp_dist     <- snp_dist
      private$.snp_params   <- snp_params
      private$.phenotypes   <- if (!is.null(phenotypes)) phenotypes else list()
      private$.mating_model <- if (!is.null(mating_model)) mating_model
    },

    #' @description Get the genotype data matrix
    #' @return Matrix or FBM object containing genotype data for all individuals
    genotype_data = function() return(private$.genotype_data),

    #' @description Get the phenotype data
    #' @return Data frame or matrix containing phenotype values for all
    #'   individuals
    phenotype_data = function() return(private$.phenotype_data),

    #' @description Get the haplotype data
    #' @return Matrix containing haplotype data for all individuals
    haplotype_data = function() return(private$.haplotype_data),

    #' @description Get the recombination rate
    #' @return Numeric value representing the recombination rate between loci
    recombination_rate = function() return(private$.recombination_rate),

    #' @description Get the number of individuals in the population
    #' @return Integer number of individuals in the population
    n_population = function() return(private$.n_population),

    #' @description Get the number of genetic loci
    #' @return Integer number of loci being simulated
    n_loci = function() return(private$.n_loci),

    #' @description Get the mating model
    #' @return MatingModel object or NULL if no mating model is set
    mating_model = function() return(private$.mating_model),

    #' @description Get the SNP minor allele frequencies
    #' @return Named numeric vector of minor allele frequencies for each SNP
    snp_mafs = function() return(private$.snp_mafs),

    #' @description Get the list of phenotypes
    #' @return Named list of Phenotype objects associated with the population
    phenotypes = function() return(private$.phenotypes),

    #' @description Get the current generation number
    #' @return Integer representing the current generation (0 for initial
    #'   generation)
    generation = function() return(private$.generation),

    #' @description Set the mating model for the population
    #' @param mating_model MatingModel object to be used for mating decisions
    #' @return No return value. Modifies the population's mating model.
    #' @details This method can only be called before the population is
    #'   compiled.
    set_mating_model = function(mating_model) {
      private$.check_compiled()
      checkmate::assert_class(mating_model, "MatingModel")
      private$.mating_model <- mating_model
    },

    #' @description Add a phenotype to the population
    #' @param phenotype Phenotype object to be added to the population
    #' @return No return value. Adds the phenotype to the population's list of
    #'   phenotypes.
    #' @details This method can only be called before the population is
    #'   compiled.
    add_phenotype = function(phenotype) {
      private$.check_compiled()
      checkmate::assert_class(phenotype, "Phenotype")
      private$.phenotypes[[phenotype$name()]] <- phenotype
    },

    #' @description Remove a phenotype from the population
    #' @param phenotype_name Character string. The name of the phenotype to
    #'   remove
    #' @return No return value. Removes the phenotype from the population's list
    #'   of phenotypes.
    #' @details This method can only be called before the population is
    #'   compiled.
    remove_phenotype = function(phenotype_name) {
      private$.check_compiled()
      phenotype_names <- names(self$phenotypes())
      checkmate::assert_choice(phenotype_name, choices = phenotype_names)
      private$.phenotypes[[phenotype_name]] <- NULL
    },

    #' @description Compile and initialise the population data
    #' @param verbose Logical. Whether to display progress messages during
    #'   compilation (default: TRUE)
    #' @return No return value. Compiles the population object and initialises
    #'   all data structures.
    #' @details This method must be called after all phenotypes and mating
    #'   models have been set up. It initialises the genotype data, generates
    #'   sibling pairs, and sets up phenotypes. Once compiled, the population
    #'   object cannot be modified.
    compile = function(verbose = TRUE) {
      private$.check_compiled()
      cli::cli_h1("Compiling Population object:")
      cli::cli_h2("Initialising genotype data")
      private$.initialise_genotypes()
      cli::cli_h2("Generating sibling pairs")
      private$.update()
      cli::cli_h2("Initialising haplotype data")
      private$.initialise_haplotypes()
      cli::cli_h2("Setting up phenotypes")
      private$.setup_phenotypes()
      private$.compiled <- TRUE
      cli::cli_alert_info("Done!")
    }
  ),

  private = list(
    # Fields
    .genotype_data      = NULL,
    .phenotype_data     = NULL,
    .haplotype_data     = NULL,
    .recombination_rate = NULL,
    .n_population       = NULL,
    .n_loci             = NULL,
    .snp_mafs           = NULL,
    .snp_dist           = NULL,
    .snp_params         = NULL,
    .summariser         = NULL,
    .phenotypes         = NULL,
    .mating_model       = NULL,
    .compiled           = FALSE,
    .generation         = 0,

    .initialise_manual = function() {
      # Initialise using parameters such as `n_population`, `n_loci`,
      # `snp_mafs`, etc.
    },

    .initialise_from_data = function() {
      # Initialise using LD matrices and SNP files
    },

    # @description Set up phenotypes for the population
    # @details This private method initialises causal SNPs and effects for each
    #   phenotype associated with the population. It randomly selects causal
    #   SNPs if not specified and computes phenotype values based on genotype
    #   data.
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
          phenotype$set_causal_snps(causal_idx)
          phenotype$set_causal_effects(causal_eff)
        }
        causal_idx  <- sort(phenotype$causal_snps())
        phenotype$set_causal_snps(causal_idx)
        phenotype$compute_data(self$genotype_data())
      }
    },

    # @description Initialise genotype data for the population
    # @param ncores Integer. Number of CPU cores to use for parallel processing
    #   (default: detected cores - 1)
    # @param seed Integer or NULL. Random seed for reproducible generation
    #   (default: NULL)
    # @return FBM object containing the initialised genotype data
    # @details This private method creates the initial genotype data matrix
    #   using parallel processing. It generates correlated genotype data for
    #   males and females separately, ensuring realistic population structure
    #   and linkage disequilibrium patterns.
    .initialise_genotypes = function(ncores = parallel::detectCores() - 1,
                                     seed = NULL) {
      # Add validation to ensure even population size
      checkmate::assert(self$n_population() %% 2 == 0,
                        .var.name = "even population size")

      n_sex      <- self$n_population() / 2
      chunk_size <- 10000

      genotypes_fbm <- bigstatsr::FBM(
        nrow        = self$n_population(),
        ncol        = self$n_loci(),
        type        = "unsigned short",
        backingfile = "genotypes_gen0"
      )

      # Create a list of chunks to fill the matrix
      start_rows <- seq(1, n_sex, by = chunk_size)
      end_rows   <- pmin(start_rows + chunk_size - 1, n_sex)
      chunk_data <- list(start = start_rows, end = end_rows)

      future::plan(future::multisession, workers = ncores)

      cli::cli_alert_info("Generating male genotype data...")
      furrr::future_walk(seq_along(chunk_data$start), ~ {
        start_row <- chunk_data$start[.x]
        end_row   <- chunk_data$end[.x]
        rows      <- end_row - start_row + 1

        male_chunk <- matrix(
          rbinom(rows * self$n_loci(),
            size = 2,
            prob = rep(self$snp_mafs(), each = rows)
          ),
          nrow = rows,
          ncol = self$n_loci()
        )

        genotypes_fbm[start_row:end_row, ] <- male_chunk
      }, .options = furrr::furrr_options(seed = TRUE))

      cli::cli_alert_info("Computing correlation matrix for genotype data...")
      bigparallelr::set_blas_ncores(ncores)

      cor_male   <- bigstatsr::big_cor(genotypes_fbm, ind.row = 1:n_sex)
      chunk_data <- list(start = start_rows + n_sex, end = end_rows + n_sex)

      cli::cli_alert_info("Generating female genotype data...")
      furrr::future_walk(seq_along(chunk_data$start), ~ {
        i1           <- chunk_data$start[.x]
        i2           <- chunk_data$end[.x]
        n_row        <- i2 - i1 + 1
        scores_chunk <- MASS::mvrnorm(
          n_row,
          mu    = rep(0, self$n_loci()),
          Sigma = cor_male[]
        )

        p0 <- qnorm((1 - self$snp_mafs())^2)
        p1 <- qnorm((1 - self$snp_mafs())^2 + 2 * self$snp_mafs() *
                      (1 - self$snp_mafs()))

        p0_matrix <- matrix(
          p0,
          nrow = nrow(scores_chunk),
          ncol = self$n_loci(),
          byrow = TRUE
        )

        p1_matrix <- matrix(
          p1,
          nrow = nrow(scores_chunk),
          ncol = self$n_loci(),
          byrow = TRUE
        )

        female_chunk <- (scores_chunk > p0_matrix) + (scores_chunk > p1_matrix)

        genotypes_fbm[i1:i2, ] <- female_chunk
      }, .options = furrr::furrr_options(seed = TRUE))

      future::plan(future::sequential)
      private$.genotype_data <- genotypes_fbm
    },

    .initialise_haplotypes = function(ncores = parallel::detectCores() - 1,
                                      seed = NULL) {
      n_pop  <- self$n_population()
      n_loci <- self$n_loci()
      haplo_fbm <- bigstatsr::FBM(
        nrow = n_pop,
        ncol = 2 * n_loci,
        type = "unsigned short",
        backingfile = "haplotypes_gen0"
      )

      chunk_size <- 10000
      start_rows <- seq(1, n_pop, by = chunk_size)
      end_rows   <- pmin(start_rows + chunk_size - 1, n_pop)
      chunk_data <- list(start = start_rows, end = end_rows)

      future::plan(future::multisession, workers = ncores)

      furrr::future_walk(seq_along(chunk_data$start), function(idx) {
        i1         <- chunk_data$start[idx]
        i2         <- chunk_data$end[idx]
        rows       <- i2 - i1 + 1
        geno_chunk <- self$genotype_data()[i1:i2, , drop = FALSE]

        # Vectorised haplotype assignment
        # For g == 0: both haplotypes 0
        # For g == 2: both haplotypes 1
        # For g == 1: randomly assign 1/0 or 0/1
        hap1 <- matrix(0L, nrow = rows, ncol = n_loci)
        hap2 <- matrix(0L, nrow = rows, ncol = n_loci)

        idx0 <- geno_chunk == 0
        idx2 <- geno_chunk == 2
        idx1 <- geno_chunk == 1

        hap1[idx2] <- 1L
        hap2[idx2] <- 1L

        if (any(idx1)) {
          # For heterozygotes, randomly assign which haplotype gets the 1
          rand    <- matrix(runif(sum(idx1)), nrow = sum(idx1))
          assign1 <- rand < 0.5
          # Assign 1 to hap1, 0 to hap2 where assign1 is TRUE, else reverse
          hap1[idx1][assign1]  <- 1L
          hap2[idx1][!assign1] <- 1L
        }

        # Combine haplotypes into one matrix for this chunk
        hap_chunk <- cbind(hap1, hap2)
        haplo_fbm[i1:i2, ] <- hap_chunk
      }, .options = furrr::furrr_options(seed = TRUE))

      future::plan(future::sequential)
      private$.haplotype_data <- haplo_fbm
    },

    # @description Update the population to the next generation
    # @param ncores Integer. Number of CPU cores to use for parallel processing
    #   (default: detected cores - 1)
    # @return No return value. Updates the population's genotype data for the
    #   next generation.
    # @details This private method simulates mating and reproduction to advance
    #   the population to the next generation. It creates gametes from parent
    #   genotypes and combines them to form offspring genotypes, using parallel
    #   processing for efficiency.
    .update = function(ncores = parallel::detectCores() - 1) {
      private$.generation <- private$.generation + 1
      genotypes_new <- bigstatsr::FBM(
        nrow = self$n_population(),
        ncol = self$n_loci(),
        type = "unsigned short",
        backingfile = paste0("genotypes_gen", self$generation())
      )

      n_sex <- self$n_population() / 2
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

        chunk_male     <- self$genotype_data()[row_idx_male, ]
        chunk_female   <- self$genotype_data()[row_idx_female, ]
        prob_male      <- rep(as.vector(chunk_male / 2), times = 2)
        prob_female    <- rep(as.vector(chunk_female / 2), times = 2)

        male_gametes <- matrix(
          rbinom(
            n_row * self$n_loci() * 2,
            size = 1,
            prob = prob_male
          ),
          nrow = n_row * 2,
          ncol = self$n_loci(),
        )

        female_gametes <- matrix(
          rbinom(
            n_row * self$n_loci() * 2,
            size = 1,
            prob = prob_female
          ),
          nrow = n_row * 2,
          ncol = self$n_loci()
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

    # @description Check if the population has been compiled
    # @return No return value. Throws an error if the population has already
    #   been compiled.
    # @details This private method is used to prevent modification of the
    #   population object after it has been compiled. It throws an error
    #   message instructing users to create a new instance to make changes.
    .check_compiled = function() {
      if (private$.compiled) cli::cli_abort("Cannot modify object after
        compilation. Create a new instance to make changes.")
    }
  )
)