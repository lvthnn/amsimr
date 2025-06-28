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
#' structure. It generates correlated genotype data, supports multiple phenotypes,
#' and can advance through generations via mating simulations. The class uses
#' parallel processing and file-backed matrices for scalability. Once compiled,
#' the object becomes immutable.
#'
#' @section Methods:
#' - `initialize()`: Set up the initial population state.
#' - `add_phenotype()`: Add a phenotype to the population.
#' - `remove_phenotype()`: Remove a phenotype from the population.
#' - `compile()`: Compile and initialise the population data.
#' - `update()`: Advance the population one generation with mate matching.
#'
#' @seealso
#' Other related classes and functions: `Phenotype`, `MatingModel`
#'
#' @export
Population <- R6::R6Class(
  "Population",

  public = list(

    #' @description Initialise a new Population object
    #' @param size Integer. The number of individuals in the population.
    #'   Must be even.
    #' @param n_loci Integer. The number of genetic loci to simulate.
    #' @param recomb_rate Numeric. The recombination rate between loci
    #'   (default: 1e-8).
    #' @param mutation_rate Numeric. The probability of mutation at loci.
    #'   (default: 1e-8)
    #' @param snp_dist Character. Distribution for SNP minor allele frequencies
    #'   ("beta" or "unif").
    #' @param snp_params Numeric vector of length 2. Parameters for the SNP
    #'   distribution.
    #' @param phenotypes Phenotype object or NULL. Optional phenotype to
    #'   associate with the population.
    #' @param mating_model MatingModel object or NULL. Optional mating model for
    #'   the population.
    #' @return A new Population object.
    initialize = function(size, n_loci, recomb_rate = 1e-8,
                          mutation_rate = 1e-8, snp_dist = "beta",
                          snp_params = c(5, 5), phenotypes = NULL,
                          mating_model = NULL) {
      checkmate::assert_count(size)
      checkmate::assert_count(n_loci)
      checkmate::assert_numeric(recomb_rate, lower = 0, upper = 1)
      checkmate::assert_numeric(mutation_rate, lower = 0, upper = 1)
      checkmate::assert_choice(snp_dist, choices = c("beta", "unif"))
      checkmate::assert_numeric(snp_params, len = 2)
      checkmate::assert_class(phenotypes, "Phenotype", null.ok = TRUE)
      checkmate::assert_class(mating_model, "MatingModel", null.ok = TRUE)

      rdist           <- get(paste0("r", snp_dist), mode = "function")
      param1          <- snp_params[[1]]
      param2          <- snp_params[[2]]
      snp_mafs        <- rdist(n_loci, param1, param2)

      checkmate::assert_numeric(snp_mafs, lower = 0, upper = 1, finite = TRUE,
                                any.missing = FALSE, all.missing = FALSE)

      # Assume there are 23 chromosomes and split SNPs evenly between them
      chrom_idx              <- seq.int(1, n_loci, by = floor(n_loci / 23))
      recomb_rate            <- rep(recomb_rate, n_loci)
      recomb_rate[chrom_idx] <- 0.5

      private$.size          <- size
      private$.n_loci        <- n_loci
      private$.recomb_rate   <- recomb_rate
      private$.mutation_rate <- mutation_rate
      private$.snp_mafs      <- snp_mafs
      private$.snp_dist      <- snp_dist
      private$.snp_params    <- snp_params
      private$.phenotypes    <- if (!is.null(phenotypes)) phenotypes else list()
      private$.mating_model  <- if (!is.null(mating_model)) mating_model
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
      private$.phenotypes[[phenotype$name]] <- phenotype
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
      phenotype_names <- names(self$phenotypes)
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
      private$.initial_update()
      cli::cli_h2("Initialising haplotype data")
      private$.initialise_haplotypes()
      cli::cli_h2("Setting up phenotypes")
      private$.setup_phenotypes()
      private$.compiled <- TRUE
      cli::cli_alert_info("Done!")
    },

    #' @description Advance the population one generation with mate matching
    #'
    #' @param ncores Integer. Number of cores to use for parallel computing.
    #' @param overwrite Logical. Should existing output files be overwritten?
    update = function(ncores = parallel::detectCores() - 1, overwrite = TRUE) {
      size   <- self$size
      n_loci <- self$n_loci()

      # Generate the mate matching
      matching <- private$.generate_mate_matching()

      # Generate recombination and mutation masks
      recomb_unif <- matrix(runif(size * n_loci), nrow = size, ncol = n_loci)
      mut_unif    <- matrix(runif(size * n_loci), nrow = size, ncol = n_loci)

      recomb_mask   <- sweep(recomb_unif, 2, self$recomb_rate, FUN = "<") + 0L
      mutation_mask <- (mut_unif < self$mutation_rate) + 0L

      # Retrieve the current haplotype data and split into paternal and maternal
      paternal  <- self$haplotypes[, 1:n_loci]
      maternal  <- self$haplotypes[, (n_loci + 1):(2 * n_loci)]
      gametes   <- (1 - recomb_mask) * paternal + recomb_mask * maternal
      gametes   <- (gametes + mutation_mask) %% 2

      new_haplotypes         <- cbind(gametes, gametes[matching, ])
      private$.hapotypes[, ] <- new_haplotypes

      # Generate offspring genotypes from the new haplotypes
      offspring_genotypes <- new_haplotypes[, 1:n_loci, drop = FALSE] +
        new_haplotypes[, (n_loci + 1):(2 * n_loci), drop = FALSE]

      private$.genotypes[, ] <- offspring_genotypes
      private$.generation   <- private$.generation + 1
    }
  ),

  private = list(
    # Fields
    .genotypes     = NULL,
    .hapotypes     = NULL,
    .recomb_rate   = NULL,
    .mutation_rate = NULL,
    .size          = NULL,
    .n_loci        = NULL,
    .snp_mafs      = NULL,
    .snp_dist      = NULL,
    .snp_params    = NULL,
    .summariser    = NULL,
    .phenotypes    = NULL,
    .mating_model  = NULL,
    .output_dir    = NULL,
    .compiled      = FALSE,
    .generation    = 0,

    # @description Set up phenotypes for the population
    .setup_phenotypes = function() {
      for (phenotype in self$phenotypes) {
        if (is.null(phenotype$causal_snps)) {
          if (phenotype$n_loci > self$n_loci) {
            cli::cli_abort(c("x" = "Number of causal SNPs specified for
                                    phenotype `{phenotype$name()}` exceeds
                                    number of SNPs modelled in the
                                    population."))
          }
          sd_eff     <- sqrt(phenotype$heritability / phenotype$n_loci)
          causal_idx <- sample(1:self$n_loci, size = phenotype$n_loci)
          causal_eff <- rnorm(phenotype$n_loci, mean = 0, sd = sd_eff)
          phenotype$causal_snps    <- causal_idx
          phenotype$causal_effects <- causal_eff
        }
        causal_idx  <- sort(phenotype$causal_snps)
        phenotype$causal_snps <- causal_idx
        phenotype$compute_values(self$genotypes)
      }
    },

    # @description Initialise genotype data for the population
    .initialise_genotypes = function(ncores = parallel::detectCores() - 1,
                                     seed = NULL) {
      checkmate::assert(self$size %% 2 == 0, .var.name = "even population size")

      n_sex      <- self$size / 2
      chunk_size <- 10000

      genotypes_fbm <- bigstatsr::FBM(
        nrow        = self$size,
        ncol        = self$n_loci,
        type        = "unsigned short",
        backingfile = "genotypes_gen0"
      )

      start_rows <- seq.int(1, n_sex, by = chunk_size)
      end_rows   <- pmin(start_rows + chunk_size - 1, n_sex)
      chunk_data <- list(start = start_rows, end = end_rows)

      future::plan(future::multisession, workers = ncores)

      cli::cli_alert_info("Generating male genotype data...")
      furrr::future_walk(seq_along(chunk_data$start), ~ {
        start_row <- chunk_data$start[.x]
        end_row   <- chunk_data$end[.x]
        rows      <- end_row - start_row + 1

        male_chunk <- matrix(
          rbinom(rows * self$n_loci,
            size = 2,
            prob = rep(self$snp_mafs, each = rows)
          ),
          nrow = rows,
          ncol = self$n_loci
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

        # TODO: Use Cholesky decomposition to generate MVN instead of
        #       MASS::mvrnorm, since the latter computes the eigendecomposition
        #       of cor_male[] for every block.

        scores_chunk <- MASS::mvrnorm(
          n_row,
          mu    = rep(0, self$n_loci),
          Sigma = cor_male[]
        )

        p0 <- qnorm((1 - self$snp_mafs)^2)
        p1 <- qnorm((1 - self$snp_mafs)^2 +
                      2 * self$snp_mafs *
                        (1 - self$snp_mafs))

        p0_matrix <- matrix(p0, nrow = nrow(scores_chunk), ncol = self$n_loci,
                            byrow = TRUE)

        p1_matrix <- matrix(p1, nrow = nrow(scores_chunk), ncol = self$n_loci,
                            byrow = TRUE)

        female_chunk <- (scores_chunk > p0_matrix) + (scores_chunk > p1_matrix)

        genotypes_fbm[i1:i2, ] <- female_chunk
      }, .options = furrr::furrr_options(seed = TRUE))

      future::plan(future::sequential)
      private$.genotypes <- genotypes_fbm
    },

    .initialise_haplotypes = function(ncores = parallel::detectCores() - 1,
                                      seed = NULL) {
      size   <- self$size
      n_loci <- self$n_loci
      haplo_fbm <- bigstatsr::FBM(
        nrow = size,
        ncol = 2 * n_loci,
        type = "unsigned short",
        backingfile = "haplotypes_gen0"
      )

      chunk_size <- 10000
      start_rows <- seq.int(1, size, by = chunk_size)
      end_rows   <- pmin(start_rows + chunk_size - 1, size)
      chunk_data <- list(start = start_rows, end = end_rows)

      future::plan(future::multisession, workers = ncores)

      furrr::future_walk(seq_along(chunk_data$start), function(idx) {
        i1         <- chunk_data$start[idx]
        i2         <- chunk_data$end[idx]
        rows       <- i2 - i1 + 1
        geno_chunk <- self$genotypes[i1:i2, , drop = FALSE]

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
          rand    <- matrix(runif(sum(idx1)), nrow = sum(idx1))
          assign1 <- rand < 0.5
          hap1[idx1][assign1]  <- 1L
          hap2[idx1][!assign1] <- 1L
        }

        hap_chunk <- cbind(hap1, hap2)
        haplo_fbm[i1:i2, ] <- hap_chunk
      }, .options = furrr::furrr_options(seed = TRUE))

      future::plan(future::sequential)
      private$.hapotypes <- haplo_fbm
    },

    .initial_update = function(ncores = parallel::detectCores() - 1) {
      private$.generation <- private$.generation + 1
      genotypes_new <- bigstatsr::FBM(
        nrow = self$size,
        ncol = self$n_loci,
        type = "unsigned short",
        backingfile = paste0("genotypes_gen", self$generation)
      )

      n_sex <- self$size / 2
      chunk_size <- 10000

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

        chunk_male     <- self$genotypes[row_idx_male, ]
        chunk_female   <- self$genotypes[row_idx_female, ]
        prob_male      <- rep(as.vector(chunk_male / 2), times = 2)
        prob_female    <- rep(as.vector(chunk_female / 2), times = 2)

        male_gametes <- matrix(
          rbinom(
            n_row * self$n_loci * 2,
            size = 1,
            prob = prob_male
          ),
          nrow = n_row * 2,
          ncol = self$n_loci
        )

        female_gametes <- matrix(
          rbinom(
            n_row * self$n_loci * 2,
            size = 1,
            prob = prob_female
          ),
          nrow = n_row * 2,
          ncol = self$n_loci
        )

        male_gametes_m   <- male_gametes[1:n_row, ]
        male_gametes_f   <- female_gametes[1:n_row, ]
        female_gametes_m <- male_gametes[(n_row + 1):(n_row * 2), ]
        female_gametes_f <- female_gametes[(n_row + 1):(n_row * 2), ]

        genotypes_new[row_idx_male, ] <- male_gametes_m + male_gametes_f
        genotypes_new[row_idx_female, ] <- female_gametes_m + female_gametes_f
      }, .options = furrr::furrr_options(seed = TRUE))

      future::plan(future::sequential)
      private$.genotypes <- genotypes_new
    },

    .check_compiled = function() {
      if (private$.compiled) cli::cli_abort("Cannot modify object after
        compilation. Create a new instance to make changes.")
    },

    .generate_mate_matching = function() {
      res <- sample.int(self$size, size = self$size)
      return(res)
    }
  ),

  active = list(
    #' @field genotypes Get the genotype data matrix.
    genotypes = function() return(private$.genotypes),

    #' @field haplotypes Get the haplotype data.
    haplotypes = function() return(private$.hapotypes),

    #' @field recomb_rate Get the recombination rate.
    recomb_rate = function() return(private$.recomb_rate),

    #' @field mutation_rate Get the mutation rate.
    mutation_rate = function() return(private$.mutation_rate),

    #' @field size Get the number of individuals in the population.
    size = function() return(private$.size),

    #' @field n_loci Get the number of loci.
    n_loci = function() return(private$.n_loci),

    #' @field mating_model Get or set the mating model.
    mating_model = function(value) {
      if (missing(value)) {
        return(private$.mating_model)
      } else {
        checkmate::assert_class(value, "MatingModel")
        private$.mating_model <- value
      }
    },

    #' @field phenotypes Get or set the list of phenotypes.
    phenotypes = function(value) {
      if (missing(value)) {
        return(private$.phenotypes)
      } else {
        checkmate::assert_list(value, types = "MatingModel",
                               names = sapply(value, function(v) v$name))
        private$.phenotypes <- value
      }
    },

    #' @field snp_mafs Get the SNP minor allele frequencies.
    snp_mafs = function() return(private$.snp_mafs),

    #' @field generation Get the current generation.
    generation = function() return(private$.generation)
  )
)
