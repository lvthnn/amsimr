#' Population class
#'
#' Represents a population for performing simulations of xAM, allowing for
#' configuration via constructor arguments.
#'
#' @importFrom R6 R6Class
#' @importFrom checkmate assert_true assert_count assert_choice assert_numeric
#' @importFrom checkmate assert_list assert_class test_count
#' @importFrom cli cli_alert_info
#' @importFrom stats var rnorm rbinom
#'
#' @export
Population <- R6Class(
  "Population",

  public = list(
    #' @field data Data frame containing SNP and phenotype data
    data         = NULL,

    #' @field n_pop Number of individuals in the population
    n_pop        = NULL,

    #' @field n_loci Number of genetic loci (SNPs) modelled in the genome
    n_loci       = NULL,

    #' @field phenotypes A list of Phenotype objects
    phenotypes   = NULL,

    #' @field maf_dist One of 'beta' or 'uniform'
    maf_dist     = NULL,

    #' @field maf_params Parameters for MAF distribution
    maf_params   = NULL,

    #' @field maf Minor allele frequencies for SNPs
    maf          = NULL,

    #' @description
    #' Instantiate a new `Population` object.
    #'
    #' @param n_pop Number of individuals in the population.
    #' @param n_loci Number of loci in the genome.
    #' @param maf_dist Minor allele frequency distribution type ("uniform" or
    #'  "beta").
    #' @param maf_params Parameters for the minor allele frequency distribution.
    #' @param phenotypes List of Phenotype objects to model in the population.
    #'
    #' @return A new `Population` object.
    initialize = function(n_pop, n_loci, maf_dist = NULL, maf_params = NULL,
                          phenotypes = NULL) {

      # Verify "n_" parameters
      assert_true(test_count(n_pop, positive = TRUE))
      if (n_pop %% 2 != 0) {
        n_pop <- n_pop - 1
        cli_alert_info("Updated parameter 'n_pop' to even value: {n_pop}")
      }
      assert_count(n_loci, positive = TRUE)

      # Verify MAF distribution
      if (is.null(maf_dist) && is.null(maf_params)) {
        maf_dist <- "beta"
        maf_params <- c(0.5, 1.5)
        cli_alert_info(c("Using default values '{maf_dist}' and '{maf_params}'",
                         "for 'maf_dist' and 'maf_params'"))
      }

      assert_choice(maf_dist, c("uniform", "beta"))

      if (maf_dist == "uniform") {
        assert_numeric(
          maf_params,
          len = 2,
          lower = 0,
          upper = 1,
          any.missing = FALSE
        )
      } else {
        assert_numeric(
          maf_params,
          len = 2,
          lower = 0,
          upper = Inf,
          finite = TRUE,
          any.missing = FALSE
        )
      }

      # Verify phenotypes object
      assert_list(phenotypes, types = "Phenotype", min.len = 1, null.ok = TRUE)

      if (is.null(phenotypes)) {
        cli_alert_info("Initialising default (empty) phenotype vector")
        self$phenotypes <- c()
      }

      # Initialise the population
      self$n_pop <- n_pop
      self$n_loci <- n_loci
      self$maf_dist <- maf_dist
      self$maf_params <- maf_params
      self$maf <- if (maf_dist == "uniform") {
        runif(n_loci, min = maf_params[1], max = maf_params[2])
      } else {
        rbeta(n_loci, shape1 = maf_params[1], shape2 = maf_params[2])
      }

      self$data <- private$initialise_population()
    },

    #' @description
    #' Add new phenotypes to model in the population
    #'
    #' @param phenotypes A list of Phenotype objects to add to the population
    add_phenotypes = function(phenotypes) {
      assert_class(phenotypes, classes = "Phenotype")
      self$phenotypes <- append(self$phenotypes, phenotypes)
      cli_alert_info("Adding new phenotypes:")
      print(self$phenotypes)
    },

    #' @description
    #' Remove phenotypes modelled in the population
    #' 
    #' @param phenotypes Names of phenotypes to remove from the population
    remove_phenotypes = function(phenotypes) {
      print("Yo!")
    }
  ),

  private = list(
    #' Initialize a population with genotypes and phenotypes
    #'
    #' @return A data frame containing the initial population
    #'
    #' @noRd
    initialise_population = function() {
      # Create sex vector (equal male/female split)
      sex <- rep(0:1, each = self$n_pop / 2)
      sibling_id <- rep(1:self$n_pop / 2, each = 2)

      # Create SNP matrix based on minor allele frequencies
      snps_paternal <- matrix(
        rbinom(
          self$n_pop * self$n_loci / 2,
          prob = rep(self$maf, each = self$n_pop / 2),
          size = 2
        ),
        nrow = self$n_pop / 2,
        ncol = self$n_loci,
        byrow = FALSE
      )

      snps_maternal <- matrix(
        rbinom(
          self$n_pop * self$n_loci / 2,
          prob = rep(self$maf, each = self$n_pop / 2),
          size = 2
        ),
        nrow = self$n_pop / 2,
        ncol = self$n_loci,
        byrow = FALSE
      )

      snp_matrix <- private$generate_snp_matrix(snps_paternal, snps_maternal)
      colnames(snp_matrix) <- paste0("rs", 1:self$n_loci)

      generate_phenotypes()
    },

    #' Generate a SNP matrix for sibling pairs given paternal and maternal SNPs
    #
    #' @param self The Population object
    #' @param snps_paternal Matrix of paternal SNP genotypes
    #' @param snps_maternal Matrix of maternal SNP genotypes
    #'
    #' @return A matrix of sibling pair genotypes
    #'
    #' @noRd
    generate_snp_matrix = function(snps_paternal, snps_maternal) {
      prob_paternal <- as.matrix(snps_paternal) / 2
      prob_maternal <- as.matrix(snps_maternal) / 2

      # For each pair of siblings
      haplotype_paternal <- matrix(
        rbinom(
          self$n_pop * self$n_loci,
          prob = rep(as.vector(prob_paternal), each = 2),
          size = 1
        ),
        nrow = self$n_pop,
        ncol = self$n_loci
      )

      haplotype_maternal <- matrix(
        rbinom(
          self$n_pop * self$n_loci,
          prob = rep(as.vector(prob_paternal), each = 2),
          size = 1
        ),
        nrow = self$n_pop,
        ncol = self$n_loci
      )

      snps_offspring <- haplotype_paternal + haplotype_maternal
      colnames(snps_offspring) <- paste0("rs", 1:self$n_loci)
      invisible(snps_offspring)
    },

    #' Generate phenotype values from genotypes and environmental effects
    #
    #' @param config Configuration list with simulation parameters
    #' @param snp_matrix Matrix of SNP genotypes (0, 1, or 2 copies of allele)
    #' @param env_variance Initial environmental variance
    #
    #' @return A matrix with one column containing the phenotype values
    #'
    #' @noRd
    generate_phenotypes = function() {
      if (is.null(self$phenotypes)) {
        cli_abort("No phenotypes defined.")
      }

      for (phenotype in phenotypes) {
        print(phenotype)
      }
    }
  )
)