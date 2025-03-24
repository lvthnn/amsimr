generate_phenotype <- function(config, snp_matrix) {
  # Extract data from the config object
  phenotype_name <- config$phenotype_name
  phenotype_causal_snps <- config$phenotype_causal_snps
  phenotype_heritability <- config$phenotype_heritability

  # Compute contribution from genetic component
  genetic_effect <- snp_matrix %*% phenotype_causal_snps
  heritability_coeff <- (1 - phenotype_heritability) / phenotype_heritability
  environment_variance <- heritability_coeff * var(genetic_effect)
  environment_effect <- rnorm(config$n_pop,
                              mean = 0,
                              sd = sqrt(environment_variance))

  # Assemble genetic and environmental effects
  phenotype <- genetic_effect + environment_effect
  colnames(phenotype) <- ifelse(!is.null(phenotype_name), phenotype_name, "X")

  return(phenotype)
}

initialise_population <- function(config) {
  # Extract data from the config object
  n_pop <- config$n_pop
  n_loci <- config$n_loci
  snp_maf <- config$snp_maf

  # Generate a vector of individual sexes
  sex <- rep(0:1, each = config$n_pop / 2)

  # Create allele count matrix for the population
  snp_matrix <- matrix(rbinom(n_pop * n_loci, size = 2,
                              prob = rep(snp_maf, each = n_pop)),
                       nrow = n_pop, ncol = n_loci, byrow = FALSE)

  colnames(snp_matrix) <- paste0("rs", 1:n_loci)

  # Compute the genetic contribution, tune heritability with random
  # normal variables
  phenotype <- generate_phenotype(config, snp_matrix)

  # Construct population data frame
  population <- cbind(sex, snp_matrix, phenotype)

  return(population)
}

produce_next_generation <- function(config, pop, matching) {
  # Retrieve some variables used in this function
  #
  # We expect `config` and `pop` to be as usual, and matching to be
  # a data frame in the form

  create_snp_matrix <- function(snps_male, snps_female) {
    # Number of offspring
    n_offspring <- config$n_pop
    n_loci <- config$n_loci

    # Each parent contributes one allele
    snps_offspring <- matrix(0, nrow = n_offspring, ncol = n_loci)

    for (i in 1:n_loci) {
      paternal_allele <- rbinom(n_offspring, 1, snps_male[, i] / 2)
      maternal_allele <- rbinom(n_offspring, 1, snps_female[, i] / 2)
      snps_offspring[, i] <- paternal_allele + maternal_allele
    }

    colnames(snps_offspring) <- paste0("rs", 1:n_loci)
    return(snps_offspring)
  }

  # Vector of sexes for the offspring
  sex <- rep(0:1, each = config$n_pop / 2)

  # Select genotype data in male and female subpopulations
  cols_select <- !(colnames(pop) %in% c("sex", config$phenotype_name))
  snps_male <- pop[matching$male_idx, cols_select]
  snps_female <- pop[matching$female_idx, cols_select]
  snp_matrix <- create_snp_matrix(snps_male, snps_female)

  # Generate phenotype
  phenotype <- generate_phenotype(config, snp_matrix)

  # Assemble the offspring population 
  offspring_pop <- data.frame(sex, snp_matrix, phenotype)

  return(offspring_pop)
}
