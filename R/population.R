generate_phenotype <- function(config, snp_matrix) {
  #' Generate phenotype values from genotypes and environmental effects
  #'
  #' @param config Configuration list with simulation parameters
  #' @param snp_matrix Matrix of SNP genotypes (0, 1, or 2 copies of alt allele)
  #'
  #' @return A matrix with one column containing the phenotype values

  # Extract data from the config object
  pheno_name <- config$phenotype_name
  pheno_causal_snps <- config$phenotype_causal_snps
  pheno_heritability <- config$phenotype_heritability
  heritability_coeff <- (1 - pheno_heritability) / pheno_heritability

  # Compute contribution from genetic component
  genetic_effect <- snp_matrix %*% pheno_causal_snps
  env_variance <- heritability_coeff * var(genetic_effect)
  env_effect <- rnorm(config$n_pop, mean = 0, sd = sqrt(env_variance))

  # Assemble genetic and environmental effects
  phenotype <- genetic_effect + env_effect
  colnames(phenotype) <- pheno_name

  return(phenotype)
}

initialise_population <- function(config) {
  #' Initialize a population with genotypes and phenotypes
  #'
  #' @param config Configuration list with simulation parameters
  #'
  #' @return A data frame containing the initial population

  # Extract data from config object
  n_pop <- config$n_pop
  n_loci <- config$n_loci

  # Create sex vector (equal male/female split)
  sex <- rep(0:1, each = n_pop / 2)

  # Create SNP matrix based on minor allele frequencies
  snp_matrix <- matrix(
    rbinom(n_pop * n_loci, size = 2, prob = rep(config$snp_maf, each = n_pop)),
    nrow = n_pop,
    ncol = n_loci,
    byrow = FALSE
  )
  colnames(snp_matrix) <- paste0("rs", 1:n_loci)

  # Compute phenotype vector
  phenotype <- generate_phenotype(config, snp_matrix)

  # Combine into population data frame
  population <- cbind(sex, snp_matrix, phenotype)

  return(population)
}

produce_next_generation <- function(config, pop, matching) {
  #' Produce the next generation of offspring through mating
  #'
  #' @param config Configuration list with simulation parameters
  #' @param pop Data frame containing the current population
  #' @param matching Data frame with male_idx and female_idx columns
  #'
  #' @return A data frame containing the offspring population

  # Helper function to create offspring SNP matrix
  create_snp_matrix <- function(snps_male, snps_female) {
    n_pop <- config$n_pop
    n_loci <- config$n_loci

    # Calculate probability matrices for inheritance
    paternal_probs <- snps_male / 2
    maternal_probs <- snps_female / 2

    # Sample alleles from each parent
    paternal_alleles <- matrix(
      rbinom(n_pop * n_loci, 1, rep(as.vector(paternal_probs), each = 2)),
      nrow = n_pop,
      ncol = n_loci
    )

    maternal_alleles <- matrix(
      rbinom(n_pop * n_loci, 1, rep(as.vector(maternal_probs), each = 2)),
      nrow = n_pop,
      ncol = n_loci
    )

    # Combine alleles from both parents
    snps_offspring <- paternal_alleles + maternal_alleles
    colnames(snps_offspring) <- paste0("rs", 1:n_loci)

    return(snps_offspring)
  }

  # Create offspring sex vector (equal male/female split)
  sex <- rep(0:1, each = config$n_pop / 2)

  # Select genotype data in male and female subpopulations
  cols_select <- !(colnames(pop) %in% c("sex", config$phenotype_name))
  snps_male <- pop[matching$male_idx, cols_select]
  snps_female <- pop[matching$female_idx, cols_select]

  # Generate SNP matrix and phenotype
  snp_matrix <- create_snp_matrix(snps_male, snps_female)
  phenotype <- generate_phenotype(config, snp_matrix)

  # Assemble offspring population
  offspring_pop <- data.frame(sex, snp_matrix, phenotype)

  return(offspring_pop)
}
