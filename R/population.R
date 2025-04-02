#' Generate a SNP matrix for sibling pairs given paternal and maternal SNPs
#'
#' @param config Configuration list with simulation parameters
#' @param snps_paternal Matrix of paternal SNP genotypes
#' @param snps_maternal Matrix of maternal SNP genotypes
#' @return A matrix of sibling pair genotypes
#'
#' @importFrom stats rbinom
#'
#' @noRd
generate_snp_matrix <- function(config, snps_paternal, snps_maternal) {
  n_pop <- config$n_pop
  n_loci <- config$n_loci

  # Calculate probability matrices for inheritance
  prob_paternal <- as.matrix(snps_paternal) / 2
  prob_maternal <- as.matrix(snps_maternal) / 2

  # For each pair of siblings
  haplotype_paternal <- matrix(
    rbinom(
      n_pop * n_loci,
      size = 1,
      prob = rep(as.vector(prob_paternal), each = 2)
    ),
    nrow = n_pop,
    ncol = n_loci
  )

  haplotype_maternal <- matrix(
    rbinom(
      n_pop * n_loci,
      size = 1,
      prob = rep(as.vector(prob_maternal), each = 2)
    ),
    nrow = n_pop,
    ncol = n_loci
  )

  snps_offspring <- haplotype_paternal + haplotype_maternal
  colnames(snps_offspring) <- paste0("rs", 1:n_loci)
  return(snps_offspring)
}

#' Generate phenotype values from genotypes and environmental effects
#'
#' @param config Configuration list with simulation parameters
#' @param snp_matrix Matrix of SNP genotypes (0, 1, or 2 copies of alt allele)
#'
#' @return A matrix with one column containing the phenotype values
#'
#' @importFrom stats var rnorm
#' @noRd
generate_phenotype <- function(config, snp_matrix) {
  # Extract data from the config object
  pheno_name <- config$phenotype_name
  pheno_causal_snps <- config$phenotype_causal_snps
  pheno_heritability <- config$phenotype_heritability
  heritability_coeff <- (1 - pheno_heritability) / pheno_heritability

  # Compute contribution from genetic component
  genetic_effect <- scale(snp_matrix) %*% pheno_causal_snps
  env_variance <- heritability_coeff * var(genetic_effect)
  env_effect <- rnorm(config$n_pop, mean = 0, sd = sqrt(env_variance))

  # Assemble genetic and environmental effects
  phenotype <- genetic_effect + env_effect
  phenotype_raw <- genetic_effect

  colnames(phenotype) <- pheno_name
  colnames(phenotype_raw) <- paste0(pheno_name, "_raw")

  return(cbind(phenotype, phenotype_raw))
}

#' Produce the next generation of offspring based on a mate matching
#'
#' @param config Configuration list with simulation parameters
#' @param population Data frame containing the current population
#' @param matching Data frame with male_idx and female_idx columns
#'
#' @return A data frame containing the offspring population
#'
#' @importFrom stats rbinom
#'
#' @noRd
generate_offspring <- function(config, population, matching) {
  # Create offspring sex vector (equal male/female split)
  n_pairs <- config$n_pop / 2
  sex <- rep(0:1, times = n_pairs)
  sibling_id <- rep(1:n_pairs, each = 2)

  # Select genotype data in male and female subpopulations
  cols_select <- !(colnames(population) %in% c("sex", config$phenotype_name))
  snps_paternal <- population[matching$male_idx, cols_select]
  snps_maternal <- population[matching$female_idx, cols_select]

  # Generate SNP matrix and phenotype
  snp_matrix <- generate_snp_matrix(config, snps_paternal, snps_maternal)
  phenotype <- generate_phenotype(config, snp_matrix)

  # Assemble offspring population
  offspring_pop <- data.frame(sex, snp_matrix, phenotype, sibling_id)

  return(offspring_pop)
}

#' Initialize a population with genotypes and phenotypes
#'
#' @param config Configuration list with simulation parameters
#'
#' @return A data frame containing the initial population
#' @importFrom stats rbinom
#'
#' @export
initialise_population <- function(config) {
  # Extract data from config object
  n_pop <- config$n_pop
  n_loci <- config$n_loci
  n_pairs <- n_pop / 2

  # Create sex vector (equal male/female split)
  sex <- rep(0:1, times = n_pop / 2)
  sibling_id <- rep(1:n_pairs, each = 2)

  # Create SNP matrix based on minor allele frequencies
  snps_paternal <- matrix(
    rbinom(
      n_pairs * n_loci,
      size = 2,
      prob = rep(config$snp_maf, each = n_pairs)
    ),
    nrow = n_pairs,
    ncol = n_loci,
    byrow = FALSE
  )

  snps_maternal <- matrix(
    rbinom(
      n_pairs * n_loci,
      size = 2,
      prob = rep(config$snp_maf, each = n_pairs)
    ),
    nrow = n_pairs,
    ncol = n_loci,
    byrow = FALSE
  )

  snp_matrix <- generate_snp_matrix(config, snps_paternal, snps_maternal)
  colnames(snp_matrix) <- paste0("rs", 1:n_loci)

  # Compute phenotype vector
  phenotype <- generate_phenotype(config, snp_matrix)

  # Combine into population data frame
  population <- as.data.frame(cbind(sex, snp_matrix, phenotype, sibling_id))

  return(population)
}
