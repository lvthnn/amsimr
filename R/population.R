#' Generate phenotype values from genotypes and environmental effects
#'
#' @param config Configuration list with simulation parameters
#' @param snp_matrix Matrix of SNP genotypes (0, 1, or 2 copies of alt allele)
#'
#' @return A matrix with one column containing the phenotype values
#'
#' @importFrom stats var rnorm
#' @NoRd
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
  colnames(phenotype) <- pheno_name

  return(phenotype)
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

  # Create sex vector (equal male/female split)
  sex <- rep(0:1, each = n_pop / 2)

  # Create SNP matrix based on minor allele frequencies
  snp_matrix <- matrix(
    rbinom(
      n_pop * n_loci,
      size = 2,
      prob = rep(config$snp_maf, each = n_pop)
    ),
    nrow = n_pop,
    ncol = n_loci,
    byrow = FALSE
  )
  colnames(snp_matrix) <- paste0("rs", 1:n_loci)

  # Compute phenotype vector
  phenotype <- generate_phenotype(config, snp_matrix)

  # Combine into population data frame
  population <- as.data.frame(cbind(sex, snp_matrix, phenotype))

  return(population)
}

#' Generate initial solution matrix for mate matching optimisation routine
#'
#' @param population Population data frame with sex column and SNP genotypes
#' @param config Configuration list with mating_model_pairs
#'
#' @return List with init_sol matrix, psi_vec, and female_swap_idx
#' @NoRd
generate_init_state <- function(config, population) {
  # Extract mating model pairs from flattened config
  mating_pairs <- config$mating_model_pairs
  male_snps <- unlist(unique(mating_pairs$male_snp))
  female_snps <- unlist(unique(mating_pairs$female_snp))

  males <- cbind(
    scale(population[population$sex == 0, male_snps, drop = FALSE]),
    id = as.integer(rownames(population[population$sex == 0, ]))
  )

  females <- cbind(
    scale(population[population$sex == 1, female_snps, drop = FALSE]),
    id = as.integer(rownames(population[population$sex == 1, ]))
  )
  females <- females[sample(nrow(females)), ]

  # Set consistent row names and prefix female columns
  rownames(males) <- rownames(females) <- 1:(config$n_pop / 2)
  colnames(females) <- paste0("f_", colnames(females))

  # Combine into initial solution matrix
  init_sol <- as.matrix(cbind(males, females))

  # Initialize psi_vec to match dimensions of the input
  snp_pairs <- matrix(0, nrow = nrow(mating_pairs), ncol = 3)
  male_cols <- match(mating_pairs[, "male_snp"], colnames(init_sol)) - 1
  female_cols <- match(paste0("f_", mating_pairs[, "female_snp"]),
                       colnames(init_sol)) - 1
  snp_pairs[, 1] <- male_cols
  snp_pairs[, 2] <- female_cols
  snp_pairs[, 3] <- as.numeric(mating_pairs$correlation)
  storage.mode(snp_pairs) <- "numeric"

  # Get indices of female columns for swapping
  female_swap_idx <- which(grepl("^f_", colnames(init_sol))) - 1

  # Return as a list
  return(list(
    init_sol = init_sol,
    snp_pairs = snp_pairs,
    female_swap_idx = female_swap_idx
  ))
}

#' Generate an approximately optimal mate matching
#'
#' @param config Configuration file with simulation parameters
#' @param population Population data frame
#' @param collect_metrics Whether to collect metrics from simulated annealing
#'  routine to evaluate performance
#'
#' @return List with init_sol matrix, psi_vec, and female_swap_idx
#' @NoRd
generate_matching <- function(config, population, collect_metrics = FALSE) {
  is_am <- config$mating_model_type == "assortative"
  sol_data <- generate_init_state(config, population)

  # Run the simulated annealing schedule to approximate the optimal matching
  if (is_am) {
    optim_sol <- optim_matching(
      sol_mat = sol_data$init_sol,
      snp_pairs = sol_data$snp_pairs,
      female_swap_idx = sol_data$female_swap_idx,
      num_iterations = 1e6,
      collect_metrics = FALSE
    )
  } else {
    optim_sol <- sol_data$init_sol
  }

  colnames(optim_sol$sol_mat) <- colnames(sol_data$init_sol)
  optim_sol$sol_mat <- as.data.frame(optim_sol$sol_mat)

  matching <- list()
  matching[["male_idx"]] <- optim_sol$sol_mat$id
  matching[["female_idx"]] <- optim_sol$sol_mat$f_id

  # Compute the attained correlation in the matching based on SNP pairs
  matching_data <- optim_sol$sol_mat

  if (is_am) {
    attained_cor <- sapply(
      seq_len(nrow(config$mating_model_pairs)),
      function(i) {
        male_snp <- config$mating_model_pairs[[i, 1]]
        female_snp <- paste0("f_", config$mating_model_pairs[[i, 2]])
        cor(matching_data[, male_snp], matching_data[, female_snp])
      }
    )
  }

  attr(matching, "attained_cor") <- attained_cor

  return(matching)
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
#' @NoRd
generate_offspring <- function(config, population, matching) {
  # Helper function to create offspring SNP matrix
  create_snp_matrix <- function(snps_male, snps_female) {
    n_pop <- config$n_pop
    n_loci <- config$n_loci

    # Calculate probability matrices for inheritance
    paternal_probs <- as.matrix(snps_male) / 2
    maternal_probs <- as.matrix(snps_female) / 2

    # Sample alleles from each parent
    paternal_alleles <- matrix(
      rbinom(
        n_pop * n_loci,
        size = 1,
        prob = rep(as.vector(paternal_probs), each = 2)
      ),
      nrow = n_pop,
      ncol = n_loci
    )

    maternal_alleles <- matrix(
      rbinom(
        n_pop * n_loci,
        size = 1,
        prob = rep(as.vector(maternal_probs), each = 2)
      ),
      nrow = n_pop,
      ncol = n_loci
    )

    # Combine alleles from both parents
    snps_offspring <- paternal_alleles + maternal_alleles
    colnames(snps_offspring) <- paste0("rs", 1:n_loci)

    return(snps_offspring)
  }

  # Create offspring sex vector (equal male/female split)
  sex <- rep(0:1, times = config$n_pop / 2)

  # Select genotype data in male and female subpopulations
  cols_select <- !(colnames(population) %in% c("sex", config$phenotype_name))
  snps_male <- population[matching$male_idx, cols_select]
  snps_female <- population[matching$female_idx, cols_select]

  # Generate SNP matrix and phenotype
  snp_matrix <- create_snp_matrix(snps_male, snps_female)
  phenotype <- generate_phenotype(config, snp_matrix)

  # Assemble offspring population
  offspring_pop <- data.frame(sex, snp_matrix, phenotype)

  return(offspring_pop)
}

#' Calculate the genetic correlation between siblings in a population
#'
#' @param config Configuration list with simulation parameters
#' @param population Data frame containing the population with siblings in
#'  consecutive pairs
#' @param phenotype_name Name of the phenotype column (optional)
#'
#' @return A list containing the genetic correlation and additional statistics
#'
#' @importFrom stats cor sd
#' @NoRd
compute_sib_cor <- function(config, population, phenotype_name) {
  # Check if population size is even (required for pairing)
  if (nrow(population) %% 2 != 0) {
    stop("Population size must be even to form sibling pairs")
  }

  # Extract genotype columns
  genotype_cols <- paste0("rs", 1:config$n_loci)
  genotypes <- population[, genotype_cols, drop = FALSE]

  # Create separate matrices for each sibling in the pairs
  sibling1_genotypes <- genotypes[seq(1, config$n_pop, by = 2), , drop = FALSE]
  sibling2_genotypes <- genotypes[seq(2, config$n_pop, by = 2), , drop = FALSE]

  # Calculate per-locus correlation
  locus_correlations <- numeric(config$n_loci)
  for (i in 1:config$n_loci) {
    locus_correlations[i] <- cor(
      sibling1_genotypes[, i],
      sibling2_genotypes[, i]
    )
  }
  names(locus_correlations) <- genotype_cols

  # Calculate average genetic correlation across all loci
  mean_genetic_correlation <- mean(locus_correlations, na.rm = TRUE)

  # Calculate phenotypic correlation if phenotype name is provided
  pheno_correlation <- NULL
  if (!missing(phenotype_name) && phenotype_name %in% colnames(population)) {
    phenos <- population[[phenotype_name]]
    sibling1_pheno <- phenos[seq(1, nrow(population), 2)]
    sibling2_pheno <- phenos[seq(2, nrow(population), 2)]
    pheno_correlation <- cor(sibling1_pheno, sibling2_pheno)
  }

  # Calculate estimated heritability if phenotypic correlation is available
  heritability_estimate <- NULL
  if (!is.null(pheno_correlation)) {
    heritability_estimate <- pheno_correlation / mean_genetic_correlation
  }

  # Calculate additional genotype statistics
  mean_genotype_similarity <- mean(
    colMeans(sibling1_genotypes == sibling2_genotypes)
  )

  # Return results
  results <- list(
    mean_genetic_correlation = mean_genetic_correlation,
    locus_correlations = locus_correlations,
    locus_correlation_sd = sd(locus_correlations, na.rm = TRUE),
    mean_genotype_similarity = mean_genotype_similarity
  )

  # Add phenotype results if available
  if (!is.null(pheno_correlation)) {
    results[paste0(config$phenotype_name, "_correlation")] <- pheno_correlation
    results$heritability_estimate <- heritability_estimate
  }

  return(results)
}
