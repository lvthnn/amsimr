initialise_population <- function(config) {
  #' Generate an initial population with initial values from config file
  #'
  #' @param sim_params Path to config file with simulation parameters
  #' @returns A data frame containing the initial population

  # Extract data from the config object
  n_pop <- config$n_pop
  n_genes <- config$n_genes
  snp_maf <- config$snp_maf

  phenotype_causal_snps <- config$phenotype_causal_snps

  # Generate a vector of individual sexes
  pop_sex <- rep(0:1, each = config$n_pop / 2)

  # Create allele count matrix for the population
  gene_matrix <- matrix(rbinom(n_pop * n_genes, size = 2,
                               prob = rep(snp_maf, each = n_pop)),
                        nrow = n_pop, ncol = n_genes, byrow = FALSE)

  # Compute the genetic contribution, tune heritability with random
  # normal variables

  phenotype_raw <- gene_matrix %*% phenotype_causal_snps
  phenotype_heritability <- config$phenotype_heritability
  alpha <- (1 - phenotype_heritability^2) / phenotype_heritability^2
}

sim_matching <- function(pop, iter = 100000, alpha = 0.9995, temp0 = 1e-9,
                         eval = FALSE) {
  #' Find an optimal matching between individuals in the population
  #'
  #' @param pop Population data frame
  #' @param iter Number of iterations simulated by Metropolis chain
  #' @param alpha Exponential decay coefficient of inverse temperature
  #' @param temp0 Initial temperature of simulated annealing procedure
  #' @param eval If TRUE, returns simulation data

  zero_idx <- function(data, cols, prefix = "") {
    if (prefix != "") cols <- paste0(prefix, cols)
    match(cols, colnames(data))
  }

  # Size of male and female populations
  n <- nrow(pop) / 2
  psi_vec <- attr(pop, "psi_vec")

  # Initial solution
  snp_cols <- grepl("sex|rs", colnames(pop))

  stdrs_male <- scale(pop[pop$sex == 0, unique(psi_vec[, 1]), drop = FALSE])
  stdrs_female <- scale(pop[pop$sex == 1, unique(psi_vec[, 2]), drop = FALSE])
  colnames(stdrs_male) <- paste0("std_", colnames(stdrs_male))
  colnames(stdrs_female) <- paste0("std_", colnames(stdrs_female))

  psi_vec[, 1] <- paste0("std_", psi_vec[, 1])
  psi_vec[, 2] <- paste0("std_", psi_vec[, 2])

  males <- cbind(
    m_id = as.integer(rownames(pop[pop$sex == 0, ])),
    pop[pop$sex == 0, snp_cols],
    stdrs_male
  )

  females <- cbind(
    id = as.integer(rownames(pop[pop$sex == 1, ])),
    pop[pop$sex == 1, snp_cols],
    stdrs_female
  )

  rownames(males) <- rownames(females) <- 1:n
  colnames(females) <- paste0("f_", colnames(females))

  init_sol <- as.matrix(cbind(males, females))

  cf_idx <- which(grepl("f_", colnames(init_sol))) - 1

  psi_vec <- as.matrix(psi_vec)
  psi_vec[, 1] <- as.integer(match(
    psi_vec[, 1],
    colnames(init_sol)
  )) - 1
  psi_vec[, 2] <- as.integer(match(
    paste0("f_", psi_vec[, 2]),
    colnames(init_sol)
  )) - 1
  psi_vec[, 3] <- as.numeric(psi_vec[, 3])
  storage.mode(psi_vec) <- "numeric"

  # Run simulated annealing algorithm
  res <- optim_matching(
    init_sol,
    psi_vec,
    cf_idx = cf_idx,
    n_iter = iter,
    alpha = alpha,
    temp0 = temp0,
    eval = eval
  )

  res$sol <- as.data.frame(res$sol)
  colnames(res$sol) <- colnames(init_sol)

  return(res)
}

generate_offspring <- function(pop, sol) {
  #' Produce a generation of offspring designated by the mate matching sol
  #'
  #' @param pop The population data frame
  #' @param sol The mate matching

  # Retrieve some variables used in this function
  n <- nrow(pop)
  p <- nrow(attr(pop, "gam_causal"))
  gam_causal <- attr(pop, "gam_causal")
  cols <- which(grepl("rs", colnames(sol)))

  pairs_geno <- split.default(
    sol[, cols],
    cut(seq_along(sol[, cols]), 2,
      labels = FALSE
    )
  )

  # Generate the SNP matrix for offspring
  snp_probs <- as.matrix(0.25 * pairs_geno[["1"]] + 0.25 * pairs_geno[["2"]])

  snps <- matrix(
    rbinom(
      n * p,
      size = 2,
      prob = rep(as.vector(snp_probs), 2)
    ),
    nrow = n,
    ncol = p
  )

  colnames(snps) <- colnames(pairs_geno[["1"]])

  # Construct parent ID columns for pedigree distance
  parents <- sol[, c("m_id", "f_id")][rep(1:(n / 2), 2), ]
  rownames(parents) <- 1:n

  # Generate phenotypes for offspring
  phenos <- snps %*% gam_causal # + matrix(rnorm(n * 2, sd = 0.5), n, 2)

  sex <- rep(c(0, 1), each = n / 2)

  offspring <- cbind(sex, snps, scale(phenos), parents) |> as.data.frame()

  return(offspring)
}
