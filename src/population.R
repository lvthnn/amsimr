population <- function(sim_params) {
  #' Generate an initial population with initial values from config file
  #'
  #' @param sim_params Path to JSON config file with simulation parameters.

  # Retrieve simulation parameters from configuration file
  config <- jsonlite::fromJSON(sim_params)

  # Extract relevant parameters
  n <- config$n_pop
  p <- config$n_gene
  psi_vec <- as.matrix(config$psi_vec)

  rownames(psi_vec) <- seq_len(nrow(psi_vec))
  colnames(psi_vec) <- c("phi_m", "phi_f", "w")

  # Simulate sex of individuals, minor allele frequencies and genotypes
  sex <- sample(rep(c(0, 1), n / 2))
  maf <- rbeta(p, 5, 5)
  gene <- matrix(rbinom(n * p, 2, rep(maf, each = n)), nrow = n, ncol = p)
  names(maf) <- colnames(gene) <- paste0("rs", 1:p)

  # Generate phenotypes, for now model as function of genotype
  snps_causal <- config$snp_causal
  phi <- matrix(nrow = n, ncol = 0)
  gam_causal <- list()

  for (pheno in names(snps_causal)) {
    snp_pheno <- snps_causal[[pheno]]
    gene_snp <- gene[, snp_pheno]
    gam_snp <- rnorm(length(snp_pheno))

    names(gam_snp) <- paste0("rs", snp_pheno)
    gam_causal[[pheno]] <- gam_snp

    phenos <- gene_snp %*% gam_snp + rnorm(n) # Phenotype generation
    colnames(phenos) <- pheno

    phi <- cbind(phi, phenos)
  }

  pop <- data.frame(sex, gene, phi)

  attr(pop, "maf_genes") <- maf
  attr(pop, "gam_causal") <- gam_causal
  attr(pop, "psi_vec") <- psi_vec

  return(pop)
}

swap_sol_old <- function(sol, swap, exch_col) {
  #' Auxiliary function to perform transposition
  #'
  #' @param sol The current state
  #' @param swap Vector of 2-swap indices to form proposal state
  #' @param exch_col Columns that specify female mate to swap

  sol[swap, exch_col] <- sol[rev(swap), exch_col]

  return(sol)
}

psi_old <- function(sol, psi_vec) {
  #' Avoidance function for a given mating permutation
  #'
  #' @param sol The state for which the density is computed
  #' @param phis Phenotype pairs to be modelled in the avoidance function

  delta <- sol[, phi_m + 1] - sol[, phi_f + 1]
  w <- matrix(as.numeric(psi_vec$w), nrow = 1)
  score <- colSums(delta^2 %*% t(w))

  return(score)
}

dpsi_old <- function(sol, swap, psi_vec, exch_col) {
  #' Compute the ratio of densities of current and proposal state.
  #'
  #' @param sol The current state
  #' @param swap Vector of 2-swap indices to form proposal state
  #' @param exch_col Columns that specify female mate to swap

  res <- psi_old(swap_sol_old(sol, swap, exch_col), psi_vec) - psi_old(sol, psi_vec)

  return(res)
}

sim_matching <- function(pop, iter = 10000, alpha = 0.9995, temp0 = 5,
                         eval = FALSE, progress = FALSE) {
  #' Find an optimal matching between individuals in the population
  #'
  #' @param pop Population data frame
  #' @param iter Number of iterations simulated by Metropolis chain
  #' @param alpha Exponential decay coefficient of inverse temperature
  #' @param temp0 Initial temperature of simulated annealing procedure
  #' @param eval If TRUE, returns simulation data

  # Size of male and female populations
  n <- nrow(pop) / 2
  psi_vec <- attr(pop, "psi_vec")
  
  # Initial solution
  males <- pop[pop$sex == 0, ]
  females <- pop[pop$sex == 1, ]
  males$id <- females$id <- 1:n
  
  colnames(females) <- paste0("f_", colnames(females))
  init_sol <- as.matrix(cbind(males, females))
  
  pheno_names <- names(attr(pop, "gam_causal"))
  cf_idx <- which(grepl("f_", colnames(init_sol))) - 1
  
  psi_vec[, 1] <- as.integer(match(psi_vec[, 1], colnames(init_sol)) - 1)
  psi_vec[, 2] <- as.integer(match(paste0("f_", psi_vec[, 2]), colnames(init_sol)) - 1)
  psi_vec[, 3] <- as.numeric(psi_vec[, 3])
  storage.mode(psi_vec) <- "numeric"
  
  # Run simulated annealing algorithm
  res <- optim_matching(init_sol, psi_vec, cf_idx = cf_idx, n_iter = iter,
                        alpha = alpha, temp0 = temp0, eval = eval, progress = progress)
  
  res$sol <- as.data.frame(res$sol)
  colnames(res$sol) <- colnames(init_sol)
  
  return(res)
}

plot_eval <- function(samps_eval) {
  #' Plot the results from the return of a `sim_matching` function call
  #'
  #' @param samps_eval The return from `sim_matching` with `eval = TRUE`

  par(mfrow = c(2, 2))
  plot(samps_eval$psi, type = "l", xlab = "Iteration", ylab = expression(psi))
  plot(samps_eval$dpsi,
    type = "l", xlab = "Iteration",
    ylab = expression(Delta * psi)
  )
  plot(samps_eval$rho, xlab = "Iteration", ylab = "Acceptance rate")
  plot(samps_eval$temp, type = "l", xlab = "Iteration", ylab = "Temperature")
}

