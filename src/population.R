population <- function(sim_params) {
  #' Generate an initial population with initial values from config file
  #'
  #' @param sim_params Path to JSON config file with simulation parameters.

  # Retrieve simulation parameters from configuration file
  config <- jsonlite::fromJSON(sim_params)

  # Extract relevant parameters
  n <- config$n_pop
  p <- config$n_gene
  psi_vec <- as.data.frame(config$psi_vec)

  rownames(psi_vec) <- seq_len(nrow(psi_vec))
  colnames(psi_vec) <- c("phi_m", "phi_f", "w")
  psi_vec$phi_f <- paste0("f_", psi_vec$phi_f)

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

swap_sol <- function(sol, swap, exch_col) {
  #' Auxiliary function to perform transposition
  #'
  #' @param sol The current state
  #' @param swap Vector of 2-swap indices to form proposal state
  #' @param exch_col Columns that specify female mate to swap

  sol[swap, exch_col] <- sol[rev(swap), exch_col]

  return(sol)
}

psi <- function(sol, psi_vec) {
  #' Avoidance function for a given mating permutation
  #'
  #' @param sol The state for which the density is computed
  #' @param phis Phenotype pairs to be modelled in the avoidance function

  delta <- sol[, psi_vec$phi_m] - sol[, psi_vec$phi_f]
  w <- as.integer(psi_vec$w)
  score <- sum(colSums(delta^2) * w)

  return(score)
}

dpsi <- function(sol, swap, psi_vec, exch_col) {
  #' Compute the ratio of densities of current and proposal state.
  #'
  #' @param sol The current state
  #' @param swap Vector of 2-swap indices to form proposal state
  #' @param exch_col Columns that specify female mate to swap

  res <- psi(swap_sol(sol, swap, exch_col), psi_vec) - psi(sol, psi_vec)

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

  # Initial sample -- change this to greedy
  males <- pop[pop$sex == 0, ]
  females <- pop[pop$sex == 1, ]
  rownames(males) <- rownames(females) <- 1:n
  names(females) <- paste0("f_", colnames(females))

  # Set up the first sample
  sol <- sol_opt <- as.matrix(cbind(males, females))
  psi_opt <- psi(sol, psi_vec)
  exch_col <- grepl("f_", names(sol))

  # Samples for exchanges
  swaps <- replicate(iter, sample(n, size = 2), simplify = FALSE)
  temp <- temp0

  # If eval is TRUE, set up vectors to store data in
  if (eval) {
    psi_eval <- numeric(iter)
    dpsi_eval <- numeric(iter)
    rho_eval <- numeric(iter)
    temp_eval <- numeric(iter)
    psi_eval[1] <- psi_opt
    dpsi_eval[1] <- 0
    rho_eval[1] <- 1
    temp_eval[1] <- temp0
  }

  # Metropolis step
  for (s in 2:iter) {
    swap <- swaps[[s]]
    psi_sol <- psi(sol, psi_vec)
    dpsi_sol <- dpsi(sol, swap, psi_vec, exch_col)
    rho_s <- min(1, exp(-dpsi_sol / temp))

    if (progress && s %% 1000 == 0) print(s)

    # Accept step
    if (runif(1) < rho_s) {
      sol <- swap_sol(sol, swap, exch_col)
      if (psi_sol < psi_opt) sol_opt <- sol
    }

    if (eval) {
      psi_eval[s] <- psi(sol, psi_vec)
      dpsi_eval[s] <- dpsi(sol, swap, psi_vec, exch_col)
      rho_eval[s] <- rho_s
    }

    temp <- alpha * temp
    temp_eval[s] <- temp
  }

  if (eval) {
    return(list(
      sol = sol,
      sol_opt = sol_opt,
      psi_opt = psi(sol_opt, psi_vec),
      dpsi = dpsi_eval,
      psi = psi_eval,
      rho = rho_eval,
      temp = temp_eval
    ))
  } else {
    return(list(
      sol = sol
    ))
  }
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

