population <- function(sim_params) {
  #' Generate an initial population with initial values from config file
  #'
  #' @param sim_params Path to JSON config file with simulation parameters.

  # Retrieve simulation parameters from configuration file
  config <- jsonlite::fromJSON(sim_params)

  # Extract relevant parameters
  n <- config$n_pop
  p <- config$n_gene
  psi_vec <- config$psi_vec
  rownames(psi_vec) <- seq_len(nrow(psi_vec))
  colnames(psi_vec) <- c("phi_m", "phi_f", "cor")

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
    gam_snp <- numeric(p)
    names(gam_snp) <- paste0("rs", 1:p)
    gam_snp[snp_pheno] <- rnorm(length(snp_pheno)) # accept snp effects from
                                                   # configuration file instead
    gam_causal[[pheno]] <- gam_snp

    phenos <- gene %*% gam_snp + rnorm(n, sd = 0.25)
    colnames(phenos) <- pheno
    phi <- cbind(phi, phenos)
  }

  phi[sex == 0, ] <- scale(phi[sex == 0, ])
  phi[sex == 1, ] <- scale(phi[sex == 1, ])

  gam_causal <- t(do.call(rbind, gam_causal) |> as.matrix())

  pop <- data.frame(sex, gene, phi, f_id = 0, m_id = 0)

  attr(pop, "maf_genes") <- maf
  attr(pop, "gam_causal") <- gam_causal
  attr(pop, "psi_vec") <- psi_vec
  attr(pop, "past_gens") <- list()

  return(pop)
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

plot_eval <- function(samps_eval) {
  #' Plot the results from the return of a `sim_matching` function call
  #'
  #' @param samps_eval The return from `sim_matching` with `eval = TRUE`
  par(mfrow = c(2, 2))
  plot(
    samps_eval$psi,
    type = "l",
    xlab = "Iteration",
    ylab = expression(psi)
  )
  plot(
    samps_eval$dpsi,
    type = "l",
    xlab = "Iteration",
    ylab = expression(Delta * psi)
  )
  plot(
    samps_eval$alpha,
    xlab = "Iteration",
    ylab = "Acceptance rate"
  )
  plot(
    samps_eval$temp,
    type = "l",
    xlab = "Iteration",
    ylab = "Temperature"
  )
}
