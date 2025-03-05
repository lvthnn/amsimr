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
    gene_snp <- gene[, snp_pheno]
    gam_snp <- rnorm(length(snp_pheno))

    names(gam_snp) <- paste0("rs", snp_pheno)
    gam_causal[[pheno]] <- gam_snp

    phenos <- gene_snp %*% gam_snp + rnorm(n, sd = 0.5)
    colnames(phenos) <- pheno

    phi <- cbind(phi, phenos)
  }

  pop <- data.frame(sex, gene, phi, family_id = 1:n)

  attr(pop, "maf_genes") <- maf
  attr(pop, "gam_causal") <- gam_causal
  attr(pop, "psi_vec") <- psi_vec

  return(pop)
}

sim_matching <- function(pop, iter = 100000, alpha = 1, temp0 = 1e-9,
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
  males <- cbind(
    scale(pop[pop$sex == 0, unique(psi_vec[, 1]), drop = FALSE]),
    id = as.integer(rownames(pop[pop$sex == 0, ]))
  )
  females <- cbind(
    scale(pop[pop$sex == 1, unique(psi_vec[, 2]), drop = FALSE]),
    id = as.integer(rownames(pop[pop$sex == 1, ]))
  )
  rownames(males) <- rownames(females) <- 1:n
  colnames(females) <- paste0("f_", colnames(females))

  init_sol <- as.matrix(cbind(males, females))

  cf_idx <- which(grepl("f_", colnames(init_sol))) - 1

  psi_vec[, 1] <- as.integer(match(psi_vec[, 1], colnames(init_sol)) - 1)
  psi_vec[, 2] <- as.integer(match(paste0("f_", psi_vec[, 2]),
                                   colnames(init_sol)) - 1)
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

generate_offspring <- function(pop, female_idx) {
  cols <- c("sex", attr(pop, "maf_genes") |> names())
  females <- pop[female_idx, cols]
  males <- pop[setdiff(1:nrow(pop), female_idx), cols]
  num_pairs <- nrow(males)
  
  offspring_boys <- data.frame()
  offspring_girls <- data.frame()
  
  for (i in 1:num_pairs) {
    male <- males[i, ]
    female <- females[i, ]
    genetic_variants_male <- male[-1]
    genetic_variants_female <- female[-1]
    
    boy_genes <- sapply(1:length(genetic_variants_male), function(j) {
      binom_male <- rbinom(1, 1, 0.5 * genetic_variants_male[[j]])
      binom_female <- rbinom(1, 1, 0.5 * genetic_variants_female[[j]])
      binom_male + binom_female
    })
    
    girl_genes <- sapply(1:length(genetic_variants_female), function(j) {
      binom_male <- rbinom(1, 1, 0.5 * genetic_variants_male[[j]])
      binom_female <- rbinom(1, 1, 0.5 * genetic_variants_female[[j]])
      binom_male + binom_female
    })
    
    names(boy_genes) <- names(girl_genes) <- names(genetic_variants_male)
   
    gam_causal <- attr(pop, "gam_causal")
    phenos <- attr(pop, "gam_causal") |> names()
    
    boy_phenos <- sapply(names(gam_causal), function(pheno) {
      gam <- gam_causal[[pheno]]
      gam_snps <- boy_genes[names(gam)]
      boy_pheno <- gam_snps %*% gam + rnorm(1, sd = 0.5)
      return(boy_pheno)
    })
    
    girl_phenos <- sapply(names(gam_causal), function(pheno) {
      gam <- gam_causal[[pheno]]
      gam_snps <- girl_genes[names(gam)]
      girl_pheno <- gam_snps %*% gam + rnorm(1, sd = 0.5)
      return(girl_pheno)
    })
    
    names(boy_phenos) <- names(girl_phenos) <- names(gam_causal)
    
    boy <- c(sex = 0, boy_genes, boy_phenos, family_id = i)
    girl <- c(sex = 1, girl_genes, girl_phenos, family_id = i)
    
    offspring_boys <- rbind(offspring_boys, boy)
    offspring_girls <- rbind(offspring_girls, girl)
  }
  
  colnames(offspring_boys) <- colnames(offspring_girls) <- c(cols, "family_id")
  next_generation <- rbind(offspring_boys, offspring_girls)
  
  return(next_generation)
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
  plot(samps_eval$alpha, xlab = "Iteration", ylab = "Acceptance rate")
  plot(samps_eval$temp, type = "l", xlab = "Iteration", ylab = "Temperature")
}
