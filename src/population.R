population <- function(params) {
  #' Generate an initial population
  #'
  #' Generate the 0th generation of a population.
  #' @param params Path to JSON config file with simulation parameters.

  # Retrieve simulation parameters from configuration file
  config <- jsonlite::fromJSON(params)
 
  # Extract relevant parameters 
  n <- config$n_pop
  p <- config$n_gene
  q <- config$n_pheno

  # Simulate sex of individuals, minor allele frequencies and genotypes
  sex <- sample(rep(c(0, 1), n / 2))
  maf <- 0.5 # @TODO: Make into beta / uniform variable
  gene <- matrix(rbinom(n * p, 2, rep(maf, each = n)), nrow = n, ncol = p)

  names(maf) <- paste0("rs", 1:p)
  colnames(gene) <- paste0("rs", 1:p)

  # Generate phenotypes, for now model as function of genotype
  snps_causal <- config$snp_causal
  phi <- matrix(nrow = n, ncol = 0)
  gam_causal <- list()
  
  for (pheno in names(snps_causal)) {
    snp_pheno <- snps_causal[[pheno]]
    gene_snp <- gene[, snp_pheno]
    gam_snp <- matrix(1)# @TODO: Make into variable
    
    names(gam_snp) <- paste0("rs", snp_pheno)
    gam_causal[[pheno]] <- gam_snp
    
    phenos <- gene_snp %*% gam_snp + rnorm(n, sd = 0.5)
    colnames(phenos) <- pheno
    
    phi <- cbind(phi, phenos)
  }
  
  pop <- data.frame(sex, gene, phi)
  sex <- pop[, "sex"] == 0
  
  # Assemble the population object
  pop <- list(
    pop = pop,
    sex = sex,
    maf_genes = maf,
    gam_causal = gam_causal,
    data = data.frame(sex, gene, phi)
  )
  
  class(pop) <- "population"
  return(pop)
}

delta_psi <- function(pop, sol, swap) {
  i <- swap[1]; j <- swap[2]
  ai <- sol[i]; aj <- sol[j]
  
  # Extract phenotype values of the original and proposed pairs
  phi_m <- pop$data[pop$sex, "phi1"][swap]
  phi_f <- pop$data[!pop$sex, "phi1"][sol[swap]]
  
  # Compute squared difference in phenotype
  return(2 * (phi_m[1] - phi_m[2]) * (phi_f[1] - phi_f[2]))
}

psi <- function(pop, sol) {
  phi_m <- pop$pop[pop$sex, "phi1"]
  phi_f <- pop$pop[!pop$sex, "phi1"]
  
  score <- 0
  for (j in 1:length(sol)) {
    score <- score + sum((phi_m[j] - phi_f[sol[j]])^2)
  }
  
  return(score)
}

sim_matching <- function(pop, iter = 500000, kappa0 = 100, eval = FALSE, lag = 1000) {
  #' Generate a partner matching for a generation of individuals.
  #'
  #' @param pop Population object
  #' @param iter Number of iterations before matching is returned
  #' @param kappa Initial value of inverse temperature parameter
  #' @param eval Return simulated values
  #' @param lag Lag interval with which to store samples if eval = TRUE
 
  # Size of male and female subpopulations
  n <- nrow(pop$data) / 2
  index <- 1
  
  # If eval = TRUE, store simulation values
  if (eval) {
    acc_prob <- numeric(iter / lag)
    psi_sols <- numeric(iter / lag)
    temp <- numeric(iter / lag)
  }
 
  # Initial sample 
  sol <- sample(1:n)
  names(sol) <- 1:n
  
  swaps <- replicate(iter, sample(1:n, size = 2), simplify = FALSE)
  phi_m <- pop$pop[pop$sex, "phi1"]
  phi_f <- pop$pop[!pop$sex, "phi1"]
  
  kappa <- kappa0 / log(1 + 100 * 1:iter)
  
  # Metropolis step
  for (s in 2:iter) {
    if (s %% 1000 == 0) print(s)
    swap <- swaps[[s]]
    
    delta <- 2 * (phi_m[swap[1]] - phi_m[swap[2]]) * (phi_f[sol[swap[1]]] - phi_f[sol[swap[2]]])
    alpha_acc <- min(1, exp(-delta / kappa[s]))

    if (runif(1) < alpha_acc) {
      sol[swap] <- sol[rev(swap)]
    }
    
    if(eval && s %% 1000 == 0) {
      acc_prob[i] <- alpha_acc
      psi_sols[i] <- psi(pop, sol) 
      temp[i] <- kappa[i]
      i <- i + 1
    }
  }
  
  if (eval) {
    return(list(
      sol = sol,
      psi_sols = psi_sols,
      acc_prob = acc_prob,
      temp = temp 
    ))
  } else {
    return(sol)
  }
}

pop <- population(params = "Github/AssocMating/config.json")

test_matching <- sim_matching(pop, iter = 1500000, kappa = 100, eval = TRUE)

par(mfrow = c(1, 3))
plot(test_matching$acc_prob)
plot(test_matching$temp, type = 'l')
plot(test_matching$psi_sol, type = 'l')
 
matching <- data.frame(names(test_matching$sol) |> as.integer(), test_matching$sol,
                       pop$pop[pop$sex, "phi1"], pop$pop[!pop$sex, "phi1"])
colnames(matching) <- c("m", "f", "phi_m", "phi_f")

plot(matching$phi_m, matching$phi_f)
