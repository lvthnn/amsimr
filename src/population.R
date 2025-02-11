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
  kappa <- config$kappa
  rules_propensity <- config$rules_propensity

  # Simulate sex of individuals, minor allele frequencies and genotypes
  sex <- sample(c(rep(0, n / 2), rep(1, n / 2)))
  maf <- runif(p, 0.1, 0.9)
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
    # gam_snp <- rnorm(length(snp_pheno), sd = 2)
    gam_snp <- c(2, 2, 2)
    
    names(gam_snp) <- paste0("rs", snp_pheno)
    gam_causal[[pheno]] <- gam_snp
    
    phenos <- gene_snp %*% gam_snp + rnorm(n, sd = 0.5)
    colnames(phenos) <- pheno
    
    phi <- cbind(phi, phenos)
  }
  
  # Assemble all of the above into the 0th generation
  pop <- list(
    kappa = kappa,
    maf_genes = maf,
    gam_causal = gam_causal,
    rules_propensity = rules_propensity,
    data = data.frame(sex, gene, phi)
  )
  
  pop_male <- pop$data[pop$data[, "sex"] == 0, ]
  pop_female <- pop$data[pop$data[, "sex"] == 1, ]
  
  # Get male/female individuals within the population
  pop$get_sex <- function(sex) {
    if (sex == 0) return(pop_male)
    else return(pop_female)
  }
  
  pop$get_geno <- function(i, rs = paste0("rs", 1:p), sex = NA) {
    pop_obj <- pop$data
    if(!missing(sex)) {
      if (sex == 0) pop_obj <- pop_male
      else pop_obj <- pop_female
    }
    return(pop_obj[i, rs])
  }
  
  pop$get_pheno <- function(i, phe = paste0("phi", 1:q), sex = NA) {
    pop_obj <- pop$data
    if(!missing(sex)) {
      if (sex == 0) pop_obj <- pop_male
      else pop_obj <- pop_female
    }
    return(pop_obj[i, phe])
  }
  
  class(pop) <- "population"

  return(pop)
}

log_propen <- function(pop, m, f) {
  #' Propensity kernel for propensity distribution on spouse pairs.
  #'
  #' @param pop Population object
  #' @param m Male individual with index no. m
  #' @param f Female individual with index no. f
  
  rules_propensity <- pop$rules_propensity
  kappa <- pop$kappa
  score <- 0
  
  for (k in 1:nrow(rules_propensity)) {
    phe1 <- rules_propensity[k, 1]
    phe2 <- rules_propensity[k, 2]
    w <- as.numeric(rules_propensity[k, 3])
    
    phi_m <- pop$get_pheno(m, phe = phe1, sex = 0)
    phi_f <- pop$get_pheno(f, phe = phe2, sex = 1)
  
    score <- score + w * (phi_m - phi_f)^2 / kappa
  }

  return(score)
}

spouse_pairs <- function(pop, n_iter = 1e5) {
  #' Generate spouse pairs using MCMC sampling.
  #'
  #' @param pop Population object
  #' @param n_iter Number of iterations before matching is returned
  
  pop_male <- pop$get_sex(0)
  pop_female <- pop$get_sex(1)
  n_male <- nrow(pop_male)
  n_female <- nrow(pop_female)
  
  curr <- sample(1:n_male)
  names(curr) <- 1:n_male

  for (s in 2:n_iter) {
    swap <- sample(curr, size = 2)
    ms <- names(curr[swap]) |> as.integer()
    fs <- curr[swap] |> as.integer()
   
    numer <- log_propen(pop, ms[1], fs[2]) + log_propen(pop, ms[2], fs[1])
    denom <- log_propen(pop, ms[1], fs[1]) + log_propen(pop, ms[2], fs[2])
    alpha_acc <- min(1, exp(numer - denom))
    
    if (runif(1) < alpha_acc) {
      curr[swap] <- curr[rev(swap)]
    }
  }
  
  return(curr)
}

pop <- population(params = "Github/AssocMating/config.json")
gen1 <- spouse_pairs(pop, n_iter = 10000)

gen1_df <- data.frame(names(gen1) |> as.integer(), gen1,
                      pop$get_sex(0)$phi1, pop$get_sex(1)$phi1)
colnames(gen1_df) <- c("m", "f", "phi_m", "phi_f")

plot(gen1_df$phi_m, gen1_df$phi_f)
