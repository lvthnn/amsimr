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
  maf <- rbeta(p, 5, 5)
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
    gam_snp <- rnorm(length(snp_pheno) + 1)
    phenos <- gam_snp[1] + gene_snp %*% gam_snp[-1] + rnorm(n)
    names(gam_snp) <- c("gam_0", paste0("rs", snp_pheno))
    colnames(phenos) <- pheno

    gam_causal[[pheno]] <- gam_snp
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
      pop_obj <- ifelse(sex == 0, pop$get_sex(0), pop$get_sex(1))
    }
    print(pop_obj)
    return(pop_obj[i, rs])
  }
  
  pop$get_pheno <- function(i, phe = paste("X", 1:q), sex = NA) {
    pop_obj <- pop$data
    if(!missing(sex)) {
      pop_obj <- ifelse(sex == 0, pop$get_sex(0), pop$get_sex(1))
    }
    return(pop_obj[i, phe])
  }
  
  class(pop) <- "population"

  return(pop)
}

alpha <- function(pop, i, j) {
  #' Potential function for propensity distribution on space of potential
  #' spouse pairs.
  #'
  #' @param pop Population object
  #' @param i Male individual with index no. i
  #' @param j Female individual with index no. j
  rules_propensity <- pop$rules_propensity
  kappa <- pop$kappa
  alpha <- 0
  
  for (k in 1:nrow(rules_propensity)) {
    phi1 <- rules_propensity[k, 1]
    phi2 <- rules_propensity[k, 2]
    w <- as.numeric(rules_propensity[k, 3])
    
    phi1_i <- pop$get_sex(0)$data[i, phi1]
    phi2_j <- pop$get_sex(1)[j, phi2]
  
    alpha <- alpha + w * (phi1_i - phi2_j)^2
  }

  return(alpha)
}

propensity_dist <- function(pop, n_samples, n_burnin) {
  pop_male <- pop$get_sex(0)
  pop_female <- pop$get_sex(1)
  n_male <- nrow(pop_male)
  n_female <- nrow(pop_female)
  
  kappa <- pop$kappa
  
  rand_state <- function() c(sample(1:n_male, 1), sample(1:n_female, 1))
  samples <- list()
  samples[[1]] <- rand_state()
  
  for (s in 2:n_samples) {
    pair_curr <- samples[[s - 1]]
    pair_prop <- rand_state()
    
    alpha_curr <- alpha(pop, pair_curr[1], pair_curr[2])
    alpha_prop <- alpha(pop, pair_prop[1], pair_prop[2])
    
    acc_ratio <- min(exp((alpha_prop - alpha_curr) / kappa), 1)
    
    if (runif(1) < acc_ratio) {
      samples[[s]] <- pair_prop
    } else {
      samples[[s]] <- pair_curr
    }
  }
  
  return(samples)
}

pop <- population(params = "Github/AssocMating/config.json")

pair_samples <- do.call(rbind, propensity_dist(pop, 1e5, 0)) |> as.data.frame()
colnames(pair_samples) <- c("i", "j")
pair_samples <- pair_samples |>
  dplyr::mutate(pairs = paste0("(", i, ",", j, ")")) |>
  dplyr::arrange(i, j)

table(pair_samples$pairs) |> prop.table() |> hist()

mating_pairs <- c()

while (length(mating_pairs) < min(pop$get_sex(0) |> nrow(), pop$get_sex(1) |> nrow())) {
  ind <- sample(1:nrow(pair_samples), size = 1)
  i_ind <- pair_samples[ind, 1]
  j_ind <- pair_samples[ind, 2]

  pair_samples <- pair_samples |> dplyr::filter(i != i_ind, j != j_ind)  
  mating_pairs <- append(mating_pairs, paste0("(", i_ind, ",", j_ind, ")")) 
}
