# Compute allele frequencies for a generation
get_snp_freq <- function(config, population) {
  rs_cols <- paste0("rs", 1:config$n_loci)
  rs_freqs <- colMeans(population[, rs_cols] / 2)
  return(rs_freqs)
}

# Compute genotype frequencies for a generation
get_genotype_freq <- function(config, population) {
  n_pop <- config$n_pop
  rs_cols <- paste0("rs", 1:config$n_loci)
  genotype_freq <- t(apply(population[, rs_cols], 2, function(x) {
    table(factor(x, levels = 0:2)) / n_pop
  }))
  return(genotype_freq)
}

# Compute allele correlation matrix for a generation
get_snp_cor <- function(config, population) {
  rs_cols <- paste0("rs", 1:config$n_loci)
  snp_matrix <- population[, rs_cols]
  return(cor(snp_matrix))
}

# Genetic correlation between siblings
get_sibling_genotype_cor <- function(config, population) {
  n_pop <- config$n_pop
  population_sibs <- population[order(population$sibling_id), ]
  snp_cols <- paste0("rs", 1:config$n_loci)

  sibs_male <- population_sibs[seq(1, n_pop, by = 2), snp_cols]
  sibs_female <- population_sibs[seq(2, n_pop, by = 2), snp_cols]

  # Compute the genetic correlation between siblings
  return(diag(cor(sibs_male, sibs_female, use = "pairwise.complete.obs")))
}

# Phenotype correlation between siblings
get_sibling_phenotype_cor <- function(config, population) {
  n_pop <- config$n_pop
  phenotype_col <- config$phenotype_name
  population_sibs <- population[order(population$sibling_id), ]
  sibs_male <- population_sibs[seq(1, n_pop, by = 2), phenotype_col]
  sibs_female <- population_sibs[seq(2, n_pop, by = 2), phenotype_col]

  # Compute the phenotype correlation between siblings
  return(cor(sibs_male, sibs_female))
}

# Heritability estimate of the phenotype
get_phenotype_heritability <- function(config, population) {
  phenotype_name <- config$phenotype_name
  phenotype <- population[, phenotype_name]
  phenotype_raw <- population[, paste0(phenotype_name, "_raw")]
  heritability <- var(phenotype_raw) / var(phenotype)

  return(heritability)
}

# ' Compute summary statistics for a simulation object
#'
#' @param generations A list of generations, each containing a population.
#' @param summarisers A named list of functions to compute summary statistics.
#'   Each function should take a population data frame as input and return a
#'   summary statistic.
#'
#' @return A list of summary statistics for each generation.
#'
#' @noRd
compute_summary_data <- function(config, population, summarisers = list()) {
  summary_results <- list()
  if (length(summarisers) == 0) return(list())

  for (summariser in names(summarisers)) {
    result <- summarisers[[summariser]](config, population)
    summary_results[[summariser]] <- result
  }

  return(summary_results)
}

#' Transform data compiled by summarisers for an amsim object.
#'
#' @param summary_data The amsim object
#'
#' @return Transformed summary data
transform_summary_data <- function(summary_data) {
  summarisers <- names(summary_data[[1]])
  transformed_data <- list()
  for (summariser in summarisers) {
    transformed_data[[summariser]] <- lapply(summary_data, function(x) {
      x[[summariser]]
    })
  }
  return(transformed_data)
}

#' Simulate one generation of the population
#'
#' Simulates one generation by performing mate matching, generating offspring,
#' and computing genetic similarity between siblings.
#'
#' @param config A list containing simulation parameters.
#' @param population A data frame representing the current population.
#'
#' @return A list containing the offspring population, SNP pair correlations
#'   from mate matching, and sibling genetic similarity.
#'
#' @noRd
simulate_generation <- function(config, population) {
  # Find mate matching and generate offspring population
  mate_matching <- generate_matching(config, population)
  population_offspring <- generate_offspring(config, population, mate_matching)

  # Attach observed SNP pair correlations to the result
  snp_pair_cors <- attr(mate_matching, "snp_pair_cors")
  attr(population_offspring, "snp_pair_cors") <- snp_pair_cors

  return(population_offspring)
}

#' Simulate a population from a config file
#'
#' @param config_path Path to configuration file with simulation parameters
#'
#' @return An amsim object containing the results of the simulation
#'
#' @export
simulate_population <- function(config_path, init_population = NULL) {
  config <- load_config(config_path)
  if (is.null(init_population)) init_population <- initialise_population(config)
  population <- init_population

  # To contain SNP correlation from annealing routine
  snp_pair_cors <- list()
  summarisers <- list(
    "snp_frequencies" = get_snp_freq,
    "snp_cor" = get_snp_cor,
    "genotype_frequencies" = get_genotype_freq,
    "sibling_genotype_cor" = get_sibling_genotype_cor,
    "sibling_phenotype_cor" = get_sibling_phenotype_cor,
    "phenotype_heritability" = get_phenotype_heritability
  )
  summary_data <- list()
  summary_data[[1]] <- compute_summary_data(config, population, summarisers)

  # Simulate the generations and capture some interesting data
  for (i in 2:(config$n_gen + 1)) {
    population_new <- simulate_generation(config, population)
    snp_pair_cors[[i]] <- attr(population_new, "snp_pair_cors")
    summary_data[[i]] <- compute_summary_data(
      config,
      population_new,
      summarisers
    )
    population <- population_new
  }

  attr(population, "class") <- "amsim"
  attr(population, "config") <- config
  attr(population, "config_path") <- config_path
  attr(population, "snp_pair_cors") <- snp_pair_cors
  attr(population, "summary_data") <- transform_summary_data(summary_data)

  return(population)
}

#' Simulate a number of population from a config file, with parallel options
#'
#' @param config_path Path to configuration file with simulation parameters
#' @param n_sims Number of populations to simulate before aggregating data
#' @param mc.cores How many cores to use for simulation (defaults to
#'  single-thread run)
#' @param initialise_population Logical specifying whether to start simulations
#'  from the same population each time
#'
#' @return A list which population, snp_pair_cors, and sib_cors objects.
#'
#' @importFrom pbmcapply pbmclapply
#'
#' @export
simulate_populations <- function(config_path, n_sims = 10, mc_cores = 1,
                                 initialise_population = FALSE) {
  config <- load_config(config_path)
  if (initialise_population) {
    init_population <- initialise_population(config)
  } else {
    init_population <- NULL
  }
  simulation_results <- pbmclapply(
    1:n_sims,
    function(i) simulate_population(config_path, init_population),
    mc.cores = mc_cores
  )
  return(simulation_results)
}

#' @export
summary.amsim <- function(object, ...) {
  # TODO Reimplement this when we have better extraction interface
}

#' @export
print.amsim <- function(object, ...) {
  config <- attr(object, "config")
  config_path <- attr(object, "config_path")

  causal_snps_idx <- which(config$phenotype_causal_snps > 0)
  causal_snps <- config$phenotype_causal_snps[causal_snps_idx]
  causal_snps_str <- paste0(names(causal_snps), " (effect = ", causal_snps, ")",
                            collapse = ", ")

  mating_model_pairs <- config$mating_model_pairs
  mating_model_pairs_str <- paste0(
    mating_model_pairs$male_snp, "-",
    mating_model_pairs$female_snp,
    " (cor = ", mating_model_pairs$correlation, ")"
  )

  cat("amsimr simulation object\n")
  cat("------------------------------------\n")
  cat("Configuration file:", normalizePath(config_path), "\n")
  cat("Random seed:", config$random_seed, "\n\n")

  cat("Number of generations simulated:", config$n_gen, "\n\n")
  cat("Summary functions:\n\n")

  cat("SNP minor allele frequencies:\n")
  print(config$snp_maf)
  cat("\n")

  cat("Phenotype name:", config$phenotype_name, "\n")
  cat("          heritability:", config$phenotype_heritability, "\n")
  cat("          causal SNPs:", causal_snps_str, "\n\n")

  cat("Mating model type:", config$mating_model_type, "\n")
  cat("             pairs:", mating_model_pairs_str, "\n")

  invisible(object)
}

#' @export
head.amsim <- function(object, ...) {
  head(object$population, ...)
}
