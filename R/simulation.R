# Compute allele frequencies for a generation
get_snp_frequencies <- function(config, population) {
  rs_cols <- paste0("rs", 1:config$n_loci)
  rs_freqs <- colMeans(population[, rs_cols] / 2)
  return(rs_freqs)
}

# Compute allele correlation matrix for a generation
get_snp_correlation <- function(config, population) {
  rs_cols <- paste0("rs", 1:config$n_loci)
  snp_matrix <- population[, rs_cols]
  return(cor(snp_matrix))
}

# Genetic correlation between siblings
get_genetic_correlation <- function(config, population) {
  n_pop <- config$n_pop
  population_sibs <- population[order(population$sibling_id), ]
  snp_cols <- paste0("rs", 1:config$n_loci)

  sibs_male <- population_sibs[seq(1, n_pop, by = 2), snp_cols]
  sibs_female <- population_sibs[seq(2, n_pop, by = 2), snp_cols]

  # Compute the genetic correlation between siblings
  return(diag(cor(sibs_male, sibs_female, use = "pairwise.complete.obs")))
}

# Phenotype correlation between siblings
get_phenotype_correlation <- function(config, population) {
  n_pop <- config$n_pop
  phenotype_col <- config$phenotype_name
  population_sibs <- population[order(population$sibling_id), ]

  sibs_male <- population_sibs[seq(1, n_pop, by = 2), phenotype_col]
  sibs_female <- population_sibs[seq(2, n_pop, by = 2), phenotype_col]

  # Compute the phenotype correlation between siblings
  return(cor(sibs_male, sibs_female))
}

# Heritability estimate of the phenotype
# TODO Find a suitable method to compute heritability estimate
get_heritability_estimate <- function(config, population) {
  return(NA)
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
compute_summary <- function(generations, summarisers = NULL) {
  # Initialize a list to store the results
  summary_results <- list()
  if (summarisers == NULL) return(list())

  # Loop through each generation
  for (i in seq_along(generations)) {
    generation <- generations[[i]]
    generation_summary <- list()

    # Loop through each summariser function
    for (name in names(summarisers)) {
      summariser_fn <- summarisers[[name]]
      generation_summary[[name]] <- summariser_fn(generation)
    }

    # Store the summary for this generation
    summary_results[[i]] <- generation_summary
  }

  return(summary_results)
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

  # Compute genetic similarity between siblings
  snp_pair_cors <- attr(mate_matching, "snp_pair_cors")
  attr(population_offspring, "snp_pair_cors") <- snp_pair_cors

  return(population_offspring)
}

#' Simulate a population from a config file
#'
#' @param config_path Path to configuration file with simulation parameters
#' @param quietly Display progress bar when running simulation
#'
#' @return A list which population, snp_pair_cors, and sib_cors objects.
#'
#' @importFrom pbapply pblapply
#'
#' @export
simulate_population <- function(config_path, progress = FALSE) {
  config <- load_config(config_path)
  init_population <- initialise_population(config)
  population <- init_population

  # Simulate the generations and capture some interesting data
  apply_fn <- ifelse(progress, lapply, pblapply)
  generations <- apply_fn(1:config$n_gen, function(i) {
    population_new <- simulate_generation(config, population)
    population <- population_new$population
    return(population)
  })

  summary_data <- compute_summary(generations, summarisers = list(
    "snp_frequencies" = get_snp_frequencies,
    "genetic_correlation" = get_genetic_correlation,
    "phenotype_correlation" = get_phenotype_correlation,
    "heritability_estimate" = get_heritability_estimate
  ))

  # Extract snp_pair_cors and sib_cors into respective lists
  snp_pair_cors <- lapply(generations, `[[`, 1)
  sib_cors <- lapply(generations, `[[`, 2)

  result <- list(
    population = population,
    snp_pair_cors = snp_pair_cors,
    sib_cors = sib_cors
  )

  attr(result, "class") <- "amsim"
  attr(result, "config") <- config
  attr(result, "config_path") <- config_path
  attr(result, "summary_data") <- summary_data

  return(result)
}

#' Simulate a number of population from a config file, with parallel options
#'
#' @param config_path Path to configuration file with simulation parameters
#' @param n_sims Number of populations to simulate before aggregating data
#' @param progress Display progress bar when running simulation
#'
#' @return A list which population, snp_pair_cors, and sib_cors objects.
#'
#' @importFrom pbapply pblapply
#'
#' @export
simulate_populations <- function(config_path, n_sims = 10, progress = FALSE) {
  return(TRUE)
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

  cat("Population size:", config$n_pop, "\n")
  cat("Number of loci:", config$n_loci, "\n")
  cat("       of generations:", config$n_gen, "\n\n")

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
