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
  snp_pairs_cor <- attr(mate_matching, "snp_pairs_cor")

  # Compute genetic similarity between siblings
  sib_cor <- compute_sib_cor(config, population_offspring)

  # Compile results into a list object
  result <- list()
  result[["population"]] <- population_offspring
  result[["snp_pairs_cor"]] <- snp_pairs_cor
  result[["sib_cor"]] <- sib_cor

  return(result)
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
simulate_pop <- function(config_path, progress = FALSE) {
  config <- load_config(config_path)
  init_population <- initialise_population(config)
  population <- init_population

  # Simulate the generations and capture some interesting data
  apply_fn <- ifelse(progress, lapply, pblapply)
  generations <- apply_fn(1:config$n_gen, function(i) {
    population_new <- simulate_generation(config, population)
    population <- population_new$population

    return(list(
      snp_pair_cor = population_new$snp_pairs_cor,
      sib_cor = population_new$sib_cor
    ))
  })

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

  return(result)
}

#' @export
summary.amsim <- function(object, ...) {
  # Extract the relevant components
  extract_data <- function(object) {
    c(
      mean_genetic_correlation = object$mean_genetic_correlation,
      mean_genotype_similarity = object$mean_genotype_similarity,
      locus_correlation_sd = object$locus_correlation_sd,
      height_correlation = object$height_correlation,
      height_heritability = object$height_heritability
    )
  }

  # Apply extraction to all list elements
  data_list <- lapply(object$sib_cors, extract_data)

  # Combine into a data frame
  df <- do.call(rbind, data_list)

  # Add locus correlations as columns
  locus_correlations <- do.call(
    rbind,
    lapply(object$sib_cors, function(x) x$locus_correlations)
  )

  # Combine locus correlations with the main dataframe
  final_df <- as.data.frame(cbind(df, locus_correlations))
  colnames(final_df)[6:ncol(final_df)] <- paste0("rs", 1:(ncol(final_df) - 5))

  return(final_df)
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
