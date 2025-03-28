simulate_pop <- function(config_path, quietly = TRUE) {
  config <- load_config(config_path)
  init_population <- initialise_population(config)

  # Prepare list of SNP correlations and siblings correlations
  snp_pair_cors <- list()
  sib_cors <- list()

  population <- init_population

  for (i in 1:config$n_gen) {
    if (!quietly) print(paste0("Simulating generation #", i))
    population_new <- simulate_generation(config, population)
    snp_pair_cors[[i]] <- population_new$snp_pairs_cor
    sib_cors[[i]] <- population_new$sib_cor

    population <- population_new$population
  }

  return(list(
    population = population,
    snp_pair_cors = snp_pair_cors,
    sib_cors = sib_cors
  ))
}

summarise_sib_cors <- function(sib_cors_list) {
  # Extract the relevant components
  extract_data <- function(x) {
    c(
      mean_genetic_correlation = x$mean_genetic_correlation,
      mean_genotype_similarity = x$mean_genotype_similarity,
      locus_correlation_sd = x$locus_correlation_sd,
      height_correlation = x$height_correlation,
      height_heritability = x$height_heritability
    )
  }

  # Apply extraction to all list elements
  data_list <- lapply(sib_cors_list, extract_data)

  # Combine into a data frame
  df <- do.call(rbind, data_list)

  # Add locus correlations as columns
  locus_correlations <- do.call(
    rbind,
    lapply(sib_cors_list, function(x) x$locus_correlations)
  )

  # Combine locus correlations with the main dataframe
  final_df <- as.data.frame(cbind(df, locus_correlations))
  colnames(final_df)[6:ncol(final_df)] <- paste0("rs", 1:(ncol(final_df) - 5))

  return(final_df)
}

