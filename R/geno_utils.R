#' Create an (n_population * n_loci) ancestral SNP matrix
#'
#' @param n_population Number of individuals in simulated population
#' @param n_loci Number of simulated genetic loci (SNPs)
#' @param snp_mafs Minor allele frequency vector of SNPs
#' @param ncores Number of cores to use for parallel processing
#' @param seed Seed for random number generation
#'
#' @noRd
.init_ancestral_snps <- function(n_population, n_loci, snp_mafs,
                                 ncores = parallel::detectCores() - 1,
                                 seed = NULL) {
  checkmate::assert_count(n_population)
  checkmate::assert_count(n_loci)
  checkmate::assert(n_population %% 2 == 0, .var.name = "even population size")
  checkmate::assert_numeric(snp_mafs, lower = 0, upper = 1, finite = TRUE,
                            any.missing = FALSE, len = n_loci)

  n_sex      <- n_population / 2
  chunk_size <- 10000

  genotypes_fbm <- bigstatsr::FBM(
    nrow        = n_population,
    ncol        = n_loci,
    type        = "unsigned short",
    backingfile = "genotypes_gen0"
  )

  # Create a list of chunks to fill the matrix
  start_rows <- seq(1, n_sex, by = chunk_size)
  end_rows   <- pmin(start_rows + chunk_size - 1, n_sex)
  chunk_data <- list(start = start_rows, end = end_rows)

  future::plan(future::multisession, workers = ncores)

  cli::cli_alert_info("Generating male genotype data...")
  furrr::future_walk(seq_along(chunk_data$start), ~ {
    start_row <- chunk_data$start[.x]
    end_row   <- chunk_data$end[.x]
    rows      <- end_row - start_row + 1

    male_chunk <- matrix(
      rbinom(rows * n_loci,
        size = 2,
        prob = rep(snp_mafs, each = rows)
      ),
      nrow = rows,
      ncol = n_loci
    )

    genotypes_fbm[start_row:end_row, ] <- male_chunk
  }, .options = furrr::furrr_options(seed = TRUE))

  cli::cli_alert_info("Computing correlation matrix for genotype data...")
  bigparallelr::set_blas_ncores(ncores)
  cor_male   <- bigstatsr::big_cor(genotypes_fbm, ind.row = 1:n_sex)
  chunk_data <- list(start = start_rows + n_sex, end = end_rows + n_sex)

  # Update chunk data for generation of female SNPs
  cli::cli_alert_info("Generating female genotype data...")
  furrr::future_walk(seq_along(chunk_data$start), ~ {
    start_row    <- chunk_data$start[.x]
    end_row      <- chunk_data$end[.x]
    n_row        <- end_row - start_row + 1
    scores_chunk <- MASS::mvrnorm(
      n_row,
      mu    = rep(0, n_loci),
      Sigma = cor_male[]
    )

    female_chunk <- matrix(0, nrow = n_row, ncol = n_loci)

    for (j in 1:n_loci) {
      scores_j          <- scores_chunk[, j]
      p_j               <- snp_mafs[j]
      p0                <- qnorm((1 - p_j)^2)
      p1                <- qnorm((1 - p_j)^2 + 2 * p_j * (1 - p_j))
      breaks            <- c(-Inf, p0, p1, Inf)
      female_chunk[, j] <- findInterval(scores_j, vec = breaks) - 1
    }
    genotypes_fbm[start_row:end_row, ] <- female_chunk
  }, .options = furrr::furrr_options(seed = TRUE))
  future::plan(future::sequential)
  return(genotypes_fbm)
}

#' Generates offspring for a Population object
#'
#' @param population A Population instance
#' @param ncores Number of cores to use for parallel processing