generate_test_cor <- function(n) {
  A <- matrix(rnorm(n^2), n, n)
  S <- A %*% t(A)

  D <- diag(1 / sqrt(diag(S)))
  R <- D %*% S %*% D

  return(R)
}

# Calibrate tetrachoric correlation given a target correlation and locus MAFs
calibrate_cor <- function(loc_cor, loc_mafs, thresh) {
  coef <- sqrt(prod(loc_mafs) * prod(1 - loc_mafs))
  z_cor <- prod(loc_mafs) + loc_cor * coef

  optim_fun <- function(rho) {
    tetra_matrix <- diag(2)
    tetra_matrix[1, 2] <- tetra_matrix[2, 1] <- rho
    pr <- mvtnorm::pmvnorm(
      lower = thresh,
      upper = Inf,
      mean = 0,
      corr = tetra_matrix
    )[1]

    return(pr - z_cor)
  }

  sol <- uniroot(optim_fun, lower = -1, upper = 1)$root

  return(sol)
}

# Use memoisation
calibrate_cor_cached <- memoise::memoise(calibrate_cor)

calibrate_cor_matrix <- function(cor_mat, loc_mafs, thresh) {
  ld_cor <- cor_mat
  indices <- which(upper.tri(cor_mat), arr.ind = TRUE)

  for (i in seq_len(nrow(indices))) {
    l <- indices[i, 1]
    r <- indices[i, 2]
    ids <- c(l, r)
    ld_cor[l, r] <- calibrate_cor_cached(
      cor_mat[l, r],
      loc_mafs[ids],
      thresh[ids]
    )
  }

  ld_cor[lower.tri(ld_cor)] <- ld_cor[upper.tri(ld_cor)]
  as.matrix(Matrix::nearPD(ld_cor, corr = TRUE)$mat)
}

# Function to preprocess LD correlation matrix for MVN thresholding
# Shouldn't this function be recombination aware?

#' @param ld_cor Founder population genetic correlation matrix
#' @param loc_mafs Vector of locus minor allele frequencies
#' @param loc_recomb Vector of locus recombination probabilities
#' @param ld_block_size Size of LD blocks in founder correlation matrix
preprocess_ld <- function(
  ld_cor,
  loc_mafs,
  loc_recomb,
  ld_block_size = 50,
  n_cores = 5
) {
  n_loc <- nrow(ld_cor)
  thresh <- qnorm(1 - loc_mafs)

  # if all loc_recomb entries equal 0.5, there should not be any founder
  # correlation between loci -- error in this case

  if (all(loc_recomb == 0.5))
    cli::cli_abort("All locus recombination probabilities are 0.5; skipping
      optimisation.")

  # create ld blocks as largest possible sequence bounded between two 0.5
  # values in loc_recomb, and divive it into ld blocks of size ld_block_size

  chrom_sep <- which(loc_recomb == 0.5)

  chrom_blocks <- lapply(seq_along(chrom_sep), function(id) {
    start <- chrom_sep[id]
    if (id == length(chrom_sep)) {
      end <- n_loc
    } else {
      end <- chrom_sep[id + 1] - 1
    }
    return(ld_cor[start:end, start:end])
  })

  ld_blocks <- pbmcapply::pbmclapply(seq_along(chrom_blocks), function(chr_id) {
    chrom_block <- chrom_blocks[[chr_id]]
    global_offset <- chrom_sep[chr_id] - 1
    n_chr <- nrow(chrom_block)

    block_ids <- ceiling(seq_len(n_chr) / ld_block_size)

    lapply(split(seq_len(n_chr), block_ids), function(idx) {
      global_idx <- idx + global_offset
      chrom_extracted <- chrom_block[idx, idx, drop = FALSE]
      chrom_calibrated <- calibrate_cor_matrix(
        chrom_extracted,
        loc_mafs[global_idx],
        thresh[global_idx]
      )
      return(chol(chrom_calibrated))
    })
  }, mc.cores = n_cores)

  ld_blocks <- unlist(ld_blocks, recursive = FALSE)

  return(ld_blocks)
}
