# Validator generator - creates vectorized validator functions
is_val <- function(check_fn) {
  validator <- function(config, x) {
    return(all(sapply(x, function(val) check_fn(config, val))))
  }
  return(validator)
}

is_val_along <- function(check_fn) {
  validator <- function(config, x) {
    for (i in seq_along(x)) {
      xi <- x[[i]]
      result <- tryCatch({
        check_fn(config, xi)
      }, error = function(e) {
        stop("Validation failed for entry #", i, ": ", e$message)
        return(FALSE)
      })

      if (!result) stop("Validation failed for entry ", i)
    }

    return(TRUE)
  }

  return(validator)
}

# Basic validator functions (now all work with both scalar and vector inputs)
is_int <- is_val(function(config, x) is.integer(x))
is_pass <- is_val(function(config, x) TRUE)
is_numeric <- is_val(function(config, x) is.numeric(x))
is_positive <- is_val(function(config, x) is.numeric(x) && x > 0)
is_positive_int <- is_val(function(config, x) as.integer(x) == x && x > 0)
is_probability <- is_val(function(config, x) is.numeric(x) && x >= 0 && x <= 1)
is_correlation <- is_val(function(config, x) is.numeric(x) && x >= -1 && x <= 1)

# Option validator
is_option <- function(config, options) {
  return(is_val(function(config, x) x %in% options))
}

# Check whether a given string is a valid RS ID
is_rs <- function(config, x) {
  valid_snps <- paste0("rs", 1:config$simulation$n_loci)
  val_fn <- is_option(config, valid_snps)

  return(val_fn(config, x))
}

# Validator for pairs in mating_model section
is_pair <- function(config, x) {
  required_keys <- c("male_snp", "female_snp", "correlation")
  missing_keys <- setdiff(names(x), required_keys)

  if (length(missing_keys) > 0) {
    stop("Pairs in 'mating_model' missing keys '", missing_keys, "'")
  }

  valid_snps <- paste0("rs", 1:config$simulation$n_loci)

  # Check correlation value
  if (!is_correlation(config, x$correlation)) {
    stop("Invalid correlation value in mating_model pair: ", x$correlation)
  }

  # Check SNP names
  if (!x$male_snp %in% valid_snps) {
    stop("Invalid male_snp in mating_model: ", x$male_snp)
  }
  if (!x$female_snp %in% valid_snps) {
    stop("Invalid female_snp in mating_model: ", x$female_snp)
  }

  return(TRUE)
}

is_pairs <- is_val_along(is_pair)

