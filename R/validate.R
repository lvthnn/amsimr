#' Create a vectorized validator function
#'
#' Wraps a check function to apply over all elements of `x`.
#'
#' @param check_fn A function of the form (config, x) returning TRUE/FALSE
#'
#' @noRd
is_val <- function(check_fn) {
  validator <- function(config, x) {
    return(all(sapply(x, function(val) check_fn(config, val))))
  }
  return(validator)
}

#' Create a validator that checks each element and reports index on failure
#'
#' Provides detailed error messaging during validation.
#'
#' @param check_fn A function of the form (config, x) returning TRUE/FALSE
#'
#' @noRd
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

#' Check if input is an integer
#'
#' @param config Parsed configuration object
#' @param x Object to validate
#'
#' @noRd
is_int <- is_val(function(config, x) is.integer(x))

#' Always passes validation
#'
#' @param config Parsed configuration object
#' @param x Object to validate
#'
#' @noRd
is_pass <- is_val(function(config, x) TRUE)

#' Check if input is numeric
#'
#' @param config Parsed configuration object
#' @param x Object to validate
#'
#' @noRd
is_numeric <- is_val(function(config, x) is.numeric(x))

#' Check if input is a positive number
#'
#' @param config Parsed configuration object
#' @param x Object to validate
#'
#' @noRd
is_positive <- is_val(function(config, x) is.numeric(x) && x > 0)

#' Check if input is a positive integer
#'
#' @param config Parsed configuration object
#' @param x Object to validate
#'
#' @noRd
is_positive_int <- is_val(function(config, x) as.integer(x) == x && x > 0)

#' Check if input is a probability (in [0, 1])
#'
#' @param config Parsed configuration object
#' @param x Object to validate
#'
#' @noRd
is_probability <- is_val(function(config, x) is.numeric(x) && x >= 0 && x <= 1)

#' Check if input is a correlation (in [-1, 1])
#' @param config Parsed configuration object
#' @param x Object to validate
#' @noRd
is_correlation <- is_val(function(config, x) is.numeric(x) && x >= -1 && x <= 1)

#' Create an option validator
#'
#' Returns a validator that checks membership in a set of options.
#'
#' @param config Parsed configuration object
#' @param options Character vector of allowed values
#'
#' @noRd
is_option <- function(config, options) {
  return(is_val(function(config, x) x %in% options))
}

#' Check if string is a valid RS ID
#'
#' Validates RS IDs based on the number of loci in config.
#'
#' @param config Parsed configuration object
#' @param x Object to validate
#'
#' @noRd
is_rs <- function(config, x) {
  valid_snps <- paste0("rs", 1:config$simulation$n_loci)
  val_fn <- is_option(config, valid_snps)

  return(val_fn(config, x))
}

#' Validate a pair entry in the mating model
#'
#' Ensures required keys and valid SNPs/correlation values.
#'
#' @param config Parsed configuration object
#' @param x A list representing a trait pair with required fields
#'
#' @noRd
is_pair <- function(config, x) {
  required_keys <- c("male_snp", "female_snp", "correlation")
  missing_keys <- setdiff(names(x), required_keys)

  if (length(missing_keys) > 0) {
    stop("Pairs in 'mating_model' missing keys '", missing_keys, "'")
  }

  valid_snps <- paste0("rs", 1:config$simulation$n_loci)

  if (!is_correlation(config, x$correlation)) {
    stop("Invalid correlation value in mating_model pair: ", x$correlation)
  }

  if (!x$male_snp %in% valid_snps) {
    stop("Invalid male_snp in mating_model: ", x$male_snp)
  }
  if (!x$female_snp %in% valid_snps) {
    stop("Invalid female_snp in mating_model: ", x$female_snp)
  }

  return(TRUE)
}

#' Validate a list of pair entries
#'
#' Applies `is_pair` to each element in the list.
#'
#' @param config Parsed configuration object
#' @param x A list of pair entries
#'
#' @noRd
is_pairs <- is_val_along(is_pair)

