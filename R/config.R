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
  valid_snps <- paste0("rs", 1:config$simulation$n_genes)
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

  valid_snps <- paste0("rs", 1:config$simulation$n_genes)

  # Check correlation value
  if (!is_correlation(config, x$correlation)) {
    stop("Invalid correlation value in mating_model pair: ", x$correlation)
  }

  # Check SNP names
  valid_snps <- paste0("rs", 1:config$simulation$n_genes)
  if (!x$male_snp %in% valid_snps) {
    stop("Invalid male_snp in mating_model: ", x$male_snp)
  }
  if (!x$female_snp %in% valid_snps) {
    stop("Invalid female_snp in mating_model: ", x$female_snp)
  }

  return(TRUE)
}

is_pairs <- is_val_along(is_pair)

# Helper function to validate key-value pair
validate_key_value <- function(config, section, key, value, rule) {
  # Check if the key and value are both valid
  key_invalid <- !is.null(rule$key_val) && !rule$key_val(config, names(value))
  val_invalid <- !is.null(rule$val) && !rule$val(config, value)

  # Throw error if either of is not according to specification
  if (key_invalid || val_invalid) {
    stop("Invalid key or value in ", key, " in ", section)
  }
}

validate_section <- function(config, section, rules) {
  #' Validate a section of the configuration against rules
  #'
  #' Validates a section of the simulation configuration file, i.e. whether
  #' a key within the section is required and how to validate it, using a
  #' validation function.
  #'
  #' @param config The configuration object
  #' @param section The section name to validate
  #' @param rules A list of validation rules
  if (!section %in% names(config)) {
    stop("Missing required section '", section, "'")
  }

  section_data <- config[[section]]

  for (key in names(rules)) {
    rule <- rules[[key]]
    key_exists <- key %in% names(section_data)

    # Check if required key exists
    if (rule$required && !key_exists) {
      stop("Missing required parameter: ", section, ".", key)
    }

    # If key exists, validate it
    if (!key_exists) next

    value <- section_data[[key]]

    validate_key_value(config, section, key, value, rule)
  }
}


validate_config <- function(config) {
  #' Validate the entire configuration file
  #'
  #' Validates that the configuration file has all required sections and
  #' that all required parameters have valid values.
  #'
  #' @param config The configuration object to validate

  # Validate simulation section
  validate_section(config, section = "simulation", list(
    n_genes = list(required = TRUE, val = is_positive_int),
    n_pop = list(required = TRUE, val = is_positive_int),
    random_seed = list(required = FALSE, val = is_positive_int)
  ))

  # Validate SNPs section
  validate_section(config, section = "snps", rules = list(
    default_maf = list(required = TRUE, val = is_probability),
    snp_maf = list(required = TRUE, key_val = is_rs, val = is_probability)
  ))

  # Validate phenotype section
  validate_section(config, section = "phenotype", rules = list(
    name = list(required = FALSE, val = is_pass),
    heritability = list(required = TRUE, val = is_probability),
    causal_snps = list(required = TRUE, key_val = is_rs, val = is_numeric)
  ))

  # Validate mating model section
  mating_types <- c("assortative", "random")

  validate_section(config, section = "mating_model", rules = list(
    type = list(required = TRUE, val = is_option(config, mating_types)),
    pairs = list(required = TRUE, val = is_pairs)
  ))
}

process_config <- function(config) {
  #' Process the configuration file
  #'
  #' Extracts and transforms data registered in configuration file
  #' to be used in the simulation.
  #'
  #' @param config The configuration object to process
  #' @return A list object with data processed from the configuration file

  # Extract some variables from the config file
  config_data <- list(
    n_genes = config$simulation[["n_genes"]],
    n_pop = config$simulation[["n_pop"]],
    random_seed = config$simulation[["random_seed"]]
  )

  # Fill in SNP MAFs with default value if unspecified
  snp_maf <- config$snps[["snp_maf"]]
  snp_ids <- paste0("rs", 1:config_data$n_genes)
  snp_default <- setdiff(snp_ids, names(snp_maf))
  snp_maf[snp_default] <- config$snps[["default_maf"]]

  config_data[["snp_maf"]] <- unlist(snp_maf)

  # Extract phenotype data into phenotype_<key>
  config_data["phenotype_name"] <- config$phenotype[["name"]]
  config_data["phenotype_heritability"] <- config$phenotype[["heritability"]]
  config_data[["phenotype_tias"]] <- unlist(config$phenotype[["causal_snps"]])

  # Extract mating_model data into mating_model_<key>. In particular, extract
  # allele pairs into a tabular / data frame format.
  mating_model_pairs <- config$mating_model$pairs
  config_data["mating_model_type"] <- config$mating_model[["type"]]
  config_data[["mating_model_pairs"]] <- do.call(rbind, mating_model_pairs)

  return(config_data)
}

load_config <- function(config_path) {
  #' Load a simulation configuration file
  #'
  #' Validates and processes data from configuration file to be used
  #' in the simulation.
  #'
  #' @param config_path
  #' @returns A list object with data processed from the configuration file

  # Read the configuration file
  config <- yaml::read_yaml(config_path)

  validate_config(config)

  # Convert into data used in simulation
  config <- process_config(config)

  return(config)
}
