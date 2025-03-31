#' Validate a configuration key-value pair
#'
#' Checks whether a key-value pair in a configuration section meets the
#' specified validation rules. Throws an error if the key or value is invalid.
#'
#' @param config A list representing the configuration.
#' @param section A character string indicating the configuration section.
#' @param key A character string specifying the key.
#' @param value The value associated with the key.
#' @param rule A list of validation functions for the key and value.
#'
#' @noRd
validate_key_value <- function(config, section, key, value, rule) {
  # Check if the key and value are both valid
  key_invalid <- !is.null(rule$key_val) && !rule$key_val(config, names(value))
  val_invalid <- !is.null(rule$val) && !rule$val(config, value)

  # Throw error if either of is not according to specification
  if (key_invalid || val_invalid) {
    stop("Invalid key or value in ", key, " in ", section)
  }
}

#' Validate a section of the configuration against rules
#'
#' Validates a section of the simulation configuration file, i.e. whether
#' a key within the section is required and how to validate it, using a
#' validation function.
#'
#' @param config The configuration object
#' @param section The section name to validate
#' @param rules A list of validation rules
#' @noRd
validate_section <- function(config, section, rules) {

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


#' Validate the entire configuration file
#'
#' Validates that the configuration file has all required sections and
#' that all required parameters have valid values.
#'
#' @param config The configuration object to validate
#' @noRd
validate_config <- function(config) {
  # Validate simulation section
  validate_section(config, section = "simulation", list(
    n_loci = list(required = TRUE, val = is_positive_int),
    n_pop = list(required = TRUE, val = is_positive_int),
    n_gen = list(required = TRUE, val = is_positive_int),
    n_iter = list(required = TRUE, val = is_positive_int),
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

#' Process the configuration file
#'
#' Extracts and transforms data registered in configuration file
#' to be used in the simulation.
#'
#' @param config The configuration object to process
#'
#' @return A list object with data processed from the configuration file
#' @noRd
process_config <- function(config) {

  # Extract some variables from the config file
  config_data <- list(
    n_loci = config$simulation[["n_loci"]],
    n_pop = config$simulation[["n_pop"]],
    n_gen = config$simulation[["n_gen"]],
    n_iter = config$simulation[["n_iter"]]
  )

  if (!is.null(config$simulation[["random_seed"]])) {
    config_data["random_seed"] <- config$simulation[["random_seed"]]
  }

  # Fill in SNP MAFs with default value if unspecified
  snp_maf <- config$snps[["snp_maf"]]
  snp_ids <- paste0("rs", 1:config_data$n_loci)
  snp_default <- setdiff(snp_ids, names(snp_maf))
  snp_maf[snp_default] <- config$snps[["default_maf"]]

  config_data[["snp_maf"]] <- unlist(snp_maf)

  # Extract phenotype data into phenotype_<key>
  config_data["phenotype_name"] <- ifelse(!is.null(config$phenotype[["name"]]),
                                          config$phenotype[["name"]], "X")
  config_data["phenotype_heritability"] <- config$phenotype[["heritability"]]

  causal_snps <- config$phenotype[["causal_snps"]]
  neutral_snps <- setdiff(snp_ids, names(causal_snps))
  causal_snps[neutral_snps] <- 0

  config_data[["phenotype_causal_snps"]] <- unlist(causal_snps)

  # Extract mating_model data into mating_model_<key>. In particular, extract
  # allele pairs into a tabular / data frame format.
  mating_model_pairs <- config$mating_model$pairs
  config_data["mating_model_type"] <- config$mating_model[["type"]]
  config_data[["mating_model_pairs"]] <- do.call(rbind, mating_model_pairs) |>
    as.data.frame()

  return(config_data)
}

#' Load a simulation configuration file
#'
#' Validates and processes data from configuration file to be used
#' in the simulation.
#'
#' @param config_path Path to configuration file
#' @returns A list object with data processed from the configuration file
#'
#' @importFrom yaml read_yaml
#'
#' @export
load_config <- function(config_path) {
  # Read the configuration file
  config <- read_yaml(config_path)

  validate_config(config)

  # Convert into data used in simulation
  config <- process_config(config)

  return(config)
}
