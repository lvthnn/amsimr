#' @export
MatingModel <- R6Class(
  "MatingModel",

  public = list(
    phenotypes = NULL,
    covariances = NULL,
    error_tol = 1e-4,

    initialize = function(phenotypes = NULL, covariances = NULL,
                          error_tol = 1e-4) {
      assert_list(phenotypes, types = "character", null.ok = TRUE)
      assert_class(covariances, classes = "matrix", null.ok = TRUE)
      assert_numeric(error_tol, null.ok = TRUE)
      invisible(self)
    },

    get_consolidated_matrix = function() {
      # collapses all phenotype pair correlation matrices into a single matrix,
      # and we take non-zero correlation entries and compress into a vector. we
      # also want a vector with indices
    }
  )
)

#' Generate initial solution matrix for mate matching optimisation routine
#'
#' @param population Population data frame with sex column and SNP genotypes
#' @param config Configuration list with mating_model_pairs
#'
#' @return List with init_sol matrix, psi_vec, and female_swap_idx
#'
#' @noRd
generate_init_state <- function(config, population) {
  # Extract mating model pairs from flattened config
  mating_pairs <- config$mating_model_pairs
  male_snps <- unlist(unique(mating_pairs$male_snp))
  female_snps <- unlist(unique(mating_pairs$female_snp))

  males <- cbind(
    scale(population[population$sex == 0, male_snps, drop = FALSE]),
    id = as.integer(rownames(population[population$sex == 0, ]))
  )

  females <- cbind(
    scale(population[population$sex == 1, female_snps, drop = FALSE]),
    id = as.integer(rownames(population[population$sex == 1, ]))
  )
  females <- females[sample(nrow(females)), ]

  # Set consistent row names and prefix female columns
  rownames(males) <- rownames(females) <- 1:(config$n_pop / 2)
  colnames(females) <- paste0("f_", colnames(females))

  # Combine into initial solution matrix
  init_sol <- as.matrix(cbind(males, females))

  # Initialize psi_vec to match dimensions of the input
  snp_pairs <- matrix(0, nrow = nrow(mating_pairs), ncol = 3)
  male_cols <- match(mating_pairs[, "male_snp"], colnames(init_sol)) - 1
  female_cols <- match(paste0("f_", mating_pairs[, "female_snp"]),
                       colnames(init_sol)) - 1
  snp_pairs[, 1] <- male_cols
  snp_pairs[, 2] <- female_cols
  snp_pairs[, 3] <- as.numeric(mating_pairs$correlation)
  storage.mode(snp_pairs) <- "numeric"

  # Get indices of female columns for swapping
  female_swap_idx <- which(grepl("^f_", colnames(init_sol))) - 1

  # Return as a list
  return(list(
    init_sol = init_sol,
    snp_pairs = snp_pairs,
    female_swap_idx = female_swap_idx
  ))
}

#' Generate an approximately optimal mate matching
#'
#' @param config Configuration file with simulation parameters
#' @param population Population data frame
#' @param diagnostics Whether to collect diagnostics from simulated annealing
#'  routine to evaluate performance
#'
#' @return List with init_sol matrix, psi_vec, and female_swap_idx
#'
#' @noRd
generate_matching <- function(config, population, collect_diagnostics = FALSE) {
  is_am <- config$mating_model_type == "assortative"
  sol_data <- generate_init_state(config, population)

  # Run the simulated annealing schedule to approximate the optimal matching
  if (is_am) {
    optim_sol <- optim_matching(
      sol_mat = sol_data$init_sol,
      snp_pairs = sol_data$snp_pairs,
      female_swap_idx = sol_data$female_swap_idx,
      num_iterations = config$n_iter,
      collect_diagnostics = collect_diagnostics
    )
  } else {
    optim_sol <- list()
    optim_sol$sol_mat <- sol_data$init_sol
  }

  colnames(optim_sol$sol_mat) <- colnames(sol_data$init_sol)
  optim_sol$sol_mat <- as.data.frame(optim_sol$sol_mat)

  matching <- list()
  matching[["male_idx"]] <- optim_sol$sol_mat$id
  matching[["female_idx"]] <- optim_sol$sol_mat$f_id
  if (is_am) {
    attr(matching, "call_params") <- optim_sol$call_params
    if (collect_diagnostics) {
      attr(matching, "diagnostics") <- list(
        energy_vals = optim_sol$energy_vals,
        energy_deltas = optim_sol$energy_deltas,
        acceptance_probs = optim_sol$acceptance_probs,
        temp_vals = optim_sol$temp_vals
      )
    }
  }

  # Compute the attained correlation in the matching based on SNP pairs
  matching_data <- optim_sol$sol_mat

  snp_pair_cors <- sapply(
    seq_len(nrow(config$mating_model_pairs)),
    function(i) {
      male_snp <- config$mating_model_pairs[[i, 1]]
      female_snp <- paste0("f_", config$mating_model_pairs[[i, 2]])
      cor(matching_data[, male_snp], matching_data[, female_snp])
    }
  )

  names(snp_pair_cors) <- paste(config$mating_model_pairs$male_snp,
                                config$mating_model_pairs$female_snp,
                                sep = "-")

  attr(matching, "snp_pair_cors") <- snp_pair_cors
  attr(matching, "class") <- "matching"

  return(matching)
}

#' @export
print.matching <- function(matching) {
  default_values <- formals(optim_matching)[c(
    "num_iterations",
    "init_temp",
    "temp_decay",
    "auto_temp_samples",
    "auto_accept_ratio",
    "collect_diagnostics",
    "quietly"
  )]
  call_params <- attr(matching, "call_params")
  is_default_value <- sapply(names(call_params), function(param) {
    call_params[[param]] == default_values[[param]]
  })

  # Print the summary of the matching object
  cat("amsimr matching object\n")
  cat("------------------------------------\n")
  cat("Number of pairs:", length(matching$male_idx), "\n")
  cat("Function call parameters:\n")
  for (param in names(call_params)) {
    param_value <- paste0("  ", param, " : ",
                          call_params[[param]],
                          ifelse(is_default_value[[param]], " (default)", ""),
                          "\n")
    cat(param_value)
  }
  cat("SNP pair correlations:\n")
  for (i in seq_along(attr(matching, "snp_pair_cors"))) {
    cat(" ", names(attr(matching, "snp_pair_cors"))[i], ":",
        attr(matching, "snp_pair_cors")[i], "\n")
  }
}

#' @export
summary.matching <- function(matching) {
  if (is.null(attr(matching, "diagnostics"))) {
    stop("Summary not supported for instance with collect_diagnostics = FALSE")
  }
  diagnostics_df <- do.call(cbind, attr(matching, "diagnostics"))
  colnames(diagnostics_df) <- names(attr(matching, "diagnostics"))
  summary(diagnostics_df)
}

#' @export
plot.matching <- function(matching) {
  if (is.null(attr(matching, "diagnostics"))) {
    stop("Plotting not supported for instance with collect_diagnostics = FALSE")
  }

  num_iterations <- attr(matching, "call_params")[["num_iterations"]]
  diagnostics <- attr(matching, "diagnostics")
  diagnostics_names <- attr(matching, "diagnostics") |> names()
  diagnostic_plot_types <- c("l", "l", "p", "l")
  names(diagnostic_plot_types) <- diagnostics_names

  par(mfrow = c(2, 2))
  for (diagnostic in diagnostics_names) {
    plot(
      x = 1:num_iterations,
      y = diagnostics[[diagnostic]],
      xlab = "Iteration",
      ylab = diagnostic,
      type = diagnostic_plot_types[diagnostic]
    )
  }
}
