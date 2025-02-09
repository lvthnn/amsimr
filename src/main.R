population <- function(params) {
  #' Generate an initial population
  #'
  #' Generate the 0th generation of a population.
  #' @param params Path to YAML config file with simulation params.

  config <- config::get(file = params)
  n <- config$n_pop
  p <- config$n_gene
  q <- config$n_vars

  sex <- rbinom(n, 1, 0.5)
  maf <- rbeta(p, 2, 5)
  gene <- matrix(rbinom(n * p, 2, rep(maf, each = n)), nrow = n, ncol = p)
  vars <- matrix(rnorm(n * q), nrow = n, ncol = q)
  colnames(gene) <- paste0("rs", 1:p)

  gen0 <- tibble::tibble(sex, data.frame(gene), data.frame(vars))
  obj <- list(params = config, maf_genes = maf, data = gen0)
  class(obj) <- "population"

  return(obj)
}

pop <- population(params = "Github/AssocMating/config.yml")