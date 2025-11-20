sim <- amsimr::Simulation$new()
sim$simulation(
  n_generations = 5,
  n_individuals = 5000,
  output_dir = "test_simulation",
  random_seed = 456
)$genome(
  n_loci = 500,
  locus_maf = 0.3,
  locus_recombination = 0.5,
  locus_mutation = 1e-8
)$phenome(
  n_phenotypes = 2,
  names = c("height", "weight"),
  n_causal_loci = c(250, 250),
  h2_genetic = c(0.6, 0.5),
  h2_environmental = c(0.25, 0.3),
  h2_vertical = c(0.15, 0.2),
  genetic_cor = diag(2),
  environmental_cor = diag(2)
)$assortative_mating(
  mate_cor = matrix(c(0.2, 0.4, 0.3, 0.5), nrow = 2),
  tol_inf = 1e-6,
  n_iterations = 1e6,
  temp_init = 1.5,
  temp_decay = 0.995
)$metrics(
  list(
    pheno_h2(),
    pheno_comp_cor("total")
  )
)

sim$run(n_replicates = 100, n_threads = 10, log_file = FALSE)
