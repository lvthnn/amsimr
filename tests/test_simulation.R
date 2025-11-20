# Test script for Simulation class
library(amsimr)

# Create a new simulation object
message("Creating Simulation object...")
sim <- amsimr::Simulation$new()

# Configure basic simulation parameters
message("Configuring simulation parameters...")
sim$simulation(
  n_generations = 10,
  n_individuals = 100,
  output_dir = tempdir(),
  random_seed = 42
)

# Verify simulation parameters
stopifnot(sim$n_generations == 10)
stopifnot(sim$n_individuals == 100)
stopifnot(sim$random_seed == 42)
message("✓ Simulation parameters configured correctly")

# Configure genome
message("Configuring genome...")
sim$genome(
  n_loci = 1000,
  locus_maf = 0.3,
  locus_recombination = 0.5,
  locus_mutation = 1e-8
)

# Verify genome parameters
stopifnot(sim$n_loci == 1000)
stopifnot(length(sim$locus_maf) == 1000)
stopifnot(all(sim$locus_maf == 0.3))
stopifnot(length(sim$locus_recombination) == 1000)
stopifnot(all(sim$locus_recombination == 0.5))
stopifnot(length(sim$locus_mutation) == 1000)
stopifnot(all(sim$locus_mutation == 1e-8))
message("✓ Genome configured correctly")

# Configure phenome
message("Configuring phenome...")
n_pheno <- 2
sim$phenome(
  n_phenotypes = n_pheno,
  names = c("trait1", "trait2"),
  n_causal_loci = c(100, 100),
  h2_genetic = c(0.5, 0.6),
  h2_environmental = c(0.3, 0.2),
  h2_vertical = c(0.2, 0.2),
  genetic_cor = diag(n_pheno),
  environmental_cor = diag(n_pheno)
)

# Verify phenome parameters
stopifnot(sim$n_phenotypes == 2)
stopifnot(all(sim$phenotype_names == c("trait1", "trait2")))
stopifnot(all(sim$n_causal_loci == c(100, 100)))
stopifnot(all(sim$h2_genetic == c(0.5, 0.6)))
stopifnot(all(sim$h2_environmental == c(0.3, 0.2)))
stopifnot(all(sim$h2_vertical == c(0.2, 0.2)))
stopifnot(length(sim$genetic_cor) == n_pheno^2)
stopifnot(length(sim$environmental_cor) == n_pheno^2)
message("✓ Phenome configured correctly")

# Test random mating
message("Configuring random mating...")
sim$random_mating()
message("✓ Random mating configured")

# Test method chaining
message("Testing method chaining...")
sim2 <- amsimr::Simulation$new()
sim2$simulation(
  n_generations = 5,
  n_individuals = 50,
  output_dir = tempdir(),
  random_seed = 123
)$genome(
  n_loci = 500,
  locus_maf = 0.25,
  locus_recombination = 0.5,
  locus_mutation = 1e-9
)$phenome(
  n_phenotypes = 1,
  names = "trait",
  n_causal_loci = 50,
  h2_genetic = 0.7,
  h2_environmental = 0.2,
  h2_vertical = 0.1,
  genetic_cor = matrix(1),
  environmental_cor = matrix(1)
)$random_mating()

stopifnot(sim2$n_generations == 5)
stopifnot(sim2$n_loci == 500)
stopifnot(sim2$n_phenotypes == 1)
message("✓ Method chaining works correctly")

# Test assortative mating configuration
message("Testing assortative mating configuration...")
sim3 <- amsimr::Simulation$new()
sim3$simulation(
  n_generations = 5,
  n_individuals = 50,
  output_dir = tempdir(),
  random_seed = 456
)$genome(
  n_loci = 500,
  locus_maf = 0.3,
  locus_recombination = 0.5,
  locus_mutation = 1e-8
)$phenome(
  n_phenotypes = 2,
  names = c("height", "weight"),
  n_causal_loci = c(50, 50),
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
)

# Verify assortative mating parameters
stopifnot(length(sim3$mate_cor) == 4)
stopifnot(sim3$tol_inf == 1e-6)
stopifnot(sim3$n_iterations == 1e6)
stopifnot(sim3$temp_init == 1.5)
stopifnot(sim3$temp_decay == 0.995)
message("✓ Assortative mating configured correctly")

message("\n=== All tests passed! ===")
