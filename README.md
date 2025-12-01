# amsimr <img src="man/figures/logo.svg" align="right" width="120" />

### Overview
**amsimr** is an R package that provides bindings to the
[amsimcpp](https://github.com/lvthnn/amsimcpp) C++ library
for forward-time genetic simulations of populations exhibiting assortative
mating. Simulations can be configured in a multitude of ways, such as locus
properties and phenotype architectures and covariance structures.

Users more familiar with Python are also encouraged to use the
[amsimpy](https://github.com/lvthnn/amsimpy) package. 

### Installation
Install the latest stable version from GitHub (requires the
[devtools](https://github.com/r-lib/devtools) package):
```
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("lvthnn/amsimr", dependencies = TRUE, build_vignettes = FALSE)
```
This installation won't include the vignettes (they take some time to build).

*Note: The package will be released on CRAN in future versions.*

### Examples

Below is an example simulation configuration:
```r
library(amsimr)

sim <- Simulation$new()

sim$simulation(
  n_generations = 10,
  n_individuals = 1000,
  output_dir = "amsim_test",
  random_seed = 42
)

sim$genome(
  n_loci = 20000,
  locus_maf = 0.5,
  locus_recombination = 0.5,
  locus_mutation = 0.5
)

sim$phenome(
  n_phenotypes = 2,
  names = c("height", "grip_strength"),
  n_causal_loci = 10000,
  h2_genetic = 0.5,
  h2_environmental = 0.5,
  h2_vertical = 0.0,
  genetic_cor = diag(2),
  environmental_cor = diag(2)
)

sim$assortative_mating(
  mate_cor = matrix(
    c(0.2, 0.3, 0.4, 0.1),
    nrow = 2,
    ncol = 2
  )
)

sim$metrics(
  pheno_h2(),
  pheno_comp_cor("genetic")
)

sim$run(
  n_replicates = 100,
  n_processes = 10,
  summarise = TRUE,
  log_file = FALSE,
  log_level = "warning"
)
```

### Citation

If you use amsimr in your research, please cite amsimcpp:

```bibtex
@software{hlynsson2024_amsimcpp,
  author = {Hlynsson, Kári},
  title = {amsimcpp: A C++ library for forward-time assortative mating simulations},
  year = {2025},
  version = {0.1.0},
  url = {https://github.com/lvthnn/amsimcpp}
}
```
