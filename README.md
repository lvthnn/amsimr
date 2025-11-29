# amsimr <img src="man/figures/logo.svg" align="right" width="120" />

### Overview
**amsimr** is an R package which exposes bindings to the
[amsimcpp](https://github.com/lvthnn/amsimcpp) C++ library
for forward-time genetic simulations of populations exhibiting assortative
mating. Simulations can be configured in a multitude of ways, such as locus
properties and phenotype architectures and covariance structures.

Users more familiar with Python are invited to use the
[amsimpy](https://github.com/lvthnn/amsimpy) package instead, which is nearly
identically implemented in terms of the user API, albeit with more sophisticated
multithreaded capabilities which are not available to the R language, resulting
in faster replicate simulations when using [user-specified metric summary
functions](vignettes/custom_metrics.Rmd).

### Getting started

If you are new to the **amsim** simulation family we recommend starting with the
tutorial vignettes and documentation available on the package site.

### Installation

* Install from CRAN:
```
install.packages("amsimr")
```
* Install latest development version from GitHub (requires the
[devtools](https://github.com/r-lib/devtools) package):
```
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("lvthnn/amsimr", dependencies = TRUE, build_vignettes = FALSE)
```
This installation won't include the vignettes (they take some time to build).

### Examples

Below is an example simulation configuration:
```r
library(amsimr)

sim <- Simulation$new()

sim$simulation(
    n_generations = 10,
    n_individuals = 1000,
    output_dir = "/tmp/sim_test",
    random_seed = 42
)

sim$genome(
    n_loci = 20000,
    locus_maf = 0.5,
    locus_recombination = 0.5,
    locus_mutation = 0.5
)
```
