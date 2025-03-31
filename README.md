# amsimr: Assortative Mating Simulation Tool in R

R package developed for my B.Sc. thesis, *Distant LD and its effects on polygenic
score construction and prediction: a simulation study*. The package implements an algorithm
to simulate non-random correlation between allele pairs.

## Installation

To install the package, clone this repository and load the package in R:

```r
# Install the package
devtools::install_github("lvthnn/amsimr")
```

## Setting up a Configuration File

To run a simulation, you will need to specify a configuration file. The required
sections are `simulation`, `snps`, `phenotype`, and `mating_model`. Shown below is
an example:

```yaml
simulation:
  n_loci: 10       # Number of genetic loci
  n_pop: 10000     # Population size
  n_iter: 10       # Number of generations to simulate
  random_seed: 42  # Random seed for reproducibility

snps:
  default_maf: 0.20  # Default minor allele frequency
  snp_maf:
    rs1: 0.15
    rs2: 0.35
    rs3: 0.05

phenotype:
  name: "height"      # Name of phenotype colum in simulated data
  heritability: 0.80  # Narrow-sense phenotype heritability
  causal_snps:
    rs1: 0.5  # Effect size of SNP rs1
    rs2: 0.3  # Effect size of SNP rs2

mating_model:
  type: "assortative"
  pairs:
    - male_snp: rs1
      female_snp: rs3
      correlation: 0.25

    - male_snp: rs1
      female_snp: rs2
      correlation: 0.75
```
You can create a template configuration file in a location of your choosing by running
```r
library(amsimr)

from_template(path = "path_to_config_file.yaml")
```
which you can succeedingly customise.

## Example Usage

A minimal example of how to use the package is shown below:

```r
# Run a simulation with the specified configuration file
sim <- simulate_pop(config = "config.yaml", progress = TRUE)
```

## Output

The simulation produces a dataset containing information about individuals across
generations, including genotype, phenotype, and mating pairs. Example output:

```r
# View simulation information
print(sim)

# View the first few individuals in the population
head(sim)

# Summarise population statistics
summary(sim)
```

## Dependencies

This package relies on the following R packages:

- `Rcpp`
- `RcppArmadillo`
- `yaml`
- `pbapply`

## Citation

If you use this package in your research, please cite:

> Hlynsson, K. (2025) *Distant LD and its effects on polygenic score construction
> and prediction: a simulation study*, B.Sc. thesis, Faculty of Physical Sciences, University of Iceland.

## License

This project is licensed under the MIT License.

## Contact

For questions or issues, please contact me at [kari.hlynsson@gmail.com].

