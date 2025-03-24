# Assortative Mating (AM) Simulation Tool

R package developed for my B.Sc. thesis, *Distant LD and its effects on polygenic
score construction and prediction: a simulation study*. The package implements an algorithm
to simulate non-random correlation between allele pairs.

## Installation

To install the package, clone this repository and load the package in R:

```sh
# Clone the repository
git clone https://github.com/lvthnn/am_simulation.git
cd am_simulation

# Load the package in R
library(devtools)
devtools::load_all()
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
  name: "height"
  heritability: 0.80  # Heritability of the trait
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

## Example Usage

A minimal example of how to use the package is shown below:

```r
# Load the package
library(am_simulation)

# Run a simulation with the specified configuration file
pop <- simulate_pop(config = "config.yaml")

# View the first few individuals in the population
head(pop)
```

## Output

The simulation produces a dataset containing information about individuals across
generations, including genotype, phenotype, and mating pairs. Example output:

```r
# Summarize the population statistics
summary(pop)

# Extract genotype frequencies
table(pop$genotype)
```

## Customizing the Simulation

The configuration file allows for various customizations:

- **Number of loci** (`n_loci`): Controls the number of SNPs simulated.
- **Population size** (`n_pop`): Determines how many individuals are included.
- **Number of generations** (`n_iter`): Defines how many generations are simulated.
- **Mating model** (`mating_model`): Specifies the assortative mating structure.

## Dependencies

This package relies on the following R packages:

- `Rcpp`
- `RcppArmadillo`
- `yaml`

Ensure dependencies are installed before running the simulation:

```
install.packages(c("yaml"))
```

## Citation

If you use this package in your research, please cite:

> Hlynsson, K. (2025). *Distant LD and its effects on polygenic score construction
> and prediction: a simulation study*.

## License

This project is licensed under the MIT License.

## Contact

For questions or issues, please contact me at [kari.hlynsson@gmail.com].

