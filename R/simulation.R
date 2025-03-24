run_simulation <- function(config_path, n_generations = 5) {
  config <- parse_config(config_path)
  init_population <- initialise_population(config)
  print(config)
}
