simulate_pop <- function(config_path, quietly = TRUE) {
  config <- load_config(config_path)
  init_population <- initialise_population(config)

  for (i in 1:config$n_gen) {
    print(paste0("Simulating generation #", i))
  }

  print(init_population)
}
