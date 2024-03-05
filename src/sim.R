library(here)
library(dplyr)
library(data.table)
library(tidyr)
library(purrr)
library(future)
library(furrr)
library(ggplot2)

source(file.path("src", "functions.R"))


## settings

n_cores <- availableCores()
n_workers <- length(unique(availableWorkers()))
cat("\ncores: ", availableCores(), "\n")
cat("\nworkers: ", unique(availableWorkers()), "\n")
if (interactive()) {
  plan(sequential)
} else {
  plan(list(
    tweak(multisession, workers = n_workers),
    tweak(multicore, workers = n_cores)
  ))
}



## sweep-00 -----
## binary classification
## between-class imbalance (different sized classes), both train and test
## balanced within-class

## Hyperparameters

pars <- 
  expand.grid(
    frac_imb1_train = c(0.1, 0.3, 0.5, 0.7, 0.9),
    frac_imb1_test = c(0.1, 0.3, 0.5, 0.7, 0.9),
    sd_error = c(1, 10, 100),
    n_resamples = c(100, 1000),
    resample_to = c(2, 5, 10) * 25  ## per unique level of conditions
  ) %>%
  filter(frac_imb1_train == frac_imb1_test) %>%
  mutate(
    n_sims = 1E2,
    n_vox = 100,
    frac_imb2_train = 0.5,
    frac_imb2_test = 0.5,
    n_trial1_train = 50,
    n_trial1_test = 50,
    n_trial2_train = 50,
    n_trial2_test = 50,
    mu_intercept = 0, 
    mu_me1 = 0,
    mu_me2 = 0,
    sd_intercept = 0,
    sd_me1 = 0,
    sd_me2 = 0,
    mu_interaction = 0,
    sd_interaction = 0,
    sd_resample = opt_smoothing(n_vox, n_trial = n_trial1_train + n_trial2_train, sd_error),
    fn = paste0(sprintf("%04s", 1:n()), ".csv")
  )

## Path to results

dir_out <- here("data", "sweep-00")
dir.create(dir_out, showWarnings = FALSE)
fwrite(pars, here("data", "hyperparameters_sweep-00.csv"))

## Run

set.seed(0)
future_pwalk(
  pars,
  simulate_methods,
  dir_out,
  .options = furrr_options(seed = TRUE)
)


