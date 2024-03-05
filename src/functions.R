source(file.path("src", "utils.R"))


simulate_methods <- function(
    n_sims, n_vox,
    mu_intercept, mu_me1, mu_me2, mu_interaction,
    sd_intercept, sd_me1, sd_me2, sd_interaction, sd_error,
    frac_imb1_train, frac_imb2_train, n_trial1_train, n_trial2_train,
    frac_imb1_test, frac_imb2_test, n_trial1_test, n_trial2_test,
    n_resamples, resample_to, sd_resample, dir_out, fn,
    methods = c("vanilla", "smoothed_boot", "boot"),
    ...
) {
  
  ## design:
  x_train <- make_design(
    frac_imb1 = frac_imb1_train, frac_imb2 = frac_imb2_train,
    n_trial1 = n_trial1_train, n_trial2 = n_trial2_train)
  x_test <- make_design(
    frac_imb1 = frac_imb1_test, frac_imb2 = frac_imb2_test, 
    n_trial1 = n_trial1_test, n_trial2 = n_trial2_test)
  idxs <- resample_idx(
    conditions = apply(x_train[, c("me1", "me2")], 1, paste, collapse = "_"),
    n_resamples = n_resamples, resample_to = resample_to, replace = TRUE)
  
  ## true voxel betas:
  mu_weights <- c(mu_intercept, mu_me1, mu_me2, mu_interaction)
  sd_weights <- c(sd_intercept, sd_me1, sd_me2, sd_interaction)
  b_true <- matrix(
    rnorm(n_vox*length(mu_weights), mean = mu_weights, sd = sd_weights),
    nrow = n_vox, byrow = TRUE)
  
  ## simulate data and fit:
  res <- future_map(
    seq_len(n_sims),
    \(sim_i) {
      y_train <- sim_data(b_true, x_train, sd_error)
      y_test  <- sim_data(b_true, x_test, sd_error)
      res_i <- data.table()
      if ("vanilla" %in% methods) {
        stat_vanilla <- cvlda(
          x_train = x_train[, "me1"], x_test = x_test[, "me1"], y_train = y_train,
          y_test = y_test)
        res_i <- cbind(res_i, stat_vanilla)
      }
      if ("smoothed_boot" %in% methods) {
        stat_smoothed <- cvlda_boot(
          x_train = x_train[, "me1"], x_test = x_test[, "me1"], y_train = y_train,
          y_test = y_test, idxs = idxs, sd_resample = sd_resample
        )
        res_i <- cbind(res_i, stat_smoothed)
      }
      if ("smoothed_boot" %in% methods) {
        stat_boot <- cvlda_boot(
          x_train = x_train[, "me1"], x_test = x_test[, "me1"], y_train = y_train,
          y_test = y_test, idxs = idxs, sd_resample = 0
        )
        res_i <- cbind(res_i, stat_boot)
      }
      res_i
    },
    .options = furrr_options(seed = TRUE)
  )
  
  d <- rbindlist(res)
  fwrite(d, file = here(dir_out, fn))
  
}


