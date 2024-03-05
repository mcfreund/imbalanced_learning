is_integer <- function(x, tol = .Machine$double.eps) {
  min(abs(c(x%%1, x%%1-1))) < tol
}


resample <- function(x, ...) {
  x[sample.int(length(x), ...)]
}


resample_idx <- function(conditions, n_resamples, resample_to = NULL, replace = FALSE) {
  stopifnot(is.character(conditions) || is.numeric(n_resamples))

  groups_list <- split(seq_along(conditions), conditions)  ## sorts list elements alphabetically
  global_min <- Reduce(min, lapply(groups_list, length))
  if (is.null(resample_to)) {
    resample_to <- global_min
  } else {
    if (!replace && resample_to > global_min) stop("When replace=FALSE, cannot draw more obs than the global min.")
  }

  idx <- replicate(n_resamples, unlist(lapply(groups_list, resample, size = resample_to, replace = replace), use.names = FALSE))
  idx <- t(idx)  ## resamples as rows, conditions as columns
  colnames(idx) <- rep(names(groups_list), each = resample_to)  ## add back names alphabetically

  idx

}


mean_diff <- function(m, x) {
  ## average m along values of x.
  ## x: a dummy variable (0 and 1 mark classes)
  ## m: numeric vector
  stopifnot(unique(x) %in% 0:1)
  x <- sign(x - 0.5)
  x[x > 0] <- x[x > 0] / sum(x[x > 0])
  x[x < 0] <- x[x < 0] / -sum(x[x < 0])
  ## apply to m:
  x %*% m
}


opt_smoothing <- function(n_vox, n_trial, sd_error) {
  ## https://link.springer.com/article/10.1007/s10618-012-0295-5
  ## equation 3
  exponent <- 1 / (n_vox + 4)
  base <- 4 / ((n_vox + 2) * n_trial)
  sd_error * base^exponent
}


make_design <- function(frac_imb1, frac_imb2, n_trial1, n_trial2) {
  
  stopifnot(is_integer(frac_imb1*n_trial1))
  stopifnot(is_integer(frac_imb2*n_trial2))
  
  ## make table of counts
  counts <- data.frame(
    var1 = c("a", "a", "b", "b"),
    var2 = c("a", "b", "a", "b"),
    n = c(
      n_trial1 * frac_imb1, n_trial1 * (1 - frac_imb1),
      n_trial2 * frac_imb2, n_trial2 * (1 - frac_imb2)
    ) %>% round
  )
  
  ## expand to observations (one obs per row)
  expanded_df <- data.frame(var1 = character(), var2 = character())
  for (i in seq_len(nrow(counts))) {
    n_times <- counts$n[i]
    expanded_rows <- data.frame(var1 = rep(counts$var1[i], n_times), var2 = rep(counts$var2[i], n_times))
    expanded_df <- rbind(expanded_df, expanded_rows)
  }
  
  ## convert to numeric design matrix
  x <- model.matrix(~ var1*var2, expanded_df)
  colnames(x) <- c("intercept", "me1", "me2", "interaction")
  attr(x, "contrasts") <- NULL
  attr(x, "assign") <- NULL
  
  x
  
}


sim_data <- function(b_true, x, sd_error) {
  n_vox <- nrow(b_true)
  n_trial <- nrow(x)
  e <- matrix(rnorm(n_vox*n_trial, sd = sd_error), ncol = n_vox)
  y <- tcrossprod(x, b_true) + e
  y
}


cvlda <- function(x_train, x_test, y_train, y_test) {
  stopifnot(unique(x_test) %in% 0:1)  ## check if dummy
  #classes <- sign(x_test - 0.5)  ## convert to 1 or -1
  ldf_train <- mean_diff(y_train, x_train)
  projs <- tcrossprod(y_test, ldf_train)  ## project test data onto training discriminant
  stdev <- sqrt(var(projs[x_test == 0]) + var(projs[x_test == 1]))  ## unpooled within-class SD
  c(mean_diff(projs, x_test) / stdev)  ## cohen's D
}


cvlda_boot <- function(x_train, x_test, y_train, y_test, idxs, sd_resample) {
  
  n_resamples <- nrow(idxs)
  resample_to <- ncol(idxs)  ## number total obs
  stats_resamp <- numeric(n_resamples)
  for (i in seq_len(n_resamples)) {
    idx <- idxs[i, ]
    stats_resamp[i] <- cvlda(
      x_train = x_train[idx],
      x_test  = x_test,
      y_train = y_train[idx, ] + rnorm(resample_to*n_vox, sd = sd_resample),
      y_test  = y_test)
  }
  stat_resamp <- mean(stats_resamp)
  
  stat_resamp
  
}

