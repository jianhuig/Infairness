# Purpose: Bootstrap variance helpers for fairness estimation.

subset_threshold <- function(threshold, idx) {
  if (length(threshold) == 1) {
    threshold
  } else {
    threshold[idx]
  }
}

stratified_bootstrap_indices <- function(indices, A) {
  if (length(indices) == 0) {
    return(integer(0))
  }

  out <- integer(0)
  for (a in sort(unique(A[indices]))) {
    idx_a <- indices[A[indices] == a]
    out <- c(out, sample(idx_a, length(idx_a), replace = TRUE))
  }

  out
}

bootstrap_sample_indices <- function(Y, A, resample = "labeled") {
  resample <- match.arg(resample, c("labeled", "all"))

  labeled_ind <- which(!is.na(Y))
  unlabeled_ind <- which(is.na(Y))

  boot_labeled <- stratified_bootstrap_indices(labeled_ind, A)
  if (resample == "all") {
    boot_unlabeled <- stratified_bootstrap_indices(unlabeled_ind, A)
  } else {
    boot_unlabeled <- unlabeled_ind
  }

  c(boot_labeled, boot_unlabeled)
}

bootstrap_variance_from_estimates <- function(template_est, est_list) {
  value_cols <- setdiff(names(template_est), "Metric")
  out <- template_est
  if (length(est_list) == 0) {
    out[, value_cols] <- NA_real_
    return(out)
  }

  arr <- array(
    NA_real_,
    dim = c(nrow(template_est), length(value_cols), length(est_list))
  )

  for (b in seq_along(est_list)) {
    arr[, , b] <- as.matrix(est_list[[b]][, value_cols, drop = FALSE])
  }

  var_mat <- apply(
    arr,
    c(1, 2),
    function(x) {
      if (all(is.na(x))) {
        NA_real_
      } else {
        stats::var(x, na.rm = TRUE)
      }
    }
  )

  out[, value_cols] <- as.data.frame(var_mat)
  out
}

bootstrap_replicate_seeds <- function(B, seed = NULL) {
  if (is.null(seed)) {
    return(rep(NA_integer_, B))
  }

  base_seed <- as.integer(seed)[1]
  as.integer((base_seed + seq_len(B) - 1L) %% .Machine$integer.max)
}

bootstrap_apply <- function(B, bootstrap_cores = 1L, FUN) {
  bootstrap_cores <- max(1L, as.integer(bootstrap_cores)[1])

  if (bootstrap_cores > 1L && .Platform$OS.type != "windows") {
    return(parallel::mclapply(seq_len(B), FUN, mc.cores = bootstrap_cores))
  }

  lapply(seq_len(B), FUN)
}

bootstrap_variance_supervised <- function(template_est,
                                          Y,
                                          S,
                                          A,
                                          threshold = 0.5,
                                          B = 200,
                                          resample = "labeled",
                                          seed = NULL,
                                          bootstrap_cores = 1L) {
  replicate_seeds <- bootstrap_replicate_seeds(B = B, seed = seed)

  est_list <- bootstrap_apply(B = B, bootstrap_cores = bootstrap_cores, FUN = function(b) {
    if (!is.na(replicate_seeds[b])) {
      set.seed(replicate_seeds[b])
    }

    boot_idx <- bootstrap_sample_indices(Y, A, resample = resample)
    Y_b <- Y[boot_idx]
    S_b <- S[boot_idx]
    A_b <- A[boot_idx]
    threshold_b <- subset_threshold(threshold, boot_idx)

    get_metric(
      Y = Y_b[!is.na(Y_b)],
      S = S_b[!is.na(Y_b)],
      A = A_b[!is.na(Y_b)],
      threshold = threshold_b
    )
  })

  bootstrap_variance_from_estimates(template_est, est_list)
}

bootstrap_variance_semisupervised <- function(template_est,
                                              Y,
                                              S,
                                              A,
                                              threshold = 0.5,
                                              X = NULL,
                                              basis = c("Poly(S)", "Poly(S)"),
                                              nknots = 3,
                                              W = NULL,
                                              k = 10,
                                              cross_fit_variance = FALSE,
                                              control_variate = FALSE,
                                              spline_ridge = TRUE,
                                              spline_gam = FALSE,
                                              ridge_unpenalize_binary_X = FALSE,
                                              spline_include_D = TRUE,
                                              log_include_D = TRUE,
                                              kernel_order = 0.45,
                                              kernel_use_ecdf = FALSE,
                                              B = 200,
                                              resample = "labeled",
                                              seed = NULL,
                                              bootstrap_cores = 1L,
                                              ...) {
  replicate_seeds <- bootstrap_replicate_seeds(B = B, seed = seed)

  est_list <- bootstrap_apply(B = B, bootstrap_cores = bootstrap_cores, FUN = function(b) {
    if (!is.na(replicate_seeds[b])) {
      set.seed(replicate_seeds[b])
    }

    boot_idx <- bootstrap_sample_indices(Y, A, resample = resample)
    Y_b <- Y[boot_idx]
    S_b <- S[boot_idx]
    A_b <- A[boot_idx]
    X_b <- if (is.null(X)) NULL else X[boot_idx, , drop = FALSE]
    W_b <- if (is.null(W) || length(W) == 1) W else W[boot_idx]
    threshold_b <- subset_threshold(threshold, boot_idx)

    SSFairness(
      Y = Y_b,
      S = S_b,
      A = A_b,
      threshold = threshold_b,
      X = X_b,
      basis = basis,
      nknots = nknots,
      W = W_b,
      k = k,
      folds = NULL,
      cross_fit_variance = cross_fit_variance,
      return_imputation_quality = FALSE,
      control_variate = control_variate,
      variance_method = "if",
      spline_ridge = spline_ridge,
      spline_gam = spline_gam,
      ridge_unpenalize_binary_X = ridge_unpenalize_binary_X,
      spline_include_D = spline_include_D,
      log_include_D = log_include_D,
      kernel_order = kernel_order,
      kernel_use_ecdf = kernel_use_ecdf,
      ...
    )$est
  })

  bootstrap_variance_from_estimates(template_est, est_list)
}
