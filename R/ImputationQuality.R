# Purpose: Imputation-quality helpers for semi-supervised estimation.

make_group_folds <- function(Y_labeled, A_labeled, k) {
  nclass <- sort(unique(A_labeled))
  folds <- setNames(vector("list", length(nclass)), as.character(nclass))

  for (a in nclass) {
    Y_a <- Y_labeled[A_labeled == a]
    k_a <- min(k, length(Y_a))

    if (k_a < 2) {
      folds[[as.character(a)]] <- rep.int(1L, length(Y_a))
      next
    }

    if (length(unique(Y_a)) > 1 && all(table(Y_a) >= 2)) {
      fold_a <- caret::createFolds(factor(Y_a), k = k_a, list = FALSE)
    } else {
      fold_a <- sample(rep(seq_len(k_a), length.out = length(Y_a)))
    }

    folds[[as.character(a)]] <- as.integer(fold_a)
  }

  folds
}

validate_group_folds <- function(folds, Y_labeled, A_labeled) {
  nclass <- sort(unique(A_labeled))

  if (!is.list(folds) || is.null(names(folds))) {
    stop("`folds` must be a named list keyed by group values.")
  }

  for (a in nclass) {
    key <- as.character(a)
    if (is.null(folds[[key]])) {
      stop("Missing fold assignments for group ", a, ".")
    }

    expected_n <- sum(A_labeled == a)
    if (length(folds[[key]]) != expected_n) {
      stop(
        "Fold assignments for group ", a,
        " must have length ", expected_n, "."
      )
    }
  }

  folds
}

make_group_gbic_splits <- function(A_labeled, folds = NULL) {
  nclass <- sort(unique(A_labeled))
  gbic_splits <- setNames(vector("list", length(nclass)), as.character(nclass))

  for (a in nclass) {
    n_a <- sum(A_labeled == a)
    if (!is.null(folds) && !is.null(folds[[as.character(a)]])) {
      fold_a <- folds[[as.character(a)]]
      gbic_splits[[as.character(a)]] <- order(fold_a, seq_along(fold_a))
    } else {
      gbic_splits[[as.character(a)]] <- seq_len(n_a)
    }
  }

  gbic_splits
}

compute_candidate_bic <- function(Y_a, design_matrix) {
  compute_candidate_bic_with_penalty(Y_a, design_matrix, penalty_rate_fun = log)
}

compute_candidate_mbic <- function(Y_a, design_matrix) {
  compute_candidate_bic_with_penalty(
    Y_a,
    design_matrix,
    penalty_rate_fun = function(n) min(n^0.1, log(n))
  )
}

compute_candidate_bic_with_penalty <- function(Y_a, design_matrix, penalty_rate_fun) {
  X <- as.matrix(design_matrix)
  y <- as.numeric(Y_a)
  n <- length(y)

  if (nrow(X) != n || n == 0L || ncol(X) == 0L) {
    return(Inf)
  }

  fit <- tryCatch(
    suppressWarnings(stats::glm.fit(x = X, y = y, family = stats::binomial())),
    error = function(e) NULL
  )
  if (is.null(fit) || !isTRUE(fit$converged)) {
    return(Inf)
  }

  beta <- fit$coefficients
  beta[is.na(beta)] <- 0
  eta <- drop(X %*% beta)
  if (any(!is.finite(eta))) {
    return(Inf)
  }

  mu <- pmin(pmax(stats::plogis(eta), .Machine$double.eps), 1 - .Machine$double.eps)
  loglik <- sum(y * log(mu) + (1 - y) * log(1 - mu))
  bic <- -2 * loglik + fit$rank * penalty_rate_fun(n)

  if (is.finite(bic)) bic else Inf
}

compute_candidate_gbic <- function(Y_a, design_matrix, split_idx) {
  # Kept as an internal compatibility wrapper for existing candidate blocks.
  compute_candidate_bic(Y_a, design_matrix)
}

compute_imputation_quality <- function(Y,
                                       S,
                                       A,
                                       threshold = 0.5,
                                       X = NULL,
                                       basis = c("Poly(S)", "Poly(S)"),
                                       nknots = 3,
                                       k = 10,
                                       folds = NULL,
                                       cross_fit = FALSE,
                                       control_variate = FALSE,
                                       variance_method = "if",
                                       bootstrap_reps = 200,
                                       bootstrap_resample = "labeled",
                                       bootstrap_seed = NULL,
                                       bootstrap_cores = 1,
                                       spline_ridge = TRUE,
                                       spline_gam = FALSE,
                                       ridge_unpenalize_binary_X = FALSE,
                                       spline_include_D = TRUE,
                                       log_include_D = TRUE,
                                       kernel_order = 0.45,
                                       kernel_use_ecdf = FALSE,
                                       alphas = NULL,
                                       ...) {
  variance_method <- match.arg(variance_method, c("if", "bootstrap"))
  labeled_ind <- which(!is.na(Y))

  # Create Indicator
  C <- ifelse(S > threshold, 1, 0)

  # Labeled data.
  Y_labeled <- Y[labeled_ind]
  S_labeled <- S[labeled_ind]
  A_labeled <- A[labeled_ind]
  C_labeled <- C[labeled_ind]


  # Unlabeled data.
  S_unlabeled <- S[-labeled_ind]
  A_unlabeled <- A[-labeled_ind]
  C_unlabeled <- C[-labeled_ind]

  if (!is.null(X)) {
    X_labeled <- X[labeled_ind, , drop = FALSE]
    X_unlabeled <- X[-labeled_ind, , drop = FALSE]
  }

  nclass <- sort(unique(A))
  if (length(basis) == 1) {
    basis <- rep(basis, length(nclass))
  }

  if (length(threshold) == 1) {
    threshold_raw <- rep(threshold, length(S))
  } else if (length(threshold) == length(S)) {
    threshold_raw <- threshold
  } else {
    stop("`threshold` must have length 1 or length equal to `S`.")
  }

  kernel_scale <- kernel_metric_scale(
    S = S,
    A = A,
    basis = basis,
    threshold = threshold,
    use_ecdf = kernel_use_ecdf
  )
  S_metric <- kernel_scale$S
  threshold_metric <- kernel_scale$threshold
  S_labeled_metric <- S_metric[labeled_ind]
  S_unlabeled_metric <- S_metric[-labeled_ind]
  threshold_labeled_metric <- threshold_metric[labeled_ind]
  threshold_unlabeled_metric <- threshold_metric[-labeled_ind]
  threshold_labeled_raw <- threshold_raw[labeled_ind]
  threshold_unlabeled_raw <- threshold_raw[-labeled_ind]
  C_metric <- ifelse(S_metric > threshold_metric, 1, 0)
  C_labeled_metric <- C_metric[labeled_ind]

  if (cross_fit) {
    if (is.null(folds)) {
      folds <- make_group_folds(Y_labeled, A_labeled, k)
    } else {
      folds <- validate_group_folds(folds, Y_labeled, A_labeled)
    }
  }
  gbic_splits <- make_group_gbic_splits(A_labeled, folds)

  m_unlabeled <- S_unlabeled
  m_labeled <- S_labeled
  alphas <- c()
  gbic_scores <- c()
  mbic_scores <- c()

  for (a in nclass) {
    basis_a <- basis[which(nclass == a)]

    if (basis_a == "Poly(S)") {
      alpha <- tryCatch(
        {
          find_alpha_glm(
            Y = Y_labeled[A_labeled == a],
            covariates_matrix = S_labeled[A_labeled == a] %>% as.matrix(),
            additional_matrix = C_labeled[A_labeled == a] %>% as.matrix()
          )
        },
        error = function(e) 1
      )
      alphas <- c(alphas, alpha)

      basis_labeled <- polynomial(S_labeled %>% as.data.frame(), alpha)
      basis_unlabeled <- polynomial(S_unlabeled %>% as.data.frame(), alpha)

      gamma <- tryCatch(
        {
          RidgeRegression(
            X = cbind(
              basis_labeled[A_labeled == a, -1],
              C_labeled[A_labeled == a]
            ) %>% as.matrix(),
            y = Y_labeled[A_labeled == a]
          )
        },
        error = function(e) {
          print("Ridge Regression produced an error")
          print(e)
        }
      )

      m_unlabeled[A_unlabeled == a] <- boot::inv.logit(
        as.matrix(
          cbind(
            1,
            basis_unlabeled[A_unlabeled == a, -1],
            C_unlabeled[A_unlabeled == a]
          )
        ) %*% gamma$coefficients
      )

      m_labeled[A_labeled == a] <- boot::inv.logit(
        as.matrix(
          cbind(
            1,
            basis_labeled[A_labeled == a, -1],
            C_labeled[A_labeled == a]
          )
        ) %*% gamma$coefficients
      )

      if (cross_fit) {
        Y_a <- Y_labeled[A_labeled == a]
        S_a <- S_labeled[A_labeled == a]
        C_a <- C_labeled[A_labeled == a]
        fold <- folds[[as.character(a)]]

        for (i in sort(unique(fold))) {
          train_id <- which(fold != i)
          test_id <- which(fold == i)
          alpha_fold <- alphas[a + 1]
          basis_train <- polynomial(S_a[train_id] %>% as.data.frame(), alpha_fold)
          basis_test <- polynomial(S_a[test_id] %>% as.data.frame(), alpha_fold)

          gamma <- tryCatch(
            {
            RidgeRegression(
              X = cbind(basis_train[, -1], C_a[train_id]) %>% as.matrix(),
              y = Y_a[train_id]
            )
          },
          error = function(e) {
            print("Ridge Regression produced an error")
            print(e)
          }
        )$coefficients

          m_labeled[which(A_labeled == a)[test_id]] <- boot::inv.logit(
            as.matrix(cbind(1, basis_test[, -1], C_a[test_id])) %*% gamma
          )
        }
      }

      design_gbic <- as.matrix(
        cbind(
          1,
          basis_labeled[A_labeled == a, -1],
          C_labeled[A_labeled == a]
        )
      )
      gbic_scores <- c(
        gbic_scores,
        compute_candidate_gbic(
          Y_labeled[A_labeled == a],
          design_gbic,
          gbic_splits[[as.character(a)]]
        )
      )
      mbic_scores <- c(
        mbic_scores,
        compute_candidate_mbic(
          Y_labeled[A_labeled == a],
          design_gbic
        )
      )
    } else if (basis_a == "Poly(S) + X") {
      alpha <- tryCatch(
        {
          find_alpha_glm(
            Y = Y_labeled[A_labeled == a],
            covariates_matrix = S_labeled[A_labeled == a] %>% as.matrix(),
            additional_matrix = cbind(
              C_labeled[A_labeled == a],
              X_labeled[A_labeled == a, ]
            ) %>% as.matrix()
          )
        },
        error = function(e) 1
      )
      alphas <- c(alphas, alpha)

      basis_labeled <- polynomial(S_labeled %>% as.data.frame(), alpha)
      basis_unlabeled <- polynomial(S_unlabeled %>% as.data.frame(), alpha)

      gamma <- tryCatch(
        {
          RidgeRegression(
            X = cbind(
              basis_labeled[A_labeled == a, -1],
              C_labeled[A_labeled == a],
              X_labeled[A_labeled == a, ]
            ) %>% as.matrix(),
            y = Y_labeled[A_labeled == a],
            penalty_factor = ridge_design_penalty_factor(
              n_basis = ncol(basis_labeled[, -1, drop = FALSE]),
              D = C_labeled[A_labeled == a],
              X = X_labeled[A_labeled == a, , drop = FALSE],
              unpenalize_binary_X = ridge_unpenalize_binary_X
            )
          )
        },
        error = function(e) {
          print("Ridge Regression produced an error")
          print(e)
        }
      )

      m_unlabeled[A_unlabeled == a] <- boot::inv.logit(
        as.matrix(
          cbind(
            1,
            basis_unlabeled[A_unlabeled == a, -1],
            C_unlabeled[A_unlabeled == a],
            X_unlabeled[A_unlabeled == a, ]
          )
        ) %*% gamma$coefficients
      )

      m_labeled[A_labeled == a] <- boot::inv.logit(
        as.matrix(
          cbind(
            1,
            basis_labeled[A_labeled == a, -1],
            C_labeled[A_labeled == a],
            X_labeled[A_labeled == a, ]
          )
        ) %*% gamma$coefficients
      )

      if (cross_fit) {
        Y_a <- Y_labeled[A_labeled == a]
        S_a <- S_labeled[A_labeled == a]
        X_a <- X_labeled[A_labeled == a, , drop = FALSE]
        C_a <- C_labeled[A_labeled == a]
        fold <- folds[[as.character(a)]]

        for (i in sort(unique(fold))) {
          train_id <- which(fold != i)
          test_id <- which(fold == i)
          alpha_fold <- alphas[a + 1]
          basis_train <- polynomial(S_a[train_id] %>% as.data.frame(), alpha_fold)
          basis_test <- polynomial(S_a[test_id] %>% as.data.frame(), alpha_fold)

          gamma <- tryCatch(
            {
            RidgeRegression(
              X = cbind(
                basis_train[, -1],
                C_a[train_id],
                X_a[train_id, ]
              ) %>% as.matrix(),
              y = Y_a[train_id],
              penalty_factor = ridge_design_penalty_factor(
                n_basis = ncol(basis_train[, -1, drop = FALSE]),
                D = C_a[train_id],
                X = X_a[train_id, , drop = FALSE],
                unpenalize_binary_X = ridge_unpenalize_binary_X
              )
            )
          },
          error = function(e) {
            print("Ridge Regression produced an error")
            print(e)
          }
        )$coefficients

          m_labeled[which(A_labeled == a)[test_id]] <- boot::inv.logit(
            as.matrix(
              cbind(1, basis_test[, -1], C_a[test_id], X_a[test_id, ])
            ) %*% gamma
          )
        }
      }

      design_gbic <- as.matrix(
        cbind(
          1,
          basis_labeled[A_labeled == a, -1],
          C_labeled[A_labeled == a],
          X_labeled[A_labeled == a, ]
        )
      )
      gbic_scores <- c(
        gbic_scores,
        compute_candidate_gbic(
          Y_labeled[A_labeled == a],
          design_gbic,
          gbic_splits[[as.character(a)]]
        )
      )
      mbic_scores <- c(
        mbic_scores,
        compute_candidate_mbic(
          Y_labeled[A_labeled == a],
          design_gbic
        )
      )
    } else if (basis_a == "Spline(S)") {
      alphas <- c(alphas, nknots)

      basis_labeled <- NaturalSplineBasis(S_labeled %>% as.matrix(), nknots, return_knots = TRUE)
      spline_knots <- attr(basis_labeled, "knots")
      basis_unlabeled <- NaturalSplineBasis(
        S_unlabeled %>% as.matrix(),
        nknots,
        knots = spline_knots
      )

      predictors_labeled <- build_spline_predictors(
        basis_matrix = basis_labeled[A_labeled == a, , drop = FALSE],
        D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D)
      )
      predictors_unlabeled <- build_spline_predictors(
        basis_matrix = basis_unlabeled[A_unlabeled == a, , drop = FALSE],
        D = spline_indicator_term(C_unlabeled[A_unlabeled == a], spline_include_D)
      )

      gamma <- tryCatch(
        {
          fit_spline_model(
            X = predictors_labeled,
            y = Y_labeled[A_labeled == a],
            use_ridge = spline_ridge,
            use_gam = spline_gam,
            S = S_labeled[A_labeled == a],
            D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D),
            num_knots = nknots,
            penalty_factor = build_spline_penalty_factor(
              basis_matrix = basis_labeled[A_labeled == a, , drop = FALSE],
              D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D),
              unpenalize_binary_X = ridge_unpenalize_binary_X
            )
          )
        },
        error = function(e) {
          print("Spline model fitting produced an error")
          print(e)
        }
      )

      m_unlabeled[A_unlabeled == a] <- predict_spline_model(
        gamma,
        predictors_unlabeled,
        S = S_unlabeled[A_unlabeled == a],
        D = spline_indicator_term(C_unlabeled[A_unlabeled == a], spline_include_D)
      )

      m_labeled[A_labeled == a] <- predict_spline_model(
        gamma,
        predictors_labeled,
        S = S_labeled[A_labeled == a],
        D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D)
      )

      if (cross_fit) {
        Y_a <- Y_labeled[A_labeled == a]
        S_a <- S_labeled[A_labeled == a]
        C_a <- C_labeled[A_labeled == a]
        fold <- folds[[as.character(a)]]

        for (i in sort(unique(fold))) {
          train_id <- which(fold != i)
          test_id <- which(fold == i)
          basis_train <- NaturalSplineBasis(S_a[train_id] %>% as.matrix(), nknots, return_knots = TRUE)
          fold_knots <- attr(basis_train, "knots")
          basis_test <- NaturalSplineBasis(
            S_a[test_id] %>% as.matrix(),
            nknots,
            knots = fold_knots
          )

          predictors_train <- build_spline_predictors(
            basis_matrix = basis_train,
            D = spline_indicator_term(C_a[train_id], spline_include_D)
          )
          predictors_test <- build_spline_predictors(
            basis_matrix = basis_test,
            D = spline_indicator_term(C_a[test_id], spline_include_D)
          )

          gamma <- tryCatch(
            {
              fit_spline_model(
                X = predictors_train,
                y = Y_a[train_id],
                use_ridge = spline_ridge,
                use_gam = spline_gam,
                S = S_a[train_id],
                D = spline_indicator_term(C_a[train_id], spline_include_D),
                num_knots = nknots,
                penalty_factor = build_spline_penalty_factor(
                  basis_matrix = basis_train,
                  D = spline_indicator_term(C_a[train_id], spline_include_D),
                  unpenalize_binary_X = ridge_unpenalize_binary_X
                )
              )
            },
            error = function(e) {
              print("Spline model fitting produced an error")
              print(e)
            }
          )

          m_labeled[which(A_labeled == a)[test_id]] <- predict_spline_model(
            gamma,
            predictors_test,
            S = S_a[test_id],
            D = spline_indicator_term(C_a[test_id], spline_include_D)
          )
        }
      }

      design_gbic <- cbind(1, predictors_labeled)
      gbic_scores <- c(
        gbic_scores,
        compute_candidate_gbic(
          Y_labeled[A_labeled == a],
          design_gbic,
          gbic_splits[[as.character(a)]]
        )
      )
      mbic_scores <- c(
        mbic_scores,
        compute_candidate_mbic(
          Y_labeled[A_labeled == a],
          design_gbic
        )
      )
    } else if (basis_a == "Spline(S) + X") {
      alphas <- c(alphas, nknots)

      basis_labeled <- NaturalSplineBasis(S_labeled %>% as.matrix(), nknots, return_knots = TRUE)
      spline_knots <- attr(basis_labeled, "knots")
      basis_unlabeled <- NaturalSplineBasis(
        S_unlabeled %>% as.matrix(),
        nknots,
        knots = spline_knots
      )

      predictors_labeled <- build_spline_predictors(
        basis_matrix = basis_labeled[A_labeled == a, , drop = FALSE],
        D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D),
        X = X_labeled[A_labeled == a, , drop = FALSE]
      )
      predictors_unlabeled <- build_spline_predictors(
        basis_matrix = basis_unlabeled[A_unlabeled == a, , drop = FALSE],
        D = spline_indicator_term(C_unlabeled[A_unlabeled == a], spline_include_D),
        X = X_unlabeled[A_unlabeled == a, , drop = FALSE]
      )

      gamma <- tryCatch(
        {
          fit_spline_model(
            X = predictors_labeled,
            y = Y_labeled[A_labeled == a],
            use_ridge = spline_ridge,
            use_gam = spline_gam,
            S = S_labeled[A_labeled == a],
            D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D),
            X_cov = X_labeled[A_labeled == a, , drop = FALSE],
            num_knots = nknots,
            penalty_factor = build_spline_penalty_factor(
              basis_matrix = basis_labeled[A_labeled == a, , drop = FALSE],
              D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D),
              X = X_labeled[A_labeled == a, , drop = FALSE],
              unpenalize_binary_X = ridge_unpenalize_binary_X
            )
          )
        },
        error = function(e) {
          print("Spline model fitting produced an error")
          print(e)
        }
      )

      m_unlabeled[A_unlabeled == a] <- predict_spline_model(
        gamma,
        predictors_unlabeled,
        S = S_unlabeled[A_unlabeled == a],
        D = spline_indicator_term(C_unlabeled[A_unlabeled == a], spline_include_D),
        X_cov = X_unlabeled[A_unlabeled == a, , drop = FALSE]
      )

      m_labeled[A_labeled == a] <- predict_spline_model(
        gamma,
        predictors_labeled,
        S = S_labeled[A_labeled == a],
        D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D),
        X_cov = X_labeled[A_labeled == a, , drop = FALSE]
      )

      if (cross_fit) {
        Y_a <- Y_labeled[A_labeled == a]
        S_a <- S_labeled[A_labeled == a]
        X_a <- X_labeled[A_labeled == a, , drop = FALSE]
        C_a <- C_labeled[A_labeled == a]
        fold <- folds[[as.character(a)]]

        for (i in sort(unique(fold))) {
          train_id <- which(fold != i)
          test_id <- which(fold == i)
          basis_train <- NaturalSplineBasis(S_a[train_id] %>% as.matrix(), nknots, return_knots = TRUE)
          fold_knots <- attr(basis_train, "knots")
          basis_test <- NaturalSplineBasis(
            S_a[test_id] %>% as.matrix(),
            nknots,
            knots = fold_knots
          )

          predictors_train <- build_spline_predictors(
            basis_matrix = basis_train,
            D = spline_indicator_term(C_a[train_id], spline_include_D),
            X = X_a[train_id, , drop = FALSE]
          )
          predictors_test <- build_spline_predictors(
            basis_matrix = basis_test,
            D = spline_indicator_term(C_a[test_id], spline_include_D),
            X = X_a[test_id, , drop = FALSE]
          )

          gamma <- tryCatch(
            {
              fit_spline_model(
                X = predictors_train,
                y = Y_a[train_id],
                use_ridge = spline_ridge,
                use_gam = spline_gam,
                S = S_a[train_id],
                D = spline_indicator_term(C_a[train_id], spline_include_D),
                X_cov = X_a[train_id, , drop = FALSE],
                num_knots = nknots,
                penalty_factor = build_spline_penalty_factor(
                  basis_matrix = basis_train,
                  D = spline_indicator_term(C_a[train_id], spline_include_D),
                  X = X_a[train_id, , drop = FALSE],
                  unpenalize_binary_X = ridge_unpenalize_binary_X
                )
              )
            },
            error = function(e) {
              print("Spline model fitting produced an error")
              print(e)
            }
          )

          m_labeled[which(A_labeled == a)[test_id]] <- predict_spline_model(
            gamma,
            predictors_test,
            S = S_a[test_id],
            D = spline_indicator_term(C_a[test_id], spline_include_D),
            X_cov = X_a[test_id, , drop = FALSE]
          )
        }
      }

      design_gbic <- cbind(1, predictors_labeled)
      gbic_scores <- c(
        gbic_scores,
        compute_candidate_gbic(
          Y_labeled[A_labeled == a],
          design_gbic,
          gbic_splits[[as.character(a)]]
        )
      )
      mbic_scores <- c(
        mbic_scores,
        compute_candidate_mbic(
          Y_labeled[A_labeled == a],
          design_gbic
        )
      )
    } else if (basis_a == "Spline Interaction") {
      alphas <- c(alphas, nknots)

      basis_labeled <- NaturalSplineBasis(S_labeled %>% as.matrix(), nknots, return_knots = TRUE)
      spline_knots <- attr(basis_labeled, "knots")
      basis_unlabeled <- NaturalSplineBasis(
        S_unlabeled %>% as.matrix(),
        nknots,
        knots = spline_knots
      )

      X_labeled_spline_int <- do.call(
        cbind,
        lapply(seq_len(ncol(X_labeled)), function(j) {
          basis_labeled * X_labeled[, j]
        })
      )
      X_unlabeled_spline_int <- do.call(
        cbind,
        lapply(seq_len(ncol(X_unlabeled)), function(j) {
          basis_unlabeled * X_unlabeled[, j]
        })
      )

      predictors_labeled <- build_spline_predictors(
        basis_matrix = basis_labeled[A_labeled == a, , drop = FALSE],
        D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D),
        X = X_labeled[A_labeled == a, , drop = FALSE],
        X_interaction = X_labeled_spline_int[A_labeled == a, , drop = FALSE]
      )
      predictors_unlabeled <- build_spline_predictors(
        basis_matrix = basis_unlabeled[A_unlabeled == a, , drop = FALSE],
        D = spline_indicator_term(C_unlabeled[A_unlabeled == a], spline_include_D),
        X = X_unlabeled[A_unlabeled == a, , drop = FALSE],
        X_interaction = X_unlabeled_spline_int[A_unlabeled == a, , drop = FALSE]
      )

      gamma <- tryCatch(
        {
          fit_spline_model(
            X = predictors_labeled,
            y = Y_labeled[A_labeled == a],
            use_ridge = spline_ridge,
            use_gam = spline_gam,
            S = S_labeled[A_labeled == a],
            D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D),
            X_cov = X_labeled[A_labeled == a, , drop = FALSE],
            interaction = TRUE,
            num_knots = nknots,
            penalty_factor = build_spline_penalty_factor(
              basis_matrix = basis_labeled[A_labeled == a, , drop = FALSE],
              D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D),
              X = X_labeled[A_labeled == a, , drop = FALSE],
              X_interaction = X_labeled_spline_int[A_labeled == a, , drop = FALSE],
              unpenalize_binary_X = ridge_unpenalize_binary_X
            )
          )
        },
        error = function(e) {
          print("Spline model fitting produced an error")
          print(e)
        }
      )

      m_unlabeled[A_unlabeled == a] <- predict_spline_model(
        gamma,
        predictors_unlabeled,
        S = S_unlabeled[A_unlabeled == a],
        D = spline_indicator_term(C_unlabeled[A_unlabeled == a], spline_include_D),
        X_cov = X_unlabeled[A_unlabeled == a, , drop = FALSE]
      )

      m_labeled[A_labeled == a] <- predict_spline_model(
        gamma,
        predictors_labeled,
        S = S_labeled[A_labeled == a],
        D = spline_indicator_term(C_labeled[A_labeled == a], spline_include_D),
        X_cov = X_labeled[A_labeled == a, , drop = FALSE]
      )

      if (cross_fit) {
        Y_a <- Y_labeled[A_labeled == a]
        S_a <- S_labeled[A_labeled == a]
        X_a <- X_labeled[A_labeled == a, , drop = FALSE]
        C_a <- C_labeled[A_labeled == a]
        fold <- folds[[as.character(a)]]

        for (i in sort(unique(fold))) {
          train_id <- which(fold != i)
          test_id <- which(fold == i)
          basis_train <- NaturalSplineBasis(
            S_a[train_id] %>% as.matrix(),
            nknots,
            return_knots = TRUE
          )
          fold_knots <- attr(basis_train, "knots")
          basis_test <- NaturalSplineBasis(
            S_a[test_id] %>% as.matrix(),
            nknots,
            knots = fold_knots
          )

          X_train_int <- do.call(
            cbind,
            lapply(seq_len(ncol(X_a)), function(j) {
              basis_train * X_a[train_id, j]
            })
          )
          X_test_int <- do.call(
            cbind,
            lapply(seq_len(ncol(X_a)), function(j) {
              basis_test * X_a[test_id, j]
            })
          )

          predictors_train <- build_spline_predictors(
            basis_matrix = basis_train,
            D = spline_indicator_term(C_a[train_id], spline_include_D),
            X = X_a[train_id, , drop = FALSE],
            X_interaction = X_train_int
          )
          predictors_test <- build_spline_predictors(
            basis_matrix = basis_test,
            D = spline_indicator_term(C_a[test_id], spline_include_D),
            X = X_a[test_id, , drop = FALSE],
            X_interaction = X_test_int
          )

          gamma <- tryCatch(
            {
              fit_spline_model(
                X = predictors_train,
                y = Y_a[train_id],
                use_ridge = spline_ridge,
                use_gam = spline_gam,
                S = S_a[train_id],
                D = spline_indicator_term(C_a[train_id], spline_include_D),
                X_cov = X_a[train_id, , drop = FALSE],
                interaction = TRUE,
                num_knots = nknots,
                penalty_factor = build_spline_penalty_factor(
                  basis_matrix = basis_train,
                  D = spline_indicator_term(C_a[train_id], spline_include_D),
                  X = X_a[train_id, , drop = FALSE],
                  X_interaction = X_train_int,
                  unpenalize_binary_X = ridge_unpenalize_binary_X
                )
              )
            },
            error = function(e) {
              print("Spline model fitting produced an error")
              print(e)
            }
          )

          m_labeled[which(A_labeled == a)[test_id]] <- predict_spline_model(
            gamma,
            predictors_test,
            S = S_a[test_id],
            D = spline_indicator_term(C_a[test_id], spline_include_D),
            X_cov = X_a[test_id, , drop = FALSE]
          )
        }
      }

      design_gbic <- cbind(1, predictors_labeled)
      gbic_scores <- c(
        gbic_scores,
        compute_candidate_gbic(
          Y_labeled[A_labeled == a],
          design_gbic,
          gbic_splits[[as.character(a)]]
        )
      )
      mbic_scores <- c(
        mbic_scores,
        compute_candidate_mbic(
          Y_labeled[A_labeled == a],
          design_gbic
        )
      )
    } else if (basis_a == "Log(S)") {
      alphas <- c(alphas, NA_real_)

      predictors_labeled <- build_log_predictors(
        S = S_labeled[A_labeled == a],
        D = if (log_include_D) C_labeled[A_labeled == a] else NULL
      )
      predictors_unlabeled <- build_log_predictors(
        S = S_unlabeled[A_unlabeled == a],
        D = if (log_include_D) C_unlabeled[A_unlabeled == a] else NULL
      )

      gamma <- tryCatch(
        {
          fit_log_basis_model(
            X = predictors_labeled,
            y = Y_labeled[A_labeled == a]
          )
        },
        error = function(e) {
          print("Log-basis model fitting produced an error")
          print(e)
        }
      )

      m_unlabeled[A_unlabeled == a] <- predict_log_basis_model(
        gamma,
        predictors_unlabeled
      )
      m_labeled[A_labeled == a] <- predict_log_basis_model(
        gamma,
        predictors_labeled
      )

      if (cross_fit) {
        Y_a <- Y_labeled[A_labeled == a]
        S_a <- S_labeled[A_labeled == a]
        C_a <- C_labeled[A_labeled == a]
        fold <- folds[[as.character(a)]]

        for (i in sort(unique(fold))) {
          train_id <- which(fold != i)
          test_id <- which(fold == i)

          predictors_train <- build_log_predictors(
            S = S_a[train_id],
            D = if (log_include_D) C_a[train_id] else NULL
          )
          predictors_test <- build_log_predictors(
            S = S_a[test_id],
            D = if (log_include_D) C_a[test_id] else NULL
          )

          gamma <- tryCatch(
            {
              fit_log_basis_model(
                X = predictors_train,
                y = Y_a[train_id]
              )
            },
            error = function(e) {
              print("Log-basis model fitting produced an error")
              print(e)
            }
          )

          m_labeled[which(A_labeled == a)[test_id]] <- predict_log_basis_model(
            gamma,
            predictors_test
          )
        }
      }

      design_gbic <- cbind(1, predictors_labeled)
      gbic_scores <- c(
        gbic_scores,
        compute_candidate_gbic(
          Y_labeled[A_labeled == a],
          design_gbic,
          gbic_splits[[as.character(a)]]
        )
      )
      mbic_scores <- c(
        mbic_scores,
        compute_candidate_mbic(
          Y_labeled[A_labeled == a],
          design_gbic
        )
      )
    } else if (basis_a == "Log(S) + X") {
      alphas <- c(alphas, NA_real_)

      predictors_labeled <- build_log_predictors(
        S = S_labeled[A_labeled == a],
        D = if (log_include_D) C_labeled[A_labeled == a] else NULL,
        X = X_labeled[A_labeled == a, , drop = FALSE]
      )
      predictors_unlabeled <- build_log_predictors(
        S = S_unlabeled[A_unlabeled == a],
        D = if (log_include_D) C_unlabeled[A_unlabeled == a] else NULL,
        X = X_unlabeled[A_unlabeled == a, , drop = FALSE]
      )

      gamma <- tryCatch(
        {
          fit_log_basis_model(
            X = predictors_labeled,
            y = Y_labeled[A_labeled == a]
          )
        },
        error = function(e) {
          print("Log-basis model fitting produced an error")
          print(e)
        }
      )

      m_unlabeled[A_unlabeled == a] <- predict_log_basis_model(
        gamma,
        predictors_unlabeled
      )
      m_labeled[A_labeled == a] <- predict_log_basis_model(
        gamma,
        predictors_labeled
      )

      if (cross_fit) {
        Y_a <- Y_labeled[A_labeled == a]
        S_a <- S_labeled[A_labeled == a]
        X_a <- X_labeled[A_labeled == a, , drop = FALSE]
        C_a <- C_labeled[A_labeled == a]
        fold <- folds[[as.character(a)]]

        for (i in sort(unique(fold))) {
          train_id <- which(fold != i)
          test_id <- which(fold == i)

          predictors_train <- build_log_predictors(
            S = S_a[train_id],
            D = if (log_include_D) C_a[train_id] else NULL,
            X = X_a[train_id, , drop = FALSE]
          )
          predictors_test <- build_log_predictors(
            S = S_a[test_id],
            D = if (log_include_D) C_a[test_id] else NULL,
            X = X_a[test_id, , drop = FALSE]
          )

          gamma <- tryCatch(
            {
              fit_log_basis_model(
                X = predictors_train,
                y = Y_a[train_id]
              )
            },
            error = function(e) {
              print("Log-basis model fitting produced an error")
              print(e)
            }
          )

          m_labeled[which(A_labeled == a)[test_id]] <- predict_log_basis_model(
            gamma,
            predictors_test
          )
        }
      }

      design_gbic <- cbind(1, predictors_labeled)
      gbic_scores <- c(
        gbic_scores,
        compute_candidate_gbic(
          Y_labeled[A_labeled == a],
          design_gbic,
          gbic_splits[[as.character(a)]]
        )
      )
      mbic_scores <- c(
        mbic_scores,
        compute_candidate_mbic(
          Y_labeled[A_labeled == a],
          design_gbic
        )
      )
    } else if (basis_a == "Log Interaction") {
      alphas <- c(alphas, NA_real_)

      predictors_labeled <- build_log_predictors(
        S = S_labeled[A_labeled == a],
        D = if (log_include_D) C_labeled[A_labeled == a] else NULL,
        X = X_labeled[A_labeled == a, , drop = FALSE],
        interaction = TRUE
      )
      predictors_unlabeled <- build_log_predictors(
        S = S_unlabeled[A_unlabeled == a],
        D = if (log_include_D) C_unlabeled[A_unlabeled == a] else NULL,
        X = X_unlabeled[A_unlabeled == a, , drop = FALSE],
        interaction = TRUE
      )

      gamma <- tryCatch(
        {
          fit_log_basis_model(
            X = predictors_labeled,
            y = Y_labeled[A_labeled == a]
          )
        },
        error = function(e) {
          print("Log-basis model fitting produced an error")
          print(e)
        }
      )

      m_unlabeled[A_unlabeled == a] <- predict_log_basis_model(
        gamma,
        predictors_unlabeled
      )
      m_labeled[A_labeled == a] <- predict_log_basis_model(
        gamma,
        predictors_labeled
      )

      if (cross_fit) {
        Y_a <- Y_labeled[A_labeled == a]
        S_a <- S_labeled[A_labeled == a]
        X_a <- X_labeled[A_labeled == a, , drop = FALSE]
        C_a <- C_labeled[A_labeled == a]
        fold <- folds[[as.character(a)]]

        for (i in sort(unique(fold))) {
          train_id <- which(fold != i)
          test_id <- which(fold == i)

          predictors_train <- build_log_predictors(
            S = S_a[train_id],
            D = if (log_include_D) C_a[train_id] else NULL,
            X = X_a[train_id, , drop = FALSE],
            interaction = TRUE
          )
          predictors_test <- build_log_predictors(
            S = S_a[test_id],
            D = if (log_include_D) C_a[test_id] else NULL,
            X = X_a[test_id, , drop = FALSE],
            interaction = TRUE
          )

          gamma <- tryCatch(
            {
              fit_log_basis_model(
                X = predictors_train,
                y = Y_a[train_id]
              )
            },
            error = function(e) {
              print("Log-basis model fitting produced an error")
              print(e)
            }
          )

          m_labeled[which(A_labeled == a)[test_id]] <- predict_log_basis_model(
            gamma,
            predictors_test
          )
        }
      }

      design_gbic <- cbind(1, predictors_labeled)
      gbic_scores <- c(
        gbic_scores,
        compute_candidate_gbic(
          Y_labeled[A_labeled == a],
          design_gbic,
          gbic_splits[[as.character(a)]]
        )
      )
      mbic_scores <- c(
        mbic_scores,
        compute_candidate_mbic(
          Y_labeled[A_labeled == a],
          design_gbic
        )
      )
    } else if (basis_a == "Interaction") {
      X_int <- S * as.matrix(X)
      X_int <- cbind(X, X_int)
      X_labeled_int <- X_int[labeled_ind, , drop = FALSE]
      X_unlabeled_int <- X_int[-labeled_ind, , drop = FALSE]

      alpha <- tryCatch(
        {
          find_alpha_glm(
            Y = Y_labeled[A_labeled == a],
            covariates_matrix = S_labeled[A_labeled == a] %>% as.matrix(),
            additional_matrix = cbind(
              C_labeled[A_labeled == a],
              X_labeled_int[A_labeled == a, ]
            ) %>% as.matrix()
          )
        },
        error = function(e) 1
      )
      alphas <- c(alphas, alpha)

      basis_labeled <- polynomial(S_labeled %>% as.data.frame(), alpha)
      basis_unlabeled <- polynomial(S_unlabeled %>% as.data.frame(), alpha)

      gamma <- tryCatch(
        {
          RidgeRegression(
            X = cbind(
              basis_labeled[A_labeled == a, -1],
              C_labeled[A_labeled == a],
              X_labeled_int[A_labeled == a, ]
            ) %>% as.matrix(),
            y = Y_labeled[A_labeled == a],
            penalty_factor = ridge_design_penalty_factor(
              n_basis = ncol(basis_labeled[, -1, drop = FALSE]),
              D = C_labeled[A_labeled == a],
              X = X_labeled[A_labeled == a, , drop = FALSE],
              X_interaction = X_labeled_int[
                A_labeled == a,
                (ncol(X_labeled) + 1):(2 * ncol(X_labeled)),
                drop = FALSE
              ],
              unpenalize_binary_X = ridge_unpenalize_binary_X
            )
          )
        },
        error = function(e) {
          print("Ridge Regression produced an error")
          print(e)
        }
      )

      m_unlabeled[A_unlabeled == a] <- boot::inv.logit(
        as.matrix(
          cbind(
            1,
            basis_unlabeled[A_unlabeled == a, -1],
            C_unlabeled[A_unlabeled == a],
            X_unlabeled_int[A_unlabeled == a, ]
          )
        ) %*% gamma$coefficients
      )

      m_labeled[A_labeled == a] <- boot::inv.logit(
        as.matrix(
          cbind(
            1,
            basis_labeled[A_labeled == a, -1],
            C_labeled[A_labeled == a],
            X_labeled_int[A_labeled == a, ]
          )
        ) %*% gamma$coefficients
      )

      if (cross_fit) {
        Y_a <- Y_labeled[A_labeled == a]
        S_a <- S_labeled[A_labeled == a]
        X_a <- X_labeled_int[A_labeled == a, , drop = FALSE]
        X_main_a <- X_labeled[A_labeled == a, , drop = FALSE]
        C_a <- C_labeled[A_labeled == a]
        fold <- folds[[as.character(a)]]

        for (i in sort(unique(fold))) {
          train_id <- which(fold != i)
          test_id <- which(fold == i)
          alpha_fold <- tryCatch(
            {
              find_alpha_glm(
                Y = Y_a[train_id],
                covariates_matrix = as.matrix(S_a[train_id]),
                additional_matrix = as.matrix(cbind(X_a[train_id, ], C_a[train_id]))
              )
            },
            error = function(e) 1
          )

          basis_train <- polynomial(S_a[train_id] %>% as.data.frame(), alpha_fold)
          basis_test <- polynomial(S_a[test_id] %>% as.data.frame(), alpha_fold)

          gamma <- tryCatch(
            {
            RidgeRegression(
              X = cbind(
                basis_train[, -1],
                C_a[train_id],
                X_a[train_id, ]
              ) %>% as.matrix(),
              y = Y_a[train_id],
              weights = NULL,
              penalty_factor = ridge_design_penalty_factor(
                n_basis = ncol(basis_train[, -1, drop = FALSE]),
                D = C_a[train_id],
                X = X_main_a[train_id, , drop = FALSE],
                X_interaction = X_a[
                  train_id,
                  (ncol(X_main_a) + 1):(2 * ncol(X_main_a)),
                  drop = FALSE
                ],
                unpenalize_binary_X = ridge_unpenalize_binary_X
              )
            )
          },
          error = function(e) {
            print("Ridge Regression produced an error")
            print(e)
          }
        )$coefficients

          m_labeled[which(A_labeled == a)[test_id]] <- boot::inv.logit(
            as.matrix(
              cbind(1, basis_test[, -1], C_a[test_id], X_a[test_id, ])
            ) %*% gamma
          )
        }
      }

      design_gbic <- as.matrix(
        cbind(
          1,
          basis_labeled[A_labeled == a, -1],
          C_labeled[A_labeled == a],
          X_labeled_int[A_labeled == a, ]
        )
      )
      gbic_scores <- c(
        gbic_scores,
        compute_candidate_gbic(
          Y_labeled[A_labeled == a],
          design_gbic,
          gbic_splits[[as.character(a)]]
        )
      )
      mbic_scores <- c(
        mbic_scores,
        compute_candidate_mbic(
          Y_labeled[A_labeled == a],
          design_gbic
        )
      )
    } else if (basis_a == "Beta") {
      alphas <- c(alphas, NA_real_)
      gbic_scores <- c(gbic_scores, Inf)
      mbic_scores <- c(mbic_scores, Inf)

      S_calibrated <- FitParametricCalibration(
        Y_labeled,
        S_labeled,
        A_labeled,
        A_val = a,
        S_unlabeled,
        method = "Beta"
      )
      m_unlabeled[A_unlabeled == a] <- S_calibrated[A_unlabeled == a]

      S_calibrated <- FitParametricCalibration(
        Y_labeled,
        S_labeled,
        A_labeled,
        A_val = a,
        S_labeled,
        method = "Beta"
      )
      m_labeled[A_labeled == a] <- S_calibrated[A_labeled == a]

      if (cross_fit) {
        Y_a <- Y_labeled[A_labeled == a]
        S_a <- S_labeled[A_labeled == a]
        fold <- folds[[as.character(a)]]

        for (i in sort(unique(fold))) {
          train_id <- which(fold != i)
          test_id <- which(fold == i)

          S_calibrated_fold <- FitParametricCalibration(
            Y_labeled = Y_a[train_id],
            S_labeled = S_a[train_id],
            A_labeled = rep(a, length(train_id)),
            A_val = a,
            S_unlabeled = S_a[test_id],
            method = "Beta"
          )

          m_labeled[which(A_labeled == a)[test_id]] <- S_calibrated_fold
        }
      }
    } else if (basis_a == "kernel") {
      alphas <- c(alphas, NA_real_)
      gbic_scores <- c(gbic_scores, Inf)
      mbic_scores <- c(mbic_scores, Inf)

      bandwidth <- kernel_bandwidth(
        S_labeled_metric[A_labeled == a],
        order = kernel_order
      )
      mhat <- npreg(
        S_labeled_metric[A_labeled == a],
        Y_labeled[A_labeled == a],
        S_unlabeled_metric,
        bandwidth
      )
      m_unlabeled[A_unlabeled == a] <- mhat[A_unlabeled == a]

      m_labeled[A_labeled == a] <- npreg(
        S_labeled_metric[A_labeled == a],
        Y_labeled[A_labeled == a],
        S_labeled_metric[A_labeled == a],
        bandwidth
      )

      if (cross_fit) {
        Y_a <- Y_labeled[A_labeled == a]
        S_a <- S_labeled_metric[A_labeled == a]
        fold <- folds[[as.character(a)]]
        for (i in sort(unique(fold))) {
          train_id <- which(fold != i)
          test_id <- which(fold == i)
          bandwidth_fold <- kernel_bandwidth(
            S_a[train_id],
            order = kernel_order
          )
          mhat_fold <- npreg(S_a[train_id], Y_a[train_id], S_a[test_id], bandwidth_fold)
          m_labeled[which(A_labeled == a)[test_id]] <- mhat_fold
        }
      }
    } else {
      stop("Basis not recognized.")
    }
  }
  
  est <- get_metric(
    Y = m_unlabeled, S = S_unlabeled,
    A = A_unlabeled, threshold = threshold_unlabeled_raw, W = NULL
  )

  # Model-selection diagnostics.
  weighted_mse <- c()
  for (a in nclass) {
    idx_a <- A_labeled == a
    y_a <- Y_labeled[idx_a]
    p_a <- m_labeled[idx_a]
    d_a <- C_labeled_metric[idx_a]

    # TPR-weighted MSE used for model selection. Use the same model-specific
    # semi-supervised TPR estimate that is passed into the influence curve.
    tpr_hat_a <- est[est$Metric == "TPR", paste0("Group", a)]
    if (is.finite(tpr_hat_a)) {
      weight_a <- (d_a - tpr_hat_a)^2
    } else {
      weight_a <- rep(1, length(y_a))
    }

    if (sum(weight_a) > 0) {
      weighted_mse_a <- mean(weight_a * (y_a - p_a)^2)
    } else {
      weighted_mse_a <- mean((y_a - p_a)^2)
    }
    weighted_mse <- c(weighted_mse, weighted_mse_a)
  }
  metrics_by_group <- cbind(
    BIC = gbic_scores,
    MBIC = mbic_scores,
    Weighted_MSE = weighted_mse
  ) %>% as.data.frame()

  rownames(metrics_by_group) <- paste0("Group_", nclass)
  
  var <- Influence_curve(
    est, Y_labeled, S_labeled, A_labeled, m_labeled, threshold_labeled_raw,
    method = "semi-supervised"
  )

  if (any(basis == "Beta")) {
    var <- beta_if_variance_update(
      var_est = var,
      est = est,
      Y_labeled = Y_labeled,
      S_labeled = S_labeled,
      A_labeled = A_labeled,
      S_unlabeled = S_unlabeled,
      A_unlabeled = A_unlabeled,
      threshold_unlabeled = threshold_unlabeled_raw,
      basis = basis
    )
  }

  if (control_variate) {
    cv_update <- apply_control_variate_update(
      est_plugin = est,
      var_plugin = var,
      Y_labeled = Y_labeled,
      S_labeled = S_labeled,
      A_labeled = A_labeled,
      m_labeled = m_labeled,
      threshold_labeled = threshold_labeled_raw,
      S_unlabeled = S_unlabeled,
      A_unlabeled = A_unlabeled,
      m_unlabeled = m_unlabeled,
      threshold_unlabeled = threshold_unlabeled_raw
    )
    est <- cv_update$est
    var <- cv_update$var
    cv_weights <- cv_update$cv_weights
  }

  if (variance_method == "bootstrap") {
    var <- bootstrap_variance_semisupervised(
      template_est = est,
      Y = Y,
      S = S,
      A = A,
      threshold = threshold,
      X = X,
      basis = basis,
      nknots = nknots,
      k = k,
      cross_fit_variance = cross_fit,
      control_variate = control_variate,
      spline_ridge = spline_ridge,
      spline_gam = spline_gam,
      ridge_unpenalize_binary_X = ridge_unpenalize_binary_X,
      spline_include_D = spline_include_D,
      log_include_D = log_include_D,
      kernel_order = kernel_order,
      kernel_use_ecdf = kernel_use_ecdf,
      B = bootstrap_reps,
      resample = bootstrap_resample,
      seed = bootstrap_seed,
      bootstrap_cores = bootstrap_cores,
      ...
    )
  }
  
  result <- list(
    est = est,
    var = var,
    alpha = alphas,
    imp_quality = metrics_by_group,
    m_labeled = m_labeled,
    m_unlabeled = m_unlabeled,
    cv_folds = folds
  )

  if (control_variate) {
    result$cv_weights <- cv_weights
  }

  return(result)
}
