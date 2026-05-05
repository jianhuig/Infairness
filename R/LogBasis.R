# Purpose: Log-basis helpers for semi-supervised estimators.

clamp_probabilities <- function(S, eps = 1e-6) {
  pmin(pmax(S, eps), 1 - eps)
}

log_basis_terms <- function(S) {
  S_safe <- clamp_probabilities(S)
  cbind(log(S_safe), log(1 - S_safe))
}

build_log_predictors <- function(S,
                                 D = NULL,
                                 X = NULL,
                                 interaction = FALSE) {
  base <- log_basis_terms(S)
  parts <- list(base)

  if (!is.null(D)) {
    parts[[length(parts) + 1]] <- as.matrix(D)
  }

  if (!is.null(X)) {
    X_mat <- as.matrix(X)
    parts[[length(parts) + 1]] <- X_mat

    if (interaction) {
      parts[[length(parts) + 1]] <- sweep(X_mat, 1, base[, 1], `*`)
      parts[[length(parts) + 1]] <- sweep(X_mat, 1, base[, 2], `*`)
    }
  }

  as.matrix(do.call(cbind, parts))
}

fit_log_basis_model <- function(X, y, weights = NULL) {
  x_design <- cbind("(Intercept)" = 1, as.matrix(X))

  if (is.null(weights)) {
    fit <- stats::glm.fit(
      x = x_design,
      y = y,
      family = stats::binomial()
    )
  } else {
    fit <- stats::glm.fit(
      x = x_design,
      y = y,
      weights = weights,
      family = stats::binomial()
    )
  }

  coefficients <- fit$coefficients
  coefficients[is.na(coefficients)] <- 0
  list(coefficients = coefficients)
}

predict_log_basis_model <- function(model, X) {
  boot::inv.logit(as.matrix(cbind(1, X)) %*% model$coefficients)
}
