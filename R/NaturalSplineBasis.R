# Updated: 2021-04-19

#' Computes the natural spline basis.
#'
#' @param X Covariate matrix.
#' @param num_knots Number of knots.
#' @param knots Optional list of knot vectors, one per column of `X`.
#' @param return_knots Logical; if `TRUE`, attach the knot locations used to the
#' returned matrix as an attribute.
#' @export
#' @return Matrix containing natural spline basis.
#' @importFrom stats coef glm quantile rbeta rbinom rlogis rnorm
#'
NaturalSplineBasis <- function(X, num_knots, knots = NULL, return_knots = FALSE) {
  X <- as.matrix(X)
  basis_parts <- vector("list", ncol(X))
  knots_used <- vector("list", ncol(X))

  for (i in 1:ncol(X)) {
    X_i <- X[, i]

    if (is.null(knots)) {
      current_knots <- quantile(X_i, seq(0, 1, length = num_knots))
      j <- 0

      while (length(unique(current_knots)) != num_knots) {
        j <- j + 1
        current_knots <- unique(quantile(X_i, seq(0, 1, length = num_knots + j)))
      }
    } else {
      current_knots <- knots[[i]]
    }
    knots_used[[i]] <- current_knots

    # Compute the natural spline basis.
    d_k <- (TruncatedCubic(X_i, current_knots[num_knots - 1])) /
      (current_knots[num_knots] - current_knots[num_knots - 1])
    evals <- sapply(1:(num_knots - 2), function(ii) {
      d_i <- (TruncatedCubic(X_i, current_knots[ii])) / (current_knots[num_knots] - current_knots[ii])
      basis.new <- d_i - d_k
    })
    evals <- matrix(evals, nrow = length(X_i))

    # Bind original variable and basis.
    basis_parts[[i]] <- cbind(X_i, evals)
  }

  # Return basis including everything.
  basis <- do.call(cbind, basis_parts)
  if (return_knots) {
    attr(basis, "knots") <- knots_used
  }

  return(basis)
}


TruncatedCubic <- function(x, knot_location) {
  return(((x > knot_location) * (x - knot_location))^3)
}

# Use a slightly more conservative ridge path for spline-expanded designs,
# which typically have wider bases than the polynomial branches.
SplineRidgeRegression <- function(X, y, weights = NULL, penalty_factor = NULL) {
  RidgeRegression(
    X = X,
    y = y,
    weights = weights,
    exponents = seq(0.8, 2.8, length.out = 120),
    penalty_factor = penalty_factor
  )
}

spline_indicator_term <- function(C, include_indicator = TRUE) {
  if (!include_indicator) {
    return(NULL)
  }

  as.matrix(C)
}

build_spline_predictors <- function(basis_matrix,
                                    D = NULL,
                                    X = NULL,
                                    X_interaction = NULL) {
  parts <- list(as.matrix(basis_matrix))

  if (!is.null(D)) {
    parts[[length(parts) + 1]] <- as.matrix(D)
  }
  if (!is.null(X)) {
    parts[[length(parts) + 1]] <- as.matrix(X)
  }
  if (!is.null(X_interaction)) {
    parts[[length(parts) + 1]] <- as.matrix(X_interaction)
  }

  as.matrix(do.call(cbind, parts))
}

build_spline_penalty_factor <- function(basis_matrix,
                                        D = NULL,
                                        X = NULL,
                                        X_interaction = NULL,
                                        unpenalize_binary_X = FALSE) {
  ridge_design_penalty_factor(
    n_basis = ncol(as.matrix(basis_matrix)),
    D = D,
    X = X,
    X_interaction = X_interaction,
    unpenalize_binary_X = unpenalize_binary_X
  )
}

spline_gam_k <- function(S, num_knots) {
  unique_n <- length(unique(as.numeric(S)))
  min(max(4L, num_knots + 1L), max(4L, unique_n))
}

spline_gam_is_binary <- function(x) {
  x_no_na <- unique(stats::na.omit(as.numeric(x)))
  length(x_no_na) <= 2L
}

build_spline_gam_data <- function(S,
                                  D = NULL,
                                  X = NULL) {
  data <- data.frame(S = as.numeric(S))

  if (!is.null(D)) {
    data$D <- as.numeric(D)
  }

  if (!is.null(X)) {
    X_df <- as.data.frame(X)
    names(X_df) <- paste0("X", seq_len(ncol(X_df)))
    data <- cbind(data, X_df)
  }

  data
}

spline_gam_formula <- function(data,
                               interaction = FALSE,
                               gam_k = 4L) {
  terms <- c(sprintf("s(S, bs = 'cs', k = %d)", gam_k))

  if ("D" %in% names(data)) {
    terms <- c(terms, "D")
  }

  x_names <- grep("^X", names(data), value = TRUE)
  if (length(x_names) > 0) {
    if (!interaction) {
      terms <- c(terms, x_names)
    } else {
      for (x_name in x_names) {
        x_values <- data[[x_name]]
        if (spline_gam_is_binary(x_values)) {
          terms <- c(
            terms,
            x_name,
            sprintf("s(S, by = %s, bs = 'cs', k = %d)", x_name, gam_k)
          )
        } else {
          x_k <- min(gam_k, max(4L, length(unique(stats::na.omit(as.numeric(x_values))))))
          terms <- c(
            terms,
            sprintf("s(%s, bs = 'cs', k = %d)", x_name, x_k),
            sprintf("ti(S, %s, bs = c('cs', 'cs'), k = c(%d, %d))", x_name, gam_k, x_k)
          )
        }
      }
    }
  }

  stats::as.formula(paste("y ~", paste(terms, collapse = " + ")))
}

fit_spline_model <- function(X,
                             y,
                             weights = NULL,
                             use_ridge = TRUE,
                             use_gam = FALSE,
                             S = NULL,
                             D = NULL,
                             X_cov = NULL,
                             interaction = FALSE,
                             num_knots = 3,
                             penalty_factor = NULL) {
  X <- as.matrix(X)

  if (use_gam) {
    gam_data <- build_spline_gam_data(S = S, D = D, X = X_cov)
    gam_data$y <- y
    gam_k <- spline_gam_k(S = S, num_knots = num_knots)
    if (is.null(weights)) {
      fit <- mgcv::gam(
        formula = spline_gam_formula(gam_data, interaction = interaction, gam_k = gam_k),
        data = gam_data,
        family = stats::binomial(),
        method = "REML",
        select = TRUE
      )
    } else {
      fit <- mgcv::gam(
        formula = spline_gam_formula(gam_data, interaction = interaction, gam_k = gam_k),
        data = gam_data,
        family = stats::binomial(),
        weights = weights,
        method = "REML",
        select = TRUE
      )
    }

    return(list(type = "gam", fit = fit))
  }

  if (use_ridge) {
    fit <- SplineRidgeRegression(
      X = X,
      y = y,
      weights = weights,
      penalty_factor = penalty_factor
    )
    return(list(type = "linear", coefficients = fit$coefficients))
  }

  x_design <- cbind("(Intercept)" = 1, X)
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
  list(type = "linear", coefficients = coefficients)
}

predict_spline_model <- function(model,
                                 X,
                                 S = NULL,
                                 D = NULL,
                                 X_cov = NULL) {
  if (identical(model$type, "gam")) {
    newdata <- build_spline_gam_data(S = S, D = D, X = X_cov)
    return(clamp_probabilities(stats::predict(model$fit, newdata = newdata, type = "response")))
  }

  boot::inv.logit(as.matrix(cbind(1, X)) %*% model$coefficients)
}
