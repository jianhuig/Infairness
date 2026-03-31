# Toy mean-estimation experiments aligned with the current SSFairnessAudit package.

library(parallel)

get_option_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  prefix <- paste0("--", name, "=")
  hit <- grep(prefix, args, value = TRUE)
  if (length(hit) == 0) {
    return(default)
  }
  sub(prefix, "", hit[1], fixed = TRUE)
}

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  match <- grep(file_arg, args, value = TRUE)
  if (length(match) == 0) {
    return(normalizePath(getwd()))
  }
  normalizePath(dirname(sub(file_arg, "", match[1])))
}

load_package_namespace <- function() {
  repo_root <- normalizePath(file.path(get_script_dir(), ".."))

  if (!"SSFairnessAudit" %in% loadedNamespaces()) {
    if (requireNamespace("pkgload", quietly = TRUE)) {
      pkgload::load_all(repo_root, quiet = TRUE, helpers = FALSE)
    } else {
      stop("Please install `pkgload` or install the `SSFairnessAudit` package first.")
    }
  }

  asNamespace("SSFairnessAudit")
}

pkg_ns <- load_package_namespace()
script_dir <- get_script_dir()
output_dir <- normalizePath(
  get_option_arg("output-dir", file.path(script_dir, "Data")),
  mustWork = FALSE
)
nsim <- as.integer(get_option_arg("nsim", "10000"))
cores <- as.integer(get_option_arg("cores", as.character(detectCores())))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

FitParametricCalibration <- get("FitParametricCalibration", envir = pkg_ns)
NaturalSplineBasis <- get("NaturalSplineBasis", envir = pkg_ns)
RidgeRegression <- get("RidgeRegression", envir = pkg_ns)
SplineRidgeRegression <- get("SplineRidgeRegression", envir = pkg_ns)
find_alpha_glm <- get("find_alpha_glm", envir = pkg_ns)
polynomial <- get("polynomial", envir = pkg_ns)
npreg <- get("npreg", envir = pkg_ns)

score_to_probability <- function(s) {
  if (all(s >= 0 & s <= 1)) {
    s
  } else {
    plogis(s)
  }
}

score_indicator <- function(s, threshold = 0.5) {
  as.integer(score_to_probability(s) > threshold)
}

estimate_oracle <- function(dat) {
  mean(dat$Y[dat$A == 1]) - mean(dat$Y[dat$A == 0])
}

estimate_supervised <- function(lab) {
  mean(lab$Y[lab$A == 1]) - mean(lab$Y[lab$A == 0])
}

estimate_naive <- function(dat) {
  mean(score_to_probability(dat$S[dat$A == 1])) - mean(score_to_probability(dat$S[dat$A == 0]))
}

estimate_beta <- function(lab, unlab) {
  s_lab <- score_to_probability(lab$S)
  s_unlab <- score_to_probability(unlab$S)

  m0 <- tryCatch(
    FitParametricCalibration(
      Y_labeled = lab$Y,
      S_labeled = s_lab,
      A_labeled = lab$A,
      A_val = 0,
      S_unlabeled = s_unlab,
      method = "Beta"
    ),
    error = function(e) NULL
  )
  m1 <- tryCatch(
    FitParametricCalibration(
      Y_labeled = lab$Y,
      S_labeled = s_lab,
      A_labeled = lab$A,
      A_val = 1,
      S_unlabeled = s_unlab,
      method = "Beta"
    ),
    error = function(e) NULL
  )

  if (is.null(m0) || is.null(m1)) {
    return(NA_real_)
  }

  mean(m1[unlab$A == 1]) - mean(m0[unlab$A == 0])
}

estimate_kernel <- function(lab, unlab) {
  bw0 <- sd(lab$S[lab$A == 0]) / (sum(lab$A == 0)^0.45)
  bw1 <- sd(lab$S[lab$A == 1]) / (sum(lab$A == 1)^0.45)

  m0 <- npreg(
    St = lab$S[lab$A == 0],
    Yt = lab$Y[lab$A == 0],
    Sv = unlab$S[unlab$A == 0],
    bw = bw0
  )
  m1 <- npreg(
    St = lab$S[lab$A == 1],
    Yt = lab$Y[lab$A == 1],
    Sv = unlab$S[unlab$A == 1],
    bw = bw1
  )

  mean(m1) - mean(m0)
}

estimate_poly <- function(lab, unlab, x_cols = NULL, threshold = 0.5) {
  out <- numeric(2)

  for (a in 0:1) {
    lab_g <- lab[lab$A == a, , drop = FALSE]
    unlab_g <- unlab[unlab$A == a, , drop = FALSE]

    c_lab <- score_indicator(lab_g$S, threshold)
    c_unlab <- score_indicator(unlab_g$S, threshold)

    additional <- as.matrix(c_lab)
    if (!is.null(x_cols)) {
      additional <- cbind(c_lab, as.matrix(lab_g[, x_cols, drop = FALSE]))
    }

    alpha <- tryCatch(
      find_alpha_glm(
        Y = lab_g$Y,
        covariates_matrix = as.matrix(lab_g$S),
        additional_matrix = additional
      ),
      error = function(e) 1
    )

    basis_lab <- polynomial(as.data.frame(lab_g$S), alpha)
    basis_unlab <- polynomial(as.data.frame(unlab_g$S), alpha)

    design_lab <- cbind(basis_lab[, -1, drop = FALSE], c_lab)
    design_unlab <- cbind(basis_unlab[, -1, drop = FALSE], c_unlab)
    if (!is.null(x_cols)) {
      design_lab <- cbind(design_lab, as.matrix(lab_g[, x_cols, drop = FALSE]))
      design_unlab <- cbind(design_unlab, as.matrix(unlab_g[, x_cols, drop = FALSE]))
    }

    fit <- RidgeRegression(X = as.matrix(design_lab), y = lab_g$Y)
    preds <- plogis(cbind(1, as.matrix(design_unlab)) %*% fit$coefficients)
    out[a + 1] <- mean(preds)
  }

  out[2] - out[1]
}

estimate_spline <- function(lab, unlab, x_cols = NULL, nknots = 3, threshold = 0.5) {
  out <- numeric(2)

  for (a in 0:1) {
    lab_g <- lab[lab$A == a, , drop = FALSE]
    unlab_g <- unlab[unlab$A == a, , drop = FALSE]

    c_lab <- score_indicator(lab_g$S, threshold)
    c_unlab <- score_indicator(unlab_g$S, threshold)

    basis_lab <- NaturalSplineBasis(as.matrix(lab_g$S), nknots, return_knots = TRUE)
    spline_knots <- attr(basis_lab, "knots")
    basis_unlab <- NaturalSplineBasis(
      as.matrix(unlab_g$S),
      nknots,
      knots = spline_knots
    )

    design_lab <- cbind(basis_lab, c_lab)
    design_unlab <- cbind(basis_unlab, c_unlab)
    if (!is.null(x_cols)) {
      design_lab <- cbind(design_lab, as.matrix(lab_g[, x_cols, drop = FALSE]))
      design_unlab <- cbind(design_unlab, as.matrix(unlab_g[, x_cols, drop = FALSE]))
    }

    fit <- SplineRidgeRegression(X = as.matrix(design_lab), y = lab_g$Y)
    preds <- plogis(cbind(1, as.matrix(design_unlab)) %*% fit$coefficients)
    out[a + 1] <- mean(preds)
  }

  out[2] - out[1]
}

run_estimators <- function(dat, label_size = 400, nknots = 3) {
  idx <- sample(nrow(dat), label_size)
  lab <- dat[idx, , drop = FALSE]
  unlab <- dat[-idx, , drop = FALSE]
  x_cols <- grep("^W", names(dat), value = TRUE)

  data.frame(
    method = c(
      "oracle",
      "supervised",
      "naive",
      "beta calibration",
      "kernel smoothing",
      "poly(S)",
      "poly(S+W)",
      "spline(S)",
      "spline(S+W)"
    ),
    estimate = c(
      estimate_oracle(dat),
      estimate_supervised(lab),
      estimate_naive(dat),
      estimate_beta(lab, unlab),
      estimate_kernel(lab, unlab),
      estimate_poly(lab, unlab, x_cols = NULL),
      estimate_poly(lab, unlab, x_cols = x_cols),
      estimate_spline(lab, unlab, x_cols = NULL, nknots = nknots),
      estimate_spline(lab, unlab, x_cols = x_cols, nknots = nknots)
    )
  )
}

draw_unit <- function(p) {
  z <- rnorm(p)
  z / sqrt(sum(z^2))
}

simulate_toy1 <- function(n = 2e4, p = 5) {
  beta_X <- draw_unit(p)
  beta_W <- draw_unit(p)

  A <- rbinom(n, 1, 0.5)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  W <- matrix(rnorm(n * p), n, p)
  colnames(W) <- paste0("W", 1:p)

  eta <- (X %*% beta_X + W %*% beta_W) / sqrt(2)
  alpha <- c(-0.5, 0.5)
  S <- plogis(alpha[A + 1] + eta)
  Y <- rbinom(n, 1, S)

  data.frame(Y, A, S, X, W)
}

simulate_toy2 <- function(n = 2e4, p = 5) {
  beta_X <- draw_unit(p)
  beta_W <- draw_unit(p)

  A <- rbinom(n, 1, 0.5)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  W <- matrix(rnorm(n * p), n, p)
  colnames(W) <- paste0("W", 1:p)

  S <- (X %*% beta_X + W %*% beta_W) / sqrt(2)
  alpha <- c(-0.5, 0.5)
  prob <- plogis(alpha[A + 1] + S)
  Y <- rbinom(n, 1, prob)

  data.frame(Y, A, S, X, W)
}

simulate_toy3 <- function(n = 2e4, p = 5) {
  beta_X <- draw_unit(p)
  beta_W <- draw_unit(p)

  A <- rbinom(n, 1, 0.5)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  W <- matrix(rnorm(n * p), n, p)
  colnames(W) <- paste0("W", 1:p)

  S <- (X %*% beta_X + W %*% beta_W) / sqrt(2)
  alpha <- c(-0.5, 0.5)
  prob <- plogis(alpha[A + 1] + S + 0.5 * S^2)
  Y <- rbinom(n, 1, prob)

  data.frame(Y, A, S, X, W)
}

simulate_toy4 <- function(n = 2e4, p = 5) {
  beta_X <- draw_unit(p)
  beta_W <- draw_unit(p)

  A <- rbinom(n, 1, 0.5)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  W <- matrix(rnorm(n * p), n, p)
  colnames(W) <- paste0("W", 1:p)

  S <- X %*% beta_X
  alpha <- c(-0.5, 0.5)
  prob <- plogis(alpha[A + 1] + S + (W %*% beta_W))
  Y <- rbinom(n, 1, prob)

  data.frame(Y, A, S, X, W)
}

run_scenario <- function(simulator, nsim = 1e4, cores = detectCores(), output_file, nknots = 3) {
  results <- mclapply(seq_len(nsim), function(i) {
    dat <- simulator()
    run_estimators(dat, nknots = nknots)
  }, mc.cores = cores)

  saveRDS(results, file = output_file)
  invisible(results)
}

run_scenario(simulate_toy1, nsim = nsim, cores = cores, output_file = file.path(output_dir, "toy1.rds"))
run_scenario(simulate_toy2, nsim = nsim, cores = cores, output_file = file.path(output_dir, "toy2.rds"))
run_scenario(simulate_toy3, nsim = nsim, cores = cores, output_file = file.path(output_dir, "toy3.rds"))
run_scenario(simulate_toy4, nsim = nsim, cores = cores, output_file = file.path(output_dir, "toy4.rds"))
