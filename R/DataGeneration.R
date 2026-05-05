# Purpose: Data generation for simulation.
# Updated: 2025-5-12

#' Data generation.
#'
#' @param n_labeled Number of labeled examples.
#' @param N_unlabeled Number of unlabeled examples.
#' @param prot_att_prevalence Prevalence of the protected attribute. For a
#' single value, this is interpreted as `P(A = 1)`. For a length-two vector,
#' entries are interpreted as the relative prevalences of groups 0 and 1.
#' @param model Indicator to generate from which model.
#' possible options: "scenario 1", "scenario 2"
#' @param rho Correlation between covariates. Default is 0.4.
#' @return Data.frame.
#' @export

resolve_binary_protected_prevalence <- function(prot_att_prevalence) {
  prevalence <- as.numeric(prot_att_prevalence)

  if (length(prevalence) == 1L) {
    p_group1 <- prevalence
  } else if (length(prevalence) == 2L) {
    if (any(prevalence < 0) || sum(prevalence) <= 0) {
      stop("`prot_att_prevalence` must contain nonnegative values with positive sum.")
    }
    p_group1 <- prevalence[2] / sum(prevalence)
  } else {
    stop("This data generator currently supports binary protected attributes only.")
  }

  if (!is.finite(p_group1) || p_group1 <= 0 || p_group1 >= 1) {
    stop("Protected-attribute prevalence for group 1 must be strictly between 0 and 1.")
  }

  p_group1
}

generate_binary_protected_attribute <- function(latent_A,
                                                prot_att_prevalence) {
  p_group1 <- resolve_binary_protected_prevalence(prot_att_prevalence)
  n <- length(latent_A)
  n_group1 <- round(n * p_group1)
  n_group1 <- min(max(n_group1, 1L), n - 1L)

  A <- integer(n)
  A[order(latent_A, decreasing = TRUE)[seq_len(n_group1)]] <- 1L
  A
}

DataGeneration <- function(n_labeled,
                           N_unlabeled,
                           prot_att_prevalence,
                           model,
                           rho = 0.4) {
  # Total sample size.
  N_total <- n_labeled + N_unlabeled
  
  # Dimension
  p <- 16
  
  # Generate Covariates
  id_matrix <- diag(p)
  ar_one_matrix <- rho^abs(row(id_matrix) - col(id_matrix))
  sigma <- 3 * ar_one_matrix
  mu <- rep(0, p)
  covariates <- MASS::mvrnorm(N_total, mu, Sigma = sigma)
  colnames(covariates) <- c(paste0("X_", 1:10), paste0("W_", 1:5), "A")
  
  # Placeholder for Y
  Y <- rep(NA, N_total)
  A <- generate_binary_protected_attribute(
    covariates[, "A"],
    prot_att_prevalence
  )
  W <- covariates[, 11:15]
  X <- covariates[, 1:10]
  
  # Placeholder for Y
  Y <- rep(NA, N_total)
  
  # Generate Y
  if (model == "scenario 1") {
    b0 <- matrix(
      c(
        -4, 1, 1, 0.5, 0.5, rep(0, 6), 0.4, 0.4, 0.4, 0, 0,
        -4, 0.9, 0.9, 0.4, 0.4, rep(0, 6), 0.3, 0.3, 0.3, 0, 0
      ),
      nrow = 2, byrow = TRUE
    )
    lin_pred <- cbind(1, X, W) %*% t(b0)
    S <- plogis(
      lin_pred + 0.3 * (X[, 2])^2 - 0.4 * (X[, 3])^3 + 0.1 * X[, 5] * X[, 6]
    )
    for (a in c(0, 1)) {
      Y[A == a] <- rbinom(sum(A == a), 1, S[A == a, (a + 1)])
    }
  }
  if (model == "scenario 2") {
    b0 <- matrix(
      c(
        1.3, 0.4, -0.3, 0.15, -0.15, rep(0, 6), 0.25, -0.2, 0.2, 0, 0,
        1.3, 0.35, -0.25, 0.2, -0.2, rep(0, 6), 0.15, -0.15, 0.2, 0, 0
      ),
      nrow = 2, byrow = TRUE
    )
    lin_pred <- cbind(1, X, W) %*% t(b0)
    for (a in c(0, 1)) {
      S <- exp(-lin_pred[, a + 1]^2)
      Y[A == a] <- rbinom(sum(A == a), 1, S[A == a])
    }
  } else if (model == "scenario 3") {
    ## ---- GLM-like index using only X1:5 and W1:3 (plus interactions / higher order) ----
    # Uses X = covariates[,1:10] and W = covariates[,11:15]
    t <- -0.5 +
      1.00*X[,1] + 0.80*X[,2] + 0.60*X[,3] + 0.50*X[,4] + 0.30*X[,5] +
      0.30*W[,1] + 0.20*W[,2] - 0.20*W[,3] +
      0.10*X[,1]*W[,1] - 0.10*X[,5]*X[,2] +          # interactions
      0.12*(X[,3]^2) - 0.08*(X[,3]^3)/3 + 0.10*(W[,1]^2)  # higher-order
    z <- as.numeric(scale(t))  # center/scale around operating region
    
    ## ---- Same shape for both groups; small coefficient tweaks only ----
    s0 <- 0.18
    # group-specific multipliers (small differences)
    b_z       <- c(3.20, 3.36)      # slope on z (≈ +5% in group 1)
    amp_up_1  <- c(7.50, 7.95)      # left bump amplitude (≈ +6%)
    amp_down  <- c(-12.50, -13.25)  # central dip amplitude (≈ +6%)
    amp_up_2  <- c(9.00, 9.45)      # right bump amplitude (≈ +5%)
    ripple_amp<- c(0.60, 0.64)      # ripple amplitude (≈ +7%)
    ripple_k  <- 4                  # same ripple frequency for both
    
    g <- numeric(N_total)
    for (a in 0:1) {
      idx <- which(A == a)
      zz  <- z[idx]
      ga <- b_z[a+1]*zz +
        amp_up_1[a+1] * exp(- (zz + 0.30)^2 / (2*s0^2)) +
        amp_down[a+1] * exp(- (zz + 0.00)^2 / (2*s0^2)) +
        amp_up_2[a+1] * exp(- (zz - 0.30)^2 / (2*s0^2)) +
        ripple_amp[a+1] * exp(-(zz/0.50)^2) * sin(2*pi*ripple_k*plogis(zz))
      g[idx] <- ga
    }
    
    ## ---- Same inverse link in both groups (cloglog), per-group prevalence ~30% ----
    p_true <- numeric(N_total)
    for (a in 0:1) {
      idx <- which(A == a)
      ga  <- g[idx]
      f <- function(d) mean(1 - exp(-exp(ga + d))) - 0.30
      delta <- uniroot(f, interval = c(-15, 15))$root
      p_true[idx] <- 1 - exp(-exp(ga + delta))
    }
    p_true <- pmin(pmax(p_true, 1e-6), 1 - 1e-6)
    
    ## ---- Sample labels ----
    Y <- rbinom(N_total, 1, p_true)
  }
  # Induce missingness.
  Y_miss <- Y
  Y_miss[sample(c(1:N_total), N_unlabeled, replace = F)] <- NA
  
  # Simulated data.
  my_data <- cbind(Y = Y, A = A, Y_miss = Y_miss, X = X, W = W)
  my_data <- data.frame(my_data)
  
  return(my_data)
}
