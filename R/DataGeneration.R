# Purpose: Data generation for simulation.
# Updated: 2023-02-05

#' Data generation.
#'
#' @param n_labeled Number of labeled examples.
#' @param N_unlabeled Number of unlabeled examples.
#' @param prot_att_prevalence Prevalence of the protected attribute
#' If the number of class is more than two, the length of this vector
#' should be equal to the number of class. If the number of class is two,
#' the length of this vector can either be one or two.
#' @param Y_prevalence Prevalence of Y.
#' @param model Indicator to generate from which model.
#' Default is logistic model.
#' possible options: "correct", "incorrect outcome", "incorrect outcome and imputation 1", and "incorrect outcome and imputation 2".
#' @param b0 Vector of parameters to generate Y.
#' @param mis_para Vector of parameters to generate S when model is "mis".
#' @param n_class Number of class in the protected attribute. Default is two
#' @param p Number of covariates. Default is 10.
#' @param rho Correlation between covariates. Default is 0.2.
#' @return Data.frame.
#' @export

DataGeneration <- function(n_labeled,
                           N_unlabeled,
                           prot_att_prevalence,
                           model = "correct",
                           b0 = NULL,
                           p = 10,
                           rho = 0.4,
                           n_class = 2) {
  # Total sample size.
  N_total <- n_labeled + N_unlabeled

  # Generate Group Membership
  A_total <- round(N_total * prot_att_prevalence)
  N_total <- sum(A_total)
  A <- sample(unlist(lapply(seq_along(A_total), function(i) rep(i, A_total[i]))), N_total,
    replace = FALSE
  )

  # Generate Covariates
  id_matrix <- diag(p)
  ar_one_matrix <- rho^abs(row(id_matrix) - col(id_matrix))
  sigma <- 3 * ar_one_matrix
  X <- cbind(1, MASS::mvrnorm(N_total, rep(0, p), Sigma = sigma))
  colnames(X) <- paste0("X_", 0:p)

  # Linear Predictor
  lin_pred <- X %*% t(b0)
  # epsilon logistic
  eps_logistic <- rlogis(N_total)
  # epsilon extreme
  eps_extreme <- evd::rgumbel(N_total, -2, 0.3)

  # Generate Y
  if (model == "correct") {
    Y <- rep(NA, N_total)
    for (a in 1:n_class) {
      Y[A == a] <- ifelse(
        lin_pred[A == a, a] + eps_logistic[A == a] > 0, 1, 0
      )
    }
  }
  if (model == "incorrect outcome") {
    Y <- rep(NA, N_total)
    S <- boot::inv.logit(lin_pred)
    # basis <- ns.basis(lin_pred, nk = 3)
    # print(head(basis))
    effect_size <- c(1, 1)
    for (a in 1:n_class) {
      Y[A == a] <- ifelse(
        lin_pred[A == a, a] +
          effect_size[a] * ns.basis(S[, a], nk = 3)[A == a, 2] +
          eps_logistic[A == a] > 0, 1, 0
      )
    }
  }
  if (model == "incorrect outcome and imputation 1") {
    Y <- rep(NA, N_total)
    for (a in 1:n_class) {
      Y[A == a] <- ifelse(
        lin_pred[A == a, a] + 0.5 * (X[, 2] * X[, 3] + X[, 2] * X[, 6] - X[, 3] * X[, 7])[A == a] +
          eps_logistic[A == a] > 0, 1, 0
      )
    }
  }
  if (model == "incorrect outcome and imputation 2") {
    Y <- rep(NA, N_total)

    for (a in 1:n_class) {
      Y[A == a] <- ifelse(
        lin_pred[A == a, a] + X[A == a, 2]^2 + X[A == a, 3]^2 + ((exp(-2 - 3 * X[, 5] - 3 * X[, 6]))[A == a]) * eps_extreme[A == a] > 0, 1, 0
      )
    }
  }


  # Induce missingness.
  Y_miss <- Y
  Y_miss[sample(c(1:N_total), N_unlabeled, replace = F)] <- NA

  # Simulated data.
  my_data <- cbind(Y = Y, A = A, Y_miss = Y_miss, X[, -1])
  my_data <- data.frame(my_data)

  return(my_data)
}
