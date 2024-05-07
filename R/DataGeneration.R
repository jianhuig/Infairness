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
#' @param model Indicator to generate from which model.
#' Default is logistic model.
#' possible options: "correct", "incorrect outcome", "incorrect outcome and imputation 1", and "incorrect outcome and imputation 2".
#' @param b0 Vector of parameters to generate Y.
#' @param p Number of covariates. Default is 10.
#' @param rho Correlation between covariates. Default is 0.2.
#' @return Data.frame.
#' @export

DataGeneration <- function(n_labeled,
                           N_unlabeled,
                           prot_att_prevalence,
                           model = "correct",
                           b0 = NULL,
                           b1 = NULL,
                           p = 10,
                           rho = 0.4) {
  # Total sample size.
  N_total <- n_labeled + N_unlabeled
  
  # Generate Group Membership
  A <- rbinom(N_total, 1, prot_att_prevalence)
  
  # Generate Covariates
  id_matrix <- diag(p)
  ar_one_matrix <- rho^abs(row(id_matrix) - col(id_matrix))
  sigma <- 3 * ar_one_matrix
  mu <- rep(0, p)
  X <- cbind(1, MASS::mvrnorm(N_total, mu, Sigma = sigma))
  colnames(X) <- paste0("X_", 0:p)
  
  # Linear Predictor
  lin_pred <- X %*% t(b0)
  
  # Generate Y
  if (model == "correct") {
    Y <- rep(NA, N_total)
    S <- plogis(lin_pred)
    for (a in c(0,1)) {
      Y[A == a] <- rbinom(sum(A == a), 1, S[A == a, (a+1)])
    }
  }
  if (model == "misspecified 1") {
    Y <- rep(NA, N_total)
    S <- plogis(lin_pred + 0.2*(X[, 2])^2 -0.1*(X[, 3])^3 +0.1* X[, 5]*X[, 6])
    for (a in c(0,1)) {
      Y[A == a] <- rbinom(sum(A == a), 1, S[A == a, (a+1)])
    }
  }
  if (model == "misspecified 2") {
    Y <- rep(NA, N_total)
    for (a in c(0,1)) {
      S <- exp(-lin_pred[,a + 1]^2)
      Y[A == a] <- rbinom(sum(A == a), 1, S[A == a])
    }
  }
  if (model == "misspecified 3") {
    Y <- rep(NA, N_total)
    lin_pred2 <- X %*% t(b1)
    for (a in c(0,1)) {
      S <- plogis(2*tanh(lin_pred[,a + 1]) + 0.3*tanh(lin_pred2[,a + 1]) - 1)
      Y[A == a] <- rbinom(sum(A == a), 1, S[A == a])
    }
  }
  if (model == "misspecified 4") {
    Y <- rep(NA, N_total)
    for (a in c(0,1)) {
      Y[A == a] <-  ifelse(lin_pred[A==a, a + 1] > 0, 1, 0) 
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
