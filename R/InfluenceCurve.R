# Purpose: Influence functions for fairness estimation
# Updated: 2024-03-01

#' @param Y Outcome in labeled dataset.
#' @param S Model score in labeled dataset.
#' @param A Group indicator in labeled dataset.
#' @param threshold Threshold for classification based on the model score.
#' Default value is 0.5.


Influence_curve <- function(pest, Y, S, A, threshold = 0.5, method) {
  class <- sort(unique(A))
  out <- pest
  if (method == "supervised") {
    for (i in class) {
      D <- 1 * (S > threshold)
      C <- 1 * (A == i)
      mu_Y <- mean(Y[A == i])
      mu_D <- mean(D[A == i])
      rho <- sum(A == i) / length(A)

      # influence curve for TPR & FNR
      out[out$Metric == "TPR", paste0("Group", i)] <- out[
        out$Metric == "FNR",
        paste0("Group", i)
      ] <- (rho)^(-2) * (mu_Y)^(-2) *
        sum(Y * (D - pest[pest$Metric == "TPR", paste0("Group", i)])**2 * C) /
        (length(Y)**2)

      # influence curve for FPR & TNR
      out[out$Metric == "FPR", paste0("Group", i)] <- out[
        out$Metric == "TNR",
        paste0("Group", i)
      ] <- (rho)^(-2) *
        (1 - mu_Y)^(-2) * sum((1 - Y) * (D - pest[
          pest$Metric == "FPR",
          paste0("Group", i)
        ])**2 * C) / (length(Y)**2)

      # influence curve for PPV and NPV
      out[out$Metric == "PPV", paste0("Group", i)] <- out[
        out$Metric == "NPV",
        paste0("Group", i)
      ] <- (rho)^(-2) * (mu_D)^(-2) *
        sum(D * (Y - pest[pest$Metric == "PPV", paste0("Group", i)])**2 * C) /
        (length(Y)**2)

      # influence curve for F1
      out[out$Metric == "F1", paste0("Group", i)] <-
        (rho)^(-2) * (mu_Y + mu_D)**(-2) *
          sum((D * (Y - pest[pest$Metric == "F1", paste0("Group", i)]) +
            Y * (D - pest[pest$Metric == "F1", paste0("Group", i)]))**2 * C) /
          (length(Y)**2)

      # influence curve for ACC
      out[out$Metric == "ACC", paste0("Group", i)] <-
        (rho)^(-2) * sum((1 - (Y - D)^2 -
          pest[pest$Metric == "ACC", paste0("Group", i)])**2 * C) /
          (length(Y)**2)

      # influence curve for BS
      out[out$Metric == "BS", paste0("Group", i)] <-
        (rho)^(-2) * sum(((S - Y)**2 - pest[
          pest$Metric == "BS",
          paste0("Group", i)
        ])**2 * C) /
          (length(Y)**2)
    }
    out[,"Delta"] <- rowSums(out[,paste0("Group", class)])
  }
  return(out)
}
