max_min_diff <- function(x) {
  return(max(x) - min(x))
}

max_min_ratio <- function(x) {
  return(max(x) / min(x))
}

max_abs_diff <- function(x) {
  return(max(x - mean(x)))
}

mad <- function(x) {
  return(mean(abs(x - mean(x))))
}

entropy <- function(x, alpha = 2) {
  k <- length(x)
  return(1 / (k * alpha * (alpha - 1)) * sum((x / mean(x))^(alpha) - 1))
}
