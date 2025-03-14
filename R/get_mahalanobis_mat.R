#' Create rank based Mahalanobis distance matrix
#'
#' @param X Matrix or dataframe to calculate distances with
#'
#' @return Rank based Mahalanobis distance matrix
#' @export
#'
#' @examples
#' get_mahalanobis_mat(mtcars)
get_mahalanobis_mat <- function(X) {
  X <- as.matrix(X)
  n <- dim(X)[1]
  k <- dim(X)[2]
  for (j in 1:k) {
    X[, j] <- rank(X[, j]) # compute on ranks
  }
  cv <- stats::cov(X)
  vuntied <- stats::var(1:n)
  rat <- sqrt(vuntied / diag(cv))
  cv <- diag(rat) %*% cv %*% diag(rat)
  out <- matrix(NA, n, n)
  icov <- MASS::ginv(cv)
  for (i in 1:n) {
    out[i, ] <- stats::mahalanobis(X, X[i, ], icov, inverted = TRUE)
  }
  out
}
