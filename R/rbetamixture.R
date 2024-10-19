#' Generate random samples from a mixture of Beta distributions
#'
#' @param n  number of observations
#' @param shape1s,shape2s  vectors of non-negative parameters of component
#' beta distributions. Length of each vector is the number of component beta
#' distributions.
#' @return a length-n vector of generated values
#' @export
rbetamixture <- function(n, shape1s = 1, shape2s = 1) {
  # Number of beta components
  n_beta1 <- length(shape1s)
  n_beta2 <- length(shape2s)
  if ((n_beta1 == 1) && (n_beta2 > 1)) {
    shape1s <- rep(shape1s, n_beta2)
    n_beta1 <- n_beta2
  }
  if ((n_beta1 > 1) && (n_beta2 == 1)) {
    shape2s <- rep(shape2s, n_beta1)
    n_beta2 <- n_beta1
  }
  if (n_beta1 != n_beta2) {
    stop("Lengths of shape1s and shape2s must equal.")
  }

  # Select beta components
  id_components <- sample(n_beta1, size = n, replace = TRUE)

  # Sample from the selected beta components
  out <- rep(-1, n)
  for (i in 1 : n) {
    id <- id_components[i]
    out[i] <- rbeta(1, shape1s[id], shape2s[id])
  }

  return(out)
}
