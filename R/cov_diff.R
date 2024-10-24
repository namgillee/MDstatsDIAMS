#' Covariance of diff between each pair of fragment ions.
#'
#' A diff for an ion \eqn{i} is defined by
#' \eqn{d_i = \bar{x}_{1i} - \bar{x}_{2i}}, where
#' \eqn{\bar{x}_{1i}} and \eqn{\bar{x}_{1i}} are the sample means for
#' conditions 1 and 2, respectively. A covariance of diff between fragment ions
#' \eqn{i} and \eqn{j} is
#'    \deqn{cov(d_i, d_j)
#'    = cov(\bar{x}_{1i}, \bar{x}_{1j}) +
#'    cov(\bar{x}_{2i}, \bar{x}_{2j})}
#'
#' @param dat_con1,dat_con2  data matrices of normalized fragment ion peak area,
#' in logarithmic scale. Each row is a replicate, each column is a fragment ion.
#' @param verbose  TRUE to print messages. Default is FALSE.
#' @return covariance matrix of mean difference
#' @examples
#'   dat_x <- cbind(1 : 4, c(1, 3, 2, 4), 4 : 1)
#'   dat_y <- cbind(4 : 1, c(4, 2, 3, 1), 1 : 4)
#'   cov_diff(dat_x, dat_y)
#' @export
cov_diff <- function(dat_con1, dat_con2, verbose = FALSE) {
  # check number of fragment ions
  if (ncol(dat_con1) != ncol(dat_con2)) {
    if (verbose) {
      print("Number of columns in dat_con1 and dat_con2 must be same.")
    }
    return(matrix())
  }

  # compute covariance matrix
  num_frag <- ncol(dat_con1)
  result <- matrix(NA, num_frag, num_frag)
  for (i in 1:num_frag) {
    for (j in i:num_frag) {
      cov_conall <- c(
        var(dat_con1[, i], dat_con1[, j], na.rm = TRUE),
        var(dat_con2[, i], dat_con2[, j], na.rm = TRUE)
      )
      num_conall <- c(
        sum((!is.na(dat_con1[, i])) & (!is.na(dat_con1[, j]))),
        sum((!is.na(dat_con2[, i])) & (!is.na(dat_con2[, j])))
      )
      if (all(!is.na(cov_conall))) {
        result[i, j] <- result[j, i] <- sum(cov_conall / num_conall)
      }
    }
  }
  return(result)
}
