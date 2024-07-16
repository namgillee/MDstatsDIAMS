#' Shrinkage-based t-test statistic
#'
#' @param dat_con1,dat_con2  data matrices of normalized fragment ion peak area,
#'   each row is a replicate, each column is a fragment ion
#' @param lambda_var_con1,lambda_var_con2  the variance shrinkage intensities
#'   for con1 and con2, respectively. If not specified (default), they are
#'   estimated by corpcor::cov.shrink function.
#' @param cov_equal  TRUE if equal covariance is assumed and a pooled
#'   covariance matrix is computed.
#' @param denom_eps  a small constant to be added to the denominator (i.e.,
#'   standard error of the numerator) of the statistic. This constant is useful
#'   for bootstrapping with a small sample size.
#' @return shrinkage-based t-test statistic
#' @examples
#'   dat_x = matrix(rnorm(12), 4, 3)
#'   dat_y = matrix(rnorm(12), 4, 3)
#'   shrinkage_t_test_statistic(dat_x, dat_y, TRUE, 0)
shrinkage_t_test_statistic <- function(
  dat_con1, dat_con2, cov_equal, denom_eps, lambda_var_con1, lambda_var_con2
) {

  result <- NaN

  if ((length(dim(dat_con1)) < 1) || (length(dim(dat_con2)) < 1)) {
    return(result)
  }

  # estimate is the sum of diff over fragment ions, and it constitutes
  # numerator of the statistic formula.
  estimate <- sum(
    apply(dat_con1, 2, mean, na.rm = TRUE) -
      apply(dat_con2, 2, mean, na.rm = TRUE),
    na.rm = TRUE
  )

  # covmat is the covariance matrix of diff between each pair of fragment ions
  covmat <- cov_diff_shrink(
    dat_con1, dat_con2, lambda_var_con1, lambda_var_con2,
    cov_equal = cov_equal, verbose = FALSE
  )

  # shrinkage statistic
  denom <- sqrt(sum(covmat, na.rm = TRUE))
  if (denom > 0) {
    result <- estimate / (denom + denom_eps)
  }

  return(result)
}


#' Shrinkage-based t-test for modified peptide using fragment ion peak area
#'
#' @param dat_con1,dat_con2  data matrices of normalized fragment ion peak area,
#'   each row is a replicate, each column is a fragment ion.
#' @param num_boot number of bootstrap replicates. Default is 200.
#' @param cov_equal  TRUE if equal covariance is assumed and a pooled
#'   covariance matrix is computed. Default is TRUE.
#' @param boot_denom_eps  a small constant to be added to the denominator (i.e.,
#'   standard error of the numerator) of the bootstrapped statistic. This
#'   constant is useful for bootstrapping with a small sample size.
#' @param verbose  TRUE to print messages. Default is FALSE.
#' @return list of statistic, df, p.value, and estimate.
#' @examples
#'   dat_x = matrix(rnorm(12), 4, 3)
#'   dat_y = matrix(rnorm(12), 4, 3)
#'   shrinkage_t_test(dat_x, dat_y)
shrinkage_t_test <- function(
  dat_con1, dat_con2, num_boot = 100, cov_equal = TRUE, boot_denom_eps = 0.5,
  verbose = FALSE
) {

  # check number of fragment ions
  if (ncol(dat_con1) != ncol(dat_con2)) {
    if (verbose) {
      print("Number of columns in dat_con1 and dat_con2 must be same.")
    }
    return(list())
  }

  # perform shrinkage t-test
  result <- list()

  result$statistic <- shrinkage_t_test_statistic(
    dat_con1, dat_con2, cov_equal, 0
  )

  shrinkage_t_statistic_wrapper <- function(dat, inds) {
    num_ions <- ncol(dat) / 2
    shrinkage_t_test_statistic(
      dat[inds, 1:num_ions],
      dat[inds, (num_ions + 1):(2 * num_ions)],
      cov_equal,
      boot_denom_eps
    )
  }
  boot_out <- boot(cbind(dat_con1, dat_con2),
                   shrinkage_t_statistic_wrapper,
                   R = num_boot)
  boot_var <- var(as.vector(boot_out$t), na.rm = TRUE)
  result$df <- ifelse(boot_var <= 1, Inf, 2 * boot_var / (boot_var - 1))
  result$p.value <- 2 * (1 - pt(abs(result$statistic), result$df))
  boot_mean <- mean(as.vector(boot_out$t), na.rm = TRUE)
  result$cv <- sqrt(boot_var) / boot_mean

  result$estimate <-
    sum(apply(dat_con1, 2, mean, na.rm = TRUE) -
          apply(dat_con2, 2, mean, na.rm = TRUE), na.rm = TRUE)

  return(result)
}
