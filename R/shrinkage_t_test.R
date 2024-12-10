#' Standard error of estimate for shrinkage-based t-test statistic.
#'
#' The estimate is the sum of diff over fragment ions, where diff for a fragment
#' ion is defined by \eqn{d_i = \bar{x}_{1i} - \bar{x}_{2i}}, and it
#' constitutes the numerator of the test statistic formula.
#' @param dat_con1,dat_con2  Data matrices of normalized fragment ion peak area
#'   in logarithmic scale. Each row is a replicate, each column is a fragment
#'   ion.
#' @param lambda_var  Variance shrinkage intensities. If not specified
#'   (default), they are estimated by corpcor::cov.shrink function.
#' @param cov_equal  TRUE if equal covariance is assumed and a pooled
#'   covariance matrix is computed.
#' @param verbose  TRUE to print messages. Default is FALSE.
#' @return standard error of the sum of mean differences
se_diff_shrink <- function(
  dat_con1, dat_con2, lambda_var, cov_equal = FALSE, verbose = FALSE
) {

  # Check number of replicates
  if (nrow(dat_con1) != nrow(dat_con2)) {
    if (verbose) {
      print("Number of rows in dat_con1 and dat_con2 must be same.")
    }
    return(matrix())
  }

  # Number of ions, I
  num_ions <- ncol(dat_con1)

  ## 1) The (2I x 2I) matrix of covariance terms between each pair of fragment
  ##    ions where replicates are equal.

  # Combine data by column
  dat <- cbind(dat_con1, dat_con2)

  # Remove rows with any missing value
  dat <- na.omit(dat)
  num_replicates <- nrow(dat)
  if (num_replicates < 3) {
    if (verbose) {
      print("Number of complete rows in the data matrix must be larger than 2")
    }
    return(matrix())
  }

  # Compute covariance matrix
  suppressWarnings(
    covmat <- corpcor::cov.shrink(
      dat, lambda.var = lambda_var, verbose = verbose
    )
  )

  # Extract sub-matrices, A1, A2, B
  cov_con1 <- covmat[1:num_ions, 1:num_ions]
  cov_con2 <- covmat[(num_ions + 1):ncol(covmat), (num_ions + 1):ncol(covmat)]
  cov_unequal_conditions <- covmat[1:num_ions, (num_ions + 1):ncol(covmat)]

  if (cov_equal) {
    cov_pooled <- (cov_con1 + cov_con2) / 2
    cov_con1 <- cov_pooled
    cov_con2 <- cov_pooled
  }

  ## 2) The covariance term between a pair of fragment ions where replicates are
  ##    unequal.

  dat_average <- (dat[, 1:num_ions] + dat[, (num_ions + 1):ncol(dat)]) / 2
  dat_con1_centered <- dat[, 1:num_ions] - dat_average
  dat_con2_centered <- dat[, (num_ions + 1):ncol(dat)] - dat_average

  cov_con1_unequal_replicates <- cov_con2_unequal_replicates <-
    matrix(0, num_ions, num_ions)
  for (shift in 1:(num_replicates - 1)) {
    idx_rows_shifted <- c((shift + 1):num_replicates, 1:shift)

    # sum_{s > 0} sum_{r} x_{r} x_{r + s}
    cov_con1_unequal_replicates <- cov_con1_unequal_replicates +
      t(dat_con1_centered) %*% dat_con1_centered[idx_rows_shifted, ]
    cov_con2_unequal_replicates <- cov_con2_unequal_replicates +
      t(dat_con2_centered) %*% dat_con2_centered[idx_rows_shifted, ]
  }

  cov_unequal_replicates <- (
    sum(cov_con1_unequal_replicates) + sum(cov_con2_unequal_replicates)
  ) / num_ions / num_ions / num_replicates / (num_replicates - 1)

  # Take an absolute value to prevent a negative value
  cov_unequal_replicates <- abs(cov_unequal_replicates)

  ## 3) standard error
  # Take an absolute value to prevent a negative value inside the sqrt()
  standard_error <- sqrt(
    abs(
      sum(cov_con1) / num_replicates + sum(cov_con2) / num_replicates -
        2 * sum(cov_unequal_conditions) / num_replicates +
        2 * num_ions * num_ions * (num_replicates - 1) / num_replicates *
          cov_unequal_replicates
    )
  )

  ## Set attributes
  attributes(standard_error) <- NULL

  attr(standard_error, "lambda") <- attr(covmat, "lambda")
  attr(standard_error, "lambda_estimated") <- attr(covmat, "lambda.estimated")
  attr(standard_error, "lambda_var") <- attr(covmat, "lambda.var")
  attr(standard_error, "lambda_var_estimated") <-
    attr(covmat, "lambda.var.estimated")

  return(standard_error)
}


#' Shrinkage-based t-test statistic
#'
#' @param dat_con1,dat_con2  Data matrices of normalized fragment ion peak area
#'   in logarithmic scale. Each row is a replicate, and each column is a
#'   fragment ion.
#' @param lambda_var  Variance shrinkage intensity. If not specified (default),
#'   they are estimated by corpcor::cov.shrink function.
#' @param cov_equal  TRUE if equal covariance is assumed and a pooled
#'   covariance matrix is computed.
#' @param denom_eps  A small constant to be added to the denominator (i.e.,
#'   standard error of the numerator) of the statistic. This constant is useful
#'   for bootstrapping with a small sample size.
#' @return shrinkage-based t-test statistic
shrinkage_t_test_statistic <- function(
  dat_con1, dat_con2, cov_equal, denom_eps, lambda_var
) {

  result <- NaN

  if ((length(dim(dat_con1)) < 1) || (length(dim(dat_con2)) < 1)) {
    return(result)
  }

  # estimate is the sum of diff over fragment ions, where diff for a fragment
  # ion is defined by \eqn{d_i = \bar{x}_{1i} - \bar{x}_{2i}}, and it
  # constitutes the numerator of the test statistic formula.
  estimate <- sum(
    apply(dat_con1, 2, mean, na.rm = TRUE) -
      apply(dat_con2, 2, mean, na.rm = TRUE),
    na.rm = TRUE
  )

  # denominator of the test statistic formula.
  denom <- se_diff_shrink(
    dat_con1, dat_con2, lambda_var,
    cov_equal = cov_equal, verbose = FALSE
  )

  # test statistic
  if (denom > 0) {
    result <- estimate / (denom + denom_eps)
  }

  return(result)
}


#' Shrinkage-based t-test for testing group mean differeneces in peptide mean
#' quantity using fragment ion peak area
#'
#' @param dat_con1,dat_con2  data matrices of normalized fragment ion peak area
#'   in logarithmic scale. Each row is a replicate, each column is a fragment
#'   ion.
#' @param num_boot number of bootstrap replicates. Default is 200.
#' @param cov_equal  TRUE if equal covariance is assumed and a pooled
#'   covariance matrix is computed. Default is TRUE.
#' @param boot_denom_eps  a small constant, boot_denom_eps * (R / 4) ^ (-3 / 2),
#'   will be added to the denominator (i.e., standard error of the numerator) of
#'   the bootstrapped statistic. This constant is useful for bootstrapping with
#'   a small sample size.
#' @param verbose  TRUE to print messages. Default is FALSE.
#' @return list of statistic, df, p.value, and estimate.
#' @examples
#'   dat_x <- cbind(1 : 4, c(1, 3, 2, 4), 4 : 1)
#'   dat_y <- cbind(4 : 1, c(4, 2, 3, 1), 1 : 4)
#'   shrinkage_t_test(dat_x, dat_y)
#' @export
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

  # check number of replicates
  if (nrow(dat_con1) != nrow(dat_con2)) {
    if (verbose) {
      print("Number of rows in dat_con1 and dat_con2 must be same.")
    }
    return(list())
  }

  # perform shrinkage t-test
  result <- list()

  result$statistic <- shrinkage_t_test_statistic(
    dat_con1, dat_con2, cov_equal, 0
  )

  shrinkage_t_statistic_wrapper <- function(dat, inds) {
    num_ions <- ceiling(ncol(dat) / 2)
    shrinkage_t_test_statistic(
      dat[inds, 1:num_ions],
      dat[inds, (num_ions + 1):ncol(dat)],
      cov_equal,
      boot_denom_eps * (nrow(dat) / 4)^(-3 / 2)
    )
  }
  boot_out <- boot::boot(
    cbind(dat_con1, dat_con2), shrinkage_t_statistic_wrapper, R = num_boot
  )
  boot_var <- var(as.vector(boot_out$t), na.rm = TRUE)
  result$df <- ifelse(boot_var <= 1, Inf, 2 * boot_var / (boot_var - 1))
  result$p.value <- 2 * (1 - pt(abs(result$statistic), result$df))
  boot_mean <- mean(as.vector(boot_out$t), na.rm = TRUE)
  result$cv <- sqrt(boot_var) / boot_mean

  result$estimate <-
    sum(
      apply(dat_con1, 2, mean, na.rm = TRUE) -
        apply(dat_con2, 2, mean, na.rm = TRUE), na.rm = TRUE
    )

  return(result)
}
