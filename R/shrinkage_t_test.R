#' Covariance of diff between each pair of fragment ions based on shrinkage.
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
#' @param dat_con1,dat_con2  data matrices of normalized fragment ion peak area
#'   in logarithmic scale. Each row is a replicate, each column is a fragment
#'   ion.
#' @param lambda_var_con1,lambda_var_con2  the variance shrinkage intensities
#'   for con1 and con2, respectively. If not specified (default), they are
#'   estimated by corpcor::cov.shrink function.
#' @param cov_equal  TRUE if equal covariance is assumed and a pooled
#'   covariance matrix is computed.
#' @param verbose  TRUE to print messages. Default is FALSE.
#' @return covariance matrix of mean difference
cov_diff_shrink <- function(
  dat_con1, dat_con2, lambda_var_con1, lambda_var_con2, cov_equal = FALSE,
  verbose = FALSE
) {

  # check number of fragment ions
  if (ncol(dat_con1) != ncol(dat_con2)) {
    if (verbose) {
      print("Number of columns in dat_con1 and dat_con2 must be same.")
    }
    return(matrix())
  }

  # remove rows with any missing value
  dat_con1 <- na.omit(dat_con1)
  num_con1 <- nrow(dat_con1)
  dat_con2 <- na.omit(dat_con2)
  num_con2 <- nrow(dat_con2)
  if ((num_con1 < 3) || (num_con2 < 3)) {
    if (verbose) {
      print("Number of rows in each data matrix must be larger than 2")
    }
    return(matrix())
  }

  # compute covariance matrix
  suppressWarnings(
    cov_con1 <- corpcor::cov.shrink(
      dat_con1, lambda.var = lambda_var_con1, verbose = verbose
    )
  )
  suppressWarnings(
    cov_con2 <- corpcor::cov.shrink(
      dat_con2, lambda.var = lambda_var_con2, verbose = verbose
    )
  )
  if (cov_equal) {
    const <- 1 / (num_con1 + num_con2 - 2) * (1 / num_con1 + 1 / num_con2)
    result <- ((num_con1 - 1) * cov_con1 + (num_con2 - 1) * cov_con2) * const
  } else {
    result <- cov_con1 / num_con1 + cov_con2 / num_con2
  }
  attributes(result) <- NULL
  dim(result) <- dim(cov_con1)

  attr(result, "lambda_con1") <- attr(cov_con1, "lambda")
  attr(result, "lambda_con1.estimated") <- attr(cov_con1, "lambda.estimated")
  attr(result, "lambda_var_con1") <- attr(cov_con1, "lambda.var")
  attr(result, "lambda_var_con1.estimated") <-
    attr(cov_con1, "lambda.var.estimated")
  attr(result, "lambda_con2") <- attr(cov_con2, "lambda")
  attr(result, "lambda_con2.estimated") <- attr(cov_con2, "lambda.estimated")
  attr(result, "lambda_var_con2") <- attr(cov_con2, "lambda.var")
  attr(result, "lambda_var_con2.estimated") <-
    attr(cov_con2, "lambda.var.estimated")
  return(result)
}


#' Shrinkage-based t-test statistic
#'
#' @param dat_con1,dat_con2  data matrices of normalized fragment ion peak area
#'   in logarithmic scale. Each row is a replicate, each column is a fragment
#'   ion.
#' @param lambda_var_con1,lambda_var_con2  the variance shrinkage intensities
#'   for con1 and con2, respectively. If not specified (default), they are
#'   estimated by corpcor::cov.shrink function.
#' @param cov_equal  TRUE if equal covariance is assumed and a pooled
#'   covariance matrix is computed.
#' @param denom_eps  a small constant to be added to the denominator (i.e.,
#'   standard error of the numerator) of the statistic. This constant is useful
#'   for bootstrapping with a small sample size.
#' @return shrinkage-based t-test statistic
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
