#' Covariance of diff between each pair of fragment ions based on shrinkage.
#'
#' A diff for an ion i is defined by d_{i} = xbar_con1{i} - xbar_con2{i}, where
#' xbar_con1 and xbar_con2 are the sample means for condition 1 and 2,
#' respectively. A covariance of diff between fragment ion i and j is
#'    cov(d_{i}, d_{j})
#'    = cov(xbar_con1{i}, xbar_con1{j}) + cov(xbar_con2{i}, xbar_con2{j})
#'
#' @param dat_con1,dat_con2  data matrices of normalized fragment ion peak area,
#'   each row is a replicate, each column is a fragment ion
#' @param lambda_var_con1,lambda_var_con2  the variance shrinkage intensities
#'   for con1 and con2, respectively. If not specified (default), they are
#'   estimated by corpcor::cov.shrink function.
#' @param cov_equal  TRUE if equal covariance is assumed and a pooled
#'   covariance matrix is computed.
#' @param verbose  TRUE to print messages. Default is FALSE.
#' @return covariance matrix of mean difference
#' @examples
#'   dat_x = matrix(rnorm(12), 4, 3)
#'   dat_y = matrix(rnorm(12), 4, 3)
#'   cov_diff_shrink(dat_x, dat_y)
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
    cov_con1 <- cov.shrink(dat_con1, lambda.var = lambda_var_con1,
                           verbose = verbose)
  )
  suppressWarnings(
    cov_con2 <- cov.shrink(dat_con2, lambda.var = lambda_var_con2,
                           verbose = verbose)
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
