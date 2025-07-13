#' Standard error of estimate for shrinkage-based t-test statistic.
#'
#' Standard error appears in the the denominator of the test statistic formula.
#' In the shrinkage-based t-test, the estimate of group mean difference is the
#' sum of \eqn{d_i} across fragment ions \eqn{i}, where
#' \eqn{d_i = \bar{x}_{1i} - \bar{x}_{2i}},
#' with \eqn{\bar{x}_{1i}} and \eqn{\bar{x}_{2i}} representing the mean
#' quantities for fragment ion \eqn{i} in groups 1 and 2, respectively.
#' @param dat_con1,dat_con2  Data matrices of normalized fragment ion peak area
#'   in logarithmic scale, \eqn{x_{1ir}^p, x_{2ir}^p}. Each row is a replicate,
#'   each column is a fragment ion.
#' @param lambda_var  Variance shrinkage intensities. If not specified
#'   (default), they are estimated by corpcor::cov.shrink function.
#' @param cov_unequal_replicates  covariance of input variables whose replicates
#'   are not equal, \eqn{Cov(x_{ci_1r_1}^p, x_{ci_2r_2}^p), r_1 \neq r_2}.
#' @param cov_equal  TRUE if equal covariance is assumed and a pooled
#'   covariance matrix is computed.
#' @param verbose  TRUE to print messages. Default is FALSE.
#' @return standard error of the sum of mean differences
se_diff_shrink <- function(
  dat_con1, dat_con2, lambda_var, cov_unequal_replicates = 0,
  cov_equal = TRUE, verbose = FALSE
) {

  # Check number of replicates
  if (nrow(dat_con1) != nrow(dat_con2)) {
    if (verbose) {
      print("Number of rows in dat_con1 and dat_con2 must be same.")
    }
    return(NaN)
  }

  # Number of ions, I
  num_ions <- ncol(dat_con1)

  ## 1) The (2I x 2I) matrix of covariance terms between each pair of fragment
  ##    ions where replicates are equal.

  # Combine data by column
  dat <- cbind(dat_con1, dat_con2)

  # Remove rows with any missing value
  dat <- stats::na.omit(dat)
  num_replicates <- nrow(dat)
  if (num_replicates < 3) {
    if (verbose) {
      print("Number of complete rows in the data matrix must be larger than 2")
    }
    return(NaN)
  }

  # Compute covariance matrix
  # Warnings are suppressed. The function corpcor::cov.shrink() raises a warning
  # when it detects "any variable with zero scale" in the input data matrix--
  # that is, when a column is constant. Nevertheless, the shrinkage estimation
  # method still returns valid covariance estimates with finite values, avoiding
  # NA or NaN entries.
  suppressWarnings(
    covmat <- corpcor::cov.shrink(
      dat, lambda.var = lambda_var, verbose = verbose
    )
  )

  # Extract sub-matrices: A1, A2, B
  cov_con1 <- covmat[1:num_ions, 1:num_ions, drop = FALSE]
  cov_con2 <- covmat[(num_ions + 1):ncol(covmat), (num_ions + 1):ncol(covmat),
                     drop = FALSE]
  cov_unequal_conditions <- covmat[1:num_ions, (num_ions + 1):ncol(covmat),
                                   drop = FALSE]

  if (cov_equal) {
    cov_pooled <- (cov_con1 + cov_con2) / 2
    cov_con1 <- cov_pooled
    cov_con2 <- cov_pooled
  }

  ## 2) standard error
  # Take an absolute value to prevent a negative value
  standard_error <- sqrt(
    abs(
      sum(cov_con1 + cov_unequal_replicates) / num_replicates +
        sum(cov_con2 + cov_unequal_replicates) / num_replicates -
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

  attr(standard_error, "a") <- sum(cov_con1 + cov_unequal_replicates) /
    num_replicates + sum(cov_con2 + cov_unequal_replicates) / num_replicates
  attr(standard_error, "a0") <- sum(diag(cov_con1 + cov_unequal_replicates)) /
    num_replicates +
    sum(diag(cov_con2 + cov_unequal_replicates)) / num_replicates
  attr(standard_error, "a_minus_d") <- sum(cov_con1) / num_replicates +
    sum(cov_con2) / num_replicates
  attr(standard_error, "a0_minus_d") <- sum(diag(cov_con1)) / num_replicates +
    sum(diag(cov_con2)) / num_replicates
  attr(standard_error, "b") <- - 2 * sum(cov_unequal_conditions) /
    num_replicates
  attr(standard_error, "b0") <- - 2 * sum(diag(cov_unequal_conditions)) /
    num_replicates
  attr(standard_error, "d") <- 2 * num_ions * num_ions * (num_replicates - 1) /
    num_replicates * cov_unequal_replicates

  return(standard_error)
}


#' Shrinkage-based t-test statistic
#'
#' @param dat_con1,dat_con2  Data matrices of normalized fragment ion peak area
#'   in logarithmic scale. Each row is a replicate, and each column is a
#'   fragment ion.
#' @param lambda_var  Variance shrinkage intensity. If not specified (default),
#'   they are estimated by corpcor::cov.shrink function.
#' @param cov_unequal_replicates  covariance of input variables whose replicates
#'   are not equal.
#' @param cov_equal  TRUE if equal covariance is assumed and a pooled
#'   covariance matrix is computed.
#' @param denom_eps  A small constant to be added to the denominator (i.e.,
#'   standard error of the numerator) of the statistic. This constant is useful
#'   for bootstrapping with a small sample size.
#' @return shrinkage-based t-test statistic
shrinkage_t_test_statistic <- function(
  dat_con1, dat_con2, cov_unequal_replicates, cov_equal, denom_eps, lambda_var
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
    cov_unequal_replicates = cov_unequal_replicates,
    cov_equal = cov_equal, verbose = FALSE
  )

  # test statistic
  if (!is.na(denom) && (denom >= 0)) {
    result <- estimate / (denom + denom_eps)
  }

  attributes(result) <- attributes(denom)

  return(result)
}


#' Shrinkage-based t-test for testing group mean differeneces in peptide mean
#' quantity using fragment ion peak area
#'
#' @param dat_con1,dat_con2  data matrices of normalized fragment ion peak area
#'   in logarithmic scale. Each row is a replicate, each column is a fragment
#'   ion.
#' @param cov_unequal_replicates  covariance of input variables whose replicates
#'   are not equal.
#' @param num_boot number of bootstrap replicates.
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
  dat_con1, dat_con2, cov_unequal_replicates = 0, num_boot = 100,
  cov_equal = TRUE, boot_denom_eps = 0.3, verbose = FALSE
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
    dat_con1, dat_con2, cov_unequal_replicates, cov_equal, 0
  )

  shrinkage_t_statistic_wrapper <- function(dat, inds) {
    num_ions <- ceiling(ncol(dat) / 2)
    shrinkage_t_test_statistic(
      dat[inds, 1:num_ions],
      dat[inds, (num_ions + 1):ncol(dat)],
      cov_unequal_replicates,
      cov_equal,
      boot_denom_eps * (nrow(dat) / 4)^(-3 / 2)
    )
  }
  boot_out <- boot::boot(
    cbind(dat_con1, dat_con2), shrinkage_t_statistic_wrapper, R = num_boot
  )
  boot_var <- stats::var(as.vector(boot_out$t), na.rm = TRUE)
  result$df <- ifelse(boot_var <= 1, Inf, 2 * boot_var / (boot_var - 1))
  result$p.value <- 2 * (1 - stats::pt(abs(result$statistic), result$df))
  boot_mean <- mean(as.vector(boot_out$t), na.rm = TRUE)
  result$cv <- sqrt(boot_var) / boot_mean

  result$estimate <-
    sum(
      apply(dat_con1, 2, mean, na.rm = TRUE) -
        apply(dat_con2, 2, mean, na.rm = TRUE), na.rm = TRUE
    )

  result$lambda <- attr(result$statistic, "lambda")
  result$lambda_var <- attr(result$statistic, "lambda_var")
  result$a <- attr(result$statistic, "a")
  result$a0 <- attr(result$statistic, "a0")
  result$a_minus_d <- attr(result$statistic, "a_minus_d")
  result$a0_minus_d <- attr(result$statistic, "a0_minus_d")
  result$b <- attr(result$statistic, "b")
  result$b0 <- attr(result$statistic, "b0")
  result$d <- attr(result$statistic, "d")

  return(result)
}


#' Wrapper function for shrinkage t-test
#'
#' @param groupdf A subset of precursor report which consists of two conditions.
#'   It contains columns named "replicate", "fragment_id", "condition", and
#'   a value column.
#' @param value_column Column name for log10-transformed fragment ion quantities
#' @param cov_unequal_replicates_column Column name for covariance of input
#'   variables whose replicates are not equal. All entries in the column is
#'   considered equal. If NULL, the covariance is set to zero.
#' @param boot_denom_eps A parameter for shrinkage t-test
#' @importFrom dplyr %>%
#' @return A data frame of the shrinkage t-test results
compute_shrink_on_group <- function(
  groupdf,
  value_column = "log10_fragment_peak_area",
  cov_unequal_replicates_column = NULL,
  boot_denom_eps = 0.3
) {
  if (length(unique(groupdf$condition)) < 2) {
    return(as.data.frame(list()))
  }

  cov_unequal_replicates <- 0
  if (!is.null(cov_unequal_replicates_column) &&
        !is.null(groupdf[[cov_unequal_replicates_column]])) {
    cov_unequal_replicates <- groupdf[[cov_unequal_replicates_column]][
      !is.na(groupdf[[cov_unequal_replicates_column]])
    ][1]
  }

  df_shrink <- groupdf %>%
    reshape::cast(
      replicate ~ fragment_id ~ condition,
      value = value_column,
      fun.aggregate = mean
    )

  dat_con1 <- as.matrix(df_shrink[, , 1])
  dat_con2 <- as.matrix(df_shrink[, , 2])
  dim(dat_con1) <- dim(dat_con2) <- dim(df_shrink)[1 : 2]

  as.data.frame(shrinkage_t_test(
    dat_con1, dat_con2, cov_unequal_replicates = cov_unequal_replicates,
    num_boot = 100, cov_equal = TRUE, boot_denom_eps = boot_denom_eps,
    verbose = FALSE
  ))
}


#' Run shrinkage t-test comparing two conditions on a standard report
#'
#' @param report Fragment ion report with the columns experiment, condition,
#'   replicate, protein_id, precursor_id, fragment_id, and fragment_peak_area.
#'   The condition column is supposed to consist of two conditions.
#' @param cov_unequal_replicates_column Column name for covariance of input
#'   variables whose replicates are not equal. If NULL, the covariance is set to
#'   zero.
#' @param boot_denom_eps A parameter for shrinkage t-test
#' @return A data frame of the test results. Columns are shrinkage_t_test
#'   results, and rows are precursors.
#' @export
compute_shrink_on_stdreport <- function(
  report,
  cov_unequal_replicates_column = "cov_unequal_replicates",
  boot_denom_eps = 0.3
) {
  report[["log10_fragment_peak_area"]] <- log10(report[["fragment_peak_area"]])
  if (!is.null(cov_unequal_replicates_column) &&
        is.null(report[[cov_unequal_replicates_column]])) {
    warning(
      paste(
        "The",
        cov_unequal_replicates_column,
        "column does not exist in the report."
      )
    )
  }

  result_shrink0 <- report %>%
    dplyr::group_by(
      .data$experiment,
      .data$protein_id,
      .data$precursor_id
    ) %>%
    dplyr::group_modify(
      ~compute_shrink_on_group(
        .x,
        value_column = "log10_fragment_peak_area",
        cov_unequal_replicates_column = cov_unequal_replicates_column,
        boot_denom_eps = boot_denom_eps
      )
    ) %>%
    as.data.frame()

  return(result_shrink0)
}
