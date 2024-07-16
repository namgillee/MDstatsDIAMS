#' Independent (two) sample t-test for modified peptide using aggregated
#' quantity
#'
#' @param quantity1,quantity2  numeric vectors of normalized peptide quantity
#' @param var_equal  TRUE to assume equal variance between conditions. Default
#'   is TRUE.
#' @param conf_level  Confidence level. Default is 0.95.
#' @param verbose  TRUE to print messages. Default is FALSE.
#' @return list of statistic, df, p.value, and estimate.
#' @examples
#' independent_t_test(rnorm(4), rnorm(4))
independent_t_test <- function(quantity1, quantity2, var_equal = TRUE,
                               conf_level = 0.95, verbose = FALSE) {
  # missing values
  is_na_x <- is.na(quantity1)
  is_na_y <- is.na(quantity2)
  if ((sum(!is_na_x) < 3) || (sum(!is_na_y) < 3)) {
    if (verbose) {
      print("Input vectors must have at least 3 entries.")
    }
    return(list())
  }

  # perform independent sample t-test
  ttest_out <- t.test(
    quantity1[!is_na_x], quantity2[!is_na_y], var.equal = var_equal,
    paired = FALSE, conf.level = conf_level
  )
  result <- list(
    statsitic = ttest_out$statistic,
    df = ttest_out$parameter,
    p.value = ttest_out$p.value,
    estimate = ttest_out$estimate[1] - ttest_out$estimate[2]
  )
  return(result)
}
