#' Paired t-test for modified peptide using fragment ion peak area
#'
#' @param quantity1,quantity2  numeric vectors of normalized fragment peak area
#' in logarithmic scale
#' @param conf_level  Confidence level. Default is 0.95.
#' @param verbose  TRUE to print messages. Default is FALSE.
#' @return list of statistic, df, p.value, and estimate.
#' @examples
#' paired_t_test(1 : 4, c(1, 3, 2, 4))
#' @export
paired_t_test <- function(
  quantity1,
  quantity2,
  conf_level = 0.95,
  verbose = FALSE
) {
  # check input length
  if (length(quantity1) != length(quantity2)) {
    print("Input vectors must have the same length.")
    return(list())
  }

  # missing values
  is_na <- is.na(quantity1) | is.na(quantity2)
  if (sum(!is_na) < 3) {
    if (verbose) {
      print("Input vectors must have at least 3 entries.")
    }
    return(list())
  }

  # perform paired t-test
  ttest_out <- t.test(
    quantity1[!is_na], quantity2[!is_na], paired = TRUE, conf.level = conf_level
  )
  result <- list(
    statistic = ttest_out$statistic,
    df = ttest_out$parameter,
    p.value = ttest_out$p.value,
    estimate = ttest_out$estimate
  )
  return(result)
}
