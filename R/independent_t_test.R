#' Independent samples t-test for precursor peptide using aggregated
#' quantity
#'
#' @param quantity1,quantity2  numeric vectors of normalized peptide quantity in
#' logarithmic scale
#' @param var_equal  TRUE to assume equal variance between conditions. Default
#'   is TRUE.
#' @param conf_level  Confidence level. Default is 0.95.
#' @param verbose  TRUE to print messages. Default is FALSE.
#' @return list of statistic, df, p.value, and estimate.
#' @examples
#' independent_t_test(1 : 4, c(1, 3, 2, 4))
#' @export
independent_t_test <- function(
  quantity1,
  quantity2,
  var_equal = TRUE,
  conf_level = 0.95,
  verbose = FALSE
) {
  # missing values
  is_na_x <- is.na(quantity1)
  is_na_y <- is.na(quantity2)
  if ((sum(!is_na_x) < 3) || (sum(!is_na_y) < 3)) {
    if (verbose) {
      print("Input vectors must have at least 3 entries.")
    }
    return(list())
  }

  # perform independent samples t-test
  ttest_out <- stats::t.test(
    quantity1[!is_na_x], quantity2[!is_na_y], var.equal = var_equal,
    paired = FALSE, conf.level = conf_level
  )
  result <- list(
    statistic = ttest_out$statistic,
    df = ttest_out$parameter,
    p.value = ttest_out$p.value,
    estimate = c(
      `mean difference` = unname(ttest_out$estimate[1] - ttest_out$estimate[2])
    )
  )
  return(result)
}


#' Wrapper function for independent samples t-test
#'
#' @param groupdf A subset of precursor report which contains two columns for
#' log10-transformed precursor quantities of conditions x and y
#' @param column_x column name for x
#' @param column y column name for y
#' @return a data frame of the independent samples t-test results
compute_indep_on_group <- function(
  groupdf,
  column_x = "Log10Quantity.x",
  column_y = "Log10Quantity.y"
) {
  as.data.frame(
    independent_t_test(
      groupdf[[column_x]],
      groupdf[[column_y]],
      verbose = FALSE
    )
  )
}


#' Run independent samples t-test comparing two conditions on a standard report
#'
#' @param report fragment ion report with the columns experiment, condition,
#'   replicate, protein_id, precursor_id, and fragment_peak_area.
#'   The condition column is supposed to consist of two conditions.
#' @return a data frame of the test results. Columns are independent_t_test
#'   results, and rows are precursors.
#' @export
compute_indep_on_stdreport <- function(report) {
  df_indep <- report %>%
    dplyr::group_by(
      .data$experiment,
      .data$protein_id,
      .data$precursor_id,
      .data$replicate,
      .data$condition
    ) %>%
    dplyr::summarise(
      log10_peptide_quantity = log10(
        sum(
          .data$fragment_peak_area,
          na.rm = TRUE
        )
      )
    ) %>%
    reshape::cast(
      experiment + protein_id + precursor_id + replicate ~ condition,
      value = "log10_peptide_quantity",
      fun.aggregate = mean,
      na.rm = TRUE
    )

  colnames(df_indep)[c(5, 6)] <- c(
    "Log10Quantity.x", "Log10Quantity.y"
  )

  result_indep0 <- df_indep %>%
    dplyr::group_by(
      .data$experiment,
      .data$protein_id,
      .data$precursor_id
    ) %>%
    dplyr::group_modify(~compute_indep_on_group(.x)) %>%
    as.data.frame()

  return(result_indep0)
}
