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
  ttest_out <- stats::t.test(
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


#' Wrapper function for paired t-test
#'
#' @param groupdf A subset of fragment ion report which contains two columns for
#' log10-transformed fragment ion quantities of conditions x and y
#' @param column_x column name for x
#' @param column y column name for y
#' @return a data frame of the paired t-test results
compute_paired_on_group <- function(
  groupdf,
  column_x = "Log10NormalizedPeakArea.x",
  column_y = "Log10NormalizedPeakArea.y"
) {
  as.data.frame(
    paired_t_test(
      groupdf[[column_x]],
      groupdf[[column_y]],
      verbose = FALSE
    )
  )
}


#' Run paired t-test comparing two conditions on a standard report
#'
#' @param report fragment ion report with the columns experiment, condition,
#'   replicate, protein_id, precursor_id, fragment_id, and fragment_peak_area.
#'   The condition column is supposed to consist of two conditions.
#' @return a data frame of the test results. Columns are paired_t_test
#'   results, and rows are precursors.
#' @export
compute_paired_on_stdreport <- function(report) {
  report[["log10_fragment_peak_area"]] <- log10(report[["fragment_peak_area"]])

  df_paired <- report %>%
    reshape::cast(
      experiment + protein_id + precursor_id + replicate + fragment_id ~
        condition,
      value = "log10_fragment_peak_area",
      fun.aggregate = mean,
      na.rm = TRUE
    )

  colnames(df_paired)[c(6, 7)] <- c(
    "Log10NormalizedPeakArea.x", "Log10NormalizedPeakArea.y"
  )

  result_paired0 <- df_paired %>%
    dplyr::group_by(
      .data$experiment,
      .data$protein_id,
      .data$precursor_id
    ) %>%
    dplyr::group_modify(~compute_paired_on_group(.x)) %>%
    as.data.frame()

  return(result_paired0)
}
