#' Run MSstatsLiP comparing two conditions on a standard report
#'
#' @param report Fragment ion report with the columns experiment, condition,
#'   replicate, protein_id, precursor_id, fragment_id, and fragment_peak_area.
#'   The condition column is supposed to consist of two conditions.
#' @return A data frame of the test results. Columns are 
#'   experiment, protein_id, precursor_id, statistic, df, p.value, estimate.
#'   And rows are precursors.
#' @importFrom dplyr %>%
#' @export
compute_msstatslip_on_stdreport <- function(report) {
  # Convert column names
  msstats_report <- convert_standard_to_ms(report)
  
  # Run MSstats processData
  summarized < MSstats::dataProcess(
    msstats_report,
    logTrans = 2,
    normalization = "equalizeMedians",
    featureSubset = "all",
    n_top_feature = 3,
    summaryMethod = "TMP",
    equalFeatureVar = TRUE,
    censoredInt = "NA",
    MBimpute = TRUE,
    use_log_file = FALSE,
    verbose = FALSE
  )
  
  # Run MSstats analysis
  model <- MSstats::groupComparison(
    "pairwise",
    summarized,
    use_log_file = FALSE,
    verbose = FALSE
  )
  
  # Convert column 
}


#' Run msqrob2 comparing two conditions on a standard report
#'
#' @param report fragment ion report with the columns experiment, condition,
#'   replicate, protein_id, precursor_id, fragment_id, and fragment_peak_area.
#'   The condition column is supposed to consist of two conditions.
#' @return a data frame of the test results. Columns are 
#'   experiment, protein_id, precursor_id, statistic, df, p.value, estimate.
#'   And rows are precursors.
#' @export
compute_msqrob2_on_stdreport <- function(report) {
}


#' Run ROTS comparing two conditions on a standard report
#'
#' @param report fragment ion report with the columns experiment, condition,
#'   replicate, protein_id, precursor_id, fragment_id, and fragment_peak_area.
#'   The condition column is supposed to consist of two conditions.
#' @return a data frame of the test results. Columns are 
#'   experiment, protein_id, precursor_id, statistic, df, p.value, estimate.
#'   And rows are precursors.
#' @export
compute_rots_on_stdreport <- function(report) {
}
