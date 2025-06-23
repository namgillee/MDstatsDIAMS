#' Run MSstats comparing two conditions on a standard report
#'
#' @param report fragment ion report with the columns experiment, condition,
#'   replicate, protein_id, precursor_id, fragment_id, and fragment_peak_area.
#'   The condition column is supposed to consist of two conditions.
#' @return a data frame of the test results. Columns are 
#'   experiment, protein_id, precursor_id, statistic, df, p.value, estimate.
#'   And rows are precursors.
#' @export
compute_msstats_on_stdreport <- function(report) {
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
