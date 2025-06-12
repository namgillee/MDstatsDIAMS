#' Convert a spectronaut fragment ion report to a standard report
#'
#' @param sn_report A spectronaut fragment ion report. Required columns are
#'   R.Condition, R.Replicate, PG.ProteinGroups, EG.ModifiedSequence,
#'   FG.Charge, EG.Qvalue, F.NormalizedPeakArea,
#'   F.FrgIon, F.FrgLossType, F.Charge, F.ExcludedFromQuantification.
#' @param filter_identified TRUE if only identified precursors remain in the
#'   report.
#' @importFrom dplyr %>%
#' @return A standard report with columns condition, replicate, experiment,
#'   protein_id, precursor_id, precursor_qvalue, fragment_id, fragment_peak_area
#' @export
convert_sn_to_standard <- function(sn_report, filter_identified = TRUE) {

  # Filter
  if (filter_identified) {
    sn_report <- sn_report[
      (sn_report$F.ExcludedFromQuantification == "False") &
        (sn_report$EG.Qvalue <= 0.01) &
        (sn_report$F.NormalizedPeakArea > 1), ]
  } else {
    sn_report <- sn_report[
      (sn_report$F.ExcludedFromQuantification == "False") &
        (sn_report$F.NormalizedPeakArea > 1), ]
  }

  # Select columns
  sn_report <- sn_report %>%
    dplyr::select(
      R.Condition, R.Replicate, PG.ProteinGroups,
      EG.ModifiedSequence, EG.Qvalue, FG.Charge,
      F.FrgIon, F.FrgLossType, F.Charge, F.NormalizedPeakArea
    )

  # Rename columns
  condition_labels <- unique(sn_report$R.Condition)
  if ("DMSO" %in% condition_labels) {
    condition_labels <- c("DMSO", condition_labels[condition_labels != "DMSO"])
  }

  sn_report <- sn_report %>%
    dplyr::mutate(
      condition = factor(.data$R.Condition, labels = condition_labels),
      experiment = "EXP",
      precursor_id = paste0(.data$EG.ModifiedSequence, ".", .data$FG.Charge),
      fragment_id = paste0(
        .data$F.FrgIon, .data$F.FrgLossType, ".", .data$F.Charge
      )
    ) %>%
    dplyr::rename(
      replicate = .data$R.Replicate,
      protein_id = .data$PG.ProteinGroups,
      precursor_qvalue = .data$EG.Qvalue,
      fragment_peak_area = .data$F.NormalizedPeakArea
    ) %>%
    dplyr::select(
      -c(R.Condition, EG.ModifiedSequence, FG.Charge,
         F.FrgIon, F.FrgLossType, F.Charge)
    )

  return(sn_report)
}
