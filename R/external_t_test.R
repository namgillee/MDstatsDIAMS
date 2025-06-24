#' Run MSstatsLiP comparing two conditions on a standard report
#'
#' @param report Fragment ion report with the columns experiment, condition,
#'   replicate, protein_id, precursor_id, fragment_id, and fragment_peak_area.
#'   The condition column is supposed to consist of two conditions.
#' @return A data frame of the test results. Columns are
#'   experiment, protein_id, peptide_sequence, statistic, df, p.value, estimate,
#'   Label, and miscellaneous columns from MSstatsLiP::groupComparisonLiP().
#'   And rows are precursors. Note that precursor_id column is missing.
#' @importFrom dplyr %>%
#' @export
compute_mslip_on_stdreport <- function(report) {
  # Convert to MSstatsLiP format
  mslip_report <- convert_standard_to_mslip(report, drop_experiment = FALSE)
  experiment_levels <- unique(report$experiment)
  test_result <- NULL

  for (exp_level in experiment_levels) {
    # MSstatsLiP dataSummarization
    # Run within suppressMessages to supress unnecessary messages arising from
    # MSstatsPTM::dataSummarizationPTM().
    suppressMessages(
      summarized <- MSstatsLiP::dataSummarizationLiP(
        list(LiP = mslip_report[mslip_report$experiment == exp_level, ]),
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
    )

    # MSstatsLiP groupComparison
    # Run within suppressMessages to supress unnecessary messages arising from
    # MSstatsPTM::dataSummarizationPTM().
    suppressMessages(
      model <- MSstatsLiP::groupComparisonLiP(
        data = summarized,
        contrast.matrix = "pairwise",
        use_log_file = FALSE,
        verbose = FALSE
      )
    )

    # Convert columns to a standard ttest result format.
    # Note: The peptide_sequence column is not standard and may be redundant.
    # The Label, SE, issue columns will remain the same.
    lip_model <- model[["LiP.Model"]]
    lip_model <- lip_model %>%
      dplyr::rename(
        protein_id = ProteinName,
        peptide_sequence = PeptideSequence,
        statistic = Tvalue,
        df = DF,
        p.value = pvalue,
        adj.p.value = adj.pvalue,
        estimate = log2FC
      ) %>%
      select(-c(FULL_PEPTIDE))

    lip_model$experiment <- exp_level

    test_result <- rbind(test_result, lip_model)
  }

  return(test_result)
}


#' Run ROTS comparing two conditions on a standard report
#'
#' @param report fragment ion report with the columns experiment, condition,
#'   replicate, protein_id, precursor_id, fragment_id, and fragment_peak_area.
#'   The condition column is supposed to consist of two conditions.
#' @return a data frame of the test results. Columns are
#'   experiment, protein_id, precursor_id, statistic, df, p.value, estimate,
#'   FDR, Label. And rows are precursors.
#' @export
compute_rots_on_stdreport <- function(report) {
  condition_levels <- levels(as.factor(report$condition))
  condition_levels <- condition_levels[condition_levels %in% report$condition]
  n_conditions <- length(condition_levels)

  test_result <- NULL

  # Summarize fragment peak area for log10 precursor quantity
  if (!("log10_precursor_quantity" %in% names(report))) {
    if (!("precursor_quantity" %in% names(report))) {
      report <- report %>%
        dplyr::group_by(
          .data$experiment,
          .data$protein_id,
          .data$precursor_id,
          .data$replicate,
          .data$condition
        ) %>%
        dplyr::mutate(
          log10_precursor_quantity = log10(
            sum(.data$fragment_peak_area, na.rm = TRUE)
          )
        )
    } else {
      report <- report %>%
        mutate(log10_precursor_quantity = log10(.data$precursor_quantity))
    }
  }

  # Precursor report
  report <- report %>% dplyr::distinct(
    .data$experiment,
    .data$protein_id,
    .data$precursor_id,
    .data$replicate,
    .data$condition,
    .keep_all = TRUE
  )

  for (con_id1 in 1:(n_conditions - 1)) {
    con1 <- condition_levels[con_id1]
    for (con_id2 in (con_id1 + 1):n_conditions) {
      con2 <- condition_levels[con_id2]

      # Convert to ROTS format
      rots_report <- report %>%
        dplyr::filter(
          (.data$condition == con1 | .data$condition == con2)
        ) %>%
        reshape::cast(
          experiment + protein_id + precursor_id ~ condition + replicate,
          value = "log10_precursor_quantity",
          fun.aggregate = mean,
          na.rm = TRUE
        )

      sample_id1 <- substr(colnames(rots_report)[-c(1 : 3)], 1, nchar(con1)) ==
        con1
      contrast <- as.numeric(sample_id1)

      # Filter rots input report
      rots_filtered <- rots_report[, -c(1 : 3)]
      id_filtered1 <- rowSums(!is.na(rots_filtered[, contrast == 1])) >= 2
      id_filtered0 <- rowSums(!is.na(rots_filtered[, contrast == 0])) >= 2
      id_filtered <- id_filtered1 & id_filtered0

      rots_model_list <- ROTS::ROTS(
        data = rots_report[id_filtered, -c(1 : 3)],
        groups = contrast,
        B = 100,
        K = 500,
        seed = 1234,
        verbose = FALSE
      )
      rots_model_df <- as.data.frame(
        rots_model_list[c("d", "pvalue", "logfc", "FDR")]
      )
      rots_model_df <- rots_model_df %>%
        dplyr::mutate(
          experiment = rots_report[id_filtered, 1],
          protein_id = rots_report[id_filtered, 2],
          precursor_id = rots_report[id_filtered, 3],
          Label = paste0(con1, " vs ", con2)
        ) %>%
        dplyr::rename(
          statistic = d,
          p.value = pvalue,
          estimate = logfc
        )

      test_result <- rbind(test_result, rots_model_df)
    }
  }

  return(test_result)
}
