## Load required libraries
library(dplyr)
library(corpcor)
library(reshape)
library(boot)


## Load required functions
source("../../R/cov_diff_shrink.R")
source("../../R/cov_diff.R")
source("../../R/independent_t_test.R")
source("../../R/paired_t_test.R")
source("../../R/shrinkage_t_test.R")


## Define wrapper functions
compute_paired_on_group <- function(groupdf) {
  as.data.frame(
    paired_t_test(
      groupdf$Log10NormalizedPeakArea.x,
      groupdf$Log10NormalizedPeakArea.y,
      verbose = FALSE
    )
  )
}

compute_indep_on_group <- function(groupdf) {
  as.data.frame(
    independent_t_test(
      groupdf$Log10Quantity.x,
      groupdf$Log10Quantity.y,
      verbose = FALSE
    )
  )
}

compute_shrink_on_group <- function(groupdf) {
  boot_denom_eps <- 0.5

  ## groupdf consists of two conditions
  df_shrink <- groupdf %>%
    cast(replicate ~ fragment_id ~ condition,
         value = "log10_fragment_peak_area",
         fun.aggregate = mean)

  dat_con1 <- as.matrix(df_shrink[, , 1])
  dat_con2 <- as.matrix(df_shrink[, , 2])

  as.data.frame(shrinkage_t_test(
    dat_con1, dat_con2, num_boot = 100, cov.equal = TRUE,
    boot_denom_eps = boot_denom_eps, verbose = FALSE
  ))
}


#' Compare performance of t-test methods
#'
#' Compute contingency tables for the pairwise comparisons: 1-2, 1-3, ..., 1-C.
#' Each pairwise comparisons are performed on each
#' experiment x protein_id x precursor_id.
#' @param report fragment ion report with the columns experiment, condition,
#' replicate, protein_id, precursor_id, precursor_quantity, fragment_id,
#' fragment_peak_area.
compute_contingency_tables <- function(report, alpha = 0.05) {
  conditions <- unique(report$condition)
  report[["log10_fragment_peak_area"]] <- log10(report[["fragment_peak_area"]])

  result_paired <- result_indep <- result_shrink <- tab_result <- list()

  for (i in 1:(length(conditions) - 1)) {
    report_twoconds <- report %>%
      filter(condition == conditions[1] | condition == conditions[i + 1])

    # Run paired t-test
    df_paired <- report_twoconds %>%
      cast(experiment + protein_id + precursor_id + replicate + fragment_id ~
             condition,
           value = "log10_fragment_peak_area",
           fun.aggregate = mean, na.rm = TRUE)

    colnames(df_paired)[c(6, 7)] <- c("Log10NormalizedPeakArea.x",
                                      "Log10NormalizedPeakArea.y")

    result_paired0 <- df_paired %>%
      group_by(experiment, protein_id, precursor_id) %>%
      group_modify(~compute_paired_on_group(.x)) %>%
      as.data.frame()

    tab_paired0 <- table(result_paired0$p.value < alpha)

    # Run independent t-test
    df_indep <- report_twoconds %>%
      group_by(experiment, protein_id, precursor_id, replicate, condition) %>%
      summarise(log10_peptide_quantity =
                  log10(sum(fragment_peak_area, na.rm = TRUE))) %>%
      cast(experiment + protein_id + precursor_id + replicate ~ condition,
           value = "log10_peptide_quantity", fun.aggregate = mean, na.rm = TRUE)

    colnames(df_indep)[c(5, 6)] <- c("Log10Quantity.x",
                                     "Log10Quantity.y")

    result_indep0 <- df_indep %>%
      group_by(experiment, protein_id, precursor_id) %>%
      group_modify(~compute_indep_on_group(.x)) %>%
      as.data.frame()

    tab_indep0 <- table(result_indep0$p.value < alpha)

    # Run shrinkage t-test
    result_shrink0 <- report_twoconds %>%
      group_by(experiment, protein_id, precursor_id) %>%
      group_modify(~compute_shrink_on_group(.x)) %>%
      as.data.frame()

    tab_shrink0 <- table(result_shrink0$p.value < alpha)

    # Collect analysis results
    result_paired[[i]] <- result_paired0
    result_indep[[i]] <- result_indep0
    result_shrink[[i]] <- result_shrink0

    # Collect contingency tables
    base_tab <- data.frame(Var1 = c(FALSE, TRUE))

    tab_result0 <- base_tab %>%
      merge(tab_paired0, by = "Var1", all = TRUE) %>%
      merge(tab_indep0, by = "Var1", all = TRUE) %>%
      merge(tab_shrink0, by = "Var1", all = TRUE)

    tab_result0[is.na(tab_result0)] <- 0

    colnames(tab_result0) <- c("", "paired", "independent", "shrink")

    tab_result[[i]] <- tab_result0
  }

  return(list(result_paired = result_paired,
              result_indep = result_indep,
              result_shrink = result_shrink,
              table = tab_result))
}
