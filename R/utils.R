library(dplyr)
library(corpcor)
library(reshape)
library(boot)

source("../../R/independent_t_test.R")
source("../../R/paired_t_test.R")
source("../../R/shrinkage_t_test.R")


#' Wrapper function for paired t-test
#'
#' @param groupdf A subset of fragment ion report which contains two columns for
#' log10-transformed fragment ion quantities of conditions 1 and 2
compute_paired_on_group <- function(groupdf) {
  as.data.frame(
    paired_t_test(
      groupdf$Log10NormalizedPeakArea.x,
      groupdf$Log10NormalizedPeakArea.y,
      verbose = FALSE
    )
  )
}


#' Wrapper function for independent samples t-test
#'
#' @param groupdf A subset of precursor report which contains two columns for
#' log10-transformed precursor quantities of conditions 1 and 2
compute_indep_on_group <- function(groupdf) {
  as.data.frame(
    independent_t_test(
      groupdf$Log10Quantity.x,
      groupdf$Log10Quantity.y,
      verbose = FALSE
    )
  )
}


#' Wrapper function for shrinkage t-test
#'
#' @param groupdf A subset of precursor report which consists of two conditions
#' @param boot_denom_eps A parameter for shrinkage t-test
compute_shrink_on_group <- function(groupdf, boot_denom_eps = 0.5) {
  df_shrink <- groupdf %>%
    cast(replicate ~ fragment_id ~ condition,
         value = "log10_fragment_peak_area",
         fun.aggregate = mean)

  dat_con1 <- as.matrix(df_shrink[, , 1])
  dat_con2 <- as.matrix(df_shrink[, , 2])

  as.data.frame(shrinkage_t_test(
    dat_con1, dat_con2, num_boot = 100, cov_equal = TRUE,
    boot_denom_eps = boot_denom_eps, verbose = FALSE
  ))
}


#' Run t-test methods
#'
#' Run t-test methods for the pairwise comparisons: 1-2, 1-3, ..., 1-C.
#' Each pairwise comparisons are performed on each
#' experiment x protein_id x precursor_id.
#' @param report fragment ion report with the columns experiment, condition,
#' replicate, protein_id, precursor_id, precursor_quantity, fragment_id,
#' fragment_peak_area.
#' @param boot_denom_eps a parameter for shrinkage t-test
#' @return list of analysis results of three t-test methods, where the result of
#' each method is a list for the pairwise comparisons.
run_ttests <- function(report, boot_denom_eps = 0.5) {
  conditions <- unique(report$condition)
  n_con <- length(conditions)
  report[["log10_fragment_peak_area"]] <- log10(report[["fragment_peak_area"]])

  result_paired <- result_indep <- result_shrink <- vector("list", n_con - 1)
  names(result_paired) <- names(result_indep) <- names(result_shrink) <-
    paste0(conditions[1], "/", conditions[2 : n_con])

  con1 <- conditions[1]
  for (con2 in conditions[2 : n_con]) {
    id <- paste0(con1, "/", con2)

    report_twoconds <- report %>%
      filter(.data$condition == con1 | .data$condition == con2)

    # Run paired t-test
    df_paired <- report_twoconds %>%
      cast(experiment + protein_id + precursor_id + replicate + fragment_id ~
             condition,
           value = "log10_fragment_peak_area",
           fun.aggregate = mean, na.rm = TRUE)

    colnames(df_paired)[c(6, 7)] <- c("Log10NormalizedPeakArea.x",
                                      "Log10NormalizedPeakArea.y")

    result_paired0 <- df_paired %>%
      group_by(.data$experiment, .data$protein_id, .data$precursor_id) %>%
      group_modify(~compute_paired_on_group(.x)) %>%
      as.data.frame()

    # Run independent t-test
    df_indep <- report_twoconds %>%
      group_by(.data$experiment, .data$protein_id, .data$precursor_id,
               .data$replicate, .data$condition) %>%
      summarise(log10_peptide_quantity =
                  log10(sum(.data$fragment_peak_area, na.rm = TRUE))) %>%
      cast(experiment + protein_id + precursor_id + replicate ~ condition,
           value = "log10_peptide_quantity", fun.aggregate = mean, na.rm = TRUE)

    colnames(df_indep)[c(5, 6)] <- c("Log10Quantity.x",
                                     "Log10Quantity.y")

    result_indep0 <- df_indep %>%
      group_by(.data$experiment, .data$protein_id, .data$precursor_id) %>%
      group_modify(~compute_indep_on_group(.x)) %>%
      as.data.frame()

    # Run shrinkage t-test
    result_shrink0 <- report_twoconds %>%
      group_by(.data$experiment, .data$protein_id, .data$precursor_id) %>%
      group_modify(
        ~compute_shrink_on_group(
          .x, boot_denom_eps = boot_denom_eps
        )
      ) %>%
      as.data.frame()

    # Collect analysis results
    result_paired[[id]] <- result_paired0
    result_indep[[id]] <- result_indep0
    result_shrink[[id]] <- result_shrink0
  }

  return(list(paired = result_paired,
              independent = result_indep,
              shrinkage = result_shrink))
}


#' Compute contingency tables
#'
#' Compute contingency tables based on the results of running t-test methods.
#' @param results_run_ttests results of run_ttest function, which is a list of
#' analysis results of each method.
#' @param alpha significance level
compute_contingency_tables <- function(results_run_ttests, alpha = 0.05) {
  id_methods <- names(results_run_ttests)
  id_comparisons <- names(results_run_ttests[[1]])

  tab_result <- vector("list", length(id_comparisons))
  names(tab_result) <- id_comparisons

  for (id in id_comparisons) {
    tab_result_comparison <- data.frame(Rejected = c(FALSE, TRUE))

    for (method in id_methods) {
      result_method <- results_run_ttests[[method]]

      # Compute contingency tables of every methods
      tab_method_comparison <- data.frame(
        table(
          result_method[[id]]$p.value < alpha
        )
      )
      names(tab_method_comparison) <- c("Rejected", method)

      # Collect contingency tables
      tab_result_comparison <- tab_result_comparison %>%
        merge(tab_method_comparison, by = "Rejected", all = TRUE)
    }

    tab_result_comparison[is.na(tab_result_comparison)] <- 0

    tab_result[[id]] <- tab_result_comparison
  }

  return(tab_result)
}


#' Generate random samples from a mixture of Beta distributions
#'
#' @param n  number of observations
#' @param shape1s,shape2s  vectors of non-negative parameters of component
#' beta distributions. Length of each vector is the number of component beta
#' distributions.
rbetamixture <- function(n, shape1s = 1, shape2s = 1) {
  # Number of beta components
  n_beta1 <- length(shape1s)
  n_beta2 <- length(shape2s)
  if ((n_beta1 == 1) && (n_beta2 > 1)) {
    shape1s <- rep(shape1s, n_beta2)
    n_beta1 <- n_beta2
  }
  if ((n_beta1 > 1) && (n_beta2 == 1)) {
    shape2s <- rep(shape2s, n_beta1)
    n_beta2 <- n_beta1
  }
  if (n_beta1 != n_beta2) {
    stop("Lengths of shape1s and shape2s must equal.")
  }

  # Select beta components
  id_components <- sample(n_beta1, size = n, replace = TRUE)

  # Sample from the selected beta components
  out <- rep(-1, n)
  for (i in 1 : n) {
    id <- id_components[i]
    out[i] <- rbeta(1, shape1s[id], shape2s[id])
  }

  return(out)
}
