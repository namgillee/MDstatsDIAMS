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


#' Compare performance of t-test methods
#'
#' Compute contingency tables for the pairwise comparisons: 1-2, 1-3, ..., 1-C.
#' Each pairwise comparisons are performed on each
#' experiment x protein_id x precursor_id.
#' @param report fragment ion report with the columns experiment, condition,
#' replicate, protein_id, precursor_id, precursor_quantity, fragment_id,
#' fragment_peak_area.
compute_contingency_tables <- function(
  report, alpha = 0.05, boot_denom_eps = 0.5
) {
  conditions <- unique(report$condition)
  report[["log10_fragment_peak_area"]] <- log10(report[["fragment_peak_area"]])

  result_paired <- result_indep <- result_shrink <- tab_result <- list()

  for (i in 1:(length(conditions) - 1)) {
    report_twoconds <- report %>%
      filter(.data$condition == conditions[1] |
               .data$condition == conditions[i + 1])

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

    tab_paired0 <- data.frame(table(result_paired0$p.value < alpha))
    names(tab_paired0) <- c("Rejected", "paired")

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

    tab_indep0 <- data.frame(table(result_indep0$p.value < alpha))
    names(tab_indep0) <- c("Rejected", "independent")

    # Run shrinkage t-test
    result_shrink0 <- report_twoconds %>%
      group_by(.data$experiment, .data$protein_id, .data$precursor_id) %>%
      group_modify(
        ~compute_shrink_on_group(
          .x, boot_denom_eps = boot_denom_eps
        )
      ) %>%
      as.data.frame()

    tab_shrink0 <- data.frame(table(result_shrink0$p.value < alpha))
    names(tab_shrink0) <- c("Rejected", "shrinkage")

    # Collect analysis results
    result_paired[[i]] <- result_paired0
    result_indep[[i]] <- result_indep0
    result_shrink[[i]] <- result_shrink0

    # Collect contingency tables
    base_tab <- data.frame(Rejected = c(FALSE, TRUE))

    tab_result0 <- base_tab %>%
      merge(tab_paired0, by = "Rejected", all = TRUE) %>%
      merge(tab_indep0, by = "Rejected", all = TRUE) %>%
      merge(tab_shrink0, by = "Rejected", all = TRUE)

    tab_result0[is.na(tab_result0)] <- 0

    tab_result[[i]] <- tab_result0
  }

  return(list(result_paired = result_paired,
              result_indep = result_indep,
              result_shrink = result_shrink,
              table = tab_result))
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
