#' Wrapper function for paired t-test
#'
#' @param groupdf A subset of fragment ion report which contains two columns for
#' log10-transformed fragment ion quantities of conditions x and y
#' @param column_x column name for x
#' @param column y column name for y
#' @return a data frame of the paired t-test results
#' @export
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


#' Wrapper function for independent samples t-test
#'
#' @param groupdf A subset of precursor report which contains two columns for
#' log10-transformed precursor quantities of conditions x and y
#' @param column_x column name for x
#' @param column y column name for y
#' @return a data frame of the independent samples t-test results
#' @export
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


#' Wrapper function for shrinkage t-test
#'
#' @param groupdf A subset of precursor report which consists of two conditions.
#'   It contains columns named "replicate", "fragment_id", "condition", and
#'   a value column.
#' @param value_column column name for log10-transformed fragment ion quantities
#' @param cov_unequal_replicates_column column name for covariance of input
#'   variables whose replicates are not equal. All entries in the column is
#'   considered equal. If NULL, the covariance is set to zero.
#' @param boot_denom_eps A parameter for shrinkage t-test
#' @return a data frame of the shrinkage t-test results
#' @export
compute_shrink_on_group <- function(
  groupdf,
  value_column = "log10_fragment_peak_area",
  cov_unequal_replicates_column = NULL,
  boot_denom_eps = 0.5
) {
  if (length(unique(groupdf$condition)) < 2) {
    return(as.data.frame(list()))
  }

  cov_unequal_replicates <- 0
  if (!is.null(cov_unequal_replicates_column) &&
        !is.null(groupdf[[cov_unequal_replicates_column]])) {
    cov_unequal_replicates <- groupdf[[cov_unequal_replicates_column]][
      !is.na(groupdf[[cov_unequal_replicates_column]])
    ][1]
  }

  df_shrink <- groupdf %>%
    reshape::cast(
      replicate ~ fragment_id ~ condition,
      value = value_column,
      fun.aggregate = mean
    )

  dat_con1 <- as.matrix(df_shrink[, , 1])
  dat_con2 <- as.matrix(df_shrink[, , 2])

  as.data.frame(shrinkage_t_test(
    dat_con1, dat_con2, cov_unequal_replicates = cov_unequal_replicates,
    num_boot = 100, cov_equal = TRUE, boot_denom_eps = boot_denom_eps,
    verbose = FALSE
  ))
}


#' Compute cov_unequal_replicates
#'
#' Compute cov_unequal_replicates from data values
#' @param x  a vector of mean precursor quantities in log10-scale
#' @return  the cov_unequal_replicates value
.comp_cov_uneq_repl <- function(x) {

  # compute the sample variance
  var_x <- var(x)

  # estimate a mixture of normal and exp-beta distributions
  fit_kde <- stats::density(x)

  idx_max_kde <- which.max(fit_kde$y)
  x_mode <- fit_kde$x[idx_max_kde]

  sigma_left <- diff(quantile(x_mode - x[x < x_mode], (c(4, 9)) / 10)) /
    diff(qnorm(0.5 + (c(4, 9)) / 20))
  sigma_right <- diff(quantile(x[x > x_mode] - x_mode, (c(4, 9)) / 10)) /
    diff(qnorm(0.5 + (c(4, 9)) / 20))

  init_sigma <- min(sigma_left, sigma_right)

  cov_unequal_replicates <- max(0, var_x - init_sigma^2)

  return(cov_unequal_replicates)
}


#' Compute and attach cov_unequal_replicates
#'
#' Compute cov_unequal_replicates by each protein, and attach it into a new
#' column. The fragment_ion_peak area is first summarized to
#' log10_peptide_quantity in precursor level. It is then used for estimating
#' a mixture of normal and expbeta distributions by each protein.
#' @param report_df  a fragment ion report with multiple conditions and the
#'   fragment_peak_area column
#' @return  the report data frame with a new column "cov_unequal_replicates"
#' @export
compute_cov_unequal_replicates <- function(report_df) {
  # compute log10 peptide quantity
  precursor_df <- report_df %>%
    dplyr::group_by(
      .data$experiment,
      .data$protein_id,
      .data$precursor_id,
      .data$replicate,
      .data$condition
    ) %>%
    dplyr::summarise(
      log10_peptide_quantity = log10(
        sum(.data$fragment_peak_area, na.rm = TRUE)
      )
    )

  # median-normalize them by each condition
  #   A negative value is not removed but ignored.
  median_log10_quantity <- median(
    precursor_df$log10_peptide_quantity, na.rm = TRUE
  )

  precursor_df_centered <- precursor_df %>%
    dplyr::group_by(
      .data$condition
    ) %>%
    dplyr::mutate(
      median_log10_quantity_by_condition =
        median(.data$log10_peptide_quantity, na.rm = TRUE)
    ) %>%
    dplyr::mutate(
      log10_peptide_quantity =
        .data$log10_peptide_quantity -
        .data$median_log10_quantity_by_condition +
        median_log10_quantity
    )

  # compute the mean of log10 peptide quantity
  mean_precursor_df_centered <- precursor_df_centered %>%
    dplyr::group_by(.data$experiment,
                    .data$protein_id,
                    .data$precursor_id,
                    .data$condition) %>%
    dplyr::summarise(
      mean_log10_peptide_quantity = mean(
        .data$log10_peptide_quantity, na.rm = TRUE
      )
    )

  # compute the cov_unequal_replicates by each protein
  cov_unequal_replicates <- .comp_cov_uneq_repl(
    mean_precursor_df_centered$mean_log10_peptide_quantity
  )

  report_df$cov_unequal_replicates <- cov_unequal_replicates

  return(report_df)
}


#' Run t-test methods
#'
#' Run t-test methods for the pairwise comparisons: 1-2, 1-3, ..., 1-C.
#' Each pairwise comparisons are performed on each
#' experiment x protein_id x precursor_id.
#' @param report fragment ion report with the columns experiment, condition,
#'   replicate, protein_id, precursor_id, fragment_id, fragment_peak_area.
#' @param boot_denom_eps a parameter for shrinkage t-test
#' @return a list of analysis results of three t-test methods, where the result
#'   of each method is a list for the pairwise comparisons.
#' @export
run_ttests <- function(report, boot_denom_eps = 0.5, base_condition = NULL) {
  conditions <- unique(report$condition)
  n_con <- length(conditions)
  report[["log10_fragment_peak_area"]] <- log10(report[["fragment_peak_area"]])

  if (is.null(base_condition)) {
    base_condition <- conditions[1]
  } else if (!(base_condition %in% conditions)) {
    print(paste(base_condition, "is not in the conditions. Using default."))
    base_condition <- conditions[1]
  }

  # Compute cov_unequal_replicates for shrinkage-t-test
  report <- compute_cov_unequal_replicates(report)

  # Run t-test methods with two conditions
  con1 <- base_condition
  con_rest <- conditions[conditions != con1]
  result_paired <- result_indep <- result_shrink <- vector("list", n_con - 1)
  names(result_paired) <- names(result_indep) <- names(result_shrink) <-
    paste0(con1, "/", con_rest)

  for (con2 in con_rest) {
    id <- paste0(con1, "/", con2)

    report_twoconds <- report %>%
      filter(.data$condition == con1 | .data$condition == con2)


    # Run paired t-test
    df_paired <- report_twoconds %>%
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


    # Run independent t-test
    df_indep <- report_twoconds %>%
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


    # Run shrinkage t-test
    result_shrink0 <- report_twoconds %>%
      dplyr::group_by(
        .data$experiment,
        .data$protein_id,
        .data$precursor_id
      ) %>%
      dplyr::group_modify(
        ~compute_shrink_on_group(
          .x,
          value_column = "log10_fragment_peak_area",
          cov_unequal_replicates_column = "cov_unequal_replicates",
          boot_denom_eps = boot_denom_eps
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
#'   analysis results of each method.
#' @param alpha significance level
#' @return list of contingency tables for every comparisons. Each element of
#'   the list is a table, whose columns are the t-test methods and the rows are
#'   TRUE/FALSE representing whether H0 is rejected or not.
#' @export
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


#' Line plot for contingency tables
#'
#' @param x horizontal coordinates of the line plots
#' @param tables list contingency tables. Lengths of x and tables must be equal.
#' @param rejected TRUE if number of rejected hypothesis is plotted in y axis.
#' @param scale_factor a numeric factor that is multiplied to data values.
#' @param add_legend If not FALSE, legend is added in the plot
#' @param legend_coord coordinate of the legend
#' @param legend_cex cex for the legend
#' @examples
#'   report <- simulate_fragment_ion_report(default_params)
#'   resu <- run_ttests(report, boot_denom_eps = 0.5)
#'   tables <- compute_contingency_tables(resu, alpha = 0.05)
#'   x <- default_params$prec_mean_condition_shift[-c(1, 2)]
#'   line_plot_contingency_tables(
#'     x, tables[-1], xlab = expression(delta),
#'     ylab = "1 - Type II error rate", cex.lab = 1.5
#'   )a
#' @export
line_plot_contingency_tables <- function(
  x, tables, rejected = TRUE, scale_factor = 1, add_legend = FALSE,
  legend_coord = "topright", legend_cex = 1, ...
) {
  name_tables <- names(tables)
  len_tables <- length(tables)
  name_methods <- colnames(tables[[1]])[-1]
  n_methods <- length(name_methods)

  values <- matrix(
    NA, len_tables, n_methods, dimnames = list(name_tables, name_methods)
  )

  for (name in name_tables) {
    a_table <- tables[[name]]
    for (method in name_methods) {
      values[name, method] <- a_table[a_table[[1]] == rejected, method]
    }
  }

  values <- values * scale_factor

  seq_hcl_colors <- rep(2 : 4, ceiling(n_methods / 3))
  line_types <- rep(1 : 2, each = ceiling(n_methods / 2))
  pch_types <- c(1, 2, 6, 3, 4, 16, 7 : 10)[1 : n_methods]
  matplot(x, values, type = "o", lty = line_types,
          col = seq_hcl_colors, pch = pch_types, ...)

  if (add_legend != FALSE) {
    legend(legend_coord, legend = name_methods, lty = line_types,
           col = seq_hcl_colors, pch = pch_types, cex = legend_cex,
           lwd = 2)
  }
}


#' bar plot for contingency tables
#'
#' @param tables list contingency tables.
#' @param rejected TRUE if number of rejected hypothesis is plotted in y axis.
#' @param scale_factor a numeric factor that is multiplied to data values.
#' @param add_legend If not FALSE, legend is added in the plot
#' @param legend_ncol ncol for the legend
#' @param legend_cex cex for the legend
#' @examples
#'   report <- simulate_fragment_ion_report(default_params)
#'   resu <- run_ttests(report, boot_denom_eps = 0.5)
#'   tables <- compute_contingency_tables(resu, alpha = 0.05)
#'   bar_plot_contingency_tables(
#'     tables[1], xlab = "Comparison", ylab = "1 - Type I error rate",
#'     cex.lab = 1.5
#'   )
#' @export
bar_plot_contingency_tables <- function(
  tables, rejected = FALSE, scale_factor = 1, add_legend = FALSE,
  legend_ncol = 2, legend_cex = 1, ...
) {
  name_tables <- names(tables)
  len_tables <- length(tables)
  name_methods <- colnames(tables[[1]])[-1]
  n_methods <- length(name_methods)

  values <- matrix(
    NA, len_tables, n_methods,
    dimnames = list(name_tables, name_methods)
  )

  for (name in name_tables) {
    a_table <- tables[[name]]
    for (method in name_methods) {
      values[name, method] <- a_table[a_table[[1]] == rejected, method]
    }
  }

  values <- values * scale_factor

  seq_hcl_colors <- colorspace::sequential_hcl(n_methods, "hawaii")
  barplot(t(values), beside = TRUE, col = seq_hcl_colors, ...)

  if (add_legend != FALSE) {
    legend_text <- name_methods
    legend(
      "top",
      fill = seq_hcl_colors,
      legend = legend_text,
      inset = c(0, -0.15),
      xpd = TRUE,
      ncol = legend_ncol,
      cex = legend_cex
    )
  }

}
