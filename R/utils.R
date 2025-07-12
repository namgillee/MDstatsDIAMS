#' Compute cov_unequal_replicates
#'
#' Compute cov_unequal_replicates from data values
#' @param x  a vector of mean precursor quantities in log10-scale
#' @return  the cov_unequal_replicates value
.comp_cov_uneq_repl <- function(x) {

  # compute the sample variance
  var_x <- stats::var(x)

  # estimate a mixture of normal and exp-beta distributions
  fit_kde <- stats::density(x)

  idx_max_kde <- which.max(fit_kde$y)
  x_mode <- fit_kde$x[idx_max_kde]

  sigma_left <- diff(stats::quantile(x_mode - x[x < x_mode], (c(4, 9)) / 10)) /
    diff(stats::qnorm(0.5 + (c(4, 9)) / 20))
  sigma_right <- diff(stats::quantile(x[x > x_mode] - x_mode, (c(4, 9)) / 10)) /
    diff(stats::qnorm(0.5 + (c(4, 9)) / 20))

  init_sigma <- min(sigma_left, sigma_right)

  cov_unequal_replicates <- max(0, var_x - init_sigma^2)
}


#' Compute cov_unequal_replicates
#'
#' Compute cov_unequal_replicates by each protein, and append it into a new
#' column. The fragment_ion_peak area is first summarized to
#' log10_peptide_quantity in precursor level. It is then used for estimating
#' a mixture of normal and expbeta distributions by each protein.
#' @param report_df  a fragment ion report with multiple conditions and the
#'   fragment_peak_area column
#' @importFrom dplyr %>%
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
  median_log10_quantity <- stats::median(
    precursor_df$log10_peptide_quantity, na.rm = TRUE
  )

  precursor_df_centered <- precursor_df %>%
    dplyr::group_by(
      .data$condition
    ) %>%
    dplyr::mutate(
      median_log10_quantity_by_condition =
        stats::median(.data$log10_peptide_quantity, na.rm = TRUE)
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
#' Run t-test methods for the pairwise comparisons:
#' cond1 / cond2, cond1 / cond3, ..., cond1 / cond_C,
#' where cond1 is called the base condition.
#' Each pairwise comparisons are performed on each
#' experiment x protein_id x precursor_id.
#' @param report Fragment ion report with the columns experiment, condition,
#'   replicate, protein_id, precursor_id, fragment_id, fragment_peak_area.
#' @param method_names A vector of characters for the method names to execute.
#'   If NULL, it is set to all the method names available. Defaults to NULL.
#' @param boot_denom_eps A parameter for shrinkage t-test.
#' @param base_condition Base condition. If NULL, the first condition among the
#'   unique list of the condition column is selected. Default to NULL.
#' @param verbose TRUE to print messages.
#' @importFrom dplyr %>%
#' @return A list of analysis results of the given methods. Each analysis result
#'   is a list of pairwise comparison results.
#' @export
run_ttests <- function(
  report, method_names = NULL, boot_denom_eps = 0.3, base_condition = NULL,
  verbose = TRUE
) {
  # List of available methods
  method_name_func <- list(
    msstatslip = compute_mslip_on_stdreport,
    rots = compute_rots_on_stdreport,
    paired = compute_paired_on_stdreport,
    independent = compute_indep_on_stdreport,
    shrinkage = compute_shrink_on_stdreport
  )

  # Set method names
  if (is.null(method_names)) {
    method_names <- names(method_name_func)
  }
  if (any(id_name_not_exist <- !(method_names %in% names(method_name_func)))) {
    warning(
      paste(
        "Method names,",
        method_names[id_name_not_exist],
        ", do not exist.",
        collapse = " "
      )
    )
    method_names <- method_names[!id_name_not_exist]
  }

  if (verbose)
    message(paste(c("Running the test methods:", method_names), collapse = " "))

  # Set base_condition
  conditions <- levels(factor(report$condition))

  if (is.null(base_condition)) {
    base_condition <- conditions[1]
  } else if (!(base_condition %in% conditions)) {
    warning(paste(base_condition, "is not in the conditions. Using default."))
    base_condition <- conditions[1]
  }

  # Compute cov_unequal_replicates for shrinkage t-test, creating a new column
  report <- compute_cov_unequal_replicates(report)

  # Run t-test methods with two conditions
  cond_rest <- conditions[conditions != base_condition]
  n_comparison <- length(cond_rest)

  ## 1. Iterate across method_names
  test_results <- vector("list", length(method_names))
  names(test_results) <- method_names
  for (method_name in method_names) {
    if (verbose)
      message(paste(method_name, "..."))

    method_func <- method_name_func[[method_name]]

    ## 2. Iterate across comparisons
    method_results <- vector("list", n_comparison)
    names(method_results) <- paste0(base_condition, "/", cond_rest)
    for (cond in cond_rest) {
      comparison <- paste0(base_condition, "/", cond)

      report_twoconds <- report %>% dplyr::filter(
        .data$condition == base_condition | .data$condition == cond
      )

      result0 <- NULL
      if (method_name != "shrinkage") {
        result0 <- method_func(report_twoconds)
      } else {
        result0 <- method_func(
          report_twoconds, boot_denom_eps = boot_denom_eps
        )
      }

      method_results[[comparison]] <- result0
    }
    test_results[[method_name]] <- method_results
  }
  return(test_results)
}


#' Compute contingency tables
#'
#' Compute contingency tables based on the results of running t-test methods.
#' @param test_results Results of run_ttest function, which is a named
#'   list of analysis results of each method.
#' @param alpha Significance level
#' @param q_value_column Name of the column having significance levels: A
#'   smaller value implies more significance.
#' @return List of contingency tables for every comparisons. Each element of
#'   the list is a table, whose columns are the t-test methods and the rows are
#'   FALSE/TRUE, where FALSE represents that null hypothesis is not rejected,
#'.  and TRUE represents rejected.
#' @export
compute_contingency_tables <- function(
  test_results, alpha = 0.05, q_value_column = "p.value"
) {
  method_names <- names(test_results)
  comparisons <- names(test_results[[1]])

  contingency_tables <- vector("list", length(comparisons))
  names(contingency_tables) <- comparisons

  for (comparison in comparisons) {
    tab_result_comparison <- data.frame(Rejected = c(FALSE, TRUE))

    for (method_name in method_names) {
      result_method <- test_results[[method_name]]

      # Compute contingency tables of every methods
      tab_method_comparison <- data.frame(
        table(result_method[[comparison]][[q_value_column]] < alpha)
      )
      names(tab_method_comparison) <- c("Rejected", method_name)

      # Collect contingency tables
      tab_result_comparison <- merge(
        tab_result_comparison,
        tab_method_comparison,
        by = "Rejected",
        all = TRUE
      )
    }

    tab_result_comparison[is.na(tab_result_comparison)] <- 0

    contingency_tables[[comparison]] <- tab_result_comparison
  }

  return(contingency_tables)
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
#' @param ... Additional arguments to graphics::matplot
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
  graphics::matplot(x, values, type = "o", lty = line_types,
                    col = seq_hcl_colors, pch = pch_types, ...)

  if (add_legend != FALSE) {
    graphics::legend(legend_coord, legend = name_methods, lty = line_types,
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
#' @param ... Additional arguments to barplot()
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
  graphics::barplot(t(values), beside = TRUE, col = seq_hcl_colors, ...)

  if (add_legend != FALSE) {
    legend_text <- name_methods
    graphics::legend(
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
