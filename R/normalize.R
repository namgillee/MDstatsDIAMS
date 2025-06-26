#' Normalize values to standard normal distribution by sliding window
#'
#' @param values_df data frame with value_column and sliding_column
#' @param value_column name of the value column
#' @param sliding_column name of the score column that sliding window runs over
#' @param window_size window size for sliding column with which an empirical
#' cdf is estimated
#' @param step_size step size for sliding column for which corresponding
#' value_column is normalized
#' @param reconstruct If TRUE, the normalized_value_column takes the robust mean
#' and robust standard deviation of the original value_column. If FALSE, the
#' mean and standard deviation are set to 0 and 1.
#' @param use_logvalues If TRUE, the value column is log10-transformed,
#' normalized, and then transformed back to the original exponential scale.
#' @param normalized_value_column name of the normalized value column
#' @param normalization_factor_column name of the normalization factor column
#' @return data frame with normalized values set in the normalized_value_column
#' @export
sw_normalize_values_df <- function(
  values_df, value_column = "precursor_quantity", sliding_column = "score",
  window_size = 100000, step_size = 50000, reconstruct = TRUE,
  use_logvalues = TRUE,
  normalized_value_column = "normalized_precursor_quantity",
  normalization_factor_column = "normalization_factor"
) {
  data_size <- nrow(values_df)

  if (window_size > data_size) {
    window_size <- data_size
  }

  # Sort by score in increasing order: sorted_values_df
  sliding_order <- order(values_df[[sliding_column]])
  sorted_values_df <- values_df[sliding_order, ]
  sorted_values_df[[normalized_value_column]] <-
    sorted_values_df[[value_column]]

  if (use_logvalues) {
    sorted_values_df[[normalized_value_column]] <-
      log10(sorted_values_df[[normalized_value_column]])
  }

  sorted_values_original <- sorted_values_df[[normalized_value_column]]

  # Update values in sorted_values_df
  for (r in seq(1, data_size, step_size)) {
    # Estimate CDF using a window
    id_window_start <- max(1, r + step_size / 2 - window_size / 2)
    id_window_end <- min(data_size, r + step_size / 2 + window_size / 2 - 1)
    y_sorted <- sort(sorted_values_original[id_window_start : id_window_end])
    cdf <- stats::stepfun(y_sorted,
                          (0 : length(y_sorted)) / (length(y_sorted) + 0.0001))

    # Transform
    id_end <- min(data_size, r + step_size - 1)
    values_original <- sorted_values_df[[normalized_value_column]][r : id_end]
    sorted_values_df[[normalized_value_column]][r : id_end] <- stats::qnorm(
      cdf(values_original)
    )

    if (reconstruct) {
      y_quantiles <- stats::quantile(
        values_original, c(0.05, 0.95), na.rm = TRUE
      )
      mn <- mean(y_quantiles)
      std <- diff(y_quantiles) / (stats::qnorm(0.95) - stats::qnorm(0.05))
      sorted_values_df[[normalized_value_column]][r : id_end] <- mn +
        std * sorted_values_df[[normalized_value_column]][r : id_end]
    }
  }

  # Sort back to the original order
  reverse_sliding_order <- order(sliding_order)
  sorted_values_df <- sorted_values_df[reverse_sliding_order, ]

  if (use_logvalues) {
    sorted_values_df[[normalized_value_column]] <-
      10 ** sorted_values_df[[normalized_value_column]]
  }

  sorted_values_df[[normalization_factor_column]] <-
    sorted_values_df[[normalized_value_column]] /
    sorted_values_df[[value_column]]

  return(sorted_values_df)
}


#' Normalize values to standard normal distribution by category
#'
#' @param values_df data frame with value_column and category_column
#' @param value_column name of the value column
#' @param category_column name of the category column
#' @param reconstruct If TRUE, the normalized_value_column is transformed to
#' take the robust mean and robust standard deviation of the value_column.
#' If FALSE, the mean and standard deviation are set to 0 and 1.
#' @param use_logvalues If TRUE, the value column is log10-transformed,
#' normalized, and then transformed back to the original exponential scale.
#' @param normalized_value_column name of the normalized value column
#' @param normalization_factor_column name of the normalization factor column
#' @return data frame with normalized values set in the new column
#' @export
ca_normalize_values_df <- function(
  values_df, value_column = "precursor_quantity", category_column = "condition",
  reconstruct = TRUE, use_logvalues = TRUE,
  normalized_value_column = "normalized_precursor_quantity",
  normalization_factor_column = "normalization_factor"
) {
  values_df[[normalized_value_column]] <- values_df[[value_column]]

  # Get a list of categories
  categories <- unique(values_df[[category_column]])

  if (use_logvalues) {
    values_df[[normalized_value_column]] <-
      log10(values_df[[normalized_value_column]])
  }

  # Update values in values_df
  for (ca in categories) {
    y <- values_df[[normalized_value_column]][
      values_df[[category_column]] == ca
    ]
    cdf <- stats::stepfun(sort(y), (0 : length(y)) / (length(y) + 0.0001))
    values_df[[normalized_value_column]][values_df[[category_column]] == ca] <-
      stats::qnorm(cdf(y))

    if (reconstruct) {
      y_quantiles <- quantile(y, c(0.05, 0.95), na.rm = TRUE)
      mn <- mean(y_quantiles)
      std <- diff(y_quantiles) / (stats::qnorm(0.95) - stats::qnorm(0.05))
      values_df[[normalized_value_column]][
        values_df[[category_column]] == ca
      ] <- mn + std *
        values_df[[normalized_value_column]][values_df[[category_column]] == ca]
    }
  }

  if (use_logvalues) {
    values_df[[normalized_value_column]] <-
      10 ** values_df[[normalized_value_column]]
  }

  values_df[[normalization_factor_column]] <-
    values_df[[normalized_value_column]] /
    values_df[[value_column]]

  return(values_df)
}
