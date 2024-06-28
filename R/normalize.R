#' Normalize values to standard normal distribution over sliding windows 
#' 
#' @param values_df data frame with value_column and sliding_column
#' @param value_column name of the value column
#' @param sliding_column name of the score column that sliding window runs over
#' @return data frame with normalized values
sw_normalize_values_df <- function(
    values_df, value_column = "value", sliding_column = "score",
    window_size = 100000, step_size = 50000
) {
  data_size <- nrow(values_df)
    
  if (window_size > data_size) {
    window_size <- data_size
  }
  
  # Sort by score in increasing order: sorted_values_df
  sliding_order <- order(values_df[[sliding_column]])
  sorted_values_df <- values_df[sliding_order, ]
  sorted_values_original <- sorted_values_df[[value_column]]
  
  # Update values in sorted_values_df
  for (r in seq(1, data_size, step_size)) {
      id_window_start <- max(1, r + step_size / 2 - window_size / 2)
      id_window_end <- min(data_size, r + step_size / 2 + window_size / 2 - 1)
      y_sorted <- sort(sorted_values_original[id_window_start : id_window_end])
      cdf <- stepfun(y_sorted, seq(0, 1, 1 / length(y_sorted)))
      
      id_end <- min(data_size, r + step_size - 1)
      sorted_values_df[[value_column]][r : id_end] <- qnorm(
        cdf(sorted_values_df[[value_column]][r : id_end])
      )
  }
  
  # Sort back to the original order
  reverse_sliding_order <- order(sliding_order)
  sorted_values_df <- sorted_values_df[reverse_sliding_order, ]
  
  return(sorted_values_df)
} 
