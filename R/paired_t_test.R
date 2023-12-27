#' Paired t-test for modified peptide using fragment ion peak area
#' 
#' @param x,y  numeric vectors of normalized fragment peak area
#' @examples 
#' paired_t_test(rnorm(10), rnorm(10))
paired_t_test <- function(x, y, var.equal = TRUE, conf.level = 0.95) {
    # check input length
    if (length(x) != length(y)) {
        print("Input vectors must have the same length.")
        return(list())
    }
  
    # missing values
    is_na <- is.na(x) | is.na(y)
    if (sum(!is_na) < 3) {
        print("Input vectors must have at least 3 entries.")
        return(list())
    }
    
    # perform paired t-test
    ttest_out <- t.test(x[!is_na], y[!is_na], var.equal=var.equal, paired=TRUE,
                        conf.level=conf.level)
    result <- list(statsitic = ttest_out$statistic, df=ttest_out$parameter, 
                   p.value=ttest_out$p.value, estimate=ttest_out$estimate)
    return(result)
}
