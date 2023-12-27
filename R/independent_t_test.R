#' Independent (two) sample t-test for modified peptide using aggregated 
#' quantity
#' 
#' @param x,y  numeric vectors of normalized peptide quantity
#' @examples 
#' independent_t_test(rnorm(4), rnorm(4))
independent_t_test <- function(x, y, var.equal = TRUE, conf.level = 0.95) {
    # missing values
    is_na_x <- is.na(x)
    is_na_y <- is.na(y)
    if ((sum(!is_na_x) < 3) | (sum(!is_na_y) < 3)) {
        print("Input vectors must have at least 3 entries.")
        return(list())
    }
  
    # perform independent sample t-test
    ttest_out <- t.test(x[!is_na_x], y[!is_na_y], var.equal=var.equal, 
                        paired=FALSE, conf.level=conf.level)
    result <- list(statsitic = ttest_out$statistic, df=ttest_out$parameter, 
                   p.value=ttest_out$p.value, estimate=ttest_out$estimate)
    return(result)
}
