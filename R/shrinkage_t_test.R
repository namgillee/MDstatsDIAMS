#' Shrinkage-based t-test statistic
#' 
#' @param dat_con1,dat_con2  data matrices of normalized fragment ion peak area, 
#'   each row is a replicate, each column is a fragment ion
#' @param lambda.var_con1,lambda.var_con2  the variance shrinkage intensities
#'   for con1 and con2, respectively. If not specified (default), they are 
#'   estimated by corpcor::cov.shrink function.
#' @return shrinkage-based t-test statistic
#' @examples 
#'   dat_x = matrix(rnorm(12), 4, 3)
#'   dat_y = matrix(rnorm(12), 4, 3)
#'   shrinkage_t_test_statistic(dat_x, dat_y)
shrinkage_t_test_statistic <- function(dat_con1, dat_con2, 
                                       lambda.var_con1, lambda.var_con2) {
  
    # estimate is the sum of diff over fragment ions, and it constitutes numerator
    # of the statistic formula.
    estimate <- sum(
        apply(dat_con1, 2, mean, na.rm = TRUE) -
            apply(dat_con2, 2, mean, na.rm = TRUE),
        na.rm = TRUE)
    
    # covmat is the covariance matrix of diff between each pair of fragment 
    # ions.
    covmat <- cov_diff_shrink(dat_con1, dat_con2, 
                              lambda.var_con1, lambda.var_con2, verbose = FALSE)
    
    # shrinkage statistic
    denom <- sqrt(sum(covmat, na.rm = TRUE))
    result <- NaN
    if (denom > 0) {
        result <- estimate / denom
    }
    
    return(result)
}


#' Shrinkage-based t-test for modified peptide using fragment ion peak area
#' 
#' @param dat_con1,dat_con2  data matrices of normalized fragment ion peak area, 
#'   each row is a replicate, each column is a fragment ion
#' @param lambda.var_con1,lambda.var_con2  the variance shrinkage intensities
#'   for con1 and con2, respectively. If not specified (default), they are 
#'   estimated by corpcor::cov.shrink function.
#' @param conf.level  Confidence level. Default is 0.95
#' @param verbose  TRUE to print messages. Default is FALSE.
#' @examples 
#'   dat_x = matrix(rnorm(12), 4, 3)
#'   dat_y = matrix(rnorm(12), 4, 3)
#'   shrinkage_t_test(dat_x, dat_y)
shrinkage_t_test <- function(dat_con1, dat_con2, 
                             lambda.var_con1, lambda.var_con2,
                             conf.level = 0.95, verbose = FALSE) {
    
    # check number of fragment ions
    if (ncol(dat_con1) != ncol(dat_con2)) {
        if (verbose) {
            print("Number of columns in dat_con1 and dat_con2 must be same.")
        }
        return(list())
    }
    
    # perform shrinkage t-test
    result <- list()
    
    result$statistic <- shrinkage_t_test_statistic(
        dat_con1, dat_con2, lambda.var_con1, lambda.var_con2,
    )
    
    shrinkage_t_test_statistic_wrapper <- function(dat, inds) {
        num_ions <- ncol(dat) / 2
        shrinkage_t_test_statistic(
            dat[inds, 1:num_ions], dat[inds, (num_ions + 1):(2 * num_ions)], 
            lambda.var_con1, lambda.var_con2)
    }
    boot_out <- boot(cbind(dat_con1, dat_con2), 
                     shrinkage_t_test_statistic_wrapper, 
                     R = 500)
    boot_var <- var(boot_out$t, na.rm = TRUE)
    result$df <- ifelse(boot_var <= 1, Inf, 2 * boot_var / (boot_var - 1))
    result$p.value <- 2 * (1 - pt(result$statistic, result$df))

    result$estimate <- sum(
        apply(dat_con1, 2, mean, na.rm = TRUE) - 
            apply(dat_con2, 2, mean, na.rm = TRUE), 
        na.rm = TRUE)

    return(result)
}
