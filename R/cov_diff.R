#' Covariance of diff between each pair of fragment ions.
#' 
#' A diff for an ion i is defined by d_{i} = xbar_con1{i} - xbar_con2{i}, where
#' xbar_con1 and xbar_con2 are the sample means for condition 1 and 2, 
#' respectively. A covariance of diff between fragment ion i and j is 
#'    cov(d_{i}, d_{j})
#'    = cov(xbar_con1{i}, xbar_con1{j}) + cov(xbar_con2{i}, xbar_con2{j})
#' 
#' @param dat_con1,dat_con2  data matrices of normalized fragment ion peak area, 
#'   each row is a replicate, each column is a fragment ion
#' @param verbose  TRUE to print messages. Default is FALSE.
#' @return covariance matrix of mean difference
#' @examples 
#'   dat_x = matrix(rnorm(12), 4, 3)
#'   dat_y = matrix(rnorm(12), 4, 3)
#'   cov_diff(dat_x, dat_y)
cov_diff <- function(dat_con1, dat_con2, verbose = FALSE) {
    # check number of fragment ions
    if (ncol(dat_con1) != ncol(dat_con2)) {
      if (verbose) {
        print("Number of columns in dat_con1 and dat_con2 must be same.")
      }
      return(matrix())
    }
  
    # compute covariance matrix
    num_frag <- ncol(dat_con1)
    result <- matrix(NA, num_frag, num_frag)
    for (i in 1:num_frag) {
        for (j in i:num_frag) {
            cov_con1 <- var(dat_con1[, i], dat_con1[, j], na.rm = TRUE)
            num_con1 <- sum((!is.na(dat_con1[, i])) & (!is.na(dat_con1[, j])))
            cov_con2 <- var(dat_con2[, i], dat_con2[, j], na.rm = TRUE)
            num_con2 <- sum((!is.na(dat_con2[, i])) & (!is.na(dat_con2[, j])))
            if ((!is.na(cov_con1)) && (!is.na(cov_con2))) {
                result[i, j] <- result[j, i] <- 
                  cov_con1 / num_con1 + cov_con2 / num_con2
            }
        }
    }
    return(result)
}
