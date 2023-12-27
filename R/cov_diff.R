#' Covariance of diff's between two fragment ions
#' 
#' cov(d_{ij}, d_{i'j'})
#' 
#' @param dat_con1,dat_con2  data matrix of normalized fragment ion peak area, 
#'   each row is a replicate, each column is a fragment ion
#' @return covariance matrix of fragment ions, ij.
#' @examples 
#'   dat_x = matrix(rnorm(9), 3, 3)
#'   dat_y = matrix(rnorm(9), 3, 3)
#'   cov_diff(dat_x, dat_y)
cov_diff <- function(dat_con1, dat_con2) {
    # check number of fragment ions
    if (ncol(dat_con1) != ncol(dat_con2)) {
        print("Number of columns in dat_con1 and dat_con2 must be same.")
        return(matrix())
    }
  
    # compute covariance matrix
    num_frag <- ncol(dat_con1)
    result <- matrix(NA, num_frag, num_frag)
    for (i in 1:num_frag) {
        for (j in i:num_frag) {
            cov_x <- var(dat_con1[, i], dat_con1[, j], na.rm = TRUE)
            cov_y <- var(dat_con2[, i], dat_con2[, j], na.rm = TRUE)
            num_x <- sum((!is.na(dat_con1[, i])) & (!is.na(dat_con1[, j])))
            num_y <- sum((!is.na(dat_con2[, i])) & (!is.na(dat_con2[, j])))
            
            result[i, j] <- result[j, i] <- cov_x / num_x + cov_y / num_y
        }
    }
    return(result)
}
