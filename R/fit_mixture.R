#' Exponential-beta density function with base 10
#'
#' @param x  non-positive values
#' @param shape1,shape2  shape parameters
#' @return  density function values
dexpbeta10 <- function(x, shape1 = 1, shape2 = 1) {
  log_f <- - lbeta(shape1, shape2) + x * (shape1 - 1) * log(10) +
    (shape2 - 1) * log1p(- 10^x) + x * log(10) + log(log(10))
  return(exp(log_f))
}


#' Fit a mixture of normal and expbeta distributions by EM algorithm
#'
#' @param x  a vector of data values
#' @param init_params  a vector of initial parameter values
#' @param max_iter maximum number of iterations
#' @param index_fixed_params  index for fixed parameters. A subset of {1,..., 4}
#' @param xlim_max an assumed maximum of x range
#' @return  a vector of fitted parameters and component weights
fit_mixture_normal_expbeta <- function(
  x, init_params, max_iter = 30, index_fixed_params = c(1), xlim_max = NULL
) {

  if (is.null(xlim_max)) {
    xlim_max <- max(x) + 0.5
  }

  fitted_params <- init_params
  compo_weights <- c(0.9, 0.1)

  for (iter in 1 : max_iter) {
    old_fitted_params <- fitted_params

    # Evaluate density functions
    weighted_density_values <- cbind(
      stats::dnorm(x, fitted_params[1], fitted_params[2]),
      dexpbeta10(x - xlim_max, fitted_params[3], fitted_params[4])
    ) %*% diag(compo_weights)

    # Posterior expectation of latent variables
    prob_compo_inclusion <- 1 / (1 + cbind(
      weighted_density_values[, 2] / weighted_density_values[, 1],
      weighted_density_values[, 1] / weighted_density_values[, 2]
    ))
    prob_compo_inclusion[is.nan(prob_compo_inclusion)] <- 0.5

    # Update component weights
    compo_weights <- apply(prob_compo_inclusion, 2, mean, na.rm = TRUE)

    # Weighted MLE for normal density
    if (!(1 %in% index_fixed_params)) {
      fitted_params[1] <- sum(x * prob_compo_inclusion[, 1]) /
        sum(prob_compo_inclusion[, 1])
    }
    if (!(2 %in% index_fixed_params)) {
      fitted_params[2] <- sqrt(
        sum((x - fitted_params[1])^2 * prob_compo_inclusion[, 1]) /
          sum(prob_compo_inclusion[, 1])
      )
    }

    # Weighted MLE for expbeta density
    fit_resu <- MASS::fitdistr(
      x = rep(
        x - xlim_max,
        round(
          prob_compo_inclusion[, 2] / sum(prob_compo_inclusion[, 2]) * 10000
        )
      ),
      densfun = dexpbeta10,
      start = list(shape1 = fitted_params[3], shape2 = fitted_params[4]),
      lower = c(0.1, 1),
      upper = c(40, Inf)
    )

    fitted_params[3] <- fit_resu$estimate[1]
    fitted_params[4] <- fit_resu$estimate[2]

    if (sum((fitted_params - old_fitted_params)^2) / sum(old_fitted_params^2) <
          1e-8) {
      break
    }
  }

  names(fitted_params) <- c("mean", "sd", "shape1", "shape2")
  names(compo_weights) <- c("component1", "component2")
  return(c(fitted_params, compo_weights, xlim_max = xlim_max))
}
