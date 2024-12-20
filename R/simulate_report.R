#' Define default parameters for sampling distributions
#' @export
default_params <- list(
  # Prefixes for run files
  experiment_prefix = "EXP",
  condition_prefix = "CON",
  fragment_prefix = "FRAG",

  # IDs
  protein_id = "prot",
  precursor_id = "prec",

  # Numbers of samples
  n_experiment = 500,
  n_condition = 4,
  n_replicate = 4,

  # Sampling of precursor mean quantity from Normal distribution
  prec_mean_mean = 5.0,
  prec_mean_std = 1,

  # Sampling of data acquisition rate from Beta distribution
  acquisition_beta_fnt = rbeta,
  acquisition_beta_shape1 = 2,
  acquisition_beta_shape2 = 10,

  # Sampling of noise from Normal distribution with the mean of zero
  noise_std = 0.1,

  # Sampling of ionization efficiency from Dirichlet distribution
  ionization_dirichlet_fnt = dirmult::rdirichlet,
  ionization_dirichlet_alpha = c(2, 2, 2),
  ionization_cor_bet_condition = 0.89,

  # Mean shift values added to the precursor mean quantity over conditions
  prec_mean_condition_shift = c(0, 0 : 2) * log10(2)
)


#' Generate a simulated fragment ion report based on a hierarchical graphical
#' model
#'
#' A fragment ion report consists of the following columns:
#'   - experiment, condition, replicate, protein_id, precursor_id,
#'     precursor_quantity, fragment_id, fragment_peak_area
#'
#' The precursor_quantity and fragment_peak_area columns are raw values, i.e.,
#' not log-transformed.
#'
#' @param params A list of parameters in the same format as default_params
#' @param seed An integer for random seed
#' @return A data frame of a fragment ion report
#' @examples
#' report <- simulate_fragment_ion_report(default_params)
#' x1 <- report[
#'   report$condition == "CON1" & report$fragment_id == "FRAG1",
#' ]$precursor_quantity
#' par(mfrow = c(1,2))
#' hist(log10(x1), main = "Histogram", xlab = "Precursor Quantity, log10",
#'      xlim = c(2, 6))
#' qqnorm(log10(x1))
#' qqline(log10(x1))
#' @export
simulate_fragment_ion_report <- function(params, seed = 100) {
  set.seed(seed)

  n_experiment <- params[["n_experiment"]]
  n_condition <- params[["n_condition"]]
  n_replicate <- params[["n_replicate"]]
  n_fragment <- length(params[["ionization_dirichlet_alpha"]])

  if (length(params[["prec_mean_condition_shift"]]) != n_condition) {
    print(
      paste(
        "!! The parameter prec_mean_condition_shift must contain",
        n_condition,
        "numbers."
      )
    )
    return(data.frame())
  }

  # Create (n_experiment * n_condition * n_replicate * n_fragment) vectors

  ## 1) Generate precursor mean values
  experiment <- rep(
    paste0(params[["experiment_prefix"]], 1:n_experiment),
    each = n_condition * n_replicate * n_fragment
  )
  condition <- rep(
    rep(
      paste0(params[["condition_prefix"]], 1:n_condition),
      each = n_replicate * n_fragment
    ),
    times = n_experiment
  )
  protein_id <- params[["protein_id"]]
  precursor_id <- params[["precursor_id"]]

  ### mu values for each experiment
  mu_values_exp <- rnorm(n_experiment,
                         mean = params[["prec_mean_mean"]],
                         sd = params[["prec_mean_std"]])
  ### mu values for each experiment * condition
  mu_values_cond <-
    rep(mu_values_exp, each = n_condition) +
    rep(params[["prec_mean_condition_shift"]], times = n_experiment)
  ### data acquisition rates for each experiment * condition
  w0_values_cond <- matrix(nrow = n_condition, ncol = n_experiment)
  rbeta_fnt <- params[["acquisition_beta_fnt"]]
  for (id_cond in 1:n_condition) {
    w0_values_cond[id_cond, ] <- rbeta_fnt(
      n_experiment,
      params[["acquisition_beta_shape1"]],
      params[["acquisition_beta_shape2"]]
    )
  }
  dim(w0_values_cond) <- NULL

  ## 2) Generate estimated precursor quantity
  ### noise for each replicate
  noise_values_rep <-
    rnorm(
      n_experiment * n_condition * n_replicate, mean = 0,
      sd = params[["noise_std"]]
    )
  ### estimated precursor quantity for each replicate
  precursor_quantity_rep <-
    rep(mu_values_cond + log10(w0_values_cond), each = n_replicate) +
    noise_values_rep

  precursor_quantity <- rep(10 ** precursor_quantity_rep, each = n_fragment)
  replicate <- rep(
    rep(1:n_replicate, each = n_fragment),
    times = n_experiment * n_condition
  )

  ## 3) Fragment peak area
  rdirichlet_fnt <- params[["ionization_dirichlet_fnt"]]
  ### generate un-correlated w1
  w1_values <- rdirichlet_fnt(
    n_replicate * n_condition * n_experiment,
    params[["ionization_dirichlet_alpha"]]
  )
  ### create correlation in w1 between conditions
  cor_w1 <- params[["ionization_cor_bet_condition"]]
  if ((cor_w1 < 0) || (cor_w1 > 1)) {
    print(
      paste(
        "In simulate_fragment_ion_report():",
        "$ionization_cor_bet_condition must be between 0 and 1.",
        "But", cor_w1, "is given. Skip creating correlation in w1."
      )
    )
  } else if ((cor_w1 > 0) && (cor_w1 <= 1)) {
    ### compute weight from correlation in w1 between conditions
    ### solve: cor_w1 = corr(X, Z), where Z = (1-weight)*Y + weight*X,
    ###   where X and Y are independent.
    weight <- ifelse(
      cor_w1 == 1 / sqrt(2),
      0.5,
      (cor_w1^2 - sqrt(cor_w1^2 * (1 - cor_w1^2))) / (2 * cor_w1^2 - 1)
    )
    ### reshape
    dim(w1_values) <- c(n_replicate * n_condition, n_experiment * n_fragment)
    w1_values <- t(w1_values) #n_experiment * n_fragment x n_replicate * n_con..
    dim(w1_values) <- c(n_experiment * n_fragment * n_replicate, n_condition)
    ### compute weighted sum of w1 values
    w1_values <- (1 - weight) * w1_values + weight * w1_values[, 1]
    ## reshape
    dim(w1_values) <- c(n_experiment * n_fragment, n_replicate * n_condition)
    w1_values <- t(w1_values) #n_replicate * n_condition x n_experiment * n_fr..
    dim(w1_values) <- c(n_replicate * n_condition * n_experiment, n_fragment)
  }
  ### reshape
  w1_values <- c(t(w1_values))

  fragment_peak_area <- precursor_quantity * w1_values
  fragment_peak_area <- pmax(fragment_peak_area, 1)

  fragment_id <- rep(paste0(params[["fragment_prefix"]], 1:n_fragment),
                     n_experiment * n_condition * n_replicate)

  return(data.frame(
    experiment, condition, replicate, protein_id, precursor_id,
    precursor_quantity, fragment_id, fragment_peak_area
  ))
}
