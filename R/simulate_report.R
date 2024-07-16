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
#' @examples
#' report <- simulate_fragment_ion_report(default_params)
#' x1 <- report[
#'     report$condition == "CON1" & report$fragment_id == "FRAG1"
#' ]$precursor_quantity
#' par(mfrow = c(1,2))
#' hist(log10(x1), main = "Histogram", xlab = "Precursor Quantity, log10",
#'      xlim = c(2, 6))
#' qqnorm(log10(x1))
#' qqline(log10(x1))

# Load packages
require(dirmult)

# Set parameters for sampling distributions
default_params <- list(
  # Prefixes for run files
  experiment_prefix = "EXP",
  condition_prefix = "CON",
  fragment_prefix = "FRAG",

  # IDs
  protein_id = "prot",
  precursor_id = "prec",

  # Numbers of samples
  n_experiment = 100,
  n_condition = 6,
  n_replicate = 4,

  # Sampling of precursor mean quantity from Normal distribution
  prec_mean_mean = 5.0,
  prec_mean_std = 0.1,

  # Sampling of data acquisition rate from Beta distribution
  acquisition_beta = c(2.0, 101.0),

  # Sampling of noise from Normal distribution with the mean of zero
  noise_std = 0.1,

  # Sampling of ionization efficiency from Dirichlet distribution
  ionization_dirichet = c(1, 1, 1)
)

default_params$prec_mean_condition_shift <-
  c(0, 0 : (default_params$n_condition - 2)) * log10(2)


#' Generate a simulated fragment ion report
simulate_fragment_ion_report <- function(params, seed = 555) {
  set.seed(seed)

  n_experiment <- params[["n_experiment"]]
  n_condition <- params[["n_condition"]]
  n_replicate <- params[["n_replicate"]]
  n_fragment <- length(params[["ionization_dirichet"]])

  if (length(params[["prec_mean_condition_shift"]]) != n_condition) {
    print(paste("!! The parameter prec_mean_condition_shift must contain",
                n_condition, "numbers."))
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
  w0_values_cond <- c()
  for (a0 in seq(params[["acquisition_beta"]][1],
                 params[["acquisition_beta"]][2],
                 length.out = n_experiment)) {
    w0_values_cond <- c(
      w0_values_cond,
      rbeta(n_condition,
            shape1 = a0,
            shape2 = (4 * log10(a0 - 1) + 1) * (a0 - 1) + 1)
    )
  }

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
  w1_values <- c(t(rdirichlet(
    n_experiment * n_condition * n_replicate,
    params[["ionization_dirichet"]]
  )))
  fragment_peak_area <- precursor_quantity * w1_values
  fragment_peak_area <- pmax(fragment_peak_area, 1)

  fragment_id <- rep(paste0(params[["fragment_prefix"]], 1:n_fragment),
                     n_experiment * n_condition * n_replicate)

  return(data.frame(
    experiment, condition, replicate, protein_id, precursor_id,
    precursor_quantity, fragment_id, fragment_peak_area
  ))
}
