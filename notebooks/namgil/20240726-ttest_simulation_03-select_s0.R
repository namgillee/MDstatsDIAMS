## Load required libraries
library(dplyr)
library(corpcor)
library(reshape)
library(boot)

## Load required functions
source("../../R/independent_t_test.R")
source("../../R/paired_t_test.R")
source("../../R/shrinkage_t_test.R")
source("../../R/simulate_report.R")
source("../../R/utils.R")


## Set initial parameters

params <- default_params

### Explicitly specify parameters for reproducibility
params$n_experiment <- 100
params$n_condition <- 4
params$prec_mean_condition_shift <-
  c(0, 0 : (params$n_condition - 2)) * log10(2)
params$prec_mean_mean <- 5.0
params$prec_mean_std <- 0.1
params$acquisition_beta_fnt <- rbetamixture
params$acquisition_beta_shape1 <- 2 : 21
params$acquisition_beta_shape2 <- 9 * (1 : 20) + 1
params$noise_std <- 0.1
params$ionization_dirichet <- c(2, 2, 2)


## Generate report

report <- simulate_fragment_ion_report(params, seed = 100)


## Run analysis for varying boot.denom.eps values

report$log10_fragment_peak_area <- log10(report$fragment_peak_area)

conditions <- unique(report$condition)
boot_consts <- seq(0.00, 1.0, 0.05)
boot_out <- lapply(boot_consts, function(boot_denom_eps) {

  print(paste0("=========== ", "Boot: ", boot_denom_eps, "==========="))

  result_shrink <- list()

  for (i in 1 : (length(conditions) - 1)) {

    print(paste0("-- ", "condition: ", i, "---"))

    report_twoconds <- report %>%
      filter(condition == conditions[1] | condition == conditions[i + 1])

    result_shrink001 <- report_twoconds %>%
      group_by(experiment, protein_id, precursor_id) %>%
      group_modify(
        ~compute_shrink_on_group(.x, boot_denom_eps = boot_denom_eps)
      ) %>%
      as.data.frame()

    result_shrink[[i]] <- result_shrink001
  }

  result_shrink
})

save(boot_consts, boot_out, file =
    "files/ttest/simulation_03/20240726_run01_boot_const_selecton.RData"
)

all_median_cv <-
  matrix(NA, length(boot_out), 3, dimnames =
    list(NULL, paste0("CON", "1-", 2 : params$n_condition))
  )

for (i in seq_along(boot_out)) {
  all_median_cv[i, 1] <- with(boot_out[[i]][[1]], - median(cv))
  all_median_cv[i, 2] <- with(boot_out[[i]][[2]], - median(cv))
  all_median_cv[i, 3] <- with(boot_out[[i]][[3]], - median(cv))
}

pdf("files/ttest/simulation_03/20240726_run01_selection_of_s0.pdf")
plot(boot_consts, all_median_cv[, 1], ylim = c(-0.6, 1), cex = 0.5,
     xlab = expression(s[0]),
     ylab = "Bootstrap CV", cex.lab = 1.5)
lines(boot_consts, all_median_cv[, 2], col = 1, lty = 1, lwd = 2)
lines(boot_consts, all_median_cv[, 3], col = 2, lty = 2, lwd = 3)
legend("topright",
  legend = c(
    expression(delta == 0),
    expression(delta == 0.3),
    expression(delta == 0.6)
  ),
  pch = c(1, -1, -1, -1),
  col = c(1, 1, 2, 4),
  lty = c(0, 1, 2, 3), lwd = 2, cex = 1.5
)
dev.off()


## Set 1% of the maximum as a threshold

z_diff <- diff(all_median_cv)
eps <- max(abs(all_median_cv)) * 0.01
print(paste("eps:", eps))
###[1]  "eps: 0.00622797562292234"

z_diff_greater_than_eps <- (- rowMeans(z_diff[, 2 : 3]) > eps)

pdf(paste0("files/ttest/simulation_03/",
           "20240726_run01_selection_of_difference_s0.pdf"))
plot(boot_consts[-1], - rowMeans(z_diff[, 2 : 3]),
     ylim = c(-0.05, 0.15),
     xlab = expression(s[0]),
     ylab = "Difference in mean bootstrap CV",
     cex.lab = 1.5)
abline(h = eps, lty = 2, lwd = 3)
dev.off()

idx_z <- min(which(!z_diff_greater_than_eps))
print(boot_consts[1 + idx_z])
###[1]  0.5
