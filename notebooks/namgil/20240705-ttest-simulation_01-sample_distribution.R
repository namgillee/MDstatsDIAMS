# Simulation

## Load required libraries
library(dplyr)
library(corpcor)
library(reshape)
library(boot)

## Load required functions
source("../../R/independent_t_test.R")
source("../../R/normalize.R")
source("../../R/paired_t_test.R")
source("../../R/shrinkage_t_test.R")
source("../../R/simulate_report.R")


## Generate mean quantities using default params

set.seed(100)
n_experiment <- 30000
mu_values_exp <- rnorm(n = n_experiment, mean = 5.0, sd = 0.1)


### Large q-value (eg, qvalue > 1e-5)
w0 <- rbetamixture(n_experiment, 1 * (1 : 20) + 1, 9 * (1 : 20) + 1)

pdf("files/ttest/simulation_01/simulated_w0_histogram_largeq.pdf")
par(mar = c(5.3, 4.6, 4.1, 2.1))
hist(w0, breaks = seq(0, 1, 0.02), freq = FALSE, 
     xlim = c(0, 1), xlab = expression(paste("u"["c"]^"p")),
     main = expression(paste("Simulated ", "u"["c"]^"p")), 
     cex.lab = 2.3, cex.main = 2.0)
dev.off()

pdf("files/ttest/simulation_01/simulated_precursor_mean_qqnorm_largeq.pdf")
par(mar = c(5.8, 4.8, 4.1, 2.1))
qqnorm(mu_values_cond + log10(w0),
       cex.lab = 2.3,
       cex.main = 2.0,
       main = "Simulated Mean Peptide Quantity", 
       cex = 0.5)
qqline(mu_values_cond + log10(w0))
dev.off()


### Small q-value (eg, qvalue <= 1e-5)
w0 <- rbetamixture(n_experiment, 1 * (1 : 20) + 1, 9 * (1 : 20) + 1) * 0.8 + 0.2

pdf("files/ttest/simulation_01/simulated_w0_histogram_smallq.pdf")
par(mar = c(5.3, 4.6, 4.1, 2.1))
hist(w0, breaks = seq(0, 1, 0.02), freq = FALSE, 
     xlim = c(0, 1), xlab = expression(paste("u"["c"]^"p")),
     main = expression(paste("Simulated ", "u"["c"]^"p")), 
     cex.lab = 2.3, cex.main = 2.0)
dev.off()

pdf("files/ttest/simulation_01/simulated_precursor_mean_qqnorm_smallq.pdf")
par(mar = c(5.8, 4.8, 4.1, 2.1))
qqnorm(mu_values_cond + log10(w0),
       cex.lab = 2.3,
       cex.main = 2.0,
       main = "Simulated Mean Peptide Quantity", 
       cex = 0.5)
qqline(mu_values_cond + log10(w0))
dev.off()


## Generate a sample report
## SAMPLE PLOT : Create a sample data report and make quantity histogram

report <- simulate_fragment_ion_report(default_params, seed = 100)

pdf(paste0("files/ttest/simulation_01/",
           "simulated_precursor_quantity_histogram_CON1.pdf"), 7, 5)
x1 <- report[report$condition == "CON1" & report$fragment_id == "FRAG1", ]
par(mfrow = c(1, 2))
hist(log10(x1$precursor_quantity),
     main = "Simulated Data",
     xlab = "Peptide Quantity, log10",
     xlim = c(3, 5),
     breaks = seq(0, 10, 0.1))
qqnorm(log10(x1$precursor_quantity))
qqline(log10(x1$precursor_quantity))
dev.off()

pdf(paste0("files/ttest/simulation_01/",
           "simulated_precursor_quantity_histogram_CONALL.pdf"), 7, 5)
xall <- report[report$fragment_id == "FRAG1", ]
par(mfrow = c(1, 2))
hist(log10(xall$precursor_quantity),
     main = "Simulated Data: All Conditions",
     xlab = "Peptide Quantity, log10",
     xlim = c(3, 7),
     breaks = seq(0, 10, 0.1))
qqnorm(log10(xall$precursor_quantity))
qqline(log10(xall$precursor_quantity))
dev.off()


## SAMPLE PLOT : Generate w0 from a mixture of Beta distribution
## The beta parameter (shape1, shape2) runs from (2, 10) to (21, 181) by
## shape2 = 9 * (shape1 - 1) + 1, which yields beta distributions with the mode
## of 0.1.

pdf("files/ttest/simulation_01/simulated_w0_histogram.pdf", 7, 5)
par(mfrow = c(1, 2))
set.seed(100)
w0 <- rbetamixture(10000, (1 : 20) + 1, 9 * (1 : 20) + 1)
hist(w0, freq = FALSE, xlab = "w0", xlim = c(0, 1), breaks = seq(0, 1, 0.02),
     main = "", cex.lab = 1.4)
hist(log10(w0), freq = FALSE, xlab = "w0, log10", xlim = c(-2.5, 0.0),
     breaks = seq(-2.5, -0.0, 0.05), main = "", cex.lab = 1.4)
lines(density(log10(w0), adjust = 2), col = 2, lwd = 2, lty = 2)
dev.off()


## CA (category) Normalization
## 1) log10-transformed precursor quantity is grouped by condition;
## 2) transform to standard normal distribution; and then
## 3) transform its mean and std by the robust mean and std

normalized_report <- ca_normalize_values_df(
  report, value_column = "precursor_quantity", category_column = "condition",
  reconstruct = TRUE, use_logvalues = TRUE,
  normalized_value_column = "normalized_precursor_quantity"
)

pdf(paste0("files/ttest/simulation_01/",
           "simulated_peptide_quantity_normalization.pdf"), 7, 5)
y <- log10(report[["precursor_quantity"]][
  report$condition == "CON1" &
    report$fragment_id == "FRAG1"
])
yn <- log10(normalized_report[["normalized_precursor_quantity"]][
  normalized_report$condition == "CON1" &
    normalized_report$fragment_id == "FRAG1"
])
hist(y, 20, freq = FALSE, ylim = c(0, 2.1), main = "",
     ylab = "Empirical Distribution",
     xlab = "Precursor Quantity, log10", cex.lab = 1.4)
lines(density(y, adjust = 1), lty = 1, col = 1, lwd = 2)
lines(density(yn, adjust = 1), type = "l", lty = 2, col = 2, lwd = 2)
legend("topright", legend = c("Before Normalization", "After Normalization"),
       lty = c(1, 2), col = c(1, 2), lwd = c(2, 2))
dev.off()
