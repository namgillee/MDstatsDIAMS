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

## Generate a sample report

report <- simulate_fragment_ion_report(default_params, seed = 100)


## SAMPLE PLOT : Create a sample data report and make quantity histogram

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


## Normalization Algorithm
# 1) Stratify data by precursor groups; They can be by
#.   - Condition (not recommended in real data),
#.   - Qvalue
# 2) Compute group mean, group std (robust std is proportional to IQR)
# 3) Normalize each group
# 4) Change its mean and std by the previous values

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
plot(density(y, adjust = 2), type = "l", lty = 1, col = 1, lwd = 2,
     ylim = c(0, 2), main = "", ylab = "Empirical Distribution",
     xlab = "Precursor Quantity, log10", cex.lab = 1.4)
lines(density(yn, adjust = 2), type = "l", lty = 2, col = 2, lwd = 2)
legend("topright", legend = c("Before Normalization", "After Normalization"),
       lty = c(1, 2), col = c(1, 2), lwd = c(2, 2))
dev.off()
