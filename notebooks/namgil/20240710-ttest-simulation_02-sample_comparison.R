# Simulation comparison of statistical testing methods

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
source("../../R/utils.R")


## Comparison on the report using default settings

report <- simulate_fragment_ion_report(default_params, seed = 100)

resu_default <- run_ttests(report, boot_denom_eps = 0.5)

tables_default <- compute_contingency_tables(resu_default, alpha = 0.05)


## Comparison on the ca-normalized report

ca_normalized_report <- ca_normalize_values_df(
  report, value_column = "precursor_quantity", category_column = "condition",
  reconstruct = TRUE, use_logvalues = TRUE,
  normalized_value_column = "normalized_precursor_quantity"
) %>%
  mutate(fragment_peak_area = normalization_factor * fragment_peak_area)

resu_ca_normalized <- run_ttests(ca_normalized_report, boot_denom_eps = 0.5)
  
tables_ca_normalized <- compute_contingency_tables(
  resu_ca_normalized, alpha = 0.05
)


## Print contingency tables

print("---------------------------------------")
print("Default setting:")
print(tables_default)

print("---------------------------------------")
print("ca-normalization:")
print(tables_ca_normalized)
