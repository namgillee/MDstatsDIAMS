# MDstatsDIAMS: Multidimensional Statistical Tools for DIA Mass Spectrometry

## Available statistical tools

1. Testing group mean differences in peptide quantities

  - ROTS
  - MSstatsLiP
  - Paired t-test
  - Independent samples t-test
  - Shrinkage-based t-test


## Examples

The package can be installed conveniently from GitHub:

```
if (!("devtools" %in% installed.packages())) 
  install.packages("devtools")

devtools::install_github("namgillee/MDstatsDIAMS")
library(MDstatsDIAMS)
```

One can generate a simulated fragment ion report using a hierarchical graphical
model, and run $t$-tests:

```
params <- default_params
params$n_experiment <- 100
report <- simulate_fragment_ion_report(params)
test_results <- run_ttests(report, boot_denom_eps = 0.3)
compute_contingency_tables(test_results)
```

Instead, one can import a Spectronaut fragment ion report, convert it into a
standard format, and run $t$-tests:

```
sample_report <- read.delim("data/sample_spectronaut_report.tsv")
report <- convert_sn_to_standard(sample_report)
test_results <- run_ttests(report, boot_denom_eps = 0.3, base_condition = "DMSO")
compute_contingency_tables(test_results)
```


### References

N. Lee, J. Kim, H. Yoo, and H. Yang. 
A shrinkage-based statistical method for testing group mean differences in 
quantitative bottom-up proteomics.
The 23th Human Proteome Organization World Congress (HUPO2024), Dresden,
Germany, October 20--24, 2024. Poster Presentation.
