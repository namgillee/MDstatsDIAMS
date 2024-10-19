# MDstatsDIAMS: Multidimensional Statistical Tools for DIA Mass Spectrometry

Available statistical tools:

1. Statistical methods for testing group mean differences in peptide quantities

  - Paired t-test
  - Independent samples t-test
  - Shrinkage-based t-test

2. Generating a simulated fragment ion report


Example:

```
library(MDstatsDIAMS)

report <- simulate_fragment_ion_report(default_params)
test_results <- run_ttests(report)
```


### References

N. Lee, J. Kim, H. Yoo, and H. Yang. 
A shrinkage-based statistical method for testing group mean differences in 
quantitative bottom-up proteomics.
The 23th Human Proteome Organization World Congress (HUPO2024), Dresden,
Germany, October 20--24, 2024. Poster Presentation.
