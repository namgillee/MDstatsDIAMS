% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shrinkage_t_test.R
\name{compute_shrink_on_stdreport}
\alias{compute_shrink_on_stdreport}
\title{Run shrinkage t-test comparing two conditions on a standard report}
\usage{
compute_shrink_on_stdreport(
  report,
  cov_unequal_replicates_column = "cov_unequal_replicates",
  boot_denom_eps = 0.3
)
}
\arguments{
\item{report}{Fragment ion report with the columns experiment, condition,
replicate, protein_id, precursor_id, fragment_id, and fragment_peak_area.
The condition column is supposed to consist of two conditions.}

\item{cov_unequal_replicates_column}{Column name for covariance of input
variables whose replicates are not equal. If NULL, the covariance is set to
zero.}

\item{boot_denom_eps}{A parameter for shrinkage t-test}
}
\value{
A data frame of the test results. Columns are shrinkage_t_test
results, and rows are precursors.
}
\description{
Run shrinkage t-test comparing two conditions on a standard report
}
