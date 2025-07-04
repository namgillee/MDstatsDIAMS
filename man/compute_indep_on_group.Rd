% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/independent_t_test.R
\name{compute_indep_on_group}
\alias{compute_indep_on_group}
\title{Wrapper function for independent samples t-test}
\usage{
compute_indep_on_group(
  groupdf,
  column_x = "Log10Quantity.x",
  column_y = "Log10Quantity.y"
)
}
\arguments{
\item{groupdf}{A subset of precursor report which contains two columns for
log10-transformed precursor quantities of conditions x and y}

\item{column_x}{column name for x}

\item{column_y}{column name for y}
}
\value{
a data frame of the independent samples t-test results
}
\description{
Wrapper function for independent samples t-test
}
