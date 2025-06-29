% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert.R
\name{convert_sk_to_standard}
\alias{convert_sk_to_standard}
\title{Convert Skyline transition report into standard report}
\usage{
convert_sk_to_standard(sk_report, annotation)
}
\arguments{
\item{sk_report}{A report in Skyline transition. Required columns are
Peptide, Protein, Fragment.Ion, Replicate, Area. The Replicate column is
unique to each Run.}

\item{annotation}{Data frame with columns Condition, Replicate, Run}
}
\value{
A standard report with columns condition, replicate, experiment,
protein_id, precursor_id, fragment_id, fragment_peak_area
}
\description{
Convert Skyline transition report into standard report.
Rows are filtered based on fragment peak area, and column names are changed.
}
