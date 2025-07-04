% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert.R
\name{convert_ms_to_standard}
\alias{convert_ms_to_standard}
\title{Convert MSstats format into standard report}
\usage{
convert_ms_to_standard(ms_report)
}
\arguments{
\item{ms_report}{A report in MSstats format. Required columns are
ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge,
Condition, BioReplicate, Run, Intensity}
}
\value{
A standard report with columns condition, replicate, experiment,
protein_id, precursor_id, fragment_id, fragment_peak_area
}
\description{
Convert MSstats format into standard report.
MSstats format is generated by running MSstatsConvert::*toMSstatsFormat().
Rows are filtered based on fragment peak area, and column names are changed.
}
