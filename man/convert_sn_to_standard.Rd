% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert.R
\name{convert_sn_to_standard}
\alias{convert_sn_to_standard}
\title{Convert spectronaut fragment ion report into standard report}
\usage{
convert_sn_to_standard(sn_report, filter_identified = TRUE)
}
\arguments{
\item{sn_report}{A spectronaut fragment ion report. Required columns are
R.Condition, R.Replicate, PG.ProteinGroups, EG.ModifiedSequence,
FG.Charge, EG.Qvalue, F.NormalizedPeakArea,
F.FrgIon, F.FrgLossType, F.Charge, F.ExcludedFromQuantification.}

\item{filter_identified}{TRUE if only identified precursors remain in the
report.}
}
\value{
A standard report with columns condition, replicate, experiment,
protein_id, precursor_id, precursor_qvalue, fragment_id, fragment_peak_area
}
\description{
Convert spectronaut fragment ion report into standard report
}
