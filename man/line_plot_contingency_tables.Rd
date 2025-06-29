% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{line_plot_contingency_tables}
\alias{line_plot_contingency_tables}
\title{Line plot for contingency tables}
\usage{
line_plot_contingency_tables(
  x,
  tables,
  rejected = TRUE,
  scale_factor = 1,
  add_legend = FALSE,
  legend_coord = "topright",
  legend_cex = 1,
  ...
)
}
\arguments{
\item{x}{horizontal coordinates of the line plots}

\item{tables}{list contingency tables. Lengths of x and tables must be equal.}

\item{rejected}{TRUE if number of rejected hypothesis is plotted in y axis.}

\item{scale_factor}{a numeric factor that is multiplied to data values.}

\item{add_legend}{If not FALSE, legend is added in the plot}

\item{legend_coord}{coordinate of the legend}

\item{legend_cex}{cex for the legend}
}
\description{
Line plot for contingency tables
}
\examples{
  report <- simulate_fragment_ion_report(default_params)
  resu <- run_ttests(report, boot_denom_eps = 0.3)
  tables <- compute_contingency_tables(resu, alpha = 0.05)
  x <- default_params$prec_mean_condition_shift[-c(1, 2)]
  line_plot_contingency_tables(
    x, tables[-1], xlab = expression(delta),
    ylab = "1 - Type II error rate", cex.lab = 1.5
  )
}
