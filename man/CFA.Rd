% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CFA.R
\name{CFA}
\alias{CFA}
\title{From the EGA's structure to CFA.}
\usage{
CFA(ega.obj, estimator, plot.CFA = TRUE, data, layout = "spring", ...)
}
\arguments{
\item{ega.obj}{An EGA object.}

\item{estimator}{The estimator used in the confirmatory factor analysis. 'WLSMV' is the estimator of choice for ordinal variables. 'ML' or 'WLS' for interval variables.}

\item{plot.CFA}{Logical. Should the CFA structure with its standardized loadings be plot?}

\item{data}{A dataframe with the variables to be used in the analysis.}

\item{...}{Arguments passed to ’cfa’ in lavaan.}
}
\description{
\code{CFA} Verifies the fit of the structure suggested by EGA using confirmatory factor analysis.
}
\examples{
ega.wmt <- EGA(data = wmt2[,7:24])
cfa.wmt <- CFA(ega.obj = ega.wmt, estimator = 'WLSMV', plot.CFA = TRUE, data = wmt2)

ega.intel <- EGA(data = intelligenceBattery[,8:66])
cfa.intel <- CFA(ega.obj = ega.intel, estimator = 'WLSMV', plot.CFA = TRUE, data = intelligenceBattery[,8:66])

\dontrun{
CFA(a, estimator = 'WLSMV', data = data, ...)
}
}
\author{
Hudson F. Golino <hfgolino at gmail.com>
}
\seealso{
\code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA and \code{\link{bootEGA}} to investigate the stability of EGA's estimation via bootstrap.
}

