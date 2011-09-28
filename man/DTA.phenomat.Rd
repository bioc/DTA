\name{DTA.phenomat}
\alias{DTA.phenomat}
\title{Create a phenomat that suits your experiment}
\description{\code{DTA.phenomat} creates a phenomat for a given experimental design, i.e. used labeling times.}
\usage{
DTA.phenomat(timepoints, timecourse = NULL)
}
\arguments{
  \item{timepoints}{The respective labeling times of the measured samples.}
  \item{timecourse}{Vector giving the order for timecourse DTA data.}
}
\value{
A matrix, containing the design of the experiment. Columns are name, fraction (U=unlabebeld, L=labeled, T=total), time and nr (=replicate number). Rows represent individual experiments.
}
\author{Bjoern Schwalb \email{schwalb@lmb.uni-muenchen.de}}
\examples{
DTA.phenomat(c(6,12))
}
\keyword{utilities}

