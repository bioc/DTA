\name{DTA.normalize}
\alias{DTA.normalize}
\title{cDTA normalization procedure.}
\description{\code{DTA.normalize} can normalize expression values from a certain species to the median of values from a reference species.}
\usage{
DTA.normalize(mat,reliable = NULL,logscale = FALSE,protocol = FALSE,center = FALSE)
}
\arguments{
\item{mat}{Expression matrix.}
\item{reliable}{The rows to be used, i.e. identifiers of the reference species to normalize on.}
\item{logscale}{Is the matrix in log-scale ?}
\item{protocol}{Should a protocol be printed ?}
\item{center}{Should the center be 0 (log-scale) or 1 (absolute scale). Otherwise the median of the medians is taken.}
}
\author{Bjoern Schwalb \email{schwalb@lmb.uni-muenchen.de}}
\examples{
### see vignette examples or reference:
### B. Schwalb, B. Zacher, S. Duemcke, D. Martin, P. Cramer, A. Tresch.
### Measurement of genome-wide RNA synthesis and decay rates with Dynamic Transcriptome Analysis (DTA/cDTA). Bioinformatics, in revision.
}
\references{
M. Sun, B. Schwalb, D. Schulz, N. Pirkl, L. Lariviere, K. Maier, A. Tresch, P. Cramer. Mutual feedback between mRNA synthesis and degradation buffers transcript levels in a eukaryote. Genome Research, in revision.
B. Schwalb, B. Zacher, S. Duemcke, D. Martin, P. Cramer, A. Tresch. Measurement of genome-wide RNA synthesis and decay rates with Dynamic Transcriptome Analysis (DTA/cDTA). Bioinformatics, in revision. 
}
\keyword{utilities}

