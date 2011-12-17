\name{Sc.phenomat.dynamic}
\docType{data}
\alias{Sc.phenomat.dynamic}
\title{Design of the Saccharomyces Cerevisiae salt stress DTA experiment from Miller et al.}
\description{
The phenotype matrix \code{Sc.phenomat.dynamic} contains information about the experimental design. It is comprised of the filename, the type of RNA fraction measured (T, U or L), the labeling time, the replicate number and an additional number indicating the timecourse order.
}
\usage{Sc.phenomat.dynamic}
\format{
The phenomat is a matrix comprised of the file name, the type of RNA fraction mesasured (T, U or L, \code{fraction} column), the labeling time (\code{time},\code{timeframe} column), the replicate number (\code{nr} column) and a number indicating the timecourse order (\code{timecourse} column). Rows in this matrix represent the individual experiments.
}
\source{
C. Miller, B. Schwalb, K. Maier, D. Schulz, S. Duemcke, B. Zacher, A. Mayer, J. Sydow, L. Marcinowski, L. Doelken, D. E. Martin, A. Tresch, and P. Cramer. Dynamic transcriptome analysis measures rates of mRNA synthesis and decay in yeast. Mol Syst Biol, 7:458, 2011.
}
\keyword{datasets}
