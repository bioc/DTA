\name{DTA.generate}
\alias{DTA.generate}
\title{Simulation of DTA experiments}
\description{DTA.generate produces the phenotype matrix and the matrix containing the simulated data according to the given parameters.}
\usage{
DTA.generate(timepoints, tnumber = NULL, plabel = NULL, nrgenes = 5000, mediantime = 12, ccl = 150, delaytime = 0, sdnoise = 0.075, nobias = FALSE, unspecific.LtoU = 0, unspec.LtoU.weighted = FALSE, unspecific.UtoL = 0, unspec.UtoL.weighted = FALSE, truehalflives = NULL, truecomplete = NULL, genenames = NULL, cDTA = FALSE)
}
\arguments{
  \item{timepoints}{Integer vector containing the labeling times for which the samples should be generated.}
  \item{tnumber}{Integer vector containing the number of uridine residues for each gene. If NULL, tnumber is sampled from an F-distribution within the function.}
  \item{plabel}{The labeling efficiency. If NULL, plabel is set to 0.005 within the function. For details, see supplemental material of Sun et al. (see references).}
  \item{nrgenes}{The number of genes the simulated experiment will have (will be cropped if it exceeds the length of tnumber).}
  \item{mediantime}{The median of the randomly drawn half-life distribution.}
  \item{ccl}{The cell cycle length (in minutes).}
  \item{delaytime}{Estimates the delay between the moment of 4sU/4tU labeling and actual incorporation of it into mRNA.}
  \item{sdnoise}{The amount of measurement noise (proportional to expression strength).}
  \item{nobias}{Should a labeling bias be added?}
  \item{unspecific.LtoU}{Proportion of labeled RNAs that unspecifically end up in the unlabeled fraction.}
  \item{unspec.LtoU.weighted}{Should unspecific proportion of labeled to unlabeled depend linearly on the length of the RNA?}
  \item{unspecific.UtoL}{Proportion of unlabeled RNAs that unspecifically end up in the labeled fraction.}
  \item{unspec.UtoL.weighted}{Should unspecific proportion of unlabeled to labeled depend linearly on the length of the RNA?}
  \item{truehalflives}{If the data should be generated using a given half-life distribution, this vector must contain the respective values for each gene.}
  \item{truecomplete}{If the data should be generated using a given expression distribution, this vector must contain the respective values for each gene.}
  \item{genenames}{An optional list of gene names.}
  \item{cDTA}{cDTA = FALSE does not rescale L and U.}
}
\value{
DTA.generate returns a list, containing the following entries
\item{phenomat}{A matrix, containing the design of the experiment as produced by \code{DTA.phenomat}.}
\item{datamat}{A matrix, containing the simulated measurements from U, L and T, according to the design given in phenomat.}
\item{tnumber}{Integer vector containing the number of uridine residues for each gene.}
\item{ccl}{The cell cycle length (in minutes).}
\item{truecomplete}{A vector, containing the true amount of total RNA.}
\item{truelambdas}{A vector, containing the true decay rates.}
\item{truemus}{A vector, containing the true synthesis rates.}
\item{truehalflives}{A vector, containing the true half-lives.}
\item{trueplabel}{The true labeling efficiency. For details, see supplemental material of Miller et al. (see references).}
\item{truear}{The true parameter ar. For details, see supplemental material of Miller et al. (see references).}
\item{truebr}{The true parameter br. For details, see supplemental material of Miller et al. (see references).}
\item{truecr}{The true parameter cr. For details, see supplemental material of Miller et al. (see references).}
\item{truecrbyar}{The true parameter cr/ar. For details, see supplemental material of Miller et al. (see references).}
\item{truecrbybr}{The true parameter cr/br. For details, see supplemental material of Miller et al. (see references).}
\item{truebrbyar}{The true parameter br/ar. For details, see supplemental material of Miller et al. (see references).}
\item{trueLasymptote}{The true parameter asymptote (labeled bias). For details, see supplemental material of Miller et al. (see references).}
\item{trueUasymptote}{The true parameter asymptote (unlabeled bias). For details, see supplemental material of Miller et al. (see references).}
}
\references{
C. Miller, B. Schwalb, K. Maier, D. Schulz, S. Duemcke, B. Zacher, A. Mayer, J. Sydow, L. Marcinowski, L. Doelken, D. E. Martin, A. Tresch, and P. Cramer. Dynamic transcriptome analysis measures rates of mRNA synthesis and decay in yeast. Mol Syst Biol, 7:458, 2011.
M. Sun, B. Schwalb, D. Schulz, N. Pirkl, L. Lariviere, K. Maier, A. Tresch, P. Cramer. Mutual feedback between mRNA synthesis and degradation buffers transcript levels in a eukaryote. Under review.
B. Schwalb, B. Zacher, S. Duemcke, D. Martin, P. Cramer, A. Tresch. Measurement of genome-wide RNA synthesis and decay rates with Dynamic Transcriptome Analysis (DTA/cDTA). Bioinformatics.
}
\author{Bjoern Schwalb \email{schwalb@lmb.uni-muenchen.de}}
\examples{
sim.object = DTA.generate(timepoints=rep(c(6,12),2))

### for control plots set 'check = TRUE' ###

res.sim = DTA.estimate(ratiomethod = "bias",simulation = TRUE,sim.object = sim.object,check = FALSE)
}
\keyword{datagen}
