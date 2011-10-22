\name{DTA.dynamic.generate}
\alias{DTA.dynamic.generate}
\title{Simulation of DTA experiments upon perturbation}
\description{DTA.dynamic.generate produces the phenotype matrix and the matrix containing the simulated data according to the given parameters.}
\usage{
DTA.dynamic.generate(duration = 60,lab.duration = 6,tnumber = NULL,plabel = NULL,nrgenes = 5000,mediantime.halflives = 12,mediantime.synthesisrates = 18,n = 1,ccl = NULL,check = TRUE,plots = FALSE,save.plots = FALSE,folder = NULL,condition = "",addformat = NULL,sdnoise = 0.075,nobias = FALSE,unspecific.LtoU = 0,unspec.LtoU.weighted = FALSE,unspecific.UtoL = 0,unspec.UtoL.weighted = FALSE,mu.values.mat = NULL,mu.breaks.mat = NULL,lambda.values.mat = NULL,lambda.breaks.mat = NULL,truehalflives = NULL,truesynthesisrates = NULL,genenames = NULL)
}
\arguments{
  \item{duration}{duration of the whole time course (min)}
  \item{lab.duration}{labeling duration for single experiments (min)}
  \item{tnumber}{Integer vector containing the number of uridine residues for each gene. If NULL, tnumber is sampled from an F-distribution within the function.}
  \item{plabel}{The labeling efficiency. If NULL, plabel is set to 0.005 within the function. For details, see the vignette.}
  \item{nrgenes}{The number of genes the simulated experiment will have (will be cropped if it exceeds the length of tnumber).}
  \item{mediantime.halflives}{the median of the half life distribution}
  \item{mediantime.synthesisrates}{the median of the synthesis rates distribution (counts/cell/cellcycle)}
  \item{n}{the number of cells N(0)}
  \item{ccl}{The cell cycle length (in minutes).}
  \item{check}{if check=TRUE, control messages will be generated}
  \item{plots}{if plots = TRUE, control plots will be plotted}
  \item{save.plots}{if save.plots = TRUE, control plots will be saved}
  \item{folder}{folder, where to save the plots}
  \item{condition}{to be added to the plotnames}
  \item{addformat}{additional fileformat for plots to be saved}
  \item{sdnoise}{The amount of measurement noise (proportional to expression strength).}
  \item{nobias}{Should a labeling bias be added?}
  \item{unspecific.LtoU}{Proportion of labeled RNAs that unspecifically end up in the unlabeled fraction.}
  \item{unspec.LtoU.weighted}{Should unspecific proportion of labeled to unlabeled depend linearly on the length of the RNA?}
  \item{unspecific.UtoL}{Proportion of unlabeled RNAs that unspecifically end up in the labeled fraction.}
  \item{unspec.UtoL.weighted}{Should unspecific proportion of unlabeled to labeled depend linearly on the length of the RNA?}
  \item{mu.values.mat}{if the data should be generated using given synthesis rates, this matrix must contain the respective values for each gene}
  \item{mu.breaks.mat}{timepoints of synthesis rate changes, this matrix must contain the respective values for each gene, only needed when mu.values.mat is given (one column less than mu.values.mat)}
  \item{lambda.values.mat}{if the data should be generated using given decay rates, this matrix must contain the respective values for each gene}
  \item{lambda.breaks.mat}{timepoints of decay rate changes, this matrix must contain the respective values for each gene, only needed when lambda.values.mat is given (one column less than lambda.values.mat)}
  \item{truehalflives}{If the data should be generated using a given half-life distribution, this vector must contain the respective values for each gene.}
  \item{truesynthesisrates}{If the data should be generated using a given synthesis rates distribution, this vector must contain the respective values for each gene}
  \item{genenames}{An optional list of gene names.}
}
\value{
DTA.dynamic.generate returns a list, containing the following entries
\item{phenomat}{A matrix, containing the design of the experiment as produced by \code{DTA.phenomat}.}
\item{datamat}{A matrix, containing the simulated measurements from U, L and T, according to the design given in phenomat.}
\item{tnumber}{Integer vector containing the number of uridine residues for each gene.}
\item{ccl}{The cell cycle length (in minutes).}
\item{truemus}{A vector, containing the true synthesis rates.}
\item{truemusaveraged}{A vector, containing the true synthesis rates, averaged over the labeling period.}
\item{truelambdas}{A vector, containing the true decay rates.}
\item{truelambdasaveraged}{A vector, containing the true decay rates, averaged over the labeling period.}
\item{truehalflives}{A vector, containing the true half-lives.}
\item{truehalflivesaveraged}{A vector, containing the true half-lives, averaged over the labeling period.}
\item{trueplabel}{The true labeling efficiency. For details, see the vignette.}
\item{truecomplete}{A vector, containing the true amount of total RNA.}
\item{truelambdas}{A vector, containing the true decay rates.}
\item{truemus}{A vector, containing the true synthesis rates.}
\item{truehalflives}{A vector, containing the true half-lives.}
\item{trueplabel}{The true labeling efficiency. For details, see the vignette.}
\item{truear}{The true parameter ar. For details, see the vignette.}
\item{truebr}{The true parameter br. For details, see the vignette.}
\item{truecr}{The true parameter cr. For details, see the vignette.}
\item{truecrbyar}{The true parameter cr/ar. For details, see the vignette.}
\item{truecrbybr}{The true parameter cr/br. For details, see the vignette.}
\item{truebrbyar}{The true parameter br/ar. For details, see the vignette.}
\item{trueLasymptote}{The true parameter asymptote (labeled bias). For details, see the vignette.}
\item{trueUasymptote}{The true parameter asymptote (unlabeled bias). For details, see the vignette.}
}
\references{
C. Miller, B. Schwalb, K. Maier, D. Schulz, S. Duemcke, B. Zacher, A. Mayer, J. Sydow, L. Marcinowski, L. Dolken, D. E. Martin, A. Tresch, and P. Cramer. Dynamic transcriptome analysis measures rates of mRNA synthesis and decay in yeast. Mol Syst Biol, 7:458, 2011.
}
\author{Bjoern Schwalb \email{schwalb@lmb.uni-muenchen.de}}
\examples{
nrgenes = 1000
truesynthesisrates = rf(nrgenes,5,5)*18
steady = rep(1,nrgenes)
shock = 1/pmax(rnorm(nrgenes,mean = 8,sd = 4),1)
induction = pmax(rnorm(nrgenes,mean = 8,sd = 4),1)
changes.mat = cbind(steady,shock,shock*induction)
mu.values.mat = changes.mat*truesynthesisrates
mu.breaks.mat = cbind(rep(12,nrgenes),rep(18,nrgenes))
truehalflives = rf(nrgenes,15,15)*12
truelambdas = log(2)/truehalflives
changes.mat = cbind(steady,shock,shock*induction,steady)
lambda.values.mat = changes.mat*truelambdas
lambda.breaks.mat = cbind(rep(12,nrgenes),rep(18,nrgenes),rep(27,nrgenes))

# it takes 2 min to build sim.object

sim.object = DTA.dynamic.generate(duration = 36,lab.duration = 6,nrgenes = nrgenes,mu.values.mat = mu.values.mat,mu.breaks.mat = mu.breaks.mat,lambda.values.mat = lambda.values.mat,lambda.breaks.mat = lambda.breaks.mat)

res = DTA.dynamic.estimate(simulation = TRUE,sim.object = sim.object,labeling = FALSE,rankpairs = FALSE,assessment = FALSE,correlation = FALSE)

dev.off()
}
\keyword{datagen}
