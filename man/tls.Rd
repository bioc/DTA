\name{tls}
\alias{tls}
\title{Weighted Total Least Square Regression}
\description{Weigthed total least square regression according to Golub and Van Loan (1980) in SIAM J.Numer.Anal Vol 17 No.6}
\usage{
tls(formula, D = NULL, T = NULL, precision = .Machine$double.eps)
}
\arguments{
  \item{formula}{An object of class formula.}
  \item{D}{Diagonal weigth matrix. Default weights are set to 1.}
  \item{T}{Diagonal weigth matrix. Default weights are set to 1.}
  \item{precision}{Smallest possible numeric value on this machine (default).}
}
\value{
\code{tls} returns a lm object.
}
\references{
Golub, G.H. and Van Loan, C.F. (1980). An analysis of the total least squares problem. SIAM J. Numer. Anal., 17:883-893.
}
\author{Sebastian Duemcke \email{duemcke@lmb.uni-muenchen.de}}
\examples{
f = 1.5 # true ratio
a = rnorm(5000)
b = f*a
a = a + rnorm(5000,sd=0.5)
b = b + rnorm(5000,sd=0.5)

coeff.tls = coef(tls(b ~ a + 0))
coeff.lm1 = coef(lm(b ~ a + 0))
coeff.lm2 = 1/coef(lm(a ~ b + 0))

heatscatter(a,b)
abline(0,coeff.lm1,col="red",pch=19,lwd=2)
abline(0,coeff.lm2,col="orange",pch=19,lwd=2)
abline(0,coeff.tls,col="green",pch=19,lwd=2)
abline(0,f,col="grey",pch=19,lwd=2,lty=2)
legend("topleft", c("Least-squares regr. (y ~ x + 0)", "Least-squares regr. (x ~ y + 0)", "Total Least-squares regr.", "True ratio"), col=c("red", "orange", "green", "grey"), lty=c(1,1,1,2), lwd=2)

results = c(coeff.tls,coeff.lm1,coeff.lm2)
names(results) = c("coeff.tls","coeff.lm1","coeff.lm2")
print(results)
}
\keyword{methods}

