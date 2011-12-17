\name{DTA.plot.it}
\alias{DTA.plot.it}
\title{Plots in any format and any quality}
\description{\code{DTA.plot.it} can save plots in any format and any quality in addition to show them in R devices}
\usage{
DTA.plot.it(filename,sw = 1,sh = 1,sres = 1,plotsfkt,ww = 7,wh = 7,pointsize = 12,dev.pointsize = 8,paper = "special",quality = 100,units = "px",bg = "white",fileformat = "jpeg",saveit = FALSE,notinR = FALSE,RStudio = FALSE,addformat = NULL)
}
\arguments{
\item{filename}{Name of the plot to be saved without the format type suffix.}
\item{sw}{Scaling factor of width. Scaled with 480px.}
\item{sh}{Scaling factor of height. Scaled with 480px.}
\item{sres}{Scaling factor of the resolution. Scaled with 72dpi.}
\item{plotsfkt}{Function of plots to be plotted.}
\item{ww}{Width of window. Needed only for plotting in R or if filformat = "pdf" or "ps". See pdf or ps.}
\item{wh}{Height of window. Needed only for plotting in R or if filformat = "pdf" or "ps". See pdf or ps.}
\item{pointsize}{The default pointsize of plotted text, interpreted as big points (1/72 inch) for plots to be saved.}
\item{dev.pointsize}{Pointsize of plotted text, interpreted as big points (1/72 inch) for display in R.}
\item{paper}{Needed only if filformat = "pdf" or "ps". See pdf or ps.}
\item{quality}{Needed only if filformat = "jpeg". See jpeg.}
\item{units}{Needed only if filformat = "jpeg", "png", "bmp" or "tiff". See corresponding function.}
\item{bg}{Backgroundcolor.}
\item{fileformat}{Save the plot as "jpeg", "png", "bmp", "tiff", "ps" or "pdf".}
\item{saveit}{Should plot be saved.}
\item{notinR}{Should plot be not plotted in R.}
\item{RStudio}{For RStudio users. Suppresses the opening of a new device, as RStudio allows only one.}
\item{addformat}{Should plot be saved additionally in another format, "jpeg", "png", "bmp", "tiff", "ps" or "pdf".}
}
\author{Bjoern Schwalb \email{schwalb@lmb.uni-muenchen.de}}
\examples{
plotsfkt = function(){
par(mfrow = c(1,2))
plot(1:10)
plot(10:1)
}
DTA.plot.it(filename = "test",plotsfkt = plotsfkt,saveit = TRUE)

dev.off()
}
\keyword{utilities}

