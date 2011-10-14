

### median center the columns of a matrix ###


mediancenter = function(mat, 				# matrix
						userows = NULL, 	# the rows to be used
						usecolumns = NULL, 	# the columns to be used
						logscale = TRUE, 	# is the matrix in log-scale ?
						protocol = TRUE, 	# should a protocol be printed ?
						center = FALSE) 	# should the center be 0 (log-scale) or 1 (absolute scale)
{
	if (is.null(userows)){if (is.null(rownames(mat))) {userows = 1:nrow(mat)} else{userows = rownames(mat)}}
	if (is.null(usecolumns)){if (is.null(colnames(mat))) {usecolumns = 1:ncol(mat)} else{usecolumns = colnames(mat)}}
	if (protocol){print("medians before:")
		print(apply(mat[userows,usecolumns],2,median))}
	medians = apply(mat[userows,usecolumns],2,median)
	medianofmedians = median(medians)
	if (center){if (logscale){medianofmedians = 0} else{medianofmedians = 1}}
	if (logscale){mat = t(t(mat) - medians + medianofmedians)} else{mat = t(t(mat) / medians * medianofmedians)}
	if (protocol){print("medians after:")
		print(apply(mat[userows,usecolumns],2,median))}
	return(mat[,usecolumns])
}

medctr = mediancenter


### bias function ###


bias = function(p,x){return(1-(1-p)^x)}


### biasdev function ###


biasdev = function(p,x){return((1-p)^x)}


### LT.plotbias function for log(L/T) ###


LT.plotbias = function(labeled,
		total,
		tnumber,
		plabel = NULL,
		asymptote = NULL,
		main = expression(paste("Labeling Bias ",l[gr]," \n")),
		xlab = expression(paste("Number of Uridine residues ",'#u'[g])),
		ylab = expression(paste("log( ",L[gr]/T[gr]," )")),
		repl = NULL,
		trueplabel = NULL,
		trueasymptote = NULL,
		upper = 500,
		lower = 500,
		ylim=c(-4,4),
		scex = 1,
		correctedlabeling = FALSE,...)
{
	tseq = min(tnumber,na.rm=TRUE):max(tnumber,na.rm=TRUE)
	if (plabel == 1){
		if (correctedlabeling){
			main = expression(paste("(corrected) Labeling Bias ",l[gr]," \n"))
		} else {
			main = expression(paste("(not assessed) Labeling Bias ",l[gr]," \n"))
		}
	}
	heatscatter(tnumber,log(labeled/total),xlab=xlab,ylab=ylab,main=main,xlim=c(0,2000),cor=FALSE,ylim=ylim,cex.main=1.25*scex,cex.lab=1*scex,cex.axis=1*scex,...)
	title(paste(" \n \n \n ( ",repl,"  p = ",signif(plabel,2),"  asymptote = ",signif(asymptote,2)," )",sep=""),col.main="darkgrey",cex.main=0.75*scex)
	points(tseq,log(bias(plabel,tseq))+asymptote,type="l",col="black",lwd=3)
	if (!is.null(trueplabel) & !is.null(trueasymptote)){points(tseq,log(bias(trueplabel,tseq))+asymptote,type="l",col="grey",lwd=3,lty=2)}
}


### UT.plotbias function for log(U/T) ###


UT.plotbias = function(unlabeled,
		total,
		tnumber,
		plabel = NULL,
		asymptote = NULL,
		main = expression(paste("Labeling Bias ",u[gr]," \n")),
		xlab = expression(paste("Number of Uridine residues ",'#u'[g])),
		ylab = expression(paste("log( ",U[gr]/T[gr]," )")),
		repl = NULL,
		trueplabel = NULL,
		trueasymptote = NULL,
		upper = 500,
		lower = 500,
		ylim=c(-4,4),
		scex = 1,
		correctedlabeling = FALSE,
		ratio = NULL,...)
{
	tseq = min(tnumber,na.rm=TRUE):max(tnumber,na.rm=TRUE)
	if (plabel == 1){
		if (correctedlabeling){
			main = expression(paste("(corrected) Labeling Bias ",u[gr]," \n"))
		} else {
			main = expression(paste("(not assessed) Labeling Bias ",u[gr]," \n"))
		}
	}
	heatscatter(tnumber,log(unlabeled/total),xlab=xlab,ylab=ylab,main=main,xlim=c(0,2000),cor=FALSE,ylim=ylim,cex.main=1.25*scex,cex.lab=1*scex,cex.axis=1*scex,...)
	title(paste(" \n \n \n ( ",repl,"  p = ",signif(plabel,2),"  asymptote = ",signif(asymptote,2)," )",sep=""),col.main="darkgrey",cex.main=0.75*scex)
	if (is.null(ratio)){ratio = 1}
	points(tseq,log(1 + ratio*biasdev(plabel,tseq))+asymptote,type="l",col="black",lwd=3)
	if (!is.null(trueplabel) & !is.null(trueasymptote)){points(tseq,log(1 + ratio*biasdev(trueplabel,tseq))+asymptote,type="l",col="grey",lwd=3,lty=2)}
}


### vlength function ###


vlength = function(x) # vector
{
	if (!is.vector(x)){print("ARGUMENT IS NOT A VECTOR !")}
	return(sqrt(sum(x^2)))
}


### vnorm function ###


vnorm = function(x) # vector
{
	if (!is.vector(x)){print("ARGUMENT IS NOT A VECTOR !")}
	return(x/vlength(x))
}


### inter function ###


inter = function(x, 	# x-vector (must contain names)
		y) 	# y-vector (must contain names)
{
	internames = intersect(names(x),names(y))
	if (length(internames) == 0) stop("VECTORS HAVE NO POINT IN COMMON !!!")
	x = x[internames]
	y = y[internames]
	if (!all(names(x) == names(y))) stop("THIS SHOULD NEVER HAPPEN !!!")
	return(list(x=x,y=y))
}


### scor pcor ###


scor = function(x, 		# x-vector
		y,...) 	# y-vector
{
	round(cor(x,y,method="spearman",...),digits=2)
}

pcor = function(x, 		# x-vector
		y,...) 	# y-vector
{
	round(cor(x,y,method="pearson",...),digits=2)
}


### quantile normalization ###


quantnorm = function(mat) # matrix
{
	ordmat = apply(mat,2,function(x){rank(x,ties.method = "first")})
	sortmat = apply(mat,2,sort)
	destin = apply(sortmat,1,mean)
	normmat = destin[ordmat]
	dim(normmat) = dim(mat)
	colnames(normmat) = colnames(mat)
	rownames(normmat) = rownames(mat)
	return(normmat)
}


### piece-wise linear synthesis rates (mu) ###


build.mu = function(d = 0,breaks = NULL,change = 12){
	if (length(breaks)+1 != length(change))stop("You need one more change than break")
	if (is.null(breaks)){d = 0}
	if (d == 0){mu = function(t){return(change[sum(breaks < t)+1])}}
	else {breaks = rep(breaks,each = 2) + rep(c(0,d),length(breaks))
		change = rep(change,each = 2)
		mu = function(t){if (odd(sum(breaks < t)+1)){return(change[sum(breaks < t)+1])}
			else {return(change[sum(breaks < t)] + (change[sum(breaks < t)+2] - change[sum(breaks < t)])*(t - breaks[sum(breaks < t)])/d)}}
	}
	return(mu)
}


### piece-wise linear decay rates (lambda) ###


build.lambda = function(d = 0,breaks = NULL,change = 0.1){
	if (length(breaks)+1 != length(change))stop("You need one more change than break")
	if (is.null(breaks)){d = 0}
	if (d == 0){lambda = function(t){return(change[sum(breaks < t)+1])}}
	else {breaks = rep(breaks,each = 2) + rep(c(0,d),length(breaks))
		change = rep(change,each = 2)
		lambda = function(t){if (odd(sum(breaks < t)+1)){return(change[sum(breaks < t)+1])}
			else {return(change[sum(breaks < t)] + (change[sum(breaks < t)+2] - change[sum(breaks < t)])*(t - breaks[sum(breaks < t)])/d)}}
	}
	return(lambda)
}


### simulate total RNA levels ###


simulate.total = function(
		mu,
		lambda,
		duration,
		lab.duration,
		n = 1,
		alpha = 0
)
{
	ncells = function(t){n*exp(alpha*t)}
	subdivisions = duration
	first = function(t){return(exp(-integrate(Vectorize(lambda),0,t,subdivisions=subdivisions,stop.on.error=FALSE)$value))}
	second = function(t){return(ncells(t)*mu(t)/(alpha + lambda(t)))}
	forthird = function(t){return(ncells(t)*mu(t)*exp(integrate(Vectorize(lambda),0,t,subdivisions=subdivisions,stop.on.error=FALSE)$value))}
	third = function(t){return(integrate(Vectorize(forthird),0,t,subdivisions=subdivisions,stop.on.error=FALSE)$value)}
	total = function(t){return(first(t)*(second(0) + third(t)))}
	total = Vectorize(total)
	sq = seq(0,duration,1)
	tsq = total(sq)
	return(tsq)
}


### TC simulate total RNA levels ###


TC.simulate.total = function(
		mu,
		lambda,
		duration,
		lab.duration,
		storage = function(t){return(0)},
		delay = 6,
		n = 1,
		ccl = NULL
)
{
	if (is.null(ccl)){alpha = 0} else {alpha = log(2)/ccl}
	ncells = function(t){n*exp(alpha*t)}
	subdivisions = duration
	syn = function(t){return(ncells(t)*mu(t))}
	sto = function(t){return(ncells(t)*storage(t))}
	syndelayed = function(t){return(ncells(t-delay)*mu(t-delay))}
	dec = function(t,initial){return(exp(-integrate(Vectorize(lambda),initial,t,subdivisions=subdivisions,stop.on.error=FALSE)$value))}
	forinta = function(t){return((syn(t-delay)+sto(t))*exp(integrate(Vectorize(lambda),delay,t,subdivisions=subdivisions,stop.on.error=FALSE)$value))}
	inta = function(t){return(integrate(Vectorize(forinta),delay,t,subdivisions=subdivisions,stop.on.error=FALSE)$value)}	
	forintb = function(t){return(syn(t-delay)*exp(integrate(Vectorize(lambda),0,t,subdivisions=subdivisions,stop.on.error=FALSE)$value))}
	intb = function(t){return(integrate(Vectorize(forintb),0,t,subdivisions=subdivisions,stop.on.error=FALSE)$value)}
	ANgrt = function(t){
		if (t < delay){
			return(integrate(Vectorize(syn),0,t,subdivisions=subdivisions,stop.on.error=FALSE)$value)
		} else {
			return(integrate(Vectorize(syn),0,t,subdivisions=subdivisions,stop.on.error=FALSE)$value-integrate(Vectorize(syndelayed),delay,t,subdivisions=subdivisions,stop.on.error=FALSE)$value)
		}
	}
	if (is.null(ccl)){CNgrtSS = function(t){return(syn(t)*delay)}} else {CNgrtSS = function(t){return((syn(t)/alpha)*(1-exp(-alpha*delay)))}}
	BNgrt = function(t){
		if (t < delay){
			return(CNgrtSS(0)-integrate(Vectorize(syndelayed),0,t,subdivisions=subdivisions,stop.on.error=FALSE)$value)
		} else {
			return(0)
		}
	}
	ACgrt = function(t){
		if (t < delay){
			return(0)
		} else {
			return(dec(t,delay)*inta(t))
		}
	}
	if (is.null(ccl)){CCgrtSS = function(t){return(((syn(0)/lambda(0))*exp(-alpha*delay)))}} else {CCgrtSS = function(t){return((syn(0)/(alpha-lambda(0)))*exp(-alpha*delay))}}
	CCgrt = function(t){return(dec(t,0)*(CCgrtSS(0) + intb(t)))}
	BCgrt = function(t){
		if (t < delay){
			return(dec(t,0)*(CCgrtSS(0) + intb(t)))
		} else {
			return(CCgrt(delay)*exp(-integrate(Vectorize(lambda),delay,t,subdivisions=subdivisions,stop.on.error=FALSE)$value))
		}
	}
	ANgrt = Vectorize(ANgrt)
	BNgrt = Vectorize(BNgrt)
	ACgrt = Vectorize(ACgrt)
	BCgrt = Vectorize(BCgrt)
	sq = seq(0,duration,1)
	tsq = ANgrt(sq) + ACgrt(sq) + BNgrt(sq) + BCgrt(sq)
	return(tsq)	
}


### sigmoidal up and downregulation ###


sigm = function(x,a,b,c){
	if (a == c){stop("a and c should differ")
	} else if (a < c){
		return((c-a)*1/(1+exp(-x+b))+a)
	} else if (a > c){
		return((a-c)*1/(1+exp(x-b))+c)
	}
}


### LT.error.progression ###


LT.error.progression = function(
		Tmean,Tsd,Lmean,Lsd,Imean = NULL,Isd = NULL,
		timepoint,alpha,
		samplesize = 10000
)
{
	TE = rnorm(samplesize,mean=Tmean,sd=Tsd)
	LE = rnorm(samplesize,mean=Lmean,sd=Lsd)
	if (!is.null(Imean) & !is.null(Isd)){TI = rnorm(samplesize,mean=Imean,sd=Isd)}
	if (!is.null(Imean) & !is.null(Isd)){dr = - alpha - ((1/timepoint)*log((exp(TE) - exp(LE))/exp(TI)))} else {dr = - alpha - ((1/timepoint)*log(1 - exp(LE)/exp(TE)))}
	dr[which(dr <= 0)] = NA
	hl = log(2)/dr
	sr = LE*(alpha + dr)/(exp(alpha*timepoint)-exp(-dr*timepoint))	
	sds = c(sd(dr,na.rm = TRUE),sd(hl,na.rm = TRUE),sd(sr,na.rm = TRUE))
	cvs = c(sd(dr,na.rm = TRUE)/mean(dr,na.rm = TRUE),sd(hl,na.rm = TRUE)/mean(hl,na.rm = TRUE),sd(sr,na.rm = TRUE)/mean(sr,na.rm = TRUE))
	return(list(sds = sds,cvs = cvs))
}


### UT.error.progression ###


UT.error.progression = function(
		Tmean,Tsd,Lmean,Lsd,Umean,Usd,Imean = NULL,Isd = NULL,
		timepoint,alpha,
		samplesize = 10000
)
{
	TE = rnorm(samplesize,mean=Tmean,sd=Tsd)
	LE = rnorm(samplesize,mean=Lmean,sd=Lsd)
	UE = rnorm(samplesize,mean=Umean,sd=Usd)
	if (!is.null(Imean) & !is.null(Isd)){TI = rnorm(samplesize,mean=Imean,sd=Isd)}
	if (!is.null(Imean) & !is.null(Isd)){dr = - alpha - ((1/timepoint)*log(exp(UE)/exp(TI)))} else {dr = - alpha - ((1/timepoint)*log(exp(UE)/exp(TE)))}
	dr[which(dr <= 0)] = NA
	hl = log(2)/dr
	sr = LE*(alpha + dr)/(exp(alpha*timepoint)-exp(-dr*timepoint))	
	sds = c(sd(dr,na.rm = TRUE),sd(hl,na.rm = TRUE),sd(sr,na.rm = TRUE))
	cvs = c(sd(dr,na.rm = TRUE)/mean(dr,na.rm = TRUE),sd(hl,na.rm = TRUE)/mean(hl,na.rm = TRUE),sd(sr,na.rm = TRUE)/mean(sr,na.rm = TRUE))
	return(list(sds = sds,cvs = cvs))
}


### BOTH.error.progression ###


BOTH.error.progression = function(
		Tmean,Tsd,Lmean,Lsd,Umean,Usd,Imean = NULL,Isd = NULL,
		timepoint,alpha,
		samplesize = 10000
)
{
	TE = rnorm(samplesize,mean=Tmean,sd=Tsd)
	LE = rnorm(samplesize,mean=Lmean,sd=Lsd)
	UE = rnorm(samplesize,mean=Umean,sd=Usd)
	if (!is.null(Imean) & !is.null(Isd)){TI = rnorm(samplesize,mean=Imean,sd=Isd)}
	if (!is.null(Imean) & !is.null(Isd)){dr = - alpha - ((1/timepoint)*log(((exp(TE) - exp(LE))/exp(TI) + exp(UE)/exp(TI))/2))} else {dr = - alpha - ((1/timepoint)*log((1 - exp(LE)/exp(TE) + exp(UE)/exp(TE))/2))}
	dr[which(dr <= 0)] = NA
	hl = log(2)/dr
	sr = LE*(alpha + dr)/(exp(alpha*timepoint)-exp(-dr*timepoint))	
	sds = c(sd(dr,na.rm = TRUE),sd(hl,na.rm = TRUE),sd(sr,na.rm = TRUE))
	cvs = c(sd(dr,na.rm = TRUE)/mean(dr,na.rm = TRUE),sd(hl,na.rm = TRUE)/mean(hl,na.rm = TRUE),sd(sr,na.rm = TRUE)/mean(sr,na.rm = TRUE))
	return(list(sds = sds,cvs = cvs))
}



