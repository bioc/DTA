

### bias function ###


bias = function(p,x)
{
	return(1-(1-p)^x)
}


### biasdev function ###


biasdev = function(p,x)
{
	return((1-p)^x)
}


### vlength function ###


vlength = function(x)
{
	return(sqrt(sum(x^2)))
}


### vnorm function ###


vnorm = function(x)
{
	return(x/vlength(x))
}


### inter function ###


inter = function(x,y)
{
	internames = intersect(names(x),names(y))
	x = x[internames]
	y = y[internames]
	return(list(x=x,y=y))
}


### scor function ###


scor = function(x,y,...)
{
	round(cor(x,y,method="spearman",...),digits=2)
}


### pcor function ###


pcor = function(x,y,...)
{
	round(cor(x,y,method="pearson",...),digits=2)
}


### median center function ###


mediancenter = function(mat, 				# matrix
		userows = NULL, 	# the rows to be used
		usecolumns = NULL, 	# the columns to be used
		logscale = TRUE, 	# is the matrix in log-scale ?
		protocol = TRUE, 	# should a protocol be printed ?
		center = FALSE	 	# should the center be 0 (log-scale) or 1 (absolute scale)
)
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


### quantnorm function ###


quantnorm = function(mat)
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


### LT.error.progression function ###


LT.error.progression = function(
		Tmean,Tsd,Lmean,Lsd,Imean = NULL,Isd = NULL,
		timepoint,alpha,
		samplesize = 10000
)
{
	TE = rnorm(samplesize,mean=Tmean,sd=Tsd)
	LE = rnorm(samplesize,mean=Lmean,sd=Lsd)
	if (!is.null(Imean) & !is.null(Isd)){TI = rnorm(samplesize,mean=Imean,sd=Isd)}
	options(warn = -1)
	if (!is.null(Imean) & !is.null(Isd)){dr = - alpha - ((1/timepoint)*log((exp(TE) - exp(LE))/exp(TI)))} else {dr = - alpha - ((1/timepoint)*log(1 - exp(LE)/exp(TE)))}
	options(warn = 0)
	dr[which(dr <= 0)] = NA
	hl = log(2)/dr
	sr = LE*(alpha + dr)/(exp(alpha*timepoint)-exp(-dr*timepoint))	
	sds = c(sd(dr,na.rm = TRUE),sd(hl,na.rm = TRUE),sd(sr,na.rm = TRUE))
	cvs = c(sd(dr,na.rm = TRUE)/mean(dr,na.rm = TRUE),sd(hl,na.rm = TRUE)/mean(hl,na.rm = TRUE),sd(sr,na.rm = TRUE)/mean(sr,na.rm = TRUE))
	return(list(sds = sds,cvs = cvs))
}


### UT.error.progression function ###


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
	options(warn = -1)
	if (!is.null(Imean) & !is.null(Isd)){dr = - alpha - ((1/timepoint)*log(exp(UE)/exp(TI)))} else {dr = - alpha - ((1/timepoint)*log(exp(UE)/exp(TE)))}
	options(warn = 0)
	dr[which(dr <= 0)] = NA
	hl = log(2)/dr
	sr = LE*(alpha + dr)/(exp(alpha*timepoint)-exp(-dr*timepoint))	
	sds = c(sd(dr,na.rm = TRUE),sd(hl,na.rm = TRUE),sd(sr,na.rm = TRUE))
	cvs = c(sd(dr,na.rm = TRUE)/mean(dr,na.rm = TRUE),sd(hl,na.rm = TRUE)/mean(hl,na.rm = TRUE),sd(sr,na.rm = TRUE)/mean(sr,na.rm = TRUE))
	return(list(sds = sds,cvs = cvs))
}


### BOTH.error.progression function ###


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
	options(warn = -1)
	if (!is.null(Imean) & !is.null(Isd)){dr = - alpha - ((1/timepoint)*log(((exp(TE) - exp(LE))/exp(TI) + exp(UE)/exp(TI))/2))} else {dr = - alpha - ((1/timepoint)*log((1 - exp(LE)/exp(TE) + exp(UE)/exp(TE))/2))}
	options(warn = 0)
	dr[which(dr <= 0)] = NA
	hl = log(2)/dr
	sr = LE*(alpha + dr)/(exp(alpha*timepoint)-exp(-dr*timepoint))	
	sds = c(sd(dr,na.rm = TRUE),sd(hl,na.rm = TRUE),sd(sr,na.rm = TRUE))
	cvs = c(sd(dr,na.rm = TRUE)/mean(dr,na.rm = TRUE),sd(hl,na.rm = TRUE)/mean(hl,na.rm = TRUE),sd(sr,na.rm = TRUE)/mean(sr,na.rm = TRUE))
	return(list(sds = sds,cvs = cvs))
}



