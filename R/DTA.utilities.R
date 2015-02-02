

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


### categorize vector to quantiles ###


toquantiles = function(t,k = 10){
	q = quantile(t,seq(0,1,1/k),na.rm = TRUE)
	dt = t
	dt[which(t <= q[2])] = 1
	for (i in 2:(k-1)){dt[which(t > q[i] & t <= q[i+1])] = i}
	dt[which(t > q[k])] = k
	return(dt)
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


### gamma.variance ###


gamma.variance = function(residuals){
	m.1 = sum(residuals)/length(residuals)
	m.2 = sum(residuals^2)/length(residuals)
	k = m.1^2/(m.2 - m.1^2)
	theta = (m.2 - m.1^2)/m.1
	v = k*theta^2
	return(c(v,k,theta))
}


### likelihood ###


likelihood = function(x,mean.x,sd.x,loess.x,gamma.var){
	prod(dnorm(x = x,mean = mean.x,sd = sd.x))*dgamma(x = sd.x, shape = (loess.x^2)/gamma.var, scale = gamma.var/loess.x)
}


### isreal function ###


is.real = function(x)
{
    if (!is.numeric(x)) stop("This is not a numeric vector")
    vec = !(is.na(x) | is.nan(x) | (x==Inf) | (x==-Inf))
    return(vec)
}


### gridfkt function for quadratic plots with log2 folds ###


gridfkt = function(lim, 			# limits for quadratic plots
lty=2, 				# linetype to be used
lwd=1,              # linewidth to be used
addlines = TRUE,    # should fold lines be added
folds = TRUE) 		# should labels be added to fold lines
{
    count = lim[2]-lim[1]+1
    sequ = lim[2]:lim[1]
    colpal = colorpalette("rdbu",count)
    abline(h=0,lty=lty,lwd=lwd)
    abline(v=0,lty=lty,lwd=lwd)
    if (addlines){for (i in 1:count){if (!(sequ[i] == 0)){abline(h=sequ[i],lty=lty,lwd=lwd,col=colpal[i])}}}
    if (folds){text(lim[1]+1,-1,"2 fold",col=colpal[which(sequ == -1)],adj=c(-.1,-.1))}
    if (folds){text(lim[1]+1,-2,"4 fold",col=colpal[which(sequ == -2)],adj=c(-.1,-.1))}
    if (folds){text(lim[1]+1,-3,"8 fold",col=colpal[which(sequ == -3)],adj=c(-.1,-.1))}
    if (addlines){for (i in 1:count){if (!(sequ[i] == 0)){abline(v=sequ[i],lty=lty,lwd=lwd,col=colpal[i])}}}
    if (folds){text(1,lim[2]-1,"2 fold",col=colpal[which(sequ == 1)],adj=c(.3,-.1))}
    if (folds){text(2,lim[2]-1,"4 fold",col=colpal[which(sequ == 2)],adj=c(.3,-.1))}
    if (folds){text(3,lim[2]-1,"8 fold",col=colpal[which(sequ == 3)],adj=c(.3,-.1))}
}

