

### DTA.dynamic.singlegenerate generates one experiment consisting of total, labeled and unlabeled measurements of a simulated time course. DATA IS RETURNED ON THE ABSOLUTE SCALE !!!! ###


DTA.dynamic.singlegenerate = function(
		mu.values,								# vector giving the synthesis rate(s) (mu)
		mu.breaks = NULL,						# vector giving the timepoints the synthesis rate changes (mu)
		lambda.values,							# vector giving the decay rate(s) (lambda)
		lambda.breaks = NULL,					# vector giving the timepoints the decay rate changes (lambda)
		n = 1,   								# the number of cells N(0)
		ccl = NULL, 							# the cell cycle length of the cells
		duration = 60,							# duration of the whole time course (min)
		lab.duration = 6,						# labeling duration for single experiments (min)
		check = TRUE, 							# if check=TRUE, control messages and plots will be generated
		plots = FALSE, 							# if plots = TRUE, control plots will be plotted
		cols = c("darkred","darkblue","black"),	# colors for plot
		save.plots = FALSE, 					# if save.plots = TRUE, control plots will be saved
		resolution = 1,							# resolution scaling factor for plotting
		folder = NULL, 							# folder, where to save the plots
		condition = "", 						# to be added to the plotnames
		addformat = NULL	 					# additional fileformat for plots to be saved
)
{
	
	### PRELIMINARIES ###
	
	if (is.null(ccl)){alpha = 0} else {alpha = log(2)/ccl}
	ncells = function(t){n*exp(alpha*t)}
	sq = seq(0,duration,1)
	sparsesq = seq(0,duration,lab.duration)
	
	### THE "TRUE" SYNTHESIS RATES (mu) ###
	
	mu = function(t){return(mu.values[sum(mu.breaks <= t)+1])}
	mu.vals = sapply(sq,mu)
	truesynthesisrate = mu.vals[sparsesq[-1]+1]
	truesynthesisrateaveraged = c()
	for (i in 1:(duration/lab.duration)){truesynthesisrateaveraged[i] = mean(mu.vals[((i*lab.duration-lab.duration):(i*lab.duration))+1])}
	
	### THE "TRUE" DECAY RATES (lambda) ###
	
	lambda = function(t){return(lambda.values[sum(lambda.breaks <= t)+1])}
	lambda.vals = sapply(sq,lambda)
	truedecayrate = lambda.vals[sparsesq[-1]+1]
	truedecayrateaveraged = c()
	for (i in 1:(duration/lab.duration)){truedecayrateaveraged[i] = mean(lambda.vals[((i*lab.duration-lab.duration):(i*lab.duration))+1])}
	piecewise.integrated.lambda = function(t){vec = sapply(c(lambda.breaks,t),function(x){min(x,t)});vec = vec - c(0,vec)[-(length(vec)+1)];return(sum(vec*lambda.values))}
	piecewise.integrated.lambda.vals = sapply(sq,piecewise.integrated.lambda)
	piecewise.integrated.lambda.vals.list = list()
	ini = 0
	for (i in 1:(duration/lab.duration)){piecewise.integrated.lambda.vals.list[[i]] = piecewise.integrated.lambda.vals[((i*lab.duration-lab.duration):(i*lab.duration))+1] - ini;ini = piecewise.integrated.lambda.vals[i*lab.duration+1]}
	for (i in 1:(duration/lab.duration)){piecewise.integrated.lambda.vals.list[[i]] = exp(-piecewise.integrated.lambda.vals.list[[i]])}
	
	### PRELIMINARY FORMULAE ###
			
	steady.state = n*mu.vals[1]/(alpha + lambda.vals[1])
	hg.b = c(sort(unique(c(mu.breaks, lambda.breaks))),duration)
	hg.a = c(0,hg.b[-length(hg.b)])
	phi.const = (sapply(hg.a,mu)*exp(sapply(hg.a,piecewise.integrated.lambda))*(exp(alpha*hg.b)-exp(alpha*hg.a + sapply(hg.a,lambda)*(hg.b-hg.a))))/(alpha - sapply(hg.a,lambda))
	phi.sums = function(t){return(sum(c(0,phi.const)[1:(sum(hg.a <= t))]))}	
	phi.partial = function(t){a = max(hg.a[hg.a <= t]);return((sapply(a,mu)*exp(sapply(a,piecewise.integrated.lambda))*(exp(alpha*t)-exp(alpha*a + sapply(a,lambda)*(t-a))))/(alpha - sapply(a,lambda)))}
	phi.vals = sapply(sq,phi.partial) + sapply(sq,phi.sums)
	
	total.vals = exp(-piecewise.integrated.lambda.vals)*(steady.state + n*phi.vals)
	unlabeled.vals = exp(-piecewise.integrated.lambda.vals)*steady.state
			
	### THE "TRUE" AMOUNT OF TOTAL RNA (Cgrt) ###
	
	truetotal = total.vals[sparsesq+1]
	
	### THE "TRUE" AMOUNT OF PRE-EXISTING RNA (Bgrt) ###
	
	unlabeled.vals.list = list()
	for (i in 1:(duration/lab.duration)){
		unlabeled.vals.list[[i]] = piecewise.integrated.lambda.vals.list[[i]]*total.vals[(i*lab.duration-lab.duration)+1]
	}
	trueunlabeled = sapply(unlabeled.vals.list,c)[lab.duration+1,]
	
	### THE "TRUE" AMOUNT OF NEWLY SYNTHESISED RNA (Agrt) ###
	
	labeled.vals = total.vals - unlabeled.vals
	labeled.vals.list = list()
	for (i in 1:(duration/lab.duration)){
		labeled.vals.list[[i]] = total.vals[((i*lab.duration-lab.duration):(i*lab.duration))+1] - unlabeled.vals.list[[i]]
	}
	truelabeled = truetotal[-1] - trueunlabeled
	
	### PLOT TIME COURSE (if plots = TRUE) ###
	
	if (plots){
		plotsfkt = function(){
			scale = 1.5
			scex = 1.75*scale
			sc1 = 1.25/scale
			sc2 = 2.5/scale
			par(mfrow=c(4,1))
			par(mar=c(5,4,4,2) + 2.5)
			plot(sq,mu.vals,col=cols[1],pch=20,xlim=c(0,duration),ylim=c(0,max(mu.vals*(4/3))),xlab="Time [min]",ylab=expression(paste("Synthesis rate  ",mu[g])),main=expression(paste("Time course of synthesis rate  ",mu[g])),cex.main=1.25*scex,cex.lab=1*scex,cex.axis=1*scex,cex=sc1*scex)
			points(sparsesq[-1]-lab.duration/2,truesynthesisrateaveraged,cex=scex*sc2,col=cols[1],lwd=2)
			abline(v=sparsesq,lty=3,col="grey")
			plot(sq,lambda.vals,col=cols[2],pch=20,xlim=c(0,duration),ylim=c(0,max(lambda.vals*(4/3))),xlab="Time [min]",ylab=expression(paste("Decay rate  ",lambda[g])),main=expression(paste("Time course of decay rate  ",lambda[g])),cex.main=1.25*scex,cex.lab=1*scex,cex.axis=1*scex,cex=sc1*scex)
			points(sparsesq[-1]-lab.duration/2,truedecayrateaveraged,cex=scex*sc2,col=cols[2],lwd=2)
			abline(v=sparsesq,lty=3,col="grey")
			plot(sq,total.vals,pch=19,xlim = c(0,duration),ylim = c(0,max(total.vals)*(4/3)),col = cols[3],xlab="Time [min]",ylab=paste("Amount of RNA per",n,"cell(s)"),main=expression(paste("Time course of T, L and U (experiment duration)")),cex.main=1.25*scex,cex.lab=1*scex,cex.axis=1*scex,cex=sc1*scex)
			points(sq,unlabeled.vals,col = cols[2],cex=sc1*scex,pch=20)
			points(sq,labeled.vals,col = cols[1],cex=sc1*scex,pch=20)
			abline(v=sparsesq,lty=3,col="grey")
			legend("topright",c("total RNA","labeled RNA","unlabeled RNA"),col=c(cols[3],cols[1],cols[2]),pch=19,cex=2,bg="white",inset=0.01)
			plot(sq,total.vals,pch=19,xlim = c(0,duration),ylim = c(0,max(total.vals)*(4/3)),col = cols[3],xlab="Time [min]",ylab=paste("Amount of RNA per",n,"cell(s)"),main=expression(paste("Time course of T, L and U (labeling durations)")),cex.main=1.25*scex,cex.lab=1*scex,cex.axis=1*scex,cex=sc1*scex)
			for (i in 1:(duration/lab.duration)){
				points((i*lab.duration-lab.duration):(i*lab.duration),unlabeled.vals.list[[i]],col = cols[2],cex=sc1*scex,pch=20)
				points((i*lab.duration-lab.duration):(i*lab.duration),labeled.vals.list[[i]],col = cols[1],cex=sc1*scex,pch=20)
			}
			points(sparsesq[-1],truetotal[-1],cex=scex*sc2,col=cols[3],lwd=2)
			points(sparsesq[-1],trueunlabeled,cex=scex*sc2,col=cols[2],lwd=2)
			points(sparsesq[-1],truelabeled,cex=scex*sc2,col=cols[1],lwd=2)
			abline(v=sparsesq,lty=3,col="grey")
			legend("topright",c("total RNA","labeled RNA","unlabeled RNA"),col=c(cols[3],cols[1],cols[2]),pch=19,cex=2,bg="white",inset=0.01)
		}
		DTA.plot.it(filename = paste(folder,"/time_course_",condition,sep=""),sw = resolution*duration/25,sh = resolution*3,sres = resolution,plotsfkt = plotsfkt,ww = 7*duration/25,wh = 21,saveit = save.plots,addformat = addformat)		
	}
	
	### RETURN RESULTS IN A LIST ###
	
	results = list()
	results[["truetotal"]] = truetotal[-1]									# the true total RNA amount (C_gr in the paper)
	results[["truelabeled"]] = truelabeled									# the true newly synthesized RNA (A_gr in the paper)
	results[["trueunlabeled"]] = trueunlabeled								# the true amount of RNA synthesized before labeling (B_gr)
	results[["truesynthesisrate"]] = truesynthesisrate						# the true synthesis rates
	results[["truedecayrate"]] = truedecayrate								# the true decay rates
	results[["truesynthesisrateaveraged"]] = truesynthesisrateaveraged		# the true (time-averaged) synthesis rates
	results[["truedecayrateaveraged"]] = truedecayrateaveraged				# the true (time-averaged) decay rates
	return(results)
}


### DTA.dynamic.generate produces the according phenotype matrix and data matrix consisting of total, labeled and unlabeled measurements of a simulated time course. DATA IS RETURNED ON THE ABSOLUTE SCALE !!!! ###


DTA.dynamic.generate = function(
		duration = 60,						# duration of the whole time course (min)
		lab.duration = 6,					# labeling duration for single experiments (min)
		tnumber = NULL,   					# the number of uridine residues for each gene
		plabel = NULL,						# the efficiency with which a Uridine position is finally Biotynilated
		nrgenes = 5000,						# the number of genes the simulated experiment will have (will be cropped if it exceeds the length of tnumber)
		mediantime.halflives = 12,			# the median of the half life distribution
		mediantime.synthesisrates = 18,		# the median of the synthesis rates distribution (counts/cell/cellcycle)
		n = 1,   							# the number of cells N(0)
		ccl = NULL, 						# the cell cycle length of the cells
		check = TRUE, 						# if check=TRUE, control messages and plots will be generated
		plots = FALSE, 						# if plots = TRUE, control plots will be plotted
		save.plots = FALSE, 				# if save.plots = TRUE, control plots will be saved
		folder = NULL, 						# folder, where to save the plots
		condition = "", 					# to be added to the plotnames
		addformat = NULL,	 				# additional fileformat for plots to be saved
		sdnoise = 0.075,					# the amount of measurement noise (proportional to expression strength)
		nobias = FALSE,						# should a labeling bias be added
		unspecific.LtoU = 0,				# proportion of labeled RNAs that unspecifically end up in the unlabeled fraction
		unspec.LtoU.weighted = FALSE,		# do unspecific proportion of labeled to unlabeled depend linearly on the length of the RNA
		unspecific.UtoL = 0,				# proportion of unlabeled RNAs that unspecifically end up in the labeled fraction
		unspec.UtoL.weighted = FALSE,		# do unspecific proportion of unlabeled to labeled depend linearly on the length of the RNA
		mu.values.mat = NULL,				# if the data should be generated using given synthesis rates, this matrix must contain the respective values for each gene 
		mu.breaks.mat = NULL,				# timepoints of synthesis rate changes, this matrix must contain the respective values for each gene, only needed when mu.values.mat is given (one column less than mu.values.mat)
		lambda.values.mat = NULL,			# if the data should be generated using given decay rates, this matrix must contain the respective values for each gene 
		lambda.breaks.mat = NULL,			# timepoints of decay rate changes, this matrix must contain the respective values for each gene, only needed when lambda.values.mat is given (one column less than lambda.values.mat)
		truehalflives = NULL,				# if the data should be generated using a given half life distribution, this vector must contain the respective values for each gene
		truesynthesisrates = NULL,			# if the data should be generated using a given synthesis rates distribution, this vector must contain the respective values for each gene
		genenames = NULL					# an optional list of gene names
)
{
	### TIMEPOINTS ###
	
	nrtimepoints = duration/lab.duration
	if (as.integer(nrtimepoints) != nrtimepoints){stop("The duration of the whole time course (duration) should be a multiple of the labeling duration for single experiments (lab.duration) !")}
	timepoints = rep(lab.duration,nrtimepoints)
	nrexperiments = length(timepoints)
	alpha = log(2)/ccl
	
	### PLABEL ###
	
	if (is.null(plabel)){plabel = rep(0.005,length(timepoints))}
	
	### GENENAMES ###
	
	if (!is.null(genenames)){
		if (length(genenames) < nrgenes){
			stop("The number of genes exceeds the length of genenames. They will not be used !")
			genenames = paste("ID",1:nrgenes,sep="-")
		}
	} else {genenames = paste("ID",1:nrgenes,sep="-")}
	
	### TNUMBER ###
	
	if (!is.null(tnumber)){
		if (length(tnumber) < nrgenes){
			stop("The number of genes exceeds the length of tnumber. They will not be used !")
			tnumber = round(rf(nrgenes,5,10)*350)
		} else {tnumber = sample(tnumber,nrgenes)}
	} else {tnumber = round(rf(nrgenes,5,10)*350)}
	names(tnumber) = genenames
	
	### HALFLIVES AND DECAYRATES ###
	
	if (!is.null(truehalflives)){
		if (length(truehalflives) < nrgenes){
			stop("The number of genes exceeds the length of truehalflives. They will not be used !")
			truehalflives = rf(nrgenes,15,15)*mediantime.halflives
		}
	} else {truehalflives = rf(nrgenes,15,15)*mediantime.halflives}
	names(truehalflives) = genenames
	truelambdas = log(2)/truehalflives
	
	### SYNTHESIS RATES ###
	
	if (!is.null(truesynthesisrates)){
		if (length(truesynthesisrates) < nrgenes){
			stop("The number of genes exceeds the length of truehalflives. They will not be used !")
			truesynthesisrates = rf(nrgenes,5,5)*mediantime.synthesisrates
		}
	} else {truesynthesisrates = rf(nrgenes,5,5)*mediantime.synthesisrates}
	names(truesynthesisrates) = genenames
	
	### CREATE PHENOMAT ###
	
	phenomat = DTA.phenomat(timepoints,timecourse = 1:nrtimepoints)
	
	### CREATE MATRICES ###
	
	valsmat = matrix(0,ncol=nrexperiments,nrow=nrgenes)
	colnames(valsmat) = phenomat[which(phenomat[,"fraction"] == "T"),"name"]
	rownames(valsmat) = genenames
	Cgrtmat = valsmat
	Agrtmat = valsmat
	Bgrtmat = valsmat
	truesynthesisratemat = valsmat
	truedecayratemat = valsmat
	truesynthesisrateaveragedmat = valsmat
	truedecayrateaveragedmat = valsmat
	
	### MU AND LAMBDA BREAKS AND CHANGES ###
	
	if (is.null(mu.values.mat)){mu.values.mat = cbind(truesynthesisrates)}
	if (is.null(lambda.values.mat)){lambda.values.mat = cbind(truelambdas)}
	
	### FILL MATRICES ###
	
	for (i in 1:nrgenes){
		res = DTA.dynamic.singlegenerate(mu.values = mu.values.mat[i,],mu.breaks = mu.breaks.mat[i,],lambda.values = lambda.values.mat[i,],lambda.breaks = lambda.breaks.mat[i,],
				n = n,ccl = ccl,duration = duration,lab.duration = lab.duration)
		Cgrtmat[i,] = res$'truetotal'
		Agrtmat[i,] = res$'truelabeled'
		Bgrtmat[i,] = res$'trueunlabeled'
		truesynthesisratemat[i,] = res$'truesynthesisrate'
		truedecayratemat[i,] = res$'truedecayrate'
		truesynthesisrateaveragedmat[i,] = res$'truesynthesisrateaveraged'
		truedecayrateaveragedmat[i,] = res$'truedecayrateaveraged'
	}
	
	### PRELIMINARIES ###
	
	truear = numeric(nrexperiments)
	truebr = numeric(nrexperiments)
	truecr = numeric(nrexperiments)
	truecrbyar = numeric(nrexperiments)
	truecrbybr = numeric(nrexperiments)
	truebrbyar = numeric(nrexperiments)
	trueLasymptote = numeric(nrexperiments)
	trueUasymptote = numeric(nrexperiments)
	
	### LABELED RNA THAT UNSPECIFICALLY ENDS UP IN UNLABELED ### 
	
	if (unspec.LtoU.weighted){weights.LtoU = tnumber/max(tnumber) - median(tnumber/max(tnumber)) + 1} else {weights.LtoU = 1}
	unspecL = Agrtmat * unspecific.LtoU * weights.LtoU
	
	### UNLABELED RNA THAT UNSPECIFICALLY ENDS UP IN LABELED ### 
	
	if (unspec.UtoL.weighted){weights.UtoL = tnumber/max(tnumber) - median(tnumber/max(tnumber)) + 1} else {weights.UtoL = 1}
	unspecU = Bgrtmat * unspecific.UtoL * weights.UtoL
	
	### THE AMOUNT OF MEASURED TOTAL RNA (Tg) ###
	
	Tgmat = apply(Cgrtmat,2,function(x){rnorm(length(x),mean=x,sd=x*sdnoise)})
	rownames(Tgmat) = genenames
	
	truecr = apply(Tgmat,2,median)/apply(Cgrtmat,2,median)
	
	### THE AMOUNT OF MEASURED LABELED RNA (Lg) ###
	
	if (nobias){labeleff = apply(as.matrix(rep(1,length(plabel))),1,function(x){bias(x,tnumber)})} else {labeleff =  apply(as.matrix(plabel),1,function(x){bias(x,tnumber)})}
	Lbetween = Agrtmat * labeleff + unspecU
	
	Lgmat = apply(Lbetween,2,function(x){rnorm(length(x),mean=x,sd=x*sdnoise)})
	Lgmat = t(t(Lgmat)/apply(Lgmat,2,median)*apply(Tgmat,2,median))
	rownames(Lgmat) = genenames
	
	truear = apply(Lgmat,2,median)/apply(Lbetween,2,median)
	
	### THE AMOUNT OF MEASURED UNLABELED RNA (Ug) ###
	
	Utrue = Cgrtmat + unspecL - Lbetween
	Ugmat = apply(Utrue,2,function(x){rnorm(length(x),mean=x,sd=x*sdnoise)})
	Ugmat = t(t(Ugmat)/apply(Ugmat,2,median)*apply(Tgmat,2,median))
	rownames(Ugmat) = genenames
	
	truebr = apply(Ugmat,2,median)/apply(Utrue,2,median)
	
	### THE TRUE PROPORTIONALITY CONSTANTS ###
	
	truecrbyar = truecr/truear
	truecrbybr = truecr/truebr
	truebrbyar = apply(Agrtmat,2,median)/apply(Bgrtmat,2,median)
	trueLasymptote = log(truecrbyar*apply(Agrtmat/Cgrtmat,2,median))
	trueUasymptote = log(truecrbybr*apply(Bgrtmat/Cgrtmat,2,median))
	
	datamat = cbind(Tgmat,Lgmat,Ugmat)
	colnames(datamat) = rownames(phenomat)
	
	### RETURN RESULTS IN A LIST ###
	
	results = list()
	results[["phenomat"]] = phenomat
	results[["datamat"]] = datamat
	results[["tnumber"]] = tnumber
	results[["ccl"]] = ccl	
	results[["truemus"]] = truesynthesisratemat
	results[["truemusaveraged"]] = truesynthesisrateaveragedmat
	results[["truelambdas"]] = truedecayratemat
	results[["truelambdasaveraged"]] = truedecayrateaveragedmat
	results[["truehalflives"]] = log(2)/truedecayratemat
	results[["truehalflivesaveraged"]] = log(2)/truedecayrateaveragedmat
	results[["trueplabel"]] = plabel
	results[["truear"]] = truear
	results[["truebr"]] = truebr
	results[["truecr"]] = truecr
	results[["truecrbyar"]] = truecrbyar
	results[["truecrbybr"]] = truecrbybr
	results[["truebrbyar"]] = truebrbyar
	results[["trueLasymptote"]] = trueLasymptote
	results[["trueUasymptote"]] = trueUasymptote
	return(results)
}



