

### singleestimate uses an experiment, given by a phenotype matrix, data matrix and the #uridines for each gene to estimate synthesis and decay rate of the genes ###


DTA.singleestimate = function(phenomat, 	# phenotype matrix, "nr" should be numbered by experiments not by replicates
		datamat, 							# data matrix, should only contain the rows of phenomat as columns
		tnumber, 							# #uridines, should have the rownames of datamat
		labelingtime,						# length of the labeling period
		ccl = NULL, 						# the cell cycle length of the cells
		mRNAs = NULL,						# estimated number of mRNAs in a cell
		reliable = NULL, 					# vector of reliable genes, which are used for regression
		mediancenter = TRUE, 				# should the L/T resp. U/T ratio of replicates be rescaled to a common median before rate extraction
		ratiomethod = "tls", 				# choose the regression method to be used, possible methods are: tls, bias, lm
		largest = 5, 						# percentage of largest residues from the first regression not to be used in the second regression
		weighted = TRUE, 					# should the regression be weighted with 1/(T^2 + median(T))
		usefractions = "LandT",				# from which fractions should the decay rate be calculated: "LandT", "UandT" or "both"
		ratio = NULL, 						# coefficient to rescale the fractions
		relevant = NULL, 					# choose the arrays to be used for halflives calculation, vector due to experiments variable 
		check = TRUE, 						# if check = TRUE, control messages and plots will be generated
		error = TRUE,						# should the measurement error be assessed
		samplesize = 1000,					# error model samplesize for resampling
		confidence.range = c(0.025,0.975),	# confidence region for error model
		bicor = TRUE, 						# should the labeling bias be corrected
		condition = "", 					# to be added to the plotnames
		timepoint = "",						# to be added to the plotnames
		upper = 700, 						# upper bound for labeling bias estimation
		lower = 500, 						# lower bound for labeling bias estimation
		save.plots = FALSE, 				# if save.plots = TRUE, control plots will be saved
		resolution = 1,						# resolution scaling factor for plotting
		notinR = FALSE,						# should plot be not plotted in R
		RStudio = FALSE,					# for RStudio users
		folder = NULL, 						# folder, where to save the plots
		fileformat = "jpeg", 				# save the plot as jpeg, png, bmp, tiff, ps or pdf
		totaloverwt = 1, 					# total mRNA over WT
		simulation = FALSE,					# simulated data via sim.object ?
		dynamic = FALSE,					# should be TRUE for timecourse data
		initials = NULL,					# initial values of total for timecourse data
		truemus = NULL,						# the true synthesis rates
		truemusaveraged = NULL,				# the true synthesis rates (averaged over labeling period)
		truelambdas = NULL,					# the true decay rates
		truelambdasaveraged = NULL,			# the true decay rates (averaged over labeling period)
		truehalflives = NULL,				# the true half-lives
		truehalflivesaveraged = NULL,		# the true half-lives (averaged over labeling period)
		trueplabel = NULL,					# the true labeling efficiency
		trueLasymptote = NULL,				# the true Lasymptote
		trueUasymptote = NULL,				# the true Uasymptote
		truecrbyar = NULL,					# the true quotient truecr/truear
		truecrbybr = NULL,					# the true quotient truecr/truebr
		truebrbyar = NULL					# the true quotient truebr/truear
)
{
	
	### PRELIMINARIES ###
	
	datamat = datamat[,rownames(phenomat)]
	datamatreliable = datamat[reliable,]
	tnumber = tnumber[rownames(datamat)]
	tnumberreliable = tnumber[reliable]
	if (is.null(ccl)){alpha = 0} else {alpha = log(2)/ccl}
	if (!is.null(initials)){alpha = 0}
	if (is.null(initials)){initialdummy = NULL} else {initialdummy = "notNULL"}
	if (is.null(ratio)){ratiodummy = NULL} else {ratiodummy = "notNULL"}
	if (!is.null(ratio)){ratiomethod = "tls"}
	experiments = unique(phenomat[,"nr"])
	nrexperiments = length(experiments)
	nrgenes = nrow(datamat)
	plabel = numeric(nrexperiments)
	labelratio = numeric(nrexperiments)
	lambda  = numeric(nrgenes)
	logquotient = numeric(nrexperiments)
	crbyarestimate = numeric(nrexperiments)
	crbybrestimate = numeric(nrexperiments)
	brbyarestimate = numeric(nrexperiments)
	planes = list()
	discards = list()
	results = list()
	calcdatamat = datamat
	calcdatamatreliable = calcdatamat[reliable,]
	
	### FIND TRIPLES ###
	
	triples = matrix(0,nrow=0,ncol=3)
	colnames(triples) = c("T","U","L")
	for (expnr in seq(experiments)){
		isT = which((phenomat[,"time"]==labelingtime) & (phenomat[,"nr"]==experiments[expnr]) & (phenomat[,"fraction"]=="T"))
		if (length(isT)==0) isU = NA
		isU = which((phenomat[,"time"]==labelingtime) & (phenomat[,"nr"]==experiments[expnr]) & (phenomat[,"fraction"]=="U"))
		if (length(isU)==0) isU = NA
		isL = which((phenomat[,"time"]==labelingtime) & (phenomat[,"nr"]==experiments[expnr]) & (phenomat[,"fraction"]=="L"))
		if (length(isL)==0) isL = NA
		triples = rbind(triples,c(isT,isU,isL))
	}
	results[["triples"]] = triples
	
	### UNLABELED CHECK ###
	
	if (any(is.na(triples[,"U"]))) {
		unlabeledfraction = FALSE
		if (is.null(ratiodummy)){stop("You have to specify the ratio of L to T (LtoTratio), as your data lacks the unlabeled fraction !")}
		ratiomethod = "tls"
		if (any(usefractions %in% c("UandT","both"))){print("usefractions is set to 'LandT' !")}
		usefractions = "LandT"
	} else {unlabeledfraction = TRUE}
	
	### ESTIMATION OF plabel L/T ###
	
	for (expnr in seq(experiments)){Tnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "T"))
		Tname = phenomat[Tnr,"name"]
		Lnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "L"))
		Lname = phenomat[Lnr,"name"]
		thistime = as.numeric(phenomat[Tnr,"time"])
		largeT = which(tnumberreliable>upper)
		asymptote = log(median(calcdatamatreliable[largeT,Lname]/calcdatamatreliable[largeT,Tname]))
		smallT = which(tnumberreliable<lower)
		
		lossfct = function(q){hilf = abs(log(calcdatamatreliable[smallT,Lname]/calcdatamatreliable[smallT,Tname]) - asymptote - log(bias(q,tnumberreliable[smallT])))
			result = median( hilf[is.double(hilf)] )
			return(result)
		}
		
		plabel[expnr] = optimize(lossfct,interval=c(0,0.05))$minimum
	}
	for (expnr in seq(experiments)){
		if (!simulation){if (check) cat("Experiment no.",expnr,", plabel: ",round(plabel[expnr],digits=4),"\n")}
		if (simulation){if (check) cat("Experiment no.",expnr,", plabel: ",round(plabel[expnr],digits=4),"(",round(trueplabel[expnr],digits=4),")","\n")}
	}
	results[["plabel"]] = plabel
	
	### ESTIMATION OF labelratio U/T ###
	
	if (unlabeledfraction){
		for (expnr in seq(experiments)){
			Tnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "T"))
			Tname = phenomat[Tnr,"name"]
			T = calcdatamat[reliable,Tname]
			Unr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "U"))
			Uname = phenomat[Unr,"name"]
			U = calcdatamat[reliable,Uname]
			
			largeT = which(tnumberreliable>upper)
			asymptote = log(median(calcdatamatreliable[largeT,Uname]/calcdatamatreliable[largeT,Tname]))
			smallT = which(tnumberreliable<lower)
			
			lossfct = function(q){hilf = abs(log(calcdatamatreliable[smallT,Uname]/calcdatamatreliable[smallT,Tname]) - asymptote - log(1 + q*biasdev(plabel[expnr],tnumberreliable[smallT])))
				result = median(hilf[is.double(hilf)])
				return(result)
			}
			
			labelratio[expnr] = optimize(lossfct,interval=c(0,5))$minimum
		}
		for (expnr in seq(experiments)){
			if (!simulation){if (check) cat("Experiment no.",expnr,", labelratio: ",round(labelratio[expnr],digits=4),"\n")}
		}
	}
	
	### ESTIMATION OF THE QUOTIENT br/ar based on plabel ###
	
	if (is.null(ratiodummy)){		
		if (ratiomethod == "bias"){
			for (expnr in seq(experiments)){
				Tnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "T"))
				Tname = phenomat[Tnr,"name"]
				T = calcdatamat[reliable,Tname]
				Unr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "U"))
				Uname = phenomat[Unr,"name"]
				U = calcdatamat[reliable,Uname]
				
				largeT = which(tnumberreliable>upper)
				asymptote = log(median(calcdatamatreliable[largeT,Uname]/calcdatamatreliable[largeT,Tname]))
				smallT = which(tnumberreliable<lower)
				
				lossfct = function(q){hilf = abs(log(calcdatamatreliable[smallT,Uname]/calcdatamatreliable[smallT,Tname]) - asymptote - log(1 + q*biasdev(plabel[expnr],tnumberreliable[smallT])))
					result = median(hilf[is.double(hilf)])
					return(result)
				}
				
				brbyarestimate[expnr] = optimize(lossfct,interval=c(0,5))$minimum
			}
			
			ratio = brbyarestimate
			for (expnr in seq(experiments)){
				if (!simulation){if (check) cat("Experiment no.",expnr,", L to U ratio: ",round(brbyarestimate[expnr],digits=2),"\n")}
				if (simulation){if (check) cat("Experiment no.",expnr,", L to U ratio: ",round(brbyarestimate[expnr],digits=2),"(",round(truebrbyar[expnr],digits=2),")","\n")}
			}
		}
		if (any(ratio <= 0)){stop("You have to specify the ratio of L to T (LtoTratio), as your data is not suitable for 'ratiomethod = bias' !")}
	}	
	
	### ESTIMATION OF THE QUOTIENTS ar/cr AND br/cr ###
	
	if (is.null(ratiodummy)){
		if (ratiomethod %in% c("tls","lm")){
			for (expnr in seq(experiments)){
				Tnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "T"))
				Tname = phenomat[Tnr,"name"]
				T = calcdatamat[reliable,Tname]
				Lnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "L"))
				Lname = phenomat[Lnr,"name"]
				L = calcdatamat[reliable,Lname]
				Unr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "U"))
				Uname = phenomat[Unr,"name"]
				U = calcdatamat[reliable,Uname]
				
				if (weighted){weights = 1/(T^2 + median(T))} else{weights = NULL}
				
				if (ratiomethod == "lm"){linmod = lm(T ~ L + U + 0,weights=weights)
					crbyarestimate[expnr] = coefficients(linmod)["L"] # /(coefficients(linmod)["L"] + coefficients(linmod)["U"])
					crbybrestimate[expnr] = coefficients(linmod)["U"] # /(coefficients(linmod)["L"] + coefficients(linmod)["U"])
				}
				
				if (ratiomethod == "tls"){linmod = tls(T ~ L + U + 0,D=diag(weights,nrow=length(weights),ncol=length(weights)))
					crbyarestimate[expnr] = coefficients(linmod)["L"] # /(coefficients(linmod)["L"] + coefficients(linmod)["U"])
					crbybrestimate[expnr] = coefficients(linmod)["U"] # /(coefficients(linmod)["L"] + coefficients(linmod)["U"])
				}
				
				TLU = cbind(T,L,U)
				normalenvec = vnorm(c(-1,crbyarestimate[expnr],crbybrestimate[expnr]))
				Lproj = c(0,1,0)-(c(0,1,0)%*%normalenvec)*normalenvec
				Uproj = c(0,0,1)-(c(0,0,1)%*%normalenvec)*normalenvec
				basis = cbind(vnorm(Lproj),vnorm(Uproj),normalenvec)
				invbasis = solve(basis)
				plane = invbasis%*%t(TLU)
				absolutes = abs(t(plane)[,3])
				nearest = paste(100-largest,"%",sep="")
				quan = quantile(absolutes,probs = seq(0, 1, 0.01))[nearest]
				discards[[expnr]] = names(absolutes[absolutes > quan])
				forreg = names(absolutes[absolutes <= quan])
				T = T[forreg]
				U = U[forreg]
				L = L[forreg]
				
				if (weighted){weights = 1/(T[forreg]^2 + median(T[forreg]))} else{weights = NULL}
				
				if (ratiomethod == "lm"){linmod = lm(T ~ L + U + 0,weights=weights)
					crbyarestimate[expnr] = coefficients(linmod)["L"] # /(coefficients(linmod)["L"] + coefficients(linmod)["U"])
					crbybrestimate[expnr] = coefficients(linmod)["U"] # /(coefficients(linmod)["L"] + coefficients(linmod)["U"])
				}
				
				if (ratiomethod == "tls"){linmod = tls(T ~ L + U + 0,D=diag(weights,nrow=length(weights),ncol=length(weights)))
					crbyarestimate[expnr] = coefficients(linmod)["L"] # /(coefficients(linmod)["L"] + coefficients(linmod)["U"])
					crbybrestimate[expnr] = coefficients(linmod)["U"] # /(coefficients(linmod)["L"] + coefficients(linmod)["U"])
				}
				
				normalenvec = vnorm(c(-1,crbyarestimate[expnr],crbybrestimate[expnr]))
				Lproj = c(0,1,0)-(c(0,1,0)%*%normalenvec)*normalenvec
				Uproj = c(0,0,1)-(c(0,0,1)%*%normalenvec)*normalenvec
				basis = cbind(vnorm(Lproj),vnorm(Uproj),normalenvec)
				invbasis = solve(basis)
				plane = invbasis%*%t(TLU)
				planes[[expnr]] = plane
			}
			
			### PLOTS OF REGRESSION ###
			
			if (check){
				plotsfkt = function(){
					par(mfrow=windowxy(nrexperiments))
					parfkt("default",nrexperiments)
					for (expnr in seq(experiments)){
						Tnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "T"))
						ylim = c(-max(abs(t(planes[[expnr]])[,3])),max(abs(t(planes[[expnr]])[,3])))
						heatscatter(log(t(planes[[expnr]])[,1]),t(planes[[expnr]])[,3],cor=FALSE,ylim=ylim,xlab="",ylab="",main="")
						points(log(t(planes[[expnr]])[,1])[discards[[expnr]]],t(planes[[expnr]])[,3][discards[[expnr]]],col = "red")
						xlab = expression(paste("log( Orthogonal projection on ",L[gr]," )"))
						ylab = expression("Normal of the plane")
						main = expression(paste("Regression plane:  ratio of  ",L[gr]/T[gr]))
						sub = paste("(",phenomat[Tnr,"name"],"  c = ",signif(crbyarestimate[expnr],digits=2),")")
						mtextfkt("default",nrexperiments,main,xlab,ylab,sub)
						abline(h = 0,col = "green",lwd=2)
						abline(h = max(t(planes[[expnr]])[,3]),col = "red",lwd=1)
						abline(h = min(t(planes[[expnr]])[,3]),col = "red",lwd=1)
					}
				}
				DTA.plot.it(filename = paste(folder,"/regression_hyperline_",condition,"_",timepoint,sep=""),sw = resolution*windowxy(nrexperiments)[2],sh = resolution*windowxy(nrexperiments)[1],sres = resolution,plotsfkt = plotsfkt,ww = 7*windowxy(nrexperiments)[2],wh = 7*windowxy(nrexperiments)[1],saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)
				
				plotsfkt = function(){
					par(mfrow=windowxy(nrexperiments))
					parfkt("default",nrexperiments)
					for (expnr in seq(experiments)){
						Tnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "T"))
						zlim = c(-max(abs(t(planes[[expnr]])[,3])),max(abs(t(planes[[expnr]])[,3])))
						s3d <- scatterplot3d(t(planes[[expnr]])[,1],t(planes[[expnr]])[,2],t(planes[[expnr]])[,3],pch=20,xlab="",ylab="",zlab="",scale.y=1,angle=40,
								mar = parfkt("scatterplot3d",nrexperiments),highlight.3d=TRUE,main="",grid=TRUE,zlim=zlim)
						
						xlab = expression(paste("Orthogonal projection on  ",L[gr]))
						ylab = expression("Normal of the plane")
						zlab = expression(paste("Orthogonal projection on  ",U[gr]))
						main = expression(paste("Regression plane 3D:  ratio of  ",L[gr]/T[gr]))
						sub = paste("(",phenomat[Tnr,"name"],"  c = ",signif(crbyarestimate[expnr],digits=2),")")
						mtextfkt("scatterplot3d",nrexperiments,main,xlab,ylab,sub,zlab)
						s3d$points3d(t(planes[[expnr]])[,1][discards[[expnr]]],t(planes[[expnr]])[,2][discards[[expnr]]],t(planes[[expnr]])[,3][discards[[expnr]]],pch=20,col = "red")
						s3d$plane3d(c(0,0,0), lty = "solid", lty.box = "solid",col = "green",lwd=2)
						s3d$plane3d(c(max(t(planes[[expnr]])[,3]),0,0), lty = "solid", lty.box = "solid",col = "red",lwd=1)
						s3d$plane3d(c(min(t(planes[[expnr]])[,3]),0,0), lty = "solid", lty.box = "solid",col = "red",lwd=1)
					}
				}
				DTA.plot.it(filename = paste(folder,"/regression_hyperplane_",condition,"_",timepoint,sep=""),sw = resolution*windowxy(nrexperiments)[2],sh = resolution*windowxy(nrexperiments)[1],sres = resolution,plotsfkt = plotsfkt,ww = 7*windowxy(nrexperiments)[2],wh = 7*windowxy(nrexperiments)[1],saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)		
			}
			
			if (ratiomethod == "lm"){ratio = crbyarestimate
				for (expnr in seq(experiments)){
					if (!simulation){if (check) cat("Experiment no.",expnr,", L to T ratio: ",round(crbyarestimate[expnr],digits=2),"\n")}
					if (simulation){if (check) cat("Experiment no.",expnr,", L to T ratio: ",round(crbyarestimate[expnr],digits=2),"(",round(truecrbyar[expnr],digits=2),")","\n")}
				}
			}	
			if (ratiomethod == "tls"){ratio = crbyarestimate
				for (expnr in seq(experiments)){
					if (!simulation){if (check) cat("Experiment no.",expnr,", L to T ratio: ",round(crbyarestimate[expnr],digits=2),"\n")}
					if (simulation){if (check) cat("Experiment no.",expnr,", L to T ratio: ",round(crbyarestimate[expnr],digits=2),"(",round(truecrbyar[expnr],digits=2),")","\n")}
				}
			}
		}	
		if (any(ratio <= 0)){stop(paste("You have to specify the ratio of L to T (LtoTratio), as your data is not suitable for 'ratiomethod = ",ratiomethod,"' !",sep=""))}
	}
	
	### MERGE METHODS FOR cr/ar, cr/br AND br/ar ###
	
	if (is.null(ratiodummy)){
		if (ratiomethod == "bias"){
			crbyarestimate = ratio/(1+ratio)
			crbybrestimate = 1-(ratio/(1+ratio))
			brbyarestimate = ratio
			results[["LtoTratio"]] = ratio/(1+ratio)
			results[["UtoTratio"]] = 1-ratio/(1+ratio)
			results[["LtoUratio"]] = brbyarestimate
		}
		else if (ratiomethod %in% c("tls","lm")){
			crbyarestimate = ratio
			crbybrestimate = 1-ratio
			brbyarestimate = ratio/(1-ratio)
			results[["LtoTratio"]] = ratio
			results[["UtoTratio"]] = 1-ratio
			results[["LtoUratio"]] = ratio/(1-ratio)
		}
	} else {
		crbyarestimate = ratio
		crbybrestimate = 1-ratio
		brbyarestimate = ratio/(1-ratio)
		results[["LtoTratio"]] = ratio
		results[["UtoTratio"]] = 1-ratio
		results[["LtoUratio"]] = ratio/(1-ratio)
	}
	
	### PLOTS OF LABELING BIAS ###
	
	if (check){
		plotsfkt = function(){
			par(mfrow=windowxy(nrexperiments))
			parfkt("default",nrexperiments)
			for (expnr in seq(experiments)){
				Tnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "T"))
				Tname = phenomat[Tnr,"name"]
				Lnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "L"))
				Lname = phenomat[Lnr,"name"]
				thistime = as.numeric(phenomat[Tnr,"time"])
				largeT = which(tnumberreliable>upper)
				asymptote = log(median( calcdatamatreliable[largeT,Lname]/calcdatamatreliable[largeT,Tname] ))
				smallT = which(tnumberreliable<lower)
				main = expression(paste("Labeling Bias ",l[gr]))
				xlab = expression(paste("Number of Uridine residues ",'#u'[g]))
				ylab = expression(paste("log( ",L[gr]/T[gr]," )"))
				sub = paste("( ",Tname,"  p = ",signif(plabel[expnr],2),"  asymptote = ",signif(asymptote,2)," )",sep="")
				tseq = min(tnumberreliable,na.rm=TRUE):max(tnumberreliable,na.rm=TRUE)
				if (!bicor){main = expression(paste("(not corrected) Labeling Bias ",l[gr]))}
				heatscatter(tnumberreliable,log(calcdatamatreliable[,Lname]/calcdatamatreliable[,Tname]),xlab="",ylab="",main="",xlim=c(0,2000),ylim=c(-4,4),cor=FALSE)
				mtextfkt("default",nrexperiments,main,xlab,ylab,sub)
				if (bicor){points(tseq,log(bias(plabel[expnr],tseq))+asymptote,type="l",col="black",lwd=3)} else {abline(h = asymptote,lwd=3)}
				if (!is.null(trueplabel[experiments[expnr]]) & !is.null(trueLasymptote[experiments[expnr]])){points(tseq,log(bias(trueplabel[experiments[expnr]],tseq))+asymptote,type="l",col="grey",lwd=3,lty=2)}
			}
		}
		DTA.plot.it(filename = paste(folder,"/estimation_bias_L_",condition,"_",timepoint,sep=""),sw = resolution*windowxy(nrexperiments)[2],sh = resolution*windowxy(nrexperiments)[1],sres = resolution,plotsfkt = plotsfkt,ww = 7*windowxy(nrexperiments)[2],wh = 7*windowxy(nrexperiments)[1],saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)    
		if (unlabeledfraction){
			if (is.null(ratiodummy) & ratiomethod == "bias"){
				plotsfkt = function(){
					par(mfrow=windowxy(nrexperiments))
					parfkt("default",nrexperiments)
					for (expnr in seq(experiments)){
						Tnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "T"))
						Tname = phenomat[Tnr,"name"]
						Unr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "U"))
						Uname = phenomat[Unr,"name"]
						thistime = as.numeric(phenomat[Tnr,"time"])
						largeT = which(tnumberreliable>upper)
						asymptote = log(median( calcdatamatreliable[largeT,Uname]/calcdatamatreliable[largeT,Tname] ))
						smallT = which(tnumberreliable<lower)
						main = expression(paste("Labeling Bias ",u[gr]))
						xlab = expression(paste("Number of Uridine residues ",'#u'[g]))
						ylab = expression(paste("log( ",U[gr]/T[gr]," )"))
						sub = paste("( ",Tname,"  p = ",signif(plabel[expnr],2),"  asymptote = ",signif(asymptote,2)," )",sep="")
						tseq = min(tnumberreliable,na.rm=TRUE):max(tnumberreliable,na.rm=TRUE)
						if (!bicor){main = expression(paste("(not corrected) Labeling Bias ",u[gr]))}
						heatscatter(tnumberreliable,log(calcdatamatreliable[,Uname]/calcdatamatreliable[,Tname]),xlab="",ylab="",main="",xlim=c(0,2000),ylim=c(-4,4),cor=FALSE)
						mtextfkt("default",nrexperiments,main,xlab,ylab,sub)
						if (bicor){points(tseq,log(1 + brbyarestimate[expnr]*biasdev(plabel[expnr],tseq))+asymptote,type="l",col="black",lwd=3)} else {abline(h = asymptote,lwd=3)}
						if (!is.null(trueplabel[experiments[expnr]]) & !is.null(trueUasymptote[experiments[expnr]])){points(tseq,log(1 + brbyarestimate[expnr]*biasdev(trueplabel[experiments[expnr]],tseq))+asymptote,type="l",col="grey",lwd=3,lty=2)}
					}
				}
				DTA.plot.it(filename = paste(folder,"/estimation_bias_U_",condition,"_",timepoint,sep=""),sw = resolution*windowxy(nrexperiments)[2],sh = resolution*windowxy(nrexperiments)[1],sres = resolution,plotsfkt = plotsfkt,ww = 7*windowxy(nrexperiments)[2],wh = 7*windowxy(nrexperiments)[1],saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)
			} else {
				plotsfkt = function(){
					par(mfrow=windowxy(nrexperiments))
					parfkt("default",nrexperiments)
					for (expnr in seq(experiments)){
						Tnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "T"))
						Tname = phenomat[Tnr,"name"]
						Unr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "U"))
						Uname = phenomat[Unr,"name"]
						thistime = as.numeric(phenomat[Tnr,"time"])
						largeT = which(tnumberreliable>upper)
						asymptote = log(median( calcdatamatreliable[largeT,Uname]/calcdatamatreliable[largeT,Tname] ))
						smallT = which(tnumberreliable<lower)
						main = expression(paste("Labeling Bias ",u[gr]))
						xlab = expression(paste("Number of Uridine residues ",'#u'[g]))
						ylab = expression(paste("log( ",U[gr]/T[gr]," )"))
						sub = paste("( ",Tname,"  p = ",signif(plabel[expnr],2),"  asymptote = ",signif(asymptote,2)," )",sep="")
						tseq = min(tnumberreliable,na.rm=TRUE):max(tnumberreliable,na.rm=TRUE)
						if (!bicor){main = expression(paste("(not corrected) Labeling Bias ",u[gr]))}
						heatscatter(tnumberreliable,log(calcdatamatreliable[,Uname]/calcdatamatreliable[,Tname]),xlab="",ylab="",main="",xlim=c(0,2000),ylim=c(-4,4),cor=FALSE)
						mtextfkt("default",nrexperiments,main,xlab,ylab,sub)
						if (bicor){points(tseq,log(1 + labelratio[expnr]*biasdev(plabel[expnr],tseq))+asymptote,type="l",col="black",lwd=3)} else {abline(h = asymptote,lwd=3)}
						if (!is.null(trueplabel[experiments[expnr]]) & !is.null(trueUasymptote[experiments[expnr]])){points(tseq,log(1 + labelratio[expnr]*biasdev(trueplabel[experiments[expnr]],tseq))+asymptote,type="l",col="grey",lwd=3,lty=2)}
					}
				}
				DTA.plot.it(filename = paste(folder,"/estimation_bias_U_",condition,"_",timepoint,sep=""),sw = resolution*windowxy(nrexperiments)[2],sh = resolution*windowxy(nrexperiments)[1],sres = resolution,plotsfkt = plotsfkt,ww = 7*windowxy(nrexperiments)[2],wh = 7*windowxy(nrexperiments)[1],saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)
			}
		}
	}
	
	### ESTIMATION OF lambda ###
	
	for (expnr in seq(experiments)){
		if (!bicor){
			plabel[expnr] = 1
			labelratio[expnr] = 1
		}
		Lnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "L"))
		Lname = phenomat[Lnr,"name"]
		calcdatamat[,Lname] = calcdatamat[,Lname]*(1/bias(plabel[expnr],tnumber))
		if (unlabeledfraction){
			Unr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "U"))
			Uname = phenomat[Unr,"name"]
			if (is.null(ratiodummy) & ratiomethod == "bias"){
				calcdatamat[,Uname] = calcdatamat[,Uname]*(1/(1 + (brbyarestimate[expnr])*biasdev(plabel[expnr],tnumber)))
			} else {
				calcdatamat[,Uname] = calcdatamat[,Uname]*(1/(1 + (labelratio[expnr])*biasdev(plabel[expnr],tnumber)))
			}
		}
	}
	
	if (is.null(ratiodummy) & ratiomethod == "bias"){calcdatamat = medctr(calcdatamat,userows = reliable,logscale = FALSE,protocol = FALSE)}
	calcdatamatreliable = calcdatamat[reliable,]
	
	### PLOTS OF LABELING BIAS CORRECTION ###
	
	if (check){
		plotsfkt = function(){
			par(mfrow=windowxy(nrexperiments))
			parfkt("default",nrexperiments)
			for (expnr in seq(experiments)){
				Tnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "T"))
				Tname = phenomat[Tnr,"name"]
				Lnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "L"))
				Lname = phenomat[Lnr,"name"]
				thistime = as.numeric(phenomat[Tnr,"time"])
				largeT = which(tnumberreliable>upper)
				asymptote = log(median( calcdatamatreliable[largeT,Lname]/calcdatamatreliable[largeT,Tname] ))
				smallT = which(tnumberreliable<lower)
				main = expression(paste("Labeling Bias ",l[gr]))
				xlab = expression(paste("Number of Uridine residues ",'#u'[g]))
				ylab = expression(paste("log( ",L[gr]/T[gr]," )"))
				sub = paste("( ",Tname,"  p = ",signif(1,2),"  asymptote = ",signif(asymptote,2)," )",sep="")
				tseq = min(tnumberreliable,na.rm=TRUE):max(tnumberreliable,na.rm=TRUE)
				main = expression(paste("(corrected) Labeling Bias ",l[gr]))
				if (!bicor){main = expression(paste("(not corrected) Labeling Bias ",l[gr]))}
				heatscatter(tnumberreliable,log(calcdatamatreliable[,Lname]/calcdatamatreliable[,Tname]),xlab="",ylab="",main="",xlim=c(0,2000),ylim=c(-4,4),cor=FALSE)
				mtextfkt("default",nrexperiments,main,xlab,ylab,sub)
				if (bicor){points(tseq,log(bias(1,tseq))+asymptote,type="l",col="black",lwd=3)} else {abline(h = asymptote,lwd=3)}
				if (!is.null(trueLasymptote[experiments[expnr]])){points(tseq,log(bias(1,tseq))+asymptote,type="l",col="grey",lwd=3,lty=2)}
			}
		}
		DTA.plot.it(filename = paste(folder,"/estimation_bias_L_",condition,"_",timepoint,"_corrected",sep=""),sw = resolution*windowxy(nrexperiments)[2],sh = resolution*windowxy(nrexperiments)[1],sres = resolution,plotsfkt = plotsfkt,ww = 7*windowxy(nrexperiments)[2],wh = 7*windowxy(nrexperiments)[1],saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)    
		if (unlabeledfraction){
			if (is.null(ratiodummy) & ratiomethod == "bias"){
				plotsfkt = function(){
					par(mfrow=windowxy(nrexperiments))
					parfkt("default",nrexperiments)
					for (expnr in seq(experiments)){
						Tnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "T"))
						Tname = phenomat[Tnr,"name"]
						Unr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "U"))
						Uname = phenomat[Unr,"name"]
						thistime = as.numeric(phenomat[Tnr,"time"])
						largeT = which(tnumberreliable>upper)
						asymptote = log(median( calcdatamatreliable[largeT,Uname]/calcdatamatreliable[largeT,Tname] ))
						smallT = which(tnumberreliable<lower)
						main = expression(paste("Labeling Bias ",u[gr]))
						xlab = expression(paste("Number of Uridine residues ",'#u'[g]))
						ylab = expression(paste("log( ",U[gr]/T[gr]," )"))
						sub = paste("( ",Tname,"  p = ",signif(1,2),"  asymptote = ",signif(asymptote,2)," )",sep="")
						tseq = min(tnumberreliable,na.rm=TRUE):max(tnumberreliable,na.rm=TRUE)
						main = expression(paste("(corrected) Labeling Bias ",u[gr]))
						if (!bicor){main = expression(paste("(not corrected) Labeling Bias ",u[gr]))}
						heatscatter(tnumberreliable,log(calcdatamatreliable[,Uname]/calcdatamatreliable[,Tname]),xlab="",ylab="",main="",xlim=c(0,2000),ylim=c(-4,4),cor=FALSE)
						mtextfkt("default",nrexperiments,main,xlab,ylab,sub)
						if (bicor){points(tseq,log(1 + 1*biasdev(1,tseq))+asymptote,type="l",col="black",lwd=3)} else {abline(h = asymptote,lwd=3)}
						if (!is.null(trueUasymptote[experiments[expnr]])){points(tseq,log(1 + 1*biasdev(1,tseq))+asymptote,type="l",col="grey",lwd=3,lty=2)}
					}
				}
				DTA.plot.it(filename = paste(folder,"/estimation_bias_U_",condition,"_",timepoint,"_corrected",sep=""),sw = resolution*windowxy(nrexperiments)[2],sh = resolution*windowxy(nrexperiments)[1],sres = resolution,plotsfkt = plotsfkt,ww = 7*windowxy(nrexperiments)[2],wh = 7*windowxy(nrexperiments)[1],saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)
			} else {
				plotsfkt = function(){
					par(mfrow=windowxy(nrexperiments))
					parfkt("default",nrexperiments)
					for (expnr in seq(experiments)){
						Tnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "T"))
						Tname = phenomat[Tnr,"name"]
						Unr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "U"))
						Uname = phenomat[Unr,"name"]
						thistime = as.numeric(phenomat[Tnr,"time"])
						largeT = which(tnumberreliable>upper)
						asymptote = log(median( calcdatamatreliable[largeT,Uname]/calcdatamatreliable[largeT,Tname] ))
						smallT = which(tnumberreliable<lower)
						main = expression(paste("Labeling Bias ",u[gr]))
						xlab = expression(paste("Number of Uridine residues ",'#u'[g]))
						ylab = expression(paste("log( ",U[gr]/T[gr]," )"))
						sub = paste("( ",Tname,"  p = ",signif(1,2),"  asymptote = ",signif(asymptote,2)," )",sep="")
						tseq = min(tnumberreliable,na.rm=TRUE):max(tnumberreliable,na.rm=TRUE)
						main = expression(paste("(corrected) Labeling Bias ",u[gr]))
						if (!bicor){main = expression(paste("(not corrected) Labeling Bias ",u[gr]))}
						heatscatter(tnumberreliable,log(calcdatamatreliable[,Uname]/calcdatamatreliable[,Tname]),xlab="",ylab="",main="",xlim=c(0,2000),ylim=c(-4,4),cor=FALSE)
						mtextfkt("default",nrexperiments,main,xlab,ylab,sub)
						if (bicor){points(tseq,log(1 + 1*biasdev(1,tseq))+asymptote,type="l",col="black",lwd=3)} else {abline(h = asymptote,lwd=3)}
						if (!is.null(trueUasymptote[experiments[expnr]])){points(tseq,log(1 + 1*biasdev(1,tseq))+asymptote,type="l",col="grey",lwd=3,lty=2)}
					}
				}
				DTA.plot.it(filename = paste(folder,"/estimation_bias_U_",condition,"_",timepoint,"_corrected",sep=""),sw = resolution*windowxy(nrexperiments)[2],sh = resolution*windowxy(nrexperiments)[1],sres = resolution,plotsfkt = plotsfkt,ww = 7*windowxy(nrexperiments)[2],wh = 7*windowxy(nrexperiments)[1],saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)
			}
		}
	}
	
	for (expnr in seq(experiments)){
		Lnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "L"))
		Lname = phenomat[Lnr,"name"]
		calcdatamat[,Lname] = calcdatamat[,Lname]*(crbyarestimate[expnr])
		if (unlabeledfraction){
			Unr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "U"))
			Uname = phenomat[Unr,"name"]
			calcdatamat[,Uname] = calcdatamat[,Uname]*(crbybrestimate[expnr])
		}
	}
	
	correcteddatamat = calcdatamat
	results[["correcteddatamat"]] = correcteddatamat
	
	calcdatamatreliable = calcdatamat[reliable,]
	
	### RELEVANT ###
	
	if (is.null(relevant)){
		relevant = sort(unique(phenomat[,"nr"]))
		phenomatrelevant = phenomat[which(phenomat[,"nr"] %in% relevant),]
		if (check){cat(labelingtime,"min experiments relevant for calculation: \n")}
		if (check){print(phenomat[which(phenomat[,"nr"] %in% relevant),])}
	}
	else {
		phenomatrelevant = phenomat[which(phenomat[,"nr"] %in% relevant),]
		if (check){cat(labelingtime,"min experiments relevant for calculation: \n")}
		if (check){print(phenomat[which(phenomat[,"nr"] %in% relevant),])}
	}
	
	### INITIALS ###
	
	if (is.null(initialdummy)){initials = calcdatamat[,phenomat[which(phenomat[,"fraction"] == "T"),"name"]]}
	
	### L/T ###
	
	LtoTmat = as.matrix(calcdatamat[,phenomat[which(phenomat[,"fraction"] == "L"),"name"]])
	if (ncol(as.matrix(LtoTmat)) == ncol(as.matrix(initials))) {LtoTmat = as.matrix(LtoTmat/initials)} else {LtoTmat = as.matrix(LtoTmat/apply(initials,1,median))}
	colnames(LtoTmat) = phenomat[which(phenomat[,"fraction"] == "T"),"name"]
	LtoTmat = cbind(LtoTmat,median = apply(as.matrix(LtoTmat[,phenomatrelevant[which(phenomatrelevant[,"fraction"] == "T"),"name"]]),1,median))
	if (mediancenter){LtoTmat = medctr(LtoTmat,reliable,logscale = FALSE,protocol = FALSE)}
	nrexp = ncol(LtoTmat)
	results[["LtoTmat"]] = LtoTmat
	LtoT = LtoTmat[,nrexp]
	results[["LtoT"]] = LtoT
	
	LTmat = as.matrix(calcdatamat[,phenomat[which(phenomat[,"fraction"] == "T"),"name"]]-calcdatamat[,phenomat[which(phenomat[,"fraction"] == "L"),"name"]])
	if (ncol(as.matrix(LTmat)) == ncol(as.matrix(initials))) {LTmat = as.matrix(LTmat/initials)} else {LTmat = as.matrix(LTmat/apply(initials,1,median))}
	colnames(LTmat) = phenomat[which(phenomat[,"fraction"] == "T"),"name"]
	LTmat = cbind(LTmat,median = apply(as.matrix(LTmat[,phenomatrelevant[which(phenomatrelevant[,"fraction"] == "T"),"name"]]),1,median))
	if (mediancenter){LTmat = medctr(LTmat,reliable,logscale = FALSE,protocol = FALSE)}
	
	### U/T ###
	
	if (unlabeledfraction){
		UtoTmat = as.matrix(calcdatamat[,phenomat[which(phenomat[,"fraction"] == "U"),"name"]])
		if (ncol(as.matrix(UtoTmat)) == ncol(as.matrix(initials))) {UtoTmat = as.matrix(UtoTmat/initials)} else {UtoTmat = as.matrix(UtoTmat/apply(initials,1,median))}
		colnames(UtoTmat) = phenomat[which(phenomat[,"fraction"] == "T"),"name"]
		UtoTmat = cbind(UtoTmat,median = apply(as.matrix(UtoTmat[,phenomatrelevant[which(phenomatrelevant[,"fraction"] == "T"),"name"]]),1,median))
		if (mediancenter){UtoTmat = medctr(UtoTmat,reliable,logscale = FALSE,protocol = FALSE)}
		nrexp = ncol(UtoTmat)
		results[["UtoTmat"]] = UtoTmat
		UtoT = UtoTmat[,nrexp]
		results[["UtoT"]] = UtoT
		
		UTmat = as.matrix(calcdatamat[,phenomat[which(phenomat[,"fraction"] == "U"),"name"]])
		if (ncol(as.matrix(UTmat)) == ncol(as.matrix(initials))) {UTmat = as.matrix(UTmat/initials)} else {UTmat = as.matrix(UTmat/apply(initials,1,median))}
		colnames(UTmat) = phenomat[which(phenomat[,"fraction"] == "T"),"name"]
		UTmat = cbind(UTmat,median = apply(as.matrix(UTmat[,phenomatrelevant[which(phenomatrelevant[,"fraction"] == "T"),"name"]]),1,median))
		if (mediancenter){UTmat = medctr(UTmat,reliable,logscale = FALSE,protocol = FALSE)}
	}
	
	if (usefractions == "LandT"){
		if (check){
			if (!is.null(initialdummy)){
				main = expression(paste("Rank heatpairs of  ",(T['gr,i']-(c['r,i']/a['r,i']*l['gr,i'])*L['gr,i'])/T['gr,i'-1]))
			} else {
				main = expression(paste("Rank heatpairs of  ",1-(c[r]/a[r]*l[gr])*L[gr]/T[gr]))
			}
			if (ncol(LTmat) > 2){
				rcex = ncol(LTmat)/2
				plotsfkt = function(){
					parfkt("rankpairs",ncol(LTmat))
					heatmat = apply(LTmat,2,rank)
					colnames(heatmat) = colnames(LTmat)
					rownames(heatmat) = rownames(LTmat)
					heatpairs(heatmat,main="",cex.labels = min(1/(max(nchar(colnames(LTmat)))/10)*3.5,2))
					sub = paste("( max median fold = ",round(max(abs(diff(apply(LTmat,2,median)))),2),")")
					mtextfkt("rankpairs",ncol(LTmat),main)
				}
				DTA.plot.it(filename = paste(folder,"/rank_heatpairs_",condition,"_",timepoint,sep=""),sw = resolution*rcex,sh = resolution*rcex,sres = resolution,plotsfkt = plotsfkt,ww = rcex*3.5,wh = rcex*3.5,saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)
			}
			
			plotsfkt = function(){
				par(mfrow=windowxy(ncol(LTmat)))
				parfkt("default",ncol(LTmat))
				for (i in 1:ncol(LTmat)){
					assess = function(x){- alpha - ((1/labelingtime)*log(x))}
					if (!is.null(initialdummy)){
						xlab = expression(paste((T['gr,i']-(c['r,i']/a['r,i']*l['gr,i'])*L['gr,i'])/T['gr,i'-1]))
						ylab = expression(paste("Decay rate  ",lambda['gr,i']))
						main = expression(paste("Limit assessment of  ",lambda['gr,i']))
						sub = paste("(",colnames(LTmat)[i],")")
					} else {
						xlab = expression(paste(1-(c[r]/a[r]*l[gr])*L[gr]/T[gr]))
						ylab = expression(paste("Decay rate  ",lambda[gr]))
						main = expression(paste("Limit assessment of  ",lambda[gr]))
						sub = paste("(",colnames(LTmat)[i],")")
					}
					plot(0,type="n",xlim=c(-2*exp(-alpha*labelingtime),3*exp(-alpha*labelingtime)),ylim=c(-0.5,max(hist(LTmat[,ncol(LTmat)],breaks=seq(floor(min(LTmat)),ceiling(max(LTmat)),0.25),plot=FALSE)$density)/2*3),xlab="",ylab="",main="")
					mtextfkt("default",ncol(LTmat),main,xlab,ylab,sub)
					hist(LTmat[,i],add=TRUE,freq=FALSE,breaks=seq(floor(min(LTmat)),ceiling(max(LTmat)),0.25),col="lightgrey")
					plot(assess,add=TRUE,xlim=c(10^-8,5),col="black",lwd=2)
					abline(h=0,col="darkred")
					abline(v=0,col="darkred",lty=2)
					abline(v=exp(-alpha*labelingtime),col="darkred",lty=2)
					text(exp(-alpha*labelingtime)/2,-0.25,sum(LTmat[,i] > 0 & LTmat[,i] < exp(-alpha*labelingtime)),cex=1.5)
					text(-exp(-alpha*labelingtime)/2,-0.25,sum(LTmat[,i] <= 0),cex=1.5)
					text(3*exp(-alpha*labelingtime)/2,-0.25,sum(LTmat[,i] >= exp(-alpha*labelingtime)),cex=1.5)
				}
			}
			DTA.plot.it(filename = paste(folder,"/range_assessment_",condition,"_",timepoint,sep=""),sw = resolution*windowxy(ncol(LTmat))[2],sh = resolution*windowxy(ncol(LTmat))[1],sres = resolution,plotsfkt = plotsfkt,ww = 7*windowxy(ncol(LTmat))[2],wh = 7*windowxy(ncol(LTmat))[1],saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)
		}
		
		nrexp = ncol(LTmat)
		
		results[["LTmat"]] = LTmat
		LT = LTmat[,nrexp]
		results[["LT"]] = LT
		
		options(warn = -1)
		LTdrmat = - alpha - ((1/labelingtime)*log(LTmat))
		options(warn = 0)
		rownames(LTdrmat) = rownames(calcdatamat)
		colnames(LTdrmat) = colnames(LTmat)
		
		LTdrmat[which(LTdrmat <= 0)] = NA
		drmat = LTdrmat
		
		if (check){
			cat("NAs produced in",labelingtime,"min experiments (LandT): \n")
			print(apply(LTdrmat,2,function(x){sum(is.na(x))}))
		}
	}
	if (usefractions == "UandT"){
		if (check){
			if (!is.null(initialdummy)){
				main = expression(paste("Rank heatpairs of  ",(c['r,i']/b['r,i']*u['gr,i'])*U['gr,i']/T['gr,i'-1]))
			} else {
				main = expression(paste("Rank heatpairs of  ",(c[r]/b[r]*u[gr])*U[gr]/T[gr]))
			}
			if (ncol(UTmat) > 2){
				rcex = ncol(UTmat)/2
				plotsfkt = function(){
					parfkt("rankpairs",ncol(UTmat))
					heatmat = apply(UTmat,2,rank)
					colnames(heatmat) = colnames(UTmat)
					rownames(heatmat) = rownames(UTmat)
					heatpairs(heatmat,main="",cex.labels = min(1/(max(nchar(colnames(UTmat)))/10)*3.5,2))
					sub = paste("( max median fold = ",round(max(abs(diff(apply(UTmat,2,median)))),2),")")
					mtextfkt("rankpairs",ncol(UTmat),main)
				}
				DTA.plot.it(filename = paste(folder,"/rank_heatpairs_",condition,"_",timepoint,sep=""),sw = resolution*rcex,sh = resolution*rcex,sres = resolution,plotsfkt = plotsfkt,ww = rcex*3.5,wh = rcex*3.5,saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)
			}
			
			plotsfkt = function(){
				par(mfrow=windowxy(ncol(UTmat)))
				parfkt("default",ncol(UTmat))
				for (i in 1:ncol(UTmat)){
					assess = function(x){- alpha - ((1/labelingtime)*log(x))}
					if (!is.null(initialdummy)){
						xlab = expression(paste((c['r,i']/b['r,i']*u['gr,i'])*U['gr,i']/T['gr,i'-1]))
						ylab = expression(paste("Decay rate  ",lambda['gr,i']))
						main = expression(paste("Limit assessment of  ",lambda['gr,i']))
						sub = paste("(",colnames(UTmat)[i],")")
					} else {
						xlab = expression(paste((c[r]/b[r]*u[gr])*U[gr]/T[gr]))
						ylab = expression(paste("Decay rate  ",lambda[gr]))
						main = expression(paste("Limit assessment of  ",lambda[gr]))
						sub = paste("(",colnames(UTmat)[i],")")
					}
					plot(0,type="n",xlim=c(-2*exp(-alpha*labelingtime),3*exp(-alpha*labelingtime)),ylim=c(-0.5,max(hist(UTmat[,ncol(UTmat)],breaks=seq(floor(min(UTmat)),ceiling(max(UTmat)),0.25),plot=FALSE)$density)/2*3),xlab="",ylab="",main="")
					mtextfkt("default",ncol(UTmat),main,xlab,ylab,sub)
					hist(UTmat[,i],add=TRUE,freq=FALSE,breaks=seq(floor(min(UTmat)),ceiling(max(UTmat)),0.25),col="lightgrey")
					plot(assess,add=TRUE,xlim=c(10^-8,5),col="black",lwd=2)
					abline(h=0,col="darkred")
					abline(v=0,col="darkred",lty=2)
					abline(v=exp(-alpha*labelingtime),col="darkred",lty=2)
					text(exp(-alpha*labelingtime)/2,-0.25,sum(UTmat[,i] > 0 & UTmat[,i] < exp(-alpha*labelingtime)),cex=1.5)
					text(-exp(-alpha*labelingtime)/2,-0.25,sum(UTmat[,i] <= 0),cex=1.5)
					text(3*exp(-alpha*labelingtime)/2,-0.25,sum(UTmat[,i] >= exp(-alpha*labelingtime)),cex=1.5)
				}
			}
			DTA.plot.it(filename = paste(folder,"/range_assessment_",condition,"_",timepoint,sep=""),sw = resolution*windowxy(ncol(UTmat))[2],sh = resolution*windowxy(ncol(UTmat))[1],sres = resolution,plotsfkt = plotsfkt,ww = 7*windowxy(ncol(UTmat))[2],wh = 7*windowxy(ncol(UTmat))[1],saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)
		}
		
		nrexp = ncol(UTmat)
		
		results[["UTmat"]] = UTmat
		UT = UTmat[,nrexp]
		results[["UT"]] = UT
		
		options(warn = -1)
		UTdrmat = - alpha - ((1/labelingtime)*log(UTmat))
		options(warn = 0)
		rownames(UTdrmat) = rownames(calcdatamat)
		colnames(UTdrmat) = colnames(UTmat)
		
		UTdrmat[which(UTdrmat <= 0)] = NA
		drmat = UTdrmat
		
		if (check){
			cat("NAs produced in",labelingtime,"min experiments (UandT): \n")
			print(apply(UTdrmat,2,function(x){sum(is.na(x))}))
		}
	}
	if (usefractions == "both"){
		
		Bmat = (LTmat + UTmat)/2
		
		if (check){
			if (!is.null(initialdummy)){
				main = expression(paste("Rank heatpairs of  ",((T['gr,i']-(c['r,i']/a['r,i']*l['gr,i'])*L['gr,i'])/T['gr,i'-1]+(c['r,i']/b['r,i']*u['gr,i'])*U['gr,i']/T['gr,i'-1])/2))
			} else {
				main = expression(paste("Rank heatpairs of  ",(1-(c[r]/a[r]*l[gr])*L[gr]/T[gr]+(c[r]/b[r]*u[gr])*U[gr]/T[gr])/2))
			}
			if (ncol(Bmat) > 2){
				rcex = ncol(Bmat)/2
				plotsfkt = function(){
					parfkt("rankpairs",ncol(Bmat))
					heatmat = apply(Bmat,2,rank)
					colnames(heatmat) = colnames(Bmat)
					rownames(heatmat) = rownames(Bmat)
					heatpairs(heatmat,main="",cex.labels = min(1/(max(nchar(colnames(Bmat)))/10)*3.5,2))
					sub = paste("( max median fold = ",round(max(abs(diff(apply(Bmat,2,median)))),2),")")
					mtextfkt("rankpairs",ncol(Bmat),main)
				}
				DTA.plot.it(filename = paste(folder,"/rank_heatpairs_",condition,"_",timepoint,sep=""),sw = resolution*rcex,sh = resolution*rcex,sres = resolution,plotsfkt = plotsfkt,ww = rcex*3.5,wh = rcex*3.5,saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)
			}
			
			plotsfkt = function(){
				par(mfrow=windowxy(ncol(Bmat)))
				parfkt("default",ncol(Bmat))
				for (i in 1:ncol(Bmat)){
					assess = function(x){- alpha - ((1/labelingtime)*log(x))}
					if (!is.null(initialdummy)){
						xlab = expression(paste(((T['gr,i']-(c['r,i']/a['r,i']*l['gr,i'])*L['gr,i'])/T['gr,i'-1]+(c['r,i']/b['r,i']*u['gr,i'])*U['gr,i']/T['gr,i'-1])/2))
						ylab = expression(paste("Decay rate  ",lambda['gr,i']))
						main = expression(paste("Limit assessment of  ",lambda['gr,i']))
						sub = paste("(",colnames(Bmat)[i],")")
					} else {
						xlab = expression(paste((1-(c[r]/a[r]*l[gr])*L[gr]/T[gr]+(c[r]/b[r]*u[gr])*U[gr]/T[gr])/2))
						ylab = expression(paste("Decay rate  ",lambda[gr]))
						main = expression(paste("Limit assessment of  ",lambda[gr]))
						sub = paste("(",colnames(Bmat)[i],")")
					}
					plot(0,type="n",xlim=c(-2*exp(-alpha*labelingtime),3*exp(-alpha*labelingtime)),ylim=c(-0.5,max(hist(Bmat[,ncol(Bmat)],breaks=seq(floor(min(Bmat)),ceiling(max(Bmat)),0.25),plot=FALSE)$density)/2*3),xlab="",ylab="",main="")
					mtextfkt("default",ncol(Bmat),main,xlab,ylab,sub)
					hist(Bmat[,i],add=TRUE,freq=FALSE,breaks=seq(floor(min(Bmat)),ceiling(max(Bmat)),0.25),col="lightgrey")
					plot(assess,add=TRUE,xlim=c(10^-8,5),col="black",lwd=2)
					abline(h=0,col="darkred")
					abline(v=0,col="darkred",lty=2)
					abline(v=exp(-alpha*labelingtime),col="darkred",lty=2)
					text(exp(-alpha*labelingtime)/2,-0.25,sum(Bmat[,i] > 0 & Bmat[,i] < exp(-alpha*labelingtime)),cex=1.5)
					text(-exp(-alpha*labelingtime)/2,-0.25,sum(Bmat[,i] <= 0),cex=1.5)
					text(3*exp(-alpha*labelingtime)/2,-0.25,sum(Bmat[,i] >= exp(-alpha*labelingtime)),cex=1.5)
				}
			}
			DTA.plot.it(filename = paste(folder,"/range_assessment_",condition,"_",timepoint,sep=""),sw = resolution*windowxy(ncol(Bmat))[2],sh = resolution*windowxy(ncol(Bmat))[1],sres = resolution,plotsfkt = plotsfkt,ww = 7*windowxy(ncol(Bmat))[2],wh = 7*windowxy(ncol(Bmat))[1],saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)
		}
		
		nrexp = ncol(Bmat)
		
		results[["Bmat"]] = Bmat
		B = Bmat[,nrexp]
		results[["B"]] = B
		
		options(warn = -1)
		Bdrmat = - alpha - ((1/labelingtime)*log(Bmat))
		options(warn = 0)
		rownames(Bdrmat) = rownames(calcdatamat)
		colnames(Bdrmat) = colnames(Bmat)
		
		Bdrmat[which(Bdrmat <= 0)] = NA
		drmat = Bdrmat
		
		if (check){
			cat("NAs produced in",labelingtime,"min experiments (both): \n")
			print(apply(Bdrmat,2,function(x){sum(is.na(x))}))
		}
	}
	
	results[["drmat"]] = drmat
	results[["hlmat"]] = log(2)/drmat
	dr = drmat[,nrexp]
	results[["dr"]] = dr
	results[["hl"]] = log(2)/dr
	
	### TE ###
	
	TEmat = cbind(calcdatamat[,which(phenomat[,"fraction"] == "T")])
	colnames(TEmat) = phenomat[which(phenomat[,"fraction"] == "T"),"name"]
	TEmat = cbind(TEmat,median = apply(as.matrix(TEmat[,phenomatrelevant[which(phenomatrelevant[,"fraction"] == "T"),"name"]]),1,median))
	colnames(TEmat) = c(phenomat[which(phenomat[,"fraction"] == "T"),"name"],"median")
	results[["TEmat"]] = TEmat
	TE = TEmat[,nrexp]
	results[["TE"]] = TE
	
	### LE ###
	
	LEmat = cbind(calcdatamat[,which(phenomat[,"fraction"] == "L")])
	colnames(LEmat) = phenomat[which(phenomat[,"fraction"] == "L"),"name"]
	LEmat = cbind(LEmat,median = apply(as.matrix(LEmat[,phenomatrelevant[which(phenomatrelevant[,"fraction"] == "L"),"name"]]),1,median))
	colnames(LEmat) = c(phenomat[which(phenomat[,"fraction"] == "L"),"name"],"median")
	results[["LEmat"]] = LEmat
	LE = LEmat[,nrexp]
	results[["LE"]] = LE
	
	### UE ###
	
	if (unlabeledfraction){
		UEmat = cbind(calcdatamat[,which(phenomat[,"fraction"] == "U")])
		colnames(UEmat) = phenomat[which(phenomat[,"fraction"] == "U"),"name"]
		UEmat = cbind(UEmat,median = apply(as.matrix(UEmat[,phenomatrelevant[which(phenomatrelevant[,"fraction"] == "U"),"name"]]),1,median))
		colnames(UEmat) = c(phenomat[which(phenomat[,"fraction"] == "U"),"name"],"median")
		results[["UEmat"]] = UEmat
		UE = UEmat[,nrexp]
		results[["UE"]] = UE
	}
	
	### ESTIMATION OF sr ###
	
	srmat = LEmat*(alpha + drmat)/(exp(alpha*labelingtime)-exp(-drmat*labelingtime))
	colnames(srmat) = c(phenomat[which(phenomat[,"fraction"] == "T"),"name"],"median")
	results[["srmat"]] = srmat
	sr = srmat[,nrexp]
	results[["sr"]] = sr
	
	if (!is.null(mRNAs)){
		Rsrmat = t(((totaloverwt*mRNAs)/apply(TEmat,2,sum))*t(srmat))*ccl
		colnames(Rsrmat) = c(phenomat[which(phenomat[,"fraction"] == "T"),"name"],"median")
		results[["Rsrmat"]] = Rsrmat
		Rsr = Rsrmat[,nrexp]
		results[["Rsr"]] = Rsr
	}
	
	### ESTIMATION OF globaldr ###
	
	globaldrmat = 1/apply(TEmat,2,sum)*apply(TEmat*drmat,2,sum)
	results[["globaldrmat"]] = globaldrmat
	globaldr = globaldrmat[nrexp]
	results[["globaldr"]] = globaldr
	
	### CORRELATION ###
	
	if (check) {
		vecmat = cbind(TE,LE,tnumber,sr,dr,(log(2)/dr))
		vecmat = vecmat[reliable,]
		compsmat = matrix(0,ncol=3,nrow=3)
		for (i in 1:3){for (j in 1:3){compsmat[i,j] = cor(vecmat[,i],vecmat[,j+3],method = "pearson",use = "na.or.complete")}}
		rownames(compsmat) = c("TE","LE","#U")
		colnames(compsmat) = c("SR","DR","HL")
		print(t(compsmat))
		
		plotsfkt = function(){
			par(mfrow=c(1,2))
			layout(matrix(data=c(1,2),nrow=1,ncol=2),widths=c(6,1)) 
			parfkt("correlationleft",1)
			image(t(t(compsmat)[3:1,]),axes=FALSE,col = colorRampPalette(c("darkred","lightgrey","darkgreen"))(500),breaks = seq(-1,1,length.out=501),main="",xlab="",ylab="")
			xlab = expression("Data inputs")
			ylab = expression("Extracted rates")
			main = expression(paste("Correlation analysis"))
			sub = paste("( pearson correlation for gene-wise medians )")
			mtextfkt("default",1,main,xlab,ylab,sub)
			axis(1,at=c(0,0.5,1),labels=c(expression(T[gr]),expression(L[gr]),expression('#u'[g])))
			axis(2,at=c(0,0.5,1),labels=c(expression(t[1/2][gr]),expression(lambda[gr]),expression(mu[gr])),las=2)
			abline(v=c(0.75,0.25),h=c(0.75,0.25),cex=0.5)
			box()
			text(c(0,0.5,1,0,0.5,1,0,0.5,1),c(0,0,0,0.5,0.5,0.5,1,1,1),round(as.vector(t(t(compsmat)[3:1,])),2))
			parfkt("correlationright",1)
			image(1,seq(-1,1,length.out=501),matrix(data=seq(-1,1,length.out=501),ncol=length(seq(-1,1,length.out=501)),nrow=1),col=colorRampPalette(c("darkred","lightgrey","darkgreen"))(500),xlab="",ylab="",axes=FALSE)
			axis(4,at=c(-1,-0.5,0,0.5,1),labels=c(-1,-0.5,0,0.5,1),las=2)
			box()
		}
		DTA.plot.it(filename = paste(folder,"/correlation_analysis_",condition,"_",timepoint,sep=""),sw = resolution,sh = resolution,sres = resolution,plotsfkt = plotsfkt,ww = 7,wh = 7,saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)
	}
	
	### ERROR MODEL ###
			
	if (error){
		TEmat = calcdatamat[,phenomatrelevant[which(phenomatrelevant[,"fraction"] == "T"),"name"]]
		LEmat = calcdatamat[,phenomatrelevant[which(phenomatrelevant[,"fraction"] == "L"),"name"]]
		if (unlabeledfraction){UEmat = calcdatamat[,phenomatrelevant[which(phenomatrelevant[,"fraction"] == "U"),"name"]]}
		if (dynamic){if (!is.null(initialdummy)){TImat = initials} else {TImat =  c()}}
		if (dynamic){if (!is.null(initialdummy)){errorpossible = (!is.null(dim(TImat)) & !is.null(dim(TEmat)))} else {errorpossible = !is.null(dim(TEmat))}} else {errorpossible = !is.null(dim(TEmat))}
		if (errorpossible){
			if (!is.null(dim(TEmat))){
				cat("Assessing error in total expression measurements ... \n")
				TEmat.norm = quantnorm(log(TEmat))
				TEmat.norm.sd = apply(TEmat.norm,1,sd)
				TEmat.norm.mean = apply(TEmat.norm,1,mean)
				TEmat.norm.loess = loess(TEmat.norm.sd ~ TEmat.norm.mean,degree = 1)
				TEmat.gamma.var = gamma.variance(TEmat.norm.loess$residuals - min(TEmat.norm.loess$residuals))
				TEmat.feat = cbind(TEmat.norm,TEmat.norm.mean,TEmat.norm.loess$fitted)
				nr.col = ncol(TEmat.norm)
				options(warn = -1)
				TEmat.norm.reg.sd = apply(TEmat.feat,1,function(x){likelihood.fct = function(q){likelihood(x[1:nr.col],x[nr.col + 1],q,x[nr.col + 2],TEmat.gamma.var[1])};optimize(likelihood.fct,interval = c(0,1),maximum = TRUE)$maximum})
				options(warn = 0)
				TEmat.norm.reg.sd[is.na(TEmat.norm.reg.sd)] = TEmat.norm.sd[is.na(TEmat.norm.reg.sd)]
				if (check) {
					plotsfkt = function(){
						par(mfrow = c(1,2))
						par(mar = c(1,4,4,2)+0.1+1)
						par(cex = 1)
						par(cex.axis = 1.5*0.75)
						min.mean = min(TEmat.norm.mean,na.rm=TRUE)
						max.mean = max(TEmat.norm.mean,na.rm=TRUE)
						max.sd = max(TEmat.norm.sd,na.rm=TRUE)
						x <- seq(min.mean,max.mean,length.out=6)
						y <- seq(0,max.sd,length.out=6)
						z <- matrix(0,6,6)
						residuals = TEmat.norm.loess$residuals - min(TEmat.norm.loess$residuals)
						res.hist = hist(residuals,plot = FALSE)
						seq.gamma = seq(0,max.sd,length.out = 1000)
						persp(x,y,z,theta = 60,phi = 30,expand = 0.5,ltheta = 120,shade = NA,ticktype = "detailed",xlab = "",ylab = "",zlab = "",zlim = c(0,max(res.hist$density)*5/4),lwd = 0.5,axes = FALSE,box = FALSE,bg = "transparent") -> persp.res
						reliable.range = range(TEmat.norm.mean[reliable],na.rm = TRUE)
						polygon(trans3d(c(reliable.range[1],reliable.range[1],reliable.range[2],reliable.range[2]),c(0,max.sd,max.sd,0),0,pmat = persp.res),col = convertcolor("lightslategray",20))
						points(trans3d(TEmat.norm.mean,TEmat.norm.sd,0,pmat = persp.res),col = "black",pch = ".")
						lines(trans3d(TEmat.norm.loess$x[order(TEmat.norm.loess$x)],TEmat.norm.loess$fitted[order(TEmat.norm.loess$x)],0,pmat = persp.res),col = "black",lwd = 2)
						offset = max.sd/25
						lines(trans3d(x,0 - offset,0,pmat = persp.res),lwd = 2)
						text(trans3d(pretty(x)[pretty(x) > min.mean & pretty(x) < max.mean],0 - 2*offset,0,pmat = persp.res),as.character(pretty(x)[pretty(x) > min.mean & pretty(x) < max.mean]))
						lines(trans3d(max.mean + offset*(max.mean-min.mean)/max.sd,y,0,pmat = persp.res),lwd = 2)
						text(trans3d(max.mean + 2*offset*(max.mean-min.mean)/max.sd,pretty(y)[pretty(y) < max.sd],0,pmat = persp.res),as.character(pretty(y)[pretty(y) < max.sd]))						
						for (j in 1:length(res.hist$density)){
							polygon(trans3d(min.mean-(max.mean-min.mean)/10,c(res.hist$breaks[j],res.hist$breaks[j],res.hist$breaks[j+1],res.hist$breaks[j+1]),c(0,res.hist$density[j],res.hist$density[j],0),pmat = persp.res),col = convertcolor("lightgray",60))	
						}
						residuals.gamma = c(0,dgamma(x = seq.gamma, shape = TEmat.gamma.var[2], scale = TEmat.gamma.var[3]),0)
						polygon(trans3d(x = min.mean-(max.mean-min.mean)/10,y = c(0,seq.gamma,max.sd),z = pmin(pmax(residuals.gamma,0),max(res.hist$density)*3),pmat = persp.res),col = "#FF000030",border = "darkred",lwd=3)
						arrow.x.red = c(0,0.8,0.8,1,0.8,0.8,0)
						arrow.y.red = c(-0.05,-0.05,-0.1,0,0.1,0.05,0.05)
						max.gamma = seq.gamma[which.min(abs(residuals.gamma-max(residuals.gamma)))]
						i = which.min(abs(TEmat.norm.loess$x-mean(TEmat.norm.mean[reliable])))
						polygon(trans3d(arrow.x.red*(TEmat.norm.loess$x[i]-min.mean)/6*4 + min.mean,(arrow.y.red/5*min(TEmat.gamma.var[2],3)*max.sd + max.gamma),min(max(residuals.gamma),max(res.hist$density))/5*4,pmat = persp.res),col = "darkred",border = "darkred")
						arrow.x.blue = c(0.3,0.8,0.8,1,0.8,0.8,0.3)-0.5
						arrow.y.blue = c(-0.05,-0.05,-0.1,0,0.1,0.05,0.05)
						rect.x.1.blue = c(0.15,0.225,0.225,0.15)-0.5
						rect.x.2.blue = c(0,0.075,0.075,0)-0.5
						rect.y.blue = c(-0.05,-0.05,0.05,0.05)
						polygon(trans3d(arrow.x.blue*0.8*(TEmat.norm.loess$x[i]-min.mean)+min.mean+(TEmat.norm.loess$x[i]-min.mean)/2,(arrow.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
						polygon(trans3d(rect.x.1.blue*0.8*(TEmat.norm.loess$x[i]-min.mean)+min.mean+(TEmat.norm.loess$x[i]-min.mean)/2,(rect.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
						polygon(trans3d(rect.x.2.blue*0.8*(TEmat.norm.loess$x[i]-min.mean)+min.mean+(TEmat.norm.loess$x[i]-min.mean)/2,(rect.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
						polygon(trans3d(x = TEmat.norm.loess$x[i],y = c(0,seq.gamma,max.sd),z = pmin(pmax(c(0,dgamma(x = seq.gamma, shape = (TEmat.norm.loess$fitted[i]^2)/TEmat.gamma.var[1], scale = TEmat.gamma.var[1]/TEmat.norm.loess$fitted[i]),0),0),max(res.hist$density)*3),pmat = persp.res),col = "#0000FF30",border = "darkblue",lwd=3)
						polygon(trans3d((-arrow.x.blue)*0.8*(max.mean-TEmat.norm.loess$x[i])+max.mean-(max.mean-TEmat.norm.loess$x[i])/2,(arrow.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
						polygon(trans3d((-rect.x.1.blue)*0.8*(max.mean-TEmat.norm.loess$x[i])+max.mean-(max.mean-TEmat.norm.loess$x[i])/2,(rect.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
						polygon(trans3d((-rect.x.2.blue)*0.8*(max.mean-TEmat.norm.loess$x[i])+max.mean-(max.mean-TEmat.norm.loess$x[i])/2,(rect.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
						xlab = expression(paste("Mean  ",T[gr]))
						ylab = expression(paste("Standard deviation  ",sigma,"(",T[gr],")"))
						main = expression(paste("Mean  ",T[gr],"  versus standard deviation  ",sigma,"(",T[gr],")"))
						sub = paste("( smoothed loess fit of ",dim(TEmat)[2]," replicates )")
						mtext(main,3,2,cex=2*0.75)
						mtext(sub,3,0,col="darkgrey",cex=1.25*0.75)
						text(trans3d(min.mean + (max.mean-min.mean)*0.6,-max.sd*0.2,-max(res.hist$density)/7,pmat = persp.res),xlab,cex = 1.65*0.75,srt = 90)
						text(trans3d(max.mean + (max.mean-min.mean)/10,max.sd*0.6,-max(res.hist$density),pmat = persp.res),ylab,cex = 1.65*0.75,srt = 0)
						legend("topright",c(expression(paste(Gamma," - distributed std. dev.  ",sigma[g])),expression(paste(Gamma," - distributed residuals  ",r[g]))),pt.bg=c("darkblue","darkred"),col=rep("black",2),bg="white",pch=21,cex=1,pt.cex=1,inset = 0.05)
						par(mar = c(1,4,4,2)+0.1+1)
						par(mai = c(1.1,1.1,1.2,0.7))
						par(cex = 1)
						par(cex.axis = 1.5*0.75)
						heatscatter(TEmat.norm.sd,TEmat.norm.reg.sd,xlab="",ylab="",main="",cor=FALSE,xlim = c(0,max(c(TEmat.norm.sd,TEmat.norm.reg.sd),na.rm = TRUE)),ylim = c(0,max(c(TEmat.norm.sd,TEmat.norm.reg.sd),na.rm = TRUE)))
						abline(0,1,lwd=1)
						xlab = expression(paste("Empirical standard deviation  ",sigma,"(",T[gr],")"))
						ylab = expression(paste("Regularized standard deviation  ",sigma^R,"(",T[gr],")"))
						main = expression(paste("Regularized versus empirical standard deviation"))
						mtext(main,3,3,cex=2*0.75)
						mtext(xlab,1,3,cex=1.65*0.75)
						mtext(ylab,2,3,cex=1.65*0.75)
					}
					DTA.plot.it(filename = paste(folder,"/error_assessment_TE_",condition,"_",timepoint,sep=""),sw = resolution*2,sh = resolution,sres = resolution,plotsfkt = plotsfkt,ww = 14,wh = 7,saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)
				}
				results[["TE.log.error"]] = cbind("log.mean" = TEmat.norm.mean,"log.sd" = TEmat.norm.reg.sd,"log.cv" = TEmat.norm.reg.sd/TEmat.norm.mean)		
			}	
			if (!is.null(dim(LEmat))){
				cat("Assessing error in labeled expression measurements ... \n")
				LEmat.norm = quantnorm(log(LEmat))
				LEmat.norm.sd = apply(LEmat.norm,1,sd)
				LEmat.norm.mean = apply(LEmat.norm,1,mean)
				LEmat.norm.loess = loess(LEmat.norm.sd ~ LEmat.norm.mean,degree = 1)
				LEmat.gamma.var = gamma.variance(LEmat.norm.loess$residuals - min(LEmat.norm.loess$residuals))
				LEmat.feat = cbind(LEmat.norm,LEmat.norm.mean,LEmat.norm.loess$fitted)
				nr.col = ncol(LEmat.norm)
				options(warn = -1)
				LEmat.norm.reg.sd = apply(LEmat.feat,1,function(x){likelihood.fct = function(q){likelihood(x[1:nr.col],x[nr.col + 1],q,x[nr.col + 2],LEmat.gamma.var[1])};optimize(likelihood.fct,interval = c(0,1),maximum = TRUE)$maximum})
				options(warn = 0)
				LEmat.norm.reg.sd[is.na(LEmat.norm.reg.sd)] = LEmat.norm.sd[is.na(LEmat.norm.reg.sd)]
				if (check) {
					plotsfkt = function(){
						par(mfrow = c(1,2))
						par(mar = c(1,4,4,2)+0.1+1)
						par(cex = 1)
						par(cex.axis = 1.5*0.75)
						min.mean = min(LEmat.norm.mean,na.rm=TRUE)
						max.mean = max(LEmat.norm.mean,na.rm=TRUE)
						max.sd = max(LEmat.norm.sd,na.rm=TRUE)
						x <- seq(min.mean,max.mean,length.out=6)
						y <- seq(0,max.sd,length.out=6)
						z <- matrix(0,6,6)
						residuals = LEmat.norm.loess$residuals - min(LEmat.norm.loess$residuals)
						res.hist = hist(residuals,plot = FALSE)
						seq.gamma = seq(0,max.sd,length.out = 1000)
						persp(x,y,z,theta = 60,phi = 30,expand = 0.5,ltheta = 120,shade = NA,ticktype = "detailed",xlab = "",ylab = "",zlab = "",zlim = c(0,max(res.hist$density)*5/4),lwd = 0.5,axes = FALSE,box = FALSE,bg = "transparent") -> persp.res
						reliable.range = range(LEmat.norm.mean[reliable],na.rm = TRUE)
						polygon(trans3d(c(reliable.range[1],reliable.range[1],reliable.range[2],reliable.range[2]),c(0,max.sd,max.sd,0),0,pmat = persp.res),col = convertcolor("lightslategray",20))
						points(trans3d(LEmat.norm.mean,LEmat.norm.sd,0,pmat = persp.res),col = "black",pch = ".")
						lines(trans3d(LEmat.norm.loess$x[order(LEmat.norm.loess$x)],LEmat.norm.loess$fitted[order(LEmat.norm.loess$x)],0,pmat = persp.res),col = "black",lwd = 2)
						offset = max.sd/25
						lines(trans3d(x,0 - offset,0,pmat = persp.res),lwd = 2)
						text(trans3d(pretty(x)[pretty(x) > min.mean & pretty(x) < max.mean],0 - 2*offset,0,pmat = persp.res),as.character(pretty(x)[pretty(x) > min.mean & pretty(x) < max.mean]))
						lines(trans3d(max.mean + offset*(max.mean-min.mean)/max.sd,y,0,pmat = persp.res),lwd = 2)
						text(trans3d(max.mean + 2*offset*(max.mean-min.mean)/max.sd,pretty(y)[pretty(y) < max.sd],0,pmat = persp.res),as.character(pretty(y)[pretty(y) < max.sd]))						
						for (j in 1:length(res.hist$density)){
							polygon(trans3d(min.mean-(max.mean-min.mean)/10,c(res.hist$breaks[j],res.hist$breaks[j],res.hist$breaks[j+1],res.hist$breaks[j+1]),c(0,res.hist$density[j],res.hist$density[j],0),pmat = persp.res),col = convertcolor("lightgray",60))	
						}
						residuals.gamma = c(0,dgamma(x = seq.gamma, shape = LEmat.gamma.var[2], scale = LEmat.gamma.var[3]),0)
						polygon(trans3d(x = min.mean-(max.mean-min.mean)/10,y = c(0,seq.gamma,max.sd),z = pmin(pmax(residuals.gamma,0),max(res.hist$density)*3),pmat = persp.res),col = "#FF000030",border = "darkred",lwd=3)
						arrow.x.red = c(0,0.8,0.8,1,0.8,0.8,0)
						arrow.y.red = c(-0.05,-0.05,-0.1,0,0.1,0.05,0.05)
						max.gamma = seq.gamma[which.min(abs(residuals.gamma-max(residuals.gamma)))]
						i = which.min(abs(LEmat.norm.loess$x-mean(LEmat.norm.mean[reliable])))
						polygon(trans3d(arrow.x.red*(LEmat.norm.loess$x[i]-min.mean)/6*4 + min.mean,(arrow.y.red/5*min(LEmat.gamma.var[2],3)*max.sd + max.gamma),min(max(residuals.gamma),max(res.hist$density))/5*4,pmat = persp.res),col = "darkred",border = "darkred")
						arrow.x.blue = c(0.3,0.8,0.8,1,0.8,0.8,0.3)-0.5
						arrow.y.blue = c(-0.05,-0.05,-0.1,0,0.1,0.05,0.05)
						rect.x.1.blue = c(0.15,0.225,0.225,0.15)-0.5
						rect.x.2.blue = c(0,0.075,0.075,0)-0.5
						rect.y.blue = c(-0.05,-0.05,0.05,0.05)
						polygon(trans3d(arrow.x.blue*0.8*(LEmat.norm.loess$x[i]-min.mean)+min.mean+(LEmat.norm.loess$x[i]-min.mean)/2,(arrow.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
						polygon(trans3d(rect.x.1.blue*0.8*(LEmat.norm.loess$x[i]-min.mean)+min.mean+(LEmat.norm.loess$x[i]-min.mean)/2,(rect.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
						polygon(trans3d(rect.x.2.blue*0.8*(LEmat.norm.loess$x[i]-min.mean)+min.mean+(LEmat.norm.loess$x[i]-min.mean)/2,(rect.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
						polygon(trans3d(x = LEmat.norm.loess$x[i],y = c(0,seq.gamma,max.sd),z = pmin(pmax(c(0,dgamma(x = seq.gamma, shape = (LEmat.norm.loess$fitted[i]^2)/LEmat.gamma.var[1], scale = LEmat.gamma.var[1]/LEmat.norm.loess$fitted[i]),0),0),max(res.hist$density)*3),pmat = persp.res),col = "#0000FF30",border = "darkblue",lwd=3)
						polygon(trans3d((-arrow.x.blue)*0.8*(max.mean-LEmat.norm.loess$x[i])+max.mean-(max.mean-LEmat.norm.loess$x[i])/2,(arrow.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
						polygon(trans3d((-rect.x.1.blue)*0.8*(max.mean-LEmat.norm.loess$x[i])+max.mean-(max.mean-LEmat.norm.loess$x[i])/2,(rect.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
						polygon(trans3d((-rect.x.2.blue)*0.8*(max.mean-LEmat.norm.loess$x[i])+max.mean-(max.mean-LEmat.norm.loess$x[i])/2,(rect.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
						xlab = expression(paste("Mean  ",L[gr]))
						ylab = expression(paste("Standard deviation  ",sigma,"(",L[gr],")"))
						main = expression(paste("Mean  ",L[gr],"  versus standard deviation  ",sigma,"(",L[gr],")"))
						sub = paste("( smoothed loess fit of ",dim(LEmat)[2]," replicates )")
						mtext(main,3,2,cex=2*0.75)
						mtext(sub,3,0,col="darkgrey",cex=1.25*0.75)
						text(trans3d(min.mean + (max.mean-min.mean)*0.6,-max.sd*0.2,-max(res.hist$density)/7,pmat = persp.res),xlab,cex = 1.65*0.75,srt = 90)
						text(trans3d(max.mean + (max.mean-min.mean)/10,max.sd*0.6,-max(res.hist$density),pmat = persp.res),ylab,cex = 1.65*0.75,srt = 0)
						legend("topright",c(expression(paste(Gamma," - distributed std. dev.  ",sigma[g])),expression(paste(Gamma," - distributed residuals  ",r[g]))),pt.bg=c("darkblue","darkred"),col=rep("black",2),bg="white",pch=21,cex=1,pt.cex=1,inset = 0.05)
						par(mar = c(1,4,4,2)+0.1+1)
						par(mai = c(1.1,1.1,1.2,0.7))
						par(cex = 1)
						par(cex.axis = 1.5*0.75)
						heatscatter(LEmat.norm.sd,LEmat.norm.reg.sd,xlab="",ylab="",main="",cor=FALSE,xlim = c(0,max(c(LEmat.norm.sd,LEmat.norm.reg.sd),na.rm = TRUE)),ylim = c(0,max(c(LEmat.norm.sd,LEmat.norm.reg.sd),na.rm = TRUE)))
						abline(0,1,lwd=1)
						xlab = expression(paste("Empirical standard deviation  ",sigma,"(",L[gr],")"))
						ylab = expression(paste("Regularized standard deviation  ",sigma^R,"(",L[gr],")"))
						main = expression(paste("Regularized versus empirical standard deviation"))
						mtext(main,3,3,cex=2*0.75)
						mtext(xlab,1,3,cex=1.65*0.75)
						mtext(ylab,2,3,cex=1.65*0.75)
					}
					DTA.plot.it(filename = paste(folder,"/error_assessment_LE_",condition,"_",timepoint,sep=""),sw = resolution*2,sh = resolution,sres = resolution,plotsfkt = plotsfkt,ww = 14,wh = 7,saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)
				}
				results[["LE.log.error"]] = cbind("log.mean" = LEmat.norm.mean,"log.sd" = LEmat.norm.reg.sd,"log.cv" = LEmat.norm.reg.sd/LEmat.norm.mean)		
			}
			if (unlabeledfraction){
				if (!is.null(dim(UEmat))){
					cat("Assessing error in unlabeled expression measurements ... \n")
					UEmat.norm = quantnorm(log(UEmat))
					UEmat.norm.sd = apply(UEmat.norm,1,sd)
					UEmat.norm.mean = apply(UEmat.norm,1,mean)
					UEmat.norm.loess = loess(UEmat.norm.sd ~ UEmat.norm.mean,degree = 1)
					UEmat.gamma.var = gamma.variance(UEmat.norm.loess$residuals - min(UEmat.norm.loess$residuals))
					UEmat.feat = cbind(UEmat.norm,UEmat.norm.mean,UEmat.norm.loess$fitted)
					nr.col = ncol(UEmat.norm)
					options(warn = -1)
					UEmat.norm.reg.sd = apply(UEmat.feat,1,function(x){likelihood.fct = function(q){likelihood(x[1:nr.col],x[nr.col + 1],q,x[nr.col + 2],UEmat.gamma.var[1])};optimize(likelihood.fct,interval = c(0,1),maximum = TRUE)$maximum})
					options(warn = 0)
					UEmat.norm.reg.sd[is.na(UEmat.norm.reg.sd)] = UEmat.norm.sd[is.na(UEmat.norm.reg.sd)]
					if (check) {
						plotsfkt = function(){
							par(mfrow = c(1,2))
							par(mar = c(1,4,4,2)+0.1+1)
							par(cex = 1)
							par(cex.axis = 1.5*0.75)
							min.mean = min(UEmat.norm.mean,na.rm=TRUE)
							max.mean = max(UEmat.norm.mean,na.rm=TRUE)
							max.sd = max(UEmat.norm.sd,na.rm=TRUE)
							x <- seq(min.mean,max.mean,length.out=6)
							y <- seq(0,max.sd,length.out=6)
							z <- matrix(0,6,6)
							residuals = UEmat.norm.loess$residuals - min(UEmat.norm.loess$residuals)
							res.hist = hist(residuals,plot = FALSE)
							seq.gamma = seq(0,max.sd,length.out = 1000)
							persp(x,y,z,theta = 60,phi = 30,expand = 0.5,ltheta = 120,shade = NA,ticktype = "detailed",xlab = "",ylab = "",zlab = "",zlim = c(0,max(res.hist$density)*5/4),lwd = 0.5,axes = FALSE,box = FALSE,bg = "transparent") -> persp.res
							reliable.range = range(UEmat.norm.mean[reliable],na.rm = TRUE)
							polygon(trans3d(c(reliable.range[1],reliable.range[1],reliable.range[2],reliable.range[2]),c(0,max.sd,max.sd,0),0,pmat = persp.res),col = convertcolor("lightslategray",20))
							points(trans3d(UEmat.norm.mean,UEmat.norm.sd,0,pmat = persp.res),col = "black",pch = ".")
							lines(trans3d(UEmat.norm.loess$x[order(UEmat.norm.loess$x)],UEmat.norm.loess$fitted[order(UEmat.norm.loess$x)],0,pmat = persp.res),col = "black",lwd = 2)
							offset = max.sd/25
							lines(trans3d(x,0 - offset,0,pmat = persp.res),lwd = 2)
							text(trans3d(pretty(x)[pretty(x) > min.mean & pretty(x) < max.mean],0 - 2*offset,0,pmat = persp.res),as.character(pretty(x)[pretty(x) > min.mean & pretty(x) < max.mean]))
							lines(trans3d(max.mean + offset*(max.mean-min.mean)/max.sd,y,0,pmat = persp.res),lwd = 2)
							text(trans3d(max.mean + 2*offset*(max.mean-min.mean)/max.sd,pretty(y)[pretty(y) < max.sd],0,pmat = persp.res),as.character(pretty(y)[pretty(y) < max.sd]))						
							for (j in 1:length(res.hist$density)){
								polygon(trans3d(min.mean-(max.mean-min.mean)/10,c(res.hist$breaks[j],res.hist$breaks[j],res.hist$breaks[j+1],res.hist$breaks[j+1]),c(0,res.hist$density[j],res.hist$density[j],0),pmat = persp.res),col = convertcolor("lightgray",60))	
							}
							residuals.gamma = c(0,dgamma(x = seq.gamma, shape = UEmat.gamma.var[2], scale = UEmat.gamma.var[3]),0)
							polygon(trans3d(x = min.mean-(max.mean-min.mean)/10,y = c(0,seq.gamma,max.sd),z = pmin(pmax(residuals.gamma,0),max(res.hist$density)*3),pmat = persp.res),col = "#FF000030",border = "darkred",lwd=3)
							arrow.x.red = c(0,0.8,0.8,1,0.8,0.8,0)
							arrow.y.red = c(-0.05,-0.05,-0.1,0,0.1,0.05,0.05)
							max.gamma = seq.gamma[which.min(abs(residuals.gamma-max(residuals.gamma)))]
							i = which.min(abs(UEmat.norm.loess$x-mean(UEmat.norm.mean[reliable])))
							polygon(trans3d(arrow.x.red*(UEmat.norm.loess$x[i]-min.mean)/6*4 + min.mean,(arrow.y.red/5*min(UEmat.gamma.var[2],3)*max.sd + max.gamma),min(max(residuals.gamma),max(res.hist$density))/5*4,pmat = persp.res),col = "darkred",border = "darkred")
							arrow.x.blue = c(0.3,0.8,0.8,1,0.8,0.8,0.3)-0.5
							arrow.y.blue = c(-0.05,-0.05,-0.1,0,0.1,0.05,0.05)
							rect.x.1.blue = c(0.15,0.225,0.225,0.15)-0.5
							rect.x.2.blue = c(0,0.075,0.075,0)-0.5
							rect.y.blue = c(-0.05,-0.05,0.05,0.05)
							polygon(trans3d(arrow.x.blue*0.8*(UEmat.norm.loess$x[i]-min.mean)+min.mean+(UEmat.norm.loess$x[i]-min.mean)/2,(arrow.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
							polygon(trans3d(rect.x.1.blue*0.8*(UEmat.norm.loess$x[i]-min.mean)+min.mean+(UEmat.norm.loess$x[i]-min.mean)/2,(rect.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
							polygon(trans3d(rect.x.2.blue*0.8*(UEmat.norm.loess$x[i]-min.mean)+min.mean+(UEmat.norm.loess$x[i]-min.mean)/2,(rect.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
							polygon(trans3d(x = UEmat.norm.loess$x[i],y = c(0,seq.gamma,max.sd),z = pmin(pmax(c(0,dgamma(x = seq.gamma, shape = (UEmat.norm.loess$fitted[i]^2)/UEmat.gamma.var[1], scale = UEmat.gamma.var[1]/UEmat.norm.loess$fitted[i]),0),0),max(res.hist$density)*3),pmat = persp.res),col = "#0000FF30",border = "darkblue",lwd=3)
							polygon(trans3d((-arrow.x.blue)*0.8*(max.mean-UEmat.norm.loess$x[i])+max.mean-(max.mean-UEmat.norm.loess$x[i])/2,(arrow.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
							polygon(trans3d((-rect.x.1.blue)*0.8*(max.mean-UEmat.norm.loess$x[i])+max.mean-(max.mean-UEmat.norm.loess$x[i])/2,(rect.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
							polygon(trans3d((-rect.x.2.blue)*0.8*(max.mean-UEmat.norm.loess$x[i])+max.mean-(max.mean-UEmat.norm.loess$x[i])/2,(rect.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
							xlab = expression(paste("Mean  ",U[gr]))
							ylab = expression(paste("Standard deviation  ",sigma,"(",U[gr],")"))
							main = expression(paste("Mean  ",U[gr],"  versus standard deviation  ",sigma,"(",U[gr],")"))
							sub = paste("( smoothed loess fit of ",dim(UEmat)[2]," replicates )")
							mtext(main,3,2,cex=2*0.75)
							mtext(sub,3,0,col="darkgrey",cex=1.25*0.75)
							text(trans3d(min.mean + (max.mean-min.mean)*0.6,-max.sd*0.2,-max(res.hist$density)/7,pmat = persp.res),xlab,cex = 1.65*0.75,srt = 90)
							text(trans3d(max.mean + (max.mean-min.mean)/10,max.sd*0.6,-max(res.hist$density),pmat = persp.res),ylab,cex = 1.65*0.75,srt = 0)
							legend("topright",c(expression(paste(Gamma," - distributed std. dev.  ",sigma[g])),expression(paste(Gamma," - distributed residuals  ",r[g]))),pt.bg=c("darkblue","darkred"),col=rep("black",2),bg="white",pch=21,cex=1,pt.cex=1,inset = 0.05)
							par(mar = c(1,4,4,2)+0.1+1)
							par(mai = c(1.1,1.1,1.2,0.7))
							par(cex = 1)
							par(cex.axis = 1.5*0.75)
							heatscatter(UEmat.norm.sd,UEmat.norm.reg.sd,xlab="",ylab="",main="",cor=FALSE,xlim = c(0,max(c(UEmat.norm.sd,UEmat.norm.reg.sd),na.rm = TRUE)),ylim = c(0,max(c(UEmat.norm.sd,UEmat.norm.reg.sd),na.rm = TRUE)))
							abline(0,1,lwd=1)
							xlab = expression(paste("Empirical standard deviation  ",sigma,"(",U[gr],")"))
							ylab = expression(paste("Regularized standard deviation  ",sigma^R,"(",U[gr],")"))
							main = expression(paste("Regularized versus empirical standard deviation"))
							mtext(main,3,3,cex=2*0.75)
							mtext(xlab,1,3,cex=1.65*0.75)
							mtext(ylab,2,3,cex=1.65*0.75)
						}
						DTA.plot.it(filename = paste(folder,"/error_assessment_UE_",condition,"_",timepoint,sep=""),sw = resolution*2,sh = resolution,sres = resolution,plotsfkt = plotsfkt,ww = 14,wh = 7,saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)
					}
					results[["UE.log.error"]] = cbind("log.mean" = UEmat.norm.mean,"log.sd" = UEmat.norm.reg.sd,"log.cv" = UEmat.norm.reg.sd/UEmat.norm.mean)		
				}
			}			
			if (dynamic){
				if (!is.null(initialdummy)){
					cat("Assessing error in total expression (initial timepoint) measurements ... \n")
					TImat.norm = quantnorm(log(TImat))
					TImat.norm.sd = apply(TImat.norm,1,sd)
					TImat.norm.mean = apply(TImat.norm,1,mean)
					TImat.norm.loess = loess(TImat.norm.sd ~ TImat.norm.mean,degree = 1)
					TImat.gamma.var = gamma.variance(TImat.norm.loess$residuals - min(TImat.norm.loess$residuals))
					TImat.feat = cbind(TImat.norm,TImat.norm.mean,TImat.norm.loess$fitted)
					nr.col = ncol(TImat.norm)
					options(warn = -1)
					TImat.norm.reg.sd = apply(TImat.feat,1,function(x){likelihood.fct = function(q){likelihood(x[1:nr.col],x[nr.col + 1],q,x[nr.col + 2],TImat.gamma.var[1])};optimize(likelihood.fct,interval = c(0,1),maximum = TRUE)$maximum})
					options(warn = 0)
					TImat.norm.reg.sd[is.na(TImat.norm.reg.sd)] = TImat.norm.sd[is.na(TImat.norm.reg.sd)]
					if (check) {
						plotsfkt = function(){
							par(mfrow = c(1,2))
							par(mar = c(1,4,4,2)+0.1+1)
							par(cex = 1)
							par(cex.axis = 1.5*0.75)
							min.mean = min(TImat.norm.mean,na.rm=TRUE)
							max.mean = max(TImat.norm.mean,na.rm=TRUE)
							max.sd = max(TImat.norm.sd,na.rm=TRUE)
							x <- seq(min.mean,max.mean,length.out=6)
							y <- seq(0,max.sd,length.out=6)
							z <- matrix(0,6,6)
							residuals = TImat.norm.loess$residuals - min(TImat.norm.loess$residuals)
							res.hist = hist(residuals,plot = FALSE)
							seq.gamma = seq(0,max.sd,length.out = 1000)
							persp(x,y,z,theta = 60,phi = 30,expand = 0.5,ltheta = 120,shade = NA,ticktype = "detailed",xlab = "",ylab = "",zlab = "",zlim = c(0,max(res.hist$density)*5/4),lwd = 0.5,axes = FALSE,box = FALSE,bg = "transparent") -> persp.res
							reliable.range = range(TImat.norm.mean[reliable],na.rm = TRUE)
							polygon(trans3d(c(reliable.range[1],reliable.range[1],reliable.range[2],reliable.range[2]),c(0,max.sd,max.sd,0),0,pmat = persp.res),col = convertcolor("lightslategray",20))
							points(trans3d(TImat.norm.mean,TImat.norm.sd,0,pmat = persp.res),col = "black",pch = ".")
							lines(trans3d(TImat.norm.loess$x[order(TImat.norm.loess$x)],TImat.norm.loess$fitted[order(TImat.norm.loess$x)],0,pmat = persp.res),col = "black",lwd = 2)
							offset = max.sd/25
							lines(trans3d(x,0 - offset,0,pmat = persp.res),lwd = 2)
							text(trans3d(pretty(x)[pretty(x) > min.mean & pretty(x) < max.mean],0 - 2*offset,0,pmat = persp.res),as.character(pretty(x)[pretty(x) > min.mean & pretty(x) < max.mean]))
							lines(trans3d(max.mean + offset*(max.mean-min.mean)/max.sd,y,0,pmat = persp.res),lwd = 2)
							text(trans3d(max.mean + 2*offset*(max.mean-min.mean)/max.sd,pretty(y)[pretty(y) < max.sd],0,pmat = persp.res),as.character(pretty(y)[pretty(y) < max.sd]))						
							for (j in 1:length(res.hist$density)){
								polygon(trans3d(min.mean-(max.mean-min.mean)/10,c(res.hist$breaks[j],res.hist$breaks[j],res.hist$breaks[j+1],res.hist$breaks[j+1]),c(0,res.hist$density[j],res.hist$density[j],0),pmat = persp.res),col = convertcolor("lightgray",60))	
							}
							residuals.gamma = c(0,dgamma(x = seq.gamma, shape = TImat.gamma.var[2], scale = TImat.gamma.var[3]),0)
							polygon(trans3d(x = min.mean-(max.mean-min.mean)/10,y = c(0,seq.gamma,max.sd),z = pmin(pmax(residuals.gamma,0),max(res.hist$density)*3),pmat = persp.res),col = "#FF000030",border = "darkred",lwd=3)
							arrow.x.red = c(0,0.8,0.8,1,0.8,0.8,0)
							arrow.y.red = c(-0.05,-0.05,-0.1,0,0.1,0.05,0.05)
							max.gamma = seq.gamma[which.min(abs(residuals.gamma-max(residuals.gamma)))]
							i = which.min(abs(TImat.norm.loess$x-mean(TImat.norm.mean[reliable])))
							polygon(trans3d(arrow.x.red*(TImat.norm.loess$x[i]-min.mean)/6*4 + min.mean,(arrow.y.red/5*min(TImat.gamma.var[2],3)*max.sd + max.gamma),min(max(residuals.gamma),max(res.hist$density))/5*4,pmat = persp.res),col = "darkred",border = "darkred")
							arrow.x.blue = c(0.3,0.8,0.8,1,0.8,0.8,0.3)-0.5
							arrow.y.blue = c(-0.05,-0.05,-0.1,0,0.1,0.05,0.05)
							rect.x.1.blue = c(0.15,0.225,0.225,0.15)-0.5
							rect.x.2.blue = c(0,0.075,0.075,0)-0.5
							rect.y.blue = c(-0.05,-0.05,0.05,0.05)
							polygon(trans3d(arrow.x.blue*0.8*(TImat.norm.loess$x[i]-min.mean)+min.mean+(TImat.norm.loess$x[i]-min.mean)/2,(arrow.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
							polygon(trans3d(rect.x.1.blue*0.8*(TImat.norm.loess$x[i]-min.mean)+min.mean+(TImat.norm.loess$x[i]-min.mean)/2,(rect.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
							polygon(trans3d(rect.x.2.blue*0.8*(TImat.norm.loess$x[i]-min.mean)+min.mean+(TImat.norm.loess$x[i]-min.mean)/2,(rect.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
							polygon(trans3d(x = TImat.norm.loess$x[i],y = c(0,seq.gamma,max.sd),z = pmin(pmax(c(0,dgamma(x = seq.gamma, shape = (TImat.norm.loess$fitted[i]^2)/TImat.gamma.var[1], scale = TImat.gamma.var[1]/TImat.norm.loess$fitted[i]),0),0),max(res.hist$density)*3),pmat = persp.res),col = "#0000FF30",border = "darkblue",lwd=3)
							polygon(trans3d((-arrow.x.blue)*0.8*(max.mean-TImat.norm.loess$x[i])+max.mean-(max.mean-TImat.norm.loess$x[i])/2,(arrow.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
							polygon(trans3d((-rect.x.1.blue)*0.8*(max.mean-TImat.norm.loess$x[i])+max.mean-(max.mean-TImat.norm.loess$x[i])/2,(rect.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
							polygon(trans3d((-rect.x.2.blue)*0.8*(max.mean-TImat.norm.loess$x[i])+max.mean-(max.mean-TImat.norm.loess$x[i])/2,(rect.y.blue/1.5*max.sd+max.sd/3),0,pmat = persp.res),col = "#0000FF30",border = "darkblue")
							xlab = expression(paste("Mean  ",T[gi]))
							ylab = expression(paste("Standard deviation  ",sigma,"(",T[gi],")"))
							main = expression(paste("Mean  ",T[gi],"  versus standard deviation  ",sigma,"(",T[gi],")"))
							sub = paste("( smoothed loess fit of ",dim(TImat)[2]," replicates )")
							mtext(main,3,2,cex=2*0.75)
							mtext(sub,3,0,col="darkgrey",cex=1.25*0.75)
							text(trans3d(min.mean + (max.mean-min.mean)*0.6,-max.sd*0.2,-max(res.hist$density)/7,pmat = persp.res),xlab,cex = 1.65*0.75,srt = 90)
							text(trans3d(max.mean + (max.mean-min.mean)/10,max.sd*0.6,-max(res.hist$density),pmat = persp.res),ylab,cex = 1.65*0.75,srt = 0)
							legend("topright",c(expression(paste(Gamma," - distributed std. dev.  ",sigma[g])),expression(paste(Gamma," - distributed residuals  ",r[g]))),pt.bg=c("darkblue","darkred"),col=rep("black",2),bg="white",pch=21,cex=1,pt.cex=1,inset = 0.05)
							par(mar = c(1,4,4,2)+0.1+1)
							par(mai = c(1.1,1.1,1.2,0.7))
							par(cex = 1)
							par(cex.axis = 1.5*0.75)
							heatscatter(TImat.norm.sd,TImat.norm.reg.sd,xlab="",ylab="",main="",cor=FALSE,xlim = c(0,max(c(TImat.norm.sd,TImat.norm.reg.sd),na.rm = TRUE)),ylim = c(0,max(c(TImat.norm.sd,TImat.norm.reg.sd),na.rm = TRUE)))
							abline(0,1,lwd=1)
							xlab = expression(paste("Empirical standard deviation  ",sigma,"(",T[gi],")"))
							ylab = expression(paste("Regularized standard deviation  ",sigma^R,"(",T[gi],")"))
							main = expression(paste("Regularized versus empirical standard deviation"))
							mtext(main,3,3,cex=2*0.75)
							mtext(xlab,1,3,cex=1.65*0.75)
							mtext(ylab,2,3,cex=1.65*0.75)
						}
						DTA.plot.it(filename = paste(folder,"/error_assessment_TI_",condition,"_",timepoint,sep=""),sw = resolution*2,sh = resolution,sres = resolution,plotsfkt = plotsfkt,ww = 14,wh = 7,saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)
					}
					results[["TI.log.error"]] = cbind("log.mean" = TImat.norm.mean,"log.sd" = TImat.norm.reg.sd,"log.cv" = TImat.norm.reg.sd/TImat.norm.mean)		
				}
			}
			if (usefractions == "LandT"){
				if (!is.null(initialdummy)){
					cat("Calculating confidence intervals ... \n")
					possible = names(dr[!is.na(dr)])
					possible.mat = cbind(TEmat.norm.mean[possible],TEmat.norm.sd[possible],LEmat.norm.mean[possible],LEmat.norm.sd[possible],TImat.norm.mean[possible],TImat.norm.sd[possible])
					conf.ints = t(apply(possible.mat,1,function(x){TE = rnorm(samplesize,mean=x[1],sd=x[2]);LE = rnorm(samplesize,mean=x[3],sd=x[4]);TI = rnorm(samplesize,mean=x[5],sd=x[6]);options(warn = -1);dr = - alpha - ((1/labelingtime)*log((exp(TE) - exp(LE))/exp(TI)));options(warn = 0);dr[which(dr <= 0)] = NA;hl = log(2)/dr;sr = LE*(alpha + dr)/(exp(alpha*labelingtime)-exp(-dr*labelingtime));TE.quantiles = quantile(exp(TE),confidence.range,na.rm = TRUE);LE.quantiles = quantile(exp(LE),confidence.range,na.rm = TRUE);TI.quantiles = quantile(exp(TI),confidence.range,na.rm = TRUE);LT.quantiles = quantile((exp(TE) - exp(LE))/exp(TI),confidence.range,na.rm = TRUE);LtoT.quantiles = quantile(exp(LE)/exp(TI),confidence.range,na.rm = TRUE);dr.quantiles = quantile(dr,confidence.range,na.rm = TRUE);hl.quantiles = quantile(hl,confidence.range,na.rm = TRUE);sr.quantiles = quantile(sr,confidence.range,na.rm = TRUE);c("TE" = paste(signif(TE.quantiles[1],3),"-",signif(TE.quantiles[2],3)),"LE" = paste(signif(LE.quantiles[1],3),"-",signif(LE.quantiles[2],3)),"TI" = paste(signif(TI.quantiles[1],3),"-",signif(TI.quantiles[2],3)),"LT" = paste(signif(LT.quantiles[1],3),"-",signif(LT.quantiles[2],3)),"LtoT" = paste(signif(LtoT.quantiles[1],3),"-",signif(LtoT.quantiles[2],3)),"dr" = paste(signif(dr.quantiles[1],3),"-",signif(dr.quantiles[2],3)),"hl" = paste(signif(hl.quantiles[1],3),"-",signif(hl.quantiles[2],3)),"sr" = paste(signif(sr.quantiles[1],3),"-",signif(sr.quantiles[2],3)))}))
					results[["TE.confidence"]] = conf.ints[,"TE"]
					results[["LE.confidence"]] = conf.ints[,"LE"]
					results[["TI.confidence"]] = conf.ints[,"TI"]
					results[["LT.confidence"]] = conf.ints[,"LT"]
					results[["LtoT.confidence"]] = conf.ints[,"LtoT"]
					results[["dr.confidence"]] = conf.ints[,"dr"]
					results[["hl.confidence"]] = conf.ints[,"hl"]
					results[["sr.confidence"]] = conf.ints[,"sr"]
				}  else {
					cat("Calculating confidence intervals ... \n")
					possible = names(dr[!is.na(dr)])
					possible.mat = cbind(TEmat.norm.mean[possible],TEmat.norm.sd[possible],LEmat.norm.mean[possible],LEmat.norm.sd[possible])
					conf.ints = t(apply(possible.mat,1,function(x){TE = rnorm(samplesize,mean=x[1],sd=x[2]);LE = rnorm(samplesize,mean=x[3],sd=x[4]);options(warn = -1);dr = - alpha - ((1/labelingtime)*log(1 - exp(LE)/exp(TE)));options(warn = 0);dr[which(dr <= 0)] = NA;hl = log(2)/dr;sr = LE*(alpha + dr)/(exp(alpha*labelingtime)-exp(-dr*labelingtime));TE.quantiles = quantile(exp(TE),confidence.range,na.rm = TRUE);LE.quantiles = quantile(exp(LE),confidence.range,na.rm = TRUE);LT.quantiles = quantile(1 - exp(LE)/exp(TE),confidence.range,na.rm = TRUE);LtoT.quantiles = quantile(exp(LE)/exp(TE),confidence.range,na.rm = TRUE);dr.quantiles = quantile(dr,confidence.range,na.rm = TRUE);hl.quantiles = quantile(hl,confidence.range,na.rm = TRUE);sr.quantiles = quantile(sr,confidence.range,na.rm = TRUE);c("TE" = paste(signif(TE.quantiles[1],3),"-",signif(TE.quantiles[2],3)),"LE" = paste(signif(LE.quantiles[1],3),"-",signif(LE.quantiles[2],3)),"LT" = paste(signif(LT.quantiles[1],3),"-",signif(LT.quantiles[2],3)),"LtoT" = paste(signif(LtoT.quantiles[1],3),"-",signif(LtoT.quantiles[2],3)),"dr" = paste(signif(dr.quantiles[1],3),"-",signif(dr.quantiles[2],3)),"hl" = paste(signif(hl.quantiles[1],3),"-",signif(hl.quantiles[2],3)),"sr" = paste(signif(sr.quantiles[1],3),"-",signif(sr.quantiles[2],3)))}))
					results[["TE.confidence"]] = conf.ints[,"TE"]
					results[["LE.confidence"]] = conf.ints[,"LE"]
					results[["LT.confidence"]] = conf.ints[,"LT"]
					results[["LtoT.confidence"]] = conf.ints[,"LtoT"]
					results[["dr.confidence"]] = conf.ints[,"dr"]
					results[["hl.confidence"]] = conf.ints[,"hl"]
					results[["sr.confidence"]] = conf.ints[,"sr"]				
				}			
			}
			if (usefractions == "UandT"){
				if (!is.null(initialdummy)){
					cat("Calculating confidence intervals ... \n")
					possible = names(dr[!is.na(dr)])
					possible.mat = cbind(TEmat.norm.mean[possible],TEmat.norm.sd[possible],LEmat.norm.mean[possible],LEmat.norm.sd[possible],UEmat.norm.mean[possible],UEmat.norm.sd[possible],TImat.norm.mean[possible],TImat.norm.sd[possible])
					conf.ints = t(apply(possible.mat,1,function(x){TE = rnorm(samplesize,mean=x[1],sd=x[2]);LE = rnorm(samplesize,mean=x[3],sd=x[4]);UE = rnorm(samplesize,mean=x[5],sd=x[6]);TI = rnorm(samplesize,mean=x[7],sd=x[8]);options(warn = -1);dr = - alpha - ((1/labelingtime)*log(exp(UE)/exp(TI)));options(warn = 0);dr[which(dr <= 0)] = NA;hl = log(2)/dr;sr = LE*(alpha + dr)/(exp(alpha*labelingtime)-exp(-dr*labelingtime));TE.quantiles = quantile(exp(TE),confidence.range,na.rm = TRUE);LE.quantiles = quantile(exp(LE),confidence.range,na.rm = TRUE);UE.quantiles = quantile(exp(UE),confidence.range,na.rm = TRUE);TI.quantiles = quantile(exp(TI),confidence.range,na.rm = TRUE);UT.quantiles = quantile(exp(UE)/exp(TI),confidence.range,na.rm = TRUE);UtoT.quantiles = quantile(exp(UE)/exp(TI),confidence.range,na.rm = TRUE);dr.quantiles = quantile(dr,confidence.range,na.rm = TRUE);hl.quantiles = quantile(hl,confidence.range,na.rm = TRUE);sr.quantiles = quantile(sr,confidence.range,na.rm = TRUE);c("TE" = paste(signif(TE.quantiles[1],3),"-",signif(TE.quantiles[2],3)),"LE" = paste(signif(LE.quantiles[1],3),"-",signif(LE.quantiles[2],3)),"UE" = paste(signif(UE.quantiles[1],3),"-",signif(UE.quantiles[2],3)),"TI" = paste(signif(TI.quantiles[1],3),"-",signif(TI.quantiles[2],3)),"UT" = paste(signif(UT.quantiles[1],3),"-",signif(UT.quantiles[2],3)),"UtoT" = paste(signif(UtoT.quantiles[1],3),"-",signif(UtoT.quantiles[2],3)),"dr" = paste(signif(dr.quantiles[1],3),"-",signif(dr.quantiles[2],3)),"hl" = paste(signif(hl.quantiles[1],3),"-",signif(hl.quantiles[2],3)),"sr" = paste(signif(sr.quantiles[1],3),"-",signif(sr.quantiles[2],3)))}))
					results[["TE.confidence"]] = conf.ints[,"TE"]
					results[["LE.confidence"]] = conf.ints[,"LE"]
					results[["UE.confidence"]] = conf.ints[,"UE"]
					results[["TI.confidence"]] = conf.ints[,"TI"]
					results[["UT.confidence"]] = conf.ints[,"UT"]
					results[["UtoT.confidence"]] = conf.ints[,"UtoT"]
					results[["dr.confidence"]] = conf.ints[,"dr"]
					results[["hl.confidence"]] = conf.ints[,"hl"]
					results[["sr.confidence"]] = conf.ints[,"sr"]				
				}  else {
					cat("Calculating confidence intervals ... \n")
					possible = names(dr[!is.na(dr)])
					possible.mat = cbind(TEmat.norm.mean[possible],TEmat.norm.sd[possible],LEmat.norm.mean[possible],LEmat.norm.sd[possible],UEmat.norm.mean[possible],UEmat.norm.sd[possible])
					conf.ints = t(apply(possible.mat,1,function(x){TE = rnorm(samplesize,mean=x[1],sd=x[2]);LE = rnorm(samplesize,mean=x[3],sd=x[4]);UE = rnorm(samplesize,mean=x[5],sd=x[6]);options(warn = -1);dr = - alpha - ((1/labelingtime)*log(exp(UE)/exp(TE)));options(warn = 0);dr[which(dr <= 0)] = NA;hl = log(2)/dr;sr = LE*(alpha + dr)/(exp(alpha*labelingtime)-exp(-dr*labelingtime));TE.quantiles = quantile(exp(TE),confidence.range,na.rm = TRUE);LE.quantiles = quantile(exp(LE),confidence.range,na.rm = TRUE);UE.quantiles = quantile(exp(UE),confidence.range,na.rm = TRUE);UT.quantiles = quantile(exp(UE)/exp(TE),confidence.range,na.rm = TRUE);UtoT.quantiles = quantile(exp(UE)/exp(TE),confidence.range,na.rm = TRUE);dr.quantiles = quantile(dr,confidence.range,na.rm = TRUE);hl.quantiles = quantile(hl,confidence.range,na.rm = TRUE);sr.quantiles = quantile(sr,confidence.range,na.rm = TRUE);c("TE" = paste(signif(TE.quantiles[1],3),"-",signif(TE.quantiles[2],3)),"LE" = paste(signif(LE.quantiles[1],3),"-",signif(LE.quantiles[2],3)),"UE" = paste(signif(UE.quantiles[1],3),"-",signif(UE.quantiles[2],3)),"UT" = paste(signif(UT.quantiles[1],3),"-",signif(UT.quantiles[2],3)),"UtoT" = paste(signif(UtoT.quantiles[1],3),"-",signif(UtoT.quantiles[2],3)),"dr" = paste(signif(dr.quantiles[1],3),"-",signif(dr.quantiles[2],3)),"hl" = paste(signif(hl.quantiles[1],3),"-",signif(hl.quantiles[2],3)),"sr" = paste(signif(sr.quantiles[1],3),"-",signif(sr.quantiles[2],3)))}))
					results[["TE.confidence"]] = conf.ints[,"TE"]
					results[["LE.confidence"]] = conf.ints[,"LE"]
					results[["UE.confidence"]] = conf.ints[,"UE"]
					results[["UT.confidence"]] = conf.ints[,"UT"]
					results[["UtoT.confidence"]] = conf.ints[,"UtoT"]
					results[["dr.confidence"]] = conf.ints[,"dr"]
					results[["hl.confidence"]] = conf.ints[,"hl"]
					results[["sr.confidence"]] = conf.ints[,"sr"]				
				}			
			}
			if (usefractions == "both"){
				if (!is.null(initialdummy)){
					cat("Calculating confidence intervals ... \n")
					possible = names(dr[!is.na(dr)])
					possible.mat = cbind(TEmat.norm.mean[possible],TEmat.norm.sd[possible],LEmat.norm.mean[possible],LEmat.norm.sd[possible],UEmat.norm.mean[possible],UEmat.norm.sd[possible],TImat.norm.mean[possible],TImat.norm.sd[possible])
					conf.ints = t(apply(possible.mat,1,function(x){TE = rnorm(samplesize,mean=x[1],sd=x[2]);LE = rnorm(samplesize,mean=x[3],sd=x[4]);UE = rnorm(samplesize,mean=x[5],sd=x[6]);TI = rnorm(samplesize,mean=x[7],sd=x[8]);options(warn = -1);dr = - alpha - ((1/labelingtime)*log(((exp(TE) - exp(LE))/exp(TI) + exp(UE)/exp(TI))/2));options(warn = 0);dr[which(dr <= 0)] = NA;hl = log(2)/dr;sr = LE*(alpha + dr)/(exp(alpha*labelingtime)-exp(-dr*labelingtime));TE.quantiles = quantile(exp(TE),confidence.range,na.rm = TRUE);LE.quantiles = quantile(exp(LE),confidence.range,na.rm = TRUE);UE.quantiles = quantile(exp(UE),confidence.range,na.rm = TRUE);TI.quantiles = quantile(exp(TI),confidence.range,na.rm = TRUE);LUT.quantiles = quantile(((exp(TE) - exp(LE))/exp(TI) + exp(UE)/exp(TI))/2,confidence.range,na.rm = TRUE);LandUtoT.quantiles = quantile((exp(LE)/exp(TI) + exp(UE)/exp(TI))/2,confidence.range,na.rm = TRUE);dr.quantiles = quantile(dr,confidence.range,na.rm = TRUE);hl.quantiles = quantile(hl,confidence.range,na.rm = TRUE);sr.quantiles = quantile(sr,confidence.range,na.rm = TRUE);c("TE" = paste(signif(TE.quantiles[1],3),"-",signif(TE.quantiles[2],3)),"LE" = paste(signif(LE.quantiles[1],3),"-",signif(LE.quantiles[2],3)),"UE" = paste(signif(UE.quantiles[1],3),"-",signif(UE.quantiles[2],3)),"TI" = paste(signif(TI.quantiles[1],3),"-",signif(TI.quantiles[2],3)),"LUT" = paste(signif(LUT.quantiles[1],3),"-",signif(LUT.quantiles[2],3)),"LandUtoT" = paste(signif(LandUtoT.quantiles[1],3),"-",signif(LandUtoT.quantiles[2],3)),"dr" = paste(signif(dr.quantiles[1],3),"-",signif(dr.quantiles[2],3)),"hl" = paste(signif(hl.quantiles[1],3),"-",signif(hl.quantiles[2],3)),"sr" = paste(signif(sr.quantiles[1],3),"-",signif(sr.quantiles[2],3)))}))
					results[["TE.confidence"]] = conf.ints[,"TE"]
					results[["LE.confidence"]] = conf.ints[,"LE"]
					results[["UE.confidence"]] = conf.ints[,"UE"]
					results[["TI.confidence"]] = conf.ints[,"TI"]
					results[["LUT.confidence"]] = conf.ints[,"LUT"]
					results[["LandUtoT.confidence"]] = conf.ints[,"LandUtoT"]
					results[["dr.confidence"]] = conf.ints[,"dr"]
					results[["hl.confidence"]] = conf.ints[,"hl"]
					results[["sr.confidence"]] = conf.ints[,"sr"]
				}  else {
					cat("Calculating confidence intervals ... \n")
					possible = names(dr[!is.na(dr)])
					possible.mat = cbind(TEmat.norm.mean[possible],TEmat.norm.sd[possible],LEmat.norm.mean[possible],LEmat.norm.sd[possible],UEmat.norm.mean[possible],UEmat.norm.sd[possible])
					conf.ints = t(apply(possible.mat,1,function(x){TE = rnorm(samplesize,mean=x[1],sd=x[2]);LE = rnorm(samplesize,mean=x[3],sd=x[4]);UE = rnorm(samplesize,mean=x[5],sd=x[6]);options(warn = -1);dr = - alpha - ((1/labelingtime)*log((1 - exp(LE)/exp(TE) + exp(UE)/exp(TE))/2));options(warn = 0);dr[which(dr <= 0)] = NA;hl = log(2)/dr;sr = LE*(alpha + dr)/(exp(alpha*labelingtime)-exp(-dr*labelingtime));TE.quantiles = quantile(exp(TE),confidence.range,na.rm = TRUE);LE.quantiles = quantile(exp(LE),confidence.range,na.rm = TRUE);UE.quantiles = quantile(exp(UE),confidence.range,na.rm = TRUE);LUT.quantiles = quantile((1 - exp(LE)/exp(TE) + exp(UE)/exp(TE))/2,confidence.range,na.rm = TRUE);LandUtoT.quantiles = quantile((exp(LE)/exp(TE) + exp(UE)/exp(TE))/2,confidence.range,na.rm = TRUE);dr.quantiles = quantile(dr,confidence.range,na.rm = TRUE);hl.quantiles = quantile(hl,confidence.range,na.rm = TRUE);sr.quantiles = quantile(sr,confidence.range,na.rm = TRUE);c("TE" = paste(signif(TE.quantiles[1],3),"-",signif(TE.quantiles[2],3)),"LE" = paste(signif(LE.quantiles[1],3),"-",signif(LE.quantiles[2],3)),"UE" = paste(signif(UE.quantiles[1],3),"-",signif(UE.quantiles[2],3)),"LUT" = paste(signif(LUT.quantiles[1],3),"-",signif(LUT.quantiles[2],3)),"LandUtoT" = paste(signif(LandUtoT.quantiles[1],3),"-",signif(LandUtoT.quantiles[2],3)),"dr" = paste(signif(dr.quantiles[1],3),"-",signif(dr.quantiles[2],3)),"hl" = paste(signif(hl.quantiles[1],3),"-",signif(hl.quantiles[2],3)),"sr" = paste(signif(sr.quantiles[1],3),"-",signif(sr.quantiles[2],3)))}))
					results[["TE.confidence"]] = conf.ints[,"TE"]
					results[["LE.confidence"]] = conf.ints[,"LE"]
					results[["UE.confidence"]] = conf.ints[,"UE"]
					results[["LUT.confidence"]] = conf.ints[,"LUT"]
					results[["LandUtoT.confidence"]] = conf.ints[,"LandUtoT"]
					results[["dr.confidence"]] = conf.ints[,"dr"]
					results[["hl.confidence"]] = conf.ints[,"hl"]
					results[["sr.confidence"]] = conf.ints[,"sr"]
				}			
			}
		} else {cat("Your data is not suitable for error assessment due to lack of replicates \n")}				
	}
	
	### SIMULATION ###
	
	if (simulation){	
		plotsfkt = function(){
			par(mfrow=c(3,3))
			#parfkt("default",9)
			if (dynamic){
				truehl = truehalflivesaveraged
				truedr = truelambdasaveraged
				truesr = truemusaveraged
			} else {
				truehl = truehalflives
				truedr = truelambdas
				truesr = truemus
			}
			
			hllims = c(0,mean(quantile(truehl,0.99,na.rm=TRUE),quantile(results[["hl"]],0.99,na.rm=TRUE)))
			hlindi = (!is.na(results[["hl"]]) & results[["hl"]] > 0 & results[["hl"]] < hllims[2])
			drlims = c(mean(quantile(truedr,0.01,na.rm=TRUE),quantile(results[["dr"]],0.01,na.rm=TRUE)),mean(quantile(truedr,0.99,na.rm=TRUE),quantile(results[["dr"]],0.99,na.rm=TRUE)))
			drindi = (!is.na(results[["dr"]]) & results[["dr"]] > 0 & results[["dr"]] < drlims[2])
			srlims = c(0,mean(quantile(truesr,0.99,na.rm=TRUE),quantile(results[["sr"]],0.99,na.rm=TRUE)))
			srindi = (!is.na(results[["sr"]]) & results[["sr"]] > 0 & results[["sr"]] < srlims[2])
			
			parfkt("default",9)
			heatscatter(truehl,results[["hl"]],main="",xlab="",ylab="",cor=FALSE,xlim=hllims,ylim=hllims)
			xlab = expression(paste("True half-life  ",t[1/2][gr]))
			ylab = expression(paste("Estimated half-life  ",t[1/2][gr]^'*'))
			main = expression(paste("Half-life  ",t[1/2][gr]))
			sub = paste("( heatscatter p-cor = ",pcor(truehl,results[["hl"]],use="na.or.complete"),")")
			mtextfkt("default",9,main,xlab,ylab,sub)
			linfactor = coefficients(tls((results[["hl"]][hlindi]) ~ (truehl[hlindi]) + 0))[1]
			abline(a=0,b=linfactor,col="black",lwd=3)
			abline(a=0,b=1,col="grey",lwd=3,lty=2)
			cat("Factor estimated/true:",linfactor,"\n")
			heatscatter(rank(truehl),rank(results[["hl"]]),main="",xlab="",ylab="",cor=FALSE)
			xlab = expression(paste("True half-life  ",t[1/2][gr],"  rank"))
			ylab = expression(paste("Estimated half-life  ",t[1/2][gr]^'*',"  rank"))
			main = expression(paste("Half-life  ",t[1/2][gr],"  rank"))
			sub = paste("( heatscatter s-cor = ",scor(rank(truehl),rank(results[["hl"]]),use="na.or.complete"),")")
			mtextfkt("default",9,main,xlab,ylab,sub)
			abline(a=0,b=1,col="black",lwd=3)
			abline(a=0,b=1,col="grey",lwd=3,lty=2)
			meanreldev = mean(abs(results[["hl"]][hlindi]-truehl[hlindi])/truehl[hlindi])
			cat("Mean relative deviation: ",signif(meanreldev,2),"\n")
			den = density(log2(results[["hl"]][hlindi]/truehl[hlindi]))
			den = cbind(den$x,den$y)
			hist(log2(results[["hl"]][hlindi]/truehl[hlindi]),xlim=c(-5,5),breaks=40,main="",xlab="",ylab="")
			xlab = expression(paste("log2( ",t[1/2][gr]^'*'/t[1/2][gr]," )"))
			ylab = expression(paste("Frequency"))
			main = expression(paste("Log-ratio  log2( ",t[1/2][gr]^'*'/t[1/2][gr]," )"))
			sub = paste("( MRD: ",signif(meanreldev,2), ",mode: ",signif(den[which(den[,2] == max(den[,2]))],2)," )")
			mtextfkt("default",9,main,xlab,ylab,sub)
			abline(v=0,col="grey",lwd=3,lty=2)
			abline(v=den[which(den[,2] == max(den[,2]))],col="black",lwd=3)
			
			heatscatter(truedr,results[["dr"]],log="xy",main="",xlab="",ylab="",cor=FALSE,xlim=drlims,ylim=drlims)
			xlab = expression(paste("True decay rate  ",lambda[gr]))
			ylab = expression(paste("Estimated decay rate  ",lambda[gr]^'*'))
			main = expression(paste("Deacy rate  ",lambda[gr]))
			sub = paste("( heatscatter p-cor = ",pcor(truedr,results[["dr"]],use="na.or.complete"),")")
			mtextfkt("default",9,main,xlab,ylab,sub)
			lindev = median(log10(results[["dr"]][drindi])-log10(truedr[drindi]))
			abline(lindev,b=1,col="black",lwd=3)
			abline(a=0,b=1,col="grey",lwd=3,lty=2)
			heatscatter(rank(truedr),rank(results[["dr"]]),main="",xlab="",ylab="",cor=FALSE)
			xlab = expression(paste("True decay rate  ",lambda[gr],"  rank"))
			ylab = expression(paste("Estimated decay  ",lambda[gr]^'*',"  rank"))
			main = expression(paste("Decay rate  ",lambda[gr],"  rank"))
			sub = paste("( heatscatter s-cor = ",scor(rank(truedr),rank(results[["dr"]]),use="na.or.complete"),")")
			mtextfkt("default",9,main,xlab,ylab,sub)
			abline(a=0,b=1,col="black",lwd=3)
			abline(a=0,b=1,col="grey",lwd=3,lty=2)
			meanreldev = mean(abs(results[["dr"]][drindi]-truedr[drindi])/truedr[drindi])
			cat("Mean relative deviation: ",signif(meanreldev,2),"\n")
			den = density(log2(results[["dr"]][drindi]/truedr[drindi]))
			den = cbind(den$x,den$y)
			hist(log2(results[["dr"]][drindi]/truedr[drindi]),xlim=c(-5,5),breaks=40,main="",xlab="",ylab="")
			xlab = expression(paste("log2( ",lambda[gr]^'*'/lambda[gr]," )"))
			ylab = expression(paste("Frequency"))
			main = expression(paste("Log-ratio  log2( ",lambda[gr]^'*'/lambda[gr]," )"))
			sub = paste("( MRD: ",signif(meanreldev,2), ",mode: ",signif(den[which(den[,2] == max(den[,2]))],2)," )")
			mtextfkt("default",9,main,xlab,ylab,sub)
			abline(v=0,col="grey",lwd=3,lty=2)
			abline(v=den[which(den[,2] == max(den[,2]))],col="black",lwd=3)
			
			heatscatter(truesr,results[["sr"]],main="",xlab="",ylab="",cor=FALSE,xlim=srlims,ylim=srlims)
			xlab = expression(paste("True synthesis rate  ",mu[gr]))
			ylab = expression(paste("Estimated synthesis rate  ",mu[gr]^'*'))
			main = expression(paste("Synthesis rate  ",mu[gr]))
			sub = paste("( heatscatter p-cor = ",pcor(truesr,results[["sr"]],use="na.or.complete"),")")
			mtextfkt("default",9,main,xlab,ylab,sub)
			linfactor = coefficients(tls((results[["sr"]][srindi]) ~ (truesr[srindi]) + 0))[1]
			abline(a=0,b=linfactor,col="black",lwd=3)
			abline(a=0,b=1,col="grey",lwd=3,lty=2)
			heatscatter(rank(truesr),rank(results[["sr"]]),main="",xlab="",ylab="",cor=FALSE)
			xlab = expression(paste("True synthesis rate  ",mu[gr],"  rank"))
			ylab = expression(paste("Estimated synthesis  ",mu[gr]^'*',"  rank"))
			main = expression(paste("Synthesis rate  ",mu[gr],"  rank"))
			sub = paste("( heatscatter s-cor = ",scor(rank(truesr),rank(results[["sr"]]),use="na.or.complete"),")")
			mtextfkt("default",9,main,xlab,ylab,sub)
			abline(a=0,b=1,col="black",lwd=3)
			abline(a=0,b=1,col="grey",lwd=3,lty=2)
			meanreldev = mean(abs(results[["sr"]][srindi]-truesr[srindi])/truesr[srindi])
			cat("Mean relative deviation: ",signif(meanreldev,2),"\n")
			den = density(log2(results[["sr"]][srindi]/truesr[srindi]))
			den = cbind(den$x,den$y)
			hist(log2(results[["sr"]][srindi]/truesr[srindi]),xlim=c(-5,5),breaks=40,main="",xlab="",ylab="")
			xlab = expression(paste("log2( ",mu[gr]^'*'/mu[gr]," )"))
			ylab = expression(paste("Frequency"))
			main = expression(paste("Log-ratio  log2( ",mu[gr]^'*'/mu[gr]," )"))
			sub = paste("( MRD: ",signif(meanreldev,2), ",mode: ",signif(den[which(den[,2] == max(den[,2]))],2)," )")
			mtextfkt("default",9,main,xlab,ylab,sub)
			abline(v=0,col="grey",lwd=3,lty=2)
			abline(v=den[which(den[,2] == max(den[,2]))],col="black",lwd=3)
		}
		DTA.plot.it(filename = paste(folder,"/simulation_",condition,"_",timepoint,sep=""),sw = resolution*3,sh = resolution*3,sres = resolution,plotsfkt = plotsfkt,ww = 21,wh = 21,saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)
	}
	
	### RETURN RESULTS IN A LIST ###
	
	return(results)
}


### estimate uses an experiment, given by a phenotype matrix, data matrix and the #uridines for each gene to estimate synthesis and decay rate of the genes ###


DTA.estimate = function(phenomat = NULL,	# phenotype matrix, "nr" should be numbered by experiments not by replicates
		datamat = NULL, 					# data matrix, should only contain the rows of phenomat as columns
		tnumber = NULL, 					# #uridines, should have the rownames of datamat
		reliable = NULL, 					# vector of reliable genes, which are used for regression
		ccl = NULL, 						# the cell cycle length of the cells
		mRNAs = NULL,						# estimated number of mRNAs in a cell
		mediancenter = TRUE, 				# should the L/T resp. U/T ratio of replicates be rescaled to a common median before rate extraction
		usefractions = "LandT",				# from which fractions should the decay rate be calculated: "LandT", "UandT" or "both"
		LtoTratio = NULL, 					# coefficient to rescale L/T
		ratiomethod = "tls", 				# choose the regression method to be used, possible methods are: tls, bias, lm
		largest = 5, 						# percentage of largest residues from the first regression not to be used in the second regression
		weighted = TRUE, 					# should the regression be weighted with 1/(T^2 + median(T))
		relevant = NULL, 					# choose the arrays to be used for halflives calculation, vector due to experiments variable 
		upper = 700, 						# upper bound for labeling bias estimation
		lower = 500, 						# lower bound for labeling bias estimation
		error = TRUE,						# should the measurement error be assessed
		samplesize = 1000,					# error model samplesize for resampling
		confidence.range = c(0.025,0.975),	# confidence region for error model
		bicor = TRUE, 						# should the labeling bias be corrected		
		check = TRUE, 						# if check = TRUE, control messages and plots will be generated
		condition = "", 					# to be added to the plotnames
		save.plots = FALSE, 				# if save.plots = TRUE, control plots will be saved
		resolution = 1,						# resolution scaling factor for plotting
		notinR = FALSE,						# should plot be not plotted in R
		RStudio = FALSE,					# for RStudio users
		folder = NULL, 						# folder, where to save the plots
		fileformat = "jpeg", 				# save the plot as jpeg, png, bmp, tiff, ps or pdf
		totaloverwt = 1, 					# total mRNA over WT
		simulation = FALSE,					# simulated data via sim.object ?
		sim.object = NULL					# simulation object created by DTA.generate
)
{
	
	### CHECK REQUIREMENTS ###
	
	if (!simulation){
		if (is.null(phenomat)){stop("You need to specify the Phenotype matrix (phenomat) !")}
		if (is.null(datamat)){stop("You need to specify the Data matrix (datamat) !")}
		if (is.null(tnumber)){stop("You need to specify the number of uridine residues per identifier (tnumber) !")}
	} else {if (is.null(sim.object)){stop("You need to specify the Simulation object (sim.object) created by DTA.generate !")}}
	if (simulation){
		phenomat = sim.object$phenomat
		datamat = sim.object$datamat
		tnumber = sim.object$tnumber
		ccl = sim.object$ccl
		truemus = sim.object$truemus
		truelambdas = sim.object$truelambdas
		truehalflives = sim.object$truehalflives
		trueplabel = sim.object$trueplabel
		names(trueplabel) = phenomat[which(phenomat[,"fraction"]=="T"),"nr"]
		trueLasymptote = sim.object$trueLasymptote
		names(trueLasymptote) = phenomat[which(phenomat[,"fraction"]=="T"),"nr"]
		trueUasymptote = sim.object$trueUasymptote
		names(trueUasymptote) = phenomat[which(phenomat[,"fraction"]=="T"),"nr"]
		truecrbyar = sim.object$truecrbyar
		names(truecrbyar) = phenomat[which(phenomat[,"fraction"]=="T"),"nr"]
		truecrbybr = sim.object$truecrbybr
		names(truecrbybr) = phenomat[which(phenomat[,"fraction"]=="T"),"nr"]
		truebrbyar = sim.object$truebrbyar
		names(truebrbyar) = phenomat[which(phenomat[,"fraction"]=="T"),"nr"]
	}
	if (!all(rownames(phenomat) %in% colnames(datamat))){stop("The rownames of the phenomat should be among the colnames of the datamat !")}
	if (length(intersect(rownames(datamat),names(tnumber))) == 0){stop("The rownames of the datamat should be the same identifiers as the names of tnumber !")}
	if (is.null(ccl)){cat("If you do not specify the Cell Cycle Length (ccl), growth is set to zero or as specified in your sim.object ! \n")}
	if (is.null(reliable)){cat("If you do not specify a vector of reliable identifiers (reliable), the parameter estimation is done on all identifiers ! \n")
		reliable = rownames(datamat)}
	if (!any(usefractions %in% c("LandT","UandT","both"))){stop("usefractions need to be 'LandT', 'UandT' or 'both' !")}
	if (!any(ratiomethod %in% c("tls","lm","bias"))){stop("ratiomethod need to be 'bias', 'tls' or 'lm' !")}
	if (is.null(LtoTratio)){cat("If you do not specify ratio of L to T (LtoTratio), it is estimated from the data ! \n")}
	if (!is.null(LtoTratio)){
		ratiomethod = "tls"
		if (length(unique(phenomat[,"time"])) != length(LtoTratio)){stop(paste("The number of specified ratios of L to T (LtoTratio) should correspond to the number of labeling durations:",length(unique(phenomat[,"time"])),"!"))}
	}
	if (save.plots){
		if (is.null(folder)){stop("You need to specify the folder, where the plots should be saved !")}
		if (!is.null(folder)){
			if (file.access(folder,0) != 0){stop("The specified folder needs to exist !")}
			if (file.access(folder,2) != 0){stop("The specified folder has no write permission !")}
		}
	}
	
	### PRELIMINARIES ###
	
	possibles = intersect(rownames(datamat),names(tnumber))
	tnumber = tnumber[possibles]
	datamat = datamat[possibles,]
	reliable = intersect(reliable,possibles)
	
	labtimes = unique(phenomat[,"time"])
	names(labtimes) = labtimes
	if (!is.null(LtoTratio)){names(LtoTratio) = labtimes}
	nrlabtimes = length(labtimes)
	
	### BUILD PHENOMATS ###
	
	phenomats = list()
	for (labtime in labtimes){
		phenomats[[labtime]] = phenomat[names(which(phenomat[,"time"] == labtime)),]
	}
	
	### BUILD DATAMATS ###
	
	datamats = list()
	for (labtime in labtimes){
		datamats[[labtime]] = datamat[,rownames(phenomats[[labtime]])]
	}
	
	### BUILD RELEVANTS ###
	
	relevants = list()
	if (is.null(relevant)){
		for (labtime in labtimes){relevants[[labtime]] = NULL}
	}
	else for (labtime in labtimes){
			relevants[[labtime]] = intersect(relevant,unique(phenomats[[labtime]][,"nr"]))
		}
	
	### BUILD LtoT RATIOS ###
	
	ratios = list()
	if (is.null(LtoTratio)){
		for (labtime in labtimes){ratios[[labtime]] = NULL}
	}
	else for (labtime in labtimes){
			if (LtoTratio[labtime] == 0){ratios[[labtime]] = NULL} else {ratios[[labtime]] = rep(LtoTratio[labtime],length(unique(phenomats[[labtime]][,"nr"])))}
		}
	
	### BUILD TRUE VALUES FROM SIMULATION ###
	
	trueplabels = list()
	trueLasymptotes = list()
	trueUasymptotes = list()
	truecrbyars = list()
	truecrbybrs = list()
	truebrbyars = list()
	if (simulation){
		for (labtime in labtimes){trueplabels[[labtime]] = trueplabel[phenomats[[labtime]][which(phenomats[[labtime]][,"fraction"]=="T"),"nr"]]}
		for (labtime in labtimes){trueLasymptotes[[labtime]] = trueLasymptote[phenomats[[labtime]][which(phenomats[[labtime]][,"fraction"]=="T"),"nr"]]}
		for (labtime in labtimes){trueUasymptotes[[labtime]] = trueUasymptote[phenomats[[labtime]][which(phenomats[[labtime]][,"fraction"]=="T"),"nr"]]}
		for (labtime in labtimes){truecrbyars[[labtime]] = truecrbyar[phenomats[[labtime]][which(phenomats[[labtime]][,"fraction"]=="T"),"nr"]]}
		for (labtime in labtimes){truecrbybrs[[labtime]] = truecrbybr[phenomats[[labtime]][which(phenomats[[labtime]][,"fraction"]=="T"),"nr"]]}
		for (labtime in labtimes){truebrbyars[[labtime]] = truebrbyar[phenomats[[labtime]][which(phenomats[[labtime]][,"fraction"]=="T"),"nr"]]}
	}
	
	### SINGLEESTIMATES ###
	
	res = list()
	for (labtime in labtimes){
		res[[labtime]] = DTA.singleestimate(phenomats[[labtime]],datamats[[labtime]],tnumber,labelingtime = as.numeric(labtimes[labtime]),ccl = ccl,mRNAs = mRNAs,
				reliable = reliable,mediancenter = mediancenter,usefractions = usefractions,ratiomethod = ratiomethod,resolution = resolution,RStudio = RStudio,
				largest = largest,weighted = weighted,ratio = ratios[[labtime]],relevant = relevants[[labtime]],check = check,totaloverwt = totaloverwt,
				error = error,samplesize = samplesize,confidence.range = confidence.range,bicor = bicor,condition = condition,timepoint = labtime,upper = upper,lower = lower,save.plots = save.plots,notinR = notinR,folder = folder,
				fileformat = fileformat,simulation = simulation,truemus = truemus,truelambdas = truelambdas,truehalflives = truehalflives,trueplabel=trueplabels[[labtime]],
				trueLasymptote=trueLasymptotes[[labtime]],trueUasymptote=trueUasymptotes[[labtime]],truecrbyar=truecrbyars[[labtime]],
				truecrbybr=truecrbybrs[[labtime]],truebrbyar=truebrbyars[[labtime]])
	}
	
	### RETURN RESULTS IN A LIST ###
	
	return(res)
}



