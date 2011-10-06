

### singleestimate uses an experiment, given by a phenotype matrix, data matrix and the #uridines for each gene to estimate synthesis and decay rate of the genes ###


DTA.singleestimate = function(phenomat, # phenotype matrix, "nr" should be numbered by experiments not by replicates
		datamat, 						# data matrix, should only contain the rows of phenomat as columns
		tnumber, 						# #uridines, should have the rownames of datamat
		labelingtime,					# length of the labeling period
		ccl = NULL, 					# the cell cycle length of the cells
		mRNAs = NULL,					# estimated number of mRNAs in a cell
		reliable = NULL, 				# vector of reliable genes, which are used for regression
		mediancenter = TRUE, 			# should the L/T resp. U/T ratio of replicates be rescaled to a common median before rate extraction
		ratiomethod = "tls", 			# choose the regression method to be used, possible methods are: tls, bias, lm
		largest = 5, 					# percentage of largest residues from the first regression not to be used in the second regression
		weighted = TRUE, 				# should the regression be weighted with 1/(T^2 + median(T))
		usefractions = "LandT",			# from which fractions should the decay rate be calculated: "LandT", "UandT" or "both"
		ratio = NULL, 					# coefficient to rescale the fractions
		relevant = NULL, 				# choose the arrays to be used for halflives calculation, vector due to experiments variable 
		check = TRUE, 					# if check=TRUE, control messages and plots will be generated
		labeling = TRUE, 				# should the labeling bias be plotted
		correctedlabeling = FALSE,		# should the corrected labeling bias be plotted
		regression = TRUE,				# should the regression results be plotted
		rankpairs = TRUE, 				# should the ranks of 1-L/T, U/T or (1-L/T+U/T)/2 be compared in heatpairs plot
		assessment = TRUE,				# should 1-L/T, U/T or (1-L/T+U/T)/2 be assessed due to limitations of the decay rate formula
		correlation = TRUE, 			# should the correlation be plotted
		error = FALSE,					# should standard deviation and coefficient of variation be calculated
		bicor = TRUE, 					# should the labeling bias be corrected
		condition = "", 				# to be added to the plotnames
		timepoint = "",					# to be added to the plotnames
		upper = 700, 					# upper bound for labeling bias estimation
		lower = 500, 					# lower bound for labeling bias estimation
		plots = FALSE, 					# if plots=TRUE, control plots will be saved
		notinR = FALSE,					# should plot be not plotted in R
		folder = NULL, 					# folder, where to save the plots
		addformat = NULL, 				# additional fileformat for plots to be saved
		totaloverwt = 1, 				# total mRNA over WT
		simulation = FALSE,				# should the simulation be plotted
		dynamic = FALSE,				# should be TRUE for timecourse data
		initials = NULL,				# initial values of total for timecourse data
		truemus = NULL,					# the true synthesis rates
		truemusaveraged = NULL,			# the true synthesis rates (averaged over labeling period)
		truelambdas = NULL,				# the true decay rates
		truelambdasaveraged = NULL,		# the true decay rates (averaged over labeling period)
		truehalflives = NULL,			# the true half-lives
		truehalflivesaveraged = NULL,	# the true half-lives (averaged over labeling period)
		trueplabel = NULL,				# the true labeling efficiency
		trueLasymptote = NULL,			# the true Lasymptote
		trueUasymptote = NULL,			# the true Uasymptote
		truecrbyar = NULL,				# the true quotient truecr/truear
		truecrbybr = NULL,				# the true quotient truecr/truebr
		truebrbyar = NULL				# the true quotient truebr/truear
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
			result = median( hilf[is.real(hilf)] )
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
				result = median(hilf[is.real(hilf)])
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
					result = median(hilf[is.real(hilf)])
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
				if (regression){
					if (any(windowxy(nrexperiments) > 2)){scex = 1/0.66} else if (all(windowxy(nrexperiments) == 2)){scex = 1/0.83} else {scex = 1}
					plotsfkt = function(){
						par(mfrow=windowxy(nrexperiments))
						par(mar=c(5,4,4,2) + 1)
						for (expnr in seq(experiments)){
							Tnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "T"))
							ylim = c(-max(abs(t(planes[[expnr]])[,3])),max(abs(t(planes[[expnr]])[,3])))
							heatscatter(log(t(planes[[expnr]])[,1]),t(planes[[expnr]])[,3],cor=FALSE,ylim=ylim,xlab=expression(paste("log( Orthogonal projection on ",L[g]," )")),ylab=expression("Normal of the plane"),main = expression(paste("Regression plane:  ratio of  ",L[g]/T[g],"\n")),cex.main=1.25*scex,cex.lab=1*scex,cex.axis=1*scex)
							points(log(t(planes[[expnr]])[,1])[discards[[expnr]]],t(planes[[expnr]])[,3][discards[[expnr]]],col = "red")
							title(paste(" \n \n \n (",phenomat[Tnr,"name"],"  c = ",signif(crbyarestimate[expnr],digits=2),")"),col.main="darkgrey",cex.main=0.75*scex)
							abline(h = 0,col = "green",lwd=2)
							abline(h = max(t(planes[[expnr]])[,3]),col = "red",lwd=1)
							abline(h = min(t(planes[[expnr]])[,3]),col = "red",lwd=1)
						}
					}
					plotit(filename = paste(folder,"/regression_hyperline_",condition,"_",timepoint,".jpg",sep=""),sw = 2*windowxy(nrexperiments)[2],sh = 2*windowxy(nrexperiments)[1],sres = 2,plotsfkt = plotsfkt,ww = 7*windowxy(nrexperiments)[2],wh = 7*windowxy(nrexperiments)[1],saveit = plots,addformat = addformat,notinR = notinR)
					
					plotsfkt = function(){
						par(mfrow=windowxy(nrexperiments))
						par(mar=c(5,4,4,2) + 1)
						for (expnr in seq(experiments)){
							Tnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "T"))
							zlim = c(-max(abs(t(planes[[expnr]])[,3])),max(abs(t(planes[[expnr]])[,3])))
							s3d <- scatterplot3d(t(planes[[expnr]])[,1],t(planes[[expnr]])[,2],t(planes[[expnr]])[,3],pch=20,xlab=expression(paste("Orthogonal projection on ",L[g])),ylab=expression(paste("Orthogonal projection on ",U[g])),zlab=expression("Normal of the plane"), scale.y=1, angle=40,
									highlight.3d=TRUE,main = expression(paste("Regression plane 3D view:  ratio of  ",L[g]/T[g],"\n")),grid=TRUE,zlim=zlim,cex.main=1.25*scex,cex.lab=1*scex,cex.axis=1*scex)
							title(paste(" \n \n \n (",phenomat[Tnr,"name"],"  c = ",signif(crbyarestimate[expnr],digits=2),")"),col.main="darkgrey",cex.main=0.75*scex)
							s3d$points3d(t(planes[[expnr]])[,1][discards[[expnr]]],t(planes[[expnr]])[,2][discards[[expnr]]],t(planes[[expnr]])[,3][discards[[expnr]]],pch=20,col = "red")
							s3d$plane3d(c(0,0,0), lty = "solid", lty.box = "solid",col = "green",lwd=2)
							s3d$plane3d(c(max(t(planes[[expnr]])[,3]),0,0), lty = "solid", lty.box = "solid",col = "red",lwd=1)
							s3d$plane3d(c(min(t(planes[[expnr]])[,3]),0,0), lty = "solid", lty.box = "solid",col = "red",lwd=1)
						}
					}
					plotit(filename = paste(folder,"/regression_hyperplane_",condition,"_",timepoint,".jpg",sep=""),sw = 2*windowxy(nrexperiments)[2],sh = 2*windowxy(nrexperiments)[1],sres = 2,plotsfkt = plotsfkt,ww = 7*windowxy(nrexperiments)[2],wh = 7*windowxy(nrexperiments)[1],saveit = plots,addformat = addformat,notinR = notinR)		
				}
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
		if (labeling){
			if (any(windowxy(nrexperiments) > 2)){scex = 1/0.66} else if (all(windowxy(nrexperiments) == 2)){scex = 1/0.83} else {scex = 1}
			plotsfkt = function(){
				par(mfrow=windowxy(nrexperiments))
				par(mar=c(5,4,4,2) + 1)
				for (expnr in seq(experiments)){
					Tnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "T"))
					Tname = phenomat[Tnr,"name"]
					Lnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "L"))
					Lname = phenomat[Lnr,"name"]
					thistime = as.numeric(phenomat[Tnr,"time"])
					largeT = which(tnumberreliable>upper)
					asymptote = log(median( calcdatamatreliable[largeT,Lname]/calcdatamatreliable[largeT,Tname] ))
					smallT = which(tnumberreliable<lower)
					LT.plotbias(labeled=calcdatamatreliable[,Lname],total=calcdatamatreliable[,Tname],tnumber=tnumberreliable,plabel=plabel[expnr],asymptote=asymptote,trueplabel=trueplabel[experiments[expnr]],trueasymptote=trueLasymptote[experiments[expnr]],repl = phenomat[Tnr,"name"],scex = scex)
				}
			}
			plotit(filename = paste(folder,"/estimation_bias_L_",condition,"_",timepoint,".jpg",sep=""),sw = 2*windowxy(nrexperiments)[2],sh = 2*windowxy(nrexperiments)[1],sres = 2,plotsfkt = plotsfkt,ww = 7*windowxy(nrexperiments)[2],wh = 7*windowxy(nrexperiments)[1],saveit = plots,addformat = addformat,notinR = notinR)    
			if (unlabeledfraction){
				if (is.null(ratiodummy) & ratiomethod == "bias"){
					plotsfkt = function(){
						par(mfrow=windowxy(nrexperiments))
						par(mar=c(5,4,4,2) + 1)
						for (expnr in seq(experiments)){
							Tnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "T"))
							Tname = phenomat[Tnr,"name"]
							Unr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "U"))
							Uname = phenomat[Unr,"name"]
							thistime = as.numeric(phenomat[Tnr,"time"])
							largeT = which(tnumberreliable>upper)
							asymptote = log(median( calcdatamatreliable[largeT,Uname]/calcdatamatreliable[largeT,Tname] ))
							smallT = which(tnumberreliable<lower)
							UT.plotbias(unlabeled=calcdatamatreliable[,Uname],total=calcdatamatreliable[,Tname],tnumber=tnumberreliable,plabel=plabel[expnr],asymptote=asymptote,trueplabel=trueplabel[experiments[expnr]],trueasymptote=trueUasymptote[experiments[expnr]],repl = phenomat[Tnr,"name"],ratio = brbyarestimate[expnr],scex = scex)
						}
					}
					plotit(filename = paste(folder,"/estimation_bias_U_",condition,"_",timepoint,".jpg",sep=""),sw = 2*windowxy(nrexperiments)[2],sh = 2*windowxy(nrexperiments)[1],sres = 2,plotsfkt = plotsfkt,ww = 7*windowxy(nrexperiments)[2],wh = 7*windowxy(nrexperiments)[1],saveit = plots,addformat = addformat,notinR = notinR)
				} else {
					plotsfkt = function(){
						par(mfrow=windowxy(nrexperiments))
						par(mar=c(5,4,4,2) + 1)
						for (expnr in seq(experiments)){
							Tnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "T"))
							Tname = phenomat[Tnr,"name"]
							Unr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "U"))
							Uname = phenomat[Unr,"name"]
							thistime = as.numeric(phenomat[Tnr,"time"])
							largeT = which(tnumberreliable>upper)
							asymptote = log(median( calcdatamatreliable[largeT,Uname]/calcdatamatreliable[largeT,Tname] ))
							smallT = which(tnumberreliable<lower)
							UT.plotbias(unlabeled=calcdatamatreliable[,Uname],total=calcdatamatreliable[,Tname],tnumber=tnumberreliable,plabel=plabel[expnr],asymptote=asymptote,trueplabel=trueplabel[experiments[expnr]],trueasymptote=trueUasymptote[experiments[expnr]],repl = phenomat[Tnr,"name"],scex = scex,ratio = labelratio[expnr])
						}
					}
					plotit(filename = paste(folder,"/estimation_bias_U_",condition,"_",timepoint,".jpg",sep=""),sw = 2*windowxy(nrexperiments)[2],sh = 2*windowxy(nrexperiments)[1],sres = 2,plotsfkt = plotsfkt,ww = 7*windowxy(nrexperiments)[2],wh = 7*windowxy(nrexperiments)[1],saveit = plots,addformat = addformat,notinR = notinR)
					
				}
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
		if (correctedlabeling){
			if (any(windowxy(nrexperiments) > 2)){scex = 1/0.66} else if (all(windowxy(nrexperiments) == 2)){scex = 1/0.83} else {scex = 1}
			plotsfkt = function(){
				par(mfrow=windowxy(nrexperiments))
				par(mar=c(5,4,4,2) + 1)
				for (expnr in seq(experiments)){
					Tnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "T"))
					Tname = phenomat[Tnr,"name"]
					Lnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "L"))
					Lname = phenomat[Lnr,"name"]
					thistime = as.numeric(phenomat[Tnr,"time"])
					largeT = which(tnumberreliable>upper)
					asymptote = log(median( calcdatamatreliable[largeT,Lname]/calcdatamatreliable[largeT,Tname] ))
					smallT = which(tnumberreliable<lower)
					LT.plotbias(labeled=calcdatamatreliable[,Lname],total=calcdatamatreliable[,Tname],tnumber=tnumberreliable,plabel=1,asymptote=asymptote,trueplabel=1,trueasymptote=trueLasymptote[experiments[expnr]],repl = phenomat[Tnr,"name"],scex = scex,correctedlabeling = correctedlabeling)
				}
			}
			plotit(filename = paste(folder,"/estimation_bias_L_",condition,"_",timepoint,"_corrected.jpg",sep=""),sw = 2*windowxy(nrexperiments)[2],sh = 2*windowxy(nrexperiments)[1],sres = 2,plotsfkt = plotsfkt,ww = 7*windowxy(nrexperiments)[2],wh = 7*windowxy(nrexperiments)[1],saveit = plots,addformat = addformat,notinR = notinR)    
			if (unlabeledfraction){
				if (is.null(ratiodummy) & ratiomethod == "bias"){
					plotsfkt = function(){
						par(mfrow=windowxy(nrexperiments))
						par(mar=c(5,4,4,2) + 1)
						for (expnr in seq(experiments)){
							Tnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "T"))
							Tname = phenomat[Tnr,"name"]
							Unr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "U"))
							Uname = phenomat[Unr,"name"]
							thistime = as.numeric(phenomat[Tnr,"time"])
							largeT = which(tnumberreliable>upper)
							asymptote = log(median( calcdatamatreliable[largeT,Uname]/calcdatamatreliable[largeT,Tname] ))
							smallT = which(tnumberreliable<lower)
							UT.plotbias(unlabeled=calcdatamatreliable[,Uname],total=calcdatamatreliable[,Tname],tnumber=tnumberreliable,plabel=1,asymptote=asymptote,trueplabel=1,trueasymptote=trueUasymptote[experiments[expnr]],repl = phenomat[Tnr,"name"],ratio = brbyarestimate[expnr],scex = scex,correctedlabeling = correctedlabeling)
						}
					}
					plotit(filename = paste(folder,"/estimation_bias_U_",condition,"_",timepoint,"_corrected.jpg",sep=""),sw = 2*windowxy(nrexperiments)[2],sh = 2*windowxy(nrexperiments)[1],sres = 2,plotsfkt = plotsfkt,ww = 7*windowxy(nrexperiments)[2],wh = 7*windowxy(nrexperiments)[1],saveit = plots,addformat = addformat,notinR = notinR)
				} else {
					plotsfkt = function(){
						par(mfrow=windowxy(nrexperiments))
						par(mar=c(5,4,4,2) + 1)
						for (expnr in seq(experiments)){
							Tnr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "T"))
							Tname = phenomat[Tnr,"name"]
							Unr = which((phenomat[,"nr"] == experiments[expnr]) & (phenomat[,"fraction"] == "U"))
							Uname = phenomat[Unr,"name"]
							thistime = as.numeric(phenomat[Tnr,"time"])
							largeT = which(tnumberreliable>upper)
							asymptote = log(median( calcdatamatreliable[largeT,Uname]/calcdatamatreliable[largeT,Tname] ))
							smallT = which(tnumberreliable<lower)
							UT.plotbias(unlabeled=calcdatamatreliable[,Uname],total=calcdatamatreliable[,Tname],tnumber=tnumberreliable,plabel=1,asymptote=asymptote,trueplabel=1,trueasymptote=trueUasymptote[experiments[expnr]],repl = phenomat[Tnr,"name"],scex = scex,correctedlabeling = correctedlabeling)
						}
					}
					plotit(filename = paste(folder,"/estimation_bias_U_",condition,"_",timepoint,"_corrected.jpg",sep=""),sw = 2*windowxy(nrexperiments)[2],sh = 2*windowxy(nrexperiments)[1],sres = 2,plotsfkt = plotsfkt,ww = 7*windowxy(nrexperiments)[2],wh = 7*windowxy(nrexperiments)[1],saveit = plots,addformat = addformat,notinR = notinR)
					
				}
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
	
	LTmat = as.matrix(calcdatamat[,phenomat[which(phenomat[,"fraction"] == "T"),"name"]]-calcdatamat[,phenomat[which(phenomat[,"fraction"] == "L"),"name"]])
	if (ncol(as.matrix(LTmat)) == ncol(as.matrix(initials))) {LTmat = as.matrix(LTmat/initials)} else {LTmat = as.matrix(LTmat/apply(initials,1,median))}
	colnames(LTmat) = phenomat[which(phenomat[,"fraction"] == "T"),"name"]
	LTmat = cbind(LTmat,median = apply(as.matrix(LTmat),1,median))
	if (mediancenter){LTmat = medctr(LTmat,reliable,logscale = FALSE,protocol = FALSE)}
	
	LtoTmat = as.matrix(calcdatamat[,phenomat[which(phenomat[,"fraction"] == "L"),"name"]])
	if (ncol(as.matrix(LtoTmat)) == ncol(as.matrix(initials))) {LtoTmat = as.matrix(LtoTmat/initials)} else {LtoTmat = as.matrix(LtoTmat/apply(initials,1,median))}
	colnames(LtoTmat) = phenomat[which(phenomat[,"fraction"] == "T"),"name"]
	LtoTmat = cbind(LtoTmat,median = apply(as.matrix(LtoTmat),1,median))
	if (mediancenter){LtoTmat = medctr(LtoTmat,reliable,logscale = FALSE,protocol = FALSE)}
	
	### U/T ###
	
	if (unlabeledfraction){
		UTmat = as.matrix(calcdatamat[,phenomat[which(phenomat[,"fraction"] == "U"),"name"]])
		if (ncol(as.matrix(UTmat)) == ncol(as.matrix(initials))) {UTmat = as.matrix(UTmat/initials)} else {UTmat = as.matrix(UTmat/apply(initials,1,median))}
		colnames(UTmat) = phenomat[which(phenomat[,"fraction"] == "T"),"name"]
		UTmat = cbind(UTmat,median = apply(as.matrix(UTmat),1,median))
		if (mediancenter){UTmat = medctr(UTmat,reliable,logscale = FALSE,protocol = FALSE)}
	}
	
	if (usefractions == "LandT"){
		if (check){
			if (rankpairs){
				if (!is.null(initialdummy)){
					main = expression(paste("Rank heatpairs of  ",(T['g,i']-c['l,i']*L['g,i'])/T['g,i'-1],"\n"))
				} else {
					main = expression(paste("Rank heatpairs of  ",1-c[l]*L[g]/T[g],"\n"))
				}
				rcex = ncol(LTmat)
				plotsfkt = function(){
					par(mar=c(5,4,4,2) + 1)
					heatmat = apply(LTmat,2,rank)
					colnames(heatmat) = colnames(LTmat)
					rownames(heatmat) = rownames(LTmat)
					heatpairs(heatmat,main = main,cex.labels = min(1/(max(nchar(colnames(LTmat)))/10)*2.5*scex,2),cex.main=1.25*scex,cex.lab=1*scex,cex.axis=1*scex)
					title(paste("( max median fold = ",round(max(abs(diff(apply(LTmat,2,median)))),2),")"),col.main="darkgrey",cex.main=0.75*scex)
				}
				plotit(filename = paste(folder,"/rank_heatpairs_",condition,"_",timepoint,".jpg",sep=""),sw = rcex,sh = rcex,sres = 2,plotsfkt = plotsfkt,ww = rcex*3.5,wh = rcex*3.5,saveit = plots,addformat = addformat,notinR = notinR)
			}
			
			if (assessment){
				if (any(windowxy(ncol(LTmat)) > 2)){scex = 1/0.66} else if (all(windowxy(ncol(LTmat)) == 2)){scex = 1/0.83} else {scex = 1}
				plotsfkt = function(){
					par(mfrow=windowxy(ncol(LTmat)))
					par(mar=c(5,4,4,2) + 1)
					for (i in 1:ncol(LTmat)){
						assess = function(x){- alpha - ((1/labelingtime)*log(x))}
						if (!is.null(initialdummy)){
							xlab = expression(paste((T['g,i']-c['l,i']*L['g,i'])/T['g,i'-1]))
							ylab=expression(paste("Decay rate  ",lambda['g,i']))
							main = expression(paste("Limit assessment of  ",(T['g,i']-c['l,i']*L['g,i'])/T['g,i'-1],"\n"))
						} else {
							xlab = expression(paste(1-c[l]*L[g]/T[g]))
							ylab=expression(paste("Decay rate  ",lambda[g]))
							main = expression(paste("Limit assessment of  ",1-c[l]*L[g]/T[g],"\n"))
						}
						plot(0,type="n",xlim=c(-2*exp(-alpha*labelingtime),3*exp(-alpha*labelingtime)),ylim=c(-0.5,max(hist(LTmat[,ncol(LTmat)],breaks=seq(floor(min(LTmat)),ceiling(max(LTmat)),0.25),plot=FALSE)$density)/2*3),xlab=xlab,ylab=ylab,main=main,cex.main=1.25*scex,cex.lab=1*scex,cex.axis=1*scex)
						title(paste(" \n \n \n (",colnames(LTmat)[i],")"),col.main="darkgrey",cex.main=0.75*scex)
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
				plotit(filename = paste(folder,"/range_assessment_",condition,"_",timepoint,".jpg",sep=""),sw = 2*windowxy(ncol(LTmat))[2],sh = 2*windowxy(ncol(LTmat))[1],sres = 2,plotsfkt = plotsfkt,ww = 7*windowxy(ncol(LTmat))[2],wh = 7*windowxy(ncol(LTmat))[1],saveit = plots,addformat = addformat,notinR = notinR)
			}
		}
		
		nrexp = ncol(LTmat)
		
		results[["LTmat"]] = LtoTmat
		LT = LtoTmat[,nrexp]
		results[["LT"]] = LT
		
		LTdrmat = - alpha - ((1/labelingtime)*log(LTmat))
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
			if (rankpairs){
				if (!is.null(initialdummy)){
					main = expression(paste("Rank heatpairs of  ",c['u,i']*U['g,i']/T['g,i'-1],"\n"))
				} else {
					main = expression(paste("Rank heatpairs of  ",c[u]*U[g]/T[g],"\n"))
				}
				rcex = ncol(UTmat)
				plotsfkt = function(){
					par(mar=c(5,4,4,2) + 1)
					heatmat = apply(UTmat,2,rank)
					colnames(heatmat) = colnames(UTmat)
					rownames(heatmat) = rownames(UTmat)
					heatpairs(heatmat,main = main,cex.labels = min(1/(max(nchar(colnames(UTmat)))/10)*2.5*scex,2),cex.main=1.25*scex,cex.lab=1*scex,cex.axis=1*scex)
					title(paste("( max median fold = ",round(max(abs(diff(apply(UTmat,2,median)))),2),")"),col.main="darkgrey",cex.main=0.75*scex)
				}
				plotit(filename = paste(folder,"/rank_heatpairs_",condition,"_",timepoint,".jpg",sep=""),sw = rcex,sh = rcex,sres = 2,plotsfkt = plotsfkt,ww = rcex*3.5,wh = rcex*3.5,saveit = plots,addformat = addformat,notinR = notinR)
			}
			
			if (assessment){
				if (any(windowxy(ncol(UTmat)) > 2)){scex = 1/0.66} else if (all(windowxy(ncol(UTmat)) == 2)){scex = 1/0.83} else {scex = 1}
				plotsfkt = function(){
					par(mfrow=windowxy(ncol(UTmat)))
					par(mar=c(5,4,4,2) + 1)
					for (i in 1:ncol(UTmat)){
						assess = function(x){- alpha - ((1/labelingtime)*log(x))}
						if (!is.null(initialdummy)){
							xlab = expression(paste(c['u,i']*U['g,i']/T['g,i'-1]))
							ylab=expression(paste("Decay rate  ",lambda['g,i']))
							main = expression(paste("Limit assessment of  ",c['u,i']*U['g,i']/T['g,i'-1],"\n"))
						} else {
							xlab = expression(paste(c[u]*U[g]/T[g]))
							ylab=expression(paste("Decay rate  ",lambda[g]))
							main = expression(paste("Limit assessment of  ",c[u]*U[g]/T[g],"\n"))
						}
						plot(0,type="n",xlim=c(-2*exp(-alpha*labelingtime),3*exp(-alpha*labelingtime)),ylim=c(-0.5,max(hist(UTmat[,ncol(UTmat)],breaks=seq(floor(min(UTmat)),ceiling(max(UTmat)),0.25),plot=FALSE)$density)/2*3),xlab=xlab,ylab=ylab,main=main,cex.main=1.25*scex,cex.lab=1*scex,cex.axis=1*scex)
						title(paste(" \n \n \n (",colnames(UTmat)[i],")"),col.main="darkgrey",cex.main=0.75*scex)
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
				plotit(filename = paste(folder,"/range_assessment_",condition,"_",timepoint,".jpg",sep=""),sw = 2*windowxy(ncol(UTmat))[2],sh = 2*windowxy(ncol(UTmat))[1],sres = 2,plotsfkt = plotsfkt,ww = 7*windowxy(ncol(UTmat))[2],wh = 7*windowxy(ncol(UTmat))[1],saveit = plots,addformat = addformat,notinR = notinR)
			}
		}
		
		nrexp = ncol(UTmat)
		
		results[["UTmat"]] = UTmat
		UT = UTmat[,nrexp]
		results[["UT"]] = UT
		
		UTdrmat = - alpha - ((1/labelingtime)*log(UTmat))
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
			if (rankpairs){
				if (!is.null(initialdummy)){
					main = expression(paste("Rank heatpairs of  ",((T['g,i']-c['l,i']*L['g,i'])/T['g,i'-1]+c['u,i']*U['g,i']/T['g,i'-1])/2,"\n"))
				} else {
					main = expression(paste("Rank heatpairs of  ",(1-c[l]*L[g]/T[g]+c[u]*U[g]/T[g])/2,"\n"))
				}
				rcex = ncol(Bmat)
				plotsfkt = function(){
					par(mar=c(5,4,4,2) + 1)
					heatmat = apply(Bmat,2,rank)
					colnames(heatmat) = colnames(Bmat)
					rownames(heatmat) = rownames(Bmat)
					heatpairs(heatmat,main = main,cex.labels = min(1/(max(nchar(colnames(Bmat)))/10)*2.5*scex,2),cex.main=1.25*scex,cex.lab=1*scex,cex.axis=1*scex)
					title(paste("( max median fold = ",round(max(abs(diff(apply(Bmat,2,median)))),2),")"),col.main="darkgrey",cex.main=0.75*scex)
				}
				plotit(filename = paste(folder,"/rank_heatpairs_",condition,"_",timepoint,".jpg",sep=""),sw = rcex,sh = rcex,sres = 2,plotsfkt = plotsfkt,ww = rcex*3.5,wh = rcex*3.5,saveit = plots,addformat = addformat,notinR = notinR)
			}
			
			if (assessment){
				if (any(windowxy(ncol(Bmat)) > 2)){scex = 1/0.66} else if (all(windowxy(ncol(Bmat)) == 2)){scex = 1/0.83} else {scex = 1}
				plotsfkt = function(){
					par(mfrow=windowxy(ncol(Bmat)))
					par(mar=c(5,4,4,2) + 1)
					for (i in 1:ncol(Bmat)){
						assess = function(x){- alpha - ((1/labelingtime)*log(x))}
						if (!is.null(initialdummy)){
							xlab = expression(paste(((T['g,i']-c['l,i']*L['g,i'])/T['g,i'-1]+c['u,i']*U['g,i']/T['g,i'-1])/2))
							ylab=expression(paste("Decay rate  ",lambda['g,i']))
							main = expression(paste("Limit assessment of  ",((T['g,i']-c['l,i']*L['g,i'])/T['g,i'-1]+c['u,i']*U['g,i']/T['g,i'-1])/2,"\n"))
						} else {
							xlab = expression(paste((1-c[l]*L[g]/T[g]+c[u]*U[g]/T[g])/2))
							ylab=expression(paste("Decay rate  ",lambda[g]))
							main = expression(paste("Limit assessment of  ",(1-c[l]*L[g]/T[g]+c[u]*U[g]/T[g])/2,"\n"))
						}
						plot(0,type="n",xlim=c(-2*exp(-alpha*labelingtime),3*exp(-alpha*labelingtime)),ylim=c(-0.5,max(hist(Bmat[,ncol(Bmat)],breaks=seq(floor(min(Bmat)),ceiling(max(Bmat)),0.25),plot=FALSE)$density)/2*3),xlab=xlab,ylab=ylab,main=main,cex.main=1.25*scex,cex.lab=1*scex,cex.axis=1*scex)
						title(paste(" \n \n \n (",colnames(Bmat)[i],")"),col.main="darkgrey",cex.main=0.75*scex)
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
				plotit(filename = paste(folder,"/range_assessment_",condition,"_",timepoint,".jpg",sep=""),sw = 2*windowxy(ncol(Bmat))[2],sh = 2*windowxy(ncol(Bmat))[1],sres = 2,plotsfkt = plotsfkt,ww = 7*windowxy(ncol(Bmat))[2],wh = 7*windowxy(ncol(Bmat))[1],saveit = plots,addformat = addformat,notinR = notinR)
			}
		}
		
		nrexp = ncol(Bmat)
		
		results[["Bmat"]] = Bmat
		B = Bmat[,nrexp]
		results[["B"]] = B
		
		Bdrmat = - alpha - ((1/labelingtime)*log(Bmat))
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
	
	if (check){
		if (correlation) {
			vecmat = cbind(TE,LE,tnumber,sr,dr,(log(2)/dr))
			vecmat = vecmat[reliable,]
			compsmat = matrix(0,ncol=3,nrow=3)
			for (i in 1:3){for (j in 1:3){compsmat[i,j] = cor(vecmat[,i],vecmat[,j+3],method = "pearson",use = "na.or.complete")}}
			rownames(compsmat) = c("TE","LE","#U")
			colnames(compsmat) = c("SR","DR","HL")
			print(t(compsmat))
			
			plotsfkt = function(){
				scex = 1/0.83
				par(mfrow=c(1,2))
				layout(matrix(data=c(1,2),nrow=1,ncol=2),widths=c(6,1)) 
				par(mar=c(5,4,4,0) + 1)
				image(t(t(compsmat)[3:1,]),axes=FALSE,col = colorRampPalette(c("darkred","lightgrey","darkgreen"))(500),breaks = seq(-1,1,length.out=501),main=expression(paste("Correlation analysis \n")),xlab=expression("Data inputs"),ylab=expression("Extracted rates"),cex.main=1.25*scex,cex.lab=1*scex,cex.axis=1*scex)
				title(paste("( pearson correlation for gene-wise medians )"),col.main="darkgrey",cex.main=0.75*scex)
				axis(1,at=c(0,0.5,1),labels=c(expression(T[g]),expression(L[g]),expression('#u')))
				axis(2,at=c(0,0.5,1),labels=c(expression(t[1/2][g]),expression(lambda[g]),expression(mu[g])),las=2)
				abline(v=c(0.75,0.25),h=c(0.75,0.25),cex=0.5)
				box()
				text(c(0,0.5,1,0,0.5,1,0,0.5,1),c(0,0,0,0.5,0.5,0.5,1,1,1),round(as.vector(t(t(compsmat)[3:1,])),2))
				par(mar = c(5,0,4,2) + 1)
				image(1,seq(-1,1,length.out=501),matrix(data=seq(-1,1,length.out=501),ncol=length(seq(-1,1,length.out=501)),nrow=1),col=colorRampPalette(c("darkred","lightgrey","darkgreen"))(500),xlab="",ylab="",axes=FALSE)
				axis(4,at=c(-1,-0.5,0,0.5,1),labels=c(-1,-0.5,0,0.5,1),las=2)
				box()
			}
			plotit(filename = paste(folder,"/correlation_analysis_",condition,"_",timepoint,".jpg",sep=""),sw = 2,sh = 2,sres = 2,plotsfkt = plotsfkt,ww = 7,wh = 7,saveit = plots,addformat = addformat,notinR = notinR)
		}
		
		### GENERATE CV PROGRESSION ###

		TEmat = cbind(calcdatamat[,which(phenomat[,"fraction"] == "T")])
		LEmat = cbind(calcdatamat[,which(phenomat[,"fraction"] == "L")])
		results[["TE.sd"]] = apply(TEmat,1,sd)
		results[["LE.sd"]] = apply(LEmat,1,sd)
		results[["TE.cv"]] = apply(TEmat,1,sd)/apply(TEmat,1,mean)
		results[["LE.cv"]] = apply(LEmat,1,sd)/apply(LEmat,1,mean)
		TEmat = log(TEmat)
		LEmat = log(LEmat)
		
		if (unlabeledfraction){
			UEmat = cbind(calcdatamat[,which(phenomat[,"fraction"] == "U")])
			results[["UE.sd"]] = apply(UEmat,1,sd)
			results[["UE.cv"]] = apply(UEmat,1,sd)/apply(UEmat,1,mean)
			UEmat = log(UEmat)
		}
		
		SDmat = matrix(NA,ncol=3,nrow=nrow(TEmat))
		rownames(SDmat) = rownames(TEmat)
		colnames(SDmat) = c("DR","HL","SR")
		CVmat = matrix(NA,ncol=3,nrow=nrow(TEmat))
		rownames(CVmat) = rownames(TEmat)
		colnames(CVmat) = c("DR","HL","SR")
		if (error){
			if (dynamic){errorpossible = (!is.null(dim(initials)) & !is.null(dim(TEmat)))} else {errorpossible = !is.null(dim(TEmat))}
			if (errorpossible){
				if (dynamic){
					if (!is.null(initialdummy)){
						TImat = initials
						results[["TI.sd"]] = apply(TImat,1,sd)
						results[["TI.cv"]] = apply(TImat,1,sd)/apply(TImat,1,mean)
						TImat = log(TImat)
						Imean = mean(TImat[i,])
						Isd = sd(TImat[i,])
					}  else {
						Imean = NULL
						Isd = NULL
					}
				} else {
					Imean = NULL
					Isd = NULL
				}
				if (usefractions == "LandT"){
					for (i in names(dr[!is.na(dr)])){
						err = LT.error.progression(Tmean = mean(TEmat[i,]),Tsd = sd(TEmat[i,]),Lmean = mean(LEmat[i,]),Lsd = sd(LEmat[i,]),timepoint = labelingtime,alpha = alpha,Imean = Imean,Isd = Isd)
						SDmat[i,] = err[["sds"]]
						CVmat[i,] = err[["cvs"]]
					}
				}
				if (usefractions == "UandT"){
					for (i in names(dr[!is.na(dr)])){
						err = UT.error.progression(Tmean = mean(TEmat[i,]),Tsd = sd(TEmat[i,]),Lmean = mean(LEmat[i,]),Lsd = sd(LEmat[i,]),Umean = mean(UEmat[i,]),Usd = sd(UEmat[i,]),timepoint = labelingtime,alpha = alpha,Imean = Imean,Isd = Isd)
						SDmat[i,] = err[["sds"]]
						CVmat[i,] = err[["cvs"]]
					}					
				}
				if (usefractions == "both"){
					for (i in names(dr[!is.na(dr)])){
						err = BOTH.error.progression(Tmean = mean(TEmat[i,]),Tsd = sd(TEmat[i,]),Lmean = mean(LEmat[i,]),Lsd = sd(LEmat[i,]),Umean = mean(UEmat[i,]),Usd = sd(UEmat[i,]),timepoint = labelingtime,alpha = alpha,Imean = Imean,Isd = Isd)
						SDmat[i,] = err[["sds"]]
						CVmat[i,] = err[["cvs"]]
					}					
				}
			} else {print("Your data is not suitable for error assessment due to lack of replicates")}	
		}
		results[["dr.sd"]] = SDmat[,"DR"]
		results[["hl.sd"]] = SDmat[,"HL"]
		results[["sr.sd"]] = SDmat[,"SR"]
		results[["dr.cv"]] = CVmat[,"DR"]
		results[["hl.cv"]] = CVmat[,"HL"]
		results[["sr.cv"]] = CVmat[,"SR"]
		
		### SIMULATION ###
		
		if (simulation){	
			plotsfkt = function(){
				scex = 1/0.66
				par(mfrow=c(3,3))
				par(mar=c(5,4,4,2) + 1.5)
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
				
				heatscatter(truehl,results[["hl"]],main=expression(paste("Half-life  ",t[1/2][g],"\n")),xlab=expression(paste("True half-life  ",t[1/2][g])),ylab=expression(paste("Estimated half-life  ",t[1/2][g]^'*')),cor=FALSE,xlim=hllims,ylim=hllims,cex.main=1.4*scex,cex.lab=1.15*scex,cex.axis=1*scex)
				title(paste(" \n \n \n ( heatscatter p-cor = ",pcor(truehl,results[["hl"]],use="na.or.complete"),")"),col.main="darkgrey",cex.main=0.75*scex)
				linfactor = coefficients(tls((results[["hl"]][hlindi]) ~ (truehl[hlindi]) + 0))[1]
				abline(a=0,b=linfactor,col="black",lwd=3)
				abline(a=0,b=1,col="grey",lwd=3,lty=2)
				cat("Factor estimated/true:",linfactor,"\n")
				heatscatter(rank(truehl),rank(results[["hl"]]),main=expression(paste("Half-life  ",t[1/2][g],"  rank")),xlab=expression(paste("True half-life  ",t[1/2][g],"  rank")),ylab=expression(paste("Estimated half-life  ",t[1/2][g]^'*',"  rank")),cor=FALSE,cex.main=1.4*scex,cex.lab=1.15*scex,cex.axis=1*scex)
				title(paste(" \n \n \n ( heatscatter s-cor = ",scor(rank(truehl),rank(results[["hl"]]),use="na.or.complete"),")"),col.main="darkgrey",cex.main=0.75*scex)
				abline(a=0,b=1,col="black",lwd=3)
				abline(a=0,b=1,col="grey",lwd=3,lty=2)
				meanreldev = mean(abs(results[["hl"]][hlindi]-truehl[hlindi])/truehl[hlindi])
				cat("Mean relative deviation: ",signif(meanreldev,2),"\n")
				den = density(log2(results[["hl"]][hlindi]/truehl[hlindi]))
				den = cbind(den$x,den$y)
				hist(log2(results[["hl"]][hlindi]/truehl[hlindi]),xlim=c(-5,5),breaks=40,main=expression(paste("Log-ratio  log2( ",t[1/2][g]^'*'/t[1/2][g]," )")),xlab=expression(paste("log2( ",t[1/2][g]^'*'/t[1/2][g]," )")),cex.main=1.4*scex,cex.lab=1.15*scex,cex.axis=1*scex)
				title(paste(" \n \n \n ( MRD: ",signif(meanreldev,2), ",mode: ",signif(den[which(den[,2] == max(den[,2]))],2)," )"),col.main="darkgrey",cex.main=0.75*scex)
				abline(v=0,col="grey",lwd=3,lty=2)
				abline(v=den[which(den[,2] == max(den[,2]))],col="black",lwd=3)
				
				heatscatter(truedr,results[["dr"]],log="xy",main=expression(paste("Deacy rate  ",lambda[g],"\n")),xlab=expression(paste("True decay rate  ",lambda[g])),ylab=expression(paste("Estimated decay rate  ",lambda[g]^'*')),cor=FALSE,xlim=drlims,ylim=drlims,cex.main=1.4*scex,cex.lab=1.15*scex,cex.axis=1*scex)
				title(paste(" \n \n \n ( heatscatter p-cor = ",pcor(truedr,results[["dr"]],use="na.or.complete"),")"),col.main="darkgrey",cex.main=0.75*scex)
				lindev = median(log10(results[["dr"]][drindi])-log10(truedr[drindi]))
				abline(lindev,b=1,col="black",lwd=3)
				abline(a=0,b=1,col="grey",lwd=3,lty=2)
				heatscatter(rank(truedr),rank(results[["dr"]]),main=expression(paste("Decay rate  ",lambda[g],"  rank")),xlab=expression(paste("True decay rate  ",lambda[g],"  rank")),ylab=expression(paste("Estimated decay  ",lambda[g]^'*',"  rank")),cor=FALSE,cex.main=1.4*scex,cex.lab=1.15*scex,cex.axis=1*scex)
				title(paste(" \n \n \n ( heatscatter s-cor = ",scor(rank(truedr),rank(results[["dr"]]),use="na.or.complete"),")"),col.main="darkgrey",cex.main=0.75*scex)
				abline(a=0,b=1,col="black",lwd=3)
				abline(a=0,b=1,col="grey",lwd=3,lty=2)
				meanreldev = mean(abs(results[["dr"]][drindi]-truedr[drindi])/truedr[drindi])
				cat("Mean relative deviation: ",signif(meanreldev,2),"\n")
				den = density(log2(results[["dr"]][drindi]/truedr[drindi]))
				den = cbind(den$x,den$y)
				hist(log2(results[["dr"]][drindi]/truedr[drindi]),xlim=c(-5,5),breaks=40,main=expression(paste("Log-ratio  log2( ",lambda[g]^'*'/lambda[g]," )")),xlab=expression(paste("log2( ",lambda[g]^'*'/lambda[g]," )")),cex.main=1.4*scex,cex.lab=1.15*scex,cex.axis=1*scex)
				title(paste(" \n \n \n ( MRD: ",signif(meanreldev,2), ",mode: ",signif(den[which(den[,2] == max(den[,2]))],2)," )"),col.main="darkgrey",cex.main=0.75*scex)
				abline(v=0,col="grey",lwd=3,lty=2)
				abline(v=den[which(den[,2] == max(den[,2]))],col="black",lwd=3)
				
				heatscatter(truesr,results[["sr"]],main=expression(paste("Synthesis rate  ",mu[g],"\n")),xlab=expression(paste("True synthesis rate  ",mu[g])),ylab=expression(paste("Estimated synthesis rate  ",mu[g]^'*')),cor=FALSE,xlim=srlims,ylim=srlims,cex.main=1.4*scex,cex.lab=1.15*scex,cex.axis=1*scex)
				title(paste(" \n \n \n ( heatscatter p-cor = ",pcor(truesr,results[["sr"]],use="na.or.complete"),")"),col.main="darkgrey",cex.main=0.75*scex)
				linfactor = coefficients(tls((results[["sr"]][srindi]) ~ (truesr[srindi]) + 0))[1]
				abline(a=0,b=linfactor,col="black",lwd=3)
				abline(a=0,b=1,col="grey",lwd=3,lty=2)
				heatscatter(rank(truesr),rank(results[["sr"]]),main=expression(paste("Synthesis rate  ",mu[g],"  rank")),xlab=expression(paste("True synthesis rate  ",mu[g],"  rank")),ylab=expression(paste("Estimated synthesis  ",mu[g]^'*',"  rank")),cor=FALSE,cex.main=1.4*scex,cex.lab=1.15*scex,cex.axis=1*scex)
				title(paste(" \n \n \n ( heatscatter s-cor = ",scor(rank(truesr),rank(results[["sr"]]),use="na.or.complete"),")"),col.main="darkgrey",cex.main=0.75*scex)
				abline(a=0,b=1,col="black",lwd=3)
				abline(a=0,b=1,col="grey",lwd=3,lty=2)
				meanreldev = mean(abs(results[["sr"]][srindi]-truesr[srindi])/truesr[srindi])
				cat("Mean relative deviation: ",signif(meanreldev,2),"\n")
				den = density(log2(results[["sr"]][srindi]/truesr[srindi]))
				den = cbind(den$x,den$y)
				hist(log2(results[["sr"]][srindi]/truesr[srindi]),xlim=c(-5,5),breaks=40,main=expression(paste("Log-ratio  log2( ",mu[g]^'*'/mu[g]," )")),xlab=expression(paste("log2( ",mu[g]^'*'/mu[g]," )")),cex.main=1.4*scex,cex.lab=1.15*scex,cex.axis=1*scex)
				title(paste(" \n \n \n ( MRD: ",signif(meanreldev,2), ",mode: ",signif(den[which(den[,2] == max(den[,2]))],2)," )"),col.main="darkgrey",cex.main=0.75*scex)
				abline(v=0,col="grey",lwd=3,lty=2)
				abline(v=den[which(den[,2] == max(den[,2]))],col="black",lwd=3)
			}
			plotit(filename = paste(folder,"/simulation_",condition,"_",timepoint,".jpg",sep=""),sw = 4,sh = 4,sres = 2,plotsfkt = plotsfkt,ww = 14,wh = 14,saveit = plots,addformat = addformat,notinR = notinR)
		}
	}
	
	### RETURN RESULTS IN A LIST ###
	
	return(results)
}


### estimate uses an experiment, given by a phenotype matrix, data matrix and the #uridines for each gene to estimate synthesis and decay rate of the genes ###


DTA.estimate = function(phenomat = NULL,# phenotype matrix, "nr" should be numbered by experiments not by replicates
		datamat = NULL, 				# data matrix, should only contain the rows of phenomat as columns
		tnumber = NULL, 				# #uridines, should have the rownames of datamat
		ccl = NULL, 					# the cell cycle length of the cells
		mRNAs = NULL,					# estimated number of mRNAs in a cell
		reliable = NULL, 				# vector of reliable genes, which are used for regression
		mediancenter = TRUE, 			# should the L/T resp. U/T ratio of replicates be rescaled to a common median before rate extraction
		usefractions = "LandT",			# from which fractions should the decay rate be calculated: "LandT", "UandT" or "both"
		LtoTratio = NULL, 				# coefficient to rescale L/T
		ratiomethod = "tls", 			# choose the regression method to be used, possible methods are: tls, bias, lm
		largest = 5, 					# percentage of largest residues from the first regression not to be used in the second regression
		weighted = TRUE, 				# should the regression be weighted with 1/(T^2 + median(T))
		relevant = NULL, 				# choose the arrays to be used for halflives calculation, vector due to experiments variable 
		check = TRUE, 					# if check=TRUE, control messages and plots will be generated
		regression = TRUE,				# should the regression results be plotted
		labeling = TRUE, 				# should the labeling bias be plotted
		correctedlabeling = FALSE,		# should the corrected labeling bias be plotted
		rankpairs = TRUE, 				# should the ranks of 1-L/T, U/T or (1-L/T+U/T)/2 be compared in heatpairs plot
		assessment = TRUE,				# should 1-L/T, U/T or (1-L/T+U/T)/2 be assessed due to limitations of the decay rate formula
		correlation = TRUE, 			# should the correlation be plotted
		error = FALSE,					# should standard deviation and coefficient of variation be calculated
		bicor = TRUE, 					# should the labeling bias be corrected
		condition = "", 				# to be added to the plotnames
		upper = 700, 					# upper bound for labeling bias estimation
		lower = 500, 					# lower bound for labeling bias estimation
		plots = FALSE, 					# if plots=TRUE, control plots will be saved
		notinR = FALSE,					# should plot be not plotted in R
		folder = NULL, 					# folder, where to save the plots
		addformat = NULL, 				# additional fileformat for plots to be saved
		totaloverwt = 1, 				# total mRNA over WT
		simulation = FALSE,				# should the simulation be plotted
		sim.object = NULL				# simulation object created by DTA.generate
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
	if (!setequal(rownames(phenomat),colnames(datamat))){stop("The rownames of the phenomat should be equal to the colnames of the datamat !")}
	if (!setequal(rownames(datamat),names(tnumber))){stop("The colnames of the datamat should be equal to the names of tnumber !")}
	if (is.null(ccl)){print("If you do not specify the Cell Cycle Length (ccl), growth is set to zero or as specified in your sim.object !")}
	if (is.null(mRNAs)){print("If you do not specify the number of mRNAs in the cells (mRNAs), the synthesis rate will be in arbitrary scale !")}
	if (is.null(reliable)){print("If you do not specify a vector of reliable identifiers (reliable), the parameter estimation is done on all identifiers !")
		reliable = rownames(datamat)}
	if (!any(usefractions %in% c("LandT","UandT","both"))){stop("usefractions need to be 'LandT', 'UandT' or 'both' !")}
	if (!any(ratiomethod %in% c("tls","lm","bias"))){stop("ratiomethod need to be 'bias', 'tls' or 'lm' !")}
	if (is.null(LtoTratio)){print("If you do not specify ratio of L to T (LtoTratio), it is estimated from the data !")}
	if (!is.null(LtoTratio)){
		ratiomethod = "tls"
		if (length(unique(phenomat[,"time"])) != length(LtoTratio)){stop(paste("The number of specified ratios of L to T (LtoTratio) should correspond to the number of labeling durations:",length(unique(phenomat[,"time"])),"!"))}
	}
	if (plots){if (is.null(folder)){stop("You need to specify the folder, where the plots should be saved !")}}
	
	### PRELIMINARIES ###
	
	datamat = datamat[,rownames(phenomat)]
	tnumber = tnumber[rownames(datamat)]
	
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
				reliable = reliable,mediancenter = mediancenter,usefractions = usefractions,ratiomethod = ratiomethod,assessment = assessment,
				largest = largest,weighted = weighted,ratio = ratios[[labtime]],relevant = relevants[[labtime]],check = check,labeling = labeling,
				totaloverwt = totaloverwt,correctedlabeling = correctedlabeling,rankpairs = rankpairs,correlation = correlation,error = error,bicor = bicor,
				condition = condition,timepoint = labtime,regression = regression,upper = upper,lower = lower,plots = plots,notinR = notinR,folder = folder,
				addformat = addformat,simulation = simulation,truemus = truemus,truelambdas = truelambdas,truehalflives = truehalflives,trueplabel=trueplabels[[labtime]],
				trueLasymptote=trueLasymptotes[[labtime]],trueUasymptote=trueUasymptotes[[labtime]],truecrbyar=truecrbyars[[labtime]],
				truecrbybr=truecrbybrs[[labtime]],truebrbyar=truebrbyars[[labtime]])
	}
	
	### RETURN RESULTS IN A LIST ###
	
	return(res)
}



