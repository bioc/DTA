

### DTA.singlegenerate generates one experiment consisting of total, labeled and unlabeled measurements of a simulated sample. DATA IS RETURNED ON THE ABSOLUTE SCALE !!!! ###


DTA.singlegenerate = function(timepoint,		# the timepoint at which the samples are taken
		tnumber = NULL,   						# the number of uridine residues for each gene
		plabel = 0.005,							# the efficiency with which a Uridine position is finally Biotynilated
		nrgenes = 5000,							# the number of genes the simulated experiment will have (will be cropped if it exceeds the length of tnumber)
		mediantime = 12,						# the median of the half life time distribution
		ccl = 150,								# the cell cycle length (in minutes)
		delaytime = 0,							# a time which estimates the delay between the moment of thio-Uridine labeling and actual incorporation of it into mRNA
		sdnoise = 0.075,						# the amount of measurement noise (proportional to expression strength)
		nobias = FALSE,							# should a labeling bias be added
		unspecific.LtoU = 0,					# proportion of labeled RNAs that unspecifically end up in the unlabeled fraction
		unspec.LtoU.weighted = FALSE,			# do unspecific proportion of labeled to unlabeled depend linearly on the length of the RNA
		unspecific.UtoL = 0,					# proportion of unlabeled RNAs that unspecifically end up in the labeled fraction
		unspec.UtoL.weighted = FALSE,			# do unspecific proportion of unlabeled to labeled depend linearly on the length of the RNA
		truehalflives = NULL,					# if the data should be generated using a given half life distribution, this vector must contain the respective values for each gene
		truecomplete = NULL,					# if the data should be generated using a given expression distribution, this vector must contain the respective values for each gene
		genenames = NULL						# an optional list of gene names
)
{
	
	### PRELIMINARIES ###	
	
	timepoint = timepoint-delaytime
	alpha = log(2)/ccl
	truelambdas = log(2)/truehalflives
	
	### THE "TRUE" AMOUNT OF TOTAL RNA (Cgrt) ###
	
	Cgr0 = truecomplete
	Cgrt = Cgr0*exp(alpha*timepoint)
	names(Cgrt) = genenames
	
	### THE "TRUE" AMOUNT OF NEWLY SYNTHESISED RNA (Agrt) ###
	
	synth = (exp(alpha*timepoint)-exp(-truelambdas*timepoint))
	Agrt = Cgr0 * synth
	names(Agrt) = genenames
	
	### THE "TRUE" AMOUNT OF PRE-EXISTING RNA (Bgrt) ###
	
	Bgrt = Cgrt - Agrt
	
	### LABELED RNA THAT UNSPECIFICALLY ENDS UP IN UNLABELED ### 
	
	if (unspec.LtoU.weighted){weights.LtoU = tnumber/max(tnumber) - median(tnumber/max(tnumber)) + 1} else {weights.LtoU = 1}
	unspecL = Agrt * unspecific.LtoU * weights.LtoU
	
	### UNLABELED RNA THAT UNSPECIFICALLY ENDS UP IN LABELED ### 
	
	if (unspec.UtoL.weighted){weights.UtoL = tnumber/max(tnumber) - median(tnumber/max(tnumber)) + 1} else {weights.UtoL = 1}
	unspecU = Bgrt * unspecific.UtoL * weights.UtoL
	
	### THE AMOUNT OF MEASURED TOTAL RNA (Tg) ###
	
	Tg = rnorm(length(Cgrt),mean=Cgrt,sd=Cgrt*sdnoise)
	names(Tg) = genenames
	
	truecr = median(Tg)/median(Cgrt)
	
	### THE AMOUNT OF MEASURED LABELED RNA (Lg) ###
	
	if (nobias){labeleff = rep(1,length(tnumber))} else {labeleff =  bias(plabel,tnumber)}
	Lbetween = Agrt * labeleff + unspecU
	
	Lg = rnorm(length(Lbetween),mean=Lbetween,sd= Lbetween*sdnoise)
	Lg = Lg/median(Lg)*median(Tg)
	names(Lg) = genenames
	
	truear = median(Lg)/median(Lbetween)
	
	### THE AMOUNT OF MEASURED UNLABELED RNA (Ug) ###
	
	Utrue = Cgrt + unspecL - Lbetween
	Ug = rnorm(length(Utrue),mean=Utrue,sd=Utrue*sdnoise)
	Ug = Ug/median(Ug)*median(Tg)
	names(Ug) = genenames
	
	truebr = median(Ug)/median(Utrue)
	
	### THE TRUE PROPORTIONALITY CONSTANTS ###
	
	truecrbyar = truecr/truear #median(Agrt)/median(Cgrt)
	truecrbybr = truecr/truebr #median(Bgrt)/median(Cgrt)
	truebrbyar = median(Agrt)/median(Bgrt) #truebr/truear
	trueLasymptote = log(truecrbyar*median(Agrt/truecomplete))
	trueUasymptote = log(truecrbybr*median(Bgrt/truecomplete))
	
	### RETURN RESULTS IN A LIST ###
	
	results = list()
	results[["truecomplete"]] = Cgrt				# the true total RNA amount (C_gr in the paper)
	results[["after"]] = Agrt						# the true newly synthesized RNA (A_gr in the paper)
	results[["before"]] = Bgrt						# the true amount of RNA synthesized before labeling (B_gr)
	results[["total"]] = Tg							# the "measured" total RNA (T_gr)
	results[["labeled"]] = Lg						# the "measured" labeled RNA (L_gr)
	results[["unlabeled"]] = Ug						# the "measured" unlabeled RNA (U_gr)
	results[["truear"]] = truear					# the proportionality constants
	results[["truebr"]] = truebr					# the proportionality constants
	results[["truecr"]] = truecr					# the proportionality constants
	results[["truecrbyar"]] = truecrbyar			# the true quotient truecr/truear
	results[["truecrbybr"]] = truecrbybr			# the true quotient truecr/truebr
	results[["truebrbyar"]] = truebrbyar			# the true quotient truebr/truear
	results[["trueLasymptote"]] = trueLasymptote	# the true Lasymptote
	results[["trueUasymptote"]] = trueUasymptote	# the true Uasymptote
	return(results)
}


### DTA.generate takes a vector of timepoints and the parameters needed for DTA.singlegenerate and produces the according phenotype matrix and data matrix ###


DTA.generate = function(timepoints, 	# the timepoints at which the samples are taken
		tnumber = NULL,   				# the number of uridine residues for each gene
		plabel = NULL,					# the efficiency with which a Uridine position is finally Biotynilated
		nrgenes = 5000,					# the number of genes the simulated experiment will have (will be cropped if it exceeds the length of tnumber)
		mediantime = 12,				# the median of the half life time distribution
		ccl = 150,						# the cell cycle length (in minutes)
		delaytime = 0,					# a time which estimates the delay between the moment of thio-Uridine labeling and actual incorporation of it into mRNA
		sdnoise = 0.075,				# the amount of measurement noise (proportional to expression strength)
		nobias = FALSE,					# should a labeling bias be added
		unspecific.LtoU = 0,			# proportion of labeled RNAs that unspecifically end up in the unlabeled fraction
		unspec.LtoU.weighted = FALSE,	# do unspecific proportion of labeled to unlabeled depend linearly on the length of the RNA
		unspecific.UtoL = 0,			# proportion of unlabeled RNAs that unspecifically end up in the labeled fraction
		unspec.UtoL.weighted = FALSE,	# do unspecific proportion of unlabeled to labeled depend linearly on the length of the RNA
		truehalflives = NULL,			# if the data should be generated using a given half life distribution, this vector must contain the respective values for each gene
		truecomplete = NULL,			# if the data should be generated using a given expression distribution, this vector must contain the respective values for each gene
		genenames = NULL				# an optional list of gene names
)
{
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
			truehalflives = rf(nrgenes,15,15)*mediantime
		}
	} else {truehalflives = rf(nrgenes,15,15)*mediantime}
	names(truehalflives) = genenames
	truelambdas = log(2)/truehalflives
	
	### THE "TRUE" AMOUNT OF TOTAL RNA (C_gr) ###
	
	if (!is.null(truecomplete)){
		if (length(truecomplete) < nrgenes){
			stop("The number of genes exceeds the length of truecomplete. They will not be used !")
			truecomplete = exp(rnorm(nrgenes,mean=8,sd=1))
		}
	} else {truecomplete = exp(rnorm(nrgenes,mean=8,sd=1))}
	names(truecomplete) = genenames
	
	### SYNTHESIS RATES ###
	
	truemus = truecomplete*(log(2)/ccl + truelambdas)
	
	### PRELIMINARIES ###
	
	nrgenes = length(genenames)
	nrexperiments = length(timepoints)
	truear = numeric(nrexperiments)
	truebr = numeric(nrexperiments)
	truecr = numeric(nrexperiments)
	truecrbyar = numeric(nrexperiments)
	truecrbybr = numeric(nrexperiments)
	truebrbyar = numeric(nrexperiments)
	trueLasymptote = numeric(nrexperiments)
	trueUasymptote = numeric(nrexperiments)
	
	### CREATE PHENOMAT ###
	
	phenomat = DTA.phenomat(timepoints)
	
	### CREATE DATAMAT ###
	
	datamat = matrix(0,ncol=nrexperiments*3,nrow=nrgenes)
	colnames(datamat) = phenomat[,"name"]
	rownames(datamat) = genenames
	
	### RESULTS ###
	
	for (nr in 1:nrexperiments){
		res = DTA.singlegenerate(timepoint=timepoints[nr],tnumber=tnumber,plabel=plabel[nr],nrgenes=nrgenes,mediantime=mediantime,ccl=ccl,
				delaytime=delaytime,sdnoise=sdnoise,nobias=nobias,unspecific.LtoU=unspecific.LtoU,unspec.LtoU.weighted=unspec.LtoU.weighted,
				unspecific.UtoL=unspecific.UtoL,unspec.UtoL.weighted=unspec.UtoL.weighted,
				truehalflives=truehalflives,truecomplete=truecomplete,genenames=genenames)		
		datamat[,nr] = res$total
		datamat[,nr+nrexperiments] = res$labeled
		datamat[,nr+2*nrexperiments] = res$unlabeled
		truear[nr] = res$truear
		truebr[nr] = res$truebr
		truecr[nr] = res$truecr
		truecrbyar[nr] = res$truecrbyar
		truecrbybr[nr] = res$truecrbybr
		truebrbyar[nr] = res$truebrbyar
		trueLasymptote[nr] = res$trueLasymptote
		trueUasymptote[nr] = res$trueUasymptote
	}
	
	### RETURN RESULTS IN A LIST ###
	
	results = list()
	results[["phenomat"]] = phenomat
	results[["datamat"]] = datamat
	results[["tnumber"]] = tnumber
	results[["ccl"]] = ccl
	results[["truecomplete"]] = truecomplete
	results[["truelambdas"]] = truelambdas
	results[["truemus"]] = truemus
	results[["truehalflives"]] = truehalflives
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


