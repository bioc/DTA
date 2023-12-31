

### estimate uses an experiment, given by a phenotype matrix, data matrix and the #uridines for each gene to estimate synthesis and decay rate of the genes ###


DTA.dynamic.estimate = function(
		phenomat = NULL,					# phenotype matrix, "nr" should be numbered by experiments not by replicates
		datamat = NULL, 					# data matrix, should only contain the rows of phenomat as columns
		tnumber = NULL, 					# #uridines, should have the rownames of datamat
		ccl = NULL, 						# the cell cycle length of the cells
		mRNAs = NULL,						# estimated number of mRNAs in a cell
		reliable = NULL, 					# vector of reliable genes, which are used for regression
		mediancenter = TRUE, 				# should the data be rescaled to a common median
		usefractions = "LandT",				# from which fractions should the decay rate be calculated: "LandT", "UandT" or "both"
		LtoTratio = NULL, 					# coefficient to rescale L/T
		ratiomethod = "tls", 				# choose the regression method to be used, possible methods are: bias, tls or lm
		largest = 5, 						# percentage of largest residues from the first regression not to be used in the second regression
		weighted = TRUE, 					# should the regression be weighted with 1/(T^2 + median(T))
		relevant = NULL, 					# choose the arrays to be used for halflives calculation, vector due to experiments variable 
		check = TRUE, 						# if check = TRUE, control messages and plots will be generated
		error = TRUE,						# should the measurement error be assessed
		samplesize = 1000,					# error model samplesize for resampling
		confidence.range = c(0.025,0.975),	# confidence region for error model
		bicor = TRUE, 						# should the labeling bias be corrected
		condition = "", 					# to be added to the plotnames
		upper = 700, 						# upper bound for labeling bias estimation
		lower = 500, 						# lower bound for labeling bias estimation
		save.plots = FALSE, 				# if save.plots = TRUE, control plots will be saved
		resolution = 1,						# resolution scaling factor for plotting
		folder = NULL, 						# folder, where to save the plots
		fileformat = "jpeg", 				# save the plot as jpeg, png, bmp, tiff, ps or pdf
		totaloverwt = 1, 					# total mRNA over WT
		sr.vs.dr.folds.lims = c(-5,5),		# limits of the folds plot
		te.vs.to.folds.lims = c(-6,6),		# limits of the folds plot
		robust = FALSE,						# if robust = TRUE, LE resp. LT is chosen instead of sr resp. dr
		clusters = "sr",					# should the dr vs sr folds be plotted with clusters, choose 'sr', 'dr' for cluster selection or 'none' to omit it
		ranktime = NULL,					# at which time should the rankgain be calculated, default is the last column
		upperquant = 0.8,					# upper quantile for cluster selection
		lowerquant = 0.6,					# lower quantile for cluster selection
		notinR = FALSE,						# should plot be not plotted in R
		RStudio = FALSE,					# for RStudio users
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
		truemusaveraged = sim.object$truemusaveraged
		truelambdas = sim.object$truelambdas
		truelambdasaveraged = sim.object$truelambdasaveraged
		truehalflives = sim.object$truehalflives
		truehalflivesaveraged = sim.object$truehalflivesaveraged
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
	if (!is.null(LtoTratio)){if (length(unique(phenomat[,"timecourse"])) != length(LtoTratio)){stop(paste("The number of specified ratios of L to T (LtoTratio) should correspond to the number of timepoints:",length(unique(phenomat[,"timecourse"])),"!"))}}
	if (save.plots){
		if (is.null(folder)){stop("You need to specify the folder, where the plots should be saved !")}
		if (!is.null(folder)){
			if (file.access(folder,0) != 0){stop("The specified folder needs to exist !")}
			if (file.access(folder,2) != 0){stop("The specified folder has no write permission !")}
		}
	}
	if (!any(clusters %in% c("sr","dr","te","to","none"))){stop("clusters need to be 'sr', 'dr', 'te', 'to' or 'none' !")}
	
	### PRELIMINARIES ###
	
	possibles = intersect(rownames(datamat),names(tnumber))
	tnumber = tnumber[possibles]
	datamat = datamat[possibles,]
	reliable = intersect(reliable,possibles)
	nrgenes = nrow(datamat)
	
	labtimes = as.character(sort(as.numeric(unique(phenomat[,"timecourse"]))))
	names(labtimes) = labtimes
	if (!is.null(LtoTratio)){names(LtoTratio) = labtimes}
	nrlabtimes = length(labtimes)
	extracttimes = as.numeric(unique(phenomat[,"time"]))*sort(as.numeric(unique(phenomat[,"timecourse"]))-1)
	
	### BUILD PHENOMATS ###
	
	phenomats = list()
	for (labtime in labtimes){
		phenomats[[labtime]] = phenomat[names(which(phenomat[,"timecourse"] == labtime)),]
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
	
	truemuslist = list()
	truemusaveragedlist = list()
	truelambdaslist = list()
	truelambdasaveragedlist = list()
	truehalfliveslist = list()
	truehalflivesaveragedlist = list()
	trueplabels = list()
	trueLasymptotes = list()
	trueUasymptotes = list()
	truecrbyars = list()
	truecrbybrs = list()
	truebrbyars = list()
	if (simulation){
		for (labtime in labtimes){truemuslist[[labtime]] = truemus[,as.numeric(labtime)]}
		for (labtime in labtimes){truemusaveragedlist[[labtime]] = truemusaveraged[,as.numeric(labtime)]}
		for (labtime in labtimes){truelambdaslist[[labtime]] = truelambdas[,as.numeric(labtime)]}
		for (labtime in labtimes){truelambdasaveragedlist[[labtime]] = truelambdasaveraged[,as.numeric(labtime)]}
		for (labtime in labtimes){truehalfliveslist[[labtime]] = truehalflives[,as.numeric(labtime)]}
		for (labtime in labtimes){truehalflivesaveragedlist[[labtime]] = truehalflivesaveraged[,as.numeric(labtime)]}
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
		if (labtime == labtimes[1]){initials = NULL} else{		
			initiallabtime = labtimes[which(labtimes == labtime)-1]
			initials = as.matrix(datamats[[initiallabtime]][,phenomats[[initiallabtime]][which(phenomats[[initiallabtime]][,"fraction"] == "T"),"name"]])
			if (check){
				cat("Current timecourse nr: ",labtime,"\n")
				cat("Previous timecourse nr: ",initiallabtime,"\n")
			}
		}
		
		res[[labtime]] = DTA.singleestimate(phenomats[[labtime]],datamats[[labtime]],tnumber,labelingtime = as.numeric(unique(phenomats[[labtime]][,"time"])),ccl = ccl,
				mRNAs = mRNAs,reliable = reliable,mediancenter = mediancenter,usefractions = usefractions,ratiomethod = ratiomethod,largest = largest,RStudio = RStudio,
				weighted = weighted,ratio = ratios[[labtime]],relevant = relevants[[labtime]],check = check,dynamic = TRUE,initials = initials,error = error,samplesize = samplesize,confidence.range = confidence.range,
				notinR = notinR,bicor = bicor,condition = condition,timepoint = labtime,upper = upper,lower = lower,save.plots = save.plots,folder = folder,fileformat = fileformat,
				simulation = simulation,truemus = truemuslist[[labtime]],truemusaveraged = truemusaveragedlist[[labtime]],totaloverwt = totaloverwt,resolution = resolution,
				truelambdas = truelambdaslist[[labtime]],truelambdasaveraged = truelambdasaveragedlist[[labtime]],truehalflives = truehalfliveslist[[labtime]],
				truehalflivesaveraged = truehalflivesaveragedlist[[labtime]],trueplabel=trueplabels[[labtime]],trueLasymptote=trueLasymptotes[[labtime]],
				trueUasymptote=trueUasymptotes[[labtime]],truecrbyar=truecrbyars[[labtime]],truecrbybr=truecrbybrs[[labtime]],truebrbyar=truebrbyars[[labtime]])
	}
	
	drmat = cbind()
	for (labtime in labtimes){drmat = cbind(drmat,res[[labtime]]$dr)}
	colnames(drmat) = labtimes
	
	srmat = cbind()
	for (labtime in labtimes){srmat = cbind(srmat,res[[labtime]]$sr)}
	colnames(srmat) = labtimes
	
	TEmat = cbind()
	for (labtime in labtimes){TEmat = cbind(TEmat,res[[labtime]]$TE)}
	colnames(TEmat) = labtimes
	
	LEmat = cbind()
	for (labtime in labtimes){LEmat = cbind(LEmat,res[[labtime]]$LE)}
	colnames(LEmat) = labtimes
	
	LTmat = cbind()
	for (labtime in labtimes){LTmat = cbind(LTmat,-res[[labtime]]$LT)}
	colnames(LTmat) = labtimes
	
	### HALF-LIFE FOLDS VS. SYNTHESIS RATE FOLDS ###
	
    ellipse.coordinates = function(r,scale = c(1,1),centre = c(0,0),level = 0.95,t = sqrt(qchisq(level,2)))
    {
        r = min(max(r,-1),1)
        d = acos(r)
        a = seq(0,2*pi,len = 100)
        matrix(c(t*scale[1]*cos(a + d/2) + centre[1],t*scale[2]*cos(a - d/2) + centre[2]),100,2,dimnames = list(NULL,c("x","y")))
    }
    
	if (check){
			if (robust) {
				srmat = LEmat
				drmat = LTmat
			}
			plotable = rownames(drmat)[!apply(is.na(drmat),1,any)]
			drfc = log2(drmat[plotable,]/drmat[plotable,1])
			srfc = log2(srmat[plotable,]/srmat[plotable,1])
			nrseries = nrlabtimes-1
			
			target.equlibrium = log2(TEmat[plotable,]/TEmat[plotable,1])
			turnover = log2((srmat[plotable,]*drmat[plotable,])/(srmat[plotable,1]*drmat[plotable,1]))
			
			if (clusters == "none"){
				plotsfkt = function(){
					par(mfrow=c(windowxy(nrseries)[1],windowxy(nrseries)[2]))
					parfkt("default",nrseries)
					for (j in 1:nrseries){
						i = inter(drfc[plotable,j+1],srfc[plotable,j+1])
						heatscatter(i$x,i$y,xlim=sr.vs.dr.folds.lims,ylim=sr.vs.dr.folds.lims,xlab="",ylab="",main="",cor=FALSE)
						xlab = expression(paste(log[2]," ( decay rate fold  ",lambda['gr,i']/lambda['gr,1']," )"))
						ylab = expression(paste(log[2]," ( synthesis rate fold  ",mu['gr,i']/mu['gr,1']," )"))
						main = paste(extracttimes[j+1],"min")
						sub = paste("( i =",j+1,")")
						mtextfkt("default",nrseries,main,xlab,ylab,sub)
						gridfkt(lim = sr.vs.dr.folds.lims)}
				}
				DTA.plot.it(filename = paste(folder,"/drfc_vs_srfc_",condition,sep=""),sw = resolution*windowxy(nrseries)[2],sh = resolution*windowxy(nrseries)[1],sres = resolution,plotsfkt = plotsfkt,ww = 7*windowxy(nrseries)[2],wh = 7*windowxy(nrseries)[1],saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)
				
				plotsfkt = function(){
					par(mfrow=c(windowxy(nrseries)[1],windowxy(nrseries)[2]))
					parfkt("default",nrseries)
					for (j in 1:nrseries){
						i = inter(target.equlibrium[plotable,j+1],turnover[plotable,j+1])
						heatscatter(i$x,i$y,xlim=te.vs.to.folds.lims,ylim=te.vs.to.folds.lims,xlab="",ylab="",main="",cor=FALSE)
						gridfkt(lim = te.vs.to.folds.lims)
						xlab = expression(paste("target equilibrium  ",mu['gr,i']/(lambda['gr,i'] + alpha)," fold"))
						ylab = expression(paste("turnover  ",mu['gr,i']*lambda['gr,i']," fold"))
						main = paste(extracttimes[j+1],"min")
						sub = paste("( i =",j+1,")")
						mtextfkt("default",nrseries,main,xlab,ylab,sub)}
				}
				DTA.plot.it(filename = paste(folder,"/target_equilibrium_vs_turnover_",condition,sep=""),sw = resolution*windowxy(nrseries)[2],sh = resolution*windowxy(nrseries)[1],sres = resolution,plotsfkt = plotsfkt,ww = 7*windowxy(nrseries)[2],wh = 7*windowxy(nrseries)[1],saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)
			} else {
				if (clusters == "sr"){
					srranks = apply(srmat[plotable,],2,rank)
					ranksnorm = srranks-srranks[,1]
				}
				if (clusters == "dr"){
					drranks = apply(drmat[plotable,],2,rank)
					ranksnorm = drranks-drranks[,1]
				}
				if (clusters == "te"){
					teranks = apply(target.equlibrium[plotable,],2,rank)
					ranksnorm = teranks-teranks[,1]
				}
				if (clusters == "to"){
					toranks = apply(turnover[plotable,],2,rank)
					ranksnorm = toranks-toranks[,1]
				}
				if (is.null(ranktime)){ranktime = dim(ranksnorm)[2]}
				rankgain = ranksnorm[,ranktime]
				upbound = quantile(abs(rankgain),upperquant)
				lowbound = quantile(abs(rankgain),lowerquant)
				up = names(rankgain[which(rankgain > upbound)])
				down = names(rankgain[which(rankgain < -upbound)])
				downeven = names(rankgain[which(rankgain <= -lowbound & rankgain >= -upbound)])
				upeven = names(rankgain[which(rankgain >= lowbound & rankgain <= upbound)])
				even = names(rankgain[which(rankgain > -lowbound & rankgain < lowbound)])
				genecluster = list()
				genecluster[["up"]] = up
				genecluster[["down"]] = down
				genecluster[["downeven"]] = downeven
				genecluster[["upeven"]] = upeven
				genecluster[["even"]] = even
				res[["genecluster"]] = genecluster
				upcol = "darkred"
				upevencol = "darkorange"
				evencol = "darkgrey"
				downevencol = "purple3"
				downcol = "darkblue"
				
				plotsfkt = function(){
					scalesd = 1
					level = 0.75
					par(mfrow=c(windowxy(nrseries)[1],windowxy(nrseries)[2]))
					parfkt("default",nrseries)
					for (j in 1:nrseries){
						i = inter(drfc[plotable,j+1],srfc[plotable,j+1])
						plot(0,0,xlim=sr.vs.dr.folds.lims,ylim=sr.vs.dr.folds.lims,col="white",xlab="",ylab="",main="")
						xlab = expression(paste(log[2]," ( decay rate fold  ",lambda['gr,i']/lambda['gr,1']," )"))
						ylab = expression(paste(log[2]," ( synthesis rate fold  ",mu['gr,i']/mu['gr,1']," )"))
						main = paste(extracttimes[j+1],"min")
						sub = paste("( i =",j+1,")")
						mtextfkt("default",nrseries,main,xlab,ylab,sub)
						gridfkt(lim = sr.vs.dr.folds.lims)
						points(i$x[intersect(plotable,genecluster[["even"]])],i$y[intersect(plotable,genecluster[["even"]])],col=evencol,pch=20)
						points(i$x[intersect(plotable,genecluster[["downeven"]])],i$y[intersect(plotable,genecluster[["downeven"]])],col=downevencol,pch=20)
						points(i$x[intersect(plotable,genecluster[["upeven"]])],i$y[intersect(plotable,genecluster[["upeven"]])],col=upevencol,pch=20)
						points(i$x[intersect(plotable,genecluster[["down"]])],i$y[intersect(plotable,genecluster[["down"]])],col=downcol,pch=20)
						points(i$x[intersect(plotable,genecluster[["up"]])],i$y[intersect(plotable,genecluster[["up"]])],col=upcol,pch=20)
						x = i$x[intersect(plotable,genecluster[["even"]])]
						y = i$y[intersect(plotable,genecluster[["even"]])]
						points(ellipse.coordinates(scor(x,y,use="na.or.complete"),scale=c(sd(x,na.rm=TRUE)*scalesd,sd(y,na.rm=TRUE)*scalesd),centre=c(mean(x,na.rm=TRUE),mean(y,na.rm=TRUE)),level=level),type = 'l',col=evencol,lwd=3)
						x = i$x[intersect(plotable,genecluster[["downeven"]])]
						y = i$y[intersect(plotable,genecluster[["downeven"]])]
						points(ellipse.coordinates(scor(x,y,use="na.or.complete"),scale=c(sd(x,na.rm=TRUE)*scalesd,sd(y,na.rm=TRUE)*scalesd),centre=c(mean(x,na.rm=TRUE),mean(y,na.rm=TRUE)),level=level),type = 'l',col=downevencol,lwd=3)
						x = i$x[intersect(plotable,genecluster[["upeven"]])]
						y = i$y[intersect(plotable,genecluster[["upeven"]])]
						points(ellipse.coordinates(scor(x,y,use="na.or.complete"),scale=c(sd(x,na.rm=TRUE)*scalesd,sd(y,na.rm=TRUE)*scalesd),centre=c(mean(x,na.rm=TRUE),mean(y,na.rm=TRUE)),level=level),type = 'l',col=upevencol,lwd=3)
						x = i$x[intersect(plotable,genecluster[["down"]])]
						y = i$y[intersect(plotable,genecluster[["down"]])]
						points(ellipse.coordinates(scor(x,y,use="na.or.complete"),scale=c(sd(x,na.rm=TRUE)*scalesd,sd(y,na.rm=TRUE)*scalesd),centre=c(mean(x,na.rm=TRUE),mean(y,na.rm=TRUE)),level=level),type = 'l',col=downcol,lwd=3)
						x = i$x[intersect(plotable,genecluster[["up"]])]
						y = i$y[intersect(plotable,genecluster[["up"]])]
						points(ellipse.coordinates(scor(x,y,use="na.or.complete"),scale=c(sd(x,na.rm=TRUE)*scalesd,sd(y,na.rm=TRUE)*scalesd),centre=c(mean(x,na.rm=TRUE),mean(y,na.rm=TRUE)),level=level),type = 'l',col=upcol,lwd=3)
						if (j == 1){legend("topleft",rev(c(paste("down (",length(intersect(plotable,genecluster[["down"]])),")",sep = ""),paste("downeven (",length(intersect(plotable,genecluster[["downeven"]])),")",sep = ""),paste("even (",length(intersect(plotable,genecluster[["even"]])),")",sep = ""),paste("upeven (",length(intersect(plotable,genecluster[["upeven"]])),")",sep = ""),paste("up (",length(intersect(plotable,genecluster[["up"]])),")",sep = ""))),pt.bg=rev(c(downcol,downevencol,evencol,upevencol,upcol)),col=c("black","black","black","black","black"),bg="white",pch=21,cex=1.25,inset=0.02)}
					}
				}
				DTA.plot.it(filename = paste(folder,"/drfc_vs_srfc_cluster_",condition,sep=""),sw = resolution*windowxy(nrseries)[2],sh = resolution*windowxy(nrseries)[1],sres = resolution,plotsfkt = plotsfkt,ww = 7*windowxy(nrseries)[2],wh = 7*windowxy(nrseries)[1],saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)
				
				plotsfkt = function(){
					scalesd = 1
					level = 0.75
					par(mfrow=c(windowxy(nrseries)[1],windowxy(nrseries)[2]))
					parfkt("default",nrseries)
					for (j in 1:nrseries){
						i = inter(target.equlibrium[plotable,j+1],turnover[plotable,j+1])
						plot(0,0,xlim=te.vs.to.folds.lims,ylim=te.vs.to.folds.lims,col="white",xlab="",ylab="",main="")
						gridfkt(lim = te.vs.to.folds.lims)
						xlab = expression(paste("target equilibrium  ",mu['gr,i']/(lambda['gr,i'] + alpha)," fold"))
						ylab = expression(paste("turnover  ",mu['gr,i']*lambda['gr,i']," fold"))
						main = paste(extracttimes[j+1],"min")
						sub = paste("( i =",j+1,")")
						mtextfkt("default",nrseries,main,xlab,ylab,sub)
						points(i$x[intersect(plotable,genecluster[["even"]])],i$y[intersect(plotable,genecluster[["even"]])],col=evencol,pch=20)
						points(i$x[intersect(plotable,genecluster[["downeven"]])],i$y[intersect(plotable,genecluster[["downeven"]])],col=downevencol,pch=20)
						points(i$x[intersect(plotable,genecluster[["upeven"]])],i$y[intersect(plotable,genecluster[["upeven"]])],col=upevencol,pch=20)
						points(i$x[intersect(plotable,genecluster[["down"]])],i$y[intersect(plotable,genecluster[["down"]])],col=downcol,pch=20)
						points(i$x[intersect(plotable,genecluster[["up"]])],i$y[intersect(plotable,genecluster[["up"]])],col=upcol,pch=20)
						x = i$x[intersect(plotable,genecluster[["even"]])]
						y = i$y[intersect(plotable,genecluster[["even"]])]
						points(ellipse.coordinates(scor(x,y,use="na.or.complete"),scale=c(sd(x,na.rm=TRUE)*scalesd,sd(y,na.rm=TRUE)*scalesd),centre=c(mean(x,na.rm=TRUE),mean(y,na.rm=TRUE)),level=level),type = 'l',col=evencol,lwd=3)
						x = i$x[intersect(plotable,genecluster[["downeven"]])]
						y = i$y[intersect(plotable,genecluster[["downeven"]])]
						points(ellipse.coordinates(scor(x,y,use="na.or.complete"),scale=c(sd(x,na.rm=TRUE)*scalesd,sd(y,na.rm=TRUE)*scalesd),centre=c(mean(x,na.rm=TRUE),mean(y,na.rm=TRUE)),level=level),type = 'l',col=downevencol,lwd=3)
						x = i$x[intersect(plotable,genecluster[["upeven"]])]
						y = i$y[intersect(plotable,genecluster[["upeven"]])]
						points(ellipse.coordinates(scor(x,y,use="na.or.complete"),scale=c(sd(x,na.rm=TRUE)*scalesd,sd(y,na.rm=TRUE)*scalesd),centre=c(mean(x,na.rm=TRUE),mean(y,na.rm=TRUE)),level=level),type = 'l',col=upevencol,lwd=3)
						x = i$x[intersect(plotable,genecluster[["down"]])]
						y = i$y[intersect(plotable,genecluster[["down"]])]
						points(ellipse.coordinates(scor(x,y,use="na.or.complete"),scale=c(sd(x,na.rm=TRUE)*scalesd,sd(y,na.rm=TRUE)*scalesd),centre=c(mean(x,na.rm=TRUE),mean(y,na.rm=TRUE)),level=level),type = 'l',col=downcol,lwd=3)
						x = i$x[intersect(plotable,genecluster[["up"]])]
						y = i$y[intersect(plotable,genecluster[["up"]])]
						points(ellipse.coordinates(scor(x,y,use="na.or.complete"),scale=c(sd(x,na.rm=TRUE)*scalesd,sd(y,na.rm=TRUE)*scalesd),centre=c(mean(x,na.rm=TRUE),mean(y,na.rm=TRUE)),level=level),type = 'l',col=upcol,lwd=3)
						if (j == 1){legend("topleft",rev(c(paste("down (",length(intersect(plotable,genecluster[["down"]])),")",sep = ""),paste("downeven (",length(intersect(plotable,genecluster[["downeven"]])),")",sep = ""),paste("even (",length(intersect(plotable,genecluster[["even"]])),")",sep = ""),paste("upeven (",length(intersect(plotable,genecluster[["upeven"]])),")",sep = ""),paste("up (",length(intersect(plotable,genecluster[["up"]])),")",sep = ""))),pt.bg=rev(c(downcol,downevencol,evencol,upevencol,upcol)),col=c("black","black","black","black","black"),bg="white",pch=21,cex=1.25,inset=0.02)}
					}
				}
				DTA.plot.it(filename = paste(folder,"/target_equilibrium_vs_turnover_cluster_",condition,sep=""),sw = resolution*windowxy(nrseries)[2],sh = resolution*windowxy(nrseries)[1],sres = resolution,plotsfkt = plotsfkt,ww = 7*windowxy(nrseries)[2],wh = 7*windowxy(nrseries)[1],saveit = save.plots,fileformat = fileformat,notinR = notinR,RStudio = RStudio)
			}
	}
	
	### RETURN RESULTS IN A LIST ###
	
	return(res)
}


