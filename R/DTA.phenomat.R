

### create a phenomat that suits your experiment ###


DTA.phenomat = function(timepoints,		# the timepoints at which the samples are taken
	timecourse = NULL
)
{
	### PRELIMINARIES ###
	
	nrexperiments = length(timepoints)
	
	### CREATE PHENOMAT ###
	
	if (!is.null(timecourse)){
		phenomat = matrix(0,ncol=5,nrow=3*length(timepoints))
		colnames(phenomat) = c("name","fraction","time","nr","timecourse")
		phenomat[,"fraction"] = rep(c("T","L","U"),rep(nrexperiments,3))
		phenomat[,"time"] = rep(timepoints,3)
		phenomat[,"nr"] = rep(1:nrexperiments,3)
		phenomat[,"timecourse"] = rep(timecourse,3)
		phenomat[,"name"] = paste(phenomat[,"fraction"],phenomat[,"time"],phenomat[,"nr"],phenomat[,"timecourse"],sep="-")
		rownames(phenomat) = phenomat[,"name"]
	} else {
		phenomat = matrix(0,ncol=4,nrow=3*length(timepoints))
		colnames(phenomat) = c("name","fraction","time","nr")
		phenomat[,"fraction"] = rep(c("T","L","U"),rep(nrexperiments,3))
		phenomat[,"time"] = rep(timepoints,3)
		phenomat[,"nr"] = rep(1:nrexperiments,3)
		phenomat[,"name"] = paste(phenomat[,"fraction"],phenomat[,"time"],phenomat[,"nr"],sep="-")
		rownames(phenomat) = phenomat[,"name"]
	}
	
	### RETURN CREATED PHENOMAT ###
	
	return(phenomat)
}



