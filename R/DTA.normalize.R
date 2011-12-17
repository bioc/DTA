

### DTA.normalize ###


DTA.normalize = function(mat,	# matrix
		reliable = NULL, 		# the rows to be used
		logscale = FALSE, 		# is the matrix in log-scale ?
		protocol = FALSE, 		# should a protocol be printed ?
		center = FALSE	 		# should the center be 0 (log-scale) or 1 (absolute scale)
	)
{
	return(medctr(mat,userows = reliable,logscale = logscale,protocol = protocol,center = center))	
}



