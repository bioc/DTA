

### fast version of mapit (sparse) ###


DTA.map.it = function(mat,		# matrix or vector
			map = NULL,			# mapping vector
			check = TRUE		# should check be printed
) 
{
	if (is.null(map)) stop("YOU NEED A MAPPING VECTOR !!!")
	mat = as.matrix(mat)
	mat = aggregate(mat,list(id = map[rownames(mat)]),function(x){median(x,na.rm = TRUE)})
	rownames(mat) = mat[,"id"]
	mat = mat[,setdiff(colnames(mat),"id"),drop = FALSE]
	mat = as.matrix(mat)
	mat = mat[,colnames(mat)]
	if(check){cat(paste(length(map)," identifiers were merged into ",length(unique(map))," identifiers in ",length(unique(map[duplicated(map)]))," occurrences \n"))}
	return(mat)
}	



