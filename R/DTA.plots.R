

### plotting parameter setting function ###


parfkt = function(parsetting,x)
{
	switch(parsetting,
			"default" = return({if (x == 1){scex = 0.75} else if (x == 2){scex = 0.85} else if (x > 2){scex = 1.15} else if (x > 6){scex = 1.75};par(mar = c(5,4,4,2)+0.1+1);par(mai = c(1.1,1.1,1.3,0.7));par(cex = 1);par(cex.axis = 1.5*scex)}),
			"scatterplot3d" = return(c(4,3,4,2)+0.1+1),
			"rankpairs" = return(par(cex.axis = 1.5)),
			"correlationleft" = return(par(mar = c(6,5,5,1)+0.1)),
			"correlationright" = return(par(mar = c(6,1,5,3)+0.1))
	)
}


### plot labels function ###


mtextfkt = function(mtextsetting,x,main,xlab = NULL,ylab = NULL,sub = NULL,zlab = NULL)
{
	switch(mtextsetting,
			"default" = return({if (x == 1){scex = 0.75} else if (x == 2){scex = 0.85} else if (x > 2){scex = 1.15} else if (x > 6){scex = 1.85};mtext(xlab,1,3,cex=1.65*scex);mtext(ylab,2,3,cex=1.65*scex);mtext(main,3,3,cex=2*scex);mtext(sub,3,1,col="darkgrey",cex=1.25*scex)}),
			"benchmark" = return({if (x == 1){scex = 0.75} else if (x == 2){scex = 0.85} else if (x > 2){scex = 1.15} else if (x > 6){scex = 1.85};mtext(xlab,1,3,cex=1.65*scex);mtext(ylab,2,3,cex=1.65*scex);mtext(main,3,3,cex=1.65*scex);mtext(sub,3,1,col="darkgrey",cex=1.25*scex)}),
			"scatterplot3d" = return({if (x == 1){scex = 0.75} else if (x == 2){scex = 0.85} else if (x > 2){scex = 1.15} else if (x > 6){scex = 1.85};mtext(xlab,1,2,cex=1.65*scex);mtext(ylab,2,2,cex=1.65*scex);mtext(zlab,4,1,adj = 0,cex=1.65*scex);mtext(main,3,3,cex=2*scex);mtext(sub,3,1,col="darkgrey",cex=1.25*scex)}),
			"rankpairs" = return({if (x == 1){scex = 0.75} else if (x == 2){scex = 0.85} else if (x > 2){scex = 1.15} else if (x > 6){scex = 1.85};mtext(main,3,2,cex=1.25*scex)})
	)
}


### plotting wrapper function ###


DTA.plot.it = function(filename, 		# name of the plot to be saved with the format type suffix
		sw = 1, 						# scaling factor of weight
		sh = 1, 						# scaling factor of height
		sres = 1, 						# scaling factor of the resolution
		plotsfkt, 						# list of plots to be plotted
		ww = 7, 						# width of window
		wh = 7, 						# height of window
		pointsize = 12,     			# the default pointsize of plotted text, interpreted as big points (1/72 inch) for plots to be saved
		dev.pointsize = 8,				# pointsize of plotted text, interpreted as big points (1/72 inch) for display in R
		paper = "special",   			# needed only if filformat = "pdf" or "ps"
		quality = 100,     				# needed only if filformat = "jpeg"
		units = "px",        			# needed only if filformat = "jpeg", "png", "bmp" or "tiff"
		bg = "white", 					# backgroundcolor
		fileformat = "jpeg", 			# save the plot as jpeg, png, bmp, tiff, ps or pdf
		saveit = FALSE, 				# should plot be saved
		notinR = FALSE, 				# should plot be not plotted in R
		RStudio = FALSE,				# for RStudio users
		addformat = NULL 				# should plot be saved additionally in another format
)
{
	pwidth = sw*480
	pheight = sh*480
	pres = sres*72
	if (saveit){
		switch(fileformat,"jpeg" = jpeg(filename = paste(filename,".",fileformat,sep=""), width = pwidth, height = pheight,units = units, pointsize = pointsize, quality = quality, bg = bg,res = pres),
				"png" = png(filename = paste(filename,".",fileformat,sep=""), width = pwidth, height = pheight,units = units, pointsize = pointsize, bg = bg,res = pres),
				"bmp" = bmp(filename = paste(filename,".",fileformat,sep=""), width = pwidth, height = pheight,units = units, pointsize = pointsize, bg = bg,res = pres),
				"tiff" = tiff(filename = paste(filename,".",fileformat,sep=""), width = pwidth, height = pheight,units = units, pointsize = pointsize, bg = bg,res = pres),
				"ps" = postscript(file = paste(filename,".",fileformat,sep=""), width = ww, height = wh, pointsize = pointsize, paper = paper),
				"pdf" = pdf(file = paste(filename,".",fileformat,sep=""), width = ww, height = wh, pointsize = pointsize, paper = paper))
		plotsfkt()
		dev.off()
	}
	if (saveit){
		if (!is.null(addformat)){
			switch(addformat,"jpeg" = jpeg(filename = paste(filename,".",addformat,sep=""), width = pwidth, height = pheight,units = units, pointsize = pointsize, quality = quality, bg = bg,res = pres),
					"png" = png(filename = paste(filename,".",addformat,sep=""), width = pwidth, height = pheight,units = units, pointsize = pointsize, bg = bg,res = pres),
					"bmp" = bmp(filename = paste(filename,".",addformat,sep=""), width = pwidth, height = pheight,units = units, pointsize = pointsize, bg = bg,res = pres),
					"tiff" = tiff(filename = paste(filename,".",addformat,sep=""), width = pwidth, height = pheight,units = units, pointsize = pointsize, bg = bg,res = pres),
					"ps" = postscript(file = paste(filename,".",addformat,sep=""), width = ww, height = wh, pointsize = pointsize, paper = paper),
					"pdf" = pdf(file = paste(filename,".",addformat,sep=""), width = ww, height = wh, pointsize = pointsize, paper = paper))
			plotsfkt()
			dev.off()
		}
	}
	if (!notinR){
		if (!RStudio){dev.new(width = ww,height = wh,pointsize = dev.pointsize)}
		plotsfkt()
	}
}



