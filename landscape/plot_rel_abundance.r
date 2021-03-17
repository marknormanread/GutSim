# Script seeks a file named 'mice-data' in the supplied directory, and draws relative abundance, absolute abundance
# and nitrogen limited landscapes based on that data. 
#
# Note that this script only deals with protein and carbohydrate intake dimensions, fat is not considered (as
# fat does not currently have any input on the model).
#
# Mark N. Read, 2018

# Clean the R environment
rm(list=ls())

require(mgcv)
require(sp)
options(error=traceback)	# More helpful stacktrace for exceptions and errors

setwd("~/Dropbox (Sydney Uni)/projects/gutsim/landscape")
directory <- '../results/20181202-adlib/tmp_test/'

drawTitle <- FALSE
drawXLabel <- FALSE
drawYLabel <- FALSE

rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"), space="Lab", interpolate="linear")


no.cols <- 256
gr <- 101  # Resolution along each axis. 


findConvex <- function(x, y, names)
{
	hull <- cbind(x,y)[chull(cbind(x,y)), ]
	px <- pretty(x)  # returns same number of items as x, but rounded and evenly spaced.
	py <- pretty(y)
	x.new <- seq(min(px),max(px),len=gr)
	y.new <- seq(min(py),max(py),len=gr)
	ingrid <- as.data.frame(expand.grid(x.new,y.new))
	Fgrid <- ingrid
	Fgrid[(point.in.polygon(ingrid[,1], ingrid[,2], hull[,1], hull[,2])==0 ), ] <- NA
	names(Fgrid) <- names
	return(Fgrid)
}


data <- read.csv(paste(directory, "/mice-data.txt", sep=""), header=1)
# Rename these columns to something shorter. 
names(data)[names(data)=="carb.intake.KJ.d"] <- "eaten.C"
names(data)[names(data)=="prot.intake.KJ.d"] <- "eaten.P"
# Create a 2D dataset at which GAM predictions will be plotted (this is very useful if the GAM was created
# using a 3D datatset, in which case a 2D plane through 3D space is needed for plotting. )
plane <- findConvex(data$eaten.P,
					          data$eaten.C,
					          c('eaten.P', 'eaten.C'))

notFlat <- function(data) 
{
	mn = min(c(unlist(data)),na.rm=TRUE)  # The minimum and maximum values enountered in the fit. 
	mx = max(c(unlist(data)),na.rm=TRUE)
	return (mn != mx)
}

# betar assumes response data lies in range 0..100, and will scale this down to 0..1. 
# Note that subsequent model predictions will also lie in range 0..1, and need to be scaled up to their original range. 
generateGamNoFat <- function(responseName, training, betar=FALSE)
{
	# k1 <- 18			# max degrees of freedom allowed in the gam predictor terms. 
	# k2 <- 12
	# fewer degrees of freedom = smoother graphs. 
	k1 <- 4  # Max degrees of freedom allowed in the gam predictor terms. 
	k2 <- 6
	if (betar)
	{
	  # Scale to range [0, 1] for use in betar GAM fitting 
	  # Note that betar does not permit values of exactly 0 or 1, and hence data is adjusted and these are never reached. 
	  # This is also unrealistic for many applications. 
	  training['gam_input'] <- training[responseName] / 100.  
	  gam_fitted <- gam(gam_input
        	             ~ s(eaten.P, k=k1,bs="tp") + s(eaten.C, k=k1, bs="tp")
        	             + s(eaten.P, eaten.C, k=k2, bs="tp"),
        	             gamma=1.0, method="REML", 
        	             family=betar(link="cauchit"), # proportional data (0,1)
        	             data=training, select=TRUE)	  
	} else {
	  gam_fitted <- gam(eval(parse(text=responseName))
          							~ s(eaten.P, k=k1,bs="tp") + s(eaten.C, k=k1, bs="tp")
          							+ s(eaten.P, eaten.C, k=k2, bs="tp"),
          							gamma=1.0, method="REML", data=training, select=TRUE)
	}
	#print(summary(gam_fitted))
	return(gam_fitted)
}


predictPlotGam <- function(
  gam_fitted, newdata, graphLoc, title, whiteCentred=NULL, betar=FALSE,
  colbar_dynamic_range=NULL  # If non-NULL, then vector of two values, (min val, max val) to map colour bar range onto.
  )
{
	# Draws a landscape based on a GAM, at the supplied data points. 
	fit <- predict(gam_fitted, newdata=newdata, type="response")
	if (betar) 
	{
	  fit <- fit * 100
	}
	# Create the plot. 
	plot2DLandscape(fit, newdata, graphLoc, title, whiteCentred=whiteCentred, colbar_dynamic_range=colbar_dynamic_range)	
	return(fit)
}


plot2DLandscape <- function(
  fit, 
  newdata, 
  graphLoc, 
  title, 
  floor_value=NULL,
  num_contour_levels=8,  # Gives some control, but does not actually force this number, nor does it cap it. 
  # whiteCentred and colbar_dynamic_range are currently mutually exclusive in useage. 
  whiteCentred=NULL,  # If non-NULL, then a value where the colour range should be centred to. 
  colbar_dynamic_range=NULL  # If non-NULL, then vector of two values, (min val, max val) to map colour bar range onto.
  )
{
  if (! is.null(whiteCentred) && ! is.null(colbar_dynamic_range)) {
    stop('Cannot both specifiy value to centre white on AND the range over which colours should be plotted; these functions are mutually exclusive.')
  }
	# Create the plot. 
	print(paste('plotting',title))
	if (! is.null(floor_value)) {
	  fit[fit < floor_value] = floor_value
	  colbar_dynamic_range[1] = floor_value
	}
	colour_map <- rgb.palette(no.cols)  # Creates a vector of colour codes
	
	mn_fit <- min(c(unlist(fit)), na.rm=TRUE)	 # the minimum and maximum values enountered in the fit. 
	mx_fit <- max(c(unlist(fit)), na.rm=TRUE)
	if (mn_fit != mx_fit)
	{
	  # Calculates range of values from the data to be mapped onto each colour.
	  # locs <- (range(unlist(fit), na.rm=TRUE) - mn) / (mx-mn)*no.cols  # DEBUG
	  # col_map_index_range <- c(0, no.cols)  
	  # Default range of colours in colour_map that are used. Conditionally modified below. 
	  col_map_index_min <- 0
	  col_map_index_max <- no.cols
	  if (! is.null(colbar_dynamic_range))
	  {
	    # Maps the colour to a value, need not be min and max in supplied fit data (e.g. to achieve similar range of colours across plots)
	    min_colour_map_value <- min(mn_fit, colbar_dynamic_range[1])
	    max_colour_map_value <- max(mx_fit, colbar_dynamic_range[2])
	    colour_map_values_range <- max_colour_map_value - min_colour_map_value
	    
	    # Where does min (or max) value in data to be plotted sit within the range of values colour_map is assigned.
	    # Calculate as proportion, then map onto indexes. 
	    col_map_index_min <- no.cols * (mn_fit - min_colour_map_value) / colour_map_values_range
	    col_map_index_max <-  no.cols * (mx_fit - min_colour_map_value) / colour_map_values_range
	  }
	  
		surf <- matrix(fit, nrow=sqrt(dim(newdata)[1]))	
		px <- pretty(newdata$eaten.P)  # the tick locations, made pretty (0  5 10 15 20 25 30)
		py <- pretty(newdata$eaten.C)	
		x.new <- seq(min(px), max(px), len=gr)  # x and y locations where coloured boxes are to be drawn on the image. 
		y.new <- seq(min(py), max(py), len=gr)
		# image creates a pixelised plot. 
		cols <- colour_map[col_map_index_min : col_map_index_max]
		if (! is.null(whiteCentred))
		{
			# The following will centre the middle color in the heatmap's range of possible colors to the specified number. 
			# The middle color is usually white here. 
			#
			# Colours are indexed by numbers (e.g. 0..256). The maximum and minimum values to be drawn on the image
			# are mapped to these extremes. To centre white on a particular value, need to adjust the indexes that the image commmand
			# can use. 
			#
			# For example, if mn=0% and mx=120%, then white will be set to be 60%. However, if white is instead to be mapped
			# to 50%, then the lowest color index should be mapped to -20%, not 0%. Some of the color index values need to 
			# be removed to compensate for mx lying further away from 50% than mn. 
			#
			# If color indexes were [0..256], then 128 represents 50%. 256 maps to 120%. 0 has been mapped to -20%, however 
			# the smallest value that the image command is going to see is 0%. Hence, range of colors that map can use must
			# start from index 37. (because 120--20 = range of 140%. mn lies 20% above lowP, hence range must start from
			# 20/140=36.6.)
			highP <- whiteCentred		# temp variables. Assigned to response variables such that white lies on 'whiteCentred'
			lowP <- whiteCentred
			mxExcess <- mx_fit - whiteCentred		# calcualte whether mn or mx lies further away from 'whiteCentred'
			mnExcess <- whiteCentred - mn_fit		
			if (mxExcess > mnExcess)			# adjust highP and lowP
			{
				highP <- mx_fit
				lowP <- whiteCentred - mxExcess
			}
			if (mnExcess > mxExcess)
			{
				lowP <- mn_fit
				highP <- whiteCentred + mnExcess
			}
			# EXPLAIN, THIS IS COMPLEX. 
			range <- highP-lowP
			lowProp <- (mn_fit-lowP) / range
			highProp <- (mx_fit-lowP) / range
			locsRange <- col_map_index_max - col_map_index_min
			lowIndex <- (locsRange * lowProp) + col_map_index_min
			highIndex <- (locsRange * highProp) + col_map_index_min
			cols <- colour_map[lowIndex:highIndex]
		}
		pdf(graphLoc)  # write pdf to this filename. 
		par(cex=1.9)  # adjust text size
		par(lwd=1.9)  # adjust line width
		xlab <- ''
		ylab <- ''
		if (drawXLabel) { xlab <- 'Protein eaten (KJ/day)' }
		if (drawYLabel) { ylab <- 'Carbohydrate eaten (KJ/day)' }
		if (! drawTitle)  { title <- ''}	
		image(x.new, y.new, surf, col=cols, xlab=xlab, ylab=ylab, main=title, axes=FALSE)	
		axis(1, lwd=1.9)  # draw the axes, adjust line width
		axis(2, lwd=1.9)
		contour(x.new, y.new,surf, add=TRUE, levels=pretty(range(mn_fit,mx_fit), num_contour_levels), labcex=1.9)  # labcex adjusts font size for contour lines. 
		dev.off()  # finish writing data to pdf file. 
	} else {  # finish writing data to pdf file. 
		print(paste('WARNING! fitted data is flat, min == max values. Cant draw landscape:',title))
	}
}


writeGamSummary <- function(gam, fileLoc)
{
	sink(file=fileLoc)
	print(summary(gam))
	sink()
}

################### RELATIVE ABUNDANCES ###########################
if (FALSE) # DEBUG
{
  if("relGrf" %in% names(data) && notFlat(data$relGrf))
  {
    gam_fitted <- generateGamNoFat('relGrf', data)
  	fit <- predictPlotGam(gam_fitted, plane, paste(directory,'/_relAb-Grf.pdf',sep=""), title='Relative abundance of Grf')
  	writeGamSummary(gam_fitted, paste(directory,'/_relAb-Grf.txt',sep=""))
  }
  if("relGrm" %in% names(data) && notFlat(data$relGrm))
  {
    gam_fitted <- generateGamNoFat('relGrm', data)
  	fit <- predictPlotGam(gam_fitted, plane, paste(directory,'/_relAb-Grm.pdf',sep=""), title='Relative abundance of Grm')
  	writeGamSummary(gam_fitted, paste(directory,'/_relAb-Grm.summary.txt',sep=""))
  }
  if("relGpf" %in% names(data) && notFlat(data$relGpf))
  {
    gam_fitted <- generateGamNoFat('relGpf', data)
  	fit <- predictPlotGam(gam_fitted, plane, paste(directory,'/_relAb-Gpf.pdf',sep=""), title='Relative abundance of Gpf')
  	writeGamSummary(gam_fitted, paste(directory,'/_relAb-Gpf.summary.txt',sep=""))
  }
  if("relGpm" %in% names(data) && notFlat(data$relGpm))
  {
    gam_fitted <- generateGamNoFat('relGpm', data)
  	fit <- predictPlotGam(gam_fitted, plane, paste(directory,'/_relAb-Gpm.pdf',sep=""), title='Relative abundance of Gpm')
  	writeGamSummary(gam_fitted, paste(directory,'/_relAb-Gpm.summary.txt',sep=""))
  }
  if("relGmm" %in% names(data) && notFlat(data$relGmm))
  {
    gam_fitted <- generateGamNoFat('relGmm', data)
  	fit <- predictPlotGam(gam_fitted, plane, paste(directory,'/_relAb-Gmm.pdf',sep=""), title='Relative abundance of Gmm')
  	writeGamSummary(gam_fitted, paste(directory,'/_relAb-Gmm.summary.txt',sep=""))
  }
  if("relGmf" %in% names(data) && notFlat(data$relGmf))
  {
    gam_fitted <- generateGamNoFat('relGmf', data)
  	fit <- predictPlotGam(gam_fitted, plane, paste(directory,'/_relAb-Gmf.pdf',sep=""), title='Relative abundance of Gmf')
  	writeGamSummary(gam_fitted, paste(directory,'/_relAb-Gmf.summary.txt',sep=""))
  }
  data$relGdietCarb <- data$relGrf + data$relGrm + data$relGpf + data$relGpm
  data$relGhostCarb <- data$relGmm + data$relGmf
  data$relGdietProt <- data$relGrf + data$relGpf + data$relGmf
  data$relGhostProt <- data$relGrm + data$relGpm + data$relGmm
  
  gam_fitted <- generateGamNoFat('relGdietCarb',data)
  fit <- predictPlotGam(gam_fitted, plane, paste(directory,'/_relAb-GdietCarb.pdf',sep=""), title='Relative abundance of dietary carb utilizers')
  writeGamSummary(gam_fitted, paste(directory,'/_relAb-GdietCarb.summary.txt',sep=""))
  
  gam_fitted <- generateGamNoFat('relGhostCarb',data)
  fit <- predictPlotGam(gam_fitted, plane, paste(directory,'/_relAb-GhostCarb.pdf',sep=""), title='Relative abundance of host carb utilizers')
  writeGamSummary(gam_fitted, paste(directory,'/_relAb-GhostCarb.summary.txt',sep=""))
  
  gam_fitted <- generateGamNoFat('relGdietProt',data)
  fit <- predictPlotGam(gam_fitted, plane, paste(directory,'/_relAb-GdietProt.pdf',sep=""), title='Relative abundance of dietary protein utilizers')
  writeGamSummary(gam_fitted, paste(directory,'/_relAb-GdietProt.summary.txt',sep=""))
  
  gam_fitted <- generateGamNoFat('relGhostProt',data)
  fit <- predictPlotGam(gam_fitted, plane, paste(directory,'/_relAb-GhostProt.pdf',sep=""), title='Relative abundance of host protein utilizers')
  writeGamSummary(gam_fitted, paste(directory,'/_relAb-GhostProt.summary.txt',sep=""))
}
#################### ABSOLUTE ABUNDANCES ###########################
# need to generate the GAMs first, and then extrapolate, so we can apply a common mapping of colours to values across all plots. 
# abs_max = max(data[, c('absGrf', 'absGrm', 'absGpf', 'absGpm', 'absGmm', 'absGmf')])
# abs_min = min(data[, c('absGrf', 'absGrm', 'absGpf', 'absGpm', 'absGmm', 'absGmf')])
rgb.palette<-colorRampPalette(c("blue","cyan","yellow", "red"), space="Lab", interpolate="linear", bias=2.5)
no.cols <- 500

min_value = 0
max_value = 0
if("absGrf" %in% names(data) && notFlat(data$absGrf))
{
  gam_fitted <- generateGamNoFat('absGrf',data)
  fit_Grf <- predict(gam_fitted, newdata=plane, type="response")
  min_value <- min(c(unlist(fit_Grf), min_value), na.rm=TRUE)
  max_value <- max(c(unlist(fit_Grf), max_value), na.rm=TRUE)
  writeGamSummary(gam_fitted, paste(directory,'/_absAb-Grf.summary.txt',sep=""))
}
if("absGrm" %in% names(data) && notFlat(data$absGrm))
{
  gam_fitted <- generateGamNoFat('absGrm',data)
  fit_Grm <- predict(gam_fitted, newdata=plane, type="response")
  min_value <- min(c(unlist(fit_Grm), min_value), na.rm=TRUE)
  max_value <- max(c(unlist(fit_Grm), max_value), na.rm=TRUE)
	writeGamSummary(gam_fitted, paste(directory,'/_absAb-Grm.summary.txt',sep=""))
}
if("absGpf" %in% names(data) && notFlat(data$absGpf))
{
  gam_fitted <- generateGamNoFat('absGpf',data)
  fit_Gpf <- predict(gam_fitted, newdata=plane, type="response")
  min_value <- min(c(unlist(fit_Gpf), min_value), na.rm=TRUE)
  max_value <- max(c(unlist(fit_Gpf), max_value), na.rm=TRUE)
	writeGamSummary(gam_fitted, paste(directory, '/_absAb-Gpf.summary.txt', sep=""))
}
if("absGpm" %in% names(data) && notFlat(data$absGpm))
{
  gam_fitted <- generateGamNoFat('absGpm',data)
  fit_Gpm <- predict(gam_fitted, newdata=plane, type="response")
  min_value <- min(c(unlist(fit_Gpm), min_value), na.rm=TRUE)
  max_value <- max(c(unlist(fit_Gpm), max_value), na.rm=TRUE)
	writeGamSummary(gam_fitted, paste(directory, '/_absAb-Gpm.summary.txt', sep=""))
}
if("absGmm" %in% names(data) && notFlat(data$absGmm))
{
  gam_fitted <- generateGamNoFat('absGmm',data)
  fit_Gmm <- predict(gam_fitted, newdata=plane, type="response")
  min_value <- min(c(unlist(fit_Gmm), min_value), na.rm=TRUE)
  max_value <- max(c(unlist(fit_Gmm), max_value), na.rm=TRUE)
	writeGamSummary(gam_fitted, paste(directory, '/_absAb-Gmm.summary.txt', sep=""))
}
if("absGmf" %in% names(data) && notFlat(data$absGmf))
{
  gam_fitted <- generateGamNoFat('absGmf',data)
  fit_Gmf <- predict(gam_fitted, newdata=plane, type="response")
  min_value <- min(c(unlist(fit_Gmf), min_value), na.rm=TRUE)
  max_value <- max(c(unlist(fit_Gmf), max_value), na.rm=TRUE)
	writeGamSummary(gam_fitted, paste(directory, '/_absAb-Gmf.summary.txt', sep=""))
}
if("absGrf" %in% names(data) && notFlat(data$absGrf))
{
  plot2DLandscape(fit_Grf, newdata=plane, 
                  graphLoc=paste0(directory,'/_absAb-Grf.pdf'), 
                  title='Absolute abundance of Grf', 
                  floor_value=0,
                  num_contour_levels=7,
                  colbar_dynamic_range=c(min_value, max_value)
  )
}
if("absGrm" %in% names(data) && notFlat(data$absGrm))
{
  plot2DLandscape(fit_Grm, newdata=plane, 
                  graphLoc=paste0(directory,'/_absAb-Grm.pdf'), 
                  title='Absolute abundance of Grm', 
                  floor_value=0,
                  colbar_dynamic_range=c(min_value, max_value)
  )
}
if("absGpf" %in% names(data) && notFlat(data$absGpf))
{
  plot2DLandscape(fit_Gpf, newdata=plane, 
                  graphLoc=paste0(directory,'/_absAb-Gpf.pdf'), 
                  title='Absolute abundance of Gpf', 
                  floor_value=0,
                  colbar_dynamic_range=c(min_value, max_value)
  )
}
if("absGpm" %in% names(data) && notFlat(data$absGpm))
{
  plot2DLandscape(fit_Gpm, newdata=plane, 
                  graphLoc=paste0(directory,'/_absAb-Gpm.pdf'), 
                  title='Absolute abundance of Gpm', 
                  floor_value=0,
                  colbar_dynamic_range=c(min_value, max_value)
  )
}
if("absGmm" %in% names(data) && notFlat(data$absGmm))
{
  plot2DLandscape(fit_Gmm, newdata=plane, 
                  graphLoc=paste0(directory,'/_absAb-Gmm.pdf'), 
                  title='Absolute abundance of Gmm', 
                  floor_value=0,
                  colbar_dynamic_range=c(min_value, max_value)
  )
}
if("absGmf" %in% names(data) && notFlat(data$absGmf))
{
  plot2DLandscape(fit_Gmf, newdata=plane, 
                  graphLoc=paste0(directory,'/_absAb-Gmf.pdf'), 
                  title='Absolute abundance of Gmf', 
                  floor_value=0,
                  colbar_dynamic_range=c(min_value, max_value)
  )
}

rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"), space="Lab", interpolate="linear")
no.cols <- 256
#################### MISCELANEOUS ###########################
if("total" %in% names(data) && notFlat(data$total))
{
  gam_fitted <- generateGamNoFat('total',data)
	fit <- predictPlotGam(gam_fitted, plane, paste(directory,'/_total.pdf',sep=""), title='Total bacteria abundance')
	writeGamSummary(gam_fitted, paste(directory,'/_total.summary.txt',sep=""))
}
#################### BACTERIAL STATES ###########################

if("dead" %in% names(data) && notFlat(data$dead))
{
  gam_fitted <- generateGamNoFat('dead',data)
	fit <- predictPlotGam(gam_fitted, plane, paste(directory,'/_b_state_dead.pdf',sep=""), title='Rel ab of Dead bacteria')
	writeGamSummary(gam_fitted, paste(directory,'/_b_state_dead.summary.txt',sep=""))
}

if("resist" %in% names(data) && notFlat(data$resist))
{
  gam_fitted <- generateGamNoFat('resist',data)
	fit <- predictPlotGam(gam_fitted, plane, paste(directory,'/_b_state_resistant.pdf',sep=""), title='Rel ab of resistant bacteria')
	writeGamSummary(gam_fitted, paste(directory,'/_b_state_resistant.summary.txt',sep=""))
}
if("stress" %in% names(data) && notFlat(data$stress))
{
  gam_fitted <- generateGamNoFat('stress',data)
	fit <- predictPlotGam(gam_fitted, plane, paste(directory,'/_b_state_stressed.pdf',sep=""), title='Rel ab of stressed bacteria')
	writeGamSummary(gam_fitted, paste(directory,'/_b_state_stressed.summary.txt',sep=""))
}
if("expon" %in% names(data) && notFlat(data$expon))
{
  gam_fitted <- generateGamNoFat('expon',data)
	fit <- predictPlotGam(gam_fitted, plane, paste(directory,'/_b_state_exponential.pdf',sep=""), title='Rel ab of exponential bacteria')
	writeGamSummary(gam_fitted, paste(directory,'/_b_state_exponential.summary.txt',sep=""))
}

#################### PROPORTIONS NITROGEN LIMITED ###########################
rgb.palette <- colorRampPalette(c("red","white","green"),space="Lab",interpolate="linear")
whiteCentre <- 50

if("nLimGrf" %in% names(data) && notFlat(data$nLimGrf))
{
  # Leaving the betar functionality in place as an example for possible future use. 
  # It concerns the link function used in the GAM. 
  # If using a gaussian family, then values <0% and >100% are permitted, which is unrealistic. 
  # Potential solution is to scale the response variable to [0, 1] and use a beta distribution (family=betar). 
  # However, that does not permit values of exactly 0 or 1, which is also unrealistic. 
  # Further, the resultant landscapes look suspiciously flat with areas of very high gradient between them, which I also
  # do not consider representative of the underlying data. 
  # No perfect solution, so sticking with gaussian for now (betar=FALSE).
  betar <- FALSE
	gam_fitted <- generateGamNoFat('nLimGrf', data, betar=betar)
	fit <- predictPlotGam(gam_fitted, plane, paste(directory, '/_nLimitAb-Grf.pdf', sep=""), betar=betar,
	                      title='Proportion of Grf N limited', whiteCentred=whiteCentre)
	writeGamSummary(gam_fitted,paste(directory, '/_nLimitAb-Grf.summary.txt', sep=""))
}
if("nLimGrm" %in% names(data) && notFlat(data$nLimGrm))
{
  gam_fitted <- generateGamNoFat('nLimGrm',data)
	fit <- predictPlotGam(gam_fitted, plane, paste(directory, '/_nLimitAb-Grm.pdf', sep=""), title='Proportion of Grm N limited',whiteCentred=whiteCentre)
	writeGamSummary(gam_fitted, paste(directory, '/_nLimitAb-Grm.summary.txt', sep=""))
}
if("nLimGpf" %in% names(data) && notFlat(data$nLimGpf))
{
  gam_fitted <- generateGamNoFat('nLimGpf',data)
	fit <- predictPlotGam(gam_fitted, plane, paste(directory, '/_nLimitAb-Gpf.pdf', sep=""), title='Proportion of Gpf N limited',whiteCentred=whiteCentre)
	writeGamSummary(gam_fitted, paste(directory, '/_nLimitAb-Gpf.summary.txt', sep=""))
}
if("nLimGpm" %in% names(data) && notFlat(data$nLimGpm))
{
  gam_fitted <- generateGamNoFat('nLimGpm',data)
	fit <- predictPlotGam(gam_fitted, plane, paste(directory, '/_nLimitAb-Gpm.pdf', sep=""), title='Proportion of Gpm N limited',whiteCentred=whiteCentre)
	writeGamSummary(gam_fitted, paste(directory, '/_nLimitAb-Gpm.summary.txt', sep=""))
}
if("nLimGmm" %in% names(data) && notFlat(data$nLimGmm))
{
  gam_fitted <- generateGamNoFat('nLimGmm',data)
	fit <- predictPlotGam(gam_fitted, plane, paste(directory, '/_nLimitAb-Gmm.pdf', sep=""), title='Proportion of Gmm N limited',whiteCentred=whiteCentre)
	writeGamSummary(gam_fitted, paste(directory, '/_nLimitAb-Gmm.summary.txt', sep=""))
}
if("nLimGmf" %in% names(data) && notFlat(data$nLimGmf))
{
  gam_fitted <- generateGamNoFat('nLimGmf',data)
	fit <- predictPlotGam(gam_fitted, plane, paste(directory, '/_nLimitAb-Gmf.pdf' ,sep=""), title='Proportion of Gmf N limited',whiteCentred=whiteCentre)
	writeGamSummary(gam_fitted, paste(directory, '/_nLimitAb-Gmf.summary.txt', sep=""))
}



