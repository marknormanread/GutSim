# Contains functions for computing generalised additive models (GAM), and predicting 
# values at particular values. 
#
# Include this in other scripts by calling "source('gam_functions.r')"



require(mgcv)
require(sp)

no.cols<-256			# number of colors
gr <- 101				# resolution along each axis. 
gamma <- 1.4
findConvex<-function(x,y,names){
	hull<-cbind(x,y)[chull(cbind(x,y)),]
	px<-pretty(x)			# returns same number of items as x, but rounded and evenly spaced.
	py<-pretty(y)
	x.new<-seq(min(px),max(px),len=gr)
	y.new<-seq(min(py),max(py),len=gr)
	ingrid<-as.data.frame(expand.grid(x.new,y.new))
	Fgrid<-ingrid
	Fgrid[(point.in.polygon(ingrid[,1], ingrid[,2], hull[,1],hull[,2])==0),]<-NA
	names(Fgrid)<-names
	return(Fgrid)
	}

generateGamNoFat <- function(responseName, training, smallK=FALSE)
{
	if (smallK)
	{
	  # Fewer degrees of freedom = smoother graphs. 
	  k1 <- 4  # Max degrees of freedom allowed in the gam predictor terms. 
	  k2 <- 6  
	} else {
    k1 <- 18			# max degrees of freedom allowed in the gam predictor terms. 
	  k2 <- 12
	}
	gam <- gam(eval(parse(text=responseName))
							~ s(eaten.P,k=k1,bs="tp") + s(eaten.C,k=k1,bs="tp")
							+ s(eaten.P,eaten.C,k=k2,bs="tp"),
							gamma=1.0,method="REML",data=training,select=TRUE)
	#print(summary(gam))
	return(gam)
}


generateGam <- function(responseName, training)
{
	k1 <- 18			# max degrees of freedom allowed in the gam predictor terms. 
	k2 <- 12
	gam <- gam(eval(parse(text=responseName))
							~ s(eaten.P,k=k1,bs='tp') + s(eaten.C,k=k1,bs='tp') + s(eaten.F,k=k1,bs='tp')
							+ s(eaten.P,eaten.C,k=k2,bs='tp') + s(eaten.C,eaten.F,k=k2,bs='tp') + s(eaten.P,eaten.F,k=k2,bs='tp')
							+ s(eaten.P,eaten.C,eaten.F,k=k2,bs='tp')
							,gamma=gamma,method="REML",data=training,select=TRUE)
	return(gam)
}

predictPlotGam <- function(gam, newdata, graphLoc, title, whiteCentred=NULL, 
						drawXLabel=TRUE, drawYLabel=TRUE, drawTitle=TRUE, 
						xdots=NULL, ydots=NULL, drawContour=TRUE)
{
	# draws a landscape based on a GAM, at the supplied data points. 
	fit <- predict(gam, newdata=newdata)
	# create the plot. 
	plot2DLandscape(fit, newdata, graphLoc, title, whiteCentred=whiteCentred, xdots=xdots, ydots=ydots, drawContour=drawContour,
	                drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)	
	
	return(fit)
}

plot2DLandscape <- function(fit, newdata, graphLoc, title, whiteCentred=NULL, 
	drawXLabel=TRUE, drawYLabel=TRUE, drawTitle=TRUE, 
	xdots=NULL, ydots=NULL,
	drawContour=TRUE)
{
	# create the plot. 
	nlev <- 8
	map <- rgb.palette(no.cols)	
	mn <- min(c(unlist(fit)),na.rm=TRUE)		# the minimum and maximum values enountered in the fit. 
	mx <- max(c(unlist(fit)),na.rm=TRUE)
	pdf(graphLoc)								# write pdf to this filename. 
	par(cex=1.9) # adjust text size
	par(lwd=1.9) # adjust line width
	if (mn != mx)
	{
		locs<-(range(unlist(fit),na.rm=TRUE)-mn)/(mx-mn)*no.cols
		surf<-matrix(fit,nrow=sqrt(dim(newdata)[1]))	
		px<-pretty(newdata$eaten.P)			# the tick locations, made pretty (0  5 10 15 20 25 30)
		py<-pretty(newdata$eaten.C)	
		x.new<-seq(min(px),max(px),len=gr)	# x and y locations where coloured boxes are to be drawn on the image. 
		y.new<-seq(min(py),max(py),len=gr)
		# image creates a pixelised plot. 
		cols<-map[locs[1]:locs[2]]
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
			highP <- whiteCentred # temp variables. Assigned to response variables such that white lies on 'whiteCentred'
			lowP <- whiteCentred
			mxExcess <- mx - whiteCentred # calcualte whether mn or mx lies further away from 'whiteCentred'
			mnExcess <- whiteCentred - mn 			
			if (mxExcess > mnExcess) # adjust highP and lowP
			{
				highP <- mx
				lowP <- whiteCentred - mxExcess
			}
			if (mnExcess > mxExcess)
			{
				lowP <- mn
				highP <- whiteCentred + mnExcess
			}
			# EXPLAIN, THIS IS COMPLEX. 
			range <- highP-lowP
			lowProp <- (mn-lowP) / range
			highProp <- (mx-lowP) / range
			locsRange <- locs[2] - locs[1]
			lowIndex <- (locsRange * lowProp) + locs[1]		
			highIndex <- (locsRange * highProp) + locs[1]
			cols <- map[lowIndex:highIndex]
		}
		xlab <- ''
		if (drawXLabel)	{ xlab <- "Protein eaten (KJ/day)" }
		ylab <- ''
		if (drawYLabel) { ylab <- "Carbohydrate eaten (KJ/day)" }
		if (! drawTitle) { title <- '' }
		image(x.new,y.new,surf,col=cols,xlab=xlab,ylab=ylab, main=title,axes=FALSE)
		axis(1, lwd=1.9)  # Draw the axes
		axis(2, lwd=1.9)
		if (drawContour)
		{
			contour(x.new, y.new, surf, add=TRUE, levels=pretty(range(mn,mx),8), labcex=1.9,
			        drawlabels=TRUE, method="flattest")  # labcex adjusts font size for contour lines. 
		}
		if (!is.null(xdots) && !is.null(ydots))
		{
			points(xdots,ydots, pch=16, cex=1.0)	
		}
		dev.off()		
	} else {					# finish writing data to pdf file. 
		print('WARNING! fitted data is flat, min value == max value in landscape. Cannot draw landscape of this.')
	}
}

writeGamSummary <- function(gam, fileLoc)
{
	sink(file=fileLoc)
	print(summary(gam))
	sink()
}