require(mgcv)
require(sp)

rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"),space="Lab",interpolate="linear")


no.cols<-256
gr<-101				# resolution along each axis. 
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

# Get hold of the data set
genus  <- read.csv("phylogenetic-genus.csv",header=1)
family <- read.csv("phylogenetic-family.csv",header=1)
class  <- read.csv("phylogenetic-class.csv",header=1)

# Pre-processing. These are the coordinates at which GAM predictions are extracted, for plotting surfaces. Can't plot surfaces
# for 3D data, so take a slice at median protein/carb/fat values. 
dfgenus <- findConvex(genus$eaten.P, genus$eaten.C, c("eaten.P","eaten.C"))
dfgenus$eaten.F <- median(genus$eaten.F)

dfgenus2 <- findConvex(genus$eaten.P, genus$eaten.F, c("eaten.P","eaten.F"))
dfgenus2$eaten.C <- median(genus$eaten.C)

dfgenus3 <- findConvex(genus$eaten.C, genus$eaten.F, c("eaten.C","eaten.F"))
dfgenus3$eaten.P <- median(genus$eaten.P)

drawResponse <- function(name, training, newdata, graphName, summaryName, protFat=FALSE, protCarb=FALSE) 
{
	# name = name of the taxa, which is the response value column in the input table. 
	# training = the dataframe containing the response data, and the predictor values. 
	# newData = dataframe describing the values at which predictions are to be calcualted for drawing a landscape. 
	# graphname = location of where landscape graph is to be written
	# summaryName = location of where summary of gam quality-of-fit should be written. 
	k1 <- 10
	k2 <- 10
	print(name)
	gam1 <- gam(eval(parse(text=name)) ~ s(eaten.P,k=k1,bs="tp")+s(eaten.C,k=k1,bs="tp")+s(eaten.F,k=k1,bs="tp") 
							+s(eaten.P,eaten.C,k=k2,bs="tp")+s(eaten.P,eaten.F,k=k2,bs="tp") +s(eaten.C,eaten.F,k=k2,bs="tp")
							+s(eaten.P,eaten.C,eaten.F,k=k2,bs="tp"), 
							gamma=1.0,method="REML",data=training,select=TRUE)

	fit1 <- predict(gam1, newdata=newdata)

	# Draw the landscape
	nlev <- 8
	map <- rgb.palette(no.cols)
	sub <- paste("(Fat: ",round(unique(newdata$eaten.F),0),")", sep="")
	if(protFat)
	{
		sub <- paste("(Carb: ",round(unique(newdata$eaten.C),0),")", sep="")
	}
	if(protCarb)
	{
		sub <- paste("(Prot: ",round(unique(newdata$eaten.P),0),")", sep="")	
	}
	mn <- min(c(unlist(fit1)),na.rm=TRUE)		# the minimum and maximum values enountered in the fit. 
	mx <- max(c(unlist(fit1)),na.rm=TRUE)

	pdf(graphName)
	par(cex=1.5)								# make text bigger on plot
	locs<-(range(unlist(fit1),na.rm=TRUE)-mn)/(mx-mn)*no.cols
	surf<-matrix(fit1,nrow=sqrt(dim(newdata)[1]))
	px<-pretty(genus$eaten.P)			# the tick locations, made pretty (0  5 10 15 20 25 30)
	py<-pretty(genus$eaten.C)
	if(protFat)
	{
		py<-pretty(genus$eaten.F)
	}	
	if(protCarb)
	{
		px<-pretty(genus$eaten.C)
		py<-pretty(genus$eaten.F)
	}
	x.new<-seq(min(px),max(px),len=gr)	# x and y locations where coloured boxes are to be drawn. 
	y.new<-seq(min(py),max(py),len=gr)
	xlab <- 'Protein eaten (KJ/day)'
	ylab <- 'Carbohydrate eaten (KJ/day)'
	if (protFat)
	{
		ylab <- 'Fat eaten (KJ/day)'
	}
	if(protCarb)
	{
		xlab <- 'Carbohydrate eaten (KJ/day)'
		ylab <- 'Fat eaten (KJ/day)'
	}
	# image creates a pixelised plot. 
	image(x.new,y.new,surf,col=map[locs[1]:locs[2]],xlab=xlab,ylab=ylab, main=name,sub=sub,axes=FALSE)
	axis(1)								# draw the axes
	axis(2)
	contour(x.new,y.new,surf,add=TRUE,levels=pretty(range(mn,mx),8),labcex=1.5)		# labcex adjusts font size for contour lines. 
	dev.off()

	print(summaryName)
	sink(file=summaryName) 
	print(summary(gam1))
	sink() 
}


######## GENUS #########
# for (j in 20:659)
# {
# 	if (mean(unlist(genus[j])) > 0.5)
# 	{
# 		name <- colnames(genus)[j]
# 		graphLoc <- paste('genus/', name, '.pdf', sep="")
# 		summaryLoc <- paste('genus/', name, '.summary.txt', sep="")
# 		drawResponse(name,genus,dfgenus,graphLoc,summaryLoc)

# 		name <- colnames(genus)[j]
# 		graphLoc <- paste('genus/', name, '-PF.pdf', sep="")
# 		drawResponse(name,genus,dfgenus2,graphLoc,summaryLoc,protFat=TRUE)

# 		name <- colnames(genus)[j]
# 		graphLoc <- paste('genus/', name, '-CF.pdf', sep="")
# 		drawResponse(name,genus,dfgenus3,graphLoc,summaryLoc,protCarb=TRUE)		
# 	}
# }

######### FAMILY #########
# for (j in 20:217)
# {
# 	if (mean(unlist(family[j])) > 0.5)
# 	{
# 		name <- colnames(family)[j]
# 		graphLoc <- paste('family/', name, '.pdf', sep='')
# 		summaryLoc <- paste('family/', name, '.summary.txt', sep='')
# 		drawResponse(name,family,dfgenus,graphLoc,summaryLoc)
		
# 		graphLoc <- paste('family/', name, '-PF.pdf', sep='')
# 		drawResponse(name,family,dfgenus2,graphLoc,summaryLoc,protFat=TRUE)		
		
# 		graphLoc <- paste('family/', name, '-CF.pdf', sep='')
# 		drawResponse(name,family,dfgenus3,graphLoc,summaryLoc,protCarb=TRUE)				
# 	}
# }

######## CLASS #########
for (j in 20:58)
{
	if (mean(unlist(class[j])) > 0.5)
	{
		name <- colnames(class)[j]
		graphLoc <- paste('class/', name, '.pdf', sep="")
		summaryLoc <- paste('class/', name, '.summary.txt', sep="")
		drawResponse(name,class,dfgenus,graphLoc,summaryLoc)

		graphLoc <- paste('class/', name, '-PF.pdf', sep="")
		drawResponse(name,class,dfgenus2,graphLoc,summaryLoc,protFat=TRUE)

		graphLoc <- paste('class/', name, '-CF.pdf', sep="")
		drawResponse(name,class,dfgenus3,graphLoc,summaryLoc,protCarb=TRUE)
	}
}





