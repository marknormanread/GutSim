# Estimates how much of an effect a particular experiment has on the guild landscapes. 
#
# This is done by taking the difference between the mouse responses under each experiment 
# (hence, it is essential that both datasets contain mice in the same locations), and 
# generating a GAM based on that. The GAM summary stats indicate whether differences
# occur across P, C or P+C dimensions. 
#
# 2016-06-22. Mark N. Read.

# Clean the R environment
rm(list=ls())

options(error=traceback)	# more helpful stacktrace for exceptions and errors
print ('starting...')

before_dir <- '../results/20181202-adlib/mice-data-const/'
after_dir  <- '../results/20190201-pre_0.1_adlib/mice-data-const/'
output_dir <- '../results/20190201-pre_0.1_adlib/adlib_vs_pre_0.1/'


drawXLabel <- FALSE
drawYLabel <- FALSE
drawTitle <- FALSE

source('gam_functions.r')  # Load functions for dealing with GAMs
dir.create(output_dir,showWarnings=FALSE)

before_data  <- read.csv(paste(before_dir, "/mice-data.txt", sep=""), header=1)
# Rename these columns to something shorter. 
names(before_data)[names(before_data) == "carb.intake.KJ.d"] <- "eaten.C"
names(before_data)[names(before_data) == "prot.intake.KJ.d"] <- "eaten.P"
# Create a 2D dataset at which GAM predictions will be plotted (this is very useful if the GAM was created
# using a 3D datatset, in which case a 2D plane through 3D space is needed for plotting. )
plane <- findConvex(before_data$eaten.P, before_data$eaten.C, c('eaten.P', 'eaten.C'))

after_data  <- read.csv(paste(after_dir, "/mice-data.txt", sep=""), header=1)
# Rename these columns to something shorter. 
names(after_data)[names(after_data) == "carb.intake.KJ.d"] <- "eaten.C"
names(after_data)[names(after_data) == "prot.intake.KJ.d"] <- "eaten.P"

# Add columns to this as we go, representing the difference between before and after for given responses.
differences <- data.frame(eaten.C=before_data$eaten.C, eaten.P=before_data$eaten.P)

# Colors should go from blue (low values) to red (high values), with white in between.
rgb.palette <- colorRampPalette(c("blue","white","red"), space="Lab", interpolate="linear")


calculateDifference <- function(response, before, after, differencesDF)
{
	print(paste('calculating difference for ', response, sep=''))
	diff <- after_data[, response] - before_data[, response]	
	differencesDF[, response] <- diff
	return(differencesDF)
}
 

#################### RELATIVE ABUNDANCES ###########################
# subtract 'before' data from 'after'
response <- 'relGrf'
if (response %in% colnames(before_data) && response %in% colnames(after_data))
{
	differences <- calculateDifference(response, before_data, after_data, differences)
	gam1 <- generateGamNoFat(response, differences)
	graphLoc <- paste(output_dir, '/_differential_relAb_Grf.pdf',sep='')
	predictPlotGam(gam1, plane, graphLoc=graphLoc, title='Grf relative abundance differential', whiteCentred=0, drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)
	writeGamSummary(gam1, paste(output_dir, '/__gamSummary_diff_relGrf.txt', sep=''))
}

response <- 'relGrm'
if (response %in% colnames(before_data) && response %in% colnames(after_data))
{
	differences <- calculateDifference(response, before_data, after_data, differences)
	gam1 <- generateGamNoFat(response, differences)
	graphLoc <- paste(output_dir, '/_differential_relAb_Grm.pdf',sep='')
	predictPlotGam(gam1, plane, graphLoc=graphLoc, title='Grm relative abundance differential', whiteCentred=0, drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)
	writeGamSummary(gam1, paste(output_dir, '/__gamSummary_diff_relGrm.txt', sep=''))
}

response <- 'relGpf'
if (response %in% colnames(before_data) && response %in% colnames(after_data))
{
	differences <- calculateDifference(response, before_data, after_data, differences)
	gam1 <- generateGamNoFat(response, differences)
	graphLoc <- paste(output_dir, '/_differential_relAb_Gpf.pdf',sep='')
	predictPlotGam(gam1, plane, graphLoc=graphLoc, title='Gpf relative abundance differential', whiteCentred=0, drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)
	writeGamSummary(gam1, paste(output_dir, '/__gamSummary_diff_relGpf.txt', sep=''))
}

response <- 'relGpm'
if (response %in% colnames(before_data) && response %in% colnames(after_data))
{
	differences <- calculateDifference(response, before_data, after_data, differences)
	gam1 <- generateGamNoFat(response, differences)
	graphLoc <- paste(output_dir, '/_differential_relAb_Gpm.pdf',sep='')
	predictPlotGam(gam1, plane, graphLoc=graphLoc, title='Gpm relative abundance differential', whiteCentred=0, drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)
	writeGamSummary(gam1, paste(output_dir, '/__gamSummary_diff_relGpm.txt', sep=''))
}

response <- 'relGmf'
if (response %in% colnames(before_data) && response %in% colnames(after_data))
{
	differences <- calculateDifference(response, before_data, after_data, differences)
	gam1 <- generateGamNoFat(response, differences)
	graphLoc <- paste(output_dir, '/_differential_relAb_Gmf.pdf',sep='')
	predictPlotGam(gam1, plane, graphLoc=graphLoc, title='Gmf relative abundance differential', whiteCentred=0, drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)
	writeGamSummary(gam1, paste(output_dir, '/__gamSummary_diff_relGmf.txt', sep=''))
}

response <- 'relGmm'
if (response %in% colnames(before_data) && response %in% colnames(after_data))
{
	differences <- calculateDifference(response, before_data, after_data, differences)
	gam1 <- generateGamNoFat(response, differences)
	graphLoc <- paste(output_dir, '/_differential_relAb_Gmm.pdf',sep='')
	predictPlotGam(gam1, plane, graphLoc=graphLoc, title='Gmm relative abundance differential', whiteCentred=0, drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)
	writeGamSummary(gam1, paste(output_dir, '/__gamSummary_diff_relGmm.txt', sep=''))
}

#################### ABSOLUTE ABUNDANCES ###########################
# subtract 'before' data from 'after'
response <- 'absGrf'
if (response %in% colnames(before_data) && response %in% colnames(after_data))
{
	differences <- calculateDifference(response, before_data, after_data, differences)
	gam1 <- generateGamNoFat(response, differences)
	graphLoc <- paste(output_dir, '/_differential_abs_Grf.pdf',sep='')
	predictPlotGam(gam1, plane, graphLoc=graphLoc, title='Grf absolute abundance differential', whiteCentred=0, drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)
	writeGamSummary(gam1, paste(output_dir, '/__gamSummary_diff_absGrf.txt', sep=''))
}

response <- 'absGrm'
if (response %in% colnames(before_data) && response %in% colnames(after_data))
{
	differences <- calculateDifference(response, before_data, after_data, differences)
	gam1 <- generateGamNoFat(response, differences)
	graphLoc <- paste(output_dir, '/_differential_abs_Grm.pdf',sep='')
	predictPlotGam(gam1, plane, graphLoc=graphLoc, title='Grm absolute abundance differential', whiteCentred=0, drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)
	writeGamSummary(gam1, paste(output_dir, '/__gamSummary_diff_absGrm.txt', sep=''))
}

response <- 'absGpf'
if (response %in% colnames(before_data) && response %in% colnames(after_data))
{
	differences <- calculateDifference(response, before_data, after_data, differences)
	gam1 <- generateGamNoFat(response, differences)
	graphLoc <- paste(output_dir, '/_differential_abs_Gpf.pdf',sep='')
	predictPlotGam(gam1, plane, graphLoc=graphLoc, title='Gpf absolute abundance differential', whiteCentred=0, drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)
	writeGamSummary(gam1, paste(output_dir, '/__gamSummary_diff_absGpf.txt', sep=''))
}

response <- 'absGpm'
if (response %in% colnames(before_data) && response %in% colnames(after_data))
{
	differences <- calculateDifference(response, before_data, after_data, differences)
	gam1 <- generateGamNoFat(response, differences)
	graphLoc <- paste(output_dir, '/_differential_abs_Gpm.pdf',sep='')
	predictPlotGam(gam1, plane, graphLoc=graphLoc, title='Gpm absolute abundance differential', whiteCentred=0, drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)
	writeGamSummary(gam1, paste(output_dir, '/__gamSummary_diff_absGpm.txt', sep=''))
}

response <- 'absGmf'
if (response %in% colnames(before_data) && response %in% colnames(after_data))
{
	differences <- calculateDifference(response, before_data, after_data, differences)
	gam1 <- generateGamNoFat(response, differences)
	graphLoc <- paste(output_dir, '/_differential_abs_Gmf.pdf',sep='')
	predictPlotGam(gam1, plane, graphLoc=graphLoc, title='Gmf absolute abundance differential', whiteCentred=0, drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)
	writeGamSummary(gam1, paste(output_dir, '/__gamSummary_diff_absGmf.txt', sep=''))
}

response <- 'absGmm'
if (response %in% colnames(before_data) && response %in% colnames(after_data))
{
	differences <- calculateDifference(response, before_data, after_data, differences)
	gam1 <- generateGamNoFat(response, differences)
	graphLoc <- paste(output_dir, '/_differential_abs_Gmm.pdf',sep='')
	predictPlotGam(gam1, plane, graphLoc=graphLoc, title='Gmm absolute abundance differential', whiteCentred=0, drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)
	writeGamSummary(gam1, paste(output_dir, '/__gamSummary_diff_absGmm.txt', sep=''))
}

#################### PROPORTION NITROGEN LIMITED ###########################
# subtract 'before' data from 'after'
response <- 'nLimGrf'
if (response %in% colnames(before_data) && response %in% colnames(after_data))
{
	differences <- calculateDifference(response, before_data, after_data, differences)
	gam1 <- generateGamNoFat(response, differences)
	graphLoc <- paste(output_dir, '/_differential_nLim_Grf.pdf',sep='')
	predictPlotGam(gam1, plane, graphLoc=graphLoc, title='Grf proportion N limited differential', whiteCentred=0, drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)
	writeGamSummary(gam1, paste(output_dir, '/__gamSummary_diff_nLimGrf.txt', sep=''))
}

response <- 'nLimGrm'
if (response %in% colnames(before_data) && response %in% colnames(after_data))
{
	differences <- calculateDifference(response, before_data, after_data, differences)
	gam1 <- generateGamNoFat(response, differences)
	graphLoc <- paste(output_dir, '/_differential_nLim_Grm.pdf',sep='')
	predictPlotGam(gam1, plane, graphLoc=graphLoc, title='Grm proportion N limited differential', whiteCentred=0, drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)
	writeGamSummary(gam1, paste(output_dir, '/__gamSummary_diff_nLimGrm.txt', sep=''))
}

response <- 'nLimGpf'
if (response %in% colnames(before_data) && response %in% colnames(after_data))
{
	differences <- calculateDifference(response, before_data, after_data, differences)
	gam1 <- generateGamNoFat(response, differences)
	graphLoc <- paste(output_dir, '/_differential_nLim_Gpf.pdf',sep='')
	predictPlotGam(gam1, plane, graphLoc=graphLoc, title='Gpf proportion N limited differential', whiteCentred=0, drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)
	writeGamSummary(gam1, paste(output_dir, '/__gamSummary_diff_nLimGpf.txt', sep=''))
}

response <- 'nLimGpm'
if (response %in% colnames(before_data) && response %in% colnames(after_data))
{
	differences <- calculateDifference(response, before_data, after_data, differences)
	gam1 <- generateGamNoFat(response, differences)
	graphLoc <- paste(output_dir, '/_differential_nLim_Gpm.pdf',sep='')
	predictPlotGam(gam1, plane, graphLoc=graphLoc, title='Gpm proportion N limited differential', whiteCentred=0, drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)
	writeGamSummary(gam1, paste(output_dir, '/__gamSummary_diff_nLimGpm.txt', sep=''))
}

response <- 'nLimGmf'
if (response %in% colnames(before_data) && response %in% colnames(after_data))
{
	differences <- calculateDifference(response, before_data, after_data, differences)
	gam1 <- generateGamNoFat(response, differences)
	graphLoc <- paste(output_dir, '/_differential_nLim_Gmf.pdf',sep='')
	predictPlotGam(gam1, plane, graphLoc=graphLoc, title='Gmf proportion N limited differential', whiteCentred=0, drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)
	writeGamSummary(gam1, paste(output_dir, '/__gamSummary_diff_nLimGmf.txt', sep=''))
}

response <- 'nLimGmm'
if (response %in% colnames(before_data) && response %in% colnames(after_data))
{
	differences <- calculateDifference(response, before_data, after_data, differences)
	gam1 <- generateGamNoFat(response, differences)
	graphLoc <- paste(output_dir, '/_differential_nLim_Gmm.pdf',sep='')
	predictPlotGam(gam1, plane, graphLoc=graphLoc, title='Gmm proportion N limited differential', whiteCentred=0, drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)
	writeGamSummary(gam1, paste(output_dir, '/__gamSummary_diff_nLimGmm.txt', sep=''))
}

#################### MISCELANEOUS ITEMS ###########################
# subtract 'before' data from 'after'
response <- 'total'
if (response %in% colnames(before_data) && response %in% colnames(after_data))
{
	differences <- calculateDifference(response, before_data, after_data, differences)
	gam1 <- generateGamNoFat(response, differences)
	graphLoc <- paste(output_dir, '/_differential_total.pdf',sep='')
	predictPlotGam(gam1, plane, graphLoc=graphLoc, title='Total load differential', whiteCentred=0, drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)
	writeGamSummary(gam1, paste(output_dir, '/__gamSummary_diff_total.txt', sep=''))
}

response <- 'dead'
if (response %in% colnames(before_data) && response %in% colnames(after_data))
{
	differences <- calculateDifference(response, before_data, after_data, differences)
	gam1 <- generateGamNoFat(response, differences)
	graphLoc <- paste(output_dir, '/_differential_dead.pdf',sep='')
	predictPlotGam(gam1, plane, graphLoc=graphLoc, title='Total dead differential', whiteCentred=0, drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)
	writeGamSummary(gam1, paste(output_dir, '/__gamSummary_diff_dead.txt', sep=''))
}

response <- 'resist'
if (response %in% colnames(before_data) && response %in% colnames(after_data))
{
	differences <- calculateDifference(response, before_data, after_data, differences)
	gam1 <- generateGamNoFat(response, differences)
	graphLoc <- paste(output_dir, '/_differential_resist.pdf',sep='')
	predictPlotGam(gam1, plane, graphLoc=graphLoc, title='Total resistant differential', whiteCentred=0, drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)
	writeGamSummary(gam1, paste(output_dir, '/__gamSummary_diff_resist.txt', sep=''))
}

response <- 'stress'
if (response %in% colnames(before_data) && response %in% colnames(after_data))
{
	differences <- calculateDifference(response, before_data, after_data, differences)
	gam1 <- generateGamNoFat(response, differences)
	graphLoc <- paste(output_dir, '/_differential_stress.pdf',sep='')
	predictPlotGam(gam1, plane, graphLoc=graphLoc, title='Total stressed differential', whiteCentred=0, drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)
	writeGamSummary(gam1, paste(output_dir, '/__gamSummary_diff_stress.txt', sep=''))
}

response <- 'expon'
if (response %in% colnames(before_data) && response %in% colnames(after_data))
{
	differences <- calculateDifference(response, before_data, after_data, differences)
	gam1 <- generateGamNoFat(response, differences)
	graphLoc <- paste(output_dir, '/_differential_exponential.pdf',sep='')
	predictPlotGam(gam1, plane, graphLoc=graphLoc, title='Total exponential differential', whiteCentred=0, drawXLabel=drawXLabel, drawYLabel=drawYLabel, drawTitle=drawTitle)
	writeGamSummary(gam1, paste(output_dir, '/__gamSummary_diff_exponential.txt', sep=''))
}