# Plots landscapes pertaining to the nutrient inputs and small intestine absorptions of simulated mice.
# Relies on a file that must be generated from the simulation. 


diets <- read.csv("../diet-analysis/nutrient_rel_abundances.csv",header=1)
out_dir <- '../diet-analysis/R/'

drawDots <- FALSE

source('gam_functions.r')

names(diets)[names(diets)=='carb.intake.KJ.d'] <- 'eaten.C'
names(diets)[names(diets)=='prot.intake.KJ.d'] <- 'eaten.P'
names(diets)[names(diets)=='fat.intake.KJ.d']  <- 'eaten.F'

plane <- findConvex(diets$eaten.P, diets$eaten.C, c('eaten.P','eaten.C'))
plane$eaten.F <- median(diets$eaten.F)

rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"),space="Lab",interpolate="linear")


###################################################
# Plot relative abundances of individual carbon/nitrogen sources as proportion of total carbon/nitrogen available in cecum.
gam <- generateGamNoFat('cecum_relCw', diets)
graphLoc <- paste(out_dir, '__cecumRelAb-Cw.pdf', sep='')
if (drawDots)
{
	fit <- predictPlotGam(gam, plane, graphLoc, 'Wheatstarch C as proportion of all carbon (%)', 
	                      xdots=diets$eaten.P, ydots=diets$eaten.C, drawContour=FALSE)
} else {
	fit <- predictPlotGam(gam, plane, graphLoc, 'Wheatstarch C as proportion of all carbon (%)')
}
writeGamSummary(gam, paste(out_dir,'__gamSummary-relCw.txt',sep=''))

gam <- generateGam('cecum_relCd', diets)
graphLoc <- paste(out_dir, '__cecumRelAb-Cd.pdf', sep='')
if (drawDots)
{
	fit <- predictPlotGam(gam, plane, graphLoc, 'Dex WS C as proportion of all carbon(%)', 
	                      xdots=diets$eaten.P, ydots=diets$eaten.C, drawContour=FALSE)
} else {
	fit <- predictPlotGam(gam, plane, graphLoc, 'Dex WS C as proportion of all carbon(%)')
}
writeGamSummary(gam, paste(out_dir,'__gamSummary-relCd.txt',sep=''))

gam <- generateGam('cecum_relCmuc', diets)
graphLoc <- paste(out_dir, '__cecumRelAb-Cmuc.pdf', sep='')
if (drawDots)
{
	fit <- predictPlotGam(gam, plane, graphLoc, 'Mucin as proportion of all carbon (%)', 
	                      xdots=diets$eaten.P, ydots=diets$eaten.C, drawContour=FALSE)
} else { 
	fit <- predictPlotGam(gam, plane, graphLoc, 'Mucin as proportion of all carbon (%)')
}
writeGamSummary(gam, paste(out_dir,'__gamSummary-relCmuc.txt',sep=''))

gam <- generateGam('cecum_relCcas', diets)
graphLoc <- paste(out_dir, '__cecumRelAb-Ccas.pdf', sep='')
if (drawDots)
{
	fit <- predictPlotGam(gam, plane, graphLoc, 'Feed protein as proportion of all carbon (%)', 
	                      xdots=diets$eaten.P, ydots=diets$eaten.C, drawContour=FALSE)
} else { 
	fit <- predictPlotGam(gam, plane, graphLoc, 'Feed protein as proportion of all carbon (%)')
}
writeGamSummary(gam, paste(out_dir,'__gamSummary-relCcas.txt',sep=''))

gam <- generateGam('cecum_relNmuc', diets)
graphLoc <- paste(out_dir, '__cecumRelAb-Nmuc.pdf', sep='')
if (drawDots)
{
	fit <- predictPlotGam(gam, plane, graphLoc, 'Mucin as proportion of all nitrogen (%)', xdots=diets$eaten.P, ydots=diets$eaten.C, drawContour=FALSE)
} else { 
	fit <- predictPlotGam(gam, plane, graphLoc, 'Mucin as proportion of all nitrogen (%)')
}
writeGamSummary(gam, paste(out_dir,'__gamSummary-relNmuc.txt',sep=''))

gam <- generateGam('cecum_relNcas', diets)
graphLoc <- paste(out_dir, '__cecumRelAb-Ncas.pdf', sep='')
if (drawDots)
{
	fit <- predictPlotGam(gam, plane, graphLoc, 'Feed protein as proportion of all nitrogen (%)', xdots=diets$eaten.P, ydots=diets$eaten.C, drawContour=FALSE)
} else { 
	fit <- predictPlotGam(gam, plane, graphLoc, 'Feed protein as proportion of all nitrogen (%)')
}
writeGamSummary(gam, paste(out_dir,'__gamSummary-relNcas.txt',sep=''))


# ###################################################
# Plot relative abundances of individual carbon/nitrogen sources as proportion of all INTAKE carbon/nitrogen.
gam <- generateGamNoFat('relCw', diets)
graphLoc <- paste(out_dir, '__relAb-Cw.pdf', sep='')
if (drawDots)
{
	fit <- predictPlotGam(gam, plane, graphLoc, 'Wheatstarch C as proportion of all carbon (%)', xdots=diets$eaten.P, ydots=diets$eaten.C, drawContour=FALSE)
} else { 
	fit <- predictPlotGam(gam, plane, graphLoc, 'Wheatstarch C as proportion of all carbon (%)')
}
writeGamSummary(gam, paste(out_dir,'__gamSummary-relCw.txt',sep=''))

gam <- generateGam('relCd', diets)
graphLoc <- paste(out_dir, '__relAb-Cd.pdf', sep='')
if (drawDots)
{
	fit <- predictPlotGam(gam, plane, graphLoc, 'Dex WS C as proportion of all carbon(%)', xdots=diets$eaten.P, ydots=diets$eaten.C, drawContour=FALSE)
} else { 
	fit <- predictPlotGam(gam, plane, graphLoc, 'Dex WS C as proportion of all carbon(%)')
}
writeGamSummary(gam, paste(out_dir,'__gamSummary-relCd.txt',sep=''))

gam <- generateGam('relCm', diets)
graphLoc <- paste(out_dir, '__relAb-Cm.pdf', sep='')
if (drawDots)
{
	fit <- predictPlotGam(gam, plane, graphLoc, 'Dex WS C as proportion of all carbon(%)', xdots=diets$eaten.P, ydots=diets$eaten.C, drawContour=FALSE)
} else { 
	fit <- predictPlotGam(gam, plane, graphLoc, 'Mucin C as proportion of all carbon (%)')
}
writeGamSummary(gam, paste(out_dir,'__gamSummary-relCm.txt',sep=''))

gam <- generateGam('relNf', diets)
graphLoc <- paste(out_dir, '__relAb-Nf.pdf', sep='')
if (drawDots)
{
	fit <- predictPlotGam(gam, plane, graphLoc, 'Dex WS C as proportion of all carbon(%)', xdots=diets$eaten.P, ydots=diets$eaten.C, drawContour=FALSE)
} else { 
	fit <- predictPlotGam(gam, plane, graphLoc, 'Feed N as proportion of all nitrogen (%)')
}
writeGamSummary(gam, paste(out_dir,'__gamSummary-relNf.txt',sep=''))

gam <- generateGam('relNm', diets)
graphLoc <- paste(out_dir, '__relAb-Nm.pdf', sep='')
if (drawDots)
{
	fit <- predictPlotGam(gam, plane, graphLoc, 'Mucin N as proportion of all nitrogen (%)', xdots=diets$eaten.P, ydots=diets$eaten.C, drawContour=FALSE)
} else { 
	fit <- predictPlotGam(gam, plane, graphLoc, 'Mucin N as proportion of all nitrogen (%)')
}
writeGamSummary(gam, paste(out_dir,'__gamSummary-relNm.txt',sep=''))


# ################################################### 
# Plot daily intakes.
gam <- generateGam('inCw', diets)
graphLoc <- paste(out_dir, '__intake-Cw.pdf', sep='')
if (drawDots)
{
	fit <- predictPlotGam(gam, plane, graphLoc, 'Wheat starch daily intake (g)', xdots=diets$eaten.P, ydots=diets$eaten.C, drawContour=FALSE)
} else { 
	fit <- predictPlotGam(gam, plane, graphLoc, 'Wheat starch daily intake (g)')
}
writeGamSummary(gam, paste(out_dir,'__gamSummary-inCw.txt',sep=''))

gam <- generateGam('inCd', diets)
graphLoc <- paste(out_dir, '__intake-Cd.pdf', sep='')
if (drawDots)
{
	fit <- predictPlotGam(gam, plane, graphLoc, 'Dextrinised wheat starch daily intake (g)', xdots=diets$eaten.P, ydots=diets$eaten.C, drawContour=FALSE)
} else { 
	fit <- predictPlotGam(gam, plane, graphLoc, 'Dextrinised wheat starch daily intake (g)')
}	
writeGamSummary(gam, paste(out_dir,'__gamSummary-inCd.txt',sep=''))

gam <- generateGam('inNf', diets)
graphLoc <- paste(out_dir, '__intake-Nf.pdf', sep='')
if (drawDots)
{
	fit <- predictPlotGam(gam, plane, graphLoc, 'Protein daily intake (g)', xdots=diets$eaten.P, ydots=diets$eaten.C, drawContour=FALSE)
} else { 
	fit <- predictPlotGam(gam, plane, graphLoc, 'Protein daily intake (g)')
}
writeGamSummary(gam, paste(out_dir,'__gamSummary-inNf.txt',sep=''))


###################################################
# Plot small intestine nutrient absorption in absolute quantities. 
gam <- generateGam('siCw', diets)
graphLoc <- paste(out_dir, '__si_absorb-Cw.pdf', sep='')
if (drawDots)
{
	fit <- predictPlotGam(gam, plane, graphLoc, 'SI wheat starch absorption (g)', xdots=diets$eaten.P, ydots=diets$eaten.C, drawContour=FALSE)
} else { 
	fit <- predictPlotGam(gam, plane, graphLoc, 'SI wheat starch absorption (g)')
}
writeGamSummary(gam, paste(out_dir,'__gamSummary-siCw.txt',sep=''))

gam <- generateGam('siCd', diets)
graphLoc <- paste(out_dir, '__si_absorb-Cd.pdf', sep='')
if (drawDots)
{
	fit <- predictPlotGam(gam, plane, graphLoc, 'SI DWS wheat starch absorption (g)', xdots=diets$eaten.P, ydots=diets$eaten.C, drawContour=FALSE)
} else { 
	fit <- predictPlotGam(gam, plane, graphLoc, 'SI DWS wheat starch absorption (g)')
}
writeGamSummary(gam, paste(out_dir,'__gamSummary-siCd.txt',sep=''))

gam <- generateGam('siNf', diets)
graphLoc <- paste(out_dir, '__si_absorb-Nf.pdf', sep='')
if (drawDots)
{
	fit <- predictPlotGam(gam, plane, graphLoc, 'SI protein absorption (g)', xdots=diets$eaten.P, ydots=diets$eaten.C, drawContour=FALSE)
} else { 
	fit <- predictPlotGam(gam, plane, graphLoc, 'SI protein absorption (g)')
}
writeGamSummary(gam, paste(out_dir,'__gamSummary-ciNf.txt',sep=''))


###################################################
# Plot small intestine nutrient absorption as a proportion of what is consumed. 
gam <- generateGam('siPCw', diets)
graphLoc <- paste(out_dir, '__si_absorb_prop-Cw.pdf', sep='')
if (drawDots)
{
	fit <- predictPlotGam(gam, plane, graphLoc, 'SI wheat starch absorption proportion', xdots=diets$eaten.P, ydots=diets$eaten.C, drawContour=FALSE)
} else { 
	fit <- predictPlotGam(gam, plane, graphLoc, 'SI wheat starch absorption proportion')
}
writeGamSummary(gam, paste(out_dir,'__gamSummary-siP-Cw.txt',sep=''))

gam <- generateGam('siPCd', diets)
graphLoc <- paste(out_dir, '__si_absorb_prop-Cd.pdf', sep='')
if (drawDots)
{
	fit <- predictPlotGam(gam, plane, graphLoc, 'SI DWS absorption proportion', xdots=diets$eaten.P, ydots=diets$eaten.C, drawContour=FALSE)
} else { 
	fit <- predictPlotGam(gam, plane, graphLoc, 'SI DWS absorption proportion')
}
writeGamSummary(gam, paste(out_dir,'__gamSummary-siP-Cd.txt',sep=''))

gam <- generateGam('siPNf', diets)
graphLoc <- paste(out_dir, '__si_absorb_prop-Nf.pdf', sep='')
if (drawDots)
{
	fit <- predictPlotGam(gam, plane, graphLoc, 'SI protein absorption proportion', xdots=diets$eaten.P, ydots=diets$eaten.C, drawContour=FALSE)
} else {
	fit <- predictPlotGam(gam, plane, graphLoc, 'SI protein absorption proportion')
}
writeGamSummary(gam, paste(out_dir,'__gamSummary-siP-Nf.txt',sep=''))


###################################################
# Plot daily nutrients available in the colon. 
gam <- generateGam('colicCw', diets)
graphLoc <- paste(out_dir, '__colic-Cw.pdf', sep='')
if (drawDots)
{
	fit <- predictPlotGam(gam, plane, graphLoc, 'Colic wheat starch availability (g)', xdots=diets$eaten.P, ydots=diets$eaten.C, drawContour=FALSE)
} else { 
	fit <- predictPlotGam(gam, plane, graphLoc, 'Colic wheat starch availability (g)')
}
writeGamSummary(gam, paste(out_dir,'__gamSummary-colicCw.txt',sep=''))

gam <- generateGam('colicCd', diets)
graphLoc <- paste(out_dir, '__colic-Cd.pdf', sep='')
if (drawDots)
{
	fit <- predictPlotGam(gam, plane, graphLoc, 'Colic DWS availability (g)', xdots=diets$eaten.P, ydots=diets$eaten.C, drawContour=FALSE)
} else { 
	fit <- predictPlotGam(gam, plane, graphLoc, 'Colic DWS availability (g)')
}
writeGamSummary(gam, paste(out_dir,'__gamSummary-colicCd.txt',sep=''))

gam <- generateGam('colicNf', diets)
graphLoc <- paste(out_dir, '__colic-Nf.pdf', sep='')
if (drawDots)
{
	fit <- predictPlotGam(gam, plane, graphLoc, 'Colic chow-protein availability (g)', xdots=diets$eaten.P, ydots=diets$eaten.C, drawContour=FALSE)
} else { 
	fit <- predictPlotGam(gam, plane, graphLoc, 'Colic chow-protein availability (g)')
}
writeGamSummary(gam, paste(out_dir,'__gamSummary-colicNf.txt',sep=''))


