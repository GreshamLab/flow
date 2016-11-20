#######################
#gating.R       
#
#started: 10/23/2016
#
#author1: D Gresham
######################

###########################################################################################################################
#This script is intended to read in .fcs files and perform manual gating for i) single cells, ii) debris and iii) fluorescence 
#
#Gating is performed with untransformed data
#
#Individual gates are saved in a file gates.Rdata for use with the Gresham Lab Flow Cytometry Analysis.Rmd pipeline
###########################################################################################################################

##To be run the first time if packages are not installed.
#source("http://bioconductor.org/biocLite.R")
#biocLite("flowViz")
#biocLite("flowCore")

#Load libraries
library(flowCore)
library(flowViz)

#Read in the data
flowData <- read.flowSet(path = ".", pattern=".fcs", alter.names=TRUE)

#Check how many cells were counted in each fcs file
fsApply(flowData, each_col, length)

##############################
#1. Generate gate for singlet cells
#this gate is defined on the basis of the relationship between forward scatter height and area
plot(flowData[[1]], c('FSC.H','FSC.A'), 
  xlim=c(0,3e6), ylim=c(0,3e6),
  main="keep clicking to surround the singles, then right click when you're done",
  smooth=F)
Agate <- locator(type="l",col="red") 
  # alternative approach could be devised with identify()
gm.1 <- as.matrix(data.frame(FSC.H=Agate$x,FSC.A=Agate$y))
pg.singlets <- polygonGate(filterId="singlets",.gate=gm.1)

#test that the singlet gate looks reasonable for the sample
xyplot(FSC.A~FSC.H,data=flowData[[1]],xlim=c(0,3e6), ylim=c(0,3e6), smooth=F, filter=pg.singlets, outline=T)

#test that the gate looks reasonable over all the samples
xyplot(FSC.A~FSC.H, data=flowData, xlim=c(0,3e6), ylim=c(0,3e6), 
       smooth=F, filter=pg.singlets, outline=T, displayFilter=TRUE,
       stat=T, pos=0.5, abs=T)

##############################
#2. Generate Gate for debris based on forward scatter and side scatter
plot(flowData[[1]], c('FSC.A','SSC.A'), xlim=c(0,3e6), ylim=c(0,3e5), smooth=F)
Bgate <- locator(10)
gm.2 <- matrix(,10,2)
colnames(gm.2) <- c('FSC.A','SSC.A')
gm.2[,1] <- Bgate$x
gm.2[,2] <- Bgate$y
pg.nondebris <- polygonGate(filterId="nonDebris",.gate=gm.2)

#test that the debris gate looks reasonable for the sample
xyplot(SSC.A ~ FSC.A, data=flowData[[1]], displayFilter=TRUE, xlim=c(0,3e6), ylim=c(0,3e6), filter=pg.nondebris, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)

#test that the decris gate looks reasonable over all the samples
xyplot(SSC.A ~ FSC.A, data=flowData, displayFilter=TRUE, xlim=c(0,3e5), ylim=c(0,3e6), filter=pg.nondebris, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)

##############################
#3. Generate Gate for non-fluorescencing cells
plot(flowData[[1]], c('FSC.A','FL1.A'), xlim=c(0,5e6), ylim=c(0,5e4), smooth=F)
Cgate <- locator(10)
gm.3 <- matrix(,10,2)
colnames(gm.3) <- c('FSC.A','FL1.A')
gm.3[,1] <- Cgate$x
gm.3[,2] <- Cgate$y
pg.nongfp <- polygonGate(filterId="GFPneg",.gate=gm.3)

#test that the non-fluorescing gate looks reasonable for the sample
xyplot(FL1.A~FSC.A,data=flowData[[1]], xlim=c(0,5e6), ylim=c(0,5e4), filter=pg.nongfp, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)

#test that the non-fluorescing gate looks reasonable over all the samples
xyplot(FL1.A~FSC.A,data=flowData, xlim=c(0,5e6), ylim=c(0,5e4), filter=pg.nongfp, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)

##############################
#4. Generate Gate for fluorescencing cells
plot(flowData[[2]], c('FSC.A','FL1.A'), xlim=c(0,5e6), ylim=c(0,5e4), smooth=F)
Dgate <- locator(10)
gm.4 <- matrix(,10,2)
colnames(gm.4) <- c('FSC.A','FL1.A')
gm.4[,1] <- Dgate$x
gm.4[,2] <- Dgate$y
pg.gfp <- polygonGate(filterId="GFPpos",.gate=gm.4)

#test that the fluorescing gate looks reasonable for the sample
xyplot(FL1.A~FSC.A,data=flowData[[1]], xlim=c(0,5e6), ylim=c(0,5e4), filter=pg.gfp, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)

#test that the fluorescing gate looks reasonable over all the samples
xyplot(FL1.A~FSC.A,data=flowData, xlim=c(0,5e6), ylim=c(0,5e4), filter=pg.gfp, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)

################################

#5. Generate Gate for high fluorescencing cells
plot(flowData[[6]], c('FSC.A','FL1.A'), xlim=c(0,5e6), ylim=c(0,5e4), smooth=F)
Egate <- locator(10)
gm.5 <- matrix(,10,2)
colnames(gm.5) <- c('FSC.A','FL1.A')
gm.5[,1] <- Egate$x
gm.5[,2] <- Egate$y
pg.hi.gfp <- polygonGate(filterId="hiGFPpos",.gate=gm.5)

#test that the high fluorescing gate looks reasonable for the sample
xyplot(FL1.A~FSC.A,data=flowData[[6]], xlim=c(0,5e6), ylim=c(0,5e4), filter=pg.hi.gfp, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)

#test that the high fluorescing gate looks reasonable over all the samples
xyplot(FL1.A~FSC.A,data=flowData, xlim=c(0,5e6), ylim=c(0,5e4), filter=pg.hi.gfp, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)


################################
##In some cases a rectangle gate may be desired, which can be generated as follows
#rectGate <- rectangleGate("FSC.A"= c(0,3.5e6),"SSC.A" = c(0,4.5e6))

#xyplot(SSC.A~FSC.A,data=flowData, xlim=c(0,5e6), ylim=c(0,5e4), filter=rectGate, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)

################################

#6. Save the gate information to an R data file

rm(list=c("flowData"))
save.image(file="gates.Rdata")


