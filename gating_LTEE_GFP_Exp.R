#######################
#gating.R       
#
#started: 01/07/2016
#
#author1: G Avecilla, S Lauer, D Gresham
######################

######################
#This script is specific for analyzing the data obtained in G Avecilla's experiment: EE_GA_glnurmix_barcode17.
#For any other purpose, the script must be modified accordingly. 

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

##Confirm the number of .fcs files in your folder. The script below is only accurate if there are 32 .fcs files
str(flowData) #it should say there are 32 observations. if not, rework the script

##############################
#1. Generate gate for singlet cells####
#this gate is defined on the basis of the relationship between forward scatter height and area
plot(flowData[[1]], c('FSC.H','FSC.A'), xlim=c(0,3e6), xaxs = "i", yaxs = "i", ylim=c(0,3e6),smooth=F)
Agate <- locator(10)
gm.1 <- matrix(,8,2)
colnames(gm.1) <- c('FSC.H','FSC.A')
gm.1[,1] <- Agate$x
gm.1[,2] <- Agate$y
pg.singlets <- polygonGate(filterId="singlets",.gate=gm.1)

polygon(Agate$x, Agate$y, border='red')

#test that the singlet gate looks reasonable for some samples (these are all 0 copy controls)
xyplot(FSC.A~FSC.H,data=flowData[[27]],xlim=c(0,3e6), ylim=c(0,3e6), smooth=F, filter=pg.singlets, outline=T)
xyplot(FSC.A~FSC.H,data=flowData[[22]],xlim=c(0,3e6), ylim=c(0,3e6), smooth=F, filter=pg.singlets, outline=T)

#a one copy control
xyplot(FSC.A~FSC.H,data=flowData[[9]],xlim=c(0,3e6), ylim=c(0,3e6), smooth=F, filter=pg.singlets, outline=T)

#a two copy control
xyplot(FSC.A~FSC.H,data=flowData[[12]],xlim=c(0,3e6), ylim=c(0,3e6), smooth=F, filter=pg.singlets, outline=T)


##############################
#2. Generate Gate for debris based on forward scatter and side scatter. ####
#This needs to be done separately for each media condition.

######GLUTAMINE#####
plot(flowData[[1]], c('FSC.A','SSC.A'), xlim=c(0,2e6), xaxs = "i", yaxs = "i", ylim=c(0,6e5), smooth=F)
Bgate <- locator(10)
gm.2 <- matrix(,8,2)
colnames(gm.2) <- c('FSC.A','SSC.A')
gm.2[,1] <- Bgate$x
gm.2[,2] <- Bgate$y
pg.nondebris.gln <- polygonGate(filterId="nonDebrisGln",.gate=gm.2)

#test that the debris gate looks reasonable for some samples (these are 1 and 2 copy controls)
xyplot(SSC.A ~ FSC.A, data=flowData[[9]], displayFilter=TRUE, xlim=c(0,2e6), ylim=c(0,1e6), filter=pg.nondebris.gln, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)
xyplot(SSC.A ~ FSC.A, data=flowData[[17]], displayFilter=TRUE, xlim=c(0,2e6), ylim=c(0,1e6), filter=pg.nondebris.gln, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)


#######UREA######
plot(flowData[[27]], c('FSC.A','SSC.A'), xlim=c(0,2e6), xaxs = "i", yaxs = "i", ylim=c(0,6e5), smooth=F)
Cgate <- locator(10)
gm.3 <- matrix(,8,2)
colnames(gm.3) <- c('FSC.A','SSC.A')
gm.3[,1] <- Cgate$x
gm.3[,2] <- Cgate$y
pg.nondebris.ur <- polygonGate(filterId="nonDebrisUr",.gate=gm.3)

#test that the debris gate looks reasonable for some samples (these are 1 and2 copy controls)
xyplot(SSC.A ~ FSC.A, data=flowData[[4]], displayFilter=TRUE, xlim=c(0,3e6), ylim=c(0,1e6), filter=pg.nondebris.ur, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)
xyplot(SSC.A ~ FSC.A, data=flowData[[12]], displayFilter=TRUE, xlim=c(0,3e6), ylim=c(0,1e6), filter=pg.nondebris.ur, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)

######GLUTAMINE and UREA mix#####
plot(flowData[[22]], c('FSC.A','SSC.A'), xlim=c(0,2e6), xaxs = "i", yaxs = "i", ylim=c(0,6e5), smooth=F)
Dgate <- locator(10)
gm.4 <- matrix(,8,2)
colnames(gm.4) <- c('FSC.A','SSC.A')
gm.4[,1] <- Dgate$x
gm.4[,2] <- Dgate$y
pg.nondebris.mix <- polygonGate(filterId="nonDebrisMix",.gate=gm.4)

#test that the debris gate looks reasonable for some samples (these are 1 and 2 copy controls)
xyplot(SSC.A ~ FSC.A, data=flowData[[30]], displayFilter=TRUE, xlim=c(0,2e6), ylim=c(0,1e6), filter=pg.nondebris.mix, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)
xyplot(SSC.A ~ FSC.A, data=flowData[[7]], displayFilter=TRUE, xlim=c(0,2e6), ylim=c(0,1e6), filter=pg.nondebris.mix, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)

#############################
#3. FLUORESCENCE####

####Generate gates for 0, 1, 2, and 3+ copies for GLUTAMINE-limitation####

##Plot the control sample that has non-fluorescing cells (0 copy)
plot(flowData[[1]], c('FSC.A','FL1.A'), xlim=c(-1e4,3e6), xaxs = "i", yaxs = "i", ylim=c(-1e4,1e5), smooth=F)
glngate.0 <- locator(10, type='l', col='red')

##Plot the control sample that has 1 copy
plot(flowData[[9]], c('FSC.A','FL1.A'), xlim=c(0,2e6), xaxs = "i", yaxs = "i", ylim=c(0,4e5), smooth=F)

##Overlay your previous gate
polygon(glngate.0$x, glngate.0$y, border='red')

##Draw a new gate for the one copy
glngate.1 <- locator(10, type='l', col='blue')

##Overlay and check the new gate
plot(flowData[[9]], c('FSC.A','FL1.A'), xlim=c(0,2e6), xaxs = "i", yaxs = "i", ylim=c(0,3e5), smooth=F)
polygon(glngate.0$x, glngate.0$y, border='red')
polygon(glngate.1$x, glngate.1$y, border='blue')

##Plot the control sample that has 2 copies
plot(flowData[[17]], c('FSC.A','FL1.A'), xlim=c(0,2e6), xaxs = "i", yaxs = "i", ylim=c(0,4e5), smooth=F)
polygon(glngate.0$x, glngate.0$y, border='red')
polygon(glngate.1$x, glngate.1$y, border='blue')
glngate.2 <- locator(10, type='l', col='green')
glngate.3 <- locator(10, type="l", col='purple')

##Check how the gates look on the sample
plot(flowData[[17]], c('FSC.A','FL1.A'), xlim=c(0,3e6), xaxs = "i", yaxs = "i", ylim=c(0,5e5), smooth=F)
polygon(glngate.0$x, glngate.0$y, border='red')
polygon(glngate.1$x, glngate.1$y, border='blue')
polygon(glngate.2$x, glngate.2$y, border='green')
polygon(glngate.3$x, glngate.3$y, border='purple')

glngm.0 <- matrix(,7,2)
glngm.1 <- matrix(,7,2)
glngm.2 <- matrix(,10,2)
glngm.3 <- matrix(,5,2)

colnames(glngm.0) <- c('FSC.A','FL1.A')
colnames(glngm.1) <- c('FSC.A','FL1.A')
colnames(glngm.2) <- c('FSC.A','FL1.A')
colnames(glngm.3) <- c('FSC.A','FL1.A')

glngm.0[,1] <- glngate.0$x
glngm.1[,1] <- glngate.1$x
glngm.2[,1] <- glngate.2$x
glngm.3[,1] <- glngate.3$x

glngm.0[,2] <- glngate.0$y
glngm.1[,2] <- glngate.1$y
glngm.2[,2] <- glngate.2$y
glngm.3[,2] <- glngate.3$y

gln.zero<- polygonGate(filterId="ZeroCopyGlutamine",.gate=glngm.0)
gln.one<- polygonGate(filterId="OneCopyGlutamine",.gate=glngm.1)
gln.two<- polygonGate(filterId="TwoCopyGlutamine",.gate=glngm.2)
gln.three<- polygonGate(filterId="ThreeCopyGlutamine",.gate=glngm.3)


#################################################################
####Generate gates for 0, 1, 2, and 3+ copies for UREA-limitation####

##Plot the control sample that has non-fluorescing cells (0 copy)
plot(flowData[[27]], c('FSC.A','FL1.A'), xlim=c(-1e4,3e6), xaxs = "i", yaxs = "i", ylim=c(-1e4,1e5), smooth=F)
urgate.0 <- locator(10, type='l', col='red')

##Plot the control sample that has 1 copy
plot(flowData[[4]], c('FSC.A','FL1.A'), xlim=c(0,3e6), xaxs = "i", yaxs = "i", ylim=c(0,5e5), smooth=F)

##Overlay your previous gate
polygon(urgate.0$x, urgate.0$y, border='red')

##Draw a new gate for the one copy
urgate.1 <- locator(10, type='l', col='blue')

##Overlay and check the new gate
plot(flowData[[4]], c('FSC.A','FL1.A'), xlim=c(0,3e6), xaxs = "i", yaxs = "i", ylim=c(0,5e5), smooth=F)
polygon(urgate.0$x, urgate.0$y, border='red')
polygon(urgate.1$x, urgate.1$y, border='blue')

##Plot the control sample that has 2 copies
plot(flowData[[12]], c('FSC.A','FL1.A'), xlim=c(0,2e6), xaxs = "i", yaxs = "i", ylim=c(0,7e5), smooth=F)
polygon(urgate.0$x, urgate.0$y, border='red')
polygon(urgate.1$x, urgate.1$y, border='blue')
urgate.2 <- locator(10, type='l', col='green')
urgate.3 <- locator(10, type="l", col='purple')

##Check how the gates look on the sample
plot(flowData[[12]], c('FSC.A','FL1.A'), xlim=c(0,3e6), xaxs = "i", yaxs = "i", ylim=c(0,5e5), smooth=F)
polygon(urgate.0$x, urgate.0$y, border='red')
polygon(urgate.1$x, urgate.1$y, border='blue')
polygon(urgate.2$x, urgate.2$y, border='green')
polygon(urgate.3$x, urgate.3$y, border='purple')

urgm.0 <- matrix(,4,2)
urgm.1 <- matrix(,8,2)
urgm.2 <- matrix(,9,2)
urgm.3 <- matrix(,6,2)

colnames(urgm.0) <- c('FSC.A','FL1.A')
colnames(urgm.1) <- c('FSC.A','FL1.A')
colnames(urgm.2) <- c('FSC.A','FL1.A')
colnames(urgm.3) <- c('FSC.A','FL1.A')

urgm.0[,1] <- urgate.0$x
urgm.1[,1] <- urgate.1$x
urgm.2[,1] <- urgate.2$x
urgm.3[,1] <- urgate.3$x

urgm.0[,2] <- urgate.0$y
urgm.1[,2] <- urgate.1$y
urgm.2[,2] <- urgate.2$y
urgm.3[,2] <- urgate.3$y

ur.zero<- polygonGate(filterId="ZeroCopyUrea",.gate=urgm.0)
ur.one<- polygonGate(filterId="OneCopyUrea",.gate=urgm.1)
ur.two<- polygonGate(filterId="TwoCopyUrea",.gate=urgm.2)
ur.three<- polygonGate(filterId="ThreeCopyUrea",.gate=urgm.3)

#################################################################
####Generate gates for 0, 1, 2, and 3+ copies for GLUTAMINE and UREA MIX####

##Plot the control sample that has non-fluorescing cells (0 copy)
plot(flowData[[22]], c('FSC.A','FL1.A'), xlim=c(-1e4,3e6), xaxs = "i", yaxs = "i", ylim=c(-1e4,1e5), smooth=F)
mixgate.0 <- locator(10, type='l', col='red')

##Plot the control sample that has 1 copy
plot(flowData[[30]], c('FSC.A','FL1.A'), xlim=c(0,2e6), xaxs = "i", yaxs = "i", ylim=c(0,4e5), smooth=F)

##Overlay your previous gate
polygon(mixgate.0$x, mixgate.0$y, border='red')

##Draw a new gate for the one copy
mixgate.1 <- locator(10, type='l', col='blue')

##Overlay and check the new gate
plot(flowData[[30]], c('FSC.A','FL1.A'), xlim=c(0,3e6), xaxs = "i", yaxs = "i", ylim=c(0,5e5), smooth=F)
polygon(mixgate.0$x, mixgate.0$y, border='red')
polygon(mixgate.1$x, mixgate.1$y, border='blue')

##Plot the control sample that has 2 copies
plot(flowData[[7]], c('FSC.A','FL1.A'), xlim=c(0,2e6), xaxs = "i", yaxs = "i", ylim=c(0,5e5), smooth=F)
polygon(mixgate.0$x, mixgate.0$y, border='red')
polygon(mixgate.1$x, mixgate.1$y, border='blue')
mixgate.2 <- locator(10, type='l', col='green')
mixgate.3 <- locator(10, type="l", col='purple')

##Check how the gates look on the sample
plot(flowData[[7]], c('FSC.A','FL1.A'), xlim=c(0,3e6), xaxs = "i", yaxs = "i", ylim=c(0,5e5), smooth=F)
polygon(mixgate.0$x, mixgate.0$y, border='red')
polygon(mixgate.1$x, mixgate.1$y, border='blue')
polygon(mixgate.2$x, mixgate.2$y, border='green')
polygon(mixgate.3$x, mixgate.3$y, border='purple')

mixgm.0 <- matrix(,4,2)
mixgm.1 <- matrix(,9,2)
mixgm.2 <- matrix(,9,2)
mixgm.3 <- matrix(,6,2)

colnames(mixgm.0) <- c('FSC.A','FL1.A')
colnames(mixgm.1) <- c('FSC.A','FL1.A')
colnames(mixgm.2) <- c('FSC.A','FL1.A')
colnames(mixgm.3) <- c('FSC.A','FL1.A')

mixgm.0[,1] <- mixgate.0$x
mixgm.1[,1] <- mixgate.1$x
mixgm.2[,1] <- mixgate.2$x
mixgm.3[,1] <- mixgate.3$x

mixgm.0[,2] <- mixgate.0$y
mixgm.1[,2] <- mixgate.1$y
mixgm.2[,2] <- mixgate.2$y
mixgm.3[,2] <- mixgate.3$y

mix.zero<- polygonGate(filterId="ZeroCopyMix",.gate=mixgm.0)
mix.one<- polygonGate(filterId="OneCopyMix",.gate=mixgm.1)
mix.two<- polygonGate(filterId="TwoCopyMix",.gate=mixgm.2)
mix.three<- polygonGate(filterId="ThreeCopyMix",.gate=mixgm.3)

#Save the gate information to an R data file
rm(list=c("flowData")) 
save.image(file="gates_sample4_word.Rdata")
