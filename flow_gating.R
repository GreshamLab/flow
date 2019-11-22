#######################
#gating.R       
#
#started: 01/07/2016
#modified: 2/14/2019
#
#author1: G Avecilla, S Lauer, D Gresham
#author2: N Brandt
######################

######################
#This script is a generic file specific for analyzing the data from flow cytometry in the FCS3 formate
#For your exact purposes, the script must be modified accordingly. 

###########################################################################################################################
#This script is intended to read in .fcs files and perform manual gating for i) single cells, ii) debris and iii) fluorescence 
#The fcs files should be in a folder with a name descriptive of the experiment and have a assocated sample sheet with the name format of samplesheet_<name>
#Should include at least the following:
#       * column1 = Well
#       * column2 = Strain

#
#Gating is performed with untransformed data
#
#Individual gates are saved in a file gates.Rdata for use with the Gresham Lab Flow Cytometry Analysis.Rmd pipeline
###########################################################################################################################

##To be run the first time if packages are not installed.
#install.packages("BiocManager")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("flowViz")
#BiocManager::install("flowCore")
#BiocManager::install("ggcyto")
#install.packages('ggforce')


#Load libraries
library(flowCore)
library(flowViz)
library(ggcyto)
library(ggforce)


#Read in the data
#Set working directory to the folder in which you have stored your .fcs files
#Read in all the fcs files in the directory.

#working directory
dir = ''

#file location
#Used for data files being outside of home directory
path.data =   ''

#set name of run to create gates for
#Should be the same name as the folder that your .fcs files are in
name <- ""

#load sample sheet
sample.sheet <- read.csv(paste(path.data,"/samplesheet_",name,".csv", sep=""))

#read in fcs files in order presented in sample sheet (based on well identifier)
files <- paste(path.data,"/",sort(factor(list.files(paste(path.data,"/", sep=""),full.names=FALSE), levels = paste(sample.sheet$Well,".fcs",sep="" ), ordered=TRUE)),sep="")
flowData <- read.ncdfFlowSet(files=files, pattern=".fcs", alter.names = TRUE)

#Ensures that the number of entries in the sample sheet match the number of flowFrames in the flowSet
sample.ind <- which(paste(sample.sheet$Well,".fcs", sep="") %in% sampleNames(flowData))
sample.sheet <- sample.sheet[sample.ind,]
sample.sheet <- sample.sheet[order(sample.sheet$Well),]

#rename sample name of flow set to make it easier to identify
sampleNames(flowData) <- paste(gsub(" ","_",sample.sheet$Strain),"_",sub(" ","_",sample.sheet$Well), sep="")

#set controls on which to gate on
#these are some of the base types of controls you may want for each channel you are taking measurements from 
neg.signal <- "THIS IS A BETTER CHANNEL"
pos.signal <- "YES"
hi.signal <- "IT IS"


##############################
#1. Generate gate for singlet cells####
#this gate is defined on the basis of the relationship between forward scatter height and area
#******Please note you may need to adjust the x and y plot limits to properly visualize your data

plot(flowData[[neg.signal]], c('FSC.H','FSC.A'), xlim=c(0,3e6), ylim=c(0,3e6),smooth=T)

singlet.gate <- locator(100, type='l', col='red')
gm.1 <- matrix(,length(singlet.gate$x),2)
colnames(gm.1) <- c('FSC.H','FSC.A')
gm.1[,1] <- singlet.gate$x
gm.1[,2] <- singlet.gate$y
pg.singlets <- polygonGate(filterId="singlets",.gate=gm.1)

ggcyto(flowData[neg.signal], aes(x = `FSC.H`, y =  `FSC.A`)) + geom_hex(bins = 512) + geom_gate(pg.singlets)

#Look at the gating on the controls
ggcyto(flowData[c(neg.signal,pos.signal,hi.signal)], aes(x = `FSC.H`, y =  `FSC.A`)) + geom_hex(bins = 512) + ggcyto_par_set(limits = list(x = c(0,1e5), y = c(0,3e4))) + geom_gate(pg.singlets) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)

#test that the singlet gate looks reasonable for All samples
for(i in 1:round(length(flowData)/4,0)){
  plot <- ggcyto(flowData, aes(x = `FSC.H`, y =  `FSC.A`)) + geom_hex(bins = 512) + ggcyto_par_set(limits = list(x = c(0,1e5), y = c(0,3e4))) + geom_gate(pg.singlets) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = i)
  print(plot)
}

#Filter out the doublets
flowData.singlets <- Subset(flowData,pg.singlets)


##############################
#2. Generate Gate for debris based on forward scatter and side scatter. ####
#This needs to be done separately for each media condition.
#******Please note you may need to adjust the x and y plot limits to properly visualize your data



plot(flowData.singlets[[neg.signal]], c('FSC.A','SSC.A'), xlim=c(0,3e6), ylim=c(0,1e6),smooth=T)

debris.gate <- locator(100, type='l', col='red')
gm.2 <- matrix(,length(debris.gate$x),2)
colnames(gm.2) <- c('FSC.A','SSC.A')
gm.2[,1] <- debris.gate$x
gm.2[,2] <- debris.gate$y
pg.nondebris <- polygonGate(filterId="nonDebris",.gate=gm.2)

#Look at the gating on the controls
ggcyto(flowData.singlets[c(neg.signal,pos.signal,hi.signal)], aes(x = `FSC.A`, y =  `SSC.A`)) + geom_hex(bins = 512) + ggcyto_par_set(limits = list(x = c(0,1e5), y = c(0,3e4))) + geom_gate(pg.nondebris) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)


#test that the singlet gate looks reasonable for All samples
for(i in 1:round(length(flowData)/4,0)){
  plot <- ggcyto(flowData.singlets, aes(x = `FSC.A`, y =  `SSC.A`)) + geom_hex(bins = 512) + ggcyto_par_set(limits = list(x = c(0,1e5), y = c(0,3e4))) + geom_gate(pg.nondebris) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = i)
  print(plot)
}

#Filter out the debris
flowData.nondebris <- Subset(flowData.singlets,pg.nondebris)


#############################
#3. FLUORESCENCE####
#******Please note you may need to adjust the x and y plot limits to properly visualize your data
# You may need to add or modify these gates depending on how you are analyzing the signal and how many channels you are working with
#This is the base code block you can use to add more gates CONTROLVARIABLE, TYPEOFSIGNAL, GATENUMBERINCODE, TYPEOFSIGNAL, should all be replaced with vairables that make sense in your code
######START
##Plot 
#     plot(flowData[[CONTROLVARIABLE]], c('FSC.A','FL1.A'), xlim=c(0,3e6), ylim=c(0,5e4),smooth=T)
##draw the gate and create a gate ploygon
#      TYPEOFSIGNAL.gate <- locator(100, type='l', col='red')
#      gm.GATENUMBERINCODE <- matrix(,length(TYPEOFSIGNAL.gate$x),2)
#      colnames(gm.GATENUMBERINCODE) <- c('FSC.A','FL1.A')
#      gm.GATENUMBERINCODE[,1] <- TYPEOFSIGNAL.gate$x
#      gm.GATENUMBERINCODE[,2] <- TYPEOFSIGNAL.gate$y
#      gate.TYPEOFSIGNAL <- polygonGate(filterId="TYPEOFSIGNAL",.gate=gm.GATENUMBERINCODE)
##Overlay and check the new gate
#     ggcyto(flowData[c(CONTROLVARIABLE)], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e4) + geom_gate(gate.TYPEOFSIGNAL)
######END

##Plot the control sample that has non-fluorescing cells
plot(flowData.nondebris[[neg.signal]], c('FSC.A','FL1.A'), xlim=c(0,3e6), ylim=c(0,5e4),smooth=T)

negsignal.gate <- locator(100, type='l', col='red')
gm.3 <- matrix(,length(negsignal.gate$x),2)
colnames(gm.3) <- c('FSC.A','FL1.A')
gm.3[,1] <- negsignal.gate$x
gm.3[,2] <- negsignal.gate$y
gate.neg <- polygonGate(filterId="zerosignal",.gate=gm.3)

##Overlay and check the new gate
ggcyto(flowData.nondebris[c(neg.signal)], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + ggcyto_par_set(limits = list(x = c(0,1e5), y = c(0,3e4))) + geom_gate(gate.neg)

#Look at the gating on all the controls
ggcyto(flowData.nondebris[c(neg.signal,pos.signal,hi.signal)], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + ggcyto_par_set(limits = list(x = c(0,1e5), y = c(0,3e4))) + geom_gate(gate.neg) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)

#Examine across all samples
for(i in 1:round(length(flowData)/4,0)){
  plot <- ggcyto(flowData.nondebris, aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + ggcyto_par_set(limits = list(x = c(0,1e5), y = c(0,3e4))) + geom_gate(gate.neg) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = i)
  print(plot)
}


##Draw a new gate for gating gate for normally flourescent cells
plot(flowData.nondebris[[pos.signal]], c('FSC.A','FL1.A'), xlim=c(0,3e6), ylim=c(0,4e4),smooth=T)
polygon(negsignal.gate)

possignal.gate <- locator(100, type='l', col='blue')
gm.4 <- matrix(,length(possignal.gate$x),2)
colnames(gm.4) <- c('FSC.A','FL1.A')
gm.4[,1] <- possignal.gate$x
gm.4[,2] <- possignal.gate$y
gate.pos <- polygonGate(filterId="positivesignal",.gate=gm.4)

##Overlay and check the new gate with the old gate
ggcyto(flowData.nondebris[pos.signal], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + geom_gate(gate.neg) + geom_gate(gate.pos) + ggcyto_par_set(limits = list(x = c(0,1e5), y = c(0,3e4)))

#Look at the gating on all the controls
ggcyto(flowData.nondebris[c(neg.signal,pos.signal,hi.signal)], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + ggcyto_par_set(limits = list(x = c(0,1e5), y = c(0,3e4))) + geom_gate(gate.neg) + geom_gate(gate.pos) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)

#Examine across all samples
for(i in 1:round(length(flowData)/4,0)){
  plot <- ggcyto(flowData.nondebris, aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + ggcyto_par_set(limits = list(x = c(0,1e5), y = c(0,3e4))) + geom_gate(gate.neg) + geom_gate(gate.pos) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = i)
  print(plot)
}
  
##Draw a new gate for gating gate for hi flourescent cells
plot(flowData.nondebris[[hi.signal]], c('FSC.A','FL1.A'), xlim=c(0,3e6), ylim=c(0,5e5),smooth=T)
polygon(negsignal.gate)
polygon(possignal.gate)

hisignal.gate <- locator(100, type='l', col='blue')
gm.5 <- matrix(,length(hisignal.gate$x),2)
colnames(gm.5) <- c('FSC.A','FL1.A')
gm.5[,1] <- hisignal.gate$x
gm.5[,2] <- hisignal.gate$y
gate.hi <- polygonGate(filterId="positivesignal",.gate=gm.5)

##Overlay and check the new gate with the old gates
ggcyto(flowData.nondebris[hi.signal], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + geom_gate(gate.neg) + geom_gate(gate.pos) + geom_gate(gate.hi) 

#Look at the gating on all the controls
ggcyto(flowData.nondebris[c(neg.signal,pos.signal,hi.signal)], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + ggcyto_par_set(limits = list(x = c(0,1e5), y = c(0,3e4))) + geom_gate(gate.neg) + geom_gate(gate.pos) + geom_gate(gate.hi) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)

#Examine across all samples
for(i in 1:round(length(flowData)/4,0)){
  plot <- ggcyto(flowData.nondebris, aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + ggcyto_par_set(limits = list(x = c(0,1e5), y = c(0,3e4))) + geom_gate(gate.neg) + geom_gate(gate.pos) + geom_gate(gate.hi) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = i)
  print(plot)
}

#Save the gate information to an R data file
#Remeber to add any additional gate variables to your save file
rm(list=c("flowData")) 
save(pg.singlets, pg.nondebris, gate.neg, gate.pos, gate.hi, file=paste(name,"_gates_",Sys.Date(),".Rdata",sep=""))
