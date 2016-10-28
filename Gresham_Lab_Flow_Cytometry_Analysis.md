# Gresham Lab Flow Core Guide
`r Sys.info()[7]`  
`r Sys.Date()`  

#1) Write a detailed description of your experiment here.  If you still see this text it means that you have not described the experiment and whatever follows is meaningless.

#2) Be sure to include a relevant title

#3) Delete text containing points 1-3 

###This code is designed for use with the Accuri flow cytometer, which is equiped with the following lasers and filters

* Blue laser (488 nm)
  + FL1 filter = 514/20nm   GFP
  + FL3 filter = 575/25nm   YFP

* Yellow/green laser (552 nm)
  + FL2 filter = 610/20nm   mCherry, dtomato
  + FL4 filter = 586/15nm   DsRed

###It is also designed for use with the CGSB Aria, which has the following lasers and printers

###In order to run this code you need to predefine your gates using the gating.R script

###This script generates a summary of results followed by quality control plots


##Step 1: Load relevant libraries 

```r
# This is a function that just makes sure you have a package, or installs it for you without prompting

requireInstall <- function(packageName,isBioconductor=F) {
  if ( !try(require(packageName,character.only=T)) ) {
    print(paste0("You don't have ",packageName," accessible, ",
      "I'm gonna install it"))
    if (isBioconductor) {
      source("http://bioconductor.org/biocLite.R")                        
      biocLite(packageName)                                                 
    } else {
      install.packages("packageName")
    }
  }
  return(1)
}

#Load libraries
requireInstall("flowCore",isBioconductor=T)
```

```
## Loading required package: flowCore
```

```
## [1] 1
```

```r
requireInstall("flowViz",isBioconductor=T)
```

```
## Loading required package: flowViz
```

```
## Loading required package: lattice
```

```
## Warning: package 'lattice' was built under R version 3.1.3
```

```
## [1] 1
```

```r
#requireInstall("flowStats")
#requireInstall("Hmisc")
requireInstall("reshape2")
```

```
## Loading required package: reshape2
```

```
## [1] 1
```

```r
requireInstall("ggplot2")
```

```
## Loading required package: ggplot2
```

```
## Warning: package 'ggplot2' was built under R version 3.1.3
```

```
## [1] 1
```

```r
#requireInstall("flowWorkspace")
#requireInstall("ggcyto", isBioconductor=T)
#requireInstall("gridExtra")
```

###Step 2: Read in .fcs files, an Rdata file containing the gates sample sheet(s) that contains four columns with 
* column1 = Well
* column2 = Strain
* column3 = Staining
* column4 = Media
* column5 = Userdefined


```r
#Read in all data for analysis. Data should be in individual directories that contain .fcs files and a corresponding sample sheet with a generic format. FCS file names should be unaltered e.g AO1.fcs, A02.fcs, ...H12.fcs 
#An abitrary number of directories can be used named dir1, dir2, dir3...with a corresponding flowData.1, flowData.2, flowData.3...and sample.sheet.1, sample.sheet.2, sample.sheet.3...

#load the Rdata file containing the gates
load("gates.Rdata") 

#Define how many directories that will be analyzed
num <- 1

#Define the directory, or directories, containing your .fcs files using absolute path names 
dir1 <- "/Users/David/Google Drive/Gresham Lab_David/flow/flow cytometry"
#dir2 <- "/Users/David/Google Drive/Gresham Lab_David/GreshamLab Flow Cytometry 2016/2015_12_02"

#Read in all the fcs files in the directory, with alter.names changing "-" to "."  
flowData.1 <- read.flowSet(path = dir1, pattern=".fcs", alter.names=TRUE)
#flowData.2 <- read.flowSet(path = dir2, pattern=".fcs", alter.names=TRUE)

#Read in the sample sheet that should be in each directory that contains the .fcs files.  
sample.sheet.1 <- read.csv(paste(dir1, "SampleSheet.csv", sep="/"))
#sample.sheet.2 <- read.csv(paste(dir2, "SampleSheet.csv", sep="/"))
```



```r
#Check how many cells were counted in each fcs file
fsApply(flowData.1, each_col, length)[1:6]
```

```
## [1] 50000 50000 50000 50000 50000 50000
```

```r
total <- fsApply(flowData.1, each_col, length)[1:6]
#fsApply(flowData.2, each_col, length)

#Print the medians of data values for each measurement
fsApply(flowData.1, each_col, median)
```

```
##            FSC.A   SSC.A   FL1.A FL2.A FL3.A FL4.A    FSC.H    SSC.H
## A01.fcs 562697.0 77701.0  2566.0   276   756   251 871409.5 108347.0
## A02.fcs 637355.5 90199.5  9941.0   558  2242   493 969225.0 122714.5
## A03.fcs 579237.0 79793.0 19631.0  1146  4984  1002 898238.0 110155.0
## A04.fcs 546466.0 80743.0 20697.5  1267  5452  1090 852802.0 113254.0
## A05.fcs 552799.5 81623.0 25457.0  1641  7192  1423 859692.0 114994.5
## A06.fcs 524191.5 78656.0 28496.5  1862  8128  1614 820285.5 111785.0
##           FL1.H FL2.H FL3.H FL4.H Width Time
## A01.fcs  2846.0   185  1594   421    58  184
## A02.fcs  9728.0   419  2453   568    60  184
## A03.fcs 18645.0   953  4677   959    58  188
## A04.fcs 19810.0  1077  5120  1042    57  184
## A05.fcs 24078.0  1412  6609  1328    58  186
## A06.fcs 27025.5  1616  7456  1508    57  188
```

```r
#fsApply(flowData.2, each_col, median)
```


###Step 3: apply filters to data and generate plots showing the effect on filtering

```r
##Subset the data by applying sequential gates##

#this filters doublets
flowData.1.1 <- Subset(flowData.1, pg.singlets) 
fsApply(flowData.1.1, each_col, length)[1:6]
```

```
## [1] 42759 41726 42816 42973 42303 42593
```

```r
singlets <- fsApply(flowData.1.1, each_col, length)[1:6]
barplot(singlets/total, ylim=c(0,1), ylab = "Proportion singlet cells", names.arg=sample.sheet.1[,1])
```

![](Gresham_Lab_Flow_Cytometry_Analysis_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

```r
#this removes debris from the singlet cells
filterData.1.1.1 <- Subset(flowData.1.1, pg.nondebris) 
fsApply(filterData.1.1.1, each_col, length)[1:6]
```

```
## [1] 42204 41196 42282 42112 41538 41535
```

```r
non.debris <- fsApply(filterData.1.1.1, each_col, length)[1:6]
barplot(non.debris/total, ylim=c(0,1), ylab = "Proportion singlet and nondebris cells", names.arg=sample.sheet.1[,1])
```

![](Gresham_Lab_Flow_Cytometry_Analysis_files/figure-html/unnamed-chunk-4-2.png)<!-- -->

```r
#non.debris is the filtered data that will be used for all subsequent analyses

#this identifies nongfp cells
gfp.neg <- Subset(filterData.1.1.1, pg.nongfp) 
fsApply(gfp.neg, each_col, length)[1:6]
```

```
## [1] 41827  7876   133   214    13     7
```

```r
non.gfp <- fsApply(gfp.neg, each_col, length)[1:6]
barplot(non.gfp/non.debris, ylim=c(0,1), ylab = "Proportion cells with no GFP", names.arg=sample.sheet.1[,1])
```

![](Gresham_Lab_Flow_Cytometry_Analysis_files/figure-html/unnamed-chunk-4-3.png)<!-- -->

```r
#this identifies gfp cells
gfp.pos <- Subset(filterData.1.1.1, pg.gfp) 
fsApply(gfp.pos, each_col, length)[1:6]
```

```
## [1] 12412 40040 33378 26725 14695  6472
```

```r
gfp.cells <- fsApply(gfp.pos, each_col, length)[1:6]
barplot(gfp.cells/non.debris, ylim=c(0,1), ylab = "Proportion cells with GFP", names.arg=sample.sheet.1[,1])
```

![](Gresham_Lab_Flow_Cytometry_Analysis_files/figure-html/unnamed-chunk-4-4.png)<!-- -->

```r
#this identifies high GFP cells
gfp.hi <- Subset(filterData.1.1.1, pg.hi.gfp) 
fsApply(gfp.hi, each_col, length)[1:6]
```

```
## [1]     0  3407 37121 37752 39674 39601
```

```r
hi.gfp.cells <- fsApply(gfp.hi, each_col, length)[1:6]
barplot(hi.gfp.cells/non.debris, ylim=c(0,1), ylab = "Proportion cells with high GFP", names.arg=sample.sheet.1[,1])
```

![](Gresham_Lab_Flow_Cytometry_Analysis_files/figure-html/unnamed-chunk-4-5.png)<!-- -->


