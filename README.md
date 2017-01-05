# flow
Scripts for processing flow cytometry and FACS data

The overall strategy for using flow cytometry analysis is:

1. Perform gating for i) doublets, ii) debris and iii) fluorescence using an R script and save gates as .Rdata file.
2. Import .Rdata file containing gates and render Rmarkdown file in order to save report of work.
