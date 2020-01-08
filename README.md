# flow
These scripts are used for processing flow cytometry and FACS data in the Gresham Lab. 
They are designed to enable quantitative analysis of cell fluorescence data. 

There are three scripts.

1. flow_gating.R

which performs gating for:

* doublets 
* debris 
* fluorescence 

using an R script and save gates as .Rdata file.
    
2. Import .Rdata file containing gates and render Rmarkdown file in order to save report of work.
