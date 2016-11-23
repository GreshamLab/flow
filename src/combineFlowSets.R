
combineFlowSets <- function(setlist) {
  if (!is.list(setlist)) { error("gimme a list yo") }
  for ( i in 1:length(setlist)) {
    if (exists("outPhenoData")) {
      tmpPhenoData <- pData(setlist[[i]])
      tmpPhenoData$name <- paste0(tmpPhenoData$name,"_",i) 
      rownames(tmpPhenoData) <- tmpPhenoData$name
      outPhenoData <- rbind(outPhenoData,tmpPhenoData)
    } else {
      tmpPhenoData <- pData(setlist[[i]])
      tmpPhenoData$name <- paste0(tmpPhenoData$name,"_",i) 
      rownames(tmpPhenoData) <- tmpPhenoData$name
      outPhenoData <- tmpPhenoData
    }
  }
  outList <- list()
  for ( i in 1:length(setlist)) {
    sampleNames(setlist[[i]]) <- 
      paste0(sampleNames(setlist[[i]]),"_",i)
    for ( j in 1:length(setlist[[i]]) ) {
      outList[[sampleNames(setlist[[i]])[j]]] <- setlist[[i]][[j]]
    }
  }
  #return(flowSet(outList,phenoData=outPhenoData))
#wtf is this environmental locks shit?
  return(outPhenoData)
}
