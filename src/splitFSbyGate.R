splitFSbyGate <- function(x,aGate=NULL) {

  if ("filter"%in%is(aGate)) {
    z <- match.call()["aGate"]
    tmp <- sub(".+\\[(.+)]]","\\1",z)
    if (exists(tmp)) { z <- get(tmp) } # might not need exists here
    gateName <- gsub("[\"]","",z)
    print(paste0("applying gate ",gateName))

    ingate  <- Subset(x,aGate)
    pData(ingate)[gateName] <- T
    sampleNames(ingate) <- paste0(sampleNames(ingate),"_in_",gateName)

    outgate <- Subset(x,!aGate)
    pData(outgate)[gateName] <- F
    sampleNames(outgate) <- paste0(sampleNames(outgate),"_out_",gateName)

    outPhenoData <- rbind(pData(ingate),pData(outgate))

    outList <- list()
    for ( i in list(ingate,outgate) ) {
      for ( j in 1:length(i)) {
        outList[[sampleNames(i)[j]]] <- i[[j]]
      }
    }

    outSet <- flowSet(outList)
    pData(outSet) <- outPhenoData
    return(outSet)   

  } else {
    print("didn't split, cuz no gate")
    return(x)
  }
}
