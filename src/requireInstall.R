
# hey I haven't tested this completely yet

requireInstall <- function(packageName,isBioconductor=F) {
  exitVal <- require(packageName,character.only=T)
  if ( !exitVal ) {
    print(paste0("You don't have ",packageName," accessible, ",
      "I'm gonna install it"))
    if (isBioconductor) {
      source("http://bioconductor.org/biocLite.R")                        
      exitVal <- biocLite(packageName,
        suppressUpdates=T,suppressAutoUpdate=T,ask=F)
    } else {
      exitVal <- install.packages(packageName,
                   repos="http://cran.us.r-project.org")
    }
  }
  if (exitVal==T) { 
    return(1) 
  } else { 
    stop(paste0("error in loading/installing pkg ",packageName)) 
  }
}
