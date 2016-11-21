
fc2df <- function(x) {
  for (i in sampleNames(x)) {
    if (exists("returnDF")) {
      returnDF <- rbind(returnDF,
        data.frame(pData(x[i]),exprs(x[[i]]),row.names=NULL)
        )
    } else {
      returnDF <- data.frame(pData(x[i]),exprs(x[[i]]),row.names=NULL)
    }
  }
  return(returnDF)
}
