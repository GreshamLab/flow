 
gator2d <- function(datar,ehx,why,events=1e5,logit=T) {
  sub <- datar[sample(1:nrow(datar),events),]
  if (logit) {
    plot(x=log10(sub[,ehx]),y=log10(sub[,why]),pch=16,cex=0.2,
        main=paste0("keep clicking to approximately define the gate",
          "\nyou want, then right click when you're done"),
        xlim=c(range(log10(sub[,ehx])[is.finite(log10(sub[,ehx]))])),
        ylim=c(range(log10(sub[,why])[is.finite(log10(sub[,why]))])))
    gp <- locator(type="l",col="red") 
    g <- as.matrix(
      data.frame(x=10^gp$x,y=10^gp$y))
  } else {
    plot(x=sub[,ehx],y=sub[,why],pch=16,cex=0.2,
        main=paste0("keep clicking to approximately define the gate",
          "\nyou want, then right click when you're done"),
        xlim=c(range(sub[,ehx])),
        ylim=c(range(sub[,why])))
    gp <- locator(type="l",col="red") 
    g <- as.matrix(
      data.frame(x=gp$x,y=gp$y))
  }
  dev.off()
  names(g) <- c(ehx,why)
  return(polygonGate(filterId=paste0("",""),.gate=g))
}
