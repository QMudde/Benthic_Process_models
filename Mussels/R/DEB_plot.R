# ==============================================================================
# ==============================================================================
# Plot function for class dt that also adds the units to the y labels.
# ==============================================================================
# ==============================================================================

plot.dtDyn <- function(x, ..., select = NULL, which = select, 
                       ylab = NULL, lty = 1, las = 1){
  
  WW <- which
  
  if (is.null(WW))
    WW <- 1:(ncol(x)-1)
  
  # find the variable, so that its units can be found
  else if (is.character(WW[1]))  
    WW <- unlist(lapply(WW, 
                        FUN=function(X) which(colnames(x)[-1] %in% X))) 
  
  if (is.null(ylab)) 
    ylab <- attributes(x)$units[WW]
  
  classx <- class(x)   
  class(x) <- classx[-which(classx == "dtDyn")]
  
  plot(x, ..., which = WW, ylab = ylab, lty = lty, las = las)  
}

# ==============================================================================

matplot.mussel <- function(x, which = NULL, 
                       ylab = NULL, lty = 1, las = 1, 
                       legend = list(x = "topleft", cex=0.5), ...){
  
  WW <- which
  ldots <- list(...)
  
  
  if (is.null(WW))
    WW <- 1:4
  
  # find the variable, so that its units can be found
  else if (is.character(WW[1]))  
    WW <- unlist(lapply(WW, 
                        FUN=function(X) which(colnames(x)[-1] %in% X))) 
  Wnames <- colnames(x)[WW+1]
  if (! length(WW)){
    ncohort <- attributes(x)$max_cohort
    WW <- paste(which[1], 1:ncohort, sep="_")
    if (is.null(ldots$main)) 
      ldots$main <- which[1]
    WW <- unlist(lapply(WW, 
                        FUN=function(X) which(colnames(x)[-1] %in% X))) 
    Wnames <- 1:ncohort
    
  }    
  
  if (is.null(ldots$xlab))
    ldots$xlab <- "times"
  if (is.null(ylab)) 
    ylab <- attributes(x)$units[WW[1]]
  
  classx <- class(x)   
  class(x) <- classx[-which(classx %in% c("dtDyn", "deSolve"))]
  
  do.call("matplot", c(alist(x = x[,1], y = x[,(WW+1)], type="l", 
          ylab = ylab, lty = lty, las = las), ldots))
  
  if (! is.null(legend)){
    if (is.null(legend$x)) legend$x <- "topleft"
     do.call("legend", 
    c(alist(legend = Wnames, col = 1: length(Wnames), lty=lty), legend))
  }
}

