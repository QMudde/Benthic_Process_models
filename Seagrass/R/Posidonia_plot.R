
# ==============================================================================
# Plot function for class dtPosidonia that also adds the units to the y labels.
# ==============================================================================

plot.dtPosidonia <- function(x, ..., select = NULL, which = select, 
                    ylab = NULL, lty = 1, las = 1){
  
  W <- which
  if (is.null(which))
    which <- 1:(ncol(x)-1)
  
  else if (is.character(which[1])){  # find the variable
    which <- unlist(lapply(which, 
                           FUN=function(X) which(colnames(x)[-1] %in% X))) 
    if(length(which) != length(W))
      stop ("Cannot find variable to plot: ", 
            paste (W[! W %in% colnames(x) [which+1]], collapse = ", "))
  
  }
  if (is.null(ylab)) 
    ylab <- attributes(x)$units[which]

  classx <- class(x)   
  class(x) <- classx[-which(classx == "dtPosidonia")]
  
  plot(x, ..., which = which, ylab = ylab, lty = lty, las = las)  
}

#setOldClass("dtPosidonia")
#setMethod("plot", "dtPosidonia", plot.dtPosidonia)

