
# ==============================================================================
# parameters, output variables, etc..
# ==============================================================================

Posidonia_get_parms <- function(out, as.vector = FALSE, which = NULL){
  
  if (missing(out))
    Parms <- .Posidonia$parms
  
  else if (length(out) == 0)     
    Parms <- .Posidonia$parms
   
  else if (inherits(out, "dtPosidonia")){
     Parms <- .Posidonia$parms
     Parms$value <- attr(out, "parms")
    }
  else stop("object 'out' not supported")
  
  if (as.vector) {
    nm           <- Parms$names 
    if ("value" %in% names(Parms))
      Parms        <- Parms$value
    else
      Parms        <- Parms$default
    names(Parms) <- nm
  }
  
  if (! is.null(which))
    Parms <- subset(Parms, subset = names %in% which)
  
  return(Parms)
}

Posidonia_get_var_or_forc <- function(out, which = NULL, DF){
  
  if (! is.null(which))
    N0D <- subset(DF, subset = names %in% which)
  
  else 
    N0D <- DF
  
  if (!missing(out)) {
    
    if (attributes(out)$model == "Posidonia_run")
      N0D <- subset(N0D, subset = names %in% colnames(out))
    
    if (inherits(out, "dtPosidonia")){
      N0D$initial_value <- out[1, N0D$names]
      N0D$mean_value    <- colMeans(out[, N0D$names])
      N0D$final_value   <- out[nrow(out), N0D$names]
      N0D$min_value     <- apply(out[ , N0D$names], MARGIN = 2, FUN = min)
      N0D$max_value     <- apply(out[ , N0D$names], MARGIN = 2, FUN = max)
    }
  }
  
  row.names(N0D)  <- NULL
  N0D
}


Posidonia_get_vars <- function(out, as.vector = FALSE, which = NULL){
  V <- Posidonia_get_var_or_forc(out, 
                                 which = which, 
                                 DF = .Posidonia$out)
  if (as.vector) V <- V$names
  V
}

  
Posidonia_get_states <- function(out, as.vector = FALSE, which = NULL){
  
  V <- Posidonia_get_var_or_forc(out, 
                                 which = which, 
                                 DF = .Posidonia$y)
  if (as.vector) {
    nm  <- V$names
    V   <- V$default
    names(V) <- nm
  }
  V
}

Posidonia_get_forcings <- function(out, as.vector = FALSE, which = NULL){
  V <- Posidonia_get_var_or_forc(out, which = which, DF = .Posidonia$forc)
  if (as.vector) V <- V$names
  V
  
}
