
# ==============================================================================
# parameters, output variables, etc..
# ==============================================================================

Macroalgae_get_parms <- function(out, as.vector = FALSE, which = NULL){
  
  if (missing(out))
    Parms <- .Macroalgae$parms
  
  else if (length(out) == 0)     
    Parms <- .Macroalgae$parms
   
  else if (inherits(out, "dtMacroalgae")){
     Parms <- .Macroalgae$parms
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

Macroalgae_get_var_or_forc <- function(out, which = NULL, DF){
  
  if (! is.null(which))
    N0D <- subset(DF, subset = names %in% which)
  else 
    N0D <- DF
  
  if (!missing(out)) {
    if (attributes(out)$model == "Macroalgae_run")
      N0D <- subset(N0D, subset = names %in% colnames(out))
    if (inherits(out, "dtMacroalgae"))
      N0D$mean_value <- colMeans(out[, N0D$names])
  }
  row.names(N0D)  <- NULL
  N0D
}

Macroalgae_get_vars <- function(out, as.vector = FALSE, which = NULL){
  V <- Macroalgae_get_var_or_forc(out, which = which, DF = .Macroalgae$out)
  if (as.vector) V <- V$names
  V
}

Macroalgae_get_states <- function(out, as.vector = FALSE, which = NULL){
  V <- Macroalgae_get_var_or_forc(out, which = which, DF = .Macroalgae$y)
  if (as.vector) {
    nm <- V$names
    V <- V$initial
    names(V) <- nm
  }
  V
}

Macroalgae_get_forcings <- function(out, as.vector = FALSE, which = NULL){
  V <- Macroalgae_get_var_or_forc(out, which = which, DF = .Macroalgae$forc)
  if (as.vector) V <- V$names
  V
}
