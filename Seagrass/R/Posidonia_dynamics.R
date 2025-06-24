
# ==============================================================================
# ==============================================================================
# single DEB model
# ==============================================================================
# ==============================================================================

Posidonia_run <- function(
                       parms         = list(),   # subset of parameters that are overruled
                       times         = 0:365, 
                       yini          = NULL, 
                       f_Temperature = 20,       # forcing functions
                       f_Light       = 100, 
                       f_Daylength   = 12,
                       f_Wind        = 1,
                       f_Piston      = 1,
                       f_SatO2       = 0.360){  
  # ---------------------------------------------------------------------------
  # parameter vector matched to default parameters
  # ---------------------------------------------------------------------------
  
  PP <- Posidonia_get_parms(as.vector = TRUE)
  if (! is.null(parms)){
    ii <- which (!names(parms) %in% names(PP))
    if (length(ii))
      stop ("names of 'parms' are not known parameters: ", 
            paste(names(parms)[ii], collapse = ", "))
    PP[names(parms)] <- unlist(parms)
  }
  
  # ---------------------------------------------------------------------------
  # Create forcing functions
  # ---------------------------------------------------------------------------

  forcings <- list(
    f_Temperature = CreateForcings(f_Temperature, "f_Temperature", times), 
    f_Light       = CreateForcings(f_Light,       "f_Light",       times),
    f_Daylength   = CreateForcings(f_Daylength,   "f_Daylength",   times),
    f_Wind        = CreateForcings(f_Wind,        "f_Wind",        times), 
    f_Piston      = CreateForcings(f_Piston,      "f_Piston",      times), 
    f_SatO2       = CreateForcings(f_SatO2,       "f_satO2",       times))
  
  # ---------------------------------------------------------------------------
  # get initial conditions
  # ---------------------------------------------------------------------------
  
  yini <- getYini_sg(yini)
  
  # ---------------------------------------------------------------------------
  # output names
  # ---------------------------------------------------------------------------
  
  # output variables, includes forcing function variables
  # variables defined in each cohort (out1D) and single values (out0D)
    
  outnames <- .Posidonia$out$names
    
  # State variable names
  ynames <- .Posidonia$y$names
  
  units  <- c(.Posidonia$y$units, .Posidonia$out$units)

  yini   <- as.double(yini)
  names(yini)  <- ynames
  
  if (length(times) > 1)  # Dynamic model
    out1 <- ode(y = yini, times = times, 
                parms = PP, forcings = forcings,
                names = ynames, nout = length(outnames), outnames = outnames, 
                func = "posidonia", initfunc = "initposidonia", 
                initforc = "initforcp", dllname = "dtPosidonia", 
                method = "vode", rtol=1e-10, atol=1e-10)
  
  else {
    out1 <- DLLfunc(y = yini, times = times, 
                    parms = PP, forcings = forcings,
                    func = "posidonia", initfunc = "initposidonia", 
                    initforc = "initforcp", dllname = "dtPosidonia", 
                    nout = length(outnames), outnames = outnames)
    names(out1$dy) <- paste("d", ynames, sep="")
  }
  
  attributes(out1)$units <- units
  attributes(out1)$max_cohort  <- 1
  attributes(out1)$dynamic_environment <- FALSE 
  attributes(out1)$parms       <- PP
  attributes(out1)$env_name    <- ".Posidonia"
  attributes(out1)$model       <- "Posidonia_run"
  class(out1) <- c("dtPosidonia", class(out1))  # dynamic digital twin application
  return(out1) 
}     ## end of Posidonia_run



# ==============================================================================
# Function to generate forcing functions
# ==============================================================================


CreateForcings <- function(forc, name, times){

  if (is.function(forc))
    forc <- cbind(times, forc(times))
  
  else if (is.vector(forc) & length(forc) == 1)
    forc <- data.frame(times, forc)
  
  else if (is.matrix(forc) | is.data.frame(forc)){
    if (min(forc[,1]) > min(times) |
        max(forc[,1]) < max(times)  )
    stop ("temporal range in ", name, 
          "forcing function does not embrace 'times'")
    forc <- forc[,1:2]
  }
  
  if (any(is.na(forc))) 
    stop ("NAs in ", name, "forcing function")
  
  forc
} 

# ==============================================================================
# Function to get all parameter values for Posidonia model
# ==============================================================================

getParms_sg <- function(parms){
  P        <- as.vector(.Posidonia$parms$default)      # default parameters
  names(P) <- .Posidonia$parms$names
  
  np  <- names(parms)
  if (any (!np %in% names(P))) 
    stop("parameter ", np[(!np %in% names(P))], "not known")
  
  P[np] <- unlist(parms)        # overrule values
  P
}

# ==============================================================================
# Function to get all initial conditions for Posidonia model
# ==============================================================================

getYini_sg <- function(yini){
  P   <- as.vector(.Posidonia$y$default)      # default initial conditions
  names(P) <- .Posidonia$y$names
  
  if (is.null(yini))
    return(P)
  
  if (inherits(yini, "dtPosidonia")) {
    yini <- Posidonia_get_states(yini)[,"final_value"]
  }
  
  np  <- names(yini)
  if (is.null(np)){
    if (length (yini) != length(P))
      stop("state variable vector 'yni' does not have the correct length \n", 
           " yini should contain: ", paste(names(P), collapse = ", "))
    np <- 1:length(yini)
    
  } else if (any (!np %in% names(P))) 
    stop("state variable ", np[(!np %in% names(P))], " in yini not known - \n",
         " yini should contain: ", paste(names(P), collapse = ", "))
  
  P[np] <- unlist(yini)        # overrule values
  P
  
}
