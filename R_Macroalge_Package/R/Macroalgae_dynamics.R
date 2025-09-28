
# ==============================================================================
# ==============================================================================
# single DEB model
# ==============================================================================
# ==============================================================================

Macroalgae_run <- function(
                       parms         = list(),   # subset of parameters that are overruled
                       times         = 0:365, 
                       yini          = NULL, 
                       f_Temperature = 20,       # forcing functions
                       f_Current     = 1, 
                       f_Light       = 100,      # W/m2
                       f_PHYTO       = 10,       # mmol/m3
                       f_ExtNO3      = 20,
                       f_ExtNH4      = 1,
                       f_ExtPO4      = 1,
                       f_ExtO2       = 200){  
  # ---------------------------------------------------------------------------
  # parameter vector matched to default parameters
  # ---------------------------------------------------------------------------
  
  PP <- Macroalgae_get_parms(as.vector = TRUE)
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
    f_Current     = CreateForcings(f_Current,     "f_Current",     times),
    f_Light       = CreateForcings(f_Light,       "f_Light",       times),
    f_PHYTO       = CreateForcings(f_PHYTO,       "f_PHYTO",       times),
    f_ExtNO3      = CreateForcings(f_ExtNO3,      "f_ExtNO3",      times),
    f_ExtNH4      = CreateForcings(f_ExtNH4,      "f_ExtNH4",      times),
    f_ExtPO4      = CreateForcings(f_ExtPO4,      "f_ExtPO4",      times),
    f_ExtO2       = CreateForcings(f_ExtO2,       "f_ExtO2",       times))
  
  # ---------------------------------------------------------------------------
  # get initial conditions
  # ---------------------------------------------------------------------------
  
  yini <- getYini_ma(yini)
  
  # ---------------------------------------------------------------------------
  # output names
  # ---------------------------------------------------------------------------
  
  # output variables, includes forcing function variables
  # variables defined in each cohort (out1D) and single values (out0D)
    
  outnames <- .Macroalgae$out$names
    
  # State variable names
  ynames <- .Macroalgae$y$names
  
  units  <- c(.Macroalgae$y$units, .Macroalgae$out$units)

  yini   <- as.double(yini)
  names(yini)  <- ynames
  
  if (length(times) > 1)  # Dynamic model
    out1 <- ode(y = yini, times = times, 
                parms = PP, forcings = forcings,
                names = ynames, nout = length(outnames), outnames = outnames, 
                func = "macroalgae", initfunc = "initmacroalgae", 
                initforc = "initforcm", dllname = "dtMacroalgae", 
                method = "lsodes", rtol=1e-12, atol=1e-12) # vode? (10?/10?)
  
  else {
    out1 <- DLLfunc(y = yini, times = times, 
                    parms = PP, forcings = forcings,
                    func = "macroalgae", initfunc = "initmacroalgae", 
                    initforc = "initforcm", dllname = "dtMacroalgae", 
                    nout = length(outnames), outnames = outnames)
    names(out1$dy) <- paste("d", ynames, sep="")
  }
  
  attributes(out1)$units <- units
  attributes(out1)$max_cohort  <- 1
  attributes(out1)$dynamic_environment <- FALSE 
  attributes(out1)$parms       <- PP
  attributes(out1)$env_name    <- ".Macroalgae"
  attributes(out1)$model       <- "Macroalgae_run"
  class(out1) <- c("dtMacroalgae", class(out1))  # dynamic digital twin application
  return(out1) 
}     ## end of Macroalgae_run



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
# Function to get all parameter values for Macroalgae model
# ==============================================================================

getParms_ma <- function(parms){
  P        <- as.vector(.Macroalgae$parms$default)      # default parameters
  names(P) <- .Macroalgae$parms$names
  
  np  <- names(parms)
  if (any (!np %in% names(P))) 
    stop("parameter ", np[(!np %in% names(P))], "not known")
  
  P[np] <- unlist(parms)        # overrule values
  P
}

# ==============================================================================
# Function to get all initial conditions for Macroalgae model
# ==============================================================================

getYini_ma <- function(yini){
  P   <- as.vector(.Macroalgae$y$initial)      # default initial conditions
  names(P) <- .Macroalgae$y$names
  
  if (is.null(yini))
    return(P)
  
  if (inherits(yini, "dtMacroalgae")) {
    yini <- Macroalgae_get_states(yini)[,"final_value"]
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
