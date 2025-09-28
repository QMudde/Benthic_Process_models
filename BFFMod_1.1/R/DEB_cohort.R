
# ==============================================================================
# ==============================================================================
# DEB model of an individual - Carbon-based
# ==============================================================================
# ==============================================================================

mussel_run_debind <- function(
                       parms       = list(),   # subset of parameters that are overruled
                       times       = 0:365, 
                       yini        = NULL, 
                       f_Phyto     = 10, 
                       f_Temp      = 20,
                       f_Detritus  = 10,
                       f_Sim       = 1,
                       f_O2        = 360,
                       f_tau       = 0.01){  # if FALSE: Algae, Poc, Sim, O2, tau are imposed
  # ---------------------------------------------------------------------------
  # parameter vector matched to default parameters
  # ---------------------------------------------------------------------------
  
  P   <- getParms_deb(parms)

  # ---------------------------------------------------------------------------
  # Create forcing functions
  # ---------------------------------------------------------------------------

  forcings <- list(
    f_Phyto    = CreateForcings(f_Phyto,     "f_Phyto",    times), 
    f_Temp     = CreateForcings(f_Temp,      "f_Temp",     times), 
    f_Detritus = CreateForcings(f_Detritus,  "f_Detritus", times), 
    f_Sim      = CreateForcings(f_Sim,       "f_Sim",      times), 
    f_O2       = CreateForcings(f_O2,        "f_O2",       times),
    f_tau      = CreateForcings(f_tau,       "f_tau",      times))
  
  # ---------------------------------------------------------------------------
  # initial conditions
  # ---------------------------------------------------------------------------
  max_cohort = 1
  
  if (is.null(yini)) {
    weight_settle <- P["weight_settle"]
    pReserve_settle <- P["pReserve_settle"]
    yini <- c(RESERVE = weight_settle*pReserve_settle, 
              STRUCT = weight_settle*(1-pReserve_settle), 
              REPROD = 0.)
  }
    
  else if (length(yini) != 3) 
    stop ("length of yini should be equal to 3: RESERVE, STRUCT, REPROD")
  
  # ---------------------------------------------------------------------------
  # names
  # ---------------------------------------------------------------------------
  
  # output variables, includes forcing function variables
  # variables defined in each cohort (out1D) and single values (out0D)
    
    outnames <- .Deb_cohort$out1D$names
    
    # remove these variables
    inot     <- which(.Deb_cohort$out1D$names %in% c("C_m2", "Settle_rate", "Mortality"))
    outnames <- outnames[-inot]
    
    # add forcing function names
    outnames <- c(outnames,  
                  .Deb_cohort$forc$names)
  
  # State variable names
    ynames <- .Deb_cohort$y1D$names[1:3]
  
    units <- .Deb_cohort$out1D$units
    units <- units[-inot]
    units <- c(units,  
                  .Deb_cohort$forc$units)
    units <- c(.Deb_cohort$y1D$units[1:3], units)
    
    yini <- as.double(yini)
    names(yini)  <- ynames
  
  # add number of cohorts (0) and dynamic interaction switch (no) to par vector
   
  PP <- as.double(c(0., 0., P)    )
  
  if (length(times) > 1)  # Dynamic model
    out1 <- ode(y = yini, times = times, func = "musselmod", method = "vode",
                parms = PP, initfunc = "initmuspar", names = ynames, 
                initforc = "initmusforc", dllname = "BFFMod", forcings = forcings,
                nout = length(outnames), outnames = outnames, 
                rtol=1e-7, atol=1e-7)
  
  
  else {
    out1 <- DLLfunc(y = yini, times = times, func = "musselmod", 
                    parms = PP, initfunc = "initmuspar", 
                    initforc = "initmusforc", dllname = "BFFMod", forcings = forcings,
                    nout = length(outnames), outnames = outnames)
    names(out1$dy) <- paste("d", ynames, sep="")
  }
  
  attributes(out1)$units <- units
  attributes(out1)$max_cohort  <- 1
  attributes(out1)$dynamic_environment <- FALSE 
  attributes(out1)$parms       <- P
  attributes(out1)$env_name    <- ".Deb_cohort"
  attributes(out1)$model       <- "mussel_run_debind"
  class(out1) <- c("dtDyn", class(out1))  # dynamic digital twin application
  return(out1) 
}     ## end of mussel_run_debind


# ==============================================================================
# ==============================================================================
# DEB model of a number of cohorts, including the interaction with the environment
# ==============================================================================
# ==============================================================================

mussel_run_debcohort <- function(
                          parms       = list(),   # subset of parameters that are overruled
                          times       = 0:365, 
                          yini        = NULL, 
                          f_Phyto     = 10, 
                          f_Temp      = 20,
                          f_Detritus  = 10,
                          f_Sim       = 1,
                          f_O2        = 360,
                          f_tau       = 0.01,
                          max_cohort  = 12,
                          dynamic_environment = TRUE){  # if FALSE: Algae, Poc, Sim, O2 are imposed
  
  # ---------------------------------------------------------------------------
  # parameter vector matched to default parameters
  # ---------------------------------------------------------------------------
  P <- getParms_deb(parms)
  
  # ---------------------------------------------------------------------------
  # Create forcing functions
  # ---------------------------------------------------------------------------
  
  forcings <- list(
    f_Phyto    = CreateForcings(f_Phyto,     "f_Phyto",    times), 
    f_Temp     = CreateForcings(f_Temp,      "f_Temp",     times), 
    f_Detritus = CreateForcings(f_Detritus,  "f_Detritus", times), 
    f_Sim      = CreateForcings(f_Sim,       "f_Sim",      times), 
    f_O2       = CreateForcings(f_O2,        "f_O2",       times),
    f_tau      = CreateForcings(f_tau,       "f_tau",      times))
  
  # ---------------------------------------------------------------------------
  # initial conditions
  # ---------------------------------------------------------------------------
  
  if (is.null(yini)) {
    
    RESERVE    = rep(0., times=max_cohort)  # [mmol.C/ind] reserve C
    STRUCT     = rep(0., times=max_cohort)  # [mmol.C/ind] structure C
    REPROD     = rep(0., times=max_cohort)  # [mmol.C/ind] carbon allocated to maturation & reproduction
    POPULATION = rep(0., times=max_cohort)  # [ind/m2] population density per cohort         
    
    # PELAGIC LARVAE 
    LARV_DENS = 1e5                          # initial seeding density of pelagic larvae [ind/m3] 
    LARV_BIOM = LARV_DENS*P[["weight_egg"]]  # biomass of pelagic larvae [mmol.C/m3]
    
    # use forcing functions for these
    PHYTO     = approx(forcings[[1]][,1], forcings[[1]][,2], xout = times[1])$y
    DETRITUS  = approx(forcings[[3]][,1], forcings[[3]][,2], xout = times[1])$y   
    Sim       = approx(forcings[[4]][,1], forcings[[4]][,2], xout = times[1])$y
    O2        = approx(forcings[[5]][,1], forcings[[5]][,2], xout = times[1])$y
    
    yini <- c(RESERVE, STRUCT, REPROD, POPULATION, 
              LARV_DENS, LARV_BIOM, 
              PHYTO, DETRITUS, Sim, O2) 
  } else {
      if (length(yini) != 4*max_cohort + 6) 
        stop ("length of yini should be equal to 4*max_cohort + 6")
      if (! dynamic_environment){  # Will stay the same
        yini["PHYTO"]    = approx(forcings[[1]][,1], forcings[[1]][,2], xout = times[1])$y
        yini["DETRITUS"] = approx(forcings[[3]][,1], forcings[[3]][,2], xout = times[1])$y
        yini["Sim"]      = approx(forcings[[4]][,1], forcings[[4]][,2], xout = times[1])$y
        yini["O2"]       = approx(forcings[[5]][,1], forcings[[5]][,2], xout = times[1])$y
        
      }
    }
  
  # ---------------------------------------------------------------------------
  # names
  # ---------------------------------------------------------------------------
  
  # output variables: 
  # variables defined in each cohort (out1D) and single values (out0D)
  if (max_cohort > 1)
    outnames <- c(paste(rep(.Deb_cohort$out1D$names, each = max_cohort), 
                        1:max_cohort, sep="_"),
                        .Deb_cohort$out0D$names, 
                        .Deb_cohort$forc$names)
  else 
    outnames <- c(.Deb_cohort$out1D$names, 
                  .Deb_cohort$out0D$names,  
                  .Deb_cohort$forc$names)
  
  units <- c(rep(.Deb_cohort$out1D$units, each = max_cohort), 
             .Deb_cohort$out0D$units, 
             .Deb_cohort$forc$units)
  units <- c(.Deb_cohort$y1D$units[1:3], units)
  
  # State variable names
  # states in each cohort (y1D) and single values (y0D)
  if (max_cohort > 1)
    ynames <- c(paste(rep(.Deb_cohort$y1D$names, each = max_cohort), 
                          1:max_cohort, sep="_"),  
                      .Deb_cohort$y0D$names)                      
  else
    ynames <- c(.Deb_cohort$y1D$names, .Deb_cohort$y0D$names)                      
  
  yini <- as.double(yini)
  names(yini)  <- ynames
  
  # add number of cohorts and dynamic interaction switch to par vector
  # 
  PP <- as.double(c(max_cohort, dynamic_environment, P)    )
  
  if (length(times) > 1)  # Dynamic model   euler method requires hini
    out1 <- ode(y = yini, times = times, func = "musselcohort", method = "euler", #lsodes # vode
                parms = PP, initfunc = "initmuspar", names = ynames, hini = 0.01,
                initforc = "initmusforc", dllname = "BFFMod", forcings = forcings,
                nout = length(outnames), outnames = outnames, maxsteps = 1e5,
                rtol=1e-7, atol=1e-7)
  
  
  else {
    out1 <- DLLfunc(y = yini, times = times, func = "musselcohort", 
                    parms = PP, initfunc = "initmuspar", 
                    initforc = "initmusforc", dllname = "BFFMod", forcings = forcings,
                    nout = length(outnames), outnames = outnames)
    colnames(out1[[1]]) <- ynames
  }
  attributes(out1)$units       <- units
  attributes(out1)$max_cohort  <- max_cohort
  attributes(out1)$spawn_seasons <- 1 + (P[["spawn_peak_2"]]>0) # 1 or 2   # NEW
  attributes(out1)$dynamic_environment <- dynamic_environment 
  attributes(out1)$parms       <- P
  attributes(out1)$env_name    <- ".Deb_cohort"
  attributes(out1)$model       <- "mussel_run_debcohort"
  class(out1) <- c("dtDyn", class(out1))  # dynamic digital twin application
  return(out1) 
}     ## end of mussel_run_debcohort


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
# Function to get all parameter values for Deb model
# ==============================================================================

getParms_deb <- function(parms){
  P   <- as.vector(.Deb_cohort$parms$default)      # default parameters
  names(P) <- .Deb_cohort$parms$names
  
  np  <- names(parms)
  if (any (!np %in% names(P))) 
    stop("parameter ", np[(!np %in% names(P))], "not known")
  
  P[np] <- unlist(parms)        # overrule values
  P
}

