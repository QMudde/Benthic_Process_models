
# ==============================================================================
# parameters, output variables, etc..
# ==============================================================================

Macroalgae_get_Cbudget <- function(out, ..., 
              which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  Macroalgae_getbudget_all(out, ..., which = which, 
              func = Macroalgae_getbudgetC_one, args = sys.call())


Macroalgae_get_Nbudget <- function(out, ..., 
                                   which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  Macroalgae_getbudget_all(out, ..., which = which, 
                           func = Macroalgae_getbudgetN_one, args = sys.call())


Macroalgae_getbudget_all <- function(
      out, ..., 
      which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat"), 
      func = Macroalgae_getbudgetC_one, args){
  
  which <- match.arg(which, 
                     choices = c("All", "Rates", "Fluxes", "Losses", "Fluxmat"))
  
  if ("All" %in% which)
    which <- c("Rates", "Fluxes", "Losses", "Fluxmat", "Delta", "dC")
  
  ALL  <- list(out, ...)
  
  NM <- unlist(lapply(args[-1], as.character))
  if (! is.null(names(NM)))
    NM <- NM[!names(NM) == "which"]
  
  budg <- func(out)  
  budgFlux    <- budg$Fluxes
  budgRates   <- budg$Rates
  budgLoss    <- budg$Losses
  budgdC      <- budg$dC
  budgDelta   <- budg$Delta
  budgFluxmat <- budg$Fluxmat
  
  if (length(ALL) > 1) {
    
    budgFlux    <- unlist(budgFlux)
    budgRates   <- unlist(budgRates)
    budgFluxmat <- list(budgFluxmat)
    for ( i in 2:length(ALL)) {
      budg <- func(ALL[[i]])
      
      budgFlux    <- cbind(budgFlux,    unlist(budg$Fluxes))
      budgRates   <- cbind(budgRates,   unlist(budg$Rates))
      budgLoss    <- c(budgLoss,        budg$Losses)
      #    budgPerturb <- cbind(unlist(budgPerturb), budg$Perturb)
      budgdC      <- cbind(budgdC,          budg$dC)
      budgDelta   <- c(budgDelta,        budg$Delta)
      budgFluxmat[[i]]  <- budg$Fluxmat
    } 
    
    cn <- rep(names(budg$Fluxes), each = 4)
    nc <- nchar(cn)
    cn <- paste (cn, c("surf", "deep", "perturb", "net"), sep="")
    
    rownames(budgFlux) <- cn
  }
  
  budg <- list()
  if ("Rates" %in% which){
    if(length(NM) > 1) colnames(budgRates) <- NM
    budg$Rates <- budgRates
  }
  if ("Fluxes" %in% which){
    if(length(NM) > 1) colnames(budgFlux) <- NM
    budg$Fluxes <- budgFlux
  }
  if ("Losses" %in% which){
    if(length(NM) > 1) names(budgLoss) <- NM
    budg$Losses <- budgLoss
  }
  if ("Fluxmat" %in% which){
    if(length(NM) > 1) names(budgFluxmat) <- NM
    budg$Fluxmat <- budgFluxmat
  }
  if ("dC" %in% which){
    if(length(NM) > 1) colnames(budgdC) <- NM
    budg$dC <- budgdC
  }
  if ("Delta" %in% which){
    if(length(NM) > 1) names(budgDelta) <- NM
    budg$Delta <- budgDelta
  }
  
  return(budg)
}

Macroalgae_getbudgetC_one <- function(out) {
  
  parms <- Macroalgae_get_parms(out, as.vector = TRUE)  
  
  dV  <- dVal(out)
  OUT <- rbind(
    Macroalgae_get_states(out)[, c("names", "mean_value")],
    Macroalgae_get_vars  (out)[, c("names", "mean_value")]
  )
  nm <- OUT$names
  OUT <- as.list(OUT[, 2])
  names(OUT) <- nm
  
  Fluxes <- data.frame(  
    F.DICtoRESERVE_C            = OUT$Photosynthesis,
    F.RESERVE_CtoDIC            = OUT$Basal_Resp + OUT$Growth_Resp,
    F.RESERVE_CtoSTRUCT_C       = OUT$Growth,
    F.RESERVE_CtoEXT            = OUT$Mortality*OUT$RC_SCr,
    F.STRUCT_CtoEXT             = OUT$Mortality,
    F.STRUCT_CtoDIC             = OUT$Regression
)
  
  Rates <- data.frame(
    Photosynthesis            = OUT$Photosynthesis,
    Respiration               = OUT$Basal_Resp + OUT$Growth_Resp,
    Growth_C                  = OUT$Growth,
    Mortality_C               = OUT$Mortality*OUT$RC_SCr + OUT$Mortality,
    Regression_C              = OUT$Regression
  )  
  
  # derivatives
    names        <- c("DIC", "RESERVE_C", "STRUCT_C", "Ext")
    # Rows = from, columns = to
    CarbonFlows  <- matrix(nrow = 4, ncol = 4, byrow = TRUE, 
                           data = with(Fluxes, c(
                             
# DIC                RESERVE_C             STRUCT_C       EXT
  0,                 F.DICtoRESERVE_C,    0,              0,               # DIC
  F.RESERVE_CtoDIC,        0,  F.RESERVE_CtoSTRUCT_C,    F.RESERVE_CtoEXT, # RESERVE_C
  F.STRUCT_CtoDIC,         0,             0,             F.STRUCT_CtoEXT,  # STRUCT_C
      0,                   0,             0,              0 )              # EXT
))
    
  rownames(CarbonFlows) <- colnames(CarbonFlows) <- names
  dC <- unlist(dV[c("RESERVE_C", "STRUCT_C", "C_CLOSING")])
  # derivatives
  return(list(Fluxes = Fluxes, Rates = Rates, dC = c(dC, sum = sum(dC)), #Perturb = Perturb,
              Losses = OUT$C_CLOSING, Fluxmat = CarbonFlows))                 
  
}

Macroalgae_getbudgetN_one <- function(out) {
  
  parms <- Macroalgae_get_parms(out, as.vector = TRUE)  
  
  dV  <- dVal(out)
  OUT <- rbind(
    Macroalgae_get_states(out)[, c("names", "mean_value")],
    Macroalgae_get_vars  (out)[, c("names", "mean_value")]
  )
  nm <- OUT$names
  OUT <- as.list(OUT[, 2])
  names(OUT) <- nm
  N_return  = OUT$Mortality*OUT$RN_SCr
  N_to_DIN  = (OUT$Mortality+OUT$Regression)*parms[["SN_SCr"]]
  
  Fluxes <- data.frame(  
    F.NO3toRESERVE_N            = OUT$NO3_upt,
    F.NH4toRESERVE_N            = OUT$NH4_upt,
    F.RESERVE_NtoSTRUCT_N       = OUT$Growth*parms[["SN_SCr"]],
    F.RESERVE_NtoNO3            = N_return*0.5,
    F.RESERVE_NtoNH4            = N_return*0.5,
    F.STRUCT_NtoNH4             = N_to_DIN*0.5,
    F.STRUCT_NtoNO3             = N_to_DIN*0.5,
    F.NO3toExt                  = (OUT$ExtNO3 - OUT$NO3)*parms[["dilution"]]*parms[["depth"]],
    F.NH4toExt                  = (OUT$ExtNH4 - OUT$NH4)*parms[["dilution"]]*parms[["depth"]]
  )
  
  Rates <- data.frame(
    NO3uptake                 = OUT$NO3_upt,
    NH4uptake                 = OUT$NH4_upt,
    Growth_N                  = OUT$Growth*parms[["SN_SCr"]],
    Mortality_N               = OUT$Mortality*(OUT$RN_SCr + parms[["SN_SCr"]]),
    Regression_N              = OUT$Regression*parms[["SN_SCr"]],
    NO3_env_exchange          = (OUT$ExtNO3 - OUT$NO3)*parms[["dilution"]],
    NH4_env_exchange          = (OUT$ExtNH4 - OUT$NH4)*parms[["dilution"]]
  )  
  
  # derivatives
  names        <- c("NO3", "NH4", "RESERVE_N", "STRUCT_N", "Ext")
  # Rows = from, columns = to
  NitrogenFlows  <- matrix(nrow = 5, ncol = 5, byrow = TRUE, 
                         data = with(Fluxes, c(
                           
   #    NO3            NH4           RESERVE_N             STRUCT_N       EXT
     0,                 0,          F.NO3toRESERVE_N,        0,       F.NO3toExt,  # NO3
     0,                 0,          F.NH4toRESERVE_N,        0,       F.NH4toExt,  # NH4
F.RESERVE_NtoNO3, F.RESERVE_NtoNH4,          0,   F.RESERVE_NtoSTRUCT_N,     0,  # RESERVE_N
 F.STRUCT_NtoNO3,  F.STRUCT_NtoNH4,          0,              0,              0,  # STRUCT_N
  0,                    0,                   0,              0,              0 ) # EXT
                         ))
  
  rownames(NitrogenFlows) <- colnames(NitrogenFlows) <- names
  dC <- unlist(dV[c("NO3", "NH4", "RESERVE_N", "STRUCT_C", "N_CLOSING")])
  dC["STRUCT_C"] <- dC["STRUCT_C"]*parms[["SN_SCr"]]
  names(dC)[4] <- "STRUCT_N"
  
  # derivatives
  return(list(Fluxes = Fluxes, Rates = Rates, dC = c(dC, sum = sum(dC)), 
              Losses = OUT$N_CLOSING, Fluxmat = NitrogenFlows))                 
  
}

# changes in state variables over the entire period
dVal <- function(out)  {# takes into account unequal timing 
  toty <- c("time", names(Macroalgae_get_states(as.vector = TRUE))) 
  
  OUT <- out[c(1, nrow(out)), toty]
  as.list((OUT[2,]-OUT[1,])/(OUT[2,1]-OUT[1,1]))
}
