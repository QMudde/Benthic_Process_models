
# ==============================================================================
# parameters, output variables, etc..
# ==============================================================================

Posidonia_get_Cbudget <- function(out, ..., 
              which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat")) 
  Posidonia_getbudget_all(out, ..., which = which, 
              func = Posidonia_getbudgetC_one, args = sys.call())



Posidonia_getbudget_all <- function(
      out, ..., 
      which = c("All", "Rates", "Fluxes", "Losses", "Fluxmat"), 
      func = Posidonia_getbudgetC_one, args){
  
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

Posidonia_getbudgetC_one <- function(out) {
  
  parms <- Posidonia_get_parms(out, as.vector = TRUE)  
  
  dV <- dVal(out)
  OUT <- rbind(
    Posidonia_get_states(out)[, c("names", "mean_value")],
    Posidonia_get_vars  (out)[, c("names", "mean_value")]
  )
  nm <- OUT$names
  OUT <- as.list(OUT[, 2])
  names(OUT) <- nm
  
  Fluxes <- data.frame(  
    F.DICtoLEAVES               = OUT$GrossPrimaryProduction - 
                                  OUT$LeafTotalRespiration,
    F.LEAVEStoDETRITUS          = OUT$LeafBreaking,
    F.LEAVEStoDEADLEAVES        = OUT$LeafMortality,
    F.ASSIMILATEStoLEAVES       = OUT$LeafAssimilation,
    F.LEAVEStoASSIMILATES       = OUT$GrossPrimaryProduction + 
                                  OUT$LeafNetReallocation,
    F.DEADLEAVEStoDETRITUS      = OUT$DeadLeafAbscission + 
                                  OUT$DeadLeafBreaking,
    F.DEADLEAVEStoDIC           = OUT$DeadLeafDecay,
    F.RHIZOMEtoDIC              = OUT$RhizomeTotalRespiration,
    F.ASSIMILATEStoRHIZOMES     = OUT$RhizomeAssimilation,
    F.RHIZOMEStoASSIMILATES     = OUT$RhizomeMobilisation,
    F.RHIZOMEStoDEADBELOWGROUND = OUT$RhizomeMortality,
    F.ROOTtoDIC                 = OUT$RootTotalRespiration,
    F.ASSIMILATEStoROOTS        = OUT$RootAssimilation,
    F.ROOTtoDEADBELOWGROUND     = OUT$RootMortality,
    F.DEADBELOWGROUNDtoDIC      = OUT$BelowgroundDecay,
    F.DETRITUStoDIC             = OUT$DetritusDecay,
    F.DETRITUStoEXPORT          = OUT$DetritusExport,
    
    F.DETRITUStoBURIAL          = OUT$Burial)
  
  Rates <- data.frame(
    GrossPrimaryProduction    = OUT$GrossPrimaryProduction,
    LeafTotalRespiration      = OUT$LeafTotalRespiration,
    LeafBreaking              = OUT$LeafBreaking,
    LeafMortality             = OUT$LeafMortality,
    LeafAssimilation          = OUT$LeafAssimilation,
    LeafNetReallocation       = OUT$LeafNetReallocation,
    
    DeadLeafAbscission        = OUT$DeadLeafAbscission,
    DeadLeafBreaking          = OUT$DeadLeafBreaking,
    DeadLeafDecay             = OUT$DeadLeafDecay,
    

    RhizomeTotalRespiration   = OUT$RhizomeTotalRespiration,
    RhizomeAssimilation       = OUT$RhizomeAssimilation,
    RhizomeAssimilation       = OUT$RhizomeMobilisation,
    RhizomeMortality          = OUT$RhizomeMortality,
    
    RootTotalRespiration      = OUT$RootTotalRespiration,
    RootAssimilation          = OUT$RootAssimilation,
    RootMortality             = OUT$RootMortality,
    
    BelowgroundDecay          = OUT$BelowgroundDecay,
    DetritusDecay             = OUT$DetritusDecay,
    DetritusExport            = OUT$DetritusExport,
    
    Burial                    = OUT$BURIED)  
  
  # derivatives
    names        <-c("DIC", "ASSIMILATES", "LEAVES", "DEADLEAVES", 
                     "RHIZOMES", "ROOTS", "DEADBELOWGROUND", "DETRITUS",
                     "BURIED")
    # Rows = from, columns = to
    CarbonFlows  <- matrix(nrow = 9, ncol = 9, byrow = TRUE, 
                           data = with(Fluxes, c(
                             
# DIC            ASSIMILATES   LEAVES        DEADLEAVES   RHIZOMES     ROOTS    DEADBELOWGROUND     DETRITUS
  0,              0,           F.DICtoLEAVES,   0,        0,            0,                  0,           0,     0, #DIC
  0,        0,  F.ASSIMILATEStoLEAVES,    0, F.ASSIMILATEStoRHIZOMES, F.ASSIMILATEStoROOTS, 0,           0,     0, #ASSIMILATES
  0,      F.LEAVEStoASSIMILATES, 0, F.LEAVEStoDEADLEAVES, 0,            0,              0,      F.LEAVEStoDETRITUS,0,    #LEAVES
  F.DEADLEAVEStoDIC, 0,        0,               0,        0,            0,              0,     F.DEADLEAVEStoDETRITUS,0, #DEADLEAVES
  F.RHIZOMEtoDIC, F.RHIZOMEStoASSIMILATES, 0,   0,        0,            0,    F.RHIZOMEStoDEADBELOWGROUND, 0,     0,     #RHIZOMES
  F.ROOTtoDIC,     0,          0,               0,        0,            0,    F.ROOTtoDEADBELOWGROUND,     0,     0,     #ROOTS
  F.DEADBELOWGROUNDtoDIC, 0,   0,               0,        0,            0,              0,                 0,     0,     #DEADBELOWGROUND
  F.DETRITUStoDIC, 0,          0,               0,        0,            0,  0,   0,  F.DETRITUStoBURIAL,     #DETRITUS
  0,               0,          0,               0,        0,            0,              0,                 0,     0)
))
    
  rownames(CarbonFlows) <- colnames(CarbonFlows) <- names
  dC <- unlist(dV[c("LEAVES", "ROOTS", "RHIZOMES", "DEADLEAVES", "DEADBELOWGROUND", "BURIED", "DETRITUS")])
  # derivatives
  return(list(Fluxes = Fluxes, Rates = Rates, dC = c(dC, sum = sum(dC)), #Perturb = Perturb,
              Losses = OUT$BURIAL, Delta = colSums(CarbonFlows)["DIC"]- OUT$BURIED, Fluxmat = CarbonFlows))                 
  
}

# changes in state variables over the entire period
dVal <- function(out)  {# takes into account unequal timing 
  toty <- c("time", names(Posidonia_get_states(as.vector = TRUE))) 
  
  OUT <- out[c(1, nrow(out)), toty]
  as.list((OUT[2,]-OUT[1,])/(OUT[2,1]-OUT[1,1]))
}
