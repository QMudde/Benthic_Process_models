
# ==============================================================================
# Mussel parameters, output variables, manipulation functions, etc..
# ==============================================================================

mussel_get_debparms <- function(out, as.vector = FALSE, which = NULL){
  if (missing(out))
    Parms <- .Deb_cohort$parms
  else if ("dtDyn" %in% class(out))
    if (as.vector) 
      Parms <- attr(out, "parms")
    else {
      Parms <- .Deb_cohort$parms
      Parms$value <- attr(out, "parms")
    }
  else stop("object 'out' not supported")
  if (! is.null(which))
    Parms <- subset(Parms, subset = names %in% which)
  return(Parms)
}

mussel_get_debforcs <- function(out, as.vector = FALSE, which = NULL){
  
  Nfrc <- .Deb_cohort$forc
  if (! is.null(which))
    Nfrc <- subset(Nfrc, subset = names %in% which)
  
  if (!missing(out)) {
    Nfrc$mean_value <- colMeans(out[, Nfrc$names])
  }
  Nfrc
}

mussel_get_deb0Dvar <- function(out, as.vector = FALSE, which = NULL, N0D){

  if (! is.null(which))
    N0D <- subset(N0D, subset = names %in% which)
  
  if (!missing(out)) {
    if (attributes(out)$model == "mussel_run_debind")
      N0D <- subset(N0D, subset = names %in% colnames(out))
    if (inherits(out, "dtDyn"))
      N0D$mean_value <- colMeans(out[, N0D$names])
  }
  N0D
}

mussel_get_deb1Dvar <- function(out, as.vector = FALSE, which = NULL, N1D){
  
  if (missing(out)) {
    if (! is.null(which))
      N1D <- subset(N1D, subset = names %in% which)
    return(N1D)
  }
  
  if (inherits(out, "dtDyn")){
    
    if (attributes(out)$model == "mussel_run_debind"){
      N1D <- subset(N1D, subset = names %in% colnames(out))
      if (! is.null(which))
        N1D <- subset(N1D, subset = names %in% which)
      
      N1D$mean_value <- colMeans(out[, N1D$names])
    } else {
      ncoh <- attributes(out)$max_cohort
      findval <- function(name)
        which(colnames(out) %in% paste(name, 1:ncoh, sep="_"))
      N1D$mean_value <- lapply(N1D$names, FUN = function(x)
        mean(colMeans(out[, findval(x)])))
      if (! is.null(which))
        N1D <- subset(N1D, subset = names %in% which)
      
    }

  }
  N1D
}

mussel_get_debvars <- function(out, as.vector = FALSE, which = NULL)
  rbind(mussel_get_deb0Dvar(out       = out, 
                        as.vector = as.vector, 
                        which     = which, 
                        N0D       = .Deb_cohort$out0D),
        mussel_get_deb1Dvar(out       = out, 
                        as.vector = as.vector, 
                        which     = which, 
                        N1D       = .Deb_cohort$out1D)
        )
  
mussel_get_debstates <- function(out, as.vector = FALSE, which = NULL){
  rbind(mussel_get_deb1Dvar(out       = out, 
                        as.vector = as.vector, 
                        which     = which, 
                        N1D       = .Deb_cohort$y1D),
        mussel_get_deb0Dvar(out       = out, 
                            as.vector = as.vector, 
                            which     = which, 
                            N0D       = .Deb_cohort$y0D)
  )
}

# ------------------------------------------------------------------
# Function to extract population structure from population model
# ------------------------------------------------------------------

mussel_pop_struct = function(model, age_classes = c(0,1,2,3), length_classes=c(0,1,2.5,5)){
  # the vectors age_classes & length_classes must be values of cutoff points of the requested sturcture
  
  stopifnot(
    "use max 5 cutoff points for structure" = (length(age_classes) < 6 & length(length_classes) < 6),
    "use min 2 cutoff points for structure" = (length(age_classes) > 1 & length(length_classes) > 1),
    "Both age and Length classes must start at 0, this is by definition the lower cutoff" = 
      (age_classes[1]==0 & length_classes[1]==0),
    class(model) == c("dtDyn","deSolve","matrix"),
    "THIS FUNCTION IS NOT UPDATED" = ((ncol(model)-26) /12) == 25
  )
  
  # define age classes based on 1 or 2 spawning seasons per year 
  df = as.data.frame(model)
  time = df$time
  doy = time%%365 + 1
  spawns = attributes(model)$spawn_seasons
  coh = attributes(model)$max_cohort
  doy_coh2 = attributes(model)$parms[["spawn_peak_2"]] - 3*attributes(model)$parms[["spawn_dev"]] #doy from which 2nd coh is defined
  yr = floor(model[,1]/365) # vector of year, with length of time
  
  new = yr%%coh + 1 # vector of cohort indexes that is newly spawned that (half) year, with length of time
  
  if(spawns == 2){
    new =   (yr%%(coh/2))*2+1+(doy>doy_coh2)
  }
  # calculate age in years rounded down for each cohort; a<-b (A becomes b); assign (a,b)
  # create empty age class vectors
  age_struct = data.frame(time = time)
  #----------------------------#
  ### EXTRACT AGE CLASS DATA ###
  #----------------------------#
  for(CLASS in 1:length(age_classes)){
    
    from <- age_classes[CLASS]
    if(CLASS<length(age_classes)){
      until <- age_classes[(CLASS+1)]
    }
    else{until <- coh} #hope this dont crash lol
    class_name <- paste0("POP_age_",from,"_",until)
    
    # assign(class_name,rep(0,nrow(model)))
    
    class <- rep(0,nrow(model))
    for(i in 1:coh){
      # index population of i-th cohort
      ii = 1+(coh*3)+i
      # different age calculation based on yearly spawning seasons
      if(spawns == 1){
        # make vector from 1 - times of ages
        ages = (new-i)%%coh # vector of ages of cohort i 
      } 
      
      if(spawns == 2){
        ages = floor( ((new-i)%%coh)/2 )
      }    
      target_age_pop <- model[,ii] *( (ages >= from) & (ages < until)) # "until" is not included in ageclass
      class <- class + target_age_pop
    }
    
    age_struct[class_name]<- class
    
  }
  
  #-------------------------------#
  ### EXTRACT LENGTH CLASS DATA ###
  #-------------------------------#
  
  for(CLASS in 1:length(length_classes)){
    
    from <- length_classes[CLASS]
    if(CLASS<length(length_classes)){
      until <- length_classes[(CLASS+1)]
    }
    else{until <- 15} #hope this dont crash lol
    class_name <- paste0("POP_length_",from,"_",until)
    class <- rep(0,nrow(model))
    for(i in 1:coh){
      # index population of i-th cohort
      ii = 1+(coh*3)+i
      # index Length (cm) of the i-th cohort
      IL = 1+(coh*4)+6+(coh*2)+i
      Lengths = model[,IL]
      target_length_pop <- model[,ii] *( (Lengths >= from) & (Lengths < until))
      class <- class + target_length_pop
    }
    
    age_struct[class_name]<- class
    
  }
  return(age_struct)
  
}


# -----------------------------------------------------------------------
# Function to Determine initial DEB conditions based on Length & Biomass 
# -----------------------------------------------------------------------

get_INI_indiv <- function(Length , AFDW = NA, DW = NA, WW = NA,
                          Rel_E_dens = 1, matured = TRUE, parms = list()){
  # get default parameters
  P<- getParms_deb(parms)
  aux<- .Deb_cohort$auxiliary
  # initial STRUCT CM = LENGTH*del.m = 4.03 * 0.33 = 1.3299
  #initial assumption: shellfish is a maximum engergy density, matured but no gametes
  ini.struct.cm = Length * P[["del.M"]]
  # struct CM to STRUCT C -> (STRUCT CM^3) * mmolC_cm3 = ini.struct.cm*3 
  ini.struct.C = (ini.struct.cm^3)*P[["mmolC_cm3"]]
  # can base R.S ratio on DEB parameters: 
  maxRS = 2190/1771 # [DEB max enegery dens [j/cm3] / J/cm3] ==> max R / S
  R.S.RATIO = Rel_E_dens*maxRS
  # in case no weight are given,
  # we can provide standard starting ratio of struct:reserve:reprod = 1:1.5:0.5 with full maturity buffer
  # these ratios are based on Karayusel ini.values for witch weight is known
  # full maturity is 3.022857 mmolC in this model
  
  # determine AFDW based on weight estimates
  afdw = AFDW 
  if(!is.na(DW)) { afdw = DW * aux[["gAFDW_gSFDW"]]}  # new formulation! # 13.89 # 15.38462
  if(!is.na(WW)) { afdw = aux[["gAFDW_gSFDW"]] * (WW /P[["gWW_gDW"]])}  # new form# 13.89 # 15.38462
  # if no weight is given, make rough estimation 
  # estimation is that struct and reprod are distributed with kappa ratio 
  # and based on maturity level, assume standard R.S ratio based on DEB pars
  
  if(is.na(afdw)) {
    #    if(matured){
    #      ini.vals = c(RESERVE = ini.struct.C * R.S.RATIO,
    #                   STRUCT = ini.struct.C,
    #                   REPROD = (pars[["mat_pub"]]))} #  
    #    else{
    #      ini.vals = c(RESERVE = ini.struct.C * R.S.RATIO,
    #                   STRUCT = ini.struct.C,
    #                   REPROD = (ini.struct.C))
    
    ini.vals = c(RESERVE = ini.struct.C * R.S.RATIO,
                 STRUCT = ini.struct.C,
                 # assume reprod ratio based on kappa carbon distribution
                 REPROD = (ini.struct.C*( 1/(P[["kappa"]]/(1-P[["kappa"]])) )))
    
  }
  else{
    # calculate rest from afwd with : AFDW.g = Biomass.C / 1000 * 12 / C_frac
    # Biomass.C = AFDW.g*C_frac/12*1000      (C frac AFDW = 0.5)
    biomass.C = afdw*41.66667  #### 0.5/12*1000
    cat("mmolC(ind) = ", biomass.C, "\n")
    # still assume struct: reserve = R.S.RATIO:1
    ini.reserve.C = ini.struct.C * R.S.RATIO
    rest.C = biomass.C - (ini.struct.C+ini.reserve.C)
    if(rest.C<0){print("ERROR, Total weight is too low for standard ratios")}
    ini.tot.repro.C <-  rest.C  #+pars[["mat_pub"]] 
    ini.vals = c(RESERVE = ini.reserve.C,
                 STRUCT =ini.struct.C,
                 REPROD = ini.tot.repro.C)
  }
  return(ini.vals)
}

## get initial conditions for a population (cohort) based on initial density
get_INI_coh <- function(Indiv_ini, Dens_ini,max_coh = 12,
                        PhyC_ini=10, DT_ini=10, 
                        Sim_ini=1, O2_ini=360){
  stopifnot(
    "max cohorts must be an integer between 2 and 20" = 
      (max_coh %in% c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))
  )
  cohNR <- max_coh-1
  
  INI_vector_coh <- c(RESERVE_ = c(Indiv_ini[[1]],rep(0,times=cohNR)),
                      
                      #RESERVE = Indiv_ini[[1]],RESERVE= rep(0,times=12), 
                      STRUCT_ = c(Indiv_ini[[2]],rep(0,times=cohNR)),
                      REPROD_ = c(Indiv_ini[[3]],rep(0,times=cohNR)),
                      POPULATION_ = c(Dens_ini,rep(0,times=cohNR)),
                      LARV_DENS=0, LARV_BIOM=0,
                      PHYTO =PhyC_ini, DETRITUS=DT_ini,Sim = Sim_ini, O2=O2_ini)
  return(INI_vector_coh) 
}

#-------------------------------------------------------------------------------
# Function to construct forcing function from time series of different unit data
#-------------------------------------------------------------------------------

make_ff <- function(day,values, recalc = "none"){
  # type = phyc, dt, O2,sim, temp
  CHL_to_PHY <- .Deb_cohort$auxiliary[["CHL_to_PHY"]]
  PHYTO_g_to_mmolC <- .Deb_cohort$auxiliary[["PHYTO_g_to_mmolC"]]
  POM_g_to_mmolC <- .Deb_cohort$auxiliary[["POM_g_to_mmolC"]]
  
  
  if(!recalc %in% c("phyc_g","pom_g","chl_mugl","none"))
    stop("type must be:'phyc_g' 'pom_g','chl_mugl','none, IF NO RECALULATION TYPE = 'none")
  
  if(recalc =="chl_mugl"){
    cat("unit recalucated from mug/l chl to mmolC phyc/m3")
    FF = approxfun(x =day,y =values*CHL_to_PHY, rule = 2)
  }
  if(recalc =="phyc_g"){
    cat("unit recalucated from gram phyc/m3 chl to mmolC phyc/m3")
    FF = approxfun(x =day,y =values*PHYTO_g_to_mmolC, rule = 2)
  }
  if(recalc =="pom_g"){
    cat("unit recalucated from gram dt/m3 to mmolC dt/m3")
    FF = approxfun(x =day,y =values*POM_g_to_mmolC, rule = 2)
  }
  if(recalc == "none"){
    cat("type=none, units are assumed corrent, no changes")
    FF = approxfun(x =day,y =values, rule = 2)
  }
  
  
  forcing_out =  cbind(day,   FF(day))
  return(forcing_out)
}

# Get auxiliry parameters

mussel_get_auxiliary <- function(out, as.vector = FALSE, which = NULL){
  
  auxP <- .Deb_cohort$auxiliary
  if (! is.null(which))
    auxP <- subset(auxP, subset = names %in% which)
  
  if (!missing(out)) {
    auxP$mean_value <- colMeans(out[, auxP$names])
  }
  auxP
}
