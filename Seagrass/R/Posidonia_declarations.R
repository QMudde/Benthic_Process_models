.Posidonia <- new.env()

# ================================================
# Parameters
# ================================================

Pparms <- c(
    
   #Physical parameters
    salinity       = 38,      # -water salinity        
    
    latitude       = 42.5833, # latitude to calculate daylength
    
    depth          = 10,     # m   mean depth of the bay
    
    C_DWT          = 0.4,    # gC/gDW required to estimate dry weight of Posidonia  
    OCr            = 1.,     # MolO2/MolC in respiration and photosynthesis
    NCr            = 0.1,    # MolN/MolC N:C ratio in (all) tissues
    
    arrheniusCt    = 1.07,   # Temperature Dependency of Process rates (repiration,..)
    
    #Photosynthetic rate determined by temperature Savva et al. 2018
    POPT           = 0.01/79.8848 , # scaling factor for T-dependence of Photosynthesis",
    TOPT           = 25.8 ,         # optimal T for Photosynthesis and N uptake",

    TMAX           = 33.8 ,    # maximal temperature above which PS and N uptake = 0
    TMIN           = 10.0 ,    # minimal temperature below which PS and N uptake = 0
    
    # Photosynthesis
    maxPhotosynthesis = 0.01,     # /day   ( 0.15 ORG)
    ksLight           = 50,       # [W/m2]  half-saturation coefficient for light
    
    k0short           = 0.05,     # m-1     Attenuation coefficient of short wavelengths
    k0long            = 0.10,     # m-1     Attenuation coefficient of long wavelengths
    plong             = 0.63,     # -       Part of PAR w/ a long wavelength
    
    ksE               = 0.0008,   # half saturation ct of epifytes covering the leaves (DIN conc)
    maxEshading       = 0.15,     # max impact Epifytes can have on carbon gain. 
    
    
    L_surf_Ratio      = 500,  # surface:C ratio, for light competition [cm2/molC] Olesen et al 2002
    max_selfshade     = 0.50, # assume at least 50% of leaves are never self shaded TBC
    ks_shade          = 5000, # [cm2] total leaf surface/m2 at which self shading is half TBC
    
    # Posidonia N uptake from soil via roots
    maxRootsNupt   = 0.00045, # [molN/molC/m2/day] max nutrient absorption from soil by roots
    ksRootsNupt    = 0.05,    # [molN] half effect value of sediment Nutrient uptake by roots
    maxRhizomeNC   = 0.3,     #  [molN/molC] max N reserves relative to rhizome biomass
    minRhizomeNC   = 0.05,    #  [molN/molC] min N reserves relative to rhizome biomass
    ksROOTS        = 45,      # [mmolC] inhibition coefficient for root uptake ~roots
    
    # Remobilization rates
    rhizomeMobilisationRate = 0.00035, #  /day     Rhizome remobilisation rate
    
    # Respiration rate
    rhizomeRespirationRate = 0.000500, #  /D    Basal respiration rate for Rhizome
    rootRespirationRate    = 0.000500, #  /D    Basal respiration rate for Roots            
    leafRespirationRate    = 0.0010,   #  /D    Basal respiration rate for leaves
    respCost               = 0.2,      #  -     Fraction of assimilation that is lost to respiration
    
    # Mortality rate and Reallocation, Abscission rate, N recycling
    belowMortalityRate     = 0.00020, # /D      Mortality rate of belowground organs    0.00045            
    leafMortalityRate      = 0.0045,  # /D     Leaves mortality rate 0.0085   0.0055
    partLeafReallocation   = 0.40,    #  -    fraction of senescence reallocated to leaves Lepoint et al. 2002
    ReallRespcost          = 0.50,    # cost of C to reallocate C and N from leaves to rhizomes  
    leafAbscissionRate     = 0.005,   # /D    Dead leaves abscission rate
    NecrosisRate           = 0.00,    # /D max necrotic mortality by temperature  0.00337
    
    #Decay rates
    belowgroundDecayRate   = 0.007,   # /D  Decay rate of dead belowground organs                   
    abovegroundDecayRate   = 0.01,    # /D  Decay rate of dead leaves                               
    detritusDecayRate      = 0.01,    # /D  Decay rate of detritus
    pBurial                = 0.3,     # [-] fraction of dead belowground C production buried
    
    #Breaking rates of leaves by wind storms and export of detritus
    leafBreakingRate     = 0.015,   # /D 
    deadleafBreakingRate = 0.0025,  # /D 
    expWind                = 2.,      # 
    ksW                    = 8.,      # m/s
    detritusExportRate     = 0.03,    #/D   during windy conditions
    mindetritusExportRate  = 0.0001,  #/D   during calm conditions  
    
    #Partitioning of growth
    shapebelow             = 0.32,    # 
    pBelowToRoot           = 0.4,     # -    fraction of belowground allocates to roots    
    pBelowGrowthMin        = 0.30,    # -    min part of assimilates allocated to belowground organs 
    
    
    #Phenological events
    TimingEndMobilRhizomes        = 270, # D   day at which mobilisation from rhizomes stops
    TimingRestartMobilRhizomes    = 366, # D   day at which mobilisation from rhizomes restarts
    TimingMaxStorageBelowground   = 365, # D   day at which maximum storage occurs
    senescenceDaylength           = 12,  # Hr  daylength at which senescence starts and stops
    
    kSenescence             = 1,  # factor for senescence scaling relative to daylength
    refDaylength            = 13, # reference daylength for senescence scaling (hours
    
    # rates that determine Nutrient concentration in water & sediment
    NexchangeRate  =  0.01    #[/d] exchange rate between water and sediment
)

.Posidonia$parms <- data.frame(names   = names(Pparms), 
                               default = Pparms)
row.names(.Posidonia$parms) <- NULL

# Description of parameters
.Posidonia$parms$description <- c(
  "salinity of water", 
  "latitude, used to estimate daylength",
  "mean depth of the bay",
  "composition required to estimate dry weight of Posidonia",
  "MolO2/MolC in respiration and photosynthesis",
  "MolN/MolC ratio in all tissues",

  "temperature Dependency of Process rates", 
  "scaling factor for temperature-dependence of Photosynthesis",
  "optimal temperature for Photosynthesis and N uptake",
  "maximal temperature above which Photosynthesis and N uptake = 0",
  "minimal temperature below which Photosynthesis and N uptake = 0",
  
  "maximum photosynthesis rate", 
  "half-saturation coefficient for light",
  
  "Attenuation coefficient of the short wavelengths", 
  "Attenuation coefficient of the long wavelengths",
  "Part of PAR with a long wavelength",
  "half saturation DIN conc, related to epifytes covering the leaves",
  "max negative impact Epifytes can have on carbon gain",
  
  "surface:C ratio, for light competition",
  "maximal fraction of leaves that are shaded",
  "total leaf surface/m2 at which self shading is half",
  
  "max nutrient absorption from soil by roots", 
  "half effect value of sediment Nutrient uptake by roots",
  "max N reserves relative to rhizome biomass",
  "min N reserves relative to rhizome biomass",
  "inhibition coefficient for root uptake by roots",
  
  "rate at which rhizome carbon is mobilised", 
 
  "basal rhizome respiration rate", 
  "basal root respiration rate", 
  "basal leaf respiration rate", 
  "mole C respired to convert of one mole glucose to biomass",

  "mortality rate of belowground tissue", 
  "mortality rate of leaves", 
  "fraction of senescence reallocated to leaves",
  "cost of C to reallocate C and N from leaves to rhizomes",
  
  "abscision rate of dead leaves",
  "max necrotic mortality by temperature", 
  
  "decay rate of belowground tissue", 
  "decay rate of above ground tissue", 
  "decay rate of detritus", 
  "fraction of produced dead belowground carbon that is sequestered",
 
  "breaking rates of leaves by wind storms",
  "breaking rates of dead leaves by wind storms", 
  "exponent of wind in breaking", 
  "half-saturation in wind speed in breaking", 
  "detritus export rate due to wind (maximum)",
  "detritus export rate in calm conditions (minimum)",
  
  "shape parameter for growth partitioning", 
  "fraction of belowground allocates to roots",
  "minimal part of assimilates allocated to belowground organs",
 
   "day at which mobilisation from rhizomes stops", 
  "day at which mobilisation from rhizomes restarts",
  "day at which storage of carbon to rhizomes and root occurs",
  "daylength (hours) at which senescence starts and stops",
  
   "factor for senescence scaling relative to daylength",
   "reference daylength for senescence scaling (hours)",
  
  "exchange rate from dissolved nutrients between water and sediment"
)
  
# Units of parameters

.Posidonia$parms$units <- c(
  salinity = "-",  latitude = "degN", depth = "m", C_DWT = "GC/GDW", 
  OCr = "molO2/molC", NCr = "molN/molC", 
  arrheniusCt =  "-", POPT = "-", TOPT = "degC", TMAX = "degC",  TMIN = "degC",  
  maxPhotosynthesis = "/day", 
  ksLight = "W/m2",  k0short = "/m", k0long = "/m", plong = "-", 
  ksE = "mol/m3", maxEshading = "-", L_surf_Ratio = "cm2/molC", 
  max_selfshade = "-", ksShade = "cm2/m2", 
  maxRootsNupt = "molN/molC/day", ksRootsNupt = "molN/m2", 
  maxRhizomeNC = "molN/molC",  minRhizomeNC = "molN/molC", ksROOTS = "molC/m2",
  rhizomeMobilisationRate = "/day", rhizomeRespirationRate = "/day", 
  rootRespirationRate = "/day", leafRespirationRate = "/day", respCost = "-", 
  belowMortalityRate = "/day", leafMortalityRate = "/day", 
  partLeafReallocation = "-", ReallRespcost = "molO2/molC", 
  leafAbscissionRate = "/day", NecrosisRate = "/day",
  belowgroundDecayRate = "/day", abovegroundDecayRate = "/day", 
  detritusDecayRate = "/day", pBurial = "-",
  leafBreakingRate =  "/day", deadleafBreakingRate = "/day", 
  expWind = "-",ksW = "m/s", 
  detritusExportRate ="/day", mindetritusExportRate = "/day", 
  shapebelow = "-", pBelowToRoot = "-", pBelowGrowthMin = "-", 
  TimingEndMobilRhizomes = "day", TimingRestartMobilRhizomes = "day", 
  TimingMaxStorageBelowground = "day",    
  senescenceDaylength = "hour", kSenescence = "-", refDaylength = "hr", 
  NexchangeRate = "m/day"
)

# ================================================
# State variables
# ================================================

States <- c(LEAVES          = 12.,      # molC/m2   initial biomass of leaves  
            ROOTS           = 32,       # molC/m2   initial biomass of roots  
            RHIZOMES        = 42,       # molC/m2   initial biomass of rhizomes 
            NRESERVES       = 11,       # molN/m2   initial Nutrient Reserve of Posidonia
            DEADLEAVES      = 5.,       # molC/m2   initial biomass of dead aboveground biomass
            DEADBELOWGROUND = 3.,       # molC/m2   initial biomass of dead belowground biomass
            BURIED          = 0.,
            DETRITUS        = 1.,       # molC/m2   initial detritus concentration
            OXYGEN          = 0.245,    # molO2/m3  initial oxygen concentration
            DIN             = 0.003,    # molN/m3   initial DIN concentration 
            DINSED          = 0.03,     # molN/m2   initial DIN in sediment 
            CLOSING         = 0.0       # molC/m3    initial value for the closing variable
            #            CLOSING         = 0.0     # molC/m3   initial value for the closing variable 
)
.Posidonia$y <- data.frame(names = 
    names(States), default = States)

.Posidonia$y$description <- 
    c("C biomass of leaves ", "C biomass of roots", 
      "C biomass of rhizomes", "Nutrient reserves",
      "dead aboveground C biomass",
      " dead belowground C biomass", "total C buried",
      "detritus biomass", 
      "pelagic oxygen concentration", "pelagic dissolved inorganic nitrogen concentration",
      "benthic dissolved inorganic nitrogen concentration", 
      "variable to close the C budget")
.Posidonia$y$units <-
    c(rep("molC/m2", times = 3), "molN/m2", rep("molC/m2", times = 4), 
      "molO2/m3", "molN/m3", "molN/m2", "molC/m2")

# ================================================
# Output variables
# ================================================

.Posidonia$out <- data.frame(names =
  
  c(      "LeavesDWT", "BelowDWT", "GrossPrimaryProduction",  
          "AssimilateProduction", "TotalRespiration", "NetO2production",  
          "BasalRespiration", "GrowthRespiration", "TotalDecay", 
          
          "LeafPrimaryProduction", 
          "LeafAssimilation",  "LeafGrowth", "LeafTotalRespiration",   
          "LeafMortality", "LeafNetReallocation", 
          "LeafReallRespiration", "LeafBreaking",   

          "RootAssimilation" , "RootGrowth", "RootTotalRespiration", 
          "RootMortality",
          
          "RhizomeAssimilation",  "RhizomeGrowth", "RhizomeTotalRespiration",            
          "RhizomeMortality",  "RhizomeMobilisation",  
          
          "DeadLeafAbscission",  "DeadLeafBreaking",  "DeadLeafDecay",
          
          "DetritusExport", "DetritusDecay",  "BelowgroundDecay", 
          "MortToDeadBelowGround", "Burial",  
          
          "RootsNupt", "NreserveToGrowth", "RespirationNRecycling",
          "LeafNRecycling", "RhizomeNRecycling", 
          "NFluxtoBottom", 
          
          "AirSeaExchange", "Piston", "Satox",
          
          "WindBreaking",  "pBelowGrowth", 
          
          "TotalC", "TotalN", "Senescence", "Remobilization", 
          
          "LossTempfac", "GrowthTempfac",
          
          "Self_shading",  "Leaf_surf", 
          
          "PSLightLim", "PSepiphyteLim" , "PSNCquotaLim",
          "PSlim", "RupNLim", "RupRootLim", "RupNCquotaLim",
          
          "Temperature", "Wind", "Light", "PAR", "Daylength" ))


.Posidonia$out$description <- 
  c("leaves dry weight", "below ground dry weight", 
    "gross primary production", "total assimilate production",
    "total respiration", "net oxygen production", 
    "basal (maintenance) respiration", "growth respiration", 
    "total decay", 
    
    "leaf (gross) photosynthesis", "leaf assimilation", 
    "leaf growth", "total leaf respiration", 
    "leaf mortality", "leaf net reallocation", 
    "leaf respiration for reallocation",  "leaf breaking (wind)",

    "root assimilation", "root growth", "total root respiration", 
    "root mortality", 
    
    "rhizome assimilation", "rhizome growth", "total rhizome respiration", 
    "rhizome mortality", "rhizome remobilisation",

    "dead leaves abscission", "dead leaves breaking", "dead leaves decay",
    
    "detritus export", "detritus decay", "total belowground decay",
    "part of belowground mortality to dead belowground",
    "burial rate (part of belowground mortality)",
    "nitrogen uptake by the roots", 
    "nitrogen from N reserve used for growth",
    "nitrogen from basal respiration stored in N reserve",
    "nitrogen from leaf reallocation stored in N reserve",
    "nitrogen from rhizome mobilization stored in N reserve",
    
    "net nitrogen flux from the water column to the bottom", 
    "air sea oxygen exchange", "piston velocity", "saturated oxygen concentration",
    
    "wind breaking", 
    "part of assimilates allocated to belowground organs",
    
    "total carbon in system", "total nitrogen in system", 
    "factor at which senescence (die-off) occurs",
    "whether Posidonia is in remobilisation phase or not", 
    
    "temperature scaling factor for loss processes", 
    "temperature scaling factor for growth processes", 
    "self-shading fraction", "total leaf surface for self-shading",
    "light limitation for photosynthesis", 
    "photosynthesis reduction due to epiphytes", 
    "photosynthesis dependence on Nstorage in rhizomes",
    "total limitation and temperature effect on photosynthesis",
    "nitrogen limitation for root N-uptake",
    "root limitation for root N-uptake",
    "root N-uptake dependence on Nstorage in rhizomes",
    
    "water temperature",
    "wind speed", "light intensity at air-water interface",
    "Photosynthetically active radiation on leaves", "length of the day")

.Posidonia$out$units <-
  c(     "gDWT/m2",  "gDWT/m2", "molC/m2/d", "molC/m2/d",
         "molC/m2/d", "molO2/m2/d", "molO2/m2/d", 
         rep("molC/m2/d", times = 27),
         
         "molN/m2/d", "molN/m2/d", "molN/m2/d",
         "molN/m2/d", "molN/m2/d", "molN/m2/d", 
         "molO2/m2/d", "m/d", "mol O2/m3",
        "-",  "-","molC/m2", "molN/m2",  "-",  
          "-", "-", "-", "-", "cm2",
         "-", "-", "-", "-", "-", "-", "-",
         "degC", "m/s", "uEinst/m2/s", "uEinst/m2/s", "hour")

.Posidonia$forc <- data.frame(names = 
   c("f_Temperature", "f_Light", "f_Wind", "f_Piston", "f_SatO2"))

.Posidonia$forc$description <- c(
   "environmental temperature", "light at the air-water interface", 
   "wind speed", "piston velocity", "oxygen saturated concentration")
.Posidonia$forc$units <- c(
   "dg C", "W/m2", "m/s", "m/d", "molO2/m3")

  
