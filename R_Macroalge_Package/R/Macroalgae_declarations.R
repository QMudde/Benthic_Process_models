.Macroalgae <- new.env()

# ================================================
# Parameters
# ================================================

Mparms <- c(

  # ------------------------------------------------ # 
  # Quota and ratios:                                #
  # ------------------------------------------------ # 
  
  pC_SDW   = 0.3,  # [gSC/gSDW] mass fraction C of Structural dry weight
  RC_SCmin = 0.05, # [molRC/molSC] minimum molar C reserve over structural C ratio

  # N quota and ratios
  SN_SCr   = 0.0355,    # [molSN/molSC] molar N:C ratio of structural tissue
  RN_SCmin = 0.0355,    # [molRN/molSC] minimal reserve N:structural C ratio  
  RN_SCmax = 0.1285625, # [molRN/molSC] maximal reserve N:structural C ratio  

  # P quota and ratios
  SP_SCr = 0.0011,    # [molSP/molSC] molar P:C ratio of structural tissue
  RP_SCmin = 0.0011,  # [molRP/molSC] minimal reserve P:structural C ratio 
  RP_SCmax = 0.0039,  # [molRP/molSC] maximal reserve P:structural C ratio  
  
  # reserve specific weight ratios 
  # reserve C galactose chains, C6H10O5; reserve N: 50%NO3, 50%NH4; reserveP: PO4
  gdw_RC = 2.25*0.012,     # [gDW/mmolRC] gram dry weight per mmol reserve carbon
  gdw_RN = 2.821429*0.014, # [gDW/mmolRN] gram dry weight per mmol reserve nitrogen  
  gdw_RP = 2.548387*0.031, # [gDW/mmolRP] gram dry weight per mmol reserve Phosphorus 
  gdw_SC = 3.894693*0.012, # [gDW/mmolSC] gram dry weight per mmol structural carbon
  
  OC_q   = 1.105,     # [mol O2/mol RC] O2:C ratio in respiration

  dw_ww  = 0.195,     # [gDW/gWW] dry:wet weight ratio
  
  
  # ------------------------------------------------ # 
  # Physiological rates: Photosynthesis & Production #
  # ------------------------------------------------ #
  
  # Photosynthesis rate = PS*(1-exp(-alpha*I/PS))*exp(-beta*I/PS) (Platt 1980)
  # PS = Pmax / ( (alpha/(alpha+beta))* (beta/(alpha+beta)^(alpha/beta)))

  pMax  = 0.258793,    # [molC/molSC/day] Gross Photosynthetic rate at reference Temperature.
  alpha = 0.002623161, # [molC/molSC/(uMol photons/m2)] Initial PI curve slope parameter
  beta  = 0.480399662, # [molC/molSC/(uMol photons/m2)] Photo inhibition parameter
  
  # dimensionless relative affinity per wavelength range to general PAR 
  # (mean must be 1)
  rho_short = 0.57, #  [-] rel affinity for short wave (blue), 400-500nm 
  rho_mid   = 1.8,  #  [-] rel affinity for midwave (green), 500-600nm
  rho_long  = 0.73, #  [-] rel affinity for longwave (red), 600-700nm
  
  # Temperature effect on photosynthesis & growth - arrhenius parameters
  t_Arh  = 2667.257,  # [dgK] Arrhenius factor  [-449.7342]
  t_A_L  = 61192.324, # [dgK] Arrhenius T factor for lower bound [20640.9346]
  t_A_H  = 17848.025, # [dgK] Arrhenius T factor for lower bound [43549.6620]
  t_LB   = 274.15,   # [dgK] lower bound of optimal temperature range [275.9967]
  t_UB   = 301.15,   # [dgK] upper bound of optimal temperature range [306.0431]
  t_ref  = 288.15,     # [dgK] reference temperature for photosynthesis
  
  ExhFmax = 0.638130226,       # [molC/molC] maximal fraction photosynthesis exudated  
  #-CALIBRATED- [0.5]
  
  # self shading 
  
  h_ALG   = 0.5,       # [m](maximal) height of Phyllophora thallus
  a_ALG = 0.00030672,  # [m2/mmolC] (total) Carbon specific shading coefficient
  
  # Impact of competing algae on light
  ww_ALG2 = 0,        # [gWW]    wet weight of competing algae
  h_ALG2  = 1,        # [m]      height of competing algae
  a_ALG2  = 0.002556, # [m2/gWW] wet weight specific shading coefficient
  
  # ----------------------------------------- # 
  # Physiological rates: Growth & Respiration #
  # ----------------------------------------- #
  
  #      = 0.0276,    # [/d] maximal specific growth of structural biomass  
  gRmax = 0.01, # [/d] maximal specific growth of structural biomass corrected
  # -CALIBRATED- 0.0282

  basalResp  = 0.01655588, # [molC/molSC/d] dark repsiration rate 
  growthResp = 0.7,        # [-] fraction of C respired for structural growth 

  # @ QUINTEN: use normal Q10 formulation instead (tRChange = log(Q10)/10)
  tRChange = 0.152748616, # speed of resp change with temperature -CALIBRATED- [0.05]
  tRRef    = 14,   # reference temperature at which dark respiration is measured (Johansen&Snoeijs 2002)

  # ------------------------------------ # 
  # Physiological rates: Nutrient uptake #
  # ------------------------------------ #
  # Kinetics for 3 different nutrient types: nitrate [NO3], ammonium [NH4] & phosphate [PO4]
  max_NO3_upt = 0.00189883, # [molN/molSC/d] max NO3 uptake rate
  ksNO3       = 9.214286,   # [mmol NO3/m3] Monod half sat constant for NO3 upt 

  max_NH4_upt = 0.01089624, # [molN/molSC/d] max NH4 uptake rate
  ksNH4       = 7.928571,   # [mmol NH4/m3] Monod half sat constant for NH4 upt 

  max_PO4_upt = 0.0001302587, #[molP/molSC/d] max PO4 uptake rate
  ksPO4       = 0.3677419,    #[mmol PO4/m3] Monod half sat constant for PO4 upt 

  # other parameters affecting Nutrient uptake
  u065        = 0.03,        # [m/s]  currentspeed where nutrient uptake=0.65*rNup  

  # --------------------------------- # 
  # External rates: Mortality factors #
  # --------------------------------- #
  
  mort_rate   = 0.0042,   # [/d] specific mortality rate -CALIBRATED- [0.008]
  nec_rate    = 0,       # [/d] heat related necrosis rate
  ero_max     = 1,       # [/d] max erosion rate per structural biomass 
  ero_inc     = 0.05, # [/(mmolC/m2)] increase of erosion as surface area increases
  # -CALIBRATED-  0.007333333 (0.022/3)

  #==========================#
  # ENVIRONMENTAL PARAMETERS #
  #==========================#
  
  # --------------------------------- # 
  # Light attenuation coefficients    # 
  # --------------------------------- # 
  
  # @ quinten: waarom niet koclear*rk_short en koclear*rk_mid ?
  # Q: done
  k0_short = 0.3*0.0384,   # [/m] light extinction coeff of 400-500nm [clear water]
  k0_mid   = 0.45*0.0384,  # [/m] light extinction coeff of 500-600nm [clear water]
  k0_long  = 2.25*0.0384,  # [/m] light extinction coeff of 600-700nm [clear water]
  
  k0_turb  = 0.1,   # [/m] PAR range attenuation for turbid water
  
  k_Phyto_short = 0.0249,      # [/(mg Chl/ m3)] extra light extinction of Chl on short waves
  k_Phyto_mid   = 0.009015,    # [/(mg Chl/ m3)] extra light extinction of Chl on mid waves
  k_Phyto_long  = 0.008642857, # [/(mg Chl/ m3)] extra light extinction of Chl on long waves

  # additional information about environment
  latitude = 45,   # [dgN] determines daylight 
  depth    = 10,   # [m] determines light that reaches bottom
  dilution = 0.25  # [/d] dilution of DIN with environemtn/day (0.25 used in Hadley et al., 2015)
)

.Macroalgae$parms <- data.frame(names   = names(Mparms), 
                                default = Mparms)
row.names(.Macroalgae$parms) <- NULL

# Description of parameters
.Macroalgae$parms$description <- c(
    "mass fraction C of Structural dry weight",
    "minimum molar reserve over structural C ratio",
    "molar N:C ratio of structural tissue",
    "minimal reserve N:structural C ratio",
    "maximal reserve N:structural C ratio",
    "molar P:C ratio of structural tissue",
    "minimal reserve P:structural C ratio",
    "maximal reserve P:structural C ratio",
    "gram dry weight per mmol reserve carbon",
    "gram dry weight per mmol reserve nitrogen",  
    "gram dry weight per mmol reserve Phosphorus", 
    "gram dry weight per mmol structural carbon",
    "O2:C ratio in respiration",
    "dry:wet weight ratio",
    "Gross Photosynthetic rate at reference Temperature.",
    "Initial PI curve slope parameter",
    "Photo inhibition parameter",
    "rel affinity for short wave (blue), 400-500nm",
    "rel affinity for midwave (green), 500-600nm",
    "rel affinity for longwave (red), 600-700nm",
    "Arrhenius factor",
    "Arrhenius T factor for lower bound",
    "Arrhenius T factor for lower bound",
    "lower bound of optimal temperature range",
    "upper bound of optimal temperature range",
    "reference temperature for photosynthesis",
    "maximal fraction photosynthesis exudated",
    "height of Phyllophora thallus",
    "(total) Carbon specific shading coefficient",
    "wet weight of competing algae",
    "height of competing algae",
    "wet weight specific shading coefficient",
    "maximal specific growth of structural biomass",
    "dark respiration rate",
    "fraction of C respired for structural growth",
    "speed of resp change with temperature",
    "reference temperature at which dark respiration is measured",
    "max NO3 uptake rate",
    "Monod half sat constant for NO3 uptake",
    "max NH4 uptake rate",
    "Monod half sat constant for NH4 uptake",
    "max PO4 uptake rate",
    "Monod half sat constant for PO4 uptake",
    "currentspeed where nutrient uptake=0.65*rNup",  
    "specific mortality rate",
    "heat related necrosis rate",
    "max erosion rate per structural biomass",
    "increase of erosion as surface area increases",
    "light extinction coeff of 400-500nm [clear water]",
    "light extinction coeff of 500-600nm [clear water]",
    "light extinction coeff of 600-700nm [clear water]",
    "PAR range extinction coefficient for turbid water",
    "Chl-specific light extinction on short waves",
    "Chl-specific light extinction on mid waves",
    "Chl-specific light extinction on long waves",
    "latitude, determines daylight", 
    "water depth, determines light that reaches bottom",
    "dilution rate of nutrients with environment")
  
# Units of parameters

.Macroalgae$parms$units <- c("gSC/gSDW",
                             "molRC/molSC",
                             "molSN/molSC",
                             "molRN/molSC",
                             "molRN/molSC",
                             "molSP/molSC",
                             "molRP/molSC",
                             "molRP/molSC",
                             "gDW/mmolRC",
                             "gDW/mmolRN",
                             "gDW/mmolRP",
                             "gDW/mmolSC",
                             "mol O2/mol RC",
                             "gDW/gWW",
                             "molC/molSC/day",
                             "molC/molSC/(uMol photons/m2)",
                             "molC/molSC/(uMol photons/m2)",
                             "-",
                             "-",
                             "-",
                             "dgK",
                             "dgK",
                             "dgK",
                             "dgK",
                             "dgK",
                             "dgK",
                             "molC/molC",
                             "m",
                             "m2/mmolC",
                             "gWW",
                             "m",
                             "m2/gWW",
                             "/d",
                             "molC/molSC/d",
                             "-",
                             "/dgC",
                             "dgC",
                             "molN/molSC/d",
                             "mmol NO3/m3",
                             "molN/molSC/d",
                             "mmol NH4/m3",
                             "molP/molSC/d",
                             "mmol PO4/m3",
                             "m/s",
                             "/d",
                             "/d",
                             "/d",
                             "/(mmolC/m2)",
                             "/m",
                             "/m",
                             "/m",
                             "/m",
                             "/(mg Chl/ m3)",
                             "/(mg Chl/ m3)",
                             "/(mg Chl/ m3)",
                             "dgN",
                             "m",
                             "/d")


# ================================================
# State variables, with initial conditions
# ================================================
States <- c(
  # MACROPHYTES structure & reserves
  
  STRUCT_C  = 50,             # [mmolC/m2] structural C of macroalgae, determining most rates. 
  RESERVE_C = 50*0.05*2,      # [mmolC/m2]  reserve C, start= structC*2*RC_SCmin   
  RESERVE_N = 50*0.0355*1.5,  # [mmolN/m2]  reserve N, start= Struct_C*RN_SCmin*1.5
  RESERVE_P = 50*0.0011*1.5,  # [mmolP/m2] reserve P [Nutrient reserves]

  # ---ENVIRONMENT of 0D box model--- #
  
  NO3     = 10,  # [mmolN/m3] nitrate concentration in water
  NH4     = 5,   # [mmolN/m3] ammonium concentration in water
  PO4     = 1,   # [mmolP/m3] phosphate concentration in water
  
  O2      = 200,  # [mmolO2/m3] oxygen concentration in water
  
  # --- Closings to check model validity ---# 
  C_CLOSING  = 0, # all C that goes in and out of the Box model [mmolC/m2]
  N_CLOSING = 0, # all N that goes in and out of the Box model [mmolN/m2]
  P_CLOSING = 0  # all P that goes in and out of the BOX model [mmolP/m,2]
)


.Macroalgae$y <- data.frame(names = 
    names(States), initial = States)

.Macroalgae$y$description <- 
    c("structural C of macroalgae", "reserve C of macroalgae", 
      "reserve N of macroalgae", "reserve P of macroalgae", 
      "nitrate concentration in water",
      "ammonium concentration in water", 
      "phosphate concentration in water",
      "oxygen concentration in water",
      "all C that goes in and out of the Box model", 
      "all N that goes in and out of the Box model",
      "all P that goes in and out of the BOX model")
.Macroalgae$y$units <-
    c(rep("mmolC/m2", times = 2), "mmolN/m2", "mmolP/m2", 
      "mmolN/m3", "mmolN/m3", "mmolP/m3", "mmolO2/m3",
      "mmolC/m2", "mmolN/m2", "mmolP/m2")
row.names(.Macroalgae$y) <- NULL

# ================================================
# Output variables
# ================================================

var_units <- c(
  RC_SCr =  "molRC/molSC",
  RN_SCr =  "molRN/molSC",
  RP_SCr =  "molRP/molSC",
  
  WW =      "gWW/m2",
  DW =      "gDW/m2",
  
  Macroalgae_C   = "mmolC/m2",
  Macroalgae_N   = "mmolN/m2",
  Macroalgae_P   = "mmolP/m2",
  
  GrossPS        = "mmolC/m2/d",
  GrossPR        = "mmolC/mmolSC/d",
  Photosynthesis = "mmolC/m2/d",
  
  ExudFrac     = "-",
  Exudation    = "mmolC/m2/d",
  
  Basal_Resp   = "mmolC/m2/d",
  Growth_Resp  = "mmolC/m2/d",
  Respiration  = "mmolC/m2/d",
  
  Respiration_loss  = "mmolC/m2/d",
  
  NO3_sat =  "-",
  NH4_sat =  "-",
  PO4_sat =  "-",
  NU      =  "-",
  
  RC_sat  =   "molRC:molSC",
  RN_sat  =   "molRN:molSC",
  RP_sat  =   "molRP:molSC",
  
  Quota_growth_fac =  "-",
  
  Photo_Eff  =   "-",
  Light_Eff  =   "-",
  
  NO3_upt_Eff =  "-",
  NH4_upt_Eff =  "-",
  PO4_upt_Eff =  "-",
  
  NO3_upt =     "mmolN/m2/d",
  NH4_upt =     "mmolN/m2/d",
  PO4_upt =     "mmolP/m2/d",
  
  Growth      =   "mmolC/m2/d",
  Growth_rate =  "/d",
  Growth_Eff  =   "-",
  
  Net_Growth_rate  =   "/d",
  
  Mortality       =   "mmolSC/m2/d",
  Necrosis        =   "mmolSC/m2/d",
  Regression      =   "mmolSC/m2/d",
  Erosion         =   "mmolSC/m2/d",
  pReg            =   "-",
  Necrosis_rate   =   "/d",
  Erosion_rate    =   "/d",
  Mortality_rate  =   "/d",
  
  SurfPAR_E =   "uMol photons/m2/s",
  E_short   =   "uM/m2/s",
  E_mid     =   "uM/m2/s",
  E_long    =   "uM/m2/s",
  E_Alg     =   "uM/m2/s",
  K_mean    =   "/m",
  
  K_MacroAlg       =  "/m",
  K_Competition    =  "/m",
  Shading_MacroAlg =  "-",
  Shading_Competition = "-",
  Shading_Eff = "-",
  TF_gain   =   "-",
  TF_resp   =   "-",
  
  Total_C   = "mmolC/m2",
  Total_N   = "mmolN/m2",
  Total_P   = "mmolP/m2",
  
  Temperature = "degC",
  Current     = "m/s",
  PHYTO       = "mmolC/m3",
  Light       = "W/m2",
  ExtNO3      = "mmolN/m3",
  ExtNH4      = "mmolN/m3",
  ExtPO4      = "mmolP/m3",
  ExtO2       = "mmolO2/m3")

var_desc <- c(  
  RC_SCr = "Carbon (C) reserve: structural C ratio (C quotum)",
  RN_SCr = "Nitrogen (N) reserve: structural C ratio (N quotum)",
  RP_SCr = "Phosphorus (P) reserve: structural C ratio (P quotum)",
  
  WW =    "wet weight macroalgae biomass",
  DW =    "dry weight macroalgae biomass",
  
  Macroalgae_C  = "total macroalgae carbon",
  Macroalgae_N  = "total macroalgae nitrogen",
  Macroalgae_P  = "total macroalgae phophorus",
  
  GrossPS          = "total gross photosynthesis",
  GrossPR          = "gross specific Photosynthesis rate at light and T",
  Photosynthesis   = "photosynthesis - exudation",
  
  ExudFrac         = "fraction of gross photosynthesis exudated",
  Exudation        = "total C exudation",
  
  Basal_Resp       = "basal respiration paid by reserve C", 
  Growth_Resp      = "C paid for growth (activity respiration)",
  Respiration      = "total respiration (-necrosis)",
  
  Respiration_loss = "total respiration (+necrosis)",
  
  NO3_sat = "Monod limitation term for NO3 uptake",
  NH4_sat = "Monod limitation term for NH4 uptake",
  PO4_sat = "Monod limitation term for PO4 uptake",
  NU      = "dependency of nutrient uptake on current velocity",
  
  RC_sat  = "N reserve saturation of nutrient uptake and growth",
  RN_sat  = "C reserve saturation of growth",
  RP_sat  = "P reserve saturation of nutrient uptake and growth",
  
  Quota_growth_fac = "limitation of growth due to quota",
  
  Photo_Eff =  "part of max Photosynthesis realised to light only",
  Light_Eff =  "part of max Photosynthesis realised due to light and temperature",
  
  NO3_upt_Eff =  "N03 uptake efficiency (relative to max uptake)", 
  NH4_upt_Eff =  "NH4 uptake efficiency (relative to max uptake)", 
  PO4_upt_Eff =  "PO4 uptake efficiency (relative to max uptake)",
  
  NO3_upt =      "total NO3 uptake",
  NH4_upt =      "total NH4 uptake",
  PO4_upt =      "total PO4 uptake",
  
  Growth      =  "total structural growth",
  Growth_rate =  "relative growth rate",
  Growth_Eff  =  "growth efficiency",
  
  Net_Growth_rate = "net structural increase (dSTRUCTC/STRUCTC)",
  
  Mortality       = "total mortality, structural C",
  Necrosis        = "total necrosis, structural C",
  Regression      = "respiration paid by structure C",
  Erosion         = "total erosition, structural C",
  pReg            = "part respriration paid by structural C",
  Necrosis_rate   = "C-specific necrosis rate",
  Erosion_rate    = "C-specific erosion rate",
  Mortality_rate  = "C-specific mortality rate",
  
  SurfPAR_E =   "irradiance on water surface",
  E_short   =   "blue-light band (400-500nm) on algae",
  E_mid     =   "green-light band (500-600nm) on algae",
  E_long    =   "red-light band (600-700nm) on algae",
  E_Alg     =   "total irradiance for photosynthesis",
  K_mean    =   "average attenuation factor for light",
  
  K_MacroAlg    = "extinction coeff due to self shading",
  K_Competition = "extinction coeff competing algae",
  Shading_MacroAlg    = "shading factor due to modeled algae",
  Shading_Competition = "shading factor due to competing algae",
  
  Shading_Eff = "part of light not shaded by algae",
  
  TF_gain =    "temperature factor for photosynthesis and growth",
  TF_resp  =   "temperature factor for respiration and necrosis",
  
  Total_C  = "total carbon in system",
  Total_N  = "total nitrogen in system",
  Total_P  = "total phophorus in system",
  
  Temperature = "temperature",
  Current     = "current speed",
  
  PHYTO       =  "phytoplankton concentration in water column ",
  Light       =  "surface irradiation",
  ExtNO3      =  "bottom water NO3 concentration",
  ExtNH4      =  "bottom water NH4 concentration",
  ExtPO4      =  "bottom water PO4 concentration",
  ExtO2       =  "bottom water O2 concentration")

.Macroalgae$out <- data.frame(names = names(var_units),
                              units = var_units,
                              description = var_desc)
row.names(.Macroalgae$out) <- NULL

.Macroalgae$forc <- data.frame(names = 
   c("f_Temperature", "f_Current", "f_Light", "f_PHYTO",
     "f_ExtNO3", "f_ExtNH4", "f_ExtPO4", "f_ExtO2"))

.Macroalgae$forc$description <- c(
   "environmental temperature", 
   "current velocity", 
   "light at the air-water interface", 
   "phytoplankton concentration (overlying water)",
   "external nitrate concentration (overlying water)", 
   "external ammonium concentration (overlying water)", 
   "external phosphate concentration (overlying water)", 
   "external oxygen concentration (overlying water)")
.Macroalgae$forc$units <- c(
   "dg C", "m/s", "W/m2", "mmolN/m3", "mmolN/m3", "mmolP/m3", 
   "mmolC/m3", "mmolO2/m3")
  
row.names(.Macroalgae$forc) <- NULL
