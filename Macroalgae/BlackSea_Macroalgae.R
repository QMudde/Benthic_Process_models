## R model of Black Sea macroalgae biomass

# ## simplified model for macroalgal growth ## #
# based on bronch & slagstad / Jiang et al., 2022

# The model captures the 2 keystone species of black sea phytobenthos: Phyllophora & Cystoceira
# studies on flora assemblages show that these 2 species are the base organims of thier phytocoenoses 

# The model describes the red seeweads from the phyllophoraceae familiy
#    # Phyllophora crispa (prev Ph. nervosa) is the target species
          # phyllophora also captures Coccotylus trucatus (prev. Ph. trucata)

# The model describes the Cystoseira genus, which dominates the shallow subtidal
#   # Cystoceira barbata (target species) (prev Gongolaria barbata)
#   # Cystoceira crinita (Ericaria crinita)

# packages:
library(deSolve)
library(FME)

# info: 2 phyllphora fields known
#  1) large (Zernov's) Phyllophora field 
# 2) small Phyllophora field (Karkinitsky Bay)



# relevant model studies 
# Vasechkina, E. (2020). Object-based modeling of marine phytoplankton and seaweeds.

#------------------------------#
#  MODEL STRUCTURE DESCRIPTION #
#------------------------------#
# The Macroalgae model is based on reserve vs uptake dynamics
# the main variables of the algae are N & C reserve and structural biomass.
# Dynamics are driven by fluxes from environmental N & C into reserve compartments
# and structural growth & respiration based on usage of those reserves




# MODEL ASSUMPTIONS : 
# 1) state variables have fixed stoichiometry (DEB - strong homeostasis)
# 2) Macroalgae blades are very flat (single cell to few cells thic) 
#         therefore a direct proportion of STRUCTURE to photosynthetic surface is assumed. 
# 3) Environmental C is unlimited, C uptake limited by light & other factors



## simplified photosynthesis model & PI curves are from Johanssen & Snoeijs 2002
# they use the hyperbolic tanget function P = Pmax * tanh(alfa*I/Pmax) + Rd
# where P = gross photosynthesis, Pmax = max photo, alfa = inital slope, Rd = dark resp
# photosynthesis expiremnts were done on young thallus pars (assume no reserve) at 14 dC.
# photonsynthesis was measured for PAR (400-700 nm)

## IRRADIANCE REACTION CURVES OF 2 GROUPS
# ----  Luning & Dring 1985 ----- #
# what part of PAR spectra (short 400-500), mid(500-600) or long(600-700) do algea react to?
# totals must be 1, these are relative photosynthesis rate spectra
Phyllo_short = c(0.05,0.05,0.05,0.05,0.1,0.16,0.39,0.5)
Phyllo_mid = c(0.5,0.52,0.59,0.55,0.6,0.55,0.39)
Phyllo_long = c(0.35,0.34,0.29,0.25,0.2,0.18,0.05,0.05)
# calculate susceptebilities to wavelength spectra of Phyllophora
mean(Phyllo_short)/sum(c(mean(Phyllo_short),mean(Phyllo_mid),mean(Phyllo_long))) # 0.19
mean(Phyllo_mid)/sum(c(mean(Phyllo_short),mean(Phyllo_mid),mean(Phyllo_long))) # 0.58
mean(Phyllo_long)/sum(c(mean(Phyllo_short),mean(Phyllo_mid),mean(Phyllo_long))) # 0.23

# same for LAMINARIA Saccharia, Proxy for Brown Weed Cystiseira
Lamin_short = c(0.42,0.54,0.49,0.52,0.55,0.56,0.52)
Lamin_mid = c(0.47,0.43,0.42,0.39,0.35,0.32,0.33)
Lamin_long = c(0.34,0.36,0.36,0.37,0.4,0.45,0.42,0.35)
# calculate susceptebilities to wavelength spectra of Phyllophora
mean(Lamin_short)/sum(c(mean(Lamin_short),mean(Lamin_mid),mean(Lamin_long))) # 0.4
mean(Lamin_mid)/sum(c(mean(Lamin_short),mean(Lamin_mid),mean(Lamin_long))) # 0.3
mean(Lamin_long)/sum(c(mean(Lamin_short),mean(Lamin_mid),mean(Lamin_long))) # 0.3


## 
# For now: all N is assumed to be NO3, (output of bamhbi, also assumed in Broch & Slagstad 2012)
##

# for now: we assume a simple temperature response, an optimum curve for photosynthesis, 
# and exponentially increasing curve for growth & respiration. 
# temperature response based on formulations from: Frieder etal., 2022? or 
# for photosynthesis : 3 pars 
# 1) T_optL (lower boundary of optimal temp range)
# 2) T_optH (higher boundary of optimal temp ragne)
# 3) T_max (upper maximum temperature for photosythesis, above this therhold, no more C fixation)
# IF T < T_optL f(t) = T/ToptL
# IF T_optL<T<T_optH f(T) = 1
# IF T_optH < T < T_max f(T) = 1-(T-T_optH)/(T_max-T_optH)
# IF T > T_max f(T) = 0

# Temperature response for growth & respiration : f(T) 1/(1+exp(-(T-T0)/Tr))
# with Tr = reference temperature for physiological rate parameter
# T0 = 
#      1/(1+exp(-(12-10)/1))
#       1/(1+exp(-ZETA*(10-THETA)))
#       1/(1+exp(-0.3*(10-10)))
# S curve f(T) = 1/(1+exp(-T1*(T-T2)))
# where T1 is the relative change of metabolic rates per degree C change
# where T2 is the reference temperature where the S curve hits exactly f(T) = 0.5

## COVERAGE biomass relationship Cystoseira branches from Orfanidis et al., 2017
# bioms = seq(0,1000) # wet weight grams biomass
# covs = (99.59 - 98.16*0.985^bioms) # allometric relationship
# covdens =  –7.47 + 1.63CVe–0.006CVe2
# plot(x = bioms, y = covs)


# Parameters 
ALG_pars <- c( 
  
  # parameter shared between species
  # for now equal recalc factors between species assumed.
  # constant recalculation factors
  gdw_SC    = 5*12     ,#[gdw/molC]   gram dry weight per mol structural carbon (BS12)
  gdw_RN    = 2.72*14  ,#[gdw/molN]   gram dry weight per mol reserve nitrogen  (BS12)
  gdw_RC    = 2.1213*12,#[gdw/molC]   gram dry weight per mol reserve carbon    (BS12)
  gdw_gww   = 0.0785   ,#[gdw/gww]    gram dry weight per gram wet weight       (BS12)
  SNCr      = 0.061  ,  #[molN/molC]   Structural NC ratio (StrucN/SrucC) (Nstruct in BS12)  
  OC_q = 1.105, # [molO2 repired per molC (reserves) respired] gerard et al., 1988
  
  
  #------------------
  # Phyllophora pars (mostly phyllophora crispa)
  #------------------
  
  # physical shape parameters 
  phy_DWC = 0.024, # Dry weigtht coupler for phyllophora structure [gDW/mmolC] (50%)
  phy_DW_WW = 0.195, # gram dw per gram ww (default ratio for saccharina = 0.08) (HERE: FUCUS distichus , Wickham et al., 2019)
  phy_COV_WWa = 0.2861739, # [ percentage cover per biomass WW for understory alga (Jung & Choi 2022)
  phy_COV_WWexp = 0.8474576, # percentage cover per biomass WW for understory alga (Jung & Choi 2022)
  # coverage = a*WW^exp in [%/m2](Jung & Choi 2022)
  phy_BCr = 0.0005, # COVERAGE biomass to seabed coverage ratio [m2/mmolC]
  #phy_selfS = 0.0003*61, # Shading, m2 shaded per mmolC [m2/mmolC]
  # phy_selfS = 1/28,
  # reserve quota
  phy_QNmin = 0.053, 
  phy_QNmax = 0.092, 
  phy_QCmin = 0.05, 
  
  # Core physiolocial rates (regarding external)
  #phy_PRmax = 
  #phy_PRmax = 0.24    , #[/d]  specific photosynth rate at T_P  (0.36 vdM18)
  # gross max photosynthesis at 15 deg at E of 10 [umol photons/m2/s] = 0.83 (mean 0.406)
  # basal repsiration = 0.145 umolO2/m2/s (both Luning & Dring 1985)
  phy_pMax = 43.7, # light saturated NET photosynthetic rate [uMol O2/kg DW/s]
  # 43.7 phyllophora crispa [muMol O2/ kg DW/sec] johannsen & snoeijs)
  phy_RD = 4.1, # dark repsiration rate [uMol O2/kgDW/s]
  phy_alpha = 0.99, # Initial PI curve slope  [(uMol O2/kgDW)/uMol photons/m2]
  phy_ExhFmax = 0.5   , #[molC/molC]   maximal fraction photosynthesis exudated (gamma in BS12) 
  
  phy_MaxNupt = 0.04062857 , #[mmolN/gSDW/d] max uptake nitrate phyllophora
  # Vmax = 23.7 ugram N/gDW/h (Wallentius 1984) 23.7/1000/14*24-> 0.04062857 [mmolN/gSDW/d] (14dgC)
  #1.69 umol/g/h = 40.56 umol/g/d = 0.04062857
  phy_ksDIN = 9.214286       , #[mmol NO3/m3] half-saturation constant for DIN uptake (nitrate)
  # ksN = 129 ugram N/liter (Wallentius 1984) -> 129/1000/14*1000 = 0.2211429 [mmolN/m3]
  phy_u065 = 0.03     , #[m/s]         currentspeed where nutrient uptake=0.65*rNup  (BS12)
  
  # Core physiolocial rates (regarding internal)
  phy_GRmax = 0.1 , #[/d]    maximal specific growth of structural biomass (0.36 vdM18) note 0.18 in BS12
  phy_RRbase = 0.02, #[/d]    basal respiration rate
  phy_growthRespFrac = 0.05, #[-] fraction of C that is repired for structural growth (NOREF)
  # (-)4.1 [muMol O2/ kg DW/sec] phyllophora crispa dark resp (johannsen & snoeijs)
  phy_MR     = 0.008, # [/d] daily, dry weight specific mortality rate (daily avg from Duarte & Ferreira +_ 0.8% mort per day (diff sp))
  
  # TEMPERATURE PARAMETERS  
  # T curve for photosynthesis & growth
  phy_T_optL = 10, #
  phy_T_optH = 18, #
  phy_T_max = 23 , #
  
  # Temperature S curve for respiration  S curve f(T) = 1/(1+exp(-T1*(T-T2)))
  phy_Tr1 = 0.2, #
  phy_TRDr = 15, # reference temperature at which dark respiration is measured
  
  
  # Effect of light on photosynthesis
  phy_rho_short = 0.19, #
  phy_rho_mid = 0.58, #
  phy_rho_long = 0.23, #
  
  # Effect of photoperiod (PPP = PhotoPeriodParameter)
  phy_PPP1 = 0.85,
  phy_PPP2 = 0.3,
  
  phy_alpha     = 0.04 , # alpha at start of PI curve (Vasechkina & Filippova)
#  phy_Isat      = 44    , # [uMol photos/m2/s]  irradiance for maximal photosynthesis (Johansen & snoeijs 2002)
  phy_PARfr     = 0.5, # fraction of irradiance used by phyllophora
  phy_Plong     = 0.65, # fraction of Phyllophora PAR that is long waves
  # Erosion (Bronch & Slagstad, 2012)
  phy_Nec_rate = 0, # heat related necrosis parameter (no info on necrosis, default = 0)
  phy_Ero_max =  1, # max erision rate per structural biomass 
  phy_Ero_inc = 0.22, # increase rate of erosion as surface area increases
  
  #-----------------
  # Cystoseira pars (mostly cystoseira barbata)
  #-----------------
  
  # physical shape parameters 
  cys_DWC = 0.024, # Dry weight - Carbon coupler for Cystoseira structure [gDW/mmolC] (50%)
  cys_DW_WW = 0.183, # gram dw per gram ww (Orfanidis et al., 2017)
  cys_WW_COVa = 99.59, # allometric ratio wet biomass- ground cover (Orifanidis et al., 2017)
  cys_WW_COVb = 98.16, # allometric ratio wet biomass- ground cover (Orifanidis et al., 2017)
  cys_WW_COVc = 0.985, # allometric ratio wet biomass- ground cover (Orifanidis et al., 2017)
  # formula for cover (max 100%) -> cover[%] = a - b*c^WW   (with WW = wet weight/biomass in grams)
  cys_WW_CanoDa = 3.08, # allometric ratio wet biomass- canopy density (Orifanidis et al., 2017)
  cys_WW_CanoDexp = 0.685, # allometric ratio wet biomass- canopy density (Orifanidis et al., 2017)
  # formula for canopy density (can exceed 100%) -> CanoD[%] = a*WW^exp (with WW = wet weight/biomass in grams)
  cys_BCr = 0.0005, # COVERAGE biomass to seabed coverage ratio [m2/mmolC]
  #cys_selfS = 0.0001*61, # Shading, m2 shaded per gSDW [m2/mmolC]
  cys_selfS = 1/50,
  # reserve quota
  cys_QNmin = 0.053, 
  cys_QNmax = 0.092, 
  cys_QCmin = 0.05, 
  
  # Core physiolocial rates (regarding external)

  # 114.2 umol/h/gDW --> 114.2*1000/60/60 == 31.72222 [umol O2/kgDW/s]  (Baghdadli et al., 1990)
  cys_pMax = 31.72222, # light saturated NET photosynthetic rate [uMol O2/kg DW/s] 

  cys_alpha     = 1.043265 , # alpha at start of PI curve (Vasechkina & Filippova) (recalc factor 0.1150235)

  cys_ExhFmax = 0.5   , #[molC/molC]   maximal fraction photosynthesis exudated (gamma in BS12) 

  cys_MaxNupt = 0.2076, # [mmolN/gSDW/d]
  # 8.65 umol/gDW/h (*24/1000 -> 0.2076) (Vmax NO3) (Vasechkina & Fillipova 2020)
  cys_ksDIN = 9.08  , #[mmol/m3]     half-saturation constant for DIN uptake
  # 9.08 umol/liter ksNo3 (Vasechkina & Fillipova 2020) 
  cys_u065 = 0.03         , #[m/s]         currentspeed where nutrient uptake=0.65*rNup  (BS12)

  # Core physiolocial rates (regarding internal)
  cys_GRmax = 0.2 , #[/d]    maximal specific growth of structural biomass (0.36 vdM18) note 0.18 in BS12
  cys_RD  = 4.432432, # dark repsiration rate [uMol O2/kgDW/s]
  cys_RRbase = 0.02, #[/d]    basal respiration rate
  cys_growthRespFrac= 0.05, #[-] fraction of C that is repired for structural growth (NOREF)

  # (-)4.1 [muMol O2/ kg DW/sec] phyllophora crispa dark resp (johannsen & snoeijs)
  cys_MR     = 0.01, # [/d] daily, dry weigth specific mortality rate (daily avg from Duarte & Ferreira +_ 0.8% mort per day (diff sp))
  
  
  # TEMPERATURE PARAMETERS  
  # T curve for photosynthesis & growth
  cys_T_optL = 10, #
  cys_T_optH = 18, #
  cys_T_max = 23 , #

  # Temperature S curve for respiration  S curve f(T) = 1/(1+exp(-T1*(T-T2)))
  cys_Tr1 = 0.2, #
  cys_TRDr = 15, # reference temperature at which dark respiration is measured

  # Effect of photoperiod (PPP = PhotoPeriodParameter)
  cys_PPP1 = 0.85,
  cys_PPP2 = 0.3,
  
  # Effect of light on photosynthesis
  cys_rho_short = 0.4, #
  cys_rho_mid = 0.3, #
  cys_rho_long = 0.3, #

  cys_Isat      = 200    , # uE/m2/s        irradiance for maximal photosynthesis 
  cys_PARfr     = 0.5,    # fraction of irradiance used by phyllophora
  cys_Plong     = 0.65,   # fraction of cystoseira PAR that is long waves
  # Erosion (Bronch & Slagstad, 2012)
  cys_Nec_rate = 0, # cystoseira heat related necrosis parameter (no info on necrosis, default = 0)
  cys_Ero_max =  1, # max erision rate per structural biomass 
  cys_Ero_inc = 0.22, # increase rate of erosion as surface area increases
  
  #==========================#
  # ENVIRONMENTAL PARAMETERS #
  #==========================#

  #-------#
  # Light # Split up in 3 wavelength regions of PAR (400 - 700 nm)
  #-------# short - 400-500 (blue) nm , mid 500-600 (green) nm, long 600-700 nm (red)

  # the 3 wavelengths have different (base) attenuation coefficients with depth
  # the 3 wavelengths have equal addition attenuation due to turbidity (plankton/detritus)
  # the second assumption is equal to bamhbi model. 

  # Light  (0.0384 for clear water) Lorenzen 1972
  # Bamhbi (0.09 for short wavelenghts) (0.23 for long wavelengths (>600 nm)) (just short/long)
  # Specific attenuation for 3 different spectra from Tait et al., 2014  (k0short = 0.235,k0mid   = 0.16,k0long  = 0.34 )
  # relative attenuations would be (short ~ 0.9591837, mid ~0.6530612, long ~1.387755)
  # specific alga susceptibility for wavelengths (Luning & Dring 1985)
  # averages are applied for the values
  
  k0Clear = 0.0384, # PAR range attenuation for clear water (Lorenzen 1972) [-m]
  k0_turb_fact = 3.5, # factor of increase for light attenuation due to CDOM, SIM & POM (equal for all PAR)
  # k0_turb_fact of 3.5 is based on offshore attenuation values from Nababan 2021

  # attenutation coefficients measured by NAbaban et al., 2021 offshore waters (short ~ 0.04, mid ~ 0.06 long ~ 0.3) 
  # relative attentuations would be (short ~ 0.3, mid ~ 0.45, long ~ 2.25)

  rk_short = 0.3, # [-] light extinction rate of 400-500nm [relative to general PAR]
  rk_mid   = 0.45,  # [-] ight extinction rate of 500-600nm [relative to general PAR]
  rk_long  = 2.25, # [-] light extinction rate of 600-700nm [relatie to general PAR]

  # absorption of wavelengths by chlorophyll in seawater (attenuation value per depth intergrated mg Chla/m2 at depth)
  # all chl A attenuation values are from Gallegos et al., 1990 (and values are average per waveband)
  # mean(c(0.0203, 0.0218, 0.0239, 0.0250, 0.0261, 0.0270, 0.0286, 0.0298,
  #  0.0298, 0.0284, 0.0268, 0.0265, 0.0268, 0.0262, 0.0249, 0.0233, 0.0217, 0.0207, 0.0207, 0.0197)) ## 400-495 ## 0.0249

  # mean(c(0.0182, 0.0166, 0.0151, 0.0138, 0.0125, 0.0112, 0.0103, 0.0097, 0.0088, 0.0083, 0.0074, 0.0067,
  # 0.0059, 0.0054, 0.0051, 0.0049, 0.0050, 0.0052, 0.0052, 0.0050)) ## 500-595 ## 0.009015
  # mean(c(0.0049, 0.0049, 0.0052, 0.0061, 0.0066, 0.0069, 0.0072, 0.0073, 0.0075, 0.0071, 0.0077, 0.0070, 
  # 0.0088, 0.0123, 0.0159, 0.0180, 0.0176, 0.0145, 0.0091, 0.0047,0.0022)) ## 600-700 ##  0.008642857
  #
  k_Phyto_short = 0.0249,# [/mg Chl/ m2] extra light extinction of phytoplankton on short waves
  k_Phyto_mid = 0.009015,# [/mg Chl/ m2] extra light extinction of phytoplankton on mid waves
  k_Phyto_long = 0.008642857,# [/mg Chl/ m2] extra light extinction of phytoplankton on long waves
  # bamhbi uses 0.0003 m2/mg chl | Lorenzen 1972 : 0.0138 (/depth intergrated mg chlorA) over all PAR
  # average plankton light attenuation = (0.0249+0.009015+0.008642857)/3 == 0.01418595
  # Kext = Kchl /(mg*m2) * Chlorophyll [mg/m3] ==> /m  

  # additional information about environment
  latitude = 45, # determines daylight 
  depth = 10, # determines light that reaches bottom
  dilution = 0.25, # dilution of DIN with environemtn/day (0.25 used in Hadley et al., 2015)
  dw_START = 1 # g/m2 dry wieght each macroalga stars with in model
  
)

# 350 g WW for canopy cover 100%..? tait et al., 2014 H. banksii

# state variables have fixed chemical compositions (strong homeostasis- DEB theory)
# start macroaglae with full reserves (based on Qmax quota from pars & starting DW)
StatVars <- c(
  # MACROPHYTES (PHY = Phyllophora [[Phyllophora crispa]])
  PHY_SDW = ALG_pars[["dw_START"]], 
  # [g/m2] Gram structural dry weight Phyllophora per m2 (constant NC ratio)
  PHY_RESERVE_C = ALG_pars[["dw_START"]] / ALG_pars[["phy_DWC"]] * (2*ALG_pars[["phy_QCmin"]]), 
  # [mmolC/m2]  reserve C, additionall to structural biomass
  PHY_RESERVE_N = ALG_pars[["dw_START"]] / ALG_pars[["phy_DWC"]] * ((ALG_pars[["phy_QNmin"]] + ALG_pars[["phy_QNmax"]])/2),  
  # [mmolN/m2] reserve N, additional to sturctural biomass * N/C
  # MACROPHYTE ( CYS = Cystoseira [[Cystoceira barbata]])
  CYS_SDW = ALG_pars[["dw_START"]], 
  # [g/m2] Gram structural dry weight Phyllophora per m2 (constant NC ratio)
  CYS_RESERVE_C = ALG_pars[["dw_START"]] / ALG_pars[["cys_DWC"]] * (2*ALG_pars[["cys_QCmin"]]), 
  # [mmolC/m2] reserve C, additionall to structural biomass
  CYS_RESERVE_N = ALG_pars[["dw_START"]] / ALG_pars[["cys_DWC"]] *  ((ALG_pars[["cys_QNmin"]] + ALG_pars[["cys_QNmax"]])/2), 
  # [mmolN/m2] reserve N, additional to sturctural biomass * N/C
  
  # ---ENVIRONMENT of 0D box model--- #
  
  DIN = 1,    # [mmolN/m3] connected to environment (can differ factor 100 in BS between regions)
  OXYGEN = 360   # [mmolO2/m3] connected to environment
)

# Forcing funcitons
# the model takes into account temperature ,photo period, light intensity , currents, Nutrients 

# daylight fucntion from Forsythe et al., 1995


fTemp   <- function(t) 15  # temperature at bottom at timestep
fLight  <- function(t) 1000 # surface irridiation at timestep [umol photons/m/s]
# recalculate irradiance with planck? N_p = I / (h * c / λ)
# NP = photons/second /m2 
# I = Irridiance (W/m2)
# h = planck constant 6.626 x 10^-34 J s
# c = speed of light 2.998 x 10^8 m/s
# λ = wavelength in meters 
IR_to_Photons = function(IR, lower = 400, upper = 700){
  # calculate Photon flux (micro mol photons / m2 /s) (unit that is used in most experiemnts)
  # from
  # Irradiance in W/m2 (unit that comes from raw data)
  # for W sunlight to photon flux, +- 4.6 muEinst/W can be used. 
  range = seq(from = lower, to = upper, by = 1)
  photonF = 0
  for(lambdas in range){
  wavelength = IR*lambdas*0.836/100
  photonF = photonF+wavelength
  }
  return(photonF)
}

PhotoP <- function(DOY,lat){
  # based on CBM model Schoolfield (1982)
  # theta = revolution angle based on DOY, phi = declination angle , daylength include twilight
  # ASSUMED: Sunrise/Sunset is when the center of the sun is even with the horizon (Forsythe et al., 1995)
  theta = 0.2163108 + 2 *atan(0.9671396 * tan(0.00860 * (DOY-186)))
  phi = asin(0.39795 * cos(theta))
  daylength_H = 24 - 24/pi*acos( (sin(lat*pi/180)*phi )/(cos(lat*pi/180)*cos(phi))    )
  return(daylength_H)
}


PhotoP(100,45)
PhotoP(1,45)
PhotoP(355/2,45)

fCurr   <- function(t) 0.1 # net currentspeed at bottom layer at timestep 
fDIN    <- function(t) 0.1 # DIN calculated (BAMHBI) for bottom layer at timestep
fPHYTO  <- function(t) 2.5 # avg mmolC/m3 phytoplankton in water column at timestep
fOX     <- function(t) 360 # external dissolved oxygen concentration mmolO2/m3

heatwave <- data.frame(times = seq(1,2*365), 
                       temp = rep(10,2*365))

(heatwave["times"]%%365 > 210 & heatwave["times"]%%365 < 230)

heatwave$heatwave <- (heatwave["times"]%%365 > 210 & heatwave["times"]%%365 < 230)
heatwave$temp2 <- 10*heatwave$heatwave + heatwave$temp

## MODEL FORMULATION

MacroAlg <- function (t, y, parms, 
                      fTemp, fCurr, fLight, fDIN, fPHYTO,fOX){
  
  with (as.list(c(y, parms)),{
    
    #====================#
    # SHARED ENVIRONMENT #
    #====================#
    
    # forcing values at time t
    # DOY = t%%365 # day of year (assume t1 is jan01)
    
    # change in daylength affects growth rate, B&S 2012. eg. a decreasing daylength decreases growth
    # this way maximum growth is less in winter than in summer
    
    Daylength = PhotoP(t%%365,latitude) # daylength based on doy and latitude
    PastDaylength = PhotoP((t-1)%%365,latitude)
    DeltaDL = Daylength-PastDaylength
    # normalized delta daylength 
    max_DeltaDL = max((PhotoP(1:365, latitude))-(PhotoP((1:365)-1, latitude)))
    normDL = DeltaDL/max_DeltaDL # normalized delta daylength relative to maximum
    
    TempC     = fTemp(t) # temp from forcing
    TempK     = TempC + 273.15 #in kelvin for calculations
    CurrentV  = fCurr(t) # current speed based on forcings
    SurfLight = fLight(t)# surface irridation from forcing [W/m2 total light spectrum]
    # here it is assumed the light input is W/ms for +- 300-2600 nm 
    extDIN       = fDIN(t) # bottom nutrient content form forcing (external)
    PHYTO     = fPHYTO(t) # average phytoplankton concentration in water column
    extOX     = fOX(t)# external Oxygen concentration at time t
    
    # can be multiplied by depth in m to approximate depth intergrated phytoplankton biomass. 
    
    

    
    #E    = E*(exp(-Light_ext_phyto*PHYTO*depth))            # shading by PHYTO in water column
    #E    = E* exp(- (phy_selfS*PHY_STRUCTURE + cys_selfS*CYS_STRUCTURE)) # shading by MacroAlg
    
    #--------------------------------------------------------------------#
    # calculate values that are required for both macroalgae dynamics    #
    #--------------------------------------------------------------------#

    # surface coverage (no additional canopy coverage for phyllophora assumed)
    WW_phy = PHY_SDW/phy_DW_WW 
    Phy_coverage = phy_COV_WWa * WW_phy^phy_COV_WWexp # % coverage with phyllophora (understory)
    # cystoseira wet weight and coverage
    WW_cys = CYS_SDW/cys_DW_WW
    Cys_coverage = cys_WW_COVa-cys_WW_COVb*cys_WW_COVc^WW_cys # % coverage with cystoseira
    # total coverage for spatial competition (logistic growth)
    tot_coverage = min(100, (Phy_coverage+Cys_coverage))
    # shading by cystoseira canopy 
    Cys_canopy_dens = cys_WW_CanoDa * WW_cys^cys_WW_CanoDexp # [% canopy density, can exceed 100% if canopy is overlapping]
    Canopy_Shading = max(0, ((Cys_canopy_dens+Phy_coverage-100)/(Cys_canopy_dens+Phy_coverage))) # assume phy coverage == phy canopy
    # specific light wavelengths of light captured by canopy not accounted for
    
    
    ####################
    ## LIGHT AT DEPTH ##
    ####################
    
    # first, Total irradiance [W/m2] (300-2600nm) is converted to PAR [W/m2] (400-700nm)
    SurfPAR_Wm2 = SurfLight * 0.43 # about 43% of total radiation spectrum is PAR (400-600nm)
    # second, irradiance [W/m2] is recalculated to photon flux [uMol photons/m2/s] via a broad estimate 
    SurfPAR_E = SurfPAR_Wm2 * 4.57 # typical conversion factor for sunlight W/m2 to umol photons in PAR range
    ## -- apply self shading as relative filter on surface par, equal effect over PAR spectrum
    # shading is caused by canopy forming macroalagae cystoseira,(not phyllophora where coverage = shading)
    PAR_E = SurfPAR_E * (1-Canopy_Shading)            #(1- (cys_selfS*CYS_SDW))
    
    # PAR will be divided into 3 bands of wavelengths short: 400-500 mid:500-600 and long 600-700
    # assume that at sea surface, these 3 bands are of equal photon flux
    E_400_500 = PAR_E*(1/3) # blue light band of wavelengths
    E_500_600 = PAR_E*(1/3) # green light band of wavelengths
    E_600_700 = PAR_E*(1/3) # red ligth band of wavelengths
    
    
    #
    k0short = k0Clear * k0_turb_fact * rk_short
    k0mid   = k0Clear * k0_turb_fact * rk_mid
    k0long  = k0Clear * k0_turb_fact * rk_long
    
    # Light based on depth, phytoplankton turbidity, all 3 wavebands attenuate.
    # 3 spectra have different background attenuation
    # all 3 spectra are equally affected by phytoplankton (turbidity) ASSUMPTION.
    # all 3 spectra are differently affected by phytoplankton light extionction
    
    # Cloern 1995: C/Chlor ratio of 30 , (factor 2.5)
    ChlA <- PHYTO/2.5 # [avg mg/m3 chlorophyll A] (or ug/l)
    
    E_short = E_400_500*(exp(-(k0short+(k_Phyto_short*ChlA))*depth))
    E_mid = E_500_600*(exp(-(k0mid+(k_Phyto_mid*ChlA))*depth))
    E_long = E_600_700*(exp(-(k0long+(k_Phyto_long*ChlA))*depth))
    
    #======================#
    # PHYLLOPHORA DYNAMICS #
    #======================#
    # constant proportions based on State Variables ,quota
    Phy_STRUC_C = PHY_SDW/phy_DWC 
    Phy_STRUC_N = Phy_STRUC_C * SNCr
    Total_Phy_C = Phy_STRUC_C + PHY_RESERVE_C
    Total_Phy_N = Phy_STRUC_N + PHY_RESERVE_N
    # current state of phyllophora (quota)
    Phy_totNC = Total_Phy_N/Total_Phy_C # total N/C ratio
    Phy_QC = PHY_RESERVE_C/Phy_STRUC_C # current energy (C) reserve quotum
    Phy_QN = PHY_RESERVE_N/Phy_STRUC_C # current nutrient (N) reserve quotum



    #-----------------------
    # UPTAKE & GROWTH TERMS  (Liebig’s law)
    #-----------------------
    
    # calculate temperature factor (TF) effect on growth and photosynthesis
    Phy_grow_TF = 0.0
    if (TempC < phy_T_optL) 
      Phy_grow_TF = TempC/phy_T_optL
    else if (TempC < phy_T_optH) 
      Phy_grow_TF = 1.0   # OPTIMAL TEMP FOR GROWTH = 10<T<18 here  (between optTL & optTH)
    else if (TempC < phy_T_max) 
      Phy_grow_TF = 1-(TempC-phy_T_optH)/(phy_T_max-phy_T_optH)
    
    Phy_grow_TF  = max(0.0, Phy_grow_TF)
    
    # Primary production (carbon from environment to reserve pool)
    # max photosynthetic rate at current temperature
    Phy_Pmax = phy_pMax*Phy_grow_TF
    
    # ---effect of actual irradiance on photosynthesis based on Phyllophora spectra ---
    E_Phyllophora = phy_rho_short*E_short + phy_rho_mid*E_mid + phy_rho_long*E_long
    
    #
    ## Classic PI curve (gross P, no RD subtraction)  
    Phy_grossPR <- Phy_Pmax * tanh(phy_alpha*E_Phyllophora/Phy_Pmax) # [uMol O2/kg DW/s]
    #
    #  1/ 1000 / 1000 * 60*60*24 == 0.0864 # uMol/kg//s -> mmol/g/d
    #
    # Calulate Gross photosynthesis in mmolC/d with DW in grams
    Phy_grossPS     = Phy_grossPR * PHY_SDW / OC_q * 0.0864 #[mmolC/d]
    
    #     Carbon exudation fraction  (surplus of C reserve is excreted constantly)
    #     calculated is a fraction of phosotynthesis that is immediately exudated.
    # formulation and parameter following Broch & Slagstad
    
  ###  Phy_ExudFrac = 1.0 - exp(-phy_ExhFmax * max(0.0, Phy_QC - phy_QCmin))
    Phy_ExudFrac = 1.0 - exp(phy_ExhFmax * min(0.0, phy_QCmin - Phy_QC))
    
    
    Phy_Exudation         =  Phy_grossPS * Phy_ExudFrac
    
    Phy_Photosynthesis    =  Phy_grossPS - Phy_Exudation
    
    #---------------------------------------------------------------------
    #  Nutrient uptake (Nutrients from environment to Nutrient reserve)
    #---------------------------------------------------------------------
    #     limited by external DIN and internal quota, dependent on current speed 
    #     if negative: excretion of N  
    
    # effect of reserveN:structureC quota (droop kinetics) on N uptake
    # described as N reserve saturation (between 0 and 1)
    Phy_RN_sat = min(1.0, (phy_QNmax-Phy_QN) / (phy_QNmax-phy_QNmin)) #[-]
    Phy_RN_sat = max(0,Phy_RN_sat) 
    
    # saturation of uptake based on seawater dissolved nutrients (between 0-1)
    Phy_DIN_sat   = DIN/(DIN+phy_ksDIN) # [-]
    
    # functional dependency of nutrient uptake on current velocity (between 0-1)
    Phy_NU = 1.0 - exp(-CurrentV/phy_u065)  #[-]           
    
    # Specific nitrogen uptake rate and total nitrogen uptake
    Phy_Nuptake_rate = phy_MaxNupt * Phy_RN_sat * Phy_DIN_sat * Phy_NU  #[mmolN/gSWD/d]          
    Phy_Nuptake      = Phy_Nuptake_rate * PHY_SDW *2              #[mmolN/m2/d] 
    # multiplied by 2: assume Nh4 uptake is equal to NO3 uptake
    
    #----------------------------------------------------------
    # STRUCTURAL GROWTH (Based on Reserve pools) & temperature
    #----------------------------------------------------------
    # Reserves get catabolised and assimilated into structural tissue 
    
    # functional dependency of growth on reserve:structure quota
    # the lowest current quotum limits the growth (C or N storage) via droop kinetics
    
    #    fdGrowthRN_SC = min(1.0, max(0.0, 1.0 - minRN_SC/RN_SC))  # for output only
    #    fdGrowthRC_SC = min(1.0, max(0.0, 1.0 - minRC_SC/RC_SC))
    
    Phy_quota_growth_lim  = (1.0 - max(phy_QNmin/Phy_QN, phy_QCmin/Phy_QC))
    Phy_quota_growth_lim  = max(0.0, Phy_quota_growth_lim) 
    Phy_quota_growth_lim  = min(1.0, Phy_quota_growth_lim)
    
    # effect of daylength (effect of photo period)
    Phy_season_effect = phy_PPP1 * (1+ sign(normDL * abs(normDL)^0.5))+phy_PPP2 
    
    
    # specific growth rate and total growth
    Phy_Growth_rate = phy_GRmax * # Phy_season_effect *
      Phy_grow_TF *Phy_quota_growth_lim         #[/d]        
    
    Phy_Growth      = Phy_Growth_rate * PHY_SDW * (100-tot_coverage)/100    #[gSDW/ind/d]
    
    #====================================================#
    # LOSS TERMS (Both losses in reserves and structure) #
    #====================================================#
    
    #--------------------------
    # Respiration
    #--------------------------
    
    # functional dependency of respiration on temperature
    # Temperature S curve for respiration  S curve f(T) = 2/(1+exp(-T1*(T-TRDr)))
    Phy_TF_Resp    = 2/(1+exp(-phy_Tr1*(TempC-phy_TRDr)))    #[-]
    
    #total mmolC needed for respiration
    Phy_Resp    = Phy_TF_Resp * phy_RD * PHY_SDW / OC_q*0.0864 #[mmolC/d] respired
    
    # respiration is paid by structure if there are not enough reserves (here: regression)
    Phy_pReg  = 1.0 - Phy_QC^3/(Phy_QC^3+phy_QCmin^3) # goes to 1 when C reserve are at minimum.
    Phy_pReg = Phy_pReg*2 # goes to 1 when C reserve are at minimum.
    Phy_pReg = min(1,Phy_pReg) 
    Phy_pReg = max(0,Phy_pReg)
    
    Phy_Regression = Phy_pReg * Phy_Resp# [mmolC/d] paid by structure, for required respiration]
    Phy_ResResp = Phy_Resp*(1-Phy_pReg) # repsiration paid by C reserves [mmolC/d]
    
    #     activity respiration and total respiration 
    Phy_Growth_Resp = phy_growthRespFrac*(Phy_Growth/phy_DWC) # mmolC/d paid for struct growth
    Phy_Total_Resp = Phy_ResResp + Phy_Growth_Resp                               
                                      
    #------------------------------------------------
    # Erosion (erosion affects structure & reserve)
    #------------------------------------------------
    
##   Erosion_rate = max_Erosion_Rate *1e-6*exp(Erosion_surf_incr*SurfaceArea)/                   
##      (1.+1e-6*(exp(Erosion_surf_incr*SurfaceArea)-1.)) # erosion [mmolC/mmolC/d/m2]  
##    Erosion      = Erosion_rate*STRUCTURE      #[mmolC/(ind or m2?)/d]
    
    #----------------------
    # Mortality & Necrosis
    #----------------------
    
    Phy_Necrosis = phy_Nec_rate * Phy_TF_Resp * (TempC > phy_T_optH) # dry weight specific necrosis mort
    Phy_Mortality = (Phy_Necrosis+phy_MR) * PHY_SDW # total mortality in gdw/d
    
    
    
    #=====================#
    # CYSTOSEIRA DYNAMICS #
    #=====================#
    
    # constant proportions based on State Variables ,quota
    
    Cys_STRUC_C = CYS_SDW/cys_DWC 
    Cys_STRUC_N = Cys_STRUC_C * SNCr
    Total_Cys_C = Cys_STRUC_C + CYS_RESERVE_C
    Total_Cys_N = Cys_STRUC_N + CYS_RESERVE_N
    # current state of cystoseira (quota)
    Cys_totNC = Total_Cys_N/Total_Cys_C # total N/C ratio
    Cys_QC = CYS_RESERVE_C/Cys_STRUC_C # current energy (C) reserve quotum
    Cys_QN = CYS_RESERVE_N/Cys_STRUC_C # current nutrient (N) reserve quotum

    
    # additional cystoseira canopy density (for shading)  
    
    #-----------------------
    # UPTAKE & GROWTH TERMS  (Liebig’s law)
    #-----------------------
    
    # calculate temperature factor (TF) effect on growth and photosynthesis
    Cys_grow_TF = 0.0
    if (TempC < cys_T_optL) 
      Cys_grow_TF = TempC/cys_T_optL
    else if (TempC < cys_T_optH) 
      Cys_grow_TF = 1.0   # OPTIMAL TEMP FOR GROWTH = 10<T<18 here  (between optTL & optTH)
    else if (TempC < cys_T_max) 
      Cys_grow_TF = 1-(TempC-cys_T_optH)/(cys_T_max-cys_T_optH)
    
    Cys_grow_TF  = max(0.0, Cys_grow_TF)
    
    # Primary production (carbon from environment to reserve pool)
    # max photosynthetic rate at current temperature
    Cys_Pmax = cys_pMax*Cys_grow_TF
    
    # ---effect of actual irradiance on photosynthesis based on Phyllophora spectra ---
    E_Cystoseira = cys_rho_short*E_short + cys_rho_mid*E_mid + cys_rho_long*E_long
    
    #
    ## Classic PI curve (gross P, no RD subtraction)  
    Cys_grossPR <- Cys_Pmax * tanh(cys_alpha*E_Cystoseira/Cys_Pmax) # [uMol O2/kg DW/s]
    #
    # OC_q / 1000 / 1000 * 60*60*24 == OC_q*0.0864
    #
    # Calulate Gross photosynthesis in mmolC/d with DW in grams
    Cys_grossPS     = Cys_grossPR * CYS_SDW / OC_q * 0.0864
    
    #     Carbon exudation fraction  (surplus of C reserve is excreted constantly)
    #     calculated is a fraction of phosotynthesis that is immediately exudated.
    # formulation and parameter following Broch & Slagstad
    
   # Cys_ExudFrac = 1.0 - exp(-cys_ExhFmax * max(0.0, Cys_QC - cys_QCmin))
    
    Cys_ExudFrac = 1.0 - exp(cys_ExhFmax * min(0.0, cys_QCmin - Cys_QC))
    
    Cys_Exudation         =  Cys_grossPS * Cys_ExudFrac
    
    Cys_Photosynthesis    =  Cys_grossPS - Cys_Exudation
    
    #---------------------------------------------------------------------
    #  Nutrient uptake (Nutrients from environment to Nutrient reserve)
    #---------------------------------------------------------------------
    #     limited by external DIN and internal quota, dependent on current speed 
    #     if negative: excretion of N  
    
    # effect of reserveN:structureC quota (droop kinetics) on N uptake
    # described as N reserve saturation (between 0 and 1)
    Cys_RN_sat = min(1.0, (cys_QNmax-Cys_QN) / (cys_QNmax-cys_QNmin)) #[-]
    Cys_RN_sat = max(0,Cys_RN_sat)
    
    # saturation of uptake based on seawater dissolved nutrients (between 0-1)
    Cys_DIN_sat   = DIN/(DIN+cys_ksDIN) # [-]
    
    # functional dependency of nutrient uptake on current velocity (between 0-1)
    Cys_NU = 1.0 - exp(-CurrentV/cys_u065)  #[-]           
    
    # Specific nitrogen uptake rate and total nitrogen uptake
    Cys_Nuptake_rate = cys_MaxNupt * Cys_RN_sat * Cys_DIN_sat * Cys_NU  #[mmolN/gSWD/d]          
    Cys_Nuptake      = Cys_Nuptake_rate * CYS_SDW *2              #[mmolN/d] 
    
    #----------------------------------------------------------
    # STRUCTURAL GROWTH (Based on Reserve pools) & temperature
    #----------------------------------------------------------
    # Reserves get catabolised and assimilated into structural tissue 
    
    # functional dependency of growth on reserve:structure quota
    # the lowest current quotum limits the growth (C or N storage) via droop kinetics
    
    #    fdGrowthRN_SC = min(1.0, max(0.0, 1.0 - minRN_SC/RN_SC))  # for output only
    #    fdGrowthRC_SC = min(1.0, max(0.0, 1.0 - minRC_SC/RC_SC))
    
    Cys_quota_growth_lim  = (1.0 - max(cys_QNmin/Cys_QN, cys_QCmin/Cys_QC))
    Cys_quota_growth_lim  = max(0.0, Cys_quota_growth_lim) 
    Cys_quota_growth_lim  = min(1.0, Cys_quota_growth_lim)
    
    # effect of daylength (effect of photo period)
    Cys_season_effect = cys_PPP1 * (1+ sign(normDL * abs(normDL)^0.5))+cys_PPP2 
    
    
    # specific growth rate and total growth
    Cys_Growth_rate = cys_GRmax *#  Cys_season_effect *
      Cys_grow_TF * Cys_quota_growth_lim         #[/d]        
    
    Cys_Growth      = Cys_Growth_rate *CYS_SDW * (100-tot_coverage)/100    #[gSDW/m2/d]
    
    #====================================================#
    # LOSS TERMS (Both losses in reserves and structure) #
    #====================================================#
    
    #--------------------------
    # Respiration
    #--------------------------
    
    # functional dependency of respiration on temperature
    # Temperature S curve for respiration  S curve f(T) = 2/(1+exp(-T1*(T-TRDr)))
    Cys_TF_Resp    = 2/(1+exp(-cys_Tr1*(TempC-cys_TRDr)))    #[-]
    
    #total mmolC needed for respiration
    Cys_Resp    = Cys_TF_Resp * cys_RD * CYS_SDW / OC_q*0.0864 #[mmolC/d] respired
    
    # respiration is paid by structure if there are not enough reserves (here: regression)
    Cys_pReg  = 1.0 - Cys_QC^3/(Cys_QC^3+cys_QCmin^3) 
    Cys_pReg = Cys_pReg*2 # goes to 1 when C reserve are at minimum.
    Cys_pReg = min(1,Cys_pReg) 
    Cys_pReg = max(0,Cys_pReg)
    
    Cys_Regression = Cys_pReg * Cys_Resp# [mmolC/d] paid by structure, for required respiration]
    Cys_ResResp = Cys_Resp*(1-Cys_pReg) # repsiration paid by C reserves [mmolC/d]
    
    #     activity respiration and total respiration 
    Cys_Growth_Resp = cys_growthRespFrac*(Cys_Growth/cys_DWC) # mmolC/d paid for struct growth
    Cys_Total_Resp = Cys_ResResp + Cys_Growth_Resp                               
    
    #------------------------------------------------
    # Erosion (erosion affects structure & reserve)
    #------------------------------------------------
    
    ##   Erosion_rate = max_Erosion_Rate *1e-6*exp(Erosion_surf_incr*SurfaceArea)/                   
    ##      (1.+1e-6*(exp(Erosion_surf_incr*SurfaceArea)-1.)) # erosion [mmolC/mmolC/d/m2]  
    ##    Erosion      = Erosion_rate*STRUCTURE      #[mmolC/(ind or m2?)/d]
    
    #----------------------
    # Mortality & Necrosis
    #----------------------
    
    Cys_Necrosis = cys_Nec_rate * Cys_TF_Resp * (TempC > cys_T_optH) # dry weight specific necrosis mort
    Cys_Mortality = (Cys_Necrosis+cys_MR) * CYS_SDW # total mortality in gdw/d
    
    #============
    # Quantities
    #============
    
    #     total dry weight of sporophyte, [mg]
##   DryWeight = gdw_SC*STRUCTURE_C + gdw_RN*RESERVE_N + gdw_RC*RESERVE_C
    
    #     total wet weight, [mg]          
##    WetWeight = DryWeight/gdw_gww
    
    #     check on mass balance
##    sumN        = DIN*depth + RESERVE_N + SN_SC*STRUCTURE_C      
    
    #========================================
    # Total exchange rates with environment
    #=======================================
    
    # totals affecting dissolved Nitrogen per m2
    Total_Nupt = Phy_Nuptake + Cys_Nuptake # [mmolN/m2/d]
    Total_N_return = (Phy_Mortality/phy_DWC*SNCr) + (Cys_Mortality/cys_DWC*SNCr) # [mmolN/m2/d mortality]
    N_env_exchage = dilution*(extDIN - DIN)
    
    # totals affecting dissolved oxygen per m2
    Total_Production = (Phy_grossPS + Cys_grossPS)*OC_q # [mmolO2/m2/d produced]
    Total_Resp = (Phy_Total_Resp+Phy_Regression + Phy_Total_Resp+Phy_Regression)*OC_q # [mmolO2/d respired]
    OX_env_exchage = dilution*(extOX - OXYGEN)
    
    # Macroalgae totals
    total_WW = WW_phy+WW_cys
    rel_Phyllophora_cov =  WW_phy/total_WW*100
    rel_Cystoseira_cov = WW_cys/total_WW*100
    rel_seaweed_coverage = total_WW/350 # 350 gww is 100% density in some seaweed species
    
    # N uptake efficiency
    Phy_Nupt_Eff= Phy_RN_sat * Phy_DIN_sat * Phy_NU
    Cys_Nupt_Eff= Cys_RN_sat * Cys_DIN_sat * Cys_NU
    
    #===========================
    # Mass Balance Equations  #no erosion
    #===========================
    
    ## --- Phyllophora ---
    dPHY_SDW   = Phy_Growth - Phy_Mortality - (Phy_Regression*phy_DWC) # [gSDW/m2/d]
    
    dPHY_RESERVE_C   = Phy_Photosynthesis - (Phy_Growth/phy_DWC) - Phy_Total_Resp # [mmolC/m2/d]
    
    dPHY_RESERVE_N   = Phy_Nuptake - (Phy_Growth/phy_DWC*SNCr)  # [mmolN/m2/d]
    
    ## --- Cystoseirsa ---
    dCYS_SDW   = Cys_Growth - Cys_Mortality - (Cys_Regression*cys_DWC) # [gSDW/m2/d]
    
    dCYS_RESERVE_C   = Cys_Photosynthesis - (Cys_Growth/cys_DWC) - Cys_Total_Resp # [mmolC/m2/d]
    
    dCYS_RESERVE_N   = Cys_Nuptake - (Cys_Growth/cys_DWC*SNCr)  # [mmolN/m2/d]
    
    ## --- Shared environment --- # dead biomass N is returned to environment directly
    ## recalculate rates per unit surface to unit volume seawater
    dDIN  = -Total_Nupt/depth + Total_N_return/depth  + N_env_exchage
    
    dOXYGEN = Total_Production/depth - Total_Resp/depth + OX_env_exchage
    
    
    
    
    list(c(dPHY_SDW, dPHY_RESERVE_C, dPHY_RESERVE_N,
           dCYS_SDW, dCYS_RESERVE_C, dCYS_RESERVE_N,
           dDIN,dOXYGEN),
         Phy_QC = Phy_QC, Phy_QN = Phy_QN, 
         Cys_QC = Cys_QC, Cys_QN = Cys_QN, 
         # production & photosynthesis outputs mmolC
         Phy_grossPS =Phy_grossPS, Cys_grossPS = Cys_grossPS,
         Phy_grossPR =Phy_grossPR , Cys_grossPR = Cys_grossPR, 
         Phy_Regression_mmC =Phy_Regression, Cys_Regression_mmC = Cys_Regression,
         Phy_Exudation = Phy_Exudation,  Cys_Exudation = Cys_Exudation,
         Phy_ExudFrac = Phy_ExudFrac, Cys_ExudFrac = Cys_ExudFrac, 
         Total_Production = Total_Production,
         
         # resp outputs
         Phy_ResResp = Phy_ResResp, Cys_ResResp = Cys_ResResp, 
         Phy_Growth_Resp = Phy_Growth_Resp, Cys_Growth_Resp = Cys_Growth_Resp, 
         Phy_Total_Resp = Phy_Total_Resp, Cys_Total_Resp = Cys_Total_Resp,
         Total_Resp = Total_Resp,
         # N dynamics
         Phy_RN_sat = Phy_RN_sat ,  Cys_RN_sat = Cys_RN_sat,
         Phy_Nuptake = Phy_Nuptake,  Cys_Nuptake = Cys_Nuptake,
         Phy_DIN_sat = Phy_DIN_sat, Cys_DIN_sat = Cys_DIN_sat, 
         Phy_NU = Phy_NU, Cys_NU = Cys_NU,
         Phy_Nupt_Eff = Phy_Nupt_Eff, Cys_Nupt_Eff = Cys_Nupt_Eff,
         # growth
         Phy_Growth = Phy_Growth, Cys_Growth = Cys_Growth,
         Phy_quota_growth_lim = Phy_quota_growth_lim,
         Cys_quota_growth_lim = Cys_quota_growth_lim,
         # mortality & regression (in gdw/day)
         Phy_Mortality = Phy_Mortality, Cys_Mortality = Cys_Mortality,
         Phy_Necrosis = Phy_Necrosis, Cys_Necrosis = Cys_Necrosis, 
         Phy_Regression = (Phy_Regression*phy_DWC), Cys_Regression = (Cys_Regression*cys_DWC),
         Phy_pReg = Phy_pReg, Cys_pReg = Cys_pReg,
         # LIGHT 
         SurfPAR_E = SurfPAR_E,
         E_short = E_short, E_mid = E_mid,E_long = E_long,
         E_Phyllophora = E_Phyllophora, 
         E_Cystoseira = E_Cystoseira,
         # WW  
         WW_phy = WW_phy, WW_cys = WW_cys, total_WW = total_WW,
         tot_coverage = tot_coverage,
         # output variables
         rel_Phyllophora_cov = rel_Phyllophora_cov,
         rel_Cystoseira_cov = rel_Cystoseira_cov,
         rel_seaweed_coverage = rel_seaweed_coverage
         
         
         )
  })
}


runtime = seq(from = 1, to = (20*365), by = 1)

# fluctuating forcings
fTemp = function(t) {(5 * sin((2*pi*(t-120)/365)) + 12)}
fPHYTO = function(t) {(5.5 * sin((2*pi*(t-120)/365)) + 6)}
fLight = function(t) {(100 * sin((2*pi*(t-120)/365)) +150)}
#fDIN    <- function(t) {(1 * sin((2*pi*(t+50)/365)) +1.1)}  
fDIN <- function(t) 10
fCurr   <- function(t) 0.1 # net currentspeed at bottom layer at timestep 
fOX     <- function(t) 200 # external dissolved oxygen concentration mmolO2/m3

#           plot(x = 1:365, y = fDIN(1:365))


system.time(
  mod <- ode(y=StatVars, func=MacroAlg, times=runtime, parms=ALG_pars,
             fTemp = fTemp, fCurr = fCurr, fLight = fLight, fPHYTO = fPHYTO, fDIN = fDIN,fOX = fOX)
)
mod
plot(mod, which = c("PHY_SDW","PHY_RESERVE_C","PHY_RESERVE_N",
                    "CYS_SDW","CYS_RESERVE_C","CYS_RESERVE_N",
                    "DIN","OXYGEN"))


# photosynthesis
plot(mod, which = c("Total_Production","Phy_grossPS","Cys_grossPS",
                    "Phy_grossPR","Cys_grossPR"))
# light attenuation & sucepetibility
plot(mod, which = c("SurfPAR_E","E_short","E_mid","E_long",
                    "E_Phyllophora","E_Cystoseira"))

# respiration
plot(mod, which = c("Total_Resp","Phy_Total_Resp","Cys_Total_Resp",
                    "Phy_Regression_mmC","Cys_Regression_mmC"))

# phyllophora carbon dynamics 
plot(mod, which = c("PHY_SDW","PHY_RESERVE_C","Phy_QC",
                    "Phy_grossPS","Phy_ExudFrac","Phy_ResResp",
                    "Phy_Growth_Resp","Phy_Growth","Phy_quota_growth_lim"))

# phyllophora N dynamics 
plot(mod, which = c("DIN","PHY_RESERVE_N","Phy_QN",
                    "Phy_RN_sat","Phy_quota_growth_lim","Phy_DIN_sat",
                    "Phy_NU","Phy_Nupt_Eff","Phy_Nuptake"))

# Cyst carbon dynamics 
plot(mod, which = c("CYS_SDW","CYS_RESERVE_C","Cys_QC",
                    "Cys_grossPS","Cys_Exudation","Cys_ResResp",
                    "Cys_Growth_Resp","Cys_Growth","Cys_quota_growth_lim"))

# Cyst N dynamics 
plot(mod, which = c("CYS_SDW","CYS_RESERVE_N","Cys_QN",
                    "Cys_RN_sat","Cys_quota_growth_lim","Cys_ResResp",
                    "Cys_Growth_Resp","Cys_Growth","Cys_Nuptake"))

# biomass growth & mortality
plot(mod, which = c("Phy_Growth","Cys_Growth",
                    "Phy_Mortality","Cys_Mortality",
                    "Phy_Necrosis","Cys_Necrosis",
                    "Phy_Regression","Cys_Regression"))

depths = c(2,5,10,20,30,40)
runtime_depths = seq(from = 1, to = (10*365), by = 1)



for(i in depths){
  P <- ALG_pars
  P[["depth"]] <- i
  fPHYTO = function(t) {(5 * sin((2*pi*(t-120)/365)) + 6)+4}
  
  assign(x = paste("mod_d_",i,sep = ""),
         value = vode(y=StatVars, func=MacroAlg, times=runtime_depths, parms=P,
                     fTemp = fTemp, fCurr = fCurr, fLight = fLight, 
                     fPHYTO = fPHYTO, fDIN = fDIN,fOX = fOX),pos = -1)

  
  
}

plot(mod_d_2, mod_d_5, mod_d_10, mod_d_20,mod_d_30,mod_d_40, 
     which = c("PHY_SDW","PHY_RESERVE_C","PHY_RESERVE_N",
               "CYS_SDW","CYS_RESERVE_C","PHY_RESERVE_N",
               "DIN","OXYGEN"), lwd = 2)

plot(mod_d_2, mod_d_5, mod_d_10, mod_d_20,mod_d_30,mod_d_40, 
     which = c("Phy_QC","Phy_QN","Cys_QC","Cys_QN",
               "Phy_Growth","Cys_Growth","Phy_quota_growth_lim",
               "Cys_quota_growth_lim"), lwd = 2)

plot(mod_d_2, mod_d_5, mod_d_10, mod_d_20,mod_d_30,mod_d_40, 
     which = c("Phy_pReg","Cys_pReg", "Phy_RN_sat","Cys_RN_sat",
               "Phy_QC","Phy_QN","Cys_QC","Cys_QN"), lwd = 2)

plot(mod_d_2, mod_d_5, mod_d_10, mod_d_20,mod_d_30,mod_d_40,
     which = c("WW_phy","WW_cys","rel_Phyllophora_cov",
                    "rel_Cystoseira_cov","tot_coverage",
               "Phy_Nuptake","Cys_Nuptake"), lwd = 2)
plot(c(1,10),c(1,10))
legend("center", col = c(1:6), legend = c("2m","5m","10m","20m","30m","40m"),
       lwd = 2, lty = c(1:6), cex = 0.5)



## PLOT LIGHT ATTENNUATION GRAPH
LIGHT_PLOT <- function(pars, depthrange, phyto, I0){
  # get light parameters (attenuation coefficients) from set  
  k0Clear <- pars[["k0Clear"]]
  k0_turb_fact <- pars[["k0_turb_fact"]]
  rk_short <- pars[["rk_short"]]
  rk_mid <- pars[["rk_mid"]]
  rk_long <- pars[["rk_long"]]
  # plankton light extinction pars
  k_Phyto_short <- pars[["k_Phyto_short"]]
  k_Phyto_mid <- pars[["k_Phyto_mid"]]
  k_Phyto_long <- pars[["k_Phyto_long"]]
  
  # apply functions from model 
  
  SurfPAR_Wm2 = I0 * 0.43 # surface total irradiance to surface PAR (43 % of light is par)
  PAR_E = SurfPAR_Wm2 * 4.57 # typical conversion factor for sunlight W/m2 to umol photons in PAR range
  # self shading by alga canopy is skipped

  # PAR will be divided into 3 bands of wavelengths short: 400-500 mid:500-600 and long 600-700
  # assume that at sea surface, these 3 bands are of equal photon flux
  E_400_500 = PAR_E*(1/3) # blue light band of wavelengths
  E_500_600 = PAR_E*(1/3) # green light band of wavelengths
  E_600_700 = PAR_E*(1/3) # red light band of wavelengths
  
  #
  k0short = k0Clear * k0_turb_fact * rk_short
  k0mid   = k0Clear * k0_turb_fact * rk_mid
  k0long  = k0Clear * k0_turb_fact * rk_long
  
  ## Calculate light intensity at depth points
  if(length(phyto) == 1){
  
  E_short = E_400_500*(exp(-(k0short+(k_Phyto_short*phyto))*depthrange))
  E_mid = E_500_600*(exp(-(k0mid+(k_Phyto_mid*phyto))*depthrange))
  E_long = E_600_700*(exp(-(k0long+(k_Phyto_long*phyto))*depthrange))
  
  # plot light extinction curves 
  plot(x = depthrange, y = E_short, type = "l", col = "blue", lwd = 2,
       ylim = c(0,I0), xlim = c(min(depthrange),max(depthrange)),
       ylab = "E [mol Photons/m2/s]", xlab = "depth")
  lines(x = depthrange, y = E_mid, col = "green", lwd = 2)
  lines(x = depthrange, y = E_long, col = "red", lwd = 2)
  legend("topright", col = c("blue","green","red"), lwd = 2, 
         legend = c("400-500nm","500-600nm","600-700nm"))
  }
  
  if(length(phyto) == 4){
    par(mfrow = c(2,2))
    for(i in phyto){
    
      E_short = E_400_500*(exp(-(k0short+(k_Phyto_short*i))*depthrange))
      E_mid = E_500_600*(exp(-(k0mid+(k_Phyto_mid*i))*depthrange))
      E_long = E_600_700*(exp(-(k0long+(k_Phyto_long*i))*depthrange))
    
      # plot light extinction curves 
      plot(x = depthrange, y = E_short, type = "l", col = "blue", lwd = 2,
           ylim = c(0,I0), xlim = c(min(depthrange),max(depthrange)),
           ylab = "E [mol Photons/m2/s]", xlab = "depth", main = paste("phyto = ",i))
      lines(x = depthrange, y = E_mid, col = "green", lwd = 2)
      lines(x = depthrange, y = E_long, col = "red", lwd = 2)
      legend("topright", col = c("blue","green","red"), lwd = 2, 
            legend = c("400-500nm","500-600nm","600-700nm"), cex = 0.6)
    }
    
    
    
  }
  
}
  
LIGHT_PLOT(pars = ALG_pars, depthrange = seq(0,40), phyto = c(0,1,10,20),I0 = 200) 

## plot for macroalgae coverage 
COVERAGE_PLOT = function(pars = ALG_pars, unit = "WW"){
  
  if(!(unit %in% c("WW","SDW","mmolC"))){stop("unit must be 'WW','SDW' or 'mmolC'")}
  
  
#  if(unit == "mmolC"){WW = }
  
  par(mfrow = c(2,2))
  
  WW = seq(0,1000)
  # phyllophora coverage (just cover, assume ground cover == canopy cover)
  
  phy_DWC = pars[["phy_DWC"]]
  phy_DW_WW = pars[["phy_DW_WW"]]
  phy_COV_WWa  = pars[["phy_COV_WWa"]]
  phy_COV_WWexp = pars[["phy_COV_WWexp"]]
  
  Phy_coverage = phy_COV_WWa * WW^phy_COV_WWexp
  
  plot(x = WW, y = Phy_coverage, type = "l", main = paste("coverage for Phyllophora in ", unit))
  
  
  
  cys_DWC = pars[["cys_DWC"]]
  cys_DW_WW = pars[["cys_DW_WW"]]
  cys_WW_COVa  = pars[["cys_WW_COVa"]]
  cys_WW_COVb = pars[["cys_WW_COVb"]]
  cys_WW_COVc  = pars[["cys_WW_COVc"]]
  cys_WW_CanoDa = pars[["cys_WW_CanoDa"]] 
  cys_WW_CanoDexp = pars[["cys_WW_CanoDexp"]]
  

  Cys_coverage = cys_WW_COVa-cys_WW_COVb*cys_WW_COVc^WW 
  
  Cys_canopy_dens = cys_WW_CanoDa * WW^cys_WW_CanoDexp
  
  Canopy_Shading = (Cys_canopy_dens+Phy_coverage-100)/(Cys_canopy_dens+Phy_coverage)
  Canopy_Shading[Canopy_Shading<0]<-0
  

  plot(x = WW, y = Cys_coverage, type = "l", 
       main = paste("ground coverage for Cystoseira in ", unit))
  plot(x = WW, y = Cys_canopy_dens, type = "l",
       main = paste("canopy density for Cystoseira in ", unit))
  plot(x = WW, y = Canopy_Shading, type = "l",
       main = paste("Macroalage self shading fraction, assuming each species ", unit))
  
}

COVERAGE_PLOT(pars = ALG_pars)
  