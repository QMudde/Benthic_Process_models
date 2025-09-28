## R model of Black Sea macroalgae biomass

# ## simplified model for macroalgal growth ## #

# The model is parametrized for the red seeweads from the phyllophoraceae familiy
#    # Phyllophora crispa (prev Ph. nervosa) is the target species
#    # Phyllophora truncata (Coccotylus trucatus) data is also used

#------------------------------#
#  MODEL STRUCTURE DESCRIPTION #
#------------------------------#
# The Macroalgae model is based on reserves and structure algae per m2 sediment
# the main variables of the algae are C, N & P reserves and structural biomass.
# Dynamics are driven by fluxes from environmental C, N & P into reserve compartments
# and structural growth & respiration based on usage of those reserves

# MODEL ASSUMPTIONS : 
# 1) state variables have fixed stoichiometry (DEB - strong homeostasis)
# 2) Macroalgae blades are very flat (single cell to few cells thick) 
#         therefore a direct proportion of STRUCTURE to photosynthetic surface is assumed. 
# 3) Environmental C is unlimited, C uptake limited by light & other factors

# Model specifics: 
# 1) this model specifically calculates irradiance at plant level for 3 distinct wavelength ranges,
# this way the attenuation and sensitivity for different wavelenghts of light is taken into account.
# 2) Photosynthesis includes photoinhibition, which can cause a very specifc depth range under certain
# attenuation circumstances 
# 3) self shading of the macroalgae biomass is taken into account, via parameters shading of competitors
# can be included in the model. 
# 4) this model is based on 3 external nutrient species (NO3, NH4 & PO4), these nutrients are implemented 
# as state variables but are added by an external environment via a dilution factor (in the parameters).
# 5) all biomass negative processes are: exudation, respiration, erosion and mortality



# ====================================================
# Phyllophora parameters (mostly phyllophora crispa)
# ====================================================
library(deSolve)

pars_Phy <- c( 
  
  # ------------------------------------------------ # 
  # Quota and ratios:                                #
  # ------------------------------------------------ # 
  
  pC_SDW   = 0.3,  # [gSC/gSDW] mass fraction C of Structural dry weight
  RC_SCmin = 0.05, # [molRC/molSC] minimum molar reserve C over structural C ratio
  
  # N quota and ratios
  SN_SCr   = 0.0355,    # [molSN/molSC] molar N:C ratio of structural tissue
  RN_SCmin = 0.0355,    # [molSN/molSC] minimal reserve N:structural C ratio  
  RN_SCmax = 0.1285625, # [molSN/molSC] maximal reserve N:structural C ratio 
  
  # P quota and ratios
  SP_SCr   = 0.0011,    # [molSP/molSC] molar P:C ratio of structural tissue
  RP_SCmin = 0.0011,    # [molSP/molSC] minimal reserve P:structural C ratio  
  RP_SCmax = 0.0039,    # [molSP/molSC] maximal reserve P:structural C ratio  
  
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
  # (mean must be 1 - sum must be 3)
  rho_short = 0.57, #  [-] rel affinity for short wave (blue), 400-500nm 
  rho_mid   = 1.8,  #  [-] rel affinity for midwave (green), 500-600nm
  rho_long  = 0.73, #  [-] rel affinity for longwave (red), 600-700nm
  
  # Temperature effect on photosynthesis & growth - arrhenius parameters
  t_Arh  = -449.7342,  # [dgK] Arrhenius factor
  t_A_L  = 20640.9346, # [dgK] Arrhenius T factor for lower bound
  t_A_H  = 43549.6620, # [dgK] Arrhenius T factor for lower bound
  t_LB   = 275.9967,   # [dgK] lower bound of optimal temperature range
  t_UB   = 306.0431,   # [dgK] upper bound of optimal temperature range
  t_ref  = 288.15,     # [dgK] reference temperature for photosynthesis
  
  ExhFmax = 0.5,       # [molC/molC] maximal fraction photosynthesis exudated 
  
  # self shading 
  
  h_ALG   = 0.5,       # [m] height of Phyllophora thallus
  a_ALG = 0.00030672,  # [m2/mmolC] (total) Carbon specific shading coefficient
  
  # Impact of competing algae on light
  ww_ALG2 = 0,        # [gWW]    wet weight of competing algae
  h_ALG2  = 1,        # [m]      height of competing algae
  a_ALG2  = 0.002556, # [m2/gWW] wet weight specific shading coefficient
  
  # ----------------------------------------- # 
  # Physiological rates: Growth & Respiration #
  # ----------------------------------------- #
  
  gRmax      = 0.0276,     # [/d] maximal specific growth of structural biomass 
  
  basalResp  = 0.01655588, # [molC/molSC/d] dark repsiration rate 
  growthResp = 0.7,        # [-] fraction of C respired for structural growth 
  
  # @ QUINTEN: use normal Q10 formulation instead (tRChange = log(Q10)/10)
  tRChange = 0.05, # speed of resp change with temperature
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
  
  mort_rate   = 0.008,   # [/d] specific mortality rate
  nec_rate    = 0,       # [/d] heat related necrosis rate
  ero_max     = 1,       # [/d] max erosion rate per structural biomass 
  ero_inc     = 0.022/3, # [/(mmolC/m2)] increase of erosion as surface area increases
  
  #==========================#
  # ENVIRONMENTAL PARAMETERS #
  #==========================#
  
  # --------------------------------- # 
  # Light attenuation coefficients    # 
  # --------------------------------- # 
  
  # @ quinten:  koclear*rk_short en koclear*rk_mid => k0_short, k0_mid, ...
  #  k0Clear = 0.0384, # [/m] PAR range attenuation for clear water (Lorenzen 1972) [-m]
  k0_short = 0.3*0.0384,   # [/m] light extinction rate of 400-500nm [relative to general PAR]
  k0_mid   = 0.45*0.0384,  # [/m] light extinction rate of 500-600nm [relative to general PAR]
  k0_long  = 2.25*0.0384,  # [/m] light extinction rate of 600-700nm [relative to general PAR]
  
  k0_turb   = 0.1,   # [/m] PAR range attenuation for turbid water
  
  k_Phyto_short = 0.0249,      # [/(mg Chl/ m2)] extra light extinction of Chl on short waves
  k_Phyto_mid   = 0.009015,    # [/(mg Chl/ m2)] extra light extinction of Chl on mid waves
  k_Phyto_long  = 0.008642857, # [/(mg Chl/ m2)] extra light extinction of Chl on long waves
  
  # additional information about environment
  latitude = 45,   # [dgN] determines daylight 
  depth    = 10,   # [m] determines light that reaches bottom
  dilution = 0.25  # [/d] dilution of DIN with environemtn/day (0.25 used in Hadley et al., 2015)
)

StatVars <- c(
  # MACROPHYTES structure & reserves
  
  # [mmolC/m2] (structural) Biomass of macroalgae, determining most rates. 
  STRUCT_C = 50, 
  # [mmolC/m2]  reserve C, [energy reserves] 
  RESERVE_C = 50*pars_Phy[["RC_SCmin"]]*2, # start at 4x min
  # [mmolN/m2]  reserve N, [nutrient reserves] 
  RESERVE_N = 50*pars_Phy[["RN_SCmin"]]*1.5,  #start at 1.5x min
  # reserve P [Nutrient reserves]
  RESERVE_P = 50*pars_Phy[["RP_SCmin"]]*1.5, # start at 1.5x min
  
  # ---ENVIRONMENT of 0D box model--- #
  
  NO3 = 10,    # [mmolN/m3] connected to environment (can differ factor 100 in BS between regions)
  NH4 = 5,     # [mmolN/m3] timeseries must be same magnetude as NO3
  PO4 = 1,     # [mmolP/m3] timesereis must be lower than NO3
  
  O2 = 200,   # [mmolO2/m3] connected to environment
  
  # --- Closings to check model validity ---# 
  C_CLOSING = 0, # all C that goes in and out of the Box model [mmolC/m2]
  N_CLOSING = 0, # all N that goes in and out of the Box model [mmolN/m2]
  P_CLOSING = 0 # all P that goes in and out of the BOX model [mmolP/m,2]
)

MacroAlg <- function (t, y, parms, 
                      f_Temperature, f_Current, f_Light, f_PHYTO, 
                      f_ExtNO3, f_ExtNH4, f_ExtPO4, f_ExtO2){
  
  with (as.list(c(y, parms)),{
    
    #---------------------------------------------------------------------------
    # forcing values at time t
    #---------------------------------------------------------------------------
    
    Temperature = f_Temperature(t)      # temperature, degrees celsius
    TempK       = Temperature + 273.15  # temperature, degrees kelvin
    Current     = f_Current(t)   # current speed, [m/s]
    Light       = f_Light(t)     # surface irradiation, [W/m2]
    ExtNO3      = f_ExtNO3(t)    # bottom water NO3 concentration [mmolN/m3]
    ExtNH4      = f_ExtNH4(t)    # bottom water NH4 concentration [mmolN/m3]
    ExtPO4      = f_ExtPO4(t)    # bottom water PO4 concentration [mmolP/m3]
    ExtO2       = f_ExtO2(t)     # bottom water O2 concentration [mmolO2/m3]

    PHYTO       = f_PHYTO(t)     # phytoplankton concentration in water column  [mmolC/m3]
    
    #---------------------------------------------------------------#
    # elemental totals, ratios and quota
    #---------------------------------------------------------------#
    
    STRUCT_N = STRUCT_C * SN_SCr         # [mmolN/m2] structural nitrogen 
    STRUCT_P = STRUCT_C * SP_SCr         # [mmolP/m2] structural phosphorous 
    Macroalgae_C  = STRUCT_C + RESERVE_C      # [mmolC/m2] total macroalgae carbon
    Macroalgae_N  = STRUCT_N + RESERVE_N      # [mmolN/m2] total macroalgae nitrogen
    Macroalgae_P  = STRUCT_P + RESERVE_P      # [mmolP/m2] total macroalgae phophorus
     
    # [[ Quota ]] (relative to structural C)
    RC_SCr  = RESERVE_C/(STRUCT_C+1e-8) # [molRC/molSC] Carbon (C) reserve quotum
    RN_SCr  = RESERVE_N/(STRUCT_C+1e-8) # [molRN/molSC] Nitrogen (N) reserve quotum
    RP_SCr  = RESERVE_P/(STRUCT_C+1e-8) # [molRP/molSC] Phosphorus (P) reserve quotum
    
    # dry and wet Weights
    
    DW       = gdw_SC * STRUCT_C + gdw_RC*RESERVE_C + 
               gdw_RN*RESERVE_N  +
               gdw_RP*RESERVE_P            # [gDW/m2] total dry weight
    
    WW       = DW / dw_ww                  # [gWW/m2] wet weight (biomass in most measurements)
    
    #---------------------------------------------------------------#
    # LIGHT AT DEPTH  Akitsu et al., 2015
    #---------------------------------------------------------------#
    
    # Total irradiance (shortwave heatflux) (300-2600nm) is converted to PAR (400-700nm)
    SurfPAR_Wm2 = Light * 0.45 # [W/m2]    
    
    # Irradiance [W/m2] is recalculated to photon flux [uMol photons/m2/s] 
    SurfPAR_E = SurfPAR_Wm2 * 4.57 
    
    # PAR will be divided into 3 bands of wavelengths 
    # short: 400-500 mid:500-600 and long 600-700
    # ASSUME:: that at sea surface, these 3 bands are of equal photon flux
    E_short_surf = SurfPAR_E*(1/3) #  blue light band of wavelengths 400-500
    E_mid_surf   = SurfPAR_E*(1/3) # green light band of wavelengths 500-600
    E_long_surf  = SurfPAR_E*(1/3) #   red ligth band of wavelengths 600-700
    
    ChlA         = PHYTO/2.5       # [mg/m3] chlorophyll A
    
    # Light based on depth, phytoplankton turbidity, all 3 wavebands attenuate.
    # 3 spectra have different background attenuation due to water (e.g. k0_short=k0clear*rk_short)
    # all 3 spectra are equally affected by non-phytoplankton turbidity ASSUMPTION. (k0_turb)
    # spectra are differently attenuated by average Chlorophyll A conc. (k_Phyto_long*ChlA)
    
    K_short = k0_short + k0_turb  + (k_Phyto_short * ChlA)
    K_mid   = k0_mid   + k0_turb  + (k_Phyto_mid   * ChlA)
    K_long  = k0_long  + k0_turb  + (k_Phyto_long  * ChlA)
    
    K_mean = (K_short + K_mid + K_long)/3 # average attenuation factor
    
    # Resulting Attenuation factors apply on light via Lambert Beer law at the provided depth
    
    # Light at middle of plant height [[h_alg*0.5]]. 
    
    E_short_canopy = E_short_surf * exp(-K_short * (depth - (h_ALG*0.5)))
    E_mid_canopy   = E_mid_surf   * exp(-K_mid   * (depth - (h_ALG*0.5)))
    E_long_canopy  = E_long_surf  * exp(-K_long  * (depth - (h_ALG*0.5)))
    
    # Self-shading
    # extinction rate due to the target Macroalga, and the competing algae
    K_MacroAlg    = a_ALG  * Macroalgae_C / (h_ALG + 1e-10) 
    K_Competition = a_ALG2 * ww_ALG2 / h_ALG2 
    
    Shading_MacroAlg   = exp(-K_MacroAlg    * h_ALG*0.5) # calculate shading at half of the height
    
    # competition over difference in height
    Shading_Competition = exp(-K_Competition * pmax(0, (h_ALG2 - (0.5*h_ALG)))) 

    E_short = E_short_canopy * Shading_MacroAlg * Shading_Competition
    E_mid   = E_mid_canopy   * Shading_MacroAlg * Shading_Competition
    E_long  = E_long_canopy  * Shading_MacroAlg * Shading_Competition
    
    # Actual irradiance on photosynthesis based on absorption spectra ---
    E_Alg  = rho_short * E_short + 
             rho_mid   * E_mid   + 
             rho_long  * E_long # [u mol photons/m2/s]

    #---------------------------------------------------------------------------
    # Arrhenius relationship of temperature for photosynthesis and growth
    #---------------------------------------------------------------------------
    
    Arrh_numerator   = exp( (t_Arh/t_ref) - (t_Arh/TempK) )
    Arrh_denominator = 1 + exp( (t_A_L/TempK) - (t_A_L/t_LB )) + 
                           exp( (t_A_H/t_UB ) - (t_A_H/TempK))
    
    TF_gain          = Arrh_numerator/Arrh_denominator
    T_atREF          = 1 / (1 + exp((t_A_L/t_ref) - (t_A_L/t_LB )) + 
                                exp((t_A_H/t_UB ) - (t_A_H/t_ref)))
    
    TF_gain          = TF_gain/T_atREF  
    
    #-----------------------
    # Photosynthesis
    
    PmaxT = pMax * TF_gain # [mmolC/mmolSC/d] max photosynthetis at current temperature
    
    # theoretical PI curve max
    PS    = PmaxT / ( (alpha / (alpha+beta)) * (beta/(alpha+beta)^(alpha/beta))) 
    a     = alpha * E_Alg/PS
    b     = beta  * E_Alg/PS
    
    # [mmolC/mmolSC/d], gross Structural C specific Photo rate at light and T
    GrossPR = PS * (1-exp(-a))*exp(-b) 
    
    Light_Eff = GrossPR / PmaxT # effectiveness of Photosynthesis due to light only
    Photo_Eff = GrossPR / pMax  # effectiveness of photosynthesis due to light and temperature
    
    # Gross photosynthesis 
    GrossPS   = GrossPR * STRUCT_C # [mmolC/m2/d] 
    
    #---------------------------------------------------------------------------
    # Carbon exudation
    #---------------------------------------------------------------------------
    
    ###  Phy_ExudFrac = 1.0 - exp(-phy_ExhFmax * max(0.0, Phy_QC - phy_QCmin))
    
    ExudFrac          = 1.0 - exp(ExhFmax * min(0.0, RC_SCmin - RC_SCr))
    Exudation         = GrossPS * ExudFrac
    Photosynthesis    = GrossPS - Exudation
    
    #---------------------------------------------------------------------------
    #  Nutrient uptake (Nutrients from environment to Nutrient reserves)[N & P]
    #---------------------------------------------------------------------------
    
    # limited by external nutrient concentrations and internal quota,
    # also dependent on current speed 
    # if flux is negative -> excretion of N or P 
    
    # effect of reserveN:structureC quota (droop kinetics) on N uptake
    # N reserve saturation of nutrient uptake and growth (between 0 and 1)
    RN_sat = min(1.0, (RN_SCmax - RN_SCr) / (RN_SCmax - RN_SCmin)) #[-]  
    RN_sat = max(0.0,  RN_sat)
    
    # effect of reserveP:structureC quota (droop kinetics) on P uptake
    # described as P reserve saturation (between 0 and 1)
    RP_sat = min(1.0, (RP_SCmax - RP_SCr) / (RP_SCmax - RP_SCmin)) #[-]
    RP_sat = max(0.0 , RP_sat)
    
    # Monod saturation formulation of uptake
    NO3_sat   = NO3 / (NO3+ksNO3) # [-] Monod limitation term for NO3 uptake
    NH4_sat   = NH4 / (NH4+ksNH4) # [-] Monod limitation term for NH4 uptake
    PO4_sat   = PO4 / (PO4+ksPO4) # [-] Monod limitation term for PO4 uptake
    
    # functional dependency of nutrient uptake on current velocity (between 0-1)
    NU        = 1.0 - exp(-Current/u065)  #[-]   
    
    # Nutrient uptake efficiency (relative to max specific uptake) dimensionless
    NO3_upt_Eff = NO3_sat * RN_sat * NU 
    NH4_upt_Eff = NH4_sat * RN_sat * NU
    PO4_upt_Eff = PO4_sat * RP_sat * NU
    
    # Specific nutrient uptake
    NO3_upt_rate = max_NO3_upt * NO3_upt_Eff # [mmolN/mmolSC/d]
    NH4_upt_rate = max_NH4_upt * NH4_upt_Eff # [mmolN/mmolSC/d]
    PO4_upt_rate = max_PO4_upt * PO4_upt_Eff # [mmolN/mmolSC/d]

    # total nutrient uptake    
    NO3_upt = NO3_upt_rate * STRUCT_C # [mmolN/m2/d] total NO3 uptake
    NH4_upt = NH4_upt_rate * STRUCT_C # [mmolN/m2/d]
    PO4_upt = PO4_upt_rate * STRUCT_C # [mmolP/m2/d]
    
    #----------------------------------------------------------
    # STRUCTURAL GROWTH (Based on Reserve pools) & temperature
    #----------------------------------------------------------
    
    # functional dependency of growth on reserve:structure quota
    # the lowest quotum limits the growth (C or N storage) via droop kinetics
    # C- limitation factor (1 is fully limiting, 0 is non-limiting)
    RC_sat = RC_SCmin/RC_SCr 
    RC_sat = min(1, RC_sat)
    RC_sat = max(0, RC_sat) 

    # acquire limiting quotum (DROOP); force between 0 and 1
    Quota_growth_fac  = (1.0 - max(RN_sat, RC_sat, RP_sat)) 
    Quota_growth_fac  = max(0.0, Quota_growth_fac) 
    Quota_growth_fac  = min(1.0, Quota_growth_fac) 

    # growth efficiency, [-]
    Growth_Eff = TF_gain * Quota_growth_fac
    
    # Relative growth rate
    Growth_rate = gRmax * Growth_Eff  # [/d] 

    # Total growth
    Growth     = Growth_rate * STRUCT_C      # [mmolC/m2/d]
    
    #====================================================#
    # LOSS TERMS (Both losses in reserves and structure) #
    #====================================================#
    
    #----------------------------------------------------------
    # Respiration
    #----------------------------------------------------------
    
    # simple exponential relationship between respiration and temperature
    TF_resp = exp(tRChange * (Temperature - tRRef))
    
    # total mmolC needed for respiration
    Resp    = TF_resp * basalResp * STRUCT_C # [mmolC/m2/d] respired
    
    # respiration is paid by structure if there are not enough reserves (here: regression)
    # part respriration paid by structural C
    pReg = 0
    if (RC_SCr <= RC_SCmin) pReg =  1. - (RC_SCr/RC_SCmin)
    
    # Total respiration that needs to be paid by structure and reserve Carbon
    
    Regression = pReg * Resp   # [mmolC/m2/d] respiration paid by structure C
    Basal_Resp = Resp*(1-pReg) # [mmolC/m2/d] respiration paid by reserve C 
    
    #     activity (growth) respiration and total respiration 
    Growth_Resp = growthResp * Growth     # [mmolC/m2/d] C paid for growth (activity respiration)
    Total_Resp = Basal_Resp + Growth_Resp # [mmolC/m2/d] total respiration                             
    
    
    #----------------------------------------------------------
    # Mortality, Erosion & Necrosis
    #----------------------------------------------------------
    # total mortality acting on plant biomass is a combination of: 
    # 1) environemtnal mortality (depth based, biomass specific)
    # 2) necrosis (temperature based, biomass specific, hard to parametrize, effect of heatwave)
    # 3) erosion  (current speed based, size dependent)
    
    # erosion formulation from Broch & Slagstad, in their model surface area is 
    # directly related to biomass, so we use biomass specific formulation instead.
    
    # in saccharina model, this scales off surface area, which is struct_C*0.0012
    
    Erosion_rate = ero_max * 1e-6 * exp(ero_inc*STRUCT_C)/
            (1. + 1e-6*( exp (ero_inc*STRUCT_C)-1.) ) # [mmolC/mmolSC/d]  C-specific erosion rate
    
    Erosion       = Erosion_rate * STRUCT_C # [mmolC/m2/d] total erosion
    
    Necrosis_rate = nec_rate * TF_resp * (TempK > t_UB) # [/d] Carbon specific necrosis (due to heat)
    Necrosis      = Necrosis_rate * STRUCT_C # [mmolC/m2/d] total necrosis
    
    Mortality_rate = (Necrosis_rate + mort_rate + Erosion_rate)  # [/d] loss rate
    Mortality      = Mortality_rate * STRUCT_C   # [mmolSC/m2/d] Total mortality, structural C

    #========================================
    # Total exchange rates with environment
    #=======================================
    
    # processes affecting dissolved Nitrogen & Phosphorus per m2
    # Nuptake [mmolN/m2/d]
    
    N_return   = (Mortality + Regression)*SN_SCr + 
                 (Mortality * RN_SCr) # [mmolN/m2/d mortality]
    
    NO3_return = N_return*0.5
    NH4_return = N_return*0.5
    
    PO4_return = (Mortality+Regression)*SP_SCr + 
                 (Mortality * RP_SCr)
    
    # processes affecting dissolved oxygen per m2
    O2_Production  = GrossPS*OC_q                 # [mmolO2/m2/d produced]
    O2_Respiration = (Total_Resp+Regression)*OC_q # [mmolO2/m2/d respired]
    
    # exchange with external environment
    NO3_env_exchange = dilution*(ExtNO3 - NO3)
    NH4_env_exchange = dilution*(ExtNH4 - NH4)
    PO4_env_exchange = dilution*(ExtPO4 - PO4)
    OX_env_exchange  = dilution*(ExtO2  - O2)
    
    # Total budgets: using Closing of C and N as a check
    
    Total_C = STRUCT_C             + 
              RESERVE_C            + 
              C_CLOSING                       # [mmolC/m2]
    Total_N = RESERVE_N            +  
              STRUCT_C * SN_SCr    + 
              (NO3 + NH4) * depth  + 
              N_CLOSING*depth                 #[mmolN/m2]
    Total_P = RESERVE_P            + 
              STRUCT_C * SP_SCr    + 
              PO4 * depth          + 
              P_CLOSING * depth               # [mmolP/m2]
    
    #===========================
    # Mass Balance Equations    
    #===========================
    
    ## --- MacroAlgae ---
    dSTRUCT_C    = (Growth                 # [mmolC/m2/d] 
                    - Mortality 
                    - Regression) 
    
    dRESERVE_C   =   (Photosynthesis         # [mmolC/m2/d]
                    - Growth 
                    - Total_Resp 
                    - Mortality*RC_SCr) 
    
    dRESERVE_N   =   (NH4_upt                # [mmolN/m2/d]
                    + NO3_upt 
                    - Growth*SN_SCr 
                    - Mortality * RN_SCr) 
    
    dRESERVE_P   =   (PO4_upt                # [mmolP/m2/d]
                    - Growth*SP_SCr 
                    - Mortality * RP_SCr)  
    
    ## dead biomass N is returned to environment directly
    ## recalculate rates per unit surface to unit volume seawater
    
    dNO3      = (- NO3_upt/depth             # [mmolN/m3/d]
                 + NO3_return/depth  
                 + NO3_env_exchange) 
    
    dNH4      = (- NH4_upt/depth             # [mmolN/m3/d]
                 + NH4_return/depth  
                 + NH4_env_exchange)
    
    dPO4      = (- PO4_upt/depth             # [mmolP/m3/d]
                 + PO4_return/depth  
                 + PO4_env_exchange)
    
    dO2       = (O2_Production/depth 
                 - O2_Respiration/depth 
                 + OX_env_exchange)
    
    # Closing (losses of mass to outside the box model)
    dC_CLOSING  = (  Total_Resp 
                   + Regression  
                   + Mortality 
                   + (Mortality*RC_SCr) 
                   - Photosynthesis)
    
    dN_CLOSING  = ( - NO3_env_exchange 
                    - NH4_env_exchange)
    
    dP_CLOSING  =   - PO4_env_exchange
    
    Net_Growth_rate = dSTRUCT_C/STRUCT_C # [/d]

    list(c(dSTRUCT_C, dRESERVE_C, dRESERVE_N, dRESERVE_P,
           dNO3, dNH4, dPO4, dO2, 
           dC_CLOSING, dN_CLOSING, dP_CLOSING),
         
         # Quota, weights, total nutrient contents
         RC_SCr = RC_SCr,   # [molRC/molSC] Carbon (C) reserve quotum
         RN_SCr = RN_SCr,   # [molRN/molSC] Nitrogen (N) reserve quotum
         RP_SCr = RP_SCr,   # [molRP/molSC] Phosphorus (P) reserve quotum
         
         # output variables & Totals
         WW = WW,           # [gWW/m2] wet weight macroalgae biomass
         DW = DW,           # [gDW/m2] dry weight macroalgae biomass
         
         Macroalgae_C = Macroalgae_C, # [mmolC/m2] total macroalgae carbon
         Macroalgae_N = Macroalgae_N, # [mmolN/m2] total macroalgae nitrogen
         Macroalgae_P = Macroalgae_P, # [mmolP/m2] total macroalgae phophorus
         
         # production & photosynthesis outputs mmolC
         GrossPS        = GrossPS, # [mmolC/m2/d]  total gross photosynthesis
         GrossPR        = GrossPR, # [mmolC/mmolSC/d] gross specific Photosynthesis rate at light and T
         Photosynthesis = Photosynthesis, # [mmolC/m2/d] photosynthesis - exudation
         
         ExudFrac       = ExudFrac, # [-] fraction of gross photosynthesis exudated
         Exudation      = Exudation, # [mmolC/m2/d] total C exudation
         
         # respiration
         Basal_Resp  = Basal_Resp,  # [mmolC/m2/d] basal respiration paid by reserve C 
         Growth_Resp = Growth_Resp, # [mmolC/m2/d] C paid for growth (activity respiration)
         Respiration = Total_Resp,  # [mmolC/m2/d] total respiration (-necrosis)
         
         Respiration_loss = Total_Resp + Regression, # [mmolC/m2/d] total respiration (+necrosis)
         
         # limitation/inhibition terms
         NO3_sat = NO3_sat, # [-] Monod limitation term for NO3 uptake
         NH4_sat = NH4_sat, # [-] Monod limitation term for NH4 uptake
         PO4_sat = PO4_sat, # [-] Monod limitation term for PO4 uptake
         NU      = NU, # [-] dependency of nutrient uptake on current velocity
         
         RC_sat  = RC_sat, # [molRC:molSC] N reserve saturation of nutrient uptake and growth 
         RN_sat  = RN_sat, # [molRN:molSC] C reserve saturation of growth
         RP_sat  = RP_sat, # [molRP:molSC] P reserve saturation of nutrient uptake and growth
         
         Quota_growth_fac = Quota_growth_fac, # [-] limitation of growth due to quota
         
         Photo_Eff = Photo_Eff, # [-] part of max Photosynthesis realised to light only
         Light_Eff = Light_Eff, # [-] part of max Photosynthesis realised due to light and temperature
         
         # Nutrient dynamics
         NO3_upt_Eff = NO3_upt_Eff, # [-] N03 uptake efficiency (relative to max uptake) 
         NH4_upt_Eff = NH4_upt_Eff, # [-] NH4 uptake efficiency (relative to max uptake) 
         PO4_upt_Eff = PO4_upt_Eff, # [-] PO4 uptake efficiency (relative to max uptake) 
         
         # total nutrient uptake    
         NO3_upt = NO3_upt,         # [mmolN/m2/d] total NO3 uptake
         NH4_upt = NH4_upt,         # [mmolN/m2/d] total NH4 uptake
         PO4_upt = PO4_upt,         # [mmolP/m2/d] total PO4 uptake
         
         # growth & limitation
         Growth      = Growth,      # [mmolC/m2/d] total structural growth
         Growth_rate = Growth_rate, # [/d] relative growth rate
         Growth_Eff  = Growth_Eff,  # [-] growth efficiency
         
         Net_Growth_rate  = Net_Growth_rate, # [/d] net structural increase (dSTRUCTC/STRUCTC)
         
         # mortality & regression 
         Mortality       = Mortality,  # [mmolSC/m2/d] total mortality, structural C
         Necrosis        = Necrosis,   # [mmolSC/m2/d] total necrosis, structural C
         Regression      = Regression, # [mmolC/m2/d] respiration paid by structure C
         Erosion         = Erosion,    # [mmolSC/m2/d] total erosition, structural C
         pReg            = pReg,       # [-] part respriration paid by structural C
         Necrosis_rate   = Necrosis_rate,  # [/d] Carbon specific necrosis rate 
         Erosion_rate    = Erosion_rate,   # [/d]  C-specific erosion rate
         Mortality_rate  = Mortality_rate, # [/d]  C-specific mortality rate
         
         # LIGHT 
         SurfPAR_E = SurfPAR_E, # [uMol photons/m2/s] irradiance on water surface
         E_short   = E_short, # [uM/m2/s] blue-light band (400-500nm) on algae
         E_mid     = E_mid,   # [uM/m2/s] green-light band (500-600nm) on algae
         E_long    = E_long,  # [uM/m2/s] red-light band (600-700nm) on algae
         E_Alg     = E_Alg,   # [uM/m2/s] total irradiance for photosynthesis
         K_mean    = K_mean,  # [/m] average attenuation factor for light
         
         K_MacroAlg    = K_MacroAlg,  # [/m] extinction coeff due to self shading
         K_Competition = K_Competition, # [/m] extinction coeff competing algae
         Shading_MacroAlg = Shading_MacroAlg, # [-] shading factor due to modeled algae
         Shading_Competition = Shading_Competition,# [-] shading factor due to competing algae
         
         Shading_Eff = (1-Shading_MacroAlg*Shading_Competition), # [-] part of light not shaded by algae
         
         # Temperature effects 
         TF_gain = TF_gain,  # [-] temperature factor for photosynthesis and growth
         TF_resp  = TF_resp, # [-] temperature factor for respiration and necrosis
         
         Total_C = Total_C, # [mmolC/m2] total carbon in system
         Total_N = Total_N, # [mmolN/m2] total nitrogen in system
         Total_P = Total_P, # [mmolP/m2] total phophorus in system
         
         # (External) forcings
         Temperature = Temperature,  # [degC] temperature
         Current     = Current,      # [m/s] current speed
         
         PHYTO     = PHYTO,     # [mmolC/m3] phytoplankton concentration in water column 
         Light     = Light,     # [W/m2] surface irradiation
         ExtNO3    = ExtNO3,    # [mmolN/m3] bottom water NO3 concentration
         ExtNH4    = ExtNH4,    # [mmolN/m3] bottom water NH4 concentration
         ExtPO4    = ExtPO4,    # [mmolP/m3] bottom water PO4 concentration
         ExtO2     = ExtO2      # [mmolO2/m3] bottom water O2 concentration
    )
  })
}

################################################################################
#  fluctuating forcings (based on Bahmbi outputs)
################################################################################

f_Temperature     <- function(t) {
          (10 * sin((2*pi*(t-120)/365)) + 14)
  }
f_PHYTO    <- function(t) {
          (5.5 * sin((2*pi*(t-120)/365)) + 6)
  }
f_Light    <- function(t) {
          (100 * sin((2*pi*(t-120)/365)) +150)
  } 

f_NO3      <- function(t) {
          (10 * sin((2*pi*(t+30)/365)) +15)
  } 
fNO3_coast <- function(t) {
          (10 * sin((2*pi*(t+50)/365)) +15)
  } 

f_NH4       <- function(t) {
           (10 * sin((2*pi*(t+30)/365)) +15)/2
  } 
fNH4_coast <- function(t) {
           (10 * sin((2*pi*(t-130)/365)) +15)
  } 

# assume PO4 is roughly 1/10 of DIN always
f_PO4       <- function(t) {
           (0.5 * sin((2*pi*(t+80)/365)) +0.6)
  } 
fPO4_coast <- function(t) {
           (fNO3_coast(t) + fNH4_coast(t))/10 
  }

f_Current      <- function(t) 0.1 # net currentspeed at bottom 

f_O2        <- function(t) 200 # dissolved oxygen conc mmolO2/m3

fOX_coast  <- function(t) {
           (160 * sin((2*pi*(t+50)/365)) +200)
} 

################################################################################
# model runs for offshore and coast
################################################################################

# run model for 10 years
runtime <- seq(from = 1, to = (10*365), by = 1)

pars    <- pars_Phy

phyMOD <- ode(y = StatVars, func = MacroAlg, times = runtime, parms = pars,
              f_Temperature = f_Temperature, f_Current = f_Current, 
              f_Light = f_Light, f_PHYTO = f_PHYTO, 
              f_ExtNO3 = f_NO3, f_ExtNH4 = f_NH4, 
              f_ExtPO4 = f_PO4, f_ExtO2 = f_O2)

phyMOD_coast <- ode(y = StatVars, func = MacroAlg, times = runtime, parms = pars,
              f_Temperature = f_Temperature, f_Current = f_Current, 
              f_Light = f_Light, f_PHYTO = f_PHYTO,
              f_ExtNO3 = fNO3_coast, f_ExtNH4 = fNH4_coast, 
              f_ExtPO4 = fPO4_coast, f_ExtO2 = fOX_coast)

plot(phyMOD, phyMOD_coast, 
     which = c("DW", "WW", "STRUCT_C", "RESERVE_C",
               "TF_gain", "TF_resp"))


plot(phyMOD, phyMOD_coast, 
     which = c("RESERVE_C", "RC_SCr", "RC_sat",
               "RESERVE_N", "RN_SCr", "RN_sat",
               "RESERVE_P", "RP_SCr", "RP_sat"))

plot(phyMOD, phyMOD_coast, 
     which = c("Shading_Eff", "Quota_growth_fac", "RC_sat",
               "Growth", "Mortality", "Respiration_loss",
               "RC_SCr", "STRUCT_C", "RESERVE_C"))

################################################################################
# what is the effect of turbidity on depth distribution ?
################################################################################

# test 1, high turbidity & HIGH PHYTO

depths   <- c(   2,    5,   10,  20,  30)

NoCysto  <- c(   0,    0,    0,   0,   0) # this is without competition

# estimate of depth spread of C. crinita (Afanasyev et al., 2017, NE BS)
Cysto    <- c(1500, 1000,  400,   0,   0) 


runtime_depths = seq(from = 1, to = (10*365), by = 1)

for(i in 1:5){
  P              <- pars_Phy
  P[["depth"]]   <- depths[i]
  P[["ww_ALG2"]] <- Cysto[i]
  
  P[["k0_turb"]]  <- 0.1   # high turbidity
  
  fPHYTO_2 = function(t) {(9 * sin((2*pi*(t-120)/365)) + 10)}
  
  assign(x = paste("mod_d_", depths[i], sep = ""),
         value = ode(y = StatVars, func = MacroAlg, 
                     times = runtime_depths, parms = P,
                     f_Temperature = f_Temperature, 
                     f_Current = f_Current, 
                     f_Light = f_Light,  f_PHYTO = fPHYTO_2 ,
                     f_ExtNO3 = f_NO3, f_ExtNH4 = f_NH4,
                     f_ExtPO4 = f_PO4, f_ExtO2 = f_O2), pos = -1)
}

plot(mod_d_2, mod_d_5, mod_d_10, mod_d_20, mod_d_30, 
     which = c("NO3", "NH4", "PO4", 
               "RN_sat", "RP_sat", "Quota_growth_fac"), lwd = 2)

plot(mod_d_2, mod_d_5, mod_d_10, mod_d_20, mod_d_30, 
     which = c("STRUCT_C", "RESERVE_C", "WW", 
               "Photo_Eff", "Quota_growth_fac"), lwd = 2)

plot(c(1,10),c(1,10))

legend("center", col = c(1 : 6), 
       legend = c("2m","5m", "10m","20m","30m"),
       lwd = 2, lty = c(1:6), cex = 1.2)

## test 2: LOW TURBIDITY & LOW PHYTO

for(i in 1:5){
  P2              <- pars_Phy
  P2[["depth"]]   <- depths[i]
  P2[["ww_ALG2"]] <- Cysto[i]
  
  P2[["k0_turb"]]  <- 0.04   # low turbidity
  
  fPHYTO_3 = function(t) {(2 * sin((2*pi*(t-120)/365)) + 6)+1}

  assign(x = paste("mod_dt_", depths[i], sep = ""),
         value = ode(y = StatVars, func = MacroAlg, 
                     times = runtime_depths, parms = P2,
                     f_Temperature = f_Temperature, f_Current = f_Current, 
                     f_Light = f_Light,  f_PHYTO = fPHYTO_3 ,
                     f_ExtNO3 = f_NO3, f_ExtNH4 = f_NH4,
                     f_ExtPO4 = f_PO4, f_ExtO2 = f_O2), pos = -1)
  
}

plot(mod_dt_2, mod_dt_5, mod_dt_10, mod_dt_20, mod_dt_30, 
     which = c("STRUCT_C", "RESERVE_C", "WW",
               "Photo_Eff", "Quota_growth_fac"), lwd = 2)
plot(c(1,10),c(1,10))
legend("center", col = c(1:6), 
       legend = c("2m", "5m", "10m", "20m", "30m"),
       lwd = 2, lty = c(1:6), cex = 1.2)

plot(mod_dt_2, mod_dt_5, mod_dt_10, mod_dt_20, mod_dt_30,
     which = c("RESERVE_C", "RC_SCr", "RC_sat",
               "RESERVE_N", "RN_SCr", "RN_sat",
               "RESERVE_P", "RP_SCr", "RP_sat"))

# less turbidity broadens the available growth depth range!

## ADDITIONAL NOTES: 
# black sea depths are known to reach  400-800 g WW/m2 (Alexandrov & Milchakova, 2021)
# -----> in MPA's range is 1-400 gWW/m2 (usually a bit lower in MPA areas)
# -----> other areas 5-850 gWW/m2


