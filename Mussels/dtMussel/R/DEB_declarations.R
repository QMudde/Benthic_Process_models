.Deb_cohort <- new.env()

.Deb_cohort$parms <- c(  #passed to Fortran surbroutines, KEEP EXACT SAME NAMES & ORDER
  
    #----------------#
    # DEB PARAMETERS #
    #----------------#
    aeff   = 0.71,    # assimilation efficiency,ESTIMATED      [-] (Sara 2011 0.88)
    kappa  = 0.5,     # fraction energy to structure  [-] Sara etal 2012 
    maint  = 0.015,   # carbon-specific maintenance rate [/d] 
    cat_s  = 0.267,   # surface-specific catabolism rate [mmolC/cm2/d] 
    eg     = 0.407,   # C respired to create structural C [mmolC/mmolC] ==> [-]
    # ==> DEB.Eg     = 1900 (Sara, 2011) 
    # energy content of STRUCT volueme unit Sara2011 --> 1350 Joule -> recalc value eg = (1900/1350)-1 = 0.4074074
    
    #---------------------#
    # Conversion factors  #
    #---------------------#
    # J_mmolC = 525 # standard values from AmP 
    
    # note on weights: there are 5 types of possible weight metrics,
    # this causes great uncertainty in literature search 
    # 1) Wet Weight inlcuding shell (live weight/Total weight), here: TW 
    # 2) Wet weight excluding shell (wet meat weight) here: WW (0.2191 used)
    # 3) Dry weight including shell (dry weight) here: not used 
    # 4) Dry weight excluding shell (dry mean weigth/ often just DW) here: DW
    # 5) Ash free dry weight (only organic compounds/pure biomass) here: AFDW
    
    # gC_gSFDW = 0.4 (smaal & vonck) expressed as % of dry meat weight
    # gAFDW_gSFDW = 0.8 # gram shell free dry weight per gram ash free dry weight
    # 0.75 is based on Palmerini & Bianchi 1994, Saraiva et al. estimated 0.85
    # ifwe  use avg: 0.8, gc_gafdw = 0.4/0.8 = 0.5 
    # gC_gAFDW = 0.5 (similar to base org ratio C:H:O:N = 1:1.8:0.5:0.15) 0.5020921
    # carbon content of AFDW is higher than SFDW because AFDW is pure organics and 
    # SFDW can contain traces of inorganic material
    # Estimation for SFDW to SFWW ratio would be : 
    
    # gAFDW/molC = 23.9 # often used conversion for pure organics
    
    # gdw/molC = 26.66666  #!!  30 with 40% Carbon value (0.09/3*1000)
    
    # to calculate gDW to cm3 Saraiva 2012 & Rosland 2009 used 0.2 gDW_cm3 becuase this 
    # was their conversion factor gDW_gWW and assumed 1gWW == 1cm3 struct. 
    # I do not agree with this assumption. 
    
    gDW_cm3   = 0.09,  # from cm3 to gram dry weight    [gDWT/cm3]
    # 0.09 is often used (monaco & mcquaid M.gallo/ Wijsman M.edulis, AmP)
    
    # gAFDW_cm3 = 0.09*0.8 (0.072)
    # mmolC_cm3 = 0.072*0.5/12*1000 = 3
    #  mmolC_cm3 = 3.375, # from cm3 to mmolC              [mmolC/cm3]  (0.09*0.45/12*1000)
    mmolC_cm3 = 3, # (0.09*0.4/12*1000)
    #  g/mmolC = 0.03
    # 3 * 12 / 1000 / 0.5 -> 0.072 gAFDW_cm3 
    # 0.072/0.8 -> 0.09 gDW_cm3
    # 
    # 
    
    gWW_gDW = 7, # gram WW (no shell) per gram DW (AFDW?!)
    #gWW_gDW = 2.62, # gram WW (no shell) per gram DW (no shell with ash) (1/0.3818281)
    #gWW_gDW = 3.022069, # calcultaed from other ratios
  #  gTW_gAFDW   = 17.24138, # from gAFDW to gram total WW (including shell!) (Palmerini & Bianchi 1994) 
    # average of 2 methods = 0.065 TW:AFDW -> 15.4 gWW/gAFDW
    # (Brigolin 2009) : from gram DW to gram WW  = 7 / gDW to TW  17 

    # gTW_gWW = 4.564126

    
    del.M     = 0.22,# from equivalent spherical length to length [-] (Sara, 2011)
    
    # FOOD parameters (food concentration [FOOD] in molC/m3) [independent of environment]
    # phy_pref  = 0.9,  # preference for algae over detritus [-]
    # sim_inhib = 1,    # half sat constant for SIM [g/m3 // or  mg/liter] (100*????)
    # k0_food   = 2.275,
    
    # energy contents affecting detritus production
    # phy_E     = 525, # [J/mmolC] (add my pet)
    # det_E     = 150, # [J/mmolC]
    # mussel_E  = 525, # [J/mmolC] (add my pet)
    
    # Oxygen limitation for feeding 1 mg/L = 1/32*1000 mmol O2/m3
    # Tang & Rissgard fig 5 M. edulis
    ks_O2    = 62.5,   # mich-menten oxygen limitation [mmol O2/m3]  (KS-crit)
    # (2 mgO2/L tang & Riisgard 2017) (M.edul)
    crit_O2  = 31.25,   # oxygen threshold conc (below this: valve shut/no filt) [mmol O2/m3]
    # (1 mgO2/L)Tang & Riisgard 2017 (fig 5)
    
    # Mortality parameters #
    e_dens_trh    = 0.78, # avg energy density at which pop mortality starts [mmolC R / mmolC S] CALIBRATED
    mort_starv = 0.031, # starvation-induced mortality rate [/d] (0.009) CALIBRATED
    mort_base   = 0.00041, # min population mortality due to environment [/day] (Perez camacho 1991)& Brigolin 2006
    # = 0.15/365 (often used mussel mortality)
    # 0.618 average NWS mortality rate Shurova, N. M., & Gomoiu, M. T. (2006)
    
    # Reproduction & Maturity parameters #
    
    mat_pub     = 3.022857, # [mmolC/ind] maturity (gonads) when reproduction starts
    spawn_peak  = 90,  # [julian day] day at which spawning optimum occurs (~april 1) # ESTIM
    spawn_peak_2 = 280, # [julian day] day at second spawning occurs (~oct 15) # ESTIM
    spawn_dev   = 5,    # [day] deviation (spread) of spawning intensity 
    spawn_eff   = 1,    # [-]  fraction of gametes that is released by spawning (1) version
    weight_egg  = 3e-06, #[mmolC/egg] # 80 ng eggs (Widdows 1991)
    # 9.246825e-07, # [mmolC/egg] recalculated from AmP ww at birth (old version)
    weight_settle   = 2.3625e-05,   # [mmolC/ind] biomass of larvae at which they start cohort
    # 2.3625e-05 = (630/1000000)*0.45/12, #[mmolC/ind]  settle afdw = 630 ng (widdows 1991) --> mmolC
    v_settle        = 1,         # [m/day], base settling rate of settling larvae
    pReserve_settle = 0.5527686,   # fraction of C into reserve for larvae that settle
    # pReserve_settle estimated so that newly settled cohorts start with max reserve density
    # 2190 /(525*3.375)  ==> 1.235979 (ratio RESERVE/STRUCT in new indiv)
    # this is often how DEB variables for new individuals are initialized
    
    # NEWPARS, EFFECT OF CROWDING
    #  substr.type    = 0.5,    # medium substrate affect settling rate of larvae to mussels
    #  crowding_start     = 60, # [cm2/m2] fraction of surface covered when crowding occurs
    #  ks_crowd           = 100, # [cm2/m2] fraction of surface covered when crowding effect = 50%
    #  crowding_mortality = 0.002 # [/d] max per capita crowding mortality per day
   
   # TEMPERATURE PARAMETERS ## (default M.edul parameters + upper/lower boundaries for M.gallo)
   
   T_a  = 5800, # Arrhenius Temperature  [[K]]
   T_r  = 293.15, # reference temperature, Kelvin to room temperature in Celcius [[K]]
   T_l  = 275.15, # lower boundary of tolerance range [[K]] (was 275) 273.15+6.7 (Jansen 2006/2009) ABT
   T_h  = 296.15, # upper boundary of tolerance range [[K]] (was 296) 273.15+27 (Jansen 2009) ABT
   T_al = 45430, # rate of decrease of LOWER boundary [[K]]
   T_ah = 31376,#, # rate of decrease of UPPER boundary [[K]]
   
   # Larvae parameters (only for cohort model)
   mort_larva    = 0.1,  # [/day] per capita density dependent morality rate larvae (here estim)
   growth_larva  = 0.13 ,  # [/day] rate of increase (80 ng to 630 ng (widdows in 16 days (lazo&pita)))
   ks_larva      = 10,      # [mmol C/m3], mich-ment half concentration of phytoplankton on larvae
   max_set_rate = 200, # max daily settlement rate of larva [ind/d] toupoint etal 2012
   ks_settlement = 2000, # competent pediveliger density at which half of max set rate is reached [ind/m3]
   # toupoint etal 2012
  
   # based on ~15k cells/ml (Pettersen et al., 2010) to reach half of max growth (at 60k cells)
   # estimated avg 15pgC/cell (Sprung, 1984b) -> 
   # could also use ~ 10 pg/cell (Llewellyn & Gibb 2000) 
   # (prev) assuming 1.1 pico molC/cell -> 1.1e-09*(15000*1000000) ->16.5 
   #new assuming 1 pico molC/cell (mix of I. galbana & diatoms) & 10kcells/ml at half growth -> 10
   #based on approximated monod curve max/2
   # was 3.75
   
   # ENVIRONMENT OF THE LOCAL WATERBOX (only for cohort model)
   
   water_renew = 1,      # water renewal rate (based on hydrology)  [/day] (simplified version = 1/depth, assuming only sinking of 1m/d?)
   O2_renew    = 1,      # O2 renewal rate  [/day] equal to water renew, can set high to diminish O2 effect
   depth  = 5,         # sediment depth [m] relevant for box volume, this is not equal across black sea bottom forcing boxes
   
   max_clear = 0.15,  # [m3/d/cm2] (Saraiva 2011 M.edul, 0.096)  ESTIMATED 
   max_phyto_filt = 1.038794, #[mmolC/cm2/d] (Rosland 2009 M.edul = 0.48) ESTIMATED (1.3)
   max_DT_filt = 2.000000,   #[mmolC/cm2/d] # ESTIMATED (1.78)
   max_spm_filt = 3.5,   #[g/cm2/d] (Saraiva 2011 M.edul)
   
   ### INGESTION/PSEUDOFAECES PARAMETERS
   rho_Phyto = 0.99,   #[-] (Saraiva 2011 M.edul)
   rho_DT = 0.6,      #[-] # ESTIMATED (prev 0.95/0.86)
   rho_SIM = 0.45,    #[-] (Saraiva 2011 M.edul)
   max.ing.Phy = 45, # [mmolC/d] (Saraiva 2011 M.edul) (prev = 9.2) (ESTIM)
   max.ing.DT = 6.8, # [mmolC/d] max DT carbon uptake # ESTIMATED (prev = 20) 
   max_SIM_ing = 0.47, # [g/d] (Saraiva 2011 M.edul = 0,15) # ESTIMATED
   
   ### ASSIMILATION/FAECEATION PARAMETERS  Parameter E ratio = 525:150 --> 3.5
   nc_tissue = 0.21, # ~~ Saraiva (Saraiva 2011 M.edul) 0.1849711 Smaal vonck 1997
   nc_Phyto = 0.151, # 0.151 ~~ redfield approx (this fluctuates heavily for phytoplankton in reality) # ESTIMATED
   #(avg, but not a constant in real life) # between 0.05 and 2 in BAMBHI model (gregoire 2008)
   nc_DT = 0.05,   # ~ 15 CN ~ Banaru et al., 2007 (black sea sediment pom nc ratio) 0.067 / 0.05 / 0.02 / 0.04 
   # NC_DT  = lowest end of POM range in black sea, since Phyto is excluded
   # max allowed
   ## shell shape and shell cover determining crowding
   shell_allo_a = -0.6456, # shell allomtric function to shell cover (intercept) (Babarro & Carrington 2013)
   shell_allo_b = 0.3978, # shell allometric function to shell cover (slope) (Babarro & Carrington 2013)
   crowd_lim = 10000, # shell cover to intialise crowding mortality ( 10000 = 100% cover) 
   ## substrate factor
   substr_factor = 1, # sediment grain size (cm diameter) (1 -> hard rock substrate)
   mort_juv = 0.03, # mort at smallest mussel size, inital thinning (Perez camacho 1991)
   mort_egg = 0.95, # inital egg mortality (fraction of eggs not survive to larva)
   # saraiva 2014 uses 0.95 for "realistic runs" 
   # Sprung (1984) suggest even higher 0.99 field mortality for eggs
   #for Q10 larva values, can set to 1 if running with bottom values not representing pelagic in 0d?
   q10_LgL = 2.7, # Q10 (lower) on larvae growth BELOW reference temp (Lazo & Pita 2012 -> 2.7)
   q10_Lm = 2,     # Q10 (mirrored) on larvae mortality around reference temp (below = 1/q10_Lm) -> 2
   length_trh = 1, # [cm] threshold above which density gets added to FIELD outputs
   recr_length = 0.5 # size at recruitment (for juvenile mortality)
)

.Deb_cohort$parms <- data.frame(names = names(.Deb_cohort$parms), 
                                default = .Deb_cohort$parms)
row.names(.Deb_cohort$parms) <- NULL

# Description of parameters
.Deb_cohort$parms$description <- c(
     "base assimilation efficiency",  "fraction catabolism energy to structure", 
     "carbon-specific maintenance rate",
     "surface-specific catabolism rate", "amount of C respired to create structural C", 
     "conversion from cm3 to gram dry weight", "conversion from cm3 to mmolC", 
     "conversion from gram Dry weight to gram wet Weight", 
     "shape factor, from equivalent spherical length to structural length", 
     "half-saturation oxygen concentration affecting metabolic rates", 
     "oxygen concentration lower threshold (below this: no more activity)", 
     "Individual average R/S ratio threshold where population starvation mortality starts", 
     "starvation-induced additional per capita mortality", 
     "basal per capita mussel mortality (environmental)",
     "maturity (gonad biomass) threshold, where reproduction (gamete production) starts", 
     "day at which (first) spawning optimum occurs",
     "day at which second spawning happens (no 2nd spawning occurs when this is 0)",
     "deviation (spread) of spawning intensity", 
     "fraction of gametes that is released by spawning (always 1 in current version)",
     "biomass of a single egg when released", 
     "biomass of average larvae at which they settle into new cohort", 
     "sinking rate of settling larvae (default = 1 = not active", 
     "fraction of C allocated into reserve for settling larvae", 
     "Arrhenius Temperature", 
     "reference temperature", 
     "lower boundary of tolerance range", 
     "upper boundary of tolerance range", 
     "rate of decrease at lower boundary",
     "rate of decrease at upper boundary",
     "base mortality of larvae", "optimal growth rate of larvae", 
     "half-saturation constant of feeding of larvae for phytoplankton",
     "max daily settlement rate of larva to a m2 of sediment",
     "competent pediveliger density at which half of max set rate is reached",
     "water renewal rate", "O2 renewal rate (set very high to ignore O2 effect in dynamic runs)", 
     "depth of water column",
     "surface specific maximum clearance rate",
     "maximum surface specific phytoplankton filtration",
     "maximum surface specific detritus filtration",
     "maximum surface specific inorganic matter filtration",
     "ingestion preference factor for phytoplankton",
     "ingestion preference factor for detritus",
     "ingestion preference factor for inorganic matter", 
     "maximum daily phytoplankton ingestion", 
     "maximum daily detritus ingestion",
     "maximum daily inorganic matter ingestion",
     "average N:C ratio mussel tissue",
     "average N:C ratio phytoplankton",
     "average N:C ratio detritus",
     "shell allomtric function to shell cover (intercept)",
     "shell allometric function to shell cover (slope)",
     "population shell surface coverage to trigger crowding mortality",
     "average grain size of substrate (>=1 -> rock (maximum value)",
     "additional juvinile mortality after settlement, represents recruitment mortality",
     "fraction of eggs that dies at spawning event",
     "Q10 factor for larval growth up to 20dgC",
     "Q10 factor for larval mortaltiy (base =20dgC)",
     "Shell length filter for comparing model output with obs",
     "Length for half juvinile mortality (recruitmentlength)") 
  
# Units of parameters
.Deb_cohort$parms$units <- c(
     "-", "-", "/d", "mmolC/cm2/d", "-",
     "gDWT/cm3", "mmolC/cm3", "gram/gram", "-", "mmolO2/m3", "mmolO2/m3",
     "/d", "/d", "/d", "mmolC/ind", "day", "day", "day", "-", "mmolC", 
     "mmolC/ind", "m/day", "-", "dgK", "dgK", "dgK", "dgK", "dgK", "dgK",
     "/d", "/d", "mmolC/m3",
     "ind/d","ind/m3",
     "/d", "/d", "m", "m3/d/cm2",
     "mmolC/cm2/d","mmolC/cm2/d","g/cm2/d",
     "-","-","-", "mmolC/d", "mmolC/d", "g/d",
     "molN/molC","molN/molC","molN/molC",
     "-","-","cm2/m2",
     "cm3","/d","-","-","-",
     "cm","cm")
  
.Deb_cohort$y1D <- data.frame(names =  #names from Fortran surbroutines, KEEP EXACT SAME NAMES & ORDER
    c("RESERVE", "STRUCT", "REPROD", "POPULATION"))
.Deb_cohort$y1D$description <- 
    c("Reserve carbon per individual", "Structural carbon per individual", 
      "Reproductive carbon (gonads + gametes) per ind", "Population density")
.Deb_cohort$y1D$units <-
    c(rep("mmolC/ind", times = 3), "ind/m2")

.Deb_cohort$y0D <- data.frame(names =  #names from Fortran surbroutines, KEEP EXACT SAME NAMES & ORDER
  c("LARV_DENS", "LARV_BIOM", "PHYTO", "DETRITUS", "SIM", "O2"))
.Deb_cohort$y0D$description <- 
  c("Density of pelagic larvae", "Biomass of pelagic larvae",
    "Phytoplankton concentration", "Detritus concentration", 
    "Suspended inorganic matter (spm) concentration", 
    "Oxygen concentration")

.Deb_cohort$y0D$units <-
  c("ind/m3", "mmolC/m3",
    rep("mmolC/m3", times = 2), "g/m3", "mmolO2/m3")

.Deb_cohort$out1D <- data.frame(names =    #names from Fortran surbroutines, KEEP EXACT SAME NAMES & ORDER
  c("Gonads", "Gametes", "Length", "C_ind", 
    "Ingestion", "Catabolism", "Struct_growth", 
    "Faeces_production", "Regression", "Assimilation",  
    "Respiration", "Spawn_rate","Reprod_growth", 
    "VB_grow","WW_ind","Pseudofaeces_C","Filtration_C","Clearance_rate",
    "Mortality", "Settle_rate", "C_m2"))


.Deb_cohort$out1D$description <- 
  c("Gonad C per individual", "Gamete C per individual", 
    "Length of individual", "Total amount of C per individual",
    "C Ingestion per individual", 
    "Catabolism per individual", "Increase of structural C per individual", 
    "faeces production per individual", "Regression (respiration of gametes) per individual", 
    "Assimilation per individual", "Respiration per individual", 
    "Spawning rate per individual (expressed in C units)",
    "Increase of reproduciton buffer per individual",
    "recalcualted Von Bertalanffy growth rate of individual",
    "Wet weight approximation of individual",
    "Carbon ejected due to Pseudofeaces flux",
    "Carbon filtered from seawater by individual",
    "Seawater filtered per day by individual",
    "Mortality of cohort", 
    "Settlement rate to new cohort", 
    "Total amount of C per cohort")
.Deb_cohort$out1D$units <-
  c(rep("mmolC/ind", times = 2), "cm", "mmolC/ind",
    rep("mmolC/ind/d", times = 9),"cm/d","gWW/ind","mmolC/ind/d","mmolC/ind/d","m3/d",
    "ind/m2/d", "ind/m2/d",  "mmolC/m2")

.Deb_cohort$out0D <- data.frame(names = #names from Fortran surbroutines, KEEP EXACT SAME NAMES & ORDER
  c("Total_C_m2", "Total_POP_m2", "Total_DW_m2", "Total_WW_m2",
    "Total_POP_m2_FIELD", "Total_WW_m2_FIELD",
    "Total_ingestion", "Total_faeces_C", "Total_respiration",
    "Total_spawn", "comp_pedi", "Larv_weight",
    "Larv_mortality_rate","Larv_growth_rate"))

.Deb_cohort$out0D$description <- 
  c("Total carbon per m2", "Total population density per m2", "Total dry weight/m2", "Tot wet weight/m2",
    "Total population density (after mesh size filter)",
    "Total population WW biomass (after mesh size filter)",
    "Total ingestion rate per m2", "Total faeces production per m2", 
    "Total respiration per m2", "Total eggs spawned per m2 (expressed in mmolC)", 
    "numbers of larva ready to settle based on size (competent pediveligers)",
    "larval mean individual biomass",
    "Larval per capita mortality","Larval growth rate")

.Deb_cohort$out0D$units <-
  c("mmolC/m2", "ind/m2", "gDW/m2","gWW/m2",
    "ind/m2","gWW/m2",
    rep("mmolC/m2/d", times = 3), 
    "mmolC/m2/d", "m/d", "mmolC/ind", "/d","/d")

.Deb_cohort$forc <- data.frame(names =  #names from Fortran surbroutines, KEEP EXACT SAME NAMES & ORDER
   c("f_Phyto", "f_Temp", "f_Detritus", "f_Sim", "f_O2"))

.Deb_cohort$forc$description <- c(
   "External phytoplankton (algal) concentration", "Environmental temperature", 
   "External detritus (POC) concentration", 
   "External concentration of total suspended matter", 
   "External oxygen concentration")
.Deb_cohort$forc$units <- c(
   "mmolC/m3", "dgC", "mmolC/m3", "g/m3", "mmolO2/m3")

# instate auxiliry variables/quota used for unitilies, not passed to FORTRAN code
.Deb_cohort$auxiliary <- c(
  
  gAFDW_gSFDW = 0.8, #gram afdw per gram shell free dry weight
  gTW_gWW = 4.564126, # gram total weight (including shell) per gram wet weight (excluding shell) for mytilus
  POM_g_to_mmolC = 37.5, # POM (g) -> POC (mmolC) ratio used in calibration procedures
  CHL_to_PHY = 2.5, # chlorophyll a (ug) -> phytoplankton (mmolC) ratio (Smaal&vonc 1995)
  PHYTO_g_to_mmolC= 37.5 # phytoplankton weight (g) -> biomass (mmolC) ratio

)

  
