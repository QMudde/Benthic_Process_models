! ==============================================================================
! ==============================================================================
!
! Cohort mussel model using DEB dynamics
!
! implementation: Quinten Mudde
!                 Karline Soetaert
! NIOZ - Yerseke
!
! subroutines that interface with R:
! initmuspar, initmusforc, musselmod
! ==============================================================================
! ==============================================================================

! ==============================================================================
! Modules with declarations
! ==============================================================================

MODULE Mussel_dim

! the dimensions
! ------------------------------------------------------------------------------
! old feeding nparms = 43, new feeding nparms = 50 (48 from R + 2 metapars) 
! (crowding + substr + juv_mort_factor + Q10s = 56+2) + 1(threshold) + larvdmort + sizemort
 INTEGER, PARAMETER :: nparms = 62,  &  ! Number of parameters passed from R
                       nforcs = 5       ! Number of forcing functions from R

! output variables are divided in those that are one value (in common block)
! and the ones defined in all cohorts                      

   ! number of outputs --> outputs from pop totals (NOT allocatable dimension vars)
 INTEGER, PARAMETER :: nout   = 14 ! output variables defined in one value
 
 INTEGER :: N_Cohort        ! to be determined from the parameter vector
 
 INTEGER :: nout_Cohort     ! output variables defined in cohorts
 INTEGER :: numOutput       ! total output variables
 
 LOGICAL :: dynamic_environment  ! include interaction with environmental variables
 LOGICAL :: includePopulation    ! include population dynamics
 
END MODULE Mussel_dim

! ==============================================================================


! ==============================================================================

MODULE Mussel_Deb

USE Mussel_dim
IMPLICIT NONE

! ------------------------------------------------------------------------------
! state variables, for all cohorts
! ------------------------------------------------------------------------------

 DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) ::   &
    RESERVE ,  & ! [mmol.C/ind] reserve C
    STRUCT,    & ! [mmol.C/ind] structure C
    REPROD,    & ! [mmol.C/ind] carbon allocated to gonads and gametes
    POPULATION   ! [ind/m2] population density per cohort         

 DOUBLE PRECISION :: LARV_DENS, LARV_BIOM, PHYTO, DETRITUS, SIM, O2 

! ------------------------------------------------------------------------------
! derivatives
! ------------------------------------------------------------------------------

 DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) ::   &
    dRESERVE ,  & ! [mmol.C/ind/d] change in reserve C
    dSTRUCT,    & ! [mmol.C/ind/d] change in structure C
    dREPROD,    & ! [mmol.C/ind/d] change in gonads and gamete C
    dPOPULATION   ! [ind/m2/d] change in population density

 DOUBLE PRECISION :: dLARV_DENS, dLARV_BIOM, dPHYTO, dDETRITUS, dSIM, dO2 

! ------------------------------------------------------------------------------
! parameters
! ------------------------------------------------------------------------------

 DOUBLE PRECISION ::   &
   
   Num          , & ! number of cohorts (if 0, then: ONLY DEB)
   Env          , & ! interaction with environment

  ! DEB parameters
   aeff         , & ! 0.88, assimilation efficiency,     [-]
   kappa        , & ! 0.7,  fraction energy to structure [-]

  ! DEB rates
   maint        , & ! maintenance rate [/d]
   cat_s        , & ! surface-specific catabolism rate [mmolC/cm2/d]
   eg           , & ! mol C respired to create one mol structural C [-] 0.4074074
  
  ! Conversion factors     
   gDW_cm3      , & ! 0.09,   from cm3 to gramDWT      [gDWT/cm3]
   mmolC_cm3    , & ! 3.375,  from cm3 to mmolC  [mmolC/cm3]   
   gWW_gDW      , & ! 7,      from gram DW to gram WW [gram/gram]  
   del_M        , & ! 0.2254, shape coeff, from equivalent spherical length to length
  
 ! FOOD parameters (food concentration [FOOD] in molC/m3) [independent of environment]
 ! removed in new feeding mechanism
  ! phy_pref     , & ! preference factor for algae in food
  ! sim_inhib    , & ! 1,         inhibition constant for SIM effect on feeding [g/m3]
  ! k0_food      , & ! 2.275,     basic michaelis menten constant for feeding [mol/m3]

 ! energy contents ! removed in new feeding mechanism
 !  phy_E        , & ! 525, #[J/mmolC] (add my pet)
  ! det_E        , & ! 150, #[J/mmolC]
   !mussel_E     , & ! 525, #[J/mmolC] (add my pet)

 ! Oxygen limitation 
   ks_o2        , & ! 20,   oxygen limitation michaelis menten coeff    [mmol O2/m3]
   crit_o2      , & ! 50,   oxygen limitation threshold (below this: 0) [mmol O2/m3]
  
 ! Mortality parameters 
   mort_starv   , & ! 0.009,  increase of mortality during starvation [/day]
   mort_base    , & ! 0.001,  population mortality due to environment [/day]
   e_dens_trh   , & ! 0.85, average indiv energy density at which pop starvation starts [mmolC R/mmolC S]

  
 ! Reproduction & Maturity parameters #
   mat_pub      , & ! 3.022857, maturity content (gonads) at puberty [mmolC/ind]
   spawn_peak   , & ! 107, day of year at which spawning peaks [day]
   spawn_peak_2 , & ! [day] of second spawning each year (default = 0 = no second spawning)
   spawn_dev    , & ! 5 deviation (spread) of spawning intensity [day]
   spawn_eff    , & ! 1, fraction of gametes released during spawning [-]  
  
   weight_egg  , & ! 9.246825e-07, # from AmP ww at birth [mmol C/egg]

 ! initial state variables for sessile stage (at metamorphosis)
   weight_settle   , & ! 2.1e-05, [mmolC/ind] biomass of larvae at which they start cohort
   v_settle        , & ! 10  [m/d]  settling velocity of larvae that settle
   pReserve_settle , & ! 0.5 

 ! Temperature dependence

   T_a,   & !  = 5800  Arrhenius Temperature  [[K]]
   T_r,   & !  = 293   reference temperature  [[K]]
   T_l,   & !  = 275   lower boundary of tolerance range [[K]]
   T_h,   & !  = 296   upper boundary of tolerance range [[K]]
   T_al,  & !  = 45430 rate of decrease of LOWER boundary [[K]]
   T_ah,  & !  = 31376 rate of decrease of UPPER boundary [[K]]

 ! Effect of crowding?
 ! Larvae parameters 

   mort_larva   , & ! 0.4385876, additional mortality of susceptible larval stage [0.05?]
   ! Bivalve larva approximate .13 mortality per day (Jørgensen, 1981)  
   growth_larva , & ! 0.184505, [/day] larval growth
   ks_larva     , & ! 
   max_set_rate , & !
   ks_settlement, & !

 ! Environmental variables

   water_renew  , & ! food conc renewal rate (default: set to 1)  [/day]
   O2_renew     , & ! O2 renewal rate  [/day]
   depth        , & ! water depth [m] relevant for box volume
   
 ! New feeding mechanism parameters 
   
   max_clear        , & ! = 0.096 [m3/d/cm2] (Saraiva, 2012) 
   max_phyto_filt   , & ! = 1.807048 [mmolC/cm2/d] (Estimated)
   max_DT_filt      , & ! = 1.807048 [mmolC/cm2/d] (Estimated)
   max_spm_filt     , & ! = 3.5 [g/cm2/d] (Saraiva 2012)
   rho_Phyto        , & ! = 0.99 [-] (Saraiva, 2012)
   rho_DT           , & ! = 0.8093233 [-] (Estimated)
   rho_SIM          , & ! = 0.45 [-] (Saraiva, 2012)
   max_Phyto_ing    , & ! = 0.65e7 (this is a strage value in source literature, recalibrated)
   max_DT_ing       , & ! = 37 [mmolC/d] (nlm) mussels make pseudo under pure DT food
   max_SIM_ing      , & ! = 0.23 [g/d] (Saraiva, 2012)
   nc_tissue        , & ! = 0.21 [molN/molC] (Saraiva, 2012)
   nc_Phyto         , & ! = 0.151 [molN/molC] (redfield ratio)
   nc_DT            , & ! = 0.04 [molN/molC] (rice et al. 1981)
   
   ! pars for crowding 
   shell_allo_a ,   & ! allometric for shell surf (linter intercept)
   shell_allo_b ,   & ! allometric for shell surf (linear slope)
   crowd_lim,       & ! (max shell cover before crowding mortality)
   substr_factor,   & ! substrate factor inpact on larvae settelment (default = 1)
   mort_juv,    & ! maximum mortality due to small size
   mort_egg,         & ! inital egg mortality fraction [-]
   q10_LgL,         & ! Q10 for larval growth Lower (<= 20dC) BASE = T_r (293k)
   q10_Lm,          & ! Q10 for larval mortality (mirrored around 20)  BASE = T_r (293k)
   length_trh,      & ! [cm] threshold to count cohort to expected field data
   recr_length        ! size at recruitment 


! The common block with all parameter declarations (order = important) 
! nr. of cohorts (overarching par) + env switch + 47 pars imported from R  +(3crowd) 
 common / musselparms / Num, Env,                                    &
    aeff, kappa, maint, cat_s, eg,                                   &
    gDW_cm3, mmolC_cm3, gWW_gDW, del_M,                              &
    ks_O2, crit_O2, e_dens_trh, mort_starv, mort_base, mat_pub,      &
    spawn_peak, spawn_peak_2, spawn_dev, spawn_eff, weight_egg,      & 
    weight_settle, v_settle, pReserve_settle,                        &
    T_a, T_r, T_l, T_h, T_al, T_ah,                                  &
    mort_larva, growth_larva, ks_larva ,                             &
    max_set_rate , ks_settlement,                                    &
    water_renew, O2_renew, depth,                                    &
    max_clear, max_phyto_filt, max_DT_filt, max_spm_filt,            & 
    rho_Phyto, rho_DT, rho_SIM, max_Phyto_ing, max_DT_ing,           &
    max_SIM_ing, nc_tissue, nc_Phyto, nc_DT,                         &     
    shell_allo_a, shell_allo_b, crowd_lim, substr_factor,            &
    mort_juv, mort_egg, q10_LgL, q10_Lm,                             &
    length_trh, recr_length
                                          

! ------------------------------------------------------------------------------
! Forcing functions
! ------------------------------------------------------------------------------

 DOUBLE PRECISION ::   &
 
    f_Phyto,    & ! algal biomass (was f.CH.a = chloro.a), [mmolC/m3]
    f_Temp,     & ! environmental temperature, [Deg_c]
    f_Detritus, & ! seawater detritus, mmolC/m3,  [mmolC/m3]
    f_Sim,      & ! suspended inorganic particulate matter [g/m3]
    f_O2          ! seawater oxygen concentration [mmol/m3]

 common / musselforcs / f_Phyto, f_Temp, f_Detritus, f_Sim, f_O2

! ------------------------------------------------------------------------------
! Output variables
! ------------------------------------------------------------------------------

 DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) ::   &
    Ingestion, Catabolism, Struct_growth, Faeces_production,         &
    Regression, Assimilation, Respiration, Mortality,                &
    Spawn_rate, Reprod_growth, Settle_rate
    
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) ::  &
   Pseudofaeces_C, FR_Phyto, FR_DT,               &
   Filtration_C, Clearance_rate

 DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) ::   &
    Gonads, Gametes, Length, C_ind, C_m2, Shell_surf, VB_grow, WW_ind  
    
 DOUBLE PRECISION :: &
    Total_C_m2, Total_POP_m2, Total_DW_m2, Total_WW_m2,              &
    Total_POP_m2_FIELD,Total_WW_m2_FIELD,                            &
    Total_ingestion,Total_FR_Phyto, Total_FR_DETRITUS,               &
    Total_psfaeces_C, Total_faeces_C, Total_respiration,             &
    Total_spawn, Total_shell_cover
    
 DOUBLE PRECISION :: comp_pedi, L_weight, larv_mort_rate, larv_growth_rate
   
 common / musselout /                                                &
    Total_C_m2, Total_POP_m2, Total_DW_m2, Total_WW_m2,              &
    Total_POP_m2_FIELD, Total_WW_m2_FIELD,                           &
    Total_ingestion, Total_faeces_C, Total_respiration,              &
    Total_spawn, comp_pedi, L_weight, larv_mort_rate, larv_growth_rate
  
! Other variables, not exported to R, but used in different subroutines

DOUBLE PRECISION ::  ox_factor, T_factor! p_phyto, p_det, (removed)

END MODULE Mussel_deb 


