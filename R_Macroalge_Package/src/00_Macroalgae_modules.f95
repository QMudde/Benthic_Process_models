!==========================================================================
!==========================================================================
! The macroalgal model in FORTRAN
!    A model for Macroalgal growth 
! author:  
! Quinten Mudde, Karline Soetaert
! 
!==========================================================================
!==========================================================================
  
!==========================================================================
! Module with common blocks
!==========================================================================
    
  MODULE macroalgae_mod

!---------------------------- 
! State variables
!---------------------------- 
    implicit none 
    
    DOUBLE PRECISION  ::     &
        STRUCT_C,     & ! [mmolC/m2] structural C of macroalgae
        RESERVE_C,    & ! [mmolC/m2] reserve C of macroalgae
        RESERVE_N,    & ! [mmolN/m2] reserve N of macroalgae
        RESERVE_P,    & ! [mmolP/m2] reserve P of macroalgae
        NO3,          & ! [mmolN/m3] nitrate concentration in water
        NH4,          & ! [mmolN/m3] ammonium concentration in water
        PO4,          & ! [mmolP/m3] phosphate concentration in water
        O2,           & ! [mmolO2/m3] oxygen concentration in water
        C_CLOSING,    & ! [mmolC/m2] all C that goes in and out of the Box model
        N_CLOSING,    & ! [mmolN/m2] all N that goes in and out of the Box model
        P_CLOSING       ! [mmolP/m2] all P that goes in and out of the BOX model

!---------------------------- 
! derivatives, 
! mmol/m2/d or mmol/m3/d
!---------------------------- 
    DOUBLE PRECISION  :: dSTRUCT_C, dRESERVE_C, dRESERVE_N, dRESERVE_P 
    DOUBLE PRECISION  :: dNO3, dNH4, dPO4, dO2 
    DOUBLE PRECISION  :: dC_CLOSING, dN_CLOSING, dP_CLOSING

!---------------------------- 
! forcing functions
!---------------------------- 
    DOUBLE PRECISION  ::    &
       f_Temperature,  & ! [dg C] environmental temperature
       f_Current,      & ! [m/s] current velocity
       f_Light,        & ! [W/m2] light at the air-water interface 
       f_PHYTO,        & ! [mmolN/m3] Phytoplankton concentration in overlying water
       f_ExtNO3,       & ! [mmolN/m3] Nitrate concentration in overlying water
       f_ExtNH4,       & ! [mmolP/m3] Ammonium concentration in overlying water
       f_ExtPO4,       & ! [mmolC/m3] Phosphate concentration in overlying water
       f_ExtO2           ! [mmolO2/m3] Oxygen concentration in overlying water       
     COMMON /myforcsM/ f_Temperature, f_Current, f_Light, f_PHYTO,              &
        f_ExtNO3, f_ExtNH4, f_ExtPO4, f_ExtO2
       
!---------------------------- 
! parameters
!---------------------------- 
    DOUBLE PRECISION  :: pC_SDW, RC_SCmin, SN_SCr, RN_SCmin, RN_SCmax,          &
       SP_SCr, RP_SCmin, RP_SCmax, gdw_RC, gdw_RN, gdw_RP, gdw_SC,              &
       OC_q,dw_ww, pMax, alpha, beta, rho_short, rho_mid, rho_long,             &
       t_Arh, t_A_L, t_A_H, t_LB, t_UB, t_ref, ExhFmax, h_ALG, a_ALG,           &
       ww_ALG2, h_ALG2, a_ALG2, gRmax, basalResp, growthResp,                   &
       tRChange, tRRef, max_NO3_upt, ksNO3, max_NH4_upt, ksNH4,                 &
       max_PO4_upt, ksPO4, u065, mort_rate, nec_rate, ero_max,                  &
       ero_inc, k0_short, k0_mid, k0_long, k0_turb,                             &
       k_Phyto_short, k_Phyto_mid, k_Phyto_long,                                &
       latitude, depth, dilution
      
    COMMON /myparmsM   /pC_SDW, RC_SCmin, SN_SCr, RN_SCmin, RN_SCmax,           &
       SP_SCr, RP_SCmin, RP_SCmax, gdw_RC, gdw_RN, gdw_RP, gdw_SC,              &
       OC_q, dw_ww, pMax, alpha, beta, rho_short, rho_mid, rho_long,            &
       t_Arh,t_A_L, t_A_H, t_LB, t_UB, t_ref, ExhFmax, h_ALG, a_ALG,            &
       ww_ALG2, h_ALG2, a_ALG2, gRmax, basalResp, growthResp,                   &
       tRChange, tRRef, max_NO3_upt, ksNO3, max_NH4_upt, ksNH4,                 &
       max_PO4_upt, ksPO4, u065, mort_rate, nec_rate, ero_max,                  &
       ero_inc, k0_short, k0_mid, k0_long, k0_turb,                             &
       k_Phyto_short, k_Phyto_mid, k_Phyto_long,                                &
       latitude, depth, dilution
       
!---------------------------- 
! output variables (66)
!---------------------------- 
    DOUBLE PRECISION  ::         & !
      
       RC_SCr,              & ! molRC/molSC    Carbon (C) reserve quotum
       RN_SCr,              & ! molRN/molSC    Nitrogen (N) reserve quotum
       RP_SCr,              & ! molRP/molSC    Phosphorus (P) reserve quotum
       WW,                  & ! gWW/m2         wet weight macroalgae biomass
       DW,                  & ! gDW/m2         dry weight macroalgae biomass
       Macroalgae_C,        & ! mmolC/m2       total macroalgae carbon
       Macroalgae_N,        & ! mmolN/m2       total macroalgae nitrogen
       Macroalgae_P,        & ! mmolP/m2       total macroalgae phophorus
       GrossPS,             & ! mmolC/m2/d     total gross photosynthesis
       GrossPR,             & ! mmolC/mmolSC/d gross specific Photosynthesis rate at light and T
       Photosynthesis,      & ! mmolC/m2/d     photosynthesis - exudation
       ExudFrac,            & ! -              fraction of gross photosynthesis exudated
       Exudation,           & ! mmolC/m2/d     total C exudation
       Basal_Resp,          & ! mmolC/m2/d     basal respiration paid by reserve C
       Growth_Resp,         & ! mmolC/m2/d     C paid for growth (activity respiration)
       Respiration,         & ! mmolC/m2/d     total respiration (-necrosis)
       Respiration_loss,    & ! mmolC/m2/d     total respiration (+necrosis)
       NO3_sat,             & ! -              Monod limitation term for NO3 uptake
       NH4_sat,             & ! -              Monod limitation term for NH4 uptake
       PO4_sat,             & ! -              Monod limitation term for PO4 uptake
       NU,                  & ! -              dependency of nutrient uptake on current velocity
       RC_sat,              & ! molRC:molSC    N reserve saturation of nutrient uptake and growth
       RN_sat,              & ! molRN:molSC    C reserve saturation of growth
       RP_sat,              & ! molRP:molSC    P reserve saturation of nutrient uptake and growth
       Quota_growth_fac,    & ! -              limitation of growth due to quota
       Photo_Eff,           & ! -              part of photosynthesis realised to light only
       Light_Eff,           & ! -              part of photosynthesis realised to light and temperature
       NO3_upt_Eff,         & ! -              N03 uptake efficiency (rel. to max uptake)
       NH4_upt_Eff,         & ! -              NH4 uptake efficiency (rel. to max uptake)
       PO4_upt_Eff,         & ! -              PO4 uptake efficiency (rel. to max uptake)
       NO3_upt,             & ! mmolN/m2/d     total NO3 uptake
       NH4_upt,             & ! mmolN/m2/d     total NH4 uptake
       PO4_upt,             & ! mmolP/m2/d     total PO4 uptake
       Growth,              & ! mmolC/m2/d     total structural growth
       Growth_rate,         & ! /d             relative growth rate
       Growth_Eff,          & ! -              growth efficiency
       Net_Growth_rate,     & ! /d             net structural increase (dSTRUCTC/STRUCTC)
       Mortality,           & ! mmolSC/m2/d    total mortality, structural C
       Necrosis,            & ! mmolSC/m2/d    total necrosis, structural C
       Regression,          & ! mmolSC/m2/d    respiration paid by structure C
       Erosion,             & ! mmolSC/m2/d    total erosition, structural C
       pReg,                & ! -              part respriration paid by structural C
       Necrosis_rate,       & ! /d             C-specific necrosis rate
       Erosion_rate,        & ! /d             C-specific erosion rate
       Mortality_rate,      & ! /d             C-specific mortality rate
       SurfPAR_E,           & ! uMolPhot/m2/s  irradiance on water surface
       E_short,             & ! uM/m2/s        blue-light band (400-500nm) on algae
       E_mid,               & ! uM/m2/s        green-light band (500-600nm) on algae
       E_long,              & ! uM/m2/s        red-light band (600-700nm) on algae
       E_Alg,               & ! uM/m2/s        total irradiance for photosynthesis
       K_mean,              & ! /m             average attenuation factor for light
       K_MacroAlg,          & ! /m             extinction coeff due to self shading
       K_Competition,       & ! /m             extinction coeff competing algae
       Shading_MacroAlg,    & ! -              shading factor due to modeled algae
       Shading_Competition, & ! -              shading factor due to competing algae
       Shading_Eff,         & ! -              part of light not shaded by algae
       TF_gain,             & ! -              temperature factor for photosynthesis and growth
       TF_resp,             & ! -              temperature factor for respiration and necrosis
       Total_C,             & ! mmolC/m2       total carbon in the system
       Total_N,             & ! mmolN/m2       total nitrogen in the system
       Total_P,             & ! mmolP/m2       total phophorus in the system
       Temperature,         & ! degC           temperature
       Current,             & ! m/s            current speed
       PHYTO,               & ! mmolC/m3       phytoplankton concentration in water column 
       Light,               & ! W/m2           surface irradiation
       ExtNO3,              & ! mmolN/m3       bottom water NO3 concentration
       ExtNH4,              & ! mmolN/m3       bottom water NH4 concentration
       ExtPO4,              & ! mmolP/m3       bottom water PO4 concentration
       ExtO2                  ! mmolO2/m3      bottom water O2 concentration


    COMMON /myoutM     /RC_SCr, RN_SCr, RP_SCr, WW, DW, Macroalgae_C,           &
       Macroalgae_N, Macroalgae_P, grossPS, grossPR, Photosynthesis,            &
       ExudFrac, Exudation, Basal_Resp, Growth_Resp, Respiration,               &
       Respiration_loss, NO3_sat,NH4_sat, PO4_sat, NU, RC_sat,                  &
       RN_sat, RP_sat, Quota_growth_fac, Photo_Eff, Light_Eff,                  &
       NO3_upt_Eff, NH4_upt_Eff, PO4_upt_Eff, NO3_upt, NH4_upt,                 &
       PO4_upt, Growth, Growth_rate, Growth_Eff, Net_Growth_rate,               &
       Mortality, Necrosis, Regression, Erosion, pReg, Necrosis_rate,           &
       Erosion_rate, Mortality_rate, SurfPAR_E, E_short, E_mid, E_long,         &
       E_Alg, K_mean, K_MacroAlg, K_Competition, Shading_MacroAlg,              &
       Shading_Competition, Shading_Eff, TF_gain, TF_resp, Total_C,             &
       Total_N, Total_P, Temperature, Current, PHYTO, Light,                    &
       ExtNO3, ExtNH4, ExtPO4, ExtO2
       

  END MODULE macroalgae_mod
      
