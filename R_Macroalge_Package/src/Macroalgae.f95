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
! initialise the common block with parameter values, 
!==========================================================================
    
       SUBROUTINE initmacroalgae (funparms)
       IMPLICIT NONE
       EXTERNAL steadyparms
  
       INTEGER,PARAMETER :: nc = 58
  
       DOUBLE PRECISION parms(nc)
       COMMON /myparmsM/parms
  
       CALL funparms(nc, parms)
  
       RETURN
       END SUBROUTINE initmacroalgae
  
!==========================================================================
! Initialise the forcing function common block  
!==========================================================================
    
       SUBROUTINE initforcm (funforcs)
       IMPLICIT NONE
       EXTERNAL steadyforcs
  
       INTEGER,PARAMETER :: N = 8 
  
       DOUBLE PRECISION forcs(N)
       COMMON /myforcsM/forcs
  
       CALL funforcs(N, forcs)
  
       RETURN
       END SUBROUTINE initforcm
  
!==========================================================================
!==========================================================================
! subroutine calculating the rate of change of
! the macroalgal model - here the quantities in the common blocks are named
!==========================================================================
!==========================================================================
    
    SUBROUTINE macroalgae (neq, t, Conc, dConc, yout, ip)
       
    USE macroalgae_mod       ! declaration of states, variables, parameters
    IMPLICIT NONE
  
!......................... declaration section.............................
    INTEGER           :: neq, ip(*), i
  
    DOUBLE PRECISION  :: t, Conc(11), dConc(11), yout(*)

!---------------------------- 
! ordinary variables 
!---------------------------- 

    DOUBLE PRECISION :: TempK, a, b, Arrh_denominator, Arrh_numerator
    DOUBLE PRECISION :: E_long_canopy,  E_long_surf,                            &
                        E_mid_canopy,   E_mid_surf,                             &
                        e_short_canopy, E_short_surf
    DOUBLE PRECISION :: ChlA, K_long, K_mid, K_short, N_return,                 &
                        NH4_env_exchange, NH4_return, NH4_upt_rate,             & 
                        NO3_env_exchange, NO3_return, NO3_upt_rate,             & 
                        PO4_env_exchange, PO4_return, PO4_upt_rate,             & 
                        O2_Production, O2_Respiration, OX_env_exchange   
    DOUBLE PRECISION :: PmaxT, PS, Resp, STRUCT_N, STRUCT_P, T_atREF,           &
                        Total_Resp, SurfPAR_Wm2  

!............................ statements ..................................
  
!     check memory allocated to output variables
      IF (ip(1) < 66)  CALL rexit("nout should be at least 66") 
  
! from Conc to state variables
      STRUCT_C    = Conc(1)
      RESERVE_C   = Conc(2)
      RESERVE_N   = Conc(3)
      RESERVE_P   = Conc(4)
      NO3         = Conc(5)
      NH4         = Conc(6)
      PO4         = Conc(7)
      O2          = Conc(8)
      C_CLOSING   = Conc(9)
      N_CLOSING   = Conc(10)
      P_CLOSING   = Conc(11)

   
!=============================!
! Initialisation              !
!=============================!
    
! set forcing functions  

    Temperature  = f_Temperature  ! temperature, degrees celsius
    Current     = f_Current       ! current speed, [m/s]
    Light       = f_Light         ! surface irradiation, [W/m2]
    ExtNO3      = f_ExtNO3        ! external water NO3 concentration [mmolN/m3]
    ExtNH4      = f_ExtNH4        ! external water NH4 concentration [mmolN/m3]
    ExtPO4      = f_ExtPO4        ! external water PO4 concentration [mmolP/m3]
    ExtO2       = f_ExtO2         ! external water O2 concentration [mmolO2/m3]
    PHYTO       = f_PHYTO         ! phytoplankton concentration in water [mmolC/m3]

    TempK     = Temperature + 273.15d0  ! temperature, degrees kelvin

    !---------------------------------------------------------------#
    ! elemental totals, ratios and quota
    !---------------------------------------------------------------#
    
    STRUCT_N     = STRUCT_C * SN_SCr      ! [mmolN/m2] structural nitrogen 
    STRUCT_P     = STRUCT_C * SP_SCr      ! [mmolP/m2] structural phosphorous 
    MacroAlgae_C = STRUCT_C + RESERVE_C   ! [mmolC/m2] total macroalgae carbon
    MacroAlgae_N = STRUCT_N + RESERVE_N   ! [mmolN/m2] total macroalgae nitrogen
    MacroAlgae_P = STRUCT_P + RESERVE_P   ! [mmolP/m2] total macroalgae phophorus
     
    ! [[ Quota ]] (relative to structural C)
    RC_SCr   = RESERVE_C/(STRUCT_C+1.d-8) ! [molRC/molSC] Carbon (C) reserve quotum
    RN_SCr   = RESERVE_N/(STRUCT_C+1.d-8) ! [molRN/molSC] Nitrogen (N) reserve quotum
    RP_SCr   = RESERVE_P/(STRUCT_C+1.d-8) ! [molRP/molSC] Phosphorus (P) reserve quotum
    
    ! dry and wet Weights
    DW       =  gdw_SC * STRUCT_C                                               &
              + gdw_RC*RESERVE_C                                                &
              + gdw_RN*RESERVE_N                                                &
              + gdw_RP*RESERVE_P           ! [gDW/m2] total dry weight
    
    WW       = DW / dw_ww                  ! [gWW/m2] wet weight (biomass in most measurements)
    
    !---------------------------------------------------------------#
    ! LIGHT AT DEPTH  Akitsu et al., 2015
    !---------------------------------------------------------------#
    
    ! Total irradiance (shortwave heatflux) (300-2600nm) is converted to PAR (400-700nm)
    SurfPAR_Wm2 = Light * 0.45d0 ! [W/m2]    
    
    ! Irradiance [W/m2] is recalculated to photon flux [uMol photons/m2/s] 
    SurfPAR_E   = SurfPAR_Wm2 * 4.57d0 
    
    ! PAR will be divided into 3 bands of wavelengths 
    ! short: 400-500 mid:500-600 and long 600-700
    ! ASSUME:: that at sea surface, these 3 bands are of equal photon flux
    E_short_surf = SurfPAR_E/3.d0 !  blue light band of wavelengths 400-500
    E_mid_surf   = SurfPAR_E/3.d0 ! green light band of wavelengths 500-600
    E_long_surf  = SurfPAR_E/3.d0 !   red ligth band of wavelengths 600-700
    
    ChlA         = PHYTO/2.5d0    ! [mg/m3] chlorophyll A
    
    ! Light based on depth, phytoplankton turbidity, all 3 wavebands attenuate.
    ! 3 spectra have different background attenuation due to water (e.g. k0clear*rk_short)
    ! all 3 spectra are equally affected by non-phytoplankton turbidity ASSUMPTION. (k0turb)
    ! spectra are differently attenuated by average Chlorophyll A conc. (k_Phyto_long*ChlA)
    
    K_short = k0_short + k0_turb  + (k_Phyto_short * ChlA)
    K_mid   = k0_mid   + k0_turb  + (k_Phyto_mid   * ChlA)
    K_long  = k0_long  + k0_turb  + (k_Phyto_long  * ChlA)
    
    K_mean = (K_short + K_mid + K_long)/3.d0 ! average attenuation factor
    
    ! Resulting Attenuation factors apply on light via Lambert Beer law at the provided depth
    
    ! Light at middle of plant height [[h_alg*0.5]]. 
    
    E_short_canopy = E_short_surf * exp(-K_short * (depth - (h_ALG*0.5d0)))
    E_mid_canopy   = E_mid_surf   * exp(-K_mid   * (depth - (h_ALG*0.5d0)))
    E_long_canopy  = E_long_surf  * exp(-K_long  * (depth - (h_ALG*0.5d0)))
    
    ! Self-shading
    ! extinction rate due to the target Macroalga, and the competing algae
    K_MacroAlg    = a_ALG  * Macroalgae_C / (h_ALG + 1e-10) 
    K_Competition = a_ALG2 * ww_ALG2 / h_ALG2 
    
    Shading_MacroAlg   = exp(-K_MacroAlg * h_ALG*0.5) ! shading at half of the height
    
    ! competition over difference in height - QUINTEN: KAN NEGATIEF WORDEN
    Shading_Competition = exp(-K_Competition *                                  &
                                max(0.d0, (h_ALG2 - (0.5d0*h_ALG)))) 

    E_short = E_short_canopy * Shading_MacroAlg * Shading_Competition
    E_mid   = E_mid_canopy   * Shading_MacroAlg * Shading_Competition
    E_long  = E_long_canopy  * Shading_MacroAlg * Shading_Competition
    
    ! Actual irradiance on photosynthesis based on absorption spectra ---
    E_Alg  =  rho_short * E_short                                               &
            + rho_mid   * E_mid                                                 & 
            + rho_long  * E_long ! [u mol photons/m2/s]

    !---------------------------------------------------------------------------
    ! Arrhenius relationship of temperature for photosynthesis and growth
    !---------------------------------------------------------------------------
    
    Arrh_numerator   = exp( (t_Arh/t_ref) - (t_Arh/TempK) )
    Arrh_denominator = 1.d0 + exp( (t_A_L/TempK) - (t_A_L/t_LB ))               &
                            + exp( (t_A_H/t_UB ) - (t_A_H/TempK))
    
    TF_gain          = Arrh_numerator/Arrh_denominator
    T_atREF          = 1.d0 / (1.d0 + exp((t_A_L/t_ref) - (t_A_L/t_LB ))        &
                                    + exp((t_A_H/t_UB ) - (t_A_H/t_ref)))
    
   ! TF_gain          = TF_gain/T_atREF ! only for rates measured at ref temp 
    
    !-----------------------
    ! Photosynthesis
    !-----------------------
    
    ! max photosynthetis at current temperature
    PmaxT = pMax * (TF_gain/T_atREF) ! [mmolC/mmolSC/d] 

    ! theoretical PI curve max
    PS    = PmaxT / ( (alpha / (alpha+beta)) * (beta/(alpha+beta)**(alpha/beta))) 
    a     = alpha * E_Alg/PS
    b     = beta  * E_Alg/PS
    
    ! [mmolC/mmolSC/d], gross Structural C specific Photo rate at light and T
    grossPR = PS * (1.d0-exp(-a))*exp(-b) 
    
    Light_Eff = grossPR / PmaxT ! effectiveness of Photosynthesis due to light only
    Photo_Eff = grossPR / pMax  ! effectiveness of photosynthesis due to light and temperature
    
    ! Gross photosynthesis 
    grossPS   = grossPR * STRUCT_C ! [mmolC/m2/d] 
    
    !---------------------------------------------------------------------------
    ! Carbon exudation
    !---------------------------------------------------------------------------
    
    ExudFrac          = 1.d0 - exp(ExhFmax * min(0.d0, RC_SCmin - RC_SCr))
    Exudation         = grossPS * ExudFrac
    Photosynthesis    = grossPS - Exudation
    
    !---------------------------------------------------------------------------
    !  Nutrient uptake (Nutrients from environment to Nutrient reserves)[N & P]
    !---------------------------------------------------------------------------
    
    ! limited by external nutrient concentrations and internal quota,
    ! also dependent on current speed 
    ! if flux is negative -> excretion of N or P 
    
    ! effect of reserveN:structureC quota (droop kinetics) on N uptake
    ! N reserve saturation of nutrient uptake and growth (between 0 and 1)
    ! low RN_sat value --> low uptake multiplyer = high internal storage
    RN_sat = min(1.d0, (RN_SCmax-RN_SCr) / (RN_SCmax-RN_SCmin)) ! [-] 
    RN_sat = max(0.d0, RN_sat)
    
    ! effect of reserveP:structureC quota (droop kinetics) on P uptake
    ! described as P reserve saturation (between 0 and 1)
    ! low RP_sat value --> low uptake multiplyer = high internal storage
    RP_sat = min(1.d0, (RP_SCmax-RP_SCr) / (RP_SCmax-RP_SCmin)) ! [-]
    RP_sat = max(0.d0 , RP_sat)
    
    ! Monod saturation formulation of uptake
    NO3_sat   = NO3 / (NO3+ksNO3) ! [-] Monod limitation term for NO3 uptake
    NH4_sat   = NH4 / (NH4+ksNH4) ! [-] Monod limitation term for NH4 uptake
    PO4_sat   = PO4 / (PO4+ksPO4) ! [-] Monod limitation term for PO4 uptake
    
    ! functional dependency of nutrient uptake on current velocity (between 0-1)
    NU        = 1.d0 - exp(-Current/u065)  ! [-]   
    
    ! Nutrient uptake efficiency (relative to max specific uptake) dimensionless
    NO3_upt_Eff = NO3_sat * RN_sat * NU 
    NH4_upt_Eff = NH4_sat * RN_sat * NU
    PO4_upt_Eff = PO4_sat * RP_sat * NU
    
    ! Specific nutrient uptake
    NO3_upt_rate = max_NO3_upt * NO3_upt_Eff ! [mmolN/mmolSC/d]
    NH4_upt_rate = max_NH4_upt * NH4_upt_Eff ! [mmolN/mmolSC/d]
    PO4_upt_rate = max_PO4_upt * PO4_upt_Eff ! [mmolN/mmolSC/d]

    ! total nutrient uptake    
    NO3_upt = NO3_upt_rate * STRUCT_C ! [mmolN/m2/d] total NO3 uptake
    NH4_upt = NH4_upt_rate * STRUCT_C ! [mmolN/m2/d]
    PO4_upt = PO4_upt_rate * STRUCT_C ! [mmolP/m2/d]
    
    !----------------------------------------------------------
    ! STRUCTURAL GROWTH (Based on Reserve pools) & temperature
    !----------------------------------------------------------
    
    ! functional dependency of growth on reserve:structure quota
    ! the lowest quotum limits the growth (C or N storage) via droop kinetics
    ! C- limitation factor (1 is fully limiting, 0 is non-limiting)
    RC_sat = RC_SCmin/RC_SCr 
    RC_sat = min(1.0d0, RC_sat)
    RC_sat = max(0.0d0, RC_sat) 

    ! acquire limiting quotum (DROOP); force between 0 and 1
    Quota_growth_fac  = (1.0 - max(RN_sat, RC_sat, RP_sat)) 
    Quota_growth_fac  = max(0.0d0, Quota_growth_fac) 
    Quota_growth_fac  = min(1.0d0, Quota_growth_fac) 

    ! growth efficiency, [-]
    Growth_Eff = TF_gain * Quota_growth_fac
    
    ! Relative growth rate
    Growth_rate = gRmax * Growth_Eff  ! [/d] 

    ! Total growth
    Growth     = Growth_rate * STRUCT_C      ! [mmolC/m2/d]
    
    !====================================================#
    ! LOSS TERMS (Both losses in reserves and structure) #
    !====================================================#
    
    !----------------------------------------------------------
    ! Respiration
    !----------------------------------------------------------
    
    ! simple exponential relationship between respiration and temperature
    TF_resp = exp(tRChange * (Temperature - tRRef))
    
    ! total mmolC needed for respiration
    Resp    = TF_resp * basalResp * STRUCT_C ! [mmolC/m2/d] respired
    
    ! respiration is paid by structure if there are not enough reserves 
    ! part respriration paid by structural C (regression)
    pReg = 0.d0
    if (RC_SCr <= RC_SCmin) pReg =  1.d0 - (RC_SCr/RC_SCmin)
    
    ! Total respiration that needs to be paid by structure and reserve Carbon
    Regression = pReg * Resp   ! [mmolC/m2/d] respiration paid by structure C
    Basal_Resp = Resp*(1.d0-pReg) ! [mmolC/m2/d] respiration paid by reserve C 
    
    !     activity (growth) respiration and total respiration 
    Growth_Resp = growthResp * Growth     ! [mmolC/m2/d] C paid for growth (activity respiration)
    Total_Resp = Basal_Resp + Growth_Resp ! [mmolC/m2/d] total respiration                             
    
    
    !----------------------------------------------------------
    ! Mortality, Erosion & Necrosis
    !----------------------------------------------------------

    Erosion_rate = ero_max * 1e-6 * exp(ero_inc*STRUCT_C)/                      &  
            (1. + 1e-6*( exp (ero_inc*STRUCT_C)-1.) )   ! [mmolC/mmolSC/d]  C-specific erosion rate
    
    Erosion       = Erosion_rate * STRUCT_C             ! [mmolC/m2/d] total erosion
    
    if (TempK > t_UB) then
       Necrosis_rate = nec_rate * TF_resp   ! [/d] Carbon specific necrosis (due to heat)
    else
       Necrosis_rate = 0.d0
    endif   
    
    Necrosis      = Necrosis_rate * STRUCT_C            ! [mmolC/m2/d] total necrosis
    
    Mortality_rate = (Necrosis_rate + mort_rate + Erosion_rate)  ! [/d] loss rate
    Mortality      = Mortality_rate * STRUCT_C   ! [mmolSC/m2/d] Total mortality

    !========================================
    ! Total exchange rates with environment
    !=======================================
    
    ! processes affecting dissolved Nitrogen & Phosphorus per m2
    ! Nuptake [mmolN/m2/d]
    
    N_return   = (Mortality + Regression) * SN_SCr                              &
                 + Mortality * RN_SCr                 ! [mmolN/m2/d mortality]
    
    NO3_return = N_return*0.5
    NH4_return = N_return*0.5
    
    PO4_return = (Mortality + Regression)*SP_SCr                                &
                 + Mortality * RP_SCr
    
    ! processes affecting dissolved oxygen per m2
    O2_Production  = grossPS*OC_q                 ! [mmolO2/m2/d]
    O2_Respiration = (Total_Resp+Regression)*OC_q ! [mmolO2/m2/d]
    
    ! exchange with external environment
    NO3_env_exchange = dilution*(ExtNO3 - NO3)
    NH4_env_exchange = dilution*(ExtNH4 - NH4)
    PO4_env_exchange = dilution*(ExtPO4 - PO4)
    OX_env_exchange  = dilution*(ExtO2  -  O2)
    
    ! Total budget: using Closing of C and N as a check
    
    Total_C = STRUCT_C                                                          & 
              + RESERVE_C                                                       & 
              + C_CLOSING                                       ! [mmolC/m2]
    Total_N = RESERVE_N                                                         &  
              + STRUCT_C * SN_SCr                                               & 
              + (NO3 + NH4) * depth                                             & 
              + N_CLOSING*depth                                 ! [mmolN/m2]
    Total_P = RESERVE_P                                                         & 
              + STRUCT_C * SP_SCr                                               & 
              + PO4 * depth                                                     &
              + P_CLOSING * depth                               ! [mmolP/m2]
    
    !===========================
    ! Mass Balance Equations    
    !===========================
    
    ! --- MacroAlgae ---
    dSTRUCT_C    =  (Growth                                                     &                
                    - Mortality                                                 &
                    - Regression)                               ! [mmolC/m2/d] 
    
    dRESERVE_C   =   (Photosynthesis                                            &
                    - Growth                                                    &
                    - Total_Resp                                                &
                    - Mortality * RC_SCr)                       ! [mmolC/m2/d]
    
    dRESERVE_N   =   (NH4_upt                                                   &
                    + NO3_upt                                                   &
                    - Growth * SN_SCr                                           &
                    - Mortality * RN_SCr)                       ! [mmolN/m2/d]
    
    dRESERVE_P   =   (PO4_upt                                                   &
                    - Growth * SP_SCr                                           &
                    - Mortality * RP_SCr)                       ! [mmolP/m2/d]
    
    ! dead biomass N is returned to environment directly
    ! recalculate rates per unit surface to unit volume seawater
    
    dNO3      = (- NO3_upt/depth                                                &
                 + NO3_return/depth                                             &
                 + NO3_env_exchange)                            ! [mmolN/m3/d]
    
    dNH4      = (- NH4_upt/depth                                                &
                 + NH4_return/depth                                             &
                 + NH4_env_exchange)                            ! [mmolN/m3/d]
    
    dPO4      = (- PO4_upt/depth                                                & 
                 + PO4_return/depth                                             &
                 + PO4_env_exchange)                            ! [mmolP/m3/d]
    
    dO2       = (O2_Production/depth                                            &
                 - O2_Respiration/depth                                         &
                 + OX_env_exchange)
    
    ! Closing (losses of mass to outside the box model)
    dC_CLOSING  = (  Total_Resp                                                 &
                   + Regression                                                 &
                   + Mortality                                                  &
                   + (Mortality*RC_SCr)                                         &
                   - Photosynthesis)
    
    dN_CLOSING  = ( - NO3_env_exchange                                          &
                    - NH4_env_exchange)
    
    dP_CLOSING  =   - PO4_env_exchange
    
    Net_Growth_rate = dSTRUCT_C/STRUCT_C ! [/d]

! The rates of change returned
      dConc(1)  = dSTRUCT_C  
      dConc(2)  = dRESERVE_C   
      dConc(3)  = dRESERVE_N 
      dConc(4)  = dRESERVE_P 
      dConc(5)  = dNO3  
      dConc(6)  = dNH4 
      dConc(7)  = dPO4
      dConc(8)  = dO2  
      dConc(9)  = dC_CLOSING
      dConc(10) = dN_CLOSING
      dConc(11) = dP_CLOSING 

      CALL getoutM(yout)
  
  RETURN
  END SUBROUTINE macroalgae
  
!==========================================================================
! put output variables in one vector
!==========================================================================
    
    SUBROUTINE getoutM(yout)
    INTEGER :: i
    INTEGER, PARAMETER :: nout = 69
    DOUBLE PRECISION :: yout(*), out(nout)
  
    COMMON /myoutM   /out
      DO i = 1, nout
       yout(i) = out (i)
      ENDDO       
  
    END SUBROUTINE getoutM
  
