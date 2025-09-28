! ==============================================================================
! ==============================================================================
!
! Cohort mussel model using DEB dynamics
!
! implementation: Quinten Mudde
!                 Karline Soetaert
! NIOZ - Yerseke
!
! Dynamics
! ==============================================================================
! ==============================================================================

! ==============================================================================
! ==============================================================================
!
! Dynamics of mussels in cohorts (dRESERVE, dSTRUCT, dREPROD)
!
! ==============================================================================
! ==============================================================================


SUBROUTINE mussel_dynamics (t)
USE Mussel_deb

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: t

DOUBLE PRECISION, EXTERNAL :: Arrhenius, stnorm  ! external functions (Dnorm vs stnorm)

! local variables
INTEGER :: I

DOUBLE PRECISION :: Surf, V,RSratio, RScat, spawn_event,        &
                    GrossdStruct, Dilution, gross_to_repro, DOY        
                    

DOUBLE PRECISION :: SPIM, Clearance_rate_max,                    &
   FR_SIM, Phyto_Ing, DT_Ing, SIM_Ing,                           &
   Pseudofaces_SIM,                                              &
   Ingestion_N, Food_NC, Assimilation_overhead,                  &
   Assimilation_N, Excretion_N, Faeces_SIM, Act_resp, AE
   
DOUBLE PRECISION :: Tot_maintenance, Smaintenance, Mmaintenance
DOUBLE PRECISION :: Cat_rate,  Maint_rate    


!-------------------------------------------------------------------------------
! day of the year (for spawning)
!-------------------------------------------------------------------------------
     DOY       = DMOD  (t , 365.d0)           ! day of year

!-------------------------------------------------------------------------------
! environmental conditions
!-------------------------------------------------------------------------------
     IF (.NOT. dynamic_environment) THEN
     
        PHYTO    = f_PHYTO
        DETRITUS = f_DETRITUS    
        SIM      = f_SIM
        O2       = f_O2 

        dPHYTO    = 0.D0
        dDETRITUS = 0.D0    
        dSIM      = 0.D0
        dO2       = 0.D0 
     ENDIF
     

!-------------------------------------------------------------------------------
!  Temperature and oxygen effects on physiological rates
!-------------------------------------------------------------------------------
    
    T_factor   = Arrhenius(f_Temp)   ! f_Temp is forcing function
    
    IF (O2 > crit_O2) THEN
      ox_factor  = (O2-crit_O2)**3 / ( (O2-crit_O2)**3 + (ks_O2-crit_O2)**3 )
    ELSE
      ox_factor  = 0.D0
    END IF  

!-------------------------------------------------------------------------------
! Environment-dependent DEB rates 
!-------------------------------------------------------------------------------
    

    
    !  basal respiration rate, 
    Maint_rate = maint * T_factor     ! [/d]   
    
    ! catabolism
    Cat_rate   = cat_s  * T_factor  * ox_factor             ! [mmolC/cm2/d]
    
    ! maintenance does not scale with O2, because at low O2 Carbon still is used in maintainance. 
    ! cabatolism scales with ox_factor, at low O2 no reserves can be mobilized, aka no new energy 
    ! can be burned due to ox limitaiton, this simulates "starvation" due to low oxygen
 
!-------------------------------------------------------------------------------
! Spawning times
!-------------------------------------------------------------------------------
    ! time around spawn peak 3*sd -> 4*sd
  IF(DOY >= (spawn_peak - 1.d0*spawn_dev)  .AND. DOY <=  (spawn_peak + 1.d0*spawn_dev)) THEN

    Spawn_event = stnorm(DOY, spawn_peak, spawn_dev)  ! fraction that spawns /day
    Spawn_event = MAX(0.d0,Spawn_event)
    
  ELSE 
    Spawn_event = 0.d0
  END IF
  
  ! repeat for 2nd spawn peak if this mechanism is active (>0)
  
  IF(spawn_peak_2 > 0) THEN
        IF(DOY >= (spawn_peak_2 - 1.d0*spawn_dev)  .AND.  &
        DOY <=  (spawn_peak_2 + 1.d0*spawn_dev)) THEN

            Spawn_event = Spawn_event +  stnorm(DOY, spawn_peak_2, spawn_dev)  
            ! fraction that spawns /day
            Spawn_event = MAX(0.d0,Spawn_event)
        END IF
         
  END IF  




    DO I = 1, N_Cohort  !Do loop for the following dynamics. I == Cohort for N cohorts nCohort
    
!-------------------------------------------------------------------------------
! size and tissue                ! Length [cm] 
!-------------------------------------------------------------------------------
     
     ! Reproductive tissue (gonads) and Gametes (eggs/sperm)
      Gonads(I)       = min(REPROD(I), mat_pub)

      Gametes(I)      = REPROD(I) - Gonads(I)

    ! individual Size 
      V            = STRUCT(I)/mmolC_cm3            ! volume, [cm3]             
      RSratio      = RESERVE(I)/(STRUCT(I)+1d-8)    ! [mol/mol]
      Surf         = V**(2.d0/3.d0)                 ! Surface, [cm2]
      Length(I)    = V**(1.d0/3.d0) / del_M
      
    ! Allometric function of shell surface based on Length
      Shell_surf(I)   = max(0.D0 ,exp(shell_allo_a + (Length(I) * shell_allo_b)))
      
      
      
!-------------------------------------------------------------------------------
! Physiological Clearance, Filtration, Ingestion, Assimilation, (Pseudo)faeces
!-------------------------------------------------------------------------------
  
  ! Surface dependent clearance rate inhibited by all suspended particles 
  
  SPIM = f_Sim ! shorthand for suspended particulate inorganic matter (causes filtration inhibition)
  
  Clearance_rate_max = max_clear * T_factor * ox_factor ! [m3/d/cm2]
  
  Clearance_rate(I) = Surf *   Clearance_rate_max / (1.d0 +             &
    (( PHYTO*Clearance_rate_max)/max_phyto_filt)    +      & ! phyto & DT separated
    ((DETRITUS*Clearance_rate_max)/max_DT_filt)    +                    &
    ((SPIM*Clearance_rate_max)/max_spm_filt))                         ! [m3/d]
 
  ! Paricle specific filtration rates based on actualised Clearing rate
 
  FR_Phyto(I) = Clearance_rate(I)*PHYTO
  FR_DT(I) = Clearance_rate(I)*DETRITUS
  FR_SIM = Clearance_rate(I)*SPIM
  
  ! Particle specific Ingestion based on binding propability (rho) and max ingestion
 
  Phyto_Ing = rho_Phyto * FR_Phyto(I) /                                              &
    (1.d0 + ((rho_phyto*FR_Phyto(I))/max_Phyto_ing) + ((rho_DT*FR_DT(I))/max_DT_ing) +  & 
    ((rho_SIM*FR_SIM)/max_SIM_ing))
    
  DT_Ing = rho_DT * FR_DT(I) /                                                       &
    (1.d0 + ((rho_phyto*FR_Phyto(I))/max_Phyto_ing) + ((rho_DT*FR_DT(I))/max_DT_ing) +  & 
    ((rho_SIM*FR_SIM)/max_SIM_ing))
  
  
    
  SIM_Ing = rho_SIM * FR_SIM /                                                       &
    (1.d0 + ((rho_phyto*FR_Phyto(I))/max_Phyto_ing) + ((rho_DT*FR_DT(I))/max_DT_ing) +  & 
    ((rho_SIM*FR_SIM)/max_SIM_ing))
    
    

  ! Pseudofaeces are filtered & not-ingested particles (clusered as detrital pallets)
  
  Pseudofaeces_C(I) = FR_Phyto(I) - Phyto_Ing + FR_DT(I)  - DT_Ing
  Pseudofaces_SIM = FR_SIM - SIM_Ing
  
  ! Ingestion, split up in Carbon and Nitrogen ingestion to calculate assimilation

  Ingestion(I) = Phyto_Ing+DT_Ing ! carbon mass
  Ingestion_N = Phyto_Ing*nc_Phyto + DT_Ing*nc_DT ! nitrogen mass
  Food_NC = Ingestion_N/(Ingestion(I) + 1d-8) !N/C ratio of all ingested organic particles
  
  ! Assimilation is based on Ingestion and N/C ratio of ingested organics
  ! N/C ratio of food closer to animal tissue results higher effective aeff (assim efficiency)
  
  AE = min(1.d0, Food_NC/nc_tissue)*aeff
  ! base assimilation efficiency (aeff) defines feeding respiration in DEB
  ! Here, we include addional assimilation penalty for the composition of the ingested food
  ! Current formulation assumes FoodNC cannot exceed nc_tissue
  ! ingested SIM only inhibits filtration and ingestion but not assimilation

  Assimilation(I) = (Ingestion(I) * AE ) 
  Assimilation_N = Assimilation(I) * nc_tissue
  
  ! Everything left is excreted
  Faeces_production(I) = Ingestion(I) - Assimilation(I) ! fecal pallets in C
  Excretion_N = Ingestion_N - Assimilation_N ! = aeff*Ing_N in most cases
  Faeces_SIM = SIM_Ing ! inorganic part of fecal pallets
  
!-------------------------------------------------------------------------------
! Respiration, catabolism [mmolC/ind/d] 
!-------------------------------------------------------------------------------
     ! part of catabolism affecting R:S ratio
      RScat         = Cat_rate * RSratio * Surf       
      
     ! Maintenance of structure and maturity (gonads))
      Smaintenance  = Maint_rate * STRUCT(I)          
      Mmaintenance  = Maint_rate * Gonads(I)              
      
      Tot_maintenance = Smaintenance + Mmaintenance 

      ! structural volume growth
      GrossdStruct = kappa*RScat - Smaintenance   
      
      ! changes second condition in if statement (games now provede all given ther is enough)
      ! this change was made because otherwisze struct_growth was always slightly neg during regresison
      
      IF (GrossdStruct < 0.d0 .AND. Gametes(I)>(-2.d0*GrossdStruct)) THEN ! maintenance paid by gametes
      ! regression == maintenance paid by GAMETES
      
        Regression(I) = -GrossdStruct ! *(1.d0 - exp(-(Gametes(I)/2.d0)**2.d0)) ! ?
      ELSE
        Regression(I) = 0.d0
      END IF  

!-------------------------------------------------------------------------------
! Production of structural C, including costs+dilution 
!-------------------------------------------------------------------------------
      Struct_growth(I) = (GrossdStruct + Regression(I)) / (1.d0 + eg + kappa*RSratio) 

      Dilution      = RSratio  * Struct_growth(I)
      
      ! catabolism, total reserves catabolised
      Catabolism(I) = RScat - Dilution       
      
!-------------------------------------------------------------------------------
! Reproduction
!-------------------------------------------------------------------------------
      Gross_to_repro  = (1.d0-kappa)*Catabolism(I) -Mmaintenance 
      
      ! Production of maturity C, including costs 
      Reprod_growth(I) = (Gross_to_repro - Regression(I))/(1.d0 + eg)
      
      ! yearly spawning empties the Gametes buffer


      Spawn_rate(I)   = Gametes(I) * spawn_eff * Spawn_event   ! [mmolC/ind/day]

      
!-------------------------------------------------------------------------------
!   MASS BALANCE EQUATIONS 
!-------------------------------------------------------------------------------
      
      dSTRUCT(I)     = Struct_growth(I)
      dRESERVE(I)    = Assimilation(I) - Catabolism(I) 
      dREPROD(I)     = Reprod_growth(I) - Spawn_rate(I)
      
!-------------------------------------------------------------------------------
! Outputs      
!-------------------------------------------------------------------------------

      ! in original DEB, the MATURITY buffer is 'complexity' of the structure
      ! here this is described by Gonads(I), can be seen as a one-time investment, 
      ! but the biomass is in STRUCT
      ! and therefore not added to weight. REPROD(I) vs Gametes(I)
      
      ! Here, we include the carbon weight of GONADS in the total weights, 
      ! as this is assumed to be live tissue
      
      ! Energy reserves & struct could have different mass densities 
      ! (g/mmolC STRUCT =/= g/mmolC RESERVE/GAMETE/REPROD )
      ! all these aspects are simplified, and not taken into account..?
      
      C_ind(I)       = RESERVE(I) + STRUCT(I) + REPROD(I) !+ Gametes(I)
      
      ! estimation of wet (meat) weight based on all carbon pools
      
      WW_ind(I)      = C_ind(I) * (gDW_cm3/mmolC_cm3) * gWW_gDW 
      
      
      Filtration_C(I) = FR_Phyto(I) + FR_DT(I)

      ! all carbon losses & overheads are assumed to be dissipated and respired
      Respiration(I) = eg*(Struct_growth(I) + Reprod_growth(I)) +    &
      Tot_maintenance ! + Act_resp
      
      VB_grow(I) = (3.d0/(STRUCT(I)+1d-8)) * Struct_growth(I) !formulation from Ben Martin 2013
      
    END DO


END SUBROUTINE mussel_dynamics

! ==============================================================================
! ==============================================================================
!
! Dynamics of mussel density in cohorts (dPOPULATION)
!
! ==============================================================================
! ==============================================================================

SUBROUTINE population_dynamics

USE Mussel_deb

IMPLICIT NONE

! local variables
INTEGER          :: I
DOUBLE PRECISION :: total_mort, overcrowd, crowdingvar, crowd_mortality, &
                    size_mortality, partial_starvation, act_mortality, RSratio2

    DO I = 1, N_Cohort
    
    ! changed env_mortality into act_mortality, 
    ! activity mortality scales with T_factor at lower temps (lower metabolism = lower moratalty)
    ! activity mortality reverse scales with T_factor at highter temps (heat mortality)
    
!-------------------------------------------------------------------------------
! Population mortality - depends on oxygen, starvation & crowding
!-------------------------------------------------------------------------------
   ! Mortality due to crowding formulation from ERSEM C based benthic overcrowding
   
   overcrowd = Total_shell_cover - crowd_lim
   crowdingvar = (overcrowd**2.d0)/(overcrowd+crowd_lim)
      
      IF (overcrowd > 0.D0) THEN
        crowd_mortality = crowdingvar/(crowdingvar+crowd_lim)
      ELSE
        crowd_mortality = 0.D0
      END IF
      
    
   
!   mort_env = mort_base + (  (mort_juv-mort_base) * (1.d0 - Length(I) / (Length(I)+recr_length))   )
   
   ! assume that 50% of individuals have less food than the average (and can therefore have starvation)
   ! this links food to population mortality and is therefore very important
   ! rs ratio below [e_dens_trh] <- partial starvation starts
   
   RSratio2      = RESERVE(I)/(STRUCT(I)+1d-8) ! this is not a global term so calculate again here
   
   
      IF (RSratio2 < e_dens_trh) THEN 
      
         partial_starvation = mort_starv * (e_dens_trh - RSratio2) / e_dens_trh 
         ! mort_starv is additional per capita mortality under complete starvation (R/S = 0)
         
      ELSE 
      
        partial_starvation = 0.d0
         
      END IF 
   
   ! assume reference mort_base is for reference temperature (T_factor = 1)
   ! ---assume additional mortality can be due to oxygen limitation---
   ! density dependence is caused by crowd mortality (individuals grow on top of eachother)
   
   IF (f_Temp > T_r .AND. T_factor < 1) THEN
   !env_mortality = mort_base*(1.d0 + (1.d0 - T_factor)) ! base + increase in heat
   act_mortality = mort_base/T_factor ! base mortality increased by heat stress
   ELSE
   act_mortality = mort_base * T_factor ! base mortality scaled by metabolic activity
   END iF
   
   ! shell size based mortality, additional per capita mortality for juveniles (unitl recruitment)
   
   size_mortality = mort_juv * (1.d0-(Length(I)**3 / (Length(I)**3 + recr_length**3)))
   
   ! mortality due to oxygen limitation removed, assumed effect via individual physiology
   
   ! sum all mortality factors for total per capita mortality
   
   total_mort = act_mortality + crowd_mortality + partial_starvation + size_mortality
   

                   
      
      ! additional layer of starvation when average animal cannot pay structural maintenance
      
      IF (Struct_growth(I) < 0.d0) THEN
         total_mort   = total_mort + mort_starv
      END IF 
      
      
      Mortality(I)  = total_mort*POPULATION(I) ![ind/m2/d]
       
      dPOPULATION(I) = -Mortality(I)
      
      ! dynamics should stop when a certain cohort is empty
      IF (POPULATION(I) < 1d-3) THEN 
      dPOPULATION(I) = 0.d0
      dSTRUCT(I) = 0.d0
      dRESERVE(I) = 0.d0
      dREPROD(I) = 0.d0
      END iF
      
      
    END DO
    
END SUBROUTINE population_dynamics

! ==============================================================================
! ==============================================================================
!
! Summed values (over all cohorts)
!
! ==============================================================================
! ==============================================================================

SUBROUTINE summedRates

USE Mussel_deb

IMPLICIT NONE

! local variables

INTEGER :: I

DOUBLE PRECISION :: Total_C_m2_FIELD



    Total_C_m2        = 0.D0
    Total_POP_m2      = 0.D0
    Total_ingestion   = 0.D0
    ! split up ingestion in phyto+DT filtration & pseudofaeces
    Total_FR_Phyto    = 0.D0
    Total_FR_DETRITUS = 0.D0
    Total_psfaeces_C = 0.D0
    Total_faeces_C      = 0.D0
    Total_respiration = 0.D0
    Total_spawn       = 0.D0
    Total_shell_cover = 0.D0
    
    Total_POP_m2_FIELD = 0.D0
    Total_C_m2_FIELD = 0.D0
    
    DO I = 1, N_Cohort
      C_m2(I)           = C_ind(I)      *POPULATION(I)  ! mol C/m2
      Total_C_m2        = Total_C_m2   + C_m2(I)
      Total_POP_m2      = Total_POP_m2 + POPULATION(I)
      Total_FR_Phyto    = Total_FR_Phyto + FR_Phyto(I)*POPULATION(I)
      Total_FR_DETRITUS = Total_FR_DETRITUS + FR_DT(I)*POPULATION(I)
      Total_psfaeces_C  = Total_psfaeces_C + Pseudofaeces_C(I)*POPULATION(I)
      Total_ingestion   = Total_ingestion       + Ingestion(I)*POPULATION(I)
      Total_faeces_C    = Total_faeces_C  + Faeces_production(I)*POPULATION(I)
      Total_respiration = Total_respiration   + Respiration(I)*POPULATION(I)
      Total_spawn       = Total_spawn          + Spawn_rate(I)*POPULATION(I)
      Total_shell_cover = Total_shell_cover + Shell_surf(I)*POPULATION(I)
      
      IF (Length(I) > length_trh) THEN
      Total_POP_m2_FIELD = Total_POP_m2_FIELD + POPULATION(I)
      Total_C_m2_FIELD = Total_C_m2_FIELD + C_m2(I)
      END IF
      
    END DO
    
    Total_DW_m2  = Total_C_m2 /mmolC_cm3*gDW_cm3
    Total_WW_m2  = Total_DW_m2 * gWW_gDW
    Total_WW_m2_FIELD = Total_C_m2_FIELD/mmolC_cm3*gDW_cm3*gWW_gDW

END SUBROUTINE summedRates


! ==============================================================================
! ==============================================================================
!
! Dynamics of larvae  (dLARV_DENS, dLARV_BIOM )
!
! ==============================================================================
! ==============================================================================

SUBROUTINE larvae_dynamics (t)
USE Mussel_deb

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: t
DOUBLE PRECISION, EXTERNAL :: Arrhenius

! local variables
INTEGER :: I, new ,Year, second_coh 

DOUBLE PRECISION :: DOY
DOUBLE PRECISION :: dL_Dens, dL_Biom, dS_Biom, dR_Biom,          &
                    L_excessW, L_DENS, recruits,                 &
                    weight_settle_full, Substrate_suitability,   &
                    Set_rate, ks_settlement_field, T_factor_l,   &
                    S_fact, Set_speed

!-------------------------------------------------------------------------------
! Cohorts affected by larval settling
!-------------------------------------------------------------------------------
     DOY       = DMOD  (t , 365.d0)           ! day of year
     Year      = FLOOR (t / 365.d0)            ! year
     
     ! make 2 spawning cohorts per year
     
    IF (DOY <  spawn_peak_2-3.d0*spawn_dev) THEN
    
           second_coh = 0 ! use as a boolean
           
      ELSE 
      
           second_coh = 1 ! use as a boolean
           
    END IF    
      
     
    IF (spawn_peak_2 > 0.d0) THEN     
     
           new = NINT(DMOD(DBLE(Year), (DBLE(N_Cohort) / 2.d0) )*2 +1 + second_coh) 
     
      ELSE 
  
           new = MOD(Year, N_Cohort) +1 ! we also put year as integer here
     
    END IF   

!-------------------------------------------------------------------------------
! Remove biomass and density from the oldest cohort that will get settlers - 
! do this for DOY < spawn
!-------------------------------------------------------------------------------
    
    IF (DOY < (spawn_peak-3.d0*spawn_dev) .AND. Year >= N_Cohort .AND. N_Cohort > 1) THEN   
    ! removal rate of 10%/day  
      dSTRUCT(new)     = - 0.1d0 * STRUCT(new)
      dRESERVE(new)    = - 0.1d0 * RESERVE(new)
      dREPROD(new)     = - 0.1d0 * REPROD(new)
      dPOPULATION(new) = - 0.1d0 * POPULATION(new)
    END IF 
    
!-------------------------------------------------------------------------------
! Larvae
!-------------------------------------------------------------------------------

!------------------------------
! Growth & Mortality of Pelagic Larva
!-----------------------------

!----- 
! Growth rate of Larva
!----

! effect of temperature (assume same arrhenius effect on adults and juviniles)
    
T_factor_l   = Arrhenius(f_Temp)   ! f_Temp is forcing function

! Larva growth : optimal larvae growth between 20 & 24 dC (Lazo & Pita, 2011) 
! can compare arrhenius curve with these physiological data
! Q10 growth analysis: 17-20 dgC Q10 = 2.7 // 17-24 dgC Q10 = 1.569 // 20-24 dgC Q10 = 0.957 = cosidered 1


IF(L_weight <= (2.d0 * weight_settle_full)) THEN

    IF(LARV_DENS > 2.d0) THEN ! growth ends before mortality for computational reasons

       
       larv_growth_rate = growth_larva*T_factor_l
       
!      IF (f_Temp >= (T_r-273.15d0)) THEN
!   
!          larv_growth_rate = growth_larva
!     ! *q10_LgU**((f_Temp - (T_r-273.15d0))/10.D0) <- removed (1q0 was 0.95, basically 1)
!       ELSE 
!      
!          larv_growth_rate = growth_larva*q10_LgL**((f_Temp - (T_r-273.15d0))/10.D0)
!      
!       END IF
  
    ELSE
     
       larv_growth_rate = 0.d0 !was 1d-5 (but this caused crashes at low dens numbers)
       
    END IF
     
ELSE
     
     larv_growth_rate = 0.d0 ! 0.d0 vs -0.1d0
     !Total_spawn = 0.d0
     
END IF

!-----------
! Mortality rate of Larva
!----------

! Larva mortality : most important factor for larva dynamics
    
 IF(LARV_DENS > 1.d0) THEN  ! was  1d-3
 
    IF (f_Temp > T_r .AND. T_factor_l < 1.d0) THEN
  
        larv_mort_rate = mort_larva/T_factor_l ! base + increase  in heat
      
    ELSE
  
        larv_mort_rate = mort_larva*T_factor_l ! decrease mortality in cooler waters
        
    END iF
    
    
    
    
!    IF (f_Temp >= (T_r-273.15d0)) THEN
!   
!       larv_mort_rate = mort_larva*q10_Lm**((f_Temp - (T_r-273.15d0))/10.d0)
!      
!    ELSE 
!      
!       larv_mort_rate = mort_larva*(1.d0/q10_Lm)**((f_Temp - (T_r-273.15d0))/10.d0)
!      
!    END IF
    
 ELSE
 
    larv_mort_rate = 1d-3! was0.d0
    
 END IF
    


!-------------
! settlement
!-------------
    
! mussels are settled form a pool of Competent Pediveligers (comp_pedi), 
! defined as larvae that have reached a certain individual biomass. 
! To determine the comp_pedi pool, we set a threshold and smoothly fill increase this value
! this smoother is needed to help the solver, if 2x the threshold is reached, all larvae are competent

! individual biomass is only calculated if a larva population exists (larv dens > 1)
! maximum larva indiviual weight is kept at 2x the threshold to ease compuations

!------ 
! Indiviudal larva weight (L_weight [mmolC/ind])
!------
    
    weight_settle_full = weight_settle*2.d0

    IF (LARV_DENS > 1.d0) THEN
       L_weight  = MIN((weight_settle_full*2.d0),(LARV_BIOM/(LARV_DENS + 1d-4)))
       ! exploding L_weight causes error
    ELSE 
       L_weight = 0.d0
    END IF
    
!---
! Competent pediveligers density (comp_pedi [ind/m3])
!---
  
    IF (L_weight > weight_settle_full) THEN
    
        comp_pedi = LARV_DENS
    
    ELSE 
 
        IF (L_weight >= weight_settle .AND. L_weight <= weight_settle_full) THEN
    
            comp_pedi = (1.d0 - (L_weight - weight_settle_full)**2.d0         &
               / ((L_weight - weight_settle_full)**2.d0                       &
              +(weight_settle*0.5d0)**2.d0 )) * LARV_DENS
    
        ELSE
    
            comp_pedi = 0.d0
       
        END IF   
    
    END IF

!---
! Settlement Rate [ind/m2/day] & effect of substrate for settlement
!---  
! Settlement rate is determined by a supply of competent/suitable pediveligers in the water column 
! the relationship used is a monod (experimental), with maximum rate and half supply value roughly 
! estimated from literature data (toupoint et al., settlement of M. edulis larva if field)
! to implement the effect of substrate (& current) we can let the substrate determine the ks value
! substrate suitability is promoted by conspecific shells or coarser substrate (very simple relationship)
! substr_factor == average grain size [cm diameter] (1 = rocky pebbles == suitable)
! shell cover is expressed in cm2, so (shellcover/10.000) is fraction cover per m2

! effect of shear stress on larva settlement? Tc can be set to 0.08 [N/m2] as critical shear stress for settlement
! ref: Daraio et al., 2010 (settlement model in a river)
! sinking speed can be calculated with stokes law: Om = -(SF/18v)*(SG-1)*g*d^(1/2)
! SF(shapefactor) = c(ab)^(-1/2) ~ 0.5-0.8 
! SG(SpecificGravity) ~ 1.2 (Daraio et al., 2010)
! d(diameter[m]) ~ 250um =  0.00025 [m] (grows with size)
! v(viscosity of water [m/s] (mu)) = 1e-6
! g(gravitational acceleration[m/s]) = 9.81
! resulting fall velocity/sinking speed = -(0.5/(18*1e-6))*(1.25-1)*9.81*0.00025^(1/2)
! setting v = 1 (?) this results in 0.00086 [m/s] -> 74 [m/d] sinking without turbulence
! the applied sinking speed falls off linearly from 74 to 0 with critical threshold of 0.08N/m2 
! we hereby assumet that in the 0D scenario, shear stress is also a proxy for upwards velocity
! stokes law for terminal (settling) velocity
! vSet = (g*d^2)/(18*mu) * (GS - 1)  === 588 [m/d sinking speed]

IF (substr_factor == 1.d0) THEN
      S_fact = 1.d0
  ELSE IF (substr_factor == 2.d0) THEN
      S_fact = 0.5d0
  ELSE IF (substr_factor == 3.d0) THEN
        S_fact = 0.1d0
  ELSE IF (substr_factor == 4.d0) THEN
        S_fact = 0.01d0
  ELSE 
        S_fact = 1.d0 ! assume that substrate is perfectly suitable if no proper par is passed
END IF

! Substrate_suitability = MIN(1.d0, (SQRT(substr_factor) + (Total_shell_cover/10000.d0) ))

Substrate_suitability = MIN(1.d0, (S_fact + (Total_shell_cover/10000.d0)) )

! PARS ==> v_settle (larva sinking velocity [m/s]) & tc (critical shear stress for settling [N/m2])
! FORC ==> TAUBOT (live shear stress forcing [N/m2])
IF (f_tau > tc) THEN
    Set_speed = 0.d0 ![m/d] (no settlement if shearstress is above critical threshold for settling)
ELSE  
    Set_speed = v_settle  * (1.d0 - (f_tau/tc)) ![m/d]
END IF

Set_rate = MIN((Set_speed * comp_pedi),max_set_rate* Substrate_suitability) ![ind/m2/d]

! now multiply suitability with max settle rate!
 
 
! the amount of supply needed to reach maximum settle rate increases with worse substrate (<1)


!  ks_settlement_field = ks_settlement / Substrate_suitability 
  
! monod equation for settlement rate
! this is then increased by the settle speed [max_set_rate, [ind/m2/d]], which can be changed in parms

!  Set_rate = max_set_rate * comp_pedi / (comp_pedi + ks_settlement_field) ! [ind/m2/d]
 
  
  

    
!---------------------------------------------------------------
! Dynamics of Larva Density [ind/m3] & Biomass [mmolC/m3] Pools
!---------------------------------------------------------------

  dLARV_DENS =    Total_spawn/2.d0/weight_egg/depth*(1.d0-mort_egg)   & !half female ratio
                - larv_mort_rate*LARV_DENS                            &
                - Set_rate/depth
                
                
! f_PHYTO/(f_PHYTO+ks_larva) do not compete for bottom phyto therefore raw forc ?
! Larva do not partake in competition of mussels (raw forcing f_phyto instead of dynamic PHYTO)
 
 IF(LARV_DENS>=1.d0) THEN
 
    dLARV_BIOM =  Total_spawn/depth/2.d0*(1.d0-mort_egg)                    & !half female ratio
                + larv_growth_rate * LARV_BIOM * f_PHYTO/(f_PHYTO+ks_larva) &  
                - larv_mort_rate*LARV_BIOM                                  &
                - Set_rate*LARV_BIOM/(LARV_DENS + 1d-4)/depth
    
 ELSE 
    IF(LARV_BIOM > weight_egg) THEN
    
       dLARV_BIOM = -0.1d0*LARV_BIOM ! reset larv biomass when density is below 1
       
    ELSE
 
       dLARV_BIOM = 0.d0
       
    END IF   
    
 END IF
 
    
    Settle_rate(:) = 0.D0
    
!-------------------------------------------------------------------------------
! New cohort dynamics
!-------------------------------------------------------------------------------

! Density increase due to settling
! recruits are proportion settlers that survive the early settlement mortality   

    Settle_rate(new) = Set_rate                  ! output variable [indivs/m2/day]
    
    
    dPOPULATION(new) = dPOPULATION(new) + Set_rate  

! Weight change due to settling, including dilution of existing 
! from dB/dt = d(N*W)/dt = N dW/dt + W dN/dt
! we have: dW/dt = (dB/dt-W dN/dt)/N

! previous: dL_Biom = recruits*LARV_BIOM

    dL_Biom = Set_rate*LARV_BIOM/(LARV_DENS + 1d-4) ! Total biomass increase
    dS_Biom = dL_Biom * (1.d0-pReserve_settle)      ! increase struct biomass
    dR_Biom = dL_Biom *  pReserve_settle            ! increase reserve biomass

    L_Dens    = POPULATION(new)                 ! Total density in new cohort
    
    
    
    
 !   IF (L_weight >= weight_settle .AND. dL_Dens > 1.D0) THEN
 
 !  we have: dW/dt = (dB/dt-W dN/dt - dB/dt-)/N   ?
    
    IF (L_Dens > 1.D0) THEN   ! (must be 1.d0!, 0 cause neg growth)
      
      dSTRUCT(new) = dSTRUCT(new)                     &
              + (dS_Biom - STRUCT(new)*Set_rate)/(L_Dens + 1d-8) 
              
      dRESERVE(new)    = dRESERVE(new)                &
              + (dR_Biom - RESERVE(new)*Set_rate)/(L_Dens + 1d-8)
              
      dREPROD(new)     = dREPROD(new)                 &
              + (0       - REPROD(new)*Set_rate) /(L_Dens + 1d-8)
    
    
    
    ELSE
      dSTRUCT(new)     = max(0.D0, dS_Biom/(Set_rate+1d-8)) ! was just dL_Dens
      dRESERVE(new)    = max(0.D0, dR_Biom/(Set_rate+1d-8))
      dREPROD(new)     = 0.D0
    END IF
    
    
   
END SUBROUTINE larvae_dynamics

! ==============================================================================
! ==============================================================================
!
! Dynamics of environmental variables (dPHYTO, dDETRITUS, dSIM, dO2)
!
! ==============================================================================
! ==============================================================================

SUBROUTINE environment_dynamics
USE Mussel_deb

IMPLICIT NONE

DOUBLE PRECISION :: renew_rate
!-------------------------------------------------------------------------------
! Environmental dynamics only when dynamic_environment is TRUE
!-------------------------------------------------------------------------------

    IF (.NOT. dynamic_environment) RETURN
    
    renew_rate = water_renew/(depth*1.*1.) ! [/d] !box volue = depth*1*1
    ! calculate relative water renewal based on daily renewal [m3/d] & box volue [m3] 
    
    !dPHYTO    = renew_rate*(f_PHYTO    - PHYTO)                    &
    !           - Total_Ingestion * p_phyto/depth
    
    dPHYTO     = renew_rate*(f_PHYTO - PHYTO)                       &
                - Total_FR_Phyto/depth
               
    !dDETRITUS = renew_rate*(f_DETRITUS - DETRITUS)                 &
    !           - Total_Ingestion * (1.d0-p_phyto)/depth
               
    dDETRITUS = renew_rate*(f_DETRITUS - DETRITUS)                  &
               - Total_FR_DETRITUS/depth
               
    dSIM      = renew_rate*(f_SIM - SIM)
    dO2       = renew_rate*(f_O2 - O2) - 1.070312 * Total_respiration/depth
    

END SUBROUTINE environment_dynamics