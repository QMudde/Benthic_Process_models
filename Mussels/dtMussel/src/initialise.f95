! ==============================================================================
! ==============================================================================
!
! Cohort mussel model using DEB dynamics
!
! implementation: Quinten Mudde
!                 Karline Soetaert
! NIOZ - Yerseke
!
! Initialisation subroutines
! ==============================================================================
! ==============================================================================

!==========================================================================
! initialise common parameter block - calls allocate function
!==========================================================================

SUBROUTINE initmuspar (steadyparms)
USE Mussel_dim
       
IMPLICIT NONE
EXTERNAL steadyparms

DOUBLE PRECISION parms(nparms)
COMMON /musselparms/parms

     CALL steadyparms(nparms, parms)
     CALL initialise_deb()

RETURN
END SUBROUTINE initmuspar

!==========================================================================
! Initialise the forcing function common block
!==========================================================================

SUBROUTINE initmusforc (steadyforcs)
USE Mussel_dim
      
IMPLICIT NONE
EXTERNAL steadyforcs

DOUBLE PRECISION forcs(nforcs)
COMMON /musselforcs/forcs

    CALL steadyforcs(nforcs, forcs)
       
RETURN
END SUBROUTINE initmusforc

! ==============================================================================
! Initialise some booleans
! ==============================================================================

SUBROUTINE initialise_deb

USE Mussel_deb
IMPLICIT NONE

INTEGER :: I, nout_0D

   dynamic_environment = FLOOR(Env + 0.1d0) > 0d0 ! interaction with environment
   
   N_Cohort    = FLOOR(Num + 0.1d0)
   nout_Cohort = 16 * N_Cohort    ! number of output variables defined in cohort
 
   includePopulation = .TRUE.
   nout_0D = nout
 
   IF (N_Cohort < 1) THEN  ! this is so when only DEB is calculated
     includePopulation = .FALSE.
     N_Cohort          = 1
     nout_Cohort       = 15 * N_Cohort   
     nout_0D           = 0
   END IF 
 
 NumOutput = nout_Cohort + nout_0D
 
 CALL allocate_deb

END SUBROUTINE initialise_deb
 
! ==============================================================================
! Allocate memory
! ==============================================================================

SUBROUTINE allocate_deb

USE Mussel_deb
IMPLICIT NONE

INTEGER :: I

! ------------------------------------------------------------------------------
! State variables
 
 IF (ALLOCATED(RESERVE))           DEALLOCATE(RESERVE)
 ALLOCATE (RESERVE        (N_Cohort))

 IF (ALLOCATED(STRUCT))            DEALLOCATE(STRUCT)
 ALLOCATE (STRUCT         (N_Cohort))

 IF (ALLOCATED(REPROD))            DEALLOCATE(REPROD)
 ALLOCATE (REPROD         (N_Cohort))

 IF (ALLOCATED(POPULATION))        DEALLOCATE(POPULATION)
 ALLOCATE (POPULATION     (N_Cohort))

! ------------------------------------------------------------------------------
! Derivatives

 IF (ALLOCATED(dRESERVE))          DEALLOCATE(dRESERVE)
 ALLOCATE (dRESERVE       (N_Cohort))

 IF (ALLOCATED(dSTRUCT))           DEALLOCATE(dSTRUCT)
 ALLOCATE (dSTRUCT        (N_Cohort))

 IF (ALLOCATED(dREPROD))           DEALLOCATE(dREPROD)
 ALLOCATE (dREPROD        (N_Cohort))

 IF (ALLOCATED(dPOPULATION))       DEALLOCATE(dPOPULATION)
 ALLOCATE (dPOPULATION    (N_Cohort))

! ------------------------------------------------------------------------------
! output variables 

 IF (ALLOCATED(Ingestion))         DEALLOCATE(Ingestion)
 ALLOCATE (Ingestion      (N_Cohort))
 
 ! 3 new allcations to calculate filtered C budged
 
 IF (ALLOCATED(Pseudofaeces_C))        DEALLOCATE(Pseudofaeces_C)
 ALLOCATE (Pseudofaeces_C     (N_Cohort))
 
 IF (ALLOCATED(FR_Phyto))        DEALLOCATE(FR_Phyto)
 ALLOCATE (FR_Phyto     (N_Cohort))
 
 IF (ALLOCATED(FR_DT))        DEALLOCATE(FR_DT)
 ALLOCATE (FR_DT     (N_Cohort))
 
 !
 
 ! New allocation to calculate total physical shell cover (for crowding morality)
 IF (ALLOCATED(Shell_surf))        DEALLOCATE(Shell_surf)
 ALLOCATE (Shell_surf     (N_Cohort))
 !
 
 ! New allocation to calculate Von Bertalanffy growth rate 
 IF (ALLOCATED(VB_grow))        DEALLOCATE(VB_grow)
 ALLOCATE (VB_grow     (N_Cohort))
 !
 
 ! New allocation to calculate WW per individual for validation
 IF (ALLOCATED(WW_ind))        DEALLOCATE(WW_ind)
 ALLOCATE (WW_ind     (N_Cohort))
 !
 
 !new
 ! New allocation to calculate WW per individual for validation
 IF (ALLOCATED(Filtration_C))        DEALLOCATE(Filtration_C)
 ALLOCATE (Filtration_C     (N_Cohort))
 
 !new
 ! New allocation to calculate WW per individual for validation
 IF (ALLOCATED(Clearance_rate))        DEALLOCATE(Clearance_rate)
 ALLOCATE (Clearance_rate     (N_Cohort))
 
 
 IF (ALLOCATED(Catabolism))        DEALLOCATE(Catabolism)
 ALLOCATE (Catabolism     (N_Cohort))
 
 IF (ALLOCATED(Struct_growth))    DEALLOCATE(Struct_growth)
 ALLOCATE (Struct_growth  (N_Cohort))
 
 IF (ALLOCATED(Faeces_production)) DEALLOCATE(Faeces_production)
 ALLOCATE (Faeces_production(N_Cohort))
 
 IF (ALLOCATED(Regression))        DEALLOCATE(Regression)
 ALLOCATE (Regression     (N_Cohort))
 
 IF (ALLOCATED(Assimilation))      DEALLOCATE(Assimilation)
 ALLOCATE (Assimilation   (N_Cohort))
 
 IF (ALLOCATED(Respiration))       DEALLOCATE(Respiration)
 ALLOCATE (Respiration    (N_Cohort))
 
 IF (ALLOCATED(Mortality))         DEALLOCATE(Mortality)
 ALLOCATE (Mortality      (N_Cohort))
 
 IF (ALLOCATED(Spawn_rate))        DEALLOCATE(Spawn_rate)
 ALLOCATE (Spawn_rate     (N_Cohort))
 
 IF (ALLOCATED(Reprod_growth))     DEALLOCATE(Reprod_growth)
 ALLOCATE (Reprod_growth  (N_Cohort))
 
 IF (ALLOCATED(Settle_rate))       DEALLOCATE(Settle_rate)
 ALLOCATE (Settle_rate    (N_Cohort))
 
 IF (ALLOCATED(Gonads))            DEALLOCATE(Gonads)
 ALLOCATE (Gonads         (N_Cohort))
 
 IF (ALLOCATED(Gametes))           DEALLOCATE(Gametes)
 ALLOCATE (Gametes        (N_Cohort))
 
 IF (ALLOCATED(Length))            DEALLOCATE(Length)
 ALLOCATE (Length         (N_Cohort))
 
 IF (ALLOCATED(C_ind))             DEALLOCATE(C_ind)
 ALLOCATE (C_ind          (N_Cohort))
 
 IF (ALLOCATED(C_m2))              DEALLOCATE(C_m2)
 ALLOCATE (C_m2           (N_Cohort))

END SUBROUTINE allocate_deb


