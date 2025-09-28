! ==============================================================================
! ==============================================================================
!
! Cohort mussel model using DEB dynamics
!
! implementation: Quinten Mudde
!                 Karline Soetaert
! NIOZ - Yerseke
!
! Main subroutines musselmod and musselcohort
! ==============================================================================
! ==============================================================================

SUBROUTINE musselmod(neq, t, y, dy, yout, ip)

USE Mussel_dim
USE Mussel_deb

IMPLICIT NONE

 INTEGER          :: neq   ! number of equations, 3 
 INTEGER          :: ip(*) ! integer vector passed from R

 DOUBLE PRECISION :: t               ! current time
 DOUBLE PRECISION :: y(neq), dy(neq) ! state variables and derivatives 
 DOUBLE PRECISION :: yout(*)         ! double precision vector passed from/to R

 INTEGER :: I

  if (ip(1) < numOutput) call rexit("nout is not large enough")
  if (neq .NE. 3) call rexit("this model should have 3 state variables")
  if (N_Cohort .NE. 1) call rexit("this model should have only one cohort")

! ------------------------------------------------------------------------------

! extract the state variables from y (only one cohort)
    RESERVE(1)     = y (1)
    STRUCT(1)      = y (2)
    REPROD(1)      = y (3)

! ------------------------------------------------------------------------------
! calculate the derivatives for 

  CALL mussel_dynamics(t)

! ------------------------------------------------------------------------------
! save the derivatives in dy

    dy (1) = dRESERVE(1)
    dy (2) = dSTRUCT(1)  
    dy (3) = dREPROD(1)  

! ------------------------------------------------------------------------------
! save the output in yout

  CALL MusselOutput(yout)
  
END SUBROUTINE MusselMod

! ==============================================================================
! Main subroutine for mussels in a cohort
! ==============================================================================

SUBROUTINE musselcohort(neq, t, y, dy, yout, ip)

USE Mussel_dim
USE Mussel_deb

IMPLICIT NONE

 INTEGER          :: neq   ! number of equations, 
 INTEGER          :: ip(*) ! integer vector passed from R

 DOUBLE PRECISION :: t               ! current time
 DOUBLE PRECISION :: y(neq), dy(neq) ! state variables and derivatives 
 DOUBLE PRECISION :: yout(*)         ! double precision vector passed from/to R

 INTEGER :: I

  if (ip(1) < numOutput) call rexit("nout is not large enough")

! ------------------------------------------------------------------------------
! Sometimes very small negative values are created - CHECK THIS  
  DO I = 1, neq   
    y(I) = MAX(y(I), 0.d0)
  END DO
  

  
! extract the state variables from y
  DO I = 1, N_Cohort    
    RESERVE(I)     = y (0*N_Cohort + I)
    STRUCT(I)      = y (1*N_Cohort + I)
    REPROD(I)      = y (2*N_Cohort + I)
    POPULATION(I)  = y (3*N_Cohort + I)  
    
!    POPULATION(I) = MAX(POPULATION(I), 0.d0) is this also necessary for population?
  END DO
  
  LARV_DENS = y (4*N_Cohort + 1) 
  LARV_BIOM = y (4*N_Cohort + 2) 
  PHYTO     = y (4*N_Cohort + 3) 
  DETRITUS  = y (4*N_Cohort + 4) 
  SIM       = y (4*N_Cohort + 5) 
  O2        = y (4*N_Cohort + 6) 

! ------------------------------------------------------------------------------
! calculate the derivatives


  CALL mussel_dynamics(t)      ! dSTRUCT, dRESERVE, dREPROD
  CALL population_dynamics     ! dPOPULATION
  
  CALL summedRates
  
  CALL larvae_dynamics(t)      ! dLARV_DENS, dLARV_BIOM 
  IF (dynamic_environment) THEN
    CALL environment_dynamics  ! dPHYTO, dDETRITUS, dSIM, dO2 
  END IF 
! ------------------------------------------------------------------------------
! save the derivatives in dy

  DO I = 1, N_Cohort    
    dy (0*N_Cohort + I) = dRESERVE(I)
    dy (1*N_Cohort + I) = dSTRUCT(I)  
    dy (2*N_Cohort + I) = dREPROD(I)  
    dy (3*N_Cohort + I) = dPOPULATION(I)         
  END DO
  
  dy (4*N_Cohort + 1) = dLARV_DENS        
  dy (4*N_Cohort + 2) = dLARV_BIOM        
  dy (4*N_Cohort + 3) = dPHYTO          
  dy (4*N_Cohort + 4) = dDETRITUS          
  dy (4*N_Cohort + 5) = dSIM          
  dy (4*N_Cohort + 6) = dO2          

! ------------------------------------------------------------------------------
! save the output in yout

  CALL MusselOutput(yout)
  
END SUBROUTINE Musselcohort

