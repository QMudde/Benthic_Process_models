! ==============================================================================
! ==============================================================================
!
! Cohort mussel model using DEB dynamics
!
! implementation: Quinten Mudde
!                 Karline Soetaert
! NIOZ - Yerseke
!
! Output
! ==============================================================================
! ==============================================================================

! ==============================================================================
! Setting the output
! ==============================================================================

SUBROUTINE MusselOutput(yout)
USE Mussel_deb

IMPLICIT NONE

 DOUBLE PRECISION :: yout(*)
 INTEGER :: II
  
  ii = 0
  
  ! save all outputs defined in n_cohort (keep the order! - will update ii)
  CALL saveOutputCohort(ii, yout, Gonads)
  CALL saveOutputCohort(ii, yout, Gametes)
  CALL saveOutputCohort(ii, yout, Length)
  CALL saveOutputCohort(ii, yout, C_ind)
  CALL saveOutputCohort(ii, yout, Ingestion)
  CALL saveOutputCohort(ii, yout, Catabolism)
  CALL saveOutputCohort(ii, yout, Struct_growth)
  CALL saveOutputCohort(ii, yout, Faeces_production)
  CALL saveOutputCohort(ii, yout, Regression)
  CALL saveOutputCohort(ii, yout, Assimilation)
  CALL saveOutputCohort(ii, yout, Respiration)
  CALL saveOutputCohort(ii, yout, Spawn_rate)
  CALL saveOutputCohort(ii, yout, Reprod_growth)
  CALL saveOutputCohort(ii, yout, VB_grow)
  CALL saveOutputCohort(ii, yout, WW_ind)
  CALL saveOutputCohort(ii, yout, Pseudofaeces_C) ! new
  CALL saveOutputCohort(ii, yout, Filtration_C) ! new
  CALL saveOutputCohort(ii, yout, Clearance_rate) ! new
  
  IF (includePopulation) THEN
    CALL saveOutputCohort(ii, yout, Mortality)
    CALL saveOutputCohort(ii, yout, Settle_rate)
    CALL saveOutputCohort(ii, yout, C_m2)
  END IF
  ! save all outputs that are summed over all cohorts
  CALL MusselOutput0D(ii, yout)
  
END SUBROUTINE MusselOutput

! ==============================================================================

SUBROUTINE MusselOutput0D(ii, yout)
USE Mussel_dim

IMPLICIT NONE
DOUBLE PRECISION :: yout(*), out(nout), forcs(nforcs)
INTEGER :: II
INTEGER :: I
common / musselout   / out
common / musselforcs / forcs

! 0-D variables
  IF (includePopulation) THEN
    DO I = 1, nout
     yout(ii + i) = out(i)
    END DO
    ii = ii + nout
  END IF
  
  
! forcing functions
  DO I = 1, nforcs
    yout(ii + i) = forcs(i)
  END DO

END SUBROUTINE MusselOutput0D

! ==============================================================================

SUBROUTINE saveOutputCohort(ii, yout, var)

 USE Mussel_dim, only :  N_Cohort

 IMPLICIT NONE
 DOUBLE PRECISION, INTENT (INOUT) :: yout(*)
 DOUBLE PRECISION, INTENT (IN)    :: var(N_Cohort)

 INTEGER, INTENT (INOUT) :: ii
 INTEGER :: I

  DO I = 1, N_Cohort
    yout(ii + I) = var(I)
  END DO

  ii = ii + N_Cohort

END SUBROUTINE saveOutputCohort

