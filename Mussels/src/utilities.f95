! ==============================================================================
! ==============================================================================
!
! Cohort mussel model using DEB dynamics
!
! implementation: Quinten Mudde
!                 Karline Soetaert
! NIOZ - Yerseke
!
! utility functions
! ==============================================================================
! ==============================================================================


! ==============================================================================
! Temperature dependence function
! ==============================================================================

DOUBLE PRECISION FUNCTION Arrhenius(Temp_C)

USE Mussel_deb

IMPLICIT NONE

 DOUBLE PRECISION :: Temp_C
 DOUBLE PRECISION :: Temp_K, Arr_base, Arr_bnd_nom, Arr_bnd_den

! Temperature in Kelvin 
 
    Temp_K      = Temp_C + 273.15d0  
    Arr_base    =       exp(T_a/T_r - T_a/Temp_K )
    
! behavior near boundaries    
    Arr_bnd_nom = 1.d0 + exp(T_al/T_r    - T_al/T_l) +                          &
                        exp(T_ah/T_h    - T_ah/T_r)
    
    Arr_bnd_den = 1.d0 + exp(T_al/Temp_K - T_al/T_l) +                          &
                        exp(T_ah/T_h    - T_ah/Temp_K)   
    
    Arrhenius = Arr_base * Arr_bnd_nom / Arr_bnd_den

END FUNCTION Arrhenius

! ==============================================================================
! Normal and Standardized normal distribution (max value at t=mean is 1)
! ==============================================================================
! Area under full dnorm curve = 1

DOUBLE PRECISION FUNCTION dnorm(t, mean, std)

IMPLICIT NONE

 DOUBLE PRECISION :: t, mean, std
 DOUBLE PRECISION, PARAMETER :: pi = 3.14159d0

   Dnorm = 1.d0/SQRT(2.d0*pi*std**2.d0) * exp(-(t-mean)**2.d0/2.d0/std**2.d0)

END FUNCTION dnorm

! ==============================================================================
! optimum point (mean) of stnorm = 1

DOUBLE PRECISION FUNCTION stnorm(t, mean, std)

IMPLICIT NONE

 DOUBLE PRECISION :: t, mean, std
 DOUBLE PRECISION, PARAMETER :: pi = 3.14159d0

   stnorm = exp(-(t-mean)**2.d0/2.d0/std**2.d0)

END FUNCTION stnorm

