!==========================================================================
!==========================================================================
! THE POSIDONIA model in FORTRAN
! Created for the FACE-iT summer School: 
!    A model for Posidonia growth in the Bay of Calvi
! author:  Marilaure Gregoire, Fabian Lenartz and Karline Soetaert
! date: June 23rd, 2017
!==========================================================================
!==========================================================================
  
!==========================================================================
! initialise the common block with parameter values, 
!==========================================================================
    
       SUBROUTINE initposidonia (funparms)
       IMPLICIT NONE
       EXTERNAL steadyparms
  
       INTEGER,PARAMETER :: nc = 57
  
       DOUBLE PRECISION parms(nc)
       COMMON /myparmsP/parms
  
       CALL funparms(nc, parms)
  
       RETURN
       END SUBROUTINE initposidonia
  
!==========================================================================
! Initialise the forcing function common block  
!==========================================================================
    
       SUBROUTINE initforcp (funforcs)
       IMPLICIT NONE
       EXTERNAL steadyforcs
  
       INTEGER,PARAMETER :: N = 6 
  
       DOUBLE PRECISION forcs(N)
       COMMON /myforcsP/forcs
  
       CALL funforcs(N, forcs)
  
       RETURN
       END SUBROUTINE initforcp
  
!==========================================================================
!==========================================================================
! subroutine calculating the rate of change of
! the Posidonia model - here the quantities in the common blocks are named
!==========================================================================
!==========================================================================
    
    SUBROUTINE posidonia (neq, t, Conc, dConc, yout, ip)
       
    USE posidonia_mod       ! declaration of states, variables, parameters
    IMPLICIT NONE
  
!......................... declaration section.............................
    INTEGER           :: neq, ip(*), i
  
    DOUBLE PRECISION  :: t, Conc(12), dConc(12), yout(*)

!---------------------------- 
! ordinary variables 
!---------------------------- 

    DOUBLE PRECISION :: DayOfYear, dTime, Force !!!! check force
    DOUBLE PRECISION :: LeafGrossReallocation, Leafnecrosis,                   &
                        LeafBasalRespiration, LeafGrowthRespiration,           &
                        LeafFall, LeafLoss
    DOUBLE PRECISION :: RhizomeBasalRespiration, RhizomeGrowthRespiration,     &
                        Rhizomenecrosis, RhizomeLoss
    DOUBLE PRECISION :: RootBasalRespiration, RootGrowthRespiration,           &
                        Rootnecrosis, RootLoss
    DOUBLE PRECISION :: maxNreserve, minNreserve, minNRhi  
    DOUBLE PRECISION :: Necroticfac, PAR0, PARtop, SedConc
    DOUBLE PRECISION :: pLeafGrowth, pRhizomeGrowth, pRootGrowth
    DOUBLE PRECISION :: breakingForce, WindLeafBreaking 
    
!............................ statements ..................................
  
!     check memory allocated to output variables
      IF (ip(1) < 64)  CALL rexit("nout should be at least 65") 
  
! from Conc to state variables
      LEAVES          = Conc(1)
      ROOTS           = Conc(2)
      RHIZOMES        = Conc(3)
      NRESERVES       = Conc(4)
      DEADLEAVES      = Conc(5)
      DEADBELOWGROUND = Conc(6)
      BURIED          = Conc(7)
      DETRITUS        = Conc(8)
      OXYGEN          = Conc(9)
      DIN             = Conc(10)
      DINSED          = Conc(11)
      CLOSING         = Conc(12)

   
!=============================!
! Initialisation              !
!=============================!
    
! set forcing functions  
    Light       = flight           ! determines photosynthesis
    Wind        = fwind            ! modulates the air-sea exchange
    Temperature = ftemperature     ! affects biological rates and air-sea exchange 
    Daylength   = fdaylength
    Piston      = fPiston
    Satox       = fSatox

! resistance of Posidonia to wind and waves : if wind > force: leaves break   
    breakingForce  = 11.5d0 + 5.d0*cos(2*3.14159*(t-50.d0)/365d0)       

! The day of the year -  used for phenology
    DayOfYear = t - INT(t/365)*365d0        ! modulus

! All rates are expressed at 18 dg - 
! LossTempfac = the factor of increase with temperature
    LossTempfac = arrheniusCt ** (Temperature-18D0)

! Growth increases with temp until it reaches a maximum and then declines

    IF (Temperature >= TMIN .and. Temperature <= TMAX) THEN
      GrowthTempfac = POPT * (((Temperature-TMAX)*(Temperature-TMIN)**2.d0)/   &
                (TOPT-TMIN) * ((TOPT-TMIN) * (Temperature-TOPT) -              &
                               (TOPT-TMAX) * (TOPT+TMIN-2*Temperature)))
    ELSE                         
      GrowthTempfac = 0.d0
    END IF 
  
!=============================!
! Respiration and mortality   !
!=============================!

! Basal respiration: biomass is respired at a constant rate, temperature dependent
    LeafBasalRespiration    = leafRespirationRate    * LossTempfac * LEAVES  
    RootBasalRespiration    = rootRespirationRate    * LossTempfac * ROOTS
    RhizomeBasalRespiration = rhizomeRespirationRate * LossTempfac * RHIZOMES
    
! respiration in sediment
    BasalRespiration       = LeafBasalRespiration +                            &
                             RootBasalRespiration +                            &
                             RhizomeBasalRespiration

! During senescence (late in the year), leaves die-off much more rapidly
     Senescence = 1.d0/(1d0 + dexp(kSenescence * (Daylength - refDaylength)))

! The mortality rates - leaves die only during the senescent phase.               
    LeafMortality     = leafMortalityRate  * LossTempfac * LEAVES * Senescence
    RootMortality     = belowMortalityRate * LossTempfac * ROOTS
    RhizomeMortality  = belowMortalityRate * LossTempfac * RHIZOMES
    
! Necrotic rates - biomass dying under heatstress

! Rate at which necrosis happens ---- NOT YET USED. 
    Necroticfac = min(1d0, 1d0 / (1d0 + dexp(-(Temperature - 28.6d0))))
    
    Leafnecrosis           = NecrosisRate * Necroticfac * LEAVES
    Rootnecrosis           = NecrosisRate * Necroticfac * ROOTS
    Rhizomenecrosis        = NecrosisRate * Necroticfac * RHIZOMES

! =========================== # 
!  Roots take N from sediment #
! =========================== #
    
! calculate max N reserves Posidonia can contain effect of RHIZOMES
    minNRhi     = RHIZOMES * NCr
    maxNreserve = RHIZOMES * maxRhizomeNC 
    minNreserve = RHIZOMES * minRhizomeNC
    
! N-uptake limited by N availability sediment, root density and NC content
    RupNLim       = DINSED/(DINSED + ksRootsNupt)
    RupRootLim    = ksROOTS/(ksROOTS + ROOTS)
    RupNCquotaLim = (maxNreserve - NRESERVES)/(maxNreserve - minNreserve)
!    RupNCquotaLim = (maxNreserve - NRESERVES)/maxNreserve
    
    RootsNupt = maxRootsNupt * GrowthTempfac  *                                &
                 RupNLim                      *                                &
                 RupRootLim                   *                                &
                 RupNCquotaLim                * ROOTS
  
! ================================ # 
!  Reallocation and remobilisation #
! ================================ #
! Reallocation of aboveground assimilates occurs only during senescence
    
    LeafGrossReallocation     =     partLeafReallocation  * LeafMortality
    LeafFall                  = (1.-partLeafReallocation) * LeafMortality
    
    LeafReallRespiration   = ReallRespcost*LeafGrossReallocation
    LeafNetReallocation    = LeafGrossReallocation - LeafReallRespiration
    
! Remobilisation from rhizomes - only when Remobilization = 1
! At certain times organic matter in the rhizomes is reallocated
     IF (Daylength >= senescenceDaylength) THEN
        Remobilization = 1.d0
        RhizomeMobilisation = rhizomeMobilisationRate * LossTempfac * RHIZOMES 
     ELSE
        Remobilization = 0.d0
        RhizomeMobilisation = 0.d0
     END IF    
     
! ================================ # 
!  Photosynthesis                  #
! ================================ #
    
! Photosynthesis depends on light, 
! light decreases with water depth, not all light is PAR, some is reflected
! assume Posidonia is 0.5 m high
    PAR0 = (1.D0-0.05)*0.46*Light ! 46 % of light in par 
    
    PARtop = PAR0 * dexp (-(   plong *k0long +                                 &
                           (1.-plong)*k0short)*(depth-0.5))  

! surfaces (cm2) that capture light for photosynthesis
    Leaf_surf      = (LEAVES + DEADLEAVES) *L_surf_Ratio
    Self_shading   = max_selfshade * Leaf_surf/(Leaf_surf+ks_shade)
 
! PAR at the leaf    
    PAR            = PARtop * (1.d0 - Self_shading)  

! Light limitation    
    PSLightLim     = PAR/(PAR + ksLight)
    
! PS is reduced if the N:C ratio of the rhizomes is very low
!    PSNCquotaLim = 1 - minNRhi / (minNRhi+NRESERVES)  ## WAS LIKE THAT
!    PSNCquotaLim      = NRESERVES/(NRESERVES + minNreserve)
    PSNCquotaLim   = (NRESERVES - minNreserve)/(maxNreserve - minNreserve)
    
! light is additionally blocked by Epifte fouling of leaves
! Epifytes decrease PS with at most maxEshading ([0-1]); 
! assume epifyte biomass ~ DIN concentration
    PSEpiphyteLim  = 1.d0 - maxEshading * DIN/(DIN+ksE)

! resulting primary production  
    LeafPrimaryProduction = maxPhotosynthesis * GrowthTempfac   *              &
                             PSLightLim                         *              &
                             PSNCquotaLim                       *              &
                             PSEpiphyteLim                      * LEAVES

    GrossPrimaryProduction = LeafPrimaryProduction

! ================================ # 
!  Assimilation                    #
! ================================ #
    
! C Assimilates available for growth  
    AssimilateProduction = LeafPrimaryProduction  +                            & 
                           LeafNetReallocation    +                            &
                           RhizomeMobilisation
    
! Assimilates are allocated to belowground and aboveground tissues. 
! The way this is distributed over the various organs is not constant.  
     IF (DayOfYear > TimingEndMobilRhizomes) THEN
        dTime = (DayOfYear                  -TimingEndMobilRhizomes)/          &
                (TimingMaxStorageBelowground-TimingEndMobilRhizomes)
                
        pBelowGrowth = min (1.D0, dTime**shapebelow + pBelowGrowthMin)
     ELSE
        pBelowGrowth = pBelowGrowthMin
     ENDIF
    
! Assimilate partitioning
    pLeafGrowth          = (1. - pBelowGrowth)
    pRootGrowth          = pBelowGrowth  * pBelowToRoot
    pRhizomeGrowth       = pBelowGrowth  * (1-pBelowToRoot)
    
    LeafAssimilation     = AssimilateProduction * pLeafGrowth
    RootAssimilation     = AssimilateProduction * pRootGrowth
    RhizomeAssimilation  = AssimilateProduction * pRhizomeGrowth
    
! ================================ # 
!  Respiration relating to growth  #
! ================================ #

    LeafGrowthRespiration    = LeafAssimilation    * respCost
    RhizomeGrowthRespiration = RhizomeAssimilation * respCost
    RootGrowthRespiration    = RootAssimilation    * respCost

    GrowthRespiration        = LeafGrowthRespiration      +                    &
                               RootGrowthRespiration      +                    & 
                               RhizomeGrowthRespiration
    
! ================================ # 
!  Nutrient recycling              #
! ================================ #
! assumed that only leaf N can be recycled (active Senescence)
    LeafNRecycling         = LeafGrossReallocation * NCr 
    RhizomeNRecycling      = RhizomeMobilisation   * NCr 
    
! N is not used in maintenance respiration - it is recycled
    RespirationNRecycling  = BasalRespiration     * NCr
  

! Total respiration (molC m-2d-1)
    LeafTotalRespiration    = LeafBasalRespiration      +                      & 
                              LeafGrowthRespiration     +                      &
                              LeafReallRespiration
                               
    RootTotalRespiration    = RootBasalRespiration    + RootGrowthRespiration 
    RhizomeTotalRespiration = RhizomeBasalRespiration + RhizomeGrowthRespiration

    TotalRespiration        = (LeafTotalRespiration    +                       & 
                               RootTotalRespiration    +                       &
                               RhizomeTotalRespiration )
    
! ================================ # 
! Growth of each organ             #
! ================================ # 
    LeafGrowth    = LeafAssimilation    - LeafGrowthRespiration
    RootGrowth    = RootAssimilation    - RootGrowthRespiration
    RhizomeGrowth = RhizomeAssimilation - RhizomeGrowthRespiration

! ================================ # 
! The decay of dead biomass        #
! ================================ # 
    
    BelowgroundDecay    = belowgroundDecayRate * LossTempfac * DEADBELOWGROUND
    DeadLeafDecay       = abovegroundDecayRate * LossTempfac * DEADLEAVES     
    DetritusDecay       = detritusDecayRate    * LossTempfac * DETRITUS       

    DeadLeafAbscission  = leafAbscissionRate   * LossTempfac * DEADLEAVES     
    
    TotalDecay          = BelowgroundDecay + DeadLeafDecay + DetritusDecay
    
 ! (simple) representation of Carbon sequestration
    Burial               =        pBurial  * (RhizomeMortality + RootMortality)
    MortToDeadBelowGround = (1.d0-pBurial) * (RhizomeMortality + RootMortality)

! ================================ # 
! Wind breaking of leaves          #
! and deadleaves                   #      
! ================================ # 
! Storms break leaves creating deadleaves and export detritus
    WindBreaking = Wind**expWind / (Wind**expWind + ksW**expWind)

! Windbreaking is 1 during a storm
    IF (Wind > breakingForce) THEN
       WindLeafBreaking = 1.d0  !((Wind - Force) > 1)
    ELSE
       WindLeafBreaking = 0.d0  !((Wind - Force) < 1)
    END IF

    LeafBreaking      = leafBreakingRate      * LEAVES     * WindLeafBreaking
    DeadLeafBreaking  = deadleafBreakingRate  * DEADLEAVES * WindBreaking
    DetritusExport    = detritusExportRate    * DETRITUS   * WindBreaking +    &
                        mindetritusExportRate * DETRITUS

! ================================ # 
! Atmospheric air-sea exchange     #
! ================================ # 
    AirSeaExchange = Piston*(Satox - OXYGEN)   ! mol O2 /m2/d
    
! ================================ # 
! Nutrient exchanges               #
! ================================ # 
    
! molN/m2/day released

! assumes thickness of sediment = 1 m
    SedConc = DINSED/1.d0
    NFluxtoBottom = (DIN - SedConc) * NexchangeRate 

!=======================================!
!=======================================!
! The rate of change of state variables !
!=======================================!
!=======================================!
    
    RootLoss       =  RootBasalRespiration    +                                &
                      RootMortality
    
    dROOTS         =  RootGrowth - RootLoss
    
    RhizomeLoss    =  RhizomeBasalRespiration +                                & 
                      RhizomeMortality        +                                &
                      RhizomeMobilisation 
    
    dRHIZOMES      =  RhizomeGrowth - RhizomeLoss
    
    LeafLoss       =  LeafBasalRespiration   +                                 &
                      LeafFall               +                                 &
                      LeafGrossReallocation  +                                 &
                      LeafBreaking  
    
    dLEAVES        =  LeafGrowth - LeafLoss
    
    ! (Leaves N uptake not considered)
    
    NreserveToGrowth = (LeafGrowth + RootGrowth + RhizomeGrowth)*NCr        

    dNRESERVES     =    RootsNupt                                              &
                      - NreserveToGrowth                                       &
                      +  RespirationNRecycling                                 &
                      +  LeafNRecycling                                        &
                      +  RhizomeNRecycling    
    

    ! LeafMortality to LeafFall
    dDEADLEAVES      =     LeafFall                                            &
                         - DeadLeafDecay                                       &
                         - DeadLeafAbscission                                  &
                         - DeadLeafBreaking
    
    dDEADBELOWGROUND = MortToDeadBelowGround                                   &
                      - BelowgroundDecay 
                      
    dBURIED          = Burial      


        
    dDETRITUS        =     DeadLeafAbscission                                  &
                          + LeafBreaking                                       &
                          + DeadLeafBreaking                                   &
                          - DetritusDecay                                      &
                          - DetritusExport 
    
    NetO2production  = GrossPrimaryProduction - TotalRespiration
                
    dOXYGEN          =   ( NetO2production                                     &
                           - TotalDecay    )  * OCr/depth                      &
                           + AirSeaExchange        /depth 
    
! Added detritusexport  
    dDIN             =    DetritusExport/depth * NCr                           &
                        - NFluxtoBottom /depth  
    
    ! note: NFluxBottom is dividided by 1 (thickness of sediment layer)
    dDINSED    =     NFluxtoBottom                                             & 
                    - RootsNupt                                                &
                    + Burial     * NCr                                         &
                    + TotalDecay * NCr  
    
! CLOSING check: it contains the losses of carbon to the outside system
    dCLOSING  =     TotalDecay                                                 &
                   + DetritusExport                                            &
                   + TotalRespiration                                          &
                   - GrossPrimaryProduction 
    
! Compute total carbon and nitrogen -  this should be constant in the model
    TotalC  =  LEAVES + ROOTS + RHIZOMES + DEADLEAVES +                        &
                DEADBELOWGROUND + DETRITUS + CLOSING + BURIED
    
    TotalN  = (LEAVES + ROOTS + RHIZOMES + DEADLEAVES +                        & 
                DEADBELOWGROUND +  DETRITUS  ) * NCr +                         &
                DIN * depth  + DINSED + NRESERVES

! The rates of change returned
      dConc(1) = dLEAVES  
      dConc(2) = dROOTS   
      dConc(3) = dRHIZOMES  
      dConc(4) = dNRESERVES  
      dConc(5) = dDEADLEAVES  
      dConc(6) = dDEADBELOWGROUND  
      dConc(7) = dBURIED
      dConc(8) = dDETRITUS  
      dConc(9) = dOXYGEN  
      dConc(10) = dDIN 
      dConc(11) = dDINSED 
      dConc(12) = dCLOSING 

      CALL getoutP(yout)
  
  RETURN
  END SUBROUTINE
  
!==========================================================================
! put output variables in one vector
!==========================================================================
    
    SUBROUTINE getoutP(yout)
    INTEGER :: i
    INTEGER, PARAMETER :: nout = 65
    DOUBLE PRECISION :: yout(*), out(nout)
  
    COMMON /myoutP   /out
      DO i = 1, nout
       yout(i) = out (i)
      ENDDO       
  
    END SUBROUTINE getoutP
  
