!==========================================================================
!==========================================================================
! THE POSIDONIA model in FORTRAN
!    A model for Posidonia growth in the Bay of Calvi
! author:  
! Karline Soetaert, Quinten Mudde, Joran de Gang, 
! Fabian Lenartz, Marilaure Gregoire
! 
!==========================================================================
!==========================================================================
  
!==========================================================================
! Module with common blocks
!==========================================================================
    
  MODULE posidonia_mod

!---------------------------- 
! State variables, mmolC/m2
!---------------------------- 

       DOUBLE PRECISION  :: LEAVES, ROOTS, RHIZOMES, NRESERVES
       DOUBLE PRECISION  :: DEADLEAVES, DEADBELOWGROUND, BURIED
       DOUBLE PRECISION  :: DETRITUS, OXYGEN, DIN, DINSED, CLOSING

!---------------------------- 
! derivatives, mmolC/m2/d
!---------------------------- 
       DOUBLE PRECISION  :: dLEAVES, dROOTS, dRHIZOMES, dNRESERVES
       DOUBLE PRECISION  :: dDEADLEAVES, dDEADBELOWGROUND, dBURIED 
       DOUBLE PRECISION  :: dDETRITUS, dOXYGEN, dDIN, dDINSED, dCLOSING

!---------------------------- 
! forcing functions
!---------------------------- 
     DOUBLE PRECISION  ::    &
       fTemperature,  & ! dg C
       fLight,        & ! W/m2
       fDaylength,    & ! hours
       fWind,         & ! m/s
       fPiston,       & ! m/d
       fSatox           ! mol/m3
       
     COMMON /myforcsP/ fTemperature, fLight, fDaylength, fWind,                &
        fPiston, fSatox
       
!---------------------------- 
! parameters
!---------------------------- 
    DOUBLE PRECISION  :: salinity, latitude, depth, C_DWT, OCr, NCr,           &
       arrheniusCt, POPT, TOPT, TMAX, TMIN, maxPhotosynthesis,                 &
       ksLight, k0short, k0long, plong, ksE, maxEshading,                      &
       L_surf_Ratio, max_selfshade, ks_shade, maxRootsNupt, ksRootsNupt,       &
       maxRhizomeNC, minRhizomeNC, ksROOTS, rhizomeMobilisationRate,           &
       rhizomeRespirationRate, rootRespirationRate, leafRespirationRate,       &
       respCost, belowMortalityRate, leafMortalityRate,                        &
       partLeafReallocation, ReallRespcost, leafAbscissionRate,                &
       NecrosisRate, belowgroundDecayRate, abovegroundDecayRate,               &
       detritusDecayRate, pBurial, leafBreakingRate, deadleafBreakingRate,     &
       expWind, ksW, detritusExportRate, mindetritusExportRate, shapebelow,    &
       pBelowToRoot, pBelowGrowthMin, TimingEndMobilRhizomes,                  &
       TimingRestartMobilRhizomes, TimingMaxStorageBelowground,                &
       senescenceDaylength, kSenescence, refDaylength, NexchangeRate
       
    COMMON /myparmsP    /salinity, latitude, depth, C_DWT, OCr, NCr,           &
       arrheniusCt, POPT, TOPT, TMAX, TMIN, maxPhotosynthesis,                 &
       ksLight, k0short, k0long, plong, ksE, maxEshading,                      &
       L_surf_Ratio, max_selfshade, ks_shade, maxRootsNupt, ksRootsNupt,       &
       maxRhizomeNC, minRhizomeNC, ksROOTS, rhizomeMobilisationRate,           &
       rhizomeRespirationRate, rootRespirationRate, leafRespirationRate,       &
       respCost, belowMortalityRate, leafMortalityRate,                        &
       partLeafReallocation, ReallRespcost, leafAbscissionRate,                &
       NecrosisRate, belowgroundDecayRate, abovegroundDecayRate,               &
       detritusDecayRate, pBurial, leafBreakingRate, deadleafBreakingRate,     &
       expWind, ksW, detritusExportRate, mindetritusExportRate, shapebelow,    &
       pBelowToRoot, pBelowGrowthMin, TimingEndMobilRhizomes,                  &
       TimingRestartMobilRhizomes, TimingMaxStorageBelowground,                &
       senescenceDaylength, kSenescence, refDaylength, NexchangeRate
       
!---------------------------- 
! output variables (64)
!---------------------------- 
    DOUBLE PRECISION  ::         & !
      
       LeavesDWT                , & ! gDWT/m2
       BelowDWT                 , & !  gDWT/m2
       GrossPrimaryProduction   , & ! molC/m2/d
       AssimilateProduction     , & ! molC/m2/d
       TotalRespiration         , & ! molC/m2/d
       NetO2production          , & ! molO2/m2/d
       BasalRespiration         , & ! molO2/m2/d
       GrowthRespiration        , & ! molC/m2/d
       TotalDecay               , & ! molC/m2/d
       LeafPrimaryProduction    , & ! molC/m2/d
       LeafAssimilation         , & ! molC/m2/d
       LeafGrowth               , & ! molC/m2/d
       LeafTotalRespiration     , & ! molC/m2/d
       LeafMortality            , & ! molC/m2/d
       LeafNetReallocation      , & ! molC/m2/d
       LeafReallRespiration     , & ! molC/m2/d
       LeafBreaking             , & ! molC/m2/d
       RootAssimilation         , & ! molC/m2/d
       RootGrowth               , & ! molC/m2/d
       RootTotalRespiration     , & ! molC/m2/d
       RootMortality            , & ! molC/m2/d
       RhizomeAssimilation      , & ! molC/m2/d
       RhizomeGrowth            , & ! molC/m2/d
       RhizomeTotalRespiration  , & ! molC/m2/d
       RhizomeMortality         , & ! molC/m2/d
       RhizomeMobilisation      , & ! molC/m2/d
       DeadLeafAbscission       , & ! molC/m2/d
       DeadLeafBreaking         , & ! molC/m2/d
       DeadLeafDecay            , & ! molC/m2/d
       DetritusExport           , & ! molC/m2/d
       DetritusDecay            , & ! molC/m2/d
       BelowgroundDecay         , & ! molC/m2/d
       MortToDeadBelowGround    , & ! molC/m2/d
       Burial                   , & ! molC/m2/d
       RootsNupt                , & ! molN/m2/d
       NreserveToGrowth         , & ! molN/m2/d
       RespirationNRecycling    , & ! molN/m2/d
       LeafNRecycling           , & ! molN/m2/d
       RhizomeNRecycling        , & ! molN/m2/d
       NFluxtoBottom            , & ! molN/m2/d
       AirSeaExchange           , & ! molO2/m2/d
       Piston                   , & ! m/d
       Satox                    , & ! mol O2/m3
       WindBreaking             , & ! -
       pBelowGrowth             , & ! -
       TotalC                   , & ! molC/m2
       TotalN                   , & ! molN/m2
       Senescence               , & ! -
       Remobilization           , & ! -
       LossTempfac              , & ! -
       GrowthTempfac            , & ! -
       Self_shading             , & ! -
       Leaf_surf                , & ! cm2
       PSLightLim               , & ! -
       PSepiphyteLim            , & ! -
       PSNCquotaLim             , & ! -
       PSlim                    , & ! -
       RupNLim                  , & ! -
       RupRootLim               , & ! -
       RupNCquotaLim            , & ! -
       Temperature              , & ! degC
       Wind                     , & ! m/s
       Light                    , & ! uEinst/m2/s
       PAR                      , & ! uEinst/m2/s
       Daylength                    ! hour


    COMMON /myoutP      /LeavesDWT, BelowDWT, GrossPrimaryProduction,          &
       AssimilateProduction, TotalRespiration, NetO2production,                &
       BasalRespiration, GrowthRespiration, TotalDecay, LeafPrimaryProduction, &
       LeafAssimilation, LeafGrowth, LeafTotalRespiration, LeafMortality,      &
       LeafNetReallocation, LeafReallRespiration, LeafBreaking,                &
       RootAssimilation, RootGrowth, RootTotalRespiration, RootMortality,      &
       RhizomeAssimilation, RhizomeGrowth, RhizomeTotalRespiration,            &
       RhizomeMortality, RhizomeMobilisation, DeadLeafAbscission,              &
       DeadLeafBreaking, DeadLeafDecay, DetritusExport, DetritusDecay,         &
       BelowgroundDecay, MortToDeadBelowGround, Burial, RootsNupt,             &
       NreserveToGrowth, RespirationNRecycling, LeafNRecycling,                &
       RhizomeNRecycling, NFluxtoBottom, AirSeaExchange, Piston, Satox,        &
       WindBreaking, pBelowGrowth, TotalC, TotalN, Senescence, Remobilization, &
       LossTempfac, GrowthTempfac, Self_shading, Leaf_surf, PSLightLim,        &
       PSepiphyteLim, PSNCquotaLim, PSlim, RupNLim, RupRootLim,                &
       RupNCquotaLim, Temperature, Wind, Light, PAR, Daylength

  END MODULE posidonia_mod
      
