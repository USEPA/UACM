MODULE module_sf_pxlsm
 
  USE module_model_constants
  USE module_sf_pxlsm_data

  INTEGER, PARAMETER   :: NSOLD   = 20
  REAL, PARAMETER      :: RD      = 287.04,   CPD   = 1004.67,             &
                          CPH2O   = 4.218E+3, CPICE = 2.106E+3,            &
                          LSUBF   = 3.335E+5, SIGMA = 5.67E-8,             &  
                          ROVCP   = RD / CPD
                          
  DOUBLE PRECISION, PARAMETER      :: CRANKP  = 0.5                    ! CRANK-NIC PARAMETER
  REAL, PARAMETER      :: RIC     = 0.25                   ! critical Richardson number
  REAL, PARAMETER      :: DENW    = 1000.0                 ! water density in KG/M3                  
  DOUBLE PRECISION, PARAMETER      :: TAUINV  = 1.0 / 86400.0          ! 1/1DAY(SEC)
  DOUBLE PRECISION, PARAMETER      :: T2TFAC  = 1.0 / 10.0             ! Bottom soil temp response factor
  DOUBLE PRECISION, PARAMETER      :: T2TFACU = 1.0 / 4.0              ! Bottom layer temp response factor (for UTK Urban)
  REAL, PARAMETER      :: PI      = 3.1415926535897932384626433 ! same as piconst in module_model_constants.F
  REAL, PARAMETER      :: RADDEG  = 180.0/PI               ! degrees per radian
  REAL, PARAMETER      :: PR0     = 0.95
  REAL, PARAMETER      :: CZO     = 0.032
  REAL, PARAMETER      :: OZO     = 1.E-4

CONTAINS
!

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
   SUBROUTINE pxlsm(U3D, V3D, PP3D, DZ8W, QV3D, T3D, TH3D, RHO,          &
                    PSFC, GSW, GLW, RAINBL, EMISS, JULIAN, XTIME, GMT,   &
                    ITIMESTEP, CURR_SECS, NSOIL, DT, ANAL_INTERVAL,      &
                    XLAND, XICE, ALBBCK, ALBEDO,                         &
                    SNOALB, SMOIS, TSLB, MAVAIL, TA2,                    &
                    QA2, QSFC, ZS, DZS, PSIH,                            &
                    LANDUSEF, SOILCTOP, SOILCBOT, VEGFRA, VEGF_PX,       &
                    ISLTYP, RA, RS, LAI, IMPERV, CANFRA, NLCAT, NSCAT,   &
                    HFX, HFX_RF, HFX_WL, ANTHFX1, ANTHFX_RF, QFX_RF,     &   !-- UTK Urban
                    LH_RF, TG_RF, T2_RF, TG_WL, T2_WL, T2_GND,           &   !-- UTK Urban
                    RADWL, RADWL_FNL, RADST, RADST_FNL, GSW_URB, GSW_WL, &   !-- UTK Urban
                    RADNET, RADNET_URB, RADNET_RF, RADNET_WL,            &   !-- UTK Urban
                    LWDN_URB, LWDN_RF, LWDN_WL,                          &   !-- UTK Urban 
                    LWUP_URB, LWUP_RF, LWUP_WL,                          &   !-- UTK Urban
                    QFX, LH, TSK, SST, ZNT, CANWAT,                      &
                    GRDFLX, SHDMIN, SHDMAX,                              &
                    SNOWC, PBLH, RMOL, UST, CAPG, DTBL,                  &
                    T2_NDG_OLD, T2_NDG_NEW,                              &
                    Q2_NDG_OLD, Q2_NDG_NEW,                              &
                    SN_NDG_OLD, SN_NDG_NEW, SNOW, SNOWH, SNOWNCV,        &
                    T2OBS, Q2OBS, PXLSM_SMOIS_INIT,                      &
                    PXLSM_SOIL_NUDGE,                                    &
                    cosz, XLAT, XLONG, SINA, COSA,                       &   !-- DDY Urban
                    LAMP, HGTBDG, lamf,                                  &
                    STDIR1,STDIR2,BLKW1,BLKW2,                           &
                    SWDDIR, SWDDIF,                                      &   !-- UTK Urban
                    pxlsm_modis_veg, ifuacm,                             &
                    LAI_PX,                                              &
                    WWLT_PX, WFC_PX, WSAT_PX,                            &
                    CLAY_PX, CSAND_PX, FMSAND_PX,                        &
                    GAB_RF, GAB_WL, WSP_CAN,                             &
                    ids,ide, jds,jde, kds,kde,                           &
                    ims,ime, jms,jme, kms,kme,                           &
                    its,ite, jts,jte, kts,kte                            )

!-------------------------------------------------------------------------
!   THIS MODULE CONTAINS THE PLEIM-XIU LAND-SURFACE MODEL (PX-LSM).
!   IT IS DESIGNED TO SIMULATE CHARACTERISTICS OF THE LAND SURFACE AND
!   VEGETATION AND EXCHANGE WITH THE PLANETARY BOUNDARY LAYER (PBL). THE
!   SOIL MOISTURE MODEL IS BASED ON THE ISBA SCHEME DEVELOPED BY NOILHAN
!   AND PLANTON (1989) AND JACQUEMIN AND NOILHAN (1990) AND INCLUDES
!   PROGNOSTIC EQUATIONS FOR SOIL MOISTURE AND SOIL TEMPERATURE IN TWO
!   LAYERS (1 CM AND 1 M) AS WELL AS CANOPY WATER CONTENT.  SURFACE
!   MOISTURE FLUXES ARE MODELED BY 3 PATHWAYS: SOIL EVAPORATION, CANOPY
!   EVAPORATION, AND VEGETATIVE EVAPOTRANSPIRATION.
!   EVAPOTRANSPIRATION DIRECTLY FROM THE ROOT ZONE SOIL LAYER IS MODELED
!   VIA A CANOPY RESISTANCE ANALOG ALGORITHM WHERE STOMATAL CONDUCTANCE
!   IS CONTROLLED BY SOLAR RADIATION, AIR TEMPERATURE, AIR HUMIDITY, AND
!   ROOT ZONE SOIL MOISTURE. REQUIRED VEGETATION CHARACTERISTICS DERIVED
!   FROM THE USGS LANDUSE DATA INCLUDE: LEAF AREA INDEX, FRACTIONAL VEGETATION
!   COVERAGE, ROUGHNESS LENGTH, AND MINIMUM STOMATAL RESISTANCE. AN INDIRECT
!   NUDGING SCHEME ADJUSTS SOIL MOISTURE ACCORDING TO DIFFERENCES BETWEEN
!   MODELED TEMPERATURE AND HUMIDITY AND ANALYSED SURFACE FIELDS.
!
! References:
! Pleim and Xiu, 1995: Development and testing of a surface flux and planetary 
!                      boundary layer model for application in mesoscale models.
!                      J. Appl. Meteoro., Vol. 34, 16-32.
! Xiu and Pleim, 2001: Development of a land surface model. Part I: Application
!                      in a mesoscale meteorological model. J. Appl. Meteoro.,
!                      Vol. 40, 192-209.
! Pleim and Xiu, 2003: Development of a land surface model. Part II: Data
!                      assimilation. J. Appl. Meteoro., Vol. 42, 1811-1822.
! 
! Pleim and Gilliam, 2009: An Indirect Data Assimilation Scheme for Deep Soil Temperature in the 
!                          Pleim-Xiu Land Surface Model. J. Appl. Meteor. Climatol., 48, 1362-1376.
!
! Gilliam and Pleim, 2010: Performance assessment of new land-surface and planetary boundary layer 
!                          physics in the WRF-ARW. Journal of Applied Meteorology and Climatology, 49, 760-774. 
! REVISION HISTORY:
!     AX     4/2005  - developed the initial WRF version based on the MM5 PX LSM
!     RG     2/2008  - Completed testing of the intial working version of PX LSM, released in WRFV3.0 early 2008
!     DW     8/2011  - Landuse specific versions of PX (USGS/NLCD/MODIS) were unified into a single code with
!                      landuse characteristics defined in module_sf_pxsfclay.F.
!     RG     12/2011 - Basic code clean, removed commented out debug statements, lined up columns, etc.
!     RG     01/2012 - Removed FIRSTIME Logic that computes PX Landuse characteristics at first time step only. This resulted 
!                      in different solutions when OpenMP was used and would not work with moving domains.
!     RG     08/2012 - Added CURR_SECS variable in argument list as replacement for PX internal CURRTIME internally comp var.
!                      This is neccessary for PX to correctly interpolate analyses for soil nudging. In this same calculation
!                      logic was added for cases where user does not specify the analysis interval, or no analysis interval is
!                      relevant as in the no PX soil nudging via namelist (pxlsm_soil_nudge = 0). Prior to this fix the default
!                      analysis interval was zero, so if not speficied a divide by zero was the result. Also, changes were made to 
!                      ensure PX LSM will work with not only MODIS and USGS, but also both the 40 and 50 class NLCD-MODIS data.
!                      Also, coupled module_sf_pxlsm_data.F was updated so landuse characteristics across datasets are more
!                      consistent. Albedo for NLCD two grassland categories were lowered from 23 and 25 to 18 and 19.
!                      For the NLCD40 and NLCD50 roughness and leaf area were made consistent between the US NLCD and
!                      outside US MODIS datasets. Prior, US boundaries created boundaries of roughness and LAI.
!
!     RG     10/2014 - Wetlands soil moisture treatment. Grid cell soil moisture cannot fall less than fraction of a grid
!                      cells wetland area * soil saturation (e.g., SMOIS of cell with 50% wetlands cannot fall below 50% of WSAT)
!                    - Both soil levels are initialized using MVAVAIL (Soil moisture availability) instead of just layer 2.
!                    - Veg Cv (heat capacity) changed from 8x10-6 to 1.2x10-5 (K-M2/J)
!                    - Alternate empirical stomatal function of PAR (F1) to better replicate photosynthesis-conductance models.
!                      The main effect is to reduce stomatal conductance for low PAR.
!                    - Snow albedo is now computed using fractional land-use weighting. Values for each land-use class
!                      are defined like other PX landuse parameters in module_sf_pxlsm_data.F. These are based on values
!                      used by NOAH LSM MODIS in VEGPARM.TBL (MAXALB), but tuned to better match satellite values in maxsnowalb 
!                      dataset. Tuning reduced the MAXALB for all forest classes from values in the 50-60% range to 30-40% range.
!                      These static values are more representative of albedo after snow has melted of fallen from trees. These
!                      values were also cross verified with http://www.globalbedo.org/global.php
!                    - USGS 28 category added as an option
!                    - Impervious surface and canopy fraction data can be used if processed (otherwise 0% so no impact)
!                      to alter surface heat capacity (See SURFPX subroutine for details) in urban areas and refine
!                      LAI and VEGF_PX estimations (see VEGLAND subroutine).
!     JP     12/2015 - Surface water vapor mixing ratio calculation added for land surface, which is passed to PX-SFCLAY 
!                      for use over all non-water and non-frozen surfaces.
!                    - PAR function and impact on transpiration modified according to Echer et al.(2015). See P-X LSM documentation
!                      for full reference. These act to reduce moisture bias near surface during PBL transition.
!     JP     11/2017 - Updated vegetation table for different land cover types. Added in WRFv4.0.
!     LR     11/2017 - Update for MODIS vegetation: many changes in soil properties. Added in WRFv4.1. 
!                      (Ran et al., 2016 JGR-atmosphere, Ran et al. 2017 in preparation)
!     JP     12/2018 - revised soil type categories (ISTI) to conform to WRF soil type input data 
!                      soil types Sand through Clay are now 1-12 rather than 1-11
!     JP     12/2021 - Added new pathway for evaporation from the ground in the vegetated fraction of the grid cell.  
!     UTK    05/2018 - Estimation of urban ground+rooftop temperature and heat fluxes. 
!                    - Estimation of evaporation rate from the urban surfaces and corresponding moisture fluxes during the precipitation events.
!                    - Estimation of anthropogenic heat fluxes from the building based upon heating degree days and cooling degree days.
!     UTK    03/2021 - Incorporation of the complete urban solar radiation model based upon Masson (2000).  
!                    - Estimation of the urban wall surface temperatures and heat fluxes.
!
!--------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------
!   ARGUMENT LIST:
!
!... Inputs:
!-- U3D          3D u-velocity interpolated to theta points (m/s)
!-- V3D          3D v-velocity interpolated to theta points (m/s)
!-- PP3D         Pressure at half levels (Pa)
!-- DZ8W         dz between full levels (m)
!-- QV3D	 3D mixing ratio
!-- T3D          Temperature (K)
!-- TH3D         Theta (K)
!-- RHO          3D dry air density (kg/m^3)

!-- PSFC         surface pressure (Pa)
!-- GSW          downward short wave flux at ground surface (W/m^2)      
!-- GLW          downward long wave flux at ground surface (W/m^2)
!-- RAINBL       Timestep rainfall 
!-- EMISS        surface emissivity (between 0 and 1)

!-- ITIMESTEP    time step number 
!-- NSOIL        number of soil layers
!-- DT           time step (second)
!-- CURR_SECS    time on model domain in seconds, universal WRF variable
!-- ANAL_INTERVAL Interval of analyses used for soil moisture and temperature nudging

!-- XLAND        land mask (1 for land, 2 for water)
!-- XICE         Sea ice
!-- ALBBCK       Background Albedo
!-- ALBEDO       surface albedo with snow cover effects
!-- SNOALB       Albedo of snow

!-- SMOIS        total soil moisture content (volumetric fraction)
!-- TSLB         soil temp (k)
!-- MAVAIL       Moisture availibility of soil
!-- TA2          2-m temperature
!-- QA2          2-m mixing ratio

!-- SVPT0        constant for saturation vapor pressure (K)
!-- SVP1         constant for saturation vapor pressure (kPa)
!-- SVP2         constant for saturation vapor pressure (dimensionless)
!-- SVP3         constant for saturation vapor pressure (K)

!-- ZS           depths of centers of soil layers
!-- DZS          thicknesses of soil layers
!-- PSIH         similarity stability function for heat

!-- LANDUSEF     Landuse fraction		
!-- SOILCTOP     Top soil fraction		
!-- SOILCBOT     Bottom soil fraction		
!-- VEGFRA       Vegetation fraction                   
!-- VEGF_PX      Veg fraction recomputed and used by PX LSM
!-- ISLTYP       Soil type

!-- RA           Aerodynamic resistence			
!-- RS           Stomatal resistence			
!-- LAI          read in Leaf area index (weighted according to fractional landuse)
!-- ZNT          roughness length
!-- QSFC         Sat. water vapor mixing ratio at the surface interface

!-- IMPERV       Fraction (percent) of grid cell that is impervious surface (concrete/road/non-veg)
!-- CANFRA       Fraction (percent) of grid cell that is covered with tree canopy

!-- NLCAT        Number of landuse categories		
!-- NSCAT        Number of soil categories		

!-- HFX          net upward heat flux at the surface (W/m^2)
!-- QFX          net upward moisture flux at the surface (kg/m^2/s)
!-- LH           net upward latent heat flux at surface (W/m^2)
!-- TSK          surface skin temperature (K)
!-- SST          sea surface temperature
!-- CANWAT       Canopy water (mm)

!-- GRDFLX       Ground heat flux
!-- SFCEVP       Evaportation from surface
!-- SHDMIN       Minimum annual vegetation fraction for each grid cell
!-- SHDMAX       Maximum annual vegetation fraction for each grid cell

!-- SNOWC        flag indicating snow coverage (1 for snow cover)
!-- PBLH         PBL height (m)
!-- RMOL         1/L Reciprocal of Monin-Obukhov length
!-- UST          u* in similarity theory (m/s)
!-- CAPG         heat capacity for soil (J/K/m^3)
!-- DTBL	       time step of boundary layer calls

!-- T2_NDG_OLD   Analysis temperature prior to current time
!-- T2_NDG_NEW   Analysis temperature ahead of current time
!-- Q2_NDG_OLD   Analysis mixing ratio prior to current time
!-- Q2_NDG_NEW   Analysis mixing ratio ahead of current time

!-- SN_NDG_OLD   Analysis snow water prior to current time
!-- SN_NDG_NEW   Analysis snow water ahead of current time
!-- SNOW         Snow water equivalent
!-- SNOWH        Physical snow depth
!-- SNOWNCV      Time step accumulated snow

!-- T2OBS             Analysis temperature interpolated from prior and next in time analysese
!-- Q2OBS             Analysis moisture interpolated from prior and next in time analysese
!-- PXLSM_SMOIS_INIT  Flag to intialize deep soil moisture to a value derived from moisture availiability.
!-- PXLSM_SOIL_NUDGE  Flag to use soil moisture and temperature nudging in the PX LSM
!                     This is typically done for the first simulation.
!-- pxlsm_modis_veg   Use MODIS vegetation option: 1 yes, 0 no
!-- ifuacm            Use Urban PX if 1, 0 no
!-- LAI_PX            LAI used for PX (m^2/m^2)
!-- WWLT_PX           Computed soil wilting point for PX (m^3/m^3)
!-- WFC_PX            Computed soil field capacity for PX (m^3/m^3)
!-- WSAT_PX           Computed soil saturation for PX (m^3/m^3) 
!-- CLAY_PX           Aggregated soil clay fraction for PX (%)
!-- CSAND_PX          Aggregated soil coarse sand fraction for PX (%)
!-- FMSAND_PX         Aggregated soil fine-medium sand fraction for PX (%)
!
!-- ids             start index for i in domain
!-- ide             end index for i in domain
!-- jds             start index for j in domain
!-- jde             end index for j in domain
!-- kds             start index for k in domain
!-- kde             end index for k in domain
!-- ims             start index for i in memory
!-- ime             end index for i in memory
!-- jms             start index for j in memory
!-- jme             end index for j in memory
!-- kms             start index for k in memory
!-- kme             end index for k in memory
!-- jts             start index for j in tile
!-- jte             end index for j in tile
!-- kts             start index for k in tile
!-- kte             end index for k in tile
!... Outputs:

     IMPLICIT NONE

!.......Arguments
! DECLARATIONS - INTEGER
   INTEGER, INTENT(IN)        :: ids,ide, jds,jde, kds,kde, &
                                 ims,ime, jms,jme, kms,kme, &
                                 its,ite, jts,jte, kts,kte 

   INTEGER, INTENT(IN)        :: NSOIL, ITIMESTEP, NLCAT, NSCAT,             &
                                 ANAL_INTERVAL, PXLSM_SMOIS_INIT, PXLSM_SOIL_NUDGE

   INTEGER, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)       :: ISLTYP

! DECLARATIONS - REAL
   REAL, INTENT(IN), OPTIONAL :: CURR_SECS

   REAL, INTENT(IN)           :: DT, DTBL 

   REAL, INTENT(IN), OPTIONAL :: JULIAN, XTIME, GMT

   REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)   :: XLAT, XLONG, SINA, COSA

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN)    :: U3D, V3D, PP3D, RHO, &
                                                                  T3D, TH3D, DZ8W, QV3D           

   REAL, DIMENSION(1:NSOIL), INTENT(IN)                        :: ZS, DZS
   REAL, DIMENSION( ims:ime, 1:NSOIL, jms:jme ), INTENT(INOUT) :: SMOIS, TSLB     

   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)          :: RA, RS, LAI, ZNT, QSFC
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)            :: GRDFLX, TSK, TA2, QA2

   REAL, DIMENSION( ims:ime, 1:NLCAT, jms:jme ), INTENT(IN)    :: LANDUSEF       
   REAL, DIMENSION( ims:ime, 1:NSCAT, jms:jme ), INTENT(IN)    :: SOILCTOP, SOILCBOT

   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)    :: PSFC, GSW, GLW, RAINBL, &                
                                                         ALBBCK, SHDMIN, SHDMAX, &
                                                         PBLH, RMOL, SNOWNCV,    &
                                                         UST, MAVAIL, SST, EMISS
                                                      
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)    :: T2_NDG_OLD, T2_NDG_NEW, &                
                                                         Q2_NDG_OLD, Q2_NDG_NEW, & 
                                                         SN_NDG_OLD, SN_NDG_NEW   

   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)   :: T2OBS, Q2OBS                
                                                     
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: CAPG, CANWAT, QFX, HFX, LH, HFX_RF, HFX_WL, &
                                                         PSIH,VEGFRA, VEGF_PX, SNOW, SNOALB,         &
                                                         SNOWH, SNOWC, ALBEDO, XLAND, XICE,          &
                                                         IMPERV, CANFRA, QFX_RF, LH_RF, ANTHFX1, ANTHFX_RF

   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)    :: SWDDIR, SWDDIF

 !-- DDY and UTK Urban
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)    :: COSZ
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: HGTBDG, LAMP, LAMF
   INTEGER, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: STDIR1, STDIR2
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)    :: BLKW1, BLKW2
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: TG_RF, RADST, RADST_FNL, GSW_URB
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: TG_WL, RADWL, RADWL_FNL, GSW_WL
   DOUBLE PRECISION,  DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) ::T2_RF,T2_WL,T2_GND
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: RADNET, RADNET_URB, RADNET_RF, RADNET_WL,   &
                                                         LWDN_URB, LWDN_RF, LWDN_WL,                 &
                                                         LWUP_URB, LWUP_RF, LWUP_WL

   INTEGER, OPTIONAL, INTENT(IN) :: pxlsm_modis_veg, ifuacm

   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(OUT) :: LAI_PX, WWLT_PX, WFC_PX, WSAT_PX,   &
                                                                 CLAY_PX, CSAND_PX, FMSAND_PX
   INTEGER  ::  KWAT

   LOGICAL :: radiation

   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)   :: GAB_RF, GAB_WL, WSP_CAN
!-------------------------------------------------------------------------
!     ---------- Local Variables --------------------------------


      !----  PARAMETERS
      INTEGER, PARAMETER :: NSTPS  = 11       ! max. soil types
      REAL, PARAMETER    :: DTPBLX = 40.0     ! Max PX timestep = 40 sec
      REAL, PARAMETER    :: EP1 = 0.6083624   ! Constant for Virtual Temperature
      
      !----  INTEGERS
      INTEGER, DIMENSION( 1: NSTPS ) :: JP
      INTEGER :: J, I, NS, NUDGE, ISTI, WEIGHT, K, CNP
      INTEGER :: NTSPS, IT

      !----  REALS
      REAL, DIMENSION( ims:ime, jms:jme ) :: XLAI, XLAIMN, RSTMIN, &
                                             XVEG, XVEGMN, XSNUP,  &
                                             XALB, XSNOALB, WETFRA
                                                  
      REAL, DIMENSION( ims:ime, jms:jme ) :: EG, ER, ETR, QST

      REAL :: SFCPRS, TA1, DENS1, DENSRF, QV1, ZLVL, SOLDN, LWDN,      &
              EMISSI, PRECIP, THETA1, VAPPRS, QSBT,                    &
              WG, W2, WR, TG, T2, USTAR, MOLX, Z0,                     &
              RAIR, CPAIR, CPAIR_RF, IFLAND, ISNOW,                    &
              ES, QSS, BETAP,                                          &
              RH2_OLD, RH2_NEW, T2_OLD, T2_NEW,                        &
              CORE, CORB, TIME_BETWEEN_ANALYSIS,                       &
              G1000, ALN10, RH2OBS, HU, SNOBS,                         &
              FWSAT, FWFC, FWWLT, FB, FCGSAT, FJP, FAS,                &
              FWRES, FC3, FCLAY, FCSAND, FFMSAND,                      &  
              FSEAS, T2I, HC_SNOW, SNOW_FRA, SNOWALB,                  &
              QST12, ZFUNC, ZF1, ZA2, QV2, DT_FDDA,                    &
              FC2R, FC1SAT, DTPBL, RAW, QSS_RF, QV_RF, QV_RF1, TA_RF,  &  !-- UTK
              THETA_RF, THETAV_RF, THETA_RF1, THETAV_RF1, DZFRF, TA_CAN, &
              THETA_RFs, THETAV_RFs, THETA_WLs, DENS_CAN, THETA_CAN,   &
              US_CAN, VS_CAN, QV_CAN, PP_CAN, CPAIR_CAN

      CHARACTER (LEN = 6) :: LAND_USE_TYPE

      !-- DDY and UTK Urban
      REAL, DIMENSION( its:ite, jts:jte )    :: URBF2D
      REAL, DIMENSION( ims:ime, jms:jme )    :: PPLN_FNL, VEGU_FNL

      REAL, DIMENSION( 0:kte )               :: ZF
      REAL, DIMENSION( kts:kte )             :: ZA, DZF
      INTEGER, DIMENSION( its:ite, jts:jte ) :: KSL

      REAL :: WR_ST, WR_VEG, WR_RF, US_RF1, US_RF, VS_RF1, VS_RF, PP_RF 
      REAL :: ZF_CNP1, ZF_CNP0, ZA_CNP2, ZA_CNP1, BB1, BB2
      REAL :: LCZtot,sumhb,sumlp,sumbk1,sumbk2
      INTEGER lcat,lcatm40 
      REAL, DIMENSION( its:ite, jts:jte ) :: ER_ST, ER_RF 
      REAL, SAVE :: STWAT(1:350,1:400), RFWAT(1:350,1:400)        ! *IMPORTANT NOTE* :- Make sure the domain of interest fits
                                                                  !                     well within the dimensions defined here 
      LOGICAL, SAVE :: FIRSTIME(1:350,1:400) = .TRUE.  
      LOGICAL, SAVE :: FIRSTIME_TG(1:350,1:400) = .TRUE.

      LOGICAL :: ISPXURB,LCZRUN

      REAL, DIMENSION(11) :: HB_LCZ, LAMP_LCZ, BLKW1_LCZ,BLKW2_LCZ, WSTT_LCZ

      DATA HB_LCZ /60., 30., 10., 30., 15., 10., 3., 10., 8.5, 10.,10./
      DATA LAMP_LCZ /0.40,0.33,0.30,0.17,0.20,0.10,0.10,0.20,0.08,0.25,0.25/
      DATA BLKW1_LCZ/70.,70.,80.,100.,80.,60.,60.,60.,60.,100.,100./
      DATA BLKW2_LCZ/180.,180.,200.,200.,200.,200.,100.,100.,100.,100.,100./
      DATA WSTT_LCZ/25.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20./
!-------------------------------------------------------------------------
!-------------------------------Executable starts here--------------------
!

      ALN10  = ALOG(10.0)
      G1000  = g * 1.0E-3          ! G/1000
      WEIGHT = 0
      ! Determine Landuse Dataset by the number of categories
      IF (NLCAT == 50) THEN
         LAND_USE_TYPE = 'NLCD50'
      ELSE IF (NLCAT == 40) THEN
         LAND_USE_TYPE = 'NLCD40'
      ELSE IF (NLCAT == 20) THEN
         LAND_USE_TYPE = 'MODIS'
      ELSE IF (NLCAT == 21) THEN
         LAND_USE_TYPE = 'MODIS'
      ELSE IF (NLCAT == 61) THEN
         LAND_USE_TYPE = 'MODIS'
      ELSE IF (NLCAT == 24) THEN
         LAND_USE_TYPE = 'USGS'
      ELSE IF (NLCAT == 28) THEN
         LAND_USE_TYPE = 'USGS28'
      ELSE IF (NLCAT == 51) THEN
         LAND_USE_TYPE = 'NLCD51'
      ELSE
         call wrf_error_fatal("Error: Unknown Land Use Category")
      END IF 
      
      IF (ITIMESTEP .EQ. 1) THEN
       CALL wrf_message( 'PX LSM will use the ' // TRIM(LAND_USE_TYPE) // ' landuse tables' )
       PRINT *, 'The analysis interval for surface soil and temp nudging = ',ANAL_INTERVAL,'sec.'
      ENDIF

      !-----------------------------------------------------------------------------------
      ! Kill WRF if user specifies soil nudging but provides no analysis interval, then provide helpful message.
      IF (ANAL_INTERVAL .LE. 0.0 .AND. PXLSM_SOIL_NUDGE .EQ. 1) THEN
       CALL wrf_message('PX LSM Error: The User specified analysis interval is zero or negative.')
       CALL wrf_message('If the PX LSM is used with soil nudging (pxlsm_soil_nudge=1) a wrfsfdda_d0* file is required.')
       CALL wrf_message('Make sure these files are present and')
       CALL wrf_message('Check the namelist to ensure sgfdda_interval_m is set to proper sfc analysis interval')
       STOP
      ENDIF                                                                    
      !-----------------------------------------------------------------------------------
      !--- Compute time relatve to old and new analysis time for timestep interpolation
      IF(PXLSM_SOIL_NUDGE .EQ. 1) THEN
        DT_FDDA = ANAL_INTERVAL * 1.0    ! Convert DT of Analysis to real
        TIME_BETWEEN_ANALYSIS = MOD(CURR_SECS,DT_FDDA)
        IF (TIME_BETWEEN_ANALYSIS .EQ. 0.0) THEN
          CORB = 1.0
          CORE = 0.0    
        ELSE
          CORE = TIME_BETWEEN_ANALYSIS  / DT_FDDA
          CORB = 1.0 - CORE
        ENDIF 
      ENDIF 

      !-------------------------------------------------------------
      !-- LSM sub-time loop too prevent dt > 40 sec
          NTSPS = INT(DT / (DTPBLX + 0.000001) + 1.0)
          DTPBL = DT / NTSPS

      !---------------------------------------------------------------------------------------
      !-- Condition to set the variable 'ISPXURB' for the domain of interest 
      !-- Urban model set for each domain according to namelist parameter ifuacm

        ISPXURB = .FALSE.
      IF (ifuacm.eq.1) ISPXURB = .TRUE.     

      !-----------------------------------------------------------------------------------
      ! Compute vegetation and land-use characteristics by land-use fraction weighting
      ! These parameters include LAI, VEGF, ZNT, ALBEDO, SNOALB, RS, etc.
      CALL VEGELAND(LANDUSEF, VEGFRA, SHDMIN, SHDMAX,         &
                    SOILCTOP, SOILCBOT, NLCAT, NSCAT,         &
                    ZNT,XLAI,XLAIMN,RSTMIN,XVEG,XVEGMN,XSNUP, &
                    XLAND, XALB,XSNOALB,WETFRA,IMPERV,CANFRA, &
                    URBF2D,ISPXURB,                           & !DDY, UTK Urban
                    ids,ide, jds,jde, kds,kde,                &
                    ims,ime, jms,jme, kms,kme,                &
                    its,ite, jts,jte, kts,kte, LAND_USE_TYPE, &
                    KWAT )
      !-----------------------------------------------------------------------------------
      !-----------------------------------------------------------------------------------
      ! Main loop over individual grid cells
      DO J = jts, jte    !-- J LOOP
        DO I = its, ite  !-- I LOOP

          !-- Estimation of height of full Model levels above ground level

          DO K = kts, kte
            DZF(K) = DZ8W(I,K,J)              
          ENDDO
          ZF(0) = 0.0


          DO K = kts, kte
            ZF(K) = DZF(K) + ZF(K-1)
            ZA(K) = 0.5 * ( ZF(K) + ZF(K-1) )
          ENDDO   
!           print *,' I,J,ZF(1),DZ8W(I,1,J),ts=',I,J,ZF(1),DZ8W(I,1,J),ITIMESTEP
          
          IFLAND = XLAND(I,J)

          ! Compute soil properties via weighting of fractional components
          IF (IFLAND .LT. 1.5 ) THEN 

          !---------------------------------------------------------
            CALL SOILPROP (SOILCBOT(I,:,J), WEIGHT,               &
                         ITIMESTEP, MAVAIL(I,J),                  &
                         PXLSM_SMOIS_INIT,                        &
                         FWSAT,FWFC,FWWLT,FCLAY,FCSAND,FFMSAND,   &
                         FB,FCGSAT,                               &  
                         FJP,FAS,FC2R,FC1SAT,FWRES, FC3, ISTI,    & 
                         SMOIS(I,1,J),  SMOIS(I,2,J)    )
          !----------------------------------------------------------
          !----------------------------------------------------------
            ISLTYP(I,J) = ISTI                      
          ELSE
            ISLTYP(I,J) = 14                    ! STATSGO type for water

            !-- aded for MODIS model
            FWWLT = 0.1
            FWFC = 1.0
            FWSAT = 1.0

            FCLAY = 0.0
            FCSAND = 0.0
            FFMSAND = 0.0
            !-- end

          ENDIF

          !-- aded for MODIS model
          WWLT_PX(I,J) = FWWLT
          WFC_PX(I,J) = FWFC
          WSAT_PX(I,J) = FWSAT

          CLAY_PX(I,J) = FCLAY * 0.01   ! percent to fraction
          CSAND_PX(I,J) = FCSAND * 0.01
          FMSAND_PX(I,J) = FFMSAND * 0.01
          !-- end

          ! initialize urban arrays for temporary hack
!          STANG_FNL(I,J) = 0.0
          VEGU_FNL(I,J)  = 0.0
          PPLN_FNL(I,J)  = 0.0
          LCZRUN = .FALSE.
          LCZtot = 0.
          sumhb = 0.0
          sumlp = 0.0
          sumbk1 = 0.0
          sumbk2 = 0.0

          IF (ISPXURB .AND. (URBF2D(I,J) .GT. 0.0)) THEN
            IF(LAND_USE_TYPE=='NLCD51')then
              STDIR1(I,J) = 05
              STDIR2(I,J) = 95               
!              print *,' nlcat,zf(1)=',i,j,URBF2D(I,J),nlcat,zf(1)
              Do lcat = 41, nlcat
                lcztot = lcztot + LANDUSEF(I,LCAT,J)
              enddo
              if(LCZtot.gt.0.01) then
                LCZRUN = .TRUE.
                Do lcat = 41, nlcat
                  lcatm40 = lcat-40
                  sumhb = sumhb + HB_LCZ(lcatm40)*LANDUSEF(I,LCAT,J)
                  sumlp = sumlp +LAMP_LCZ(lcatm40)*LANDUSEF(I,LCAT,J)
                  sumbk1 = sumbk1 + BLKW1_LCZ(lcatm40)*LANDUSEF(I,LCAT,J)
                  sumbk2 = sumbk2 + BLKW2_LCZ(lcatm40)*LANDUSEF(I,LCAT,J)
                ENDDO   
                HGTBDG(I,J) = sumhb/lcztot
                LAMP(I,J)   = sumlp/lcztot
                blkw1(i,J) =  sumbk1/lcztot
                blkw2(i,J) =  sumbk2/lcztot
              else
               HGTBDG(I,J) = (10.*LANDUSEF(I,24,J)+ 20.*LANDUSEF(I,25,J)+ 50.*LANDUSEF(I,26,J))/URBF2D(I,J)
               LAMP(I,J)   = (0.15*LANDUSEF(I,24,J)+ 0.25*LANDUSEF(I,25,J)+ 0.4*LANDUSEF(I,26,J))/URBF2D(I,J)
               blkw1(i,J) = (80.*LANDUSEF(I,24,J)+ 70.*LANDUSEF(I,25,J)+ 70.*LANDUSEF(I,26,J))/URBF2D(I,J)
               blkw2(i,J)  = (200*LANDUSEF(I,24,J)+ 180*LANDUSEF(I,25,J)+ 180.*LANDUSEF(I,26,J))/URBF2D(I,J)
!               WSTT(I,J) = (20.*LANDUSEF(I,24,J)+ 20.*LANDUSEF(I,25,J)+ 25.*LANDUSEF(I,26,J))/URBF2D(I,J)
    !             print *,'i.j.URBF2D,lcztot =',i,j,URBF2D(I,J),LCZtot
              endif
              BB1 = blkw1(i,J) +20.  ! WST = 20. which is subtract in SURFPX
              BB2 = blkw2(i,J) +20.
            ELSE   ! No LCZ
              IF(HGTBDG(I,J) .LT. 0.1) THEN
                HGTBDG(I,J) = (10.*LANDUSEF(I,24,J)+ 20.*LANDUSEF(I,25,J)+ 50.*LANDUSEF(I,26,J))/URBF2D(I,J)
                LAMP(I,J)   = (0.15*LANDUSEF(I,24,J)+ 0.25*LANDUSEF(I,25,J)+ 0.4*LANDUSEF(I,26,J))/URBF2D(I,J)
                LAMF(I,J)   = (0.08*LANDUSEF(I,24,J)+ 0.1*LANDUSEF(I,25,J)+ 0.4*LANDUSEF(I,26,J))/URBF2D(I,J)
              ENDIF
              BB1 = blkw1(i,J)*(1-LANDUSEF(I,KWAT,J))
              BB2 = blkw2(i,J)*(1-LANDUSEF(I,KWAT,J))
            ENDIF
            PPLN_FNL(I,J)  = HGTBDG(I,J) * LAMP(I,J) * 4444.4 ! p/km2
            VEGU_FNL(I,J)  = (0.15*LANDUSEF(I,24,J)+ 0.05*LANDUSEF(I,25,J)+ 0.01*LANDUSEF(I,26,J))/URBF2D(I,J)
          ENDIF

          !--  Variables Sub. SURFPX needs
          SFCPRS = PSFC(i,j) / 1000.0                 ! surface pressure in cb

          TA1    = T3D(i,1,j)                         ! air temperature at first layer
          DENS1  = RHO(I,1,J)                         ! air density at first layer

          QV1    = QV3D(i,1,j)                        ! water vapor mixing ratio at first layer
          QV2    = QV3D(i,2,j)

          ZLVL   = 0.5 * DZ8W(i,1,j)                  ! thickness of lowest half level
          ZF1    = DZ8W(i,1,j)
          ZA2    = ZF1 + 0.5 * DZ8W(i,2,j)

          LWDN   = GLW(I,J)                           ! longwave radiation
          EMISSI = EMISS(I,J)                         ! emissivity
          PRECIP = MAX ( 1.0E-3*RAINBL(i,j)/DTBL,0.0) ! accumulated precip. rate during DT (=dtpbl)
                                                      ! convert RAINBL from mm to m for PXLSM
          WR     = 1.0E-3*CANWAT(I,J)                 ! convert CANWAT from mm to m for PXLSM
          THETA1 = TH3D(i,1,j)                        ! potential temp at first layer

          SNOBS  = SNOW(I,J)                          ! Set snow cover to existing model value
                                                      !   this is overwritten below if snow analysis is availiable
                                                      !   otherwise snow cover remains constant through simulation

          IF (ITIMESTEP.EQ.1) THEN
            
            STWAT(I,J) = CANWAT(I,J)
            RFWAT(I,J) = CANWAT(I,J) 
            TG_RF(I,J) = TSLB(I,1,J)                ! Initialization of the surface and bottom layer temperatures 
            T2_RF(I,J) = TSLB(I,2,J)                !  for the Roof and Walls with that of soil layers    
            TG_WL(I,J) = TSLB(I,1,J)
            T2_WL(I,J) = TSLB(I,2,J)
            T2_GND(I,J) = TSLB(I,2,J)
          ENDIF     

          CNP=-1
          !-- Estimation of the urban canopy height level (Model-level)
          IF (ISPXURB .AND. (URBF2D(I,J) .GT. 0.0)) THEN 
            DO K = kts, kte
              IF ((HGTBDG(I,J) .GE. ZF(K-1)) .AND. (HGTBDG(I,J) .LE. ZF(K))) THEN 
                  KSL(I,J) = K
              ENDIF
            ENDDO
            if(ksl(i,j).lt.1) Then
                  PRINT *, ' HGTBDG,ksl=',HGTBDG(I,J),KSL(i,j)
                 PRINT *, ' ZF0-4=',ZF(0) ,ZF(1),ZF(2),ZF(3),ZF(4)
                 PRINT *, ' DZF1-4=',DZ8w(i,1,j),DZ8w(i,2,j),DZF(3),DZF(4)
            ENDIF
          ENDIF
          
          IF (ISPXURB .AND. (URBF2D(I,J) .GT. 0.0)) THEN    ! For UTK-Urban 
            cnp = ksl(i,j)
            if(cnp.lt.1) PRINT *, ' URBF2D(I,J),CNP,HGTBDG(I,J)=',i,j,URBF2D(I,J),CNP,HGTBDG(I,J),KSL(i,j)
            WR_VEG    = (1.0E-3*CANWAT(I,J))
            WR_ST     = (1.0E-3*STWAT(I,J))     
            WR_RF     = (1.0E-3*RFWAT(I,J))
            DENSRF    = RHO(i,CNP,j)       
            QV_RF     = QV3D(i,CNP,j)                 ! water-vapor mixing ratio at urban canopy layer
            QV_RF1    = QV3D(i,CNP+1,j)               ! water-vapor mixing ratio at one layer above the urban canopy layer
            TA_RF     = T3D(i,CNP,j)                  ! air temperature at roof layer
            THETA_RF  = TH3D(i,CNP,j)                 ! potential temp at urban canopy layer
            THETA_RF1 = TH3D(i,CNP+1,j)               ! potential temp at one layer above the urban canopy layer
            DZFRF     = DZF(CNP)                      ! vertical layer thickness

            US_RF1 = U3D(i,CNP+1,j)                   ! Horz. wind components at one layer above the urban canopy level
            VS_RF1 = V3D(i,CNP+1,j)
            US_RF  = U3D(i,CNP,j)                     ! Horz. wind components at urban canopy level
            VS_RF  = V3D(i,CNP,j)        
            PP_RF  = PP3D(i,CNP,j) / 1000.0           ! Pressure at urban canopy level in cb (kPa) 

            THETAV_RF  = THETA_RF * (1.0 + EP1 * QV_RF)
            THETAV_RF1 = THETA_RF1 * (1.0 + EP1 * QV_RF1)  

            THETA_RFs  = TG_RF(I,J) * ((100.0 / PP_RF)**0.285)     ! Conversion of Roof- surface temp. to Pot. temp.
            THETAV_RFs = THETA_RFs * (1.0 + EP1 * QV_RF)

            ZF_CNP1 = ZF(CNP)
            ZF_CNP0 = ZF(CNP-1)

            ZA_CNP2 = ZA(CNP+1) 
            ZA_CNP1 = ZA(CNP)

            !-- Urban canopy calculations for walls
            DENS_CAN  = 0.0
            THETA_CAN = 0.0
            TA_CAN    = 0.0
            US_CAN    = 0.0
            VS_CAN    = 0.0
            WSP_CAN(I,J)   = 0.0 
            QV_CAN    = 0.0
            PP_CAN    = 0.0
            DO K = kts,CNP
              DENS_CAN  = DENS_CAN + RHO(i,K,j)  
              THETA_CAN = THETA_CAN + TH3D(i,K,j)
              TA_CAN    = TA_CAN + T3D(i,K,j)
              US_CAN    = U3D(i,K,j)
              VS_CAN    = V3D(i,K,j)
              WSP_CAN(I,J)   = WSP_CAN(I,J) + MAX(SQRT(US_CAN**2 + VS_CAN**2), 0.01) 
              QV_CAN    = QV_CAN + QV3D(i,K,j)
              PP_CAN    = PP_CAN + (PP3D(i,K,j) / 1000.0)          ! Air pressure in kPa
            ENDDO
            
            DENS_CAN  = DENS_CAN/CNP                  ! Average air density within the urban canopy
            THETA_CAN = THETA_CAN/CNP                 ! Average air potential temperature (K) within the urban canopy
            TA_CAN    = TA_CAN/CNP                    ! Average air temperature (K) within the urban canopy
            WSP_CAN(I,J)   = WSP_CAN(I,J)/CNP                   ! Average horz. wind speed within the urban canopy
            QV_CAN    = QV_CAN/CNP                    ! Average water-vapor mixing ratio within the urban canopy
            PP_CAN    = PP_CAN/CNP                    ! Average air pressure within the urban canopy

            THETA_WLs = TG_WL(I,J) * ((100.0 / PP_CAN)**0.285)     ! Conversion of Walls- surface temp. to Pot. temp.
          ENDIF  ! urban


          IF(PXLSM_SOIL_NUDGE .EQ. 1) THEN                   
            !-- 2 m Temp and RH for Nudging
            T2_OLD = T2_NDG_OLD(I,J)
            T2_NEW     = T2_NDG_NEW(I,J)
            VAPPRS     = SVP1 * EXP(SVP2 * (T2_OLD - SVPT0) / ( T2_OLD - SVP3))
            QSBT       = EP_2 * VAPPRS / (SFCPRS - VAPPRS)          
            RH2_OLD    = Q2_NDG_OLD(I,J) / QSBT
            VAPPRS     = SVP1 * EXP(SVP2 * (T2_NEW - SVPT0) / (T2_NEW - SVP3))
            QSBT       = EP_2 * VAPPRS / (SFCPRS - VAPPRS)          
            RH2_NEW    = Q2_NDG_NEW(I,J) / QSBT
            RH2OBS     = CORB * RH2_OLD +  CORE * RH2_NEW  
            T2OBS(I,J) = CORB * T2_OLD +  CORE * T2_NEW
            Q2OBS(I,J) = CORB * Q2_NDG_OLD(I,J) +  CORE * Q2_NDG_NEW(I,J)
            SNOBS = CORB * SN_NDG_OLD(I,J) +  CORE * SN_NDG_NEW(I,J)  
          ENDIF

          USTAR  = MAX(UST(I,J),0.005)
          
          IF (IFLAND .GE. 1.5) THEN ! if over water                            
            ZNT(I,J) = CZO * UST(I,J) * UST(I,J) / G + OZO
          ENDIF                                                              
          
          Z0       = ZNT(I,J)
          CPAIR    = CPD * (1.0 + 0.84 * QV1)            ! J/(K KG)

          IF (ISPXURB .AND. (URBF2D(I,J) .GT. 0.0)) THEN   ! For UTK-Urban     
            CPAIR_RF  = CPD * (1.0 + 0.84 * QV_RF)             ! J/(K KG)
            CPAIR_CAN = CPD * (1.0 + 0.84 * QV_CAN)            ! J/(K KG) 
          ENDIF

          ! Set WRF Snow albedo to PX snow albedo based on fractional landuse
          ! Snow albedo for each landuse class is defined in module_sf_pxlsm_data.F
          SNOALB(I,J) = XSNOALB(I,J)
          !-------------------------------------------------------------
          ! Compute fractional snow area and snow albedo
          CALL PXSNOW (ITIMESTEP, SNOBS, SNOWNCV(I,J), SNOW(I,J),  &
                       SNOWH(I,J), XSNUP(I,J),  XALB(i,j),         &
                       SNOALB(I,J),VEGF_PX(I,J), SHDMIN(I,J),      &
                       HC_SNOW, SNOW_FRA, SNOWC(I,J),  ALBEDO(I,J) ) 
          !-------------------------------------------------------------
                                 
          !-------------------------------------------------------------
          ! Sea Ice from analysis and water cells that are very cold, but more than 50% water
          ! are converted to ice/snow for more reasonable treatment.
          IF( (XICE(I,J).GE.0.5) .OR.   &
              (SST(I,J).LE.270.0.AND.XLAND(I,J).GE.1.50) ) THEN
            XLAND(I,J)   = 1.0
            IFLAND       = 1.0
            ZNT(I,J)     = 0.001  
            SMOIS(I,1,J) = 1.0 
            SMOIS(I,2,J) = 1.0 
            XICE(I,J)    = 1.0
            ALBEDO(I,J)  = 0.7
            SNOWC(I,J)   = 1.0
            SNOW_FRA     = 1.0
            VEGF_PX(I,J) = 0.0
            LAI_PX(I,J)  = 0.0
            LAI(I,J)     = 0.0
            FCGSAT       = 3.670
          ENDIF
          !-------------------------------------------------------------

          !-------------------------------------------------------------
          !-- Note that when IFGROW = 0 is selected in Vegeland then max and min           
          !-- LAI and Veg are the same                                                     
          T2I = TSLB(I,2,J)                                                       
          FSEAS = AMAX1(1.0 - 0.015625 * (290.0 - T2I) ** 2,0.0)           
          IF (T2I .GE. 290.0) FSEAS = 1.0                                          
          
          !get PX table vegetation
          LAI_PX(I,J)    = XLAIMN(I,J) + FSEAS*(XLAI(I,J) - XLAIMN(I,J))
          VEGF_PX(I,J)   = XVEGMN(I,J) + FSEAS*(XVEG(I,J) - XVEGMN(I,J))                       
          
!... use MODIS LAI and VEGFRA from wrflowinp
          IF ( pxlsm_modis_veg .EQ. 1 ) THEN

                ! get LAI for vegetated area
                IF ( VEGFRA(I,J) .GT. 0.0 )  THEN
                   LAI_PX(I,J) = LAI(I,J) / ( VEGFRA(I,J) / 100.0)
                ELSE
                   LAI_PX(I,J) = 0.0
                ENDIF

                VEGF_PX(I,J) = VEGFRA(I,J) / 100.0

                !vegF is just for the land
                IF ( LANDUSEF(I,KWAT,J) .LT. 1.0 )  THEN
                   VEGF_PX(I,J) = VEGF_PX(I,J) / (1.0 - LANDUSEF(I,KWAT,J))
                ELSE
                   VEGF_PX(I,J) = 0.0
                ENDIF

          ENDIF

          LAI_PX(I,J) = MIN(LAI_PX(I,J), 8.0)
          LAI_PX(I,J) = MAX(LAI_PX(I,J), 0.0001)

          VEGF_PX(I,J) = MIN(VEGF_PX(I,J), 1.0)
          VEGF_PX(I,J) = MAX(VEGF_PX(I,J), 0.0001)

!... END OF MODIS LAI and FPAR


          ! Ensure veg algorithms not used for water 
          IF (IFLAND .GE. 1.5) THEN                        
             VEGF_PX(I,J) = 0.0
             LAI_PX(I,J) = 0.0         
          ENDIF                                                                   
          !-------------------------------------------------------------


          SOLDN = GSW(I,J) / (1.0 - ALBEDO(I,J))     ! downward shortwave radiaton

          ISNOW = SNOWC(I,J)

          NUDGE = PXLSM_SOIL_NUDGE
          IF ( J .LE. 2 .OR. J .GE. (jde-1) ) NUDGE = 0
          IF ( I .LE. 2 .OR. I .GE. (ide-1) ) NUDGE = 0
          
          IF ( RMOL(I,J) .GT. 0.0 )  THEN
              MOLX = AMIN1(1/RMOL(I,J),1000.0)
          ELSE IF ( RMOL(I,J) .LT. 0.0 ) THEN
              MOLX = AMAX1(1/RMOL(I,J),-1000.0)
          ELSE
              MOLX = 1000
          ENDIF   
 
          ZFUNC = ZF1 * (1.0 - ZF1 / MAX(100.,PBLH(I,J))) ** 2
          QST12 = KARMAN * ZFUNC*(QV2-QV1) / (ZA2-ZLVL)
               

          !-------------------------------------------------------------
          !-- LSM sub-time loop too prevent dt > 40 sec       
          NTSPS = INT(DT / (DTPBLX + 0.000001) + 1.0)                                
          DTPBL = DT / NTSPS                                                       

          DO IT = 1, NTSPS                                                     
          
            !... SATURATION VAPOR PRESSURE (MB) 
            IF ( TSLB(I,1,J) .LE. SVPT0 ) THEN        ! For ground that is below freezing
              ES = SVP1 * EXP(22.514 - 6.15E3 / TSLB(I,1,J))      ! cb
            ELSE
              ES = SVP1 * EXP(SVP2 * (TSLB(I,1,J) - SVPT0) /  (TSLB(I,1,J) - SVP3))
            ENDIF
            QSS  = ES * 0.622 / (SFCPRS - ES)

            !...Saturation vapor pressure and saturation mixing ratio at Roof for Urban grid cell
            IF (ISPXURB .AND. (URBF2D(I,J) .GT. 0.0)) THEN 
              ES     = SVP1 * EXP(SVP2 * (TG_RF(I,J) - SVPT0) / (TG_RF(I,J) - SVP3)) 
              QSS_RF = ES * 0.622 / (PP_RF - ES) 
            ENDIF
          
            !... beta method, Lee & Pielke (JAM,May1992)
            BETAP = 1.0
            IF (IFLAND .LT. 1.5 .AND. ISNOW .LT. 0.5 .AND. SMOIS(I,1,J) .LE. FWFC) THEN       
              BETAP = 0.25 * (1.0 - COS(SMOIS(I,1,J) / FWFC * PI)) ** 2
            ENDIF
            
            !-------------------------------------------------------------------------
            CALL SURFPX(ITIMESTEP,DTPBL, ISPXURB, IFLAND, SNOWC(I,J),  NUDGE, XICE(I,J),          & !in
                        JULIAN, XTIME, GMT, SOLDN,  GSW(I,J), LWDN,   EMISSI, ZLVL,     & !in
                        MOLX, Z0, USTAR, ZA_CNP2, ZA_CNP1, ZF_CNP1, ZF_CNP0, DZFRF,     & !in
                        SFCPRS, DENS1, DENSRF, QV1, QV_RF, QSS, QSS_RF, TA1, TA_RF,     & !in --UTK-Urban
                        THETA1, THETA_RF, PRECIP, THETAV_RF, THETA_RF1, THETAV_RF1,     & !in --UTK-Urban
                        THETA_RFs, THETAV_RFs, THETA_WLs,                               & !in --UTK-Urban
                        CPAIR, CPAIR_RF, PSIH(I,J), US_RF1, US_RF, VS_RF1, VS_RF,       & !in --UTK-Urban
                        DENS_CAN, THETA_CAN, TA_CAN, WSP_CAN(I,J), CPAIR_CAN,           & !in --UTK-Urban
                        RH2OBS, T2OBS(I,J),                                             & !in
                        VEGF_PX(I,J), ISTI, LAI_PX(I,J), IMPERV(I,J), CANFRA(I,J),      & !in
                        BETAP, RSTMIN(I,J), HC_SNOW, SNOW_FRA, WETFRA(I,J),             & !in          
                        FWWLT, FWFC, FWRES, FCGSAT,  FWSAT, FB,                         & !in
                        FC1SAT, FC2R, FAS, FJP, FC3, DZS(1), DZS(2), QST12,             & !in
                        URBF2D(I,J), COSZ(I,J), CURR_SECS, ALBEDO(I,J),                 & !-- DDY and UTK Urban in
                        XLAT(I,J), XLONG(I,J), SINA(I,J), COSA(I,J),                    & !in --Urban
                        RADNET(I,J), GRDFLX(I,J), HFX(I,J), QFX(I,J), LH(I,J),          & !out
                        HFX_RF(I,J), HFX_WL(I,J), ANTHFX1(I,J), ANTHFX_RF(I,J), QFX_RF(I,J), LH_RF(I,J),   & !out --UTK-Urban
                        EG(I,J), ER(I,J), ETR(I,J), ER_ST(I,J), ER_RF(I,J),                                & !out --UTK-Urban
                        QST(I,J), CAPG(I,J), RS(I,J), RA(I,J),                                             & !out
                        TSLB(I,1,J), T2_GND(I,J), TG_RF(I,J), T2_RF(I,J), TG_WL(I,J), T2_WL(I,J),          & !out --UTK-Urban
                        RADWL(I,J), RADWL_FNL(I,J), RADST(I,J), RADST_FNL(I,J), GSW_URB(I,J), GSW_WL(I,J), & !out --UTK-Urban
                        RADNET_URB(I,J), RADNET_RF(I,J), RADNET_WL(I,J), LWDN_URB(I,J), LWDN_RF(I,J),      & !out --UTK-Urban
                        LWDN_WL(I,J), LWUP_URB(I,J), LWUP_RF(I,J), LWUP_WL(I,J),                           & !out --UTK-Urban
                        SMOIS(I,1,J), SMOIS(I,2,J), WR, WR_ST, WR_VEG, WR_RF,                              & !out --UTK-Urban
                        TA2(I,J), QA2(I,J), LAND_USE_TYPE, I, J,                                           &
                        HGTBDG(I,J), LAMP(I,J), LAMF(I,J), STDIR1(I,J), STDIR2(I,J), BB1, BB2,             & !    --UTK Urban
                        PPLN_FNL(I,J), VEGU_FNL(I,J), KSL(I,J), SWDDIR(I,J), SWDDIF(I,J),                  &
                        GAB_RF(I,J),GAB_WL(I,J)                                                            ) !    --UTK Urban
            !-------------------------------------------------------------------------
          
          END DO                        ! Time internal PX time loop   

          IF (IFLAND .GE. 1.5) THEN                        
            TSK(I,J)     = SST(I,J)     ! Skin temp set to sea surface temperature for open water
            GRDFLX(I,J)  = 0.0
          ELSE
            TSK(I,J)  = TSLB(I,1,J)     ! Skin temp set to 1 cm soil temperature in PX for now
          ENDIF
          IF (ISPXURB .AND. (URBF2D(I,J) .GT. 0.0)) THEN  ! For UTK-Urban
             CANWAT(I,J) =  WR_VEG * 1000.    ! convert WR back to mm for CANWAT
             STWAT(I,J)  =  WR_ST * 1000.     ! convert WR back to mm for STWAT       
             RFWAT(I,J)  =  WR_RF * 1000.     ! convert WR back to mm for RFWAT
          ELSE
             CANWAT(I,J) =  WR * 1000.        ! convert WR back to mm for CANWAT
          ENDIF
          RAW        = RA(I,J) + 4.503 / USTAR
          QSFC(I,J)  = QFX(I,J) * RAW / DENS1 + QV1
          
        ENDDO !  END MAIN I LOOP
      ENDDO   !  END MAIN J LOOP
      
!------------------------------------------------------
   END SUBROUTINE pxlsm
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
      SUBROUTINE VEGELAND(landusef, vegfra,                            & 
                          shdmin, shdmax,                              &                          
                          soilctop, soilcbot, nlcat, nscat,znt, xlai,  &
                          xlaimn, rstmin, xveg, xvegmn, xsnup, xland,  &
                          xalb, xsnoalb, wetfra, imperv, canfra,       &
                          urbf2d, ispxurb,                             & !-- DDY, UTK
                          ids,ide, jds,jde, kds,kde,                   &
                          ims,ime, jms,jme, kms,kme,                   &
                          its,ite, jts,jte, kts,kte,                   &
                          LAND_USE_TYPE, KWAT_OUT                      )
!-------------------------------------------------------------------------
!           
!   CALLED FROM Sub. bl_init in module_physics.init.F
!           
!   THIS PROGRAM PROCESSES THE USGS LANDUSE DATA
!   WHICH HAS BEEN GRIDDED BY THE WPS SYSTEM
!   AND PRODUCES 2-D FIELDS OF LU RELATED PARAMETERS
!   FOR USE IN THE PX SURFACE MODEL
!
!
! REVISION HISTORY:
!  AX      Oct 2004 - developed WRF version based on VEGELAND in the MM5 PX LSM
!  RG      Dec 2006 - revised for WRFV2.1.2
!  JP      Dec 2007 - revised for WRFV3 - removed IFGROW options
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
                              
      IMPLICIT NONE
!... 
      INTEGER, INTENT(IN)        ::  ids,ide, jds,jde, kds,kde,  &
                                     ims,ime, jms,jme, kms,kme,  &
                                     its,ite, jts,jte, kts,kte
                                                                            
      INTEGER, INTENT(IN)        ::  NSCAT, NLCAT

      REAL, DIMENSION( ims:ime, 1:NLCAT, jms:jme ), INTENT(IN) :: LANDUSEF                        
      REAL, DIMENSION( ims:ime, 1:NSCAT, jms:jme ), INTENT(IN) :: SOILCTOP, SOILCBOT
                                                                   
      REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)    :: VEGFRA, SHDMIN, SHDMAX 

      REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: ZNT, IMPERV, CANFRA, XLAND

      REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)   :: XLAI, XLAIMN, RSTMIN, XALB, &
                                                            XVEG, XVEGMN, XSNUP,        &
                                                            WETFRA, XSNOALB 

      INTEGER, INTENT(OUT)            :: KWAT_OUT

      CHARACTER (LEN = 6), INTENT(IN) :: LAND_USE_TYPE

      !-UTK Urban 
      LOGICAL, INTENT(IN) :: ISPXURB
      REAL, PARAMETER     :: Z0ST = 30.0            ! Urban street surface roughness length (cm)

                                     
!... local variables

      INTEGER :: ITF, JTF, K, J, I
      REAL    :: SUMLAI, SUMLMN, SUMRSI, SUMLZ0, SUMVEG, SUMVMN, &
                 ALAI, VEGF, SUMSNUP, SUMALB, SUMSNOALB
                 
      REAL    :: VFMX, VFMN, VSEAS, FAREA, FWAT, ZNOTC, FCAN, FIMP, FORFRA, EXTFOR

      REAL, DIMENSION( NLCAT ) :: LAIMX, LAIMN, Z0, VEG, VEGMN, SNUP, ALB, SNOALB

      REAL, PARAMETER :: ZNOTCMN = 5.0  ! CM, MIN Zo FOR CROPS
      REAL, PARAMETER :: ZNOTCMX = 15.0 ! CM, MAX Zo FOR CROPS

      REAL, SAVE, DIMENSION(:), POINTER  :: RSMIN, Z00, VEG0, VEGMN0, LAI0, &
                                            LAIMN0, SNUP0, ALBF, SNOALBF

      !---- INITIALIZE PARAMETERS
      INTEGER, SAVE :: KWAT, LIMIT1, LIMIT2

      !-- DDY urban
      REAL, DIMENSION(its:ite, jts:jte), INTENT(INOUT) :: URBF2D
      REAL :: LUTOT, URBF
      INTEGER NLCATX
      

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------

      NLCATX = NLCAT
     ! Initialize LU characteristics by LU Dataset
      IF (LAND_USE_TYPE == 'USGS') THEN
         KWAT = 16
         RSMIN  => RSMIN_USGS
         Z00    => Z00_USGS  
         VEG0   => VEG0_USGS 
         VEGMN0 => VEGMN0_USGS 
         LAI0   => LAI0_USGS 
         LAIMN0 => LAIMN0_USGS 
         SNUP0  => SNUP0_USGS 
         ALBF   => ALBF_USGS
         SNOALBF=> SNOALB_USGS
         LIMIT1 = 2
         LIMIT1 = 6
      ELSE IF (LAND_USE_TYPE == 'USGS28') THEN
         KWAT = 16
         RSMIN  => RSMIN_USGS28
         Z00    => Z00_USGS28  
         VEG0   => VEG0_USGS28 
         VEGMN0 => VEGMN0_USGS28 
         LAI0   => LAI0_USGS28 
         LAIMN0 => LAIMN0_USGS28 
         SNUP0  => SNUP0_USGS28 
         ALBF   => ALBF_USGS28
         SNOALBF=> SNOALB_USGS28
         LIMIT1 = 2
         LIMIT1 = 6
      ELSE IF (LAND_USE_TYPE == 'NLCD50') THEN
         KWAT = 1
         RSMIN  => RSMIN_NLCD50
         Z00    => Z00_NLCD50  
         VEG0   => VEG0_NLCD50 
         VEGMN0 => VEGMN0_NLCD50 
         LAI0   => LAI0_NLCD50 
         LAIMN0 => LAIMN0_NLCD50 
         SNUP0  => SNUP0_NLCD50 
         ALBF   => ALBF_NLCD50
         SNOALBF=> SNOALB_NLCD50
         LIMIT1 = 20
         LIMIT1 = 43
      ELSE IF (LAND_USE_TYPE == 'NLCD40') THEN
         KWAT = 17
         RSMIN  => RSMIN_NLCD40
         Z00    => Z00_NLCD40  
         VEG0   => VEG0_NLCD40 
         VEGMN0 => VEGMN0_NLCD40 
         LAI0   => LAI0_NLCD40 
         LAIMN0 => LAIMN0_NLCD40 
         SNUP0  => SNUP0_NLCD40 
         ALBF   => ALBF_NLCD40
         SNOALBF=> SNOALB_NLCD40
         LIMIT1 = 20
         LIMIT1 = 43
      ELSE IF (LAND_USE_TYPE == 'NLCD51') THEN
         KWAT = 17
         RSMIN  => RSMIN_NLCD51
         Z00    => Z00_NLCD51  
         VEG0   => VEG0_NLCD51 
         VEGMN0 => VEGMN0_NLCD51 
         LAI0   => LAI0_NLCD51 
         LAIMN0 => LAIMN0_NLCD51 
         SNUP0  => SNUP0_NLCD51 
         ALBF   => ALBF_NLCD51
         SNOALBF=> SNOALB_NLCD51
         LIMIT1 = 20
         LIMIT1 = 43
         NLCATX  = 40
      ELSE IF (LAND_USE_TYPE == 'MODIS') THEN
         KWAT = 17
         RSMIN  => RSMIN_MODIS
         Z00    => Z00_MODIS  
         VEG0   => VEG0_MODIS 
         VEGMN0 => VEGMN0_MODIS 
         LAI0   => LAI0_MODIS 
         LAIMN0 => LAIMN0_MODIS 
         SNUP0  => SNUP0_MODIS 
         ALBF   => ALBF_MODIS
         SNOALBF=> SNOALB_MODIS
         LIMIT1 = 12
         LIMIT1 = 14
      END IF

      KWAT_OUT = KWAT
      !--------------------------------------------------------------------
      DO J = jts,jte
        DO I = its,ite
          XLAI(I,J)   = 0.0                                                        
          XLAIMN(I,J) = 0.0                                                     
          RSTMIN(I,J) = 9999.0                                                    
          XVEG(I,J)   = 0.0                                                       
          XVEGMN(I,J) = 0.0
          XSNUP(I,J)  = 0.0                                                            
          XALB(I,J)   = 0.0                                                            
          XSNOALB(I,J)= 0.0                                                            

          ! Code that may be needed in case these arrays are not intialized
          ! with zero by real.exe when not defined by GEOGRID or present 
          ! in met_em* files
          !IMPERV(I,J) = AMAX1(0.0001,IMPERV(I,J))                                                           
          !CANFRA(I,J) = AMAX1(0.0001,CANFRA(I,J))                                                           

        ENDDO   ! END LOOP THROUGH GRID CELLS
      ENDDO     ! END LOOP THROUGH GRID CELLS
      !--------------------------------------------------------------------

      DO J = jts,jte
        DO I = its,ite
          !-- DDY urban, calculate urban fraction
          URBF2D(I,J) = -1.
          LUTOT = 0.0
          DO K = 1,NLCATX
            LUTOT = LUTOT + LANDUSEF(I,K,J)
          ENDDO
          IF (LUTOT .GT. 0.0) THEN
            URBF = (LANDUSEF(I,24,J)+LANDUSEF(I,25,J)+LANDUSEF(I,26,J))/LUTOT  !-Only for NLCD40 here.
            IF(URBF.GT.0.05) URBF2D(I,J) = URBF
          ENDIF

          !-- Initialize 2 and 3-D veg parameters to be caculated
          DO K = 1, NLCAT
            LAIMX(K)  = LAI0(K)                                                       
            LAIMN(K)  = LAIMN0(K)                                                   
            Z0(K)     = Z00(K)                                                   
            VEG(K)    = VEG0(K)                                              
            VEGMN(K)  = VEGMN0(K)
            SNUP(K)   = SNUP0(K)                                            
            ALB(K)    = ALBF(K)                                            
            SNOALB(K) = SNOALBF(K)                                            
          ENDDO                                              

          !-- Set Urban street roughness length here to 'Z0ST' (in cm.) for the domain of interest  !--UTK urban
          IF (ISPXURB .AND. (URBF2D(I,J) .GT. 0.0)) THEN
            IF (LAND_USE_TYPE == 'USGS') THEN           !--For USGS (24-Cat) landuse here. 
              Z0(1)  = Z0ST 
            ELSE IF (LAND_USE_TYPE == 'USGS28') THEN
              Z0(1)  = Z0ST
            ELSE  !  NLCD40
              Z0(24) = Z0ST
              Z0(25) = Z0ST
              Z0(26) = Z0ST
            ENDIF
          ENDIF

          !--  INITIALIZE SUMS
          SUMLAI    = 0.0
          SUMLMN    = 0.0
          SUMRSI    = 0.0
          SUMLZ0    = 0.0
          SUMVEG    = 0.0
          SUMVMN    = 0.0
          ALAI      = 0.0
          SUMSNUP   = 0.0
          SUMALB    = 0.0
          SUMSNOALB = 0.0

          !--   ESTIMATE CROP EMERGANCE DATE FROM VEGFRAC
          VFMX   = SHDMAX(I,J)
          VFMN   = SHDMIN(I,J)
          VEGF   = VEGFRA(I,J) 
                
          !-- Computations for VEGETATION CELLS ONLY
          IF(VFMX.GT.0.0.AND.LANDUSEF(I,KWAT,J).LT.1.00) THEN
            VSEAS = VEGF/VFMX
            IF(VSEAS.GT.1.0.OR.VSEAS.LT.0.0) THEN
              VSEAS = MIN(VSEAS,1.0)
              VSEAS = MAX(VSEAS,0.0)
            ENDIF

            ZNOTC = ZNOTCMN * (1-VSEAS) + ZNOTCMX * VSEAS  ! Zo FOR CROPS
            DO K = 1, NLCATX
              IF (LAND_USE_TYPE == 'MODIS') THEN
                !-- USE THE VEGFRAC DATA ONLY FOR CROPS
                IF (K.EQ.12 .OR. K.EQ.14) THEN
                   LAIMX(K) = LAIMN0(K) * (1-VSEAS) + LAI0(K) * VSEAS
                   LAIMN(K) = LAIMX(K)
                   VEG(K)   = VEGMN0(K) * (1-VSEAS) + VEG0(K) * VSEAS
                   VEGMN(K) = VEG(K)
                   !-- SEASONALLY VARY Zo FOR MODIS DryCrop (k=12)
                   IF (K .EQ. 12) THEN
                      Z0(K) = ZNOTC
                   !-- CrGrM (k=14) USE AVG WITH GRASS AND FOREST
                   ELSE IF (K .EQ.14) THEN
                      Z0(K)  = 0.5 * (ZNOTC + Z00(K))
                   ENDIF
                ENDIF
              ELSE IF (LAND_USE_TYPE == 'NLCD50') THEN
                !-- USE THE VEGFRAC DATA ONLY FOR CROPS
                IF (K.EQ.20 .OR. K.EQ.43 .OR. K.EQ.45) THEN
                   LAIMX(K) = LAIMN0(K) * (1-VSEAS) + LAI0(K) * VSEAS
                   LAIMN(K) = LAIMX(K)
                   VEG(K)   = VEGMN0(K) * (1-VSEAS) + VEG0(K) * VSEAS
                   VEGMN(K) = VEG(K)
                   !-- SEASONALLY VARY Zo FOR DryCrop (k=20) OR Irigated Crop (k=43) OR  Mix Crop (k=4)
                   IF (K.EQ.20 .OR. K.EQ.43) THEN
                      Z0(K) = ZNOTC
                   !-- CrNatM (k=45) USE AVG WITH GRASS AND FOREST
                  ELSE IF (K.EQ.45) THEN
                      Z0(K)  = 0.5 * (ZNOTC + Z00(K))
                  ENDIF
                ENDIF
              ELSE IF (LAND_USE_TYPE == 'NLCD40'.or.LAND_USE_TYPE == 'NLCD51') THEN
                !-- USE THE VEGFRAC DATA ONLY FOR CROPS
                IF (K.EQ.12 .OR. K.EQ.14 .OR. K.EQ.38) THEN
                   LAIMX(K) = LAIMN0(K) * (1-VSEAS) + LAI0(K) * VSEAS
                   LAIMN(K) = LAIMX(K)
                   VEG(K)   = VEGMN0(K) * (1-VSEAS) + VEG0(K) * VSEAS
                   VEGMN(K) = VEG(K)
                   !-- SEASONALLY VARY Zo FOR Crop (k=12 for MODIS or 38 for NLCD)
                   IF (K.EQ.12 .OR. K.EQ.38) THEN
                      Z0(K) = ZNOTC
                   !-- CrNatM (k=14) USE AVG WITH GRASS AND FOREST
                   ELSE IF (K.EQ.14) THEN
                      Z0(K)  = 0.5 * (ZNOTC + Z00(K))
                  ENDIF
                ENDIF
              ELSE IF (LAND_USE_TYPE == 'USGS') THEN
                !-- USE THE VEGFRAC DATA ONLY FOR CROPS 
                IF (K .GE. 2 .AND. K .LE. 6) THEN
                   LAIMX(K) = LAIMN0(K) * (1-VSEAS) + LAI0(K) * VSEAS
                   LAIMN(K) = LAIMX(K)
                   VEG(K)   = VEGMN0(K) * (1-VSEAS) + VEG0(K) * VSEAS
                   VEGMN(K) = VEG(K)
                   !-- SEASONALLY VARY Zo FOR DryCrop (k=2) OR Irigated Crop (k=3) OR  Mix Crop (k=4)
                   IF (K .GE. 2 .AND. K .LE. 4) THEN
                      Z0(K) = ZNOTC
                   !-- CrGrM (k=5) or CrWdM (k=6) USE AVG WITH GRASS AND FOREST
                   ELSE IF (K .GE.5 .AND. K .LE. 6) THEN
                      Z0(K)  = 0.5 * (ZNOTC + Z00(K))
                   ENDIF
                ENDIF
              ELSE IF (LAND_USE_TYPE == 'USGS28') THEN
                !-- USE THE VEGFRAC DATA ONLY FOR CROPS 
                IF (K .GE. 2 .AND. K .LE. 6) THEN
                   LAIMX(K) = LAIMN0(K) * (1-VSEAS) + LAI0(K) * VSEAS
                   LAIMN(K) = LAIMX(K)
                   VEG(K)   = VEGMN0(K) * (1-VSEAS) + VEG0(K) * VSEAS
                   VEGMN(K) = VEG(K)
                   !-- SEASONALLY VARY Zo FOR DryCrop (k=2) OR Irigated Crop (k=3) OR  Mix Crop (k=4)
                   IF (K .GE. 2 .AND. K .LE. 4) THEN
                      Z0(K) = ZNOTC
                   !-- CrGrM (k=5) or CrWdM (k=6) USE AVG WITH GRASS AND FOREST
                   ELSE IF (K .GE.5 .AND. K .LE. 6) THEN
                      Z0(K)  = 0.5 * (ZNOTC + Z00(K))
                   ENDIF
                ENDIF

              END IF

            ENDDO

          ENDIF       !-- IF cell is vegetation

          !-------------------------------------
          !-- LOOP THROUGH LANDUSE Fraction and compute totals
          DO K = 1, NLCATX 
            FAREA     = LANDUSEF(I,K,J)
            SUMLAI    = SUMLAI + LAIMX(K) * FAREA
            SUMLMN    = SUMLMN + LAIMN(K) * FAREA
            ALAI      = ALAI + FAREA
            SUMRSI    = SUMRSI + FAREA * LAIMX(K) / RSMIN(K)
            SUMLZ0    = SUMLZ0 + FAREA * ALOG(Z0(K))
            SUMVEG    = SUMVEG + FAREA * VEG(K)
            SUMVMN    = SUMVMN + FAREA * VEGMN(K)
            SUMSNUP   = SUMSNUP+ FAREA * SNUP(K)
            SUMALB    = SUMALB + FAREA * ALB(K)
            SUMSNOALB = SUMSNOALB + FAREA * SNOALB(K)
          ENDDO

          FWAT = LANDUSEF(I,KWAT,J)
          !-- CHECK FOR WATER
          IF (FWAT .GE. 0.50) THEN        ! Changed WRFV3.7
            XLAI(I,J)    = LAIMX(KWAT)
            XLAIMN(I,J)  = LAIMN(KWAT)
            RSTMIN(I,J)  = RSMIN(KWAT)
            ZNT(I,J)     = Z0(KWAT)
            XVEG(I,J)    = VEG(KWAT)
            XVEGMN(I,J)  = VEGMN(KWAT)
            XSNUP(I,J)   = SNUP(KWAT)
            XALB(I,J)    = ALB(KWAT)
            XSNOALB(I,J) = SNOALB(KWAT)
          ELSE
            IF (FWAT .GT. 0.10) THEN
              ALAI   = ALAI - FWAT
              SUMLZ0 = SUMLZ0 - FWAT * ALOG(Z0(KWAT))
            ENDIF
            XLAI(I,J)    = SUMLAI / ALAI
            XLAIMN(I,J)  = SUMLMN / ALAI
            RSTMIN(I,J)  = SUMLAI / SUMRSI
            ZNT(I,J)     = EXP(SUMLZ0/ALAI)
            XVEG(I,J)    = SUMVEG / ALAI
            XVEGMN(I,J)  = SUMVMN / ALAI
            XSNUP(I,J)   = SUMSNUP
            XALB(I,J)    = SUMALB
            XSNOALB(I,J) = SUMSNOALB
          ENDIF

          IF (FWAT .GT. 0.50) THEN
            ZNT(I,J)     = Z0(KWAT)
            XALB(I,J)    = ALB(KWAT)
            XSNOALB(I,J) = SNOALB(KWAT)
          ENDIF

          !-- Compute wetlands fraction for proper MMLUIN data set
          !-- Note: if LU categories change, these hard coded indicies must be updated
          IF (LAND_USE_TYPE == 'USGS') THEN
            WETFRA(I,J) = LANDUSEF(I,17,J)+LANDUSEF(I,18,J)
          ELSE IF (LAND_USE_TYPE == 'USGS28') THEN
            WETFRA(I,J) = LANDUSEF(I,17,J)+LANDUSEF(I,18,J)
          ELSE IF (LAND_USE_TYPE == 'NLCD50') THEN
            WETFRA(I,J) = LANDUSEF(I,22,J)+LANDUSEF(I,23,J)+LANDUSEF(I,27,J)+LANDUSEF(I,28,J)+LANDUSEF(I,42,J)
          ELSE IF (LAND_USE_TYPE == 'NLCD40'.or.LAND_USE_TYPE == 'NLCD51') THEN
            WETFRA(I,J) = LANDUSEF(I,39,J)+LANDUSEF(I,40,J)+LANDUSEF(I,11,J)
          ELSE IF (LAND_USE_TYPE == 'MODIS') THEN
            WETFRA(I,J) = LANDUSEF(I,11,J)
          END IF

          ZNT(I,J)     = ZNT(I,J) * 0.01           !CONVERT TO M
          XVEG(I,J)    = XVEG(I,J) * 0.01          !CONVERT TO FRAC
          XVEGMN(I,J)  = XVEGMN(I,J) * 0.01
          XLAND(I,J)   = 1.0 + FWAT
          XALB(I,J)    = XALB(I,J) * 0.01
          XSNOALB(I,J) = XSNOALB(I,J) * 0.01
        
          !-------Adjustment according to CANFRA and IMPERV fo NLCD40 only -----------
          FIMP = IMPERV(I,J) * 0.01
          FCAN = CANFRA(I,J) * 0.01
          IF (LAND_USE_TYPE == 'NLCD40'.or.LAND_USE_TYPE == 'NLCD51') THEN
             XVEG(I,J)   = MIN(XVEG(I,J),1.0-FIMP)
             XVEGMN(I,J) = MIN(XVEGMN(I,J),1.0-FIMP)
             XVEG(I,J)   = MAX(XVEG(I,J),FCAN)
             XVEGMN(I,J) = MAX(XVEGMN(I,J),FCAN)
           
             FORFRA = LANDUSEF(I,39,J)+LANDUSEF(I,30,J)+LANDUSEF(I,29,J)+LANDUSEF(I,28,J)
             EXTFOR =  FCAN - FORFRA
             IF (EXTFOR.GE.0.01) THEN
                XLAI(I,J)   = LAIMX(30) * EXTFOR + XLAI(I,J) * (1 - EXTFOR)
                XLAIMN(I,J) = LAIMN(30) * EXTFOR + XLAIMN(I,J) * (1 - EXTFOR)
             ENDIF
          ENDIF
          !--------------------------------------------------------------------------
          IF(XLAND(I,J).GT.1.5) URBF2D(I,J) = -1.   ! set urban frac to -1 for water
        ENDDO     ! END LOOP THROUGH GRID CELLS
      ENDDO       ! END LOOP THROUGH GRID CELLS
      !--------------------------------------------------------------------

  END SUBROUTINE vegeland

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      SUBROUTINE SURFPX(ITIMESTEP,DTPBL, ISPXURB, IFLAND, ISNOW, NUDGEX, XICE1,                & !in
                        JULIAN, XTIME, GMT, SOLDN, GSW, LWDN, EMISSI, Z1,            & !in
                        MOL, ZNT, UST, ZA_CNP2, ZA_CNP1, ZF_CNP1, ZF_CNP0, DZF,      & !in UTK Urban
                        PSURF, DENS1, DENSRF,                                        & !in UTK Urban
                        QV1, QV_RF, QSS, QSS_RF, TA1, TA_RF,                         & !in UTK Urban
                        THETA1, THETA_RF, PRECIP, THETAV_RF, THETA_RF1, THETAV_RF1,  & !in UTK Urban
                        THETA_RFs, THETAV_RFs, THETA_WLs, CPAIR, CPAIR_RF, PSIH,     & !in UTK Urban
                        US_RF1, US_RF, VS_RF1, VS_RF,                           & !in UTK Urban
                        DENS_CAN, THETA_CAN, TA_CAN, WSP_CAN, CPAIR_CAN,        & !in UTK Urban
                        RH2OBS, T2OBS,                                          & !in
                        VEGFRC, ISTI, LAI, IMPERV, CANFRA,                      & !in
                        BETAP, RSTMIN, HC_SNOW, SNOW_FRA, WETFRA,               & !in
                        WWLT, WFC, WRES, CGSAT, WSAT, B,                        & !in
                        C1SAT, C2R, AS, JP, C3, DS1, DS2, QST12,                & !in
                        URBF, COSZ, CURR_SECS, ALB,                             & !-- DDY and UTK Urban in
                        XLATIJ, XLONIJ, SINAIJ, COSAIJ,                         & !in --Urban
                        RADNET, GRDFLX, HFX, QFX, LH,                           & !out
                        HFX_RF, HFX_WL, ANTHFX1, ANTHFX_RF, QFX_RF, LH_RF,      & !out UTK Urban
                        EG, ER, ETR, ER_ST, ER_RF,                              & !out UTK Urban
                        QST, CAPG, RS, RA, TG, T2, TG_RF, T2_RF, TG_WL, T2_WL,  &
                        RADWL, RADWL_FNL, RADST, RADST_FNL, GSW_URB, GSW_WL,    & !out UTK Urban
                        RADNET_URB, RADNET_RF, RADNET_WL, LWDN_URB, LWDN_RF,    & !out UTK Urban
                        LWDN_WL, LWUP_URB, LWUP_RF, LWUP_WL,                    & !out UTK Urban  
                        WG, W2, WR, WR_ST, WR_VEG, WR_RF,                       & !out
                        TA2, QA2, LAND_USE_TYPE, I, J,                          & !out
                        HGTBDG, LAMP, LAMF, STDIR1, STDIR2, BB1, BB2,           & !-- UTK Urban in
                        PPLN, VEGU,KSL, SWDDIR, SWDDIF,GAB_RF,GAB_WL            ) !-- UTK Urban in

!------------------------------------------------------------------------------
! 
!  FUNCTION:
!    THIS SUBROUTINE COMPUTES SOIL MOISTURE AND TEMPERATURE TENDENCIES
!    BY SOLVING THE PROGNOSTIC EQUATIONS IN PX95.
!
!  SUBROUTINES CALLED:
!    SUB. QFLUX   compute the soil and canopy evaporation, and transpiration
!    SUB. SMASS   compute nudging coefficients for soil moisture and temp nudging
!
!  ARGUMENTS:
!    DTPBL:   TIME STEP OF THE MINOR LOOP FOR THE LAND-SURFACE/PBL MODEL
!    IFLAND:  INDEX WHICH INDICATES THE TYPE OF SURFACE,=1,LAND;=2,SEA
!    ISNOW:   SNOW (=1) OR NOT (=0)
!    NUDGE:   SWITCH FOR SOIL MOISTURE NUDGING
!    SOLDN:   SHORT-WAVE RADIATION
!    LWDN:    LONG-WAVE RADIATION
!    EMISSI:  EMISSIVITY
!    Z1:      HEIGHT OF THE LOWEST HALF LAYER
!    MOL:     MONIN-OBUKOV LENGH (M)
!    ZNT:     ROUGHNESS LENGTH (M)
!    UST:     FRICTION VELOCITY (M/S)
!    TST:     Turbulent moisture scale
!    RA:      AERODYNAMIC RESISTENCE
!    PSURF:   P AT SURFACE (CB)
!    DENS1:   AIR DENSITY AT THE FIRST HALF LAYER
!    QV1:     Air humidity at first half layer
!    QSS:     Saturation mixing ratio at ground
!    TA1:     Air temperature at first half layer
!    THETA1:  Potential temperature at first half layer
!    PRECIP:  Precipitation rate in m/s
!    STBOLT:  STEFAN BOLTZMANN'S CONSTANT
!    KARMAN:  VON KARMAN CONSTANT
!    CPAIR:   Specific heat of moist air (M^2 S^-2 K^-1)
!    TAUINV:  1/1DAY(SEC)
!    VEGFRC:  Vegetation coverage
!    ISTI:    soil type
!    LAI:     Leaf area index
!    IMPERV:  Percentage of IMPERVIOUS Fraction 
!    CANFRA:  Percentage of Canopy/Tree Fraction 
!    BETAP:   Coefficient for bare soil evaporation
!    WETFRA:  Fraction of Wetlands area 
!    THZ1OB:  Observed TEMP FROM ANAL OF OBS TEMPS AT SCREEN HT
!    RHOBS:   Observed relative humidity at SCREEN HT
!    RSTMIN   Minimum Stomatol resistence
!... Outputs from SURFPX
!    RADNET: Net radiation
!    HFX:    SENSIBLE HEAT FLUX (W / M^2)
!    QFX:    TOTAL EVAP FLUX (KG / M^2 S)
!    GRDFLX: Ground heat flux (W / M^2)
!    QST:    Turbulent moisture scale
!    CAPG:   THERMAL CAPACITY OF GROUND SLAB (J/M^2/K)
!    RS:     Surface resistence
!    RA:     Surface aerodynamic resistence
!    EG:     evaporation from ground (bare soil)
!    ER:     evaporation from canopy
!    ETR:    transpiration from vegetation
!    TA2:    diagnostic 2-m temperature from surface layer and lsm
!
!... Updated variables in this subroutine
!    TG:     Soil temperature at first soil layer
!    T2:     Soil temperature in root zone
!    WG:     Soil moisture at first soil layer
!    W2:     Soil moisture at root zone
!    WR:     Canopy water content in m
!
!  REVISION HISTORY:
!     AX     02/2005 - developed WRF version based on the MM5 PX LSM
!   
!     UTK    05/2018 - Estimation of Urban ground+rooftop surface temperatures and heat fluxes.
!                    - Estimation of evaporation rate from the urban surfaces and corresponding 
!                      moisture fluxes during the precipitation events.
!                    - Estimation of anthropogenic heat fluxes from the buildings based upon 
!                      the heating degree days and cooling degree days approach. 
!
!     UTK    03/2021 - Incorporation of the complete urban solar radiation model based upon Masson (2000).
!                    - Estimation of the urban wall surface temperatures and heat fluxes.     
!
!
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      IMPLICIT NONE

!.......Arguments

!.. Integer
      INTEGER, INTENT(IN) :: ISTI, NUDGEX, I, J,ITIMESTEP

!... Real
      REAL, INTENT(IN)    :: DTPBL, DS1, DS2
      REAL, INTENT(IN)    :: IFLAND, ISNOW, XICE1
      REAL, INTENT(IN)    :: JULIAN, XTIME, GMT
      REAL, INTENT(IN)    :: SOLDN, GSW, LWDN, EMISSI, Z1
      REAL, INTENT(IN)    :: ZNT
      REAL, INTENT(IN)    :: PSURF, DENS1, QV1, QSS, TA1, THETA1, PRECIP
      REAL, INTENT(IN)    :: CPAIR, CPAIR_RF
      REAL, INTENT(IN)    :: VEGFRC, LAI, IMPERV, CANFRA
      REAL, INTENT(IN)    :: RSTMIN, HC_SNOW, SNOW_FRA, WETFRA 
      REAL, INTENT(IN)    :: WWLT, WFC, WRES, CGSAT, WSAT, B, C1SAT, C2R, AS, JP, C3
      REAL, INTENT(IN)    :: RH2OBS, T2OBS
      REAL, INTENT(IN)    :: QST12

      REAL, INTENT(OUT)   :: EG, ER, ETR
      REAL, INTENT(OUT)   :: QST, CAPG, RS, TA2, QA2

      REAL, INTENT(INOUT) :: WG, W2, WR, UST, RA, BETAP , TG
      DOUBLE PRECISION, INTENT(INOUT) :: T2
      REAL, INTENT(INOUT) :: GRDFLX, QFX, HFX, LH, PSIH, MOL

      CHARACTER (LEN = 6), INTENT(IN) :: LAND_USE_TYPE

      !-- UTK Urban
      REAL, INTENT(IN)    :: QV_RF, QSS_RF, DENSRF, THETA_RF
      REAL, INTENT(IN)    :: ZA_CNP2, ZA_CNP1, DZF, TA_RF
      REAL, INTENT(IN)    :: THETAV_RF, THETA_RF1, THETAV_RF1
      REAL, INTENT(IN)    :: ZF_CNP1, ZF_CNP0, THETA_RFs, THETAV_RFs, THETA_WLs
      REAL, INTENT(IN)    :: US_RF1, US_RF, VS_RF1, VS_RF
      REAL, INTENT(IN)    :: DENS_CAN, THETA_CAN, TA_CAN, WSP_CAN, CPAIR_CAN
      REAL, INTENT(OUT)   :: ER_ST, ER_RF
      REAL, INTENT(INOUT) :: WR_ST, WR_VEG, WR_RF
      REAL, INTENT(INOUT) :: HFX_RF, HFX_WL, ANTHFX1, ANTHFX_RF, QFX_RF, LH_RF
      REAL, INTENT(INOUT) :: TG_RF, TG_WL
      DOUBLE PRECISION, INTENT(INOUT) :: T2_RF, T2_WL
      REAL, INTENT(INOUT) :: RADWL, RADWL_FNL, RADST, RADST_FNL, GSW_URB, GSW_WL
      REAL, INTENT(INOUT) :: RADNET, RADNET_URB, RADNET_RF
      REAL, INTENT(INOUT) :: LWDN_URB, LWDN_RF, LWUP_URB, LWUP_RF
      REAL, INTENT(INOUT) :: RADNET_WL, LWDN_WL, LWUP_WL 
      REAL, INTENT(OUT)   :: GAB_RF,GAB_WL

      LOGICAL, INTENT(IN) :: ISPXURB


!... Local Variables

!... Real
      REAL :: HF, LV, CQ4, WETSAT, SM2, LV_RF 
      REAL :: RAH, RAW, ET, ET_RF, W2CG, CG, CT, CT_RF, CT_WL, SOILFLX
      REAL :: CPOT, THETAG, THETAG_RF
      REAL :: ZOL, ZOBOL, ZNTOL, Y, Y0, PSIH15, YNT
      REAL :: WGNUDG, W2NUDG, ALPH1, ALPH2, BET1, BET2
      REAL :: CQ1, CQ2, CQ3, COEFFNP1, COEFFN, TSNEW, TSHLF
      REAL :: CQ1_G, CQ2_G, COEFFNP1_G, COEFFN_G, TSNEW_G, TSHLF_G
      REAL :: CQ1_RF, CQ2_RF, COEFFNP1_RF, COEFFN_RF, TSNEW_RF, TSHLF_RF 
      DOUBLE PRECISION ::T2NEW_RF, T2NEW_WL, T2NEW_G, T2NEW, T2U
      REAL :: COEFFNP1_WL, COEFFN_WL, TSNEW_WL, TSHLF_WL
      REAL :: ROFF, WRMAX, PC, DWR, PNET, TENDWR, WRNEW
      REAL :: COF1, CFNP1WR, CFNWR, PG, FASS
      REAL :: TENDW2, W2NEW, W2HLF, W2REL, C1, C2, WEQ, CFNP1, CFN, WGNEW
      REAL :: ALN10, TST
      REAL :: QST1, PHIH, PSIOB
      REAL :: T2NUD, T2NUDF
      REAL :: VAPPRS, QSBT, RH2MOD, IMF, VEGF, SOILF
      REAL :: RSOIL, LDRY, DP
      REAL :: C1MAX, ZZA, ZZB, ZDEL, ZLY, ZA, ZB, ZY2
      REAL :: Rinc, Hcan
      REAL :: CQ3BG, CQ3VW, CQ3VG, SIGG

      !-- DDY and UTK Urban
      REAL, INTENT(IN) :: URBF, COSZ, XLATIJ, XLONIJ, SINAIJ, COSAIJ, CURR_SECS, ALB
      REAL :: SXLONG, DA, XT24, TLOCTM, ARGU, CSZA
      REAL :: WSTT, BBDG, TANZ
      REAL :: AZMA    ! Solar azimuth angle 
      REAL :: HRA     ! Hour angle in degrees
      REAL :: DECLN   ! Declination angle
      REAL :: RFFA    ! Roof fraction area
      REAL :: WLFA    ! Wall fraction area
      REAL :: STVEGFA ! Street + Vegetation fraction area 
      REAL :: STFA    ! Street fraction area
      REAL :: VEGFA   ! Vegetation fraction area
      REAL :: EOT     ! Equation of time in minutes
      REAL :: RESWL   ! Aerodynamic resistance for walls
      REAL :: ETAC    ! Eta. coefficient for walls

      REAL :: WRMAX_VEG, WRMAX_ST, WRNEW_VEG, WRNEW_ST
      REAL :: WRNEW_RF, WRMAX_RF, CQ4_RF, CQ4_ST
      REAL :: WSP1, AA, BB, WSPD, RIN, FH2
      REAL :: RZ0RT, CF, FHH2, RBPW
      REAL :: MR, MW                      ! Sum of the reflections against the road and wall
      REAL :: RR, RW                      ! Reflected radiation parts from the road and wall
      REAL :: CON1, CON2, HNX, CNX, GFX   ! For Anthro flux estimation  
      REAL :: STANGX
      INTEGER :: IDIR   
      REAL :: FWALL, RBPW1, RBPW2 
      REAL :: time
      integer ihr
      
      !-- DDY, UTK Urban morphological parameters
      REAL, INTENT(IN) :: LAMF            ! Frontal area density
      REAL, INTENT(IN) :: LAMP            ! Plan area density
      REAL, INTENT(IN) :: HGTBDG          ! Average building height
      INTEGER, INTENT(IN) :: STDIR1       ! Dominant street angle (canyon orientation)
      INTEGER, INTENT(IN) :: STDIR2       ! 2nd Dominant street angle (canyon orientation)
      REAL, INTENT(INOUT) :: BB1             ! Block width between streets in STDIR1 direction
      REAL, INTENT(INOUT) :: BB2             ! Block width between streets in STDIR2 direction
      REAL, INTENT(IN) :: PPLN            ! Population density  
      REAL, INTENT(IN) :: VEGU            ! veg frac in Urban part
      REAL :: SKVFRD                      ! Sky view factor for canyon road
      REAL :: SKVFWL                      ! Sky view factor for canyon one wall
      REAL, PARAMETER :: RDALB = 0.2      ! Albedo for the concrete road surface
      REAL, PARAMETER :: WLALB = 0.4      ! Albedo for the concrete building wall surface
      REAL, PARAMETER :: SNALB = 0.8      ! Albedo for the snow covered surface
      REAL, PARAMETER :: RDEMISSI = 0.9   ! Emissivity for the concrete road surface 
      REAL, PARAMETER :: WLEMISSI = 0.9   ! Emissivity for the concrete building wall surface
      REAL, PARAMETER :: RFEMISSI = 0.9   ! Emissivity for the concrete building roof surface
      REAL, PARAMETER :: SNEMISSI = 0.99  ! Emissivity for the snow surface      
      INTEGER, INTENT(IN) :: KSL          ! Urban canopy level K-index 
      REAL, INTENT(IN) :: SWDDIR          ! Direct SW radiation from goddard scheme
      REAL, INTENT(IN) :: SWDDIF          ! Diffused SW radiation from goddard scheme


!... Parameters
      REAL :: ZOBS, GAMAH, BETAH, SIGF, BH, CT_SNOW, CT_IMPERV
      REAL :: RDEMISSIBR                  ! Emissivity for the road surface together with snow
      REAL :: RFEMISSIBR                  ! Emissivity for the rooftop surface together with snow
      REAL :: RDALBBR                     ! Albedo for for the road surface together with snow
      REAL :: TGU
      REAL :: SWDIFRD, SWDIFWL, CHWL

      REAL, PARAMETER :: CV = 1.2E-5   ! K-M2/J Note: Update from 8E-6 10/14 Jon Pleim

      PARAMETER (ZOBS  = 1.5)    ! height for observed screen temp., (m)
      PARAMETER (BH    = 15.7)
      PARAMETER (GAMAH = 16.)    !11.6)
      PARAMETER (BETAH = 5.0)    !8.21)
      PARAMETER (SIGF  = 0.5)    ! rain interception see LSM (can be 0-1)
      REAL, PARAMETER  :: DWAT   = 0.2178  ! [cm^2 / s] at 273.15K

      !--------------------------------------------------------------------
      !  OLD PX legacy value from MM5 ... unknown origin PARAMETER (CT_SNOW  = 5.54E-5)   
      ! New value of CT_SNOW calibrated using multilayer soil model where csnow=6.9E5 J/(m3 K) 
      ! from NCAR (WRFV3.2 -WRFV3.6.1) CSM PARAMETER (CT_SNOW  = 2.0E-5) 
      PARAMETER (CT_SNOW  = 2.0E-5)  

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! CS for concrete/asphalt/road material (http://www.engineeringtoolbox.com/specific-heat-capacity-d_391.html)
      !                                        cs_imperv = 920 J kg-1 K-1 
      ! CS for asphalt and concrete from 0.1785 to 0.20 cal g-1 K-1 (http://pages.towson.edu/morgan/files/Impervious.pdf)
      ! the values above translate to ~750 to 850 J kg-1 K-1
      !
      ! CAPG used for WRF urban physics ranges from 1.0E6 for roof and building walls to 1.4E6 J m-3 K-1 for roads/urban ground
      ! Using these values to back out CS along with 0.15 m thickness 1.4E6 J m-3 K-1 * 0.15 m = 2.1E5 J m-2 K-1
      ! inverse of the value above gives CT_IMPERV value of 1/0.000021 = 4.762E-6 K m2 J-1

      ! The middle value will be used for now. 850 J kg-1 K-1. This needs to be converted from J/K per kg to area using 
      ! the approxiate concrete/asphalt density and layer thickness or represenative thickness. For starters (12/2011)
      ! well not use the PX first layer thickness, but representative thickness of most roads/parkinglots/buildings.
      ! for now welll use 6 inches or about 15 cm or 0.15 m. Density of concrete (normal) from 
      ! Dorf, Richard. Engineering Handbook. New York: CRC Press, 1996. is ~2400 kg m-3.
      ! Using these values 850 J kg-1 K-1 * 2400 kg m-3 * 0.15 m = 3.06E5 J m-2 K-1 or in CT form (inverse) 3.268E-6 K m2 J-1
      ! If you look at the range of possible values considering density differences of concrete/asphalt/etc
      ! Values can range from 2.5 to 6.0 E-6                                 
      PARAMETER (CT_IMPERV  = 3.268E-6)  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !---UTK Urban
      REAL, PARAMETER :: CT_STREET = 3.268E-6            !--UTK Urban (******CHANGE IT******)
      REAL, PARAMETER :: SIGFF     = 0.7                 ! Rain interception for the street
      REAL, PARAMETER :: CT_ROOF   = 8.6E-6   ! 3.268E-6 
      REAL, PARAMETER :: CT_WALL   = 6.7E-6   ! 3.268E-6
      REAL, PARAMETER :: B1        = 9.4
      REAL, PARAMETER :: B2        = B1/2.
      REAL, PARAMETER :: CSTH      = 5.3
      REAL, PARAMETER :: Z0RT      = 0.05                ! Building roof roughness length (m)
      REAL, PARAMETER :: KARMAN    = 0.4 
      REAL, PARAMETER :: PPLNNYC   = 10430.   ! NYC pop den Sailor et al 2015
      REAL, PARAMETER :: Qfmax     = 63.      ! NYC
      REAL, PARAMETER :: Fveh      = 0.39
      REAL, DIMENSION(24) :: ProHr    ! Summer
      DATA PROHR/0.25,0.23,0.25,0.21,0.22,0.29,0.53,0.82,0.87,0.80,0.80,0.84,    &
                 0.89,0.89,0.93,1.00,0.90,0.78,0.56,0.48,0.44,0.41,0.36,0.30/ 
      REAL ProFac,LHR
      INTEGER LastHr, NextHr
      REAL BW1,BW2

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

        ALN10  = ALOG(10.0)
        RADNET = GSW - (EMISSI *(STBOLT *TG **4 - LWDN)) ! NET RADIATION (Only a diagnostic term and not used in model further)

        CPOT   = (100.0 / PSURF) ** ROVCP                ! rcp is global constant(module_model_constants)
        THETAG = TG * CPOT

        IF (ISPXURB .AND. (URBF .GT. 0.0)) THEN
          THETAG_RF = TG_RF * CPOT
        ENDIF

        ZOL   = Z1 / MOL                                                       
        ZOBOL = ZOBS / MOL                                                      
        ZNTOL = ZNT / MOL                              
        RZ0RT = 1.0 / Z0RT


        IF (ISPXURB .AND. (URBF .GT. 0.0)) THEN
          !-- Building (Urban Morphology) Parameters 
!         BBDG = 60. !(LAMP / LAMF) * HGTBDG                 ! Indv. Building Width
!          WSTT = (1.0 / SQRT(LAMP) - 1.0) * BBDG        ! Indv. Street Width
!          WSTT = MIN(WSTT, 40.)                         ! max street width
!          WSTT = MAX(WSTT, 10.)                         ! min street width
          WSTT = 20.0
          BW1 = max(BB1-WSTT,10.)
          BW2 = max(BB2-WSTT,10.)
!          RBPW = 1.0 / (BBDG + WSTT)     
          RBPW1   = 1.0 / (BW1 + WSTT)
          RBPW2   = 1.0 / (BW2 + WSTT)    
!         RBPW = 1.0 / (BBDG + WSTT)
        ENDIF


        !-----------------------------------------------------------------------------------------
        IF (MOL .LT. 0.0) THEN
          Y      = ( 1.0 - GAMAH * ZOL ) ** 0.5
          Y0     = ( 1.0 - GAMAH * ZOBOL ) ** 0.5
          YNT    = ( 1.0 - GAMAH * ZNTOL ) ** 0.5
          PSIH15 =  2.0 * ALOG((Y + 1.0) / (Y0 + 1.0))
          PSIH   =  2.0 * ALOG((Y + 1.0) / (YNT + 1.0))
          PSIOB  =  2.0 * ALOG((Y0 + 1.0) / (YNT + 1.0))
          PHIH = 1.0 / Y
        ELSE
          IF((ZOL-ZNTOL).LE.1.0) THEN                                         
            PSIH = -BETAH*(ZOL-ZNTOL)                                        
          ELSE                                                               
            PSIH = 1.-BETAH-(ZOL-ZNTOL)                                      
          ENDIF             
          IF ((ZOBOL - ZNTOL) .LE. 1.0) THEN
            PSIOB = -BETAH * (ZOBOL - ZNTOL)
          ELSE
            PSIOB = 1.0 - BETAH - (ZOBOL - ZNTOL)
          ENDIF
          PSIH15 =  PSIH - PSIOB
          IF(ZOL.LE.1.0) THEN
            PHIH = 1.0 + BETAH * ZOL
          ELSE
            PHIH = BETAH + ZOL
          ENDIF
        ENDIF

        !-----------------------------------------------------------------------------------------
        !-- ADD RA AND RB FOR HEAT AND MOISTURE
        !... RB FOR HEAT = 5 /UST
        !... RB FOR WATER VAPOR =  5*(0.599/0.709)^2/3 /UST = 4.503/UST
        RA  = PR0 * ( ALOG(Z1 / ZNT) - PSIH ) / (KARMAN*UST)
        RAH = RA + 5.0 / UST
        RAW = RA + 4.503 / UST
        IF (IFLAND .LT. 1.5.AND. XICE1.LT.0.5) THEN
          LDRY  = 1.75 * DS1 * (EXP((1. - WG/WSAT)**5) - 1.)/1.718  ! 1.75 cm is the layer thickness used by S&Z09
          DP    = DWAT * 1.E-4 * WSAT**2 * (1. - WRES/WSAT)**(2. + 3./B)
          RSOIL = LDRY / DP
        ELSE
          RSOIL = 0.0
        ENDIF

        !--------------------------------------------------------------------
        ! Compute soil moisture layer 2 that considers fraction of saturated 
        ! wetlands. If 100% of cell is wetland, soil moisture can be no lower
        ! than full soil saturation. If half wetland, no less than half saturated 
        IF (IFLAND .LT. 1.5) THEN
          WETSAT = 1.00 * WSAT                           ! Wetlands soil moisture
          SM2    = (WETFRA * WETSAT)                     ! 
          W2     = AMAX1(SM2, W2)                        ! In case that W2 > Field capacity (heavy precip), use wetter W2
        ENDIF

        !--------------------------------------------------------------------------
        !--Computation of Roof fraction, street fraction and Veg fraction area for
        !---an urban grid cell
        !-------------------------------------------
        ! The street fraction area is calculated by subtracting urban plan area index from 1.0.
        ! And consequently, the vegetation fraction area and soil fraction area are obtained too.
        IF (ISPXURB .AND. (URBF .GT. 0.0)) THEN  
          RFFA    = LAMP
!          STVEGFA = 1.0 - LAMP
!          VEGFA   = VEGFRC
!          STFA    = 1.0 - VEGFA
          VEGFA   = VEGU
          STFA    = MAX(1.0 - LAMP - VEGU,0.0)
!          WLFA    = LAMF          
        ENDIF

        !--------------------------------------------------------------------
        ! In-canopy resistance  - Erisman et al (1994)
        Hcan  = ZNT * 10.0
        Rinc  = 14.0 * LAI* Hcan / UST

        !--------------------------------------------------------------------
        !--  COMPUTE MOISTURE FLUX
        CALL QFLUX(URBF, DTPBL, ISPXURB, DENS1, DENSRF,                   &
                   QV1, QV_RF, TA1, SOLDN, RAW, QSS, QSS_RF,              &
                   VEGFRC, ISNOW,  ISTI, IFLAND, LAI, BETAP,              &
                   WG, W2, WR, WR_ST, WR_VEG, WR_RF, STFA, VEGFA, RFFA,   &
                   RSTMIN, WWLT, WFC, RSOIL, RINC, SIGG,                  &
                   EG, ER, ETR, ER_ST, ER_RF,                             &
                   CQ4, CQ4_RF, CQ4_ST, RS, FASS                          )
        !--------------------------------------------------------------------

        !--------------------------------------------------------------------
        !    Compute Total evaporation (ET) from various modes of moisture flux
        !-----From Urban-veg grid in D4
        IF (ISPXURB .AND. (URBF .GT. 0.0)) THEN
          ET    = EG*(1-IMPERV*0.01) + ER_ST + ER + ETR
          ET_RF = ER_RF
        ELSE  
          !-----(Default PX)
          ET  = EG*(1-IMPERV*0.01) + ER + ETR
        ENDIF

        QST = -ET / (DENS1 * UST)
 
        LV  = 2.83E6                        !-- LATENT HEAT OF SUBLIMATION AT 0C FROM STULL(1988)
        IF ((ISNOW .LT. 0.5) .AND. (TG .GT. 273.15))                            &
                    LV = (2.501 - 0.00237 * (TG - 273.15)) * 1.E6  !-- FROM STULL(1988) in J/KG
        
        QFX = ET
        LH  = LV * QFX                      !-- latent heat flux


        IF (ISPXURB .AND. (URBF .GT. 0.0)) THEN
          LV_RF  = (2.501 - 0.00237 * (TG_RF - 273.15)) * 1.E6         !-- FROM STULL(1988) in J/KG
          QFX_RF = ET_RF                                 
          LH_RF  = LV_RF * QFX_RF           !-- latent heat flux at roof  
        ENDIF
        !-----------------------------------------------------------------------------------------


        !-----------------------------------------------------------------------------------------
        ! Anthropogenic Heat Flux Estimation from Urban Grid Cells ---------------
        !----------------------------------------------------------------------------------------
        IF (ISPXURB .AND. (URBF .GT. 0.0)) THEN
        ! This HRA calculation is based on subroutine calc_coszen in module_radiation_driver.F
           DA     = 6.2831853071795862 * (JULIAN-1)/365.
           EOT    = (0.000075 + 0.001868 * COS(DA) - 0.032077 * SIN(DA) &
                    -0.014615 * COS(2*DA) - 0.04089 * SIN(2*DA)) * (229.18)
           XT24   = MOD(XTIME,1440.) + EOT
           TLOCTM = GMT + XT24/60. + XLONIJ/15.
           HRA    = 15. * (TLOCTM - 12.)   ! Hour angle in degrees

        !-----------------------------------------------------------------------------------------
        ! Anthropogenic Heat Flux Estimation from Urban Grid Cells ---------------
        !  References: Sailor et al. 2015
        !----------------------------------------------------------------------------------------
          ANTHFX1   = 0.0
          ANTHFX_RF = 0.0
          if(TLOCTM.lt.1.0) then
            LHR = TLOCTM + 24.
          else
            LHR = TLOCTM
          endif
          LastHr = int(LHR) 
          NextHr = LastHr + 1
          IF(LHR.LT.24.) Then
            ProFac = ProHr(LastHr) + (LHR-float(LastHr)) * (ProHr(NextHr)-ProHr(LastHr))
          ELSE 
            ProFac = ProHr(24) + (LHR-24.) * (ProHr(1)-ProHr(24))
          ENDIF
          ANTHFX1   =   0.5* (PPLN / PPLNNYC * (1.-Fveh) +  Fveh) * ProFac * Qfmax
          ANTHFX_RF = ANTHFX1
!          if(HGTBDG.gt.120.0) then
!             print *,' LHR,ProFac,PPLN,ANTHFX1=',LHR,ProFac,PPLN,ANTHFX1
!          endif
        ENDIF


        !-----------------------------------------------------------------------------------------
        ! Surface Sensible Heat Flux  ( for Urban/Non-Urban ) ------------------------
        !-----------------------------------------------------------------------------------------
        TST = (THETA1 - THETAG ) / (UST*RAH)
        HF  = UST * TST           

        GFX = 0.0
        HFX = 0.0 
        IF (ISPXURB .AND. (URBF .GT. 0.0)) THEN       ! Urban Ground Heat Flux Estimation 
           GFX = AMAX1(-DENS1 * CPAIR * HF, -250.0)  
           HFX = GFX + ANTHFX1 * URBF   ! using -250. from MM5                   !-- UTK Urban 
        ELSE    ! Non-Urban Ground Heat Flux           
           HFX = AMAX1(-DENS1 * CPAIR * HF, -250.0)   ! using -250. from MM5
        ENDIF


        !-----------------------------------------------------------------------------------------
        ! Urban Rooftop Heat Flux Estimation ----------------------------
        !-----------------------------------------------------------------------------------------
        HFX_RF = 0.0 
        GAB_RF = 0.0
        IF (ISPXURB .AND. (URBF .GT. 0.0)) THEN  
           !-- Richardson Number Estimation 
           WSP1 = MAX(SQRT(US_RF**2 + VS_RF**2), 0.01)
           AA   = THETAV_RF1 - THETAV_RF
!           AA   = THETAV_RF - THETAV_RFs
           BB   = 2 * 9.81 / (THETAV_RF1 + THETAV_RF) 
!           BB   = 9.81 / THETAV_RFs

           WSPD = (US_RF1 - US_RF)**2 + (VS_RF1 - VS_RF)**2      
           WSPD = MAX(WSPD,0.1)
           RIN  = ABS(ZA_CNP2 - ZA_CNP1) * AA * BB / WSPD

           !-- Determination of Stability conditions for (Louis,1979) functions FH 
           IF (RIN .GE. 0.0) THEN
             !- Stable condition
             AA  = 1.0 + B2 * RIN
             AA  = AA * AA
             FH2 = 1.0 / AA
           ELSE
             !- Unstable condition
             AA  = DZF * 0.5 * RZ0RT
             BB  = KARMAN / LOG(AA)
             CF  = BB * BB * B1 * SQRT(AA)
             AA  = B1 * RIN
             BB  = CF * SQRT(-RIN)
             FH2 = 1.0 - AA / (1.0 + CSTH * BB)
           ENDIF

           !-- Urban Rooftop heat flux estimation (Martilli et al.,2002)
           AA   = KARMAN / LOG(DZF * 0.5 * RZ0RT)
           AA   = AA * AA * WSP1
           FHH2 = AA * FH2 / 0.74 
           GAB_RF = FHH2 * URBF * RFFA   ! m/s
           HFX_RF = - DENSRF * CPAIR_RF * GAB_RF * (TA_RF - TG_RF) + ANTHFX_RF *URBF   !-- UTK Urban 
        ENDIF  


        !-----------------------------------------------------------------------------------------
        ! Urban Wall Heat Flux Estimation ----------------------------
        !-----------------------------------------------------------------------------------------
        RESWL  = 0.0
        ETAC   = 0.0
        HFX_WL = 0.0
        GAB_WL = 0.0
        IF (ISPXURB .AND. (URBF .GT. 0.0)) THEN
           ! Forced convection coefficient for walls [Wm^-2/deg.K] (Clarke,1985)
           ETAC   = 5.678*(1.09 + 0.23*(WSP_CAN/0.3048))
!           CHWL   = WSP_CAN * DENS_CAN*CPAIR_CAN/ 35.   ! Ra+Rb = 35./U where Zo=0.1
           ! Urban Walls heat flux estimation (Martilli et al.,2002)
           Fwall = 2. * HGTBDG * ((RBPW2+RBPW1) - WSTT*RBPW2*RBPW1)/(1-lamp)
           GAB_WL = ETAC * URBF * Fwall 
           HFX_WL = - GAB_WL * (THETA_CAN - THETA_WLs)                  !-- UTK Urban
           GAB_WL = GAB_WL / (DENS_CAN* CPAIR_CAN) 
            time = itimestep*dtpbl  !s
            ihr = mod(int(time),3600)
!           if(HGTBDG.gt.120.0.and.ihr.eq.0) then
!                print *,'FWALL,URBF,RFFA=',FWALL,URBF,RFFA
!                print *,'HGTBDG,time,GAB_WL,GAB_RF=',HGTBDG,time,GAB_WL,GAB_RF,wsp_can
!           endif
!           HFX_WL = - CHWL * (THETA_CAN - THETA_WLs) * URBF * (2*HGTBDG/WSTT) 
!           IF(HGTBDG.GT.40.0.AND.SWDDIR.GT.200.)print*,' ETAC,CHWL,FHH2=',ETAC,CHWL,FHH2
        ENDIF 
        !-----------------------------------------------------------------------------------------


        !-----------------------------------------------------------------------------------------
        ! Compute the diagnosed 2m Q and T consistent with PX LSM
        QST1 = 0.5*(QST+QST12/PHIH)
        TA2  = (THETAG + TST * (PR0 / KARMAN * (ALOG(ZOBS / ZNT) - PSIOB)+5.))/CPOT
        QA2  = QV1 - QST1 * PR0/ KARMAN *  (ALOG(Z1 / ZOBS) - PSIH15)

        IF (QA2.LE.0.0) QA2 = QV1

        !--  Relative humidity
        VAPPRS = SVP1 * EXP(SVP2 * (TA2 - SVPT0) / (TA2 - SVP3))
        QSBT   = EP_2 * VAPPRS / (PSURF - VAPPRS)
        RH2MOD = QA2 / QSBT
        !-----------------------------------------------------------------------------------------

        !-----------------------------------------------------------------------------------------
        !-----ESTIMATION OF CT, CAPG, GRDFLX for Urban grid cell-----------   !---UTK Urban
        CT    = 0.0
        CT_RF = 0.0
        CT_WL = 0.0
        IF (ISPXURB .AND. (URBF .GT. 0.0)) THEN
!          CT    = 1./(STFA/CT_STREET + VEGFA/CV)
!          CT    = 1./(SNOW_FRA/CT_SNOW + (1-SNOW_FRA)/CT)
!          CAPG  = 1.0/CT 

          CT_RF = CT_ROOF
          CT_RF = 1./(SNOW_FRA/CT_SNOW + (1-SNOW_FRA)/CT_RF)
 
          CT_WL = CT_WALL

          RDEMISSIBR = (1.0-SNOW_FRA)*RDEMISSI + SNOW_FRA*SNEMISSI
          RFEMISSIBR = (1.0-SNOW_FRA)*RFEMISSI + SNOW_FRA*SNEMISSI  

          RDALBBR    = (1.0-SNOW_FRA)*RDALB + SNOW_FRA*SNALB          

          SOILFLX    = 2.0 * PI * TAUINV * (TG - T2)     !-- It's just an output, not used  
!          GRDFLX     = SOILFLX / CT                      !-- It's just an output, not used 
        ENDIF

        IF (IFLAND .LT. 1.5) THEN
          W2CG = AMAX1(W2,WWLT)
          CG   = CGSAT * 1.0E-6 * (WSAT/ W2CG) **    &  
                 (0.5 * B / ALN10)
          ! IMPERVIOUS weighting scheme -- Subtract highly accurate impervious fraction from cell
          ! remainder is split between ground and vegetation. CT is a weighted fractional average.
          ! Snow CT is then applied for final heat capacity
          IMF   = AMAX1(0.0, IMPERV/100.0)
          VEGF  = (1.0 - IMF) * VEGFRC
          SOILF = (1.0 - IMF) * (1.0 - VEGFRC)
          CT    = 1./( IMF/CT_IMPERV + VEGF/CV + SOILF/CG)                        
          CT    = 1./(SNOW_FRA/CT_SNOW + (1 - SNOW_FRA)/CT)
          CAPG  = 1.0/CT          

          SOILFLX = 2.0 * PI * TAUINV * (TG - T2)     !-- It's just an output, not used
          GRDFLX  = SOILFLX / CT                      !-- It's just an output, not used
        ENDIF

        !-----------------------------------------------------------------------------------------

        !--------------------------------------------------------------------
        !-- ASSIMILATION --- COMPUTE SOIL MOISTURE NUDGING FROM TA2 and RH2               
        !-------COMPUTE ASSIMILATION COEFFICIENTS FOR ALL I                            

        T2NUD = 0.0
        IF (IFLAND .LT. 1.5) THEN
          IF (NUDGEX .EQ. 0) THEN                                          !-- NO NUDGING CASE                
            WGNUDG = 0.0                                                        
            W2NUDG = 0.0    
            T2NUD  = 0.0                                    
          ELSE                                                               !-- NUDGING CASE        
            
            CALL SMASS (ISTI,  FASS,  SOLDN,   VEGFRC, RA, WWLT, WFC,   &                       
                        ALPH1, ALPH2, BET1, BET2, T2NUDF)                             

            !--COMPUTE SOIL MOISTURE NUDGING
            WGNUDG = ALPH1 * (T2OBS - TA2) + ALPH2 * (RH2OBS - RH2MOD) * 100  !NUDGING W2 FOR NON-VEG
            W2NUDG = BET1  * (T2OBS - TA2) + BET2  * (RH2OBS - RH2MOD) * 100
            IF (W2 .GE. WFC)  W2NUDG = AMIN1(W2NUDG,0.0)
            IF (W2 .LE. WWLT) W2NUDG = AMAX1(W2NUDG,0.0)
            IF (W2 .GE. WFC)  WGNUDG = AMIN1(WGNUDG,0.0)
            IF (W2 .LE. WWLT) WGNUDG = AMAX1(WGNUDG,0.0)
            T2NUD = T2NUDF * (T2OBS - TA2)
          ENDIF
        ENDIF
        !-----------------------------------------------------------------------------------------


        !------------------
        !-- Initialization
        !- For Urban Roof
        CQ1_RF      = 0.0
        CQ2_RF      = 0.0
        COEFFNP1_RF = 0.0
        COEFFN_RF   = 0.0
        TSNEW_RF    = 0.0
        TSHLF_RF    = 0.0 
        T2NEW_RF    = 0.0d0

        !- For Urban Ground
        CQ1_G       = 0.0
        CQ2_G       = 0.0
        COEFFNP1_G  = 0.0
        COEFFN_G    = 0.0
        TSNEW_G     = 0.0
        TSHLF_G     = 0.0
        T2NEW_G     = 0.0d0

        !- For Urban Walls
        COEFFNP1_WL = 0.0
        COEFFN_WL   = 0.0
        TSNEW_WL    = 0.0
        TSHLF_WL    = 0.0
        T2NEW_WL    = 0.0d0

        RADST       = 0.0
        RADST_FNL   = 0.0
        GSW_URB     = 0.0            
        RADWL       = 0.0
        RADWL_FNL   = 0.0  
        GSW_WL      = 0.0

        !- Default PX Land
        CQ1         = 0.0
        CQ2         = 0.0
        CQ3         = 0.0
        COEFFNP1    = 0.0
        COEFFN      = 0.0
        TSNEW       = 0.0
        TSHLF       = 0.0
        T2NEW       = 0.0d0

        !-----------------------------------------------------------------------------------------
        !-- Compute new values for TS,T2,WG,W2 and WR. No change over ice or water (XLAND > 1)
        !-------ESTIMATION of Urban Surface Temperatures----------    !----UTK Urban
        IF (ISPXURB .AND. (URBF .GT. 0.0)) THEN     ! For Urban Land Surface
           !-- Computation of Urban Roof Temperature, Solving by CRANK-NIC Implicit method:
           !-- TENDTS=CT*(RADNET-HFX-QFX)
           !-- Calculate the coefficients for implicit calculation of ROOF TEMP.
           CQ1_RF      = (1.0 - 0.622 * LV_RF * CRANKP / (r_d * TG_RF)) * QSS_RF
           CQ2_RF      = 0.622 * LV_RF * QSS_RF * CRANKP / (r_d * TG_RF * TG_RF)


           COEFFNP1_RF = 1.0 + DTPBL * CRANKP * (4.0 * RFEMISSIBR * STBOLT * TG_RF ** 3   &
                        * CT_RF + DENSRF * CPAIR_RF * FHH2 * CT_RF + 2.0 * PI             &
                        * TAUINV) + DTPBL * (CT_RF * LV_RF * CQ2_RF * CQ4_RF)

           COEFFN_RF   = CT_RF * (GSW + RFEMISSIBR * (STBOLT * (4.0 * CRANKP - 1.0)       &
                        * TG_RF*TG_RF*TG_RF*TG_RF + LWDN)                                 & !NET RAD
                        + DENSRF * CPAIR_RF * FHH2 * (TA_RF - (1.0 - CRANKP) * TG_RF)     &
                        - LV_RF * (CQ4_RF * (CQ1_RF - QV_RF)))                            & !SFC HEAT FLUX 
                        - 2.0 * PI * TAUINV * ((1.0 - CRANKP) * TG_RF - T2_RF)  


           TSNEW_RF    = (TG_RF + DTPBL * COEFFN_RF) / COEFFNP1_RF
           !-- FOR SNOW COVERED SURFACE TEMPERATURE IS NOT MORE THAN ZERO
           IF (XICE1 .GT. 0.5) TSNEW_RF = AMIN1(TSNEW_RF,273.15)                          ! Re-added Jan 2010 to keep ice surface at or below freezing (J. Pleim)
           TSHLF_RF = 0.5 * (TSNEW_RF + TG_RF)
           T2NEW_RF = (T2_RF + DBLE(DTPBL) * TAUINV * T2TFACU * (TSHLF_RF - (1 - CRANKP) * T2_RF)   &
                    + DTPBL*T2NUD)  &              ! Added deep temperature nudging 
                    / (1.0 + DBLE(DTPBL) * TAUINV * T2TFACU * CRANKP)
           !-- REPLACE OLD with NEW Value
           TG_RF = TSNEW_RF
           T2_RF = T2NEW_RF

           !-- Diagnostic Outputs
           LWDN_RF   = RFEMISSIBR * LWDN 
           LWUP_RF   = RFEMISSIBR * STBOLT * (TG_RF**4)
           RADNET_RF = GSW - (RFEMISSIBR * (STBOLT * (TG_RF**4) - LWDN))   
 


           !--------------------------------------------------------------------------------------------------------------
           !-- UTK Urban, correction of shortwave radiation that reaches urban ground and wall surface (Radiation Model)
           !--------------------------------------------------------------------------------------------------------------

           !-- Urban morphology parameters estimation 
!           BBDG   = (LAMP / LAMF) * HGTBDG
!           WSTT   = (1.0 / SQRT(LAMP) - 1.0) * BBDG
           SKVFRD = SQRT((HGTBDG/WSTT)**2 + 1.0) - (HGTBDG/WSTT) 
           SKVFWL = 0.5*( (HGTBDG/WSTT) + 1.0 - SQRT((HGTBDG/WSTT)**2 + 1.0) ) / (HGTBDG/WSTT)

           !-- Initialize 
           RADST     = 0.0
           RADST_FNL = 0.0
           GSW_URB   = 0.0 
           RR        = 0.0
           MW        = 0.0
           
           RADWL     = 0.0
           RADWL_FNL = 0.0  
           GSW_WL    = 0.0
           RW        = 0.0
           MR        = 0.0 

           DECLN = 0.0
           EOT   = 0.0
!           HRA   = 0.0 
           AZMA  = 0.0

           IF (SOLDN .GT. 0.0 .AND. COSZ .GT. 0.001) THEN               
              TANZ = ABS((SQRT(1.0-COSZ*COSZ))/COSZ)

              ! This DECLN calculation is based on subroutine radconst in module_radiation_driver.F
              IF (JULIAN .GE. 80.) THEN
                 SXLONG = DPD * (JULIAN - 80.)
              ELSE              
                 SXLONG = DPD * (JULIAN + 285.)
              ENDIF
              DECLN = ASIND(SIND(23.45)*SIND(SXLONG))   ! Declination angle in degrees


              ! This AZMA calculation is based on subroutine toposhad in module_radiation_driver.F
              CSZA = SIND(XLATIJ) * SIND(DECLN) + COSD(XLATIJ) * COSD(DECLN) * COSD(HRA)
              ARGU = (CSZA * SIND(XLATIJ) - SIND(DECLN)) / (SIN(ACOS(CSZA))*COSD(XLATIJ))
              IF (ARGU .GT. 1.) ARGU = 1.
              IF (ARGU .LT. -1.) ARGU = -1.
              AZMA = SIGN(ACOS(ARGU),SIND(HRA)) + PI  ! azimuth angle of the sun
              IF (COSAIJ .GE. 0.) THEN
                 AZMA = AZMA + ASIN(SINAIJ)    ! rotation towards WRF grid 
              ELSE
                 AZMA = AZMA + PI - ASIN(SINAIJ)
              ENDIF
              AZMA = AZMA * RADDEG    ! Solar azimuth angle in degrees 
               
              IF (HRA .GT. 0.0) AZMA = 360.0-AZMA
             RADST = 0.0
             RADWL = 0.0
             DO idir=1,2   ! 2 street dir
             IF(idir.eq.1) Then
                STANGX = float(STDIR1)
             Else
                STANGX = float(STDIR2)
                if (STANGX .GT. 360.0) STANGX = STANGX - 360.
             EndIf
              IF (HRA .LE. 0.0) THEN       ! Morning hours   
                IF (STANGX .GE. 0.0 .AND. STANGX .LE. 180.0) THEN
                  IF ( TANZ .LT. (WSTT/(ABS(SIND(STANGX-AZMA))*HGTBDG)) )THEN    ! Radiation reaching on streets
                    RADST = RADST +SWDDIR*(1.0 - (HGTBDG/WSTT)*TANZ*ABS(SIND(STANGX-AZMA)))        ! Direct solar flux received by the road
                    RADWL = RADWL +SWDDIR*(0.5*TANZ*ABS(SIND(STANGX-AZMA)))                        ! Direct solar flux received by the walls
                  ELSE 
                    RADST = RADST +0.0
                    RADWL = RADWL +SWDDIR*(0.5*(WSTT/HGTBDG))
                  ENDIF
                ENDIF

                IF (STANGX .GT. 180.0 .AND. STANGX .LE. 360.0) THEN
                  IF ( TANZ .LT. (WSTT/(ABS(SIND(ABS(180.0-STANGX)-AZMA))*HGTBDG)) )THEN  
                    RADST = RADST +SWDDIR*(1.0 - (HGTBDG/WSTT)*TANZ*ABS(SIND(ABS(180.0-STANGX)-AZMA)))
                    RADWL = RADWL +SWDDIR*(0.5*TANZ*ABS(SIND(ABS(180.0-STANGX)-AZMA)))
                  ELSE 
                    RADST = RADST +0.0
                    RADWL = RADWL +SWDDIR*(0.5*(WSTT/HGTBDG))      
                  ENDIF
                ENDIF

              ELSE IF (HRA .GT. 0.0) THEN  ! Afternoon hours 
                IF (STANGX .GE. 0.0 .AND. STANGX .LE. 180.0) THEN
                  IF ( TANZ .LT. (WSTT/(ABS(SIND(ABS(180.0+STANGX)-AZMA))*HGTBDG)) )THEN                 
                    RADST = RADST +SWDDIR*(1.0 - (HGTBDG/WSTT)*TANZ*ABS(SIND(ABS(180.0+STANGX)-AZMA)))
                    RADWL = RADWL +SWDDIR*(0.5*TANZ*ABS(SIND(ABS(180.0+STANGX)-AZMA)))
                  ELSE
                    RADST = RADST +0.0
                    RADWL = RADWL +SWDDIR*(0.5*(WSTT/HGTBDG))
                  ENDIF
                ENDIF 

                IF (STANGX .GT. 180.0 .AND. STANGX .LE. 360.0) THEN
                  IF ( TANZ .LT. (WSTT/(ABS(SIND(STANGX-AZMA))*HGTBDG)) )THEN      
                    RADST = RADST +SWDDIR*(1.0 - (HGTBDG/WSTT)*TANZ*ABS(SIND(STANGX-AZMA)))
                    RADWL = RADWL +SWDDIR*(0.5*TANZ*ABS(SIND(STANGX-AZMA)))
                  ELSE 
                    RADST = RADST +0.0
                    RADWL = RADWL +SWDDIR*(0.5*(WSTT/HGTBDG))    
                  ENDIF    
                ENDIF 
              ENDIF 
!              if(HGTBDG.GT.100.0)print *,' px-HGTBDG,stdir1,BB1',HGTBDG,stdir1,BB1
             ENDDO       ! 2 street dir
              RADST = RADST * 0.5
              RADWL = RADWL * 0.5
              SWDIFWL = SWDDIF * 2.*SKVFWL     ! diffuse SW x Skyview for 2 walls
              SWDIFRD = SWDDIF * SKVFRD     ! diffuse SW x Skyview for road
              RR = (RDALBBR*RADST + RDALBBR*SWDIFRD)                        ! First reflected part from the road            
              RW = (WLALB*RADWL + WLALB*SWDIFWL)                            ! First reflected part from the wall 

              MW = ( RW + SKVFWL*WLALB*RR ) / ( 1.0 - (1.0 - 2.0*SKVFWL)*WLALB + (1.0 - SKVFRD)*SKVFWL*RDALBBR*WLALB )                                      ! Multiple reflections (infinite times) from the walls 
              MR = ( RR + (1.0 - SKVFRD)*RDALBBR*(RW + (SKVFWL*WLALB*RR)) ) / ( 1.0 - (1.0 - 2.0*SKVFWL)*WLALB + (1.0 - SKVFRD)*SKVFWL*RDALBBR*WLALB )      ! Multiple reflections (infinite times) from the road

              RADST_FNL = (1.0 - RDALBBR)*RADST + (1.0 - RDALBBR)*SWDIFRD + (1.0 - RDALBBR)*(1.0 - SKVFRD)*MW                                            ! Total SW radiation absorbed by the road surface
              RADWL_FNL = (1.0 - WLALB)*RADWL + (1.0 - WLALB)*SWDIFWL + (1.0 - WLALB)*(1.0 - 2.0*SKVFWL)*MW + (1.0 - WLALB)*SKVFWL*MR                    ! Total SW radiation absorbed by the walls surface

              GSW_URB = URBF*RADST_FNL + (1.0-URBF)*SOLDN*(1.0-ALB)        ! Total SW radiation absorbed at the ground in an urban grid cell
              GSW_WL  = RADWL_FNL                                          ! Total SW radiation absorbed at the walls in an urban grid cell
!            IF(HGTBDG.GT.40.0.AND.SWDDIR.GT.200.)print*,' HGTBDG,SWDDIR=',HGTBDG,SWDDIR,RADWL
           ENDIF


           !------------------------------------------------------------------------------------------
           !-- Computation of Urban Ground Temperature, Solving by CRANK-NIC Implicit method:
           !-- TENDTS=CT*(RADNET-HFX-QFX)
           !-- Calculate the coefficients for implicit calculation of Urban Ground (STREET+VEG) TEMP.
           CQ1_G      = (1.0 - 0.622 * LV * CRANKP / (r_d * TG)) * QSS
           CQ2_G      = 0.622 * LV * QSS * CRANKP / (r_d * TG * TG)

           COEFFNP1_G = 1.0 + DTPBL * CRANKP * (4.0 * RDEMISSIBR * STBOLT * TG**3 * CT                &
                       - 4.0 * RDEMISSIBR * (1.0 - WLEMISSI) * (1.0 - SKVFRD) * SKVFWL * STBOLT       &
                       * ((1.0 - SNOW_FRA) * RDEMISSIBR * TG**3 + SNOW_FRA * SNEMISSI * TG**3) * CT   &
                       + DENS1 * CPAIR / RAH * CPOT * CT + 2.0 * PI                                   &
                       * TAUINV) + DTPBL * (CT * LV * CQ2_G * CQ4_ST)

           COEFFN_G   = CT * (GSW_URB + RDEMISSIBR * STBOLT * (4.0 * CRANKP - 1.0)            &
                       * TG*TG*TG*TG + RDEMISSIBR * SKVFRD * LWDN                                            & 
                       + RDEMISSIBR * WLEMISSI * (1.0 - SKVFRD) * STBOLT * TG_WL**4           & 
                       + RDEMISSIBR * (1.0 - WLEMISSI) * (1.0 - SKVFRD) * SKVFWL * LWDN                      &
                       + RDEMISSIBR * WLEMISSI * (1.0 - WLEMISSI) * (1.0 - SKVFRD) * (1.0 - 2.0*SKVFWL)      &      
                       * STBOLT * TG_WL**4                                                    &
                       - RDEMISSIBR * (1.0 - WLEMISSI) * (1.0 - SKVFRD) * SKVFWL * STBOLT     &
                       * ((1.0 - SNOW_FRA) * RDEMISSIBR * (4.0 * CRANKP - 1.0) * TG**4        & 
                       + SNOW_FRA * SNEMISSI * (4.0 * CRANKP - 1.0) * TG**4)                  & !NET RAD
                       + DENS1 * CPAIR / RAH * (THETA1 - (1.0 - CRANKP) * THETAG)             &
                       - LV * (CQ4_ST * (CQ1_G - QV1)))      & !SFC HEAT FLUX
                       - 2.0 * PI * TAUINV * ((1.0 - CRANKP) * TG - T2)                         !SOIL FLUX

           TSNEW_G    = (TG + DTPBL * COEFFN_G) / COEFFNP1_G
           !-- FOR SNOW COVERED SURFACE TEMPERATURE IS NOT MORE THAN ZERO
           IF (XICE1 .GT. 0.5) TSNEW_G = AMIN1(TSNEW_G,273.15)                         ! Re-added Jan 2010 to keep ice surface at or below freezing (J. Pleim)
           TSHLF_G = 0.5 * (TSNEW_G + TG)
           T2NEW_G = (T2 + DBLE(DTPBL) * TAUINV * T2TFACU * (TSHLF_G - (1 - CRANKP) * T2)               &
                    + DTPBL*T2NUD)  &              ! Added deep temperature nudging 
                    / (1.0 + DBLE(DTPBL) * TAUINV * T2TFACU * CRANKP)
           !-- REPLACE OLD with NEW Value
           TGU = TSNEW_G
           T2U = T2NEW_G

           !-- Diagnostic Outputs
           LWDN_URB   = RDEMISSIBR * SKVFRD * LWDN + RDEMISSIBR * (1.0 - WLEMISSI) * (1.0 - SKVFRD) * SKVFWL * LWDN
           LWUP_URB   = RDEMISSIBR * STBOLT *(TGU**4) + RDEMISSIBR * WLEMISSI * (1.0 - SKVFRD) * STBOLT * TG_WL**4                     &
                        + RDEMISSIBR * WLEMISSI * (1.0 - WLEMISSI) * (1.0 - SKVFRD) * (1.0 - 2.0*SKVFWL) * STBOLT * TG_WL**4          &
                        + RDEMISSIBR * (1.0 - WLEMISSI) * (1.0 - SKVFRD) * SKVFWL * STBOLT * ((1.0 - SNOW_FRA) * RDEMISSIBR * TGU**4   &
                         + SNOW_FRA * SNEMISSI * TGU**4)
           RADNET_URB = GSW_URB + RDEMISSIBR * SKVFRD * LWDN + RDEMISSIBR * (1.0 - WLEMISSI) * (1.0 - SKVFRD)                        &
                        * SKVFWL * LWDN - RDEMISSIBR * STBOLT * (TGU**4)                                                             &
                        + RDEMISSIBR * WLEMISSI * (1.0 - SKVFRD) * STBOLT * TG_WL**4 + RDEMISSIBR * WLEMISSI * (1.0 - WLEMISSI)       &
                         * (1.0 - SKVFRD) * (1.0 - 2.0*SKVFWL) * STBOLT * TG_WL**4    &
                        + RDEMISSIBR * (1.0 - WLEMISSI) * (1.0 - SKVFRD) * SKVFWL * STBOLT * ((1.0 - SNOW_FRA) * RDEMISSIBR * TGU**4   &
                         + SNOW_FRA * SNEMISSI * TGU**4) 


           !-----------------------------------------------------------------------------------
           !-- Computation of Urban Walls Temperature, Solving by CRANK-NIC Implicit method:
           !-- TENDTS=CT*(RADNET-HFX-QFX)
           !-- Calculate the coefficients for implicit calculation of WALLS TEMP.

           COEFFNP1_WL = 1.0 + DTPBL * CRANKP * (4.0 * WLEMISSI * STBOLT * TG_WL**3 * CT_WL                             &
                        - 4.0 * WLEMISSI**2 * (1.0 - 2.0*SKVFWL) * STBOLT * TG_WL**3 * CT_WL                            &
                        - 4.0 * WLEMISSI**2 * (1.0 - WLEMISSI) * (1.0 - 2.0*SKVFWL)**2 * STBOLT * TG_WL**3 * CT_WL      &
                        - 4.0 * WLEMISSI**2 * (1.0 - RDEMISSIBR) * SKVFWL * (1.0 - SKVFRD) * STBOLT * TG_WL**3 * CT_WL  &
                        + ETAC * CT_WL + 2.0 * PI * TAUINV)   

           COEFFN_WL   = CT_WL * (GSW_WL + WLEMISSI * STBOLT * (4.0 * CRANKP - 1.0)               &
                        * TG_WL*TG_WL*TG_WL*TG_WL + WLEMISSI * SKVFWL * LWDN                            &
                        + WLEMISSI * SKVFWL * STBOLT * ((1.0 - SNOW_FRA) * RDEMISSI               &   
                        * TGU**4 + SNOW_FRA * SNEMISSI * TGU**4)                                    &
                        - WLEMISSI**2 * (1.0 - 2.0*SKVFWL) * STBOLT * (4.0 * CRANKP - 1.0) * TG_WL**4   &
                        + WLEMISSI * (1.0 - RDEMISSIBR) * SKVFWL * SKVFRD * LWDN                        &
                        + WLEMISSI * (1.0 - WLEMISSI) * SKVFWL * (1.0 - 2.0*SKVFWL) * LWDN              &
                        - WLEMISSI**2 * (1.0 - WLEMISSI) * (1.0 - 2.0*SKVFWL)**2                  &
                        * STBOLT * (4.0 * CRANKP - 1.0) * TG_WL**4                                      &  
                        - WLEMISSI**2 * (1.0 - RDEMISSIBR) * SKVFWL * (1.0 - SKVFRD)              &
                        * STBOLT * (4.0 * CRANKP - 1.0) * TG_WL**4                                      &
                        + WLEMISSI * (1.0 - WLEMISSI) * SKVFWL * (1.0 - 2.0*SKVFWL)               &
                        * STBOLT * ((1.0 - SNOW_FRA) * RDEMISSI * TGU**4 + SNOW_FRA * SNEMISSI * TGU**4)  & !NET RAD    
                        + ETAC * (TA_CAN - (1.0 - CRANKP) * TG_WL))                               & !SFC HEAT FLUX
                        - 2.0 * PI * TAUINV * ((1.0 - CRANKP) * TG_WL - T2_WL)


           TSNEW_WL    = (TG_WL + DTPBL * COEFFN_WL) / COEFFNP1_WL
           !-- FOR SNOW COVERED SURFACE TEMPERATURE IS NOT MORE THAN ZERO
           IF (XICE1 .GT. 0.5) TSNEW_WL = AMIN1(TSNEW_WL,273.15)                         ! Re-added Jan 2010 to keep ice surface at or below freezing (J. Pleim)
           TSHLF_WL = 0.5 * (TSNEW_WL + TG_WL)
           T2NEW_WL = (T2_WL + DBLE(DTPBL) * TAUINV * T2TFACU * (TSHLF_WL - (1 - CRANKP) * T2_WL)           &
                    + DTPBL*T2NUD)  &              ! Added deep temperature nudging
                    / (1.0 + DBLE(DTPBL) * TAUINV * T2TFACU * CRANKP)
           !-- REPLACE OLD with NEW Value
           TG_WL = TSNEW_WL
           T2_WL = T2NEW_WL 

           !-- Diagnostic Outputs
           LWDN_WL   = WLEMISSI * SKVFWL * LWDN + WLEMISSI * (1.0 - RDEMISSIBR) * SKVFWL * SKVFRD * LWDN    &
                       + WLEMISSI * (1.0 - WLEMISSI) * SKVFWL * (1.0 - 2.0*SKVFWL) * LWDN 
           LWUP_WL   = WLEMISSI * STBOLT * (TG_WL**4) + WLEMISSI * SKVFWL * STBOLT * ((1.0 - SNOW_FRA) * RDEMISSI      &
                        * TGU**4 + SNOW_FRA * SNEMISSI * TGU**4)   &
                       + WLEMISSI**2 * (1.0 - 2.0*SKVFWL) * STBOLT * TG_WL**4 + WLEMISSI**2 * (1.0 - WLEMISSI) * (1.0 - 2.0*SKVFWL)**2 * STBOLT * TG_WL**4     &
                       + WLEMISSI**2 * (1.0 - RDEMISSIBR) * SKVFWL * (1.0 - SKVFRD) * STBOLT * TG_WL**4     &
                       + WLEMISSI * (1.0 - WLEMISSI) * SKVFWL * (1.0 - 2.0*SKVFWL) * STBOLT * ((1.0 - SNOW_FRA) * RDEMISSI * TGU**4 + SNOW_FRA * SNEMISSI * TGU**4)   
           RADNET_WL = GSW + WLEMISSI * SKVFWL * LWDN + WLEMISSI * (1.0 - RDEMISSIBR) * SKVFWL * SKVFRD * LWDN         &
                       + WLEMISSI * (1.0 - WLEMISSI) * SKVFWL * (1.0 - 2.0*SKVFWL) * LWDN    &
                       - WLEMISSI * STBOLT * (TG_WL**4) + WLEMISSI * STBOLT * (TG_WL**4) + WLEMISSI * SKVFWL * STBOLT * ((1.0 - SNOW_FRA) * RDEMISSI           &
                        * TGU**4 + SNOW_FRA * SNEMISSI * TGU**4)   &
                       + WLEMISSI**2 * (1.0 - 2.0*SKVFWL) * STBOLT * TG_WL**4 + WLEMISSI**2 * (1.0 - WLEMISSI) * (1.0 - 2.0*SKVFWL)**2 * STBOLT * TG_WL**4     &
                       + WLEMISSI**2 * (1.0 - RDEMISSIBR) * SKVFWL * (1.0 - SKVFRD) * STBOLT * TG_WL**4     &
                       + WLEMISSI * (1.0 - WLEMISSI) * SKVFWL * (1.0 - 2.0*SKVFWL) * STBOLT * ((1.0 - SNOW_FRA) * RDEMISSI * TGU**4 + SNOW_FRA * SNEMISSI * TGU**4)  


        ENDIF   ! Urban

        IF (IFLAND .LT. 1.5) THEN         ! For Non-urban ground surface
          !-- SOLVE BY CRANK-NIC --   TENDTS=CT*(RADNET-HFX-QFX)-SOILFLX 
          !-- Calculate the coefficients for implicit calculation of TG
          CQ1      = (1.0 - 0.622 * LV * CRANKP / (r_d * TG)) * QSS
          CQ2      = 0.622 * LV * QSS * CRANKP / (r_d * TG * TG)
          CQ3BG    = DENS1 * (1.0 - VEGFRC) / (RAW + RSOIL)
          CQ3VW    = DENS1 * VEGFRC*SIGG / RAW
          CQ3VG    = DENS1 * VEGFRC / (RAW + RSOIL + RINC)
          CQ3      = CQ3BG + CQ3VW + CQ3VG

          COEFFNP1 = 1.0 + DTPBL * CRANKP * (4.0 * EMISSI *  STBOLT * TG ** 3    &
                     * CT + DENS1 * CPAIR / RAH * CPOT * CT + 2.0 * PI           &
                     * TAUINV) + DTPBL * (CT * LV * CQ2 * (CQ3 + CQ4))
          COEFFN   = CT * (GSW + EMISSI * (STBOLT * (4.0 * CRANKP - 1.0)         &
                     * TG*TG*TG*TG + LWDN)                                       & !NET RAD
                     + DENS1 * CPAIR / RAH * (THETA1 - (1.0 - CRANKP) * THETAG)  &
                     - LV * (CQ3 * (CQ1 - QV1) + CQ4 * (CQ1 - QV1)))             & !SFC HEAT FLUX 
                     - 2.0 * PI * TAUINV * ((1.0 - CRANKP) * TG - T2)              !SOIL FLUX
          TSNEW    = (TG + DTPBL * COEFFN) / COEFFNP1
          !-- FOR SNOW COVERED SURFACE TEMPERATURE IS NOT MORE THAN ZERO
          IF (XICE1 .GT. 0.5) TSNEW = AMIN1(TSNEW,273.15)                         ! Re-added Jan 2010 to keep ice surface at or below freezing (J. Pleim)
          TSHLF = 0.5 * (TSNEW + TG)
          T2NEW = (T2 + DBLE(DTPBL) * TAUINV * T2TFAC * (TSHLF - (1 - CRANKP) * T2)     &
                   + DTPBL*T2NUD)  &              ! Added deep temperature nudging 
                   / (1.0 + DBLE(DTPBL) * TAUINV * T2TFAC * CRANKP)
          !-- REPLACE OLD with NEW Value
          IF (ISPXURB .AND. (URBF .GT. 0.0)) THEN  
            TG = (1.-MAX(URBF,0.0))*TSNEW + URBF*TGU
            T2 = (1.d0-MAX(DBLE(URBF),0.0d0))*T2NEW + DBLE(URBF)*T2U 
          ELSE
            TG = TSNEW
            T2 = T2NEW
!          if(I.eq.18.and.j.eq.131)  PRINT *, 'T2,URBF,T2NEW,TSHLF,I,J,ISPXURB,ITIMESTEP',      &
!                               T2,URBF,T2NEW,TSHLF,I,J,ISPXURB,ITIMESTEP
          ENDIF
        ENDIF
        !-----------------------------------------------------------------------------------------

        !-----------------------------------------------------------------------------------------
        ! Compute new subsurface soil and canopy moisture values DENS1. No change required over ocean.
        !--For Urban-vegetation grid -----------------------   !---UTK-Urban
        IF (IFLAND .LT. 1.5.AND. XICE1.LT.0.5) THEN     !------------(1) 
          IF (ISPXURB .AND. (URBF .GT. 0.0)) THEN    ! For Urban Land 
             !-- Compute WR_ST
             ROFF = 0.0    
             WRMAX_ST = 0.247E-3     ! (meter) From Nakayoshi et al., 2009

             IF (WRMAX_ST .GT. 0.0) THEN
               !--  PC is precip. intercepted by street.(M/S)
               PC    = STFA * SIGFF * PRECIP
               DWR   = (WRMAX_ST - WR_ST) / DTPBL                       ! the tendency to reach max.
               PNET  = PC - ER_ST/ DENW                                 ! residual of precip. and evap.
               IF (PNET .GT. DWR) THEN
                 ROFF = PNET - DWR
                 PC   = PC - ROFF
               ENDIF
               IF (QSS .LT. QV1) THEN
                 TENDWR   = PC - ER_ST / DENW
                 WRNEW_ST = WR_ST + DTPBL * TENDWR
               ELSE
                 COF1     = DENS1 / DENW * STFA * (QSS - QV1) / RAW
                 !-- using delta=wr/wrmax
                 CFNP1WR  = 1.0 + DTPBL * COF1 * CRANKP / WRMAX_ST  
                 CFNWR    = PC - COF1 * (1.0 - CRANKP) * WR_ST / WRMAX_ST
                 WRNEW_ST = (WR_ST + DTPBL * CFNWR) / CFNP1WR
               ENDIF
             ELSE
               PC = 0.0
               WRNEW_ST = 0.0  
             ENDIF 
             !------------------------------
             !.. new value             
             WR_ST = AMIN1(WRMAX_ST,WRNEW_ST)
     
             !----------------------------------------------  

             !-- Compute WR_RF
             ROFF = 0.0    
             WRMAX_RF = 0.239E-3     ! (meter) From Nakayoshi et al., 2009
 
             IF (WRMAX_RF .GT. 0.0) THEN
               !--  PC is precip. intercepted by roof.(M/S)
               PC    = RFFA * SIGFF * PRECIP
               DWR   = (WRMAX_RF - WR_RF) / DTPBL                    ! the tendency to reach max.
               PNET  = PC - ER_RF/ DENW                              ! residual of precip. and evap.
               IF (PNET .GT. DWR) THEN
                 ROFF = PNET - DWR
                 PC   = PC - ROFF
               ENDIF
               IF (QSS_RF .LT. QV_RF) THEN      
                 TENDWR   = PC - ER_RF / DENW
                 WRNEW_RF = WR_RF + DTPBL * TENDWR
               ELSE
                 COF1      = DENS1 / DENW * RFFA * (QSS_RF - QV_RF) / RAW
                 !-- using delta=wr/wrmax
                 CFNP1WR   = 1.0 + DTPBL * COF1 * CRANKP / WRMAX_RF  
                 CFNWR     = PC - COF1 * (1.0 - CRANKP) * WR_RF / WRMAX_RF
                 WRNEW_RF  = (WR_RF + DTPBL * CFNWR) / CFNP1WR
               ENDIF
             ELSE
               PC = 0.0
               WRNEW_RF = 0.0  
             ENDIF 
             !------------------------------
             !.. new value             
             WR_RF = AMIN1(WRMAX_RF,WRNEW_RF)

             !----------------------------------------------
    
          ENDIF   ! Urban Part
 
          !-- Compute WR (Default PX)
          ROFF  = 0.0
          WRMAX = 0.2E-3 * VEGFRC * LAI                     ! max. WRMAX IN m
          
          IF (WRMAX.GT.0.0) THEN
            !--  PC is precip. intercepted by veg.(M/S)
            PC    = VEGFRC * SIGF * PRECIP
            DWR   = (WRMAX - WR) / DTPBL                       ! the tendency to reach max.
            PNET  = PC - ER/ DENW                              ! residual of precip. and evap.
            IF (PNET .GT. DWR) THEN
              ROFF = PNET - DWR
              PC   = PC - ROFF
            ENDIF
            IF (QSS .LT. QV1) THEN
              TENDWR = PC - ER / DENW
              WRNEW  = WR + DTPBL * TENDWR
            ELSE
              COF1    = DENS1 / DENW * VEGFRC * (QSS - QV1) / RAW
              !-- using delta=wr/wrmax
              CFNP1WR = 1.0 + DTPBL * COF1 * CRANKP / WRMAX  
              CFNWR   = PC - COF1 * (1.0 - CRANKP) * WR / WRMAX
              WRNEW   = (WR + DTPBL * CFNWR) / CFNP1WR
            ENDIF
          ELSE
            PC    = 0.0
            WRNEW = 0.0  
          ENDIF
          !---------------------------------------------
          !-- Compute W2
          PG     = DENW * (PRECIP - PC)               ! PG is precip. reaching soil (PC already including ROFF)
          TENDW2 = 1.0 / (DENW * DS2) * (PG - EG - ETR) &
                   - C3/DS2 * TAUINV * AMAX1(0.0,(W2 - WFC)) &
                   + (W2NUDG + WGNUDG) / DS2                    ! NUDGING
          W2NEW  = W2 + DTPBL * TENDW2
          W2NEW  = AMIN1(W2NEW,WSAT)
          W2NEW  = AMAX1(W2NEW,WRES) !0.05)  !Limei 08/02/2017
          W2HLF  = 0.5 * (W2 + W2NEW)
          !.. new values
          W2     = W2NEW
          WR     = AMIN1(WRMAX,WRNEW)
        ENDIF       !------------------(1)

        !-----------------------------------------------------------------------------------------

        !-----------------------------------------------------------------------------------------
        ! Compute new surface soil moisture values (WG).
        IF (IFLAND .LT. 1.5 .AND. XICE1 .LT. 0.5) THEN  ! over ocean no change to wg w2,wr
          !--    FOR SNOW COVERED SURFACE, ASSUME SURFACE IS SATURATED AND
          !      WG AND W2 ARE NOT CHANGED
          IF (ISNOW .GT.0.5) THEN
            WG = WSAT
          ELSE
            W2REL = W2HLF / WSAT
            IF (WG .GT. WWLT) THEN
              C1 = DS1*C1SAT * (WSAT / WG) ** (0.5 * B + 1.0) 
            ELSE   ! revise C1 for wg < wilting point Noilhan & Mahfouf (1996)
              ZY2   = C1SAT * (WSAT / WWLT) ** (0.5 * B + 1.0)
              C1MAX = (1.19*WWLT - 5.09)*TG - 146.*WWLT + 1786.
              C1MAX = MAX(MAX(C1MAX,ZY2),10.)
!* Giard-Bazile formulation (resolution of a second order equation)
!
              ZLY =   LOG( C1MAX/10.)
              ZZA = - LOG( ZY2  /10.)
              ZZB = 2. * WWLT * ZLY
              ZDEL = 4. * (ZLY+ZZA) * ZLY * WWLT**2
              ZA  = (-ZZB+SQRT(ZDEL)) / (2.*ZZA)
              ZB  = ZA**2 / ZLY
              C1 = DS1*C1MAX * EXP(-(WG-ZA)**2/ZB)
            ENDIF
            C2 = C2R * W2HLF / (WSAT - W2HLF + 1.E-11) 
            IF (W2HLF .GE. WSAT) THEN
              WEQ = WSAT
            ELSE
              WEQ = W2HLF - AS * WSAT * W2REL ** JP *         &
                   (1.0 - W2REL ** (8.0 * JP))
            ENDIF

            !.... The diffusion method in Sakaguchi and Zeng [2009] (JGR-Atmos.)
            CFNP1 = 1.0 + DTPBL * C2 * TAUINV * CRANKP
            CFN   = C1 / (DENW * DS1) * (PG - EG) - C2 * TAUINV *               &
                    ((1.0 - CRANKP) * WG - WEQ) + WGNUDG/ DS1

            WGNEW = AMAX1((WG + DTPBL * CFN) / CFNP1, WRES)
            !-- NEW VALUES
            WG    = AMIN1(WGNEW, WSAT)
          ENDIF                  !endif for ISNOW
        ENDIF           !----------------(2)
            
      END SUBROUTINE surfpx
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
      SUBROUTINE QFLUX(URBF, DTPBL, ISPXURB, DENS1, DENSRF,                  & !in
                       QV1, QV_RF, TA1, RG, RAW, QSS, QSS_RF,                & !in
                       VEGFRC, ISNOW, ISTI, IFLAND, LAI, BETAP,              & !in
                       WG, W2, WR, WR_ST, WR_VEG, WR_RF, STFA, VEGFA, RFFA,  & !in and UTK-in
                       RSTMIN, WWLT, WFC, RSOIL, RINC, SIGG,                 & !in
                       EG, ER, ETR, ER_ST, ER_RF,                            & !out and UTK-out
                       CQ4, CQ4_RF, CQ4_ST, RS, FASS                         ) !out

!-------------------------------------------------------------------------
!                                                              
!  FUNCTION:                                                          
!    THIS SUBROUTINE COMPUTES EVAPORATION FROM BARE SOIL (EG) AND FROM
!    THE WET PART OF CANOPY (ER) AND TRANSPIRATION FROM THE DRY PART OF
!    CANOPY (ETR).                       
!                                        
!  REVISION HISTORY: 
!    A. Xiu       Oct  2004  - adapted from the PX LSM in MM5 for the WRF system                            
!    R. Gilliam   Dec  2006  - Completed WRF V2.1.2 implementation                                 
!                                                               
!-------------------------------------------------------------------------
!  QFLUX ARGUMENT LIST:                                                        
!-------------------------------------------------------------------------
! INPUTS:
!-- DENS1        air density at first layer
!-- QV1          air humidity at first layer
!-- TA1          air temperature at first layer
!-- RG           shortwave radition reaching the ground
!-- RAW          RA+RB for moisture
!-- QSS          saturation mixing ratio at ground
!-- VEGFRC       vegetation coverage
!-- ISNOW        if snow on the ground 
!-- ISTI         soil type
!-- IFLAND       if land (=1) or water (=2)
!-- LAI          leaf area index
!-- BETAP        
!-- WG           soil moisture at first soil layer
!-- W2           soil moisture at root zone
!-- WR           Canopy water
!
!  OUTPUTS FROM QFLUX:
!-- EG           evaporation from ground (bare soil)
!-- ER           evaporation from canopy
!-- ETR          transpiration from vegetation
!-- CQ4          CQ4 is used for the implicit calculation of TG in SURFACE
!-- RS           surface resistence
!-- FASS         parameter for soil moisture nudging
!-- SIGG         near ground-first layer moisture difference for evap or dew (1)
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      IMPLICIT NONE

! DECLARATIONS - INTEGER
      INTEGER, INTENT(IN) :: ISTI

! DECLARATIONS - REAL
      REAL, INTENT(IN)    :: ISNOW, IFLAND
      REAL, INTENT(IN)    :: DENS1, QV1, TA1, RG, RAW, QSS,            &
                             VEGFRC, LAI,                              &
                             WG, W2, WR, RSTMIN, QV_RF, QSS_RF, DENSRF
      REAL, INTENT(IN)    :: WR_ST, WR_VEG, WR_RF                  ! UTK Urban in
      REAL, INTENT(INOUT) :: BETAP, RSOIL
      REAL, INTENT(IN)    :: WWLT, WFC, RINC
      REAL, INTENT(IN)    :: URBF, DTPBL, STFA, VEGFA, RFFA        ! UTK Urban in
      REAL, INTENT(OUT)   :: EG, ER, ETR, CQ4, RS, FASS, SIGG
      REAL, INTENT(OUT)   :: ER_ST,  ER_RF, CQ4_RF, CQ4_ST   ! UTK out
      LOGICAL, INTENT(IN) :: ISPXURB

!... Local Variables

!... Real
      REAL :: WRMAX, DELTA, RADL, RADF, W2AVAIL, W2MXAV
      REAL :: FTOT, F1, F2, F3,  F4
      REAL :: FSHELT, GS, GA, FX
      REAL :: PAR, F1MAX
      REAL :: WRMAX_ST, WRMAX_VEG, WRMAX_RF
      REAL :: SIGG1, SIGG2, SIGG3

!... Parameters
      REAL, PARAMETER :: RSMAX = 5000.0                ! s/m
      REAL, PARAMETER :: FTMIN = 0.0000001             ! m/s
      REAL, PARAMETER :: F3MIN = 0.25

!... for water surface, no canopy evaporation and transpiration

      !-- Initialize
      ER  = 0.0
      ETR = 0.0
      CQ4 = 0.0

      ER_ST  = 0.0
      ER_RF  = 0.0
      CQ4_RF = 0.0           ! Coeff. for urban roof water in Crank-Nicolson
      CQ4_ST = 0.0           ! Coeff. for urban street water in Crank-Nicolson

!... GROUND EVAPORATION (DEPOSITION)        [From soil]

      IF (QSS .LT. QV1) RSOIL = 0.0

      EG = DENS1 * (QSS - QV1) * ((1.0 - VEGFRC) / (RAW + RSOIL) + VEGFRC / (RAW + RINC + RSOIL))


!!----------------------------------------------------------------------------
!... Evaporation from Urban Street , Veg., Roof Canopy    !-----UTK Urban    
      IF (IFLAND .LT. 1.5 .AND. VEGFRC .GT. 0.0) THEN     !----------------(3) 
        IF (ISPXURB .AND. (URBF .GT. 0.0)) THEN      ! For Urban Land 
          WRMAX_ST = 0.247E-3      ! in (m) unit from Nakayoshi et al., 2009
          IF (WR_ST .LE. 0.0) THEN
             DELTA = 0.0  
          ELSE
             DELTA = WR_ST / WRMAX_ST
          ENDIF

          IF (QSS .GE. QV1) THEN
            SIGG1 = DELTA 
          ELSE
            SIGG1 = 1.0
          END IF
     
          ER_ST  = DENS1 * STFA * SIGG1 * (QSS - QV1) / RAW
          CQ4_ST = DENS1 * SIGG1 / RAW     


          ! From Roof Surface
          WRMAX_RF = 0.239E-3      ! in (m) unit from Nakayoshi et al., 2009
          IF (WR_RF .LE. 0.0) THEN
            DELTA = 0.0  
          ELSE
            DELTA = WR_RF / WRMAX_RF
          ENDIF

          IF (QSS_RF .GE. QV_RF) THEN
            SIGG3 = DELTA 
          ELSE
            SIGG3 = 1.0
          END IF
     
          ER_RF  = DENSRF * RFFA * SIGG3 * (QSS_RF - QV_RF) / RAW
          CQ4_RF = DENSRF * SIGG3 / RAW


!... VEG CANOPY only  
     
        ENDIF   ! FUrban Part

        WRMAX = 0.2E-3 * VEGFRC * LAI   ! in unit m
        IF (WR .LE. 0.0) THEN
          DELTA = 0.0
        ELSE
          DELTA = WR / WRMAX           ! referred to SiB model
        ENDIF
        
        IF (QSS .GE. QV1) THEN
          SIGG = DELTA
        ELSE
          SIGG = 1.0
        ENDIF

        ER = DENS1 * VEGFRC * SIGG * (QSS - QV1) / RAW  

      ENDIF        !---------------(3)
      !!---------------------------------------------------------------------

      !!---------------------------------------------------------------------
      !-- TRANSPIRATION from Vegetation
      !!---------------------------------------------------------------------
      IF (IFLAND .LT. 1.5 .AND. VEGFRC .GT. 0.0)  THEN      !--------------(4)
        
        !!!-RADIATION
        IF (RSTMIN .GT. 130.0) THEN
          F1MAX = 1. - 0.02 * LAI    !Echer2015  Trees
        ELSE
          F1MAX = 1. - 0.07 * LAI    !Echer2015 crops/grass
        ENDIF
        PAR = 0.45 * RG * 4.566  ! converted from W/m2 to umoles/m2/s  (1/.219)  Echer2015
        F1 = F1MAX * (1.0 - exp(-0.0017 * PAR))   !Echer2015
        F1 = AMAX1(F1, RSTMIN / RSMAX)

        !!!-SOIL MOISTURE
        W2AVAIL = W2 - WWLT
        W2MXAV  = WFC - WWLT
        F2      = 1.0 / (1.0 + EXP(-5.0 * ( W2AVAIL / W2MXAV -               &
                  (W2MXAV / 3.0 + WWLT))))    ! according JP, 9/94

        !-AIR TEMP
        !... according to Avissar (1985) and AX 7/95
        IF (TA1 .LE. 302.15) THEN
          F4 = 1.0 / (1.0 + EXP(-0.41 * (TA1 - 282.05)))
        ELSE
          F4 = 1.0 / (1.0 + EXP(0.5 * (TA1 - 314.0)))
        ENDIF

        FTOT = LAI * F1 * F2 * F4
      ENDIF       !---------------(4)


      IF (IFLAND .LT. 1.5 .AND. VEGFRC .GT. 0.0)  THEN      !---------------(5)
        FSHELT = 1.0   ! go back to NP89
        GS     = FTOT / (RSTMIN * FSHELT)
        GA     = 1.0 / RAW
        !-- Compute humidity effect according to RH at leaf surf
        F3 = 0.5 * (GS - GA + SQRT(GA * GA + GA * GS *                       &
             (4.0 * QV1 / QSS - 2.0) + GS * GS)) / GS
        F3 = AMIN1(AMAX1(F3,F3MIN),1.0)
        RS = 1.0 / (GS * F3)
        
        !--- Compute Assimilation factor for soil moisture nudging - jp 12/94
        !--  Note that the 30 coef is to result in order 1 factor for max
        IF (RG .LT. 0.00001) THEN   ! do not nudge during night
          FX = 0.0
        ELSE
          FX = 30.0 * F1 * F4 * LAI / (RSTMIN * FSHELT)
        ENDIF

        FASS = FX

        ETR = DENS1 * VEGFRC * (1.0 - SIGG) * (QSS - QV1) / (RAW + RS)
        !..... CQ4 is used for the implicit calculation of TG in SURFACE
        CQ4 = DENS1 * VEGFRC * ((1.0 - SIGG) / (RAW + RS) + SIGG / RAW)
      ENDIF      !---------------(5)

    END SUBROUTINE qflux
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
      SUBROUTINE SMASS (ISTI,  FASS,  RG,   VEGFRC, RA,              & !in                      
                        WWLT, WFC,                                   & !in
                        ALPH1, ALPH2, BET1, BET2, T2NUDF )             !out  
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
      IMPLICIT NONE
!------------------------------------------------------------------------------------------
!     SMASS COMPUTES SOIL MOISTURE NUDGING COEFFICIENTS
!------------------------------------------------------------------------------------------
!
!.........Arguments
      INTEGER, PARAMETER             :: NSCAT    = 16     ! max. soil types
      
      INTEGER,                   INTENT(IN)   :: ISTI 
      REAL,                      INTENT(IN)   :: FASS, RG, VEGFRC, RA
      REAL,                      INTENT(IN)   :: WWLT, WFC
      REAL,                      INTENT(OUT)  :: ALPH1, ALPH2, BET1, BET2, T2NUDF


!........Local variables

!... Real
      REAL    :: FBET, FALPH, FRA, FTEXT
      
!... Parameters
      REAL, PARAMETER  :: A1MAX = -10.E-5, A2MAX = 1.E-5  ! m/K, m for 6hr period
      REAL, PARAMETER  :: B1MAX = -10.E-3, B2MAX = 1.E-3  ! m/K, m (Bouttier et al 1993)
      REAL, PARAMETER  :: TASSI = 4.6296E-5               ! 1/6hr in 1/sec
      REAL, PARAMETER  :: RAMIN = 10.0                    ! 0.1 s/cm
      REAL, PARAMETER  :: WFCX  = 0.243                   ! middle of WFC range 
      REAL, PARAMETER  :: WWLTX = 0.169                   ! middle of WWLT range 

!
      FBET  = FASS
      FALPH = RG / 1370.0
!--TEXTURE FACTOR NORMALIZED BY LOAM (IST=5)
      FRA   = RAMIN / RA            ! scale by aerodynamic resistance
      FTEXT = TASSI * (WWLT + WFC) / (WWLTX + WFCX) * FRA
!        write(6,*) ' ftot, fbet=',ftot, fbet,' ftext=',ftext/tassi
!
      ALPH1 = A1MAX * FALPH * (1.0 - VEGFRC) * FTEXT
      ALPH2 = A2MAX * FALPH * (1.0 - VEGFRC) * FTEXT
      BET1  = B1MAX * FBET  *        VEGFRC  * FTEXT
      BET2  = B2MAX * FBET  *        VEGFRC  * FTEXT
      T2NUDF = 1.0E-5 * ( VEGFRC*MAX((1.0 - 5.0 * FALPH),0.0) + (1-VEGFRC) )  !  T2 Nudging at night and day for non-veg frac - jp 10/30/14

    END SUBROUTINE smass
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

      SUBROUTINE SOILPROP (SOILCBOT,WEIGHT, ITIMESTEP, MAVAIL, &  ! IN
                           PXLSM_SMOIS_INIT,                   &  ! IN
                           FWSAT,FWFC,FWWLT,FCLAY,FCSAND,      &  ! OUT
                           FFMSAND,FB,FCGSAT,                  &  ! OUT
                           FJP,FAS,FC2R,FC1SAT,FWRES,FC3,ISTI, &  ! OUT
                           WG, W2    )                            ! OUT
     !------------------------------------------------------------------------
     !     SOILPROP COMPUTES SOIL PARAMETERS FOR BOTH BOTTOM AND TOP LAYERS
     !              USING FRACTIONAL SOIL TYPE. A HARD CODED OPTION IS AVAILIABLE
     !              TO COMPUTE THE SOIL PARAMETERS USING FRACTIONAL INFORMATION, OR
     !              TO JUST USE SOIL PARAMETERS OF THE DOMINANT SOIL TYPE
     !------------------------------------------------------------------------
     !-- SOIL PARAMETERS ARE SPECIFIED BY SOIL TYPE:
     !   #  SOIL TYPE  WSAT  WFC  WWLT    B   CGSAT   JP   AS   C2R  C1SAT WRES
     !   _  _________  ____  ___  ____  ____  _____   ___  ___  ___  _____ ____
     !   1  SAND       .395 .135  .068  4.05  3.222    4  .387  3.9  .082  0.020
     !   2  LOAMY SAND .410 .150  .075  4.38  3.057    4  .404  3.7  .098  0.035
     !   3  SANDY LOAM .435 .195  .114  4.90  3.560    4  .219  1.8  .132  0.041
     !   4  SILT LOAM  .485 .255  .179  5.30  4.418    6  .105  0.8  .153  0.015
     !   5  SILT       .485 .255  .179  5.30  4.418    6  .105  0.8  .153  0.015
     !   6  LOAM       .451 .240  .155  5.39  4.111    6  .148  0.8  .191  0.027
     !   7  SND CLY LM .420 .255  .175  7.12  3.670    6  .135  0.8  .213  0.068
     !   8  SLT CLY LM .477 .322  .218  7.75  3.593    8  .127  0.4  .385  0.040
     !   9  CLAY LOAM  .476 .325  .250  8.52  3.995   10  .084  0.6  .227  0.075
     !  10  SANDY CLAY .426 .310  .219 10.40  3.058    8  .139  0.3  .421  0.109
     !  11  SILTY CLAY .482 .370  .283 10.40  3.729   10  .075  0.3  .375  0.056
     !  12  CLAY       .482 .367  .286 11.40  3.600   12  .083  0.3  .342  0.090
     !------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
      IMPLICIT NONE

     !.........Arguments
      INTEGER, PARAMETER             :: NSCAT    = 16     ! max. soil types
      INTEGER, PARAMETER             :: NSCATMIN = 16     ! min soil types

      INTEGER,                      INTENT(IN)   :: WEIGHT, ITIMESTEP, PXLSM_SMOIS_INIT
      REAL,                         INTENT(IN)   :: MAVAIL
      REAL,     DIMENSION(1:NSCAT), INTENT(IN)   :: SOILCBOT
      REAL,                         INTENT(OUT)  :: FWSAT,FWFC,FWWLT,FCLAY,       &
                                                    FCSAND,FFMSAND,FB,FCGSAT,     &
                                                    FJP,FAS,FC2R,FC1SAT,FWRES,FC3
      REAL,                         INTENT(INOUT)  :: W2, WG


      INTEGER,                      INTENT(OUT)  :: ISTI
!........Local variables
      CHARACTER*4, AVCLASS
      CHARACTER*4, DIMENSION( 1: NSCAT ) ::  TEXID 
!... Integer
      INTEGER:: S
!... Real
      REAL:: TFRACBOT, CFRAC, SUMCSND, SUMFMSND, SUMCLY,        &
             AVS, AVCS, AVFMS, AVC, AVSLT,                      &
             SSMPOT, DSMPOT         ! saturated and air-dry soil matric potential

      REAL,    DIMENSION( 1: NSCAT ) :: WSAT, WFC, WWLT, B, CGSAT, AS,  &
                                        JP, C2R, C1SAT, WRES

      REAL,    DIMENSION( 1: NSCATMIN ) ::  CSAND,FMSAND, CLAY
!.......... DATA statement for SOIL PARAMETERS for the 11 soil types
!...........Follow Menut et al., 2013 JGR
      DATA CSAND  /46.0,40.0,29.0, 0.0, 0.0,                        &
                    0.0,29.0, 0.0, 0.0, 0.0,                        &
                    0.0, 0.0, 0.0, 0.0,46.0,                        &
                    0.0/
      DATA FMSAND /46.0,40.0,29.0,17.0,10.0,                        &
                   43.0,29.0,10.0,32.0,52.0,                        &
                    6.0,22.0,43.0,43.0,46.0,                        &
                   32.0/

      DATA CLAY   / 3.0, 4.0,10.0,13.0, 5.0,                        &
                   18.0,27.0,34.0,34.0,42.0,                        &
                   47.0,58.0,18.0,18.0, 3.0,                        &
                   34.0/

      DATA TEXID/'Sand','Lsan','Sloa','Sill','Silt',                &
                 'Loam','Sclo','Sicl','Cllo','Sacl',                &
                 'Sicy','Clay','Ormt','Wate','Bedr',                &
                 'Othe'/                                         
!-- Removed soil parameter lookup table data since these parameters are now computed analytically (Noilhan and Mahfouf 1996)

      DSMPOT  =  -1.0E7          ! mm air-dry water matric potential
!   
!-------------------------------Exicutable starts here--------------------

!     Compute soil characteristics by sand (coarse and fine-medium) and clay fraction
    
      CFRAC      = 0.0                                                               
      SUMCSND    = 0.0
      SUMFMSND   = 0.0
      SUMCLY     = 0.0                                                               
      TFRACBOT   = 0.0                                                              
                                                                                
      DO S = 1,NSCAT                                                        

        TFRACBOT  = TFRACBOT  + SOILCBOT(S)
        SUMCSND   = SUMCSND   + CSAND(S)  * SOILCBOT(S)
        SUMFMSND  = SUMFMSND  + FMSAND(S) * SOILCBOT(S)
        SUMCLY    = SUMCLY    + CLAY(S)   * SOILCBOT(S)
        
        IF(SOILCBOT(S).GE.CFRAC) THEN   ! Find Dominant Category and fraction
          ISTI=S
          CFRAC=SOILCBOT(S)
        ENDIF

      ENDDO

         
      IF(TFRACBOT.GT.0.001) THEN                                                  
        AVCS  = SUMCSND / TFRACBOT
        AVFMS = SUMFMSND / TFRACBOT
        AVS   = AVCS + AVFMS                                               

        AVC   = SUMCLY / TFRACBOT                                           
        AVSLT = 100.0 - AVS - AVC                                   
        
        IF(AVS.GT.(85.+ 0.5*AVC)) THEN
          AVCLASS= 'Sand'                                                   
          ISTI = 1
        ELSE IF(AVS.GT.(70.+ AVC)) THEN
          AVCLASS= 'Lsan'                                             
          ISTI = 2
        ELSE IF((AVC.LT.20..AND.AVS.GT.52.) &
                .OR.(AVC.LE.7.5.AND.AVSLT.LT.50.)) THEN                
          AVCLASS= 'Sloa'                                             
          ISTI = 3
        ELSE IF(AVC.LT.35..AND.AVS.GT.45..AND.AVSLT.LT.28.) THEN      
          AVCLASS= 'Sclo'                                                  
          ISTI = 7
        ELSE IF(AVC.GE.35..AND.AVS.GT.45.) THEN                           
          AVCLASS = 'Sacl'                                                  
          ISTI = 10
        ELSE IF(AVC.LT.27.0.AND.AVSLT.LT.50.) THEN                        
          AVCLASS= 'Loam'                                                
          ISTI = 6
        ELSE IF(AVC.LT.12..AND.AVSLT.GT.80.) THEN                
          AVCLASS = 'Silt'                                         
          ISTI = 5
        ELSE IF(AVC.LT.27.) THEN                                                 
          AVCLASS = 'Sill'                                                       
          ISTI = 4
        ELSE IF(AVC.LT.40..AND.AVS.GT.20.) THEN                                  
          AVCLASS = 'Cllo'                                                       
          ISTI = 9 
        ELSE IF(AVC.LT.40.) THEN                                                 
          AVCLASS = 'Sicl'                                                       
          ISTI = 8       
        ELSE IF(AVSLT.GE.40.) THEN                                               
          AVCLASS = 'Sicy'                                                       
          ISTI = 11            
        ELSE                                                                     
          AVCLASS = 'Clay'                                                       
          ISTI = 12                        
        ENDIF                                                                          
      ELSE
        ! set no soil to 9 - clay loam
        ISTI = 9
        AVCLASS = TEXID(ISTI)  

        AVCS  = CSAND(ISTI)
        AVFMS = FMSAND(ISTI)
        AVS   = AVCS + AVFMS

        AVC   = CLAY(ISTI)
        AVSLT = 100.0 - AVS - AVC

      ENDIF

      FCSAND  = AVCS
      FFMSAND = AVFMS
      FCLAY   = AVC
  
      ! Continous formulation of secondary soil parmeters (Noilhan and Mahfouf 1996)
      FWSAT = (-1.08 * AVS + 494.305) * 1.0E-3
      FWWLT = 37.1342E-3 * SQRT(AVC)
      FWFC  = 89.0467E-3 * AVC**0.3496
      FB    = 0.137 * AVC + 3.501
      FCGSAT= -1.557E-2 * AVS - 1.441E-2 * AVC + 4.7021
      FC1SAT= (5.58 * AVC + 84.88) * 1.0E-2
      FC2R  = 13.815 * AVC**(-0.954)
      FC3   = 5.327 * AVC **(-1.043)
      FAS   = 732.42E-3 * AVC **(-0.539)
      FJP   = 0.134 * AVC + 3.4
      FWRES = 0.00123 * AVC - 0.00066 * AVSLT + 0.0405  !J. Pleim fitted function
      FWRES = AMAX1(FWRES, 0.01)                        !L. Ran set minimum
    
    ! Compute W2 using soil moisture availiability if pxlsm_smois_init (in namelist) is not zero
    IF (ITIMESTEP .EQ. 1 .AND. PXLSM_SMOIS_INIT .GT. 0) THEN
       WG = FWWLT + (0.5*(FWSAT+FWFC) - FWWLT)  * MAVAIL 
       W2 = FWWLT + (0.5*(FWSAT+FWFC) - FWWLT)  * MAVAIL 
    ENDIF                                                                    
    
    END SUBROUTINE soilprop
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
      SUBROUTINE PXSNOW (ITIMESTEP, ASNOW, CSNOW, SNOW,       &
                         SNOWH, SNUP,                         &
                         ALB, SNOALB, SHDFAC, SHDMIN,         & 
                         HC_SNOW, SNOW_FRA, SNOWC, SNOWALB)                   
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!     Pleim-Xiu LSM Simple Snow model
!---------------------------------------------------------------------------------------------------
!     ITIMESTEP -- Model time step index
!     ASNOW --  Analyzed snow water equivalent (mm)
!     CSNOW --  Current snow water equivalent (mm)
!     SNOW  --  Snow water equivalent in model (mm)
!     SNOWH --  Physical snow depth (m)
!     SNUP --   Physical snow depth (landuse dependent) where when below, snow cover is fractional
!
!     HC_SNOW -- Heat capacity of snow as a function of depth
!     SNOW_FRA-- Factional snow area 
!     IFSNOW  -- Snow flag
!     SNOWALB -- Snow albedo
!---------------------------------------------------------------------------------------------------

      IMPLICIT NONE

!---  Arguments
      REAL, PARAMETER  :: W2SN_CONV   =   10.0
      REAL, PARAMETER  :: CS_SNOWPACK = 2092.0
      REAL, PARAMETER  :: SALP        =    2.6
!-----------------------------------------------
      INTEGER,       INTENT(IN)    :: ITIMESTEP
      REAL,          INTENT(IN)    :: ASNOW, CSNOW, SNUP, ALB, SNOALB, SHDFAC, SHDMIN
      REAL,          INTENT(INOUT) :: SNOW, SNOWH, SNOWC
      REAL,          INTENT(OUT)   :: HC_SNOW, SNOW_FRA, SNOWALB
!------------------------------------------------------------------------------------

                                                   
!-----------------------------------------------
!    Local variables
     !... Real
     REAL:: CONV_WAT2SNOW, CSNOWW, RHO_SNOWPACK,   &
            LIQSN_RATIO, SNEQV, RSNOW
!-----------------------------------------------
             
     SNEQV=ASNOW*0.001              ! Snow water in meters
     RHO_SNOWPACK = 450                   ! Snowpack density (kg/m3), this should be computed in the future
     LIQSN_RATIO  = DENW/RHO_SNOWPACK     ! Ratio of water density and snowpack density
     
     CONV_WAT2SNOW = LIQSN_RATIO/1000     ! Conversion factor for snow liquid equiv. (mm) to snowpack (m)
      
     SNOW = ASNOW                         ! Set snow water to analysis value for now, implement a nudging scheme later
     SNOWH= SNOW * CONV_WAT2SNOW          ! Conversion of snow water (mm) to snow depth (m)     
     
     
     ! Is snow cover now present. The limit is 0.45 mm of water eqivalent or about 2 inches snow depth 
     SNOWC = 0.0
     IF (SNOWH .GT. 0.005) SNOWC = 1.0  
     
     HC_SNOW = RHO_SNOWPACK * CS_SNOWPACK * SNOWH
     
      IF (SNEQV < SNUP) THEN                                                  
         RSNOW = SNEQV / SNUP                                                    
         SNOW_FRA = 1. - ( EXP ( - SALP * RSNOW) - RSNOW * EXP ( - SALP))          
      ELSE                                                                       
         SNOW_FRA = 1.0                                                            
      END IF   
      
      SNOWC  = SNOW_FRA

      SNOWALB = ALB + SNOWC*(SNOALB-ALB)

    END SUBROUTINE pxsnow
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

END MODULE module_sf_pxlsm

