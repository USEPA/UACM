!WRF:MODEL_LAYER:PHYSICS
!
MODULE module_bl_acm

!  USE module_model_constants

  REAL, PARAMETER      :: RIC    = 0.25                ! critical Richardson number
  REAL, PARAMETER      :: CRANKP = 0.5                 ! CRANK-NIC PARAMETER

CONTAINS

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE ACMPBL(XTIME,    DTPBL,    U3D,   V3D,                  &
                     PP3D,     DZ8W,     TH3D,  T3D,                  &
                     QV3D,     QC3D,     QI3D,  RR3D,                 &
#if (WRF_CHEM == 1)
                     CHEM3D,   VD3D,     NCHEM,                       &  ! For WRF-Chem
                     KDVEL, NDVEL, NUM_VERT_MIX,                      &  ! For WRF-Chem
#endif
                     UST,      HFX,      QFX,   TSK,                  &
                     PSFC,     EP1,      G,                           &
                     ROVCP,    RD,       CPD,                         &
                     PBLH,     KPBL2D,   EXCH_H, EXCH_M, REGIME,      &
                     GZ1OZ0,   WSPD,     PSIM, MUT, RMOL,             &
                     RUBLTEN,  RVBLTEN,  RTHBLTEN,                    &
                     RQVBLTEN, RQCBLTEN, RQIBLTEN,                    &
                     NLCAT,    LANDUSEF,                              &  ! DDY urban
                     COSZ,     SWDOWN, IFUACM,                        &  ! DDY urban
                     HFX_RF,   HFX_WL,   QFX_RF,                      &  ! UTK urban
                     LAMP, BUILD_HEIGHT, LAMF,                        &
                     STDIR1, STDIR2, BLKW1, BLKW2,                    &
                     UST_WL,UST_RF,UST_ST,                            &      
                     ids,ide, jds,jde, kds,kde,                       &
                     ims,ime, jms,jme, kms,kme,                       &
                     its,ite, jts,jte, kts,kte)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!   THIS MODULE COMPUTES VERTICAL MIXING IN AND ABOVE THE PBL ACCORDING TO 
!   THE ASYMMETRICAL CONVECTIVE MODEL, VERSION 2  (ACM2), WHICH IS A COMBINED 
!   LOCAL NON-LOCAL CLOSURE SCHEME BASED ON THE ORIGINAL ACM (PLEIM AND CHANG 1992)
!
!   REFERENCES: 
!   Pleim (2007) A combined local and non-local closure model for the atmospheric
!                boundary layer.  Part1: Model description and testing.  
!                JAMC, 46, 1383-1395
!   Pleim (2007) A combined local and non-local closure model for the atmospheric
!                boundary layer.  Part2: Application and evaluation in a mesoscale
!                meteorology model.  JAMC, 46, 1396-1409
!
!  REVISION HISTORY:
!     AX        3/2005   - developed WRF version based on the MM5 PX LSM
!     RG and JP 7/2006   - Finished WRF adaptation
!     JP 12/2011 12/2011 - ACM2 modified so it's not dependent on first layer thickness. 
!     JP        3/2013   - WRFChem version. Mixing of chemical species are added
!     JP 12/2014         - Km and Kh updated. New Richardson number stability function.
!                          Reduces day and night 2-m temperature bias.
!     JP 12/2014         - Minimum PBL height bug caused PBLH lower than first level thickness in very stable 
!                          conditions. Fixed so PBLH minimum is layer 1 thickness.
!     JP 12/2015         - RMOL updated after initial calculation in the PX-SFCLAY.
!     DDY 1/2016         - Inclusion of the urban momentum drag
!     HF        5/2016   - MPAS version  
!     JP        8/2017   - Z-coord version from MPAS version for hybrid coords  
!     UTK 5/2018         - Inclusion of the urban ground and rooftop heat fluxes
!                        - Inclusion of the urban moisture flux
!     JP        10/2020  - Fixed a bug where KSRC is now an array dimensioned ( its:ite )
!                        - Modified the PBL height algorithm for stable conditions so
!                          that the Richardson number is computed using windspeed in layer
!                          k rather than wind speed difference between layer k and ksrc
!                        - Added EXCH_M to argument list 
!
!**********************************************************************
!   ARGUMENT LIST:
!
!... Inputs:
!-- XTIME           Time since simulation start (min)
!-- DTPBL           PBL time step
!-- U3D             3D u-velocity interpolated to theta points (m/s)
!-- V3D             3D v-velocity interpolated to theta points (m/s)
!-- PP3D            Pressure at half levels (Pa)
!-- DZ8W            dz between full levels (m)
!-- TH3D            Potential Temperature (K)
!-- T3D             Temperature (K)
!-- QV3D            3D water vapor mixing ratio (Kg/Kg)
!-- QC3D            3D cloud mixing ratio (Kg/Kg)
!-- QI3D            3D ice mixing ratio (Kg/Kg)
!-- RR3D            3D dry air density (kg/m^3)
!-- CHEM3D          Chemical species mixing ratios (ppm)  Optional for WRFChem   
!-- VD3D            Dry deposition velocity (m/s)         Optional for WRFChem
!-- NCHEM           Number of chemical species            Optional for WRFChem
!-- UST             Friction Velocity (m/s)
!-- HFX		    Upward heat flux at the surface (w/m^2)
!-- QFX		    Upward moisture flux at the surface (Kg/m^2/s)
!-- TSK             Surface temperature (K)
!-- PSFC            Pressure at the surface (Pa)
!-- EP1             Constant for virtual temperature (r_v/r_d-1) (dimensionless)
!-- G               Gravity (m/s^2)
!-- ROVCP           r/cp
!-- RD              gas constant for dry air (j/kg/k)
!-- CPD             heat capacity at constant pressure for dry air (j/kg/k)
!-- GZ1OZ0          log(z/z0) where z0 is roughness length
!-- WSPD            wind speed at lowest model level (m/s)
!-- PSIM            similarity stability function for momentum
!-- MUT             Total Mu : Psfc - Ptop
 
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
!
!... Outputs: 
!-- PBLH            PBL height (m)
!-- KPBL2D          K index for PBL layer
!-- REGIME          Flag indicating PBL regime (stable, unstable, etc.)
!-- RUBLTEN         U tendency due to PBL parameterization (m/s^2)
!-- RVBLTEN         V tendency due to PBL parameterization (m/s^2)
!-- RTHBLTEN        Theta tendency due to PBL parameterization (K/s)
!-- RQVBLTEN        Qv tendency due to PBL parameterization (kg/kg/s)
!-- RQCBLTEN        Qc tendency due to PBL parameterization (kg/kg/s)
!-- RQIBLTEN        Qi tendency due to PBL parameterization (kg/kg/s)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
     IMPLICIT NONE

!.......Arguments
! DECLARATIONS - INTEGER
    INTEGER, INTENT(IN) :: ids,ide, jds,jde, kds,kde,       &
                           ims,ime, jms,jme, kms,kme,       &
                           its,ite, jts,jte, kts,kte, XTIME

! DECLARATIONS - REAL
    REAL, INTENT(IN)    :: DTPBL, EP1, G, ROVCP, RD, CPD

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                              &
             INTENT(IN)                                :: U3D, V3D,            &
                                                          PP3D, DZ8W, T3D,     &
                                                          QV3D, QC3D, QI3D,    &
                                                          RR3D, TH3D

    REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)    :: PSIM, GZ1OZ0,        &
                                                          HFX, QFX, TSK,       &
                                                          PSFC, WSPD, MUT

    REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: PBLH, REGIME,        & 
                                                          UST, RMOL

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                              &
             INTENT(INOUT)                             :: RUBLTEN, RVBLTEN,    &
                                                          RTHBLTEN, RQVBLTEN,  &
                                                          RQCBLTEN, RQIBLTEN

    real, dimension( ims:ime, kms:kme, jms:jme ),                              &
             intent(inout)                             :: exch_h, exch_m

    INTEGER, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: KPBL2D
    real, dimension( ims:ime, kms:kme, jms:jme ), INTENT(OUT) :: UST_WL
    real, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: UST_RF, UST_ST

#if (WRF_CHEM == 1)
!... Chem
    INTEGER, INTENT(IN) :: nchem, kdvel, ndvel, num_vert_mix
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme, nchem ), INTENT(INOUT) :: CHEM3D
    REAL, DIMENSION( ims:ime, kdvel, jms:jme, ndvel ), INTENT(IN)      :: VD3D
#endif

!-- DDY-UTK urban
    INTEGER, INTENT(IN) :: NLCAT      
    REAL, DIMENSION( ims:ime , 1:NLCAT, jms:jme ), INTENT(IN) :: LANDUSEF
    REAL, DIMENSION(ims:ime,jms:jme),INTENT(IN)               :: COSZ,SWDOWN
    REAL, DIMENSION(ims:ime,jms:jme),INTENT(IN)               :: HFX_RF, HFX_WL, QFX_RF
    REAL, DIMENSION(ims:ime,jms:jme),INTENT(IN)               :: LAMP, BUILD_HEIGHT, LAMF
    INTEGER, DIMENSION(ims:ime,jms:jme),INTENT(IN)            :: STDIR1, STDIR2
    REAL, DIMENSION(ims:ime,jms:jme),INTENT(IN)               :: BLKW1, BLKW2
    INTEGER, INTENT(IN   )                                    :: IFUACM
!... Local Variables

!... Integer
      INTEGER :: I, J, K, L
!... Real
      REAL  RDT
      REAL, PARAMETER :: KARMAN = 0.4

#if (WRF_CHEM == 1)
!... Chem
    REAL,    DIMENSION( ims:ime, kms:kme, nchem ) :: CHEM2D
    REAL,    DIMENSION( ims:ime, kdvel, ndvel ) :: VD2D
#endif


!--------------------Executables start here--------------------------
!--------------------------------------------------------------------

   
   RDT = 1.0 / DTPBL

!==================================================

   DO j = jts,jte   
#if (WRF_CHEM == 1)
      DO L = 1, nchem
      DO K = kms,kme
      DO I = ims, ime
        CHEM2D(i,k,l) = chem3d(i,k,j,l)
      ENDDO
      ENDDO
      ENDDO

      DO L = 1, ndvel
      DO K = 1, kdvel
      DO I = ims, ime
        VD2D(i,k,l) = VD3D(i,k,j,l)
      ENDDO
      ENDDO
      ENDDO
#endif
      CALL ACM2D(j=j,xtime=xtime,dtpbl=dtpbl                       &
              ,us=u3d(ims,kms,j),vs=v3d(ims,kms,j)                 &
              ,theta=th3d(ims,kms,j),tt=t3d(ims,kms,j)             &
              ,qvs=qv3d(ims,kms,j),qcs=qc3d(ims,kms,j)             &
              ,qis=qi3d(ims,kms,j)                                 &
#if (WRF_CHEM == 1)
              ,chem=chem2d                                         &
              ,vd=vd2d                                             &
              ,nchem=nchem,kdvel=kdvel,ndvel=ndvel                 &
              ,num_vert_mix=num_vert_mix                           &
#endif
              ,dzf=DZ8W(ims,kms,j)                                 &
              ,densx=RR3D(ims,kms,j)                               &
              ,utnp=rublten(ims,kms,j),vtnp=rvblten(ims,kms,j)     &
              ,ttnp=rthblten(ims,kms,j),qvtnp=rqvblten(ims,kms,j)  &
              ,qctnp=rqcblten(ims,kms,j),qitnp=rqiblten(ims,kms,j) &
              ,cpd=cpd,g=g,rovcp=rovcp,rd=rd,rdt=rdt               &
              ,psfcpa=psfc(ims,j),ust=ust(ims,j)                   &
              ,pbl=pblh(ims,j)                                     &
              ,exch_hx=exch_h(ims,kms,j)                           &
              ,exch_mx=exch_m(ims,kms,j)                           &
              ,regime=regime(ims,j),psim=psim(ims,j)               &
              ,hfx=hfx(ims,j),hfx_rf=hfx_rf(ims,j)                 & ! UTK urban
              ,hfx_wl=hfx_wl(ims,j),qfx_rf=qfx_rf(ims,j)           & ! UTK urban
              ,qfx=qfx(ims,j)                                      & ! UTK urban
              ,tg=tsk(ims,j),gz1oz0=gz1oz0(ims,j)                  &
              ,wspd=wspd(ims,j) ,klpbl=kpbl2d(ims,j)               &
              ,mut=mut(ims,j), rmol=rmol(ims,j)                    &
              ,ep1=ep1,karman=karman                               &
              ,nlcat=nlcat,landusef=landusef(ims,1,j)              & ! DDY urban
              ,cosz=cosz(ims,j),swdown=swdown(ims,j)               & ! DDY urban
              ,lamp=lamp(ims,j), HGTBDG=BUILD_HEIGHT(ims,j), lamf=lamf(ims,j) &
              ,stdir1=stdir1(ims,j), stdir2=stdir2(ims,j)          &
              ,blkw1=blkw1(ims,j), blkw2=blkw2(ims,j)              &
              ,UST_WL=UST_WL(ims,kms,j),UST_RF=UST_RF(ims,j)       &
              ,UST_ST=UST_ST(ims,j),IFUACM=ifuacm                  &
              ,ids=ids,ide=ide, jds=jds,jde=jde, kds=kds,kde=kde   &
              ,ims=ims,ime=ime, jms=jms,jme=jme, kms=kms,kme=kme   &
              ,its=its,ite=ite, jts=jts,jte=jte, kts=kts,kte=kte   )
#if (WRF_CHEM == 1)
      DO L = 1, nchem
      DO I = ims, ime
        chem3d(i,kms:kme,j,l) = CHEM2D(i,kms:kme,l)
      ENDDO
      ENDDO
#endif
   ENDDO

   END SUBROUTINE ACMPBL
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE ACM2D(j, xtime, dtpbl                                &
              ,us, vs, theta, tt, qvs, qcs, qis                    &
#if (WRF_CHEM == 1)
              ,chem, vd, nchem, kdvel, ndvel                       &
              ,num_vert_mix                                        &
#endif
              ,dzf, densx, utnp, vtnp, ttnp, qvtnp, qctnp, qitnp   &
              ,cpd, g, rovcp, rd, rdt, psfcpa, ust                 &
              ,pbl, exch_hx, exch_mx, regime, psim                 &
              ,hfx, hfx_rf, hfx_wl, qfx_rf, qfx, tg, gz1oz0        &
              ,wspd, klpbl                                         &
              ,mut, rmol                                           &
              ,ep1, karman                                         &
              ,nlcat, landusef                                     & ! DDY urban 
              ,cosz, swdown                                        & ! DDY urban
              ,lamp, HGTBDG, lamf                                  &
              ,stdir1, stdir2, blkw1, blkw2                        &
              ,UST_WL,UST_RF,UST_ST,IFUACM                         &
              ,ids,ide, jds,jde, kds,kde                           &
              ,ims,ime, jms,jme, kms,kme                           &
              ,its,ite, jts,jte, kts,kte                           )

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      IMPLICIT NONE

!.......Arguments

!... Real
      REAL,                                INTENT(IN)    :: DTPBL, G, RD, ep1, karman, CPD, ROVCP, RDT
      REAL, DIMENSION( ims:ime ),          INTENT(INOUT) :: PBL, UST
      
      REAL, DIMENSION( ims:ime, kms:kme ), INTENT(IN)    :: US, VS, THETA, TT,   &
                                                            QVS, QCS, QIS, DENSX
      REAL, DIMENSION( ims:ime, kms:kme ), intent(in)    :: DZF
      REAL, DIMENSION( ims:ime, kms:kme ), intent(inout) :: utnp,  &
                                                            vtnp,  &
                                                            ttnp,  &
                                                            qvtnp, &
                                                            qctnp, &
                                                            qitnp
      real, dimension( ims:ime ), intent(in)     :: psfcpa
      real, dimension( ims:ime ), intent(in)     :: tg
      real, dimension( ims:ime ), intent(inout)  :: regime, rmol
      real, dimension( ims:ime ), intent(in)     :: wspd, psim, gz1oz0
      real, dimension( ims:ime ), intent(in)     :: hfx, qfx, hfx_rf, hfx_wl, qfx_rf
      real, dimension( ims:ime ), intent(in)     :: mut
      real, dimension( ims:ime, kms:kme ),                         &
               intent(inout)                     :: exch_hx, exch_mx
      real, dimension( ims:ime, kms:kme ), INTENT(OUT)    :: UST_WL
      REAL, DIMENSION( ims:ime ), INTENT(OUT)    :: UST_RF, UST_ST
!... Integer
      INTEGER, DIMENSION( ims:ime ), INTENT(OUT) :: KLPBL
      INTEGER, INTENT(IN)                        :: XTIME
      INTEGER, INTENT(IN   )                     :: IFUACM
      integer, intent(in)                        :: ids,ide, jds,jde, kds,kde, &
                                                    ims,ime, jms,jme, kms,kme, &
                                                    its,ite, jts,jte, kts,kte, j
#if (WRF_CHEM == 1)
!....Chem
      INTEGER,  INTENT(IN)      ::      NCHEM, KDVEL, NDVEL, NUM_VERT_MIX
      REAL, DIMENSION( ims:ime, kms:kme, NCHEM ), INTENT(INOUT)  :: CHEM
      REAL, DIMENSION( ims:ime, KDVEL, NDVEL ), INTENT(IN)       :: VD
#endif

!-- DDY urban
      INTEGER, INTENT(IN)                              :: nlcat
      REAL, DIMENSION( ims:ime, 1:nlcat ), INTENT(IN)  :: landusef
      REAL, DIMENSION( ims:ime ), INTENT(IN)           :: COSZ, SWDOWN
      REAL, DIMENSION( ims:ime ), INTENT(IN)           :: lamp, HGTBDG, lamf
      INTEGER, DIMENSION( ims:ime ), INTENT(IN)        :: stdir1, stdir2
      REAL, DIMENSION( ims:ime ), INTENT(IN)           :: blkw1, blkw2
      logical :: restart

!--------------------------------------------------------------------
!--Local 
      INTEGER :: I, K, M
      INTEGER :: KPBLHT
      INTEGER, DIMENSION( its:ite ) :: KPBLH, NOCONV, KSRC, KSL

!... Real
      REAL    ::  TVCON, WSS, TCONV, TH1, TOG, DTMP, WSSQ
      REAL    ::  ZH1,UH1,VH1                                   ! NEW FOR V3.7
      REAL    ::  psix, THV1
      REAL, DIMENSION( its:ite )          :: URBF
      REAL, DIMENSION( its:ite,kts:kte, 1:2 )  :: IMU, EXU, IMV, EXV
      REAL, DIMENSION( its:ite,kts:kte )       :: MIXLKZ
      REAL, DIMENSION( its:ite )          :: FINT, PSTAR, CPAIR, QVS_CAN,        &
                                             CPAIR_RF, CPAIR_CAN, DENSX_CAN
      REAL, DIMENSION( its:ite, kts:kte ) :: THETAV, RIB,                        &
                                             EDDYZ, EDDYZM, UX, VX, THETAX,      &
                                             QVX, QCX, QIX, ZA, DZHI, DZFI, DZH
      REAL, DIMENSION( its:ite, 0:kte )   :: ZF
      REAL, DIMENSION( its:ite )          :: WST, TST, QST, USTM, TSTV,          &
                                             TMPFX_RF, TMPFX_WL, MOISFX_RF
      REAL, DIMENSION( its:ite )          :: MOL
      REAL :: FINTT, ZMIX, UMIX, VMIX
      REAL :: TMPFX, TMPVTCON, TMPP, TMPTHS, TMPTH1, TMPVCONV, WS1, DTH
      REAL :: A, TST12, RL, ZFUNC, DENSF
      REAL       :: TFXT, QFXT      ! For UTK urban


!... Integer
      INTEGER :: KL, jtf, ktf, itf, KMIX, CNP
!...

     !-- DDY-UTK urban
      real :: LUTOT, GAMDIS
      integer :: NL
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      LOGICAL :: eof
      INTEGER :: stat, c2, r
      REAL, DIMENSION( its:ite )          :: UST_CAN, USTT
      LOGICAL :: ISPXURB

!...
      character*512 :: message


!---------Executables start here------------------
!-------------------------------------------------


      !--------------------------------------------------------------------------------------------------------------
      !-- Condition to set the variable 'ISPXURB' for the domain of interest(Domain-4 here)
      !-- DDY-UTK urban (*IMPORTANT NOTE*:-set the time step here for nested domain accord. to the main time step)

      if (ifuacm .eq. 1) then
          ispxurb = .true.
      else
          ispxurb = .false.
      endif

!-----initialize vertical tendencies and

      DO i = its,ite
        DO k = kts,kte
          utnp(i,k) = 0.
          vtnp(i,k) = 0.
          ttnp(i,k) = 0.
        ENDDO
      ENDDO

      DO k = kts,kte
        DO i = its,ite
          qvtnp(i,k) = 0.
        ENDDO
      ENDDO

      DO k = kts,kte
        DO i = its,ite
          qctnp(i,k) = 0.
          qitnp(i,k) = 0.
        ENDDO
      ENDDO

      !-DDY urban

      !-UTK urban
      !-- Urban Roof and Wall fluxes initialization
      DO i = its,ite
        TMPFX_RF(i)  = 0.0
        TMPFX_WL(i)  = 0.0
        MOISFX_RF(i) = 0.0
      ENDDO

      IF (ISPXURB ) THEN
        DO i = its,ite
          !--DDY urban calculate urban fractions for USGS-landuse only (NOW)

          URBF(I) = 0.0
          LUTOT   = 0.0
          DO NL = 1,NLCAT
            LUTOT = LUTOT + LANDUSEF(I,NL)
          ENDDO
          IF (LUTOT .GT. 0.0) THEN
!            URBF(I) = LANDUSEF(I,1)/LUTOT
            URBF(I) = (LANDUSEF(I,24)+LANDUSEF(I,25)+LANDUSEF(I,26))/LUTOT  !-Only for NLCD40 here. 
          ENDIF
          IF(URBF(I).LT.0.05) URBF(I)=-1.0
        ENDDO
      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!- Compute the height of full- and half-sigma level above ground level
!- UTK urban
!- Estimation of Urban Roof level Temp. & Mois. fluxes for ACM module

        DO I = its,ite
          ZF(I,0) = 0.0
        ENDDO
        KSL = 0

        DO K = kts,kte
          DO I = its,ite
            ZF(I,K) = DZF(I,K) + ZF(I,K-1)
            ZA(I,K) = 0.5 * (ZF(I,K) + ZF(I,K-1))
            !--DDY-UTK urban  ! Determine urban canopy layer
            IF (ispxurb .AND. URBF(I) .GT. 0.0)  THEN
!              IF ((h(I,j) .GE. ZF(I,K-1)).AND.(h(I,j) .LE. ZF(I,K))) THEN
              IF ((HGTBDG(I) .GT. ZF(I,K-1)).AND.(HGTBDG(I) .LE. ZF(I,K))) THEN  ! temp hack to run without building data
                GAMDIS = (HGTBDG(I) - ZF(I,K-1)) / DZF(I,K)
                CNP = K
                KSL(I) = CNP
                CPAIR_RF(I)  = CPD * (1.0 + 0.84 * QVS(I,CNP))
                TMPFX_RF(I)  = HFX_RF(I) / (CPAIR_RF(I) * DENSX(I,CNP))
                MOISFX_RF(I) = QFX_RF(I) / DENSX(I,CNP)
  
                QVS_CAN(I)   = 0.0
                DENSX_CAN(I) = 0.0
                UST_CAN(I)   = 0.0
                DO M = kts,CNP-1
                  QVS_CAN(I)   = QVS_CAN(I) + QVS(I,M) * DZF(I,M)  
                  DENSX_CAN(I) = DENSX_CAN(I) + DENSX(I,M) * DZF(I,M)
                  UST_CAN(I) = UST_CAN(I) + UST_WL(I,M)**2 * DZF(I,M)
                ENDDO
                QVS_CAN(I)   = QVS_CAN(I) + QVS(I,CNP) * DZF(I,CNP)* GAMDIS
                DENSX_CAN(I) = DENSX_CAN(I) + DENSX(I,CNP) * DZF(I,CNP)* GAMDIS
                UST_CAN(I) = UST_CAN(I) + UST_WL(I,CNP)**2 * DZF(I,CNP)
                QVS_CAN(I)   = QVS_CAN(I)/HGTBDG(I)   !ZF(I,CNP)
                DENSX_CAN(I) = DENSX_CAN(I)/HGTBDG(I)   !ZF(I,CNP)
                UST_CAN(I) = sqrt(UST_CAN(I)/HGTBDG(I))        !ZF(I,CNP) Fix on July 9, 2024
                CPAIR_CAN(I) = CPD * (1.0 + 0.84 * QVS_CAN(I))    
                TMPFX_WL(I)  = HFX_WL(I) / (CPAIR_CAN(I) * DENSX_CAN(I))
              ENDIF
            ENDIF
          ENDDO
        ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Compute Micromet Scaling variables, not availiable in WRF for ACM
     DO I = its,ite
        CPAIR(I)  = CPD * (1.0 + 0.84 * QVS(I,1))                    ! J/(K KG)
        TMPFX     = HFX(I)  / (cpair(i) * DENSX(I,1))
        TMPVTCON  = 1.0 + EP1 * QVS(I,1)                             ! COnversion factor for virtual temperature
        WS1       = SQRT(US(I,1)**2 + VS(I,1)**2)                    ! Level 1 wind speed
        TST(I)    = -TMPFX / UST(I)
        QST(I)    = -QFX(I) / (UST(I)*DENSX(I,1))
        USTM(I)   = UST(I) * WS1 / wspd(i)
        THV1      = TMPVTCON * THETA(I,1) 
        IF (ispxurb .AND. URBF(I) .GT. 0.0.and.HGTBDG(I) .GT.0.0) THEN
          UST(I) = MAX(UST(I),UST_ST(I))   
!          TSTT(I) = -(TMPFX+TMPFX_RF(I)) / UST(I) 
          TFXT = -(TMPFX+TMPFX_RF(I) + TMPFX_WL(I))
          QFXT = -(QFX(I)+QFX_RF(I)) / (DENSX(I,1))
          TSTV(I) = TFXT*TMPVTCON + THV1*EP1*QFXT
          USTT(I) = sqrt(UST_CAN(I)**2 + UST(I)**2 + UST_RF(I)**2)

          TSTV(I) = (TFXT*TMPVTCON + THV1*EP1*QFXT)/USTT(I)
        ELSE
          TSTV(I) = TST(I)*TMPVTCON + THV1*EP1*QST(I)
          USTT(I) = UST(i)
        ENDIF
        IF (ABS(TSTV(I)) .LT. 1.0E-6) THEN
          TSTV(I) = SIGN(1.0E-6,TSTV(I))
        ENDIF
        MOL(I)    = THV1 * USTT(I)**2/(KARMAN*G*TSTV(I))
        RMOL(I)   = 1.0/MOL(I)
        WST(I)    = USTT(I) * (PBL(I)/(KARMAN*ABS(MOL(I)))) ** 0.333333       
        PSTAR(I)  = MUT(I)/1000.                                     ! P* in cb 
     ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!... Compute PBL height

!... compute the height of full- and half-sigma level above ground level

     DO K = kts,kte
       DO I = its,ite
         DZH(I,K) = ZF(I,K) - ZF(I,K-1)
         DZHI(I,K)= 1./DZH(I,K)
       ENDDO
     ENDDO

     DO K = kts,kte-1
       DO I = its,ite
         DZFI(I,K) = 1./(ZA(I,K+1)-ZA(I,K))
       ENDDO
     ENDDO
     DO I = its,ite
       DZFI(I,KTE) = DZFI(I,KTE-1)
     ENDDO

     DO K = kts,kte
       DO I = its,ite
         TVCON       = 1.0 + EP1 * QVS(I,K)
         THETAV(I,K) = THETA(I,K) * TVCON
       ENDDO
     ENDDO


!...  COMPUTE PBL WHERE RICHARDSON NUMBER = RIC (0.25) HOLTSLAG ET AL 1990  
     DO 100 I = its,ite
       DO K = 1,kte
         RIB(I,K) = 0.0
         KSRC(I) = K
         IF (ZF(I,K) .GT. 30.0) GO TO 69
       ENDDO
69     CONTINUE 

       TH1 = 0.0
       ZH1 = 0.0
       UH1 = 0.0
       VH1 = 0.0
       DO K = 1,KSRC(I)
         TH1 = TH1 + THETAV(I,K)  
         ZH1 = ZH1 + ZA(I,K)
         UH1 = UH1 + US(I,K)
         VH1 = VH1 + VS(I,K)
       ENDDO  
       TH1 = TH1/KSRC(I)
       ZH1 = ZH1/KSRC(I)
       UH1 = UH1/KSRC(I)
       VH1 = VH1/KSRC(I)
       IF(MOL(I).LT.0.0 .AND. XTIME.GT.1) then
         WSS   = (USTT(I) ** 3 + 0.6 * WST(I) ** 3) ** 0.33333
         TCONV = -8.5 * USTT(I) * TSTV(I) / WSS
         TH1   = TH1 + TCONV
       ENDIF

99     KMIX = KSRC(I)
       DO K = KSRC(I),kte
         DTMP   = THETAV(I,K) - TH1
         IF (DTMP.LT.0.0) KMIX = K
       ENDDO
       IF(KMIX.GT.KSRC(I)) THEN
         FINTT = (TH1 - THETAV(I,KMIX)) / (THETAV(I,KMIX+1)               &
               - THETAV(I,KMIX))
         ZMIX = FINTT * (ZA(I,KMIX+1)-ZA(I,KMIX)) + ZA(I,KMIX)
         UMIX = FINTT * (US(I,KMIX+1)-US(I,KMIX)) + US(I,KMIX)
         VMIX = FINTT * (VS(I,KMIX+1)-VS(I,KMIX)) + VS(I,KMIX)
       ELSE
         ZMIX = ZH1
         UMIX = UH1
         VMIX = VH1
       ENDIF
       DO K = KMIX,kte
         DTMP   = THETAV(I,K) - TH1
         TOG = 0.5 * (THETAV(I,K) + TH1) / G
         WSSQ = (US(I,K)-UMIX)**2                                     &
              + (VS(I,K)-VMIX)**2
         IF (KMIX .EQ. KSRC(I)) THEN
               WSSQ = US(I,K)**2 + VS(I,K)**2
               WSSQ = WSSQ + 100.*USTT(I)*USTT(I) 
         ENDIF
         WSSQ = MAX( WSSQ, 0.1 )
         RIB(I,K) = ABS(ZA(I,K)-ZMIX) * DTMP / (TOG * WSSQ)
         IF (RIB(I,K) .GE. RIC) GO TO 201
       ENDDO


       write (message,*)' RIBX never exceeds RIC, RIB(i,kte) = ',rib(i,5),        &
               ' THETAV(i,1) = ',thetav(i,1),' MOL=',mol(i),            &
               ' TCONV = ',TCONV,' WST = ',WST(I),                      &
               ' KMIX = ',kmix,' USTT = ',USTT(I),                       &
               ' TST = ',TST(I),' U,V = ',US(I,1),VS(I,1),              &
               ' I,J=',I,J,'UF=',URBF(I),' TSTV=',TSTV(I)
       CALL wrf_error_fatal ( message )

201    CONTINUE

       KPBLH(I) = K

100  CONTINUE

     DO I = its,ite
       IF (KPBLH(I) .GT. KSRC(I)) THEN
!---------INTERPOLATE BETWEEN LEVELS -- jp 7/93
         FINT(I) = (RIC - RIB(I,KPBLH(I)-1)) / (RIB(I,KPBLH(I)) -       &
                    RIB(I,KPBLH(I)-1))
         IF (FINT(I) .GT. 0.5) THEN
           KPBLHT  = KPBLH(I)
           FINT(I) = FINT(I) - 0.5
         ELSE
           KPBLHT  = KPBLH(I) - 1
           FINT(I) = FINT(I) + 0.5
         ENDIF
         PBL(I)  = FINT(I) * (ZF(I,KPBLHT) - ZF(I,KPBLHT-1)) +          &
                     ZF(I,KPBLHT-1)
         KLPBL(I) = KPBLHT
       ELSE
         KLPBL(I) = KSRC(I)
         PBL(I)    = ZA(I,KSRC(I))                                                  
       ENDIF
     ENDDO

     DO I = its,ite       
       NOCONV(I) = 0
! Check for CBL and identify conv. vs. non-conv cells
       IF (PBL(I) / MOL(I) .LT. -0.02 .AND. KLPBL(I) .GT. 3        &
           .AND. THETAV(I,1) .GT. THETAV(I,2) .AND. XTIME .GT. 1) THEN
          NOCONV(I)   = 1
          REGIME(I) = 4.0                     ! FREE CONVECTIVE - ACM
       ENDIF
     ENDDO

!--- DDY urban, to call ACMURB subroutine
     restart = .false.
     MIXLKZ = -1.0
     if(ispxurb) then
       CALL ACMURB(COSZ, SWDOWN, CPD,                               &
                   DZF, NLCAT, LANDUSEF, DENSX,                     &
                   THETAV, ZA,                                      &
                   US, VS, THETA, XTIME, DTPBL, restart,            &
                   MIXLKZ, IMU, EXU, IMV, EXV,                      & !--DDY urban 
                   UST_WL,UST_RF,UST_ST,                            &
                   ispxurb,lamp, HGTBDG, lamf,                      &
                   stdir1,stdir2,blkw1,blkw2,                       &
                   ims,ime, kms,kme, jms,jme, its,ite, kts,kte,j)
     endif

!... Calculate Kz
     CALL EDDYX(DTPBL, ZF,  ZA, MOL, PBL, UST, USTT, UST_WL,UST_RF,UST_ST,  &
                US,    VS,  TT,  THETAV, DENSX, PSTAR,              &
                QVS,   QCS, QIS, G, RD, CPAIR,                      &
                EDDYZ, EDDYZM,                                      &
                MIXLKZ, ispxurb, KSL,                               & !--DDY urban
                its,ite, kts,kte, ims,ime, kms,kme)

     CALL ACM(DTPBL, PSTAR,  NOCONV, ZF, DZH, DZHI, J,              &
              KLPBL, PBL,   DZFI, MOL,  UST,                        &
              TST, QST,  USTM,   EDDYZ,  DENSX, DENSX_CAN,          &
              THETA,  QVS,    QCS,    QIS,                          &
              THETAX, QVX,    QCX,    QIX,                          &
#if (WRF_CHEM == 1)
              CHEM,  VD,  NCHEM, KDVEL, NDVEL, NUM_VERT_MIX,        &
#endif
              TMPFX_RF, TMPFX_WL, MOISFX_RF, XTIME, ispxurb,        & !DDY-UTK urban
              NLCAT, LANDUSEF, LAMP, HGTBDG, LAMF,                  &
              ids,ide, jds,jde, kds,kde,                            &
              ims,ime, jms,jme, kms,kme,                            &
              its,ite, jts,jte, kts,kte)

     CALL ACMM(DTPBL, PSTAR,  NOCONV, ZF, DZH, DZHI, J,             &
               KLPBL, PBL,   DZFI, MOL,  UST,                       &
               TST, QST,  USTM,  EDDYZM, DENSX,                     &
               US,    VS,                                           &
               UX,    VX,                                           &
               IMU, EXU, IMV, EXV, URBF, XTIME, ispxurb,            & !DDY urban
               ids,ide, jds,jde, kds,kde,                           &
               ims,ime, jms,jme, kms,kme,                           &
               its,ite, jts,jte, kts,kte)

!.. Load exch_h for use in CCN activation

     DO K = kts,kte-1
       DO I = its,ite
         exch_hx(I,K) = EDDYZ(I,K)
         exch_mx(I,K) = EDDYZM(I,K)
       ENDDO
     ENDDO

!... Calculate tendency due to PBL parameterization

     DO K = kts,kte
       DO I = its,ite
         UTNP(I,K)  = UTNP(I,K) + (UX(I,K) - US(I,K)) * RDT
         VTNP(I,K)  = VTNP(I,K) + (VX(I,K) - VS(I,K)) * RDT
         TTNP(I,K)  = TTNP(I,K) + (THETAX(I,K) - THETA(I,K)) * RDT
         QVTNP(I,K) = QVTNP(I,K) + (QVX(I,K) - QVS(I,K)) * RDT
         QCTNP(I,K) = QCTNP(I,K) + (QCX(I,K) - QCS(I,K)) * RDT
         QITNP(I,K) = QITNP(I,K) + (QIX(I,K) - QIS(I,K)) * RDT
       ENDDO
     ENDDO

   END SUBROUTINE ACM2D
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE ACMINIT(RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,           &
                      RQCBLTEN,RQIBLTEN,P_QI,P_FIRST_SCALAR,       &
                      restart, allowed_to_read ,                   &
                      ids, ide, jds, jde, kds, kde,                &
                      ims, ime, jms, jme, kms, kme,                &
                      its, ite, jts, jte, kts, kte                 )
!-----------------------------------------------------------------------
!
!    This subroutine is for preparing ACM PBL variables. 
!    Called from module_physics_init.F
!
!  REVISION HISTORY:
!     AX     3/2005 - Originally developed
!-----------------------------------------------------------------------
!   ARGUMENT LIST:
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   IMPLICIT NONE
!
   LOGICAL, INTENT(IN)          :: restart , allowed_to_read

   INTEGER, INTENT(IN)          :: ids, ide, jds, jde, kds, kde,   &
                                   ims, ime, jms, jme, kms, kme,   &
                                   its, ite, jts, jte, kts, kte

   INTEGER, INTENT(IN)          :: P_QI,P_FIRST_SCALAR

!   REAL, DIMENSION( kms:kme ), INTENT(IN)  :: SHALF
   REAL, DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(OUT) ::  &
                                                         RUBLTEN,  &
                                                         RVBLTEN,  &
                                                         RTHBLTEN, &
                                                         RQVBLTEN, &
                                                         RQCBLTEN, & 
                                                         RQIBLTEN

!... Local Variables
   INTEGER :: i, j, k, itf, jtf, ktf

!
   jtf=min0(jte,jde-1)
   ktf=min0(kte,kde-1)
   itf=min0(ite,ide-1)

   IF(.not.restart)THEN
     DO j=jts,jtf
     DO k=kts,ktf
     DO i=its,itf
        RUBLTEN(i,k,j)=0.
        RVBLTEN(i,k,j)=0.
        RTHBLTEN(i,k,j)=0.
        RQVBLTEN(i,k,j)=0.
        RQCBLTEN(i,k,j)=0.
     ENDDO
     ENDDO
     ENDDO
   ENDIF

   IF (P_QI .ge. P_FIRST_SCALAR .and. .not.restart) THEN
      DO j=jts,jtf
      DO k=kts,ktf
      DO i=its,itf
         RQIBLTEN(i,k,j)=0.
      ENDDO
      ENDDO
      ENDDO
   ENDIF


   END SUBROUTINE acminit
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!-------------------------------------------------------------------          
   SUBROUTINE EDDYX(DTPBL, ZF,  ZA, MOL, PBL, UST,USTT, UST_WL,UST_RF,UST_ST, &
                    US,    VS,  TT,  THETAV, DENSX, PSTAR,            &
                    QVS,   QCS, QIS, G, RD, CPAIR,                    &
                    EDDYZ, EDDYZM,                                    &
                    MIXLKZ,   ispxurb,  KSL,                          & !-- DDY urban
                    its,ite, kts,kte, ims,ime, kms,kme )


!**********************************************************************
!   Two methods for computing Kz:
!   1.  Boundary scaling similar to Holtslag and Boville (1993)
!   2.  Local Kz computed as function of local Richardson # and vertical 
!       wind shear, similar to LIU & CARROLL (1996)
!
!**********************************************************************
!
!-- DTPBL           time step of the minor loop for the land-surface/pbl model
!-- ZF              height of full level
!-- ZA              height of half level
!-- MOL             Monin-Obukhov length in 1D form
!-- PBL             PBL height in 1D form
!-- UST             friction velocity U* in 1D form (m/s)
!-- US              U wind 
!-- VS              V wind
!-- TT              temperature
!-- THETAV          potential virtual temperature
!-- DENSX           dry air density (kg/m^3)
!-- PSTAR           P*=Psfc-Ptop
!-- QVS             water vapor mixing ratio (Kg/Kg)
!-- QCS             cloud mixing ratio (Kg/Kg)
!-- QIS             ice mixing ratio (Kg/Kg)
!-- G               gravity
!-- RD              gas constant for dry air (j/kg/k)
!-- CPAIR           specific heat of moist air (M^2 S^-2 K^-1)
!-- EDDYZ           eddy diffusivity for heat KZ
!-- EDDYZM          eddy diffusivity for momentum KM
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      IMPLICIT NONE

!.......Arguments
  
!... Integer
      INTEGER, INTENT(IN)   :: its, ite, kts, kte, ims, ime, kms, kme
!... Real
      REAL, DIMENSION( ims:ime ),          INTENT(IN)    :: PBL, UST, UST_RF, UST_ST
      REAL,                                INTENT(IN)    :: DTPBL, G, RD
      REAL, DIMENSION( its:ite ),          INTENT(IN)    :: MOL, PSTAR, CPAIR, USTT

      REAL, DIMENSION( ims:ime, kms:kme ), INTENT(IN)    :: US,VS, TT,   &
                                                            QVS, QCS, QIS, DENSX, UST_WL
      REAL, DIMENSION( its:ite, kts:kte ), INTENT(IN)    :: ZA, THETAV
      REAL, DIMENSION( its:ite, 0:kte ),   INTENT(IN)    :: ZF
      
      REAL, DIMENSION( its:ite, kts:kte ), INTENT(OUT)   :: EDDYZ,EDDYZM

      !-- DDY urban 
      REAL, DIMENSION( its:ite, kts:kte ), INTENT(INOUT) :: MIXLKZ
      INTEGER, DIMENSION( its:ite ),       INTENT(IN)    :: KSL

!.......Local variables

!... Integer
      INTEGER  :: ILX, KL, KLM, K, I

!... Real
      REAL     :: ZOVL, PHIH, WT, ZSOL, ZFUNC, DZF, SS, GOTH, EDYZ, DZZ
      REAL     :: RI, QMEAN, TMEAN, XLV, ALPH, CHI, ZK, SQL, DENSF, KZO
      REAL     :: FH, FM
      REAL     :: WM, EDYZM, PHIM, MIXWL, USTK
      REAL, DIMENSION( its:ite ) :: USTK2
!... Parameters
      REAL, PARAMETER :: RV     = 461.5
      REAL, PARAMETER :: RC     = 0.25
      REAL, PARAMETER :: RLAM   = 80.0
      REAL, PARAMETER :: GAMH   = 16.0  ! Dyer74 !15.0  !  Holtslag and Boville (1993)
      REAL, PARAMETER :: GAMM   = 16.0  ! Dyer74
      REAL, PARAMETER :: BETAH  = 5.0   ! Holtslag and Boville (1993)  BETAM = BETAH
      REAL, PARAMETER :: KARMAN = 0.4
      REAL, PARAMETER :: P      = 2.0   ! ZFUNC exponent
      REAL, PARAMETER :: EDYZ0  = 0.01  ! New Min Kz
      REAL, PARAMETER :: PR     = 0.8   ! Prandtl #
      REAL, PARAMETER :: WSTT   = 20.0   !m
      INTEGER, PARAMETER :: imvdif = 1
      logical, intent(in) :: ispxurb
!
      ILX = ite 
      KL  = kte
      KLM = kte - 1

      USTK2 = 0.0
      
      DO K = kts,KLM
        DO I = its,ILX
          EDYZ = 0.0
          ZOVL = 0.0
          DZF  = ZA(I,K+1) - ZA(I,K)
          DZZ  = ZF(I,K) -ZF(I,K-1)
          KZO = EDYZ0
          IF(ispxurb) THEN
            IF(K.LT.KSL(I)) THEN
               USTK2(I) = USTK2(I) + UST_WL(I,K)**2 * DZZ     
               USTK = SQRT( USTK2(I) / ZF(I,K) + UST(I)**2)      
            ELSE 
               USTK = USTT(I)
            ENDIF
          ELSE
            USTK = UST(I)
          ENDIF
          MIXLKZ(I,K) = KARMAN * ZF(I,K)
          MIXWL = 0.0
          IF(ispxurb) MIXWL = 1./(1./ZF(I,K) + 1./WSTT)     
!--------------------------------------------------------------------------
          IF (ZF(I,K) .LT. PBL(I)) THEN
            ZOVL = ZF(I,K) / MOL(I)
            IF (ZOVL .LT. 0.0) THEN
              IF (ZF(I,K) .LT. 0.1 * PBL(I)) THEN
                PHIH = 1.0 / SQRT(1.0 - GAMH * ZOVL)
                PHIM = (1.0 - GAMM * ZOVL)**(-0.25)
                WT   = USTK / PHIH
                WM   = USTK / PHIM
              ELSE
                ZSOL = 0.1 * PBL(I) / MOL(I)
                PHIH = 1.0 / SQRT(1.0 - GAMH * ZSOL)
                PHIM = (1.0 - GAMM * ZSOL)**(-0.25)
                WT   = USTK / PHIH
                WM   = USTK / PHIM
              ENDIF
            ELSE IF (ZOVL .LT. 1.0) THEN
              PHIH = 1.0 + BETAH * ZOVL
              WT   = USTK / PHIH
              WM   = WT
            ELSE
              PHIH = BETAH + ZOVL
              WT   = USTK / PHIH
              WM   = WT
            ENDIF
            !-- DDY, mixing length profile

            ZFUNC = MIXLKZ(I,K)* (1.0 - ZF(I,K) / PBL(I)) ** P

            EDYZ  = WT * ZFUNC
            EDYZM = WM * ZFUNC
          ENDIF
!--------------------------------------------------------------------------
          SS   = ((US(I,K+1) - US(I,K)) ** 2 + (VS(I,K+1) - VS(I,K)) ** 2)   &
                  / (DZF * DZF) + 1.0E-9
          GOTH = 2.0 * G / (THETAV(I,K+1) + THETAV(I,K))
          RI   = GOTH * (THETAV(I,K+1) - THETAV(I,K)) / (DZF * SS)
!--------------------------------------------------------------------------
!         Adjustment to vert diff in Moist air
          IF(imvdif.eq.1)then
            IF ((QCS(I,K)+QIS(I,K)) .GT. 0.01E-3 .OR. (QCS(I,K+1)+             &
                 QIS(I,K+1)) .GT. 0.01E-3) THEN
              QMEAN = 0.5 * (QVS(I,K) + QVS(I,K+1))
              TMEAN = 0.5 * (TT(I,K) + TT(I,K+1))
              XLV   = (2.501 - 0.00237 * (TMEAN - 273.15)) * 1.E6
              ALPH  =  XLV * QMEAN / RD / TMEAN
              CHI   =  XLV * XLV * QMEAN / CPAIR(I) / RV / TMEAN / TMEAN
              RI    = (1.0 + ALPH) * (RI -G * G / SS / TMEAN / CPAIR(I) *       &
                      ((CHI - ALPH) / (1.0 + CHI)))
            ENDIF
          ENDIF
!--------------------------------------------------------------------------
            
          !-- DDY (01/2016)
          SQL = (MIXLKZ(I,K) * RLAM / (RLAM + MIXLKZ(I,K))) ** 2
            
          IF (RI .GE. 0.0) THEN
            FH = 1./(1.+10.*RI+50.*RI**2+5000.*RI**4)+0.0012  !pleim5
            FM = PR*FH + 0.00104

            if(K.EQ.KSL(I)) then
              EDDYZ(I,K) = KZO + SQRT(SS) * FH * SQL + sqrt(UST_WL(I,K)**2+UST_RF(I)**2) * MIXWL
              EDDYZM(I,K) = KZO + SQRT(SS) * FM * SQL + sqrt(UST_WL(I,K)**2+UST_RF(I)**2) * MIXWL
            else
              EDDYZ(I,K) = KZO + SQRT(SS) * FH * SQL + UST_WL(I,K) * MIXWL
              EDDYZM(I,K) = KZO + SQRT(SS) * FM * SQL + UST_WL(I,K) * MIXWL
            endif
          ELSE
            if(K.EQ.KSL(I)) then
              EDDYZ(I,K) = KZO + SQRT(SS*(1.0-25.0*RI))*SQL+sqrt(UST_WL(I,K)**2+UST_RF(I)**2)*MIXWL
            else
              EDDYZ(I,K) = KZO + SQRT(SS * (1.0 - 25.0 * RI)) * SQL + UST_WL(I,K)*MIXWL
            endif
            EDDYZM(I,K) = EDDYZ(I,K) * PR
          ENDIF
  
          IF(EDYZ.GT.EDDYZ(I,K)) THEN
            EDDYZ(I,K) = EDYZ
            EDDYZM(I,K) = MIN(EDYZM,EDYZ*0.8)  !PR
          ENDIF

          EDDYZ(I,K) = MIN(1000.0,EDDYZ(I,K))
          EDDYZ(I,K) = MAX(KZO,EDDYZ(I,K))
          EDDYZM(I,K) = MIN(1000.0,EDDYZM(I,K))
          EDDYZM(I,K) = MAX(KZO,EDDYZM(I,K))

        ENDDO             ! for I loop
      ENDDO               ! for k loop
!
      DO I = its,ILX
        EDDYZ(I,KL) = 0.0
        EDDYZM(I,KL) = 0.0
      ENDDO

   END SUBROUTINE EDDYX
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!-------------------------------------------------------------------          
   SUBROUTINE ACM(DTPBL, PSTAR,  NOCONV, ZF, DZH, DZHI, JX,       &
                  KLPBL, PBL,   DZFI, MOL,  UST,                  &
                  TST, QST,  USTM,   EDDYZ,  DENSX, DENSX_CAN,    &
                  THETA,  QVS,    QCS,    QIS,                    &
                  THETAX, QVX,    QCX,    QIX,                    &
#if (WRF_CHEM == 1)
                  CHEM,  VD, NCHEM, KDVEL, NDVEL,                 &
                  NUM_VERT_MIX,                                   &
#endif
                  TMPFX_RF, TMPFX_WL, MOISFX_RF, XTIME, ispxurb,  & !DDY-UTK urban
                  NLCAT, LANDUSEF, LAMP, HGTBDG, LAMF,            &
                  ids,ide, jds,jde, kds,kde,                      &
                  ims,ime, jms,jme, kms,kme,                      &
                  its,ite, jts,jte, kts,kte)
!**********************************************************************
!   PBL model called the Asymmetric Convective Model, Version 2 (ACM2) 
!   -- See top of module for summary and references
!
!---- REVISION HISTORY:
!   AX     3/2005 - developed WRF version based on ACM2 in the MM5 PX LSM
!   JP and RG 8/2006 - updates
!   JP     3/2013 - Chem additions
!
!**********************************************************************
!  ARGUMENTS:
!-- DTPBL           PBL time step
!-- PSTAR           Psurf - Ptop in cb
!-- NOCONV          If free convection =0, no; =1, yes
!-- ZF              Height for full layer
!-- DZH             Layer thickness --> ZF(K) - ZF(K-1)
!-- DZHI            Inverse of layer thickness
!-- JX              N-S index
!-- KLPBL           PBL level at K index
!-- PBL             PBL height in m
!-- DZFI            Inverse layer thickness --> 1/(Z(K+1)-Z(K))
!-- MOL             Monin-Obukhov length in 1D form
!-- UST             U* in 1D form
!-- TST             Theta* in 1D form
!-- QST             Q* in 1D form
!-- USTM            U* for computation of momemtum flux 
!-- EDDYZ           eddy diffusivity KZ
!-- DENSX           dry air density (kg/m^3)
!-- US              U wind 
!-- VS              V wind
!-- THETA           potential temperature
!-- QVS             water vapor mixing ratio (Kg/Kg)
!-- QCS             cloud mixing ratio (Kg/Kg)
!-- QIS             ice mixing ratio (Kg/Kg)
!-- UX              new U wind 
!-- VX              new V wind
!-- THETAX          new potential temperature
!-- QVX             new water vapor mixing ratio (Kg/Kg)
!-- QCX             new cloud mixing ratio (Kg/Kg)
!-- QIX             new ice mixing ratio (Kg/Kg)
!-- CHEM            Chemical species mixing ratios (ppm)  WRFChem   
!-- VD              Dry deposition velocity (m/s)         WRFChem
!-- NCHEM           Number of chemical species            WRFChem
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      IMPLICIT NONE

!.......Arguments

!... Integer
      INTEGER,  INTENT(IN)      ::      ids,ide, jds,jde, kds,kde, &
                                        ims,ime, jms,jme, kms,kme, &
                                        its,ite, jts,jte, kts,kte, &
                                        JX, XTIME
      INTEGER,  DIMENSION( its:ite ), INTENT(IN)  :: NOCONV
      INTEGER,  DIMENSION( ims:ime ), INTENT(IN)  :: KLPBL

!... Real
      REAL, DIMENSION( ims:ime ),          INTENT(IN)    :: PBL, UST
      REAL,                                INTENT(IN)    :: DTPBL
      REAL, DIMENSION( its:ite ),          INTENT(IN)    :: PSTAR,    &
                                                            MOL, TST, &
                                                            QST, USTM
      REAL, DIMENSION( ims:ime ),          INTENT(IN)    :: LAMP, HGTBDG, LAMF
      REAL, DIMENSION( its:ite, kts:kte ), INTENT(IN)    :: DZHI, DZH, DZFI
      REAL, DIMENSION( its:ite, 0:kte ),   INTENT(IN)    :: ZF
      REAL, DIMENSION( its:ite, kts:kte ), INTENT(INOUT) :: EDDYZ
      REAL, DIMENSION( ims:ime, kms:kme ), INTENT(IN)    :: THETA,    &
                                                            QVS, QCS, QIS, DENSX
      REAL, DIMENSION( its:ite, kts:kte ), INTENT(OUT)   :: THETAX,   &
                                                            QVX, QCX, QIX
#if (WRF_CHEM == 1)
!......Chem
      INTEGER,  INTENT(IN)      ::   NCHEM, KDVEL, NDVEL, NUM_VERT_MIX
      REAL, DIMENSION( ims:ime, kms:kme, nchem ), INTENT(INOUT) :: CHEM   
      REAL, DIMENSION( ims:ime, KDVEL, NDVEL ), INTENT(IN)      :: VD
#endif

!-- UTK urban
      REAL, DIMENSION( its:ite ), INTENT(IN)     :: DENSX_CAN, TMPFX_RF, TMPFX_WL, MOISFX_RF
      INTEGER, INTENT(IN) :: NLCAT
      REAL, DIMENSION( ims:ime,1:NLCAT ), INTENT(IN)  :: LANDUSEF
      INTEGER, DIMENSION( its:ite )  :: KSL, URBF
      logical, intent(in) :: ispxurb
      INTEGER :: CNP

      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      LOGICAL :: eof
      INTEGER :: stat, c2, r
      REAL :: BBDG, WSTT, RBPW, RFFRC, LUTOT

!.......Local variables

!... Parameters
      INTEGER, PARAMETER :: NSP   = 4 
!
      REAL,    PARAMETER :: XX    = 0.5          ! FACTOR APPLIED TO CONV MIXING TIME STEP
      REAL,    PARAMETER :: KARMAN = 0.4

!... Integer
      INTEGER :: ILX, KL, KLM, I, K, NSPX, NLP, NL, JJ, L,LL
      INTEGER :: KCBLMX
      INTEGER, DIMENSION( its:ite ) :: KCBL

!... Real
      REAL                               :: MBMAX, HOVL, MEDDY, MBAR
      REAL                               :: EKZ, RZ, FM, WSPD, DTS, DTRAT, F1
      REAL, DIMENSION( its:ite )         :: PSTARI, FSACM, DTLIM
      REAL, DIMENSION( kts:kte, its:ite) :: MBARKS, MDWN, VRF, NFRF
      REAL, DIMENSION( kts:kte )         :: XPLUS, XMINUS
      REAL                               :: DELC
      REAL, DIMENSION( kts:kte )                :: AI, BI, CI, EI
      REAL, ALLOCATABLE, DIMENSION( : , : )     :: DI, UI    
      REAL, ALLOCATABLE, DIMENSION( : , : )     :: FS
      REAL, ALLOCATABLE, DIMENSION( : , : , : ) :: VCI
      REAL, ALLOCATABLE, DIMENSION( : , : )     :: FS_RF     ! Urban roof fluxes
      REAL, ALLOCATABLE, DIMENSION( : , : )     :: FS_WL     ! Urban wall fluxes

      CHARACTER*80 :: message
      REAL CVOL
!
!--Start Executables ----

!-----

      ILX = ite
      KL  = kte
      KLM = kte - 1
      NSPX = NSP
#if (WRF_CHEM == 1)
      NSPX = NSPX + NUM_VERT_MIX 
#endif

      KCBLMX = 0
      MBMAX  = 0.0
!...Allocate species variables
      ALLOCATE (DI( 1:NSPX,kts:kte ))       
      ALLOCATE (UI( 1:NSPX,kts:kte ))  
      ALLOCATE (FS( 1:NSPX,its:ite )) 
      ALLOCATE (VCI( 1:NSPX,its:ite,kts:kte ))
      ALLOCATE (FS_RF( 1:NSPX,its:ite ))
      ALLOCATE (FS_WL( 1:NSPX,its:ite ))

!---COMPUTE ACM MIXING RATE
      DO I = its,ILX
        DTLIM(I)  = DTPBL
        PSTARI(I) = 1.0 / PSTAR(I)
        KCBL(I)   = 1
        FSACM(I)  = 0.0

        !--DDY urban calculate urban fractions for USGS-landuse only (NOW)
        URBF(I) = 0.0
        LUTOT   = 0.0
        DO NL = 1,NLCAT
          LUTOT = LUTOT + LANDUSEF(I,NL)
        ENDDO
        IF (LUTOT .GT. 0.0) THEN
!          URBF(I) = LANDUSEF(I,1)/LUTOT
          URBF(I) = (LANDUSEF(I,24)+LANDUSEF(I,25)+LANDUSEF(I,26))/LUTOT  !-Only for NLCD40 here. 
        ENDIF
        IF(URBF(I).LT.0.05) URBF(I)=-1.0
        !-- Estimation of the Urban canopy height level (Model-level)
        IF (ISPXURB .AND. (URBF(I) .GT. 0.0)) THEN   !temp hack for testing
          DO K = kts,kte
            IF (HGTBDG(I) .LE. ZF(I,K) .AND. HGTBDG(I) .GE. ZF(I,K-1)) THEN
              KSL(I) = K
            ENDIF
          ENDDO
        ENDIF
! Adjustment for volume in urban canopies
         DO K = kts,kte    
           VRF(K,I) =1.0
           NFRF(K,I) = 1.0 
        ENDDO
        IF (ispxurb .AND. URBF(I) .GT. 0.0) THEN
          CNP = KSL(I)   
          RFFRC = LAMP(I)*URBF(I)
!          FS(1,I) = (1.-RFFRC)*FS(1,I)
!          FS(2,I) = (1.-RFFRC)*FS(2,I)
          DO K = kts,kte  
            IF( K.LT.CNP) then
               VRF(K,I) = 1. -RFFRC
               NFRF(K,I) = 1. -RFFRC
            ELSE
               VRF(K,I) =1.0
               NFRF(K,I) = 1.0 
            ENDIF 
          ENDDO
          VRF(cnp,i) = 1. -(HGTBDG(I)-ZF(I,cnp-1))*DZHI(I,cnp)*RFFRC
!          NFRF(cnp)=VRF(cnp)
        ENDIF

        IF (NOCONV(I) .EQ. 1) THEN
          KCBL(I) = KLPBL(I)

!-------MBARKS IS UPWARD MIXING RATE; MDWN IS DOWNWARD MIXING RATE
!--New couple ACM & EDDY-------------------------------------------------------------
          HOVL     = -PBL(I) / MOL(I)
          FSACM(I) = 1./(1.+((KARMAN/(HOVL))**0.3333)/(0.72*KARMAN))
        IF (ispxurb .AND. URBF(I) .GT. 0.0) THEN
          IF(KCBL(I).gt.CNP)THEN 
            CVOL = (PBL(I) - ZF(i,1)) - RFFRC*(HGTBDG(I)-ZF(i,1))
          ELSEIF (KCBL(I).eq.CNP) Then
            CVOL = DZH(I,CNP)*VRF(CNP,I) + (ZF(i,cnp-1) - ZF(i,1)) * (1. -RFFRC)
          ELSE
            CVOL = (PBL(I) - ZF(i,1)) * (1. -RFFRC)
          endif
        else
          CVOL = (PBL(I) - ZF(i,1))
        endif
          MEDDY    = EDDYZ(I,1) * DZFI(i,1)/VRF(1,I) / CVOL
          MBAR     = MEDDY * FSACM(I)
          DO K = kts,KCBL(I)-1
            EDDYZ(I,K) = EDDYZ(I,K) * (1.0 - FSACM(I))
          ENDDO

          MBMAX = AMAX1(MBMAX,MBAR)

          DO K = kts+1,KCBL(I)
            MBARKS(K,I) = MBAR
            MDWN(K,I)   = MBAR * (PBL(I) - ZF(i,K-1)) * DZHI(i,K)
          ENDDO
          MBARKS(1,I) = MBAR
          MBARKS(KCBL(I),I) = MDWN(KCBL(I),I)
          MDWN(KCBL(I)+1,I) = 0.0
        ENDIF

      ENDDO                              ! end of I loop

      DO K = kts,KLM
        DO I = its,ILX
          EKZ = EDDYZ(I,K) * DZFI(i,K) * DZHI(i,K)
          DTLIM(I) = AMIN1(0.75 / EKZ,DTLIM(I))
        ENDDO
      ENDDO
       
      DO I = its,ILX 
        IF (NOCONV(I) .EQ. 1) THEN
          KCBLMX = AMAX0(KLPBL(I),KCBLMX)
          RZ     = (ZF(i,KCBL(I)) - ZF(i,1)) * DZHI(i,1)
          DTLIM(I)  = AMIN1(XX / (MBARKS(1,I) * RZ),DTLIM(I))
        ENDIF
      ENDDO

      DO K = kts,KL
        DO I = its,ILX
          VCI(1,I,K) = THETA(I,K)
          VCI(2,I,K) = QVS(I,K)
          ! -- Also mix cloud water and ice IF necessary
          ! IF (IMOISTX.NE.1.AND.IMOISTX.NE.3) THEN  !!! Check other PBL models
          VCI(3,I,K) = QCS(I,K)
          VCI(4,I,K) = QIS(I,K)
#if (WRF_CHEM == 1)
          DO L = NSP+1, NSPX
             VCI(L,I,K) = CHEM(I,K,L-NSP)
          ENDDO
#endif
        ENDDO
      ENDDO

      DO I = its,ILX
        FS(1,I) = -UST(I) * TST(I)
        FS(2,I) = -UST(I) * QST(I)
        FS(3,I) = 0.0
        FS(4,I) = 0.0                      ! SURFACE FLUXES OF CLOUD WATER AND ICE = 0
#if (WRF_CHEM == 1)
        DO L = NSP+1, NSPX
          FS(L,I) = -VD(I,1,L-NSP) * CHEM(I,1,L-NSP)
        ENDDO
#endif
      ENDDO

      !-- UTK urban
      DO I = its,ILX
        IF (ispxurb .AND. URBF(I) .GT. 0.0) THEN
          FS_RF(1,I) = TMPFX_RF(I) !* DENSX(I,CNP) * PSTARI(I)
          FS_RF(2,I) = MOISFX_RF(I) !* DENSX(I,CNP) * PSTARI(I)
          FS_RF(3,I) = 0.0
          FS_RF(4,I) = 0.0               ! ROOF FLUXES OF CLOUD WATER AND ICE = 0

          FS_WL(1,I) = TMPFX_WL(I) !* DENSX_CAN(I) * PSTARI(I)
          FS_WL(2,I) = 0.0
          FS_WL(3,I) = 0.0
          FS_WL(4,I) = 0.0               ! WALL FLUXES OF MOISTURE, CLOUD WATER AND ICE = 0 
        ENDIF
      ENDDO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I = its,ILX      
        CNP = KSL(I)
        NLP   = INT(DTPBL / DTLIM(I) + 1.0)
        DTS   = (DTPBL / NLP)
        DTRAT = DTS / DTPBL

        DO NL = 1,NLP           ! LOOP OVER SUB TIME LOOP              

!-- COMPUTE ARRAY ELEMENTS THAT ARE INDEPENDANT OF SPECIES
          DO K = kts,kte
            AI(K) = 0.0
            BI(K) = 0.0
            CI(K) = 0.0
            EI(K) = 0.0
          ENDDO

          DO K = 2, KCBL(I)
            EI(K-1) = -CRANKP * MDWN(K,I) * DTS * DZH(i,K) * DZHI(i,K-1) !* NFRF(k-1) / VRF(k-1)
            BI(K)   = 1.0 + CRANKP * MDWN(K,I) * DTS !* NFRF(k-1) / VRF(k)
            AI(K)   = -CRANKP * MBARKS(K,I) * DTS !* NFRF(1) / VRF(k)
          ENDDO

          EI(1) = EI(1) -EDDYZ(I,1) * CRANKP * DZHI(i,1) * DZFI(i,1) * DTS * NFRF(1,I) / VRF(1,I)
          AI(2) = AI(2) -EDDYZ(I,1) * CRANKP * DZHI(i,2) * DZFI(i,1) * DTS * NFRF(1,I) / VRF(2,I)

          DO K =  KCBL(I)+1, KL
            BI(K) = 1.0
          ENDDO

          DO K = 2,KL
            XPLUS(K)  = EDDYZ(I,K) * DZHI(i,K) * DZFI(i,K) * DTS * NFRF(K,I) / VRF(K,I)
            XMINUS(K) = EDDYZ(I,K-1) * DZHI(i,K) * DZFI(i,K-1) * DTS * NFRF(k-1,I) / VRF(K,I)
            CI(K)     = - XMINUS(K) * CRANKP
            EI(K)     = EI(K) - XPLUS(K) * CRANKP
            BI(K)     = BI(K) + XPLUS(K) * CRANKP + XMINUS(K) * CRANKP
          ENDDO

          IF (NOCONV(I) .EQ. 1) THEN
            BI(1) = 1.0 + CRANKP * MBARKS(1,I) * (PBL(I)-ZF(i,1)) * DTS *DZHI(i,1) &
                        + EDDYZ(I,1) * CRANKP * DZHI(i,1) * DZFI(i,1) * DTS * NFRF(1,I) / VRF(1,I)
          ELSE
            BI(1) = 1.0 + EDDYZ(I,1) * CRANKP * DZHI(i,1) * DZFI(i,1) * DTS * NFRF(1,I) / VRF(1,I)
          ENDIF


          DO K = 1,KL
            DO L = 1,NSPX                    
              DI(L,K) = 0.0
            ENDDO
          ENDDO
!
!**   COMPUTE TENDENCY OF CBL CONCENTRATIONS - SEMI-IMPLICIT SOLUTION
          DO K = 2,KCBL(I)
            DO L = 1,NSPX                    
              DELC = DTS * (MBARKS(K,I) * VCI(L,I,1) - MDWN(K,I) *       &
                 VCI(L,I,K)  + DZH(i,K+1) * DZHI(i,K) *                      &
                        MDWN(K+1,I) * VCI(L,I,K+1))    !* NFRF(k) / VRF(k)
              DI(L,K)   = VCI(L,I,K) + (1.0 - CRANKP) * DELC
            ENDDO
          ENDDO

          DO K = KCBL(I)+1, KL
            DO L = 1,NSPX                    
              DI(L,K) = VCI(L,I,K)
            ENDDO
          ENDDO

          DO K = 2,KL
            IF (K .EQ. KL) THEN
              DO L = 1,NSPX                    
                DI(L,K) = DI(L,K)  - (1.0 - CRANKP) * XMINUS(K) *           &
                          (VCI(L,I,K) - VCI(L,I,K-1))
              ENDDO
            ELSE
              DO L = 1,NSPX                    
                DI(L,K) = DI(L,K) + (1.0 - CRANKP) * XPLUS(K) *             &
                          (VCI(L,I,K+1) - VCI(L,I,K))  -                    &
                          (1.0 - CRANKP) * XMINUS(K) *                      &
                          (VCI(L,I,K) - VCI(L,I,K-1))
              ENDDO
            ENDIF
          ENDDO

          IF (NOCONV(I) .EQ. 1) THEN
            DO L = 1,NSPX                    
              DI(L,1) = VCI(L,I,1) + (FS(L,I) - (1.0 - CRANKP)              &
                        * (MBARKS(1,I) *                                    &
                        (PBL(I) - ZF(i,1)) * VCI(L,I,1)  -                   &
                        MDWN(2,I) * VCI(L,I,2) * DZH(i,2))) * DZHI(i,1) * DTS
            ENDDO
          ELSE
            DO L = 1,NSPX                    
              DI(L,1) = VCI(L,I,1) + FS(L,I) * DZHI(i,1) * DTS
            ENDDO
          ENDIF
          DO L = 1,NSPX                    
            DI(L,1) = DI(L,1) + (1.0 - CRANKP) * EDDYZ(I,1) * DZHI(i,1)     &
                    * DZFI(i,1) * DTS * (VCI(L,I,2) - VCI(L,I,1)) * NFRF(1,I) / VRF(1,I)
          ENDDO

          !-- UTK Urban --addition of rooftop and wall fluxes
          IF (ispxurb .AND. URBF(I) .GT. 0.0) THEN
            DO L = 1,NSPX
              DI(L,CNP) = DI(L,CNP) - FS_RF(L,I) * DZHI(i,CNP) * DTS
            ENDDO

            DO L = 1,NSPX
              DO K = 1,CNP
                DI(L,K) = DI(L,K) + (FS_WL(L,I) / HGTBDG(I)) * DTS
              ENDDO
            ENDDO  
          ENDIF

          IF ( NOCONV(I) .EQ. 1 ) THEN
            CALL MATRIX (AI, BI, CI, DI, EI, UI, KL, NSPX)
          ELSE
            CALL TRI (CI, BI, EI, DI, UI, KL, NSPX)
          END IF

!
!-- COMPUTE NEW THETAV AND Q
          DO K = 1,KL
            DO L = 1,NSPX                    
              VCI(L,I,K) = UI(L,K)
            ENDDO
          ENDDO

        ENDDO                   ! END SUB TIME LOOP
      ENDDO                     ! END I LOOP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
      DO K = kts,KL
        DO I = its,ILX
          THETAX(I,K) = VCI(1,I,K)
          QVX(I,K)    = VCI(2,I,K)
          QCX(I,K)    = VCI(3,I,K)
          QIX(I,K)    = VCI(4,I,K)
#if (WRF_CHEM == 1)
          DO LL = 7, NSPX
             CHEM(I,K,LL-NSP) = VCI(LL,I,K) 
          ENDDO
#endif
      ENDDO
      ENDDO

      DEALLOCATE (DI)       
      DEALLOCATE (UI)  
      DEALLOCATE (FS)
      DEALLOCATE (VCI)
      DEALLOCATE (FS_RF)
      DEALLOCATE (FS_WL)

   END SUBROUTINE ACM
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-------------------------------------------------------------------          
   SUBROUTINE ACMM (DTPBL, PSTAR,  NOCONV, ZF, DZH, DZHI, JX,      &
                   KLPBL, PBL, DZFI, MOL, UST,                     &
                   TST, QST,  USTM, EDDYZM, DENSX,                 &
                   US,    VS,                                      &
                   UX,    VX,                                      &
                   IMU, EXU, IMV, EXV, URBF, XTIME, ispxurb,       &    !DDY urban
                   ids,ide, jds,jde, kds,kde,                      &
                   ims,ime, jms,jme, kms,kme,                      &
                   its,ite, jts,jte, kts,kte)
!**********************************************************************
!   PBL model called the Asymmetric Convective Model, Version 2 (ACM2) 
!   -- See top of module for summary and references
!
!---- REVISION HISTORY:
!   AX     3/2005 - developed WRF version based on ACM2 in the MM5 PX LSM
!   JP and RG 8/2006 - updates
!
!**********************************************************************
!  ARGUMENTS:
!-- DTPBL           PBL time step
!-- PSTAR           Psurf - Ptop in cb
!-- NOCONV          If free convection =0, no; =1, yes
!-- ZF              Height for full layer
!-- DZH             Layer thickness --> ZF(K) - ZF(K-1)
!-- DZHI            Inverse of layer thickness
!-- JX              N-S index
!-- KLPBL           PBL level at K index
!-- PBL             PBL height in m
!-- DZFI            Inverse layer thickness --> 1/(Z(K+1)-Z(K))
!-- MOL             Monin-Obukhov length in 1D form
!-- UST             U* in 1D form
!-- TST             Theta* in 1D form
!-- QST             Q* in 1D form
!-- USTM            U* for computation of momemtum flux 
!-- EDDYZM          eddy diffusivity for momentum KM
!-- DENSX           dry air density (kg/m^3)
!-- US              U wind 
!-- VS              V wind
!-- THETA           potential temperature
!-- QVS             water vapor mixing ratio (Kg/Kg)
!-- QCS             cloud mixing ratio (Kg/Kg)
!-- QIS             ice mixing ratio (Kg/Kg)
!-- UX              new U wind 
!-- VX              new V wind
!-- THETAX          new potential temperature
!-- QVX             new water vapor mixing ratio (Kg/Kg)
!-- QCX             new cloud mixing ratio (Kg/Kg)
!-- QIX             new ice mixing ratio (Kg/Kg)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      IMPLICIT NONE

!.......Arguments

!... Integer
      INTEGER,  INTENT(IN)      ::      ids,ide, jds,jde, kds,kde, &
                                        ims,ime, jms,jme, kms,kme, &
                                        its,ite, jts,jte, kts,kte, &
                                        JX, XTIME
      INTEGER,  DIMENSION( its:ite ), INTENT(IN)  :: NOCONV
      INTEGER,  DIMENSION( ims:ime ), INTENT(IN)  :: KLPBL

!... Real
      REAL, DIMENSION( ims:ime ),          INTENT(IN)    :: PBL, UST
      REAL,                                INTENT(IN)    :: DTPBL
      REAL, DIMENSION( its:ite ),          INTENT(IN)    :: PSTAR,    &
                                                            MOL, TST, &
                                                            QST, USTM
      REAL, DIMENSION( its:ite, kts:kte ), INTENT(IN)    :: DZHI, DZH, DZFI
      REAL, DIMENSION( its:ite, 0:kte ),   INTENT(IN)    :: ZF
      REAL, DIMENSION( its:ite, kts:kte ), INTENT(INOUT) :: EDDYZM
      REAL, DIMENSION( ims:ime, kms:kme ), INTENT(IN)    :: US, VS,   &
                                                            DENSX
      REAL, DIMENSION( its:ite, kts:kte ), INTENT(OUT)   :: UX, VX

     !-- DDY urban
      REAL, DIMENSION( its:ite, kts:kte, 1:2 ), INTENT(IN)    :: IMU, EXU, IMV, EXV

!.......Local variables

!... Parameters
      INTEGER, PARAMETER :: NSP   = 2
!
      REAL,    PARAMETER :: XX    = 0.5          ! FACTOR APPLIED TO CONV MIXING TIME STEP
      REAL,    PARAMETER :: KARMAN = 0.4

!... Integer
      INTEGER :: ILX, KL, KLM, I, K, NSPX, NLP, NL, JJ, L, IDR
      INTEGER :: KCBLMX
      INTEGER, DIMENSION( its:ite ) :: KCBL

!... Real
      REAL                                :: MBMAX, HOVL, MEDDY, MBAR
      REAL                                :: EKZ, RZ, FM, WSPD, DTS, DTRAT, F1
      REAL, DIMENSION( its:ite )          :: PSTARI, FSACM, DTLIM
      REAL, DIMENSION( kts:kte, its:ite ) :: MBARKS, MDWN
      REAL, DIMENSION( 1:NSP, its:ite )   :: FS
      REAL, DIMENSION( kts:kte )          :: XPLUS, XMINUS
      REAL                                :: DELC
      REAL, DIMENSION( 1:NSP,its:ite,kts:kte  ) :: VCI

      REAL, DIMENSION( kts:kte )                :: AI, BI, CI, EI !, Y
      REAL, DIMENSION( 1:NSP,kts:kte )          :: DI, UI    

     !-- DDY urban
      REAL, DIMENSION( kts:kte )                :: BI2      
      REAL, DIMENSION( 1:1,kts:kte )            :: DI2, UI2
      REAL, DIMENSION( its:ite ), INTENT(IN)    :: URBF
      INTEGER                                   :: NSP2  
      logical, intent(in)                       :: ispxurb
      REAL                                      :: CO

!
!--Start Executable ----

      NSP2 = 1 !-- DDY
      ILX = ite
      KL  = kte
      KLM = kte - 1

      KCBLMX = 0
      MBMAX  = 0.0

!---COMPUTE ACM MIXING RATE
      DO I = its, ILX
        DTLIM(I)  = DTPBL
        PSTARI(I) = 1.0 / PSTAR(I)
        KCBL(I)   = 1
        FSACM(I)  = 0.0

        IF (NOCONV(I) .EQ. 1) THEN
          KCBL(I) = KLPBL(I)

!-------MBARKS IS UPWARD MIXING RATE; MDWN IS DOWNWARD MIXING RATE
!--New couple ACM & EDDY-------------------------------------------------------------
          HOVL     = -PBL(I) / MOL(I)
          FSACM(I) = 1./(1.+((KARMAN/(HOVL))**0.3333)/(0.72*KARMAN))
          MEDDY    = EDDYZM(I,1) * DZFI(i,1) / (PBL(I) - ZF(i,1))
          MBAR     = MEDDY * FSACM(I)
          DO K = kts,KCBL(I)-1
            EDDYZM(I,K) = EDDYZM(I,K) * (1.0 - FSACM(I))
          ENDDO

          MBMAX = AMAX1(MBMAX,MBAR)
          DO K = kts+1,KCBL(I)
            MBARKS(K,I) = MBAR
            MDWN(K,I)   = MBAR * (PBL(I) - ZF(i,K-1)) * DZHI(i,K)
          ENDDO
          MBARKS(1,I) = MBAR
          MBARKS(KCBL(I),I) = MDWN(KCBL(I),I)
          MDWN(KCBL(I)+1,I) = 0.0
        ENDIF
      ENDDO                              ! end of I loop

      DO K = kts,KLM
        DO I = its,ILX
          EKZ   = EDDYZM(I,K) * DZFI(i,K) * DZHI(i,K)
          DTLIM(I) = AMIN1(0.75 / EKZ,DTLIM(I))
        ENDDO
      ENDDO
       
      DO I = its,ILX 
        IF (NOCONV(I) .EQ. 1) THEN
          KCBLMX = AMAX0(KLPBL(I),KCBLMX)
          RZ     = (ZF(i,KCBL(I)) - ZF(i,1)) * DZHI(i,1)
          DTLIM(I)  = AMIN1(XX / (MBARKS(1,I) * RZ),DTLIM(I))
        ENDIF
      ENDDO

      DO K = kts,KL
        DO I = its,ILX
          VCI(1,I,K) = US(I,K)
          VCI(2,I,K) = VS(I,K)
        ENDDO
      ENDDO

      NSPX = 2

      DO I = its,ILX
        WSPD    = SQRT(US(I,1) * US(I,1) + VS(I,1) * VS(I,1)) + 1.E-9
        FM      = -USTM(I) * USTM(I) / WSPD
        IF (URBF(I) .GT. 0.0 .and. ispxurb) THEN
          FM = FM * (1 - URBF(I))
        ENDIF
!        FM = FM * DENSX(I,1) * PSTARI(I)
        FS(1,I) = FM * US(I,1) !/ WSPD
        FS(2,I) = FM * VS(I,1) !/ WSPD
      ENDDO

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I = its,ILX      

        NLP   = INT(DTPBL / DTLIM(I) + 1.0)
        DTS   = (DTPBL / NLP)
        DTRAT = DTS / DTPBL

        DO NL = 1,NLP           ! LOOP OVER SUB TIME LOOP              

!-- COMPUTE ARRAY ELEMENTS THAT ARE INDEPENDANT OF SPECIES
          DO K = kts,KL
            AI(K) = 0.0
            BI(K) = 0.0
            CI(K) = 0.0
            EI(K) = 0.0
          ENDDO

          DO K = 2, KCBL(I)
            EI(K-1) = -CRANKP * MDWN(K,I) * DTS * DZH(i,K) * DZHI(i,K-1)
            BI(K)   = 1.0 + CRANKP * MDWN(K,I) * DTS
            AI(K)   = -CRANKP * MBARKS(K,I) * DTS
          ENDDO

          EI(1) = EI(1) -EDDYZM(I,1) * CRANKP * DZHI(i,1) * DZFI(i,1) * DTS
          AI(2) = AI(2) -EDDYZM(I,1) * CRANKP * DZHI(i,2) * DZFI(i,1) * DTS

          DO K =  KCBL(I)+1, KL
            BI(K) = 1.0
          ENDDO

          DO K = 2,KL
            XPLUS(K)  = EDDYZM(I,K) * DZHI(i,K) * DZFI(i,K) * DTS
            XMINUS(K) = EDDYZM(I,K-1) * DZHI(i,K) * DZFI(i,K-1) * DTS
            CI(K)     = - XMINUS(K) * CRANKP
            EI(K)     = EI(K) - XPLUS(K) * CRANKP
            BI(K)     = BI(K) + XPLUS(K) * CRANKP + XMINUS(K) * CRANKP
          ENDDO

          IF (NOCONV(I) .EQ. 1) THEN
            BI(1) = 1.0 + CRANKP * MBARKS(1,I) * (PBL(I) - ZF(i,1)) * DTS   &
                  * DZHI(i,1) + EDDYZM(I,1) * DZHI(i,1) * DZFI(i,1) * CRANKP * DTS
          ELSE
            BI(1) = 1.0  + EDDYZM(I,1) * DZHI(i,1) * DZFI(i,1) * CRANKP * DTS
          ENDIF


          DO K = 1,KL
            DO L = 1,NSPX                    
              DI(L,K) = 0.0
            ENDDO
          ENDDO
!
!**   COMPUTE TENDENCY OF CBL CONCENTRATIONS - SEMI-IMPLICIT SOLUTION
          DO K = 2,KCBL(I)
            DO L = 1,NSPX                    
              DELC = DTS * (MBARKS(K,I) * VCI(L,I,1) - MDWN(K,I) *          &
                 VCI(L,I,K) + DZH(i,K+1) * DZHI(i,K) *                      &
                        MDWN(K+1,I) * VCI(L,I,K+1))
              DI(L,K)   = VCI(L,I,K) + (1.0 - CRANKP) * DELC
            ENDDO
          ENDDO

          DO K = KCBL(I)+1, KL
            DO L = 1,NSPX                    
              DI(L,K) = VCI(L,I,K)
            ENDDO
          ENDDO

          DO K = 2,KL
            IF (K .EQ. KL) THEN
              DO L = 1,NSPX                    
                DI(L,K) = DI(L,K)  - (1.0 - CRANKP) * XMINUS(K) *           &
                          (VCI(L,I,K) - VCI(L,I,K-1))
              ENDDO
            ELSE
              DO L = 1,NSPX                    
                DI(L,K) = DI(L,K) + (1.0 - CRANKP) * XPLUS(K) *             &
                          (VCI(L,I,K+1) - VCI(L,I,K))  -                    &
                          (1.0 - CRANKP) * XMINUS(K) *                      &
                          (VCI(L,I,K) - VCI(L,I,K-1))
              ENDDO
            ENDIF
          ENDDO

          IF (NOCONV(I) .EQ. 1) THEN
            DO L = 1,NSPX                    
              DI(L,1) = VCI(L,I,1) + (FS(L,I) - (1.0 - CRANKP)              &
                        * (MBARKS(1,I) *                                    &
                        (PBL(I) - ZF(i,1)) * VCI(L,I,1) -                   &
                        MDWN(2,I) * VCI(L,I,2) * DZH(i,2))) * DZHI(i,1) * DTS
            ENDDO
          ELSE
            DO L = 1,NSPX                    
              DI(L,1) = VCI(L,I,1) + FS(L,I) * DZHI(i,1) * DTS
            ENDDO
          ENDIF
          DO L = 1,NSPX                    
            DI(L,1) = DI(L,1) + (1.0 - CRANKP) * EDDYZM(I,1) * DZHI(i,1)    &
                    * DZFI(i,1) * DTS * (VCI(L,I,2) - VCI(L,I,1))
          ENDDO

          !-- DDY urban
          IF (ispxurb .AND. URBF(I) .GT. 0.0) THEN
          DO IDR = 1,2
            DO L = 1,NSPX
              DO K = 1,KL
                CO = DZHI(I,K) * DTS
                IF (L .EQ. 1) THEN
                  BI2(K)   = BI(K) - IMU(I,K,IDR) * CO
                  DI2(1,K) = DI(1,K) + EXU(I,K,IDR) * CO
                ELSE
                  BI2(K)   = BI(K) - IMV(I,K,IDR) * CO
                  DI2(1,K) = DI(2,K) + EXV(I,K,IDR) * CO
                ENDIF
              ENDDO

              IF ( NOCONV(I) .EQ. 1 ) THEN
                CALL MATRIX (AI, BI2, CI, DI2, EI, UI2, KL, NSP2)
              ELSE
                CALL TRI (CI, BI2, EI, DI2, UI2, KL, NSP2)
              ENDIF

              !-- COMPUTE NEW U and V
            if (IDR.eq.1) then
              DO K = 1,KL                
                VCI(L,I,K) = UI2(1,K)
              ENDDO
            else
              DO K = 1,KL                
                VCI(L,I,K) = 0.5 * (VCI(L,I,K) + UI2(1,K))
              ENDDO
            endif
            ENDDO 
          ENDDO   !end 2dir
          ELSE
            IF ( NOCONV(I) .EQ. 1 ) THEN
              CALL MATRIX (AI, BI, CI, DI, EI, UI, KL, NSPX)
            ELSE
              CALL TRI (CI, BI, EI, DI, UI, KL, NSPX)
            ENDIF

!
!-- COMPUTE NEW THETAV AND Q
            DO K = 1,KL
              DO L = 1,NSPX                    
                VCI(L,I,K) = UI(L,K)
              ENDDO
            ENDDO
          ENDIF !-- END DDY URBAN

        ENDDO                   ! END SUB TIME LOOP
      ENDDO                     ! END I LOOP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
      DO K = kts,KL
        DO I = its,ILX
          UX(I,K) = VCI(1,I,K)
          VX(I,K) = VCI(2,I,K)
        ENDDO
      ENDDO

   END SUBROUTINE ACMM
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   SUBROUTINE MATRIX(A,B,C,D,E,X,KL,NSP)
   
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   IMPLICIT NONE
!
!-- Bordered band diagonal matrix solver for ACM2

!-- ACM2 Matrix is in this form:
!   B1 E1
!   A2 B2 E2
!   A3 C3 B3 E3
!   A4    C4 B4 E4
!   A5       C5 B5 E5
!   A6          C6 B6

!--Upper Matrix is
!  U11 U12
!      U22 U23
!          U33 U34
!              U44 U45
!                  U55 U56
!                      U66

!--Lower Matrix is:
!  1
! L21  1
! L31 L32  1
! L41 L42 L43  1
! L51 L52 L53 L54  1
! L61 L62 L63 L64 L65 1
!---------------------------------------------------------
!...Arguments
      INTEGER, INTENT(IN)   :: KL
      INTEGER, INTENT(IN)   :: NSP
      REAL A(KL),B(KL),E(KL)
      REAL C(KL),D(NSP,KL),X(NSP,KL)

!...Locals
      REAL Y(NSP,KL),AIJ,SUM
      REAL L(KL,KL),UII(KL),UIIP1(KL),RUII(KL)
      INTEGER I,J,V

!-- Define Upper and Lower matrices
      L(1,1) = 1.
      UII(1) = B(1)
      RUII(1) = 1./UII(1)
      DO I = 2, KL
        L(I,I) = 1.
        L(I,1) = A(I)/B(1)
        UIIP1(I-1)=E(I-1)
        IF(I.GE.3) THEN
          DO J = 2,I-1
            IF(I.EQ.J+1) THEN
              AIJ = C(I)
            ELSE
              AIJ = 0.
            ENDIF
            L(I,J) = (AIJ-L(I,J-1)*E(J-1))/      &
                     (B(J)-L(J,J-1)*E(J-1))
          ENDDO
        ENDIF
      ENDDO
  
      DO I = 2,KL
        UII(I) = B(I)-L(I,I-1)*E(I-1)
        RUII(I) = 1./UII(I)
      ENDDO
  
!-- Forward sub for Ly=d
      DO V= 1, NSP
        Y(V,1) = D(V,1)
        DO I=2,KL
          SUM = D(V,I)
          DO J=1,I-1
            SUM = SUM - L(I,J)*Y(V,J)
          ENDDO
          Y(V,I) = SUM
        ENDDO
      ENDDO

!-- Back sub for Ux=y

      DO V= 1, NSP
        X(V,KL) = Y(V,KL)*RUII(KL)
      ENDDO
      DO I = KL-1,1,-1
        DO V= 1, NSP
         X(V,I) = (Y(V,I)-UIIP1(I)*X(V,I+1))*RUII(I)
        ENDDO
      ENDDO

   END SUBROUTINE MATRIX


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE TRI(L, D, U, B, X, KL, NSP)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!  FUNCTION:
!    Solves tridiagonal system by Thomas algorithm. 
!   The associated tri-diagonal system is stored in 3 arrays
!   D : diagonal
!   L : sub-diagonal
!   U : super-diagonal
!   B : right hand side function
!   X : return solution from tridiagonal solver

!     [ D(1) U(1) 0    0    0 ...       0     ]
!     [ L(2) D(2) U(2) 0    0 ...       .     ]
!     [ 0    L(3) D(3) U(3) 0 ...       .     ]
!     [ .       .     .     .           .     ] X(i) = B(i)
!     [ .             .     .     .     0     ]
!     [ .                   .     .     .     ]
!     [ 0                           L(n) D(n) ]

!-----------------------------------------------------------------------

      IMPLICIT NONE

! Arguments:

      INTEGER, INTENT(IN)   :: KL
      INTEGER, INTENT(IN)   :: NSP

      REAL        L( KL )               ! subdiagonal
      REAL        D(KL)   ! diagonal
      REAL        U( KL )               ! superdiagonal
      REAL        B(NSP,KL )   ! R.H. side
      REAL        X( NSP,KL )   ! solution

! Local Variables:

      REAL        GAM( KL )
      REAL        BET
      INTEGER     V, K

! Decomposition and forward substitution:
      BET = 1.0 / D( 1 )
      DO V = 1, NSP
         X( V,1 ) = BET * B(V,1 )
      ENDDO

      DO K = 2, KL
        GAM(K ) = BET * U( K-1 )
        BET = 1.0 / ( D( K ) - L( K ) * GAM( K ) )
        DO V = 1, NSP
          X( V, K ) = BET * ( B( V,K ) - L( K ) * X( V,K-1 ) )
        ENDDO
      ENDDO

! Back-substitution:

      DO K = KL - 1, 1, -1
        DO V = 1, NSP
          X( V,K ) = X( V,K ) - GAM( K+1 ) * X( V,K+1 )
        ENDDO
      ENDDO
     
  END SUBROUTINE TRI


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-- DDY urban subroutines
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE ACMURB(COSZ, SWDOWN, CPD,                         &
                        DZF, NLCAT, LANDUSEF, DENSX,               &
                        THETAV, ZA,                                &
                        US, VS, THETA, XTIME, DTPBL, restart,      &
                        MIXLKZ, IMU, EXU, IMV, EXV,                &
                        UST_WL,UST_RF,UST_ST,                      &
                        ispxurb,lamp, HGTBDG, lamf,                &
                        stdir1,stdir2,blkw1,blkw2,                 &
                        ims,ime,kms,kme,jms,jme,its,ite,kts,kte,j)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!--   The subroutine ACMURB is created to account for urban effects 
!--   on wind speed and potential temperature in the ACM2 PBL scheme
!--   
!--   ACMURB subroutine will output, 
!--   (1) MIXLKZ, veritcal mixing length profile, used in EDDYX
!--   (2) IMU, EXU, IMV and EXV, the drag terms on momentum due to urban
!--       used in ACMM
!--   
!--   ACMURB will also account for building surfaces and street temperature
!--   which will be used in calculating the potential temperature drag terms
!--   
!--   REVISION HISTORY: 
!--     DDY    07/2016   first version
!--     UTK    05/2018   second version 
!--
!--   ARGUMENTS LIST:
!--
!--   INPUT VARIABLES:
!--   COSZ      cosine of solar zenith angle
!--   SWDOWN    downward shortwave radiation
!--   ZF        height of each vertical layer
!--   DZF       vertical layer thickness
!--   NLCAT     number of landuse categories
!--   LANDUSEF  landuse fraction
!--   DENSX     density of air at vertical layer height
!--   U, V, THETA         state variables wind speed and potential temperature
!--   ims,ime,kms,kme,its,ite,kts,kte,j   dimensions of arrays
!--
!--
!--   OUTPUT VARIABLES:
!--   MIXLKZ    mixing length profile
!--   IMU       implicit component of U-momentum drag
!--   EXU       explicit component of U-momentum drag
!--   IMV       implicit component of V-momentum drag
!--   EXV       explicit component of V-momentum drag
!--
!--   URBAN PARAMETERS:
!--   LAMF      frontal area density
!--   LAMP      plan area density
!--   HGTBDG    average building height
!--   STANG     street angle
!--   Z0RT      building roof roughness length
!--   Z0ST      street roughness length
!--   ALB       albedo of building surfaces
!--   EMISS     emissivity of building surfaces
!--   SURFDX    building surface thickness
!--
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: XTIME,ims,ime,kms,kme,jms,jme,its,  &
                             ite,kts,kte,j
      LOGICAL, INTENT(IN) :: restart
      REAL, INTENT(IN)    :: DTPBL, CPD
      REAL, DIMENSION( its:ite,kts:kte ), INTENT(IN)  :: THETAV, ZA 
      REAL, DIMENSION( its:ite,kts:kte ), INTENT(OUT) :: MIXLKZ
      REAL, DIMENSION( its:ite,kts:kte, 1:2 ), INTENT(OUT) ::         &
                             IMU,EXU,IMV,EXV
      REAL, DIMENSION( ims:ime,kms:kme ), INTENT(OUT) :: UST_WL
      REAL, DIMENSION( ims:ime ), INTENT(OUT)         :: UST_RF, UST_ST
      REAL, DIMENSION( ims:ime ), INTENT(IN)          :: SWDOWN, COSZ
      REAL, DIMENSION( ims:ime,kms:kme ), INTENT(IN)  ::         &
                             DZF,DENSX,US,VS,THETA
      INTEGER, INTENT(IN) :: NLCAT
      REAL, DIMENSION( ims:ime,1:NLCAT ), INTENT(IN)  :: LANDUSEF

      
      REAL, PARAMETER :: PI = 3.141592653
      !-- urban morphological parameters
      REAL, DIMENSION( ims:ime ), INTENT(IN) :: LAMF   ! frontal area density
      REAL, DIMENSION( ims:ime ), INTENT(IN) :: LAMP   ! plan area density
      REAL, DIMENSION( ims:ime ), INTENT(IN) :: HGTBDG ! average building height
      INTEGER, DIMENSION( ims:ime ), INTENT(IN) :: STDIR1
      INTEGER, DIMENSION( ims:ime ), INTENT(IN) :: STDIR2
      REAL, DIMENSION( ims:ime ), INTENT(IN) :: BLKW1
      REAL, DIMENSION( ims:ime ), INTENT(IN) :: BLKW2
      REAL, DIMENSION( its:ite ) :: STANG  ! street angle (true bearing [deg])

      !-- physical properties
      REAL, PARAMETER :: ALB = 0.2, EMISS = 0.9, SURFDX = 0.25,     &
                         CPCON = 850*2400*SURFDX, KDIFF = 6.7E-5
      !-- constants
      REAL, PARAMETER :: KARMAN = 0.4, SBCON = 5.67E-8
      REAL, PARAMETER :: Z0RT   = 0.05        ! building roof roughness length in meters
      REAL, PARAMETER :: Z0ST   = 0.3         ! street roughness length in meters
      !-- constants
      REAL, PARAMETER :: B1 = 9.4,  B2 = B1/2, CSTM = 7.4,  CSTH = 5.3 
      REAL, PARAMETER :: CN = 0.5
      INTEGER, PARAMETER :: KWAT = 17
      real fwat,RBPW1,RBPW2, BB1, BB2

      INTEGER :: NLEV, I, K, K2, ZCU1, ZCU2, NL, NSP, KL
      REAL    :: DZBDG, TW, ACOEF, LCAN, ZW, ZC, LM, LUTOT, LMU
      REAL    :: GAMDEN, GAMDIS, MORT1, WSP1, CF
      REAL    :: FMHU1, FMHV1, FMV1, FMV2, FM
      REAL    :: CO1, CO2, AA, C, BB, CO3, FMHU2, FMHV2, FM2
      REAL    :: SHAZCU, SHAZ, ZR, SW
      REAL    :: RADRT, RADST, BBDG, WSTT,W, RIN, WSPD
      LOGICAL :: ISMIXL, ISDRAGM
      REAL, DIMENSION( kts:kte ) :: VL, RADZ
      REAL, DIMENSION( 0:kte )   :: Z1, Z2
      REAL, DIMENSION(kts:kte+2) ::& ! D, L, U, DIF,               &
                    SWF, LWF, HTF,TAIR
      REAL, DIMENSION( kts:kte+2 ) :: F,Y
      REAL, DIMENSION( its:ite )    :: URBF
      INTEGER, DIMENSION( its:ite ) :: KSL
      REAL :: RDZBDG, RZWMC, TANZ
      REAL, SAVE ::  RZ0ST, RZ0RT
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      LOGICAL :: eof
      logical,intent(in) :: ispxurb
      INTEGER :: stat, c2, r
      REAL, SAVE :: h(1:350,1:400), lf(1:350,1:400), lp(1:350,1:400), st(1:350,1:400)
      REAL :: data(5)   ! *IMPORTANT NOTE*:- Change the array dim. accord. to dim. of data 
!      REAL, SAVE :: H2, LF2, LP2, STANG2
!      REAL, DIMENSION( ims:ime ) :: STANG2  ! street angle (true bearing [rad])
      REAL :: STANGX
      INTEGER :: IDIR,ihr
! testing
      real UST_ST_U,UST_ST_V,UST_RF_U,UST_RF_V,UST_WL_U,UST_WL_V
      real time,USTWL2
      LOGICAL :: LCZRUN = .TRUE.
!
!------------------Executables start here-----------------------
!---------------------------------------------------------------


      RZ0ST   = 1/Z0ST
      RZ0RT   = 1/Z0RT

      ISMIXL  = .FALSE.
      ISDRAGM = .TRUE.

      !-- initialization of urban momentum drag terms
      Z2  = 0.0
      IMU = 0.0
      IMV = 0.0
      EXU = 0.0
      EXV = 0.0

      DO I = its,ite
        DO K = kts,kte
          Z2(K)       = Z2(K-1) + DZF(I,K)
          MIXLKZ(I,K) = KARMAN * Z2(K)
        ENDDO
      ENDDO

      if (ispxurb) then
        UST_WL =0.0
        DO I = its,ite

        !-DDY urban

          STANG(I)  = 0.0

          !--DDY urban calculate urban fractions for USGS-landuse only (NOW)
          URBF(I) = 0.0
          LUTOT   = 0.0
          DO NL = 1,  40   !NLCAT
            LUTOT = LUTOT + LANDUSEF(I,NL)
          ENDDO
          IF (LUTOT .GT. 0.0) THEN
!            URBF(I) = LANDUSEF(I,1)/LUTOT
            URBF(I) = (LANDUSEF(I,24)+LANDUSEF(I,25)+LANDUSEF(I,26))/LUTOT  !-Only for NLCD40 here. 
            FWAT = LANDUSEF(I,KWAT)
          ENDIF
          IF(URBF(I).LT.0.05.OR.FWAT.GT.0.5) URBF(I)=-1.0
          IF( NLCAT.EQ.40) then 
            LCZRUN = .FALSE.
          ELSE IF ( NLCAT.EQ.51) then
            LCZRUN = .TRUE.
          ELSE
            STOP ' Must be NLCD40 or NLCD 51'
          Endif
          !-- Estimation of the Urban canopy height level (Model-level)
          IF (ISPXURB .AND. (URBF(I) .GT. 0.0)) THEN   !temp hack for testing
            !-- convert urban morphology input parameters into building width and street width
!            BBDG   = 60.  !(LAMP(I) / LAMF(I)) * HGTBDG(I)
            WSTT   = 20.  !(1.0 / SQRT(LAMP(I)) - 1.0) * BBDG
            IF(LCZRUN) then
              BB1 = blkw1(i) - WSTT
              BB2 = blkw2(i) - WSTT
            else
              BB1 = max(blkw1(i)*(1-fwat) - WSTT,10.)
              BB2 = max(blkw2(i)*(1-fwat) - WSTT,10.)
            endif
            RBPW1   = 1.0 / (BB1 + WSTT)
            RBPW2   = 1.0 / (BB2 + WSTT)
!            STANG2(I) = STANG(I)*PI/180.
            !-- calculate each level height                   
!            if(BB1.gt.444..or.BB2.gt.444.) print *, ' BB1, BB2=', BB1,BB2
            Z2(0) = 0.0
            DO K = kts,kte
              Z2(K) = Z2(K-1) + DZF(I,K)
              IF (HGTBDG(I) .LE. Z2(K) .AND. HGTBDG(I) .GE. Z2(K-1)) THEN
                KSL(I) = K
              ENDIF
            ENDDO
            IF (ISMIXL) THEN  
              !-- parameters required for mixing length change
              ACOEF = 9.6 * LAMF(I)
              LCAN  = HGTBDG(I) * ( 1.2 * LAMF(I) * (1 - EXP(-2 * ACOEF) ) &
                      / 4 / ACOEF ** 3) ** 0.5
              ZW    = 2 * HGTBDG(I)
              ZC    = HGTBDG(I)
              RZWMC = 1/(ZW - ZC)
              DO K = kts,kte
                !-- calculate mixing length
                IF (Z2(K) .LT. ZC) THEN
                  LM = LCAN
                ELSEIF (Z2(K) .GT. ZW ) THEN
                  LM = KARMAN * Z2(K)
                ELSE
                  LM = LCAN + (KARMAN * ZW - LCAN) * &
                       (Z2(K) - ZC) * RZWMC
                ENDIF
                !-- calculate weighted mixing length with urban fraction
                MIXLKZ(I,K) = KARMAN * Z2(K) + URBF(I) * &
                              (LM  - KARMAN * Z2(K))
              ENDDO
            ENDIF

            DO K = kts,KSL(I)
              IF(ISDRAGM) THEN
                GAMDEN = 0.
                GAMDIS = 0.
                IF (K .EQ. KSL(I)) THEN 
                  !! get building distribution parameters
                  !! at building top
                  GAMDEN = 1.
                  GAMDIS = (HGTBDG(I) - Z2(K-1)) / (Z2(K) - Z2(K-1))
                ELSE
                  !! lies between building blocks
                  GAMDEN = 0.
                  GAMDIS = 1.
                ENDIF
              ENDIF

              IF (ISDRAGM) THEN
                IF (GAMDIS .GT. 0.0) THEN 
                  WSP1  = MAX(SQRT(US(I,K) * US(I,K) + VS(I,K) * VS(I,K) ),0.01) 

                  AA    = THETAV(I,K+1) - THETAV(I,K)
                  BB    = 2 * 9.81 / (THETAV(I,K+1) + THETAV(I,K))
                  WSPD  = (US(I,K+1)-US(I,K))*(US(I,K+1)-US(I,K)) &
                          +(VS(I,K+1)-VS(I,K))*(VS(I,K+1)-VS(I,K))
                  WSPD  = MAX(WSPD,0.1)
                  RIN   = ABS(ZA(I,K+1)-ZA(I,K)) * AA *BB / WSPD

                  IF (RIN .GE. 0.0) THEN
                    !  stable condition
                    AA  = 1 + B2 * RIN
                    AA  = AA * AA
                    FM  = 1 / AA
                    FM2 = FM
                  ELSE
                    ! unstable condition
                    IF (K .EQ. kts) THEN
                      AA = DZF(I,K) * 0.5 * RZ0ST
                      BB = KARMAN / LOG(AA)
                      CF = BB * BB * B1 * SQRT(AA)
                      AA = B1 * RIN
                      BB = CF * SQRT(-RIN)
                      FM = 1 - AA / (1 + CSTM * BB)
                    ENDIF
                    IF (K .EQ. KSL(I)) THEN
                      AA = DZF(I,K) * 0.5 * RZ0RT
                      BB = KARMAN / LOG(AA)
                      CF = BB * BB * B1 * SQRT(AA)
                      AA = B1 * RIN
                      BB = CF * SQRT(-RIN)
                      FM2 = 1 - AA / (1 + CSTM * BB)
                    ENDIF
                  ENDIF
                  FMHU1 = 0.0
                  FMHV1 = 0.0
                  FMHU2 = 0.0
                  FMHV2 = 0.0
                  IF (K .EQ. kts) THEN
                    AA    = KARMAN / LOG(DZF(I,K) * 0.5 * RZ0ST)
                    AA    = - AA * AA * WSP1
                    FMHU1 = AA * FM 
                    FMHV1 = FMHU1
                  ENDIF
                  IF (K .EQ. KSL(I)) THEN 
                    AA    = KARMAN / LOG(DZF(I,K) * 0.5 * RZ0RT)
                    AA    = - AA * AA * WSP1
                    FMHU2 = AA * FM2 
                    FMHV2 = FMHU2
                  ENDIF

!                  VL(K) = 1. - GAMDIS * BBDG * RBPW
!                  AA    = URBF(I) / VL(K) * RBPW
!                  IF (VL(K) .NE. 0.0) AA    = URBF(I) / VL(K) * RBPW  ! temp hack
!                  CO1   = WSTT * AA
                  CO1   = URBF(I) * (1. - LAMP(I))/(1-gamdis*LAMP(I))
                  CO3   = URBF(I) * LAMP(I)

                USTWL2 = 0.0
                DO idir = 1,2  ! 2 street dir
!                  CO2   = DZF(I,K) * GAMDIS * AA
!                  CO2   = DZF(I,K) * GAMDIS * URBF(I) / ((WSTT+BBDG) - BBDG**2 * RBPW)
!                  CO2   = 2. * DZF(I,K) * GAMDIS * URBF(I) * (1. - WSTT*RBPW1)*RBPW2
!                  CO3   = BBDG * AA
                IF(idir.eq.1) Then
                   STANGX = float(stdir1(i))
                   CO2   = DZF(I,K) * GAMDIS * URBF(I) * (1. - WSTT*RBPW2)*RBPW1/(1. - LAMP(I))
                Else
                   STANGX = float(stdir2(i))
                   CO2   = DZF(I,K) * GAMDIS * URBF(I) * (1. - WSTT*RBPW1)*RBPW2/(1. - LAMP(I))
                EndIf
                STANGX = STANGx*PI/180.
                  MORT1 = ABS(COS(STANGX) * US(I,K) - SIN(STANGX) * VS(I,K))
                  FMV1  = - 0.4 * MORT1 
                  FMV2  = - 0.4 * MORT1 
                  IMU(I,K,idir) = FMV1 * COS(STANGX) * COS(STANGX) * CO2
                  EXU(I,K,idir) = -FMV1 * SIN(STANGX) * COS(STANGX) * VS(I,K) * CO2
                  EXV(I,K,idir) = -FMV2 * SIN(STANGX) * COS(STANGX) * US(I,K) * CO2
                  IMV(I,K,idir) = FMV2 * SIN(STANGX) * SIN(STANGX) * CO2
!                  if(HGTBDG(I).gt.100.0.and.k.eq.2) print *,'HGTBDG,STANGX,EXU,EXV=',HGTBDG(I),STANGX,EXU(I,K),EXV(I,K)
                  IF (K .EQ. kts) THEN
                    IMU(I,K,idir) = IMU(I,K,idir) + FMHU1 * CO1
                    IMV(I,K,idir) = IMV(I,K,idir) + FMHV1 * CO1
                    IF(idir.eq.1) Then
                      UST_ST_U = abs(FMHU1* US(I,KTS)* CO1) 
                      UST_ST_V = abs(FMHV1* VS(I,KTS)* CO1)
                    ENDIF

                  ENDIF

                  IF (K .EQ. KSL(I)) THEN
                    IMU(I,K,idir) = IMU(I,K,idir) + FMHU2 * CO3
                    IMV(I,K,idir) = IMV(I,K,idir) + FMHV2 * CO3
                    IF(idir.eq.1) Then
                      UST_RF_U = abs(FMHU2* US(I,KSL(i))* CO3)
                      UST_RF_V = abs(FMHV2* VS(I,KSL(i))* CO3)
                    ENDIF

                  ENDIF

                  IF (XTIME.EQ.1) THEN
                    IMU(I,K,idir) = 0.
                    EXU(I,K,idir) = 0.
                    EXV(I,K,idir) = 0.
                    IMV(I,K,idir) = 0.
                  ENDIF
                  
                  USTWL2 = USTWL2 + 0.4*MORT1*MORT1 * CO2  !  sqrt(k Uort**2 fwl))
                ENDDO ! 2 street dir
                UST_WL(i,k) = sqrt(USTWL2)
                ENDIF !-- end if gamdis

              ENDIF !-- end if drag flags
            ENDDO !-- end K loop

!             DO K = kts,KSL(I)
!               UST_WL_U=0.
!               UST_WL_v=0.
!               DO idir = 1,2
!                 IF (K .EQ. kts) THEN
!                    UST_WL_U=UST_WL_U+abs((IMU(I,K,idir)-FMHU1 * CO1)*US(I,k))+abs(EXU(I,K,idir))
!                    UST_WL_V=UST_WL_V+abs((IMV(I,K,idir)-FMHV1 * CO1)*VS(I,k))+abs(EXV(I,K,idir))
!                 else
!                 UST_WL_U=UST_WL_U+abs(IMU(I,K,idir)*US(I,k))+abs(EXU(I,K,idir))
!                 UST_WL_V=UST_WL_V+abs(IMV(I,K,idir)*VS(I,k))+abs(EXV(I,K,idir))
!                 endif
!               ENDDO
!               UST_WL(i,k) = sqrt(UST_WL_U + UST_WL_V)
!             ENDDO
!            USTOU_WL(I) = sqrt(UST_WL_U**2 + UST_WL_V**2)/HGTBDG(I)
            UST_RF(I)   = sqrt(UST_RF_U + UST_RF_V)
            UST_ST(I)   = sqrt(UST_ST_U + UST_ST_V)
            time = xtime*dtpbl  !s
            ihr = mod(int(time),3600)
!            if(HGTBDG(I).gt.120.0.and.ihr.eq.0) then
!             print *,'HGTBDG,time,ksl=',HGTBDG(I),time,ksl(i)
!             print *,'UST_WL_U,UST_WL_V,UST_WL=',UST_WL_U,UST_WL_V,UST_WL(I,2)
!             print *,'UST_ST_U,UST_ST_V=',UST_ST_U,UST_ST_V,US(I,KTS),VS(I,KTS)
!             print *,'UST_RF_U,UST_RF_V,UST_RF=',UST_RF_U,UST_RF_V,UST_RF(I)
!            ENDIF
          ENDIF !-- end if urban fraction .gt. 0.0
        ENDDO !-- end I loop
      endif

      END SUBROUTINE ACMURB

END MODULE module_bl_acm
