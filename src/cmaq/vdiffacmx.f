!--------------------------------------------------------------------!
!  The Community Multiscale Air Quality (CMAQ) system software is in     !
!  continuous development by various groups and is based on information  !
!  from these groups: Federal Government employees, contractors working  !
!  within a United States Government contract, and non-Federal sources   !
!  including research institutions.  These groups give the Government    !
!  permission to use, prepare derivative works of, and distribute copies !
!  of their work in the CMAQ system to the public and to permit others   !
!  to do so.  The United States Environmental Protection Agency          !
!  therefore grants similar permission to use the CMAQ system software,  !
!  but users are requested to provide copies of derivative works or      !
!  products designed to operate in the CMAQ system to the United States  !
!  Government without restrictions as to use by others.  Software        !
!  that is used with the CMAQ system but distributed under the GNU       !
!  General Public License or the GNU Lesser General Public License is    !
!  subject to their copyright restrictions.                              !
!------------------------------------------------------------------------!

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE VDIFFACMX ( DTSEC, SEDDY, DDEP, ICMP, CNGRD )

C-----------------------------------------------------------------------
C Asymmetric Convective Model v2 (ACM2/ACM1) -- Pleim(2006/2014)
C Function:
C   calculates vertical diffusion

C Subroutines and Functions Called:
C   SEC2TIME, TIME2SEC, WRITE3, NEXTIME,
C   M3EXIT, EDDYX, TRI, MATRIX, PA_UPDATE_EMIS, PA_UPDATE_DDEP

C Revision History:
C   Analogous to VDIFFACM2
C 11 Apr 13 J.Young: fix double adjustment of conc for DDBF in heterogeneous HONO
C           if-then-else clauses; eliminate some white space
C 13 May 13 J.Young: access met data from VDIFF_MET module
C                    change CRANKP to THBAR, CRANKQ to THETA
C 25 May 13 J.Young: re-do the acm/eddy algorithm for computational efficiency
C 30 Apr 14 J.Young: switch THBAR and THETA
C  2 May 14 J.Pleim, J.Young: replace the banded tridiagonal matrix solver for the
C           convective PBL, with the ACM1 matrix solver followed by the tridiagonal
C           matrix solver
C   30 May 14 J.Young: split vdiff calculation out of vdiff proc.
C   07 Nov 14 J.Bash: Updated for the ASX_DATA_MOD shared data module. 
C   02 Nov 2018: L.Zhou, S.Napelenok: isam implementation
C   May 2019  J.Pleim Changed from sigma coords to Z coords for compatability w/ MPAS and WRF
C   12 Dec 19 S.L.Napelenok: ddm-3d implementation for version 5.3.1
C   15 Jun 21 J. Pleim: implemented HONO fix for dry depsotion flux
C-----------------------------------------------------------------------

      USE CGRID_SPCS          ! CGRID mechanism species
      USE GRID_CONF
      USE DESID_VARS, ONLY : VDEMIS_DIFF,DESID_LAYS
      USE DESID_PARAM_MODULE, ONLY : DESID_N_SRM
      USE DEPV_DEFN
      USE ASX_DATA_MOD
      USE VDIFF_MAP
      USE UTILIO_DEFN
!      USE BIDI_MOD
!      USE LSM_MOD, ONLY: N_LUFRAC
      USE VDIFF_DIAG, NLPCR => NLPCR_MEAN
      USE HGRD_DEFN,only : COLSX_PE, ROWSX_PE
      USE BDSNP_MOD, ONLY: GET_N_DEP



      IMPLICIT NONE

!      INCLUDE SUBST_FILES_ID  ! file name parameters

      CHARACTER( 120 ) :: XMSG = ' '

C Arguments:
      REAL, INTENT( IN )    :: DTSEC                ! model time step in seconds
C--- SEDDY is strictly an input, but it gets modified here
      REAL, INTENT( INOUT ) :: SEDDY    ( :,:,: )   ! flipped EDDYV
      REAL, INTENT( INOUT ) :: DDEP     ( :,:,: )   ! ddep accumulator
      REAL, INTENT( INOUT ) :: ICMP     ( :,:,: )   ! component flux accumlator 
      REAL, INTENT( INOUT ) :: CNGRD    ( :,:,:,: ) ! cgrid replacement

C Parameters:

C explicit, THETA = 0, implicit, THETA = 1     ! Crank-Nicholson: THETA = 0.5
      REAL, PARAMETER :: THETA = 0.5,
     &                   THBAR = 1.0 - THETA

      REAL, PARAMETER :: EPS = 1.0E-06

C External Functions: None

C Local Variables:

      CHARACTER( 16 ), SAVE :: PNAME = 'VDIFFACMX'

      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      LOGICAL, SAVE :: SPECLOG = .TRUE.             ! For BDSNP
      LOGICAL, SAVE :: ispxurb = .TRUE.
!      LOGICAL, SAVE :: ispxurb = .FALSE.
      REAL, ALLOCATABLE, SAVE :: DD_FAC     ( : )   ! combined subexpression
      REAL, ALLOCATABLE, SAVE :: DDBF       ( : )   ! secondary DDEP
      REAl, ALLOCATABLE, SAVE :: CMPF       ( : )   ! intermediate CMP
      REAL, ALLOCATABLE, SAVE :: CONC       ( :,: ) ! secondary CGRID_DATA expression
      REAL, ALLOCATABLE, SAVE :: EMIS       ( :,: ) ! emissions subexpression
      REAL        DTDENS1                       ! DT * layer 1 air density
      REAL URBF, RFFRC
      REAL     :: VRF   ( NLAYS ) 
      REAL     :: NFRF  ( NLAYS )    

C ACM Local Variables
      REAL     :: EDDY  ( NLAYS )               ! local converted eddyv
      REAL        MEDDY                         ! ACM2 intermediate var
      REAL        MBAR                          ! ACM2 mixing rate (S-1)
      REAL     :: MBARKS( NLAYS )               ! by layer
      REAL     :: MDWN  ( NLAYS )               ! ACM down mix rate
      REAL     :: MFAC  ( NLAYS )               ! intermediate loop factor
      REAL     :: AA    ( NLAYS )               ! matrix column one
      REAL     :: BB1   ( NLAYS )               ! diagonal for MATRIX1
      REAL     :: BB2   ( NLAYS )               ! diagonal for TRI
      REAL     :: CC    ( NLAYS )               ! subdiagonal
      REAL     :: EE1   ( NLAYS )               ! superdiagonal for MATRIX1
      REAL     :: EE2   ( NLAYS )               ! superdiagonal for TRI
      REAL, ALLOCATABLE, SAVE :: DD ( :,: )     ! R.H.S
      REAL, ALLOCATABLE, SAVE :: UU ( :,: )     ! returned solution
      REAL        DFACP, DFACQ
      REAL     :: DFSP( NLAYS ), DFSQ( NLAYS )  ! intermediate loop factors
      REAL        DELC, DELP, RP, RQ
      REAL     :: LFAC1( NLAYS )                ! intermediate factor for CONVT
      REAL     :: LFAC2( NLAYS )                ! intermediate factor for CONVT
      REAL     :: LFAC3( NLAYS )                ! intermediate factor for eddy
      REAL     :: LFAC4( NLAYS )                ! intermediate factor for eddy
      REAL, ALLOCATABLE, SAVE :: DEPVCR     ( : )   ! dep vel in one cell
                                                    ! one cell for each landuse category
      REAL, ALLOCATABLE, SAVE :: EFAC1 ( : )
      REAL, ALLOCATABLE, SAVE :: EFAC2 ( : )
      REAL, ALLOCATABLE, SAVE :: EFAC1_wl ( :,: )
      REAL, ALLOCATABLE, SAVE :: EFAC2_wl ( :,: )
      REAL, ALLOCATABLE, SAVE :: EFAC1_rf ( : )
      REAL, ALLOCATABLE, SAVE :: EFAC2_rf ( : )
      REAL, ALLOCATABLE, SAVE :: POL   ( : )    ! prodn/lossrate = PLDV/DEPV
      REAL        PLDV_HONO                     ! PLDV for HONO
      REAL        DEPV_NO2                      ! dep vel of NO2
      REAL        DEPV_HNO3                     ! dep vel of HNO3
      REAL        FNL                           ! ACM2 Variable
      INTEGER     NLP, NL, LCBL
      INTEGER, SAVE :: NO2_HIT = 0, HONO_HIT = 0, HNO3_HIT = 0, NO2_MAP= 0, HONO_MAP = 0, HNO3_MAP = 0
      INTEGER, SAVE :: O3_HIT = 0, O3_MAP = 0
      INTEGER, SAVE :: NH3_HIT = 0
      REAL        DTLIM, DTS, DTACM, RZ

      INTEGER     ASTAT
      INTEGER     C, R, L, S, V, I, J           ! loop induction variables
      INTEGER     MDATE, MTIME                  ! internal simulation date&time
      INTEGER     LAYRF                         ! Layer of roof for urban
      REAL        HGTBDG, CVOL
!--Local Arrays for Z-coord implimentation
      REAL     :: DZH   ( NLAYS )               ! ZF(L) - ZF(L-1)
      REAL     :: DZHI  ( NLAYS )               ! 1/DZH
      REAL     :: DZFI  ( NLAYS )               ! ZH(L+1) - ZH(L)
      integer  gl_c, gl_r


      INTERFACE
         SUBROUTINE MATRIX1 ( KL, A, B, E, D, X )
            INTEGER,        INTENT( IN )  :: KL
            REAL,           INTENT( IN )  :: A( : ), B( : ), E( : )
            REAL,           INTENT( IN )  :: D( :,: )
            REAL,           INTENT( OUT ) :: X( :,: )
         END SUBROUTINE MATRIX1
         SUBROUTINE TRI ( L, D, U, B, X )
            REAL,           INTENT( IN )  :: L( : ), D( : ), U( : )
            REAL,           INTENT( IN )  :: B( :,: )
            REAL,           INTENT( OUT ) :: X( :,: )
         END SUBROUTINE TRI
      END INTERFACE

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN

         FIRSTIME = .FALSE.

         MDATE = 0; MTIME = 0

C set auxiliary depv arrays

         ALLOCATE ( DD_FAC( N_SPC_DEPV  ),
     &              DDBF  ( N_SPC_DEPV ),
     &              DEPVCR( N_SPC_DEPV ),
     &              EFAC1 ( N_SPC_DEPV ),
     &              EFAC2 ( N_SPC_DEPV ),
     &              EFAC1_wl ( N_SPC_DEPV,NLAYS ),
     &              EFAC2_wl ( N_SPC_DEPV,NLAYS ),
     &              EFAC1_rf ( N_SPC_DEPV ),
     &              EFAC2_rf ( N_SPC_DEPV ),
     &              POL   ( N_SPC_DEPV ), STAT = ASTAT )
         IF ( ASTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating DD_FAC, DDBF, DEPVCR, EFAC1, EFAC2, or POL'
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
         END IF

         ALLOCATE ( CMPF( LCMP ), STAT = ASTAT )
         IF ( ASTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating CMPF'
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
         END IF

         ALLOCATE ( CONC( N_SPC_DIFF,NLAYS ),
     &              EMIS( N_SPC_DIFF,NLAYS ), STAT = ASTAT )
         IF ( ASTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating CONC or EMIS'
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
         END IF
         CONC = 0.0; EMIS = 0.0   ! array assignment

         ALLOCATE ( DD( N_SPC_DIFF,NLAYS ),
     &              UU( N_SPC_DIFF,NLAYS ), STAT = ASTAT )
         IF ( ASTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating DD or UU'
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
         END IF
         DD = 0.0; UU = 0.0   ! array assignment


         HONO_HIT = 0; HNO3_HIT = 0; NO2_HIT  = 0; NH3_HIT = 0
         HONO_MAP = 0; HNO3_MAP = 0; NO2_MAP  = 0
         DO V = 1, N_SPC_DEPV
            IF ( DV2DF_SPC( V ) .EQ. 'NO2' ) THEN
               NO2_HIT = V
               NO2_MAP = DV2DF( V )
            ELSE IF ( DV2DF_SPC( V ) .EQ. 'HONO' ) THEN
               HONO_HIT = V
               HONO_MAP = DV2DF( V )
            ELSE IF ( DV2DF_SPC( V ) .EQ. 'HNO3' ) THEN
               HNO3_HIT = V
               HNO3_MAP = DV2DF( V )
            ELSE IF ( DV2DF_SPC( V ) .EQ. 'NH3' ) THEN
               NH3_HIT = V
            ELSE IF ( DV2DF_SPC( V ) .EQ. 'O3' ) THEN
               O3_HIT = V
               O3_MAP = DV2DF( V )
            END IF
         END DO


 
      END IF   !  if Firstime

C ------------------------------------------- Row, Col LOOPS -----------

      DO 345 R = 1, NROWS
      DO 344 C = 1, NCOLS
         DZH(1)  =  Met_Data%ZF( C,R,1 )
         DZHI(1) =  1./DZH(1)
         DO L = 2, NLAYS
            DZH(L)  =  Met_Data%ZF( C,R,L ) - Met_Data%ZF( C,R,L-1 ) 
            DZHI(L) =  1./DZH(L)
         ENDDO
         DO L = 1, NLAYS - 1
            DZFI(L) = 1. / ( Met_Data%ZH( C,R,L+1 ) - Met_Data%ZH( C,R,L ) )
         ENDDO
         DZFI(NLAYS) = DZFI(NLAYS-1)

         URBF = (GRID_DATA%lufrac(c,r,24) + 
     &               GRID_DATA%lufrac(c,r,25) +
     &               GRID_DATA%lufrac(c,r,26))     !-Only for NLCD40 here.         
         IF (ISPXURB .AND. (URBF .GT. 0.001)) then
           HGTBDG =  GRID_DATA%build_height(C,R)
           LAYRF = 1
           DO L = 2, NLAYS - 1 
             IF ((HGTBDG .GE. Met_Data%ZF( C,R,L-1 )).AND.(HGTBDG .LE.Met_Data%ZF( C,R,L ))) THEN 
              LAYRF = L
             ENDIF
           ENDDO
!          if (hgtbdg.gt.200.) print *,'LAYRF,HGTBDG,URBF=',
!     &         layrf,HGTBDG,URBF         
         ENDIF

! Adjustment for volume in urban canopies
        DO L = 1, NLAYS   
           VRF(L) =1.0
           NFRF(L) = 1.0 
        ENDDO

        IF (ispxurb .AND. (URBF .GT. 0.001).and.(HGTBDG.gt.1.0)) THEN 
          RFFRC = URBF*GRID_DATA%build_area_fraction(C,R)
          if(RFFRC.gt.0.9)print *,'RFFRC,URBF,LAMP=',RFFRC,URBF,GRID_DATA%build_area_fraction(C,R)
!          FS(1,I) = (1.-RFFRC)*FS(1,I)
!          FS(2,I) = (1.-RFFRC)*FS(2,I)
          DO L = 1, NLAYS
            IF( L.LT.Layrf) then
               VRF(L) = 1. -RFFRC
               NFRF(L) = 1. -RFFRC
            ELSE
               VRF(L) =1.0
               NFRF(L) = 1.0 
            ENDIF
          ENDDO
          IF (LAYRF.EQ.1) THEN
             VRF(LAYRF) = 1. -HGTBDG*DZHI(LAYRF)*RFFRC
          ELSE
             VRF(LAYRF) = 1. -(HGTBDG-Met_Data%ZF( C,R,LAYRF-1 ))*DZHI(LAYRF)*RFFRC
          ENDIF
!          if (hgtbdg.gt.200.) print *,'VRF,LAYRF,HGTBDG,zf-1,DZHI,RFFRC=',
!     &         VRF(LAYRF),layrf,HGTBDG,Met_Data%ZF( C,R,LAYRF-1 ),DZHI(LAYRF),RFFRC
        ENDIF

C for ACM time step
         DTLIM = DTSEC

C dt = .75 dzf*dzh / Kz
         DO L = 1, NLAYS - 1
            DTLIM = MIN( DTLIM, 0.5 / ( SEDDY( L,C,R ) * DZHI(L)*DZFI(L) ) )
         END DO
         MBARKS = 0.0   ! array assignment
         MDWN = 0.0     ! array assignment

C conjoin ACM & EDDY ---------------------------------------------------

         MBAR = 0.0
         FNL = 0.0

         IF ( Met_Data%CONVCT( C,R ) ) THEN   ! Do ACM for this column
            CVOL = (Met_Data%PBL( C,R ) - Met_Data%ZF(C,R,1))
            LCBL = Met_Data%LPBL( C,R )
            IF (ispxurb .AND. (URBF .GT. 0.001).and.(HGTBDG.gt.1.0)) THEN
               IF(LCBL.gt.LAYRF)THEN 
                 CVOL = CVOL - RFFRC*(HGTBDG-Met_Data%ZF( C,R,1))
               ELSEIF (LCBL.eq.LAYRF) Then
                 CVOL = DZH(LAYRF)*VRF(LAYRF) + (Met_Data%ZF( C,R,LAYRF-1 ) - 
     &                  Met_Data%ZF( C,R,1)) * (1. -RFFRC)
               ELSE
                 CVOL = CVOL * (1. -RFFRC)
               endif
            endif


            MEDDY = SEDDY( 1,C,R ) * DZFI(1) / CVOL/VRF(1)
            FNL = 1.0 / ( 1.0 + ( ( KARMAN / ( -Met_Data%HOL( C,R ) ) ) ** 0.3333 )
     &                / ( 0.72 * KARMAN ) )
            MBAR = MEDDY * FNL
            IF ( MEDDY .LT. EPS ) THEN
               gl_c = c + COLSX_PE(1,mype+1) -1
               gl_r = r + ROWSX_PE(1,mype+1) -1
               WRITE( LOGDEV,* ) ' Warning --- MEDDY < 1e-6 s-1'
               WRITE( LOGDEV,* ) ' SEDDY, MEDDY, FNL, HOL = ',
     &                             SEDDY( 1,C,R ), MEDDY, FNL, Met_Data%HOL( C,R )
               XMSG = '*** ACM fails ***'
               WRITE( LOGDEV,*)' c,r=', gl_c,gl_r,' pbl,ust=',Met_Data%PBL( C,R ),Met_Data%USTAR( C,R )
!               CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT2 )
            END IF

            IF ( ( FNL .LE. 0.0 ) .OR.            ! never gonna happen for CONVCT
     &           ( LCBL .GE. NLAYS-1 ) .OR.       ! .GT. never gonna happen
     &           ( Met_Data%HOL( C,R ) .GT. -0.00001 ) )   ! never gonna happen
     &         WRITE( LOGDEV,1015 ) LCBL, MBAR, FNL, SEDDY( 1,C,R ), Met_Data%HOL( C,R )
1015           FORMAT( ' LCBL, MBAR, FNL, SEDDY1, HOL:', I3, 1X, 4(1PE13.5) )

            DO L = 2, LCBL 
               SEDDY( L,C,R ) = ( 1.0 - FNL ) * SEDDY( L,C,R  )
               MBARKS( L ) = MBAR
               MDWN( L )   = MBAR * (Met_Data%PBL( C,R ) - Met_Data%ZF(C,R,L-1)) * DZHI(L)
            END DO
            SEDDY( 1,C,R ) = ( 1.0 - FNL ) * SEDDY( 1,C,R  )
            MBARKS(1) = MBAR
            MBARKS(LCBL) = MDWN(LCBL)
            MDWN(LCBL+1) = 0.0

C Modify Timestep for ACM2
            RZ     = (Met_Data%ZF(C,R,LCBL) - Met_Data%ZF(C,R,1)) * DZHI(1)
            DTACM  = 1.0 / ( MBAR * RZ )
            DTLIM  = MIN( 0.75 * DTACM, DTLIM )
         ELSE
            LCBL = 1
         END IF

C-----------------------------------------------------------------------
         DTLIM= DTLIM   !/5.

         NLP = INT( DTSEC / DTLIM + 0.99 )
         IF ( VDIFFDIAG ) NLPCR( C,R ) = REAL( NLP )
         DTS = DTSEC / REAL( NLP )
         if(dts.lt.0.0.or.dts.gt.300.)print *,'dts,dtsec,nlp,dtlim,MBAR=',dts,dtsec,nlp,dtlim,MBAR
         DTDENS1 = DTS * Met_Data%DENS1( C,R )
         DFACP = THETA * DTS
         DFACQ = THBAR * DTS


         DO L = 1, NLAYS
            DO V = 1, N_SPC_DIFF
               CONC( V,L ) = CNGRD( DIFF_MAP( V ),L,C,R )
            END DO
         END DO


 
         EMIS = 0.0      ! array assignment
         IF ( DESID_N_SRM .GE. 1 ) 
     &        EMIS( :,1:DESID_LAYS ) = DTS * VDEMIS_DIFF( :,:,C,R )
                
         DO L = 1, NLAYS
            DFSP( L ) = DFACP * DZHI( L )
            DFSQ( L ) = DFACQ * DZHI( L )
            EDDY( L ) = SEDDY( L,C,R ) * DZFI(L)
         END DO
         RP = DFACP * Met_Data%RDEPVHT( C,R )
         RQ = DFACQ * Met_Data%RDEPVHT( C,R )
         DO V = 1, N_SPC_DEPV
            DDBF( V )   = DDEP( V,C,R )
            DEPVCR( V ) = DEPV( V,C,R )
            DD_FAC( V ) = DTDENS1 * DD_CONV( V ) * DEPVCR( V )
            EFAC1 ( V ) = EXP( -DEPVCR( V ) * RP )
            EFAC2 ( V ) = EXP( -DEPVCR( V ) * RQ )
            POL   ( V ) = PLDV( V,C,R ) / DEPVCR( V )
         ENDDO
         EFAC1_wl= 1.0
         EFAC2_wl = 1.0
         EFAC1_rf = 1.0
         EFAC2_rf = 1.0 
         IF (ispxurb .AND. (URBF .GT. 0.001).and.(HGTBDG.gt.1.0)) THEN 
            DO V = 1, N_SPC_DEPV
               DO L=1,NLAYS
                  IF(DEPV_wl( V,C,R,L ).lt.1.e-6) EXIT
                  EFAC1_wl(V,L) = EXP( -DEPV_wl( V,C,R,L ) * DFSP( L ) )
                  EFAC2_wl(V,L) = EXP( -DEPV_wl( V,C,R,L ) * DFSQ( L ) )
!                  if(EFAC1_wl(V,L).lt.0.9.or.EFAC1_wl(V,L).gt.1.0)print *,'EFAC1_wl,depv_wl,HGTBDG,L= ',
!     &            EFAC1_wl(v,L),depv_wl( V,C,R,L ),HGTBDG,L
               ENDDO 
               IF(DEPV_rf( V,C,R ).gt.1.e-6) then
                 EFAC1_rf(V) = EXP( -DEPV_rf( V,C,R ) * DFSP( Layrf ) )
                 EFAC2_rf(V) = EXP( -DEPV_rf( V,C,R ) * DFSQ( Layrf ) )
!                 print *,'V,C,R,EFAC1_rf,depv_rf,HGTBDG= ',V,C,R,
!     &            EFAC1_rf(v),depv_rf( V,C,R ),HGTBDG
               endif
            ENDDO
         endif

         PLDV_HONO = PLDV( HONO_HIT,C,R )


C These don`t change in the NLP sub-time step loop:---------------------
         DO L = 1, NLAYS
            AA ( L ) = 0.0
            BB1( L ) = 0.0
            EE1( L ) = 0.0
            CC ( L ) = 0.0
            EE2( L ) = 0.0
            BB2( L ) = 0.0
         END DO
         IF ( Met_Data%CONVCT( C,R ) ) THEN
            L = 1
            DELP = Met_Data%PBL( C,R ) - Met_Data%ZF( C,R,L )
            BB1( L ) = 1.0 + DELP * DFSP( L ) * MBARKS( L ) !* NFRF(L) / VRF(L)
            LFAC1( L ) = DFSQ( L ) * DELP * MBARKS( L ) !* NFRF(L) / VRF(L)
            LFAC2( L ) = DFSQ( L ) * MDWN( L+1 ) * DZH( L+1 ) !* NFRF(L) / VRF(L)
            DO L = 2, LCBL
               AA ( L ) = -DFACP * MBARKS( L ) !* NFRF(1) / VRF(L)
               BB1( L ) = 1.0 + DFACP * MDWN( L ) !* NFRF(L-1) / VRF(L)
               EE1( L ) = -DFSP( L-1 ) * DZH( L ) * MDWN( L ) !* NFRF(L-1) / VRF(L-1)
               MFAC( L ) = DZH( L+1 ) * DZHI( L ) * MDWN( L+1 )
            END DO
         END IF

         DO L = 1, NLAYS
            EE2( L ) = - DFSP( L ) * EDDY( L ) * NFRF(L) / VRF(L)
            LFAC3( L ) = DFSQ( L ) * EDDY( L ) * NFRF(L) / VRF(L)
         END DO

         BB2( 1 ) = 1.0 - EE2( 1 )
         DO L = 2, NLAYS
            CC ( L ) = - DFSP( L ) * EDDY( L-1 ) * NFRF(l-1) / VRF(l)
            BB2( L ) = 1.0 - CC( L ) - EE2( L ) 
            LFAC4( L ) = DFSQ( L ) * EDDY( L-1 ) * NFRF(l-1) / VRF(l)
         END DO

         DO 301 NL = 1, NLP      ! loop over sub time

            DO V = 1, N_SPC_DEPV

C --------- HET HONO RX -----------------

C Use special treatment for HNO3
C HNO3 produced via the heterogeneous reaction sticks on surfaces and
C is accounted as depositional loss; calculate increased deposition loss
               IF ( V .EQ. HNO3_HIT ) THEN
                  S = HNO3_MAP
                  CONC( S,1 ) = POL( V ) + ( CONC( S,1 ) - POL( V ) ) * EFAC1( V )
                  DEPV_HNO3 = DEPVCR( V ) + PLDV_HONO / CONC( NO2_MAP,1 )
                  DD_FAC( V ) = DTDENS1 * DD_CONV( V ) * DEPV_HNO3
                  DDBF( V ) = DDBF( V ) + THETA * DD_FAC( V ) * CONC( S,1 )

C Use special treatment for NO2
C Loss of NO2 via the heterogeneous reaction is accounted for as an additional
C depositional loss. Add the loss of NO2 via the heterogeneous reaction
C to the regular deposition velocity (increased dep. vel.).  This will
C reduce the NO2 conc. in the atmosphere without affecting the depositional loss.
               ELSE IF ( V .EQ. NO2_HIT ) THEN
                  S = NO2_MAP
                  DEPV_NO2 = DEPVCR( V ) + 2.0 * PLDV_HONO / CONC( S,1 )
                  EFAC1 ( V ) = EXP( -DEPV_NO2 * RP )
                  EFAC2 ( V ) = EXP( -DEPV_NO2 * RQ )
                  POL   ( V ) = PLDV( V,C,R ) / DEPV_NO2
                  CONC( S,1 ) = POL( V ) + ( CONC( S,1 ) - POL( V ) ) * EFAC1( V )
                  DDBF( V ) = DDBF( V ) + THETA * DD_FAC( V ) * CONC( S,1 )

               ELSE IF ( V .EQ. HONO_HIT ) THEN
                  S = HONO_MAP
                  CONC( S,1 ) = POL( V ) + ( CONC( S,1 ) - POL( V ) ) * EFAC1( V )
                  DDBF( V ) = DDBF( V ) + THETA * DD_FAC( V ) * CONC( S,1 )
  
C --------- END of HET HONO RX ----------

               ELSE

C Pass selected N species to the BDSNP Soil NO emissions scheme

                  IF ( MGN_ONLN_DEP ) THEN

                    IF(SPECLOG) then
                      IF( V .eq. N_SPC_DEPV)  THEN
                       SPECLOG = .false. ! no need to do any species more than once
                       WRITE( LOGDEV,*) 'BDSNP Species list complete', speclog
                      END IF
                    END IF

                    IF ( (INDEX(TRIM( DV2DF_SPC( V ) ), 'NH3') .NE. 0) .OR.
     &                 (INDEX(TRIM( DV2DF_SPC( V ) ), 'NH4') .NE. 0) .OR.         
     &                 (INDEX(TRIM( DV2DF_SPC( V ) ), 'HNO3').NE. 0) .OR.
     &                 (INDEX(TRIM( DV2DF_SPC( V ) ), 'NO3') .NE. 0) .OR.
     &                 (INDEX(TRIM( DV2DF_SPC( V ) ), 'NO2') .NE. 0) .OR.
     &                 (INDEX(TRIM( DV2DF_SPC( V ) ), 'PAN') .NE. 0)) THEN


                      IF( SPECLOG ) THEN !write species each time it is used
                        WRITE( LOGDEV,*) 'BDSNP Dry Species Used:', TRIM(DV2DF_SPC( V ) ), V, N_SPC_DEPV
                      END IF

                      IF ( ( DDBF(V)- DDEP( V,C,R) ) .LT. 0.0 ) THEN !negative error checking

                       XMSG = 'Negative Deposition'
!                  WRITE( LOGDEV,*) 'BDSNP Negative Deposition vdiff, variable:', 
!     &            TRIM( DV2DF_SPC( V )), ( DDBF(V)- DDEP( V,C,R) ), C, R                  
!                      CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
                       CALL GET_N_DEP (DV2DF_SPC( V ), 0/
     &                           DTSEC, C, R )
                      else
                       CALL GET_N_DEP (DV2DF_SPC( V ), ( DDBF(V)- DDEP( V,C,R) )/
     &                           DTSEC, C, R )
                      END IF !end negative error checking

                     
                    END IF !end species check
                  
                  END IF !end BDSNP check

                  S = DV2DF( V )
                  CONC( S,1 ) = POL( V ) + ( CONC( S,1 ) - POL( V ) ) * EFAC1( V )
                  DDBF( V ) = DDBF( V ) + THETA * ( DD_FAC( V ) * CONC( S,1 ) 
C Add evasion as negative dep flux                  
     &                      - DTDENS1 * DD_CONV( V ) * PLDV( V,C,R ) )         
               END IF
               IF (ISPXURB .AND. (URBF .GT. 0.001)) then
                 DO L = 1, layrf
                 CONC( S,L ) =  CONC( S,L ) * EFAC1_wl( V,L )  
                 ENDDO
                 L = layrf
                 CONC( S,L ) =  CONC( S,L ) * EFAC1_rf( V )
               ENDIF

            END DO
   
            DO L = 1, NLAYS
               DO V = 1, N_SPC_DIFF
                  DD( V,L ) = 0.0
                  UU( V,L ) = 0.0
               END DO
            END DO


C Compute tendency of CBL concentrations - semi-implicit solution
C Set MATRIX1 elements A (col 1), B (diag.), E (superdiag.) and D (RHS)

            IF ( Met_Data%CONVCT( C,R ) ) THEN

               L = 1
               DO V = 1, N_SPC_DIFF
                  DD( V,L ) = CONC( V,L )
     &                      - LFAC1( L ) * CONC( V,L )
     &                      + LFAC2( L ) * CONC( V,L+1 ) 
               END DO

               DO L = 2, LCBL
                  DO V = 1, N_SPC_DIFF
                     DELC = MBARKS( L ) * CONC( V,1 ) !* NFRF(1) / VRF(L)
     &                    -   MDWN( L ) * CONC( V,L ) !* NFRF(L-1) / VRF(L)
     &                    +   MFAC( L ) * CONC( V,L+1 ) !* NFRF(L) / VRF(L)
                     DD( V,L ) = CONC( V,L ) + DFACQ * DELC
                  END DO


               END DO

               CALL MATRIX1 ( LCBL, AA, BB1, EE1, DD, UU )



C update conc
               DO L = 1, LCBL
                  DO V = 1, N_SPC_DIFF
                     CONC( V,L ) = UU( V,L )
                  END DO
               END DO

C reinitialize for TRI solver
               DO L = 1, NLAYS
                  DO V = 1, N_SPC_DIFF
                     DD( V,L ) = 0.0
                     UU( V,L ) = 0.0
                  END DO
               END DO

            END IF

            L = 1
            DO V = 1, N_SPC_DIFF
               DD( V,L ) = CONC( V,L )
     &                   + LFAC3( L ) * ( CONC( V,L+1 ) - CONC( V,L ) )
     &                   + EMIS( V,L ) / VRF(L)
            END DO


            DO L = 2, NLAYS-1
               DO V = 1, N_SPC_DIFF
                  DD( V,L ) = CONC( V,L )
     &                      + LFAC3( L ) * ( CONC( V,L+1 ) - CONC( V,L ) )
     &                      - LFAC4( L ) * ( CONC( V,L ) - CONC( V,L-1 ) )
     &                      + EMIS( V,L ) / VRF(L)
               END DO
            END DO

            L = NLAYS
            DO V = 1, N_SPC_DIFF
               DD( V,L ) = CONC( V,L )
     &                   - LFAC4( L ) * ( CONC( V,L ) - CONC( V,L-1 ) )
            END DO
            DO L = 1, NLAYS
               DO V = 1, N_SPC_DIFF
                  IF(DD(V,L)/CONC( V,L ).LT.0.5) then
                    print *,'C,R,L,V,DD,CONC,EDDY( L)=',C,R,L,V,DD( V,L ),CONC( V,L ),EDDY( L)
                    print *,'DEPV_wl,DEPV_rf=',DEPV_wl( 7,C,R,L ),DEPV_rf( 7,C,R )
                  ENDIF
               ENDDO
            ENDDO

            CALL TRI ( CC, BB2, EE2, DD, UU )

C Load into CGRID
            DO L = 1, NLAYS
               DO V = 1, N_SPC_DIFF
                  CONC( V,L ) = UU( V,L )
               END DO
            END DO
            






            DO V = 1, N_SPC_DEPV

C --------- HET HONO RX -----------------

               IF ( V .EQ. HNO3_HIT ) THEN
                  S = HNO3_MAP
                  CONC( S,1 ) = POL( V ) + ( CONC( S,1 ) - POL( V ) ) * EFAC2( V )
                  DDBF( V ) = DDBF( V ) + THBAR * DD_FAC( V ) * CONC( S,1 )


               ELSE IF ( V .EQ. NO2_HIT ) THEN
                  S = NO2_MAP
                  CONC( S,1 ) = POL( V ) + ( CONC( S,1 ) - POL( V ) ) * EFAC2( V )
                  DDBF( V ) = DDBF( V ) + THBAR * DD_FAC( V ) * CONC( S,1 )


               ELSE  IF ( V .EQ. HONO_HIT ) THEN
                  S = HONO_MAP
                  CONC( S,1 ) = POL( V ) + ( CONC( S,1 ) - POL( V ) ) * EFAC2( V )
                  DDBF( V ) = DDBF( V ) + THBAR * DD_FAC( V ) * CONC( S,1 )


C --------- END of HET HONO RX ----------

               ELSE
                  S = DV2DF( V )
                  CONC( S,1 ) = POL( V ) + ( CONC( S,1 ) - POL( V ) ) * EFAC2( V )
                  DDBF( V ) = DDBF( V ) + THBAR * ( DD_FAC( V ) * CONC( S,1 )
C Add evasion as negative dep flux                      
     &                      - DTDENS1 * DD_CONV( V ) * PLDV( V,C,R ) )
                 
               END IF
               IF (ISPXURB .AND. (URBF .GT. 0.001)) then
                 DO L = 1, layrf
                    IF(DEPV_wl( V,C,R,L ).lt.1.e-6) EXIT
                    CONC( S,L ) =  CONC( S,L ) * EFAC2_wl( V,L )  
                    DDBF( V ) = DDBF( V ) + CONC( S,l ) * DTS*Met_Data%DENS( C,R,L )
     &                        * DD_CONV( V ) * DEPV_wl( V,C,R,L )
                 ENDDO
                 IF(DEPV_rf( V,C,R ).gt.1.e-6) then
                   L = layrf
                   CONC( S,L ) =  CONC( S,L ) * EFAC2_rf( V )    
                   DDBF( V ) = DDBF( V ) + CONC( S,l ) * DTS*Met_Data%DENS( C,R,L )
     &                       * DD_CONV( V ) * DEPV_rf( V,C,R )
                 endif            
               ENDIF                  

            END DO
    

301      CONTINUE                 ! end sub time loop


         DO L = 1, NLAYS
            DO V = 1, N_SPC_DIFF
               CNGRD( DIFF_MAP( V ),L,C,R ) = CONC( V,L )
            END DO
         END DO

         DO V = 1, N_SPC_DEPV
            DDEP( V,C,R ) = DDBF( V )
         END DO




 
344   CONTINUE         !  end loop on col C
345   CONTINUE         !  end loop on row R

      RETURN
      END
