
!------------------------------------------------------------------------!
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
      SUBROUTINE INIT_CCTM_MODULES ( CGRID, JDATE, JTIME, TSTEP, ASTEP )

      USE RUNTIME_VARS, ONLY : COLUMN_MODEL

      IMPLICIT NONE

C Arguments:
      REAL, POINTER :: CGRID( :,:,:,: )

      INTEGER, INTENT(IN)    :: JDATE        ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN)    :: JTIME        ! current model time, coded HHMMSS
      INTEGER, INTENT(IN)    :: TSTEP( 3 )   ! time step vector (HHMMSS)
                                             ! TSTEP(1) = local output step
                                             ! TSTEP(2) = sciproc sync. step (chem)
                                             ! TSTEP(3) = twoway model time step w.r.t. wrf time
                                             !            step and wrf/cmaq call frequency

      INTEGER, INTENT(IN)    :: ASTEP( : )   ! layer advection time step
      

      INTERFACE
#ifndef mpas
         SUBROUTINE HADV ( CGRID, JDATE, JTIME, TSTEP, ASTEP )
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 3 )
            INTEGER, INTENT( IN )     :: ASTEP( : )
         END SUBROUTINE HADV
         SUBROUTINE ZADV ( CGRID, JDATE, JTIME, TSTEP )
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 3 )
         END SUBROUTINE ZADV
         SUBROUTINE COUPLE ( CGRID, JDATE, JTIME, TSTEP )
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 3 )
         END SUBROUTINE COUPLE
         SUBROUTINE DECOUPLE ( CGRID, JDATE, JTIME, TSTEP )
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 3 )
         END SUBROUTINE DECOUPLE
         SUBROUTINE HDIFF ( CGRID, JDATE, JTIME, TSTEP )
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 3 )
         END SUBROUTINE HDIFF
#endif
         SUBROUTINE VDIFF ( CGRID, JDATE, JTIME, TSTEP )
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 3 )
         END SUBROUTINE VDIFF
         SUBROUTINE PHOT ( CGRID, JDATE, JTIME, TSTEP )
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( : )
         END SUBROUTINE PHOT
         SUBROUTINE CLDPROC ( CGRID, JDATE, JTIME, TSTEP )
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 3 )
         END SUBROUTINE CLDPROC
         SUBROUTINE CHEM ( CGRID, JDATE, JTIME, TSTEP )
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 3 )
         END SUBROUTINE CHEM
         SUBROUTINE PVO3 ( CGRID, JDATE, JTIME )
            IMPLICIT NONE
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
         END SUBROUTINE PVO3
         SUBROUTINE AERO ( CGRID, JDATE, JTIME, TSTEP )
            REAL, POINTER             :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )     :: JDATE, JTIME
            INTEGER, INTENT( IN )     :: TSTEP( 3 )
         END SUBROUTINE AERO

      END INTERFACE

      ! Calculate Vertical Diffusion, Dry Deposition, and Emissions
      CALL VDIFF ( CGRID, JDATE, JTIME, TSTEP )

      IF ( .NOT. COLUMN_MODEL ) THEN

          !CALL HADV ( CGRID, JDATE, JTIME, TSTEP, ASTEP )
          !CALL ZADV ( CGRID, JDATE, JTIME, TSTEP )
          !CALL HDIFF ( CGRID, JDATE, JTIME, TSTEP )

      END IF
          
      CALL PHOT ( CGRID, JDATE, JTIME, TSTEP )    
      CALL CLDPROC ( CGRID, JDATE, JTIME, TSTEP )
      !CALL CHEM ( CGRID, JDATE, JTIME, TSTEP )
      !CALL AERO ( CGRID, JDATE, JTIME, TSTEP )

      ! Write CGRID state file for subsequent runs
      CALL WR_CGRID ( CGRID, JDATE, JTIME, TSTEP( 1 ) )
#ifdef isam
      CALL WR_SA_CGRID ( JDATE, JTIME, TSTEP( 1 ) )
#endif
#ifdef sens
      CALL WR_SENGRID ( JDATE, JTIME, TSTEP( 1 )  )
#endif

      ! Open ISAM and DDM Concentration Files
#ifdef isam
      ! Write out ISAM source-resolved concentrations
      CALL WR_AVG_SA ( JDATE, JTIME, TSTEP( 1 ) )
#endif


      END SUBROUTINE INIT_CCTM_MODULES

