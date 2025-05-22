
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

      MODULE MIO_ASCII

      IMPLICIT NONE

      INTEGER, SAVE :: N_MIO_IN, N_MIO_OUT

      TYPE MIO_FILE_DATA_TYPE
          CHARACTER(16)              :: FILE_LOGICAL
          INTEGER                    :: NLAYS
          INTEGER                    :: BLEV
          INTEGER                    :: ELEV
          INTEGER                    :: NVARS
          CHARACTER(32), ALLOCATABLE :: VARNAMES(:)
          INTEGER,       ALLOCATABLE :: VARTYPES(:)
          CHARACTER(16), ALLOCATABLE :: VARUNITS(:)
          CHARACTER(256),ALLOCATABLE :: VARDESC (:)
          INTEGER,       ALLOCATABLE :: NDIMS(:)
          LOGICAL,       ALLOCATABLE :: L_TSTEP(:)
          LOGICAL,       ALLOCATABLE :: L_LAY(:)
          LOGICAL,       ALLOCATABLE :: L_COL(:)
          LOGICAL,       ALLOCATABLE :: L_ROW(:)
          LOGICAL,       ALLOCATABLE :: L_VEXT(:)
      END TYPE MIO_FILE_DATA_TYPE

      TYPE( MIO_FILE_DATA_TYPE ) :: MIO_TXT_FILE_DATA( 200 )

      INTEGER :: NDIMS3D(2048)
      LOGICAL :: L_TSTEP(2048)
      LOGICAL :: L_LAY(2048)
      LOGICAL :: L_COL(2048)
      LOGICAL :: L_ROW(2048)
      LOGICAL :: L_VEXT(2048)

      CONTAINS

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE INIT_MIO_FILE_DATA

          N_MIO_IN  = 0
          N_MIO_OUT = 0

      END SUBROUTINE INIT_MIO_FILE_DATA

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE LOAD_MIO_FILE ( FILE_LOGICAL, NLAYS, BLEV, ELEV, VARNAMES, 
     &           VARTYPES, VARUNITS, VARDESC, NDIMS, 
     &           L_TSTEP, L_LAY, L_COL, L_ROW, L_VEXT )


      IMPLICIT NONE

      CHARACTER(*),  INTENT( IN )  :: FILE_LOGICAL
      INTEGER,       INTENT( IN )  :: NLAYS
      INTEGER,       INTENT( IN )  :: BLEV
      INTEGER,       INTENT( IN )  :: ELEV
      CHARACTER(*),  INTENT( IN )  :: VARNAMES( : )
      INTEGER,       INTENT( IN )  :: VARTYPES( : )
      CHARACTER(16), INTENT( IN )  :: VARUNITS( : )
      CHARACTER(*),  INTENT( IN )  :: VARDESC ( : )
      INTEGER,       INTENT( IN )  :: NDIMS( : )
      LOGICAL,       INTENT( IN )  :: L_TSTEP(:), L_LAY(:), 
     &                                L_COL(:), L_ROW(:), L_VEXT(:)

      INTEGER NVARS
      INTEGER IOS

      NVARS = SIZE( VARTYPES )
      N_MIO_OUT = N_MIO_OUT + 1

      ! Set high-level parameters for this output file in the MIO file
      ! structure
      MIO_TXT_FILE_DATA( N_MIO_OUT )%FILE_LOGICAL = FILE_LOGICAL
      MIO_TXT_FILE_DATA( N_MIO_OUT )%NLAYS = NLAYS
      MIO_TXT_FILE_DATA( N_MIO_OUT )%BLEV  = BLEV
      MIO_TXT_FILE_DATA( N_MIO_OUT )%ELEV  = ELEV
      MIO_TXT_FILE_DATA( N_MIO_OUT )%NVARS = NVARS

      ! Allocate MIO File Structure for this Output File
      ALLOCATE( MIO_TXT_FILE_DATA( N_MIO_OUT )%VARNAMES( NVARS ) ,
     &          MIO_TXT_FILE_DATA( N_MIO_OUT )%VARTYPES( NVARS ) ,
     &          MIO_TXT_FILE_DATA( N_MIO_OUT )%VARUNITS( NVARS ) ,
     &          MIO_TXT_FILE_DATA( N_MIO_OUT )%VARDESC ( NVARS ) ,
     &          MIO_TXT_FILE_DATA( N_MIO_OUT )%NDIMS   ( NVARS ) ,
     &          MIO_TXT_FILE_DATA( N_MIO_OUT )%L_TSTEP ( NVARS ) ,
     &          MIO_TXT_FILE_DATA( N_MIO_OUT )%L_LAY   ( NVARS ) ,
     &          MIO_TXT_FILE_DATA( N_MIO_OUT )%L_COL   ( NVARS ) ,
     &          MIO_TXT_FILE_DATA( N_MIO_OUT )%L_ROW   ( NVARS ) ,
     &          MIO_TXT_FILE_DATA( N_MIO_OUT )%L_VEXT  ( NVARS ) ,
     &          STAT = IOS )
      CALL CHECKMEM( IOS, 'MIO_TXT_FILE_DATA','LOAD_MIO_FILE' )

      ! Populate MIO File Structure
      MIO_TXT_FILE_DATA( N_MIO_OUT )%VARNAMES = VARNAMES
      MIO_TXT_FILE_DATA( N_MIO_OUT )%VARTYPES = VARTYPES
      MIO_TXT_FILE_DATA( N_MIO_OUT )%VARUNITS = VARUNITS
      MIO_TXT_FILE_DATA( N_MIO_OUT )%VARDESC  = VARDESC
      MIO_TXT_FILE_DATA( N_MIO_OUT )%NDIMS    = NDIMS
      MIO_TXT_FILE_DATA( N_MIO_OUT )%L_TSTEP  = L_TSTEP
      MIO_TXT_FILE_DATA( N_MIO_OUT )%L_LAY    = L_LAY
      MIO_TXT_FILE_DATA( N_MIO_OUT )%L_COL    = L_COL
      MIO_TXT_FILE_DATA( N_MIO_OUT )%L_ROW    = L_ROW
      MIO_TXT_FILE_DATA( N_MIO_OUT )%L_VEXT   = L_VEXT

      END SUBROUTINE LOAD_MIO_FILE

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE WRITE_MIO_ASCII ( FNAME )

C   This subroutine will create a mock-up of the MIO input file which
C   specifies the number of input and output files, their filenames, and
C   the variables on each file, along with metadata for each variable.
C   The subroutine is expected to be used to demonstrate that the
C   updates made to driver.F to initialize the CMAQ system are
C   successfully making available all of the information needed by MIO.
C   After these updates are confirmed, approved, and merged, it is
C   expected that this routine will either be removed, or it may be kept
C   as an optional diagnostic for MIO characterization. It should not be
C   used in default operations of CMAQ.
C-----------------------------------------------------------------------

      USE UTILIO_DEFN 
      USE RUNTIME_VARS, ONLY : OUTDIR, APPL_NAME

      IMPLICIT NONE

      CHARACTER(*)   :: FNAME
      INTEGER        :: IOUT, IFILE, IVAR
      CHARACTER(10)  :: VT
      CHARACTER(100) :: DIMSTRING
      CHARACTER(400) :: FPATH

C-----------------------------------------------------------------------


      FPATH = TRIM(OUTDIR) // '/' // TRIM(FNAME) // '_' //
     &        TRIM(APPL_NAME) // '.txt'
      IOUT = JUNIT()
      OPEN( IOUT,  FILE = FPATH, STATUS = 'UNKNOWN' )

      WRITE ( IOUT, '(A3, I0, A5, I0)' ) 'in ' , N_MIO_IN, 
     &                                  ' out ', N_MIO_OUT

      ! Write Output File Metadata
      DO IFILE = 1,N_MIO_OUT
         WRITE ( IOUT,* ) 
         WRITE ( IOUT, '(A)' ) MIO_TXT_FILE_DATA( IFILE )%FILE_LOGICAL
         WRITE ( IOUT, '(I0)' ) MIO_TXT_FILE_DATA( IFILE )%NVARS

         DO IVAR = 1,MIO_TXT_FILE_DATA( IFILE )%NVARS 
           ! Find string for vairbale type
#ifndef mpas
             IF ( MIO_TXT_FILE_DATA( IFILE )%VARTYPES( IVAR ) .EQ. M3REAL ) VT = 'real'
#else 
           ! Assume real
             VT = 'real'
#endif

           ! Make string for dimensions
           DIMSTRING = ''
           IF ( MIO_TXT_FILE_DATA( IFILE )%L_TSTEP(IVAR) ) DIMSTRING = TRIM(DIMSTRING) // ' TSTEP'
           IF ( MIO_TXT_FILE_DATA( IFILE )%L_LAY(IVAR) )   DIMSTRING = TRIM(DIMSTRING) // ' LAY'
           IF ( MIO_TXT_FILE_DATA( IFILE )%L_ROW(IVAR) )   DIMSTRING = TRIM(DIMSTRING) // ' ROW'
           IF ( MIO_TXT_FILE_DATA( IFILE )%L_COL(IVAR) )   DIMSTRING = TRIM(DIMSTRING) // ' COL'
           IF ( MIO_TXT_FILE_DATA( IFILE )%L_VEXT(IVAR) )  DIMSTRING = TRIM(DIMSTRING) // ' VEXT'

           WRITE ( IOUT, '(A,1x,A,1x,I0,A,1x,A1,A,A1,1x,A1,A,A1)' ) 
     &             TRIM( MIO_TXT_FILE_DATA( IFILE )%VARNAMES( IVAR ) ), TRIM(VT),
     &             MIO_TXT_FILE_DATA( IFILE )%NDIMS( IVAR ),TRIM(DIMSTRING),
     &             '"',TRIM(MIO_TXT_FILE_DATA( IFILE )%VARUNITS( IVAR )),'"',
     &             '"',TRIM(MIO_TXT_FILE_DATA( IFILE )%VARDESC ( IVAR )),'"'
         END DO
      END DO

      CLOSE( IOUT )

      END SUBROUTINE WRITE_MIO_ASCII

      END MODULE MIO_ASCII
