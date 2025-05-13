!------------------------------------------------------------------------------!
!  The Community Multiscale Air Quality (CMAQ) system software is in           !
!  continuous development by various groups and is based on information        !
!  from these groups: Federal Government employees, contractors working        !
!  within a United States Government contract, and non-Federal sources         !
!  including research institutions.  These groups give the Government          !
!  permission to use, prepare derivative works of, and distribute copies       !
!  of their work in the CMAQ system to the public and to permit others         !
!  to do so.  The United States Environmental Protection Agency                !
!  therefore grants similar permission to use the CMAQ system software,        !
!  but users are requested to provide copies of derivative works or            !
!  products designed to operate in the CMAQ system to the United States        !
!  Government without restrictions as to use by others.  Software              !
!  that is used with the CMAQ system but distributed under the GNU             !
!  General Public License or the GNU Lesser General Public License is          !
!  subject to their copyright restrictions.                                    !
!------------------------------------------------------------------------------!

MODULE const

!------------------------------------------------------------------------------! 
! Purpose: Fundamental constants for air quality modeling. Includes functions
!          for calculating saturation vapor pressure and its reciprocal.
!------------------------------------------------------------------------------! 
 
!------------------------------------------------------------------------------! 
! References:
!   CRC Handbook of Chemistry and Physics (76th Ed), CRC Press, 1995 
!
!   Hobbs, P.V., Basic Physical Chemistry for the Atmospheric Sciences,
!     Cambridge Univ. Press, 206 pp, 1995.  
!
!   Snyder, J.P., Map Projections-A Working Manual, U.S. Geological Survey
!     Paper 1395 U.S.GPO, Washington, DC, 1987.
!
!   Stull, R.B., An Introduction to Boundary Layer Meteorology, Kluwer, 
!     Dordrecht, 1988
!
!   NIST, The International System of Units (SI). Newell, D.B. and 
!     Tiesinga, E., eds. NIST Special Publication 330, 2019.
!     doi: 10.6028/nist.sp.330-2019.
!------------------------------------------------------------------------------! 

  IMPLICIT NONE

! Geometric Constants:
  REAL,      PARAMETER :: PI = 3.14159265
  REAL( 8 ), PARAMETER :: DPI = 3.14159265358979324D0
  Real( 8 ), Parameter :: f6dpi = 6.0D0 / dpi
  Real,      Parameter :: f6pi = 6.0 / pi
  Real( 8 ), Parameter :: dpi6 = dpi / 6.0D0
  Real,      Parameter :: pi6 = pi / 6.0

! pi/180 [ rad/deg ]
  REAL,      PARAMETER :: PI180  = PI / 180.0

! Geodetic Constants:
 
  ! radius of the earth [ m ]
  ! -- radius of sphere having same surface area as Clarke ellipsoid of 1866 
  !    ( Source: Snyder, 1987)
 !REAL,      PARAMETER :: REARTH = 6370997.0
  REAL,      PARAMETER :: REARTH = 6370000.0    ! default Re in MM5 and WRF
 
  ! length of a sidereal day [ sec ]  (Source: CRC76 p. 14-6 )
  REAL,      PARAMETER :: SIDAY = 86164.09
 
  ! mean gravitational acceleration [ m/sec**2 ]
  ! -- mean of polar and equatorial values  (Source: CRC76 p. 14-6)
  REAL,      PARAMETER :: GRAV = 9.80622

  ! latitude degrees to meters
  REAL,      PARAMETER :: DG2M = REARTH * PI180

  ! Solar constant [ W/m**2 ] (Source: CRC76 p. 14-2)
  ! UNUSED IN CMAQ
  !   REAL, PARAMETER :: SOLCNST = 1373.0

! Fundamental Constants: (Source: CRC76, pp. 1-1 to 1-6 except where 
!                         specified otherwise)

  ! Avogadro's Constant [ number/mol ] 
  ! Exact definition. Source: NIST 2019 
  REAL,      PARAMETER :: AVO  = 6.02214076E23 
  REAL( 8 ), PARAMETER :: DAVO = 6.02214076D23

  ! Boltzmann Constant [ J K-1 ]
  ! Exact definition. Source: NIST 2019
  REAL,      PARAMETER :: KBOLTZ   = 1.380649E-23 
  REAL( 8 ), PARAMETER :: DKBOLTZ  = 1.380649D-23 

  ! Exact definition. Source: NIST 2019
  REAL,      PARAMETER :: RGASUNIV  =  AVO * KBOLTZ
  REAL( 8 ), PARAMETER :: DRGASUNIV = DAVO * DKBOLTZ

  ! standard atmosphere  [ Pa ]
  REAL,      PARAMETER :: STDATMPA = 101325.0

  ! Standard Temperature [ K ]
  REAL,      PARAMETER :: STDTEMP = 273.15

  ! Stefan-Boltzmann constant [ W/(m**2 K**4) ]
  ! UNUSED IN CMAQ
  ! REAL, PARAMETER :: STFBLZ = 5.67051E-8

! Non-MKS

  ! Molar volume at STP [ L/mol ] Non MKS units 
   REAL,     PARAMETER :: MOLVOL = 22.41410

! Atmospheric Constants: 

! mean molecular weight for dry air [ g/mol ]
! FSB: 78.06% N2, 21% O2, and 0.943% Ar on a mole 
! fraction basis (Source: Hobbs, pp. 69-70)
   REAL,     PARAMETER :: MWAIR = 28.9628
   Real,     Parameter :: inv_mwair   = 1.0E3 / mwair  ! [ mol/kg ]

! dry-air gas constant [ J / kg-K ]
   REAL,     PARAMETER :: RDGAS = 1.0E3 * RGASUNIV / MWAIR   ! 287.07548994

! mean molecular weight for water vapor [ g/mol ]
   REAL,     PARAMETER :: MWWAT = 18.0153

! gas constant for water vapor [ J/kg-K ]
   REAL,     PARAMETER :: RWVAP = 1.0E3 * RGASUNIV / MWWAT   ! 461.52492604

! FSB NOTE: CPD, CVD, CPWVAP and CVWVAP are calculated assuming dry air and
! water vapor are classical ideal gases, i.e. vibration does not contribute
! to internal energy.

! specific heat of dry air at constant pressure [ J/kg-K ]
   REAL,     PARAMETER :: CPD = 7.0 * RDGAS / 2.0            ! 1004.7642148 

! specific heat of dry air at constant volume [ J/kg-K ]
   REAL,     PARAMETER :: CVD = 5.0 * RDGAS / 2.0            ! 717.68872485

! specific heat for water vapor at constant pressure [ J/kg-K ]
   REAL,     PARAMETER :: CPWVAP = 4.0 * RWVAP               ! 1846.0997042

! specific heat for water vapor at constant volume [ J/kg-K ]
   REAL,     PARAMETER :: CVWVAP = 3.0 * RWVAP               ! 1384.5747781

! vapor press of water at 0 C [ Pa ] Source: CRC76 pp. 6-15
   REAL,     PARAMETER :: VP0 = 611.29

! latent heat of vaporization of water at 0 C [ J/kg ] 
! Source: Stull, p. 641
   REAL,     PARAMETER :: LV0 = 2.501E6

! Rate of change of latent heat of vaporization with
! respect to temperature [ J/kg-K ]
! Source: Stull, p. 641
   REAL,     PARAMETER :: DLVDT = 2370.0

! latent heat of fusion of water at 0 C [ J/kg ]
! Source: Stull, p. 641
   REAL,     PARAMETER :: LF0 = 3.34E5

   CONTAINS

! Function to compute saturation vapor pressure of water as a function 
!  of air temperature [Pa]
      REAL FUNCTION ESATL( TT )
        REAL, INTENT(IN) :: TT  ! air temperature in K

        REAL, PARAMETER :: AL = 610.94
        REAL, PARAMETER :: BL = 17.625
        REAL, PARAMETER :: CL = 243.04

! *** values of AL, BL, and CL are from:
!     Alduchov and Eskridge, "Improved Magnus Form Approximations of
!                            Saturation Vapor Pressure,"
!                            Jour. of Applied Meteorology, vol. 35,
!                            pp 601-609, April, 1996.

         ESATL = AL * EXP( BL * ( TT - 273.15 ) / ( TT - 273.15 + CL ) )
       END FUNCTION ESATL

! Function to compute reciprocal of saturation vapor pressure of water
!   as a function of air temperature. [Pa-1] 
      REAL( 8 ) FUNCTION INV_ESATL( TT )
        REAL( 8 ), INTENT( IN ) :: TT
! *** values of AL, BL, and CL are from:
!     Alduchov and Eskridge, "Improved Magnus Form Approximations of
!                            Saturation Vapor Pressure,"
!                            Jour. of Applied Meteorology, vol. 35,
!                            pp 601-609, April, 1996.
        REAL( 8 ), PARAMETER :: AL = 610.94D0 
        REAL( 8 ), PARAMETER :: BL = 17.625D0
        REAL( 8 ), PARAMETER :: CL = 243.04D0
        REAL( 8 ), PARAMETER :: DL = 1.0D0 / AL
        INV_ESATL = DL * DEXP( BL * ( 273.15D0 - TT ) / ( TT - 273.15D0 + CL ) )
      END FUNCTION INV_ESATL

END MODULE const

