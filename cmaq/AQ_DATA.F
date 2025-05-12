
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
      Module aq_data

C  Defines gas and aerosol surrogate data.
C  this module replace the AQ_PARAMS.EXT include file
C
C  Revision History:
C     First version was coded in May 2010 by Steve Howard with
C     Prakash Bhave, Jeff Young, and Sergey Napelenok.
C     Feb 2011: Steve Howard, Shawn Roselle, Prakash Bhave, Bill Hutzell -
C     Mar 2011: Sergey Napelenok - AE6 update for Ca/Mg/K/Fe/Mn 
C               indexing, density, variable names
C     Mar 2011: Shawn Roselle, replaced I/O API include files with UTILIO_DEFN
C     Sep 2011: B.Hutzell: adapted for mercury and HAP mechanisms
C     Jun 2012: J.Young: remove full character blank padding for GNU Fortran (GCC) 4.1.2
C     Oct 2014: B.Hutzell:forcing to skip mapping if aerosol species
C     index equals the values for H3OP species as for aerosol water
C     species 
C-----------------------------------------------------------------------
      Use aero_data, Only: n_aerolist, n_aerospc, conmin, n_mode, 
     &                     iait, iacc, icor
      Use runtime_vars
#ifdef mpas
      Use util_module, only : index1
#endif

      Implicit None

      Integer, Parameter :: n_aq_gas = 15   ! number of gas-phase species for AQCHEM
      Integer, Parameter :: n_aq_aer = 31   ! number of aerosol surrogate species for AQCHEM

      Integer, Parameter :: n_aq_liqs = 64  ! number of liquid-phase species in AQCHEM
      Integer, Parameter :: mxsrg = 70  ! maximum number of species for surrogate

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...Gas surrogates species 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      Type sur_gas_type
         Character( 16 ) :: name               ! surrogate name
         Real            :: backgnd            ! background value [ppm]
      End Type sur_gas_type

      Type ( sur_gas_type ), Parameter :: srggas( n_aq_gas ) = (/
C                          Name       Background
C                       ----------    ----------
     & sur_gas_type( 'SO2       ',        0.0 ),
     & sur_gas_type( 'HNO3      ',        0.0 ),
     & sur_gas_type( 'N2O5      ',        0.0 ),
     & sur_gas_type( 'CO2       ',      340.0 ),
     & sur_gas_type( 'NH3       ',        0.0 ),
     & sur_gas_type( 'H2O2      ',        0.0 ),
     & sur_gas_type( 'O3        ',        0.0 ),
     & sur_gas_type( 'FOA       ',        0.0 ),
     & sur_gas_type( 'MHP       ',        0.0 ),
     & sur_gas_type( 'PAA       ',        0.0 ),
     & sur_gas_type( 'H2SO4     ',        0.0 ),
     & sur_gas_type( 'HCL       ',        0.0 ),
     & sur_gas_type( 'GLY       ',        0.0 ),
     & sur_gas_type( 'MGLY      ',        0.0 ),
     & sur_gas_type( 'HO        ',        0.0 ) / )

C...Pointers for the AQCHEM array GAS

      Integer, Parameter :: lso2    = 1  ! Sulfur Dioxide
      Integer, Parameter :: lhno3   = 2  ! Nitric Acid
      Integer, Parameter :: ln2o5   = 3  ! Dinitrogen Pentoxide
      Integer, Parameter :: lco2    = 4  ! Carbon Dioxide
      Integer, Parameter :: lnh3    = 5  ! Ammonia
      Integer, Parameter :: lh2o2   = 6  ! Hydrogen Peroxide
      Integer, Parameter :: lo3     = 7  ! Ozone
      Integer, Parameter :: lfoa    = 8  ! Formic Acid
      Integer, Parameter :: lmhp    = 9  ! Methyl Hydrogen Peroxide
      Integer, Parameter :: lpaa    = 10 ! Peroxyacetic Acid
      Integer, Parameter :: lh2so4  = 11 ! Sulfuric Acid
      Integer, Parameter :: lhcl    = 12 ! Hydrogen Chloride
      Integer, Parameter :: lgly    = 13 ! Glyoxal
      Integer, Parameter :: lmgly   = 14 ! Methylglyoxal
      Integer, Parameter :: lho     = 15 ! OH, hydroxyl radical

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...Aerosol surrogates species 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      Character(32) :: srgaer_name( n_aq_aer,n_mode )
      Type sur_aer_type
         Character( 16 ) :: name     ! bulk surrogate name
         Real            :: backgnd  ! background value [ug m-3]
         Real            :: density  ! density [kg m-3] 
         Real            :: molwt    ! molecular weight [g mol-1]
      End Type sur_aer_type

      Type ( sur_aer_type ), Parameter  :: srgaer( n_aq_aer ) = (/

C                      Name     Background  Density  MolWt
C                    ----------  ----------  ------  -----
     & sur_aer_type('SO4       ',   0.0,     1800.0,  96.0 ),
     & sur_aer_type('NH4       ',   0.0,     1800.0,  18.0 ),
     & sur_aer_type('NO3       ',   0.0,     1800.0,  62.0 ),
     & sur_aer_type('SOA       ',   0.0,     1400.0, 200.0 ),
     & sur_aer_type('AORGC     ',   0.0,     1400.0, 177.0 ),
     & sur_aer_type('POA       ',   0.0,     1400.0, 220.0 ),
     & sur_aer_type('EC        ',   0.0,     2200.0,  12.0 ),
     & sur_aer_type('PRI       ',   0.0,     1400.0, 218.0 ),
     & sur_aer_type('NA        ',   0.0,     2200.0,  23.0 ),
     & sur_aer_type('CL        ',   0.0,     2200.0,  35.5 ),
     & sur_aer_type('CA        ',   0.0,     2200.0,  40.1 ),
     & sur_aer_type('MG        ',   0.0,     2200.0,  24.3 ),
     & sur_aer_type('K         ',   0.0,     2200.0,  39.1 ),
     & sur_aer_type('FE        ',   0.0,     2200.0,  55.8 ),
     & sur_aer_type('MN        ',   0.0,     2200.0,  54.9 ),
     & sur_aer_type('SOIL      ',   0.0,     2600.0, 100.0 ),
     & sur_aer_type('ANTH      ',   0.0,     2200.0, 100.0 ),
     & sur_aer_type('SEAS      ',   0.0,     2200.0,  23.75),
     & sur_aer_type('TRACER    ',   0.0,     2200.0,  60.0 ),
     & sur_aer_type('PHG       ',   0.0,     2200.0, 200.5 ),
     & sur_aer_type('CACO3     ',   0.0,     2200.0, 100.08),
     & sur_aer_type('MGCO3     ',   0.0,     2200.0,  84.31),
     & sur_aer_type('A3FE      ',   0.0,     2200.0,  55.85),
     & sur_aer_type('B2MN      ',   0.0,     2200.0,  54.93),
     & sur_aer_type('SO4_AQH2O2',   0.0,     1800.0,  96.0 ),
     & sur_aer_type('SO4_AQO3  ',   0.0,     1800.0,  96.0 ),
     & sur_aer_type('SO4_AQFEMN',   0.0,     1800.0,  96.0 ),
     & sur_aer_type('SO4_AQMHP ',   0.0,     1800.0,  96.0 ),
     & sur_aer_type('SO4_AQPAA ',   0.0,     1800.0,  96.0 ),
     & sur_aer_type('NUM       ',   0.0,        1.0,   1.0 ),
     & sur_aer_type('SRF       ',   0.0,        1.0,   1.0 )/)

C...Pointers for the AQCHEM array AEROSOL
      Integer, Parameter :: lso4       = 1  !  Sulfate
      Integer, Parameter :: lnh4       = 2  !  Ammonium
      Integer, Parameter :: lno3       = 3  !  Nitrate
      Integer, Parameter :: lsoa       = 4  !  SOA
      Integer, Parameter :: lorgc      = 5  !  SOA (anthropogenic organic oligomers)
      Integer, Parameter :: lpoa       = 6  !  Primary organic aerosol
      Integer, Parameter :: lec        = 7  !  Elemental carbon
      Integer, Parameter :: lpri       = 8  !  Primary aerosol (Aitken and Accumulation)
      Integer, Parameter :: lna        = 9  !  Sodium
      Integer, Parameter :: lcl        = 10 !  Chloride ion
      Integer, Parameter :: lca        = 11 !  Ca in Accumulation mode (AE6) SLN 16March2011
      Integer, Parameter :: lmg        = 12 !  Mg in Accumulation mode (AE6) SLN 16March2011
      Integer, Parameter :: lk         = 13 !  K in Accumulation mode (AE6)  SLN 16March2011
      Integer, Parameter :: lfe        = 14 !  Fe in Accumulation mode (AE6) SLN 22March2011
      Integer, Parameter :: lmn        = 15 !  Mn in Accumulation mode (AE6) SLN 22March2011
      Integer, Parameter :: lsoil      = 16 !  SOIL in Coarse mode (AE6)     SLN 16March2011
      Integer, Parameter :: lanth      = 17 !  CORS in Coarse mode (AE6)     SLN 16March2011
      Integer, Parameter :: lseas      = 18 !  SEAS in Coarse mode (AE6)     SLN 16March2011
      Integer, Parameter :: ltracer    = 19 !  Tracer aerosols (Aitken)
      Integer, Parameter :: lphg       = 20 !  Mercuric aerosols (Aitken and Accumulation)
      Integer, Parameter :: lcaco3     = 21 !  Calcium Carbonate aerosol (place holder)
      Integer, Parameter :: lmgco3     = 22 !  Magnesium Carbonate aerosol (place holder)
      Integer, Parameter :: la3fe      = 23 !  Iron aerosol (place holder)
      Integer, Parameter :: lb2mn      = 24 !  Manganese aerosol (place holder)
      Integer, Parameter :: lso4aqh2o2 = 25 ! Sulfate from aq H2O2 rxn
      Integer, Parameter :: lso4aqo3   = 26 ! Sulfate from aq O3 rxn
      Integer, Parameter :: lso4aqfemn = 27 ! Sulfate from aq rxn catalyzed by FE&MN
      Integer, Parameter :: lso4aqmhp  = 28 ! Sulfate from aq MHP rxn
      Integer, Parameter :: lso4aqpaa  = 29 ! Sulfate from aq PAA rxn
      Integer, Parameter :: lnum       = 30 !  Number
      Integer, Parameter :: lsrf       = 31 !  Surface area

C...Pointers for the AQCHEM arrays LIQUID and WETDEP
      Integer, Parameter :: lacl        = 01  ! Hydrogen ion
      Integer, Parameter :: lnh4accl    = 02  ! Ammonium(attributed to accum. mode)
      Integer, Parameter :: lcacorl     = 03  ! coarse mode Calcium 
      Integer, Parameter :: lnaaccl     = 04  ! Sodium
      Integer, Parameter :: lohl        = 05  ! Hydroxyl radical
      Integer, Parameter :: lso4accl    = 06  ! Sulfate (attributed to accumulation mode)
      Integer, Parameter :: lhso4accl   = 07  ! bisulfate (attributed to accumulation mode)
      Integer, Parameter :: lso3l       = 08  ! sulfite
      Integer, Parameter :: lhso3l      = 09  ! bisulfite
      Integer, Parameter :: lso2l       = 10  ! sulfur dioxide
      Integer, Parameter :: lco3l       = 11  ! carbonate
      Integer, Parameter :: lhco3l      = 12  ! bicarbonate
      Integer, Parameter :: lco2l       = 13  ! carbon dioxide
      Integer, Parameter :: lno3accl    = 14  ! nitrate(attributed to accumulation mode)
      Integer, Parameter :: lnh3l       = 15  ! ammonia
      Integer, Parameter :: lclaccl     = 16  ! chloride ion (attributed to accumulation mode)
      Integer, Parameter :: lh2o2l      = 17  ! hydrogen peroxide
      Integer, Parameter :: lo3l        = 18  ! ozone
      Integer, Parameter :: lfecorl     = 19  ! iron
      Integer, Parameter :: lmncorl     = 20  ! Manganese
      Integer, Parameter :: lal         = 21  ! generalized anion associated with iron
      Integer, Parameter :: lfoal       = 22  ! Formic acid
      Integer, Parameter :: lhco2l      = 23  ! HCOO- ion
      Integer, Parameter :: lmhpl       = 24  ! Methyl hydrogen peroxide
      Integer, Parameter :: lpaal       = 25  ! Peroxyacidic acid
      Integer, Parameter :: lhcll       = 26  ! Hydrogen chloride
      Integer, Parameter :: lpriml      = 27  ! primary aerosol
      Integer, Parameter :: lmgcorl     = 28  ! coarse mode Magnesium
      Integer, Parameter :: lkcorl      = 29  ! coarse mode Potassium
      Integer, Parameter :: lbl         = 30  ! generalized anion associated with manganese
      Integer, Parameter :: lhno3l      = 31  ! nitric acid
!     Integer, Parameter :: lprimcorl   = 32  ! coarse-mode primary aerosol
      Integer, Parameter :: lnumcorl    = 32  ! coarse-mode number
      Integer, Parameter :: lts6corl    = 33  ! sulfate (attributed to coarse mode)
      Integer, Parameter :: lnacorl     = 34  ! sodium (attributed to coarse mode)
      Integer, Parameter :: lclcorl     = 35  ! chloride ion (attributed to coarse mode)
      Integer, Parameter :: lno3corl    = 36  ! nitrate (attributed to coarse mode)
      Integer, Parameter :: lnh4corl    = 37  ! ammonium (attributed to coarse mode)
      Integer, Parameter :: lsoal       = 38  ! secondary organic aerosols
      Integer, Parameter :: lorgcl      = 39  ! in-cloud soa
      integer, parameter :: lpoal       = 40  ! primary organic aerosols
      Integer, Parameter :: lecl        = 41  ! elemental carbon
      Integer, Parameter :: lglyl       = 42  ! glyoxal
      Integer, Parameter :: lmglyl      = 43  ! methylglyoxal
      Integer, Parameter :: lcaaccl     = 44  ! accumulation mode Calcium (AE6)   SLN 16March2011
      Integer, Parameter :: lmgaccl     = 45  ! accumulation mode Magnesium (AE6) SLN 16March2011
      Integer, Parameter :: lkaccl      = 46  ! accumulation mode Potassium (AE6) SLN 16March2011
      Integer, Parameter :: lsoilcl     = 47  ! coarse mode SOIL (AE6)   SLN 16March2011
      Integer, Parameter :: lanthcl     = 48  ! coarse mode CORS (lol) (AE6)   SLN 16March2011
      Integer, Parameter :: lseascl     = 49  ! coarse mode SEAS (AE6) SLN 16March2011
      Integer, Parameter :: lfeaccl     = 50  ! accumulation mode Iron (AE6)      SLN 22March2011
      Integer, Parameter :: lmnaccl     = 51  ! accumulation mode Manganese (AE6) SLN 22March2011
C for TX version
      integer, parameter :: ltracerl    = 52 ! TRACER
      integer, parameter :: ltracercorl = 53 ! Coarse-mode TRACER
      integer, parameter :: lphgfinel   = 54 ! Fine-mode Mercury Tracer
      integer, parameter :: lphgcorl    = 55 ! Coarse-mode Mercury Tracer
C for sulfur tracking model
      Integer, Parameter :: lts6aqh2o2l = 56 ! total S(VI) from aq H2O2 rxn
      Integer, Parameter :: lts6aqo3l   = 57 ! total S(VI) from aq O3 rxn
      Integer, Parameter :: lts6aqfemnl = 58 ! total S(VI) from aq FE&MN catalyzed rxn
      Integer, Parameter :: lts6aqmhpl  = 59 ! total S(VI) from aq MHP rxn
      Integer, Parameter :: lts6aqpaal  = 60 ! total S(VI) from aq PAA rxn
      
      Integer, Parameter :: ltracil = 61
      Integer, Parameter :: lpoail = 62
      Integer, Parameter :: lpriil = 63
      Integer, Parameter :: lsoail = 64

      Integer, Save :: n_aqsrg_gas( n_aq_gas )             ! Store Number of CGRID species pointing to each surrogate
      Integer, Save :: n_aqsrg_aer( n_aq_aer,n_mode )      ! Store Number of CGRID species pointing to each surrogate
      Integer, Allocatable, Save :: MAP_CGRIDtoAQGAS( : )  ! Map from CGRID to Surrogates
      Integer, Allocatable, Save :: MAP_CGRIDtoAQAERO( : ) ! Map from CGRID to Surrogates
      Integer, Allocatable, Save :: MAP_CGRIDtoAQMODE( : ) ! Map from CGRID to Surrogates
      Integer, Allocatable, Save :: MAP_ACCtoAIT( : )      ! Map Accumulation to Aitken Mode Species
      Real( 8 ),    Save         :: FRACAIT( N_AQ_AER )    ! Store fraction of Aitken mode in each aerosol surrogate

      Real( 8 )   :: jh2o2_hydrometeors    ! H2O2 photolysis rate in condensed water [sec-1]
      Real( 8 )   :: jhno3_hydrometeors    ! HNO3 photolysis rate in condensed water [sec-1]

      Logical, Private, Save :: mapped = .False.

      Contains
 
C------------------------------------------------------------------------
      Subroutine aq_data_init()

C  Defines mapping from CGRID for species concentration and moments.
 
C  Revision History:
C     First version was coded in May 2010 by Steve Howard with
C     Prakash Bhave, Jeff Young, and Sergey Napelenok.
C     Jan 2017: Additional error check added to avoid array out of
C     bounds when additional aerosol species are added. HOT Pye.
C-----------------------------------------------------------------------

      Use cgrid_spcs          ! CGRID mechanism species
      Use utilio_defn
      Use aero_data 
      Use phot_mod, Only: init_phot_shared  ! set-up shared photolysis arrays
      
      Implicit None

C Functions:
      Integer, External :: indexn

C Local Variables:
      Character( 300 ) :: xmsg
      Character( 10 ) :: bkground
      Integer :: igas
      Integer :: iaer
      Integer :: isrg
      Integer :: spc
      Integer :: i, im, im0, iaero, iaero0
      Integer :: n
      Integer :: var
      Logical :: new, Jexist
      Integer :: lsrg( mxsrg )
      Integer :: ios

C Local Variable used to determine what optional TX to add to surrogates
      !Logical :: found_all_modes
      !Integer :: n_missing_tx
      !Integer :: n_which_txadd( nadd_tx )
      CHARACTER(16) :: NAMEI, NAMEJ
      INTEGER JAIT, JACC, JAQ_ACC, JAQ_AIT

      Logical, Save :: Firstcall = .True.
      
C-----------------------------------------------------------------------


      If ( Firstcall ) Then
         Firstcall = .False.
      Else
         Return
      End if
 
      Call init_phot_shared()
      
C...Allocate Maps from CGRID to Aqueous Chemistry Vectors
      Allocate ( MAP_CGRIDtoAQGAS( n_cgrid_spc ), stat = ios )
      If ( ios .Ne. 0 ) Then
         xmsg = '*** Error allocating MAP_CGRIDtoAQGAS'
         Call m3exit( 'aq_data_init', 0, 0, xmsg, xstat3 )
      End If

      Allocate ( MAP_CGRIDtoAQAERO( n_cgrid_spc ), stat = ios )
      If ( ios .Ne. 0 ) Then
         xmsg = '*** Error allocating MAP_CGRIDtoAQAERO'
         Call m3exit( 'aq_data_init', 0, 0, xmsg, xstat3 )
      End If

      Allocate ( MAP_CGRIDtoAQMODE( n_cgrid_spc ), stat = ios )
      If ( ios .Ne. 0 ) Then
         xmsg = '*** Error allocating MAP_CGRIDtoAQMODE'
         Call m3exit( 'aq_data_init', 0, 0, xmsg, xstat3 )
      End If

      Allocate ( MAP_ACCtoAIT( n_cgrid_spc ), stat = ios )
      If ( ios .Ne. 0 ) Then
         xmsg = '*** Error allocating MAP_ACCtoAIT'
         Call m3exit( 'aq_data_init', 0, 0, xmsg, xstat3 )
      End If

C...Map aero species
      Call map_aero()

C...Initialize aerosol surrogate species names 
      Do spc = 1, n_aq_aer
        Do im = 1,n_mode
          srgaer_name( spc,im ) = trim(srgaer( spc )%name) // '_' // trim(modename( im ))
        End Do
      End Do

C...Create Map Directly from CGRID to the gas and aerosol species in AQCHEM
      MAP_CGRIDtoAQGAS  = 0
      MAP_CGRIDtoAQAERO = 0
      MAP_CGRIDtoAQMODE = 0
      N_AQSRG_GAS = 0
      N_AQSRG_AER = 0
            
      DO I = 1,N_SPC_AQ
         SPC = MAP_AQtoCGRID( I ) 
         IF ( CGRID_MASK_GAS( SPC ) .OR. CGRID_MASK_NR( SPC ) .OR.
     &        CGRID_MASK_TRAC( SPC ) ) THEN
           ! Map Gas-Phase Species
           igas = index1( aq_spc( I ), n_aq_gas, SRGGAS(:)%name )
           if ( igas .eq. 0 ) then
             Write ( xmsg, '(A,A,A,A,A,A)' ),'Aqueous Chem Surrogate ',trim( aq_spc(I) ),
     &             ' for species ',trim( CGRID_NAME( SPC )),' is not available as a gas in the ',
     &             ' Aq Chem Solver. Please Correct or Delete.'
             Call m3exit( 'AQ_INIT', 0, 0, xmsg, xstat3 )
           else
             MAP_CGRIDtoAQGAS( SPC ) = IGAS
             N_AQSRG_GAS( IGAS ) = N_AQSRG_GAS( IGAS ) + 1
           end if
         ELSEIF ( CGRID_MASK_AERO( SPC ) ) THEN
           ! Map Aerosol Species
           ! Loop over modes; find the AE surrogate name on the local
           ! mode-dependent cloud chemistry list; there will only be one
           ! match across all modes.
           iaero = 0
           do im0 = 1,n_mode
             iaero0 = index1( aq_spc( I ), n_aq_aer, SRGAER_name(:,im0) )
             if ( iaero0 .ne. 0 ) then
                 iaero = iaero0
                 im = im0
             end if
           end do
           if ( iaero .eq. 0 ) then
             Write ( xmsg, '(A,A,A,A,A,A)' ),'Aqueous Chem Surrogate ',trim( aq_spc(I) ),
     &             ' for species ',trim( CGRID_NAME( SPC )),' is not available as an ',
     &             'aerosol in the Aq Chem Solver. Please Correct or Delete.'
             Call m3exit( 'AQ_INIT', 0, 0, xmsg, xstat3 )
           else
             MAP_CGRIDtoAQAERO( SPC ) = IAERO
             MAP_CGRIDtoAQMODE( SPC ) = IM
             N_AQSRG_AER( IAERO,IM ) = N_AQSRG_AER( IAERO,IM ) + 1
           end if
        END IF
      END DO

      mapped = .True.
      
C...Match Accumulation mode surrogate indices to their Aitken mode counterparts
      MAP_ACCtoAIT = 0
      DO IAERO = 1,N_AEROSPC
         ! Find the CGRID Indeces for the Accumulation and Aitken 
         ! mode species
         JACC = AEROSPC_MAP( IAERO, IACC )
         JAIT = AEROSPC_MAP( IAERO, IAIT )

         ! Determine if these species have aqueous chemistry surrogates
         IF ( JACC .GT. 0 .AND. JAIT .GT. 0 ) THEN 
            JAQ_ACC = MAP_CGRIDtoAQAERO( JACC )
            JAQ_AIT = MAP_CGRIDtoAQAERO( JAIT )

            ! If both species are mapped and both have aqueous-chemistry
            ! surrogates that are the same, then link them.
            IF ( JAQ_ACC .GT. 0 .AND. JAQ_ACC .EQ. JAQ_AIT ) 
     &         MAP_ACCtoAIT( JACC ) = JAIT
         END IF
      END DO

      Write( logdev,'( / 5x, a / )' ) ' --- AQCHEM Gas and Aerosol Species Mapped ---'

      Return
      End Subroutine aq_data_init

C-----------------------------------------------------------------------
      Subroutine set_aeconcmin( aeconcmin )
C-----------------------------------------------------------------------
C
C  Function:
C       set minimum concentrations for aerosol species
C       routine called from convcld_acm
C
C  Preconditions:
C       Dates and times represented YYYYDDD:HHMMSS.
C
C  Revision History:
C       5/05 copied code from aero_depv to begin subroutine
C       4/08 J.Kelly: added code to set min for coarse surface area
C       21 Jun 10 J.Young: convert for Namelist redesign
C-----------------------------------------------------------------------

      Use cgrid_spcs            ! CGRID mechanism species
      Use aero_data
      Use utilio_defn
      Use const, only : mwair
      Implicit None

C..Arguments:
      Real, Intent( Out ) :: aeconcmin( : ) ! array of minimum concentrations

C..Local Variables:

      Integer m        ! mode loop counter
      Integer n        ! aerospc loop counter
      Integer spc      ! species loop counter


      ! compute minimums for numbers and surfaces

      Do m = 1, n_mode
         spc = index1( aeromode( m )%num_name, n_ae_spc, ae_spc )
         aeconcmin( spc ) = aeromode_minNum( m ) * mwair * 1.0E-3

         spc = index1( aeromode( m )%srf_name, n_ae_spc, ae_spc )
         aeconcmin( spc ) = aeromode_minM2( m ) * mwair * 1.0E-3
      End Do

      ! compute minimums for aerosol species
      Do m = 1, n_mode
         Do n = 1, n_aerospc
            If ( .NOT. aero_missing(n,m) ) Then
               spc = index1( aerospc( n )%name( m ), n_ae_spc, ae_spc )
               If ( spc .Ne. 0 ) Then
                  aeconcmin( spc ) = aerospc( n )%min_conc( m ) * mwair * 1.0E-9 / ae_molwt( spc )
               End If               
            End If
         End Do
      End Do 

      Return
      End Subroutine set_aeconcmin

      End Module aq_data
