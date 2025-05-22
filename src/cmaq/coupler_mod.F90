!===============================================================================
! Purpose:  To provide a platform for exchanging data between met and CMAQ 
!           models in the WRF-CMAQ and MPAS-CMAQ coupled model systems.
!
! Revised:  31 Aug 2023  Original version.  (D. Wong)
!===============================================================================

! Note: In a generic form, a coupler serves as a middle man to supply a piece
!       of data from model X to model Y or vice versa (diagram below):
!
!                  storing               retrieving
!       var A    ----------->  coupler  ----------->  var B
!       model X  <-----------           <-----------  model Y
!                 retrieving               storing
!
!       In terms of memory requirement, it needs space to hold A in model X,
!       B in model Y, and A and B in coupler.
!
!       In some instances, retrieving var B takes place in multiple subourintes
!       in model Y. This means multiple copies of B are required, i.e. increase
!       of memory footprint. To mitigate this potential large memory footprint
!       scenario, implementation of the generic form of the coupler was relaxed
!       by accessing the coupler's data structure directly. A set of indices
!       were determined beforehand to facilitate such access.
!
!       In MPAS, spatial data is stored in 1D array (NROWS = 1) but such data is
!       stored in coupler's 2D data structure. Similarly, three dimensional 
!       volumn data in MPAS is stored in a 2D array (NROWS = 1) but such data is 
!       stored in coupler's 3D data structure. In addition, 2D MPAS array 
!       defined as (nlays, ncols) is converted to (ncols, nrows, nlays) 
!       structure, where nrows = 1.
!
!       Variable naming convention:
!         time dependency ( i / d )        -- time independent vs time dependent
!         grep type ( r / e )              -- regular cross vs extended dot grid
!         transfer direction ( m2a / a2m ) -- from met to CMAQ or CMAQ to met
!
!       For MPAS-CMAQ coupled model: 
!         * there are m met variables (2D and 3D together) and all CMAQ 
!           concentration fields, n variables, are transferring from MPAS to 
!           CMAQ. Only those m met variables' name are stored and sorted. 
!         * All n variables are communicated back to MPAS from CMAQ at once so 
!           coupler_a2m_vname does not bear any significance.
!         * All n variables are communicated between MPAS and CMAQ using
!           coupler_m2a_dr_3d_data only without coupler_a2m_3d_data to save
!           space.

module coupler_module

  implicit none

  integer, parameter :: wrf_cmaq  = 0
  integer, parameter :: mpas_cmaq = 1
  integer, parameter :: mta       = 0
  integer, parameter :: atm       = 1
  integer, parameter :: vname_max_str_len = 32

  integer, parameter :: n_so4 = 3

! various data types
  integer, parameter :: ir2 = 1   ! time independent regular  met      2d data
  integer, parameter :: ie2 = 2   ! time independent extended met      2d data
  integer, parameter :: dr2 = 3   ! time   dependent regular  met      2d data
  integer, parameter :: dr3 = 4   ! time   dependent regular  met      3d data
  integer, parameter :: de3 = 5   ! time   dependent extended met      3d data
  integer, parameter :: fb3 = 6   ! time independent regular  feedback 3d data

! main coupler data structure
! time independent 2d data (:,:,:)     - first two are spatial dimensions and
!                                        last one is variable
! time   dependent 2d data (:,:,:,:)   - first two are spatial dimensions, third
!                                        one variable, and last time for time
!                                        circular buffer, 0:1
! time   dependent 3d data (:,:,:,:,:) - first three are spatial dimensions,
!                                        fourth one variable, and last time for
!                                        time circular buffer, 0:1 for WRF-CMAQ
!                                        and 1 for MPAS-CMAQ since it dose not 
!                                        require time interpolation

! coupler_m2a_rd_3d_data contains data transfer from met model to CMAQ model
! namely met fields information. However on the MPAS-CMAQ side, it also contains
! conc fields. Hence the data arrangement is conc fields first and then met
! fields. The conc fields are in and out at once always but the met fields are
! being accessed individually. With this consideration, variable 
! coupler_a2m_nvars indicates the number of met variables from met model to CMAQ
! model and it is the same as coupler_m2a_mnvars (for MPAS-CMAQ only). With this
! reasoning, coupler_m2a_vname only contains met variable names only. Variable
! coupler_m2a_ncvars indicates the number of conc variables in the coupled model
! (for MPAS-CMAQ only).

  character (vname_max_str_len), allocatable :: coupler_m2a_vname(:)              ! var names from met to CMAQ model
  character (vname_max_str_len), allocatable :: coupler_a2m_vname(:)              ! var names from CMAQ to met model
  integer, allocatable                       :: coupler_m2a_data_type(:)          ! data type from met to AQ  model
  integer, allocatable                       :: coupler_a2m_data_type(:)          ! data type from AQ  to met model
  integer, allocatable                       :: coupler_m2a_nlays(:)
  integer                                    :: coupler_a2m_nlays
  integer, allocatable                       :: coupler_m2a_loc(:)
  integer, allocatable                       :: coupler_a2m_loc(:)
  integer, allocatable                       :: data_m2a_head(:)
  integer, allocatable                       :: data_m2a_tail(:)
  real, allocatable                          :: coupler_m2a_ir_2d_data(:,:,:)     ! time indep. regular  2d met data
  real, allocatable                          :: coupler_m2a_ie_2d_data(:,:,:)     ! time indep. extended 2d met data
  real, allocatable                          :: coupler_m2a_dr_2d_data(:,:,:,:)   ! time   dep. regular  2d met data
  real, allocatable                          :: coupler_m2a_dr_3d_data(:,:,:,:,:) ! time   dep. regular  3d met data
  real, allocatable                          :: coupler_m2a_de_3d_data(:,:,:,:,:) ! time   dep. extended 3d met data
  real, allocatable                          :: coupler_a2m_3d_data(:,:,:,:)

  integer :: coupler_m2a_nvars             ! # of vars from met to CMAQ model
  integer :: coupler_a2m_nvars             ! # of vars from CMAQ to met modol
  integer :: coupler_m2a_nmvars            ! # of met  vars from met to CMAQ model, MPAS only
  integer :: coupler_m2a_ncvars            ! # of conc vars from met to CMAQ model, MPAS only

  integer :: coupler_ncols, coupler_nrows

  integer :: so4_ind(n_so4)                ! location of different mode of aerosol SO4 

  integer :: coupled_model                 ! 0 - WRF-CMAQ, 1 - MPAS-CMAQ

  integer :: coupled_model_tstep
  integer :: coupled_model_stdate
  integer :: coupled_model_sttime
  real    :: coupled_model_xcent
  real    :: coupled_model_ycent

  logical :: coupler_initialized = .false.

! for WRF-CMAQ model:
! met data
  integer, parameter :: n_gridcro2d_var = 7
  integer, parameter :: n_griddot2d_var = 1
  integer, parameter :: n_metcro3d_var  = 18
  integer, parameter :: n_metdot3d_var  = 4
  integer, parameter :: n_metcro2d_var  = 44

  character (len = vname_max_str_len), parameter :: gridcro2d_vlist(n_gridcro2d_var) = &
    (/ 'LAT   ', 'LON   ', 'MSFX2 ', 'HT    ', 'LWMASK',                               &
       'PURB  ', 'DLUSE '                                                    /)

  ! for checking purposes only
  character (len = vname_max_str_len), parameter :: gridcro2d_units(n_gridcro2d_var) = &
    (/ 'DEGREES ', 'DEGREES ', '(M/M)**2', 'M       ', '-       ',                     &
       'PERCENT ', 'CATEGORY'                                                /)

  character (len = vname_max_str_len), parameter :: griddot2d_vlist(n_griddot2d_var) = &
    (/ 'MSFD2'                                                               /)

  ! for checking purposes only
  character (len = vname_max_str_len), parameter :: griddot2d_units(n_griddot2d_var) = &
    (/ '(M/M)**2'                                                            /)

  character (len = vname_max_str_len), parameter :: metcro3d_vlist(n_metcro3d_var) =   &
    (/ 'JACOBF ', 'JACOBM ', 'DENSA_J', 'TA     ', 'QV     ',                          &
       'QC     ', 'QR     ', 'QI     ', 'QS     ', 'QG     ',                          &
       'PRES   ', 'DENS   ', 'ZH     ', 'ZF     ', 'UWIND  ',                          &
       'VWIND  ', 'PV     ', 'UST_WL '                                       /)

  ! for checking purposes only
  character (len = vname_max_str_len), parameter :: metcro3d_units(n_metcro3d_var) =   &
    (/ 'M               ', 'M               ',                                         &
       'KG/M**2         ', 'K               ',                                         &
       'KG/KG           ', 'KG/KG           ',                                         &
       'KG/KG           ', 'KG/KG           ',                                         &
       'KG/KG           ', 'KG/KG           ',                                         &
       'Pa              ', 'KG/M**3         ',                                         &
       'M               ', 'M               ',                                         &
       'M/S             ', 'M/S             ',                                         &
       'M^2*K/KG/S * E-6', 'M/S             '                                /)

  character (len = vname_max_str_len), parameter :: metdot3d_vlist(n_metdot3d_var) =   &
    (/ 'UWINDC ', 'VWINDC ', 'UHAT_JD', 'VHAT_JD'                            /)

  ! for checking purposes only
  character (len = vname_max_str_len), parameter :: metdot3d_units(n_metdot3d_var) =   &
    (/ 'M/S     ', 'M/S     ', 'KG/(M*S)', 'KG/(M*S)'                        /)

  character (len = vname_max_str_len), parameter :: metcro2d_vlist(n_metcro2d_var) =   &
    (/ 'PRSFC    ', 'USTAR    ', 'WSTAR    ', 'PBL      ', 'ZRUF     ',                &
       'MOLI     ', 'HFX      ', 'RA       ', 'RS       ', 'WSPD10   ',                &
       'GSW      ', 'RGRND    ', 'RNA      ', 'RCA      ', 'CFRAC    ',                &
       'CLDT     ', 'CLDB     ', 'WBAR     ', 'SNOCOV   ', 'VEG      ',                &
       'TEMP2    ', 'WR       ', 'TEMPG    ', 'LAI      ', 'SLTYP    ',                &
       'Q2       ', 'SEAICE   ', 'SOIM1    ', 'SOIM2    ', 'SOIT1    ',                &
       'SOIT2    ', 'LH       ', 'WWLT_PX  ', 'WFC_PX   ', 'WSAT_PX  ',                &
       'CLAY_PX  ', 'CSAND_PX ', 'FMSAND_PX', 'BUILD_AREA_FRACTION', 'BUILD_HEIGHT',   &
       'BLKW1    ', 'BLKW2    ', 'UST_RF   ', 'GAB_RF   '                    /)

  ! for checking purposes only
  character (len = vname_max_str_len), parameter :: metcro2d_units(n_metcro2d_var) =   &
    (/ 'Pascal    ', 'M/S       ', 'M/S       ', 'M         ', 'M         ',           &
       '1/M       ', 'WATTS/M**2', 'S/M       ', 'S/M       ', 'M/S       ',           &
       'WATTS/M**2', 'WATTS/M**2', 'CM        ', 'CM        ', 'FRACTION  ',           &
       'M         ', 'M         ', 'G/M**3    ', 'NODIM     ', 'NO UNIT   ',           &
       'K         ', 'M         ', 'K         ', 'AREA/AREA ', '-         ',           &
       'KG/KG     ', 'FRACTION  ', 'M**3/M**3 ', 'M**3/M**3 ', 'K         ',           &
       'K         ', 'WATTS/M**2', 'M**3/M**3 ', 'M**3/M**3 ', 'M**3/M**3 ',           &
       'FRACTION  ', 'FRACTION  ', 'FRACTION  ', 'FRACTION  ', 'M         ',           &
       'M         ', 'M         ', 'M/S       ', 'M/S       '                /)

! feedback data
    ! water soluble
    integer :: num_ws_spc(3)

    integer, allocatable :: ws_spc_index(:,:)

    ! water insoluble
    integer :: num_wi_spc(3)

    integer, allocatable :: wi_spc_index(:,:)

    ! elmental carbon
    integer, parameter :: num_ec_spc = 2

    integer :: ec_spc_index(num_ec_spc)

    character (len = vname_max_str_len), parameter :: ec_spc(num_ec_spc) = &
      (/ 'AECI', 'AECJ'                                        /)

    ! sea salt
    integer, parameter :: num_ss_spc = 5

    integer :: ss_spc_index(num_ss_spc)

    character (len = vname_max_str_len), parameter :: ss_spc(num_ss_spc) = &
      (/ 'ANAJ   ', 'ACLJ   ', 'ACLK   ', 'ASO4K  ', 'ASEACAT' /)

    ! water
    integer, parameter :: num_h2o_spc = 4

    integer :: h2o_spc_index(num_h2o_spc)

    character (len = vname_max_str_len), parameter :: h2o_spc(num_h2o_spc) = &
      (/ 'AH2OI   ', 'AH2OJ   ', 'AH2OK   ', 'AORGH2OJ'        /)

    integer, parameter :: n_cmaq_dfb_vars = 22   ! number of aerosol direct feedback (dfb) variables

    character (16), parameter :: cmaq_dfb_vlist(n_cmaq_dfb_vars) =    &
      (/ 'ws_1            ', 'ws_2            ', 'ws_3            ',  &
         'is_1            ', 'is_2            ', 'is_3            ',  &
         'ec_1            ', 'ec_2            ', 'ec_3            ',  &
         'ss_1            ', 'ss_2            ', 'ss_3            ',  &
         'h2o_1           ', 'h2o_2           ', 'h2o_3           ',  &
         'diam_1          ', 'diam_2          ', 'diam_3          ',  &
         'sd_1            ', 'sd_2            ', 'sd_3            ',  &
         'O3              '                                           &
      /)

! for 2d variable
!        integer, parameter :: n2d_data_old   = 35

!        integer, parameter :: cfrac2dr_ind =  1
!        integer, parameter :: chlo_ind     =  2
!        integer, parameter :: dms_ind      =  3
!        integer, parameter :: hfx_ind      =  4
!        integer, parameter :: ht_ind       =  5
!        integer, parameter :: lai_ind      =  6
!        integer, parameter :: lat_ind      =  7
!        integer, parameter :: lh_ind       =  8
!        integer, parameter :: lon_ind      =  9
!        integer, parameter :: lwmask_ind   = 10
!        integer, parameter :: open_ind     = 11
!        integer, parameter :: pbl_ind      = 12
!        integer, parameter :: prsfc_ind    = 13
!        integer, parameter :: purb_ind     = 14
!        integer, parameter :: q2_ind       = 15
!        integer, parameter :: ra_ind       = 16     ! aerodynamic resistance
!        integer, parameter :: rainc_ind    = 17     ! time-step convective precipitation
!        integer, parameter :: rgrnd_ind    = 18
!        integer, parameter :: rainnc_ind   = 19     ! time-step total grid-scale precipitation
!        integer, parameter :: rs_ind       = 20     ! surface resistance
!        integer, parameter :: seaice_ind   = 21
!        integer, parameter :: sltyp_ind    = 22
!        integer, parameter :: snocov_ind   = 23
!        integer, parameter :: soit1_ind    = 24
!        integer, parameter :: surf_ind     = 25
!        integer, parameter :: temp2_ind    = 26
!        integer, parameter :: tempg_ind    = 27
!        integer, parameter :: ustar_ind    = 28
!        integer, parameter :: vegpx_ind    = 29
!        integer, parameter :: canwat_ind   = 30 
!        integer, parameter :: wspd10_ind   = 31     ! 2m wind speed
!        integer, parameter :: znt_ind      = 32
!        integer, parameter :: cellArea_ind = 33     ! cell area, m**2
!        integer, parameter :: cfrac2dt_ind = 34
!        integer, parameter :: rmol_ind     = 35

! alphabetical order and upper case letter goes first than lower case letter
!        character (20), parameter :: vname_2d_old(n2d_data_old) =                &   ! in ascending order
!           (/ 'CFRAC   ',  'CHLO    ',  'DMS     ',  'HFX     ',  'HT      ',    &
!              'LAI     ',  'LAT     ',  'LH      ',  'LON     ',  'LWMASK  ',    &
!              'OPEN    ',  'PBL     ',  'PRSFC   ',  'PURB    ',  'Q2      ',    &
!              'RA      ',  'RC      ',  'RGRND   ',  'RN      ',  'RS      ',    &
!              'SEAICE  ',  'SLTYP   ',  'SNOCOV  ',  'SOIT1   ',  'SURF    ',    &
!              'TEMP2   ',  'TEMPG   ',  'USTAR   ',  'VEG     ',  'WR      ',    &
!              'WSPD10  ',  'ZRUF    ',  'cellArea',  'cfrac2dt',  'rmol    '  /)

! for 3d variable
!        integer, parameter :: n3d_data_old     = 18

!        integer, parameter :: cfrac3d_ind    =  1
!        integer, parameter :: dens_ind       =  2
!        integer, parameter :: densa_j_ind    =  3
!        integer, parameter :: qc_ind         =  5
!        integer, parameter :: qg_ind         =  6
!        integer, parameter :: qi_ind         =  7
!        integer, parameter :: qr_ind         =  8
!        integer, parameter :: qs_ind         =  9
!        integer, parameter :: cldfracwcu_ind = 15
!        integer, parameter :: eddy_ind       = 16
!        integer, parameter :: qc_cu_ind      = 17
!        integer, parameter :: qi_cu_ind      = 18

    real, allocatable :: smois_data(:,:,:)                ! surface layer

    integer, private :: num_land_cat   ! for WRF-CMAQ only

    character (16) :: mminlu_wrf

! for MPAS-CMAQ model:
! met data

! 2d variable
    integer, parameter :: n2d_data = 32

! Note: RA       - aerodynamic resistance
!       RC       - time-step convective precipitation
!       RN       - time-step total grid-scale precipitation
!       RS       - surface resistance
!       WSPD10   - 2m wind speed
!       cellArea - cell area, m**2
    character (vname_max_str_len), parameter :: vname_2d (n2d_data) =        &
      (/ 'CFRAC   ',  'HFX     ',  'HT      ',  'LAI     ',  'LAT     ',     &
         'LH      ',  'LON     ',  'LWMASK  ',  'PBL     ',  'PRSFC   ',     &
         'PURB    ',  'Q2      ',  'RA      ',  'RC      ',  'RGRND   ',     &
         'RN      ',  'RS      ',  'SEAICE  ',  'SLTYP   ',  'SNOCOV  ',     &
         'SOIT1   ',  'SOIT2   ',  'TEMP2   ',  'TEMPG   ',  'USTAR   ',     &
         'VEG     ',  'WR      ',  'WSPD10  ',  'ZRUF    ',  'cellArea',     &
         'cfrac2dt',  'rmol    '                                         /)

! 3d variable
    integer, parameter :: n3d_data = 18

! Note: 1. cellht means cell thickness
!          cellvol means cell volumn
!          inv means inverse/reciprical
!          mlvl means mid level
!       2. currently cldfracwcu is not used

    character (vname_max_str_len), parameter :: vname_3d (n3d_data) =                                &
      (/ 'CFRAC_3D     ',  'DENS         ',  'DENSA_J      ',  'PRES         ',  'QC           ',    &
         'QG           ',  'QI           ',  'QR           ',  'QS           ',  'QV           ',    &
         'TA           ',  'WSPD         ',  'ZF           ',  'ZH           ',  'cldfracwcu   ',    &
         'eddy         ',  'qc_cu        ',  'qi_cu        '                                      /)

  character (16) :: mminlu_mpas

  character (19) :: ctm_out_clock

  real, allocatable :: lufrac_data(:,:),                 &
                       coupler_ocean(:,:),               &
                       coupler_szone(:,:),               &
                       cell_area(:,:),                   &
                       cell_vol(:,:,:),                  &     ! cell volume!
                       inv_cell_vol(:,:,:),              &     ! reciprical of cell volume
                       cell_ht(:,:,:),                   &     ! cell thickness full level
                       inv_cell_ht(:,:,:),               &     ! reciprical of cell thickness
                       inv_mlvl_cell_ht(:,:,:)                 ! reciprical of cell thickness mid level

  character(1000) :: namelist

  integer  :: my_gc_adj
  integer  :: my_ae_adj
  integer  :: my_nr_adj

  integer :: my_emis_buffer_ind, my_emis_tstep, zf_loc, zh_loc,      &
             cfrac3d_loc, eddy_loc, wspd_loc, temp_loc, dens_loc,    &
             qv_loc, pres_loc, qc_loc, qr_loc, qg_loc, qi_loc,       &
             qs_loc, qc_cu_loc, qi_cu_loc, densa_j_loc, ustar_loc,   &
             hfx_loc, pbl_loc

  logical :: mpas_cmaq_last_step = .false.
  logical :: mpas_diag

  private :: quicksort_vname, find_n_items_in_namelist, binary_search

  interface coupler_data_storing
    module procedure coupler_1d_rdata_storing,      &
                     coupler_1d_idata_storing,      &
                     coupler_2d_rdata_storing,      &
                     coupler_3d_rdata_storing
  end interface

  interface coupler_data_retrieving
    module procedure coupler_1d_rdata_retrieving,   &
                     coupler_1d_idata_retrieving,   &
                     coupler_2d_rdata_retrieving,   &
                     coupler_2d_idata_retrieving,   &
                     coupler_3d_rdata_retrieving
  end interface

  contains

! -------------------------------------------------------------------------
  subroutine coupler_init (dim1, dim2, dim3, num_land_cat)

    use get_env_module

    integer, intent(in) :: dim1, dim2, dim3
    integer, optional, intent(in) :: num_land_cat

    integer :: stat, start, end,       &
               my_n_gc,                &
               my_n_ae,                &
               my_n_nr,                &
               nvars, n, nn

#ifdef twoway
    coupled_model = wrf_cmaq
#else
    coupled_model = mpas_cmaq
#endif

    coupler_a2m_nlays = dim3

    if (coupled_model == wrf_cmaq) then
       nvars =   n_gridcro2d_var      &
               + num_land_cat         &
               + n_griddot2d_var      &
               + n_metcro3d_var       &
               + n_metdot3d_var       &
               + n_metcro2d_var

       coupler_m2a_nvars = nvars

       coupler_a2m_nvars = n_cmaq_dfb_vars

       allocate (coupler_m2a_vname(coupler_m2a_nvars),         &
                 coupler_m2a_data_type(coupler_m2a_nvars),     &
                 coupler_m2a_nlays(coupler_m2a_nvars),         &
                 coupler_m2a_loc(coupler_m2a_nvars),           &
                 data_m2a_head(coupler_m2a_nvars),             &
                 data_m2a_tail(coupler_m2a_nvars),             &
                 coupler_a2m_vname(coupler_a2m_nvars),         &
                 coupler_a2m_data_type(coupler_a2m_nvars),     &
                 coupler_a2m_loc(coupler_a2m_nvars),           &
                 stat=stat)

       coupler_a2m_vname = cmaq_dfb_vlist
       do n = 1, coupler_a2m_nvars
          coupler_a2m_loc(n)       = n
          coupler_a2m_data_type(n) = fb3
       end do

       data_m2a_head = -1
       data_m2a_tail = -1

       nn = 0
       do n = 1, n_gridcro2d_var
          nn = nn + 1
          coupler_m2a_vname(nn) = 'c_' // gridcro2d_vlist(n)
          coupler_m2a_data_type(nn) = ir2
          coupler_m2a_nlays(nn)     = 1
          coupler_m2a_loc(nn)       = n
       end do

       do n = 1, num_land_cat
          nn = nn + 1
          write (coupler_m2a_vname(nn), '(a9, i2.2)') 'c_LUFRAC_', n
          coupler_m2a_data_type(nn) = ir2
          coupler_m2a_nlays(nn)     = 1
          coupler_m2a_loc(nn)       = nn
       end do

       do n = 1, n_griddot2d_var
          nn = nn + 1
          coupler_m2a_vname(nn) = 'c_' // griddot2d_vlist(n)
          coupler_m2a_data_type(nn) = ie2
          coupler_m2a_nlays(nn)     = 1
          coupler_m2a_loc(nn)       = n
       end do

       do n = 1, n_metcro2d_var
          nn = nn + 1
          coupler_m2a_vname(nn) = 'c_' // metcro2d_vlist(n)
          coupler_m2a_data_type(nn) = dr2
          coupler_m2a_nlays(nn)     = 1
          coupler_m2a_loc(nn)       = n
       end do

       do n = 1, n_metcro3d_var
          nn = nn + 1
          coupler_m2a_vname(nn) = 'c_' // metcro3d_vlist(n)
          coupler_m2a_data_type(nn) = dr3
          coupler_m2a_nlays(nn)     = dim3
          coupler_m2a_loc(nn)       = n
       end do

       do n = 1, n_metdot3d_var
          nn = nn + 1
          coupler_m2a_vname(nn) = 'c_' // metdot3d_vlist(n)
          coupler_m2a_data_type(nn) = de3
          coupler_m2a_nlays(nn)     = dim3
          coupler_m2a_loc(nn)       = n
       end do

       allocate (coupler_m2a_ir_2d_data(dim1,   dim2,         n_gridcro2d_var+num_land_cat),  &  ! for gridcro2d
                 coupler_m2a_ie_2d_data(dim1+1, dim2+1,       n_griddot2d_var),               &  ! for griddot2d +1 for dot
                 coupler_m2a_dr_2d_data(dim1+2, dim2+2,       n_metcro2d_var, 0:1),           &  ! for  metcro2d +2 for ghost cells
                 coupler_m2a_dr_3d_data(dim1+2, dim2+2, dim3, n_metcro3d_var, 0:1),           &  ! for  metcro3d +2 for ghost cells
                 coupler_m2a_de_3d_data(dim1+3, dim2+3, dim3, n_metdot3d_var, 0:1),           &  ! for  metdot3d +2 for ghost cells +1 for dot
                 coupler_a2m_3d_data(dim1,dim2,dim3,coupler_a2m_nvars),                       &
                 stat=stat)

       call quicksort_vname (coupler_m2a_vname, coupler_m2a_data_type, &
                             coupler_m2a_loc, 1, coupler_m2a_nvars,    &
                             coupler_m2a_nlays)

       call quicksort_vname (coupler_a2m_vname, coupler_a2m_data_type, &
                             coupler_a2m_loc, 1, coupler_a2m_nvars)

    else if (coupled_model == mpas_cmaq) then

       call get_env (namelist, 'gc_matrix_nml', ' ')
       call find_n_items_in_namelist (namelist, my_n_gc)
       call get_env (namelist, 'ae_matrix_nml', ' ')
       call find_n_items_in_namelist (namelist, my_n_ae, so4_ind)
       call get_env (namelist, 'nr_matrix_nml', ' ')
       call find_n_items_in_namelist (namelist, my_n_nr)

       so4_ind = so4_ind + my_n_gc

       coupler_a2m_nvars = my_n_gc + my_n_ae + my_n_nr

       coupler_m2a_ncvars = coupler_a2m_nvars

       nvars = n2d_data + n3d_data
       coupler_m2a_nmvars = nvars
       coupler_m2a_nvars  = nvars

       allocate (coupler_m2a_vname(coupler_m2a_nvars),                                        &
                 coupler_a2m_vname(coupler_a2m_nvars),                                        &
                 coupler_m2a_data_type(coupler_m2a_nvars),                                    &
                 coupler_m2a_nlays(coupler_m2a_nvars),                                        &
                 coupler_m2a_loc(coupler_m2a_nvars),                                          &
                 coupler_m2a_dr_2d_data(dim1, dim2, n2d_data, 1),                             &
                 coupler_m2a_dr_3d_data(dim1, dim2, dim3, n3d_data+coupler_m2a_ncvars, 1),    &
!                coupler_a2m_3d_data(dim1,dim2,dim3,coupler_m2a_ncvars),                      &
                 stat=stat)

       coupler_a2m_vname      = 'fb'
       coupler_m2a_dr_2d_data = 0.0
       coupler_m2a_dr_3d_data = 0.0
!      coupler_a2m_3d_data    = 0.0

       nn = 0
       do n = 1, n2d_data
          nn = nn + 1
          coupler_m2a_vname(nn)     = 'c_' // vname_2d(n)
          coupler_m2a_data_type(nn) = dr2
          coupler_m2a_nlays(nn)     = 1
          coupler_m2a_loc(nn)       = n
       end do

       do n = 1, n3d_data
          nn = nn + 1
          coupler_m2a_vname(nn)     = 'c_' // vname_3d(n)
          coupler_m2a_data_type(nn) = dr3
          coupler_m2a_nlays(nn)     = dim3
          coupler_m2a_loc(nn)       = n + coupler_m2a_ncvars
       end do

       call quicksort_vname (coupler_m2a_vname, coupler_m2a_data_type, &
                             coupler_m2a_loc, 1, coupler_m2a_nmvars,   &
                             coupler_m2a_nlays)

       wspd_loc = coupler_search_name ('WSPD', coupler_m2a_vname, coupler_m2a_nvars, mta)
       if (wspd_loc < 0) then
          write (6, *) ' Abort: cannot find variable WSPD '
          stop
       end if

       zf_loc = coupler_search_name ('ZF', coupler_m2a_vname, coupler_m2a_nvars, mta)
       if (zf_loc < 0) then
          write (6, *) ' Abort: cannot find variable ZF '
          stop
       end if

       zh_loc = coupler_search_name ('ZH', coupler_m2a_vname, coupler_m2a_nvars, mta)
       if (zh_loc < 0) then
          write (6, *) ' Abort: cannot find variable ZH '
          stop
       end if

       cfrac3d_loc = coupler_search_name ('CFRAC_3D', coupler_m2a_vname, coupler_m2a_nvars, mta)
       if (cfrac3d_loc < 0) then
          write (6, *) ' Abort: cannot find variable CFRAC_3D '
          stop
       end if

       eddy_loc = coupler_search_name ('eddy', coupler_m2a_vname, coupler_m2a_nvars, mta)
       if (eddy_loc < 0) then
          write (6, *) ' Abort: cannot find variable eddy '
          stop
       end if

       temp_loc = coupler_search_name ('TA', coupler_m2a_vname, coupler_m2a_nvars, mta)
       if (temp_loc < 0) then
          write (6, *) ' Abort: cannot find variable TA '
          stop
       end if

       dens_loc = coupler_search_name ('DENS', coupler_m2a_vname, coupler_m2a_nvars, mta)
       if (dens_loc < 0) then
          write (6, *) ' Abort: cannot find variable DENS '
          stop
       end if

       qv_loc = coupler_search_name ('QV', coupler_m2a_vname, coupler_m2a_nvars, mta)
       if (qv_loc < 0) then
          write (6, *) ' Abort: cannot find variable QV '
          stop
       end if

       pres_loc = coupler_search_name ('PRES', coupler_m2a_vname, coupler_m2a_nvars, mta)
       if (pres_loc < 0) then
          write (6, *) ' Abort: cannot find variable PRES '
          stop
       end if

       qc_loc = coupler_search_name ('QC', coupler_m2a_vname, coupler_m2a_nvars, mta)
       if (qc_loc < 0) then
          write (6, *) ' Abort: cannot find variable QC '
          stop
       end if

       qr_loc = coupler_search_name ('QR', coupler_m2a_vname, coupler_m2a_nvars, mta)
       if (qr_loc < 0) then
          write (6, *) ' Abort: cannot find variable QR '
          stop
       end if

       qg_loc = coupler_search_name ('QG', coupler_m2a_vname, coupler_m2a_nvars, mta)
       if (qg_loc < 0) then
          write (6, *) ' Abort: cannot find variable QG '
          stop
       end if

       qi_loc = coupler_search_name ('QI', coupler_m2a_vname, coupler_m2a_nvars, mta)
       if (qi_loc < 0) then
          write (6, *) ' Abort: cannot find variable QI '
          stop
       end if

       qs_loc = coupler_search_name ('QS', coupler_m2a_vname, coupler_m2a_nvars, mta)
       if (qs_loc < 0) then
          write (6, *) ' Abort: cannot find variable QS '
          stop
       end if

       qc_cu_loc = coupler_search_name ('qc_cu', coupler_m2a_vname, coupler_m2a_nvars, mta)
       if (qc_cu_loc < 0) then
          write (6, *) ' Abort: cannot find variable qc_cu '
          stop
       end if

       qi_cu_loc = coupler_search_name ('qi_cu', coupler_m2a_vname, coupler_m2a_nvars, mta)
       if (qi_cu_loc < 0) then
          write (6, *) ' Abort: cannot find variable qi_cu '
          stop
       end if

       densa_j_loc = coupler_search_name ('DENSA_J', coupler_m2a_vname, coupler_m2a_nvars, mta)
       if (densa_j_loc < 0) then
          write (6, *) ' Abort: cannot find variable DENSA_J '
          stop
       end if

       ustar_loc = coupler_search_name ('USTAR', coupler_m2a_vname, coupler_m2a_nvars, mta)
       if (ustar_loc < 0) then
          write (6, *) ' Abort: cannot find variable USTAR '
          stop
       end if

       hfx_loc = coupler_search_name ('HFX', coupler_m2a_vname, coupler_m2a_nvars, mta)
       if (hfx_loc < 0) then
          write (6, *) ' Abort: cannot find variable HFX '
          stop
       end if

       pbl_loc = coupler_search_name ('PBL', coupler_m2a_vname, coupler_m2a_nvars, mta)
       if (pbl_loc < 0) then
          write (6, *) ' Abort: cannot find variable PBL '
          stop
       end if

    end if

    coupler_initialized = .true.

  end subroutine coupler_init

! -----------------------------------------------------------------------------
  subroutine coupler_1d_rdata_storing (vname, data)

    character (*), intent(in) :: vname
    real, intent(in)          :: data(:)

    integer :: list_loc, var_loc, var_type, transfer_type, tsize  ! tsize is need for MPAS where data storage is larger

    var_loc = coupler_search_name (vname, coupler_m2a_vname, coupler_m2a_nvars, mta, list_loc)
    if (var_loc < 0) then
       write (6, *) ' Abort in coupler storing: variable does not exist ', vname
       stop
    else
       transfer_type = mta
    end if

    if (coupled_model == wrf_cmaq) then
       var_type = coupler_m2a_data_type(list_loc)
    else
       if (transfer_type == mta) then
          tsize = size(coupler_m2a_dr_2d_data,1)
          coupler_m2a_dr_2d_data(:,1,var_loc,1) = data(1:tsize)
       end if
    end if

  end subroutine coupler_1d_rdata_storing

! -----------------------------------------------------------------------------
  subroutine coupler_1d_idata_storing (vname, data)

    character (*), intent(in) :: vname
    integer, intent(in)       :: data(:)

    integer :: list_loc, var_loc, var_type, transfer_type, tsize

    var_loc = coupler_search_name (vname, coupler_m2a_vname, coupler_m2a_nvars, mta, list_loc)
    if (var_loc < 0) then
       write (6, *) ' Abort in coupler storing: variable does not exist ', vname
       stop
    else
       transfer_type = mta
    end if

    if (coupled_model == wrf_cmaq) then
       var_type = coupler_m2a_data_type(list_loc)
    else
       if (transfer_type == mta) then
          tsize = size(coupler_m2a_dr_2d_data,1)
          coupler_m2a_dr_2d_data(:,1,var_loc,1) = real(data(1:tsize))
       end if
    end if

  end subroutine coupler_1d_idata_storing

! -----------------------------------------------------------------------------
  subroutine coupler_2d_rdata_storing (vname, data)

    character (*), intent(in) :: vname
    real, intent(in)          :: data(:,:)

    integer :: list_loc, var_loc, var_type, transfer_type, tsize1, tsize2, k, i, lhead

    var_loc = coupler_search_name (vname, coupler_m2a_vname, coupler_m2a_nvars, mta, list_loc)
    if (var_loc < 0) then
       write (6, *) ' Abort in coupler storing: variable does not exist ', vname
       stop
    else
       transfer_type = mta
    end if

    if (coupled_model == wrf_cmaq) then
       var_type = coupler_m2a_data_type(list_loc)
       if (var_type == ir2) then
          coupler_m2a_ir_2d_data(:,:,var_loc) = data 
       else if (var_type == ie2) then
          coupler_m2a_ie_2d_data(:,:,var_loc) = data 
       else if (var_type == dr2) then
          if (data_m2a_head(list_loc) == -1) then
             coupler_m2a_dr_2d_data(:,:,var_loc,0) = data 
             data_m2a_head(list_loc) = 0
          else if (data_m2a_tail(list_loc) == -1) then
             coupler_m2a_dr_2d_data(:,:,var_loc,1) = data 
             data_m2a_tail(list_loc) = 1
          else
             lhead = data_m2a_head(list_loc)
             coupler_m2a_dr_2d_data(:,:,var_loc,lhead) = data 
          end if
       else
          write (6, *) ' Abort in coupler_2d_rdata_storing: Unknown data type'
          stop
       end if
    else
       if (transfer_type == mta) then
          tsize1 = size(coupler_m2a_dr_3d_data,1)
          tsize2 = size(coupler_m2a_dr_3d_data,3)
          do k = 1, tsize2
             do i = 1, tsize1
                coupler_m2a_dr_3d_data(i,1,k,var_loc,1) = data(k, i)
             end do
          end do
       end if
    end if

  end subroutine coupler_2d_rdata_storing

! -----------------------------------------------------------------------------
  subroutine coupler_3d_rdata_storing (vname, data)

    character (*), intent(in) :: vname
    real, intent(in)          :: data(:,:,:)

    integer :: list_loc, var_loc, var_type, transfer_type, lhead

    var_loc = coupler_search_name (vname, coupler_m2a_vname, coupler_m2a_nvars, mta, list_loc)
    if (var_loc < 0) then
       var_loc = coupler_search_name (vname, coupler_a2m_vname, coupler_a2m_nvars, atm, list_loc)
       if (var_loc < 0) then
          write (6, *) ' Abort in coupler storing: variable does not exist ', vname
          stop
       else
          transfer_type = atm
       end if
    else
       transfer_type = mta
    end if

    if (coupled_model == wrf_cmaq) then
       if (transfer_type == mta) then
          var_type = coupler_m2a_data_type(list_loc)
          if (var_type == dr3) then
             if (data_m2a_head(list_loc) == -1) then
                coupler_m2a_dr_3d_data(:,:,:,var_loc,0) = data 
                data_m2a_head(list_loc) = 0
             else if (data_m2a_tail(list_loc) == -1) then
                coupler_m2a_dr_3d_data(:,:,:,var_loc,1) = data 
                data_m2a_tail(list_loc) = 1
             else
                lhead = data_m2a_head(list_loc)
                coupler_m2a_dr_3d_data(:,:,:,var_loc,lhead) = data 
             end if
          else if (var_type == de3) then
             if (data_m2a_head(list_loc) == -1) then
                coupler_m2a_de_3d_data(:,:,:,var_loc,0) = data 
                data_m2a_head(list_loc) = 0
             else if (data_m2a_tail(list_loc) == -1) then
                coupler_m2a_de_3d_data(:,:,:,var_loc,1) = data 
                data_m2a_tail(list_loc) = 1
             else
                lhead = data_m2a_head(list_loc)
                coupler_m2a_de_3d_data(:,:,:,var_loc,lhead) = data 
             end if
          else
             write (6, *) ' Abort in coupler_3d_rdata_storing: Unknown data type'
             stop
          end if
       else
          coupler_a2m_3d_data(:,:,:,list_loc) = data 
       end if
    else
       if (transfer_type == mta) then
          coupler_m2a_dr_3d_data(:,1,:,var_loc,1) = data(:,1,:)
       end if
    end if

  end subroutine coupler_3d_rdata_storing

! -----------------------------------------------------------------------------
  subroutine coupler_1d_rdata_retrieving (vname, data)

    character (*), intent(in) :: vname
    real, intent(out)         :: data(:)

    integer :: list_loc, var_loc, var_type, transfer_type, lhead, dsize, d1, d2, d3

    var_loc = coupler_search_name (vname, coupler_m2a_vname, coupler_m2a_nvars, mta, list_loc)
    if (var_loc < 0) then
       write (6, *) ' Abort in coupler retrieving: variable does not exist ', vname
       stop
    else
       transfer_type = mta
    end if

    if (coupled_model == wrf_cmaq) then
       var_type = coupler_m2a_data_type(list_loc)
       lhead = data_m2a_head(list_loc)
       dsize = size(data)
       if (var_type == dr2) then
          d1    = size(coupler_m2a_dr_2d_data, 1)
          d2    = size(coupler_m2a_dr_2d_data, 2)
          data = reshape (coupler_m2a_dr_2d_data(2:d1-1,2:d2-1,var_loc,lhead), (/ dsize /))
       else if (var_type == dr3) then
          data = reshape (coupler_m2a_dr_3d_data(:,:,:,var_loc,lhead), (/ dsize /))
       else if (var_type == de3) then
          d1    = size(coupler_m2a_de_3d_data, 1)
          d2    = size(coupler_m2a_de_3d_data, 2)
          data = reshape (coupler_m2a_de_3d_data(2:d1-1,2:d2-1,:,var_loc,lhead), (/ dsize /))
       else
          write (6, *) ' Abort in coupler_1d_rdata_retrieving: Unknown data type'
          stop
       end if
       data_m2a_head(list_loc) = data_m2a_tail(list_loc)
       data_m2a_tail(list_loc) = lhead
    else
       if (transfer_type == mta) then
          data = coupler_m2a_dr_2d_data(:,1,var_loc,1)
       end if
    end if

  end subroutine coupler_1d_rdata_retrieving

! -----------------------------------------------------------------------------
  subroutine coupler_1d_idata_retrieving (vname, data)

    character (*), intent(in) :: vname
    integer, intent(out)      :: data(:)

    integer :: list_loc, var_loc, var_type, transfer_type

    var_loc = coupler_search_name (vname, coupler_m2a_vname, coupler_m2a_nvars, mta, list_loc)
    if (var_loc < 0) then
       write (6, *) ' Abort in coupler retrieving: variable does not exist ', vname
       stop
    else
       transfer_type = mta
    end if

    if (coupled_model == wrf_cmaq) then
       var_type = coupler_m2a_data_type(list_loc)
    else
       if (transfer_type == mta) then
          data = int(coupler_m2a_dr_2d_data(:,1,var_loc,1))
       end if
    end if

  end subroutine coupler_1d_idata_retrieving

! -----------------------------------------------------------------------------
  subroutine coupler_2d_rdata_retrieving (vname, data)

    character (*), intent(in) :: vname
    real, intent(out)         :: data(:,:)

    integer :: list_loc, var_loc, var_type, transfer_type, lhead, d1, d2

    var_loc = coupler_search_name (vname, coupler_m2a_vname, coupler_m2a_nvars, mta, list_loc)
    if (var_loc < 0) then
       write (6, *) ' Abort in coupler retrieving: variable does not exist ', vname
       stop
    else
       transfer_type = mta
    end if

    if (coupled_model == wrf_cmaq) then
       var_type = coupler_m2a_data_type(list_loc)
       if (var_type == ir2) then
          data = coupler_m2a_ir_2d_data(:,:,var_loc)
       else if (var_type == ie2) then
          data = coupler_m2a_ie_2d_data(:,:,var_loc)
       else if (var_type == dr2) then
          lhead = data_m2a_head(list_loc)
          d1 = size(coupler_m2a_dr_2d_data,1)
          d2 = size(coupler_m2a_dr_2d_data,2)
          data = coupler_m2a_dr_2d_data(2:d1-1,2:d2-1,var_loc,lhead)
          data_m2a_head(list_loc) = data_m2a_tail(list_loc)
          data_m2a_tail(list_loc) = lhead
       else
          write (6, *) ' Abort in coupler_2d_rdata_retrieving: Unknown data type'
          stop
       end if
    else
       if (transfer_type == mta) then
          data(:,1) = coupler_m2a_dr_2d_data(:,1,var_loc,1)
       end if
    end if

  end subroutine coupler_2d_rdata_retrieving

! -----------------------------------------------------------------------------
  subroutine coupler_2d_idata_retrieving (vname, data)

    character (*), intent(in) :: vname
    integer, intent(out)      :: data(:,:)

    integer :: list_loc, var_loc, var_type, transfer_type

    var_loc = coupler_search_name (vname, coupler_m2a_vname, coupler_m2a_nvars, mta, list_loc)
    if (var_loc < 0) then
       write (6, *) ' Abort in coupler retrieving: variable does not exist ', vname
       stop
    else
       transfer_type = mta
    end if

    if (coupled_model == wrf_cmaq) then
       var_type = coupler_m2a_data_type(list_loc)
    else
       if (transfer_type == mta) then
          data(:,1) = int(coupler_m2a_dr_2d_data(:,1,var_loc,1))
       end if
    end if

  end subroutine coupler_2d_idata_retrieving

! -----------------------------------------------------------------------------
  subroutine coupler_3d_rdata_retrieving (vname, data)

    character (*), intent(in) :: vname
    real, intent(out)         :: data(:,:,:)

    integer :: list_loc, var_loc, var_type, transfer_type

    var_loc = coupler_search_name (vname, coupler_m2a_vname, coupler_m2a_nvars, mta, list_loc)
    if (var_loc < 0) then
       var_loc = coupler_search_name (vname, coupler_a2m_vname, coupler_a2m_nvars, atm, list_loc)
       if (var_loc < 0) then
          write (6, *) ' Abort in coupler retrieving: variable does not exist ', vname
          stop
       else
          transfer_type = atm
       end if
    else
       transfer_type = mta
    end if

    if (coupled_model == wrf_cmaq) then
       if (transfer_type == atm) then
          data = coupler_a2m_3d_data(:,:,:,list_loc)
       else
          var_type = coupler_m2a_data_type(list_loc)
       end if
    else
       if (transfer_type == mta) then
          data(:,1,:) = coupler_m2a_dr_3d_data(:,1,:,var_loc,1)
       end if
    end if

  end subroutine coupler_3d_rdata_retrieving

! -----------------------------------------------------------------------------
  subroutine find_n_items_in_namelist (file, n, so4_ind)

    character (*), intent(in)      :: file
    integer, intent(out)           :: n
    integer, intent(out), optional :: so4_ind(3)

    integer, parameter :: iunit = 7
    integer :: loc_n, stat, n_so4
    logical :: eof
    character (200) :: line

    open (unit = iunit, file=file, status='old')

    n_so4 = 0
    loc_n = 0
    eof = .false.
    do while (.not. eof)
       read (iunit, '(a200)', iostat=stat) line
       if (stat .ne. 0) then
          eof = .true.
       else
          if (line(1:1) == "'") then
! for 5.4: AE speicies does not provide full name in the *.nml file and only 
!          indicates under Aitken, Accum, and Coarse modes
             if (present(so4_ind)) then
                if (line(2:5) == 'ASO4') then
                   if (line(21:21) == 'T') then
                      loc_n = loc_n + 1
                      n_so4 = n_so4 + 1
                      so4_ind(n_so4) = loc_n
                   end if
                   if (line(29:29) == 'T') then
                      loc_n = loc_n + 1
                      n_so4 = n_so4 + 1
                      so4_ind(n_so4) = loc_n
                   end if
                   if (line(36:36) == 'T') then
                      loc_n = loc_n + 1
                      n_so4 = n_so4 + 1
                      so4_ind(n_so4) = loc_n
                   end if
                else
                   if (line(21:21) == 'T') loc_n = loc_n + 1
                   if (line(29:29) == 'T') loc_n = loc_n + 1
                   if (line(36:36) == 'T') loc_n = loc_n + 1
                end if
             else
                loc_n = loc_n + 1
             end if
          end if
       end if
    end do

    n = loc_n

    close (iunit)

  end subroutine find_n_items_in_namelist

! -------------------------------------------------------------------------
  recursive subroutine quicksort_vname (name, dtype, loc, begin, end, nlays)

    character (*), intent(out)       :: name(:)
    integer, intent(inout)           :: dtype(:)
    integer, intent(inout), optional :: nlays(:)
    integer, intent(inout)           :: loc(:)
    integer, intent(in)              :: begin, end

    integer        :: i, j, temp
    character (50) :: str1, str2
    logical        :: done

    str1 = name( (begin + end) / 2 )
    i = begin
    j = end
    done = .false.
    do while (.not. done)
       do while (name(i) < str1)
          i = i + 1
       end do
       do while (str1 < name(j))
          j = j - 1
       end do
       if (i .ge. j) then
          done = .true.
       else
          str2 = name(i)
          name(i) = name(j)
          name(j) = str2

          temp     = dtype(i)
          dtype(i) = dtype(j)
          dtype(j) = temp

          if (present(nlays)) then
             temp     = nlays(i)
             nlays(i) = nlays(j)
             nlays(j) = temp
          end if

          temp   = loc(i)
          loc(i) = loc(j)
          loc(j) = temp

          i = i + 1
          j = j - 1
       end if
    end do
    if (present(nlays)) then
       if (begin < i-1) call quicksort_vname(name, dtype, loc, begin, i-1, nlays)
       if (j+1 < end)   call quicksort_vname(name, dtype, loc, j+1, end, nlays)
    else
       if (begin < i-1) call quicksort_vname(name, dtype, loc, begin, i-1)
       if (j+1 < end)   call quicksort_vname(name, dtype, loc, j+1, end)
    end if

  end subroutine quicksort_vname

! -------------------------------------------------------------------------
  recursive integer function binary_search (vname, vname_list, lower, upper) result (pos)

    character (*), intent(in) :: vname
    character (*), intent(in) :: vname_list(:)
    integer, intent(in)       :: lower, upper 

    integer :: loc_pos, middle

    if (lower > upper) then
       loc_pos = -1
    else
       middle = (lower + upper) / 2
       if (vname == vname_list(middle)) then
          loc_pos = middle
       else if (vname > vname_list(middle)) then
          loc_pos = binary_search (vname, vname_list, middle+1, upper)
       else
          loc_pos = binary_search (vname, vname_list, lower, middle-1)
       end if
    end if

    pos = loc_pos

  end function binary_search

! -------------------------------------------------------------------------
  integer function coupler_search_name (vname, vname_list, nvars, transfer_type, list_loc)

    character (*), intent(in)      :: vname
    character (*), intent(in)      :: vname_list(:)
    integer, intent(in)            :: nvars, transfer_type
    integer, intent(out), optional :: list_loc

    integer :: temp_loc
    character (vname_max_str_len) :: loc_vname

    if (transfer_type == mta) then
       loc_vname = 'c_' // vname
    else
       loc_vname = vname
    end if

    temp_loc = binary_search (loc_vname, vname_list, 1, nvars)

    if (temp_loc > 0) then
       if (transfer_type == mta) then
          coupler_search_name = coupler_m2a_loc(temp_loc)
       else
          coupler_search_name = coupler_a2m_loc(temp_loc)
       end if
    else
       coupler_search_name = -1
    end if

    if (present(list_loc)) then
       list_loc = temp_loc
    end if

  end function coupler_search_name

end module coupler_module
