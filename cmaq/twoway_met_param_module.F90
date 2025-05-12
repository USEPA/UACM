!===============================================================================
! Purpose:  Define general information in each of the buffer files
!
! Revised:  April 2007  Original version.  David Wong
!           July, 16 2013  David Wong -- corrected the unit for RA and RS
!                                     -- added DLUSE in GRIDCRO2D
!           Jan, 11 2016  David Wong -- added a new variable PV
!           Mar 04, 2019  Gilliam and Wong -- added new metcro2d variables 
!                         according to new PX implementation in WRFv4.1+
!           Aug 01, 2019  Wong -- added two new metdot3d variables, UWIND and
!                                 VWIND (wind component of the mass point)
!           Jan 10, 2024  David Wong -- Incorporated unified coupler implmentation
!===============================================================================

  module twoway_met_param_module

    use coupler_module, only : n_gridcro2d_var, &
                               n_griddot2d_var, &
                               n_metcro3d_var,  &
                               n_metdot3d_var,  &
                               n_metcro2d_var

    INTEGER, PARAMETER :: max_nvars       = 1000

end module twoway_met_param_module
