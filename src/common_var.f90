!================================================================================
module common_var
! 
  integer :: mpisize, mpirank, mpierror
  double precision :: dt, tmax, xmin, xmax, vmin, vmax, dx, dv, kT, x0now, dtint_prop
  integer :: nx,nt,nttot, dtmult, n_tprop, prop_order
  integer, allocatable :: ibeg_tprop(:),iend_tprop(:)
  double precision :: init_taug, opt_temp1,opt_temp2, target_acc
  integer :: type_Langevin, ratio_Langevin_MD, ntraj_Langevin, ntraj_MD
  integer :: fit_F, fit_gamm, fit_taug, fit_mass, fix_mass0, weight_nll
  integer :: iopt,opt_niter, type_error, pos_dep_gamma,pos_dep_mass, use_velocity, ntraj_prop, bin_index
  logical :: print_traj,test_propagator
  double precision, allocatable :: colvar(:,:), x0(:)
  double precision, allocatable :: Pref(:,:,:), Pmod(:,:,:)
  double precision :: mass0,gamm,taug, bin_width
  double precision :: dxgrid,dxgrid2
  integer, parameter :: ngrid=1000 ! this is arbitrary...
  integer, parameter :: nbins=20 ! this is arbitrary...
  double precision :: prof_F(ngrid),prof_force(ngrid),prof_g(ngrid),prof_m(ngrid)
  double precision :: max_Gaussian_h(3),max_Gaussian_w(3)
  double precision, DIMENSION(20) :: inv_hist_traj, hist_traj
  double precision, dimension(20) :: lower_bound, upper_bound
  double precision :: term1, term2
!
end module common_var
!================================================================================
