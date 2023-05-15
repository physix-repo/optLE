!================================================================================
subroutine read_input
!
  use common_var
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! README_begin
!
!  The input file should be called "input",
!  and inside it, the order of the lines is fixed.
!
!  input file format: keyword value (fixed order, keyword is to check)
!
!colvar_file	   colvar        (format: t x, multiple trajectoires appended)	
!dt        	     0.01        (time step in colvar file)
!dtmult    	     10          (time step for integration dtint = dt/dtmult)
!xmin            -3.0            (better to take xmin and xmax large)
!xmax             3.0
!xbins	          40             (used for x and v histograms)
!kT               5.0
!init_tau         400.
!opt_niter        500            (optimization steps, 0 = just print one traj and stop)
!opt_temp         1e-4 1e-5 0.05    (annealing T: initial, final, target acceptance)
!type_Langevin      1            (0 = overdamped, 1 = std, 2 = GLE)
!ratio_Langevin_MD  10           (ratio between number of Langevin traj and MD traj)
!fit_F              1
!fit_gamma          1
!fit_tau            0
!fit_mass           0
!type_error         1            (1 = -logLik(KL), 2 = RMSD, 3 = -logLik(analyt.propagator (-3=test&exit)), 4 = -logLik(num.propagator))
!pos_dep_gamma      1            (0 = fixed gamma, 1 = position-dependent gamma)
!pos_dep_mass       1            (0 = fixed mass , 1 = position-dependent mass )
!max_Gaussian_h   10.  5.   1.   (max height of Gaussians added to profiles F,g,m)
!max_Gaussian_w   0.05 0.05 0.05 (max width of Gaussians added, units of (xmax-xmin))
!fix_mass0          1            (keep fixed the value of the mass at the TS)
!use_velocity       1            (0 = optimize P(x,t), 1 = optimize P(x,v,t)
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! README_end
  implicit none
  character :: keyword*20
  character :: line*100
  double precision :: tmp1,tmp2,x
  integer :: i,j,it
  !
  !----------------------------------------------- read input file
  write(*,*) ""
  write(*,*) "----------------------- parsing input file ------------------------"
  write(*,*) ""
  open(inputfile_id,file="input",status="old")
  !
  read(inputfile_id,*) keyword,colvar_file
  if (trim(keyword)/="colvar_file")    call error("input: expected keyword colvar_file")
  write(*,*) keyword,colvar_file
  !
  read(inputfile_id,*) keyword,dt
  if (trim(keyword)/="dt")             call error("input: expected keyword dt")
  write(*,*) keyword,dt
  !
  read(inputfile_id,*) keyword,dtmult
  if (trim(keyword)/="dtmult")         call error("input: expected keyword dtmult")
  write(*,*) keyword,dtmult
  !
  read(inputfile_id,*) keyword,xmin
  if (trim(keyword)/="xmin")           call error("input: expected keyword xmin")
  write(*,*) keyword,xmin
  !
  read(inputfile_id,*) keyword,xmax
  if (trim(keyword)/="xmax")           call error("input: expected keyword xmax")
  write(*,*) keyword,xmax
  !
  read(inputfile_id,*) keyword,nx
  if (trim(keyword)/="xbins")          call error("input: expected keyword xbins")
  write(*,*) keyword,nx
  !
  read(inputfile_id,*) keyword,kT
  if (trim(keyword)/="kT")             call error("input: expected keyword kT")
  write(*,*) keyword,kT
  !
  read(inputfile_id,*) keyword,init_taug
  if (trim(keyword)/="init_tau")     call error("input: expected keyword init_tau")
  write(*,*) keyword,init_taug
  !
  read(inputfile_id,*) keyword,opt_niter
  if (trim(keyword)/="opt_niter") call error("input: expected keyword opt_niter")
  write(*,*) keyword,opt_niter
  !
  read(inputfile_id,*) keyword,opt_temp1,opt_temp2,target_acc
  if (trim(keyword)/="opt_temp") call error("input: expected keyword opt_temp")
  write(*,*) keyword,opt_temp1,opt_temp2,target_acc
  !
  read(inputfile_id,*) keyword,type_Langevin
  if (trim(keyword)/="type_Langevin")  call error("input: expected keyword type_Langevin")
  write(*,*) keyword,type_Langevin
  !
  read(inputfile_id,*) keyword,ratio_Langevin_MD
  if (trim(keyword)/="ratio_Langevin_MD") call error("input: expected keyword ratio_Langevin_MD")
  write(*,*) keyword,ratio_Langevin_MD
  !
  read(inputfile_id,*) keyword,fit_F
  if (trim(keyword)/="fit_F") call error("input: expected keyword fit_F")
  write(*,*) keyword,fit_F
  !
  read(inputfile_id,*) keyword,fit_gamm
  if (trim(keyword)/="fit_gamma") call error("input: expected keyword fit_gamma")
  write(*,*) keyword,fit_gamm
  !
  read(inputfile_id,*) keyword,fit_taug
  if (trim(keyword)/="fit_tau") call error("input: expected keyword fit_tau")
  write(*,*) keyword,fit_taug
  !
  read(inputfile_id,*) keyword,fit_mass
  if (trim(keyword)/="fit_mass") call error("input: expected keyword fit_mass")
  write(*,*) keyword,fit_mass
  !
  read(inputfile_id,*) keyword,type_error
  if (trim(keyword)/="type_error") call error("input: expected keyword type_error")
  write(*,*) keyword,type_error
  !
  read(inputfile_id,*) keyword,pos_dep_gamma
  if (trim(keyword)/="pos_dep_gamma") call error("input: expected keyword pos_dep_gamma")
  write(*,*) keyword,pos_dep_gamma
  !
  read(inputfile_id,*) keyword,pos_dep_mass
  if (trim(keyword)/="pos_dep_mass") call error("input: expected keyword pos_dep_mass")
  write(*,*) keyword,pos_dep_mass
  !
  read(inputfile_id,*) keyword,max_Gaussian_h(1:3)
  if (trim(keyword)/="max_Gaussian_h") call error("input: expected keyword max_Gaussian_h")
  write(*,*) keyword,max_Gaussian_h(1:3)
  !
  read(inputfile_id,*) keyword,max_Gaussian_w(1:3)
  if (trim(keyword)/="max_Gaussian_w") call error("input: expected keyword max_Gaussian_w")
  write(*,*) keyword,max_Gaussian_w(1:3)
  !
  read(inputfile_id,*) keyword,fix_mass0
  if (trim(keyword)/="fix_mass0") call error("input: expected keyword fix_mass0")
  write(*,*) keyword,pos_dep_mass
  !
  read(inputfile_id,*) keyword,use_velocity
  if (trim(keyword)/="use_velocity") call error("input: expected keyword use_velocity")
  write(*,*) keyword,use_velocity
  !
  read(inputfile_id,*) keyword,dtint_prop
  if (trim(keyword)/="dtint_prop") call error("input: expected keyword dtint_prop")
  write(*,*) keyword,dtint_prop
  !
  read(inputfile_id,*) keyword,ntraj_prop
  if (trim(keyword)/="ntraj_prop") call error("input: expected keyword ntraj_prop")
  write(*,*) keyword,ntraj_prop
  !
  close(inputfile_id)
  !
  write(*,*) ""
  if (type_Langevin.eq.0) then
    write(*,*) "type_Langevin = 0: using overdamped Langevin equation"
  elseif (type_Langevin.eq.1) then
    write(*,*) "type_Langevin = 1: using standard Langevin equation (memoryless)"
    if (fit_taug.eq.1) then
      write(*,*) "(overriding fit_tau=1: setting fit_tau=0)"
      fit_taug=0
    endif
  elseif (type_Langevin.eq.2) then
    write(*,*) "type_Langevin = 2: using generalized Langevin equation (exponential memory)"
  else  
    call error("only type_Langevin = 0 1 2 are implemented")
  endif
  !
  write(*,*) ""
  if (pos_dep_gamma.eq.0) then
    write(*,*) "pos_dep_gamma = 0: using constant gamma (indepedent of x)"
  elseif (pos_dep_gamma.eq.1) then
    write(*,*) "pos_dep_gamma = 1: using position-dependent gamma"
  else  
    call error("only pos_dep_gamma = 0 1 are implemented")
  endif
  !
  write(*,*) ""
  if (pos_dep_mass.eq.0) then
    write(*,*) "pos_dep_mass = 0: using constant mass (indepedent of x)"
  elseif (pos_dep_mass.eq.1) then
    write(*,*) "pos_dep_mass = 1: using position-dependent mass"
  else  
    call error("only pos_dep_mass = 0 1 are implemented")
  endif
  !
  write(*,*) ""
  !
  if (fit_F   .ne.0.and.fit_F   .ne.1) call error("fit_F     must be 0 or 1")
  if (fit_gamm.ne.0.and.fit_gamm.ne.1) call error("fit_gamma must be 0 or 1")
  if (fit_taug.ne.0.and.fit_taug.ne.1) call error("fit_tau   must be 0 or 1")
  if (fit_mass.ne.0.and.fit_mass.ne.1) call error("fit_mass  must be 0 or 1")
  !
  test_propagator=.false.
  write(*,*) ""
  if (type_error.eq.1) then
    write(*,*) "type_error = 1: using -log(Likelihood) (equiv. to Kullback-Leibler divergence)"
  elseif (type_error.eq.2) then
    write(*,*) "type_error = 2: using RMSD of prob. distributions"
  elseif (type_error.eq.3.or.type_error.eq.-3) then
    write(*,*) "type_error = 3: using -log(Likelihood) (from analytical propagator)"
    if (type_error.eq.-3) then
      test_propagator=.true.
      type_error=3
      write(*,*) "WARNING: testing the analytical propagator against the numerical one and exiting !"
    endif
    if (type_Langevin.gt.1) then
      call error("this error is implemented only for overdamped and standard Langevin dynamics")
    endif
  else  
    call error("only type_error = 1 2 3 4 are implemented")
  endif
  !
  write(*,*) ""
  if (use_velocity.eq.1) then
    write(*,*) "optimization of p(x,v,t)"
  else
    write(*,*) "optimization of p(x,t) : velocity is not used"
  endif
  !
  write(*,*) ""
  if (target_acc.gt.0.d0.and.target_acc.lt.1.d0) then
    write(*,*) "adapting opt_temp at each step to target acceptance =",target_acc
  else
    write(*,*) "opt_temp will be linearly changed between initial and final one"
  endif
  !
  dxgrid=(xmax-xmin)/dble(ngrid-1) 
  dxgrid2=2.d0*dxgrid
  !
  !----------------------------- read colvar_file
  write(*,*) ""
  write(*,*) "----------------------- reading colvar_file -----------------------"
  write(*,*) ""
  open(intraj_id,file=colvar_file,status="old")
  tmax=-1.d0
  i=0
  ntraj_MD=0
  do
    read(intraj_id,'(A100)',end=102) line
    i=i+1
    if (line(1:1).ne."#") then
      read (line(:),*,err=101) tmp1,tmp2
      if (dabs(tmp1).lt.1.d-7) ntraj_MD=ntraj_MD+1
      if (tmp1.gt.tmax) tmax=tmp1
    endif
  enddo
  101 continue
      write(*,*) "problem reading file at line ",i
      call error("format of colvar_file: t x (two columns)") 
  102 continue
  rewind(intraj_id)
  allocate(x0(ntraj_MD)) ! here we store the shooting points (they can be different)
  nttot=i                ! total lines in colvar, i.e. total number of time frames
  nt=nint(tmax/dt)+1     ! number of time frames per trajectory, assuming equal durations 
  allocate(colvar(nttot,3)) ! time q dq/dt
  it=0
  do j=1,nttot
    read(intraj_id,*) colvar(j,1),colvar(j,2)
    if (colvar(j,1).lt.dt/10.) then ! fill x0 with points x(t=0)
      ! note: the syntax in if() here above is a complicated way of testing if t=0 ...
      it=it+1
      x0(it)=colvar(j,2)
!      write(*,*) "x0=",x0(it) ! DEBUG
    endif
    ! TODO can we allow points outside, skipping them?
    if (colvar(j,2).gt.xmax) call error("position larger than xmax !")
    if (colvar(j,2).lt.xmin) call error("position smaller than xmin !")
  enddo
  close(intraj_id)
  ! compute numerical velocity dq/dt
  write(111,'(A)') "# time position numerical_velocity"
  colvar(:,3)=0.d0
  do j=2,nttot-1
    if (colvar(j+1,1).gt.colvar(j-1,1)) then ! avoid time discontinuity
      colvar(j,3) = ( colvar(j+1,2)-colvar(j-1,2) )/(2.d0*dt)
    endif
    open(rewrite_intraj_id, file=trim(colvar_file)//".trajectories")
    write(rewrite_intraj_id,'(F18.8,2F18.10)') colvar(j,1),colvar(j,2),colvar(j,3)
    close(rewrite_intraj_id)
  enddo
  write(*,*) "written numerical velocities in fort.111"
  ! store indeces of q,dq/dt points for propagator estimation
  ! note: when using type_err=3 (propagator-based likelihood)
  !       we use q,dq/dt skipping every dtmult, where
  !       i_prop contains the starting point of every trajectory
  allocate(ibeg_tprop(nttot/dtmult+1),iend_tprop(nttot/dtmult+1))
  n_tprop=1
  ibeg_tprop(1)=2
  iend_tprop(1)=2
  write(*,*) "lines at beginning and end of trajectories with stride dtmult:"
  do j=3,nttot-1
    if (colvar(j,1).lt.colvar(j-1,1)) then
      n_tprop=n_tprop+1
      ibeg_tprop(n_tprop)=j+1 ! skip t=0 to allow calculation of velocity
      iend_tprop(n_tprop)=j+1 
    else
      ! skip last time of each traj to allow calculation of velocity
      if (((j-iend_tprop(n_tprop)).eq.dtmult).and.(colvar(j+1,1).gt.colvar(j,1))) then
        iend_tprop(n_tprop)=j
      endif
    endif
  enddo
  do j=1,n_tprop
    if ((j.le.5).or.(j.ge.n_tprop-4)) then
      write(*,'(i6,i9,i9)') j,ibeg_tprop(j),iend_tprop(j)
    else
      if (j.eq.6) write(*,*) "..."
    endif
  enddo
  !
  ntraj_Langevin=ntraj_MD*ratio_Langevin_MD
  write(*,'(A,I9,I6,F12.3)') " total time frames, time frames per traj, tmax = ",nttot,nt,tmax
  write(*,'(A,I9)')          " number of MD shooting trajectories = ",ntraj_MD
!  if (ntraj_MD.ne.nint((1.d0*nttot)/nt)) call error("number of MD traj must be equal to total frames / frames per traj")
  if (ntraj_MD.ne.nint((1.d0*nttot)/nt)) write(*,*) "WARNING: number of MD traj must be equal to total frames / frames per traj"
  !----------------------------- read RESTART
  ! if (restart.eq.1) then
    write(*,*) ""
    write(*,*) "----------------------- reading RESTART ---------------------------"
    write(*,*) ""
    !
    write(*,*) "reading F(x) gamma(x) mass(x) profiles from file RESTART"
    open(restart_id,file="RESTART",status="old") ! TODO: check for errors in format and ngrid...
    read(restart_id,*) line
    do i=1,ngrid
      read(restart_id,*) x,prof_F(i),tmp1,prof_g(i),prof_m(i) ! x,F(with units),gamma,mass
      if (i.eq.1    .and.abs(x-xmin)>dxgrid) call error("xmin is different from first position in RESTART")
      if (i.eq.ngrid.and.abs(x-xmax)>dxgrid) call error("xmax is different from last  position in RESTART")
    enddo
    close(restart_id) 
    !
    call update_prof_force ! important, every time prof_F is changed!
    !
    if (fix_mass0.eq.1) then
      mass0=prof_m(int((x0(1)-xmin)/dxgrid)+1)
      write(*,*) "  mass0 (at initial point of first MD trajectory) =",mass0 
    endif
    !
  ! endif
!
end subroutine read_input
!================================================================================
subroutine init_Pref
!
  use common_var
  !
  implicit none
  integer i,ix,iv,it,nout
  double precision :: v

  !
  allocate(Pref(nx,nx,0:nt),Pmod(nx,nx,0:nt))
  !
  write(*,*) ""
  write(*,*) "----------------------- computing Pref(x,v,t) from colvar_file ------"
  write(*,*) ""
  nout=0
  Pref=0.D0
  dx=(xmax-xmin)/dble(nx-1)
  !
  ! determining velocity range:
  vmin= 1.d15
  vmax=-1.d15
  do i=1,nttot
    it=nint(dble(nt)*colvar(i,1)/tmax)
    if (it<nt) then
      v=((colvar(i+1,2)-colvar(i,2))/dt) 
      if (v>vmax) vmax=v
      if (v<vmin) vmin=v
    endif
  enddo
  vmin=vmin*1.5D0 ! TODO: this is arbitrary !!
  vmax=vmax*1.5D0 ! TODO: this is arbitrary !!
  dv=(vmax-vmin)/dble(nx-1)
  write(*,'(A,2E16.6)') " setting vmin vmax =",vmin,vmax
  !
  ! computing Pref:
  do i=1,nttot
    it=nint(dble(nt)*colvar(i,1)/tmax)
    if (it<0.or.it>nt) then
      write(*,*) "at line ",i," of colvar_file time is out of bounds:",colvar(i,1),it
      call error("error in colvar_file")
    endif
    !
    ix=int((colvar(i,2)-xmin)/dx)+1
    if (use_velocity==1) then
      if (it<nt) then
        v=((colvar(i+1,2)-colvar(i,2))/dt) 
        iv=int((v-vmin)/dv)+1
      else
        iv=-1
      endif
    else
      iv=1
    endif
    !
    if (ix<1.or.ix>nx.or.iv<1.or.iv>nx) then
      nout=nout+1
    else
      Pref(ix,iv,it)=Pref(ix,iv,it)+1.d0
    endif
  enddo
  !
  write(*,*) "points out of x,v,t bounds =",nout
  if (nout.eq.nttot) call error("all points are out of bounds !")
  !
  write(*,*) "normalizing Pref by the number of initial points =",sum(Pref(:,:,0))
  Pref=Pref/sum(Pref(:,:,0))
  !
#ifndef DEBUG
  if (type_error.ne.3) then
    write(*,*) "writing Pref file"
    open(pref_id,file=trim(colvar_file)//".Pref",status="unknown")
    do ix=1,nx
      do iv=1,nx
        do it=0,nt
          write(pref_id,'(4E11.3)') xmin+(dble(ix)-1.d0)*dx,vmin+(dble(iv)-1.d0)*dv,dble(it)*dt,Pref(ix,iv,it)
        enddo
      enddo
    enddo
    close(pref_id)
  endif
#endif
!
end subroutine init_Pref
!================================================================================
!   subroutine init_Pref_onlyxt
!   !
!     use common_var
!     !
!     implicit none
!     integer i,ix,it,nout
!     !
!     allocate(Pref(nx,0:nt),Pmod(nx,0:nt))
!     !
!     write(*,*) ""
!     write(*,*) "----------------------- computing Pref(x,t) from colvar_file ------"
!     write(*,*) ""
!     nout=0
!     Pref=0.D0
!     dx=(xmax-xmin)/dble(nx-1)
!     do i=1,nttot
!       it=nint(dble(nt)*colvar(i,1)/tmax)
!       if (it<0.or.it>nt) then
!         write(*,*) "at line ",i," of colvar_file time is out of bounds:",colvar(i,1),it
!         call error("error in colvar_file")
!       endif
!       ix=int((colvar(i,2)-xmin)/dx)+1
!       if (ix<1.or.ix>nx) then
!         nout=nout+1
!       else
!         Pref(ix,it)=Pref(ix,it)+1.d0
!       endif
!     enddo
!     !
!     write(*,*) "points out of x,t bounds =",nout
!     if (nout.eq.nttot) call error("all points are out of bounds !")
!     write(*,*) "normalizing Pref by the number of initial points =",sum(Pref(:,0))
!     Pref=Pref/sum(Pref(:,0))
!     !
!     write(*,*) "writing Pref file"
!     open(pref_id,file="Pref",status="unknown")
!     do ix=1,nx
!       do it=0,nt
!         write(pref_id,'(3F15.6)') xmin+(dble(ix)-1.d0)*dx,dble(it)*dt,Pref(ix,it)
!       enddo
!       write(pref_id,*) ""
!     enddo
!     close(pref_id)
!   !
!   end subroutine init_Pref_onlyxt
!================================================================================
!   subroutine init_Qref
!   !
!     use common_var
!     !
!     implicit none
!     integer :: i,iv,it
!     double precision :: v
!     !
!     allocate(Qref(nx,0:nt-1),Qmod(nx,0:nt-1))
!     !
!     write(*,*) ""
!     write(*,*) "----------------------- computing Qref(v,t) from colvar_file ------"
!     write(*,*) ""
!     !
!     Qref=0.d0
!     vmin= 1.d15
!     vmax=-1.d15
!     !
!     do i=1,nttot
!       it=nint(dble(nt)*colvar(i,1)/tmax)
!       if (it<nt) then
!         v=((colvar(i+1,2)-colvar(i,2))/dt) 
!         if (v>vmax) vmax=v
!         if (v<vmin) vmin=v
!       endif
!     enddo
!     dv=(vmax-vmin)/dble(nx-1)
!     !
!     do i=1,nttot
!       it=nint(dble(nt)*colvar(i,1)/tmax)
!       if (it<nt) then
!         v=((colvar(i+1,2)-colvar(i,2))/dt) 
!         iv=int((v-vmin)/dv)+1
!         Qref(iv,it)=Qref(iv,it)+1.d0
!       endif
!     enddo
!     !
!     write(*,*) "normalizing Qref by the number of initial points =",sum(Qref(:,0))
!     Qref=Qref/sum(Qref(:,0))
!     !
!     write(*,*) "writing Qref file"
!     open(qref_id,file="Qref",status="unknown")
!     ! NOTE THE INVERSION OF INDECES COMPARED TO Pref FILE
!     do it=0,nt-1
!       do iv=1,nx
!         write(qref_id,'(3F15.6)') vmin+(dble(iv)-1.d0)*dv,dble(it)*dt,Qref(iv,it)
!       enddo
!       write(qref_id,*) ""
!     enddo
!     close(qref_id)
!   !
!   end subroutine init_Qref
!================================================================================
subroutine init_mass
!
  use common_var
  !
  implicit none
  double precision :: v,vv,vv1,vv2
  integer :: i,it,nave,nave1,nave2
  !
  write(*,*) ""
  write(*,*) "----------------------- estimating mass -------------------------"
  write(*,*) ""
  vv=0.D0
  nave=0
  do i=1,nttot
    if (dabs(colvar(i,1)-0).lt.1.d-10) then
      nave=nave+1
      vv=vv+((colvar(i+1,2)-colvar(i,2))/dt)**2
    endif
  enddo
  !
  vv=vv/dble(nave)
  mass0=kT/vv
  write(*,*) ""
  write(*,*) "*** using velocities at the barrier top"
  write(*,*) "number of initial velocities =",nave
  write(*,*) "initial <v^2> =",vv
  write(*,*) "mass0 = kT/<v^2> =",mass0
  !
  write(*,*) ""
  write(*,*) "*** testing velocities at the minima bottom (last 30% of each traj.)"
  vv1=0.D0
  vv2=0.D0
  nave1=0
  nave2=0
  do i=1,nttot-1
    it=nint(dble(nt)*colvar(i,1)/tmax)
    if (colvar(i,1)>tmax*0.7d0.and.it<nt-1) then
    !test: if (colvar(i,1)>tmax*0.1d0.and.colvar(i,1)<tmax*0.2d0.and.it<nt) then
      v=((colvar(i+1,2)-colvar(i,2))/dt)
      if (colvar(i,2)<x0(1)) then
        vv1=vv1+v*v
        nave1=nave1+1
      else
        vv2=vv2+v*v
        nave2=nave2+1
      endif
    endif
  enddo
  write(*,*) "number of velocities =",nave1,nave2
  if (nave1>0) vv1=vv1/dble(nave1)
  if (nave2>0) vv2=vv2/dble(nave2)
  write(*,*) "<v^2> =",vv1,vv2
  write(*,*) "approx. masses = kT/<v^2> =",kT/vv1,kT/vv2 
  write(*,*) "WARNING: if trajectories are not equilibrated to the T given in input, these 2 masses are meaningless"
  !
  ! if (restart.eq.0) then
  !   write(*,*) "intialization: mass(x) = mass0"
  !   prof_m(:)   = mass0
  !   ! 
  !   ! TODO: init non-uniform mass profile using TS and minima values of mass ?
  ! else
    write(*,*) "mass read from RESTART discarding previous analysis"
    mass0=prof_m(int((x0(1)-xmin)/dxgrid)+1)
    write(*,*) "  mass0 (at initial point of first MD trajectory) =",mass0
  ! endif
  !
!
end subroutine init_mass
!================================================================================
subroutine init_potential
!
  use common_var
  !
  implicit none
  !double precision :: bottom1,bottom2,tmp
  !integer :: nave1,nave2,i,j,ibest,itype
  !double precision :: x,xx,dbest,control_p(2,5)
  !character :: line*80
  !
  !if (restart.ne.1) then
  !  !
  !  write(*,*) ""
  !  write(*,*) "----------------------- initializing potential ------------------"
  !  write(*,*) ""
  !  !
  !  write(*,*) "setting 5 points to control the potential: xmin, bottom1, x0, bottom2, xmax"
  !  !
  !  write(*,*) "estimated bottom1 and bottom2 from last 10% of trajectories:"
  !  nave1=0
  !  nave2=0
  !  bottom1=0.D0
  !  bottom2=0.D0
  !  do i=1,nttot
  !    if (colvar(i,1)>tmax*0.9D0) then
  !      if (colvar(i,2)<x0) then
  !        bottom1=bottom1+colvar(i,2)
  !        nave1=nave1+1
  !      else
  !        bottom2=bottom2+colvar(i,2)
  !        nave2=nave2+1
  !      endif    
  !    endif
  !  enddo
  !  if (nave1==0) call error("not enough points to the left  of x0 !")
  !  if (nave2==0) call error("not enough points to the right of x0 !")
  !  bottom1=bottom1/dble(nave1)
  !  bottom2=bottom2/dble(nave2)
  !  write(*,*) bottom1,bottom2
  !  !
  !  ! the following vertical values (kT units) are somehow arbitrary...
  !   control_p(1,1)=xmin          ! left border
  !   control_p(2,1)=0.5D0*init_barrier
  !  control_p(1,2)=bottom1       ! left minimum
  !  control_p(2,2)=0.D0
  !   control_p(1,3)=x0            ! saddle (origin of shootings)
  !   control_p(2,3)=init_barrier
  !  control_p(1,4)=bottom2       ! right minimum
  !  control_p(2,4)=0.D0
  !   control_p(1,5)=xmax          ! right border
  !   control_p(2,5)=0.5D0*init_barrier
  !  !
  !  call interpolate_potential_simplecos(control_p)
  !  !
  !endif ! restart
!
end subroutine init_potential
!================================================================================
subroutine init_friction
!
  use common_var
  !
  implicit none
  integer :: i,j,it=0
  double precision :: nave(2),ave(2),var(2),tau(2),Dtmp(2),told
  double precision, allocatable :: corr(:,:)
  !
  write(*,*) ""
  write(*,*) "----------------------- initializing friction ---------------------"
  write(*,*) ""
  !
  write(*,*) "analysis of the last 25% of trajectories (Hummer New J Phys 2005):"
  write(*,*) "(data below are for left and right basin)"
  !
  nave=0
  ave=0.d0
  var=0.d0
  do i=1,nttot
    if (colvar(i,1)>tmax*0.75d0) then
      if (colvar(i,2)<x0(1)) then
        j=1
      else
        j=2
      endif
      nave(j)=nave(j)+1  
      ave(j)=ave(j)+colvar(i,2)
      var(j)=var(j)+colvar(i,2)**2
    endif
  enddo
  write(*,*) "number of samples in basins left and right =",nave(1),nave(2)
  do j=1,2
    ave(j)=ave(j)/dble(nave(j))
    var(j)=var(j)/dble(nave(j))-ave(j)**2
  enddo
  write(*,*) "variance of the coordinate =",var(:)
  !
  allocate(corr(0:nt,2))
  corr=0.d0
  told=1.d9
  do i=1,nttot
    if (colvar(i,1)>tmax*0.75d0) then
      if (colvar(i,2)<x0(1)) then
        j=1
      else
        j=2
      endif
      if (colvar(i,1)<told) it=-1
      told=colvar(i,1)
      it=it+1
      if (it>nt) call error("problem in init_friction")
      corr(it,j)=corr(it,j)+(colvar(i-it,2)-ave(j))*(colvar(i,2)-ave(j))
    endif
  enddo
  do j=1,2
    corr(:,j)=corr(:,j)/(nave(j)*var(j))
    tau(j)=0.d0
    do i=0,it
      tau(j)=tau(j)+i*dt*corr(i,j)
    enddo
    tau(j)=tau(j)/sum(corr(:,j))
  enddo
  open(corrfunc_id,file=trim(colvar_file)//".corr_func",status="unknown")
  do i=0,it
    write(corrfunc_id,'(3F14.5)') i*dt,corr(i,1),corr(i,2)
  enddo
  close(corrfunc_id)
  write(*,*) "autocorrelation time       =",tau(:)
  do i=1,2
    if (tau(i)<0.) then
      write(*,*) "WARNING: negative tau, changing sign..."
      tau(i)=-tau(i)
    endif
  enddo
  write(*,*) "(correlation functions written in corr_func)"
  !
  Dtmp(1)=var(1)/tau(1)
  Dtmp(2)=var(2)/tau(2)
  write(*,*) "D=kT/m*gamma=var/tau =",Dtmp(:)
  write(*,*) "gamma                =",kT/(mass0*Dtmp(:))
  deallocate(corr)
  !
  ! if (restart.eq.0) then
  !   if (init_gamm.lt.0.D0) then
  !     gamm=0.5d0*sum(kT/(mass0*Dtmp(:)))
  !     gamm=dsqrt(dabs(gamm))
  !     write(*,*) "average gamma (sqrt) to be used =",gamm
  !   else
  !     gamm=init_gamm
  !     write(*,*) "discarding the previous analysis: user-provided gamma =",gamm
  !   endif
  !   !
  !   write(*,*) "intialization: gamma(x) = gamma"
  !   prof_g(:) = gamm
  ! else
    write(*,*) "gamma read from RESTART (discarding previous analysis)"
  ! endif
  !
  if (type_Langevin.ne.2) then
    taug=0.d0
    write(*,*) "ignoring tau (memory-less Langevin equation)"
  elseif (type_Langevin.eq.2) then
    if (init_taug.lt.0.D0) then
      ! TODO this choice is not correct...
      taug=0.5d0*(tau(1)+tau(2))
      write(*,*) "average tau to be used =",taug
    else
      taug=init_taug
      write(*,*) "discarding the previous analysis: user-provided tau =",taug
    endif
  endif
!
end subroutine init_friction
!================================================================================

