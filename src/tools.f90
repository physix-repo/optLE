!================================================================================
subroutine error(line)
!
  implicit none
  character :: line*(*)
  write(*,*) "######## ",line
  write(*,*) "######## EXITING"
  stop
!
end subroutine error
!================================================================================
subroutine update_prof_force
!
  use common_var
  !
  implicit none
  integer :: i
  !
  do i=1,ngrid-1
    prof_force(i)=(prof_F(i)-prof_F(i+1))/dxgrid
  enddo
  prof_force(ngrid)=prof_force(ngrid-1)
!
end subroutine update_prof_force
!================================================================================
subroutine print_profiles
!
  use common_var
  !
  implicit none
  integer :: i
  double precision :: x
  !
  write(*,*) ""
  write(*,*) "writing F,gamma,mass in PROFILES"
  write(*,*) ""
  open(77,file="PROFILES",status="unknown")
  write(77,'(A,E13.4)') "# x F F/kT gamma mass; taug=",taug
  do i=1,ngrid
    x=xmin+dble(i-1)*dxgrid
    write(77,'(5E15.6)') x,prof_F(i)-minval(prof_F(:)),(prof_F(i)-minval(prof_F(:)))/kT,&
     prof_g(i),prof_m(i)
  enddo
  close(77)
!
end subroutine print_profiles
!================================================================================
subroutine compute_error(error,type_err)
!
  use common_var
  !
  implicit none
  integer :: ix,iv,it,i,igrid,nave,type_err, prop_order
  double precision :: error,eps,newnorm
  double precision :: dq,ddq,diff(ngrid),delta, M1,M11
  double precision :: a(ngrid),da(ngrid),dda(ngrid) ! drift       and derivatives
  double precision :: D(ngrid),dD(ngrid),ddD(ngrid) ! diff.coeff. and derivatives
  double precision :: q2,q1,q0,f1,g1,sigma2,tau
  double precision, parameter :: pi=3.14159265359d0
  !
  if (type_err.eq.1) then
    !
    ! Kullback-Leibler divergence: err=sum(-Pref*log(Pmod/Pref))
    ! log likelihood: L = sum(Pref*log(Pmod))
    ! I use the latter, simpler, and shift a little bit to include empty regions ("smoothing"?)
    eps=(0.001D0)/dble(nx*nx*(nt+1))  ! altering total norm by 0.001 (arbitrary, but small)
    newnorm=1.D0+0.001D0
    ! TODO improve using Gaussian kernels...
    error = -sum( Pref(:,:,:)*log((Pmod(:,:,:)+eps)/newnorm) )
  elseif (type_err.eq.2) then
    !
    ! root mean square deviation
    if (use_velocity.eq.1) then
      error = dsqrt( sum((Pref(:,:,:)-Pmod(:,:,:))**2)/dble(nx*nx*(nt+1)) )
    else
      error = dsqrt( sum((Pref(:,1,:)-Pmod(:,1,:))**2)/dble(nx*(nt+1)) )
    endif
    !
  elseif (type_err.eq.3) then
    !
    if (type_Langevin.eq.0) then
      !
      ! log likelihood from short-time overdamped propagator
      prop_order = 2   ! order of the propagator: TODO choose in the input
      !
      if (prop_order.eq.1) then
        D(:)=kT/(prof_m(:)*prof_g(:))
        do igrid=1,ngrid-1
          dD(igrid)=( D(igrid+1)-D(igrid) )/dxgrid
        enddo
        dD(ngrid)=dD(ngrid-1)
        !diff(:)=4.d0*D(:)*dt
        a(:)=( D(:)*prof_force(:)/kT + dD(:) )
        error=0.d0
        nave=0
        do i=2,nttot
          if (colvar(i,1)-colvar(i-1,1)>0.d0) then ! (avoid time discontinuity)
            dq    = colvar(i,2)-colvar(i-1,2)
            igrid = int((colvar(i-1,2)-xmin)/dxgrid)+1
            M1    = a(igrid)*dt ! we leave q out
            M11   = 2.d0*D(igrid)*dt
            error = error + 0.5d0*log(2.d0*pi*M11) + (dq-M1)**2 / (2.d0*M11)
            !error = error + 0.5d0*log(diff(igrid)*pi) + (dq-a(igrid)*dt)**2/(diff(igrid))
            nave  = nave+1
          endif 
        enddo
        error = error/dble(nave)
      endif
      !
      if (prop_order.eq.2) then
        D(:)=kT/(prof_m(:)*prof_g(:))
        do igrid=1,ngrid-1
          dD(igrid)=( D(igrid+1)-D(igrid) )/dxgrid
        enddo
        dD(ngrid)=dD(ngrid-1)
        do igrid=1,ngrid-1
          ddD(igrid)=( dD(igrid+1)-dD(igrid) )/dxgrid
        enddo
        ddD(ngrid)=ddD(ngrid-1)
        a(:)=( D(:)*prof_force(:)/kT + dD(:) )
        do igrid=1,ngrid-1
          da(igrid)=( a(igrid+1)-a(igrid) )/dxgrid
        enddo
        da(ngrid)=da(ngrid-1)
        do igrid=1,ngrid-1
          dda(igrid)=( da(igrid+1)-da(igrid) )/dxgrid
        enddo
        dda(ngrid)=dda(ngrid-1)
        error=0.d0
        nave=0
        do i=2,nttot
          if (colvar(i,1)-colvar(i-1,1)>0.d0) then ! (avoid time discontinuity)
            dq    = colvar(i,2)-colvar(i-1,2)
            igrid = int((colvar(i-1,2)-xmin)/dxgrid)+1 ! pre-point
            M1    = a(igrid)*dt+0.5d0*( a(igrid)*da(igrid)+dda(igrid)*D(igrid) )*dt*dt ! we leave q out
            M11   = 2.d0*D(igrid)*dt+( a(igrid)*dD(igrid)+2.d0*da(igrid)*D(igrid)+D(igrid)*ddD(igrid) )*dt*dt
            error = error + 0.5d0*log(2.d0*pi*M11) + (dq-M1)**2 / (2.d0*M11)
            nave  = nave+1
          endif 
        enddo
        error = error/dble(nave)
      endif
      !
    elseif (type_Langevin.eq.1) then
      !
      ! log likelihood from short-time standard (Kramers) propagator
      ! TODO: this is is a bad surrogate, by inverting the integrator, replace with Drozdov!
      error=0.d0
      nave=0
      tau=dt/dble(dtmult)
      !
      do i=3,nt_interp
        if (interp_colvar(i,1)-interp_colvar(i-2,1)>0.d0) then ! (avoid time discontinuity)
          q2     = interp_colvar(i,2)
          q1     = interp_colvar(i-1,2)
          q0     = interp_colvar(i-2,2)
          ddq    = ( q2-2.d0*q1+q0 )/tau
          dq     = q1-q0 
          igrid  = int((q1-xmin)/dxgrid)+1
          f1     = prof_force(igrid)
          m1     = prof_m(igrid)
          g1     = prof_g(igrid)
          delta  = ( ddq - tau*f1/m1 + g1*dq )
          sigma2 = tau*2.d0*kT*g1/m1
          delta  = delta*delta / sigma2
          error  = error + 0.5d0*log(2.d0*pi*sigma2) + 0.5d0*delta ! missing the normalization here
          nave   = nave+1
        endif 
      enddo
      error = error/dble(nave)
      !
    endif
    !
  endif
  ! TODO: add also Kolmogorov-Smirnov statistic on P(x,t) (discarding velocities...)
!
end subroutine compute_error
!======================================================================
subroutine compute_typical_min_error(minerr1,minerr2)
!
  use common_var
  !
  double precision :: minerr1,minerr2
  double precision :: Ptmp(nx,nx,0:nt)
  !
  ! store the correct Pref (from MD data)
  Ptmp=Pref
  ! compute 1st Pmod (a particular realization depending on random sequence)
  call compute_Pmod ! run ntraj_Langevin simulations and compute Pmod
  ! pretend Pref is Pmod
  Pref=Pmod
  ! compute 2nd Pmod (a particular realization depending on random sequence)
  call compute_Pmod ! run ntraj_Langevin simulations and compute Pmod
  ! compute errors between the two Pmod obtained from identical parameters
  call compute_error(minerr1,1)  
  call compute_error(minerr2,2)  
  write(*,'(A,2E11.4)') " typical min error (between two realizations of the same model) = ",&
   minerr1,minerr2
  ! store back the original Pref
  Pref=Ptmp
!
end subroutine compute_typical_min_error
!======================================================================
function lcg(s)
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  use iso_fortran_env, only: int64
  integer :: lcg
  integer(int64) :: s
  if (s == 0) then
     s = 104729
  else
     s = mod(s, 4294967296_int64)
  end if
  s = mod(s * 279470273_int64, 4294967291_int64)
  lcg = int(mod(s, int(huge(0), int64)), kind(0))
end function lcg
!======================================================================
subroutine init_random_seed()
  use iso_fortran_env, only: int64
  implicit none
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid, lcg
  integer(int64) :: t
  double precision :: G

  ! this code is taken from https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
  write(*,*) ""
  write(*,*) "initializing random number generator" 
  call random_seed(size = n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
  open(newunit=un, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
     read(un) seed
     close(un)
  else
     ! Fallback to XOR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.
     call system_clock(t)
     if (t == 0) then
        call date_and_time(values=dt)
        t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24_int64 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
     end if
     pid = getpid()
     t = ieor(t, int(pid, kind(t)))
     do i = 1, n
        seed(i) = lcg(t)
     end do
  end if
  call random_seed(put=seed)
  !
  !!!!!!!!!!!!! test Gaussian random numbers:
  do i=1,10000
    call noise(G)
    write(222,'(F12.6)') G
  enddo
  !
end subroutine init_random_seed
!======================================================================
