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
  do i=2,ngrid-1
    prof_force(i)=(prof_F(i-1)-prof_F(i+1))/dxgrid2
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
subroutine compute_error(error,type_err,iprintGauss)
!
  use common_var
  !
  implicit none
  integer :: ix,iv,it,i,j,ig,nave,type_err,iprintGauss, prop_order
  double precision :: error,eps,newnorm
  double precision :: dq,ddq,diff(ngrid),delta
  double precision :: a(ngrid),da(ngrid),dda(ngrid)       ! drift       and derivatives
  double precision :: D(ngrid),dD(ngrid),ddD(ngrid)       ! diff.coeff. and derivatives
  !double precision :: b(ngrid),db(ngrid),ddb(ngrid)       ! drift       and derivatives
  double precision :: b,db,ddb                            ! drift       and derivatives
  double precision :: phi(ngrid),dphi(ngrid),ddphi(ngrid) ! force/m     and derivatives
  double precision :: gam(ngrid),dgam(ngrid),ddgam(ngrid) ! gamma       and derivatives
  double precision :: q,v,v1,v2,q3,q2,q1,q0,f1,g1,sigma2,tau 
  double precision :: g,dg,ddt,ddt2,ddt3,tmp, qq,vv
  double precision :: Mq,Mqq,Mv,Mvv,Mqv,detM
  double precision :: Lqq,Lqv,Lvv
  double precision :: minuslogL,err_prop_L,err_prop_KS,dKS,pKS,p_traj
  integer, parameter :: ntraj=100 ! TODO: ask in input
  integer :: itraj
  double precision :: qlast(ntraj),qtmp,etmp,erftmp
  double precision, parameter :: pi=3.14159265359d0, pi2=pi*pi, sqrt2=sqrt(2.)
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
    ddt =dt*dtmult
    ddt2=ddt*ddt
    ddt3=ddt2*ddt
    !
    if (iprintGauss.eq.1) then
      ! here we write deviations from model: it should be a normalized Gaussian
      open(123,file="colvar_disp_scaled_from_prop",status="unknown") 
    endif
    !
    if (type_Langevin.eq.0) then
      !
      ! log likelihood from short-time overdamped propagator
      prop_order = 2   ! order of the propagator: TODO choose in the input
      !
      if (test_propagator) then
        open(60,file="err_prop",status="unknown")
        write(60,'(A)') "# propagator error: avg scaled Likelihood and Kolmogorov-Smirnov (dist, prob)"
        write(60,'(A)') "# note: for the first, the ideal value is 0.5*[log(2pi)+1] = 1.419"
        open(61,file="shooting_disp_scaled_from_prop",status="unknown")
        write(61,'(A)') "# error of the propagator: scaled displacements from Langevin shootings"
        write(61,'(A)') "# ideally, they should be distributed as a Gaussian N(0,1)"
        err_prop_L  = 0.D0
        err_prop_KS = 0.D0
        write(*,*) "TESTING THE PROPAGATOR AND EXITING ..."
      endif
      !
      D(:)=kT/(prof_m(:)*prof_g(:))
      do ig=2,ngrid-1
        dD(ig)=( D(ig+1)-D(ig-1) )/dxgrid2
      enddo
      dD(1)    =dD(2)
      dD(ngrid)=dD(ngrid-1)
      a(:)=( D(:)*prof_force(:)/kT + dD(:) )
      !
      if (prop_order.eq.1) then
        error=0.d0
        nave=0
        do i=1,n_tprop
          do j=ibeg_tprop(i),iend_tprop(i)-dtmult,dtmult
             q    = colvar(j,2)
            q2    = colvar(j+dtmult,2)
            dq    = q2-q
            ig    = int((q-xmin)/dxgrid)+1 ! pre-point
            Mq    = a(ig)*ddt ! we leave q out
            Mqq   = 2.d0*D(ig)*ddt
            minuslogL = 0.5d0*log(2.d0*pi*Mqq) + (dq-Mq)**2 / (2.d0*Mqq)
            error = error + minuslogL
            nave  = nave+1
            if (test_propagator) then 
               call prop_via_traj(q,q2,dt*dtmult,ntraj,qlast)
               ! compute likelihood of propagator given Langevin shootings on scaled data:
               qlast(:) = (qlast(:)-q-Mq)/sqrt(Mqq)
               etmp = 0.5d0*log(2.d0*pi) + &
                sum( qlast(1:ntraj)**2 ) / (ntraj*2.d0) 
               write(60,'(E14.6,$)') etmp
               err_prop_L = err_prop_L + etmp
               ! compare CDF of qlast with Gaussian CDF:
               call quicksort(qlast,1,ntraj)
               write(61,'(F6.3)') qlast(:)
               ! compute Kolmogorov-Smirnov statistic for normal distribution:
               call KSone_normal(qlast,ntraj,dKS,pKS)
               write(60,'(2F11.6)') dKS,pKS
               err_prop_KS = err_prop_KS + dKS
            endif
          enddo
        enddo
        error = error/dble(nave)
        if (test_propagator) then 
          write(60,'(A,2E14.6)') "# average errors (Likelihood, Kolmogorov-Smirnov) = ", &
           err_prop_L/dble(nave),err_prop_KS/dble(nave)
          close(60)
          close(61)
          stop
        endif
      endif
      !
      if (prop_order.eq.2) then
        do ig=2,ngrid-1
          ddD(ig)=( dD(ig+1)-dD(ig-1) )/dxgrid2
        enddo
        ddD(1)    =ddD(2)
        ddD(ngrid)=ddD(ngrid-1)
        do ig=2,ngrid-1
          da(ig)=( a(ig+1)-a(ig-1) )/dxgrid2
        enddo
        da(1)    =da(2)
        da(ngrid)=da(ngrid-1)
        do ig=2,ngrid-1
          dda(ig)=( da(ig+1)-da(ig-1) )/dxgrid2
        enddo
        dda(1)    =dda(2)
        dda(ngrid)=dda(ngrid-1)
        !
        error=0.d0
        nave=0
        do i=1,n_tprop
          do j=ibeg_tprop(i),iend_tprop(i)-dtmult,dtmult
             q    = colvar(j,2)
            q2    = colvar(j+dtmult,2)
            dq    = q2-q
            ig    = int((q-xmin)/dxgrid)+1 ! pre-point
            Mq    = a(ig)*ddt+0.5d0*( a(ig)*da(ig)+dda(ig)*D(ig) )*ddt2 ! we leave q out
            Mqq   = 2.d0*D(ig)*ddt+( a(ig)*dD(ig)+2.d0*da(ig)*D(ig)+D(ig)*ddD(ig) )*ddt2
            qq    = dq-Mq
            minuslogL = 0.5d0*log(2.d0*pi*Mqq) + qq*qq / (2.d0*Mqq)
            error = error + minuslogL
            nave  = nave+1
            if (test_propagator) then 
               call prop_via_traj(q,q2,dt*dtmult,ntraj,qlast)
               ! compute likelihood of propagator given Langevin shootings on scaled data:
               qlast(:) = (qlast(:)-q-Mq)/sqrt(Mqq)
               etmp = 0.5d0*log(2.d0*pi) + &
                sum( qlast(1:ntraj)**2 ) / (ntraj*2.d0) 
               write(60,'(E14.6,$)') etmp
               err_prop_L = err_prop_L + etmp
               ! compare CDF of qlast with Gaussian CDF:
               call quicksort(qlast,1,ntraj)
               write(61,'(F6.3)') qlast(:)
               ! compute Kolmogorov-Smirnov statistic for normal distribution:
               call KSone_normal(qlast,ntraj,dKS,pKS)
               write(60,'(2F11.6)') dKS,pKS
               err_prop_KS = err_prop_KS + dKS
            endif
            !
            if (iprintGauss.eq.1) then
              write(123,'(f9.5)') sqrt(1.d0/Mqq)*qq
            endif
          enddo
        enddo
        error = error/dble(nave)
        if (test_propagator) then 
          write(60,'(A,2E14.6)') "# average errors (Likelihood, Kolmogorov-Smirnov) = ", &
           err_prop_L/dble(nave),err_prop_KS/dble(nave)
          close(60)
          close(61)
          stop
        endif
      endif
      !
    elseif (type_Langevin.eq.1) then
      !
      ! log likelihood from short-time standard (Kramers) propagator, third order:
      ! note that b=phi-gamma*v, phi=force/m, D=gamma/m*beta, db/dv=-gamma
      !
      D(:)=prof_g(:)*kT/prof_m(:)
      do ig=2,ngrid-1
        dD(ig)=( D(ig+1)-D(ig+1) )/dxgrid2
      enddo
      dD(1)    =dD(2)
      dD(ngrid)=dD(ngrid-1)
      do ig=2,ngrid-1
        ddD(ig)=( dD(ig+1)-dD(ig+1) )/dxgrid2
      enddo
      ddD(1)    =ddD(2)
      ddD(ngrid)=ddD(ngrid-1)
      !
      phi(:)=prof_force(:)/prof_m(:)
      do ig=2,ngrid-1
        dphi(ig)=( phi(ig+1)-phi(ig-1) )/dxgrid2
      enddo
      dphi(1)    =dphi(2)
      dphi(ngrid)=dphi(ngrid-1)
      do ig=2,ngrid-1
        ddphi(ig)=( dphi(ig+1)-dphi(ig-1) )/dxgrid2
      enddo
      ddphi(1)    =ddphi(2)
      ddphi(ngrid)=ddphi(ngrid-1)
      !
      gam(:)=prof_g(:)
      do ig=2,ngrid-1
        dgam(ig)=( gam(ig+1)-gam(ig-1) )/dxgrid2
      enddo
      dgam(1)    =dgam(2)
      dgam(ngrid)=dgam(ngrid-1)
      do ig=2,ngrid-1
        ddgam(ig)=( dgam(ig+1)-dgam(ig-1) )/dxgrid2
      enddo
      ddgam(1)    =ddgam(2)
      ddgam(ngrid)=ddgam(ngrid-1)
      !
      error=0.d0
      nave=0
      do i=1,n_tprop
        do j=ibeg_tprop(i),iend_tprop(i)-dtmult,dtmult
          q      = colvar(j,2)
          q2     = colvar(j+dtmult,2)
          v      = colvar(j,3)
          v2     = colvar(j+dtmult,3)
          ig     = int((q-xmin)/dxgrid)+1
          b      =   phi(ig)  -gam(ig)*v
          db     =  dphi(ig) -dgam(ig)*v
          ddb    = ddphi(ig)-ddgam(ig)*v
          g      = gam(ig)
          dg     = dgam(ig)
          !
          tmp = (db*v-b*g)
          Mq  = q + v*ddt + b*ddt2/2.d0 + tmp*ddt3/6.d0
          Mqq = 2.d0*D(ig)*ddt3/3.d0
          Mv  = v + b*ddt + tmp*ddt2/2.d0 &
              + ((ddb*v-db*g-2.d0*b*dg)*v+b*db+b*g*g-2.d0*D(ig)*dg)*ddt3/6.d0
          Mvv = 2.d0*D(ig)*ddt + (dD(ig)*v-2.d0*D(ig)*g)*ddt2 &
              + ((ddD(ig)*v-2.d0*dD(ig)*g-4.d0*D(ig)*dg)*v+dD(ig)*b+2.d0*D(ig)*(db+2.d0*g*g))*ddt3/3.d0
          Mqv = D(ig)*ddt2 + (dD(ig)*v-3.d0*D(ig)*g)*ddt3/3.d0
          detM = Mqq*Mvv-Mqv*Mqv 
          !
          qq = (q2-Mq)
          vv = (v2-Mv)
          minuslogL = 0.5d0*log(4.d0*pi2*detM) &
                + ( 0.5d0*Mvv*qq*qq + 0.5d0*Mqq*vv*vv - Mqv*qq*vv )/detM
          error = error + minuslogL
          nave  = nave+1
          !
          if (iprintGauss.eq.1) then
            ! we need Cholesky decomposition to get vector 
            !  of Gaussian numbers with zero mean and unit variance:
            Lqq = dsqrt(Mvv)
            Lqv = -Mqv/Lqq
            Lvv = dsqrt(Mqq-Mqv*Mqv/Mvv)
            write(123,'(2f9.5)') (Lqq*qq+Lqv*vv)/dsqrt(detM), Lqv*vv/dsqrt(detM)
          endif
        enddo
      enddo
      error = error/dble(nave)
      !
    endif
    !
    if (iprintGauss.eq.1) then
      close(123)
    endif
    !
  endif
  ! TODO: add also Kolmogorov-Smirnov statistic on P(x,t) (discarding velocities...)
!
end subroutine compute_error
!======================================================================
subroutine prop_via_traj(q1,q2,dt_data,ntraj,qlast)
!
  use common_var
  !
  double precision :: q1,q2,dt_data,p_traj,qlast(ntraj)
  double precision :: dtint
  integer :: ntraj ! important parameter: number of shootings
  integer :: i,it,igrid1,igrid2,nstep
  double precision, allocatable :: G_vec(:)
  double precision :: q,G,noisefac,mynoise,tmp
  double precision :: D_over_kT(ngrid),sqrt_2_D_over_dt(ngrid),dD(ngrid)
  double precision :: kde_sigma
  double precision, parameter :: pi=3.14159265359d0
  
  !---------------------------------------------------------------------------------- 
  ! integrate overdamped Langevin traj. to estimate propagator q1,0 -> q2,dt_data
  ! 
  ! D = kT / mg  (gamma is an inverse time)
  ! we use Milstein integrator (better for pos-dep D, see Kloeden-Platen 10.3):
  !   q' = q + ( (D/kT)*F + 0.5*dD/dq )*dt + sqrt(2*D*dt)*G + 0.5*dD/dq*dt*G^2
  ! instead of Euler-Maruyama:
  !   q' = q + ( (D/kT)*F + dD/dq )*dt + sqrt(2*D*dt)*G
  ! with G a Gaussian random number with <G>=0, <G^2>=1
  !----------------------------------------------------------------------------------
  dtint=dt/100. ! TODO arbitrary: rather, give in input and check it is good on-the-fly!
  nstep = dt_data/dtint
  D_over_kT(:)   = 1.d0/(prof_m(:)*prof_g(:))
  sqrt_2_D_over_dt(:) = sqrt(2.d0*kT*D_over_kT(:)/dtint)
  do igrid1=2,ngrid-1
    dD(igrid1)=kT*(D_over_kT(igrid1+1)-D_over_kT(igrid1-1))/dxgrid2
  enddo
  dD(1)     = dD(2)
  dD(ngrid) = dD(ngrid-1)
  allocate(G_vec(nstep+1))
  !
  do it = 1,ntraj ! loop on ntraj
    call noise_vector(G_vec,nstep) 
    q = q1
    !
    do i = 1,nstep ! main loop
      ig = int((q-xmin)/dxgrid)+1  
      if (ig<1) then ! the test is expensive, but the shootings can go beyond borders!
        ig=1             
      else
        if (ig>ngrid) ig=ngrid 
      endif
      G = G_vec(i)
      tmp = 0.5d0*dD(ig)
      mynoise  = ( sqrt_2_D_over_dt(ig) + tmp*G )*G                  ! Milstein 
      q = q + ( D_over_kT(ig)*prof_force(ig) + tmp + mynoise )*dtint ! Milstein
    enddo ! main loop
    qlast(it) = q
  enddo ! loop on ntraj
!!!  !
!!!  ! kernel density estimation of the transition probability q1->q2
!!!  ! I use Silverman's rule to guess a good sigma of the Gaussian kernels
!!!  ! (see wikipedia, https://archive.org/details/densityestimatio00silv_0/page/45)
!!!  ! which is suitable for prob in the form of a Gaussian (maybe a strong approx?)
!!!  !
!!!  kde_sigma = 1.06d0*dsqrt( sum(qlast(:)**2)/ntraj - (sum(qlast(:))/ntraj)**2 )
!!!  kde_sigma = kde_sigma / (ntraj**(0.2d0))
!!!  p_traj = 0.d0
!!!  do it = 1,ntraj
!!!    p_traj = p_traj + dexp( -0.5d0*( (q2-qlast(it))/kde_sigma )**2 )
!!!  enddo
!!!  p_traj = p_traj / ( ntraj*dsqrt(2.d0*pi)*kde_sigma )
!
end subroutine prop_via_traj
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
  call compute_error(minerr1,1,0)  
  call compute_error(minerr2,2,0)  
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
  !do i=1,10000
   ! call noise(G)
    !write(222,'(F12.6)') G
  !enddo
  !
end subroutine init_random_seed
!======================================================================
subroutine noise_vector(G_vec,n)
! this generates a vector of Gaussian-distributed random numbers
!
implicit none
integer :: i,j,n
double precision :: G_vec(n+1) ! to take care of possible odd n
integer :: iset
double precision :: r,fac,gset,rsq,v1,v2
!
do i=1,(n+1)/2 ! to take care of possible odd n
1 call random_number(r)
  v1=2.*r-1.
  call random_number(r)
  v2=2.*r-1.
  rsq=v1**2+v2**2
  if (rsq.ge.1..or.rsq.eq.0.) goto 1
  fac=sqrt(-2.*log(rsq)/rsq)
  j=i*2
  G_vec(j-1)=v1*fac
  G_vec(j)  =v2*fac
enddo
!
end subroutine noise_vector
!================================================================================
! quicksort.f -*-f90-*-
! Author: t-nissie
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!
recursive subroutine quicksort(a, first, last)
  implicit none
  real*8  a(*), x, t
  integer first, last
  integer i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, first, i-1)
  if (j+1 < last)  call quicksort(a, j+1, last)
end subroutine quicksort
!================================================================================
