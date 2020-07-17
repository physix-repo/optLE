subroutine Langevin_traj_overdamped
!
  use common_var
  !
  implicit none
  integer :: ix,iv,it,i,igrid1,igrid2,nstep
  double precision :: dtint,x,t,v,G,force,mforce,xold,xnew,vnew, mass
  double precision :: noisefac,mynoise, tmp 
  double precision :: D_over_kT(ngrid),sqrt_2_D_over_dt(ngrid),dD_over_dx(ngrid)
  !---------------------------------------------------------------------------------- 
  ! D = kT / mg  (gamma is an inverse time)
  ! As integrator we can use Euler-Maruyama:
  !   q' = q + ( (D/kT)*F + dD/dq )*dt + sqrt(2*D*dt)*G
  ! with G a Gaussian random number with <G>=0, <G^2>=1
  ! or we can use Milstein (better for pos-dep D, see Kloeden-Platen 10.3):
  !   q' = q + ( (D/kT)*F + 0.5*dD/dq )*dt + sqrt(2*D*dt)*G + 0.5*dD/dq*dt*G^2
  !---------------------------------------------------------------------------------- 
  !TODO: store noise sequence (to avoid calling many times "noise")
  dtint    = dt/dble(dtmult) ! TODO: do automatic test to see if it is sufficiently accurate...
  nstep    = nint(tmax/dtint) 
  t        = 0.D0
  x        = x0now     ! initial position: saddle point 
  v        = 0.D0      ! useless in overdamped eq
  ! igrid1 = int((x0now-xmin)/dxgrid)+1
  D_over_kT(:)   = 1.d0/(prof_m(:)*prof_g(:))
  sqrt_2_D_over_dt(:) = sqrt(2.d0*kT*D_over_kT(:)/dtint)
  do igrid1=1,ngrid-1
    dD_over_dx(igrid1)=kT*(D_over_kT(igrid1+1)-D_over_kT(igrid1))/dxgrid
  enddo
  dD_over_dx(ngrid)=dD_over_dx(ngrid-1)
  !
#ifdef DEBUG
  write(111,'(A)') "# t,x,v  (with mforce)"
  write(111,'(E16.8,2E14.5)') t,x,v ! DEBUG traj
#endif
  !
  xold=x
  !!!!!!!!!!!!!!!! run one trajectory
  do i = 1,nstep ! main loop
    !
    t = t+dtint
    !
    igrid1=int((x-xmin)/dxgrid)+1
    !igrid2=min(igrid1+1,ngrid)
    !if (igrid1>=1.and.igrid2<=ngrid) then ! TODO avoid this to speed-up
      force  = prof_force(igrid1) ! dim = m*s/t**2
      !!!if (mod(i,dtmult).eq.0) write(112,'(4E14.5)') force,mforce,mass,gamm !DDDDD
    !endif
    ! 
    call noise(G) ! TODO: store in memory vector of G to accelerate?
    ! mynoise  = sqrt_2_D_over_dt(igrid1)*G ! Euler-Maruyama
    mynoise  = ( sqrt_2_D_over_dt(igrid1) + 0.5d0*dD_over_dx(igrid1)*G )*G ! Milstein
    ! 
    xnew = x + (D_over_kT(igrid1)*force+0.5d0*dD_over_dx(igrid1)+mynoise)*dtint 
    x = xnew
    ! TODO: here is arbitrary, to avoid the cost of "if"
    x=max(x,xmin)
    x=min(x,xmax)
!    v=max(v,vmin) ! careful: you need wide enough boundaries...
!    v=min(v,vmax)
    !
#ifdef DEBUG
    if (mod(i,dtmult).eq.0) write(111,'(E16.8,E14.5)') t,x ! DEBUG traj
#endif
    !
    !it=nint(t/dt) ! important to use nint!
    if (mod(i,dtmult)==0) then 
!AAA      ix=int((x-xmin)/dx)+1
      ix=int((xold-xmin)/dx)+1
      !iv=int(((x-xold)/dt-vmin)/dv)+1
      !iv=min(iv,nx)
      !iv=max(iv,1)
      xold=x
    ! if ((ix>=1.and.ix<=nx).and.(it>=0.and.it<=nt)) then ! TODO avoid this to speed-up
    !  if (iv>=1.and.iv<=nx) then ! TODO avoid this to speed-up
!AAA        it=nint(t/dt) 
        it=nint(t/dt)-1 
        Pmod(ix,1,it) = Pmod(ix,1,it) + 1.D0 ! no velocities here
    !  endif
      !debug write(*,*) ix,it,t,x
    !endif
    endif
    ! 
  enddo ! main loop
!
end subroutine Langevin_traj_overdamped
!================================================================================
subroutine Langevin_traj_std ! TODO TODO TODO replace with Haynes JCP 101 7811 1994
!
  use common_var
  !
  implicit none
  integer :: ix,iv,it,i,igrid1,igrid2,nstep
  double precision :: dtint,x,t,v,G,force,mforce,xold,xnew,vnew, mass
  double precision :: noisefac,mynoise, tmp
  !---------------------------------------------------------------------------------- 
  ! We follow Risken: m*a = F - m*g*v + R, <R(0)R(t)> = 2*m*g*k*T*delta(t), D = k*T/m*g
  ! in this way gamma is an inverse time.
  ! As integrator we use Euler-Maruyama:
  !   x' = x + v*dt
  !   v' = v + (F/m)*dt - g*v*dt + sqrt(dt*2*k*T*g/m)*G
  ! with G a Gaussian random number with <G>=0, <G^2>=1
  ! (e.g., see Vanden-Eijnden & Ciccotti, Chem Phys Lett 429, 310 (2006)) 
  !---------------------------------------------------------------------------------- 
  igrid1 = int((x0now-xmin)/dxgrid)+1
  gamm   = prof_g(igrid1)
  mass   = prof_m(igrid1)
  !
  dtint    = dt/dble(dtmult) ! TODO: do automatic test to see if it is sufficiently accurate...
  nstep    = nint(tmax/dtint) 
  t        = 0.D0
  x        = x0now     ! initial position: saddle point 
  call noise(G)
  v        = dsqrt(kT/mass)*G ! initial velocity: Boltzmann distribution
  tmp      = 2.D0*kT*dtint
  !
#ifdef DEBUG
  write(111,'(A)') "# t,x,v  (with mforce)"
  write(111,'(E16.8,2E14.5)') t,x,v ! DEBUG traj
#endif
  !traj(0,1)=x
  !traj(0,2)=v
  !traj(0,3:5)=0.d0
  !
  xold=x
!AAA  ix=int((x-xmin)/dx)+1
!AAA  iv=int((v-vmin)/dv)+1
!AAA  it=nint(t/dt) ! important to use nint!
!AAA  if ((ix>=1.and.ix<=nx).and.(iv>=1.and.iv<=nx).and.(it>=0.and.it<=nt)) then
!AAA    Pmod(ix,iv,it) = Pmod(ix,iv,it) + 1.D0
!AAA    !debug write(*,*) ix,it,t,x
!AAA  endif
  !!!!!!!!!!!!!!!! run one trajectory
  do i = 1,nstep ! main loop
    !
    t = t+dtint
    !
    igrid1=int((x-xmin)/dxgrid)+1
    igrid2=min(igrid1+1,ngrid)
    !if (igrid1>=1.and.igrid2<=ngrid) then ! TODO avoid this to speed-up
      force  = prof_force(igrid1) ! dim = m*s/t**2
      mforce = (prof_m(igrid1)-prof_m(igrid2))/dxgrid ! dim = m*s/t**2
      mforce = 0.5d0*mforce*v*v
      gamm   = prof_g(igrid1)
      mass   = prof_m(igrid1)
      noisefac = dsqrt(tmp*gamm/mass) ! dim = s/t
      !!!if (mod(i,dtmult).eq.0) write(112,'(4E14.5)') force,mforce,mass,gamm !DDDDD
    !endif
    ! 
    call noise(G) ! TODO: store in memory vector of G to accelerate?
    mynoise = noisefac*G ! dim = s/t
    ! 
    xnew = x + dtint*v
    !vnew = v + dtint*((force+mforce)/mass -gamm*v) + mynoise ! TODO: put back mforce
    vnew = v + dtint*((force)/mass -gamm*v) + mynoise
    !write(114,*) (vnew-v-dtint*((force)/mass -gamm*v))/noisefac ! DEBUG
    x = xnew
    v = vnew
    ! TODO: here is arbitrary, to avoid the cost of "if"
    x=max(x,xmin)
    x=min(x,xmax)
!    v=max(v,vmin) ! careful: you need wide enough boundaries...
!    v=min(v,vmax)
    !
#ifdef DEBUG
    if (mod(i,dtmult).eq.0) write(111,'(E16.8,2E14.5)') t,x,v ! DEBUG traj
#endif
    !
    !it=nint(t/dt) ! important to use nint!
    if (mod(i,dtmult)==0) then  ! avoid using double loop on time
!AAA      ix=int((x-xmin)/dx)+1
      ix=int((xold-xmin)/dx)+1
      iv=int(((x-xold)/dt-vmin)/dv)+1
      iv=min(iv,nx)
      iv=max(iv,1)
      xold=x
    ! if ((ix>=1.and.ix<=nx).and.(it>=0.and.it<=nt)) then ! TODO avoid this to speed-up
    !  if (iv>=1.and.iv<=nx) then ! TODO avoid this to speed-up
!AAA        it=nint(t/dt) 
        it=nint(t/dt)-1 
        Pmod(ix,iv,it) = Pmod(ix,iv,it) + 1.D0
    !  endif
      !debug write(*,*) ix,it,t,x
    !endif
    endif
    ! 
  enddo ! main loop
!
end subroutine Langevin_traj_std
!================================================================================
subroutine Langevin_traj_GLEexp ! TODO TODO TODO replace with Haynes JCP 101 7811 1994
!---------------------------------------------------------------------------------- 
! We adopt the following formulation:
!  m*a = -dF/dx -int_0^t dt' g(x(t))*eta(t-t')*g(x(t'))*p(t') + g(x(t))*R
!  <R(t)R(t')>_x0 = m*kT*eta(t-t')  (x0 = TS, g(x0)=1)
!  eta(t-t') = eta0*exp(-(t-t')/tau) 
! We use the integrator in Haynes JCP 101 7811 1994:
!  dx = (p/m)*dt
!  dp = dt*(-dF/dx - 0.5*dm/dx*v**2 + z*g)
!  dz = dt*(-eta0*p*g - z/tau) + sqrt((eta0/tau)*m*kT*dt)*G
!  TODO : how to initialize z ???
!---------------------------------------------------------------------------------- 
!
  use common_var
  !
  implicit none
  integer :: ix,iv,it,i,igrid1,igrid2,nstep
  double precision :: dtint,x,t,v,p,z,f,R,xold,xnew,pnew,znew,vnew,fnew,Rnew, mass
  double precision :: eta0,G,force,mforce,facnoise
  !double precision :: traj(nstep,5)
  igrid1    = int((x0now-xmin)/dxgrid)+1
  eta0      = prof_g(igrid1)**2      ! during integration we use this eta0 ...
  prof_g(:) = prof_g(:)/dsqrt(eta0)  ! ... and we normalize gamma = 1 at shooting point (at the end we reset)
  !write(777,'(e12.6)') prof_g(:) ! XXXXXX
  facnoise  = dsqrt((eta0/taug)*kT*dtint) 
  mass      = prof_m(igrid1)
  !
  dtint    = dt/dble(dtmult) ! TODO: do automatic test to see if it is sufficiently accurate...
  nstep    = nint(tmax/dtint)
  !!! init !!!!
  t        = 0.D0
  x        = x0now      ! initial position: saddle point 
  call noise(G)
  v        = dsqrt(kT/mass)*G ! initial velocity: Boltzmann distribution TODO check wrt COLVAR !!!
  p        = mass*v
  z        = 0.D0 ! TODO how to correctly initialize z ???
  !!!!!!!!!!!!!!
#ifdef DEBUG
  !write(111,'(A)') "# t,x,v,force,friction,noise" ! DEBUG traj
  write(111,'(A)') "# t,x,v,p,z (with mforce)"
  write(111,'(E16.8,5E14.5)') t,x,v,p,z ! DEBUG traj
#endif
  !
  xold=x
!AAA  ix=int((x-xmin)/dx)+1
!AAA  iv=int((v-vmin)/dv)+1
!AAA  it=nint(t/dt) ! important to use nint!
!AAA  if ((ix>=1.and.ix<=nx).and.(iv>=1.and.iv<=nx).and.(it>=0.and.it<=nt)) then
!AAA    Pmod(ix,iv,it) = Pmod(ix,iv,it) + 1.D0
!AAA    !debug write(*,*) ix,it,t,x
!AAA  endif
  !!!!!!!!!!!!!!!! run one trajectory
  do i = 1,nstep ! main loop
    !
    t = t+dtint
    !
    igrid1=int((x-xmin)/dxgrid)+1
    !igrid2=min(igrid1+1,ngrid)
    !if (igrid1>=1.and.igrid2<=ngrid) then ! we avoid this to speed-up
      force = prof_force(igrid1)
      mass  = prof_m(igrid1)
      v     = p/mass
      ! mforce = 0.5*(mass-prof_m(igrid2))*v*v/dxgrid   ! (with -) we still need /dxgrid  
      mforce = 0.d0 ! TODO: check if we need to put back...
      gamm  = prof_g(igrid1)
    !endif
    ! 
    call noise(G)
    ! 
    xnew = x + dtint*v
    pnew = p + dtint*((force + mforce) + z*gamm)
    znew = z - dtint*(eta0*p*gamm + z/taug) + facnoise*dsqrt(mass)*G
    ! TODO you can avoid computing dsqrt by storing also vector sqrtm(:) = dsqrt(prof_m(:)) ...
    x = xnew
    p = pnew
    z = znew
    ! TODO: here is arbitrary, to avoid the cost of "if"
    x=max(x,xmin)
    x=min(x,xmax)
!    v=max(v,vmin) ! careful: you need wide enough boundaries...
!    v=min(v,vmax)
    !
#ifdef DEBUG
    if (mod(i,dtmult).eq.0) write(111,'(E16.8,5E14.5)') t,x,p/mass,p,z ! DEBUG traj
#endif
    !
    if (mod(i,dtmult)==0) then 
!AAA      ix=int((x-xmin)/dx)+1
      ix=int((xold-xmin)/dx)+1
      iv=int(((x-xold)/dt-vmin)/dv)+1
      iv=min(iv,nx)
      iv=max(iv,1)
      xold=x
    ! if ((ix>=1.and.ix<=nx).and.(it>=0.and.it<=nt)) then ! TODO avoid this to speed-up
    !  if (iv>=1.and.iv<=nx) then ! TODO avoid this to speed-up
!AAA        it=nint(t/dt) 
        it=nint(t/dt)-1 
        Pmod(ix,iv,it) = Pmod(ix,iv,it) + 1.D0
    !  endif
      !debug write(*,*) ix,it,t,x
    !endif
    endif
    ! 
  enddo ! main loop
  prof_g(:) = prof_g(:)*dsqrt(eta0)  ! here we reset gamma to the definition including eta0
!
end subroutine Langevin_traj_GLEexp
!======================================================================
subroutine noise(dw)
! this generates random number according to a Gaussian distribution
!
implicit none
double precision :: dw
integer :: iset
double precision :: r,fac,gset,rsq,v1,v2
save :: iset,gset
data iset/0/
!
if (iset.eq.0) then
1 call random_number(r)
  v1=2.*r-1.
  call random_number(r)
  v2=2.*r-1.
  rsq=v1**2+v2**2
  if(rsq.ge.1..or.rsq.eq.0.)goto 1
  fac=sqrt(-2.*log(rsq)/rsq)
  gset=v1*fac
  dw=v2*fac
  iset=1
else
  dw=gset
  iset=0
endif
!
end subroutine noise
!================================================================================
!!!!!!!! OLD STUFF !!!!!!!!!!
!================================================================================
!  subroutine Langevin_traj_std_fixedgamma
!  !
!    use common_var
!    !
!    implicit none
!    integer :: ix,iv,it,i,igrid1,igrid2,nstep
!    double precision :: dtint,x,t,v,G,force,xnew,vnew
!    double precision :: noisefac,mynoise
!    !---------------------------------------------------------------------------------- 
!    ! We follow Risken: m*a = F - m*g*v + R, <R(0)R(t)> = 2*m*g*k*T*delta(t), D = k*T/m*g
!    ! in this way gamma is an inverse time.
!    ! As integrator we use Euler-Maruyama:
!    !   x' = x + v*dt
!    !   v' = v + F*dt - g*v*dt + sqrt(dt*2*k*T*g/m)*G
!    ! with G a Gaussian random number with <G>=0, <G^2>=1
!    ! (e.g., see Vanden-Eijnden & Ciccotti, Chem Phys Lett 429, 310 (2006)) 
!    !---------------------------------------------------------------------------------- 
!    dtint    = dt ! TODO: how to choose? 
!    noisefac = dsqrt(2.D0*kT*gamm*dtint/mass) ! dim = s/t
!    nstep    = nt ! nint(tmax/dtint) !
!    t        = 0.D0
!    x        = x0        ! initial position: saddle point 
!    call noise(G)
!    v        = dsqrt(kT/mass)*G ! initial velocity: Boltzmann distribution
!    !
!  #ifdef DEBUG
!    write(111,'(3F12.4)') t,x,v ! DEBUG traj
!  #endif
!    !traj(0,1)=x
!    !traj(0,2)=v
!    !traj(0,3:5)=0.d0
!    !
!    ix=int((x-xmin)/dx)+1
!    iv=int((v-vmin)/dv)+1
!    it=nint(t/dt) ! important to use nint!
!    if ((ix>=1.and.ix<=nx).and.(iv>=1.and.iv<=nx).and.(it>=0.and.it<=nt)) then
!      Pmod(ix,iv,it) = Pmod(ix,iv,it) + 1.D0
!      !debug write(*,*) ix,it,t,x
!    endif
!    !!!!!!!!!!!!!!!! run one trajectory
!    do i = 1,nstep ! main loop
!      !
!      t = t+dtint
!      !
!      igrid1=int((x-xmin)/dxgrid)+1
!      igrid2=igrid1+1
!      !if (igrid1>=1.and.igrid2<=ngrid) then ! TODO avoid this to speed-up
!        force=(prof_F(igrid1)-prof_F(igrid2))/dxgrid ! dim = m*s/t**2
!      !endif
!      ! 
!      call noise(G)
!      mynoise = noisefac*G ! dim = s/t
!      ! 
!      xnew = x + dtint*v
!      vnew = v + dtint*(force/mass -gamm*v) + mynoise
!      x = xnew
!      v = vnew
!      ! TODO: here is arbitrary, to avoid the cost of "if"
!      x=max(x,xmin)
!      x=min(x,xmax)
!      v=max(v,vmin)
!      v=min(v,vmax)
!      !
!  #ifdef DEBUG
!      write(111,'(3F12.4)') t,x,v ! DEBUG traj
!  #endif
!      !
!      ix=int((x-xmin)/dx)+1
!      iv=int((v-vmin)/dv)+1
!      !it=nint(t/dt) ! important to use nint!
!      it=i ! in general it could be better to use dtint < dt, but for now we use the same
!      ! if ((ix>=1.and.ix<=nx).and.(it>=0.and.it<=nt)) then ! TODO avoid this to speed-up
!      !  if (iv>=1.and.iv<=nx) then ! TODO avoid this to speed-up
!          Pmod(ix,iv,it) = Pmod(ix,iv,it) + 1.D0
!      !  endif
!        !debug write(*,*) ix,it,t,x
!      !endif
!      ! 
!    enddo ! main loop
!  !
!  end subroutine Langevin_traj_std_fixedgamma
!================================================================================
!old!======================================================================
!old!  subroutine Langevin_traj_GLEexp_fixedgamma
!old!  !
!old!    use common_var
!old!    !
!old!    implicit none
!old!    integer :: ix,iv,it,i,igrid1,igrid2,nstep
!old!    double precision :: dtint,x,t,v,f,R,G,force,xnew,vnew,fnew,Rnew
!old!    double precision :: fac1,fac2,fac3
!old!    !double precision :: traj(nstep,5)
!old!    !---------------------------------------------------------------------------------- 
!old!    ! We adopt the following notation:
!old!    !  m*a = F -m*int_0^inf(dt' K(t')*v(t-t')) + R(t)     (K=Gamma)
!old!    !  <R(0)*R(t)> = m*k*T*K(t), K(t) = (g/tau)*exp(-t/tau)
!old!    !  int_0^inf(dt K(t)) = g  (gamma in memory-less Langevin equation, an inverse time)
!old!    ! We use the integrator in Fox et al., Phys Rev A 38, 5938 (1988):
!old!    !  (see pdf...) 
!old!    !---------------------------------------------------------------------------------- 
!old!    dtint    = dt ! TODO: how to choose? 
!old!    nstep    = nt ! nint(tmax/dtint) !
!old!    fac1     = mass*gamm*dtint/taug
!old!    fac2     = dsqrt((mass*gamm*kT/taug)*(1.d0-dexp(-2.d0*dtint/taug)))
!old!    fac3     = dexp(-dtint/taug)
!old!    !!! init !!!!
!old!    t        = 0.D0
!old!    x        = x0        ! initial position: saddle point 
!old!    call noise(G)
!old!    v        = dsqrt(kT/mass)*G ! initial velocity: Boltzmann distribution
!old!    f        = -fac1*v
!old!    call noise(G)
!old!    R        = dsqrt(mass*gamm*kT/taug)*G
!old!    !!!!!!!!!!!!!!
!old!  #ifdef DEBUG
!old!    write(111,'(A)') "# t,x,v,force,friction,noise" ! DEBUG traj
!old!    write(111,'(6E14.5)') t,x,v,force,f,R ! DEBUG traj
!old!  #endif
!old!    !
!old!    ix=int((x-xmin)/dx)+1
!old!    iv=int((v-vmin)/dv)+1
!old!    it=nint(t/dt) ! important to use nint!
!old!    if ((ix>=1.and.ix<=nx).and.(iv>=1.and.iv<=nx).and.(it>=0.and.it<=nt)) then
!old!      Pmod(ix,iv,it) = Pmod(ix,iv,it) + 1.D0
!old!      !debug write(*,*) ix,it,t,x
!old!    endif
!old!    !!!!!!!!!!!!!!!! run one trajectory
!old!    do i = 1,nstep ! main loop
!old!      !
!old!      t = t+dtint
!old!      !
!old!      igrid1=int((x-xmin)/dxgrid)+1
!old!      igrid2=igrid1+1
!old!      !if (igrid1>=1.and.igrid2<=ngrid) then ! TODO avoid this to speed-up
!old!        force=(prof_F(igrid1)-prof_F(igrid2))/dxgrid ! dim = m*s/t**2
!old!      !endif
!old!      ! 
!old!      call noise(G)
!old!      ! 
!old!      xnew = x + dtint*v
!old!      vnew = v + dtint*(force + f + R)/mass
!old!      fnew = f*fac3 - fac1*v
!old!      Rnew = R*fac3 + fac2*G
!old!      x = xnew
!old!      v = vnew
!old!      f = fnew
!old!      R = Rnew
!old!      ! TODO: here is arbitrary, to avoid the cost of "if"
!old!      x=max(x,xmin)
!old!      x=min(x,xmax)
!old!      v=max(v,vmin)
!old!      v=min(v,vmax)
!old!      !
!old!  #ifdef DEBUG
!old!      write(111,'(6E14.5)') t,x,v,force,f,R ! DEBUG traj
!old!  #endif
!old!      !
!old!      ix=int((x-xmin)/dx)+1
!old!      iv=int((v-vmin)/dv)+1
!old!      !it=nint(t/dt) ! important to use nint!
!old!      it=i ! in general it could be better to use dtint < dt, but for now we use the same
!old!      ! if ((ix>=1.and.ix<=nx).and.(it>=0.and.it<=nt)) then ! TODO avoid this to speed-up
!old!      !  if (iv>=1.and.iv<=nx) then ! TODO avoid this to speed-up
!old!          Pmod(ix,iv,it) = Pmod(ix,iv,it) + 1.D0
!old!      !  endif
!old!        !debug write(*,*) ix,it,t,x
!old!      !endif
!old!      ! 
!old!    enddo ! main loop
!old!  !
!old!  end subroutine Langevin_traj_GLEexp_fixedgamma
!old!======================================================================
!Fox!================================================================================
!Foxsubroutine Langevin_traj_GLEexp
!Fox!
!Fox  use common_var
!Fox  !
!Fox  implicit none
!Fox  integer :: ix,iv,it,i,igrid1,igrid2,nstep
!Fox  double precision :: dtint,x,t,v,f,R,G,force,mforce,xold,xnew,vnew,fnew,Rnew
!Fox  double precision :: fac1,fac2,fac3,tmp
!Fox  !double precision :: traj(nstep,5)
!Fox  !---------------------------------------------------------------------------------- 
!Fox  ! We adopt the following notation:
!Fox  !  m*a = F -m*int_0^inf(dt' K(t')*v(t-t')) + R(t)     (K=Gamma)
!Fox  !  <R(0)*R(t)> = m*k*T*K(t), K(t) = (g/tau)*exp(-t/tau)
!Fox  !  int_0^inf(dt K(t)) = g  (gamma in memory-less Langevin equation, an inverse time)
!Fox  ! We use the integrator in Fox et al., Phys Rev A 38, 5938 (1988):
!Fox  !  (see pdf...) 
!Fox  !---------------------------------------------------------------------------------- 
!Fox  igrid1 = int((x0-xmin)/dxgrid)+1
!Fox  gamm   = prof_g(igrid1)
!Fox  mass   = prof_m(igrid1)
!Fox  !
!Fox  dtint    = dt/dble(dtmult) ! TODO: do automatic test to see if it is sufficiently accurate...
!Fox  nstep    = nint(tmax/dtint)
!Fox  fac3     = dexp(-dtint/taug)
!Fox  tmp      = (kT/taug)*(1.d0-dexp(-2.d0*dtint/taug))
!Fox  !!! init !!!!
!Fox  t        = 0.D0
!Fox  x        = x0        ! initial position: saddle point 
!Fox  call noise(G)
!Fox  v        = dsqrt(kT/mass)*G ! initial velocity: Boltzmann distribution
!Fox  f        = -fac1*v
!Fox  call noise(G)
!Fox  R        = dsqrt(mass*gamm*kT/taug)*G
!Fox  !!!!!!!!!!!!!!
!Fox#ifdef DEBUG
!Fox  !write(111,'(A)') "# t,x,v,force,friction,noise" ! DEBUG traj
!Fox  write(111,'(A)') "# with mforce!!!"
!Fox  write(111,'(E16.8,5E14.5)') t,x,v,force,f,R ! DEBUG traj
!Fox#endif
!Fox  !
!Fox  xold=x
!Fox!AAA  ix=int((x-xmin)/dx)+1
!Fox!AAA  iv=int((v-vmin)/dv)+1
!Fox!AAA  it=nint(t/dt) ! important to use nint!
!Fox!AAA  if ((ix>=1.and.ix<=nx).and.(iv>=1.and.iv<=nx).and.(it>=0.and.it<=nt)) then
!Fox!AAA    Pmod(ix,iv,it) = Pmod(ix,iv,it) + 1.D0
!Fox!AAA    !debug write(*,*) ix,it,t,x
!Fox!AAA  endif
!Fox  !!!!!!!!!!!!!!!! run one trajectory
!Fox  do i = 1,nstep ! main loop
!Fox    !
!Fox    t = t+dtint
!Fox    !
!Fox    igrid1=int((x-xmin)/dxgrid)+1
!Fox    igrid2=min(igrid1+1,ngrid)
!Fox    !if (igrid1>=1.and.igrid2<=ngrid) then ! TODO avoid this to speed-up
!Fox      force = (prof_F(igrid1)-prof_F(igrid2))/dxgrid ! dim = m*s/t**2
!Fox      mforce = (prof_m(igrid1)-prof_m(igrid2))/dxgrid ! dim = m*s/t**2
!Fox      mforce = 0.5d0*mforce*v*v
!Fox      gamm  = prof_g(igrid1)
!Fox      mass  = prof_m(igrid1)
!Fox    !endif
!Fox    ! 
!Fox    call noise(G)
!Fox    ! 
!Fox    xnew = x + dtint*v
!Fox    vnew = v + dtint*(force+mforce + f + R)/mass
!Fox    fac1 = mass*gamm*dtint/taug
!Fox    fac2 = dsqrt(mass*gamm*tmp)
!Fox    fnew = f*fac3 - fac1*v
!Fox    Rnew = R*fac3 + fac2*G
!Fox    x = xnew
!Fox    v = vnew
!Fox    f = fnew
!Fox    R = Rnew
!Fox    ! TODO: here is arbitrary, to avoid the cost of "if"
!Fox    x=max(x,xmin)
!Fox    x=min(x,xmax)
!Fox!    v=max(v,vmin) ! careful: you need wide enough boundaries...
!Fox!    v=min(v,vmax)
!Fox    !
!Fox#ifdef DEBUG
!Fox    if (mod(i,dtmult).eq.0) write(111,'(E16.8,5E14.5)') t,x,v,force,f,R ! DEBUG traj
!Fox#endif
!Fox    !
!Fox    if (mod(i,dtmult)==0) then 
!Fox!AAA      ix=int((x-xmin)/dx)+1
!Fox      ix=int((xold-xmin)/dx)+1
!Fox      iv=int(((x-xold)/dt-vmin)/dv)+1
!Fox      iv=min(iv,nx)
!Fox      iv=max(iv,1)
!Fox      xold=x
!Fox    ! if ((ix>=1.and.ix<=nx).and.(it>=0.and.it<=nt)) then ! TODO avoid this to speed-up
!Fox    !  if (iv>=1.and.iv<=nx) then ! TODO avoid this to speed-up
!Fox!AAA        it=nint(t/dt) 
!Fox        it=nint(t/dt)-1 
!Fox        Pmod(ix,iv,it) = Pmod(ix,iv,it) + 1.D0
!Fox    !  endif
!Fox      !debug write(*,*) ix,it,t,x
!Fox    !endif
!Fox    endif
!Fox    ! 
!Fox  enddo ! main loop
!Fox!
!Foxend subroutine Langevin_traj_GLEexp
!Fox!======================================================================
