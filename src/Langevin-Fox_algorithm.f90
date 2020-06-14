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
!        force=(prof(1,igrid1)-prof(1,igrid2))/dxgrid ! dim = m*s/t**2
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
subroutine Langevin_traj_std
!
  use common_var
  !
  implicit none
  integer :: ix,iv,it,i,igrid1,igrid2,nstep
  double precision :: dtint,x,t,v,G,force,mforce,xold,xnew,vnew
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
  igrid1 = int((x0-xmin)/dxgrid)+1
  gamm   = prof(2,igrid1)
  mass   = prof(3,igrid1)
  !
  dtint    = dt/dble(dtmult) ! TODO: do automatic test to see if it is sufficiently accurate...
  nstep    = nint(tmax/dtint) 
  t        = 0.D0
  x        = x0        ! initial position: saddle point 
  call noise(G)
  v        = dsqrt(kT/mass)*G ! initial velocity: Boltzmann distribution
  tmp      = 2.D0*kT*dtint
  !
#ifdef DEBUG
  write(111,'(A)') "# with mforce!!!"
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
      force  = (prof(1,igrid1)-prof(1,igrid2))/dxgrid ! dim = m*s/t**2
      mforce = (prof(3,igrid1)-prof(3,igrid2))/dxgrid ! dim = m*s/t**2
      mforce = 0.5d0*mforce*v*v
      gamm   = prof(2,igrid1)
      mass   = prof(3,igrid1)
      noisefac = dsqrt(tmp*gamm/mass) ! dim = s/t
      !!!if (mod(i,dtmult).eq.0) write(112,'(4E14.5)') force,mforce,mass,gamm !DDDDD
    !endif
    ! 
    call noise(G) ! TODO: store in memory vector of G to accelerate?
    mynoise = noisefac*G ! dim = s/t
    ! 
    xnew = x + dtint*v
    vnew = v + dtint*((force+mforce)/mass -gamm*v) + mynoise
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
!
end subroutine Langevin_traj_std
!======================================================================
!  subroutine Langevin_traj_GLEexp_fixedgamma
!  !
!    use common_var
!    !
!    implicit none
!    integer :: ix,iv,it,i,igrid1,igrid2,nstep
!    double precision :: dtint,x,t,v,f,R,G,force,xnew,vnew,fnew,Rnew
!    double precision :: fac1,fac2,fac3
!    !double precision :: traj(nstep,5)
!    !---------------------------------------------------------------------------------- 
!    ! We adopt the following notation:
!    !  m*a = F -m*int_0^inf(dt' K(t')*v(t-t')) + R(t)     (K=Gamma)
!    !  <R(0)*R(t)> = m*k*T*K(t), K(t) = (g/tau)*exp(-t/tau)
!    !  int_0^inf(dt K(t)) = g  (gamma in memory-less Langevin equation, an inverse time)
!    ! We use the integrator in Fox et al., Phys Rev A 38, 5938 (1988):
!    !  (see pdf...) 
!    !---------------------------------------------------------------------------------- 
!    dtint    = dt ! TODO: how to choose? 
!    nstep    = nt ! nint(tmax/dtint) !
!    fac1     = mass*gamm*dtint/taug
!    fac2     = dsqrt((mass*gamm*kT/taug)*(1.d0-dexp(-2.d0*dtint/taug)))
!    fac3     = dexp(-dtint/taug)
!    !!! init !!!!
!    t        = 0.D0
!    x        = x0        ! initial position: saddle point 
!    call noise(G)
!    v        = dsqrt(kT/mass)*G ! initial velocity: Boltzmann distribution
!    f        = -fac1*v
!    call noise(G)
!    R        = dsqrt(mass*gamm*kT/taug)*G
!    !!!!!!!!!!!!!!
!  #ifdef DEBUG
!    write(111,'(A)') "# t,x,v,force,friction,noise" ! DEBUG traj
!    write(111,'(6E14.5)') t,x,v,force,f,R ! DEBUG traj
!  #endif
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
!        force=(prof(1,igrid1)-prof(1,igrid2))/dxgrid ! dim = m*s/t**2
!      !endif
!      ! 
!      call noise(G)
!      ! 
!      xnew = x + dtint*v
!      vnew = v + dtint*(force + f + R)/mass
!      fnew = f*fac3 - fac1*v
!      Rnew = R*fac3 + fac2*G
!      x = xnew
!      v = vnew
!      f = fnew
!      R = Rnew
!      ! TODO: here is arbitrary, to avoid the cost of "if"
!      x=max(x,xmin)
!      x=min(x,xmax)
!      v=max(v,vmin)
!      v=min(v,vmax)
!      !
!  #ifdef DEBUG
!      write(111,'(6E14.5)') t,x,v,force,f,R ! DEBUG traj
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
!  end subroutine Langevin_traj_GLEexp_fixedgamma
!======================================================================
subroutine Langevin_traj_GLEexp
!
  use common_var
  !
  implicit none
  integer :: ix,iv,it,i,igrid1,igrid2,nstep
  double precision :: dtint,x,t,v,f,R,G,force,mforce,xold,xnew,vnew,fnew,Rnew
  double precision :: fac1,fac2,fac3,tmp
  !double precision :: traj(nstep,5)
  !---------------------------------------------------------------------------------- 
  ! We adopt the following notation:
  !  m*a = F -m*int_0^inf(dt' K(t')*v(t-t')) + R(t)     (K=Gamma)
  !  <R(0)*R(t)> = m*k*T*K(t), K(t) = (g/tau)*exp(-t/tau)
  !  int_0^inf(dt K(t)) = g  (gamma in memory-less Langevin equation, an inverse time)
  ! We use the integrator in Fox et al., Phys Rev A 38, 5938 (1988):
  !  (see pdf...) 
  !---------------------------------------------------------------------------------- 
  igrid1 = int((x0-xmin)/dxgrid)+1
  gamm   = prof(2,igrid1)
  mass   = prof(3,igrid1)
  !
  dtint    = dt/dble(dtmult) ! TODO: do automatic test to see if it is sufficiently accurate...
  nstep    = nint(tmax/dtint)
  fac3     = dexp(-dtint/taug)
  tmp      = (kT/taug)*(1.d0-dexp(-2.d0*dtint/taug))
  !!! init !!!!
  t        = 0.D0
  x        = x0        ! initial position: saddle point 
  call noise(G)
  v        = dsqrt(kT/mass)*G ! initial velocity: Boltzmann distribution
  f        = -fac1*v
  call noise(G)
  R        = dsqrt(mass*gamm*kT/taug)*G
  !!!!!!!!!!!!!!
#ifdef DEBUG
  !write(111,'(A)') "# t,x,v,force,friction,noise" ! DEBUG traj
  write(111,'(A)') "# with mforce!!!"
  write(111,'(E16.8,5E14.5)') t,x,v,force,f,R ! DEBUG traj
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
    igrid2=min(igrid1+1,ngrid)
    !if (igrid1>=1.and.igrid2<=ngrid) then ! TODO avoid this to speed-up
      force = (prof(1,igrid1)-prof(1,igrid2))/dxgrid ! dim = m*s/t**2
      mforce = (prof(3,igrid1)-prof(3,igrid2))/dxgrid ! dim = m*s/t**2
      mforce = 0.5d0*mforce*v*v
      gamm  = prof(2,igrid1)
      mass  = prof(3,igrid1)
    !endif
    ! 
    call noise(G)
    ! 
    xnew = x + dtint*v
    vnew = v + dtint*(force+mforce + f + R)/mass
    fac1 = mass*gamm*dtint/taug
    fac2 = dsqrt(mass*gamm*tmp)
    fnew = f*fac3 - fac1*v
    Rnew = R*fac3 + fac2*G
    x = xnew
    v = vnew
    f = fnew
    R = Rnew
    ! TODO: here is arbitrary, to avoid the cost of "if"
    x=max(x,xmin)
    x=min(x,xmax)
!    v=max(v,vmin) ! careful: you need wide enough boundaries...
!    v=min(v,vmax)
    !
#ifdef DEBUG
    if (mod(i,dtmult).eq.0) write(111,'(E16.8,5E14.5)') t,x,v,force,f,R ! DEBUG traj
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
