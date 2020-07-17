!================================================================================
subroutine optimize_Pmod
!
  use common_var
  !
  implicit none
  integer :: i,j,iopt,ix,iv,it,imove,nacc,nacc_last,iu,ifit,jfit
  double precision :: err1,err2,minerr1,minerr2, opt_temp 
  double precision :: err,err_old,err_best,gamm_old,gamm_best
  double precision :: mass_old,mass_best,taug_old,taug_best
  double precision :: r,rr,rrr,x,deltae
  double precision :: scale_t,scale_f,scale_g,scale_m
  double precision :: rand_taug, rand_global, rand_global_scaleall
  double precision :: prob_taug, prob_global, prob_global_scaleall
  integer, parameter :: nbestF=50
  double precision :: bestF(nbestF+2,0:ngrid),deltaf
  character :: acc*2,mmmm*4
  logical :: gaussian
  double precision :: prof_F_old(ngrid),prof_g_old(ngrid),prof_m_old(ngrid)
  double precision :: prof_F_best(ngrid),prof_g_best(ngrid),prof_m_best(ngrid)
  !
#ifdef DEBUG
  write(*,*) "DEBUG: writing out one trajectory in fort.111, without changing parameters"
  call compute_Pmod
  stop
#endif
  write(*,*) ""
  write(*,*) "----------------------- optimizing Langevin model -----------------"
  write(*,*) ""
  !
  ! initial error:
  print_traj=.false.
  !
  if (type_error==1.or.type_error==2) then
    call compute_typical_min_error(minerr1,minerr2)
  endif
  !
  call compute_Pmod ! run ntraj_Langevin simulations and compute Pmod
  call compute_error(err_old,type_error)
  write(*,'(A,E13.4)') " initial error = ",err_old
  !
  prof_F_old  = prof_F ! note: prof_F has units, is not divided by kT
  prof_g_old  = prof_g
  prof_m_old  = prof_m
  prof_F_best = prof_F
  prof_g_best = prof_g
  prof_m_best = prof_m
  err_best    = err_old
  taug_old    = taug
  taug_best   = taug
  !
  bestF(:,:)=0.d0
  bestF(:,0)=1.d10
  !
  scale_f = 0.3d0
  scale_g = 0.3d0
  scale_m = 0.3d0
  scale_t = 0.1d0
  write(*,'(A,4F9.6)') " maximum scaling for F,gamma,m,taug random moves =", &
    scale_f,scale_g,scale_m,scale_t
  !
  prob_taug            = 0.1d0 ! formerly 0.02
  prob_global          = 0.20d0
  prob_global_scaleall = 0.80d0
  if (fit_taug.eq.1) write(*,'(A,F9.6)') " probability of attempting to move taug =",prob_taug
  if (fit_F+fit_gamm+fit_mass>0) then
    write(*,'(A,F9.6)') &
     " probability of attempting global vs individual F,gamma,m moves =",prob_global
    write(*,'(A,F9.6)') &
     " probability of rigidly scaling whole profiles in global moves  =",prob_global_scaleall
  endif
  !
  nacc=0         ! total n of accepted moves
  nacc_last=0    ! n of accepted moves in the last 100 steps
  !
  write(*,*) "########################################################################"
  write(*,*) "MAIN LOOP: iopt, err(-logL,RMSD), acc, gamm, taug, mass, moved_par, optT"
  write(*,*) "########################################################################"
  !
  gaussian=.true.
  opt_temp=opt_temp1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do iopt=1,opt_niter ! opt loop
    !
    ! linear interpolation of T between T1 and T2
    !opt_temp=opt_temp1+(opt_temp2-opt_temp1)*dble(iopt)/dble(opt_niter)
    opt_temp=opt_temp+(opt_temp2-opt_temp)/dble(opt_niter-iopt)

    ! reset parameters to last accepted values
    prof_F = prof_F_old 
    prof_g = prof_g_old 
    prof_m = prof_m_old 
    taug   = taug_old          ! time constant of exp. memory function
    call update_prof_force ! important, every time prof_F is changed!

    ! now move randomly the parameters
    mmmm="____"

    call random_number(rand_taug)
    if (fit_taug.eq.1.and.rand_taug<prob_taug) then  ! change taug (a rare move...)
      call scaling_factor(gaussian,scale_t,r)
      taug=taug_old*r
      write(mmmm(3:3),'(A)') "t"
    else               ! change F, gamma, m

      call random_number(rand_global) ! this allows simultaneous moves of F,gamma,m
      call random_number(rand_global_scaleall) ! scaling of full profiles, not single control points
      call random_number(rr)
      ifit=int(rr*(fit_F+fit_gamm+fit_mass))+1
      jfit=1

      if (fit_F.eq.1) then
        if ((ifit.eq.jfit).or.(rand_global<prob_global)) then
          if (rand_global_scaleall>prob_global_scaleall) then            ! local move
            call add_gaussian(1,max_Gaussian_h(1),max_Gaussian_w(1))
          else                        ! global move
            call scaling_factor(gaussian,scale_f,r)
            prof_F=prof_F*r
          endif
          write(mmmm(1:1),'(A)') "f"
          call update_prof_force ! important, every time prof_F is changed!
        endif
        jfit=jfit+1
      endif
      !
      if (fit_gamm.eq.1) then
        if ((ifit.eq.jfit).or.(rand_global<prob_global)) then
          if (pos_dep_gamma.eq.0) then
            call scaling_factor(gaussian,scale_g,r)
            prof_g=prof_g*r
          else
            if (rand_global_scaleall>prob_global_scaleall) then            ! local move
              call add_gaussian(2,max_Gaussian_h(2),max_Gaussian_w(2))
            else                        ! global move
              call scaling_factor(gaussian,scale_g,r)
              prof_g=prof_g*r
            endif
          endif
          write(mmmm(2:2),'(A)') "g"
        endif
        jfit=jfit+1
      endif
      !
      if (fit_mass.eq.1) then
        if ((ifit.eq.jfit).or.(rand_global<prob_global)) then
          if (pos_dep_mass.eq.0) then
            call scaling_factor(gaussian,scale_m,r)
            prof_m=prof_m*r
          else
            if (rand_global_scaleall>prob_global_scaleall) then            ! local move
              call add_gaussian(3,max_Gaussian_h(3),max_Gaussian_w(3))
            else                        ! global move
              call scaling_factor(gaussian,scale_m,r)
              prof_m=prof_m*r
            endif
          endif
          write(mmmm(4:4),'(A)') "m"
        endif
        jfit=jfit+1
      endif

    endif
    !
    if (fix_mass0.eq.1) prof_m = mass0*prof_m/prof_m(int((x0-xmin)/dxgrid)+1)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call compute_Pmod ! run ntraj_Langevin simulations and compute Pmod
    !
    ! TODO: you can save time by computing only the necessary type of error!
    if (type_error==1.or.type_error==2) then 
      call compute_error(err1,1)
      call compute_error(err2,2)
      if (type_error.eq.1) err=err1
      if (type_error.eq.2) err=err2
    elseif (type_error==3) then
      call compute_error(err,3)
      err1=err
      err2=0.d0
    endif
    !
    ! store the best F profiles in memory
    do i=1,nbestF
      if (err<bestF(i,0)) then
        bestF(i,0)=err
        bestF(i,1:ngrid)=prof_F(1:ngrid)
        exit
      endif
    enddo
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! acceptance test
    acc="N "
    deltae=err-err_old
    call random_number(r)
    ! OLD VERSION WITH RELATIVE ERROR: if (r<dexp(-(err-err_old)/(err_old*opt_temp))) then
    if (r<dexp(-deltae/opt_temp)) then ! accepted
      err_old =err
      taug_old=taug
      prof_F_old  = prof_F
      prof_g_old  = prof_g
      prof_m_old  = prof_m
      acc=" Y"
      nacc=nacc+1
      nacc_last=nacc_last+1
      !
      if (err<err_best) then
        ! write profiles only if err is the best so far
        err_best =err
        taug_best=taug
        prof_F_best  = prof_F
        prof_g_best  = prof_g
        prof_m_best  = prof_m
        !
        iu=5000000+iopt
        open(iu,status="unknown")
        write(iu,'(A,2E13.5)') "# x F F/kT gamma mass ; taug, err=",taug,err
        do i=1,ngrid
          x=xmin+dble(i-1)*dxgrid
          write(iu,'(5E13.4)') x,prof_F(i)-minval(prof_F(:)), &
           (prof_F(i)-minval(prof_F(:)))/kT,prof_g(i),prof_m(i)
        enddo
        close(iu)
        !
        !open(33,file="Pmod",status="unknown")
        !do ix=1,nx
        !  do iv=1,nx
        !    do it=0,nt
        !      write(33,'(4F11.6)') xmin+(dble(ix)-1.d0)*dx,vmin+(dble(iv)-1.d0)*dv,dble(it)*dt,Pmod(ix,iv,it)
        !    enddo
        !  enddo
        !enddo
        !close(33)
      endif
      !
    endif
    !
    if (acc.eq." Y") then
      write(*,'(I8,2x,2E13.5,3X,A,3E12.4,3X,A4,3X,E11.4)') &
       iopt,err1,err2,acc,sum(prof_g(:))/ngrid,taug,sum(prof_m(:))/ngrid,mmmm,opt_temp
    endif
    !
    if (mod(iopt,100).eq.0) then
      if (target_acc.gt.0.d0.and.target_acc.lt.1.d0) then
        if (nacc_last.lt.(target_acc*70.d0))  opt_temp=opt_temp*1.2d0
        if (nacc_last.gt.(target_acc*130.d0)) opt_temp=opt_temp*0.8d0
        write(*,'(A,I3,A,I3,A,E11.4)') &
         "### acc.moves/100= ",nacc_last," target= ",nint(target_acc*100)," newT= ",opt_temp
        nacc_last=0
      endif
    endif
    !
    call flush()
  enddo ! opt loop
  !
  write(*,*) ""
  write(*,*) "* * * * * * * * * *  final results  * * * * * * * * * *"
  write(*,*) ""
  prof_F = prof_F_best 
  prof_g = prof_g_best 
  prof_m = prof_m_best 
  write(*,'(A,F14.8,3E12.4)') " best error, average mass, average gamma, tau =",&
   err_best,sum(prof_m(:))/ngrid,sum(prof_g(:))/ngrid,taug_best
  write(*,'(A,F14.8)') " fraction of accepted moves = ",dble(nacc)/dble(opt_niter)
  !
  if (fix_mass0.eq.1) prof_m(:) = mass0*prof_m(:)/prof_m(int((x0-xmin)/dxgrid)+1)
  !
  call print_profiles
  !
  write(*,*) "errors corresponding to the",nbestF," best free energy profiles:"
  write(*,'(10E13.5)') bestF(1:nbestF,0)
  bestF(nbestF+1,1:ngrid)=0.d0 ! here we store the average
  do i=1,nbestF
    x=sum(bestF(i,1:ngrid))/dble(ngrid)
    bestF(i,1:ngrid)=bestF(i,1:ngrid)-x
  enddo
  do i=2,ngrid
    deltaf=sum(bestF(1:nbestF,i)-bestF(1:nbestF,i-1))/dble(nbestF)
    bestF(nbestF+1,i)=bestF(nbestF+1,i-1)+deltaf ! average profile, as integral of average derivative
  enddo
  bestF(nbestF+1,1:ngrid)=bestF(nbestF+1,1:ngrid)-sum(bestF(nbestF+1,1:ngrid))/dble(ngrid)
  do i=2,ngrid
    bestF(nbestF+2,i)=sqrt(sum((bestF(1:nbestF,i)-bestF(nbestF+1,i))**2)/dble(nbestF)) ! RMSD
  enddo
  bestF(nbestF+1,1:ngrid)=bestF(nbestF+1,1:ngrid)-minval(bestF(nbestF+1,1:ngrid))
  open(77,file="AVERF",status="unknown")
  write(77,'(A)') "# x averF RMSD averF(kT) RMSD(kT)"
  x=0.d0
  do i=1,ngrid
    x=xmin+dble(i-1)*dxgrid
    write(77,'(5F16.8)') x,bestF(nbestF+1,i),bestF(nbestF+2,i),bestF(nbestF+1,i)/kT,bestF(nbestF+2,i)/kT
  enddo
  close(77)
  write(*,*) "written average of the ",nbestF," best free energy profiles in AVERF"
  !
  print_traj=.true. ! TODO this is now disabled...
  !
  call compute_Pmod ! run ntraj_Langevin simulations and compute Pmod
  !write(*,*) "final trajectories written in fort.111"
  if (type_error.ne.3) then
    open(33,file="Pmod",status="unknown")
    open(34,file="Pdiff",status="unknown")
    do ix=1,nx
      do iv=1,nx
        do it=0,nt
          write(33,'(4E11.3)') xmin+(dble(ix)-1.d0)*dx,vmin+(dble(iv)-1.d0)*dv,dble(it)*dt,Pmod(ix,iv,it)
          write(34,'(4E11.3)') xmin+(dble(ix)-1.d0)*dx,vmin+(dble(iv)-1.d0)*dv,dble(it)*dt,Pref(ix,iv,it)-Pmod(ix,iv,it)
        enddo
      enddo
    enddo
    close(33)
    close(34)
    write(*,*) "written Pmod and Pdiff"
  endif
  !
  ! clean memory
  deallocate(colvar)
  deallocate(Pref,Pmod)
!
end subroutine optimize_Pmod
!================================================================================
subroutine compute_Pmod
!
  use common_var
  !
  implicit none
  integer :: itraj_MD,imult,ix,it
  !
  Pmod=0.D0
  !
  if (type_Langevin.eq.0) then
    if (type_error.ne.3) then 
      do itraj_MD=1,ntraj_MD
        x0now=x0(itraj_MD)
        do imult=1,ratio_Langevin_MD
          call Langevin_traj_overdamped   ! this updates Pmod with one traj
        enddo
      enddo
    else
      ! likelihood of MD trajectory from overdamped propagator
      continue 
    endif
  elseif (type_Langevin.eq.1) then
    do itraj_MD=1,ntraj_MD
      x0now=x0(itraj_MD)
      do imult=1,ratio_Langevin_MD
        call Langevin_traj_std          ! this updates Pmod with one traj
      enddo
    enddo
  elseif (type_Langevin.eq.2) then
    do itraj_MD=1,ntraj_MD
      x0now=x0(itraj_MD)
      do imult=1,ratio_Langevin_MD
        call Langevin_traj_GLEexp       ! this updates Pmod with one traj
      enddo
    enddo
  endif
  !
  Pmod=Pmod/sum(Pmod(:,:,0))
!
end subroutine compute_Pmod
!================================================================================
subroutine scaling_factor(gaussian,s,r)
!
  implicit none
  logical :: gaussian
  double precision :: s,r
  !
  if (gaussian) then
    call noise(r)
    if (r<-3.d0) r=-3.d0  ! prevent extreme changes 
    if (r> 3.d0) r= 3.d0
    r=1.d0+r*s
  else
    call random_number(r)
    r=(1.d0-s)+r*s*2.d0
  endif
!
end subroutine scaling_factor
!================================================================================
subroutine add_gaussian(iprof,max_height,max_width)
!
  use common_var
  !
  implicit none
  integer :: iprof,i,deltai,pos,imin,imax
  double precision :: r,width,height,max_height,max_width,profmin
  double precision, parameter :: shift=dexp(-4.d0*4.d0/2.d0) ! about 0.0003
  !
  ! --- set Gaussian parameters
  call random_number(r)
  pos    = int(r*ngrid)+1
  call random_number(r)
  r = max(max_width,r) ! arbitrary, to avoid Gaussians smaller than xrange/25
  width  = ngrid*r
  call random_number(r)
  height = max_height*(2.d0*r-1.d0) ! positive or negative
  !
  ! --- add Gaussian (truncated and shifted, for continuity) to profile
  deltai = int(width*4.d0/dxgrid)
  imin   = max(1,pos-deltai)
  imax   = min(ngrid,pos+deltai)
  if (iprof==1) then ! free energy
    do i=imin,imax
      prof_F(i) = prof_F(i)+height*(dexp(-0.5d0*((pos-i)/width)**2)-shift)
    enddo
  endif
  if (iprof==2) then ! gamma
    profmin = minval(prof_g(imin:imax))
    if (profmin+height<=0.d0) height=-0.5d0*profmin ! avoid negative gamma
    do i=imin,imax
      prof_g(i) = prof_g(i)+height*(dexp(-0.5d0*((pos-i)/width)**2)-shift)
    enddo
  endif
  if (iprof==3) then ! mass
    profmin = minval(prof_m(imin:imax))
    if (profmin+height<=0.d0) height=-0.5d0*profmin ! avoid negative mass
    do i=imin,imax
      prof_m(i) = prof_m(i)+height*(dexp(-0.5d0*((pos-i)/width)**2)-shift)
    enddo
  endif
!
end subroutine add_gaussian
!================================================================================
