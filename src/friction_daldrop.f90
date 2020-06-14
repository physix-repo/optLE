!================================================================================
subroutine friction_daldrop
!
  use common_var
  !
  implicit none
  double precision :: tmp1,tmp2
  integer :: ntcorr, ntcorrave, ncorrave=10, i,j,igrid1,igrid2
  double precision :: v0,vt,at,ft
  double precision, allocatable :: corrvv(:),corrva(:),corrvf(:), cv(:,:)
  double precision :: corraa0,corraf0, dtcorr
  !
!TODO: set constant mass and variable free eneegy based on RESTART, determine constant gamma
!!!  open(111,file="colvar_daldrop",status="old",err=333)
!!!  ! if the file exists, we estimate friction as in Daldrop-PNAS-18
!!!  write(*,*) ""
!!!  write(*,*) "------------ estimating friction as in Daldrop-PNAS-18 ------------------"
!!!  write(*,*) ""
!!!  !
!!!  ntcorr=0
!!!  do
!!!    read(111,'(A100)',end=332) tmp1,tmp2
!!!    ntcorr=ntcorr+1
!!!  enddo
!!!  332 continue
!!!  rewind(111)
!!!  ntcorrave=(ntcorr-3)/ncorrave ! to estimate v and a without troubles...
!!!  allocate(cv(ntcorr,2))
!!!  do j=1,ntcorr
!!!    read(111,*) cv(j,1),cv(j,2)
!!!  enddo
!!!  close(111)
!!!  !
!!!  allocate(corrvv(0:ntcorr),corrva(0:ntcorr),corrvf(0:ntcorr))
!!!  corrvv=0.d0
!!!  corrva=0.d0
!!!  corrvf=0.d0
!!!  corraa0=0.d0
!!!  corraf0=0.d0
!!!  dtcorr=cv(2,1)-cv(1,1)
!!!  !
!!!  write(*,'(A,I9)')    " time points in colvar_daldrop         = ",ntcorr
!!!  write(*,'(A,I9)')    " time points for correlation functions = ",ntcorrave
!!!  write(*,'(A,I9)')    " number of averages                    = ",ncorrave
!!!  write(*,'(A,F12.6)') " time step                             = ",dtcorr
!!!  !
!!!  do i=1,ncorrave
!!!    t1=2+(i-1)*ntcorrave
!!!    do j=0,ntcorrave
!!!      t2=t1+j
!!!      ! computing the correlation functions:
!!!      v0=( cv(t1+1,2)-cv(t1,2) )/dt
!!!      vt=( cv(t2+1,2)-cv(t2,2) )/dt
!!!      at=( cv(t2-1,2)-2.d0*cv(t2,2)+cv(t2+1,2) )/(dt*dt)
!!!      igrid1=int((x-xmin)/dxgrid)+1
!!!      igrid2=igrid1+1
!!!      ft = (prof(1,igrid1)-prof(1,igrid2))/dxgrid ! dim = m*s/t**2
!!!      !gamm  = prof(2,igrid1)
!!!      !mass  = prof(3,igrid1)
!!!      corrvv(j)=corrvv(j)+v0*vt
!!!      corrva(j)=corrva(j)+v0*at
!!!      corrvf(j)=corrvf(j)+v0*ft
!!!      endif
!!!    endif
!!!  enddo
!!!  !
!!!  corrvv(:)=corrvv(:)/ncorrave
!!!  corrva(:)=corrva(:)/ncorrave
!!!  corrvf(:)=corrvf(:)/ncorrave
!!!  !
!!!  ! TODO: finish... what about position dependence of gamma and mass ?
!!!  !
!!!  deallocate(corrvv,corrva,corrvf)
!!!  !
!!!  ! nothing is done if file "colvar_daldrop" does not exist
!!!  333 continue
!
end subroutine friction_daldrop
!================================================================================
