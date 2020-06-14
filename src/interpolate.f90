!================================================================================
subroutine interpolate_potential_simplecos(control_p)
!
! just paste together cos curves so that F' is flat in the two minima and saddle
! no problem if points are non-uniform
!
  use common_var
  !
  implicit none
  double precision :: control_p(2,5)
  integer :: i
  double precision :: x,s
  double precision, parameter :: pi=4.D0*datan(1.D0)
  !
  prof_F(:)=0.d0
  do i=1,ngrid
    x=xmin+dble(i-1)*dxgrid
    if (x<control_p(1,2)) then
      s=(x-control_p(1,1))/(control_p(1,2)-control_p(1,1))
      prof_F(i)=control_p(2,2)+(control_p(2,1)-control_p(2,2))*(dcos((s+1.d0)*pi*0.5d0)+1.d0)
    elseif (x<control_p(1,3)) then
      s=(x-control_p(1,2))/(control_p(1,3)-control_p(1,2))
      prof_F(i)=control_p(2,2)+(control_p(2,3)-control_p(2,2))*0.5d0*(dcos((s+1.d0)*pi)+1.d0)
    elseif (x<control_p(1,4)) then
      s=(x-control_p(1,3))/(control_p(1,4)-control_p(1,3))
      prof_F(i)=control_p(2,4)+(control_p(2,3)-control_p(2,4))*0.5d0*(dcos(s*pi)+1.d0)
    else     
      s=(x-control_p(1,4))/(control_p(1,5)-control_p(1,4))
      prof_F(i)=control_p(2,4)+(control_p(2,5)-control_p(2,4))*(dcos((s+2.d0)*pi*0.5d0)+1.d0)
    endif
  enddo
  !
  prof_F(:)=prof_F(:)*kT
!
end subroutine interpolate_potential_simplecos
! !================================================================================
! subroutine interpolate_prefsplines(ic)
! !
!   use common_var
!   !
!   implicit none
!   integer :: ic ! 1 = free energy, 2 = friction, 3 = mass
!   integer :: i,j,jleft,n
!   double precision :: x,ff,a,a2,a3
!   double precision, allocatable :: c(:)
!   !
!   ! storing (prefiltered?) coefficients, extending at borders by linear extrapolation
!   n=ncontrol_p(ic)
!   allocate(c(0:n+2))
!   c(1:n)=control_p(ic,2,1:n)
!   c(0)=control_p(ic,2,1)+(control_p(ic,2,1)-control_p(ic,2,2))
!   c(n+1)=control_p(ic,2,n)+(control_p(ic,2,n)-control_p(ic,2,n-1))
!   c(n+2)=control_p(ic,2,n)+2.d0*(control_p(ic,2,n)-control_p(ic,2,n-1))
!   !
!   ! fill free energy or friction or mass profile on grid
!   prof(ic,:)=0.d0
!   ! 
!   do i=1,ngrid
!     x=xmin+dble(i-1)*dxgrid
!     ! find interval
!     do j=1,ncontrol_p(ic)-1
!       if (x>=control_p(ic,1,j).and.x<control_p(ic,1,j+1)) then
!         a=(x-control_p(ic,1,j))/(control_p(ic,1,j+1)-control_p(ic,1,j)) ! fractional pos btw nodes
!         a2=a*a
!         a3=a*a*a
!         jleft=j
!       endif
!     enddo
!     ! interpolate 
!     ff=0.d0
!     ff=ff+c(jleft-1)*(-a3+3.d0*a2-3.d0*a+1.d0)
!     ff=ff+c(jleft)*(3.d0*a3-6.d0*a2+4.d0)
!     ff=ff+c(jleft+1)*(-3.d0*a3+3.d0*a2+3.d0*a+1.d0)
!     ff=ff+c(jleft+2)*a3
!     ff=ff/6.D0
!     !
!     prof(ic,i)=ff
!   enddo
!   !
!   deallocate(c)
!   !
!   if (ic.eq.1) prof(1,:)=prof(1,:)*kT
!   !
! !!    if (x<control_p(2,1)) then
! !!      y1=control_p(1,2)+(control_p(1,2)-control_p(2,2)) ! extrapolate left
! !!      y2=control_p(1,2)
! !!      y3=control_p(2,2)
! !!      xfrac=(x-control_p(1,1))/(control_p(2,1)-control_p(1,1))
! !!    elseif (x<control_p(3,1)) then
! !!      y1=control_p(1,2)
! !!      y2=control_p(2,2)
! !!      y3=control_p(3,2)
! !!      xfrac=(x-control_p(2,1))/(control_p(3,1)-control_p(2,1))
! !!    elseif (x<control_p(4,1)) then
! !!      y1=control_p(2,2)
! !!      y2=control_p(3,2)
! !!      y3=control_p(4,2)
! !!      xfrac=(x-control_p(3,1))/(control_p(4,1)-control_p(3,1))
! !!    else
! !!      y1=control_p(3,2)
! !!      y2=control_p(4,2)
! !!      y3=control_p(5,2) ! ??? extrapolate right ???
! !!      xfrac=(x-control_p(4,1))/(control_p(5,1)-control_p(4,1))
! !!    endif
! !!    call uniform_quadratic_B_spline(xfrac,y1,y2,y3,y,yp)
! !!    prof(1,i)=y
! !!  enddo
! !
! end subroutine interpolate_prefsplines
! !================================================================================
! !================================================================================
! ! subroutine potential_Bsplines(x,F)
! ! !
! ! real*8 :: x,F
! ! !
! ! do i=1,n
! !   do j=0,9
! !     tval=dble(j)*0.1d0
! !     ! TODO: at the borders, add one point left and right extrapolating the neighbors...
! !     call uniform_quadratic_B_spline (tval,y(i-1),y(i),y(i+1),yval,ypval)
! !     write(103,*) tval+dble(i)-0.5d0,yval,ypval
! !   enddo
! ! enddo
! ! !
! ! end subroutine potential_Bsplines
! !================================================================================
! subroutine uniform_quadratic_B_spline(t,y1,y2,y3,y,yp)
! !                       [  1 -2  1 ] [ y_i-1 ]
! ! S(t) = [t**2 t 1] 0.5 [ -2  2  0 ] [ y_i   ]    with t in [0,1]
! !                       [  1  1  0 ] [ y_i+1 ]
! ! = 0.5 ( t**2 ( y_i-1 - 2 y_i + y_i+1 )
! !          + t ( -2 y_i-1 + 2 y_i )
! !            + y_i-1 + y_i               ) 
! ! or, as a scalar product:
! ! S(t) = [ 0.5*t**2-t+0.5 , -t**2+t+0.5 , 0.5*t**2 ] [ y_i-1 , y_i , y_i+1 ]
! ! S'(t)= [ t-1 , -2*t+1 , t ] [ y_i-1 , y_i , y_i+1 ]
! ! note: the first derivative is a continuous piece-wise linear function
! ! therefore the second derivative is a discontinuous piece-wise constant function.
! ! t = 0.5 is the central node, and at that point:
! ! S = 0.125*y_i-1 + 0.75*y_i + 0.125*y_i+1
! ! S'= -0.5*y_i-1 + 0.5*y_i+1  
! implicit none
! real*8 :: t,y1,y2,y3,y,yp
! real*8 :: ftt,ft
! ftt=y1-2.d0*y2+y3
! ft=2.d0*(-y1+y2)
! y = 0.5d0*( t*t*ftt + t*ft + y1+y2 )
! yp = t*ftt + 0.5d0*ft
! end subroutine uniform_quadratic_B_spline
! !================================================================================
! subroutine prefilter_splines(ic)
! !
!   use common_var
!   !
!   implicit none
!   double precision :: pole,lambda,s,zk, c(1000), ctmp(1000), delta
!   integer :: i,n,no,ic
!   !
!   pole=(sqrt(3.d0)-2.d0)
!   lambda=6.d0
!   ! artificial: linear extrapolation to zero at left and right
!   no=ncontrol_p(ic)
!   n=ncontrol_p(ic)+6
!   c(4:no+3)=control_p(ic,2,1:no) ! true data
!   c(3)=c(4)*2.d0/3.d0
!   c(2)=c(4)*1.d0/3.d0
!   c(1)=0.d0
!   c(no+4)=c(no+3)*2.d0/3.d0
!   c(no+5)=c(no+3)*1.d0/3.d0
!   c(no+6)=0.d0
!   ctmp(1:no+6)=c(1:no+6)
!   
!   ! causal initialization
!   ! note: here boundary conditions are "f(i<1)=f(i>n)=0"
!   ! (in reality it interpolates well at the borders only if f(i)=f(n)=0)
!   !*** c(1)=0.d0
!   ! boundary conditions: f(i<1)=f(i), f(i>n)=f(n) (Ruijters and Thevenaz...)
!   c(1)=0.d0
!   do i=1,n
!     c(1)=c(1)+ctmp(i)*(pole**i+pole**(2*n-(i-1)))
!   enddo
!   c(1)=c(1)/(1.D0-pole**(2*n))
!   c(1)=(c(1)+ctmp(1))*lambda
!   ! causal filtering
!   do i=2,n
!     c(i)=lambda*c(i)+pole*c(i-1)
!   enddo
!   ! anticausal initialization
!   !***  c(n)=c(n)*pole/(pole**2-1.d0)
!   c(n)=c(n)*pole/(pole-1.d0)
!   ! anticausal filtering
!   do i=n-1,1,-1
!     c(i)=pole*(c(i+1)-c(i))
!   enddo
!   !
!   control_p(ic,2,1:no)=c(4:no+3)
! !
! end subroutine prefilter_splines
! !================================================================================

