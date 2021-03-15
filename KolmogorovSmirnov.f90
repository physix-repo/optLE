!--------------------------------------------------------
! adapted from Numerical Recipes in Fortran, 2nd edition
!--------------------------------------------------------
SUBROUTINE KSone_normal(data1,n,d,prob)
  INTEGER n
  double precision :: d,data1(n),prob !,func
  double precision, parameter :: sqrt2=1.4142135623731
  !EXTERNAL func
  !USES probKS,sort
  ! Given an array data1(1:n) , and given a user-supplied function of a single variable func
  ! which is a cumulative distribution function ranging from 0 (for smallest values of its argu-
  ! ment) to 1 (for largest values of its argument), this routine returns the K–S statistic d , and
  ! the significance level prob . Small values of prob show that the cumulative distribution
  ! function of data is significantly different from func . The array data1 is modified by being
  ! sorted into ascending order.
  INTEGER j
  double precision :: dt,en,ff,fn,fo,probKS
  !call sort(n,data1) 
  !call quicksort(data1,1,n)
  en=n
  d=0.
  fo=0.                               ! Data’s c.d.f. before the next step.
  do j=1,n                            ! Loop over the sorted data points.
    fn=j/en                             ! Data’s c.d.f. after this step.
    !ff=func(data1(j))                   ! Compare to the user-supplied function.
    ff=0.5*(erf(data1(j)/sqrt2)+1)       ! here we use the normal law CDF:
              ! 0.5*(erf(x)+1) is the CDF of exp(-x*x)/sqrt(pi) with sigma=1/sqrt(2)
              ! 0.5*(erf(x/sqrt(2))+1) is the CDF of exp(-x*x/2)/sqrt(2pi) with sigma=1
    dt=max(abs(fo-ff),abs(fn-ff))       ! Maximum distance.
    if (dt.gt.d) d=dt
    fo=fn
  enddo 
  en=sqrt(en)
  prob=probKS((en+0.12+0.11/en)*d)    ! Compute significance.
  return
END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE KStwo(data1,n1,data2,n2,d,prob)
  INTEGER n1,n2
  double precision :: d,prob,data1(n1),data2(n2)
  ! USES probKS,sort
  ! Given an array data1(1:n1) , and an array data2(1:n2) , this routine returns the K-S
  ! statistic d , and the significance level prob for the null hypothesis that the data sets
  ! are drawn from the same distribution. Small values of prob show that the cumulative
  ! distribution function of data1 is significantly different from that of data2 . The arrays
  ! data1 and data2 are modified by being sorted into ascending order.
  INTEGER j1,j2
  double precision :: d1,d2,dt,en1,en2,en,fn1,fn2,probKS
  !call sort(n1,data1) 
  !call sort(n2,data2) 
  call quicksort(data1,1,n1)
  call quicksort(data2,1,n2)
  en1=n1
  en2=n2
  j1=1    ! Next value of data1 to be processed.
  j2=1    ! Ditto, data2.
  fn1=0.
  fn2=0.
  d=0.
1 if (j1.le.n1.and.j2.le.n2) then        ! If we are not done...
    d1=data1(j1)
    d2=data2(j2)
  if(d1.le.d2)then                     ! Next step is in data1.
    fn1=j1/en1
    j1=j1+1
  endif
  if(d2.le.d1)then                     ! Next step is in data2.
    fn2=j2/en2
    j2=j2+1
  endif
  dt=abs(fn2-fn1)
  if(dt.gt.d)d=dt
    goto 1
  endif
  en=sqrt(en1*en2/(en1+en2))
  prob=probKS((en+0.12+0.11/en)*d)     ! Compute significance.
  return
END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION probKS(alam)
  double precision :: probKS,alam,EPS1,EPS2
  PARAMETER (EPS1=0.001, EPS2=1.e-8)
  ! Kolmogorov-Smirnov probability function.
  INTEGER j
  double precision :: a2,fac,term,termbf
  !
  a2=-2.*alam**2
  fac=2.
  probKS=0.
  termbf=0.                     ! Previous term in sum.
  do j=1,100
    term=fac*exp(a2*j**2)
    probKS=probKS+term
    if (abs(term).le.EPS1*termbf.or.abs(term).le.EPS2*probKS) return
    fac=-fac                      ! Alternating signs in sum.
    termbf=abs(term)
  enddo
  probKS=1.                     ! Get here only by failing to converge.
  return
END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
