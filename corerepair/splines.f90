  module splines

  use kinds, only: dp

  implicit none

  type spline_struct
    integer :: ideg
    integer :: n_knot
    real(dp),dimension(:),pointer :: t_knot => null()
    real(dp),dimension(:),pointer :: coeff => null()
  end type spline_struct

  private
  public :: curev, splint, splder, splev, fourco, curfit

  public :: spline_struct, destroy_spline_struct, create_spline_struct
  public :: spline_fit, spline_eval, spline_deriv, spline_integral

  interface spline_eval
    module procedure spline_eval_scalar
    module procedure spline_eval_array
  end interface spline_eval


  contains

  subroutine create_spline_struct( spline, n_knot, ideg, ierr )

  type(spline_struct) :: spline
  integer :: n_knot, ideg, ierr

  call destroy_spline_struct( spline )

  spline%n_knot = n_knot
  spline%ideg = ideg
  allocate( spline%t_knot(n_knot), stat=ierr )
  allocate( spline%coeff(n_knot),  stat=ierr )

  end subroutine create_spline_struct

  subroutine destroy_spline_struct( spline )

  use pointers

  type(spline_struct) :: spline
  integer :: ierr

  call pointer_destroy( spline%t_knot )
  call pointer_destroy( spline%coeff )

  !if( associated( spline%t_knot ) ) deallocate( spline%t_knot, stat=ierr )
  !nullify( spline%t_knot )
  !if( ierr/=0 ) then
  !  write(*,*) 't_knot'
  !  write(*,*) 'ierr= ', ierr
  !  write(*,*) 'associated: ', associated(spline%t_knot)
  !  stop
  !endif
  !if( associated( spline%coeff ) ) deallocate( spline%coeff, stat=ierr )
  !nullify( spline%coeff )
  !if( ierr/=0 ) then
  !  write(*,*) 'coeff'
  !  write(*,*) 'ierr= ', ierr
  !  write(*,*) 'associated: ', associated(spline%t_knot)
  !  stop
  !endif

  end subroutine destroy_spline_struct

! ====================================================================== 
! Wrapper higher level routines from D. Prendergast
! ====================================================================== 

  subroutine spline_fit( ideg, n_data, x_data, y_data, spline, ier, &
                         xb, xe )

  ! Does a straight interpolative fit passing through all data points
  ! of the degree specified by ideg. This version stops splining at the
  ! ends of the abcissa set x_data (no extrapolation) unless optional
  ! variables are set

  integer,intent(in) :: ideg, n_data
  real(dp),intent(in) :: x_data(n_data)
  real(dp),intent(in) :: y_data(n_data)
  type(spline_struct),intent(out) :: spline
  integer,intent(out) :: ier
  real(dp),intent(in),optional :: xb, xe

  integer :: nest, lwrk
  real(dp),allocatable :: weight(:), wrk(:)
  integer,allocatable :: iwrk(:)
  integer :: iopt
  real(dp) :: s, fp
  real(dp) :: xb_, xe_

  integer :: i, ixb, ixe, n_data_
  logical :: errflag

  iopt=0
  s=0.d0

  ixb=1
  if( .not.(present(xb)) ) then
    xb_=minval(x_data)
  else
    xb_=xb
    ! if the data includes abcissa before xb then adjust
    if( x_data(1) < xb_ ) then
      ixb = 2
      do while( x_data(ixb) < xb_ .and. ixb < n_data )
        ixb=ixb+1
      enddo
    endif
  endif
  if( ixb == n_data ) then
    write(*,*) 'spline_fit error : x_data all less than requested starting point'
    ier=1
    return
  endif

  ixe=n_data
  if( .not.(present(xe)) ) then 
    xe_=maxval(x_data)
  else
    xe_=xe
    ! if the data includes abcissa after xe then adjust
    if( x_data(n_data) > xe_ ) then
      ixe = n_data-1
      do while( x_data(ixe) > xe_ .and. ixe > 1 )
        ixe=ixe-1
      enddo
    endif
  endif
  if( ixe == 1 ) then
    write(*,*) 'spline_fit error : x_data all greater than requested end point'
    ier=1
    return
  endif
                                                                                
  ! redefine number of useful points
  n_data_ = ixe - ixb + 1
  if( n_data_ < 1 ) then
    write(*,*) 'spline_fit error : number of useful points < 1 ', n_data_
    ier = 1
    return
  endif

  ! initialize some variables
  nest=n_data_+ideg+1
  lwrk=n_data_*(ideg+1)+nest*(7+3*ideg)
  allocate( weight(n_data_), wrk(lwrk), iwrk(nest) )
  weight=1.d0

  call create_spline_struct( spline, nest, ideg, ier )

  if( ier /= 0 ) then
    write(*,*) 'spline_fit error : unable to allocate space for splines'
    return
  endif

  call curfit( iopt,n_data_,x_data(ixb),y_data(ixb), &
               weight,xb_,xe_,ideg,s,nest, &
               spline%n_knot,spline%t_knot,spline%coeff,           &
               fp,wrk,lwrk,iwrk,ier )

  if( ier == 10 ) then
    write(*,*) 'curfit error ', ier
    write(*,*) 'checking various causes...'

    write(*,'(2x,a)',advance='no') '-1<=iopt<=1 : '
    if( iopt >= -1 .and. iopt <= 1 ) then
      write(*,*) 'ok'
    else
      write(*,*) 'error'
      write(*,*) ' iopt = ', iopt
    endif

    write(*,'(2x,a)',advance='no') '1<=k<=5 : '
    if( ideg >= 1 .and. ideg <= 5 ) then
      write(*,*) 'ok'
    else
      write(*,*) 'error'
      write(*,*) ' k = ', ideg
    endif

    write(*,'(2x,a)',advance='no') 'm>k : '
    if( n_data_ > ideg ) then
      write(*,*) 'ok'
    else
      write(*,*) 'error'
      write(*,*) ' m = ', n_data_
      write(*,*) ' k = ', ideg
    endif

    write(*,'(2x,a)',advance='no') 'nest>2*k+2 : '
    if( nest > 2*ideg+2 ) then
      write(*,*) 'ok'
    else
      write(*,*) 'error'
      write(*,*) ' nest = ', nest
      write(*,*) ' k = ', ideg
    endif

    write(*,'(2x,a)',advance='no') 'w(i)>0,i=1,2,...,m : '
    if( all( weight > 0.d0 ) ) then
      write(*,*) 'ok'
    else
      write(*,*) 'error'
      do i=1,n_data_
        if( weight(i) > 0.d0 ) write(*,*) ' w(',i,') = ', weight(i)
      enddo
    endif
                                                                                
    write(*,'(2x,a)',advance='no') 'xb<=x(1)<x(2)<...<x(m)<=xe : '
    errflag=.false.
    if( xb_ > x_data(ixb) ) errflag=.true.
    if( xe_ < x_data(ixe) ) errflag=.true.
    do i=ixb+1,ixe
      if( x_data(i-1) >= x_data(i) ) errflag=.true.
    enddo
    if( errflag ) then
      write(*,*) 'error'
    else
      write(*,*) 'ok'
    endif
    if( xb_ > x_data(ixb) ) write(*,*) 'error xb>x(1)', xb_, x_data(ixb)
    if( xe_ < x_data(ixe) ) write(*,*) 'error xe<x(m)', xe_, x_data(ixe)
    do i=ixb+1,ixe
      if( x_data(i-1) >= x_data(i) ) write(*,*) 'error x(',i-1,')>=x(',i,')', x_data(i-1), x_data(i)
    enddo

    write(*,'(2x,a)',advance='no') 'lwrk>=(k+1)*m+nest*(7+3*k) : '
    if( lwrk >= (ideg+1)*n_data_+nest*(7+3*ideg) ) then
      write(*,*) 'ok'
    else
      write(*,*) 'error'
      write(*,*) ' lwrk = ', lwrk
      write(*,*) ' k = ', ideg
      write(*,*) ' m = ', n_data_
      write(*,*) ' nest = ', nest
    endif

    write(*,'(2x,a)',advance='no') 'iopt=0, s=0, nest >= m+k+1 : '
    if( iopt == 0 .and. s == 0.d0 .and. nest >= n_data_+ideg+1 ) then
      write(*,*) 'ok'
    else
      write(*,*) 'error'
      write(*,*) ' iopt = ', iopt
      write(*,*) ' s = ', s
      write(*,*) ' nest = ', nest
      write(*,*) ' m = ', n_data_
      write(*,*) ' k = ', ideg
    endif
     
  endif

  deallocate( weight, wrk, iwrk )

  end subroutine spline_fit


  subroutine spline_eval_array( spline, n_data, x_data, y_data, ier )

  type(spline_struct),intent(in) :: spline
  integer,intent(in) :: n_data
  real(dp),intent(in) :: x_data(n_data)
  real(dp),intent(out) :: y_data(n_data)
  integer,intent(out) :: ier

  call splev( spline%t_knot, spline%n_knot, spline%coeff, spline%ideg, &
              x_data, y_data, n_data, ier )

  end subroutine spline_eval_array

  subroutine spline_eval_scalar( spline, x_data, y_data, ier )

  type(spline_struct),intent(in) :: spline
  real(dp),intent(in) :: x_data
  real(dp),intent(out) :: y_data
  integer,intent(out) :: ier

  integer :: n
  real(dp) :: x(1), y(1)

  n=1
  x(1) = x_data
  call splev( spline%t_knot, spline%n_knot, spline%coeff, spline%ideg, &
              x, y, n, ier )
  y_data = y(1)

  end subroutine spline_eval_scalar


  subroutine spline_deriv( spline, ideriv, n_data, x_data, y_data, ier )

  type(spline_struct),intent(in) :: spline
  integer,intent(in) :: ideriv
  integer,intent(in) :: n_data
  real(dp),intent(in) :: x_data(n_data)
  real(dp),intent(out) :: y_data(n_data)
  integer,intent(out) :: ier

  real(dp),allocatable :: wrk(:)

  allocate( wrk(spline%n_knot) )

  call splder( spline%t_knot, spline%n_knot, spline%coeff, spline%ideg, &
               ideriv, x_data, y_data, n_data, wrk, ier )

  deallocate( wrk )

  end subroutine spline_deriv


  function spline_integral( spline, a, b )

  real(dp) :: spline_integral
  type(spline_struct),intent(in) :: spline
  real(dp),intent(in) :: a, b

  real(dp),allocatable :: wrk(:)

  allocate( wrk(spline%n_knot) )

  spline_integral = splint( spline%t_knot, spline%n_knot, spline%coeff, &
                            spline%ideg, a, b, wrk )

  deallocate( wrk )

  end function spline_integral


! ====================================================================== 
! Lower level routines from dierckx
! ====================================================================== 

      subroutine curev(idim,t,n,c,nc,k,u,m,x,mx,ier)
!  subroutine curev evaluates in a number of points u(i),i=1,2,...,m
!  a spline curve s(u) of degree k and dimension idim, given in its
!  b-spline representation.
!
!  calling sequence:
!     call curev(idim,t,n,c,nc,k,u,m,x,mx,ier)
!
!  input parameters:
!    idim : integer, giving the dimension of the spline curve.
!    t    : array,length n, which contains the position of the knots.
!    n    : integer, giving the total number of knots of s(u).
!    c    : array,length nc, which contains the b-spline coefficients.
!    nc   : integer, giving the total number of coefficients of s(u).
!    k    : integer, giving the degree of s(u).
!    u    : array,length m, which contains the points where s(u) must
!           be evaluated.
!    m    : integer, giving the number of points where s(u) must be
!           evaluated.
!    mx   : integer, giving the dimension of the array x. mx >= m*idim
!
!  output parameters:
!    x    : array,length mx,giving the value of s(u) at the different
!           points. x(idim*(i-1)+j) will contain the j-th coordinate
!           of the i-th point on the curve.
!    ier  : error flag
!      ier = 0 : normal return
!      ier =10 : invalid input data (see restrictions)
!
!  restrictions:
!    m >= 1
!    mx >= m*idim
!    t(k+1) <= u(i) <= u(i+1) <= t(n-k) , i=1,2,...,m-1.
!
!  other subroutines required: fpbspl.
!
!  references :
!    de boor c : on calculating with b-splines, j. approximation theory
!                6 (1972) 50-62.
!    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
!                applics 10 (1972) 134-149.
!    dierckx p. : curve and surface fitting with splines, monographs on
!                 numerical analysis, oxford university press, 1993.
!
!  author :
!    p.dierckx
!    dept. computer science, k.u.leuven
!    celestijnenlaan 200a, b-3001 heverlee, belgium.
!    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  latest update : march 1987
!
!  ..scalar arguments..
      integer idim,n,nc,k,m,mx,ier
!  ..array arguments..
      real(dp) t(n),c(nc),u(m),x(mx)
!  ..local scalars..
      integer i,j,jj,j1,k1,l,ll,l1,mm,nk1
      real(dp) arg,sp,tb,te
!  ..local array..
      real(dp) h(6)
!  ..
!  before starting computations a data check is made. if the input data
!  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(m-1) 100,30,10
  10  do 20 i=2,m
        if(u(i).lt.u(i-1)) go to 100
  20  continue
  30  if(mx.lt.(m*idim)) go to 100
      ier = 0
!  fetch tb and te, the boundaries of the approximation interval.
      k1 = k+1
      nk1 = n-k1
      tb = t(k1)
      te = t(nk1+1)
      l = k1
      l1 = l+1
!  main loop for the different points.
      mm = 0
      do 80 i=1,m
!  fetch a new u-value arg.
        arg = u(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
!  search for knot interval t(l) <= arg < t(l+1)
  40    if(arg.lt.t(l1) .or. l.eq.nk1) go to 50
        l = l1
        l1 = l+1
        go to 40
!  evaluate the non-zero b-splines at arg.
  50    call fpbspl(t,n,k,arg,l,h)
!  find the value of s(u) at u=arg.
        ll = l-k1
        do 70 j1=1,idim
          jj = ll
          sp = 0.
          do 60 j=1,k1
            jj = jj+1
            sp = sp+c(jj)*h(j)
  60      continue
          mm = mm+1
          x(mm) = sp
          ll = ll+n
  70    continue
  80  continue
 100  return
      end subroutine curev
      subroutine fpbspl(t,n,k,x,l,h)
!  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
!  degree k at t(l) <= x < t(l+1) using the stable recurrence
!  relation of de boor and cox.
!  ..
!  ..scalar arguments..
      real(dp) x
      integer n,k,l
!  ..array arguments..
      real(dp) t(n),h(6)
!  ..local scalars..
      real(dp) f,one
      integer i,j,li,lj
!  ..local arrays..
      real(dp) hh(5)
!  ..
      one = 0.1e+01
      h(1) = one
      do 20 j=1,k
        do 10 i=1,j
          hh(i) = h(i)
  10    continue
        h(1) = 0.
        do 20 i=1,j
          li = l+i
          lj = li-j
          f = hh(i)/(t(li)-t(lj))
          h(i) = h(i)+f*(t(li)-x)
          h(i+1) = f*(x-t(lj))
  20  continue
      return
      end subroutine fpbspl
      subroutine curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp, &
        wrk,lwrk,iwrk,ier)
!  given the set of data points (x(i),y(i)) and the set of positive
!  numbers w(i),i=1,2,...,m,subroutine curfit determines a smooth spline
!  approximation of degree k on the interval xb <= x <= xe.
!  if iopt=-1 curfit calculates the weighted least-squares spline
!  according to a given set of knots.
!  if iopt>=0 the number of knots of the spline s(x) and the position
!  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
!  ness of s(x) is then achieved by minimalizing the discontinuity
!  jumps of the k-th derivative of s(x) at the knots t(j),j=k+2,k+3,...,
!  n-k-1. the amount of smoothness is determined by the condition that
!  f(p)=sum((w(i)*(y(i)-s(x(i))))**2) be <= s, with s a given non-
!  negative constant, called the smoothing factor.
!  the fit s(x) is given in the b-spline representation (b-spline coef-
!  ficients c(j),j=1,2,...,n-k-1) and can be evaluated by means of
!  subroutine splev.
!
!  calling sequence:
!     call curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,
!    * lwrk,iwrk,ier)
!
!  parameters:
!   iopt  : integer flag. on entry iopt must specify whether a weighted
!           least-squares spline (iopt=-1) or a smoothing spline (iopt=
!           0 or 1) must be determined. if iopt=0 the routine will start
!           with an initial set of knots t(i)=xb, t(i+k+1)=xe, i=1,2,...
!           k+1. if iopt=1 the routine will continue with the knots
!           found at the last call of the routine.
!           attention: a call with iopt=1 must always be immediately
!           preceded by another call with iopt=1 or iopt=0.
!           unchanged on exit.
!   m     : integer. on entry m must specify the number of data points.
!           m > k. unchanged on exit.
!   x     : real(dp) array of dimension at least (m). before entry, x(i)
!           must be set to the i-th value of the independent variable x,
!           for i=1,2,...,m. these values must be supplied in strictly
!           ascending order. unchanged on exit.
!   y     : real(dp) array of dimension at least (m). before entry, y(i)
!           must be set to the i-th value of the dependent variable y,
!           for i=1,2,...,m. unchanged on exit.
!   w     : real(dp) array of dimension at least (m). before entry, w(i)
!           must be set to the i-th value in the set of weights. the
!           w(i) must be strictly positive. unchanged on exit.
!           see also further comments.
!   xb,xe : real(dp) values. on entry xb and xe must specify the boundaries
!           of the approximation interval. xb<=x(1), xe>=x(m).
!           unchanged on exit.
!   k     : integer. on entry k must specify the degree of the spline.
!           1<=k<=5. it is recommended to use cubic splines (k=3).
!           the user is strongly dissuaded from choosing k even,together
!           with a small s-value. unchanged on exit.
!   s     : real(dp).on entry (in case iopt>=0) s must specify the smoothing
!           factor. s >=0. unchanged on exit.
!           for advice on the choice of s see further comments.
!   nest  : integer. on entry nest must contain an over-estimate of the
!           total number of knots of the spline returned, to indicate
!           the storage space available to the routine. nest >=2*k+2.
!           in most practical situation nest=m/2 will be sufficient.
!           always large enough is  nest=m+k+1, the number of knots
!           needed for interpolation (s=0). unchanged on exit.
!   n     : integer.
!           unless ier =10 (in case iopt >=0), n will contain the
!           total number of knots of the spline approximation returned.
!           if the computation mode iopt=1 is used this value of n
!           should be left unchanged between subsequent calls.
!           in case iopt=-1, the value of n must be specified on entry.
!   t     : real(dp) array of dimension at least (nest).
!           on succesful exit, this array will contain the knots of the
!           spline,i.e. the position of the interior knots t(k+2),t(k+3)
!           ...,t(n-k-1) as well as the position of the additional knots
!           t(1)=t(2)=...=t(k+1)=xb and t(n-k)=...=t(n)=xe needed for
!           the b-spline representation.
!           if the computation mode iopt=1 is used, the values of t(1),
!           t(2),...,t(n) should be left unchanged between subsequent
!           calls. if the computation mode iopt=-1 is used, the values
!           t(k+2),...,t(n-k-1) must be supplied by the user, before
!           entry. see also the restrictions (ier=10).
!   c     : real(dp) array of dimension at least (nest).
!           on succesful exit, this array will contain the coefficients
!           c(1),c(2),..,c(n-k-1) in the b-spline representation of s(x)
!   fp    : real(dp). unless ier=10, fp contains the weighted sum of
!           squared residuals of the spline approximation returned.
!   wrk   : real(dp) array of dimension at least (m*(k+1)+nest*(7+3*k)).
!           used as working space. if the computation mode iopt=1 is
!           used, the values wrk(1),...,wrk(n) should be left unchanged
!           between subsequent calls.
!   lwrk  : integer. on entry,lwrk must specify the actual dimension of
!           the array wrk as declared in the calling (sub)program.lwrk
!           must not be too small (see wrk). unchanged on exit.
!   iwrk  : integer array of dimension at least (nest).
!           used as working space. if the computation mode iopt=1 is
!           used,the values iwrk(1),...,iwrk(n) should be left unchanged
!           between subsequent calls.
!   ier   : integer. unless the routine detects an error, ier contains a
!           non-positive value on exit, i.e.
!    ier=0  : normal return. the spline returned has a residual sum of
!             squares fp such that abs(fp-s)/s <= tol with tol a relat-
!             ive tolerance set to 0.001 by the program.
!    ier=-1 : normal return. the spline returned is an interpolating
!             spline (fp=0).
!    ier=-2 : normal return. the spline returned is the weighted least-
!             squares polynomial of degree k. in this extreme case fp
!             gives the upper bound fp0 for the smoothing factor s.
!    ier=1  : error. the required storage space exceeds the available
!             storage space, as specified by the parameter nest.
!             probably causes : nest too small. if nest is already
!             large (say nest > m/2), it may also indicate that s is
!             too small
!             the approximation returned is the weighted least-squares
!             spline according to the knots t(1),t(2),...,t(n). (n=nest)
!             the parameter fp gives the corresponding weighted sum of
!             squared residuals (fp>s).
!    ier=2  : error. a theoretically impossible result was found during
!             the iteration proces for finding a smoothing spline with
!             fp = s. probably causes : s too small.
!             there is an approximation returned but the corresponding
!             weighted sum of squared residuals does not satisfy the
!             condition abs(fp-s)/s < tol.
!    ier=3  : error. the maximal number of iterations maxit (set to 20
!             by the program) allowed for finding a smoothing spline
!             with fp=s has been reached. probably causes : s too small
!             there is an approximation returned but the corresponding
!             weighted sum of squared residuals does not satisfy the
!             condition abs(fp-s)/s < tol.
!    ier=10 : error. on entry, the input data are controlled on validity
!             the following restrictions must be satisfied.
!             -1<=iopt<=1, 1<=k<=5, m>k, nest>2*k+2, w(i)>0,i=1,2,...,m
!             xb<=x(1)<x(2)<...<x(m)<=xe, lwrk>=(k+1)*m+nest*(7+3*k)
!             if iopt=-1: 2*k+2<=n<=min(nest,m+k+1)
!                         xb<t(k+2)<t(k+3)<...<t(n-k-1)<xe
!                       the schoenberg-whitney conditions, i.e. there
!                       must be a subset of data points xx(j) such that
!                         t(j) < xx(j) < t(j+k+1), j=1,2,...,n-k-1
!             if iopt>=0: s>=0
!                         if s=0 : nest >= m+k+1
!             if one of these conditions is found to be violated,control
!             is immediately repassed to the calling program. in that
!             case there is no approximation returned.
!
!  further comments:
!   by means of the parameter s, the user can control the tradeoff
!   between closeness of fit and smoothness of fit of the approximation.
!   if s is too large, the spline will be too smooth and signal will be
!   lost ; if s is too small the spline will pick up too much noise. in
!   the extreme cases the program will return an interpolating spline if
!   s=0 and the weighted least-squares polynomial of degree k if s is
!   very large. between these extremes, a properly chosen s will result
!   in a good compromise between closeness of fit and smoothness of fit.
!   to decide whether an approximation, corresponding to a certain s is
!   satisfactory the user is highly recommended to inspect the fits
!   graphically.
!   recommended values for s depend on the weights w(i). if these are
!   taken as 1/d(i) with d(i) an estimate of the standard deviation of
!   y(i), a good s-value should be found in the range (m-sqrt(2*m),m+
!   sqrt(2*m)). if nothing is known about the statistical error in y(i)
!   each w(i) can be set equal to one and s determined by trial and
!   error, taking account of the comments above. the best is then to
!   start with a very large value of s ( to determine the least-squares
!   polynomial and the corresponding upper bound fp0 for s) and then to
!   progressively decrease the value of s ( say by a factor 10 in the
!   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
!   approximation shows more detail) to obtain closer fits.
!   to economize the search for a good s-value the program provides with
!   different modes of computation. at the first call of the routine, or
!   whenever he wants to restart with the initial set of knots the user
!   must set iopt=0.
!   if iopt=1 the program will continue with the set of knots found at
!   the last call of the routine. this will save a lot of computation
!   time if curfit is called repeatedly for different values of s.
!   the number of knots of the spline returned and their location will
!   depend on the value of s and on the complexity of the shape of the
!   function underlying the data. but, if the computation mode iopt=1
!   is used, the knots returned may also depend on the s-values at
!   previous calls (if these were smaller). therefore, if after a number
!   of trials with different s-values and iopt=1, the user can finally
!   accept a fit as satisfactory, it may be worthwhile for him to call
!   curfit once more with the selected value for s but now with iopt=0.
!   indeed, curfit may then return an approximation of the same quality
!   of fit but with fewer knots and therefore better if data reduction
!   is also an important objective for the user.
!
!  other subroutines required:
!    fpback,fpbspl,fpchec,fpcurf,fpdisc,fpgivs,fpknot,fprati,fprota
!
!  references:
!   dierckx p. : an algorithm for smoothing, differentiation and integ-
!                ration of experimental data using spline functions,
!                j.comp.appl.maths 1 (1975) 165-184.
!   dierckx p. : a fast algorithm for smoothing data on a rectangular
!                grid while using spline functions, siam j.numer.anal.
!                19 (1982) 1286-1304.
!   dierckx p. : an improved algorithm for curve fitting with spline
!                functions, report tw54, dept. computer science,k.u.
!                leuven, 1981.
!   dierckx p. : curve and surface fitting with splines, monographs on
!                numerical analysis, oxford university press, 1993.
!
!  author:
!    p.dierckx
!    dept. computer science, k.u. leuven
!    celestijnenlaan 200a, b-3001 heverlee, belgium.
!    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  creation date : may 1979
!  latest update : march 1987
!
!  ..
!  ..scalar arguments..
      real(dp) xb,xe,s,fp
      integer iopt,m,k,nest,n,lwrk,ier
!  ..array arguments..
      real(dp) x(m),y(m),w(m),t(nest),c(nest),wrk(lwrk)
      integer iwrk(nest)
!  ..local scalars..
      real(dp) tol
      integer i,ia,ib,ifp,ig,iq,iz,j,k1,k2,lwest,maxit,nmin
!  ..
!  we set up the parameters tol and maxit
      maxit = 20
      tol = 0.1e-02
!  before starting computations a data check is made. if the input data
!  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(k.le.0 .or. k.gt.5) go to 50
      k1 = k+1
      k2 = k1+1
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 50
      nmin = 2*k1
      if(m.lt.k1 .or. nest.lt.nmin) go to 50
      lwest = m*k1+nest*(7+3*k)
      if(lwrk.lt.lwest) go to 50
      if(xb.gt.x(1) .or. xe.lt.x(m) .or. w(1).le.0.) go to 50
      do 10 i=2,m
         if(x(i-1).ge.x(i) .or. w(i).le.0.) go to 50
  10  continue
      if(iopt.ge.0) go to 30
      if(n.lt.nmin .or. n.gt.nest) go to 50
      j = n
      do 20 i=1,k1
         t(i) = xb
         t(j) = xe
         j = j-1
  20  continue
      call fpchec(x,m,t,n,k,ier)
      if(ier) 50,40,50
  30  if(s.lt.0.) go to 50
      if(s.eq.0. .and. nest.lt.(m+k1)) go to 50
      ier = 0
! we partition the working space and determine the spline approximation.
  40  ifp = 1
      iz = ifp+nest
      ia = iz+nest
      ib = ia+nest*k1
      ig = ib+nest*k2
      iq = ig+nest*k2
      call fpcurf(iopt,x,y,w,m,xb,xe,k,s,nest,tol,maxit,k1,k2,n,t,c,fp, &
       wrk(ifp),wrk(iz),wrk(ia),wrk(ib),wrk(ig),wrk(iq),iwrk,ier)
  50  return
      end subroutine curfit
      subroutine fpback(a,z,n,k,c,nest)
!  subroutine fpback calculates the solution of the system of
!  equations a*c = z with a a n x n upper triangular matrix
!  of bandwidth k.
!  ..
!  ..scalar arguments..
      integer n,k,nest
!  ..array arguments..
      real(dp) a(nest,k),z(n),c(n)
!  ..local scalars..
      real(dp) store
      integer i,i1,j,k1,l,m
!  ..
      k1 = k-1
      c(n) = z(n)/a(n,1)
      i = n-1
      if(i.eq.0) go to 30
      do 20 j=2,n
        store = z(i)
        i1 = k1
        if(j.le.k1) i1 = j-1
        m = i
        do 10 l=1,i1
          m = m+1
          store = store-c(m)*a(i,l+1)
  10    continue
        c(i) = store/a(i,1)
        i = i-1
  20  continue
  30  return
      end subroutine fpback
      subroutine fpchec(x,m,t,n,k,ier)
!  subroutine fpchec verifies the number and the position of the knots
!  t(j),j=1,2,...,n of a spline of degree k, in relation to the number
!  and the position of the data points x(i),i=1,2,...,m. if all of the
!  following conditions are fulfilled, the error parameter ier is set
!  to zero. if one of the conditions is violated ier is set to ten.
!      1) k+1 <= n-k-1 <= m
!      2) t(1) <= t(2) <= ... <= t(k+1)
!         t(n-k) <= t(n-k+1) <= ... <= t(n)
!      3) t(k+1) < t(k+2) < ... < t(n-k)
!      4) t(k+1) <= x(i) <= t(n-k)
!      5) the conditions specified by schoenberg and whitney must hold
!         for at least one subset of data points, i.e. there must be a
!         subset of data points y(j) such that
!             t(j) < y(j) < t(j+k+1), j=1,2,...,n-k-1
!  ..
!  ..scalar arguments..
      integer m,n,k,ier
!  ..array arguments..
      real(dp) x(m),t(n)
!  ..local scalars..
      integer i,j,k1,k2,l,nk1,nk2,nk3
      real(dp) tj,tl
!  ..
      k1 = k+1
      k2 = k1+1
      nk1 = n-k1
      nk2 = nk1+1
      ier = 10
!  check condition no 1
      if(nk1.lt.k1 .or. nk1.gt.m) go to 80
!  check condition no 2
      j = n
      do 20 i=1,k
        if(t(i).gt.t(i+1)) go to 80
        if(t(j).lt.t(j-1)) go to 80
        j = j-1
  20  continue
!  check condition no 3
      do 30 i=k2,nk2
        if(t(i).le.t(i-1)) go to 80
  30  continue
!  check condition no 4
      if(x(1).lt.t(k1) .or. x(m).gt.t(nk2)) go to 80
!  check condition no 5
      if(x(1).ge.t(k2) .or. x(m).le.t(nk1)) go to 80
      i = 1
      l = k2
      nk3 = nk1-1
      if(nk3.lt.2) go to 70
      do 60 j=2,nk3
        tj = t(j)
        l = l+1
        tl = t(l)
  40    i = i+1
        if(i.ge.m) go to 80
        if(x(i).le.tj) go to 40
        if(x(i).ge.tl) go to 80
  60  continue
  70  ier = 0
  80  return
      end subroutine fpchec
      subroutine fpcurf(iopt,x,y,w,m,xb,xe,k,s,nest,tol,maxit,k1,k2, &
       n,t,c,fp,fpint,z,a,b,g,q,nrdata,ier)
!  ..
!  ..scalar arguments..
      real(dp) xb,xe,s,tol,fp
      integer iopt,m,k,nest,maxit,k1,k2,n,ier
!  ..array arguments..
      real(dp) x(m),y(m),w(m),t(nest),c(nest),fpint(nest), &
     & z(nest),a(nest,k1),b(nest,k2),g(nest,k2),q(m,k1)
      integer nrdata(nest)
!  ..local scalars..
      real(dp) acc,con1,con4,con9,cos,half,fpart,fpms,fpold,fp0,f1,f2,f3, &
       one,p,pinv,piv,p1,p2,p3,rn,sin,store,term,wi,xi,yi
      integer i,ich1,ich3,it,iter,i1,i2,i3,j,k3,l,l0, &
       mk1,new,nk1,nmax,nmin,nplus,npl1,nrint,n8
!  ..local arrays..
      real(dp) h(7)
!  ..function references
! Edit DGP
      real(dp) abs
!      real(dp) abs,fprati
      integer max0,min0
!  ..subroutine references..
!    fpback,fpbspl,fpgivs,fpdisc,fpknot,fprota
!  ..
!  set constants
      one = 0.1e+01
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
      half = 0.5e0
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  part 1: determination of the number of knots and their position     c
!  **************************************************************      c
!  given a set of knots we compute the least-squares spline sinf(x),   c
!  and the corresponding sum of squared residuals fp=f(p=inf).         c
!  if iopt=-1 sinf(x) is the requested approximation.                  c
!  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
!    if fp <=s we will continue with the current set of knots.         c
!    if fp > s we will increase the number of knots and compute the    c
!       corresponding least-squares spline until finally fp<=s.        c
!    the initial choice of knots depends on the value of s and iopt.   c
!    if s=0 we have spline interpolation; in that case the number of   c
!    knots equals nmax = m+k+1.                                        c
!    if s > 0 and                                                      c
!      iopt=0 we first compute the least-squares polynomial of         c
!      degree k; n = nmin = 2*k+2                                      c
!      iopt=1 we start with the set of knots found at the last         c
!      call of the routine, except for the case that s > fp0; then     c
!      we compute directly the least-squares polynomial of degree k.   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  determine nmin, the number of knots for polynomial approximation.
      nmin = 2*k1
      if(iopt.lt.0) go to 60
!  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
!  determine nmax, the number of knots for spline interpolation.
      nmax = m+k1
      if(s.gt.0.) go to 45
!  if s=0, s(x) is an interpolating spline.
!  test whether the required storage space exceeds the available one.
      n = nmax
      if(nmax.gt.nest) go to 420
!  find the position of the interior knots in case of interpolation.
  10  mk1 = m-k1
      if(mk1.eq.0) go to 60
      k3 = k/2
      i = k2
      j = k3+2
      if(k3*2.eq.k) go to 30
      do 20 l=1,mk1
        t(i) = x(j)
        i = i+1
        j = j+1
  20  continue
      go to 60
  30  do 40 l=1,mk1
        t(i) = (x(j)+x(j-1))*half
        i = i+1
        j = j+1
  40  continue
      go to 60
!  if s>0 our initial choice of knots depends on the value of iopt.
!  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
!  polynomial of degree k which is a spline without interior knots.
!  if iopt=1 and fp0>s we start computing the least squares spline
!  according to the set of knots found at the last call of the routine.
  45  if(iopt.eq.0) go to 50
      if(n.eq.nmin) go to 50
      fp0 = fpint(n)
      fpold = fpint(n-1)
      nplus = nrdata(n)
      if(fp0.gt.s) go to 60
  50  n = nmin
      fpold = 0.
      nplus = 0
      nrdata(1) = m-2
!  main loop for the different sets of knots. m is a save upper bound
!  for the number of trials.
  60  do 200 iter = 1,m
        if(n.eq.nmin) ier = -2
!  find nrint, tne number of knot intervals.
        nrint = n-nmin+1
!  find the position of the additional knots which are needed for
!  the b-spline representation of s(x).
        nk1 = n-k1
        i = n
        do 70 j=1,k1
          t(j) = xb
          t(i) = xe
          i = i-1
  70    continue
!  compute the b-spline coefficients of the least-squares spline
!  sinf(x). the observation matrix a is built up row by row and
!  reduced to upper triangular form by givens transformations.
!  at the same time fp=f(p=inf) is computed.
        fp = 0.
!  initialize the observation matrix a.
        do 80 i=1,nk1
          z(i) = 0.
          do 80 j=1,k1
            a(i,j) = 0.
  80    continue
        l = k1
        do 130 it=1,m
!  fetch the current data point x(it),y(it).
          xi = x(it)
          wi = w(it)
          yi = y(it)*wi
!  search for knot interval t(l) <= xi < t(l+1).
  85      if(xi.lt.t(l+1) .or. l.eq.nk1) go to 90
          l = l+1
          go to 85
!  evaluate the (k+1) non-zero b-splines at xi and store them in q.
  90      call fpbspl(t,n,k,xi,l,h)
          do 95 i=1,k1
            q(it,i) = h(i)
            h(i) = h(i)*wi
  95      continue
!  rotate the new row of the observation matrix into triangle.
          j = l-k1
          do 110 i=1,k1
            j = j+1
            piv = h(i)
            if(piv.eq.0.) go to 110
!  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(j,1),cos,sin)
!  transformations to right hand side.
            call fprota(cos,sin,yi,z(j))
            if(i.eq.k1) go to 120
            i2 = 1
            i3 = i+1
            do 100 i1 = i3,k1
              i2 = i2+1
!  transformations to left hand side.
              call fprota(cos,sin,h(i1),a(j,i2))
 100        continue
 110      continue
!  add contribution of this row to the sum of squares of residual
!  right hand sides.
 120      fp = fp+yi**2
 130    continue
        if(ier.eq.(-2)) fp0 = fp
        fpint(n) = fp0
        fpint(n-1) = fpold
        nrdata(n) = nplus
!  backward substitution to obtain the b-spline coefficients.
        call fpback(a,z,nk1,k1,c,nest)
!  test whether the approximation sinf(x) is an acceptable solution.
        if(iopt.lt.0) go to 440
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
!  if f(p=inf) < s accept the choice of knots.
        if(fpms.lt.0.) go to 250
!  if n = nmax, sinf(x) is an interpolating spline.
        if(n.eq.nmax) go to 430
!  increase the number of knots.
!  if n=nest we cannot increase the number of knots because of
!  the storage capacity limitation.
        if(n.eq.nest) go to 420
!  determine the number of knots nplus we are going to add.
        if(ier.eq.0) go to 140
        nplus = 1
        ier = 0
        go to 150
 140    npl1 = nplus*2
        rn = nplus
        if(fpold-fp.gt.acc) npl1 = rn*fpms/(fpold-fp)
        nplus = min0(nplus*2,max0(npl1,nplus/2,1))
 150    fpold = fp
!  compute the sum((w(i)*(y(i)-s(x(i))))**2) for each knot interval
!  t(j+k) <= x(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = 0.
        i = 1
        l = k2
        new = 0
        do 180 it=1,m
          if(x(it).lt.t(l) .or. l.gt.nk1) go to 160
          new = 1
          l = l+1
 160      term = 0.
          l0 = l-k2
          do 170 j=1,k1
            l0 = l0+1
            term = term+c(l0)*q(it,j)
 170      continue
          term = (w(it)*(term-y(it)))**2
          fpart = fpart+term
          if(new.eq.0) go to 180
          store = term*half
          fpint(i) = fpart-store
          i = i+1
          fpart = store
          new = 0
 180    continue
        fpint(nrint) = fpart
        do 190 l=1,nplus
!  add a new knot.
          call fpknot(x,m,t,n,fpint,nrdata,nrint,nest,1)
!  if n=nmax we locate the knots as for interpolation.
          if(n.eq.nmax) go to 10
!  test whether we cannot further increase the number of knots.
          if(n.eq.nest) go to 200
 190    continue
!  restart the computations with the new set of knots.
 200  continue
!  test whether the least-squares kth degree polynomial is a solution
!  of our approximation problem.
 250  if(ier.eq.(-2)) go to 440
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  part 2: determination of the smoothing spline sp(x).                c
!  ***************************************************                 c
!  we have determined the number of knots and their position.          c
!  we now compute the b-spline coefficients of the smoothing spline    c
!  sp(x). the observation matrix a is extended by the rows of matrix   c
!  b expressing that the kth derivative discontinuities of sp(x) at    c
!  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
!  ponding weights of these additional rows are set to 1/p.            c
!  iteratively we then have to determine the value of p such that      c
!  f(p)=sum((w(i)*(y(i)-sp(x(i))))**2) be = s. we already know that    c
!  the least-squares kth degree polynomial corresponds to p=0, and     c
!  that the least-squares spline corresponds to p=infinity. the        c
!  iteration process which is proposed here, makes use of rational     c
!  interpolation. since f(p) is a convex and strictly decreasing       c
!  function of p, it can be approximated by a rational function        c
!  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
!  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
!  to calculate the new value of p such that r(p)=s. convergence is    c
!  guaranteed by taking f1>0 and f3<0.                                 c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  evaluate the discontinuity jump of the kth derivative of the
!  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
      call fpdisc(t,n,k2,b,nest)
!  initial value for p.
      p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = 0.
      do 255 i=1,nk1
         p = p+a(i,1)
 255  continue
      rn = nk1
      p = rn/p
      ich1 = 0
      ich3 = 0
      n8 = n-nmin
!  iteration process to find the root of f(p) = s.
      do 360 iter=1,maxit
!  the rows of matrix b with weight 1/p are rotated into the
!  triangularised observation matrix a which is stored in g.
        pinv = one/p
        do 260 i=1,nk1
          c(i) = z(i)
          g(i,k2) = 0.
          do 260 j=1,k1
            g(i,j) = a(i,j)
 260    continue
        do 300 it=1,n8
!  the row of matrix b is rotated into triangle by givens transformation
          do 270 i=1,k2
            h(i) = b(it,i)*pinv
 270      continue
          yi = 0.
          do 290 j=it,nk1
            piv = h(1)
!  calculate the parameters of the givens transformation.
            call fpgivs(piv,g(j,1),cos,sin)
!  transformations to right hand side.
            call fprota(cos,sin,yi,c(j))
            if(j.eq.nk1) go to 300
            i2 = k1
            if(j.gt.n8) i2 = nk1-j
            do 280 i=1,i2
!  transformations to left hand side.
              i1 = i+1
              call fprota(cos,sin,h(i1),g(j,i1))
              h(i) = h(i1)
 280        continue
            h(i2+1) = 0.
 290      continue
 300    continue
!  backward substitution to obtain the b-spline coefficients.
        call fpback(g,c,nk1,k2,c,nest)
!  computation of f(p).
        fp = 0.
        l = k2
        do 330 it=1,m
          if(x(it).lt.t(l) .or. l.gt.nk1) go to 310
          l = l+1
 310      l0 = l-k2
          term = 0.
          do 320 j=1,k1
            l0 = l0+1
            term = term+c(l0)*q(it,j)
 320      continue
          fp = fp+(w(it)*(term-y(it)))**2
 330    continue
!  test whether the approximation sp(x) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
!  test whether the maximal number of iterations is reached.
        if(iter.eq.maxit) go to 400
!  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 340
        if((f2-f3).gt.acc) go to 335
!  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p=p1*con9 + p2*con1
        go to 360
 335    if(f2.lt.0.) ich3=1
 340    if(ich1.ne.0) go to 350
        if((f1-f2).gt.acc) go to 345
!  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 360
        if(p.ge.p3) p = p2*con1 + p3*con9
        go to 360
 345    if(f2.gt.0.) ich1=1
!  test whether the iteration process proceeds as theoretically
!  expected.
 350    if(f2.ge.f1 .or. f2.le.f3) go to 410
!  find the new value for p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 360  continue
!  error codes and messages.
 400  ier = 3
      go to 440
 410  ier = 2
      go to 440
 420  ier = 1
      go to 440
 430  ier = -1
 440  return
      end subroutine fpcurf
      subroutine fpdisc(t,n,k2,b,nest)
!  subroutine fpdisc calculates the discontinuity jumps of the kth
!  derivative of the b-splines of degree k at the knots t(k+2)..t(n-k-1)
!  ..scalar arguments..
      integer n,k2,nest
!  ..array arguments..
      real(dp) t(n),b(nest,k2)
!  ..local scalars..
      real(dp) an,fac,prod
      integer i,ik,j,jk,k,k1,l,lj,lk,lmk,lp,nk1,nrint
!  ..local array..
      real(dp) h(12)
!  ..
      k1 = k2-1
      k = k1-1
      nk1 = n-k1
      nrint = nk1-k
      an = nrint
      fac = an/(t(nk1+1)-t(k1))
      do 40 l=k2,nk1
        lmk = l-k1
        do 10 j=1,k1
          ik = j+k1
          lj = l+j
          lk = lj-k2
          h(j) = t(l)-t(lk)
          h(ik) = t(l)-t(lj)
  10    continue
        lp = lmk
        do 30 j=1,k2
          jk = j
          prod = h(j)
          do 20 i=1,k
            jk = jk+1
            prod = prod*h(jk)*fac
  20      continue
          lk = lp+k1
          b(lmk,j) = (t(lk)-t(lp))/prod
          lp = lp+1
  30    continue
  40  continue
      return
      end subroutine fpdisc
      subroutine fpgivs(piv,ww,cos,sin)
!  subroutine fpgivs calculates the parameters of a givens
!  transformation .
!  ..
!  ..scalar arguments..
      real(dp) piv,ww,cos,sin
!  ..local scalars..
      real(dp) dd,one,store
!  ..function references..
      real(dp) abs,sqrt
!  ..
      one = 0.1e+01
      store = abs(piv)
      if(store.ge.ww) dd = store*sqrt(one+(ww/piv)**2)
      if(store.lt.ww) dd = ww*sqrt(one+(piv/ww)**2)
      cos = ww/dd
      sin = piv/dd
      ww = dd
      return
      end subroutine fpgivs
      subroutine fpknot(x,m,t,n,fpint,nrdata,nrint,nest,istart)
!  subroutine fpknot locates an additional knot for a spline of degree
!  k and adjusts the corresponding parameters,i.e.
!    t     : the position of the knots.
!    n     : the number of knots.
!    nrint : the number of knotintervals.
!    fpint : the sum of squares of residual right hand sides
!            for each knot interval.
!    nrdata: the number of data points inside each knot interval.
!  istart indicates that the smallest data point at which the new knot
!  may be added is x(istart+1)
!  ..
!  ..scalar arguments..
      integer m,n,nrint,nest,istart
!  ..array arguments..
      real(dp) x(m),t(nest),fpint(nest)
      integer nrdata(nest)
!  ..local scalars..
      real(dp) an,am,fpmax
      integer ihalf,j,jbegin,jj,jk,jpoint,k,maxbeg,maxpt, &
       next,nrx,number
!  ..
      k = (n-nrint-1)/2
!  search for knot interval t(number+k) <= x <= t(number+k+1) where
!  fpint(number) is maximal on the condition that nrdata(number)
!  not equals zero.
      fpmax = 0.
      jbegin = istart
      do 20 j=1,nrint
        jpoint = nrdata(j)
        if(fpmax.ge.fpint(j) .or. jpoint.eq.0) go to 10
        fpmax = fpint(j)
        number = j
        maxpt = jpoint
        maxbeg = jbegin
  10    jbegin = jbegin+jpoint+1
  20  continue
!  let coincide the new knot t(number+k+1) with a data point x(nrx)
!  inside the old knot interval t(number+k) <= x <= t(number+k+1).
      ihalf = maxpt/2+1
      nrx = maxbeg+ihalf
      next = number+1
      if(next.gt.nrint) go to 40
!  adjust the different parameters.
      do 30 j=next,nrint
        jj = next+nrint-j
        fpint(jj+1) = fpint(jj)
        nrdata(jj+1) = nrdata(jj)
        jk = jj+k
        t(jk+1) = t(jk)
  30  continue
  40  nrdata(number) = ihalf-1
      nrdata(next) = maxpt-ihalf
      am = maxpt
      an = nrdata(number)
      fpint(number) = fpmax*an/am
      an = nrdata(next)
      fpint(next) = fpmax*an/am
      jk = next+k
      t(jk) = x(nrx)
      n = n+1
      nrint = nrint+1
      return
      end subroutine fpknot
      real(dp) function fprati(p1,f1,p2,f2,p3,f3)
!  given three points (p1,f1),(p2,f2) and (p3,f3), function fprati
!  gives the value of p such that the rational interpolating function
!  of the form r(p) = (u*p+v)/(p+w) equals zero at p.
!  ..
!  ..scalar arguments..
      real(dp) p1,f1,p2,f2,p3,f3
!  ..local scalars..
      real(dp) h1,h2,h3,p
!  ..
      if(p3.gt.0.) go to 10
!  value of p in case p3 = infinity.
      p = (p1*(f1-f3)*f2-p2*(f2-f3)*f1)/((f1-f2)*f3)
      go to 20
!  value of p in case p3 ^= infinity.
  10  h1 = f1*(f2-f3)
      h2 = f2*(f3-f1)
      h3 = f3*(f1-f2)
      p = -(p1*p2*h3+p2*p3*h1+p3*p1*h2)/(p1*h1+p2*h2+p3*h3)
!  adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
  20  if(f2.lt.0.) go to 30
      p1 = p2
      f1 = f2
      go to 40
  30  p3 = p2
      f3 = f2
  40  fprati = p
      return
      end function fprati
      subroutine fprota(cos,sin,a,b)
!  subroutine fprota applies a givens rotation to a and b.
!  ..
!  ..scalar arguments..
      real(dp) cos,sin,a,b
! ..local scalars..
      real(dp) stor1,stor2
!  ..
      stor1 = a
      stor2 = b
      b = cos*stor2+sin*stor1
      a = cos*stor1-sin*stor2
      return
      end subroutine fprota
      subroutine fourco(t,n,c,alfa,m,ress,resc,wrk1,wrk2,ier)
!  subroutine fourco calculates the integrals
!                    /t(n-3)
!    ress(i) =      !        s(x)*sin(alfa(i)*x) dx    and
!              t(4)/
!                    /t(n-3)
!    resc(i) =      !        s(x)*cos(alfa(i)*x) dx, i=1,...,m,
!              t(4)/
!  where s(x) denotes a cubic spline which is given in its
!  b-spline representation.
!
!  calling sequence:
!     call fourco(t,n,c,alfa,m,ress,resc,wrk1,wrk2,ier)
!
!  input parameters:
!    t    : real(dp) array,length n, containing the knots of s(x).
!    n    : integer, containing the total number of knots. n>=10.
!    c    : real(dp) array,length n, containing the b-spline coefficients.
!    alfa : real(dp) array,length m, containing the parameters alfa(i).
!    m    : integer, specifying the number of integrals to be computed.
!    wrk1 : real(dp) array,length n. used as working space
!    wrk2 : real(dp) array,length n. used as working space
!
!  output parameters:
!    ress : real(dp) array,length m, containing the integrals ress(i).
!    resc : real(dp) array,length m, containing the integrals resc(i).
!    ier  : error flag:
!      ier=0 : normal return.
!      ier=10: invalid input data (see restrictions).
!
!  restrictions:
!    n >= 10
!    t(4) < t(5) < ... < t(n-4) < t(n-3).
!    t(1) <= t(2) <= t(3) <= t(4).
!    t(n-3) <= t(n-2) <= t(n-1) <= t(n).
!
!  other subroutines required: fpbfou,fpcsin
!
!  references :
!    dierckx p. : calculation of fouriercoefficients of discrete
!                 functions using cubic splines. j. computational
!                 and applied mathematics 3 (1977) 207-209.
!    dierckx p. : curve and surface fitting with splines, monographs on
!                 numerical analysis, oxford university press, 1993.
!
!  author :
!    p.dierckx
!    dept. computer science, k.u.leuven
!    celestijnenlaan 200a, b-3001 heverlee, belgium.
!    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  latest update : march 1987
!
!  ..scalar arguments..
      integer n,m,ier
!  ..array arguments..
      real(dp) t(n),c(n),wrk1(n),wrk2(n),alfa(m),ress(m),resc(m)
!  ..local scalars..
      integer i,j,n4
      real(dp) rs,rc
!  ..
      n4 = n-4
!  before starting computations a data check is made. in the input data
!  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(n.lt.10) go to 50
      j = n
      do 10 i=1,3
        if(t(i).gt.t(i+1)) go to 50
        if(t(j).lt.t(j-1)) go to 50
        j = j-1
  10  continue
      do 20 i=4,n4
        if(t(i).ge.t(i+1)) go to 50
  20  continue
      ier = 0
!  main loop for the different alfa(i).
      do 40 i=1,m
!  calculate the integrals
!    wrk1(j) = integral(nj,4(x)*sin(alfa*x))    and
!    wrk2(j) = integral(nj,4(x)*cos(alfa*x)),  j=1,2,...,n-4,
!  where nj,4(x) denotes the normalised cubic b-spline defined on the
!  knots t(j),t(j+1),...,t(j+4).
         call fpbfou(t,n,alfa(i),wrk1,wrk2)
!  calculate the integrals ress(i) and resc(i).
         rs = 0.
         rc = 0.
         do 30 j=1,n4
            rs = rs+c(j)*wrk1(j)
            rc = rc+c(j)*wrk2(j)
  30     continue
         ress(i) = rs
         resc(i) = rc
  40  continue
  50  return
      end subroutine fourco
      subroutine fpbfou(t,n,par,ress,resc)
!  subroutine fpbfou calculates the integrals
!                    /t(n-3)
!    ress(j) =      !        nj,4(x)*sin(par*x) dx    and
!              t(4)/
!                    /t(n-3)
!    resc(j) =      !        nj,4(x)*cos(par*x) dx ,  j=1,2,...n-4
!              t(4)/
!  where nj,4(x) denotes the cubic b-spline defined on the knots
!  t(j),t(j+1),...,t(j+4).
!
!  calling sequence:
!     call fpbfou(t,n,par,ress,resc)
!
!  input parameters:
!    t    : real(dp) array,length n, containing the knots.
!    n    : integer, containing the number of knots.
!    par  : real, containing the value of the parameter par.
!
!  output parameters:
!    ress  : real(dp) array,length n, containing the integrals ress(j).
!    resc  : real(dp) array,length n, containing the integrals resc(j).
!
!  restrictions:
!    n >= 10, t(4) < t(5) < ... < t(n-4) < t(n-3).
!  ..
!  ..scalar arguments..
      integer n
      real(dp) par
!  ..array arguments..
      real(dp) t(n),ress(n),resc(n)
!  ..local scalars..
      integer i,ic,ipj,is,j,jj,jp1,jp4,k,li,lj,ll,nmj,nm3,nm7
      real(dp) ak,beta,con1,con2,c1,c2,delta,eps,fac,f1,f2,f3,one,quart, &
       sign,six,s1,s2,term
!  ..local arrays..
      real(dp) co(5),si(5),hs(5),hc(5),rs(3),rc(3)
!  ..function references..
      real(dp) cos,sin,abs
!  ..
!  initialization.
      one = 0.1e+01
      six = 0.6e+01
      eps = 0.1e-07
      quart = 0.25e0
      con1 = 0.5e-01
      con2 = 0.12e+03
      nm3 = n-3
      nm7 = n-7
      if(par.ne.0.) term = six/par
      beta = par*t(4)
      co(1) = cos(beta)
      si(1) = sin(beta)
!  calculate the integrals ress(j) and resc(j), j=1,2,3 by setting up
!  a divided difference table.
      do 30 j=1,3
        jp1 = j+1
        jp4 = j+4
        beta = par*t(jp4)
        co(jp1) = cos(beta)
        si(jp1) = sin(beta)
        call fpcsin(t(4),t(jp4),par,si(1),co(1),si(jp1),co(jp1), &
        rs(j),rc(j))
        i = 5-j
        hs(i) = 0.
        hc(i) = 0.
        do 10 jj=1,j
          ipj = i+jj
          hs(ipj) = rs(jj)
          hc(ipj) = rc(jj)
  10    continue
        do 20 jj=1,3
          if(i.lt.jj) i = jj
          k = 5
          li = jp4
          do 20 ll=i,4
            lj = li-jj
            fac = t(li)-t(lj)
            hs(k) = (hs(k)-hs(k-1))/fac
            hc(k) = (hc(k)-hc(k-1))/fac
            k = k-1
            li = li-1
  20    continue
        ress(j) = hs(5)-hs(4)
        resc(j) = hc(5)-hc(4)
  30  continue
      if(nm7.lt.4) go to 160
!  calculate the integrals ress(j) and resc(j),j=4,5,...,n-7.
      do 150 j=4,nm7
        jp4 = j+4
        beta = par*t(jp4)
        co(5) = cos(beta)
        si(5) = sin(beta)
        delta = t(jp4)-t(j)
!  the way of computing ress(j) and resc(j) depends on the value of
!  beta = par*(t(j+4)-t(j)).
        beta = delta*par
        if(abs(beta).le.one) go to 60
!  if !beta! > 1 the integrals are calculated by setting up a divided
!  difference table.
        do 40 k=1,5
          hs(k) = si(k)
          hc(k) = co(k)
  40    continue
        do 50 jj=1,3
          k = 5
          li = jp4
          do 50 ll=jj,4
            lj = li-jj
            fac = par*(t(li)-t(lj))
            hs(k) = (hs(k)-hs(k-1))/fac
            hc(k) = (hc(k)-hc(k-1))/fac
            k = k-1
            li = li-1
  50    continue
        s2 = (hs(5)-hs(4))*term
        c2 = (hc(5)-hc(4))*term
        go to 130
!  if !beta! <= 1 the integrals are calculated by evaluating a series
!  expansion.
  60    f3 = 0.
        do 70 i=1,4
          ipj = i+j
          hs(i) = par*(t(ipj)-t(j))
          hc(i) = hs(i)
          f3 = f3+hs(i)
  70    continue
        f3 = f3*con1
        c1 = quart
        s1 = f3
        if(abs(f3).le.eps) go to 120
        sign = one
        fac = con2
        k = 5
        is = 0
        do 110 ic=1,20
          k = k+1
          ak = k
          fac = fac*ak
          f1 = 0.
          f3 = 0.
          do 80 i=1,4
            f1 = f1+hc(i)
            f2 = f1*hs(i)
            hc(i) = f2
            f3 = f3+f2
  80      continue
          f3 = f3*six/fac
          if(is.eq.0) go to 90
          is = 0
          s1 = s1+f3*sign
          go to 100
  90      sign = -sign
          is = 1
          c1 = c1+f3*sign
 100      if(abs(f3).le.eps) go to 120
 110    continue
 120    s2 = delta*(co(1)*s1+si(1)*c1)
        c2 = delta*(co(1)*c1-si(1)*s1)
 130    ress(j) = s2
        resc(j) = c2
        do 140 i=1,4
          co(i) = co(i+1)
          si(i) = si(i+1)
 140    continue
 150  continue
!  calculate the integrals ress(j) and resc(j),j=n-6,n-5,n-4 by setting
!  up a divided difference table.
 160  do 190 j=1,3
        nmj = nm3-j
        i = 5-j
        call fpcsin(t(nm3),t(nmj),par,si(4),co(4),si(i-1),co(i-1), &
        rs(j),rc(j))
        hs(i) = 0.
        hc(i) = 0.
        do 170 jj=1,j
          ipj = i+jj
          hc(ipj) = rc(jj)
          hs(ipj) = rs(jj)
 170    continue
        do 180 jj=1,3
          if(i.lt.jj) i = jj
          k = 5
          li = nmj
          do 180 ll=i,4
            lj = li+jj
            fac = t(lj)-t(li)
            hs(k) = (hs(k-1)-hs(k))/fac
            hc(k) = (hc(k-1)-hc(k))/fac
            k = k-1
            li = li+1
 180    continue
        ress(nmj) = hs(4)-hs(5)
        resc(nmj) = hc(4)-hc(5)
 190  continue
      return
      end subroutine fpbfou
      subroutine fpcsin(a,b,par,sia,coa,sib,cob,ress,resc)
!  fpcsin calculates the integrals ress=integral((b-x)**3*sin(par*x))
!  and resc=integral((b-x)**3*cos(par*x)) over the interval (a,b),
!  given sia=sin(par*a),coa=cos(par*a),sib=sin(par*b) and cob=cos(par*b)
!  ..
!  ..scalar arguments..
      real(dp) a,b,par,sia,coa,sib,cob,ress,resc
!  ..local scalars..
      integer i,j
      real(dp) ab,ab4,ai,alfa,beta,b2,b4,eps,fac,f1,f2,one,quart,six, &
       three,two
!  ..function references..
      real(dp) abs
!  ..
      one = 0.1e+01
      two = 0.2e+01
      three = 0.3e+01
      six = 0.6e+01
      quart = 0.25e+0
      eps = 0.1e-09
      ab = b-a
      ab4 = ab**4
      alfa = ab*par
! the way of calculating the integrals ress and resc depends on
! the value of alfa = (b-a)*par.
      if(abs(alfa).le.one) go to 100
! integration by parts.
      beta = one/alfa
      b2 = beta**2
      b4 = six*b2**2
      f1 = three*b2*(one-two*b2)
      f2 = beta*(one-six*b2)
      ress = ab4*(coa*f2+sia*f1+sib*b4)
      resc = ab4*(coa*f1-sia*f2+cob*b4)
      go to 400
! ress and resc are found by evaluating a series expansion.
 100  fac = quart
      f1 = fac
      f2 = 0.
      i = 4
      do 200 j=1,5
        i = i+1
        ai = i
        fac = fac*alfa/ai
        f2 = f2+fac
        if(abs(fac).le.eps) go to 300
        i = i+1
        ai = i
        fac = -fac*alfa/ai
        f1 = f1+fac
        if(abs(fac).le.eps) go to 300
 200  continue
 300  ress = ab4*(coa*f2+sia*f1)
      resc = ab4*(coa*f1-sia*f2)
 400  return
      end subroutine fpcsin
      subroutine splder(t,n,c,k,nu,x,y,m,wrk,ier)
!  subroutine splder evaluates in a number of points x(i),i=1,2,...,m
!  the derivative of order nu of a spline s(x) of degree k,given in
!  its b-spline representation.
!
!  calling sequence:
!     call splder(t,n,c,k,nu,x,y,m,wrk,ier)
!
!  input parameters:
!    t    : array,length n, which contains the position of the knots.
!    n    : integer, giving the total number of knots of s(x).
!    c    : array,length n, which contains the b-spline coefficients.
!    k    : integer, giving the degree of s(x).
!    nu   : integer, specifying the order of the derivative. 0<=nu<=k
!    x    : array,length m, which contains the points where the deriv-
!           ative of s(x) must be evaluated.
!    m    : integer, giving the number of points where the derivative
!           of s(x) must be evaluated
!    wrk  : real(dp) array of dimension n. used as working space.
!
!  output parameters:
!    y    : array,length m, giving the value of the derivative of s(x)
!           at the different points.
!    ier  : error flag
!      ier = 0 : normal return
!      ier =10 : invalid input data (see restrictions)
!
!  restrictions:
!    0 <= nu <= k
!    m >= 1
!    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
!
!  other subroutines required: fpbspl
!
!  references :
!    de boor c : on calculating with b-splines, j. approximation theory
!                6 (1972) 50-62.
!    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
!                applics 10 (1972) 134-149.
!   dierckx p. : curve and surface fitting with splines, monographs on
!                numerical analysis, oxford university press, 1993.
!
!  author :
!    p.dierckx
!    dept. computer science, k.u.leuven
!    celestijnenlaan 200a, b-3001 heverlee, belgium.
!    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  latest update : march 1987
!
!  ..scalar arguments..
      integer n,k,nu,m,ier
!  ..array arguments..
      real(dp) t(n),c(n),x(m),y(m),wrk(n)
!  ..local scalars..
      integer i,j,kk,k1,k2,l,ll,l1,l2,nk1,nk2,nn
      real(dp) ak,arg,fac,sp,tb,te
!  ..local arrays ..
      real(dp) h(6)
!  before starting computations a data check is made. if the input data
!  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(nu.lt.0 .or. nu.gt.k) go to 200
      if(m-1) 200,30,10
  10  do 20 i=2,m
        if(x(i).lt.x(i-1)) go to 200
  20  continue
  30  ier = 0
!  fetch tb and te, the boundaries of the approximation interval.
      k1 = k+1
      nk1 = n-k1
      tb = t(k1)
      te = t(nk1+1)
!  the derivative of order nu of a spline of degree k is a spline of
!  degree k-nu,the b-spline coefficients wrk(i) of which can be found
!  using the recurrence scheme of de boor.
      l = 1
      kk = k
      nn = n
      do 40 i=1,nk1
         wrk(i) = c(i)
  40  continue
      if(nu.eq.0) go to 100
      nk2 = nk1
      do 60 j=1,nu
         ak = kk
         nk2 = nk2-1
         l1 = l
         do 50 i=1,nk2
            l1 = l1+1
            l2 = l1+kk
            fac = t(l2)-t(l1)
            if(fac.le.0.) go to 50
            wrk(i) = ak*(wrk(i+1)-wrk(i))/fac
  50     continue
         l = l+1
         kk = kk-1
  60  continue
      if(kk.ne.0) go to 100
!  if nu=k the derivative is a piecewise constant function
      j = 1
      do 90 i=1,m
         arg = x(i)
  70     if(arg.lt.t(l+1) .or. l.eq.nk1) go to 80
         l = l+1
         j = j+1
         go to 70
  80     y(i) = wrk(j)
  90  continue
      go to 200
 100  l = k1
      l1 = l+1
      k2 = k1-nu
!  main loop for the different points.
      do 180 i=1,m
!  fetch a new x-value arg.
        arg = x(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
!  search for knot interval t(l) <= arg < t(l+1)
 140    if(arg.lt.t(l1) .or. l.eq.nk1) go to 150
        l = l1
        l1 = l+1
        go to 140
!  evaluate the non-zero b-splines of degree k-nu at arg.
 150    call fpbspl(t,n,kk,arg,l,h)
!  find the value of the derivative at x=arg.
        sp = 0.
        ll = l-k1
        do 160 j=1,k2
          ll = ll+1
          sp = sp+wrk(ll)*h(j)
 160    continue
        y(i) = sp
 180  continue
 200  return
      end subroutine splder
      subroutine splev(t,n,c,k,x,y,m,ier)
!  subroutine splev evaluates in a number of points x(i),i=1,2,...,m
!  a spline s(x) of degree k, given in its b-spline representation.
!
!  calling sequence:
!     call splev(t,n,c,k,x,y,m,ier)
!
!  input parameters:
!    t    : array,length n, which contains the position of the knots.
!    n    : integer, giving the total number of knots of s(x).
!    c    : array,length n, which contains the b-spline coefficients.
!    k    : integer, giving the degree of s(x).
!    x    : array,length m, which contains the points where s(x) must
!           be evaluated.
!    m    : integer, giving the number of points where s(x) must be
!           evaluated.
!
!  output parameter:
!    y    : array,length m, giving the value of s(x) at the different
!           points.
!    ier  : error flag
!      ier = 0 : normal return
!      ier =10 : invalid input data (see restrictions)
!
!  restrictions:
!    m >= 1
!    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
!
!  other subroutines required: fpbspl.
!
!  references :
!    de boor c  : on calculating with b-splines, j. approximation theory
!                 6 (1972) 50-62.
!    cox m.g.   : the numerical evaluation of b-splines, j. inst. maths
!                 applics 10 (1972) 134-149.
!    dierckx p. : curve and surface fitting with splines, monographs on
!                 numerical analysis, oxford university press, 1993.
!
!  author :
!    p.dierckx
!    dept. computer science, k.u.leuven
!    celestijnenlaan 200a, b-3001 heverlee, belgium.
!    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  latest update : march 1987
!
!  ..scalar arguments..
      integer n,k,m,ier
!  ..array arguments..
      real(dp) t(n),c(n),x(m),y(m)
!  ..local scalars..
      integer i,j,k1,l,ll,l1,nk1
      real(dp) arg,sp,tb,te
!  ..local array..
      real(dp) h(6)
!  ..
!  before starting computations a data check is made. if the input data
!  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(m-1) 100,30,10
  10  do 20 i=2,m
        if(x(i).lt.x(i-1)) go to 100
  20  continue
  30  ier = 0
!  fetch tb and te, the boundaries of the approximation interval.
      k1 = k+1
      nk1 = n-k1
      tb = t(k1)
      te = t(nk1+1)
      l = k1
      l1 = l+1
!  main loop for the different points.
      do 80 i=1,m
!  fetch a new x-value arg.
        arg = x(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
!  search for knot interval t(l) <= arg < t(l+1)
  40    if(arg.lt.t(l1) .or. l.eq.nk1) go to 50
        l = l1
        l1 = l+1
        go to 40
!  evaluate the non-zero b-splines at arg.
  50    call fpbspl(t,n,k,arg,l,h)
!  find the value of s(x) at x=arg.
        sp = 0.
        ll = l-k1
        do 60 j=1,k1
          ll = ll+1
          sp = sp+c(ll)*h(j)
  60    continue
        y(i) = sp
  80  continue
 100  return
      end subroutine splev
      subroutine fpintb(t,n,bint,nk1,x,y)
!  subroutine fpintb calculates integrals of the normalized b-splines
!  nj,k+1(x) of degree k, defined on the set of knots t(j),j=1,2,...n.
!  it makes use of the formulae of gaffney for the calculation of
!  indefinite integrals of b-splines.
!
!  calling sequence:
!     call fpintb(t,n,bint,nk1,x,y)
!
!  input parameters:
!    t    : real(dp) array,length n, containing the position of the knots.
!    n    : integer value, giving the number of knots.
!    nk1  : integer value, giving the number of b-splines of degree k,
!           defined on the set of knots ,i.e. nk1 = n-k-1.
!    x,y  : real(dp) values, containing the end points of the integration
!           interval.
!  output parameter:
!    bint : array,length nk1, containing the integrals of the b-splines.
!  ..
!  ..scalars arguments..
      integer n,nk1
      real(dp) x,y
!  ..array arguments..
      real(dp) t(n),bint(nk1)
!  ..local scalars..
      integer i,ia,ib,it,j,j1,k,k1,l,li,lj,lk,l0,min
      real(dp) a,ak,arg,b,f,one
!  ..local arrays..
      real(dp) aint(6),h(6),h1(6)
!  initialization.
      one = 0.1e+01
      k1 = n-nk1
      ak = k1
      k = k1-1
      do 10 i=1,nk1
        bint(i) = 0.
  10  continue
!  the integration limits are arranged in increasing order.
      a = x
      b = y
      min = 0
      if(a-b) 30,160,20
  20  a = y
      b = x
      min = 1
  30  if(a.lt.t(k1)) a = t(k1)
      if(b.gt.t(nk1+1)) b = t(nk1+1)
!  using the expression of gaffney for the indefinite integral of a
!  b-spline we find that
!  bint(j) = (t(j+k+1)-t(j))*(res(j,b)-res(j,a))/(k+1)
!    where for t(l) <= x < t(l+1)
!    res(j,x) = 0, j=1,2,...,l-k-1
!             = 1, j=l+1,l+2,...,nk1
!             = aint(j+k-l+1), j=l-k,l-k+1,...,l
!               = sumi((x-t(j+i))*nj+i,k+1-i(x)/(t(j+k+1)-t(j+i)))
!                 i=0,1,...,k
      l = k1
      l0 = l+1
!  set arg = a.
      arg = a
      do 90 it=1,2
!  search for the knot interval t(l) <= arg < t(l+1).
  40    if(arg.lt.t(l0) .or. l.eq.nk1) go to 50
        l = l0
        l0 = l+1
        go to 40
!  calculation of aint(j), j=1,2,...,k+1.
!  initialization.
  50    do 55 j=1,k1
          aint(j) = 0.
  55    continue
        aint(1) = (arg-t(l))/(t(l+1)-t(l))
        h1(1) = one
        do 70 j=1,k
!  evaluation of the non-zero b-splines of degree j at arg,i.e.
!    h(i+1) = nl-j+i,j(arg), i=0,1,...,j.
          h(1) = 0.
          do 60 i=1,j
            li = l+i
            lj = li-j
            f = h1(i)/(t(li)-t(lj))
            h(i) = h(i)+f*(t(li)-arg)
            h(i+1) = f*(arg-t(lj))
  60      continue
!  updating of the integrals aint.
          j1 = j+1
          do 70 i=1,j1
            li = l+i
            lj = li-j1
            aint(i) = aint(i)+h(i)*(arg-t(lj))/(t(li)-t(lj))
            h1(i) = h(i)
  70    continue
        if(it.eq.2) go to 100
!  updating of the integrals bint
        lk = l-k
        ia = lk
        do 80 i=1,k1
          bint(lk) = -aint(i)
          lk = lk+1
  80    continue
!  set arg = b.
        arg = b
  90  continue
!  updating of the integrals bint.
 100  lk = l-k
      ib = lk-1
      do 110 i=1,k1
        bint(lk) = bint(lk)+aint(i)
        lk = lk+1
 110  continue
      if(ib.lt.ia) go to 130
      do 120 i=ia,ib
        bint(i) = bint(i)+one
 120  continue
!  the scaling factors are taken into account.
 130  f = one/ak
      do 140 i=1,nk1
        j = i+k1
        bint(i) = bint(i)*(t(j)-t(i))*f
 140  continue
!  the order of the integration limits is taken into account.
      if(min.eq.0) go to 160
      do 150 i=1,nk1
        bint(i) = -bint(i)
 150  continue
 160  return
      end subroutine fpintb
      real(dp) function splint(t,n,c,k,a,b,wrk)
!  function splint calculates the integral of a spline function s(x)
!  of degree k, which is given in its normalized b-spline representation
!
!  calling sequence:
!     aint = splint(t,n,c,k,a,b,wrk)
!
!  input parameters:
!    t    : array,length n,which contains the position of the knots
!           of s(x).
!    n    : integer, giving the total number of knots of s(x).
!    c    : array,length n, containing the b-spline coefficients.
!    k    : integer, giving the degree of s(x).
!    a,b  : real(dp) values, containing the end points of the integration
!           interval. s(x) is considered to be identically zero outside
!           the interval (t(k+1),t(n-k)).
!
!  output parameter:
!    aint : real(dp), containing the integral of s(x) between a and b.
!    wrk  : real(dp) array, length n.  used as working space
!           on output, wrk will contain the integrals of the normalized
!           b-splines defined on the set of knots.
!
!  other subroutines required: fpintb.
!
!  references :
!    gaffney p.w. : the calculation of indefinite integrals of b-splines
!                   j. inst. maths applics 17 (1976) 37-41.
!    dierckx p. : curve and surface fitting with splines, monographs on
!                 numerical analysis, oxford university press, 1993.
!
!  author :
!    p.dierckx
!    dept. computer science, k.u.leuven
!    celestijnenlaan 200a, b-3001 heverlee, belgium.
!    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  latest update : march 1987
!
!  ..scalar arguments..
      real(dp) a,b
      integer n,k
!  ..array arguments..
      real(dp) t(n),c(n),wrk(n)
!  ..local scalars..
      integer i,nk1
!  ..
      nk1 = n-k-1
!  calculate the integrals wrk(i) of the normalized b-splines
!  ni,k+1(x), i=1,2,...nk1.
      call fpintb(t,n,wrk,nk1,a,b)
!  calculate the integral of s(x).
      splint = 0.
      do 10 i=1,nk1
        splint = splint+c(i)*wrk(i)
  10  continue
      return
      end function splint

  end module splines
