    function intradialprod( ngrid, r, f, g )

    use kinds, only : dp
    use splines

    real(dp) :: intradialprod
    integer,intent(in) :: ngrid
    real(dp),intent(in) :: r(ngrid)
    real(dp),intent(in) :: f(ngrid)
    real(dp),intent(in) :: g(ngrid)

    real(dp) :: p(ngrid)

    ! Spline variables
    integer :: iopt,n_knot,ier,ideg,nest,lwrk
    integer,allocatable :: iwrk(:)
    real(dp),allocatable :: weight(:),t_knot(:),coeff(:),wrk(:)
    real(dp) :: xb,xe,s,fp

    p = f*g

    ! Initialize spline dimensions
    ideg=5
    nest=ngrid+ideg+1
    lwrk=ngrid*(ideg+1)+nest*(7+3*ideg)
    allocate( weight(ngrid),t_knot(nest),coeff(nest),wrk(lwrk),iwrk(nest) )

    iopt=0
    weight=1.d0
    xb=0.d0
    xe=r(ngrid)
    s=0.d0
    n_knot=nest

    call curfit(iopt,ngrid,r,p,weight,xb,xe,         &
                ideg,s,nest,n_knot,t_knot,coeff,fp, &
                wrk,lwrk,iwrk,ier)

    intradialprod = splint(t_knot,n_knot,coeff,ideg,xb,xe,wrk)

    deallocate( weight,t_knot,coeff,wrk,iwrk )

    end function intradialprod

