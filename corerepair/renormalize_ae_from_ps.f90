  subroutine renormalize_ae_from_ps( ae_n, ae_r, ae_wave, &
                                     ps_n, ps_r, ps_wave, &
                                     rmax, tnormcnsv, alpha )

! Given input all-electron and pseudo waves of the same angular momentum
! renormalize the all-electron wave to be consistent with the pseudo wave
! Consistent is defined as
!  (i) norm-conserving : ae_wave = ae_wave * alpha sth
!      alpha**2 * integral_0^R ae_wave**2 = integral_0^R ps_wave**2
!  (ii) otherwise : alpha * ae_wave(R) = ps_wave(R)
!  R an input parameter

  use kinds
  use splines

  implicit none

  integer,parameter :: ideg = 5
  real(dp),parameter :: zero = 0.d0

  integer :: ae_n, ps_n
  real(dp) :: ae_r(ae_n), ps_r(ps_n)
  real(dp) :: ae_wave(ae_n), ps_wave(ps_n)
  real(dp) :: rmax
  logical :: tnormcnsv
  real(dp) :: alpha

  real(dp) :: f(max(ae_n,ps_n))
  type(spline_struct) :: f_spl
  integer :: ierr
  real(dp) :: ae_norm, ps_norm
  real(dp) :: ae_rmax, ps_rmax

  if( tnormcnsv ) then
  ! norm conserving
    f(1:ae_n) = ae_wave(:)**2
    call spline_fit( ideg, ae_n, ae_r, f, f_spl, ierr, zero, rmax )
    if( ierr > 0 ) then
      write(0,*) 'renormalize_ae_from_ps : error 1 in spline_fit', ierr
      stop 1
    endif
    ae_norm = spline_integral( f_spl, zero, rmax )

    f(1:ps_n) = ps_wave(:)**2
    call spline_fit( ideg, ps_n, ps_r, f, f_spl, ierr, zero, rmax )
    if( ierr > 0 ) then
      write(0,*) 'renormalize_ae_from_ps : error 2 in spline_fit', ierr
      stop 1
    endif
    ps_norm = spline_integral( f_spl, zero, rmax )

    alpha = sqrt(ps_norm/ae_norm)
    write(*,*) 'ps norm = ', ps_norm
    write(*,*) 'ae norm = ', ae_norm
    write(*,*) ' factor = ', alpha
  else
    f(1:ae_n) = ae_wave(:)
    call spline_fit( ideg, ae_n, ae_r, f, f_spl, ierr, zero, rmax )
    if( ierr > 0 ) then
      write(0,*) 'renormalize_ae_from_ps : error 3 in spline_fit', ierr
      stop 1
    endif
    call spline_eval( f_spl, rmax, ae_rmax, ierr )

    f(1:ps_n) = ps_wave(:)
    call spline_fit( ideg, ps_n, ps_r, f, f_spl, ierr, zero, rmax )
    if( ierr > 0 ) then
      write(0,*) 'renormalize_ae_from_ps : error 4 in spline_fit', ierr
      stop 1
    endif
    call spline_eval( f_spl, rmax, ps_rmax, ierr )

    alpha = ps_rmax/ae_rmax
    write(*,*) 'ps(rmax) = ', ps_rmax
    write(*,*) 'ae(rmax) = ', ae_rmax
    write(*,*) '  factor = ', alpha
  endif

  return

  end subroutine renormalize_ae_from_ps
