  module fermi

  implicit none

  contains
! ======================================================================
  pure function fermifunc( e, ef, kT )
! ======================================================================
! Compute the Fermi-Dirac distribution function
  use kinds, only : dp
  real(dp),intent(in) :: e, ef, kT
  real(dp) :: fermifunc

  real(dp) :: de

  de=(e-ef)/kT
  if( abs(kT) < 1.d-10) then
    if( de > 0 ) then
      fermifunc = 0.d0
    else if( de < 0 ) then
      fermifunc = 1.d0
    else
      fermifunc = 0.5d0
    endif
  else
    if( de > 23.0d0 ) then
      fermifunc = 0.d0
    else if( de < -23.0d0 ) then
      fermifunc = 1.d0
    else
      fermifunc = 1.d0 / ( exp( de ) + 1 )
    endif
  endif
  end function fermifunc

! ======================================================================
  pure function fermideriv( e, ef, kT )
! ======================================================================
! Compute the Fermi-Dirac distribution function
  use kinds, only : dp
  real(dp),intent(in) :: e, ef, kT
  real(dp) :: fermideriv

  real(dp) :: de

  fermideriv=0.d0
  de=(e-ef)/kT
  if( abs(kT) >= 1.d-10) then
    if( de >= -23.0d0 .and. de <= 23.0d0 ) then
      fermideriv = -1.d0 / ( exp( de ) + 1 )**2 * exp( de )/kT
    endif
  endif
  end function fermideriv

! ======================================================================
  subroutine fermi_energy( nelec, nk, nbnd, wk, eigvalk, kT, smearing, ef )
! ======================================================================

  ! assumes parallel distribution of eigenvalues and k-points
  use kinds, only : dp
  use mp, only : mp_sum
  use mp_world, only : world_comm
  use io_global, only : stdout
  use shirley_constants, only : maxchar, rytoev, kelvin2rydberg

  real(dp),intent(in) :: nelec
  integer,intent(in) :: nk, nbnd
  real(dp),intent(in) :: wk(nk)
  real(dp),intent(in) :: eigvalk(nbnd,nk)
  real(dp),intent(in) :: kT
  character(maxchar),intent(in) :: smearing
  real(dp),intent(out) :: ef

  integer :: nbndf
  integer :: i, ik
  real(dp) :: nelecf, def
  real(dp) :: wktot, wkfac
  real(dp) :: wg(nbnd,nk)
  integer :: niter
  
  if( trim(smearing) == 'fermi-dirac' ) then
  else if( trim(smearing) == 'gaussian' ) then
  else if( trim(smearing) == 'lorentzian' ) then
  else
    call errore('fermi_energy','unrecognized smearing type: '//trim(smearing),1)
  endif
  write(stdout,*) ' smearing = ', trim(smearing)

  ! first guess
  ! assume spin-unpolarized
  ! normalize weights to 2.0
  wktot = sum( wk )
  call mp_sum( wktot, world_comm )
  wkfac = 2.d0 / wktot
  nbndf = int(nelec*0.5d0)
  ef = sum(eigvalk(nbndf,:)*wk(:)*wkfac)*0.5d0
  call mp_sum(ef, world_comm)

  ! first estimate of integral over bands to get nelecf
  forall( i=1:nbnd, ik=1:nk ) &
    wg(i,ik) = fermi_smearing( eigvalk(i,ik), ef, kT, smearing ) * wk(ik)*wkfac
  nelecf = sum(wg)
  call mp_sum( nelecf, world_comm )

  def = 1.d0
  if( nelecf > nelec ) def = -1.d0
  write(stdout,'(4a25)') 'ef', 'def', 'nelecf', 'nelec'
  write(stdout,'(4f16.8)') ef*rytoev, def, nelecf, nelec

  niter=0
  do while( abs(def) > 1.d-10 .and. abs(nelecf-nelec) > 1.d-10 .and. niter < 1000 )
    ! bisection
    if( def*(nelecf-nelec) > 0.d0 ) def = -0.5d0*def
    ef = ef + def
    forall( i=1:nbnd, ik=1:nk ) &
      wg(i,ik) = fermi_smearing( eigvalk(i,ik), ef, kT, smearing ) * wk(ik)*wkfac
    nelecf = sum(wg)
    call mp_sum( nelecf, world_comm )

    write(stdout,'(4f16.8)') ef*rytoev, def, nelecf, nelec

    niter=niter+1
  enddo

  if( niter>= 1000 ) call errore('fermi_energy','Fermi energy did not converge',1)

  end subroutine fermi_energy


! ======================================================================
  pure function fermi_smearing( eigval, ef, kT, smearing )
! ======================================================================

  use kinds, only : dp
  use shirley_constants, only : maxchar, sqrt2, pi

  real(dp),intent(in) :: eigval, ef, kT
  character(maxchar),intent(in) :: smearing
  real(dp) :: fermi_smearing
#ifdef __PGI
  real(dp),external :: erf
#endif

  
  if( trim(smearing) == 'fermi-dirac' ) then

    fermi_smearing = fermifunc( eigval, ef, kT )

  else if( trim(smearing) == 'gaussian' ) then

    fermi_smearing = 0.5d0 + 0.5d0*erf( (ef-eigval)/(sqrt2*kT) )

  else if( trim(smearing) == 'lorentzian' ) then

    fermi_smearing = 0.5d0 + atan( (ef-eigval)/kT ) / pi

  endif

  end function fermi_smearing


  end module fermi
