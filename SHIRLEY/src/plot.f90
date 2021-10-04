  module plot_module

  use kinds, only : dp

  implicit none

  public

  type plotspec_type
    integer :: nplot
    integer,pointer :: ink(:,:)
  end type plotspec_type

  real(dp) :: alat
  real(dp) :: at(3,3)
  integer :: nat, ntyp
  real(dp),allocatable :: tau(:,:)
  character(3),allocatable :: atm(:)
  integer,allocatable :: ityp(:)
  integer :: nr1, nr2, nr3, nbnd

  contains

! ----------------------------------------------------------------------
  subroutine plotspec_read( plotspec, iunit )
! ----------------------------------------------------------------------
  type(plotspec_type) :: plotspec
  integer,intent(in) :: iunit
  integer :: i

  read(iunit,*) ! heading
  read(iunit,*) plotspec%nplot
  allocate( plotspec%ink(2,plotspec%nplot) )
  do i=1,plotspec%nplot
    read(iunit,*) plotspec%ink(1:2,i)
  enddo

  end subroutine plotspec_read


! ----------------------------------------------------------------------
  subroutine plotspec_bcast( plotspec, root, comm )
! ----------------------------------------------------------------------
  use mp, only : mp_bcast, mp_rank

  type(plotspec_type) :: plotspec
  integer,intent(in) :: root, comm

  integer :: mpime 

  mpime = mp_rank( comm )

  call mp_bcast( plotspec%nplot, root, comm )
  if( mpime /= root ) then
    allocate( plotspec%ink(2,plotspec%nplot) )
  endif
  call mp_bcast( plotspec%ink, root, comm )

  end subroutine plotspec_bcast
  

! ----------------------------------------------------------------------
  subroutine plotspec_scatter( plotspec, root, comm )
! ----------------------------------------------------------------------

  use mp, only : mp_bcast, mp_rank
  use mp_scatt, only: mp_scatter_size, mp_scatter

  type(plotspec_type) :: plotspec
  integer,intent(in) :: root, comm

  integer :: mpime
  type(plotspec_type) :: plotspec_l

  mpime = mp_rank( comm )

  ! broadcast size
  call mp_bcast( plotspec%nplot, root, comm )

  ! find local size
  call mp_scatter_size( plotspec%nplot, plotspec_l%nplot, root, comm )

  ! make local temporary space 
  allocate( plotspec_l%ink(2,plotspec_l%nplot) )

  ! scatter root's plotspec to local plotspec_l
  call mp_scatter( plotspec%ink, plotspec_l%ink, root, comm )

  ! don't forget to deallocate the kptlist on root
  if( mpime==root ) then
    if( associated(plotspec%ink) ) deallocate(plotspec%ink)
  endif

  ! dump contents of plotspec_l into resized plotspec
  allocate( plotspec%ink(2,plotspec_l%nplot) )
  plotspec = plotspec_l
  
  end subroutine plotspec_scatter


! ----------------------------------------------------------------------
  subroutine plot_basis_read_header( iunplt )
! ----------------------------------------------------------------------

  integer,intent(in) :: iunplt

  read(iunplt) alat
  read(iunplt) at

  read(iunplt) nat, ntyp
  if( allocated(tau) ) deallocate(tau)
  if( allocated(atm) ) deallocate(atm)
  if( allocated(ityp) ) deallocate(ityp)
  allocate( tau(3,nat), atm(ntyp), ityp(nat) )
  read(iunplt) tau
  read(iunplt) atm
  read(iunplt) ityp

  read(iunplt) nr1, nr2, nr3
  read(iunplt) nbnd

  end subroutine plot_basis_read_header


! ----------------------------------------------------------------------
  subroutine plot_basis_write_header( iunplt )
! ----------------------------------------------------------------------

  integer,intent(in) :: iunplt

  write(iunplt) alat
  write(iunplt) at

  write(iunplt) nat, ntyp
  write(iunplt) tau
  write(iunplt) atm
  write(iunplt) ityp

  write(iunplt) nr1, nr2, nr3
  write(iunplt) nbnd

  end subroutine plot_basis_write_header


! ----------------------------------------------------------------------
  subroutine plot_basis_save_header( alat_, at_, nat_, ntyp_, &
               tau_, atm_, ityp_, nr1_, nr2_, nr3_, nbnd_ )
! ----------------------------------------------------------------------

  real(dp),intent(in) :: alat_, at_(3,3)
  integer,intent(in) :: nat_, ntyp_
  real(dp),intent(in) :: tau_(3,*)
  character(3),intent(in) :: atm_(*)
  integer,intent(in) :: ityp_(*)
  integer,intent(in) :: nr1_, nr2_, nr3_, nbnd_

  alat = alat_
  at = at_
  nat = nat_
  ntyp = ntyp_
  if( allocated(tau) ) deallocate( tau )
  if( allocated(atm) ) deallocate( atm )
  if( allocated(ityp) ) deallocate( ityp )
  allocate( tau(3,nat), atm(ntyp), ityp(nat) )
  tau = tau_(1:3,1:nat)
  atm = atm_(1:ntyp)
  ityp = ityp_(1:nat)

  nr1=nr1_
  nr2=nr2_
  nr3=nr3_
  nbnd=nbnd_

  end subroutine plot_basis_save_header

  end module plot_module
