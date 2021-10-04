! ----------------------------------------------------------------------
  module corerepair_module
! ----------------------------------------------------------------------
  use kinds, only : dp

  implicit none
  private
  public :: core_type, corerepair_type
  public :: read_corerepair, bcast_corerepair

  type core_type
    integer :: species
    integer :: atom
    character(255) :: filename
    character(255) :: label
    integer :: nwfc1, nwfc2
    integer,pointer :: lwfc1(:), lwfc2(:)
    integer :: nproj1, nproj2
    integer :: nnonzero
    complex(dp),pointer :: matrix(:,:,:)
  end type core_type

  type corerepair_type
    integer :: nspecies, natom
    integer :: ncore
    type(core_type),pointer :: core(:)    
  end type corerepair_type

  contains

! ----------------------------------------------------------------------
  subroutine read_corerepair( iunit, corerep )
! ----------------------------------------------------------------------
  integer,intent(in) :: iunit
  type(corerepair_type) :: corerep

  integer :: i, iuncore, ierr
  integer,external :: freeunit

  read(iunit,*) ! heading
  read(iunit,*) corerep%nspecies, corerep%natom
  read(iunit,*) corerep%ncore
  allocate( corerep%core(corerep%ncore) )
  do i=1,corerep%ncore
    read(iunit,*) corerep%core(i)%species, corerep%core(i)%atom, &
                  corerep%core(i)%filename

    iuncore=freeunit()
    open(iuncore,file=trim(corerep%core(i)%filename),form='formatted', &
         iostat=ierr)
    if( ierr/=0 ) call errore( 'read_corerepair','unable to open file '// &
                               trim(corerep%core(i)%filename), 1 )

    call read_core( iuncore, corerep%core(i) )
    close(iuncore)
  enddo

  end subroutine read_corerepair

! ----------------------------------------------------------------------
  subroutine read_core( iunit, core )
! ----------------------------------------------------------------------
  integer,intent(in) :: iunit
  type(core_type) :: core

  integer :: i
  integer :: ip1, ip2, ixyz
  real(dp) :: cR, cI

  read(iunit,*) core%label
  read(iunit,*) core%nwfc1, core%nwfc2

  allocate( core%lwfc1(core%nwfc1) )
  allocate( core%lwfc2(core%nwfc2) )
  read(iunit,*) core%lwfc1(1:core%nwfc1)
  read(iunit,*) core%lwfc2(1:core%nwfc2)

  core%nproj1 = sum(2*core%lwfc1(1:core%nwfc1)+1)
  core%nproj2 = sum(2*core%lwfc2(1:core%nwfc2)+1)

  allocate( core%matrix(core%nproj1,core%nproj2,9) ) ! Yufeng: ixyz 1 to 9
  core%matrix = 0.d0
  read(iunit,*) core%nnonzero
  do i=1,core%nnonzero
    read(iunit,*) ip1, ip2, ixyz, cR, cI
    core%matrix(ip1,ip2,ixyz) = cmplx(cR,cI,kind=dp)
  enddo
  
  end subroutine read_core

! ----------------------------------------------------------------------
  subroutine bcast_corerepair( corerep, root, comm )
! ----------------------------------------------------------------------
  use mp, only : mp_bcast, mp_rank
  type(corerepair_type) :: corerep
  integer,intent(in) :: root, comm

  integer :: i
  integer :: mpime

  mpime = mp_rank( comm )

  call mp_bcast( corerep%nspecies, root, comm )
  call mp_bcast( corerep%natom,    root, comm )
  call mp_bcast( corerep%ncore,    root, comm )

  if( mpime/=root ) allocate( corerep%core(corerep%ncore) )
  do i=1,corerep%ncore
    call bcast_core( corerep%core(i), root, comm )
  enddo

  end subroutine bcast_corerepair

! ----------------------------------------------------------------------
  subroutine bcast_core( core, root, comm )
! ----------------------------------------------------------------------
  use mp, only : mp_bcast, mp_rank
  type(core_type) :: core
  integer,intent(in) :: root, comm

  integer :: i
  integer :: mpime

  mpime = mp_rank( comm )

  call mp_bcast( core%species, root, comm )
  call mp_bcast( core%atom, root, comm )
  call mp_bcast( core%filename, root, comm )

  call mp_bcast( core%label, root, comm )

  call mp_bcast( core%nwfc1, root, comm )
  call mp_bcast( core%nwfc2, root, comm )
  if( mpime/=root ) allocate(core%lwfc1(core%nwfc1))
  if( mpime/=root ) allocate(core%lwfc2(core%nwfc2))
  call mp_bcast( core%lwfc1, root, comm )
  call mp_bcast( core%lwfc2, root, comm )

  call mp_bcast( core%nproj1, root, comm )
  call mp_bcast( core%nproj2, root, comm )
  call mp_bcast( core%nnonzero, root, comm )
  if( mpime/=root ) allocate(core%matrix(core%nproj1,core%nproj2,9) ) ! Yufeng: ixyz 1 to 9
  call mp_bcast( core%matrix, root, comm )

  end subroutine bcast_core


  end module corerepair_module
