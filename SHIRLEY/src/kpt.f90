  module kpt_module

  use kinds, only : dp
  use io_global, only : stdout

  implicit none

  private
  public :: kpt_read, kpt_scatter, kpt_bcast
  public :: kpt_type, kptparam_type, bandstructure_type, &
            kptlist_type, kpttetra_type, kpt_units
!  public :: kptpool_type

  real(dp),parameter :: one=1.d0
  real(dp),parameter :: two=2.d0
  real(dp),parameter :: zero=0.d0
  real(dp),parameter :: half=0.5d0

  type kptparam_type
    logical :: cartesian
    character(255) :: grid_type
    integer :: nktot
    integer :: nkgrid(3), ikgrid(3)
    character(255) :: shift_type
    real(dp) :: shift(3)
  end type kptparam_type

  type bandstructure_type
    integer :: nksp, kres, nk
    real(dp),pointer :: xksp(:,:)
    character(255),pointer :: labelsp(:)
    real(dp),pointer :: kpathlensp(:)
    real(dp),pointer :: kpathlen(:)
  end type bandstructure_type

  type kptlist_type
    integer :: nk
    integer :: nspin
    real(dp),pointer :: kvec(:,:)
    real(dp),pointer :: wk(:)
  end type kptlist_type

  type kpttetra_type
    integer :: ntetra
    integer,pointer :: tetra(:,:)
  end type kpttetra_type

  type kpt_type
    type(kptparam_type) param
    type(bandstructure_type) bandstructure
    type(kptlist_type) list
    type(kpttetra_type) tetra
  end type kpt_type

  type kptpool_type
    integer :: npool            ! number of pools
    integer :: my_pool          ! my pool index
    integer,pointer :: my_displ ! displacement of pool
    integer,pointer :: my_nk    ! number of k in pool
    integer :: mpiroot          ! MPI root for this pool
    integer :: mpime            ! MPI pid  for this pool
  end type kptpool_type

  contains


! ----------------------------------------------------------------------
  subroutine kpt_type_init( kpt )
! ----------------------------------------------------------------------

  type(kpt_type) :: kpt

  call bandstructure_type_init( kpt%bandstructure )
  call kptlist_type_init( kpt%list )
  call kpttetra_type_init( kpt%tetra )

  end subroutine kpt_type_init


! ----------------------------------------------------------------------
  subroutine bandstructure_type_init( bandstructure )
! ----------------------------------------------------------------------

  type(bandstructure_type) :: bandstructure

  nullify( bandstructure%xksp, &
           bandstructure%labelsp, &
           bandstructure%kpathlensp, &
           bandstructure%kpathlen )

  end subroutine bandstructure_type_init


! ----------------------------------------------------------------------
  subroutine kptlist_type_init( kptlist )
! ----------------------------------------------------------------------

  type(kptlist_type) :: kptlist

  nullify( kptlist%kvec, &
           kptlist%wk )

  end subroutine kptlist_type_init

! ----------------------------------------------------------------------
  subroutine kpttetra_type_init( kpttetra )
! ----------------------------------------------------------------------

  type(kpttetra_type) :: kpttetra

  nullify( kpttetra%tetra )

  end subroutine kpttetra_type_init


! ----------------------------------------------------------------------
  subroutine kpt_read( input_unit, kpt, nspin )
! ----------------------------------------------------------------------

  integer,intent(in) :: input_unit
  type(kpt_type),intent(out) :: kpt
  integer,optional :: nspin

  read(input_unit,*) ! heading

  call kpt_type_init( kpt )

  if( present(nspin) ) then
    kpt%list%nspin = nspin
  else
    kpt%list%nspin = 1
  endif

  read(input_unit,*) kpt%param%grid_type ! type of coordinates
  kpt%param%grid_type = adjustl(kpt%param%grid_type)

  kpt%param%cartesian = .false.
  if( trim(kpt%param%grid_type) == 'tpiba' .or. &
      trim(kpt%param%grid_type) == 'crystal' ) then

    call kpt_read_userlist( input_unit, kpt%list )

    if( trim(kpt%param%grid_type) == 'tpiba' ) kpt%param%cartesian = .true.

  else if( trim(kpt%param%grid_type) == 'automatic' ) then

    call kpt_read_automatic( input_unit, kpt%param, kpt%list )

  else if( trim(kpt%param%grid_type) == 'fermisurface' ) then

    call kpt_read_fermisurface( input_unit, kpt%param, kpt%list )

  else if( trim(kpt%param%grid_type) == 'shifted' ) then

    call kpt_read_shifted( input_unit, kpt%param, kpt%list )

  elseif( trim(kpt%param%grid_type) == 'bandstructure' ) then

    kpt%param%cartesian = .true.
    call kpt_read_bandstructure( input_unit, kpt%bandstructure, kpt%list )

  elseif( trim(kpt%param%grid_type) == 'random' ) then

    kpt%param%cartesian = .false.
    call kpt_read_random( input_unit, kpt%list )

  else if( trim(kpt%param%grid_type) == 'timereversed' ) then

    call kpt_read_timereversed( input_unit, kpt%param, kpt%list )

  else if( trim(kpt%param%grid_type) == 'tetrahedra' ) then

    call kpt_read_timereversed( input_unit, kpt%param, kpt%list )
    call tetrahedra( kpt%param, kpt%list, kpt%tetra, .true. )

  else

    write(stdout,*) ' kpoints flag unrecognized'
    write(stdout,*) ' should be: tpiba, crystal, bandstructure'
    write(stdout,*) '            fermisurface or automatic'
    call errore('kpt_read','stopping',1)

  endif

  ! store total number of k-points in parameters aswell
  kpt%param%nktot = kpt%list%nk

  end subroutine kpt_read


! ----------------------------------------------------------------------
  subroutine kpt_read_userlist( input_unit, kptlist )
! ----------------------------------------------------------------------

  integer,intent(in) :: input_unit
  type(kptlist_type),intent(out) :: kptlist

  integer :: ik

  read(input_unit,*) kptlist%nk

  allocate( kptlist%kvec(3,kptlist%nk), kptlist%wk(kptlist%nk) )
  do ik=1,kptlist%nk
    read(input_unit,*) kptlist%kvec(1:3,ik), kptlist%wk(ik)
  enddo

  end subroutine kpt_read_userlist


! ----------------------------------------------------------------------
  subroutine kpt_read_automatic( input_unit, kptparam, kptlist )
! ----------------------------------------------------------------------

  integer,intent(in) :: input_unit
  type(kptparam_type),intent(out) :: kptparam
  type(kptlist_type),intent(out) :: kptlist

  integer :: i,j,k

  read(input_unit,*) kptparam%nkgrid, kptparam%ikgrid

  if( any(kptparam%ikgrid > 1) .or. any(kptparam%ikgrid < 0) ) &
    call errore('kpt_module','grid shift values should be 0 or 1',201)

  if( any(kptparam%nkgrid < 1) ) &
    call errore('kpt_module','grid values should be positive',202)

  kptlist%nk = product(kptparam%nkgrid)
  allocate( kptlist%kvec(3,kptlist%nk), kptlist%wk(kptlist%nk) )

  kptlist%nk=0
  do i=1,kptparam%nkgrid(1)
    do j=1,kptparam%nkgrid(2)
      do k=1,kptparam%nkgrid(3)
        kptlist%nk=kptlist%nk+1
        !  xk are the components of the complete grid in crystal coords
        kptlist%kvec(1,kptlist%nk) = DBLE(i-1)/kptparam%nkgrid(1)         &
          + DBLE(kptparam%ikgrid(1))/2/kptparam%nkgrid(1)
        kptlist%kvec(2,kptlist%nk) = DBLE(j-1)/kptparam%nkgrid(2)         &
          + DBLE(kptparam%ikgrid(2))/2/kptparam%nkgrid(2)
        kptlist%kvec(3,kptlist%nk) = DBLE(k-1)/kptparam%nkgrid(3)         &
          + DBLE(kptparam%ikgrid(3))/2/kptparam%nkgrid(3)
      enddo
    enddo
  enddo
  kptlist%wk = two/dble(kptlist%nspin)/dble(kptlist%nk)

  end subroutine kpt_read_automatic


! ----------------------------------------------------------------------
  subroutine kpt_read_fermisurface( input_unit, kptparam, kptlist )
! ----------------------------------------------------------------------

  integer,intent(in) :: input_unit
  type(kptparam_type),intent(out) :: kptparam
  type(kptlist_type),intent(out) :: kptlist

  integer :: i,j,k

  read(input_unit,*) kptparam%nkgrid, kptparam%ikgrid

  if( any(kptparam%ikgrid > 1) .or. any(kptparam%ikgrid < 0) ) &
    call errore('kpt_module','grid shift values should be 0 or 1',201)

  if( any(kptparam%nkgrid < 1) ) &
    call errore('kpt_module','grid values should be positive',202)

  kptlist%nk = (kptparam%nkgrid(1)+1) &
             * (kptparam%nkgrid(2)+1) &
             * (kptparam%nkgrid(3)+1) 
  allocate( kptlist%kvec(3,kptlist%nk), kptlist%wk(kptlist%nk) )

  kptlist%nk=0
  do i=1,kptparam%nkgrid(1)+1
    do j=1,kptparam%nkgrid(2)+1
      do k=1,kptparam%nkgrid(3)+1
        kptlist%nk=kptlist%nk+1
        !  xk are the components of the complete grid in crystal coords
        kptlist%kvec(1,kptlist%nk) = DBLE(i-1)/kptparam%nkgrid(1)         &
          - DBLE(kptparam%ikgrid(1))/2/kptparam%nkgrid(1)
        kptlist%kvec(2,kptlist%nk) = DBLE(j-1)/kptparam%nkgrid(2)         &
          - DBLE(kptparam%ikgrid(2))/2/kptparam%nkgrid(2)
        kptlist%kvec(3,kptlist%nk) = DBLE(k-1)/kptparam%nkgrid(3)         &
          - DBLE(kptparam%ikgrid(3))/2/kptparam%nkgrid(3)
      enddo
    enddo
  enddo
  kptlist%wk = two/dble(kptlist%nspin)/dble(kptlist%nk)

  end subroutine kpt_read_fermisurface


! ----------------------------------------------------------------------
  subroutine kpt_read_shifted( input_unit, kptparam, kptlist )
! ----------------------------------------------------------------------

  integer,intent(in) :: input_unit
  type(kptparam_type),intent(out) :: kptparam
  type(kptlist_type),intent(out) :: kptlist

  integer :: i,j,k
  real(dp) :: trn(3)

  read(input_unit,*) kptparam%nkgrid, kptparam%shift_type

  if( trim(kptparam%shift_type) == 'random' .or. &
      trim(kptparam%shift_type) == 'user') then

    if( any(kptparam%nkgrid < 1) ) &
      call errore('kpt_module','grid values should be positive',202)

    kptlist%nk = product(kptparam%nkgrid)
    kptparam%ikgrid = 0 ! no half-shift here
    allocate( kptlist%kvec(3,kptlist%nk), kptlist%wk(kptlist%nk) )

    kptlist%nk=0
    do i=1,kptparam%nkgrid(1)
      do j=1,kptparam%nkgrid(2)
        do k=1,kptparam%nkgrid(3)
          kptlist%nk=kptlist%nk+1
          !  kvec are the components of the complete grid in crystal coords
          kptlist%kvec(1,kptlist%nk) = DBLE(i-1)/kptparam%nkgrid(1)
          kptlist%kvec(2,kptlist%nk) = DBLE(j-1)/kptparam%nkgrid(2)
          kptlist%kvec(3,kptlist%nk) = DBLE(k-1)/kptparam%nkgrid(3)
        enddo
      enddo
    enddo
    ! weights
    kptlist%wk = two/dble(kptlist%nspin)/dble(kptlist%nk)

    ! now shift 
    if( trim(kptparam%shift_type) == 'random' ) then
      write(stdout,*) ' k-point grid - adopting random shift:'
      call random_seed ! initialize random number generator
      ! generate translation
      call random_number( trn )
      trn = trn * two - one ! stretch to [-1:1]

      kptparam%shift = trn

    else ! user-defined
      write(stdout,*) ' k-point grid - adopting user shift:'
      read(input_unit,*) (kptparam%shift(i), i=1,3)
    endif

    ! report
    write(stdout,*) ' to generate the same k-grid, use the following input:'
    write(stdout,*) ' ====================================================='
    write(stdout,*) 'K_POINTS'
    write(stdout,*) ' shifted'
    write(stdout,'(3i6,a)') kptparam%nkgrid, ' user'
    write(stdout,'(12e12.5)')  kptparam%shift
    write(stdout,*) ' ====================================================='

    ! shift
    forall( i=1:3, j=1:kptlist%nk ) &
      kptlist%kvec(i,j) = kptlist%kvec(i,j) + trn(i)

    ! reduce back to BZ [0,1)^3
    forall(i=1:3, j=1:kptlist%nk) &
      kptlist%kvec(i,j) = kptlist%kvec(i,j) - floor(kptlist%kvec(i,j))

    ! report
    write(stdout,*) ' k-points:'
    do i=1,kptlist%nk
      write(stdout,'(2x,i6,2x,3e12.5,2x,e12.5)') i, kptlist%kvec(1:3,i), kptlist%wk(i)
    enddo
    write(stdout,*)

  else

    call errore('kpt_module','rotateshift_type must be random or user',201)

  endif
 
  end subroutine kpt_read_shifted


! ----------------------------------------------------------------------
  subroutine kpt_read_bandstructure( input_unit, bs, kptlist )
! ----------------------------------------------------------------------

  integer,intent(in) :: input_unit
  type(bandstructure_type),intent(out) :: bs
  type(kptlist_type),intent(out) :: kptlist

  integer :: i, j
  real(dp) :: lambda, dk, dkvec(3)

  read(input_unit,*) bs%nksp, bs%kres

  allocate( bs%xksp(3,bs%nksp), bs%labelsp(bs%nksp) )

  do i=1,bs%nksp
    read(input_unit,*) bs%xksp(:,i), bs%labelsp(i)
  enddo

  kptlist%nk=1
  do i=2,bs%nksp
    kptlist%nk=kptlist%nk+bs%kres
  enddo
  allocate( kptlist%kvec(3,kptlist%nk), kptlist%wk(kptlist%nk) )
  kptlist%nk=1
  kptlist%kvec(:,kptlist%nk) = bs%xksp(:,1)
  do i=2,bs%nksp
    do j=1,bs%kres
      lambda=dble(j)/dble(bs%kres)
      kptlist%nk=kptlist%nk+1
      kptlist%kvec(:,kptlist%nk) = (one-lambda)*bs%xksp(:,i-1) &
                                 + lambda*bs%xksp(:,i)
    enddo
  enddo
  kptlist%wk = two/dble(kptlist%nspin) ! kind of irrelevant for band structure

! now make the path length
  allocate( bs%kpathlensp(bs%nksp) )
  bs%kpathlensp(1) = 0.d0
  do i=2,bs%nksp
    dkvec = bs%xksp(:,i) - bs%xksp(:,i-1)
    dk = sqrt( dot_product( dkvec, dkvec ) )
    bs%kpathlensp(i) = bs%kpathlensp(i-1) + dk
  enddo

! and also for the list
  bs%nk=kptlist%nk
  allocate( bs%kpathlen(bs%nk) )
  bs%kpathlen(1) = 0.d0
  do i=2,bs%nk
    dkvec = kptlist%kvec(:,i) - kptlist%kvec(:,i-1)
    dk = sqrt( dot_product( dkvec, dkvec ) )
    bs%kpathlen(i) = bs%kpathlen(i-1) + dk
  enddo

  end subroutine kpt_read_bandstructure


! ----------------------------------------------------------------------
  subroutine kpt_read_random( input_unit, kptlist )
! ----------------------------------------------------------------------

  integer,intent(in) :: input_unit
  type(kptlist_type),intent(out) :: kptlist

  integer :: ik

  read(input_unit,*) kptlist%nk

  call random_seed

  allocate( kptlist%kvec(3,kptlist%nk), kptlist%wk(kptlist%nk) )
  do ik=1,kptlist%nk
    call random_number( kptlist%kvec(1:3,ik) )
  enddo
  kptlist%wk = 2.d0/dble(kptlist%nspin)/dble(kptlist%nk)

  ! report
  write(stdout,*) 'K_POINTS'
  write(stdout,*) '  random'
  write(stdout,*) kptlist%nk
  do ik=1,kptlist%nk
    write(stdout,'(2x,4f12.5)') kptlist%kvec(1:3,ik), kptlist%wk(ik)
  enddo

  end subroutine kpt_read_random


! ----------------------------------------------------------------------
  subroutine kpt_scatter( kptlist, kptparam, root, comm )
! ----------------------------------------------------------------------

  ! *************************************************************************
  ! This routine should be rewritten to reflect that kpathlen has been moved
  ! to the bandstructure part of kpt and is unnecessary for kptlist
  ! - davegp
  ! *************************************************************************

  use mp, only : mp_bcast, mp_barrier, mp_rank, mp_size
  use mp_scatt, only : mp_scatter_size, mp_scatter
  use hamq_pool, only : rootpool, mypool, cross_pool_comm, &
                        mypoolroot, mypoolid, intra_pool_comm

  type(kptlist_type),intent(inout) :: kptlist
  type(kptparam_type),intent(inout) :: kptparam

  integer,intent(in) :: root, comm

  integer :: mpime

  type(kptlist_type) :: kptlist_l
  integer :: i, ik
  logical :: havekpath

  mpime = mp_rank( comm )

  ! broadcast parameters
  call mp_bcast( kptparam%cartesian, root, comm )
  call mp_bcast( kptparam%grid_type, root, comm )
  call mp_bcast( kptparam%nktot,     root, comm )
  call mp_bcast( kptparam%nkgrid,    root, comm )
  call mp_bcast( kptparam%ikgrid,    root, comm )
  call mp_bcast( kptparam%shift_type,root, comm )
  call mp_bcast( kptparam%shift,     root, comm )

  havekpath = .false.
  if( trim(kptparam%grid_type) == 'bandstructure' ) then
    havekpath = .true.
  endif

  ! find local size
  if( mypoolid==mypoolroot ) &
    call mp_scatter_size( kptparam%nktot, kptlist_l%nk, rootpool, cross_pool_comm )
  call mp_bcast( kptlist_l%nk, mypoolroot, intra_pool_comm )

  ! make local temporary space kptlist_l
  allocate( kptlist_l%kvec(3,kptlist_l%nk), &
            kptlist_l%wk(kptlist_l%nk) )

  ! scatter root's kptlist to local kptlist_l
  if( mypoolid==mypoolroot ) &
    call mp_scatter( kptlist%kvec, kptlist_l%kvec, rootpool, cross_pool_comm )
  call mp_bcast( kptlist_l%kvec, mypoolroot, intra_pool_comm )

  if( mypoolid==mypoolroot ) &
    call mp_scatter( kptlist%wk,   kptlist_l%wk, rootpool, cross_pool_comm )
  call mp_bcast( kptlist_l%wk, mypoolroot, intra_pool_comm )

  ! don't forget to deallocate the kptlist on root
  if( mpime==root ) then
    if( associated(kptlist%kvec) ) deallocate(kptlist%kvec)
    if( associated(kptlist%wk) ) deallocate(kptlist%wk)
  endif

  ! dump contents of kptlist_l into resized kptlist
  allocate( kptlist%kvec(3,kptlist_l%nk),  &
            kptlist%wk(kptlist_l%nk) )
  
  ! report
  write(stdout,'(a,i6)') ' # k-points globally = ', kptparam%nktot
  do i=1,mp_size( comm )
    if( mpime==i-1 ) then
      write(*,'(a,i4,a,i6)') ' # k-points locally on processor ', i, &
                      ' = ', kptlist%nk
      do ik=1,kptlist%nk
        write(*,'(4x,i4,3(2x,e12.6),4x,e12.6)') &
          ik, kptlist%kvec(:,ik), kptlist%wk(ik)
      enddo
    endif
    call mp_barrier( comm )
  enddo

  end subroutine kpt_scatter


! ----------------------------------------------------------------------
  subroutine kpt_bcast( kpt, root, comm )
! ----------------------------------------------------------------------

  type(kpt_type) :: kpt
  integer,intent(in) :: root, comm

  call kptparam_bcast( kpt%param, root, comm )
  if( trim(kpt%param%grid_type) == 'bandstructure' ) then
    call kptbandstructure_bcast( kpt%bandstructure, root, comm )
  endif
  call kptlist_bcast( kpt%list, root, comm )
  if( trim(kpt%param%grid_type) == 'tetrahedra' ) then
    call kpttetra_bcast( kpt%tetra, root, comm )
  endif

  end subroutine kpt_bcast


! ----------------------------------------------------------------------
  subroutine kptparam_bcast( kptparam, root, comm ) 
! ----------------------------------------------------------------------

  use mp, only : mp_bcast

  type(kptparam_type) :: kptparam
  integer,intent(in) :: root, comm

  ! broadcast parameters
  call mp_bcast( kptparam%cartesian, root, comm )
  call mp_bcast( kptparam%grid_type, root, comm )
  call mp_bcast( kptparam%nktot,     root, comm )
  call mp_bcast( kptparam%nkgrid,    root, comm )
  call mp_bcast( kptparam%ikgrid,    root, comm )
  call mp_bcast( kptparam%shift_type,root, comm )
  call mp_bcast( kptparam%shift,     root, comm )

  end subroutine kptparam_bcast


! ----------------------------------------------------------------------
  subroutine kptbandstructure_bcast( kptbs, root, comm ) 
! ----------------------------------------------------------------------

  use mp, only : mp_bcast, mp_rank

  type(bandstructure_type) kptbs
  integer,intent(in) :: root, comm

  logical :: isassoc(4)
  integer :: mpime

  mpime = mp_rank( comm )

  call mp_bcast( kptbs%nksp, root, comm )
  call mp_bcast( kptbs%kres, root, comm )
  call mp_bcast( kptbs%nk, root, comm )
  isassoc = .false.
  if( mpime == root ) then
    if( associated(kptbs%xksp) )       isassoc(1)=.true.
    if( associated(kptbs%labelsp) )    isassoc(2)=.true.
    if( associated(kptbs%kpathlensp) ) isassoc(3)=.true.
    if( associated(kptbs%kpathlen) )   isassoc(4)=.true.
  endif
  call mp_bcast( isassoc, root, comm )
  if( mpime /= root ) then
    if( isassoc(1) ) allocate( kptbs%xksp(3,kptbs%nksp) )
    if( isassoc(2) ) allocate( kptbs%labelsp(kptbs%nksp) )
    if( isassoc(3) ) allocate( kptbs%kpathlensp(kptbs%nksp) )
    if( isassoc(4) ) allocate( kptbs%kpathlen(kptbs%nk) )
  endif
  call mp_bcast( kptbs%xksp, root, comm )
  call mp_bcast( kptbs%labelsp, root, comm )
  call mp_bcast( kptbs%kpathlensp, root, comm )
  call mp_bcast( kptbs%kpathlen, root, comm )

  end subroutine kptbandstructure_bcast

 
! ----------------------------------------------------------------------
  subroutine kptlist_bcast( kptlist, root, comm ) 
! ----------------------------------------------------------------------

  use mp, only : mp_bcast, mp_rank

  type(kptlist_type),intent(inout) :: kptlist
  integer,intent(in) :: root, comm

  integer :: mpime
  logical :: isassoc(2)

  mpime = mp_rank( comm )

  call mp_bcast( kptlist%nk,    root, comm )
  call mp_bcast( kptlist%nspin, root, comm )
  isassoc = .false.
  if( mpime == root ) then
    if( associated(kptlist%kvec) )     isassoc(1)=.true.
    if( associated(kptlist%wk) )       isassoc(2)=.true.
  endif
  call mp_bcast( isassoc, root, comm )
  if( mpime /= root ) then
    if( isassoc(1) ) allocate( kptlist%kvec(3,kptlist%nk) )
    if( isassoc(2) ) allocate( kptlist%wk(kptlist%nk) )
  endif
  if( isassoc(1) ) call mp_bcast( kptlist%kvec, root, comm )
  if( isassoc(2) ) call mp_bcast( kptlist%wk, root, comm )

  end subroutine kptlist_bcast


! ----------------------------------------------------------------------
  subroutine kpttetra_bcast( kpttetra, root, comm ) 
! ----------------------------------------------------------------------

  use mp, only : mp_bcast, mp_rank

  type(kpttetra_type),intent(inout) :: kpttetra
  integer,intent(in) :: root, comm

  integer :: mpime
  logical :: isassoc(1)

  mpime = mp_rank( comm )

  call mp_bcast( kpttetra%ntetra, root, comm )
  isassoc = .false.
  if( mpime == root ) then
    if( associated(kpttetra%tetra) )     isassoc(1)=.true.
  endif
  call mp_bcast( isassoc, root, comm )
  if( mpime /= root ) then
    if( isassoc(1) ) allocate( kpttetra%tetra(4,kpttetra%ntetra) )
  endif
  call mp_bcast( kpttetra%tetra, root, comm )

  end subroutine kpttetra_bcast


!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine kpt_read_timereversed( input_unit, kptparam, kptlist )
!-----------------------------------------------------------------------
! ----------------------------------------------------------------------
!  subroutine kpt_read_automatic( input_unit, kptparam, kptlist )
! ----------------------------------------------------------------------
!
!  Automatic generation of a uniform grid of k-points
!
  integer,intent(in) :: input_unit
  type(kptparam_type),intent(out) :: kptparam
  type(kptlist_type),intent(out) :: kptlist

  real(DP), parameter :: eps=1.0e-5
  real(dp) :: epsk(3)
  real(DP) xkr(3), fact, xx, yy, zz
  real(DP), allocatable:: xkg(:,:), wkk(:)
  integer nkr, i,j,k, n, nk
  integer, allocatable :: equiv(:)
  logical :: in_the_list
  !
  read(input_unit,*) kptparam%nkgrid, kptparam%ikgrid

  if( any(kptparam%ikgrid > 1) .or. any(kptparam%ikgrid < 0) ) &
    call errore('kpt_module','grid shift values should be 0 or 1',201)

  if( any(kptparam%nkgrid < 1) ) &
    call errore('kpt_module','grid values should be positive',202)

  nkr = product(kptparam%nkgrid)
  epsk = eps / dble(nkr)
  allocate (xkg( 3,nkr),wkk(nkr))    
  allocate (equiv( nkr))    
  !
  do i=1,kptparam%nkgrid(1)
     do j=1,kptparam%nkgrid(2)
        do k=1,kptparam%nkgrid(3)
           !  this is nothing but consecutive ordering
           n = (k-1) + (j-1)*kptparam%nkgrid(3) + (i-1)*kptparam%nkgrid(2)*kptparam%nkgrid(3) + 1
           !  xkg are the components of the complete grid in crystal axis
           xkg(1,n) = DBLE(i-1)/kptparam%nkgrid(1) + DBLE(kptparam%ikgrid(1))/2/kptparam%nkgrid(1)
           xkg(2,n) = DBLE(j-1)/kptparam%nkgrid(2) + DBLE(kptparam%ikgrid(2))/2/kptparam%nkgrid(2)
           xkg(3,n) = DBLE(k-1)/kptparam%nkgrid(3) + DBLE(kptparam%ikgrid(3))/2/kptparam%nkgrid(3)
        end do
     end do
  end do

  !  equiv(nk) =nk : k-point nk is not equivalent to any previous k-point
  !  equiv(nk)!=nk : k-point nk is equivalent to k-point equiv(nk)

  do nk=1,nkr
     equiv(nk)=nk
  end do

  do nk=1,nkr
     !  check if this k-point has already been found equivalent to another
     if (equiv(nk).eq.nk) then
        wkk(nk)   = 1.0
        !  check if there are equivalent k-point to this in the list
        !  (excepted those previously found to be equivalent to another)
        !  check both k and -k
        xkr(:) = xkg(:,nk)
        xkr(:) = xkr(:) - nint( xkr(:) )
        xx = xkr(1)*kptparam%nkgrid(1) - 0.5d0*kptparam%ikgrid(1)
        yy = xkr(2)*kptparam%nkgrid(2) - 0.5d0*kptparam%ikgrid(2)
        zz = xkr(3)*kptparam%nkgrid(3) - 0.5d0*kptparam%ikgrid(3)
        in_the_list=abs(xx-nint(xx)).le.epsk(1).and.abs(yy-nint(yy)).le.epsk(2) &
                                               .and.abs(zz-nint(zz)).le.epsk(3) 
        if (in_the_list) then
           i = mod ( nint ( xkr(1)*kptparam%nkgrid(1) - 0.5 * kptparam%ikgrid(1) &
             + 2*kptparam%nkgrid(1)), kptparam%nkgrid(1) ) + 1
           j = mod ( nint ( xkr(2)*kptparam%nkgrid(2) - 0.5 * kptparam%ikgrid(2) &
             + 2*kptparam%nkgrid(2)), kptparam%nkgrid(2) ) + 1
           k = mod ( nint ( xkr(3)*kptparam%nkgrid(3) - 0.5 * kptparam%ikgrid(3) &
             + 2*kptparam%nkgrid(3)), kptparam%nkgrid(3) ) + 1
           n = (k-1) + (j-1)*kptparam%nkgrid(3) + (i-1)*kptparam%nkgrid(2)*kptparam%nkgrid(3) + 1
           if (n.gt.nk .and. equiv(n).eq.n) then
              equiv(n) = nk
              wkk(nk)=wkk(nk)+1.0
           else
              if (equiv(n).ne.nk.or.n.lt.nk) call errore('kpoint_grid', &
              'something wrong in the checking algorithm',2)
           end if
        end if
        ! now check -k
        xx =-xkr(1)*kptparam%nkgrid(1) - 0.5d0*kptparam%ikgrid(1)
        yy =-xkr(2)*kptparam%nkgrid(2) - 0.5d0*kptparam%ikgrid(2)
        zz =-xkr(3)*kptparam%nkgrid(3) - 0.5d0*kptparam%ikgrid(3)
        in_the_list=abs(xx-nint(xx)).le.epsk(1).and.abs(yy-nint(yy)).le.epsk(2) &
                                               .and.abs(zz-nint(zz)).le.epsk(3) 
        if (in_the_list) then
           i = mod ( nint (-xkr(1)*kptparam%nkgrid(1) - 0.5 * kptparam%ikgrid(1) &
             + 2*kptparam%nkgrid(1)), kptparam%nkgrid(1) ) + 1
           j = mod ( nint (-xkr(2)*kptparam%nkgrid(2) - 0.5 * kptparam%ikgrid(2) &
             + 2*kptparam%nkgrid(2)), kptparam%nkgrid(2) ) + 1
           k = mod ( nint (-xkr(3)*kptparam%nkgrid(3) - 0.5 * kptparam%ikgrid(3) &
             + 2*kptparam%nkgrid(3)), kptparam%nkgrid(3) ) + 1
           n = (k-1) + (j-1)*kptparam%nkgrid(3) + (i-1)*kptparam%nkgrid(2)*kptparam%nkgrid(3) + 1
           if (n.gt.nk .and. equiv(n).eq.n) then
              equiv(n) = nk
              wkk(nk)=wkk(nk)+1.0
           else
              if (equiv(n).ne.nk.or.n.lt.nk) call errore('kpoint_grid', &
              'something wrong in the checking algorithm',2)
           end if
        end if
     endif
  end do

  !  count time-reversed points and order them

  kptlist%nk=0
  do nk=1,nkr
     if (equiv(nk).eq.nk) then
        kptlist%nk=kptlist%nk+1
     endif
  enddo
  allocate( kptlist%kvec(3,kptlist%nk), kptlist%wk(kptlist%nk) )
  kptlist%nk=0
  fact=0.0
  do nk=1,nkr
     if (equiv(nk).eq.nk) then
        kptlist%nk=kptlist%nk+1
        kptlist%wk(kptlist%nk) = wkk(nk)
        fact    = fact+kptlist%wk(kptlist%nk)
        !  bring back into to the first BZ
        do i=1,3
           kptlist%kvec(i,kptlist%nk) = xkg(i,nk)-nint(xkg(i,nk))
        end do
     end if
  end do
  !  normalize weights to 2.d0 (including spin)
  do nk=1,kptlist%nk
     kptlist%wk(nk) = kptlist%wk(nk)/fact*(2.d0/dble(kptlist%nspin))
  end do

  deallocate(equiv)
  deallocate(xkg,wkk)

  return
  end subroutine kpt_read_timereversed

  !-----------------------------------------------------------------------
  subroutine tetrahedra ( kptparam, kptlist, kpttetra, minus_q )
  !-----------------------------------------------------------------------
  !
  ! Tetrahedron method according to P. E. Bloechl et al, PRB49, 16223 (1994)
  !
  ! INPUT:
  type(kptparam_type),intent(in) :: kptparam
  type(kptlist_type),intent(in) :: kptlist
  logical,intent(in) :: minus_q  ! true if we have time-reversal
  ! OUTPUT:
  type(kpttetra_type),intent(out) :: kpttetra
  ! LOCAL:
  real(DP) :: xkr(3), deltap(3), deltam(3)
  real(DP), parameter:: eps=1.0d-5
  real(DP) :: epsk(3)
  real(DP), allocatable :: xkg(:,:)
  integer :: nkr, i,j,k, n, nk, ip1,jp1,kp1, &
       n1,n2,n3,n4,n5,n6,n7,n8
  integer, allocatable:: equiv(:)
  !
  ! Re-generate a uniform grid of k-points xkg
  !
  nkr=product(kptparam%nkgrid)
  epsk = eps / dble(nkr)
  kpttetra%ntetra=6*nkr
  allocate( kpttetra%tetra(4,kpttetra%ntetra) )

  allocate (xkg( 3,nkr))    
  allocate (equiv( nkr))    
!
  do i=1,kptparam%nkgrid(1)
     do j=1,kptparam%nkgrid(2)
        do k=1,kptparam%nkgrid(3)
           !  this is nothing but consecutive ordering
           n = (k-1) + (j-1)*kptparam%nkgrid(3) + (i-1)*kptparam%nkgrid(2)*kptparam%nkgrid(3) + 1
           !  xkg are the components of the complete grid in crystal axis
           xkg(1,n) = DBLE(i-1)/kptparam%nkgrid(1) + DBLE(kptparam%ikgrid(1))/2/kptparam%nkgrid(1)
           xkg(2,n) = DBLE(j-1)/kptparam%nkgrid(2) + DBLE(kptparam%ikgrid(2))/2/kptparam%nkgrid(2)
           xkg(3,n) = DBLE(k-1)/kptparam%nkgrid(3) + DBLE(kptparam%ikgrid(3))/2/kptparam%nkgrid(3)
        end do
     end do
  end do

  !  locate k-points of the uniform grid in the list of irreducible k-points
  !  that was previously calculated

  do nk=1,nkr
     do n=1,kptlist%nk
        xkr = kptlist%kvec(:,n)
        !  xkr is the n-th k-point 
        do i=1,3
           deltap(i) = xkr(i)-xkg(i,nk) - nint (xkr(i)-xkg(i,nk) )
           deltam(i) = xkr(i)+xkg(i,nk) - nint (xkr(i)+xkg(i,nk) )
        end do
        !  deltap is the difference vector, brought back in the first BZ
        !  deltam is the same but with k => -k (for time reversal)
        if ( sqrt ( deltap(1)**2 + &
                    deltap(2)**2 + &
                    deltap(3)**2 ) .lt. eps .or. ( minus_q .and. &
             sqrt ( deltam(1)**2 +  &
                    deltam(2)**2 +  &
                    deltam(3)**2 ) .lt. eps ) ) then
           !  equivalent irreducible k-point found
           equiv(nk) = n
           go to 15
        end if
     end do
     !  equivalent irreducible k-point found - something wrong
     call errore('tetrahedra','cannot locate  k point',nk)
15   continue
  end do

  do n=1,kptlist%nk
     do nk=1,nkr
        if (equiv(nk).eq.n) go to 20
     end do
     !  this failure of the algorithm may indicate that the displaced grid
     !  (with kptparam%ikgrid(1),kptparam%ikgrid(2),kptparam%ikgrid(3).ne.0) does not have the full symmetry of the lattice
     call errore('tetrahedra','cannot remap grid on k-point list',n)
20   continue
  end do

  !  construct tetrahedra

  do i=1,kptparam%nkgrid(1)
     do j=1,kptparam%nkgrid(2)
        do k=1,kptparam%nkgrid(3)
           !  n1-n8 are the indices of k-point 1-8 forming a cube
           ip1 = mod(i,kptparam%nkgrid(1))+1
           jp1 = mod(j,kptparam%nkgrid(2))+1
           kp1 = mod(k,kptparam%nkgrid(3))+1
           n1 = (  k-1) + (  j-1)*kptparam%nkgrid(3) + (  i-1)*kptparam%nkgrid(2)*kptparam%nkgrid(3) + 1
           n2 = (  k-1) + (  j-1)*kptparam%nkgrid(3) + (ip1-1)*kptparam%nkgrid(2)*kptparam%nkgrid(3) + 1
           n3 = (  k-1) + (jp1-1)*kptparam%nkgrid(3) + (  i-1)*kptparam%nkgrid(2)*kptparam%nkgrid(3) + 1
           n4 = (  k-1) + (jp1-1)*kptparam%nkgrid(3) + (ip1-1)*kptparam%nkgrid(2)*kptparam%nkgrid(3) + 1
           n5 = (kp1-1) + (  j-1)*kptparam%nkgrid(3) + (  i-1)*kptparam%nkgrid(2)*kptparam%nkgrid(3) + 1
           n6 = (kp1-1) + (  j-1)*kptparam%nkgrid(3) + (ip1-1)*kptparam%nkgrid(2)*kptparam%nkgrid(3) + 1
           n7 = (kp1-1) + (jp1-1)*kptparam%nkgrid(3) + (  i-1)*kptparam%nkgrid(2)*kptparam%nkgrid(3) + 1
           n8 = (kp1-1) + (jp1-1)*kptparam%nkgrid(3) + (ip1-1)*kptparam%nkgrid(2)*kptparam%nkgrid(3) + 1
           !  there are 6 tetrahedra per cube (and kptparam%nkgrid(1)*kptparam%nkgrid(2)*kptparam%nkgrid(3) cubes)
           n  = 6 * ( (k-1) + (j-1)*kptparam%nkgrid(3) + (i-1)*kptparam%nkgrid(3)*kptparam%nkgrid(2) )

           kpttetra%tetra (1,n+1) = equiv(n1)
           kpttetra%tetra (2,n+1) = equiv(n2)
           kpttetra%tetra (3,n+1) = equiv(n3)
           kpttetra%tetra (4,n+1) = equiv(n6)

           kpttetra%tetra (1,n+2) = equiv(n2)
           kpttetra%tetra (2,n+2) = equiv(n3)
           kpttetra%tetra (3,n+2) = equiv(n4)
           kpttetra%tetra (4,n+2) = equiv(n6)

           kpttetra%tetra (1,n+3) = equiv(n1)
           kpttetra%tetra (2,n+3) = equiv(n3)
           kpttetra%tetra (3,n+3) = equiv(n5)
           kpttetra%tetra (4,n+3) = equiv(n6)

           kpttetra%tetra (1,n+4) = equiv(n3)
           kpttetra%tetra (2,n+4) = equiv(n4)
           kpttetra%tetra (3,n+4) = equiv(n6)
           kpttetra%tetra (4,n+4) = equiv(n8)

           kpttetra%tetra (1,n+5) = equiv(n3)
           kpttetra%tetra (2,n+5) = equiv(n6)
           kpttetra%tetra (3,n+5) = equiv(n7)
           kpttetra%tetra (4,n+5) = equiv(n8)

           kpttetra%tetra (1,n+6) = equiv(n3)
           kpttetra%tetra (2,n+6) = equiv(n5)
           kpttetra%tetra (3,n+6) = equiv(n6)
           kpttetra%tetra (4,n+6) = equiv(n7)
        end do
     end do
  end do

  !  check

  do n=1,kpttetra%ntetra
     do i=1,4
        if ( kpttetra%tetra(i,n).lt.1 .or. kpttetra%tetra(i,n).gt.kptlist%nk ) &
             call errore ('tetrahedra','something wrong',n)
     end do
  end do

  deallocate(equiv)
  deallocate(xkg)

  return
  end subroutine tetrahedra


! ----------------------------------------------------------------------
  subroutine kpt_units( kpt, tpiba )
! ----------------------------------------------------------------------

  type(kpt_type) :: kpt
  real(dp),intent(in) :: tpiba

  if( trim(kpt%param%grid_type) == 'tpiba' ) then
    kpt%list%kvec = kpt%list%kvec * tpiba
  endif

  end subroutine kpt_units


  end module kpt_module
