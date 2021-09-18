! ----------------------------------------------------------------------
  module scalapack_module
! ----------------------------------------------------------------------

! An interface to scalapack routines
! designed to make distributed matrix inversion easier

  use kinds, only : dp
  use io_global, only : stdout
  implicit none
  public

! parameters
  integer,parameter :: DLEN_ = 9
  integer,parameter :: DTYPE_ = 1, &
                       CTXT_  = 2, &
                       M_     = 3, &
                       N_     = 4, &
                       MB_    = 5, &
                       NB_    = 6, &
                       RSRC_  = 7, &
                       CSRC_  = 8, &
                       LLD_   = 9 
! useful variables
  integer :: nprocs, nprow, npcol
  integer :: context
  integer :: iam, myrow, mycol

  integer :: nblock
  integer :: desca(DLEN_)
  
  integer,external :: numroc, ilcm, iceil

  contains


! ----------------------------------------------------------------------
  subroutine scalapack_init
! ----------------------------------------------------------------------

  call blacs_pinfo( iam, nprocs )

  call scalapack_defgrid( nprocs, nprow, npcol )

  call blacs_get( -1, 0, context )
  call blacs_gridinit( context, 'c', nprow, npcol )
  call blacs_gridinfo( context, nprow, npcol, myrow, mycol )

  end subroutine scalapack_init


! ----------------------------------------------------------------------
  subroutine init_scalapack_grid( nprow, npcol, context, comm )
! ----------------------------------------------------------------------

  use parallel_include
  use mp, only : mp_sum, mp_barrier

  integer,intent(in) :: nprow, npcol
  integer,intent(out) :: context
  integer,intent(in),optional :: comm

  integer :: cm, np, me, i, j, ierr 
  integer :: mpime
  integer,allocatable :: usermap(:,:)

  cm=mpi_comm_world
  if( present(comm) ) cm=comm

  call mpi_comm_size( cm, np, ierr )
  call mpi_comm_rank( cm, me, ierr )
  call mpi_comm_rank( MPI_COMM_WORLD, mpime, ierr )

  if( np /= nprow*npcol ) call errore('init_scalapack_grid','process grid dimensions do not match comm size',1)

  ! define the usermap
  allocate( usermap(nprow,npcol) )
  usermap = 0
  i = me / nprow + 1
  j = mod( me, nprow ) + 1
  usermap(j,i) = mpime
  call mp_sum( usermap, cm )

  ! get a fresh context
  call blacs_get( 0, 0, context )

  ! define the subgrid according to usermap
  call blacs_gridmap( context, usermap, nprow, &
                      nprow, npcol )

  write(stdout,*) ' creating a new BLACS context'
  write(stdout,*) ' new context info: ', mpime, me, cm, context
  write(stdout,*) mpime, ' usermap: ', usermap

  deallocate( usermap )

  return

  end subroutine init_scalapack_grid


! ----------------------------------------------------------------------
  subroutine scalapack_distrib( nr, nc, nr_l, nc_l )
! ----------------------------------------------------------------------
  integer,intent(in) :: nr, nc
  integer,intent(out) :: nr_l, nc_l
  integer :: info

  call scalapack_blocksize( nblock, 32, nr, nc, nprow, npcol )
  ! padding constants
  nr_l = numroc( nr, nblock, myrow, 0, nprow )
  nc_l = numroc( nc, nblock, mycol, 0, npcol )
  call descinit( desca, nr, nc, nblock, nblock, 0, 0, context, max(1,nr_l), info )
  
  end subroutine scalapack_distrib


! ----------------------------------------------------------------------
  subroutine scalapack_localindex( i, j, li, lj, islocal )
! ----------------------------------------------------------------------
  integer,intent(in) :: i, j
  integer,intent(out) :: li, lj 
  logical,intent(out) :: islocal
  integer :: prow, pcol

  call infog2l( i, j, desca, nprow, npcol, myrow, mycol, li, lj, prow, pcol )
  islocal = .false.
  if( myrow==prow .AND. mycol==pcol ) islocal = .true.
  
  end subroutine scalapack_localindex


! ----------------------------------------------------------------------
  subroutine scalapack_invert( n, a_l )
! ----------------------------------------------------------------------

  integer :: n
  complex(dp) :: a_l(*)
  integer,allocatable :: ipiv(:), iwork(:)
  complex(dp),allocatable :: work(:)
  integer :: lcm, nr_l, nc_l
  integer,save :: lwork, liwork
  integer :: info
  integer,save :: first_time=0

  ! workspace dimensions
  nr_l = numroc( n, nblock, myrow, 0, nprow )
  nc_l = numroc( n, nblock, mycol, 0, npcol )
!  lwork = nr_l * nblock
!  if( nprow==npcol ) then
!    liwork = nc_l + nblock
!  else
!    lcm = ilcm( nprow, npcol )
!    liwork = nc_l + max( iceil( iceil(nr_l,nblock) , (lcm/nprow) ), nblock )
!  endif

  call blacs_barrier( context, 'All' )

  allocate( ipiv(nr_l+nblock) )
  call pzgetrf( n, n, a_l, 1, 1, desca, ipiv, info )

  ! workspace query
  if( first_time /= n )  then
    lwork=-1
    liwork = -1
    allocate( work(1), iwork(1) )
    call pzgetri( n, a_l, 1, 1, desca, ipiv, &
                  work, lwork, iwork, liwork, info )
    lwork = work(1)
    liwork = iwork(1)
    deallocate( work, iwork )
    first_time = n
  endif

  allocate( work(lwork), iwork(liwork), stat=info )
  call errore('scalapack_invert','unable to allocate workspace',abs(info))

  call pzgetri( n, a_l, 1, 1, desca, ipiv, &
                work, lwork, iwork, liwork, info )

  deallocate( ipiv, work, iwork )

  end subroutine scalapack_invert


! ----------------------------------------------------------------------
  subroutine scalapack_diag( n, a_l, eig, z_l )
! ----------------------------------------------------------------------

  integer :: n
  complex(dp) :: a_l(*), z_l(*)
  real(dp) :: eig(n)
  complex(dp),allocatable :: work(:)
  real(dp),allocatable :: rwork(:)
  integer,allocatable :: iwork(:)
  integer,save :: lwork, lrwork, liwork
  integer :: info
  integer,save :: first_time=0
  integer :: descz(DLEN_)

  integer :: nb, np0, nq0, np, nq

  integer :: ia, ja
  real(dp) :: trace
  complex(dp),external :: pzlatra

  call blacs_barrier( context, 'All' )

  ! check trace
  ia = 1
  ja = 1
  trace = pzlatra( n, a_l, ia, ja, desca )
  write(*,*) ' scalapack_diag: trace (matrix estimate) = ', trace

  ! assume same descriptor for z as a
  descz = desca

  ! workspace query
  if( first_time /= n )  then
!    write(*,*) 'first_time'
    lwork=-1
    lrwork = -1
    liwork = -1
    allocate( work(1), rwork(1), iwork(1) )
    call pzheevd( 'V', 'U', n, a_l, 1, 1, desca, eig, z_l, 1, 1, descz, &
                   work, lwork, rwork, lrwork, iwork, liwork, info )
    lwork = work(1)
    lrwork = rwork(1)
    liwork = iwork(1)
    deallocate( work, rwork, iwork )
!    first_time = n
!    write(*,*) 'first_time', lwork, lrwork, liwork
!    write(*,*) ' n = ', n
!    nb = desca( MB_ )
!    !np0 = numroc( n, nb, 0, 0, nprow )
!    np0 = numroc( max(n,nb,2), nb, 0, 0, nprow )
!    nq0 = numroc( max(n,nb,2), nb, 0, 0, npcol )
!    !lwork = (np0+nq0+nb)*nb+3*n+n**2
!    !lrwork = 2*n + 2*n-2
!    lwork = n + (np0+nq0+nb)*nb
!    np = numroc( n, nb, myrow, iarow, nprow )
!    nq = numroc( n, nb, mycol, iacol, npcol )
!    lrwork = 1+9*n+3*np*nq
!    liwork = 7*n+8*npcol+2
!    write(*,*) ' lwork >= ', lwork
!    write(*,*) ' lrwork >= ', lrwork
  endif

  allocate( work(lwork), rwork(lrwork), iwork(liwork), stat=info )
  call errore('scalapack_diag','unable to allocate workspace',abs(info))

  !write(*,*) ' pzheev '
  !call pzheev( 'V', 'U', n, a_l, 1, 1, desca, eig, z_l, 1, 1, descz, &
  !              work, lwork, rwork, lrwork, info )
  !write(*,*) ' pzheevd '
  call pzheevd( 'V', 'U', n, a_l, 1, 1, desca, eig, z_l, 1, 1, descz, &
                work, lwork, rwork, lrwork, iwork, liwork, info )
  !write(*,*) ' done '

  deallocate( work, rwork, iwork )

  ! check trace
  trace = sum(eig)
  write(stdout,'(a,f16.8)') ' scalapack_diag: trace (eigval estimate) = ', trace
  end subroutine scalapack_diag


! ----------------------------------------------------------------------
  subroutine scalapack_defgrid( nprocs, nprow, npcol )
! ----------------------------------------------------------------------

  integer,intent(in) :: nprocs
  integer,intent(out) :: nprow, npcol
  integer :: sqside

  ! try to make a square grid
  sqside = sqrt( dble(nprocs) )
  nprow = max( sqside, 1 )
  npcol = nprocs / nprow
  do while( nprocs - nprow*npcol > 0 ) 
    nprow=nprow+1
    npcol = nprocs / nprow
  enddo
  
  end subroutine scalapack_defgrid


! ----------------------------------------------------------------------
  subroutine scalapack_blocksize( nb, nbpref, nr, nc, nprow, npcol )
! ----------------------------------------------------------------------
  integer,intent(out) :: nb
  integer,intent(in) :: nbpref, nr, nc, nprow, npcol

  nb = min (nr / nprow, nc / npcol)
  if (nbpref.gt.0) then
     nb = min (nb, nbpref)
  endif
  nb = min (nb, nr, nc)
  if (min(nr,nc).lt.10) nb = 1

  end subroutine scalapack_blocksize


  end module scalapack_module
