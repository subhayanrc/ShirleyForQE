  module hamq_pool

! This module is used to enable a pool distribution of the diagonalization
! and storage workload for hamq_shirley and diag modules
! For a total number of processors nproc, we receive as input nproc_per_pool and
! divide things up accordingly. k-points are distributed across pools, while
! diagonalization and storage are distributed within a pool over the available
! nproc_per_pool processors.

! Hamiltonian matrices are distributed in the ScaLAPACK block-cyclic fashion to
! enable efficient diagonalization using ScaLAPACK routines.
! The projectors are distributed according to atoms, since the non-local
! potential is naturally decomposed into a sum of atomic contributions.

  use parallel_include
  use scalapack_module
  use mp_world, only : mpime, nproc, root, world_comm
  use mp_global, only : mp_sum
  use io_global, only : stdout
  implicit none

! keep things private by default
  private
  public :: hamq_pool_init, &
            read_pool_matrix, read_pool_projs, &
            diag_pool_matrix, hamq_pool_io
  public :: hamq_pool_unpack

  public :: create_cyclic_desc, local_cyclic_dims, &
            cyclic_localindex, &
            set_cyclic_diagonal, &
            create_striped_desc, local_striped_dim, &
            create_striped_desc2, local_striped_dim2, &
            distrib_projs_byatom

  public :: nproc_per_pool, npool, natom_global, nproj_global, &
            index_betaq_global, index_global2local, &
            context_global, context_cyclic, context_striped, &
            mypoolid, mypoolroot, intra_pool_comm, &
            mypool, rootpool, cross_pool_comm, &
            type_atom_global, nbasis_row, nbasis_col, &
            desc_cyclic, desc_striped, desc_striped2, &
            ndiag_cyclic, diag_cyclic, &
            nbasis_stripe, nbasis_stripe_max, nb_stripe

  public :: natomproj_global, index_ldaUq_global, &
            distrib_ldaU_byatom

  interface read_pool_matrix
    module procedure read_pool_matrix_z2, read_pool_matrix_z3
  end interface

  interface diag_pool_matrix
    module procedure diag_pool_matrix_divide_and_conquer
    module procedure diag_pool_matrix_generalized_expert
  end interface

! mpi variables
  integer :: nproc_per_pool, npool
  integer :: mypoolid, mypoolroot, intra_pool_comm
  integer :: mypool, rootpool, cross_pool_comm

! scalapack variables
  integer :: context_global
  integer :: context_cyclic
  integer :: context_striped

  integer :: desc_cyclic(DLEN_)
  integer :: desc_striped(DLEN_)
  integer :: desc_striped2(DLEN_)


! copies of variables that are kept as global versions
  integer :: nbasis_row, nbasis_col
  integer :: natom_global
  integer,allocatable :: type_atom_global(:)
  integer,allocatable :: index_atom_global(:)
  integer :: nproj_global
  integer :: nproj_nl_global
  integer :: nproj_offset, atom_offset
  integer,allocatable :: index_betaq_global(:,:)
  integer,allocatable :: index_global2local(:)
  integer,allocatable :: index_nlproj_betaq_global(:)

  integer :: natomproj_global
  integer,allocatable :: index_ldaUq_global(:,:)

  integer :: ndiag_cyclic
  integer,allocatable :: diag_cyclic(:,:)
  integer :: nbasis_stripe, nbasis_stripe_max, nb_stripe

  integer :: istat(MPI_STATUS_SIZE)


  contains


! ---------------------------------------------------------------------- 
  subroutine hamq_pool_init
! ---------------------------------------------------------------------- 

  use mp, only : mp_bcast, mp_barrier, mp_sum

  character(20) :: cproc_per_pool
  integer :: narg, ierr, id, pid, nproc_per_pool_, npool_
  integer :: ndims, dims(2)
  logical :: periods(2), reorder, remain_dims(2)
  integer :: cartesian_comm
  integer :: nprow, npcol
  integer :: ipool, context, i, j
  integer,allocatable :: usermap(:,:)
#ifdef __PGI
  integer,external :: iargc
#else
  integer,intrinsic :: iargc
#endif

  if( mpime==root ) then
    narg=iargc()
    if( narg /= 1 ) then
      write(stdout,*) ' usage: [exe] nproc_per_pool'
      write(stdout,*) ' assuming nproc_per_pool=1'
      nproc_per_pool=1
    else
      call getarg( 1, cproc_per_pool )
      read(cproc_per_pool,*) nproc_per_pool
    endif
  endif
  call mp_bcast(nproc_per_pool, root, world_comm)

  if( nproc < nproc_per_pool ) nproc_per_pool=nproc

  npool = nproc / nproc_per_pool
  if( mod(nproc,npool) /= 0 ) then
    write(stdout,*) ' error: nproc_per_pool should divide nproc '
    stop
  endif

  if( mpime == root ) then
    write(stdout,*) ' nproc = ', nproc
    write(stdout,*) ' npool = ', npool
    write(stdout,*) ' nproc_per_pool = ', nproc_per_pool
    write(stdout,*)
  endif

  ! new code to organize pools as a Cartesian grid
  ! note that procs will be distributed in C-ordering 
  ndims=2
  dims=(/npool,nproc_per_pool/)
  periods=(/.false.,.false./)
  reorder=.false.
  call mpi_cart_create( mpi_comm_world, ndims, dims, periods, reorder, &
                        cartesian_comm, ierr )

  ! now subdivide this grid to make other communicators

  ! by columns - communicates within a pool only
  ! we use this communicator a lot
  remain_dims=(/.false.,.true./)
  call mpi_cart_sub( cartesian_comm, remain_dims, intra_pool_comm, ierr )

  ! by rows - communicates between pools by rank only - called a slice
  ! this comm is not used so much
  remain_dims=(/.true.,.false./)
  call mpi_cart_sub( cartesian_comm, remain_dims, cross_pool_comm, ierr )
  
  call mpi_comm_size( intra_pool_comm, nproc_per_pool_, ierr )
  if( nproc_per_pool_ /= nproc_per_pool ) &
    call errore('hamq_pool_init','problem with intra_pool_comm',1)
  call mpi_comm_size( cross_pool_comm, npool_, ierr )
  if( npool_ /= npool ) &
    call errore('hamq_pool_init','problem with cross_pool_comm',1)

  call mpi_comm_rank( intra_pool_comm, mypoolid, ierr )
  call mpi_comm_rank( cross_pool_comm, mypool, ierr )

  mypoolroot=0
  rootpool=0

  ipool=mypoolid+mypool*nproc_per_pool
  call mp_sum(ipool,intra_pool_comm)

  if( mpime==root ) then
    write(stdout,*)
    write(stdout,*) ' pool distribution of processes'
    write(stdout,*)
    write(stdout,'(2x,4a18)') 'mpime', 'mypoolid', 'mypool', 'nproc_per_pool'
    call flush(6)
  endif
  call mp_barrier( world_comm )
  do id=0,nproc-1
    if( mpime == id ) then
      write(stdout,'(2x,4i18)') mpime, mypoolid, mypool, nproc_per_pool
      call flush(6)
    endif
    call mp_barrier( world_comm )
  enddo
  if( mpime==0 ) then
    write(stdout,*)
    call flush(6)
  endif
  call mp_barrier( world_comm )


  ! now make ScaLAPACK distributions for each pool
  call blacs_get( -1, 0, context_global )
  call blacs_gridinit( context_global, 'R', 1, nproc )

  ! try to make a square process grid
  ndims=2
  dims=(0,0)
  call mpi_dims_create( nproc_per_pool, ndims, dims, ierr )

  nprow=dims(2)
  npcol=dims(1)

  allocate( usermap(nprow,npcol) )
  usermap = 0
  i = mypoolid / nprow + 1
  j = mod( mypoolid, nprow ) + 1
  usermap(j,i) = mpime
  call mp_sum( usermap, intra_pool_comm )

  ! get a fresh context
  call blacs_get( -1, 0, context )
  ! define the subgrid according to usermap
  call blacs_gridmap( context, usermap, nprow, &
                      nprow, npcol )
  context_cyclic = context
  deallocate( usermap )


  nprow=nproc_per_pool
  npcol=1

  id=0
  allocate( usermap(nprow,npcol) )
  usermap = 0
  i = mypoolid / nprow + 1
  j = mod( mypoolid, nprow ) + 1
  usermap(j,i) = mpime
  call mp_sum( usermap, intra_pool_comm )

  ! get a fresh context
  call blacs_get( -1, 0, context )
  ! define the subgrid according to usermap
  call blacs_gridmap( context, usermap, nprow, &
                      nprow, npcol )
  context_striped = context

  deallocate( usermap )

  end subroutine hamq_pool_init


! ---------------------------------------------------------------------- 
  subroutine hamq_pool_io( iunit )
! ---------------------------------------------------------------------- 

  integer,intent(in) :: iunit
  logical :: io_open
  character(11) :: io_form

  ! I/O
  inquire(iunit,opened=io_open,form=io_form)
  if( io_open ) then
    if( mypoolid /= mypoolroot ) then
      close(iunit,status='delete')
      open(iunit,file='/dev/null',form=io_form)
    else
      write(stdout,*) mpime, ' root of pool ', mypool, ' handling output'
      call flush(6)
    endif
  endif

  end subroutine hamq_pool_io


! ---------------------------------------------------------------------- 
  subroutine create_cyclic_desc( n )
! ---------------------------------------------------------------------- 

  ! distribute nxn in block-cyclic fashion

  use mp, only : mp_barrier

  integer,intent(in) :: n

  integer :: id, ierr
  integer :: nblock
  ! these will be complex double precision arrays - so block sizes of 32
  integer :: nblock_pref=32
  integer :: nprow, npcol, myrow, mycol
  integer :: nr, nc


  call blacs_gridinfo( context_cyclic, nprow, npcol, myrow, mycol )

  ! find block size - prefered size is 32
  call scalapack_blocksize( nblock, nblock_pref, n, n, nprow, npcol )
  write(stdout,*) ' scalapack blocksize ', nblock
  
  nr = numroc( n, nblock, myrow, 0, nprow )
  nc = numroc( n, nblock, mycol, 0, npcol )

  call descinit( desc_cyclic, n, n, nblock, nblock, 0, 0, &
                 context_cyclic, nr, ierr )

  ! since we'll need access to diagonal elements
  call set_cyclic_diagonal

  ! report
  if( mpime==0 ) then
    write(stdout,*) ' nprow = ', nprow, ' npcol = ', npcol
    write(stdout,*) ' memory distribution of matrices n x n'
    write(stdout,'(2x,8a18)') 'mpime', 'mypool', 'mypoolid', 'myrow', 'mycol', &
                         'n', 'nr', 'nc'
    call flush(6)
  endif
  call mp_barrier( world_comm )
  do id=0,nproc-1
    if( mpime == id ) then
      write(stdout,'(2x,8i18)') mpime, mypool, mypoolid, myrow, mycol, &
                           n, nr, nc
      call flush(6)
    endif
    call mp_barrier( world_comm )
  enddo
  if( mpime==0 ) then
    write(stdout,*)
    call flush(6)
  endif
  call mp_barrier( world_comm )

  end subroutine create_cyclic_desc


! ---------------------------------------------------------------------- 
  subroutine local_cyclic_dims( nr, nc )
! ---------------------------------------------------------------------- 

  ! get local dimensions of nxn block-cyclic distribution

  integer,intent(out) :: nr, nc

  integer :: n, nblock
  integer :: nprow, npcol, myrow, mycol


  call blacs_gridinfo( context_cyclic, nprow, npcol, myrow, mycol )
  n=desc_cyclic(M_)
  nblock=desc_cyclic(MB_)
  nr = numroc( n, nblock, myrow, 0, nprow )
  nblock=desc_cyclic(NB_)
  nc = numroc( n, nblock, mycol, 0, npcol )

  end subroutine local_cyclic_dims


! ---------------------------------------------------------------------- 
  subroutine create_striped_desc( n )
! ---------------------------------------------------------------------- 

  ! distribute nxn in striped fashion

  use mp, only : mp_barrier

  integer,intent(in) :: n

  integer :: id, ierr
  integer :: nblock
  integer :: nprow, npcol, myrow, mycol
  integer :: np, nr, nc


  ! find block size
  call mpi_comm_size( intra_pool_comm, np, ierr )
  nblock = max( 0, (n-1)/np ) + 1
  
  call blacs_gridinfo( context_striped, nprow, npcol, myrow, mycol )
  nr = numroc( n, nblock, myrow, 0, nprow )
  nc = numroc( n, n, mycol, 0, npcol )

  call descinit( desc_striped, n, n, nblock, n, 0, 0, &
                 context_striped, max(1,nr), ierr )

  if( mpime==0 ) then
    write(stdout,*)
    write(stdout,*) ' memory distribution of projectors (by nbasis)'
    write(stdout,'(2x,6a18)') 'mpime', 'mypool', 'mypoolid', &
                  'n', 'nr', 'nc'
    call flush(6)
  endif
  call mp_barrier( world_comm )
  do id=0,nproc-1
    if( mpime == id ) then
      write(stdout,'(2x,6i18)') mpime, mypool, mypoolid, &
                           n, nr, nc
      call flush(6)
    endif
    call mp_barrier( world_comm )
  enddo
  if( mpime==0 ) then
    write(stdout,*)
    call flush(6)
  endif
  call mp_barrier( world_comm )

  end subroutine create_striped_desc


! ---------------------------------------------------------------------- 
  subroutine local_striped_dim( nr )
! ---------------------------------------------------------------------- 

  ! get local dimension of nxm striped distribution

  integer,intent(out) :: nr

  integer :: n, nblock
  integer :: nprow, npcol, myrow, mycol


  call blacs_gridinfo( context_striped, nprow, npcol, myrow, mycol )
  n=desc_striped(M_)
  nblock=desc_striped(MB_)
  nr = numroc( n, nblock, myrow, 0, nprow )

  end subroutine local_striped_dim


! ---------------------------------------------------------------------- 
  subroutine create_striped_desc2( n )
! ---------------------------------------------------------------------- 

  ! distribute nxn in striped fashion

  use mp, only : mp_barrier

  integer,intent(in) :: n

  integer :: id, ierr
  integer :: nblock
  integer :: nprow, npcol, myrow, mycol
  integer :: np, nr, nc


  ! find block size
  call mpi_comm_size( intra_pool_comm, np, ierr )
  nblock = max( 0, (n-1)/np ) + 1
  
  call blacs_gridinfo( context_striped, nprow, npcol, myrow, mycol )
  nr = numroc( n, nblock, myrow, 0, nprow )
  nc = numroc( n, n, mycol, 0, npcol )

  call descinit( desc_striped2, n, n, nblock, n, 0, 0, &
                 context_striped, max(1,nr), ierr )

  if( mpime==0 ) then
    write(stdout,*)
    write(stdout,*) ' memory distribution of projectors (by nbasis)'
    write(stdout,'(2x,6a18)') 'mpime', 'mypool', 'mypoolid', &
                  'n', 'nr', 'nc'
    call flush(6)
  endif
  call mp_barrier( world_comm )
  do id=0,nproc-1
    if( mpime == id ) then
      write(stdout,'(2x,6i18)') mpime, mypool, mypoolid, &
                           n, nr, nc
      call flush(6)
    endif
    call mp_barrier( world_comm )
  enddo
  if( mpime==0 ) then
    write(stdout,*)
    call flush(6)
  endif
  call mp_barrier( world_comm )

  end subroutine create_striped_desc2


! ---------------------------------------------------------------------- 
  subroutine local_striped_dim2( nr )
! ---------------------------------------------------------------------- 

  ! get local dimension of nxm striped distribution

  integer,intent(out) :: nr

  integer :: n, nblock
  integer :: nprow, npcol, myrow, mycol


  call blacs_gridinfo( context_striped, nprow, npcol, myrow, mycol )
  n=desc_striped2(M_)
  nblock=desc_striped2(MB_)
  nr = numroc( n, nblock, myrow, 0, nprow )

  end subroutine local_striped_dim2


! ---------------------------------------------------------------------- 
  subroutine set_cyclic_diagonal
! ---------------------------------------------------------------------- 

  integer :: n, i, j, li, lj
  logical :: islocal
  integer :: id

  ! assumes desc_cyclic has already been defined
  n=desc_cyclic(M_)

  if( nproc_per_pool > 1 ) then
    ! define diagonal elements
    ndiag_cyclic=0
    do i=1,n
      call cyclic_localindex( i, i, li, lj, islocal )
      if( islocal ) then
        ndiag_cyclic=ndiag_cyclic+1
      endif
    enddo
    allocate( diag_cyclic(2,ndiag_cyclic) )
    ndiag_cyclic=0
    do i=1,n
      call cyclic_localindex( i, i, li, lj, islocal )
      if( islocal ) then
        ndiag_cyclic=ndiag_cyclic+1
        diag_cyclic(:,ndiag_cyclic) = (/li,lj/)
      endif
    enddo
  else
    ndiag_cyclic=n
    allocate( diag_cyclic(2,n) )
    forall(i=1:n) diag_cyclic(:,i)=(/i,i/)
  endif

!  do id=0,nproc-1
!    if( mpime == id ) then
!      write(*,*)
!      write(*,*) ' mpime = ', mpime
!      write(*,*) ' mypool = ', mypool
!      write(*,*) ' mypoolid = ', mypoolid
!      write(*,*)
!      write(*,*) '   ndiag_cyclic = ', ndiag_cyclic
!      if( ndiag_cyclic>0 ) write(*,*) '    diag_cyclic = ', diag_cyclic
!    endif
!    call mp_barrier()
!  enddo

  end subroutine set_cyclic_diagonal


! ---------------------------------------------------------------------- 
  subroutine hamq_pool_unpack( n, ham, hmat )
! ---------------------------------------------------------------------- 

  use mp_world, only : mpime, root, nproc, world_comm
  use mp, only : mp_max, mp_sum
  use mp_scatt, only: mp_scatter_size, mp_scatter_displ

  integer,intent(in) :: n 
  complex(dp),intent(in) :: ham(:)
  complex(dp),intent(out) :: hmat(:,:)

  complex(dp),parameter :: zero=cmplx(0.d0, 0.d0)

  integer :: nprow, npcol, myrow, mycol
  integer :: ij, i, j, n_l, n_lmx, displ, ierr
  integer :: ip, li, lj, ijl, prow, pcol
  integer,allocatable :: lengths(:), displs(:)


  call blacs_gridinfo( context_cyclic, nprow, npcol, myrow, mycol )

  ! collect info on dimensions
  allocate( lengths(nproc), displs(nproc) )
  call mp_scatter_size( (n*(n+1))/2, n_l, root )
  call mp_scatter_displ( (n*(n+1))/2, displ, root )

  n_lmx=size(ham)
  call mp_max( n_lmx, world_comm )

  displs=0
  displs(mpime+1)=displ
  call mp_sum( displs, world_comm )

  lengths=0
  lengths(mpime+1)=n_l
  call mp_sum( lengths, world_comm )

  ! distribute according to pool
  hmat=zero
  do ip=1,nproc
    i=0 ; j=0
    do ij=1,(n*(n+1))/2
      i=i+1
      if( i>j ) then
        j=j+1; i=1
      endif
      if( ij <= displs(ip) .or. ij > displs(ip)+lengths(ip) ) cycle

      ijl=ij-displs(ip)
      call infog2l( i, j, desc_cyclic, nprow, npcol, &
                    myrow, mycol, li, lj, prow, pcol )
      if( myrow==prow .and. mycol==pcol ) then
        hmat(li,lj) = ham(ijl)
      endif
    enddo
  enddo

  end subroutine hamq_pool_unpack


! ---------------------------------------------------------------------- 
  subroutine read_pool_matrix_z2( unit, mat )
! ---------------------------------------------------------------------- 

  use kinds, only : dp
  use mp, only : mp_bcast
  

  integer,intent(in) :: unit
  complex(dp),intent(out) :: mat(:,:)

  complex(dp),allocatable :: mat_tmp(:,:)
  integer :: nbasis, desc_root(DLEN_), ierr

  nbasis=desc_cyclic(M_)
  if( mypoolid == mypoolroot ) allocate( mat_tmp(nbasis,nbasis) )

  call descinit( desc_root, nbasis, nbasis, nbasis, nbasis, 0, 0, &
                 context_cyclic, max(1,nbasis), ierr )

  ! read on root and send to mypoolroot's
  if( mpime==root ) read(unit) mat_tmp
  if( mypoolid == mypoolroot ) call mp_bcast( mat_tmp, root, cross_pool_comm )

  ! redistribute within pool
  call PZGEMR2D( nbasis, nbasis, mat_tmp, 1, 1, desc_root, &
                 mat, 1, 1, desc_cyclic, desc_cyclic(CTXT_) )

  if( mypoolid == mypoolroot ) deallocate( mat_tmp )

  end subroutine read_pool_matrix_z2


! ---------------------------------------------------------------------- 
  subroutine read_pool_matrix_z3( unit, mat )
! ---------------------------------------------------------------------- 

  use kinds, only : dp
  use mp, only : mp_bcast

  integer,intent(in) :: unit
  complex(dp),intent(out) :: mat(:,:,:)

  complex(dp),allocatable :: mat_tmp(:,:,:)
  integer :: nbasis, desc_root(DLEN_), ierr, i

  nbasis=desc_cyclic(M_)
  if( mypoolid == mypoolroot ) allocate( mat_tmp(nbasis,nbasis,size(mat,3)) )

  call descinit( desc_root, nbasis, nbasis, nbasis, nbasis, 0, 0, &
                 context_cyclic, max(1,nbasis), ierr )

  ! read on root and send to mypoolroot's
  if( mpime==root ) read(unit) mat_tmp
  if( mypoolid == mypoolroot ) call mp_bcast( mat_tmp, root, cross_pool_comm )

  do i=1,size(mat,3)
  ! redistribute within pool
  call PZGEMR2D( nbasis, nbasis, mat_tmp(:,:,i), 1, 1, desc_root, &
                 mat(:,:,i), 1, 1, desc_cyclic, desc_cyclic(CTXT_) )
  enddo

  if( mypoolid == mypoolroot ) deallocate( mat_tmp )

  end subroutine read_pool_matrix_z3


! ---------------------------------------------------------------------- 
  subroutine read_pool_projs( unit, proj, nproj_type, type_atom )
! ---------------------------------------------------------------------- 

  use kinds, only : dp
  use mp, only : mp_bcast

  integer,intent(in) :: unit
  ! proj dimensions are (coefs, projs, basis)
  real(dp),intent(out) :: proj(:,:,:)
  integer,intent(in) :: nproj_type(:), type_atom(:)

  real(dp),allocatable :: proj_tmp(:,:,:), p(:)
  integer :: ierr, i, j, k, j_l, dest
  integer :: iatom, iatom_global
  integer :: it, ip, ip_global, np, np_global, nproj_type_max
  integer :: p_size, id, p_size1, tag

  if( mpime==0 ) write(stdout,*) ' read_pool_projs '

  if( mypoolid == mypoolroot ) then
    allocate( proj_tmp(size(proj,1),nproj_global,size(proj,3)) )
  endif
  nproj_type_max=maxval(nproj_type)
  p_size=size(proj,1)*nproj_type_max*size(proj,3)
  allocate( p(p_size) )

  ! read on root and send to mypoolroot's
  if( mpime==root ) read(unit) proj_tmp
  if( mypoolid == mypoolroot ) call mp_bcast( proj_tmp, root, cross_pool_comm )

  if( mpime==0 ) write(stdout,*) ' root -> mypoolroot '

  ! redistribute within pool
  ip=0
  ip_global=0
  do it=1,size(nproj_type)
  iatom=0
  do iatom_global=1,natom_global
    dest=mod(iatom_global-1,nproc_per_pool)
    if( dest==mypoolid ) then
      iatom=iatom+1
!      it = type_atom(iatom)
      np = nproj_type(it)
    endif
    if( mypoolid==mypoolroot ) then
!      it = type_atom_global(iatom_global)
      np_global = nproj_type(it)
    endif

    if( type_atom_global(iatom_global) == it ) then

!    if( mypoolid==mypoolroot ) write(*,*) ' iatom = ', iatom_global, 'dest = ', dest, ' pool = ', mypool
    if( mypoolid==mypoolroot ) then
      ! if sending to itself - just copy
      if( dest==mypoolid ) then
        proj(:,ip+1:ip+np,:) = proj_tmp(:,ip_global+1:ip_global+np,:)
      ! otherwise package and send to dest
      else
        p_size=size(proj,1)*np_global*size(proj,3)
        p(1:p_size) = reshape( proj_tmp(:,ip_global+1:ip_global+np_global,:), (/ p_size /) )
        call mpi_send( p, p_size, MPI_DOUBLE_PRECISION, dest, &
                       iatom_global-1, intra_pool_comm, ierr )
      endif
    else
    ! so, mypoolid is not mypoolroot
      if( dest==mypoolid ) then
        p_size=size(proj,1)*np*size(proj,3)
        call mpi_recv( p, p_size, MPI_DOUBLE_PRECISION, mypoolroot, &
                       iatom_global-1, intra_pool_comm, istat, ierr )
        proj(:,ip+1:ip+np,:) = reshape( p(1:p_size), (/ size(proj,1), np, size(proj,3) /) )
      endif
    endif

    if( dest==mypoolid ) ip = ip + np
    if( mypoolid==mypoolroot ) ip_global = ip_global + np_global

    endif
  enddo
  enddo

  if( mypoolid == mypoolroot ) then
    deallocate( proj_tmp )
  endif
  deallocate( p )

  end subroutine read_pool_projs

! ----------------------------------------------------------------------
  subroutine cyclic_localindex( i, j, li, lj, islocal )
! ----------------------------------------------------------------------
  integer,intent(in) :: i, j
  integer,intent(out) :: li, lj 
  logical,intent(out) :: islocal
  integer :: nprow, npcol, myrow, mycol
  integer :: prow, pcol

  if( nproc_per_pool > 1 ) then

  call blacs_gridinfo( context_cyclic, nprow, npcol, myrow, mycol )

  call infog2l( i, j, desc_cyclic, nprow, npcol, &
                myrow, mycol, li, lj, prow, pcol )
  islocal = .false.
  if( myrow==prow .AND. mycol==pcol ) islocal = .true.

  else

    li=i
    lj=j
    islocal=.true.

  endif
  
  end subroutine cyclic_localindex


! ----------------------------------------------------------------------
  subroutine distrib_projs_byatom( natom, ntype, type_atom, nproj_type, &
                                   nproj_type_nl, nproj_type_max, nproj, nproj_nl, &
                                   index_betaq, index_nlproj_betaq )
! ----------------------------------------------------------------------

  use mp, only : mp_barrier

!  use hamq_shirley, only : natom, type_atom, nproj_type, nproj_type_nl, &
!                           nproj, nproj_nl 

  integer,intent(inout) :: natom
  integer,intent(in)    :: ntype
  integer,intent(inout) :: type_atom(natom)
  integer,intent(in)    :: nproj_type(ntype), &
                           nproj_type_nl(ntype), nproj_type_max
  integer,intent(inout) :: nproj, nproj_nl, &
                           index_betaq(nproj_type_max,natom), index_nlproj_betaq(nproj_nl)

  integer :: iatom, iatom_l, id, it, ip, i, ierr
  integer :: ikb, ikbnl, ikb_global, ikbnl_global
  logical :: lda_plus_u

  ! copy
  natom_global = natom
  allocate( type_atom_global(natom_global) )
  type_atom_global = type_atom
  nproj_global = nproj
  nproj_nl_global = nproj_nl

  allocate( index_betaq_global(nproj_type_max,natom_global), &
            index_nlproj_betaq_global(nproj_nl_global) )
  index_betaq_global = index_betaq
  index_nlproj_betaq_global = index_nlproj_betaq

  ! if we only have one proc per pool jump out
  ! Yufeng Liang: for general considerations, we still continue (to calculate index_global2local ... )
  ! if( nproc_per_pool == 1 ) return

  nproj = 0
  nproj_nl = 0
  nproj_offset = 0
  natom = 0
  do iatom=1,natom_global
    it    = type_atom_global(iatom)
    ! Calculate the global offset for the projectors in the current pool
    if( mod(iatom-1,nproc_per_pool) < mypoolid ) &
      nproj_offset = nproj_offset + nproj_type(it)

    if( mod(iatom-1,nproc_per_pool) /= mypoolid ) cycle

    natom = natom + 1
    nproj    = nproj    + nproj_type(it)
    nproj_nl = nproj_nl + nproj_type_nl(it)
  enddo

  type_atom=0
  allocate( index_atom_global(natom) )
  ! Yufeng Liang
  ! index_betaq is the distributed version of index_betaq_global
  ! example:
  ! atoms as listed in the pw input file (ATOMIC_POSITIONS):
  ! 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
  ! C, N, O, H, H, N, P, H, C, O, C
  ! order of atom in index_betaq_global is species-major:
  ! 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
  ! C, C, C, H, H, H, N, N, O, O, P
  !
  ! Assume there are 3 pools:
  ! 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1
  ! order of atom in each index_betaq (sorted in each pool):
  !    0, 1, 2, 3
  ! 0: C, H, O, P  
  ! 1: C, H, H, N 
  ! 2: C, N, O
  ! 
  ! If these pools are combined (with MPI_AlltoAll), the order will be:
  ! 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
  ! C, H, O, P, C, H, H, N, C, N, O
  !
  ! Then index_global2local is used to map a given atomic list as in 
  ! ATOMIC_POSITIONS to this particular order. Using species-major order
  ! is not a good idea for two different scfs.

  allocate( index_global2local(nproj_global) )
  index_global2local = 0

  natom = 0
  do iatom=1,natom_global
    if( mod(iatom-1,nproc_per_pool) /= mypoolid ) cycle
    natom = natom + 1
    index_atom_global(natom) = iatom
    type_atom(natom) = type_atom_global(iatom)
  enddo
#ifdef __SHIRLEY_DEBUG
  ! print *, " mypoolid = ", mypoolid, " index_atom_global = ", index_atom_global(1 : natom), " type_atom = ", type_atom(1 : natom)
  write(stdout, *) " nproj_global = ", nproj_global
  print *, " mypoolid = ", mypoolid, " nproj_offset = ", nproj_offset
#endif
  index_betaq=0
  index_nlproj_betaq=0

  ikb=0
  ikbnl=0
  ikb_global=0
  ikbnl_global=0
  do it=1,ntype
    iatom_l=0
    atom_offset = 0
    do iatom=1,natom_global
      ! local atomic index
      if( mod(iatom-1,nproc_per_pool) == mypoolid ) iatom_l=iatom_l+1
      ! is this the right type?
      if( type_atom_global(iatom) == it ) then

        if( mod(iatom-1,nproc_per_pool) == mypoolid ) then
          do ip=1,nproj_type(it)
            index_betaq(ip,iatom_l) = ikb+ip
            ! local index as on this pool
            if(nproj_offset + ikb + ip > nproj_global) then ! test if out of range
              print *, " seg fault: mypoolid = ", mypoolid, " nproj_offset = ", nproj_offset, "ikb = ", ikb, "ip = ", ip
            else
              ! map from species-major order
              ! index_global2local(ikb_global + ip) = nproj_offset + ikb + ip
              ! map from ATOMIC_POSITIONS
              index_global2local(atom_offset + ip) = nproj_offset + ikb + ip
            end if
          enddo

          index_nlproj_betaq(ikbnl+1:ikbnl+nproj_type_nl(it)) = &
            index_nlproj_betaq_global(ikbnl_global+1:ikbnl_global+nproj_type_nl(it)) &
            - (ikb_global-ikb)

          ikb = ikb + nproj_type(it)
          ikbnl = ikbnl + nproj_type_nl(it)
        endif

        ! global index as sorted in species-major order
        ikb_global = ikb_global + nproj_type(it)
        ikbnl_global = ikbnl_global + nproj_type_nl(it)

      endif
      atom_offset = atom_offset + nproj_type( type_atom_global(iatom) )
    enddo ! iatom
  enddo ! it

  ! Yufeng Liang: collect all indices here
  call mp_sum( index_global2local, intra_pool_comm )

  if( mpime==0 ) then
    write(stdout,*)
    write(stdout,*) ' memory distribution of projectors (by atom)'
    write(stdout,'(2x,9a18)') 'mpime', 'mypool', 'mypoolid', &
                  'natom_global', 'natom', 'nproj_global', 'nproj', &
                  'nproj_nl_global', 'nproj_nl'
    call flush(6)
  endif
  call mp_barrier( world_comm )
  do id=0,nproc-1
    if( mpime == id ) then
      write(stdout,'(2x,9i18)') mpime, mypool, mypoolid, &
                           natom_global, natom, nproj_global, nproj, &
                           nproj_nl_global, nproj_nl
      call flush(6)
    endif
    call mp_barrier( world_comm )
  enddo
  if( mpime==0 ) then
    write(stdout,*)
    call flush(6)
  endif
  call mp_barrier( world_comm )

  end subroutine distrib_projs_byatom


! ----------------------------------------------------------------------
  subroutine distrib_ldaU_byatom( ntype, &
                                  natomproj, Hubbard_lmax, Hubbard_l, &
                                  Hubbard_U, Hubbard_alpha, index_ldaUq )
! ----------------------------------------------------------------------

  use mp, only : mp_barrier

!  use hamq_shirley, only : natom, type_atom, nproj_type, nproj_type_nl, &
!                           nproj, nproj_nl 

  integer,intent(in)    :: ntype
  integer,intent(inout) :: natomproj
  integer,intent(in)    :: Hubbard_lmax, Hubbard_l(ntype)
  real(dp),intent(in)   :: Hubbard_U(ntype), Hubbard_alpha(ntype)
  integer,intent(inout) :: index_ldaUq(2*Hubbard_lmax+1,natom_global)

  integer :: iatom, iatom_l, id, it, ip, i, ierr
  integer :: ikb

  ! note that distrib_projs_byatom must have been called first
  ! to define natom_global etc.

  natomproj_global = natomproj
  allocate( index_ldaUq_global(2*Hubbard_lmax+1,natom_global) )
  index_ldaUq_global = index_ldaUq

  ! if we only have one proc per pool jump out
  if( nproc_per_pool == 1 ) return

  natomproj = 0
  do iatom=1,natom_global
    if( mod(iatom-1,nproc_per_pool) /= mypoolid ) cycle

    it        = type_atom_global(iatom)
    if( Hubbard_U(it) /= 0.d0 .or. Hubbard_alpha(it) /= 0.d0 ) then
      natomproj = natomproj + 2*Hubbard_l(it)+1
    endif
  enddo

  index_ldaUq=0

  ikb=0
  do it=1,ntype
    iatom_l=0
    do iatom=1,natom_global
      ! local atomic index
      if( mod(iatom-1,nproc_per_pool) == mypoolid ) iatom_l=iatom_l+1
      ! is this the right type?
      if( type_atom_global(iatom) == it ) then

        if( mod(iatom-1,nproc_per_pool) == mypoolid ) then
          if( Hubbard_U(it) /= 0.d0 .or. Hubbard_alpha(it) /= 0.d0 ) then
            do ip=1,2*Hubbard_l(it)+1
              index_ldaUq(ip,iatom_l) = ikb+ip
            enddo

            ikb = ikb + 2*Hubbard_l(it)+1
          endif
        endif

      endif
    enddo
  enddo

  end subroutine distrib_ldaU_byatom


! ---------------------------------------------------------------------- 
  subroutine diag_pool_matrix_divide_and_conquer( a, eig, z, info )
! ---------------------------------------------------------------------- 

  complex(dp) :: a(:,:), z(:,:)
  real(dp) :: eig(*)
  integer :: info

  integer :: n
  complex(dp),allocatable :: work(:)
  real(dp),allocatable :: rwork(:)
  integer,allocatable :: iwork(:)
  integer,save :: lwork, lrwork, liwork
  integer,save :: work_space=0

  integer :: nb, np0, nq0, np, nq

  call blacs_barrier( context_cyclic, 'All' )

  n=desc_cyclic(M_)

  ! workspace query
  if( work_space /= n )  then
    !write(*,*) ' pzheevd workspace'
    lwork=-1
    lrwork = -1
    liwork = -1
    allocate( work(1), rwork(1), iwork(1) )
    call pzheevd( 'V', 'U', n, a, 1, 1, desc_cyclic, eig, z, 1, 1, desc_cyclic, &
                  work, lwork, rwork, lrwork, iwork, liwork, info )
    lwork = work(1)
    lrwork = rwork(1)
    liwork = iwork(1)
    deallocate( work, rwork, iwork )
    work_space=n
    !write(*,*) ' lwork ', lwork, ' lrwork ', lrwork, ' liwork ', liwork
    !write(*,*) ' pzheevd workspace done'
  endif

  allocate( work(lwork), rwork(lrwork), iwork(liwork), stat=info )
  !if( info /= 0 ) then
  !  write(stdout,*) lwork, lrwork, liwork
  !endif
  call errore('diag_pool_matrix_divide_and_conquer','unable to allocate workspace',abs(info))

  if( mpime==root ) write(stdout,*) ' pzheevd diag'
  call pzheevd( 'V', 'U', n, a, 1, 1, desc_cyclic, eig, z, 1, 1, desc_cyclic, &
                work, lwork, rwork, lrwork, iwork, liwork, info )
  if( mpime==root ) write(stdout,*) ' pzheevd done'

  deallocate( work, rwork, iwork )

  end subroutine diag_pool_matrix_divide_and_conquer


! ---------------------------------------------------------------------- 
  subroutine diag_pool_matrix_generalized_expert( a, b, eig, z, info, l_eig, u_eig, neig_found )
! ---------------------------------------------------------------------- 

  use mp, only : mp_bcast

  integer :: ierr
  complex(dp) :: a(:,:), b(:,:), z(:,:)
  real(dp) :: eig(:)
  integer :: info
  integer,optional :: l_eig, u_eig, neig_found

  integer :: n
  real(dp) :: vl, vu, abstol, orfac
  integer :: il, iu, n_found, nv_found
  character :: eigval_range
  complex(dp),allocatable :: a_copy(:,:), b_copy(:,:)
  complex(dp),allocatable :: work(:)
  real(dp),allocatable :: rwork(:)
  integer,allocatable :: iwork(:)
  integer,save :: lwork, lrwork, liwork
  integer,allocatable :: ifail(:), icluster(:)
  real(dp),allocatable :: gap(:)
  logical,save :: first_time=.true.
  integer,save :: work_space=0

  integer :: nprow, npcol, myrow, mycol

  integer :: nb, np0, nq0, np, nq
  ! Added a save here to keep track of how much extra temp space is required
  ! beyond the crappy (under)estimate provided by PZHEGVX (clearly broken)
  integer,save :: memfac
  integer :: i, j

  call blacs_gridinfo( context_cyclic, nprow, npcol, myrow, mycol )

  n=desc_cyclic(M_)

  allocate( ifail(n), icluster(2*nprow*npcol), gap(nprow*npcol) )

  if( mypoolid == mypoolroot ) then
  if( first_time ) then
    memfac=max(2,nproc_per_pool/2)
    first_time=.false.
  endif
  endif
  call mp_bcast( memfac, mypoolroot, intra_pool_comm )

  ! keep copies of a,b
  allocate( a_copy(size(a,1),size(a,2)), b_copy(size(b,1),size(b,2)), stat=ierr )
  if( ierr/=0 ) call errore('diag_pool_matrix_generalized_expert', &
    'unable to allocate space for copies of matrices',abs(ierr))
  a_copy = a
  b_copy = b

  call pool_diag
  do while( info==2 )
    write(stdout,*) ' insufficient workspace'
    memfac=memfac+1
    work_space=0
    a=a_copy
    b=b_copy
    write(stdout,*) ' increasing memfac to ', memfac
    call pool_diag
  enddo

  neig_found = n_found

  deallocate( a_copy, b_copy )

  deallocate( ifail, icluster, gap )

  contains

    subroutine pool_diag

    call blacs_barrier( context_cyclic, 'All' )

    abstol=-1.d0
    orfac=-1.d0 ! sets a default value of 1.d-3 for reorthogonalization
    il=1
    iu=n
    eigval_range='A'
    if( present(l_eig) ) il=l_eig
    if( present(u_eig) ) iu=u_eig
    if( il /= 1 .or. iu /= n ) eigval_range='I'

    ! workspace query
    if( work_space /= n )  then
      lwork=-1
      lrwork = -1
      liwork = -1
      allocate( work(1), rwork(3), iwork(1) )
      call PZHEGVX( 1, 'V', eigval_range, 'U', n, a, 1, 1, &
                    desc_cyclic, b, 1, 1, desc_cyclic, vl, vu, il, iu, &
                    abstol, n_found, nv_found, eig, orfac, z, 1, 1, desc_cyclic, &
                    work, lwork, rwork, lrwork, iwork, liwork, &
                    ifail, icluster, gap, info )

      lwork = work(1)
      lrwork = rwork(1)
      liwork = iwork(1)
      !write(stdout,*) ' lwork ', work(1), ' lrwork ', rwork(1), ' liwork ', iwork(1)
      deallocate( work, rwork, iwork )
      work_space=n

      lwork=lwork*memfac
      lrwork=lrwork*memfac
      liwork=liwork*memfac
      if( mpime==root ) write(stdout,*) ' including extra factor of', memfac, ' in workspace'
      !write(stdout,*) ' lwork ', lwork, ' lrwork ', lrwork, ' liwork ', liwork
    endif

    allocate( work(lwork), rwork(max(3,lrwork)), iwork(liwork), stat=info )
    call errore('diag_pool_matrix_generalized_expert','unable to allocate workspace',abs(info))

    if( mpime==root ) write(stdout,*) ' pzhegvx diag'
    call PZHEGVX( 1, 'V', eigval_range, 'U', n, a, 1, 1, &
                  desc_cyclic, b, 1, 1, desc_cyclic, vl, vu, il, iu, &
                  abstol, n_found, nv_found, eig, orfac, z, 1, 1, desc_cyclic, &
                  work, lwork, rwork, lrwork, iwork, liwork, &
                  ifail, icluster, gap, info )
    if( mpime==root ) write(stdout,*) ' pzhegvx done'

    deallocate( work, rwork, iwork )

    call blacs_barrier( context_cyclic, 'All' )

    if( info /= 0 ) then
      write(stdout,*) ' PZHEGVX: info = ', info
      if( mod(info/16,2)/=0 ) write(stdout,*) ' overlap matrix not positive definite'

!      write(*,*) ' ifail = ', ifail
!      write(*,*) ' icluster = ', icluster
!      write(*,*) ' gap = ', gap
    endif

    end subroutine pool_diag

  end subroutine diag_pool_matrix_generalized_expert


  end module hamq_pool
