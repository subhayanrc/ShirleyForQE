  program shirley_qdiag

  ! stand-alone utility to read in the Hamiltonian in the
  ! optimal Shirley basis and solve it for a given input
  ! q-point list

  ! David Prendergast, UCB, Nov 2007

  ! now parallelized
  use kinds, only : dp
  use parallel_include
  use hamq_shirley
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_world, only : nproc, mpime, root, world_comm, mp_world_end
  use mp, only : mp_bcast, mp_sum
  use mpio
  use kpt_module
  use corerepair_module
  use shirley_input_module
  use diag_module
  use hamq_pool, only : nproc_per_pool, npool, &
                        mypool, rootpool, mypoolid, mypoolroot, &
                        cross_pool_comm, intra_pool_comm, &
                        desc_cyclic, desc_striped, context_cyclic, &
                        local_striped_dim, cyclic_localindex, &
                        context_global, nproj_global, &
                        natom_global, type_atom_global
!  use hamq_pool, only : nbasis_global

  implicit none

  integer,external :: freeunit

  REAL(DP), PARAMETER :: rytoev=13.6058d0
  complex(dp),parameter :: ONE=(1.d0,0.d0)
  complex(dp),parameter :: ZERO=(0.d0,0.d0)

  character(255) :: eigval_file, eigvec_file, info_file, proj_file
  logical :: ex

  integer(kind=MPI_OFFSET_KIND) :: offset
  integer :: fheigval, fheigvec, fhproj
  integer :: status(MPI_STATUS_SIZE)
  integer :: ierr

  character(255) :: fmtstr

  real(dp) :: lambda

  real(dp) :: kvec(3)
  real(dp) :: qpathlen, cqvec(3), dqvec(3)
  character(255) :: ic

  integer :: ik
  integer :: i,j,k
  integer :: nbasis_subset

  integer :: nk, nbnd, ispin
  real(dp) :: nelec_, alat, volume, at(3,3), bg(3,3), tpiba
  integer :: nspin
  real(dp) :: fermi_energy
  complex(dp),allocatable :: ztmp(:,:)

  integer :: nbasis_stripe, stripe_displ
  integer,allocatable :: stripe_sizes(:)
  complex(dp),allocatable :: eigvec_s(:,:)

  complex(dp),allocatable :: beta_nk(:,:)
  integer :: iatom_global, it, dest, offbeta

  integer :: iuninf
  namelist /info/ nk, nbnd, nelec, alat, volume, &
                  at, bg, tpiba, fermi_energy, nspin, lda_plus_u

  call shirley_input

  if( band_subset(1) < 1 ) band_subset(1)=1
!  if( band_subset(2) < band_subset(1) .or. band_subset(2) > nbasis ) band_subset(2)=nbasis_global
  if( band_subset(2) < band_subset(1) .or. band_subset(2) > nbasis ) band_subset(2)=nbasis
  nbasis_subset = band_subset(2)-band_subset(1)+1
  write(stdout,*) ' band_subset = ', band_subset

  call diagx_init( band_subset(1), band_subset(2) )

  call dump_system( nelec_, alat, volume, at, bg, tpiba, nspin, lda_plus_u )
  write(stdout,*) ' nspin = ', nspin
  write(stdout,*) ' lda_plus_u = ', lda_plus_u

  ! band structure is given in units of tpiba
  if( trim(kpt%param%grid_type) == 'bandstructure' ) then
    kpt%list%kvec = kpt%list%kvec*tpiba
    kpt%bandstructure%kpathlen = kpt%bandstructure%kpathlen*tpiba
    kpt%bandstructure%kpathlensp = kpt%bandstructure%kpathlensp*tpiba
    kpt%bandstructure%xksp = kpt%bandstructure%xksp*tpiba
  endif
 
  ! store total number of k-points
  nk=kpt%list%nk

  if( mpime==root ) then
    info_file=trim(outfile)//'.info'
    nbnd=nbasis_subset  ! important difference

    fermi_energy = efermi

    write(stdout,*) ' Fermi energy = ', fermi_energy

    iuninf=freeunit()
    open(iuninf,file=info_file,form='formatted')
    write(iuninf,nml=info)
    write(iuninf,*) kpt%param%cartesian
    write(iuninf,*) trim(kpt%param%grid_type)
    write(iuninf,*) kpt%list%wk
    write(iuninf,*) kpt%list%kvec
    if( trim(kpt%param%grid_type) == 'tetrahedra' ) then
      write(iuninf,*) kpt%tetra%ntetra
      write(iuninf,*) kpt%tetra%tetra
    else if( trim(kpt%param%grid_type) == 'bandstructure' ) then
      write(iuninf,*) kpt%bandstructure%nksp
      do i=1,kpt%bandstructure%nksp
        write(iuninf,*) kpt%bandstructure%kpathlensp(i), &
                        trim(kpt%bandstructure%labelsp(i))
      enddo
      write(iuninf,*) kpt%bandstructure%kpathlen
    endif
    close(iuninf)

    write(stdout,*) info_file
    
  endif

  ! MPI-IO
  if( mypoolid==mypoolroot ) then
    eigval_file=trim(outfile)//'.eigval'
    inquire(file=trim(eigval_file),exist=ex)
    if( ex ) then
      ! delete pre-existing files
      fheigval=freeunit()
      open(fheigval,file=trim(eigval_file),form='unformatted')
      close(fheigval,status='delete')
    endif

    call mp_file_open_dp( eigval_file, fheigval, rootpool, cross_pool_comm )
  endif
  ! open eigvec file on all pool procs
  if( eigvec_output ) then
    eigvec_file=trim(outfile)//'.eigvec'
    call mp_file_open_dp( eigvec_file, fheigvec, rootpool, intra_pool_comm )
  endif
  ! open proj file on all pool procs
  if( proj_output ) then
    proj_file=trim(outfile)//'.proj'
    call mp_file_open_dp( proj_file, fhproj, rootpool, intra_pool_comm )
  endif

  do ispin=1,nspin

  do ik=1,kpt%list%nk
! ======================================================================

    if( mod(ik-1,npool)/=mypool ) cycle

    write(stdout,'(a,3f12.5,i6,a,i6,a,i3)') ' k-point ', &
      kpt%list%kvec(1:3,ik), ik, ' of ', kpt%list%nk, &
                                  ' on node ', mpime

    ! build the Hamiltonian for this q-point
    call diag_build_hamk( kpt%list%kvec(1:3,ik), kpt%param%cartesian, ispin )

    if( kinetic_only ) then
      allocate( ztmp(nbasis,nbasis) )
      call diag_build_kink( kpt%list%kvec(1:3,ik), kpt%param%cartesian, ztmp )
      forall( i=1:nbasis ) eigval(i)=ztmp(i,i)
      deallocate( ztmp )
    else if( local_only ) then
      allocate( ztmp(nbasis,nbasis) )
      call diag_build_vlock( kpt%list%kvec(1:3,ik), kpt%param%cartesian, ispin, ztmp )
      forall( i=1:nbasis ) eigval(i)=ztmp(i,i)
      deallocate( ztmp )
    else if( nonlocal_only ) then
      allocate( ztmp(nbasis,nbasis) )
      call diag_build_vnlk( kpt%list%kvec(1:3,ik), kpt%param%cartesian, ispin, ztmp )
      forall( i=1:nbasis ) eigval(i)=ztmp(i,i)
      deallocate( ztmp )
    else if( smatrix_only ) then
      allocate( ztmp(nbasis,nbasis) )
      call diag_build_sk( kpt%list%kvec(1:3,ik), kpt%param%cartesian, ztmp )
      forall( i=1:nbasis ) eigval(i)=ztmp(i,i)
      deallocate( ztmp )
    else

    call diagx_ham

    endif

    write(stdout,*) ik, eigval(band_subset(1):band_subset(2))*rytoev

    if( mypoolid==mypoolroot ) then
      ! write eigenvalues
      offset = ((ispin-1)*kpt%list%nk + ik-1)*nbasis_subset
      call mpi_file_write_at( fheigval, offset, &
                              eigval(band_subset(1):band_subset(2)), nbasis_subset, &
                              MPI_DOUBLE_PRECISION, status, ierr )
    endif

    ! dump eigenvectors
    if( eigvec_output ) then
      if( nproc_per_pool > 1 ) then
        call local_striped_dim( nbasis_stripe )
        allocate( eigvec_s(nbasis_stripe,nbasis) )

        ! redistribute eigenbasis
        call PZGEMR2D(nbasis,nbasis, &
                      eigvec,1,1,desc_cyclic, &
                      eigvec_s,1,1,desc_striped,context_cyclic)
                   
        allocate( stripe_sizes(nproc_per_pool) )
        stripe_sizes = 0
        stripe_sizes(mypoolid+1) = nbasis_stripe
        call mp_sum( stripe_sizes, intra_pool_comm )
        stripe_displ=0
        do i=1,mypoolid,1
          stripe_displ=stripe_displ+stripe_sizes(i)
        enddo
        deallocate( stripe_sizes )
        ! write eigenvectors
        offset = ((ispin-1)*kpt%list%nk + ik-1)*(nbasis_stripe*nbasis)*2
        offset = offset + stripe_displ
        call mpi_file_write_at( fheigvec, offset, &
                                eigvec_s, nbasis_stripe*nbasis*2, &
                                MPI_DOUBLE_PRECISION, status, ierr )
        deallocate( eigvec_s )
      else
        ! write eigenvectors
        offset = ((ispin-1)*kpt%list%nk + ik-1)*(nbasis*nbasis)*2
        call mpi_file_write_at( fheigvec, offset, &
                                eigvec, (nbasis*nbasis)*2, &
                                MPI_DOUBLE_PRECISION, status, ierr )
      endif
    endif

    ! construct <nk|beta> for all
    if( proj_output ) then
      allocate( beta_nk(nproj,nbasis) )
      if( nproc_per_pool > 1 ) then

        write(mpime+100,*) ' betaq and nproj = ', size(betaq,1), nproj_global, nproj
        ! expand projectors in the eigenbasis
        CALL ZGEMM( 'N', 'N', nproj, nbasis, nbasis, one, &
                    betaq, nproj, &
                    eigvec, nbasis, &
                    zero, beta_nk, nproj )

        offset = ((ispin-1)*kpt%list%nk + ik-1)*nproj_global*nbasis*2
        offbeta = 0
        do iatom_global=1,natom_global
          dest=mod(iatom_global-1,nproc_per_pool)
          it=type_atom_global(iatom_global)
          if( dest==mypoolid ) then
            ! write projectors
            call mpi_file_write_at( fhproj, offset, &
                                    beta_nk(1+offbeta,1), nproj_type(it)*nbasis*2, &
                                    MPI_DOUBLE_PRECISION, status, ierr )
            offbeta=offbeta+nproj_type(it)
          endif
          offset = offset + nproj_type(it)
        enddo

      else

        ! transform all projectors to eigenbasis for this k
        CALL ZGEMM( 'N', 'N', nproj, nbasis, nbasis, one, &
                    betaq, nproj, &
                    eigvec, nbasis, &
                    zero, beta_nk, nproj )
        ! write projectors
        offset = ((ispin-1)*kpt%list%nk + ik-1)*nproj*nbasis*2
        call mpi_file_write_at( fhproj, offset, &
                                beta_nk, nproj*nbasis*2, &
                                MPI_DOUBLE_PRECISION, status, ierr )
      endif
      deallocate( beta_nk )
    endif

  enddo ! ik

  enddo ! ispin

  ! close MPI-IO files
  if( mypoolid==mypoolroot ) then
    ! close binary dump
    call mpi_file_close( fheigval, ierr )
    if( ierr/=0 ) &
      call errore('shirley_qdiagp','problem closing eigval file',abs(ierr))

    if( eigvec_output ) then
      ! close eigvec file
      call mpi_file_close( fheigvec, ierr )
      if( ierr/=0 ) &
        call errore('shirley_qdiagp','problem closing eigvec file',abs(ierr))
    endif
  endif

  write(stdout,*) ' end shirley_qdiag'
  call mp_world_end( )
  
  end program shirley_qdiag
