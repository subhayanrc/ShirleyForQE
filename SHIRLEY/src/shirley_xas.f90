program shirley_xas

  ! stand-alone utility to read in the Hamiltonian in the
  ! optimal Shirley basis and produce x-ray absorption spectra
  ! by averaging over a k-point grid

  ! David Prendergast, UCB, Jul 2007

  use kinds, only : dp
  use parallel_include
  use hamq_shirley
!  use diag_shirley
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_world, only : nproc, mpime, root, world_comm, mp_world_end
  use mp, only : mp_bcast, mp_barrier, mp_sum
  use mp_scatt, only : mp_scatter_size, mp_scatter
  use mpio
  use scalapack_module
  use kpt_module
  use corerepair_module
  use shirley_input_module, nelec_in=>nelec
  use diag_module
  use hamq_pool, only : nproc_per_pool, npool, &
                        mypool, rootpool, mypoolid, mypoolroot, &
                        cross_pool_comm, intra_pool_comm, &
                        desc_cyclic, desc_striped, context_cyclic, &
                        local_striped_dim, cyclic_localindex, &
                        create_striped_desc2, local_striped_dim2, &
                        context_global, nproj_global, &
                        natom_global, type_atom_global, index_global2local

  implicit none

  REAL(DP), PARAMETER :: rytoev=13.6058d0
  complex(dp),parameter :: ONE=(1.d0,0.d0)
  complex(dp),parameter :: ZERO=(0.d0,0.d0)


  character(255) :: eigval_file, info_file, haminf_file, xmat_file, eigvec_file, proj_file
  character(255) :: str_err

  integer(kind=MPI_OFFSET_KIND) :: offset, offset_sk
  integer :: fh, fheigval, fheigvec, fhxmat, fhproj
  integer :: status(MPI_STATUS_SIZE)
  integer :: ierr, len_err, ierr_str

  integer :: group_ranks(1), group_size, ionode_group, ionode_comm
  integer :: wcomm, wgroup
  integer :: nbasis_stripe
  integer :: nbasis_stripe2

  complex(dp),allocatable :: ztmp(:,:), ztmp_local(:,:)
  complex(dp),allocatable :: beta_nk(:,:), beta_nk_c(:,:), betaq_send(:,:), betaq_s(:,:), betaq_ss(:,:)
  real(dp) :: fermi_energy

  complex(dp),allocatable :: posn(:,:,:), posn_l(:,:)

  integer :: iuninf

  integer,external :: freeunit
  integer,external :: ilaenv

  logical :: islocal

  integer :: ik, i, j, nbnd, nk, ispin
  integer :: icore, ibasis, ncp, ncm, ixyz, core_root, icore_atom
  integer :: iatom, it, np, npm, iatom_global
  integer :: dest, offbeta

  complex(dp), allocatable :: eigvec_s(:, :)
  ! local # of basis and projectors
  integer, allocatable :: stripe_sizes(:), stripe_displs(:), proj_sizes(:), proj_displs(:)
  integer :: stripe_displ, stripe_lo, stripe_hi, proj_displ

  integer :: desc1(DLEN_), desc2(DLEN_)

  real(dp) :: nelec, alat, volume, at(3,3), bg(3,3), tpiba
  integer :: nspin

  namelist /info/ nk, nbnd, ncp, nelec, nbasis, alat, volume, &
                  at, bg, tpiba, fermi_energy, nspin, lda_plus_u


  call shirley_input

  ! ======================================================================
  ! print out header and info

  call dump_system( nelec, alat, volume, at, bg, tpiba, nspin, lda_plus_u )
  write(stdout,*)
  write(stdout,*) ' System details: '
  write(stdout,*) '       nelec = ', nelec
  write(stdout,*) '        alat = ', alat
  write(stdout,*) ' cell volume = ', volume
  write(stdout,*) '       nspin = ', nspin
  write(stdout,*)

  ! fix band_subset if undefined
  if( band_subset(1) < 1 ) band_subset(1) = 1
  if( band_subset(2) < band_subset(1) .or. &
      band_subset(2) > nbasis ) band_subset(2) = nbasis

  if( any( band_subset > 0 ) ) then
    write(stdout,*) ' band_subset = ', band_subset
    call diagx_init( band_subset(1), band_subset(2) )
  else
    call diag_init
  endif

  write(stdout,*) ' shirley_xas'
  write(stdout,*)

  ! store total number of k-points
  nk = kpt % list % nk

  ! output information
  if( mpime == root ) then

    info_file = trim(outfile)//'.info'

    if( any( band_subset > 0 ) ) then
      nbnd = band_subset(2) - band_subset(1) + 1
    else
      nbnd = nbasis
    endif

    fermi_energy = efermi

    write(stdout,*) ' Fermi energy = ', fermi_energy

    iuninf = freeunit()
    open( iuninf, file = info_file, form = 'formatted' )

    ncp = 1 ! ncp is undefined here. Need to consider multiple core levels ?
    write(iuninf, nml = info)
    write(iuninf, *) kpt % param % cartesian
    write(iuninf, *) trim(kpt % param % grid_type)
    write(iuninf, *) kpt % list % wk
    write(iuninf, *) kpt % list % kvec
    if( trim( kpt % param % grid_type ) == 'tetrahedra' ) then
      write(iuninf, *) kpt % tetra % ntetra
      write(iuninf, *) kpt % tetra % tetra
    endif

    close(iuninf)
    
    haminf_file=trim(outfile)//'.haminf'

    fhhaminf=freeunit()
    open(fhhaminf,file=trim(haminf_file),form='formatted')

    call write_haminf()

    close(fhhaminf)

  endif

  call mp_bcast( nbnd, root, world_comm )

  ! make another descriptor for matrix multiplies determined by splitting
  ! data across those bands we want
  call create_striped_desc2( nbnd )

  ! ======================================================================
  ! MPI-IO

  if ( mypoolid == mypoolroot ) then

    eigval_file = trim(outfile)//'.eigval'
    call mp_file_open_dp( eigval_file, fheigval, rootpool, cross_pool_comm )

  endif

  if ( corerep % ncore > 0 ) then

    xmat_file = trim(outfile)//'.xmat'
    call mp_file_open_dp( xmat_file, fhxmat, rootpool, cross_pool_comm )
    write(stdout, *) ' corerep % ncore = ', corerep % ncore
  endif

  ! all pool processes involved in writing eigvec
  if( eigvec_output ) then

    eigvec_file = trim(outfile)//'.eigvec'
    call mp_file_open_dp( eigvec_file, fheigvec, rootpool, intra_pool_comm )

  endif

  if( proj_output ) then

    proj_file=trim(outfile)//'.proj'
    call mp_file_open_dp( proj_file, fhproj, rootpool, intra_pool_comm )

  endif

  ! ======================================================================

  ! allocations for variables

  allocate( beta_nk(nproj_global, nbnd) )
 
  stripe_displ = 0
  proj_displ = 0
  nbasis_stripe = nbasis
  nbasis_stripe2 = nbasis
  ! find out the stripes: in future this block should be rewritten in OOP
  if( nproc_per_pool > 1 ) then

    ! split up nbasis (nbnd)
    call local_striped_dim( nbasis_stripe )
    call local_striped_dim2( nbasis_stripe2 )


    ! find stripe sizes and distribute
    allocate( stripe_sizes(nproc_per_pool), stripe_displs(nproc_per_pool) )
    stripe_sizes  = 0
    stripe_displs = 0
  
    stripe_sizes(mypoolid + 1) = nbasis_stripe
    call mp_sum( stripe_sizes, intra_pool_comm )
  
    ! calculate the current stripe displacement
    do i = 1, mypoolid
      stripe_displ = stripe_displ + stripe_sizes(i)
    enddo

    stripe_displs(mypoolid + 1) = stripe_displ
    call mp_sum( stripe_displs, intra_pool_comm )

    ! split up nproj_global
    allocate( proj_sizes(nproc_per_pool), proj_displs(nproc_per_pool) )
    proj_sizes = 0
    proj_displs = 0
    
    proj_sizes(mypoolid + 1) = nproj
    call mp_sum( proj_sizes, intra_pool_comm )
    
    do i = 1, mypoolid
      proj_displ = proj_displ + proj_sizes(i)
    enddo

    proj_displs(mypoolid + 1) = proj_displ
    call mp_sum( proj_displs, intra_pool_comm )

  else
    allocate( proj_sizes(1), proj_displs(1) )
    proj_sizes = 0
    proj_displs = 0
    proj_sizes(1) = nproj
  end if

  allocate( eigvec_s(nbasis_stripe, nbasis) ) 

#ifdef __SHIRLEY_DEBUG
  ! print *, " mypoolid = ", mypoolid, " index_betaq = ", index_betaq
#endif

! ======================================================================
  do ispin = 1, nspin
  do ik = 1, kpt % list % nk
! ======================================================================

    if ( mod(ik - 1, npool) /= mypool ) cycle

    write(stdout, '(a,3f12.5,a,f12.5,i6,a,i6,a,a)') &
    ' k-point ', kpt%list%kvec(1:3,ik), &
     ' weight ',       kpt%list%wk(ik), ik,   &
      ' of ', kpt%list%nk, ' on node ', trim(nodenumber)

    ! calculate the common offset due to spins and k-points
    offset_sk =  ( ispin - 1 ) * (kpt % list % nk) + ik - 1 

    ! build the Hamiltonian for this q-point
    ! local contribution

    call diag_build_hamk( kpt % list % kvec(1 : 3, ik), kpt % param % cartesian, ispin )

    call mp_barrier( intra_pool_comm )

    ! diagonalize matrix to find eigenvalues (eigval) and vectors (eigvec)

    if( any( band_subset > 0 ) ) then

      call diagx_ham

    else

      call diag_ham

    endif

    call mp_bcast( neig_found, mypoolroot, intra_pool_comm )

    if ( neig_found /= nbnd ) then

      write(stdout, *) ' only found ', neig_found, ' eigenvalues of a total ', nbnd
      call errore('shirley_xas','problem diagonalizing hamiltonian',1)

    endif

    ! ======================================================================

    ! write eigenvalues: eigval
    if( mypoolid == mypoolroot ) then
      offset = offset_sk * neig_found
      call mpi_file_write_at( fheigval, offset, &
                               eigval, neig_found, &
                               MPI_DOUBLE_PRECISION, status, ierr )

#ifdef __SHIRLEY_DEBUG
      ! write the eigenvalue for ik = 0
      if ( ik == 1 ) then
        write(stdout, *) 'Output eigval to eigval.dat for debug. '
        open(unit = fh, file = 'eigval.dat', form = 'formatted')
        write(fh, *) eigval
        close(fh)
      endif
#endif
    endif ! write eigenvalues

    ! ======================================================================

    ! redistribute eigenbasis
    if( nproc_per_pool > 1 ) then

        call PZGEMR2D(nbasis, nbasis, &
                      eigvec, 1, 1, desc_cyclic, &
                      eigvec_s, 1, 1, desc_striped, context_cyclic)

    else
        eigvec_s = eigvec
    end if
 
    ! write eigenvectors: eigvec, or < B_i | nk >
    if ( eigvec_output ) then

        offset = (offset_sk * nbasis +  stripe_displ) * nbnd * 2

        ! write TRANSPOSE(eigvec_s)
        call mpi_file_write_at( fheigvec, offset, &
                                transpose(eigvec_s), nbasis_stripe * nbnd * 2, &
                                MPI_DOUBLE_PRECISION, status, ierr )

#ifdef __SHIRLEY_DEBUG
        if ( mypoolid == mypoolroot .and. ik == 1 ) then          
          write(stdout, *) 'Output eigvec_s to eigvec.dat for debug. '
          open(unit = fh, file = 'eigvec.dat', form = 'formatted')
          write(fh, *) eigvec_s(1, :)
          close(fh)
        endif
#endif
    endif ! write eigenvectors

    ! ======================================================================
    ! Process betaq

    ! betaq
    ! < beta_l | B_i > (nproj * nbasis) = 
    ! < beta_1 | B_1 >   < beta_1 | B_2 >    ...
    ! < beta_2 | B_1 >   < beta_2 | B_2 >    ...
    ! ...

    if ( nproc_per_pool > 1 ) then
      ! Convert betaq from row stripes into column stripes
      ! so as to find beta_nk

      ! betaq:                 betaq_col:
      ! 0 0 0 0                0 1 2 3
      ! 1 1 1 1          =>    0 1 2 3
      ! 2 2 2 2                0 1 2 3
      ! nproj * nbasis         nproj_global * nbasis_stripe

      ! define sending buffer
      allocate( betaq_send( nproj, nbasis ) )

      ! transpose each sending block
      do i = 1, nproc_per_pool

        stripe_lo = stripe_displs(i) + 1
        stripe_hi = stripe_displs(i) + stripe_sizes(i)

        betaq_send( :, stripe_lo : stripe_hi ) = &
        reshape(transpose(betaq( :, stripe_lo : stripe_hi )), &
                (/ nproj, stripe_sizes(i) /))
      enddo

      allocate( betaq_ss( nbasis_stripe, nproj_global ) )

      call MPI_AlltoAllv( betaq_send, &
                          nproj * stripe_sizes * 2, nproj * stripe_displs * 2, MPI_DOUBLE_PRECISION, &
                          betaq_ss, &
                          nbasis_stripe * proj_sizes * 2, nbasis_stripe * proj_displs * 2, MPI_DOUBLE_PRECISION, &
                          intra_pool_comm, ierr )

      ! Now we transpose it back
      allocate( betaq_s( nproj_global, nbasis_stripe ) )
      betaq_s = transpose( betaq_ss )

      deallocate( betaq_send )

    endif

    ! ======================================================================

    ! construct beta_nk, or < beta | nk > for * ALL * projectors

    ! < beta | nk > (nproj * neig_found) =
    ! < beta_1 | 1 k >    < beta_1 | 2 k >   ...
    ! < beta_2 | 1 k >    < beta_2 | 2 k >   ...
    ! ...

    ! Calculate \sum_i < beta_l | B_i > < B_i | nk >
    if ( nproc_per_pool > 1 ) then

      if ( nbasis_stripe > 0 ) then
        ! i iterates over the local stripe
        CALL ZGEMM( 'N', 'N', nproj_global, neig_found, nbasis_stripe, one, &
                    betaq_s, nproj_global, &
                    eigvec_s, nbasis_stripe, &
                    zero, beta_nk, nproj_global )
      else

        beta_nk = zero

      endif

      CALL mp_sum( beta_nk, intra_pool_comm )

      deallocate( betaq_s, betaq_ss )

    else if ( nproc_per_pool == 1 ) then

      ! transform all projectors to eigenbasis for this k
      CALL ZGEMM( 'N', 'N', nproj, neig_found, nbasis, one, &
                  betaq, nproj, &
                  eigvec, nbasis, &
                  zero, beta_nk, nproj )

    endif

    ! ======================================================================

    ! write projectors
    if( proj_output .and. mypoolid == mypoolroot ) then

      offset = ( offset_sk * nproj_global ) * neig_found * 2

      ! write beta_nk
      call mpi_file_write_at( fhproj, offset, &
                              beta_nk(index_global2local, :), nproj_global * neig_found * 2, &
                              MPI_DOUBLE_PRECISION, status, ierr )

#ifdef __SHIRLEY_DEBUG
      if ( ik == 1 ) then          
        write(stdout, *) 'Output projectors to proj.dat for debug. '
        open(unit = fh, file = 'proj.dat', form = 'formatted')
        write(fh, *) beta_nk
        close(fh)
      endif
      write(stdout, *) ' index_global2local = ', index_global2local
#endif

    endif ! write projectors < beta | nk >

    ! ======================================================================
    do icore = 1, corerep % ncore
    ! ======================================================================

      icore_atom = corerep % core(icore) % atom
      core_root = mod(icore_atom - 1, nproc_per_pool)

      ! check if < beta | nk > for icore is processed on the current proc
      if ( core_root /= mypoolid ) cycle 

      ncp = corerep % core(icore) % nproj2      ! # of core levels
      allocate( posn(nbnd, ncp, 9) )            ! < nk | r | phi_c > for all bands; 3 --> 9 Yufeng

      ! This process will have access to the corresponding beta
      ! convert icore_atom to local index
      icore_atom = (icore_atom - 1) / nproc_per_pool + 1
      it = type_atom(icore_atom)
      np = nproj_type(it)                      ! # of valence projectors (beta_l)

      !  beta_nk_c < beta | nk > for icore, select the relevant projectors
      allocate( beta_nk_c( np, nbnd ) )

      ! Note that index_betaq is only used for projectors on the same proc but not global
      forall (i = 1 : np) &
         beta_nk_c(i, :) = beta_nk(index_betaq(i, icore_atom) + proj_displ, :)

#ifdef __SHIRLEY_DEBUG
      if ( ik == 1 ) then
        write(stdout, *) 'Output beta_nk_c to beta_nk_c.dat for debug. '
        open(unit = fh, file = 'beta_nk_c.dat', form = 'formatted')
        write(fh, '(1F20.14)') abs(beta_nk_c(:, 7))
        close(fh)
      endif  
#endif

      allocate(posn_l(nbnd, ncp))

      ! Maybe we can parallelize this over nbnd in future
      ! < nk | r | phi_c > = 
      ! \sum_l  < nk | beta_l > < psi_l | r | phi_c >
      do ixyz = 1, 9 ! Yufeng: ixyz 1 to 9

        CALL ZGEMM( 'C', 'N', nbnd, ncp, np, one, &
                  beta_nk_c, np, &
                  corerep % core(icore) % matrix(1, 1, ixyz), &
                  corerep % core(icore) % nproj1, &
                  zero, posn_l, nbnd)

        posn(:, :, ixyz) = posn_l

      enddo ! ixyz

      ! write xas
      ! 3 --> 9 Yufeng
      offset = (offset_sk * corerep % ncore + icore - 1) * (neig_found * ncp * 9) * 2 
      call mpi_file_write_at( fhxmat, offset, &
                              posn, neig_found * ncp * 9 * 2, &
                              MPI_DOUBLE_PRECISION, status, ierr )

      deallocate(beta_nk_c, posn_l, posn)

    enddo ! loop over icore

! ======================================================================
  enddo ! loop over k-points ik
  enddo ! loop over spin ispin
! ======================================================================

  ! deallocate space
  if( nproc_per_pool > 1 ) then 
    deallocate( eigvec_s )
    deallocate( stripe_sizes, stripe_displs )
    deallocate( proj_sizes, proj_displs )
  end if

  deallocate( beta_nk ) ! You can deallocate an array even if its length is zero.

  if( mypoolid == mypoolroot ) then

    ! close binary dump
    call mpi_file_close( fheigval, ierr )
    if( ierr /= 0 ) &
      call errore('shirley_xas', 'problem closing eigval file', abs(ierr))

  endif

  ! close xmat file
  if( corerep % ncore > 0 ) call mpi_file_close( fhxmat, ierr )
  if( ierr /= 0 ) &
    call errore('shirley_xas', 'problem closing xmat file', abs(ierr))

  ! close eigvec output file
  if( eigvec_output ) then
    call mpi_file_close( fheigvec, ierr )
    if( ierr /= 0 ) &
      call errore('shirley_xas', 'problem closing eigvec file', abs(ierr))
  endif

  if( proj_output ) then
    call mpi_file_close( fhproj, ierr )
    if( ierr /= 0 ) &
      call errore('shirley_xas', 'problem closing proj file', abs(ierr))
  endif

  999 continue
  write(stdout,*) ' waiting for other nodes'
  call mp_barrier( world_comm )
  call mp_world_end
  write(stdout,*) ' end shirley_xas'
  stop
  
  end program shirley_xas
