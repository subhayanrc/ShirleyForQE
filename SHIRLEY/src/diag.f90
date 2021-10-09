! ----------------------------------------------------------------------
  module diag_module
! ----------------------------------------------------------------------

  use kinds, only : dp

  implicit none
  private
  public :: neig, eigval, eigvec, betaq, usq, smatrix
  public :: diag_init, diag_build_hamk, diag_ham, diag_build_dhamk
  public :: diag_build_dham
  public :: diagx_init, diagx_ham
  public :: diag_build_kink, diag_build_vlock, diag_build_vnlk, diag_build_sk

  public :: neig_found
  
  integer :: neig, neig_found
  integer :: il, iu
  real(dp),allocatable :: eigval(:)
  complex(dp),allocatable :: eigvec(:,:)
  complex(dp),allocatable :: betaq(:,:)
  complex(dp),allocatable :: vbetaq(:,:)

  complex(dp),allocatable :: beta_block(:,:), vnl(:,:), ztmp(:,:)
  complex(dp),allocatable :: usq(:,:), sbetaq(:,:), smatrix(:,:)
  integer :: npm

  complex(dp),allocatable :: atomq(:,:)
  complex(dp),allocatable :: vatomq(:,:)
  complex(dp),allocatable :: atom_block(:,:), vhU(:,:)

  integer :: lwork
  complex(dp),allocatable :: work(:)
  real(dp),allocatable :: rwork(:)
  integer,allocatable :: iwork(:), ifail(:)
  complex(dp),allocatable :: atmp(:,:)

  contains

! ----------------------------------------------------------------------
  subroutine diag_init
! ----------------------------------------------------------------------

  use hamq_shirley, only : nbasis, nproj, nproj_nl, nproj_type, ntype, ncpp, natomproj, Hubbard_lmax
  use hamq_pool, only : nproc_per_pool, local_cyclic_dims

  integer :: nblock, nr, nc, ierr
  integer,external :: ILAENV

  call local_cyclic_dims( nr, nc )
  allocate( eigval(nbasis), eigvec(nr,nc), smatrix(nr,nc), &
            betaq(nproj,nbasis), vbetaq(nproj,nbasis), &
            sbetaq(nproj,nbasis), &
            atomq(natomproj,nbasis), vatomq(natomproj,nbasis), &
            stat=ierr )
  if( ierr/=0 ) call errore('diag_init','problem allocating space',1)

  if( nproc_per_pool == 1 ) then
  ! is this even necessary?
  ! tmp space for diagonalization - with what?
  nblock = ILAENV( 1, 'ZHETRD', 'U', nbasis, - 1, - 1, - 1 )
  IF ( nblock < 1 ) nblock = MAX( 1, nbasis )
  IF ( nblock == 1 .OR. nblock >= nbasis ) THEN
     lwork = 2 * nbasis - 1
  ELSE
     lwork = ( nblock + 1 ) * nbasis
  END IF
  allocate( work(lwork), rwork(3*nbasis-2) )
  endif

  if( nproc_per_pool > 1 ) allocate( atmp(nr,nc) )

  npm = maxval( nproj_type(1:ntype) )  
  allocate( beta_block(npm,nbasis), &
            atom_block(2*Hubbard_lmax+1,nbasis), &
            vhU(2*Hubbard_lmax+1,2*Hubbard_lmax+1), &
            vnl(npm,npm), usq(npm,npm), &
            ztmp(npm,nbasis) )

  neig = nbasis

  end subroutine diag_init

! ----------------------------------------------------------------------
  subroutine diagx_init( il_, iu_ )
! ----------------------------------------------------------------------

  use hamq_shirley, only : nbasis, nproj, nproj_nl, nproj_type, ntype, natomproj, Hubbard_lmax
  use hamq_pool, only : nproc_per_pool, local_cyclic_dims

  integer,intent(in) :: il_, iu_

  integer :: nblock, nr, nc, ierr
  integer,external :: ILAENV

  call local_cyclic_dims( nr, nc )

  ! indices of eigenvectors to find
  il=il_
  iu=iu_

  allocate( eigval(nbasis), eigvec(nr,nc) )
  allocate( smatrix(nr,nc) )
  allocate( betaq(nproj,nbasis) )
  allocate( vbetaq(nproj,nbasis) )
  allocate( sbetaq(nproj,nbasis) )
  allocate( atomq(natomproj,nbasis) )
  allocate( vatomq(natomproj,nbasis) )

  if( nproc_per_pool == 1 ) then
  ! tmp space for diagonalization
  nblock = ILAENV( 1, 'ZHETRD', 'U', nbasis, - 1, - 1, - 1 )
  nblock = max( nblock, ILAENV( 1, 'ZUNMTR', 'U', nbasis, - 1, - 1, - 1 ) )
  IF ( nblock < 1 ) nblock = MAX( 1, nbasis )
  IF ( nblock == 1 .OR. nblock >= nbasis ) THEN
     lwork = 2 * nbasis
  ELSE
     lwork = ( nblock + 1 ) * nbasis
  END IF
  allocate( work(lwork) )
  allocate( rwork(7*nbasis), iwork(5*nbasis), ifail(nbasis) )
  endif

  allocate( atmp(nr,nc) )

  npm = maxval( nproj_type(1:ntype) )  
  allocate( beta_block(npm,nbasis), &
            atom_block(2*Hubbard_lmax+1,nbasis), &
            vhU(2*Hubbard_lmax+1,2*Hubbard_lmax+1), &
            vnl(npm,npm), usq(npm,npm), &
            ztmp(npm,nbasis) )

  neig = nbasis

  end subroutine diagx_init


! ----------------------------------------------------------------------
  subroutine diag_build_hamk( kvec, cartesian, ispin )
! ----------------------------------------------------------------------

  use hamq_shirley, only : natom, type_atom, nproj_type_nl, nbasis, &
                           nproj, nproj_nl, &
                           index_betaq, &
                           index_nlproj_type, ncpp, &
                           build_hamq_local, build_hamq_projs, nproj_type, ntype, &
                           lda_plus_u, build_hamq_atomprojs, &
                           natomproj, Hubbard_lmax, &
                           Hubbard_l, Hubbard_U, Hubbard_alpha, &
                           index_ldaUq
  use shirley_input_module, only : vnl_atom, usq_atom, vhU_atom
  use scalapack_module, only : CTXT_, MB_, DLEN_
  use hamq_pool, only : intra_pool_comm, natom_global, &
                        nproj_global, nproc_per_pool, &
                        context_cyclic, mypoolid, context_striped, &
                        desc_cyclic, desc_striped, &
                        ndiag_cyclic, diag_cyclic, &
                        local_cyclic_dims, local_striped_dim, &
                        cyclic_localindex, &
                        npool, mypool, context_global
  use mp, only : mp_sum, mp_max, mp_barrier, mp_bcast
  use mp_world, only : mpime, root
  use io_global, only : stdout

  complex(dp),parameter :: ONE=(1.d0,0.d0)
  complex(dp),parameter :: ZERO=(0.d0,0.d0)

  real(dp),intent(in) :: kvec(3)
  logical,intent(in) :: cartesian
  integer,intent(in) :: ispin

  integer :: iatom, it, np, i, j
  integer :: il, jl
  logical :: islocal
  complex(dp),external :: ZDOTC
  complex(dp) :: hamij, sij
  complex(dp),allocatable :: nlham(:,:), tmp_stripe(:,:), nlham_s(:,:), smatrix_s(:,:)

  integer :: iproc, rdest, cdest
  integer :: nbasis_stripe, nbasis_stripe_max
  integer :: nr, nc
  integer :: ipool, ierr
  integer :: nprow, npcol, myrow, mycol
  integer :: desc1(DLEN_), desc2(DLEN_)

  integer :: nstripe, istripe
  integer :: ioffset

  !if( mpime==root ) write(stdout,*) ' building Hamiltonian'

  ! build the Hamiltonian for this k-point
  ! local contribution
  call build_hamq_local( kvec, cartesian, ispin, eigvec )

  ! non-local contribution

  ! projectors
  call build_hamq_projs( kvec, cartesian, betaq )

  ! atomic projectors for LDA+U
  if( lda_plus_u ) call build_hamq_atomprojs( kvec, cartesian, atomq )

  ! potential summed atom by atom
  if( nproj_global > 0 ) then

    call local_striped_dim( nbasis_stripe )
    nbasis_stripe_max=nbasis_stripe
    call mp_max( nbasis_stripe_max, intra_pool_comm )

    allocate( nlham(nbasis_stripe_max,nbasis), &
              nlham_s(nbasis_stripe,nbasis), &
              smatrix_s(nbasis_stripe,nbasis), &
              tmp_stripe(nproj,nbasis_stripe_max) )

    ! important
    vbetaq=ZERO
    if( .not. ncpp ) sbetaq=ZERO

    do iatom=1,natom
      it = type_atom(iatom)
      np = nproj_type_nl(it)

      if( np > 0 ) then

      forall( i=1:np, j=1:nbasis ) &
        beta_block(i,j) = &
          betaq(index_betaq(index_nlproj_type(i,it),iatom),j)

      forall( i=1:np, j=1:np ) &
        vnl(i,j) = cmplx( vnl_atom(iatom,ispin)%matrix( index_nlproj_type(i,it), &
                                                  index_nlproj_type(j,it) ) )
      CALL ZGEMM( 'N', 'N', np, nbasis, np, ONE, &
                  vnl, npm, &
                  beta_block, npm, &
        ZERO, vbetaq(index_betaq(index_nlproj_type(1,it),iatom),1), nproj )


      if( .not. ncpp ) then
        ! build the S-matrix
        forall( i=1:np, j=1:np ) &
          usq(i,j) = cmplx( usq_atom(iatom)%matrix( index_nlproj_type(i,it), &
                                                    index_nlproj_type(j,it) ) )

          CALL ZGEMM( 'N', 'N', np, nbasis, np, ONE, &
                      usq, npm, &
                      beta_block, npm, &
            ZERO, sbetaq(index_betaq(index_nlproj_type(1,it),iatom),1), nproj )
      endif ! .not. ncpp
      endif ! np > 0
    enddo ! loop over atoms

    ! offset for striping
    istripe=0

    ! loop over processors in the pool
    do iproc=0,nproc_per_pool-1

      ! all processors need to know the size of the stripe for process iproc
      if( mypoolid==iproc ) nstripe=nbasis_stripe
      call mp_bcast( nstripe, iproc, intra_pool_comm )

      ! nothing to do if no atoms on this proc (nproj should be zero too)
      if( nproj==0 ) then
        nlham=ZERO
      ! otherwise
      else
        ! load tmp_stripe for the current proc
        forall( i=1:nproj, j=1:nstripe ) &
          tmp_stripe(i,j) = betaq(i,istripe+j)

        ! make the product beta . (v . beta)
        CALL ZGEMM( 'C', 'N', nstripe, nbasis, nproj, ONE, &
                    tmp_stripe, nproj, &
                    vbetaq, nproj, &
                    ZERO, nlham, nbasis_stripe_max )
        ! pad with zeros
        if( nbasis_stripe_max > nstripe ) &
          forall( i=nstripe+1:nbasis_stripe_max, j=1:nbasis ) &
            nlham(i,j)=ZERO
      endif

      ! sum over procs to accumulate over atoms and save to correct process
      ! coordinate in the grid
      rdest=iproc
      cdest=0

      call ZGSUM2D(context_striped,'A',' ',nbasis_stripe_max,nbasis, &
                   nlham,nbasis_stripe_max,rdest,cdest)
!      call mp_sum( nlham, intra_pool_comm )

      ! copy to correct sized array
      if( mypoolid==iproc) then
        forall(i=1:nbasis_stripe, j=1:nbasis ) nlham_s(i,j) = nlham(i,j)
      endif

      if( .not. ncpp ) then
        ! nothing to do if no atoms on this proc (nproj should be zero too)
        if( nproj==0 ) then
          nlham=ZERO
        ! otherwise
        else
          ! make the product beta . (s . beta)
          CALL ZGEMM( 'C', 'N', nstripe, nbasis, nproj, ONE, &
                      tmp_stripe, nproj, &
                      sbetaq, nproj, &
                      ZERO, nlham, nbasis_stripe_max )
          ! pad with zeros
          if( nbasis_stripe_max > nstripe ) &
            forall( i=nstripe+1:nbasis_stripe_max, j=1:nbasis ) &
              nlham(i,j)=ZERO
        endif

        ! sum over procs to accumulate over atoms and save to correct process
        ! coordinate in the grid
        rdest=iproc
        cdest=0
        call ZGSUM2D(context_striped,'A',' ',nbasis_stripe_max,nbasis, &
                     nlham,nbasis_stripe_max,rdest,cdest)
!        call mp_sum( nlham, intra_pool_comm )

        ! copy to correct sized array
        if( mypoolid==iproc) then
          forall(i=1:nbasis_stripe, j=1:nbasis ) smatrix_s(i,j) = nlham(i,j)
        endif
      endif  ! .not. ncpp

      ! update offset for striping
      istripe=istripe+nstripe

    enddo ! loop over procs

    deallocate( tmp_stripe ) 

    ! LDA+U
    if( lda_plus_u ) then

      allocate( tmp_stripe(natomproj,nbasis_stripe_max) )

      vatomq=ZERO

      ioffset=0

      do iatom=1,natom
        it = type_atom(iatom)

        if( Hubbard_U(it) /= 0.d0 .or. Hubbard_alpha(it) /= 0.d0 ) then
          np = 2*Hubbard_l(it)+1

          if( np > 0 ) then

            forall( i=1:np, j=1:nbasis ) &
              atom_block(i,j) = &
                atomq(ioffset+i,j)

            forall( i=1:np, j=1:np ) &
              vhU(i,j) = cmplx( vhU_atom(iatom,ispin)%matrix(i,j) )

            CALL ZGEMM( 'N', 'N', np, nbasis, np, ONE, &
                        vhU, (2*Hubbard_lmax+1), &
                        atom_block, (2*Hubbard_lmax+1), &
                        ZERO, vatomq(ioffset+1,1), natomproj )
            ioffset=ioffset+np
          endif ! np > 0
        endif ! U,alpha/=0
      enddo ! loop over atoms

      ! offset for striping
      istripe=0

      ! loop over processors in the pool
      do iproc=0,nproc_per_pool-1

        ! all processors need to know the size of the stripe for process iproc
        if( mypoolid==iproc ) nstripe=nbasis_stripe
        call mp_bcast( nstripe, iproc, intra_pool_comm )

        ! nothing to do if no atom projections on this proc (nproj should be zero too)
        if( natomproj==0 ) then
          nlham=ZERO
        ! otherwise
        else
          ! load tmp_stripe for the current proc
          forall( i=1:natomproj, j=1:nstripe ) &
            tmp_stripe(i,j) = atomq(i,istripe+j)

          ! make the product beta . (v . beta)
          CALL ZGEMM( 'C', 'N', nstripe, nbasis, natomproj, ONE, &
                      tmp_stripe, natomproj, &
                      vatomq, natomproj, &
                      ZERO, nlham, nbasis_stripe_max )
          ! pad with zeros
          if( nbasis_stripe_max > nstripe ) &
            forall( i=nstripe+1:nbasis_stripe_max, j=1:nbasis ) &
              nlham(i,j)=ZERO
        endif

        ! sum over procs to accumulate over atoms and save to correct process
        ! coordinate in the grid
        rdest=iproc
        cdest=0
        call ZGSUM2D(context_striped,'A',' ',nbasis_stripe_max,nbasis, &
                     nlham,nbasis_stripe_max,rdest,cdest)

        ! add to correct sized array
        if( mypoolid==iproc) then
          forall(i=1:nbasis_stripe, j=1:nbasis ) &
            nlham_s(i,j) = nlham_s(i,j) + nlham(i,j)
        endif

        ! update offset for striping
        istripe=istripe+nstripe

      enddo ! loop over procs

      deallocate( tmp_stripe )

    endif ! lda_plus_u

    ! non-local Hamiltonian computed

    ! re-use the array nlham to redistribute the non-local hamiltonian
    deallocate( nlham )
    call local_cyclic_dims( nr, nc )
    allocate( nlham(nr,nc) )

    if( nproc_per_pool > 1 ) then

      if( mpime==root ) write(stdout,*) ' redistribute non-local Hamiltonian'
      ! by now every process has its own stripe of the non-local hamiltonian
      ! we need to redistribute this data in a block-cyclic fashion
      call PZGEMR2D(nbasis,nbasis, &
                    nlham_s,1,1,desc_striped, &
                    nlham,1,1,desc_cyclic,desc_cyclic(CTXT_))

    else

      nlham = nlham_s

    endif

    eigvec = eigvec + nlham
    deallocate( nlham )

    if( .not. ncpp ) then
      if( nproc_per_pool > 1 ) then

      smatrix=ZERO
      call PZGEMR2D(nbasis,nbasis, &
                    smatrix_s,1,1,desc_striped, &
                    smatrix,1,1,desc_cyclic,desc_cyclic(CTXT_))

      else

      smatrix = smatrix_s

      endif

      ! add identity 
      if( ndiag_cyclic > 0 ) &
        forall(i=1:ndiag_cyclic) &
          smatrix(diag_cyclic(1,i),diag_cyclic(2,i)) &
            = smatrix(diag_cyclic(1,i),diag_cyclic(2,i)) + ONE

    endif

    deallocate( nlham_s, smatrix_s )

    !if( mpime==root ) write(stdout,*) ' done with nonlocal '

  endif ! if nproj > 0

  !if( mpime==root ) write(stdout,*) ' done building Hamiltonian'
  return
  end subroutine diag_build_hamk


! ----------------------------------------------------------------------
  subroutine diag_build_dham( ideriv, kvec, cartesian, ispin, dham )
! ----------------------------------------------------------------------

  ! this computes the derivative in the crystal coord (dimensionless) dirn

  use hamq_shirley, only : natom, type_atom, nproj, nproj_type_nl, nbasis, &
                           index_betaq, index_nlproj_type, &
                           build_dhamq_local, build_hamq_projs, build_dhamq_projs
  use shirley_input_module, only : vnl_atom
  use mp_world, only : mpime

  complex(dp),parameter :: ONE=(1.d0,0.d0)
  complex(dp),parameter :: ZERO=(0.d0,0.d0)

  integer,intent(in) :: ideriv
  real(dp),intent(in) :: kvec(3)
  logical,intent(in) :: cartesian
  integer,intent(in) :: ispin
  complex(dp),intent(out) :: dham(nbasis,nbasis)

  complex(dp) :: dbetaq(nproj,nbasis), dbeta_block(npm,nbasis)
  integer :: iatom, it, np, i, j

  ! assumes that eigvec has already been determined by diagonalization of Ham
  ! with call to diag_hamk

  ! build the Hamiltonian for this k-point
  ! local contribution
  call build_dhamq_local( ideriv, kvec, cartesian, dham )

  ! non-local contribution
  ! projectors
  ! I assume that diag_build_hamk was just called, so betaq is already loaded 
  ! call build_hamq_projs( kvec, cartesian, betaq )

  ! the derivative remain to be computed
  call build_dhamq_projs( ideriv, kvec, cartesian, dbetaq )

  ! potential summed atom by atom
  do iatom=1,natom
    it = type_atom(iatom)
    np = nproj_type_nl(it)

    forall( i=1:np, j=1:nbasis ) &
      beta_block(i,j) = &
        betaq(index_betaq(index_nlproj_type(i,it),iatom),j)

    forall( i=1:np, j=1:nbasis ) &
      dbeta_block(i,j) = &
        dbetaq(index_betaq(index_nlproj_type(i,it),iatom),j)

    forall( i=1:np, j=1:np ) &
      vnl(i,j) = cmplx( vnl_atom(iatom,ispin)%matrix( index_nlproj_type(i,it), &
                                                index_nlproj_type(j,it) ) )

    ! dbeta V beta
    CALL ZGEMM( 'N', 'N', np, nbasis, np, ONE, &
                vnl, npm, &
                beta_block, npm, &
                ZERO, ztmp, npm )
    CALL ZGEMM( 'C', 'N', nbasis, nbasis, np, ONE, &
                dbeta_block, npm, &
                ztmp, npm, &
                ONE, dham, nbasis )

    ! beta V dbeta
    CALL ZGEMM( 'N', 'N', np, nbasis, np, ONE, &
                vnl, npm, &
                dbeta_block, npm, &
                ZERO, ztmp, npm )
    CALL ZGEMM( 'C', 'N', nbasis, nbasis, np, ONE, &
                beta_block, npm, &
                ztmp, npm, &
                ONE, dham, nbasis )
  enddo

  end subroutine diag_build_dham


! ----------------------------------------------------------------------
  subroutine diag_build_dhamk( ideriv, kvec, cartesian, ispin, dhamk )
! ----------------------------------------------------------------------

  use hamq_shirley, only : natom, type_atom, nproj, nproj_type_nl, nbasis, &
                           index_betaq, index_nlproj_type, &
                           build_dhamq_local, build_hamq_projs, build_dhamq_projs

  complex(dp),parameter :: ONE=(1.d0,0.d0)
  complex(dp),parameter :: ZERO=(0.d0,0.d0)

  integer,intent(in) :: ideriv
  real(dp),intent(in) :: kvec(3)
  logical,intent(in) :: cartesian
  integer,intent(in) :: ispin
  complex(dp),intent(out) :: dhamk(nbasis,nbasis)

  !complex(dp) :: etmp(nbasis,nbasis)
  complex(dp),allocatable :: etmp(:,:)

  allocate( etmp(nbasis,nbasis) )

  ! assumes that eigvec has already been determined by diagonalization of Ham
  ! with call to diag_hamk

  ! build the Hamiltonian derivative
  call diag_build_dham( ideriv, kvec, cartesian, ispin, dhamk )

  ! now transform to eigenbasis
  ! Z dH Z
  CALL ZGEMM( 'N', 'N', nbasis, neig_found, nbasis, ONE, &
              dhamk, nbasis, &
              eigvec, nbasis, &
              ZERO, etmp, nbasis )

  CALL ZGEMM( 'C', 'N', neig_found, nbasis, nbasis, ONE, &
              eigvec, nbasis, &
              etmp, nbasis, &
              ZERO, dhamk, nbasis )

  deallocate( etmp )

  end subroutine diag_build_dhamk

! ----------------------------------------------------------------------
  subroutine diag_build_kink( kvec, cartesian, kink )
! ----------------------------------------------------------------------

  use hamq_shirley, only : nbasis, &
                           build_hamq_kin

  complex(dp),parameter :: ONE=(1.d0,0.d0)
  complex(dp),parameter :: ZERO=(0.d0,0.d0)

  real(dp),intent(in) :: kvec(3)
  logical,intent(in) :: cartesian
  complex(dp),intent(out) :: kink(nbasis,nbasis)

  complex(dp),allocatable :: etmp(:,:)

  allocate( etmp(nbasis,nbasis) )

  ! assumes that eigvec has already been determined by diagonalization of Ham
  ! with call to diag_hamk

  ! build the Kinetic contribution to the Hamiltonian
  call build_hamq_kin( kvec, cartesian, kink )

  ! now transform to eigenbasis
  ! Z dH Z
  CALL ZGEMM( 'N', 'N', nbasis, neig_found, nbasis, ONE, &
              kink, nbasis, &
              eigvec, nbasis, &
              ZERO, etmp, nbasis )

  CALL ZGEMM( 'C', 'N', neig_found, nbasis, nbasis, ONE, &
              eigvec, nbasis, &
              etmp, nbasis, &
              ZERO, kink, nbasis )

  deallocate( etmp )

  end subroutine diag_build_kink

! ----------------------------------------------------------------------
  subroutine diag_build_vlock( kvec, cartesian, ispin, vlock )
! ----------------------------------------------------------------------

  use hamq_shirley, only : nbasis, &
                           build_hamq_vloc

  complex(dp),parameter :: ONE=(1.d0,0.d0)
  complex(dp),parameter :: ZERO=(0.d0,0.d0)

  real(dp),intent(in) :: kvec(3)
  logical,intent(in) :: cartesian
  integer,intent(in) :: ispin
  complex(dp),intent(out) :: vlock(nbasis,nbasis)

  !complex(dp) :: etmp(nbasis,nbasis)
  complex(dp),allocatable :: etmp(:,:)

  allocate( etmp(nbasis,nbasis) )

  ! assumes that eigvec has already been determined by diagonalization of Ham
  ! with call to diag_hamk

  ! build the Local Potential contribution to the Hamiltonian
  call build_hamq_vloc( kvec, cartesian, ispin, vlock )

  ! now transform to eigenbasis
  ! Z dH Z
  CALL ZGEMM( 'N', 'N', nbasis, neig_found, nbasis, ONE, &
              vlock, nbasis, &
              eigvec, nbasis, &
              ZERO, etmp, nbasis )

  CALL ZGEMM( 'C', 'N', neig_found, nbasis, nbasis, ONE, &
              eigvec, nbasis, &
              etmp, nbasis, &
              ZERO, vlock, nbasis )

  deallocate( etmp )

  end subroutine diag_build_vlock

! ----------------------------------------------------------------------
  subroutine diag_build_vnlk( kvec, cartesian, ispin, vnlk )
! ----------------------------------------------------------------------

  use hamq_shirley, only : natom, type_atom, nproj, nproj_type_nl, nbasis, &
                           index_betaq, index_nlproj_type, &
                           build_hamq_projs
  use shirley_input_module, only : vnl_atom

  complex(dp),parameter :: ONE=(1.d0,0.d0)
  complex(dp),parameter :: ZERO=(0.d0,0.d0)

  real(dp),intent(in) :: kvec(3)
  logical,intent(in) :: cartesian
  integer,intent(in) :: ispin
  complex(dp),intent(out) :: vnlk(nbasis,nbasis)

  complex(dp),allocatable :: etmp(:,:)
  integer :: iatom, it, np, i, j

  allocate( etmp(nbasis,nbasis) )

  ! assumes that eigvec has already been determined by diagonalization of Ham
  ! with call to diag_hamk

  vnlk=ZERO

  ! projectors
  call build_hamq_projs( kvec, cartesian, betaq )
  ! potential summed atom by atom
  do iatom=1,natom
    it = type_atom(iatom)
    np = nproj_type_nl(it)

    forall( i=1:np, j=1:nbasis ) &
      beta_block(i,j) = &
        betaq(index_betaq(index_nlproj_type(i,it),iatom),j)

    forall( i=1:np, j=1:np ) &
      vnl(i,j) = cmplx( vnl_atom(iatom,ispin)%matrix( index_nlproj_type(i,it), &
                                                index_nlproj_type(j,it) ) )

    CALL ZGEMM( 'N', 'N', np, nbasis, np, ONE, &
                vnl, npm, &
                beta_block, npm, &
                ZERO, ztmp, npm )
    CALL ZGEMM( 'C', 'N', nbasis, nbasis, np, ONE, &
                beta_block, npm, &
                ztmp, npm, &
                ONE, vnlk, nbasis )
  enddo

  ! now transform to eigenbasis
  ! Z dH Z
  CALL ZGEMM( 'N', 'N', nbasis, neig_found, nbasis, ONE, &
              vnlk, nbasis, &
              eigvec, nbasis, &
              ZERO, etmp, nbasis )

  CALL ZGEMM( 'C', 'N', neig_found, nbasis, nbasis, ONE, &
              eigvec, nbasis, &
              etmp, nbasis, &
              ZERO, vnlk, nbasis )

  deallocate( etmp )

  end subroutine diag_build_vnlk

! ----------------------------------------------------------------------
  subroutine diag_build_sk( kvec, cartesian, sk )
! ----------------------------------------------------------------------

  use hamq_shirley, only : nbasis

  complex(dp),parameter :: ONE=(1.d0,0.d0)
  complex(dp),parameter :: ZERO=(0.d0,0.d0)

  real(dp),intent(in) :: kvec(3)
  logical,intent(in) :: cartesian
  complex(dp),intent(out) :: sk(nbasis,nbasis)

  !complex(dp) :: etmp(nbasis,nbasis)
  complex(dp),allocatable :: etmp(:,:)

  allocate( etmp(nbasis,nbasis) )

  ! assumes that eigvec has already been determined by diagonalization of Ham
  ! with call to diag_hamk

  sk = smatrix

  ! now transform to eigenbasis
  ! Z dH Z
  CALL ZGEMM( 'N', 'N', nbasis, neig_found, nbasis, ONE, &
              sk, nbasis, &
              eigvec, nbasis, &
              ZERO, etmp, nbasis )

  CALL ZGEMM( 'C', 'N', neig_found, nbasis, nbasis, ONE, &
              eigvec, nbasis, &
              etmp, nbasis, &
              ZERO, sk, nbasis )

  deallocate( etmp )

  end subroutine diag_build_sk

! ----------------------------------------------------------------------
  subroutine diag_ham
! ----------------------------------------------------------------------

  use mp_world, only : mpime
  use hamq_shirley, only : ncpp, nbasis
  use hamq_pool, only : nproc_per_pool, diag_pool_matrix
  use io_global, only : stdout

  integer :: ierr
  integer :: i,j
  logical :: islocal

  ierr=0
  ! diagonalize matrix to find eigenvalues and vectors
  if( ncpp ) then

    if( nproc_per_pool>1 ) then
      atmp=eigvec
      call diag_pool_matrix( atmp, eigval, eigvec, ierr )
    else
      CALL ZHEEV( 'V', 'U', neig, eigvec, neig, eigval, &
                  work, lwork, rwork, ierr )
    endif

  else
  ! ultrasoft case
  ! Hx=eSx
    if( nproc_per_pool>1 ) then
      atmp=eigvec
      call diag_pool_matrix( atmp, smatrix, eigval, eigvec, ierr )
    else
      CALL ZHEGV( 1, 'V', 'U', neig, eigvec, neig, smatrix, neig, &
                  eigval, work, lwork, rwork, ierr )
    endif

  endif

  call errore('diag_ham','error in diagonalization of hamiltonian',abs(ierr))
  neig_found = neig

  end subroutine diag_ham


! ----------------------------------------------------------------------
  subroutine diagx_ham
! ----------------------------------------------------------------------

  use mp_world, only : mpime
  use hamq_shirley, only : ncpp, nbasis
  use hamq_pool, only : nproc_per_pool, diag_pool_matrix
  use io_global, only : stdout

  integer :: ierr
  real(dp) :: vl, vu, abstol
  integer :: i, j
  logical :: islocal


  abstol=-1.d0
  atmp = eigvec

  if( ncpp ) then
    if( nproc_per_pool>1 ) then
      call diag_pool_matrix( atmp, eigval, eigvec, ierr )
      neig_found=neig
    else
      ! diagonalize matrix to find eigenvalues and vectors
      CALL ZHEEVX( 'V', 'I', 'U', neig, atmp, neig, &
                   vl, vu, il, iu, abstol, neig_found, eigval, eigvec, neig, &
                   work, lwork, rwork, iwork, ifail, ierr )
    endif
  else
  ! ultrasoft case
  ! Hx=eSx
    if( nproc_per_pool>1 ) then
      call diag_pool_matrix( atmp, smatrix, eigval, eigvec, ierr, il, iu, neig_found )
    else
      CALL ZHEGVX( 1, 'V', 'I', 'U', neig, atmp, neig, smatrix, neig, &
                   vl, vu, il, iu, abstol, neig_found, eigval, eigvec, neig, &
                   work, lwork, rwork, iwork, ifail, ierr )
    endif
  endif
  call errore('diagx_ham','error in diagonalization of hamiltonian',abs(ierr))

  end subroutine diagx_ham


  end module diag_module

