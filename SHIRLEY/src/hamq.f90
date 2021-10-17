! ---------------------------------------------------------------------- 
  subroutine hamq( )
! ---------------------------------------------------------------------- 

  USE parallel_include
  USE io_global,  ONLY : stdout
  USE constants, ONLY : rytoev, fpi, e2, tpi
  USE ions_base, ONLY : ntyp=>nsp, nat, ityp
  USE cell_base
  USE gvect  
  use gvecw,     only : gcutw, ecutwfc
  USE klist, ONLY: xk, nks, nkstot, nelec, ngk, igk_k
  USE wvfct
  use control_flags, only : gamma_only
  USE uspp, ONLY: nkb, vkb, dvan, okvan, deeq, nhtol, qq_at 
! davegp
  USE us, ONLY: dq, tab
  USE cellmd, ONLY : cell_factor
! davegp
  USE scf,      ONLY : vrs, vltot, rho, v
  USE uspp_param, ONLY : upf, nhm, nh, newpseudo
  USE becmod,   ONLY: becp, calbec, &
                      allocate_bec_type, deallocate_bec_type, bec_type
  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc, seqopn
  USE wavefunctions, ONLY: evc, psic
  use buffers, only : open_buffer, close_buffer, get_buffer
  USE lsda_mod, ONLY : nspin, lsda
  use ldaU, only : lda_plus_u_=>lda_plus_u, &
                   Hubbard_l_ => Hubbard_l, &
                   Hubbard_lmax_ => Hubbard_lmax, &
                   Hubbard_U_ => Hubbard_U, &
                   Hubbard_alpha_ => Hubbard_alpha, &
                   oatwfc, nwfcU, offsetU
  USE noncollin_module,     ONLY : npol
  use basis, only : natomwfc, swfcatom
  use ener, only : ehart, etxc, vtxc
  use mp,        only : mp_max, mp_sum, mp_barrier, mp_bcast, mp_root_sum
  use mp_scatt,  only :  mp_scatter_size, mp_scatter
  use mp_world,  ONLY : mpime, root, nproc, world_comm
  use mp_pools,  only : intra_pool_comm
  USE mp_wave,   ONLY : mergewf
  use fft_base,  only : dffts
  use fft_interfaces, only : fwfft, invfft

  use hamq_shirley
  use shirley_ham_input, only : nkgrid, ikgrid, ksplord, &
                                local_channel, debug, &
                                band_subset

  implicit none

  real(dp),parameter :: eps=1.d-10
  complex(dp),parameter :: zero=cmplx(0.d0,0.d0)

  real(dp),allocatable :: xk_x(:), xk_y(:), xk_z(:), xk_cart(:,:)
  integer,allocatable :: kmap(:,:)
  integer :: kxord, kyord, kzord, nk3(3)
  integer :: nkr

  integer :: ik, ibnd, jbnd, ispin
  integer,allocatable :: ibnd_indx(:)
  real(dp),allocatable :: norm(:)
  real(dp) :: v0, charge, fac

  complex(dp),allocatable :: ham(:), hamtmp(:)
  integer(kind=MPI_OFFSET_KIND) :: fposnloc
  integer(kind=MPI_OFFSET_KIND) :: fposnprj

  integer :: nbnd_l, nbnd_lmx, ip, ij, ijl, nbnd_ip
  integer,allocatable :: lengths(:), starts(:)

  type(matrix_list),allocatable :: vnl_atom(:,:)
  type(matrix_list),allocatable :: usq_atom(:)
  type(matrix_list),allocatable :: vhU_atom(:,:)
  complex(dp),allocatable :: becp_tmp(:), becp_l(:,:)
  complex(dp),allocatable :: swfcatom_ldaU(:,:)
  integer :: na, nt, ih, jh

  integer :: nkb_orig
  integer :: m1, m2
  type(bec_type) :: proj

  integer :: ierr, ixyz, ig, i, j, k, n, ikb
  complex(dp),allocatable :: jtmp(:)
  real(dp),allocatable :: gtmp(:)
  real(dp) :: qvec(3), tqvec(3)
  real(dp) :: atp(3,3), bgp(3,3)
  integer :: iunhq, iunvnl, igwx
  logical :: exst

  REAL(DP), EXTERNAL :: erf
  COMPLEX(DP), EXTERNAL :: ZDOTU
  integer,external :: freeunit

  integer,external :: n_plane_waves

  WRITE( stdout, '(/5x,"Calling hamq .... ",/)')

  ! ======================================================================
  ! sort out which band subset we will work with
  ! ======================================================================
  if( band_subset(1) > band_subset(2) ) then
    i=band_subset(2)
    band_subset(2) = band_subset(1)
    band_subset(1) = i
  endif
  if( band_subset(1)>=nbnd .or. band_subset(1)<=0 ) band_subset(1) = 1
  if( band_subset(2)>=nbnd .or. band_subset(2)<=0 ) band_subset(2) = nbnd

  if( band_subset(2)-band_subset(1)+1 < nbnd ) then
    write(stdout,*) ' Requested band subset:', band_subset(1), &
                    ' ... ', band_subset(2)
    write(stdout,*) ' Reducing total number of bands from ', nbnd, &
                    ' to ', band_subset(2)-band_subset(1)+1
  endif

  nbnd = band_subset(2)-band_subset(1)+1
  allocate( ibnd_indx(nbnd) )
  do i=1,nbnd
    ibnd_indx(i) = band_subset(1)+i-1
  enddo

  ! I think this shouldn't be allowed - I'm pretty sure that the optimal basis is not real-valued
  ! I'm not sure if this has been implemented everywhere.
  ! Check this in the future
  IF ( gamma_only ) THEN
     WRITE( stdout, '(5x,"gamma-point specific algorithms are used",/)')
  END IF
  !
  call summary
  !
  ! ======================================================================
  ! be sure that number of planewaves is set up correctly in ngk(:)
  ! try to use ngk(ik) instead of npw from now on
  ! ======================================================================
  npwx = n_plane_waves (gcutw, nkstot, xk, g, ngm )
  write(stdout,*) '     ecutwfc = ', ecutwfc
  write(stdout,*) '      tpiba2 = ', tpiba2
  write(stdout,*) ' nks, nkstot = ', nks, nkstot
  write(stdout,*) '          xk = ', xk(1:3,1:nkstot)
  write(stdout,*) '         npw = ', ngk(1:nks)
  write(stdout,*) '        npwx = ', npwx
  write(stdout,*) '         nkb = ', nkb

  ! ======================================================================
  ! open files for output of Hamiltonian
  ! ======================================================================
  call open_hamq( trim(tmp_dir)//trim(prefix) )


  ! ======================================================================
  ! allocate space for Hamiltonian
  ! ======================================================================

  ! scatter upper triangular
  call mp_scatter_size( (nbnd*(nbnd+1))/2, nbnd_l, root )
  ! collect info on dimensions
  nbnd_lmx=nbnd_l
  call mp_max( nbnd_lmx, world_comm )

  allocate( lengths(nproc), starts(nproc) )
  call mpi_allgather( nbnd_l, 1, MPI_INTEGER, lengths, 1, MPI_INTEGER, &
                      MPI_COMM_WORLD, ierr )
  starts(1)=0
  do ip=2,nproc
    starts(ip)=starts(ip-1)+lengths(ip-1)
  enddo

  allocate( ham(nbnd_l), hamtmp(nbnd_lmx), jtmp(npwx), gtmp(npwx), stat=ierr )
  if( ierr/=0 ) then
    call errore('hamq','problem allocating space for shirley local hamiltonian',1)
  endif


  ! set up output for Hamiltonian module hamq_shirley
  call init_stdout( stdout )

  ! set up system variables
  call init_system( nelec, alat, omega, at, bg, tpiba, nspin, lda_plus_u_ )
 

  ! set up transformation matrices from crystal to Cartesian and back
  atp=transpose(at)/tpiba
  bgp=bg*tpiba
  !call init_qtrans( transpose(at)/tpiba, bg*tpiba )
  call init_qtrans( atp, bgp )

! N.B. I'm not allocating space
! instead initialize nbasis here
  nbasis = nbnd

  jtmp = zero
  gtmp = zero


  CALL gk_sort( xk(1,1), ngm, g, ecutwfc / tpiba2, npw, igk_k(1,1), g2kin )
  g2kin = g2kin * tpiba2

  ! load basis functions
  write(stdout,*)
  write(stdout,*) ' load wave function'
  CALL open_buffer ( iunwfc, 'wfc', nwordwfc, 1, exst )
  CALL get_buffer( evc, nwordwfc, iunwfc, 1 )
  CALL close_buffer ( iunwfc, 'KEEP' )

  ! report norms
  allocate( norm(nbnd) )
  do ibnd=1,nbnd
    norm(ibnd) = dot_product( evc(:,ibnd_indx(ibnd)), evc(:,ibnd_indx(ibnd)) )
  enddo
  call mp_sum( norm, intra_pool_comm )
  do ibnd=1,nbnd
    if( abs(norm(ibnd)-1.d0) > eps ) then
      write(stdout,'(a,i6,a,f14.10)') ' band ', ibnd, ' norm = ', norm(ibnd)
      call errore('hamq','wave function norm is not 1',1)
    endif
  enddo
  deallocate( norm )
 
  write(stdout,*)
  write(stdout,*) ' construct hamiltonian:'
  write(stdout,*)

  ! initialize file position for local write
  fposnloc=0

  ! ======================================================================
  ! kinetic energy - constant term
  ! ======================================================================
  write(stdout,*) ' 1. kinetic energy - constant term'
  do ip=1,nproc
    i=0 ; j=0
    do ij=1,(nbnd*(nbnd+1))/2
      i=i+1
      if( i>j ) then
        j=j+1; i=1
      endif
      if( ij <= starts(ip) .or. ij > starts(ip)+lengths(ip) ) cycle

      ijl=ij-starts(ip)
      if( ijl==1 .or. i==1 ) then
        jtmp(1:npw) = cmplx(g2kin(1:npw),kind=dp) * conjg(evc(1:npw,ibnd_indx(j)))
      endif
      hamtmp(ijl) = ZDOTU( npw, evc(1,ibnd_indx(i)), 1, jtmp, 1 )
    enddo

    call mp_sum( hamtmp, intra_pool_comm )
    if( mpime==ip-1 ) then
      ham = conjg( hamtmp(1:nbnd_l) )
    endif
  enddo

  call write_hamloc( fposnloc, nbnd, ham )


  ! ======================================================================
  ! kinetic energy - linear term
  ! ======================================================================
  write(stdout,*) ' 2. kinetic energy - linear term'
  do ixyz=1,3
    gtmp(1:npw) = g(ixyz,igk_k(1:npw,1))*tpiba

    do ip=1,nproc
      i=0 ; j=0
      do ij=1,(nbnd*(nbnd+1))/2
        i=i+1
        if( i>j ) then
          j=j+1; i=1
        endif
        if( ij <= starts(ip) .or. ij > starts(ip)+lengths(ip) ) cycle

        ijl=ij-starts(ip)
        if( ijl==1 .or. i==1 ) then
          jtmp(1:npw) = cmplx(gtmp(1:npw),kind=dp) * conjg(evc(1:npw,ibnd_indx(j)))
        endif
        hamtmp(ijl) = ZDOTU( npw, evc(1,ibnd_indx(i)), 1, jtmp, 1 )
      enddo

      call mp_sum( hamtmp, intra_pool_comm )
      if( mpime==ip-1 ) then
        ham = conjg( hamtmp(1:nbnd_l) )
      endif
    enddo

    call write_hamloc( fposnloc, nbnd, ham )
  enddo

  ! ======================================================================
  ! local potential
  ! ======================================================================
  do ispin=1,nspin

  if( nspin == 1 ) then
    write(stdout,*) ' 3. local potential'
  else
    write(stdout,*) ' 3. local potential: spin channel ', ispin
  endif

  do ip=1,nproc
    i=0 ; j=0
    do ij=1,(nbnd*(nbnd+1))/2
      i=i+1
      if( i>j ) then
        j=j+1; i=1
      endif
      if( ij <= starts(ip) .or. ij > starts(ip)+lengths(ip) ) cycle

      ijl=ij-starts(ip)
      if( ijl==1 .or. i==1 ) then
        !
        CALL start_clock( 'firstfft' )
        !
        psic(:) = ( 0.D0, 0.D0 )
        psic(dffts%nl(igk_k(1:npw,1))) = evc(1:npw,ibnd_indx(j))
        CALL invfft ('Wave', psic, dffts)
        !
        CALL stop_clock( 'firstfft' )
        !
        ! ... product with the potential vrs = (vltot+v) on the smooth grid
        !
        psic(1:dffts%nnr) = psic(1:dffts%nnr) * vrs(1:dffts%nnr,ispin)
        !
        ! ... back to reciprocal space
        !
        CALL start_clock( 'secondfft' )
        !
        CALL fwfft ('Wave', psic, dffts)
        !
        ! ... store with correct ordering
        !
        jtmp(1:npw) = conjg(psic(dffts%nl(igk_k(1:npw,1))))
        !
        CALL stop_clock( 'secondfft' )
        !
      endif
      hamtmp(ijl) = ZDOTU( npw, evc(1,ibnd_indx(i)), 1, jtmp, 1 )
    enddo

    call mp_sum( hamtmp, intra_pool_comm )
    if( mpime==ip-1 ) then
      ham = conjg( hamtmp(1:nbnd_l) )
    endif
  enddo

  call write_hamloc( fposnloc, nbnd, ham )

  enddo ! ispin

  ! ======================================================================
  ! deallocate space
  ! ======================================================================
  deallocate( ham, gtmp, stat=ierr )
  deallocate( ibnd_indx, jtmp )


  ! just in case we have a local potential only
  ! then we have no need for projectors, overlaps, etc
  ncpp=.true.
  ! if we have projectors
  if( nkb > 0 ) then
    ! ======================================================================
    ! nonlocal projections
    ! ======================================================================
    write(stdout,*) ' 4. nonlocal projections'

    ! the index mapping from becp or vkb to atoms and projections 
    call update_atomic_proj( nkb, nat, ntyp, ityp, nhm, nh, nhtol, local_channel )

    ! the atomic matrix for the nonlocal potential, i.e. D
    allocate( vnl_atom(nat,nspin) )

    ! if any of the pseudopotentials are Vanderbilt style then we need S
    ! also let hamq_shirley know through ncpp=.false.
    if( okvan ) then
      allocate( usq_atom(nat) )

      ncpp=.false.
      write(stdout,*)
      write(stdout,*) ' ====================================================='
      write(stdout,*) ' Note that I have found Vanderbilt pseudopotentials'
      write(stdout,*) ' I am setting ncpp = ', ncpp
      write(stdout,*) ' This means we solve a generalized eigenproblem H-eS=0'
      write(stdout,*) ' ====================================================='
      write(stdout,*)
    else
      ! use this flag to reduce the diagonalization cost for NCPP's
      ncpp=.true.
      write(stdout,*)
      write(stdout,*) ' ====================================================='
      write(stdout,*) ' Note that I can find no Vanderbilt pseudopotentials'
      write(stdout,*) ' I am setting ncpp = ', ncpp
      write(stdout,*) ' This means we solve a simple eigenproblem H-eI=0'
      write(stdout,*) ' ====================================================='
      write(stdout,*)
    endif

    if( lda_plus_u ) then
      ! the atomic matrix for the Hubbard potential
      allocate( vhU_atom(nat,nspin) )
    endif


    do ispin=1,nspin

      if( nspin == 1 ) then
        write(stdout,*) '    atomic matrices'
      else
        write(stdout,*) '    atomic matrices : spin channel ', ispin
      endif

      do na=1,nat

        nt = ityp(na)
        allocate( vnl_atom(na,ispin)%matrix( nh(nt), nh(nt) ) ) 

        vnl_atom(na,ispin)%matrix = 0.d0
        if( .not. upf(nt)%tvanp ) then
          if( .not. newpseudo(nt) ) then
            write(stdout,*) ' inverting definition of deeq'
            ! if ncpp we know that deeq is diagonal, independent of atomic position,
            ! and needs to be inverted to be consistent with USPPs
            ! I won't reduce the size from nat to ntyp for now to make for easier
            ! code hopping between US and NC cases.
            do ih=1,nh(nt)
              if( abs(deeq(ih,ih,na,ispin)) < 1.d-12 ) then
                call errore('hamq','seems like a diagonal element of deeq for this pseudopotential is zero',1)
              endif
              vnl_atom(na,ispin)%matrix(ih,ih) = 1.d0 / deeq(ih,ih,na,ispin)
            enddo
          else
            ! if we are using the newpseudo (atomic ld1.x) then we may have
            ! non-zero off-diagonal elements of deeq
            do jh=1,nh(nt)
              do ih=1,nh(nt)
                if( abs(deeq(ih,jh,na,ispin)) < 1.d-12 ) cycle
                vnl_atom(na,ispin)%matrix(ih,jh)=deeq(ih,jh,na,ispin)
              enddo
            enddo
          endif
        else
          do jh=1,nh(nt)
            do ih=1,nh(nt)
              vnl_atom(na,ispin)%matrix(ih,jh)=deeq(ih,jh,na,ispin)
            enddo
          enddo
        endif

      enddo ! nat
    enddo ! ispin

    if( okvan ) then
      do na=1,nat
        nt = ityp(na)
        allocate( usq_atom(na)%matrix( nh(nt), nh(nt) ) )
        usq_atom(na)%matrix = 0.d0
      
        if( upf(nt)%tvanp ) then
          do jh=1,nh(nt)
            do ih=1,nh(nt)
              usq_atom(na)%matrix(ih,jh)=qq_at(ih,jh,na)
            enddo
          enddo
        endif
      enddo
    endif

    if( lda_plus_u ) then
      call update_atomic_proj_ldaU( natomwfc, nat, ntyp, ityp, Hubbard_lmax_, Hubbard_l_, Hubbard_U_, Hubbard_alpha_ )

      ! transfer the Hubbard potential atomic matrices
      do na=1,nat
        nt = ityp(na)
        IF (Hubbard_U(nt).NE.0.d0 .OR. Hubbard_alpha(nt).NE.0.d0) THEN
          DO ispin = 1, nspin
            allocate( vhU_atom(na,ispin)%matrix( 2*Hubbard_l(nt)+1, 2*Hubbard_l(nt)+1 ) ) 
            DO m2 = 1, 2 * Hubbard_l(nt) + 1
              DO m1 = 1, 2 * Hubbard_l(nt) + 1
                vhU_atom(na,ispin)%matrix(m1,m2) = v%ns(m1,m2,ispin,na)
              ENDDO
            ENDDO
          ENDDO
        else
          do ispin=1,nspin
            allocate( vhU_atom(na,ispin)%matrix(0,0) )
          enddo 
        ENDIF
      enddo
    endif

    ! split over bands to reduce memory overhead for projectors
    call mp_scatter_size( nbnd, nbnd_l, root )
    write(stdout,*) ' number of basis functions ', nbnd, ' split over processors leaving locally ', nbnd_l, ' basis functions'
    allocate( becp_l(nkb,nbnd_l) )

    ! k-point grid - uniform only
    write(stdout,*) ' generate nonlocal projectors on the following k-point grid:'
    write(stdout,'(2(2x,a,3i6))') 'nkgrid =', nkgrid(1:3), 'ikgrid =', ikgrid(1:3)

    call init_vnl_kgrid( nkgrid, ikgrid, nkb, nbnd_l, nkr )

    allocate( xk_cart(3,nkr) )
    call get_vnl_kgrid_cart( xk_cart )
    xk_cart = xk_cart / tpiba


    ! k-point loop for non-local potential

    call allocate_bec_type( nkb, ceiling(dble(nbnd)/dble(nproc)), becp )

    do ik=1,nkr

      write(stdout,*) ' load non-local potential for ik = ', ik
      write(stdout,'(x,a,3f12.5)') ' xk_cart = ', xk_cart(1:3,ik)

      call init_us_2_shirley( npw, igk_k(1,1), xk_cart(1:3,ik), vkb )

      becp_l = zero
      ibnd=0
      do ip=1,nproc
        if( mpime==(ip-1) ) nbnd_ip=nbnd_l
        call mp_bcast( nbnd_ip, (ip-1), intra_pool_comm )

        CALL ZGEMM( 'C', 'N', nkb, nbnd_ip, npw, (1.0_DP,0.0_DP), &
           vkb, npwx, evc(1,band_subset(1)+ibnd), npwx, (0.0_DP,0.0_DP), becp%k, nkb )
        ! collect results at process ip-1
        call mp_root_sum( becp%k(:,1:nbnd_ip), becp_l, ip-1, world_comm )
        ibnd=ibnd+nbnd_ip
      enddo

      ! if we have a norm-conserving pseudopotential then a different convention
      ! has been adopted for the normalization of the projector |beta> 
      ! which is actually equal to |chi>=dV|phi> rather than the 
      ! normalized version D|chi> where D is 1/<phi|dV|phi>
      ! This should have been fixed in the conversion to UPF format but it isn't
      !   ( note: becp(j,n) = <beta_j|nk> )
      if( .not. all(upf(1:ntyp)%tvanp) ) then
        do na=1,nat
          nt=ityp(na)
          if( .not. upf(nt)%tvanp .and. .not. newpseudo(nt) ) then
            ! note that the matrix for vnl_atom is diagonal for ncpp
            do ih=1,nh(nt)
              j = index_betaq(ih,na)
              becp_l(j,1:nbnd_l) = becp_l(j,1:nbnd_l) / vnl_atom(na,1)%matrix(ih,ih)
            enddo
          endif
        enddo
      endif

      ! store the projectors for this k-point
      call init_vnl_k( nkb, nbnd_l, ik, becp_l )

    enddo  ! ik
    write(stdout,*) 'init_vnl_k for all k'

    ! ask for default spline orders
    kxord=ksplord(1)
    kyord=ksplord(2)
    kzord=ksplord(3)
    write(stdout,*) ' set up splines with the following orders'
    write(stdout,*) kxord, kyord, kzord

    call init_vnl_spline( kxord, kyord, kzord )
    write(stdout,*) 'init_vnl_spline'

    fposnprj=0
    call write_hamprj( fposnprj )
    write(stdout,*) 'write_hamprj'

    ! ======================================================================
    ! deallocate space
    ! ======================================================================
    deallocate( becp_l )
    call deallocate_bec_type( becp )

  endif  ! if nkb > 0

  ! dump the Hamiltonian dimensions and pseudopotential matrices to disk
  if( mpime==root ) then
    iunhq = freeunit()
    call seqopn(iunhq,'hamq','unformatted',exst)

    iunvnl = freeunit()
    call seqopn(iunvnl,'vnla','unformatted',exst)

    write(stdout,*)
    write(stdout,*) ' 5. dumping Hamiltonian dimensions to file ', trim(prefix) // '.hamq'
    call write_hamq( iunhq )

    write(stdout,*)
    write(stdout,*) ' 6. dumping atomic contribution to non-local potential to file ', trim(prefix) // '.vnla'
    do ispin=1,nspin
      write(stdout,*) ' ispin = ', ispin
      call write_nloper( iunvnl, vnl_atom(1:nat,ispin) )
    enddo
    if( okvan ) then
      call write_nloper( iunvnl, usq_atom )
    endif

    if( lda_plus_u ) then
      write(stdout,*) ' 7. dumping atomic contribution to LDA+U to file ', trim(prefix) // '.vnla'
 
      do ispin=1,nspin
        write(stdout,*) ' ispin = ', ispin
        call write_nloper( iunvnl, vhU_atom(1:nat,ispin) )
      enddo
    endif

    close( iunhq )
    close( iunvnl )
  endif

  ! new  
  ! if LDA+U repeat the whole procedure to accumulate the atomic wave functions in the optimal basis
  ! new
  if( lda_plus_u ) then

    write(stdout,*) ' 8. atomic projections for LDA+U'

    ! davegp
    ALLOCATE( swfcatom( npwx*npol, natomwfc) )
    !IF ( lda_plus_u .AND. (U_projection.NE.'pseudo') ) &
    !   ALLOCATE( wfcU(npwx*npol, nwfcU) )

    write(stdout,*) '  natomwfc = ', natomwfc
    write(stdout,*) ' natomproj = ', natomproj

    call allocate_bec_type( natomproj, ceiling(dble(nbnd)/dble(nproc)), proj )
    allocate( becp_l(natomproj,nbnd_l) )
    call resize_vnl( nkr, natomproj, nbnd_l )

    ! I need to reset the number of projectors for the dimensions of LDA+U
    ! if we want to use init_vnl_k and init_vnl_spline
    nproj = natomproj

    ! make a temporary array for swfcatom
    allocate( swfcatom_ldaU(size(swfcatom,1),natomproj) )

    do ik=1,nkr

      write(stdout,*) ' load Hubbard projections for ik = ', ik
      write(stdout,'(x,a,3f12.5)') ' xk_cart = ', xk_cart(1:3,ik)

      call orthoatwfc_shirley( npw, igk_k(1,1), xk_cart(1:3,ik) )

      ! reduce swfcatom to only those elements necessary for LDA+U (i.e., with U neq 0)
      ikb=0
      do na=1,nat
        nt=ityp(na)
        if( Hubbard_U(nt) /= 0.d0 .or. Hubbard_alpha(nt) /= 0.d0 ) then
          do m1=1,2*Hubbard_l(nt)+1
            swfcatom_ldaU(:,ikb+m1) = swfcatom(:,oatwfc(na)+m1)
          enddo
          ikb=ikb + 2*Hubbard_l(nt)+1
        endif
      enddo

      ! distribute the computation of the projections in the eigenbasis (like for vkb above)
      becp_l = zero
      ibnd=0
      do ip=1,nproc
        if( mpime==(ip-1) ) nbnd_ip=nbnd_l
        call mp_bcast( nbnd_ip, (ip-1), intra_pool_comm )

        CALL ZGEMM( 'C', 'N', ikb, nbnd_ip, npw, (1.0_DP,0.0_DP), &
           swfcatom_ldaU, size(swfcatom_ldaU,1), evc(1,band_subset(1)+ibnd), size(evc,1), (0.0_DP,0.0_DP), proj%k, size(proj%k,1) )
        ! collect results at process ip-1
        call mp_root_sum( proj%k(1:ikb,1:nbnd_ip), becp_l, ip-1, world_comm )
        ibnd=ibnd+nbnd_ip
      enddo

      ! store the projectors for this k-point
      call init_vnl_k( natomproj, nbnd_l, ik, becp_l )

    enddo

    ! ask for default spline orders
    kxord=ksplord(1)
    kyord=ksplord(2)
    kzord=ksplord(3)
    write(stdout,*) ' set up splines with the following orders'
    write(stdout,*) kxord, kyord, kzord

    call init_vnl_spline( kxord, kyord, kzord )

    ! append the projectors
    call write_hamprj()

    deallocate( swfcatom_ldaU )
    deallocate( becp_l )
    call deallocate_bec_type( proj )
    
    deallocate( swfcatom )

  endif  ! if lda_plus_u


  ! close out the projector file
  call close_hamq()

  write(stdout,*) ' Done with hamq'
  return

  end subroutine hamq
! ---------------------------------------------------------------------- 

