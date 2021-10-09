  program efermi

  use parallel_include
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_world, only : nproc, mpime, world_comm, mp_world_end
  use mp, only : mp_bcast, mp_barrier, mp_sum, mp_max, mp_min
  use fermi, only : fermifunc, fermideriv, fermi_energy_sub=>fermi_energy
  use mpio

  implicit none

  real(dp),parameter :: evtory=7.3498649d-2
  real(dp),parameter :: rytoev=1.d0/evtory
  real(dp),parameter :: kelvin2rydberg=6.333630d-6

  character(len=3) :: nodenumber

  integer :: narg

  real(dp) :: kT

  integer :: ik, i, j, k, ij, nij, ideriv, niter, ibnd

  real(dp),allocatable :: eigval(:,:), wg(:,:)
  real(dp),allocatable :: wk(:), kvec(:,:)
  logical :: cartesian
  real(dp) :: ef, cbm, vbm
  integer :: nk_loc
  real(dp),allocatable :: wk_loc(:)

  character(255) :: filename
  character(255) :: ckT
  character(255) :: cnelec

  integer :: iuninf, ierr
  integer :: fheigval, fhdhdk
  character(255) :: eigval_file
  integer :: status(MPI_STATUS_SIZE)
  integer(kind=MPI_OFFSET_KIND) :: offset
  integer :: reqeigval

  integer,external :: freeunit
#ifdef __PGI
  integer,external :: iargc
#endif

  character(255) :: smearing, grid_type

  integer :: nk, nbnd, ncp, nspin, nbasis
  logical :: lda_plus_u
  real(dp) :: nelec, alat, volume, &
              at(3,3), bg(3,3), tpiba, fermi_energy
  namelist /info/ nk, nbnd, ncp, nelec, nbasis, alat, volume, &
                  at, bg, tpiba, fermi_energy, nspin, lda_plus_u
  real(dp) :: nelec_in


  
  ! initialize mpi
  CALL start_shirley (nodenumber)

  if( ionode ) then
    narg = iargc()
  endif
  call mp_bcast( narg, ionode_id, world_comm )

  if( narg /= 2 .and. narg /= 3 ) then
    write(stdout,*) ' usage: efermi temp filename [nelec]'
    call mp_world_end
    stop
  endif

  if( ionode ) then

  call getarg( 1, ckT )
  call getarg( 2, filename )
  if( narg==3 ) then
    call getarg( 3, cnelec )
    read(cnelec,*) nelec_in
  endif

  read(ckT,*) kT

  kT = kT*kelvin2rydberg

  endif

  call mp_bcast( kT, ionode_id, world_comm )
  call mp_bcast( filename, ionode_id, world_comm )

  if( ionode ) then
    iuninf=freeunit()
    open(iuninf,file=trim(filename)//'.info',form='formatted')
    read(iuninf,nml=info)
    allocate( wk(nk), kvec(1:3,nk) )
    read(iuninf,*) cartesian
    read(iuninf,*) grid_type
    read(iuninf,*) wk
    read(iuninf,*) kvec
    close(iuninf)
    write(stdout,*) ' reading info from ', trim(filename)
    write(stdout,*) ' fermi_energy (stored) = ', fermi_energy

    if( narg==3 ) then
      write(stdout,*) ' number of electrons modified based on input:'
      write(stdout,*) ' nelec from Hamiltonian = ', nelec
      nelec=nelec_in
      write(stdout,*) ' nelec from input       = ', nelec
    endif
  endif

  call mp_bcast( nk, ionode_id, world_comm )
  call mp_bcast( nbnd, ionode_id, world_comm )
  call mp_bcast( nelec, ionode_id, world_comm )
  call mp_bcast( alat, ionode_id, world_comm )
  call mp_bcast( volume, ionode_id, world_comm )
  call mp_bcast( at, ionode_id, world_comm )
  call mp_bcast( bg, ionode_id, world_comm )
  call mp_bcast( tpiba, ionode_id, world_comm )
  call mp_bcast( fermi_energy, ionode_id, world_comm )
  call mp_bcast( nspin, ionode_id, world_comm )
  call mp_bcast( lda_plus_u, ionode_id, world_comm )

  if( .not. ionode ) allocate( wk(nk), kvec(3,nk) )
  call mp_bcast( cartesian, ionode_id, world_comm )
  call mp_bcast( wk, ionode_id, world_comm )
  call mp_bcast( kvec, ionode_id, world_comm )

  ef = fermi_energy / rytoev

  ! MPI-IO
  eigval_file=trim(filename)//'.eigval'
  
  call mp_file_open_dp( eigval_file, fheigval, ionode_id, world_comm )

  write(stdout,*) '    running ...'

  nk_loc=0
  do ik=1,nk
    if( mod(ik,nproc)/=mpime ) cycle
    nk_loc = nk_loc + 1
  enddo

! debug
!  write(mpime+100,*) nk_loc
! debug

  allocate( eigval(nbnd,nk_loc), &
            wg(nbnd,nk_loc), wk_loc(nk_loc), stat=ierr )
  if( ierr/=0 ) call errore('allocation error',ierr)

  ! distribute weights
  nk_loc=0
  do ik=1,nk
    if( mod(ik,nproc)/=mpime ) cycle

    nk_loc = nk_loc + 1
    wk_loc(nk_loc) = wk(ik)
  enddo

  nk_loc=0
  do ik=1,nk
    if( mod(ik,nproc)/=mpime ) cycle

    nk_loc = nk_loc + 1

    write(stdout,*) ' reading ik= ', ik, ' of ', nk, ' wt= ', &
                    wk_loc(nk_loc),  kvec(1:3,ik)

    ! read eigenvalues
    offset = (ik-1)*nbnd
    call mpi_file_iread_at( fheigval, offset, &
                           eigval(1,nk_loc), nbnd, &
                           MPI_DOUBLE_PRECISION, reqeigval, ierr )
  enddo

  ! close - otherwise I don't know if I have the data or not
  do ik=1,min(nk,nproc)
    if( mod(ik,nproc)/=mpime ) cycle
    call mp_wait( reqeigval, status, ierr )
  enddo
  call mpi_file_close( fheigval, ierr )

! debug
!  nk_loc=0
!  do ik=1,nk
!    if( mod(ik,nproc)/=mpime ) cycle
!    nk_loc = nk_loc + 1
!    if( cartesian ) then
!      write(mpime+100,'(i,3f)') ik, matmul( transpose(at), kvec(:,ik) )/tpiba
!      write(mpime+100,'(i,3f)') ik, kvec(1:3,ik)/tpiba
!    else
!      write(mpime+100,'(i,3f)') ik, kvec(1:3,ik)
!      write(mpime+100,'(i,3f)') ik, matmul( bg, kvec(:,ik) )
!    endif
!    write(mpime+100,*) eigval(1:nbnd,nk_loc)*rytoev
!  enddo
! debug

  write(stdout,*) '    input Fermi energy ...'
  write(stdout,*) ' Fermi Energy = ', ef*rytoev, ' eV'
  write(stdout,*) ' temperature = ', kT*rytoev, ' eV'
  write(stdout,*) ' temperature = ', kT/kelvin2rydberg, ' K'
  write(stdout,*) ' nelec = ', nelec

  ! determine fermi energy
  smearing='fermi-dirac'
  call fermi_energy_sub( nelec, nk_loc, nbnd, wk_loc, eigval, kT, smearing, ef )

  ! find VBM and CBM if present
  cbm= 9.d+20
  vbm=-9.d+20
  do ik=1,nk_loc
    do ibnd=1,nbnd
      if( eigval(ibnd,ik) >= ef ) exit
    enddo
    cbm = min( cbm, eigval(min(ibnd,nbnd),ik) )
    vbm = max( vbm, eigval(ibnd-1,ik) )
  enddo
  call mp_max( vbm, world_comm )
  call mp_min( cbm, world_comm )

  write(stdout,*) '          VBM = ', vbm*rytoev, ' eV'
  write(stdout,*) '          CBM = ', cbm*rytoev, ' eV'
  write(stdout,*) ' Fermi Energy = ', ef*rytoev, ' eV'

  call mp_barrier( world_comm )
  call mp_world_end
  stop

  end program efermi

