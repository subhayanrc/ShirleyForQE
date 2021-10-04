  module shirley_input_module

! All the variables associated with initializing shirley 
! and routines that generate them.

!#include "f_defs.h"

  use kinds, only : dp
  use shirley_constants, only : maxchar
  use kpt_module
  use corerepair_module
  use plot_module, only : plotspec_type, plotspec_read, plotspec_bcast
  use hamq_shirley

  implicit none

  integer,parameter :: stdin =5

  type(kpt_type) :: kpt, qpt
  type(corerepair_type) :: corerep
  type(plotspec_type) :: plotspec
  type(matrix_list),allocatable :: vnl_atom(:,:)
  type(matrix_list),allocatable :: usq_atom(:)
  type(matrix_list),allocatable :: vhU_atom(:,:)

  character(maxchar) :: prefix
  character(maxchar) :: outdir

  character(maxchar) :: hamqfile
  character(maxchar) :: vnlfile
  character(maxchar) :: pwifile
  character(maxchar) :: pwmfile
  character(maxchar) :: pltfile
  character(maxchar) :: outfile
  character(maxchar) :: pseudo_dir

  integer :: iunhq, iunvnl, iunout, iunplt
  logical :: iunout_opened
  character(maxchar) :: iunout_name
  character(maxchar) :: tmpfile
  real(dp) :: epsilon_cutoff
  real(dp) :: efermi, temperature, nelec
  character(maxchar) :: smearing
  logical :: readqpts, readcorerep, readplotspec
  character(maxchar) :: trifile
  character(maxchar) :: int0file
  character(maxchar) :: vxcfile

  character(maxchar) :: epsfile
  character(maxchar) :: eps00file
  character(maxchar) :: epsijfile
  character(maxchar) :: epsi0file
  character(maxchar) :: eps0jfile
  character(maxchar) :: epsqfile

  real(dp) :: deltaq(3), epstol
  integer :: band_subset(2)
  integer :: nener_sigma
  real(dp) :: dener_sigma

  real(dp) :: memfac

  integer :: niter_equil_mc, niter_samp_mc
  real(dp) :: mc_step_radius
  real(dp) :: spec_min, spec_max, spec_sigma
  character(maxchar) :: spec_smearing

  real(dp) :: omega_resp, eta_resp
  character(maxchar) :: plot_style

  integer :: nspline_fit
  character(maxchar) :: spline_fit(10)

  integer :: niter_stieltjes
  integer :: nener_stieltjes
  real(dp) :: eta_stieltjes, emin_stieltjes, emax_stieltjes

  logical :: kinetic_only, local_only, nonlocal_only, smatrix_only
  logical :: eigvec_output
  logical :: proj_output

  logical :: debug

  character(len=3) :: nodenumber

  contains

  subroutine shirley_input

  use io_global, only : stdout, ionode, ionode_id
  use mp, only : mp_bcast
  use mp_world, only : mpime, world_comm
  !use pwmat_module, only : pwi_read, pwm_open
  use hamq_pool

  integer :: ierr

  character(maxchar) :: outfilenode, fileprefix
  integer :: iatom, it, np, icore, ispin

  integer,external :: freeunit

  real(dp) :: nelec_, alat, omega, at(3,3), bg(3,3), tpiba
  integer :: nspin
  logical :: lda_plus_u

  namelist / input / prefix, outdir, &
                     hamqfile, vnlfile, outfile, pwifile, pwmfile, &
                     pltfile, readplotspec, &
                     epsilon_cutoff, efermi, temperature, nelec, smearing, &
                     tmpfile, readqpts, &
                     readcorerep, pseudo_dir, &
                     trifile, int0file, vxcfile, &
                     epsfile, eps00file, epsijfile, epsi0file, eps0jfile, epsqfile, &
                     deltaq, epstol, &
                     band_subset, nener_sigma, dener_sigma, &
                     niter_equil_mc, niter_samp_mc, mc_step_radius, &
                     spec_min, spec_max, spec_sigma, spec_smearing, &
                     nspline_fit, spline_fit, &
                     omega_resp, eta_resp, plot_style, &
                     memfac, &
                     niter_stieltjes, nener_stieltjes, eta_stieltjes, &
                     emin_stieltjes, emax_stieltjes, &
                     kinetic_only, local_only, nonlocal_only, smatrix_only, &
                     eigvec_output, &
                     proj_output, &
                     debug



! Initialize for shirley_epsilon

  ! Initialize MPI
  call start_shirley( nodenumber )

  ! Modify stdout
  call init_stdout( stdout )

  ! read name list
  if( ionode ) then
    call input_from_file ()

    ! default value
    prefix = 'pwscf_opt'
    outdir = './'
    debug = .false.
    epsilon_cutoff = 0.d0
    efermi = 0.d0
    temperature = 0.d0
    nelec = -1.d0
    smearing = 'fermi-dirac'
    spec_smearing = 'lorentzian'
    tmpfile = 'tmpfile'
    readqpts = .false.
    readcorerep = .false.
    pwifile = ''
    pwmfile = ''
    outfile = 'shirley'
    pltfile = ''
    band_subset = 0
    nener_sigma = 3
    dener_sigma = 1.0
    nspline_fit = 0
    spline_fit = ''
    omega_resp = 0.d0
    eta_resp = 0.1d0
    plot_style = 'cube'
    memfac = 0.5d0
    kinetic_only=.false.
    local_only=.false.
    nonlocal_only=.false.
    smatrix_only=.false.
    eigvec_output=.false.
    proj_output=.false.
    ierr=0
    read(stdin,nml=input,iostat=ierr)
    if( ierr /= 0 ) &
      call errore('shirley_input','problem reading namelist &input',abs(ierr))
  endif
  ! broadcast some input
  call mp_bcast( debug, ionode_id, world_comm )

  call mp_bcast( prefix, ionode_id, world_comm )
  call mp_bcast( outdir, ionode_id, world_comm )

  call mp_bcast( pltfile, ionode_id, world_comm )
  call mp_bcast( pwmfile, ionode_id, world_comm )
  call mp_bcast( outfile, ionode_id, world_comm )
  call mp_bcast( tmpfile, ionode_id, world_comm )
  call mp_bcast( epsilon_cutoff, ionode_id, world_comm )

  call mp_bcast( efermi, ionode_id, world_comm )
  call mp_bcast( temperature, ionode_id, world_comm )
  call mp_bcast( nelec, ionode_id, world_comm )
  call mp_bcast( smearing, ionode_id, world_comm )

  call mp_bcast( readqpts, ionode_id, world_comm )
  call mp_bcast( readcorerep, ionode_id, world_comm )
  call mp_bcast( readplotspec, ionode_id, world_comm )
  call mp_bcast( pseudo_dir, ionode_id, world_comm )
  call mp_bcast( trifile, ionode_id, world_comm )
  call mp_bcast( int0file, ionode_id, world_comm )
  call mp_bcast( vxcfile, ionode_id, world_comm )

  call mp_bcast( epsfile, ionode_id, world_comm )
  call mp_bcast( eps00file, ionode_id, world_comm )
  call mp_bcast( epsijfile, ionode_id, world_comm )
  call mp_bcast( eps0jfile, ionode_id, world_comm )
  call mp_bcast( epsi0file, ionode_id, world_comm )
  call mp_bcast( epsqfile, ionode_id, world_comm )

  call mp_bcast( deltaq, ionode_id, world_comm )
  call mp_bcast( epstol, ionode_id, world_comm )
  call mp_bcast( band_subset, ionode_id, world_comm )
  call mp_bcast( nener_sigma, ionode_id, world_comm )
  call mp_bcast( dener_sigma, ionode_id, world_comm )

  call mp_bcast( niter_equil_mc, ionode_id, world_comm )
  call mp_bcast( niter_samp_mc, ionode_id, world_comm )
  call mp_bcast( mc_step_radius, ionode_id, world_comm )

  call mp_bcast( spec_min, ionode_id, world_comm )
  call mp_bcast( spec_max, ionode_id, world_comm )
  call mp_bcast( spec_sigma, ionode_id, world_comm )
  call mp_bcast( spec_smearing, ionode_id, world_comm )

  call mp_bcast( omega_resp, ionode_id, world_comm )
  call mp_bcast( eta_resp, ionode_id, world_comm )
  call mp_bcast( plot_style, ionode_id, world_comm )

  call mp_bcast( nspline_fit, ionode_id, world_comm )
  call mp_bcast( spline_fit, ionode_id, world_comm )

  call mp_bcast( niter_stieltjes, ionode_id, world_comm )
  call mp_bcast( nener_stieltjes, ionode_id, world_comm )
  call mp_bcast( eta_stieltjes, ionode_id, world_comm )

  call mp_bcast( memfac, ionode_id, world_comm )

  call mp_bcast( kinetic_only, ionode_id, world_comm )
  call mp_bcast( local_only, ionode_id, world_comm )
  call mp_bcast( nonlocal_only, ionode_id, world_comm )
  call mp_bcast( smatrix_only, ionode_id, world_comm )

  call mp_bcast( eigvec_output, ionode_id, world_comm )
  call mp_bcast( proj_output, ionode_id, world_comm )

  ! check if we have pools
  call hamq_pool_init

  ! ------------------------------------------------------------------------
  ! open files
  ! ------------------------------------------------------------------------
  fileprefix=trim(outdir)//'/'//trim(prefix)
  call open_hamq( fileprefix )

  if( ionode ) then
  ! open file for Hamiltonian input
    iunhq = freeunit()
    hamqfile=trim(fileprefix)//'.hamq'
    open(iunhq,file=trim(hamqfile),form='unformatted',iostat=ierr)
    if( ierr /= 0 ) &
      call errore('shirley_input','problem opening hamq file '//trim(hamqfile),101)

  ! open file for Non-local potential input
    iunvnl = freeunit()
    vnlfile=trim(fileprefix)//'.vnla'
    open(iunvnl,file=trim(vnlfile),form='unformatted',iostat=ierr)
    if( ierr /= 0 ) &
      call errore('shirley_input','problem opening vnla file '//trim(vnlfile),102)
  endif

  ! now start reading files
  write(stdout,*) ' reading hamiltonian ...'
  call read_hamq( iunhq )

  call close_hamq

  call dump_system( nelec_, alat, omega, at, bg, tpiba, nspin, &
                    lda_plus_u )
  write(stdout,*) ' nspin = ', nspin
  write(stdout,*) ' lda_plus_u = ', lda_plus_u
  write(stdout,*) ' natomproj = ', natomproj

  ! read atomic matrix elements for non-local potential
  allocate( vnl_atom(natom,nspin) )
  do ispin=1,nspin
  do iatom=1,natom
    it = type_atom(iatom)
    np = nproj_type(it)
    allocate( vnl_atom(iatom,ispin)%matrix(np,np) )
  enddo
  call read_nloper( iunvnl, vnl_atom(1:natom,ispin) )
  enddo

  if( .not. ncpp ) then
    if( mpime==ionode_id ) write(stdout,*) ' Not a norm-conserving pseudpotential'
    allocate( usq_atom(natom) )
    do iatom=1,natom
      it = type_atom(iatom)
      np = nproj_type(it)
      allocate( usq_atom(iatom)%matrix(np,np) )
    enddo
    call read_nloper( iunvnl, usq_atom )
  endif

  if( lda_plus_u ) then
    if( mpime==ionode_id ) write(stdout,*) ' LDA+U Hamiltonian'
    
    allocate( vhU_atom(natom,nspin) )
    do ispin=1,nspin
      do iatom=1,natom
        it = type_atom(iatom)
        if( Hubbard_U(it) /= 0.d0 .or. Hubbard_alpha(it) /= 0.d0 ) then
          np = 2 * Hubbard_l(it) + 1
          allocate( vhU_atom(iatom,ispin)%matrix(np,np) )
        else
          allocate( vhU_atom(iatom,ispin)%matrix(0,0) )
        endif
      enddo
      call read_nloper( iunvnl, vhU_atom(1:natom,ispin) )
    enddo
  endif


  ! input blocks
  if( ionode ) then
    ! read k-points
    call kpt_read( stdin, kpt, nspin )
    write(stdout,*) 'done reading k-points'

    if( readplotspec ) then
      ! read bands to plot
      call plotspec_read( plotspec, stdin )
      write(stdout,*) 'done reading plot specification'
    endif

    if( readqpts ) then
      ! read q-points
      call kpt_read( stdin, qpt )
      write(stdout,*) 'done reading q-points'
    endif

    if( readcorerep ) then
      ! read atomic matrix elements for core-valence position
      call read_corerepair( stdin, corerep )
      write(stdout,*) 'done reading corerepair'
    endif
  endif
  
  ! broadcast input
  call kpt_bcast( kpt, ionode_id, world_comm )
  if( readplotspec ) then
    call plotspec_bcast( plotspec, ionode_id, world_comm )
  endif
  if( readqpts ) then
    call kpt_bcast( qpt, ionode_id, world_comm )
  endif
  if( readcorerep ) then
    call bcast_corerepair( corerep, ionode_id, world_comm )
  endif


  if( nelec < 0.d0 ) then
    nelec = nelec_
    write(stdout,*) ' no nelec provided - using hamq value instead'
  endif
  write(stdout,*) ' nelec = ', nelec

  ! check k-point units
  call kpt_units( kpt, tpiba )
  if( readqpts ) call kpt_units( qpt, tpiba )

  if( readcorerep ) then
    if( nproc_per_pool > 1 ) then

      ! checks on core-repair
      do icore=1,corerep%ncore
        write(stdout,*) corerep%core(icore)%atom
        write(stdout,*) type_atom_global(corerep%core(icore)%atom)
        write(stdout,*) corerep%core(icore)%species
        write(stdout,*) nproj_type(corerep%core(icore)%species)
        write(stdout,*) corerep%core(icore)%nproj1,corerep%core(icore)%nproj2
      enddo

      if( corerep%nspecies /= ntype ) &
        call errore('shirley_input','core-repair has incorrect no. of species', &
                    abs(corerep%nspecies))
      if( corerep%natom /= natom_global ) &
        call errore('shirley_input','core-repair has incorrect no. of atoms', &
                    abs(corerep%natom))
      do icore=1,corerep%ncore
        iatom = corerep%core(icore)%atom
        it = type_atom_global(iatom)
        np = nproj_type(corerep%core(icore)%species)
        if( it /= corerep%core(icore)%species ) then
          call errore('shirley_input','species mismatch with core-repair',icore)
        endif
        if( np /= corerep%core(icore)%nproj1 .and. &
            np /= corerep%core(icore)%nproj2 ) then
          write(stdout,*) ' no. valence projections from core-repair = ', &
                          corerep%core(icore)%nproj1, &
                          corerep%core(icore)%nproj2
          write(stdout,*) ' no. valence projections for species ', &
                          corerep%core(icore)%species, ' = ', np
          call errore('shirley_input','projector mismatch with core-repair',icore)
        endif
      enddo

    else

      ! checks on core-repair
      write(stdout,*) corerep%nspecies,ntype
      write(stdout,*) corerep%natom,natom
      write(stdout,*) corerep%ncore
      do icore=1,corerep%ncore
        write(stdout,*) corerep%core(icore)%atom
        write(stdout,*) type_atom(corerep%core(icore)%atom)
        write(stdout,*) corerep%core(icore)%species
        write(stdout,*) nproj_type(corerep%core(icore)%species)
        write(stdout,*) corerep%core(icore)%nproj1,corerep%core(icore)%nproj2
      enddo

      if( corerep%nspecies /= ntype ) &
        call errore('shirley_input','core-repair has incorrect no. of species', &
                    abs(corerep%nspecies))
      if( corerep%natom /= natom ) &
        call errore('shirley_input','core-repair has incorrect no. of atoms', &
                    abs(corerep%natom))
      do icore=1,corerep%ncore
        iatom = corerep%core(icore)%atom
        it = type_atom(iatom)
        np = nproj_type(corerep%core(icore)%species)
        if( it /= corerep%core(icore)%species ) then
          call errore('shirley_input','species mismatch with core-repair',icore)
        endif
        if( np /= corerep%core(icore)%nproj1 .and. &
            np /= corerep%core(icore)%nproj2 ) then
          write(stdout,*) ' no. valence projections from core-repair = ', &
                          corerep%core(icore)%nproj1, &
                          corerep%core(icore)%nproj2
          write(stdout,*) ' no. valence projections for species ', &
                          corerep%core(icore)%species, ' = ', np
          call errore('shirley_input','projector mismatch with core-repair',icore)
        endif
      enddo
    endif
  endif

  if( readplotspec ) then
    iunplt = freeunit()
    open(iunplt,file=trim(pltfile),form='unformatted',iostat=ierr)
    if( ierr /= 0 ) &
      call errore('shirley_input','problem opening plot file '//trim(pltfile),104)
  endif

  !if( pwifile /= '' .and. pwmfile /= '' ) then
  !  ! pwmatrix elements
  !  call pwi_read( pwifile, ionode_id, mpime )
  !  call pwm_open( pwmfile )
  !endif

  return
  end subroutine shirley_input


  end module shirley_input_module
