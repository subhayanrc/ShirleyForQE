  module shirley_ham_input

  use kinds, only : dp

  implicit none

  character(255) :: calculation
  character(255) :: outdir
  ! band subset to read and process
  integer :: band_subset(2)
  ! three-center band subset to read and process
  integer :: triples_subset(2)
  ! k-point grid and offset
  integer :: nkgrid(3), ikgrid(3)
  ! k-point spline orders
  integer :: ksplord(3)
  ! pseudo potential info
  logical :: ncpp
  logical :: updatepp
  character(255) :: pseudo_dir
  integer :: nspecies
  integer,allocatable :: local_channel(:)
  character(255),allocatable :: pseudo_file(:)
  ! matrix element info
  logical :: pwmtxel, elphmtxel
  integer :: pwmtxel_method
  real(dp) :: pwmtxel_cutoff
  integer :: nrdv1, nrdv2, nrdv3, nspindv, npedv
  character(255) :: dvfile, elphfile, dynfile
  ! plotting
  integer :: nplot
  character(255) :: plot_prefix, plot_style
  integer :: plot_center_atom
  ! flags
  logical :: reduced_io
  ! debug flag
  logical :: debug
  ! spin
  integer :: nspin_ham
  !
  !
  contains

  !----------------------------------------------------------------------- 
    subroutine get_input
  !----------------------------------------------------------------------- 
  !
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE io_files,   ONLY : nd_nmbr, prefix, tmp_dir 
  USE mp_world,   ONLY : world_comm
  USE mp,         ONLY : mp_bcast, mp_barrier      
  !
  INTEGER :: ios
  integer :: i
  integer :: ibuf5(5)
  ! 
  NAMELIST / input / calculation, outdir, prefix, &
                     band_subset, ksplord, updatepp, &
                     pseudo_dir, &
                     debug, ncpp, &
                     pwmtxel_method, pwmtxel_cutoff, &
                     triples_subset, reduced_io, &
                     nrdv1, nrdv2, nrdv3, nspindv, npedv, &
                     dvfile, elphfile, dynfile, &
                     nplot, plot_prefix, plot_style, &
                     plot_center_atom, &
                     nspin_ham
  !
  !   set default values for variables in namelist 
  ! 
  calculation = 'ham'
  prefix = 'pwscf_opt' 
  outdir = './' 
  ksplord = 0
  ncpp = .false.
  updatepp = .false.
  pseudo_dir = './' 
  band_subset = 0
  triples_subset = 0
  pwmtxel_method = 1
  pwmtxel_cutoff = 0.d0
  debug = .false.
  reduced_io = .false.
  nrdv1=0 ; nrdv2=0 ; nrdv3=0 ; nspindv=0 ; npedv=0
  dvfile = ' '
  elphfile = ' '
  dynfile = ' '
  nplot=0
  plot_prefix = trim(prefix)
  plot_style = 'cube'
  plot_center_atom = 0
  nspin_ham = 1 ! default
  !
  IF ( ionode )  THEN  
     !
     CALL input_from_file ( )
     !
     READ (5, input, err = 200, iostat = ios) 
200  CALL errore ('shirley_ham', 'reading input namelist', ABS (ios) ) 
     ! 
     call append_missing_slash( outdir )
     tmp_dir = TRIM(outdir) 
     !
     if( ncpp ) write(stdout,*) ' ncpp = .T. This is a norm-conserving pseudopotential'
     if( nspin_ham > 1 ) write(stdout,*) ' This is a spin-polarized Hamiltonian: nspin = ', nspin_ham
     !
     write(stdout,*)
     !
     ! read k-point list and shift integers
     read(5,*) ! heading
     read(5,*) nkgrid(1:3), ikgrid(1:3)
     write(stdout,*) 'k-point grid for non-local potential:'
     write(stdout,*) nkgrid
     write(stdout,*) ikgrid
     !
     if( updatepp ) then
     ! read updates to pseudopotential files
       read(5,*) ! heading
       read(5,*) nspecies
       allocate( local_channel(nspecies), pseudo_file(nspecies) )
       write(stdout,*) ' updated potentials:'
       do i=1,nspecies
         read(5,*) local_channel(i), pseudo_file(i)
         write(stdout,*) ' l_local = ', local_channel(i), ' pseudo = ', trim(pseudo_file(i))
       enddo
     endif
     !
     ibuf5 = (/ nrdv1, nrdv2, nrdv3, nspindv, npedv /)
  END IF 
  ! 
  ! ... Broadcast variables 
  ! 
  CALL mp_bcast( calculation,  ionode_id, world_comm ) 
  CALL mp_bcast( tmp_dir,  ionode_id, world_comm ) 
  CALL mp_bcast( prefix,   ionode_id, world_comm ) 
  CALL mp_bcast( ksplord,  ionode_id, world_comm )
  CALL mp_bcast( ncpp, ionode_id, world_comm )
  CALL mp_bcast( updatepp, ionode_id, world_comm )
  CALL mp_bcast( pwmtxel_method, ionode_id, world_comm )
  CALL mp_bcast( pwmtxel_cutoff, ionode_id, world_comm )
  CALL mp_bcast( debug, ionode_id, world_comm )
  CALL mp_bcast( reduced_io, ionode_id, world_comm )
  CALL mp_bcast( band_subset, ionode_id, world_comm )
  CALL mp_bcast( triples_subset, ionode_id, world_comm )
  ! 
  CALL mp_bcast( nplot, ionode_id, world_comm )
  CALL mp_bcast( plot_prefix, ionode_id, world_comm )
  CALL mp_bcast( plot_style, ionode_id, world_comm )
  CALL mp_bcast( plot_center_atom, ionode_id, world_comm )
  ! 
  CALL mp_bcast( nkgrid,   ionode_id, world_comm )
  CALL mp_bcast( ikgrid,   ionode_id, world_comm )
  !
  CALL mp_bcast( nspin_ham, ionode_id, world_comm )
  !
  if( updatepp ) then
    CALL mp_bcast( pseudo_dir, ionode_id, world_comm )
    CALL mp_bcast( nspecies, ionode_id, world_comm )
    if( .not. ionode ) allocate( local_channel(nspecies), pseudo_file(nspecies) )
    CALL mp_bcast( local_channel, ionode_id, world_comm )
    CALL mp_bcast( pseudo_file, ionode_id, world_comm )
  endif

  ! elphmtxel
  call mp_bcast( ibuf5, ionode_id, world_comm )
  call mp_bcast( dvfile, ionode_id, world_comm )
  call mp_bcast( elphfile, ionode_id, world_comm )
  call mp_bcast( dynfile, ionode_id, world_comm )
  nrdv1   = ibuf5(1)
  nrdv2   = ibuf5(2)
  nrdv3   = ibuf5(3)
  nspindv = ibuf5(4)
  npedv   = ibuf5(5)

  return
  end subroutine get_input
! ---------------------------------------------------------------------- 

  end module shirley_ham_input

