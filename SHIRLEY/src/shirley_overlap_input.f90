  module shirley_overlap_input

  use kinds, only : dp

  implicit none

  character(256) :: outdir1, outdir2
  CHARACTER(len=256) :: prefix1, prefix2
  ! flags
  logical :: reduced_io
  ! debug flag
  logical :: debug
  ! spin
  integer :: nspin_ham
  !
  ! nbnd
  integer :: nbnd1, nbnd2
  contains

  !----------------------------------------------------------------------- 
    subroutine get_input
  !----------------------------------------------------------------------- 
  !
  USE io_global,  ONLY : stdout, ionode, ionode_id
  !USE io_files,   ONLY : tmp_dir 
  USE mp,         ONLY : mp_bcast, mp_barrier      
  USE mp_world,   ONLY : world_comm
  !
  INTEGER :: ios
  integer :: i
  integer :: ibuf5(5)
  ! 
  NAMELIST / input / outdir1, outdir2, prefix1, prefix2, &
                     debug, reduced_io, nspin_ham, &
                     nbnd1, nbnd2
  !
  !   set default values for variables in namelist 
  ! 
  prefix1 = 'pwscf_opt1' 
  prefix2 = 'pwscf_opt2' 
  outdir1 = './' 
  outdir2 = './' 
  debug = .false.
  reduced_io = .false.
  nspin_ham = 1 ! default
  nbnd1 = 0
  nbnd2 = 0
  !
  IF ( ionode )  THEN  
     !
     CALL input_from_file ( )
     !
     READ (5, input, err = 200, iostat = ios) 
200  CALL errore ('shirley_overlap', 'reading input namelist', ABS (ios) ) 
     ! 
     call append_missing_slash( outdir1 )
     call append_missing_slash( outdir2 )
     !tmp_dir = TRIM(outdir) 
     !
     if( nspin_ham > 1 ) write(stdout,*) ' This is a spin-polarized Hamiltonian: nspin = ', nspin_ham
     !
     write(stdout,*)
     !
  END IF 
  ! 
  ! ... Broadcast variables 
  ! 
  CALL mp_bcast( outdir1,  ionode_id, world_comm ) 
  CALL mp_bcast( outdir2,  ionode_id, world_comm ) 
  CALL mp_bcast( prefix1,   ionode_id, world_comm ) 
  CALL mp_bcast( prefix2,   ionode_id, world_comm ) 
  CALL mp_bcast( debug, ionode_id, world_comm )
  CALL mp_bcast( reduced_io, ionode_id, world_comm )
  ! 
  CALL mp_bcast( nspin_ham, ionode_id, world_comm )
  CALL mp_bcast( nbnd1, ionode_id, world_comm )
  CALL mp_bcast( nbnd2, ionode_id, world_comm )
  !
  return
  end subroutine get_input
! ---------------------------------------------------------------------- 

  end module shirley_overlap_input

