  module shirley_basis_input

  use kinds, only : dp

  implicit none

  character(255) :: outdir, prefix
  ! k-point grid and pools
  integer :: band_subset(2)
  integer :: ndim(3)
  real(dp) :: trace_tol, ecut_in
  logical :: expand_kpts, debug
  !
  NAMELIST / input / outdir, prefix, trace_tol, ndim, band_subset, &
                     expand_kpts, ecut_in, debug

  contains

!----------------------------------------------------------------------- 
  subroutine get_input
!----------------------------------------------------------------------- 
  !
  use mp, only : mp_bcast
  USE mp_images, ONLY : intra_image_comm
  use io_global, only : ionode, ionode_id
  ! 
  INTEGER :: ios, l
  !
  !   set default values for variables in namelist 
  ! 
  prefix = 'pwscf' 
  outdir = './' 
  trace_tol = -1.d0  ! negative values mean no truncation
  ndim = (/ 1, 1, 1 /)
  band_subset = 0
  expand_kpts = .true.
  ecut_in = 0.d0
  debug = .false.
  !
  IF ( ionode )  THEN  
     !
     CALL input_from_file ( )
     !
     READ (5, nml=input, err = 200, iostat = ios) 
200  CALL errore ('get_input', 'reading inputpp namelist', ABS (ios) ) 
     ! 
     if( any(ndim > 1 .or. ndim < 0 ) ) &
       call errore('shirley_basis','input parameter ndim must be 0 or 1',1)
     !
     call append_missing_slash( outdir )
     !
  END IF 
  ! 
  ! ... Broadcast variables 
  ! 
  CALL mp_bcast( outdir, ionode_id, intra_image_comm ) 
  CALL mp_bcast( prefix,  ionode_id, intra_image_comm ) 
  CALL mp_bcast( trace_tol,  ionode_id, intra_image_comm ) 
  CALL mp_bcast( ndim,  ionode_id, intra_image_comm ) 
  CALL mp_bcast( band_subset,  ionode_id, intra_image_comm ) 
  CALL mp_bcast( expand_kpts,  ionode_id, intra_image_comm ) 
  CALL mp_bcast( ecut_in,  ionode_id, intra_image_comm ) 
  CALL mp_bcast( debug,  ionode_id, intra_image_comm ) 
  ! 
  return
  end subroutine get_input
! ---------------------------------------------------------------------- 

  end module shirley_basis_input
