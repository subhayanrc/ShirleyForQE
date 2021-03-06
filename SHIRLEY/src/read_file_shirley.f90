!
! Copyright (C) 2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! davegp - copied from PW/src/read_file_new.f90
!----------------------------------------------------------------------------
SUBROUTINE read_file_shirley( nspin_ham )
  !----------------------------------------------------------------------------
  !
  ! Read data produced by pw.x or cp.x - new xml file and binary files
  ! Wrapper routine for backwards compatibility
  !
  USE io_files,             ONLY : nwordwfc, iunwfc, prefix, tmp_dir, wfc_dir
  USE io_global,            ONLY : stdout, ionode
  USE buffers,              ONLY : open_buffer, close_buffer
  USE wvfct,                ONLY : nbnd, npwx
  USE noncollin_module,     ONLY : npol
  USE paw_variables,        ONLY : okpaw, ddd_PAW
  USE paw_onecenter,        ONLY : paw_potential
  USE uspp,                 ONLY : becsum
  USE scf,                  ONLY : rho
  USE realus,               ONLY : betapointlist, &
                                   init_realspace_vars,real_space
  USE dfunct,               ONLY : newd
  USE ldaU,                 ONLY : lda_plus_u, U_projection
  USE pw_restart_new,       ONLY : read_collected_to_evc
  USE control_flags,        ONLY : twfcollect
  USE io_files,             ONLY : tmp_dir, prefix, postfix
  USE control_flags,        ONLY : io_level
  USE klist,                ONLY : init_igk
  USE gvect,                ONLY : ngm, g
  USE gvecw,                ONLY : gcutw
  !
  USE lsda_mod,             ONLY : nspin
  !
  IMPLICIT NONE 

  integer, intent(in) :: nspin_ham

  INTEGER :: ierr
  LOGICAL :: exst
  CHARACTER( LEN=256 )  :: dirname
  !
  !
  ierr = 0 
  !
  ! ... Read the contents of the xml data file
  !
  dirname = TRIM( tmp_dir ) // TRIM( prefix ) // postfix
  IF ( ionode ) WRITE( stdout, '(/,5x,A,/,5x,A)') &
     'Reading data from directory:', TRIM( dirname )
  !
  CALL read_xml_file_shirley ( nspin_ham )
  !
  ! ... Open unit iunwfc, for Kohn-Sham orbitals - we assume that wfcs
  ! ... have been written to tmp_dir, not to a different directory!
  ! ... io_level = 1 so that a real file is opened
  !
  wfc_dir = tmp_dir
  nwordwfc = nbnd*npwx*npol
  io_level = 1
  CALL open_buffer ( iunwfc, 'wfc', nwordwfc, io_level, exst )
  !
  ! ... Allocate and compute k+G indices and number of plane waves
  ! ... FIXME: should be read from file, not re-computed
  !
  CALL init_igk ( npwx, ngm, g, gcutw ) 
  !
  ! davegp - momentarily set the variable nspin=1 for reading optimal basis
  if( nspin_ham > 1 ) then
    nspin=1
  endif
  !
  IF ( twfcollect )  CALL read_collected_to_evc ( TRIM ( dirname )) 
  !
  ! now reset nspin
  if( nspin_ham > 1 ) then
    nspin=nspin_ham
  endif
  !
  ! ... Assorted initialization: pseudopotentials, PAW
  ! ... Not sure which ones (if any) should be done here
  !
  CALL init_us_1()
  !
  IF (lda_plus_u .AND. (U_projection == 'pseudo')) CALL init_q_aeps()
  !
  IF (okpaw) THEN
     becsum = rho%bec
     CALL PAW_potential(rho%bec, ddd_PAW)
  ENDIF 
  !
  IF ( real_space ) THEN
    CALL betapointlist()
    CALL init_realspace_vars()
    IF( ionode ) WRITE(stdout,'(5x,"Real space initialisation completed")')
  ENDIF
  CALL newd()
  !
  CALL close_buffer  ( iunwfc, 'KEEP' )
  !
END SUBROUTINE read_file_shirley
!
!----------------------------------------------------------------------------
SUBROUTINE read_xml_file_shirley ( nspin_ham )
  !----------------------------------------------------------------------------
  !
  ! ... This routine allocates space for all quantities already computed
  ! ... in the pwscf program and reads them from the data file.
  ! ... All quantities that are initialized in subroutine "setup" when
  ! ... starting from scratch should be initialized here when restarting
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, nsp, ityp, tau, extfor
  USE cell_base,            ONLY : tpiba2, alat,omega, at, bg, ibrav
  USE force_mod,            ONLY : force
  USE klist,                ONLY : nkstot, nks, xk, wk
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE wvfct,                ONLY : nbnd, nbndx, et, wg
  USE symm_base,            ONLY : irt, d1, d2, d3, checkallsym, nsym
  USE extfield,             ONLY : forcefield, tefield, gate, forcegate
  USE cellmd,               ONLY : cell_factor, lmovecell
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE fft_types,            ONLY : fft_type_allocate
  USE recvec_subs,          ONLY : ggen, ggens
  USE gvect,                ONLY : gg, ngm, g, gcutm, mill, ngm_g, ig_l2g, &
                                   eigts1, eigts2, eigts3, gstart
  USE fft_base,             ONLY : dfftp, dffts
  USE gvecs,                ONLY : ngms, gcutms 
  USE spin_orb,             ONLY : lspinorb, domag
  USE scf,                  ONLY : rho, rho_core, rhog_core, v
  USE wavefunctions, ONLY : psic
  USE vlocal,               ONLY : strf
  USE io_files,             ONLY : tmp_dir, prefix, iunpun, nwordwfc, iunwfc
  USE noncollin_module,     ONLY : noncolin, npol, nspin_lsda, nspin_mag, nspin_gga
  USE pw_restart_new,       ONLY :  pw_readschema_file, init_vars_from_schema 
  USE qes_types_module,     ONLY :  output_type, parallel_info_type, general_info_type, input_type
  USE qes_libs_module,      ONLY :  qes_reset_output, qes_reset_input, qes_reset_general_info, qes_reset_parallel_info 
  USE io_rho_xml,           ONLY : read_scf
  USE fft_rho,              ONLY : rho_g2r
  USE read_pseudo_mod,      ONLY : readpp
  USE uspp,                 ONLY : becsum
  USE uspp_param,           ONLY : upf
  USE paw_variables,        ONLY : okpaw, ddd_PAW
  USE paw_init,             ONLY : paw_init_onecenter, allocate_paw_internals
  USE ldaU,                 ONLY : lda_plus_u, eth, init_lda_plus_u
  USE control_flags,        ONLY : gamma_only
  USE funct,                ONLY : get_inlc, get_dft_name
  USE kernel_table,         ONLY : initialize_kernel_table
  USE esm,                  ONLY : do_comp_esm, esm_init
  USE mp_bands,             ONLY : intra_bgrp_comm, nyfft
  USE Coul_cut_2D,          ONLY : do_cutoff_2D, cutoff_fact 
#if defined(__BEOWULF)
  USE io_global,             ONLY : ionode, ionode_id
  USE bcast_qes_types_module,ONLY : qes_bcast 
  USE mp_images,             ONLY : intra_image_comm
#endif
  !
  use io_global, only : stdout
  !
  IMPLICIT NONE

  integer,intent(in) :: nspin_ham

  INTEGER  :: i, is, ik, ibnd, nb, nt, ios, isym, ierr, inlc
  REAL(DP) :: rdum(1,1), ehart, etxc, vtxc, etotefield, charge
  REAL(DP) :: sr(3,3,48)
  CHARACTER(LEN=20) dft_name
  TYPE ( output_type)                   :: output_obj 
  TYPE (parallel_info_type)             :: parinfo_obj
  TYPE (general_info_type )             :: geninfo_obj
  TYPE (input_type)                     :: input_obj
  !
  !
#if defined(__BEOWULF)
   IF (ionode) THEN
      CALL pw_readschema_file ( ierr, output_obj, parinfo_obj, geninfo_obj, input_obj)
      IF ( ierr /= 0 ) CALL errore ( 'read_schema', 'unable to read xml file', ierr ) 
   END IF
   CALL qes_bcast(output_obj, ionode_id, intra_image_comm)
   CALL qes_bcast(parinfo_obj, ionode_id, intra_image_comm)
   CALL qes_bcast(geninfo_obj, ionode_id, intra_image_comm) 
   CALL qes_bcast(input_obj, ionode_id, intra_image_comm)
#else
  CALL pw_readschema_file ( ierr, output_obj, parinfo_obj, geninfo_obj, input_obj)
  IF ( ierr /= 0 ) CALL errore ( 'read_schema', 'unable to read xml file', ierr ) 
#endif
  ! ... first we get the version of the qexml file
  !     if not already read
  !
  ! ... here we read the variables that dimension the system
  !
  CALL init_vars_from_schema( 'dim',   ierr , output_obj, parinfo_obj, geninfo_obj )
  CALL errore( 'read_xml_file ', 'problem reading file ' // &
             & TRIM( tmp_dir ) // TRIM( prefix ) // '.save', ierr )
  !
  CALL init_vars_from_schema( 'boundary_conditions',   ierr , output_obj, parinfo_obj, geninfo_obj )
  CALL errore( 'read_xml_file ', 'problem reading file ' // &
             & TRIM( tmp_dir ) // TRIM( prefix ) // '.save', ierr )
  ! ... allocate space for atomic positions, symmetries, forces
  !
  IF ( nat < 0 ) CALL errore( 'read_xml_file', 'wrong number of atoms', 1 )
  !
  ! ... allocation
  !
  ALLOCATE( ityp( nat ) )
  ALLOCATE( tau(    3, nat ) )
  ALLOCATE( force(  3, nat ) )
  ALLOCATE( extfor(  3, nat ) )
  !
  IF ( tefield ) ALLOCATE( forcefield( 3, nat ) )
  IF ( gate ) ALLOCATE( forcegate( 3, nat ) ) ! TB
  !
  ALLOCATE( irt( 48, nat ) )
  !
  CALL set_dimensions()
  CALL fft_type_allocate ( dfftp, at, bg, gcutm, intra_bgrp_comm, nyfft=nyfft )
  CALL fft_type_allocate ( dffts, at, bg, gcutms, intra_bgrp_comm, nyfft=nyfft )
  !
  ! ... check whether LSDA
  !
  IF ( lsda ) THEN
     !
     nspin = 2
     npol  = 1
     !
  ELSE IF ( noncolin ) THEN
     !
     nspin        = 4
     npol         = 2
     current_spin = 1
     !
  ELSE
     !
     nspin        = 1
     npol         = 1
     current_spin = 1
     !
  END IF
  !
  if (cell_factor == 0.d0) cell_factor = 1.D0
  !
  ! ... allocate memory for eigenvalues and weights (read from file)
  !
  nbndx = nbnd
  ALLOCATE( et( nbnd, nkstot ) , wg( nbnd, nkstot ) )
  !
  ! ... here we read all the variables defining the system
  !
  CALL init_vars_from_schema ( 'nowave', ierr, output_obj, parinfo_obj, geninfo_obj, input_obj )
  !
  ! ... distribute across pools k-points and related variables.
  ! ... nks is defined by the following routine as the number 
  ! ... of k-points in the current pool
  !
  CALL divide_et_impera( nkstot, xk, wk, isk, nks )
  !
  CALL poolscatter( nbnd, nkstot, et, nks, et )
  CALL poolscatter( nbnd, nkstot, wg, nks, wg )
  !
  ! ... check on symmetry
  !
  IF (nat > 0) CALL checkallsym( nat, tau, ityp)
  !
  ! davegp - major edit with respect to original
  !
  if( nspin_ham /= 0 ) then
  !
  write(stdout,*) ' nspin_ham = ', nspin_ham
  !
  ! davegp - adjusting spin now to true value: nspin_ham
  !          The input basis funcitons have no spin index, but the Hamiltonian does
  !
  if( nspin_ham > 1 ) then
    nspin = nspin_ham
    write(stdout,*)
    write(stdout,*) ' resetting the spin dimension according to input '
    write(stdout,*) ' nspin = ', nspin
    write(stdout,*)
  endif
  !
    SELECT CASE( nspin )
  CASE( 1 )
     !
     lsda = .false.
     IF ( noncolin ) nspin = 4
     !
  CASE( 2 )
     !
     lsda = .true.
     IF ( noncolin ) &
        CALL errore( 'iosys', &
                     'noncolin .and. nspin==2 are conflicting flags', 1 )
     !
  CASE( 4 )
     !
     lsda = .false.
     noncolin = .true.
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys', 'wrong input value for nspin', 1 )
     !
  END SELECT
  !
  ! ... check whether LSDA
  !
  IF ( lsda ) THEN
     !
     nspin = 2
     npol  = 1
     !
  ELSE IF ( noncolin ) THEN
     !
     nspin        = 4
     npol         = 2
     current_spin = 1
     !
  ELSE
     !
     nspin        = 1
     npol         = 1
     current_spin = 1
     !
  END IF
  !
  ! very important check
  !
  write(stdout,*) '        nspin = ', nspin
  write(stdout,*) '         npol = ', npol
  write(stdout,*) ' current_spin = ', current_spin
  !
  endif ! end of if-block for nspin_ham /= 0
  !
  !  Set the different spin indices
  !
  nspin_mag  = nspin
  nspin_lsda = nspin
  nspin_gga  = nspin
  IF (nspin==4) THEN
     nspin_lsda=1
     IF (domag) THEN
        nspin_gga=2
     ELSE
        nspin_gga=1
        nspin_mag=1
     ENDIF
  ENDIF
  !
  ! ... read pseudopotentials
  !
  CALL init_vars_from_schema ( 'pseudo', ierr, output_obj, parinfo_obj, geninfo_obj ) 
  !
  dft_name = get_dft_name () ! already set, should not be set again
  CALL readpp ( dft_name )
  !
  ! ... read the vdw kernel table if needed
  !
  inlc = get_inlc()
  if (inlc > 0 ) then
      call initialize_kernel_table(inlc)
  endif
  !
  okpaw = ANY ( upf(1:nsp)%tpawp )
  !
  IF ( .NOT. lspinorb ) CALL average_pp ( nsp )
  !
  ! ... allocate memory for G- and R-space fft arrays
  !
  CALL pre_init()
  CALL data_structure ( gamma_only )
  CALL allocate_fft()
  CALL ggen ( dfftp, gamma_only, at, bg, gcutm, ngm_g, ngm, &
       g, gg, mill, ig_l2g, gstart ) 
  CALL ggens( dffts, gamma_only, at, g, gg, mill, gcutms, ngms ) 
  IF (do_comp_esm) THEN
     CALL init_vars_from_schema ( 'esm', ierr, output_obj, parinfo_obj, geninfo_obj ) 
     CALL esm_init()
  END IF
  CALL gshells ( lmovecell ) 
  !
  ! ... allocate the potential and wavefunctions
  !
  CALL allocate_locpot()
!  CALL allocate_nlpot()
! davegp - this should work and also extend the size of tab_at
  CALL reallocate_nlpot()
  IF (okpaw) THEN
     CALL allocate_paw_internals()
     CALL paw_init_onecenter()
     CALL d_matrix(d1,d2,d3)
  ENDIF
  !
  IF ( lda_plus_u ) THEN
     CALL init_lda_plus_u ( upf(1:nsp)%psd, noncolin )
     CALL init_at_1()
  ENDIF
  !
  CALL allocate_wfc()
  !
  ! ... read the charge density
  !
  CALL read_scf( rho, nspin, gamma_only )
  ! FIXME: for compatibility. rho was previously read and written in real space
  ! FIXME: now it is in G space - to be removed together with old format
  CALL rho_g2r ( dfftp, rho%of_g, rho%of_r )
  !
  ! ... re-calculate the local part of the pseudopotential vltot
  ! ... and the core correction charge (if any) - This is done here
  ! ... for compatibility with the previous version of read_file
  !
  !2D calculations: re-initialize cutoff fact before calculating potentials
  IF(do_cutoff_2D) CALL cutoff_fact()
  !
  CALL init_vloc()
  CALL struc_fact( nat, tau, nsp, ityp, ngm, g, bg, dfftp%nr1, dfftp%nr2, &
                   dfftp%nr3, strf, eigts1, eigts2, eigts3 )
  CALL setlocal()
  CALL set_rhoc()
  !
  ! ... bring rho to G-space
  !
  DO is = 1, nspin
     !
     psic(:) = rho%of_r(:,is)
     CALL fwfft ('Rho', psic, dfftp)
     rho%of_g(:,is) = psic(dfftp%nl(:))
     !
  END DO
  !
  ! ... recalculate the potential
  !
  CALL v_of_rho( rho, rho_core, rhog_core, &
                 ehart, etxc, vtxc, eth, etotefield, charge, v )
  !
  !
  CALL qes_reset_output ( output_obj )  
  CALL qes_reset_general_info ( geninfo_obj ) 
  CALL qes_reset_parallel_info ( parinfo_obj ) 
  IF ( TRIM(input_obj%tagname) == "input") CALL qes_reset_input ( input_obj) 
  ! 
  RETURN
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE set_dimensions()
      !------------------------------------------------------------------------
      !
      USE constants, ONLY : pi, eps8
      USE cell_base, ONLY : alat, tpiba, tpiba2
      USE gvect,     ONLY : ecutrho, gcutm
      USE gvecs,     ONLY : gcutms, dual, doublegrid
      USE gvecw,     ONLY : gcutw, ecutwfc
      !
      !
      ! ... Set the units in real and reciprocal space
      !
      tpiba  = 2.D0 * pi / alat
      tpiba2 = tpiba**2
      !
      ! ... Compute the cut-off of the G vectors
      !
      gcutw =        ecutwfc / tpiba2
      gcutm = dual * ecutwfc / tpiba2
      ecutrho=dual * ecutwfc
      !
      doublegrid = ( dual > 4.0_dp + eps8 )
      IF ( doublegrid ) THEN
         gcutms = 4.D0 * ecutwfc / tpiba2
      ELSE
         gcutms = gcutm
      END IF
      !
    END SUBROUTINE set_dimensions
    !
  END SUBROUTINE read_xml_file_shirley
