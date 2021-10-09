!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! 
!----------------------------------------------------------------------- 
PROGRAM shirley_overlap 
!----------------------------------------------------------------------- 

  ! Step 1 

  ! David Prendergast
  ! TMF, Sep 20, 2016

  ! Generate the overlap matrix between two Shirley optimal basis sets
  !
  ! The assumption is that the output of two independent shirley_basis.x
  ! calculations provides two optimal basis sets for systems with exactly
  ! the same unit cell and plane-wave kinetic energy cut-offs. Imagine that
  ! each basis set is as follows:
  !   <G|B1_i> and <G|B2_j>
  ! It is not necessary that the number of B1 and B2 functions be the same
  ! but we do require that the number of G vectors is identical in order to
  ! enable overlap matrix elements as dot_products
  !  S_ij = dot_product( B1_i, B2_j ) = sum_G <G|B1_i>* . <G|B2_j>

  ! Step 2

  ! Yufeng Liang
  ! TMF, Jan 2017
  !
  ! Exploit the eigenvectors < B_i | nk > obtained with shirley_xas to find:
  !                                ~
  !                         < nk | n'k > ^ PS
  !                                 ~       ~     ~
  ! = sum_{ij} < nk | B_j > < B_j | B_i > < B_i | n'k >
  !
  ! < B_j | nk > input from eigvec file in the GS directory
  !         ~
  ! < B_i | n'k >      from the excited-atom directory

  ! Step 3
  ! 
  ! Exploit the projectors < beta_l | nk > obtained with shirley_xas to find:
  !                                       ~        ~
  ! \sum_{Ill'} <nk|beta^I_l> Q^I_{ll'} <beta^I_l'|n'k>
  !
  ! <beta^I_l|nk> obtained from the proj file
  !                         ~                         ~
  ! Q^I_{ll'} = <phi^I_l | phi^I_l'>^AE - <phi^I_l | phi^I_l'>^PS
  ! obtained from the preprocessed file Qij.dat.

  ! Summing up contribution from Step 2 and 3, we can get
  !                              ~
  !                       < nk | n'k > ^ AE
  !

  USE kinds,                 ONLY : DP 
  USE parameters,            ONLY : ntypx, npk, lmaxx
  USE io_global,             ONLY : stdout, ionode
  USE io_files,              ONLY : nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc, wfc_dir, postfix
  USE control_flags,         ONLY : lscf, gamma_only, io_level, twfcollect
  USE basis,                 ONLY : starting_pot, starting_wfc, natomwfc

  USE mp_world,              ONLY : world_comm
  USE mp_global,             ONLY : npool, mp_sum
  USE buffers,               ONLY : open_buffer, close_buffer, get_buffer

  USE klist,                 ONLY : nkstot, nks, ngk
  USE noncollin_module,      ONLY : npol
  USE wavefunctions,         ONLY : evc
  USE wvfct,                 ONLY : nbnd, nbndx, npwx, npw

  USE shirley_overlap_input, ONLY : get_input, prefix1, prefix2, &
                                    outdir1, outdir2, nspin_ham, &
                                    nbnd1, nbnd2

  IMPLICIT NONE

  LOGICAL  :: exst
  character(256) :: fout_olp
  integer :: iunolp = 100
  complex(dp),allocatable :: evc1(:,:)
  complex(dp),allocatable :: overlap(:,:)
  complex(dp),parameter :: one=(1.d0,0.d0)
  complex(dp),parameter :: zero=(0.d0,0.d0)
  integer :: npw1, npw2


  CALL start_shirley (nd_nmbr) 

  if( npool /= 1 ) then
    call errore('shirley_overlap', 'number of pools should be 1', abs(npool))
  endif

  call start_clock( 'shirley' )

  IF ( ionode ) THEN

     WRITE( UNIT = stdout, &
            FMT = '(/5X,"Ultrasoft (Vanderbilt) Pseudopotentials")')

     WRITE( unit = stdout, FMT = 9010 ) & 
         ntypx, npk, lmaxx

9010 FORMAT( /,5X,'Current dimensions of program PWSCF are:', &
           & /,5X,'Max number of different atomic species (ntypx) = ',I2,&
           & /,5X,'Max number of k-points (npk) = ',I6,&
           & /,5X,'Max angular momentum in pseudopotentials (lmaxx) = ',i2)

  END IF  

  ! ======================================================================
  ! read input and distribute
  ! ======================================================================

  call get_input

  ! this provides us with two prefix variables (prefix1 and prefix2) and
  ! one outdir where we will dump the overlap

  write(stdout,*) ' generating overlap matrix between optimal bases at'
  write(stdout,*) '   ', trim(outdir1)//trim(prefix1)
  write(stdout,*) '   ', trim(outdir2)//trim(prefix2)
  write(stdout,*)

  ! ======================================================================
  ! Setup flags for a nonselfconsistent run
  ! ======================================================================
 
  lscf = .false.
  starting_pot = 'file'
  starting_wfc = 'file'

  ! ======================================================================
  !   Now allocate space for pwscf variables, read and check them. 
  ! ======================================================================
  ! would this work? just changing the prefix and reading again?
  ! do we need openfil?
  prefix = prefix1
  tmp_dir = trim(outdir1)
  !
  ! ======================================================================
  !   Now allocate space for pwscf variables, read and check them. 
  ! ======================================================================
  call read_file_shirley( nspin_ham )
  !
  write(stdout,*) '           nks = ', nks
  write(stdout,*) '          npwx = ', npwx
  write(stdout,*) '           npw = ', ngk(1:nks)
  write(stdout,*) '          npol = ', npol
  write(stdout,*) '          nbnd = ', nbnd
  write(stdout,*) '     size(evc) = ', size(evc,1), size(evc,2)
  !
  CALL open_buffer ( iunwfc, 'wfc', nwordwfc, io_level, exst )
  CALL get_buffer( evc, nwordwfc, iunwfc, 1 )
  CALL close_buffer ( iunwfc, 'DELETE' )

  ! ======================================================================
  ! copy evc to evc1
  ! ======================================================================
  allocate( evc1(size(evc,1),size(evc,2)) )
  evc1 = evc
  npw1 = ngk(1)
  if ( nbnd1 == 0 ) then 
    nbnd1 = nbnd
  end if
  write(stdout,*) ' nbnd = ', nbnd1, ' npw = ', npw1

  ! ======================================================================
  ! reset
  ! ======================================================================
  deallocate( evc )
  npwx=0
  nbnd=0
  npol=0

  call clean_pw( .true. )

  ! ======================================================================
  ! now read second file
  ! ======================================================================
  prefix=prefix2
  tmp_dir=trim(outdir2)
  !
  call read_file_shirley( nspin_ham )
  !
  write(stdout,*) '           nks = ', nks
  write(stdout,*) '          npwx = ', npwx
  write(stdout,*) '           npw = ', ngk(1:nks)
  write(stdout,*) '          npol = ', npol
  write(stdout,*) '          nbnd = ', nbnd
  write(stdout,*) '     size(evc) = ', size(evc,1), size(evc,2)
  !
  write(stdout,*) ' reading basis 2'
  CALL open_buffer ( iunwfc, 'wfc', nwordwfc, io_level, exst )
  CALL get_buffer( evc, nwordwfc, iunwfc, 1 )
  CALL close_buffer ( iunwfc, 'DELETE' )

  npw2 = ngk(1)
  if ( nbnd2 == 0 ) then
    nbnd2=nbnd
  end if
  write(stdout,*) ' nbnd = ', nbnd2, ' npw = ', npw2

  ! ======================================================================
  ! check size consistency
  ! ======================================================================
  write(stdout,*) ' number of plane-waves ', npw1, npw2, size(evc1,1), size(evc,1)
  if( npw1 /= npw2 ) &
    call errore('shirley_overlap', 'inconsistent number of plane waves', abs(npw1-npw2))

  ! ======================================================================
  ! make the overlap
  ! ======================================================================
  allocate( overlap(nbnd1,nbnd2) )
  write(stdout,*) ' computing overlap'
  call zgemm( 'C', 'N', nbnd1, nbnd2, npw1, one, evc1, size(evc1,1), evc, size(evc,1), &
              zero, overlap, nbnd1 )
  ! reduce
  call mp_sum( overlap, world_comm )

  ! ======================================================================
  ! dump overlap to file as ASCII and binary
  ! ======================================================================
  fout_olp = 'overlap.dat'
  open(iunolp, file=trim(fout_olp), form='formatted')
  write(iunolp,'(2f20.15)') overlap
  close(iunolp)

  fout_olp = 'overlap'
  open(iunolp, file=trim(fout_olp), form='unformatted')
  write(iunolp) overlap
  close(iunolp)
  ! ======================================================================
  ! stop
  ! ======================================================================
  call stop_clock( 'shirley' )
  call stop_shirley

END PROGRAM shirley_overlap 
!
