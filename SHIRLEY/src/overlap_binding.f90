!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! 
!----------------------------------------------------------------------- 
PROGRAM overlap_binding 
!----------------------------------------------------------------------- 

  ! David Prendergast
  ! TMF, OCT 18, 2021

  ! Generate the overlap matrix between two SCF calculations and estimate
  ! binding energies of the unoccupied space based on the Bethe-Salpeter
  ! equation, as follows:
  !
  ! E_b(f) = sum_c A_cf* (e_c) A_cf / sum_c A_cf* A_cf - e~_f
  !
  ! The code reads two SCFs - the first is assumed to be the initial state
  ! containing e_c above, the second is the final state containing e~_f.
  ! The overlap matrix elements can be used to define the exciton 
  ! eigenvectors: A_cf = <c|f~> between initial orbital c and final 
  ! orbital f~. 
  !
  ! written based on shirley_overlap.x but using read_file instead of
  ! read_file_shirley.

  USE kinds,                 ONLY : DP 
  USE parameters,            ONLY : ntypx, npk, lmaxx
  USE constants,             ONLY : rytoev
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
  USE wvfct,                 ONLY : wg, et
  USE ener,                  ONLY : ef
  USE wvfct,                 ONLY : nbnd, nbndx, npwx, npw

  USE shirley_overlap_input, ONLY : get_input, prefix1, prefix2, &
                                    outdir1, outdir2, nspin_ham, &
                                    nbnd1, nbnd2

  IMPLICIT NONE

  LOGICAL  :: exst
  character(256) :: fout_olp
  integer :: iunolp = 100
  complex(dp),allocatable :: evc1(:,:)
  real(dp),allocatable :: et1(:,:)
  complex(dp),allocatable :: overlap(:,:)
  complex(dp),parameter :: one=(1.d0,0.d0)
  complex(dp),parameter :: zero=(0.d0,0.d0)
  integer :: npw1, npw2
  integer :: i, ifermi, ifermi1
  real(dp) :: ebi, o, osum, norm


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
  call read_file( )
  !
  write(stdout,*) '           nks = ', nks
  write(stdout,*) '          npwx = ', npwx
  write(stdout,*) '           npw = ', ngk(1:nks)
  write(stdout,*) '          npol = ', npol
  write(stdout,*) '          nbnd = ', nbnd
  write(stdout,*) '     size(evc) = ', size(evc,1), size(evc,2)
  !
  write(stdout,*) '      size(et) = ', size(et,1), size(et,2)
  write(stdout,*) '   eigenvalues = ', et(1:nbnd,1)*rytoev
  !
  ! find Fermi band index
  do i=1,nbnd
    if( et(i,1) > ef ) exit
  enddo
  ifermi = i
  !
  CALL open_buffer ( iunwfc, 'wfc', nwordwfc, io_level, exst )
  CALL get_buffer( evc, nwordwfc, iunwfc, 1 )
  CALL close_buffer ( iunwfc, 'DELETE' )

  ! ======================================================================
  ! copy evc to evc1
  ! ======================================================================
  allocate( evc1(size(evc,1),size(evc,2)) )
  evc1 = evc
  allocate( et1(size(et,1),size(et,2)) )
  et1 = et
  ifermi1=ifermi
  npw1 = ngk(1)
  if ( nbnd1 == 0 ) then 
    nbnd1 = nbnd
  end if
  write(stdout,*) ' nbnd = ', nbnd1, ' npw = ', npw1

  ! ======================================================================
  ! reset
  ! ======================================================================
  !deallocate( evc )
  call clean_pw( .true. )

  npwx=0
  nbnd=0
  npol=0

  ! ======================================================================
  ! now read second file
  ! ======================================================================
  prefix=prefix2
  tmp_dir=trim(outdir2)
  !
  call read_file( )
  !
  write(stdout,*) '           nks = ', nks
  write(stdout,*) '          npwx = ', npwx
  write(stdout,*) '           npw = ', ngk(1:nks)
  write(stdout,*) '          npol = ', npol
  write(stdout,*) '          nbnd = ', nbnd
  write(stdout,*) '     size(evc) = ', size(evc,1), size(evc,2)
  !
  write(stdout,*) '      size(et) = ', size(et,1), size(et,2)
  write(stdout,*) '   eigenvalues = ', et(1:nbnd,1)*rytoev
  !
  ! find Fermi band index
  do i=1,nbnd
    if( et(i,1) > ef ) exit
  enddo
  ifermi = i
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

  write(stdout,*) ' unoccupied indices: ', ifermi1, ifermi
  write(stdout,'(a6,a6,4a16)') 'band', 'occ', 'SCF2 eigval', &
                  'binding energy', 'exciton CB norm', 'full norm'
  do i=ifermi,nbnd
    norm = dot_product( overlap(:,i), overlap(:,i) )
    osum=sum( overlap(ifermi1:nbnd1,i)*conjg(overlap(ifermi1:nbnd1,i)) )
    ebi = sum( et1(ifermi1:nbnd1,1) * overlap(ifermi1:nbnd1,i)*conjg(overlap(ifermi1:nbnd1,i)) )/osum - et(i,1)
    write(stdout,'(i6,f6.2,4f16.8)') i, wg(i,1), et(i,1)*rytoev, ebi*rytoev, osum, norm
  enddo 
  !
  goto 999
  ! I don't need to write the overlap
  ! ======================================================================
  ! dump overlap to file as ASCII and binary
  ! ======================================================================
  !fout_olp = 'overlap.dat'
  !open(iunolp, file=trim(fout_olp), form='formatted')
  !write(iunolp,'(2f20.15)') overlap
  !close(iunolp)

  !fout_olp = 'overlap'
  !open(iunolp, file=trim(fout_olp), form='unformatted')
  !write(iunolp) overlap
  !close(iunolp)
  ! ======================================================================
  ! stop
  ! ======================================================================
 999 continue
  call stop_clock( 'shirley' )
  call stop_shirley

END PROGRAM overlap_binding 
!
