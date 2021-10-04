!----------------------------------------------------------------------- 
SUBROUTINE update_pseudo( nspecies, pseudo_dir_, pseudo_file )
  !----------------------------------------------------------------------- 
  ! 
  USE io_global,  ONLY : stdout, ionode
  use ions_base, only: ntyp=>nsp
!  USE atom 
!  USE basis,     ONLY : natomwfc
!  USE cell_base 
!  USE constants, ONLY: rytoev, eps4
!  USE gvect 
!  USE klist, ONLY: xk, nks, nkstot, nelec
!  USE ldaU 
!  USE lsda_mod, ONLY: nspin, isk, current_spin
!  USE symme, ONLY: nsym, irt 
!  USE wvfct 
!  USE uspp, ONLY: nkb, vkb
!  USE uspp_param, ONLY : lmaxkb, nbeta, lll, lmaxq, nhm, tvanp, nh
!  USE becmod,   ONLY: becp, rbecp
!  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc 
  USE io_files,   ONLY : pseudo_dir, psfile
  use dfunct,     only : newd
  USE funct,      ONLY : get_dft_name
  use read_pseudo_mod, only : readpp
!  USE wavefunctions_module, ONLY: evc 
  !
  !
  IMPLICIT NONE
  !
  integer,intent(in) :: nspecies
  character (len=255), intent(in) :: pseudo_dir_, pseudo_file(nspecies)
  character (len=20) :: input_dft
  !
  integer :: nt
  !
  write(stdout,*) 'update_pseudo'
  !
  ! This check is important, otherwise there will be a search for
  ! pseudopotential files that don't exist
  if( ntyp /= nspecies ) &
    call errore('update_pseudo','inconsistent number of atomic species',ntyp)
  !
  ! davegp
  ! update any previous definition of pseudo_dir
  pseudo_dir = trim(adjustl(pseudo_dir_))
  write(stdout,*) ' pseudo_dir = ', trim(pseudo_dir)
  do nt=1,nspecies
    psfile(nt) = trim(pseudo_file(nt))
    write(stdout,*) 'psfile ', nt, ' : ', trim(psfile(nt))
  enddo
  input_dft = get_dft_name ()
  call readpp ( input_dft )

!  write(stdout,*) ' just out of readpp'
!  do nt=1,ntyp
!    write(stdout,*) '    type = ', nt
!    write(stdout,*) '   nbeta = ', nbeta(nt)
!    write(stdout,*) '   lmaxkb = ', lmaxkb
!    write(stdout,*) '   lll = ', lll(1:nbeta(nt),nt)
!    write(stdout,*) '   lmaxq = ', lmaxq
!    write(stdout,*) '   nhm = ', nhm
!    write(stdout,*) '   nh = ', nh(nt)
!    write(stdout,*)
!  enddo

  call reallocate_nlpot

!  write(stdout,*) ' just out of reallocate_nlpot'
!  do nt=1,ntyp
!    write(stdout,*) '    type = ', nt
!    write(stdout,*) '   nbeta = ', nbeta(nt)
!    write(stdout,*) '   lmaxkb = ', lmaxkb
!    write(stdout,*) '   lll = ', lll(1:nbeta(nt),nt)
!    write(stdout,*) '   lmaxq = ', lmaxq
!    write(stdout,*) '   nhm = ', nhm
!    write(stdout,*) '   nh = ', nh(nt)
!    write(stdout,*)
!  enddo

  call init_us_1
  CALL newd()
  !
end subroutine update_pseudo
!
