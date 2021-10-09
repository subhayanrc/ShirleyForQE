!----------------------------------------------------------------------- 
SUBROUTINE update_pseudo( nspecies, pseudo_dir_, pseudo_file )
  !----------------------------------------------------------------------- 
  ! 
  USE io_global,  ONLY : stdout
  use ions_base, only: ntyp=>nsp
  USE io_files,   ONLY : pseudo_dir, psfile
  use dfunct,     only : newd
  USE funct,      ONLY : get_dft_name
  use read_pseudo_mod, only : readpp
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

  ! now read the pseudos
  input_dft = get_dft_name ()
  call readpp ( input_dft )

  ! update sizes of data structures based on input pseudos
  call pre_init ()

  ! based on updated sizes, reallocate some arrays that depend on them
  call reallocate_nlpot

  ! repopulate those arrays with data from the pseudos
  call init_us_1
  CALL newd()
  !
end subroutine update_pseudo
!
