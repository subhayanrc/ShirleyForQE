!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! 
!----------------------------------------------------------------------- 
PROGRAM shirley_basis 
  !----------------------------------------------------------------------- 
  ! 
  ! David Prendergast
  ! UCB, Dec 2006
  !
  ! Generates the optimal basis set for Shirley Brillouin zone interpolation
  ! of electronic wave functions for use in the parametrized Shirley Hamiltonian
  !
  USE parameters,         ONLY : ntypx, npk, lmaxx
  USE io_global,  ONLY : stdout, ionode
  USE io_files,   ONLY : prefix, tmp_dir 
  USE mp,         ONLY : mp_rank      
  USE mp_world,   ONLY : world_comm
  USE mp_pools,   ONLY : npool
  USE mp_bands,   ONLY : ntask_groups
  USE basis,                ONLY : starting_pot, starting_wfc
  !
  use shirley_basis_input, only : get_input, &
                                  prefix_input=>prefix, &
                                  tmp_dir_input=>outdir, &
                                  debug
  !
  implicit none
  !
  CALL start_shirley () 
  !
  if( npool /= 1 ) then
    call errore('shirley_basis','number of pools should be 1',abs(npool))
  endif
  !
  if( ntask_groups /= 1 ) then
    call errore('shirley_basis','number of task groups should be 1',abs(ntask_groups))
  endif
  !
  IF ( ionode ) THEN
     !     
     WRITE( UNIT = stdout, &
            FMT = '(/5X,"Ultrasoft (Vanderbilt) Pseudopotentials")')
     !     
     WRITE( unit = stdout, FMT = 9010 ) & 
         ntypx, npk, lmaxx
     !     
9010 FORMAT( /,5X,'Current dimensions of program PWSCF are:', &
           & /,5X,'Max number of different atomic species (ntypx) = ',I2,&
           & /,5X,'Max number of k-points (npk) = ',I6,&
           & /,5X,'Max angular momentum in pseudopotentials (lmaxx) = ',i2)
  !

  END IF  
  !
  ! read input and distribute
  call get_input
  !
  if( debug .and. .not. ionode ) stdout = 500+mp_rank(world_comm)
  !
  prefix = prefix_input
  tmp_dir = tmp_dir_input
  !
  call read_file()
  !
  call openfil
  !
  CALL wfcinit()
  !
  ! build optimal basis by constructing and diagonalizing overlap
  !
  call build_optimal_basis()
  !
  ! stop
  call stop_shirley


END PROGRAM shirley_basis 
!
