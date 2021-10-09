!
! Copyright (C) 2001-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! davegp - copied from PW/src/allocate_fft
!-----------------------------------------------------------------------
SUBROUTINE deallocate_fft
  !-----------------------------------------------------------------------
  !
  !     This routine *de*allocates memory for FFT-related arrays - IMPORTANT:
  !     routine "data_structure" must be called before it in order to
  !     set the proper dimensions and grid distribution across processors
  !     these dimensions
  !
  USE scf,       ONLY : rho, v, vnew, vltot, vrs, rho_core, rhog_core, &
                        kedtau, destroy_scf_type
  USE noncollin_module, ONLY : pointlist, factlist, r_loc
  USE wavefunctions, ONLY : psic, psic_nc
  IMPLICIT NONE
  !
  !     *DE*Allocate memory for all kind of stuff.
  !
  CALL destroy_scf_type(rho)
  CALL destroy_scf_type(v)
  CALL destroy_scf_type(vnew)
  if( allocated(vltot) ) deallocate(vltot)
  if( allocated(rho_core) ) deallocate(rho_core)
  if( allocated(kedtau) ) deallocate(kedtau)
  if( allocated(rhog_core) ) deallocate(rhog_core)
  if( allocated(psic) ) deallocate(psic)
  if( allocated(vrs) ) deallocate(vrs)

  if( allocated(psic_nc) ) deallocate(psic_nc)
  if( allocated(pointlist) ) deallocate(pointlist)
  if( allocated(factlist) ) deallocate(factlist)
  if( allocated(r_loc) ) deallocate(r_loc)

  RETURN
END SUBROUTINE deallocate_fft
