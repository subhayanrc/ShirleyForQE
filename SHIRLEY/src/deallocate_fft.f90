!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine deallocate_fft
  !-----------------------------------------------------------------------
  !     This routine deallocates the data structure associated to the FFT
  !
!#include "f_defs.h"
  USE gvect,     ONLY : deallocate_gvect
  USE scf,       ONLY : rho, v, vltot, vrs, vnew, rho_core, rhog_core, kedtau
  use scf,       only : destroy_scf_type
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : pointlist, factlist, r_loc, &
      report, i_cons, noncolin
  USE wavefunctions, ONLY : psic, psic_nc
  implicit none
  !
  !     DeAllocate memory for all kind of stuff.
  !
  call deallocate_gvect()
  call destroy_scf_type( rho )
  call destroy_scf_type( v )
  call destroy_scf_type( vnew )
  if( allocated( vltot ) ) deallocate (vltot)    
  if( allocated( vrs ) ) deallocate (vrs)    
  if( allocated( rho_core ) ) deallocate (rho_core)    
  if( allocated( kedtau ) ) deallocate (kedtau)    
  if( allocated( rhog_core ) ) deallocate (rhog_core)    

  if( allocated( psic ) ) deallocate (psic)    
  if (noncolin) then
     if( allocated( psic_nc ) ) deallocate (psic_nc)    
     if (((report.ne.0).or.(i_cons.ne.0)).and.(noncolin)) then
!
! In order to print out local quantities, integrated around the atoms,
! we need the following variables
!
        if( allocated( pointlist ) ) deallocate(pointlist)
        if( allocated( factlist ) ) deallocate(factlist)
        if( allocated( r_loc ) ) deallocate(r_loc)

     endif
  endif

  return
end subroutine deallocate_fft

