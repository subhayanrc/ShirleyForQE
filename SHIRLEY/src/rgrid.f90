!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#include "f_defs.h"
!
!---------------------------------------------------------------------------
SUBROUTINE rgrid( r )
  !---------------------------------------------------------------------------
  !
  ! modeled after PW/src/compute_dip.f90 compute_el_dip()
  !
  USE kinds,      ONLY : DP
  USE fft_base,   ONLY : dffts
  !
  IMPLICIT NONE
  !
  REAL(DP)     :: r(dffts%nnr,3)
  !
  INTEGER  :: i, j, k, j0, k0, ir, idx
  REAL(DP) :: inv_nr1, inv_nr2, inv_nr3
  !
  !
  inv_nr1 = 1.D0 / DBLE( dffts%nr1 )
  inv_nr2 = 1.D0 / DBLE( dffts%nr2 )
  inv_nr3 = 1.D0 / DBLE( dffts%nr3 )
  !
  j0 = dffts%my_i0r2p ; k0 = dffts%my_i0r3p
  DO ir = 1, dffts%nr1x*dffts%my_nr2p*dffts%my_nr3p
     !
     ! ... three dimensional indexes
     !
     idx = ir - 1
     k   = idx / (dffts%nr1x*dffts%my_nr2p)
     idx = idx - (dffts%nr1x*dffts%my_nr2p)*k
     k   = k + k0
     j   = idx / dffts%nr1x
     idx = idx - dffts%nr1x*j
     j   = j + j0
     i   = idx
     !
     ! lattice coordinates
     !
     r(ir,1) = DBLE( i )*inv_nr1
     r(ir,2) = DBLE( j )*inv_nr2
     r(ir,3) = DBLE( k )*inv_nr3
     !
  END DO
  !
  RETURN
  !
END SUBROUTINE rgrid
!
