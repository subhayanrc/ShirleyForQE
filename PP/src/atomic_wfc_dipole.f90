!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE atomic_wfc_dipole (ik, wfcatom, natom, atom_index, &
                              nl_tab, l_tab, tab)
  !-----------------------------------------------------------------------
  !
  ! This routine is modeled after atomic_wfc in PW/src
  ! modified to take as input an updated table (tab) of atomic dipole
  ! for a range of l values
  !
  ! This routine computes the superposition of atomic wavefunctions
  ! for k-point "ik" - output in "wfcatom"
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : tpi, fpi, pi
  USE cell_base,  ONLY : omega, tpiba
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau
  USE basis,      ONLY : natomwfc
  USE gvect,      ONLY : mill, eigts1, eigts2, eigts3, g
  USE klist,      ONLY : xk, igk_k, ngk
  USE wvfct,      ONLY : npwx
  USE us,         ONLY : tab_at, dq, nqx
  USE uspp_param, ONLY : upf
  USE noncollin_module, ONLY : noncolin, npol, angle1, angle2
  USE spin_orb,   ONLY : lspinorb, rot_ylm, fcoef, lmaxx, domag, &
                         starting_spin_angle
  USE mp_bands,   ONLY : inter_bgrp_comm
  USE mp,         ONLY : mp_sum
  !
  implicit none
  !
  integer, intent(in) :: ik
  complex(DP), intent(out) :: wfcatom (npwx, npol, natomwfc)
  integer,intent(in) :: natom, atom_index(natom)
  integer,intent(in) :: nl_tab, l_tab(3)
  real(dp),intent(in) :: tab(nqx,3)
  !
  integer :: iatom, ltab
  !
  integer :: n_starting_wfc, lmax_wfc, nt, l, nb, na, m, lm, ig, iig, &
             i0, i1, i2, i3, npw
  real(DP), allocatable :: qg(:), ylm (:,:), chiq (:,:), gk (:,:)
  complex(DP), allocatable :: sk (:), aux(:)
  complex(DP) :: kphase, lphase
  real(DP) :: arg, px, ux, vx, wx
  integer :: ig_start, ig_end

  call start_clock ('atomic_wfc_dipole')

  ! calculate max angular momentum required in wavefunctions
  lmax_wfc = MAX ( 0, MAXVAL (l_tab(1:nl_tab)) )
  !
  npw = ngk(ik)
  !
  allocate ( ylm (npw,(lmax_wfc+1)**2), chiq(npw,nl_tab), &
             sk(npw), gk(3,npw), qg(npw) )
  !
  do ig = 1, npw
     iig = igk_k (ig,ik)
     gk (1,ig) = xk(1, ik) + g(1,iig)
     gk (2,ig) = xk(2, ik) + g(2,iig)
     gk (3,ig) = xk(3, ik) + g(3,iig)
     qg(ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  enddo
  !
  !  ylm = spherical harmonics
  !
  call ylmr2 ((lmax_wfc+1)**2, npw, gk, qg, ylm)

  ! from now to the end of the routine the ig loops are distributed across bgrp
  call divide(inter_bgrp_comm,npw,ig_start,ig_end)
  !
  ! set now q=|k+G| in atomic units
  !
  do ig = ig_start, ig_end
     qg(ig) = sqrt(qg(ig))*tpiba
  enddo
  !
  n_starting_wfc = 0
  !
  ! chiq = radial fourier transform of atomic orbitals chi
  !
  do ltab = 1, nl_tab
           do ig = ig_start, ig_end
              px = qg (ig) / dq - int (qg (ig) / dq)
              ux = 1.d0 - px
              vx = 2.d0 - px
              wx = 3.d0 - px
              i0 = INT( qg (ig) / dq ) + 1
              i1 = i0 + 1
              i2 = i0 + 2
              i3 = i0 + 3
              chiq (ig, ltab) = &
                     tab (i0, ltab) * ux * vx * wx / 6.d0 + &
                     tab (i1, ltab) * px * vx * wx / 2.d0 - &
                     tab (i2, ltab) * px * ux * wx / 2.d0 + &
                     tab (i3, ltab) * px * ux * vx / 6.d0
           enddo
  enddo

  deallocate (qg, gk)
  allocate ( aux(npw) )
  !
  wfcatom(:,:,:) = (0.0_dp, 0.0_dp)
  !
  do iatom = 1, natom
     na = atom_index(iatom)
     arg = (xk(1,ik)*tau(1,na) + xk(2,ik)*tau(2,na) + xk(3,ik)*tau(3,na)) * tpi
     kphase = CMPLX(cos (arg), - sin (arg) ,kind=DP)
     !
     !     sk is the structure factor
     !
     do ig = ig_start, ig_end
        iig = igk_k (ig,ik)
        sk (ig) = kphase * eigts1 (mill (1,iig), na) * &
                           eigts2 (mill (2,iig), na) * &
                           eigts3 (mill (3,iig), na)
     enddo
     !
     do ltab=1,nl_tab
           l=l_tab(ltab)
           lphase = (0.d0,1.d0)**l
           !
           !  the factor i^l MUST BE PRESENT in order to produce
           !  wavefunctions for k=0 that are real in real space
           !
           IF ( noncolin ) THEN
              !
              IF ( upf(nt)%has_so ) THEN
                 !
                 IF (starting_spin_angle.OR..not.domag) THEN
!                    call atomic_wfc_dipole_so ( )
                 ELSE
!                    call atomic_wfc_dipole_so_mag ( )
                 ENDIF
                 !
              ELSE
                 !
!                 call atomic_wfc_dipole_nc ( )
                 !
              ENDIF
              !
           ELSE
              !
              call atomic_wfc_dipole__ ( )
              !
           END IF
           !
     END DO
     !
  END DO

  if (n_starting_wfc /= natomwfc) call errore ('atomic_wfc', &
       'internal error: some wfcs were lost ', 1)

  deallocate(aux, sk, chiq, ylm)

  ! collect results across bgrp
  call mp_sum(wfcatom, inter_bgrp_comm)

  call stop_clock ('atomic_wfc_dipole')
  return

CONTAINS

   SUBROUTINE atomic_wfc_dipole__( )
   !
   ! ... LSDA or nonmagnetic case
   !
   DO m = 1, 2 * l + 1
      lm = l**2 + m
      n_starting_wfc = n_starting_wfc + 1
      if (n_starting_wfc > natomwfc) call errore &
         ('atomic_wfc_dipole__', 'internal error: too many wfcs', 1)
      !
      do ig = ig_start, ig_end
         wfcatom (ig, 1, n_starting_wfc) = lphase * &
            sk (ig) * ylm (ig, lm) * chiq (ig, ltab)
      ENDDO
      !
   END DO
   !
   END SUBROUTINE atomic_wfc_dipole__
   !
END SUBROUTINE atomic_wfc_dipole
