!
! Copyright (C) 2001-2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE check_energies()
  !----------------------------------------------------------------------------
  !
  ! ... this is a wrapper to specific calls
  !
  ! ... internal procedures :
  !
  ! ... c_bands_gamma()   : for gamma sampling of the BZ (optimized algorithms)
  ! ... c_bands_k()       : for arbitrary BZ sampling (general algorithm)
  ! ... test_exit_cond()  : the test on the iterative diagonalization
  !
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : eps4, rytoev
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : iunigk, nwordatwfc, iunat, iunwfc, nwordwfc, iunefield
  USE cell_base,            ONLY : tpiba2, at
  USE klist,                ONLY : nkstot, nks, wk, xk, nelec
  USE uspp,                 ONLY : vkb, nkb, okvan
  USE gvect,                ONLY : g, gstart, ngm
  USE wvfct,                ONLY : g2kin, wg, nbndx, et, nbnd, npwx, igk, &
                                   npw, current_k, ecutwfc, ecfixed, qcutz, q2sigma
  USE control_flags,        ONLY : diis_ndim, istep, ethr, lscf, max_cg_iter, &
                                   isolve, io_level, tr2
  USE ldaU,                 ONLY : lda_plus_u, swfcatom
  USE scf,                  ONLY : vltot
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE noncollin_module,     ONLY : noncolin, npol
  USE wavefunctions,        ONLY : evc
  USE bp,                   ONLY : lelfield, evcel
  USE mp_global,            ONLY : intra_pool_comm
  use mp, only : mp_sum
  !
  !!!use haydock
  !
  !
  IMPLICIT NONE
  !
  ! ... First the I/O variables
  !
  INTEGER :: ik_, iter
    ! k-point already done
    ! current iterations
  REAL(DP) :: dr2
    ! current accuracy of self-consistency
  !
  ! ... local variables
  !
  REAL(DP) :: avg_iter, cg_iter, v_of_0
    ! average number of iterations
    ! number of iterations in Conjugate-Gradient
    ! the average of the potential
  INTEGER :: ik, ig, ibnd, dav_iter, diis_iter, ntry, notconv
    ! counter on k points
    ! counter on G vectors
    ! counter on bands
    ! number of iterations in Davidson
    ! number of iterations in DIIS
    ! number or repeated call to diagonalization in case of non convergence
    ! number of notconverged elements
  LOGICAL :: lrot
    ! .TRUE. if the wfc have already be rotated
  INTEGER, ALLOCATABLE :: btype(:)
    ! type of band: valence (1) or conduction (0)  
  !
  ! ... external functions
  !
  REAL(DP), EXTERNAL :: erf
    ! error function  
  !
  !
       !-----------------------------------------------------------------------
       !
       ! ... This routine is a driver for the diagonalization routines of the
       ! ... total Hamiltonian at each k-point.
       ! ... It reads the Hamiltonian and an initial guess of the wavefunctions
       ! ... from a file and computes initialization quantities for the
       ! ... diagonalization routines.
       ! ... There are three types of iterative diagonalization:
       ! ... a) Davidson algorithm (all-band)
       ! ... b) Conjugate Gradient (band-by-band)
       ! ... c) DIIS algorithm
       !
       !
       !
       ! ... here the local variables
       !
       INTEGER :: ipol
       complex(dp) :: hpsi(npwx), eval(nbnd), ovrlp
       !
       !
       ! ... For each k point diagonalizes the hamiltonian
       !
       call hinit0
       call potinit
       k_loop: DO ik = 1, nks
          !
          current_k = ik
          !
          IF ( lsda ) current_spin = isk(ik)
          !
          ! ... sets the kinetic energy
          !
          call gk_sort (xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
          g2kin = g2kin * tpiba2
          !
          ! ... various initializations
          !
          IF ( nkb > 0 ) &
             CALL init_us_2( npw, igk, xk(1,ik), vkb )
          !
          ! ... read in wavefunctions from the previous iteration
          !
             IF ( nks > 1 .OR. ( io_level > 1 ) .OR. lelfield) &
                CALL davcio( evc, 2*nwordwfc, iunwfc, ik, -1 )
          !   
          !
          write(stdout,'(i6,3f12.5)') ik, xk(:,ik)
          write(stdout,'(i6,3f12.5)') ik, matmul( transpose(at), xk(:,ik) )
          do ibnd=1,nbnd
            hpsi = 0.d0
            call h_psi( npwx, npw, 1, evc(:,ibnd), hpsi )
            eval(ibnd) = dot_product( evc(1:npw,ibnd), hpsi(1:npw) )
            ovrlp = dot_product( evc(1:npw,ibnd), evc(1:npw,ibnd) )
#ifdef __MPI
            call mp_sum( eval(ibnd), intra_pool_comm )
#endif
            write(stdout,'(2i6,2f12.5)') ik, ibnd, eval(ibnd)/real(ovrlp)*rytoev
          enddo
          !
       END DO k_loop
       !
       !
       RETURN
       !
     !     
END SUBROUTINE check_energies
