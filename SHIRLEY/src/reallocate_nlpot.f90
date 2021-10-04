!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine reallocate_nlpot
  !-----------------------------------------------------------------------
  !
  ! This routine computes the dimension of the Hamiltonian matrix and
  ! allocates arrays containing the non-local part of the pseudopotential
  !
  ! It computes the following global quantities:
  !
  !     ngk           !  number of plane waves (for each k point)
  !     npwx          !  maximum number of plane waves
  !     nqx           !  number of points of the interpolation table
  !     nqxq          !  as above, for q-function interpolation table
  !
  !
  USE ions_base,        ONLY : nat, nsp, ityp
  USE cell_base,        ONLY : tpiba2
  USE cellmd,           ONLY : cell_factor
  USE gvect,            ONLY : ngm, gcutm, g
  use gvecw,            only : ecutwfc
  USE klist,            ONLY : xk, wk, ngk, nks, qnorm, igk_k
  USE lsda_mod,         ONLY : nspin
  USE ldaU,             ONLY : Hubbard_lmax
  USE scf,              ONLY : rho
  USE noncollin_module, ONLY : noncolin
  USE wvfct,            ONLY : npwx, npw, g2kin
  USE us,               ONLY : qrad, tab, tab_d2y, tab_at, dq, nqx, &
                               nqxq, spline_ps
  USE uspp,             ONLY : indv, nhtol, nhtolm, ijtoh, qq_at, dvan, deeq, vkb, &
                               nkb, nkbus, nhtoj, becsum, qq_so, dvan_so, deeq_nc
  USE uspp_param,       ONLY : upf, lmaxq, lmaxkb, nh, nhm, nbetam
  USE spin_orb,         ONLY : lspinorb, fcoef
  USE paw_variables,    ONLY : okpaw
  !
  implicit none
  !
  !    a few local variables
  !
  integer :: nwfcm  
  ! counters on atom type, atoms, beta functions
  !
  !   calculate number of PWs for all kpoints
  !
  if( allocated( ngk ) ) deallocate( ngk ); allocate (ngk( nks ))
  !
  call n_plane_waves (ecutwfc, tpiba2, nks, xk, g, ngm, npwx, ngk)
  !
  !   igk_k relates the index of PW k+G to index in the list of G vector
  !
  if( allocated( igk_k ) ) deallocate( igk_k ) ; allocate (igk_k( npwx, nks ) )
  if( allocated( g2kin ) ) deallocate( g2kin ); allocate ( g2kin ( npwx ) )    
  !
  ! Note: computation of the number of beta functions for
  ! each atomic type and the maximum number of beta functions
  ! and the number of beta functions of the solid has been
  ! moved to init_run.f90 : pre_init()
  !
  if( allocated(indv) ) deallocate(indv); allocate (indv( nhm, nsp))    
  if( allocated(nhtol) ) deallocate(nhtol); allocate (nhtol(nhm, nsp))    
  if( allocated(nhtolm) ) deallocate(nhtolm); allocate (nhtolm(nhm, nsp))    
  if( allocated(nhtoj) ) deallocate(nhtoj); allocate (nhtoj(nhm, nsp))    
  if( allocated(ijtoh) ) deallocate(ijtoh); allocate (ijtoh(nhm, nhm, nsp))
  if( allocated(deeq) ) deallocate(deeq); allocate (deeq( nhm, nhm, nat, nspin))    
  if (noncolin) then
     if( allocated(deeq_nc) ) deallocate(deeq_nc); allocate (deeq_nc( nhm, nhm, nat, nspin))    
  endif
  if( allocated(qq_at) ) deallocate(qq_at); allocate (qq_at(   nhm, nhm, nsp))    
  if (lspinorb) then
    if( allocated(qq_so) ) deallocate(qq_so); allocate (qq_so(nhm, nhm, 4, nsp))    
    if( allocated(dvan_so) ) deallocate(dvan_so); allocate (dvan_so( nhm, nhm, nspin, nsp))    
    if( allocated(fcoef) ) deallocate(fcoef); allocate (fcoef(nhm,nhm,2,2,nsp))
  else
    if( allocated(dvan) ) deallocate(dvan); allocate (dvan( nhm, nhm, nsp))    
  endif
  !
  ! This routine is called also by the phonon code, in which case it should
  ! allocate an array that includes q+G vectors up to |q+G|_max <= |Gmax|+|q|
  !
  nqxq = INT( ( (sqrt(gcutm) + qnorm ) / dq + 4) * cell_factor )
  lmaxq = 2*lmaxkb+1
  !
  if (lmaxq > 0) then
    if( allocated(qrad) ) deallocate(qrad); allocate (qrad( nqxq, nbetam*(nbetam+1)/2, lmaxq, nsp))    
  endif
  if (nkb > 0) then
    if( allocated(vkb) ) deallocate(vkb); allocate (vkb( npwx,  nkb))    
  endif
  if( allocated(becsum) ) deallocate(becsum); allocate (becsum( nhm * (nhm + 1)/2, nat, nspin))    
  !
  ! Calculate dimensions for array tab (including a possible factor
  ! coming from cell contraction during variable cell relaxation/MD)
  !
  nqx = INT( (sqrt (ecutwfc) / dq + 4) * cell_factor )

  if( allocated(tab) ) deallocate(tab); allocate (tab( nqx , nbetam , nsp))

  ! d2y is for the cubic splines
  if (spline_ps) then
    if( allocated(tab_d2y) ) deallocate(tab_d2y); allocate (tab_d2y( nqx , nbetam , nsp))
  endif

  nwfcm = MAXVAL ( upf(1:nsp)%nwfc )
  if( allocated(tab_at) ) deallocate(tab_at); allocate (tab_at( nqx , nwfcm , nsp))

  return
end subroutine reallocate_nlpot

