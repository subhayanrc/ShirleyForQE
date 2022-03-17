!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine init_at_dipole( l, nr, r, u, nl_tab, l_tab, tab )
  !-----------------------------------------------------------------------
  !
  ! This routine computes a table with the radial Fourier transform 
  ! of the atomic wavefunctions.
  !
  USE kinds,      ONLY : dp
  USE constants,  ONLY : fpi
  USE cell_base,  ONLY : omega
  USE us,         ONLY : nqx, dq
  USE mp_bands,   ONLY : intra_bgrp_comm
  USE mp,         ONLY : mp_sum
  use io_global,  only : stdout
  !
  implicit none
  !
  integer,intent(in) :: l, nr
  real(dp),intent(in) :: r(nr), u(nr)
  integer,intent(out) :: nl_tab, l_tab(3)
  real(dp),intent(out) :: tab(nqx,3)
  !
  integer :: ltab, i, lmax, lmin, lp
  !
  integer :: iq, ir, startq, lastq, ndm
  !
  real(DP), allocatable :: rab(:), aux (:), vchi (:)
  real(DP) :: vqint, pref, q
  real(dp) :: y, p

  call start_clock ('init_at_dipole')

  write(stdout,*) ' init_at_dipole'
  write(stdout,*) ' l, nr = ', l, nr

  allocate (rab(nr),aux(nr),vchi(nr))

  ! build rab from what we know of r (assuming log-radial)
  ! r(i)=((exp((i-1)*dx)-1)*exp(xmin)/zmesh
  ! rab(i)=[r(i)+exp(xmin)/zmesh]*dx
  ! y=exp(dx)=[r(i+1)-r(i)]/r(i) = (r(3)-r(2))/r(2)
  ! p=exp(xmin)/zmesh=(r(3)-r(2))/(y*y-y)
  ! dx=log(y)
  p = r(3)-r(2)
  y = p/r(2)
  p = p/(y*y-y)
  y = log(y)
  forall( i=1:nr ) rab(i) = (r(i)+p)*y
  !
  ! chiq = radial fourier transform of atomic orbitals chi
  !
  pref = fpi/sqrt(3.d0*omega)
  ! factor of sqrt(3) comes from conversion of dipole to spherical Harmonic
  ! r = Y_1m / sqrt(3)
  ! necessary to compute correctly lda+U projections)
  call divide (intra_bgrp_comm, nqx, startq, lastq)
  ! tab should have size (nqx,3)
  tab(:,:) = 0.d0

  if( l<0 ) call errore('init_at_dipole','l value less than zero',abs(l))
  lmin=max(l-1,0)
  lmax=l+1
  ltab=0
  write(stdout,*) l, lmin, lmax, ltab
  do lp=lmin,lmax
    ltab=ltab+1
    l_tab(ltab)=lp
    write(stdout,*) l_tab(ltab)
    do iq = startq, lastq
      q = dq * (iq - 1)
      call sph_bes (nr, r, q, lp, aux)
      do ir = 1, nr
        vchi(ir) = u(ir) * aux(ir) * r(ir) * r(ir)
      enddo
      call simpson (nr, vchi, rab, vqint)
      tab (iq, ltab) = vqint * pref
    enddo
  enddo
  nl_tab=ltab
  call mp_sum ( tab, intra_bgrp_comm )

  deallocate(rab, aux ,vchi)

  write(stdout,*) ' init_at_dipole'

  call stop_clock ('init_at_dipole')
  return

end subroutine init_at_dipole

