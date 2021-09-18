  module wfc_shirley

  ! this module contains details of the Fourier space wave functions
  ! used in expanding the optimal Shirley basis set

  ! David Prendergast

  use kinds, only : dp

  USE wavefunctions, ONLY : evc
  use shirley_basis_input, only : debug

  implicit none

  public

  ! additional array used to store wave functions from different k-points
  complex(dp),allocatable :: evc_2(:,:)
  integer :: nwordwfc_gamma, iunwfc_gamma

  ! G-space details at the Gamma-point
  integer :: npw_gamma, nbnd_gamma
  integer :: igwx, igwx_gamma
  real(dp) :: ecutwfc_gamma, gcutw_gamma
  integer,allocatable :: igk_gamma(:)
  real(dp),allocatable :: g2kin_gamma(:)
  integer,allocatable :: igk_l2g(:,:)
  integer,allocatable :: igk_l2g_gamma(:)

  ! inverse mapped to unmapped
  type kmap
    integer :: nmap
    real(dp),pointer :: xmap(:,:)
  end type kmap
  type(kmap),allocatable :: xk_map(:)

  complex(dp),parameter :: zero = cmplx(0.d0,0.d0)
  complex(dp),parameter :: one  = cmplx(1.d0,0.d0)


  contains


! ----------------------------------------------------------------------
  subroutine expand_Gspace( ecut_in )
! ----------------------------------------------------------------------

  ! determines the largest G-vector index over all k-points and expands
  ! the G-space cut-off to include this G-vector at the Gamma-point

  use gvect,     only : ngm, gcutm, g, ig_l2g
  use cell_base, only : bg, tpiba2
  use klist,     only : nkstot, nks, xk, wk, ngk, igk_k
  use wvfct,     only : npw, g2kin, npwx
  use gvecw,     only : gcutw, ecutwfc
  use mp,        only : mp_max
  USE mp_global, ONLY : intra_pool_comm
  use io_global, only : stdout
  use fft_base,               only : dffts, dfftp

  integer,external :: n_plane_waves

  real(dp),intent(in) :: ecut_in

  integer :: ik, i, npw_store
  integer :: ikm, imap
  real(dp) :: gcut
  real(dp) :: xk_gamma(3,1)
  integer :: ngk_gamma(1)
  real(dp),allocatable :: xk_exp(:,:)
  integer :: ierr
  real(dp) :: g2, g2max
  real(dp) :: bglen


  ! summary of input details
  call summary

  ! old cut-off in Gspace units
  gcut = gcutw

  ! note that in case the lattice vectors have not been normalized
  bglen=dot_product(bg(:,1),bg(:,1))
 
  ! max G index
  call fix_shirley_planewaves
  write(stdout,*) 'max G-index over procs igwx = ', igwx

  ! ---------------------------------------------
  ! increase npwx to the max possible
  npwx = ngm
  ! ---------------------------------------------
  ! reallocate these arrays larger than necessary
  if( allocated( igk_gamma ) ) deallocate( igk_gamma )
  if( allocated( g2kin_gamma ) ) deallocate( g2kin_gamma )
  allocate( igk_gamma(npwx), g2kin_gamma(npwx) )

  if( allocated(igk_l2g_gamma) ) deallocate( igk_l2g_gamma )
  allocate( igk_l2g_gamma(npwx) )


  write(stdout,*)
  write(stdout,*) ' expand_Gspace -'
  write(stdout,*) '   increasing cut-off to include G-vectors from all k-points'
  write(stdout,*) '     at the Gamma-point'
  write(stdout,*)
  write(stdout,'(2x,a6,2a12,2a10,a10)') ' iter', ' max-G-index', ' max-k-index', ' new Ecut', ' old Ecut', ' npwx'

  ! Gamma-point vector
  xk_gamma = 0.d0

  i = 0
  igwx_gamma = 0
  do while( igwx_gamma < igwx .and. gcut <= gcutm )
    call gk_sort (xk_gamma(1,1), ngm, g, gcut, npw_gamma, igk_gamma, g2kin_gamma)
    call gk_l2gmap (ngm, ig_l2g, npw_gamma, igk_gamma, igk_l2g_gamma)

    ! remember where we started in size
    if( igwx_gamma==0 ) npw_store = npw_gamma

    ! for the Gamma-point the ordering of ig_l2g is correct
    ! or is it?
    igwx_gamma = maxval( igk_l2g_gamma(1:npw_gamma) )
    call mp_max( igwx_gamma, intra_pool_comm )

    i=i+1
    ! new cut-off
    gcutw_gamma = gcut
    ecutwfc_gamma = gcut * tpiba2

    ! report
    write(stdout,'(2x,i6,2i12,2f10.5,i10)') &
      i, igwx_gamma, igwx, ecutwfc_gamma, ecutwfc, npw_gamma

    gcut = gcut + 1.d0*bglen
  enddo
  gcut = gcut - 1.d0*bglen ! correct for additional increment
  write(stdout,*)
  write(stdout,*) ' all-encompassing cut-off      = ', ecutwfc_gamma, npw_gamma
  write(stdout,*)

  if( gcut > gcutm ) call errore('expand_Gspace','expanded cut-off too large',1)

  if( debug ) then
    write(stdout,*) ' dfftp%nnr(before) = ', dfftp%nnr
    write(stdout,*) ' dffts%nnr(before) = ', dffts%nnr
    write(stdout,*) '       ngm(before) = ', ngm
  endif

  do i=1,20
    write(stdout,'(2i6,3f7.3,"-",3f7.3)') i, igk_gamma(i), g(1:3,igk_gamma(i)), g(1:3,i)
  enddo

  ! copy new cut-off
  gcutw = gcutw_gamma
  ecutwfc = ecutwfc_gamma
  ! update FFT grids based on the new cut-off
  write(stdout,*) ' update FFT grids '
  call update_shirley_fft

  if( debug ) then
    write(stdout,*) ' dfftp%nnr(after) = ', dfftp%nnr
    write(stdout,*) ' dffts%nnr(after) = ', dffts%nnr
    write(stdout,*) '       ngm(after) = ', ngm
  endif

  write(stdout,*) 'local number of G-vectors before FFT update: ', npw_store

  npw_gamma = n_plane_waves (gcutw, 1, xk_gamma, g, ngm)

  write(stdout,*) 'local number of G-vectors after FFT update: ', npw_gamma

  do i=1,20
    write(stdout,'(2i6,3f7.3,"-",3f7.3)') i, igk_gamma(i), g(1:3,igk_gamma(i)), g(1:3,i)
  enddo

  ! allocate the correct amount of space and reinitialize indexing arrays
  if( allocated( igk_gamma ) ) deallocate( igk_gamma )
  if( allocated( g2kin_gamma ) ) deallocate( g2kin_gamma )
  if( allocated(igk_l2g_gamma) ) deallocate( igk_l2g_gamma )
  allocate( igk_gamma(npw_gamma), g2kin_gamma(npw_gamma), igk_l2g_gamma(npw_gamma) )
  ! ------------------------------
  ! set npwx - very important
  ! gk_sort uses npwx to size incoming arrays igk, g2kin
  npwx = npw_gamma
  ! ------------------------------
  call gk_sort (xk_gamma(1,1), ngm, g, gcutw_gamma, npw_gamma, igk_gamma, g2kin_gamma)
  call gk_l2gmap (ngm, ig_l2g, npw_gamma, igk_gamma, igk_l2g_gamma)

  write(stdout,*) ' all-encompassing cut-off = ', ecutwfc_gamma
  write(stdout,*)

  ! now restore commonly used variables: npwx, igk, g2kin, igk_l2g
  call fix_shirley_planewaves
  write(stdout,*) 'leaving expand_Gspace'
  return
  end subroutine expand_Gspace


! ----------------------------------------------------------------------
  subroutine update_shirley_fft( )
! ----------------------------------------------------------------------
  !
  ! update all fft information based on the new cut-off ecutwfc_gamma
  !
  use gvecw,                  only : gcutw, ecutwfc
  use gvect
  !use fft_types,              only : realspace_grids_init
  use gvecs,                  only : dual, doublegrid, gcutms, ngms
  use recvec_subs,            only : ggen, ggens
  use cell_base,              only : tpiba2, at, bg
  use control_flags,          only : gamma_only, smallmem
  use fft_base,               only : dffts, dfftp
  use io_global, only : stdout
  USE scf,       ONLY : rho


  !
  ! deallocate space used for FFT's and g-vectors
  !
  call deallocate_fft
  !
  ! ... Compute the cut-off of the G vectors
  !
  gcutm = dual * gcutw
  !
  doublegrid = ( dual > 4.D0 )
  !
  IF ( doublegrid ) THEN
     !
     gcutms = 4.D0 * gcutw
     !
  ELSE
     !
     gcutms = gcutm
     !
  END IF
  !
  ! ... calculate dimensions of the FFT grid
  !
  !call realspace_grids_init( at, bg, gcutm, gcutms )
  !
  ! ... determine the data structure for fft arrays
  !
  CALL data_structure( gamma_only )
  !
  ! ... print a summary and a memory estimate before starting allocating
  !
  CALL summary()
  CALL memory_report()
  !
  ! ... allocate memory for G- and R-space fft arrays
  !
  CALL allocate_fft()
  !
  ! ... generate reciprocal-lattice vectors and fft indices
  !
  IF( smallmem ) THEN
     CALL ggen( dfftp, gamma_only, at, bg, gcutm, ngm_g, ngm, &
          g, gg, mill, ig_l2g, gstart, no_global_sort = .TRUE. )
  ELSE
     CALL ggen( dfftp, gamma_only, at, bg, gcutm, ngm_g, ngm, &
       g, gg, mill, ig_l2g, gstart )
  END IF
  CALL ggens( dffts, gamma_only, at, g, gg, mill, gcutms, ngms )
  if (gamma_only) THEN
     ! ... Solvers need to know gstart
     call export_gstart_2_solvers(gstart)
  END IF
  !!
  !! ... allocate memory for G- and R-space fft arrays
  !!
  !call allocate_fft
  !rho%of_g = 0.d0
  !rho%of_r = 0.d0
  !!
  !! ... generate reciprocal-lattice vectors and fft indices
  !!
  !call ggen ( gamma_only, at, bg )
  !
  call summary
  !
  if( debug ) then
    write(stdout,*) ' update_shirley_fft grid dimensions:'
    write(stdout,*) ' dfftp%nnr = ', dfftp%nnr
    write(stdout,*) '      ngm  = ', ngm
  endif
  !
  return
  end subroutine update_shirley_fft


! ----------------------------------------------------------------------
  subroutine init_shirley_wfc_2( )
! ----------------------------------------------------------------------

  ! Initialize the memory and disk for shirley wave functions

  integer :: ierr

  allocate( evc_2(npw_gamma,nbnd_gamma), stat=ierr )
  if( ierr/=0 ) call errore('init_shirley_wfc','unable to allocate space for shirley wave function #2',1)

  return
  end subroutine init_shirley_wfc_2


! ----------------------------------------------------------------------
  subroutine write_shirley_wfc( ik )
! ----------------------------------------------------------------------

  integer,intent(in) :: ik

  call davcio( evc, nwordwfc_gamma, iunwfc_gamma, ik, +1 )

  return
  end subroutine write_shirley_wfc


! ---------------------------------------------------------------------- 
  subroutine read_shirley_wfc_1( ik )
! ---------------------------------------------------------------------- 
  use gvect,     only : g
  integer,intent(in) :: ik
  integer,save :: ik_current=0

  integer :: i, ibnd

  if( ik/=ik_current ) then
    call davcio( evc, nwordwfc_gamma, iunwfc_gamma, ik, -1 )
    ik_current = ik
  endif

  return
  end subroutine read_shirley_wfc_1


! ---------------------------------------------------------------------- 
  subroutine read_shirley_wfc_2( jk )
! ---------------------------------------------------------------------- 
  use gvect,     only : g
  integer,intent(in) :: jk
  integer,save :: jk_current=0

  integer :: i, ibnd

  if( jk/=jk_current ) then
    call davcio( evc_2, nwordwfc_gamma, iunwfc_gamma, jk, -1 )
    jk_current = jk
  endif

  return
  end subroutine read_shirley_wfc_2


! ---------------------------------------------------------------------- 
  subroutine delete_shirley_wfc
! ---------------------------------------------------------------------- 

  close(iunwfc_gamma,status='delete')

  return
  end subroutine delete_shirley_wfc


! ---------------------------------------------------------------------- 
  function dotprod_shirley_wfc( i, j )
! ---------------------------------------------------------------------- 
  USE mp,        ONLY : mp_sum
  USE mp_global, ONLY : intra_pool_comm
  complex(dp) :: dotprod_shirley_wfc
  integer,intent(in) :: i, j
  complex(dp) :: ZDOTC
  complex(dp) :: dpshwfc

  dpshwfc = ZDOTC( size(evc_2,1), evc_2(1,i), 1, evc(1,j), 1 )
  call mp_sum( dpshwfc, intra_pool_comm )
  dotprod_shirley_wfc = dpshwfc

  return
  end function dotprod_shirley_wfc


! ---------------------------------------------------------------------- 
  subroutine overlap_shirley_wfc( S )
! ---------------------------------------------------------------------- 
  complex(dp) :: S(:,:)
  complex(dp),parameter :: zero = cmplx(0.d0,0.d0)
  complex(dp),parameter :: one  = cmplx(1.d0,0.d0)

  call ZGEMM('C','N',size(evc_2,2),size(evc,2),size(evc_2,1),one, &
             evc_2,size(evc_2,1),evc,size(evc,1),zero,S,size(S,1)) 

  return
  end subroutine overlap_shirley_wfc


! ---------------------------------------------------------------------- 
  subroutine fix_shirley_planewaves
! ---------------------------------------------------------------------- 

  ! note that nkstot should be 1 now and ecutwfc should be larger

  ! be sure that number of planewaves is set up correctly in ngk(:)
  ! and that igk_l2g is set up for indexing wave function components

  use io_global, only : stdout
  use cell_base, only : tpiba2
  use klist,     only : xk, nks, nkstot, ngk, igk_k
  USE gvect,     ONLY : ngm, g, ig_l2g
  USE wvfct,     ONLY : npw, npwx, g2kin
  use gvecw,     only : gcutw, ecutwfc
  use mp,        only : mp_max
  USE mp_global, ONLY : intra_pool_comm

  integer,external :: n_plane_waves

  integer :: ik

  if( debug) then
    write(stdout,*)
    write(stdout,*) ' fix_shirley_planewaves:'
    write(stdout,*) ' establishing correct dimensions for planewaves'
    write(stdout,*) ' nkstot  = ', nkstot
    write(stdout,*) ' nks     = ', nks   
    write(stdout,*) ' ecutwfc = ', ecutwfc   
    write(stdout,*) ' gcutw   = ', gcutw, ecutwfc/tpiba2
    write(stdout,*) ' ngm = ', ngm   
    write(stdout,*) ' npw = ', npw   
    write(stdout,*)
  endif

  npwx = n_plane_waves (gcutw, nkstot, xk, g, ngm )
  if( allocated( igk_k ) ) deallocate( igk_k )
  allocate( igk_k(npwx,nkstot) )
  if( allocated( g2kin ) ) deallocate( g2kin )
  allocate( g2kin(npwx) )
  if( allocated( igk_l2g ) ) deallocate( igk_l2g )
  allocate( igk_l2g(npwx,nkstot) )
  igk_l2g=0
  do ik=1,nkstot
    call gk_sort (xk(1,ik), ngm, g, gcutw, npw, igk_k(1,ik), g2kin)
    call gk_l2gmap (ngm, ig_l2g(1), npw, igk_k(1,ik), igk_l2g(1,ik))
  enddo
  igwx=maxval(igk_l2g)
  call mp_max( igwx, intra_pool_comm )

  return
  end subroutine fix_shirley_planewaves


! ---------------------------------------------------------------------- 
  subroutine map_shirley_bz( ndim, band_subset, nkstot_mapped )
! ---------------------------------------------------------------------- 

  ! This is an important mapping of the available wave functions from
  ! their original k-points to k in [0,1]
  ! Note that we include the zone boundary at 1 aswell
  ! This means for example that the Gamma-point will be mapped
  ! to itself (000) and also (001), (010), (011), (100), (101), (110), (111)
  ! If we have a system of reduced dimensions (e.g. 1D-periodic) then the
  ! possible mappings is reduced, e.g. (000) and (001) for 1D along z

  ! I need to add inversion symmetry to this routine

  use constants, only : tpi
  use klist,     only : xk, nks, nkstot, ngk, wk, igk_k
  use cell_base, only : ibrav, at, bg, tpiba2, tpiba, alat, omega
  use symm_base, only : s, sname, ftau, nsym, t_rev, irt, ftau, invsym
  USE gvect, ONLY : ngm, ngm_g, gcutm, g, ig_l2g, gg
!  use grid_dimensions, only : nrxx
  use gvecs, only : dual, doublegrid, gcutms
!  use smooth_grid_dimensions, only : nrxxs
  USE wvfct,         ONLY : nbnd, npwx, npw, et, wg, nbndx, g2kin
  use gvecw,         only : gcutw, ecutwfc
  USE wavefunctions, ONLY : psic
  use noncollin_module, only : npol
  USE io_files,      ONLY : iunwfc, nwordwfc, prefix, diropn
  use io_global, only : stdout
  use control_flags, only: twfcollect
  USE mp_global, ONLY : me_pool, nproc_pool, intra_pool_comm, root_pool
  USE mp_wave,   ONLY : mergewf, splitwf
  USE mp,        ONLY : mp_max, mp_barrier, mp_sum
  use fft_base,  only : dffts, dfftp
  use fft_interfaces, only : fwfft, invfft

  use shirley_basis_input, only : ecut_in


  real(dp),parameter :: eps=1.d-12
  complex(dp),parameter :: iota=cmplx(0.d0,1.d0)

  integer,intent(in) :: ndim(3)
  integer,intent(in) :: band_subset(2)
  ! output the expanded and mapped set of k-points in [0,1]^3
  integer,intent(out) :: nkstot_mapped

  integer :: ik, isym, i, j
  integer :: imap
  integer :: ixyz

  integer :: nkstot_in
  integer :: nkstot_orig, npw_orig, npwx_orig
  real(dp),allocatable :: xk_orig(:,:)
  integer,allocatable :: igk_orig(:,:)
  integer,allocatable :: igk_l2g_orig(:,:)
  integer,allocatable :: ngk_orig(:)
  real(dp),allocatable :: g2kin_orig(:,:)

  integer :: igwx

  real(dp),allocatable :: r(:,:)
  complex(dp),allocatable :: expikr(:)
  integer :: ibnd, jbnd, ikmapped
  complex(dp),allocatable :: evc_tmp(:,:)
  integer :: nbnd_subset

  logical :: exst, opend
  integer,external :: freeunit

  complex(dp) :: norm
  complex(dp),external :: ZDOTC
  complex(dp) :: dotprod, dotprodw
  real(dp) :: dg(3)

  complex(dp),allocatable :: wtmp(:), evc1(:)
  integer :: ierr

  write(stdout,*)
  write(stdout,*) ' map_shirley_bz '
  write(stdout,*) ' Map wave functions on input k-point set to all equivalent'
  write(stdout,*) ' points in the positive BZ [0,1]^3'
  write(stdout,*)

  if( debug ) then
    ! check incoming vectors
    ibnd=1; ik=1
    write(stdout,*) ' npw = ', npw
    write(stdout,*) ' ngk(ik) = ', ngk(ik), ik
    write(stdout,*) '    wave function coefficients of band ', ibnd
    do i=1,min(50,ngk(ik))
      write(stdout,'(i6,3f7.3,2f16.8,i6)') i, g(1:3,igk_k(i,ik)), evc(i,ibnd), igk_k(i,ik)
    enddo
    write(stdout,*)
  endif

  ! copy original k-point set information
  nkstot_orig = nkstot
  npwx_orig = npwx
  allocate( xk_orig(3,nkstot_orig) )
  xk_orig = xk(1:3,1:nkstot_orig)
  allocate( igk_orig(npwx_orig,nkstot_orig) )
  allocate( g2kin_orig(npwx_orig,nkstot_orig) )
  allocate( ngk_orig(nkstot_orig) )
  ngk_orig = ngk(1:nkstot_orig)
  allocate( igk_l2g_orig(npwx_orig,nkstot_orig) )
  ! I no longer have access to igk_l2g - need to regenerate
  if( allocated(igk_l2g) ) deallocate( igk_l2g )
  allocate( igk_l2g( npwx, nks ) )
  igk_l2g = 0
  do ik=1,nks
    call gk_sort (xk_orig(1,ik), ngm, g, gcutw, npw, igk_orig(1,ik), g2kin_orig)
    call gk_l2gmap( ngm, ig_l2g, npw, igk_orig(1,ik), igk_l2g(1,ik) )
  enddo
  igk_l2g_orig = igk_l2g


  ! expand kpoints using symmetry, inversion and translation
  ! note that currently only translation is included
  call kpoint_expand( ndim )


  ! update k-points based on expansion from kpoint_expand
  nkstot = sum(xk_map(:)%nmap)
  nks = nkstot
  write(stdout,*) ' new number of mapped k-points = ', nkstot
  ikmapped=0
  do ik=1,nkstot_orig
    do imap=1,xk_map(ik)%nmap
      ikmapped=ikmapped+1
      xk(:,ikmapped) = matmul( bg, xk_map(ik)%xmap(:,imap) ) + xk_orig(:,ik)
      wk(ikmapped) = 0.d0 ! I'm zeroing this just in case
    enddo
  enddo
  if( allocated( ngk ) ) deallocate( ngk )
  allocate( ngk(nkstot) )


  ! expand G-space size to include G's from all k-points
  ! also updates FFT grids
  call expand_Gspace( ecut_in )

  ! store new larger k-point set size
  nkstot_mapped = nkstot


  ! new output file for mapped and reordered  wave functions
  iunwfc_gamma = freeunit()
  nbnd_subset = band_subset(2)-band_subset(1)+1
  nwordwfc_gamma = 2 * nbnd_subset * npw_gamma * npol
  call diropn( iunwfc_gamma, 'wfcr', nwordwfc_gamma, exst )

  ! reallocate evc for output of final wave functions
  if( allocated(evc) ) deallocate(evc)
  allocate( evc(npw_gamma,nbnd_subset) )

  ! make real-space grid for r, and expikr
  write(stdout,*) ' rgrid '
  allocate( r(dffts%nnr,3), expikr(dffts%nnr) )
  call rgrid( r )

  ! temporary arrays for evc
  allocate( evc_tmp(1:npwx_orig,1:nbnd) )
 
  ! begin mapping wave functions
  ikmapped = 0
  do ik=1,nkstot_orig
    
    write(stdout,'(a,i6,a,3f8.3)') ' for original k-point ', ik, ' = ', xk_orig(1:3,ik)
!    if( debug ) then
      write(stdout,*) '    loading wave functions from file unit ', iunwfc
!   endif

    ! load wave functions for ik
    call davcio( evc_tmp, 2*nwordwfc, iunwfc, ik, -1 )

    write(stdout,*) 'mapping k-point ', ik

    ! make phases
    do imap=1,xk_map(ik)%nmap
      
      ikmapped=ikmapped+1

      ! determine igk and igk_l2g
      write(stdout,*) '   mapped k-point ', ikmapped
      CALL gk_sort (xk(1,ikmapped), ngm, g, gcutw, npw, igk_k(1,ikmapped), g2kin)
      call gk_l2gmap (ngm, ig_l2g(1), npw, igk_k(1,ikmapped), igk_l2g(1,ikmapped))
      ngk(ikmapped)=npw

      ! allocate a 1D eigenvector
      allocate( evc1( max(npw,ngk_orig(ik),npw_gamma) ) )

      igwx=max( maxval(igk_l2g(1:npw,ikmapped)),          &
                maxval(igk_l2g_orig(1:ngk_orig(ik),ik)) )
      call mp_max( igwx, intra_pool_comm )
      if( allocated(wtmp) ) deallocate(wtmp)
      allocate( wtmp(igwx), stat=ierr )
      if( ierr /= 0 ) &
        call errore('map_shirley_bz','unable to allocate temporary space',igwx)

      ! e^(-i 2 pi xk_map.r ) - r and xk in crystal units
      expikr(:) = matmul( r, xk_map(ik)%xmap(:,imap) )
      expikr(:)= exp( ( - iota * tpi ) * expikr(:) )



      ! loop over bands to map and re-order wave function coefficients
      do ibnd=1,nbnd_subset
        !
        jbnd = ibnd-1+band_subset(1)
        !
        ! tmp copy of evc_tmp
        evc1(1:ngk_orig(ik))=evc_tmp(1:ngk_orig(ik),jbnd)

        ! re-order based on new g-vector ordering with larger inclusive cut-off
        ! zero the temp array
        !if( debug ) write(stdout,*) ' reorder original k-point ', ik, ' band ', jbnd
        write(stdout,*) ' reorder original k-point ', ik, ' band ', jbnd
        wtmp = zero

        ! merge wave function from evc into wtmp
        call mergewf(evc1(:), wtmp, ngk_orig(ik), igk_l2g_orig(:,ik), me_pool, nproc_pool, root_pool, intra_pool_comm)

        if( debug ) then
          ! check
          dotprod = dot_product( evc1(1:ngk_orig(ik)), evc1(1:ngk_orig(ik)) )
#ifdef __MPI
          call mp_sum( dotprod, intra_pool_comm )
#endif
          dotprodw = dot_product(wtmp, wtmp)
          if( abs(real(dotprodw-dotprod)) > 1.d-12 ) then
            write(stdout,'(a,i6,3e12.5)') 'warning: bad merged norm for band ', &
              jbnd, real(dotprod), real(dotprodw), real(dotprodw-dotprod) 
          endif
        endif

        ! zero the former wave function
        evc1(:) = zero

        ! split wave function from wtmp into evc
        call splitwf(evc1(:), wtmp, npw, igk_l2g(:,ikmapped), me_pool, nproc_pool, root_pool, intra_pool_comm)

        if( debug ) then
        if( ibnd==1 ) then
          write(stdout,*) '    wave function coefficients of band ', jbnd
          do i=1,min(50,npw)
            write(stdout,'(i6,3f7.3,2f16.8,i6)') i, g(1:3,igk_k(i,ikmapped)), evc1(i), igk_k(i,ikmapped)
          enddo
          write(stdout,*)
        endif
        endif

        ! now FFT
        psic(:) = ( 0.0D0, 0.0D0 )
        psic(dffts%nl(igk_k(1:npw,ikmapped))) = evc1(1:npw)

!        if( debug ) then
!          ! check norm
!          norm = ZDOTC( ngk_orig(ik), evc_tmp(1,jbnd), 1, evc_tmp(1,jbnd), 1 )
!          call mp_sum(norm, intra_pool_comm)
!          write(stdout,*) ' mapping ', ikmapped, jbnd, norm
!        endif

        CALL start_clock( 'firstfft' )
        !
        CALL invfft ('Wave', psic, dffts)
        !
        CALL stop_clock( 'firstfft' )
        !
        if( debug ) then
          ! check norm
          norm = ZDOTC( size(psic), psic(1), 1, psic(1), 1 )
          call mp_sum(norm, intra_pool_comm)
          norm = norm / dble(dffts%nr1*dffts%nr2*dffts%nr3)
          write(stdout,*) ' ifft ', ikmapped, jbnd, norm
        endif

        !
        ! ... product with the phase factor expikr =  exp(i xmap.r) 
        !
        psic(1:dffts%nnr) = psic(1:dffts%nnr) * expikr(1:dffts%nnr)

        if( debug ) then
          ! check norm
          norm = ZDOTC( size(psic), psic(1), 1, psic(1), 1 )
          call mp_sum(norm, intra_pool_comm)
          norm = norm / dble(dffts%nr1*dffts%nr2*dffts%nr3)
          write(stdout,*) ' expikr ', ikmapped, jbnd, norm
        endif
        !
        ! ... back to reciprocal space
        !
        CALL start_clock( 'secondfft' )
        !
        call fwfft ('Wave', psic, dffts)
        !
        CALL stop_clock( 'secondfft' )
        !
        ! ... store with correct ordering
        !
        evc1(1:npw) = psic(dffts%nl(igk_k(1:npw,ikmapped)))

        if( debug ) then
          ! check norm
          norm = ZDOTC( size(evc1,1), evc1(1), 1, evc1(1), 1 )
          call mp_sum(norm, intra_pool_comm)
          write(stdout,*) ' fft ', ikmapped, jbnd, norm
        endif

        !

        if( debug ) then
          if( ibnd==1 ) then
            write(stdout,'(a,i6,a,3f8.3)') ' for mapped k-point ', ikmapped, ' = ', xk(1:3,ikmapped)

            write(stdout,*) '    wave function coefficients of band ', ibnd
            do i=1,min(50,ngk(ikmapped))
              write(stdout,'(i6,3f7.3,2f16.8)') i, g(1:3,igk_k(i,ikmapped)), evc1(i)
            enddo
            write(stdout,*) '    ...'
            write(stdout,*) ' size igk = ', size(igk_k,1), ngk(ikmapped)
            do i=max(1,ngk(ikmapped)-10),ngk(ikmapped)
              write(stdout,'(i6,3f7.3,2f16.8)') i, g(1:3,igk_k(i,ikmapped)), evc1(i)
            enddo
            write(stdout,*)
          endif
        endif
       
        !
        ! re-order as for Gamma-point with larger inclusive cut-off
        ! zero the temp array
        !if( debug ) write(stdout,*) ' reorder mapped k-point ', ikmapped, ' band ', ibnd
        write(stdout,*) ' reorder mapped k-point ', ikmapped, ' band ', ibnd
        wtmp = zero

        ! merge wave function from evc into wtmp
        call mergewf(evc1, wtmp, ngk(ikmapped), igk_l2g(:,ikmapped), me_pool, nproc_pool, root_pool, intra_pool_comm)

        if( debug ) then
          ! check
          dotprod = dot_product( evc1(1:ngk(ikmapped)), evc1(1:ngk(ikmapped)) )
#ifdef __MPI
          call mp_sum( dotprod, intra_pool_comm )
#endif
          dotprodw = dot_product(wtmp, wtmp)
          if( abs(real(dotprodw-dotprod)) > 1.d-12 ) then
            write(stdout,'(a,i6,3e12.5)') 'warning: bad merged norm for band ', &
              jbnd, real(dotprod), real(dotprodw), real(dotprodw-dotprod) 
          endif
        endif

        ! zero the output wave function
        evc(:,ibnd) = zero

        ! split wave function from wtmp into evc
        call splitwf(evc(:,ibnd), wtmp, npw_gamma, igk_l2g_gamma(:), me_pool, nproc_pool, root_pool, intra_pool_comm)

        if( debug ) then
          ! check
          dotprod = dot_product( evc(1:npw_gamma,ibnd), evc(1:npw_gamma,ibnd) )
#ifdef __MPI
          call mp_sum( dotprod, intra_pool_comm )
#endif
          if( abs(real(dotprodw-dotprod)) > 1.d-12 ) then
            write(stdout,'(a,i6,3e12.5)') 'warning: bad split norm for band ', &
              jbnd, real(dotprod), real(dotprodw), real(dotprodw-dotprod) 
          endif
        endif

        if( debug ) then
          if( ibnd==1 ) then
            write(stdout,'(a,i6,a,3f16.8)') ' for reordered k-point ', ikmapped, ' = ', xk(1:3,ikmapped)

            write(stdout,*) '    wave function coefficients of band ', ibnd
            do i=1,min(50,npw_gamma)
              write(stdout,'(i6,3f7.3,2f16.8)') i, g(1:3,igk_gamma(i)), evc(i,ibnd)
            enddo
            write(stdout,*) '    ...'
            do i=max(1,npw_gamma-10),npw_gamma
              write(stdout,'(i6,3f7.3,2f16.8)') i, g(1:3,igk_gamma(i)), evc(i,ibnd)
            enddo
            write(stdout,*)
          endif
        endif
       

      enddo

      ! write to file
      write(stdout,*) '     dumping mapped k-point ', ikmapped
      call davcio( evc, nwordwfc_gamma, iunwfc_gamma, ikmapped, +1 )

      if( debug ) then
        ! check norm of mapped wave functions - evc
        do ibnd=1,nbnd_subset
          norm = ZDOTC( size(evc,1), evc(1,ibnd), 1, evc(1,ibnd), 1 )
          call mp_sum(norm, intra_pool_comm)
          write(stdout,'(a,i6,a,i6,a,2f16.8)') 'map_shirley_bz (mapped): norm for k ', ikmapped, ' n ', ibnd, ' = ', norm
        enddo
      endif

      ! free memory
      deallocate( evc1 )
      deallocate( wtmp )

    enddo
    
  enddo

  if( allocated(evc_tmp) ) deallocate( evc_tmp )
  if( allocated(r) ) deallocate( r )
  if( allocated(expikr) ) deallocate( expikr )
  if( allocated(xk_orig) ) deallocate( xk_orig )
  if( allocated(igk_orig) ) deallocate( igk_orig )


  ! modify the wave function variables for future reference
  inquire(unit=iunwfc,opened=opend)
  if( opend ) then
    if( twfcollect ) then
      close( iunwfc, status='delete' )
    else
      close( iunwfc )
    endif
  endif
  ! switch wave function unit to wfcr's
  iunwfc = iunwfc_gamma
  ! new convention for nwordwfc in v4.3
  nwordwfc = nwordwfc_gamma/2
  nbnd = nbnd_subset
  nbndx = nbnd


  ! from this point on we assume only one k-point k=0
  ! we set nkstot and nks to 1
  nkstot = 1
  nks = 1
  xk(1:3,1) = 0.d0  ! N.B. Important to set this to Gamma point
                    !      It is stored in the save file
  wk(1) = 2.d0      ! the weight is two for non-spin-polarized
  ! spin should have no meaning for the optimal basis - it should work for both
  ! spin-channels, so really wk is just a dummy variable now
  gcutw = gcutw_gamma
  ecutwfc = ecutwfc_gamma
  ! copy number of bands
  nbnd_gamma = nbnd_subset

  write(stdout,*) ' done with map_shirley_bz'

  end subroutine map_shirley_bz


  end module wfc_shirley
