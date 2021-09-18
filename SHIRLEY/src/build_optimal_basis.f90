  subroutine build_optimal_basis()

  ! build optimal basis set by expanding existing periodic functions for each k-point to fill the unit cube
  ! and then diagonalize their overlap matrix

  ! this routine reads existing wavefunctions from file and rewrites the
  ! orthogonalized versions in the ordering of the Gamma-point (NB)

  USE kinds,         ONLY : dp
  USE io_global,     ONLY : stdout
  USE wvfct,         ONLY : nbnd, npwx, et, wg, nbndx
  USE wavefunctions, ONLY : evc
  use noncollin_module, only : npol
  USE io_files,      ONLY : iunwfc, nwordwfc, prefix, tmp_dir, diropn
  use input_parameters, only : pseudo_dir
  use mp_global, only : kunit
  USE lsda_mod,             ONLY : nspin, lsda
  USE noncollin_module,     ONLY : noncolin
  use control_flags, only: twfcollect

  USE basis_shirley
  USE wfc_shirley, only : expand_Gspace, &
                          init_shirley_wfc_2, &
                          write_shirley_wfc, &
                          read_shirley_wfc_1, read_shirley_wfc_2, &
                          dotprod_shirley_wfc, npw_gamma, ecutwfc_gamma, &
                          map_shirley_bz, &
                          delete_shirley_wfc, fix_shirley_planewaves
  use pot_rho_G
  use shirley_basis_input, only : trace_tol, ndim, band_subset, &
                                  expand_kpts, ecut_in, debug

  implicit none

  complex(dp),parameter :: zero=(0.d0,0.d0)
  complex(dp),parameter :: one =(1.d0,0.d0)

  integer :: nkstot_mapped
  integer :: ik, ibnd
  integer :: jk, jbnd
  integer :: i, j, k

  real(dp) :: invnorm

  integer :: ierr, kunittmp, iunpun
  logical :: exst, opend

  complex(dp),allocatable :: ztmp(:)

  integer :: ios

  complex(dp),external :: ZDOTC
  complex(dp) :: norm

  write(stdout,*)
  write(stdout,*) ' total number of bands = ', nbnd
  if( band_subset(1) == 0 ) band_subset(1) = 1
  if( band_subset(2) == 0 ) band_subset(2) = nbnd
  if( band_subset(1) < 1 .or. band_subset(1) > nbnd .or. &
      band_subset(1) < 1 .or. band_subset(1) > nbnd ) then
    call errore('build_optimal_basis','band subset is screwy',-1)
    band_subset(1) = 1
    band_subset(2) = nbnd
  endif
  if( band_subset(1) > band_subset(2) ) then
    i = band_subset(1)
    band_subset(1) = band_subset(2)
    band_subset(2) = i
  endif
  write(stdout,*) ' band subset = ', band_subset(1), ' ... ', band_subset(2)
  write(stdout,*) ' number of bands in subset = ', band_subset(2)-band_subset(1)+1

  ! read the potential and density from file and transform to G-space
  call read_pot_rho_G

  ! map wave functions to [0,1] 
  ! (Note we are including 1, so the Gamma point will get mapped irregardless)
  ! this routine also redefines nbnd based on band_subset and some terms
  ! related to wfc storage
  !if( expand_kpts ) call map_shirley_bz( ndim, band_subset )
  ! need to rewrite for the case where we don't expand
  call map_shirley_bz( ndim, band_subset, nkstot_mapped )
  ! from now on nks=1 and nkstot=1

  ! allocate another shirley wavefunction array
  call init_shirley_wfc_2

  ! initialize matrices for Shirley basis
  call basis_init( nkstot_mapped, nbnd )

  ! construct overlap
  call construct_overlap( nkstot_mapped, nbnd )

  ! diagonalize overlap
  call diagonalize_overlap

  ! truncate the basis
  call truncate_basis( trace_tol )

  ! expand new basis
  call expand_basis( npw_gamma )

  ! don't need reordered wave functions anymore
  call delete_shirley_wfc

  ! change dimensions etc of klist and wvfct
  ! Do I need this?
  call fix_shirley_planewaves

  ! new band dimensions based on truncation
  nbnd = nbasis_trunc
  nbndx = nbnd

  ! et and wg
  if( allocated( et ) ) deallocate( et )
  if( allocated( wg ) ) deallocate( wg )
  allocate( et(nbnd,1), wg(nbnd,1) )
  et = 0.d0
  wg = 1.d0


  ! note that we should copy the old pseudopotentials from the .save directory to the new _opt directory
  pseudo_dir = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
  write(stdout,*) ' changing pseudo_dir to ', trim(pseudo_dir)
  ! prepare for save to new _opt directory
  prefix = trim(prefix) // '_opt'
  write(stdout,*) ' changing prefix to ', trim(prefix)

  if( debug ) then
    write(stdout,*) ' dfftp info'
    write(stdout,*) dfftp%nr1, dfftp%nr2, dfftp%nr3, &
                    dfftp%nr1x, dfftp%nr2x, dfftp%nnp
  endif

  ! transform the potential and density from G-space and write to file
  write(stdout,*) ' dump charge and potential to file using updated FFT grid'
  ! be careful not to mess with nspin before this point, since it
  ! determines how the charge is written
  ! nspin=2 also writes spin-polarization
  call write_pot_rho_G

  nwordwfc = nbndx * npwx * npol
  write(stdout,*) 'nbndx = ', nbndx, ' npwx = ', npwx, ' npol = ', npol
  inquire(unit=iunwfc,opened=opend)
  if( opend ) then
    if( twfcollect ) then
      close( iunwfc, status='delete', iostat=ios )
    else
      close( iunwfc, iostat=ios )
    endif
    if( ios/=0 ) call errore('build_optimal_basis','unable to close unit',-1*iunwfc)
  endif
  call diropn( iunwfc, 'wfc', 2*nwordwfc, exst )
  write(stdout,*) ' nwordwfc = ', nwordwfc
  write(stdout,*) ' size(evc) = ', size(evc,1), size(evc,2), size(evc)
  call davcio( evc, 2*nwordwfc, iunwfc, 1, +1 )
    
  write(stdout,*) ' about to save'

  kunittmp = kunit

  ! quick switch off of spin before writing
  ! the optimal basis should be independent of spin
  if( lsda .or. (nspin > 1) ) then
    write(stdout,*) ' Note: spin is not used for the optimal basis'
    lsda=.false.
    nspin=1
  endif
  ! I don't yet have a noncolin implementation
  noncolin=.false.

  call punch( 'all' )

  write(stdout,*) ' done'
  return
  end subroutine build_optimal_basis

