!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM do_projwfc
  !-----------------------------------------------------------------------
  !
  ! This program is based on projwfc.f90
  ! The key difference is that this program aims to find the entire
  ! projection of a given angular momentum channel on a specific atomic site.
  ! Unlike projwfc, this code does not rely on radial atomic functions to
  ! define the projection, since these typically only reflect a subspace of
  ! the entire l projection.
  !
  ! For example, if we have a F atom whose pseudopotential is well represented
  ! using only 2s and 2p projectors, then we may neglect any higher principal
  ! quantum number components (3s, 3p, etc.) of a given l projection, say l=1.
  ! Typically, this is not a problem for neutral ground state atoms, where 
  ! the higher quantum number orbitals are very high in energy, or unbound.
  ! However, for positively charged states or core-excited atoms, there will be
  ! an entire Rydberg series of bound states that site below the vacuum 
  ! potential and could easily contribute to electronic structure not far from
  ! the band edges.
  !
  ! My edit is to remove the use of atomwfc in the projection and to supply 
  ! an atomic wave function for the ultimate matrix element.
  ! To start with, the 1s dipole matrix element r phi_1s(r) can be rewritten
  ! as a product of the 1s radial function and the spherical harmonics from
  ! l=1 (selection rule). This will be the projection we care about for one
  ! specific atom.
  !
  ! old notes:
  !
  ! projects wavefunctions onto orthogonalized atomic wavefunctions,
  ! projected DOS
  !
  ! See files INPUT_PROJWFC.* in Doc/ directory for usage
  ! IMPORTANT: since v.5 namelist name is &projwfc and no longer &inputpp
  !
  USE parameters, ONLY : npk
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE constants,  ONLY : rytoev
  USE kinds,      ONLY : DP
  USE klist,      ONLY : nks, nkstot, xk, degauss, ngauss, lgauss, ltetra
  USE io_files,   ONLY : nd_nmbr, prefix, tmp_dir
  USE noncollin_module, ONLY : noncolin
  USE mp,         ONLY : mp_bcast
  USE spin_orb,   ONLY: lforcet
  USE mp_world,   ONLY : world_comm
  USE mp_global,  ONLY : mp_startup, nproc_ortho, nproc_pool, nproc_pool_file
  USE environment,ONLY : environment_start, environment_end
  USE wvfct,      ONLY : et, nbnd
  USE basis,      ONLY : natomwfc
  USE control_flags, ONLY: twfcollect
  USE paw_variables, ONLY : okpaw
  ! following modules needed for generation of tetrahedra
  USE ktetra,     ONLY : tetra, tetra_type, opt_tetra_init
  USE symm_base,  ONLY : nsym, s, time_reversal, t_rev
  USE cell_base,  ONLY : at, bg
  USE start_k,    ONLY : k1, k2, k3, nk1, nk2, nk3
  USE lsda_mod,   ONLY : lsda
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CHARACTER (len=256) :: filpdos, filproj, outdir
  REAL (DP), allocatable :: xk_collect(:,:)
  REAL (DP) :: Emin, Emax, DeltaE, degauss1, ef_0
  INTEGER :: nks2, ngauss1, ios
  LOGICAL :: lwrite_overlaps, lbinary_data
  LOGICAL :: lsym, kresolveddos, tdosinboxes, plotboxes, pawproj
  INTEGER, PARAMETER :: N_MAX_BOXES = 999
  INTEGER :: n_proj_boxes, irmin(3,N_MAX_BOXES), irmax(3,N_MAX_BOXES)
  LOGICAL :: lgww  !if .true. use GW QP energies from file bands.dat
  integer :: i
  !
  integer :: natom_center
  integer,allocatable :: atom_center_index(:)
  integer :: nr_center, l_center
  real(dp), allocatable :: r_center(:), u_center(:)
  !
  NAMELIST / projwfc / outdir, prefix, ngauss, degauss, lsym, &
             Emin, Emax, DeltaE, filpdos, filproj, lgww, &
             kresolveddos, tdosinboxes, n_proj_boxes, irmin, irmax, plotboxes, &
             lwrite_overlaps, lbinary_data, pawproj, lforcet, ef_0
  !
  ! initialise environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'PROJWFC' )
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  filproj= ' '
  filpdos= ' '
  Emin   =-1000000.d0
  Emax   =+1000000.d0
  DeltaE = 0.01d0
  ngauss = 0
  lsym   = .true.
  degauss= 0.d0
  lgww   = .false.
  pawproj= .false.
  lwrite_overlaps   = .false.
  lbinary_data = .false.
  kresolveddos = .false.
  tdosinboxes = .false.
  plotboxes   = .false.
  n_proj_boxes= 1
  irmin(:,:)  = 1
  irmax(:,:)  = 0
  !
  ios = 0
  !

  ef_0 = 0.d0
  lforcet = .false.


  IF ( ionode )  THEN
     !
     CALL input_from_file ( )
     !
     READ (5, projwfc, iostat = ios)
     !
     tmp_dir = trimcheck (outdir)
     ! save the value of degauss and ngauss: they are read from file
     degauss1=degauss
     ngauss1 = ngauss
     !
     ! now read the atom center info
     !
     write(stdout,*) ' Reading details of atomic centers '
     !
     read(5,*) natom_center
     allocate( atom_center_index(natom_center) )
     read(5,*) atom_center_index
     read(5,*)
     read(5,*) nr_center, l_center
     allocate( r_center(nr_center), u_center(nr_center) )
     do i=1,nr_center
       read(5,*) r_center(i), u_center(i)
     enddo
     !
     write(stdout,*) ' natom_center = ', natom_center
     write(stdout,*) ' atom_center_index = ', atom_center_index
     write(stdout,*) ' nr_center, l_center = ', nr_center, l_center
     !
  ENDIF
  !
  CALL mp_bcast (ios, ionode_id, world_comm )
  IF (ios /= 0) CALL errore ('do_projwfc', 'reading projwfc namelist', abs (ios) )
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir,   ionode_id, world_comm )
  CALL mp_bcast( prefix,    ionode_id, world_comm )
  CALL mp_bcast( filproj,   ionode_id, world_comm )
  CALL mp_bcast( ngauss1,   ionode_id, world_comm )
  CALL mp_bcast( degauss1,  ionode_id, world_comm )
  CALL mp_bcast( DeltaE,    ionode_id, world_comm )
  CALL mp_bcast( lsym,      ionode_id, world_comm )
  CALL mp_bcast( Emin,      ionode_id, world_comm )
  CALL mp_bcast( Emax,      ionode_id, world_comm )
  CALL mp_bcast( lwrite_overlaps, ionode_id, world_comm )
  CALL mp_bcast( lbinary_data,    ionode_id, world_comm )
  CALL mp_bcast( lgww,      ionode_id, world_comm )
  CALL mp_bcast( pawproj,   ionode_id, world_comm )
  CALL mp_bcast( tdosinboxes,     ionode_id, world_comm )
  CALL mp_bcast( n_proj_boxes,    ionode_id, world_comm )
  CALL mp_bcast( irmin,     ionode_id, world_comm )
  CALL mp_bcast( irmax,     ionode_id, world_comm )
  CALL mp_bcast( ef_0, ionode_id, world_comm )
  CALL mp_bcast( lforcet, ionode_id, world_comm )

  ! broadcast atom center info
  CALL mp_bcast( natom_center,   ionode_id, world_comm )
  if( .not. allocated(atom_center_index) ) allocate(atom_center_index(natom_center))
  CALL mp_bcast( atom_center_index,   ionode_id, world_comm )
  CALL mp_bcast( nr_center,   ionode_id, world_comm )
  CALL mp_bcast( l_center,   ionode_id, world_comm )
  if( .not. allocated(r_center) ) allocate(r_center(nr_center))
  if( .not. allocated(u_center) ) allocate(u_center(nr_center))
  CALL mp_bcast( r_center,   ionode_id, world_comm )
  CALL mp_bcast( u_center,   ionode_id, world_comm )

  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file ( )
  !
  IF (pawproj) THEN
    IF ( .NOT. okpaw ) CALL errore ('projwfc','option pawproj only for PAW',1)
    IF ( noncolin )  CALL errore ('projwfc','option pawproj and noncolinear spin not implemented',2)
  END IF
  !
  IF (nproc_pool /= nproc_pool_file .and. .not. twfcollect)  &
     CALL errore('projwfc',&
     'pw.x run with a different number of procs/pools. Use wf_collect=.true.',1)
  !
  CALL openfil_pp ( )
  !
  !   Tetrahedron method
  !
  IF ( ltetra ) THEN
     !
     ! info on tetrahedra is no longer saved to file and must be rebuilt
     ! workaround for old xml file, to be removed when the old xml file is
     IF(ALLOCATED(tetra)) DEALLOCATE(tetra)
     !
     ! in the lsda case, only the first half of the k points
     ! are needed in the input of "tetrahedra"
     !
     IF ( lsda ) THEN
        nks2 = nkstot / 2
     ELSE
        nks2 = nkstot
     END IF
     IF(tetra_type < 2) THEN
        ! use linear tetrahedron for both tetra_type=0 and tetra_type=1
        IF ( tetra_type == 0 ) tetra_type = 1
        WRITE( stdout,'(/5x,"Linear tetrahedron method (read from file) ")')
     ELSE
        WRITE( stdout,'(/5x,"Optimized tetrahedron method (read from file) ")')
     END IF
     !
     ! not sure this is needed
     !
     ALLOCATE(xk_collect(3,nkstot))
     CALL poolcollect(3, nks, xk, nkstot, xk_collect)
     !
     CALL opt_tetra_init(nsym, s, time_reversal, t_rev, at, bg, npk, k1,k2,k3, &
          &              nk1, nk2, nk3, nks2, xk_collect, 1)
     !
     DEALLOCATE(xk_collect)
     !
  ELSE IF (degauss1/=0.d0) THEN
     degauss=degauss1
     ngauss =ngauss1
     WRITE( stdout,'(/5x,"Gaussian broadening (read from input): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
     lgauss=.true.
  ELSE IF (lgauss) THEN
     WRITE( stdout,'(/5x,"Gaussian broadening (read from file): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
  ELSE
     degauss=DeltaE/rytoev
     ngauss =0
     WRITE( stdout,'(/5x,"Gaussian broadening (default values): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
     lgauss=.true.
  ENDIF
  !
  IF ( filpdos == ' ') filpdos = prefix
  !


IF ( lforcet ) THEN
!    CALL projwave_nc(filproj,lsym,lwrite_overlaps,lbinary_data,ef_0)
ELSE
  IF ( tdosinboxes ) THEN
!     CALL projwave_boxes (filpdos, filproj, n_proj_boxes, irmin, irmax, plotboxes)
  ELSE IF ( pawproj ) THEN
!     CALL projwave_paw (filproj)
  ELSE
     IF ( natomwfc <= 0 ) CALL errore &
        ('do_projwfc', 'Cannot project on zero atomic wavefunctions!', 1)
     IF (noncolin) THEN
!        CALL projwave_nc(filproj, lsym, lwrite_overlaps, lbinary_data,ef_0)
     ELSE
        IF( nproc_ortho > 1 ) THEN
!           CALL pprojwave (filproj, lsym, lwrite_overlaps, lbinary_data )
        ELSE
           CALL lprojwave (natom_center, atom_center_index, &
                           l_center, nr_center, r_center, u_center, &
                           filproj, lwrite_overlaps, lbinary_data)
        ENDIF
     ENDIF
  ENDIF
  !
  IF ( ionode ) THEN
     IF ( tdosinboxes ) THEN
!        CALL partialdos_boxes (Emin, Emax, DeltaE, kresolveddos, filpdos, n_proj_boxes)
     ELSE IF ( lsym .OR. kresolveddos ) THEN
        IF (noncolin) THEN
!           CALL partialdos_nc (Emin, Emax, DeltaE, kresolveddos, filpdos)
        ELSE
           CALL partialdos (Emin, Emax, DeltaE, kresolveddos, filpdos)
        ENDIF
     ENDIF
  ENDIF
ENDIF


  !
  CALL environment_end ( 'PROJWFC' )
  !
  CALL stop_pp
  !
END PROGRAM do_projwfc

!-----------------------------------------------------------------------
SUBROUTINE lprojwave( natom, atom_index, l, nr, r, u, filproj, lwrite_ovp, lbinary )
  !-----------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout, ionode
  USE run_info, ONLY: title
  USE ions_base, ONLY : zv, tau, nat, ntyp => nsp, ityp, atm
  USE basis,     ONLY : natomwfc
  USE cell_base
  USE constants, ONLY: rytoev, eps4
  USE gvect
  USE gvecs,   ONLY: dual
  USE gvecw,   ONLY: ecutwfc
  USE fft_base, ONLY : dfftp
  USE klist, ONLY: xk, nks, nkstot, nelec, ngk, igk_k
  USE lsda_mod, ONLY: nspin, isk, current_spin
  USE symm_base, ONLY: nsym, irt, d1, d2, d3
  USE wvfct, ONLY: npwx, nbnd, et, wg
  USE control_flags, ONLY: gamma_only
  USE uspp, ONLY: nkb, vkb
  USE us,         ONLY : nqx
  USE becmod,   ONLY: bec_type, becp, calbec, allocate_bec_type, deallocate_bec_type
  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc
  USE wavefunctions, ONLY: evc
  !
  USE projections
  !
  IMPLICIT NONE
  !
  ! atoms
  !
  integer,intent(in) :: natom, atom_index(natom)
  !
  ! atomic wave
  !
  integer,intent(in) :: l, nr
  real(dp),intent(in) :: r(nr), u(nr)
  !
  CHARACTER (len=*) :: filproj
  LOGICAL           :: lwrite_ovp, lbinary
  !
  integer :: nl_tab
  integer,allocatable :: l_tab(:)
  real(dp),allocatable :: tab(:,:)
  !
  integer :: ltab, lp
  !
  INTEGER :: npw, ik, ibnd, i, j, k, na, nb, nt, isym, n,  m, m1, nwfc,&
       nwfc1, lmax_wfc, is, iunproj
  REAL(DP), ALLOCATABLE :: e (:)
  COMPLEX(DP), ALLOCATABLE :: wfcatom (:,:)
  COMPLEX(DP), ALLOCATABLE :: work(:,:),work1(:), proj0(:,:)
  ! Some workspace for k-point calculation ...
  REAL   (DP), ALLOCATABLE :: rwork1(:),rproj0(:,:)
  ! ... or for gamma-point.
  REAL(DP), ALLOCATABLE :: charges(:,:,:), charges_lm(:,:,:,:), proj1 (:)
  REAL(DP) :: psum
  INTEGER  :: nksinit, nkslast
  CHARACTER(len=256) :: filename
  INTEGER, ALLOCATABLE :: idx(:)
  !
  !
  WRITE( stdout, '(/5x,"Calling lprojwave .... ")')
  IF ( gamma_only ) THEN
     WRITE( stdout, '(5x,"gamma-point specific algorithms are used")')
  ENDIF
  !
  allocate( l_tab(3), tab(nqx,3) )
  CALL init_at_dipole( l, nr, r, u, nl_tab, l_tab, tab )

  ! redefine natomwfc
  natomwfc=0
  do ltab=1,nl_tab
    lp=l_tab(ltab)
    natomwfc=natomwfc+2*lp+1
    write(stdout,*) ltab, lp, 2*lp+1, natomwfc
  enddo

  write(stdout,*) ' redefined natomwfc = ', natomwfc

  CALL fill_nlmchi_dipole ( natom, atom_index, nl_tab, l_tab, natomwfc, nwfc, lmax_wfc )
  !
  ! allocations

  ALLOCATE( proj (natomwfc, nbnd, nkstot) )
  ALLOCATE( proj_aux (natomwfc, nbnd, nkstot) )
  proj      = 0.d0
  proj_aux  = (0.d0, 0.d0)
  !
  ALLOCATE(wfcatom (npwx, natomwfc) )
  !
  !    loop on k points
  !
  Write (stdout,*)'ispin,  ik,  ibnd,  atomic_wfc,  real proj.,    imag. proj.,    abs proj.'
  DO ik = 1, nks

     ! <psi_i|
     npw = ngk(ik)
     CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)

     ! |phi_j>
     CALL atomic_wfc_dipole (ik, wfcatom, natom, atom_index, nl_tab, l_tab, tab )
     !
     ! make the projection <psi_i| phi_j>
     !
     IF ( gamma_only ) THEN
        !
        ALLOCATE(rproj0(natomwfc,nbnd), rwork1 (nbnd) )
        CALL calbec ( npw, wfcatom, evc, rproj0)
        !
        proj_aux(:,:,ik) = cmplx( rproj0(:,:), 0.0_dp, kind=dp )
        !
     ELSE
        !
        ALLOCATE(proj0(natomwfc,nbnd), work1 (nbnd) )
        CALL calbec ( npw, wfcatom, evc, proj0)
        !
        proj_aux(:,:,ik) = proj0(:,:)
        !
     ENDIF
     !
        IF ( gamma_only ) THEN
           DO nwfc=1,natomwfc
              DO ibnd=1,nbnd
                 proj(nwfc,ibnd,ik)=abs(rproj0(nwfc,ibnd))**2
              ENDDO
           ENDDO
        ELSE
           DO nwfc=1,natomwfc
              DO ibnd=1,nbnd
                 proj(nwfc,ibnd,ik)=abs(proj0(nwfc,ibnd))**2
                 !SUBHAYAN
                 if (nspin .eq. 1) then
                  Write (stdout,'(A1,1x,I5,1x,I6,1x,I6,3(4x,F8.6))')'1',ik,ibnd,nwfc,real(proj0(nwfc,ibnd)),aimag(proj0(nwfc,ibnd)),proj(nwfc,ibnd,ik)
                 endif
                 if (nspin .eq. 2) then
                  !write(stdout,*)"nks is",nks
                  if (ik .le. int(nks*1.0/2.0)) Write (stdout,'(A1,1x,I5,1x,I6,1x,I6,4x,F8.6,4x,F8.6,4x,F8.6)')'1',ik,ibnd,nwfc,real(proj0(nwfc,ibnd)),aimag(proj0(nwfc,ibnd)),proj(nwfc,ibnd,ik)
                  if (ik .gt. int(nks*1.0/2.0)) Write (stdout,'(A1,1x,I5,1x,I6,1x,I6,4x,F8.6,4x,F8.6,4x,F8.6)')'2',(ik-int(nks*1.0/2.0)),ibnd,nwfc,real(proj0(nwfc,ibnd)),aimag(proj0(nwfc,ibnd)),proj(nwfc,ibnd,ik)
                 endif
              ENDDO
           ENDDO
        ENDIF
     IF ( gamma_only ) THEN
        DEALLOCATE (rwork1)
        DEALLOCATE (rproj0)
     ELSE
        DEALLOCATE (work1)
        DEALLOCATE (proj0)
     ENDIF
     ! on k-points
  ENDDO
  write(stdout,*)'DONE PRINTING COMPLEX PROJECTIONS'
  !
  DEALLOCATE (wfcatom)
  !
  !   vectors et and proj are distributed across the pools
  !   collect data for all k-points to the first pool
  !
  CALL poolrecover (et,       nbnd, nkstot, nks)
  CALL poolrecover (proj,     nbnd * natomwfc, nkstot, nks)
  CALL poolrecover (proj_aux, 2 * nbnd * natomwfc, nkstot, nks)
  !
  IF ( ionode ) THEN
     !
     ! write on the file filproj
     !
     IF (filproj/=' ') THEN
        DO is=1,nspin
           IF (nspin==2) THEN
              IF (is==1) filename=trim(filproj)//'.projwfc_up'
              IF (is==2) filename=trim(filproj)//'.projwfc_down'
              nksinit=(nkstot/2)*(is-1)+1
              nkslast=(nkstot/2)*is
           ELSE
              filename=trim(filproj)//'.projwfc_up'
              nksinit=1
              nkslast=nkstot
           ENDIF
           iunproj=33
           CALL write_io_header(filename, iunproj, title, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, &
                dfftp%nr1, dfftp%nr2, dfftp%nr3, nat, ntyp, ibrav, celldm, at, gcutm, dual,   &
                ecutwfc, nkstot/nspin, nbnd, natomwfc)
           DO nwfc = 1, natomwfc
              WRITE(iunproj,'(2i5,1x,a4,1x,a2,1x,3i5)') &
                  nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
                  nlmchi(nwfc)%els, nlmchi(nwfc)%n, nlmchi(nwfc)%l, nlmchi(nwfc)%m
              DO ik=nksinit,nkslast
                 DO ibnd=1,nbnd
                   WRITE(iunproj,'(2i8,f20.10)') ik,ibnd, &
                                                 abs(proj(nwfc,ibnd,ik))
                 ENDDO
              ENDDO
           ENDDO
           CLOSE(iunproj)
        ENDDO
     ENDIF

     !
     ! write projections to file using iotk
     !
     CALL write_proj( "atomic_proj", lbinary, proj_aux, lwrite_ovp, ovps_aux )
     !
     DEALLOCATE( proj_aux )

     !
     ! write on the standard output file
     !
     WRITE( stdout,'(/5x,"Atomic states used for projection")')
     WRITE( stdout,'( 5x,"(read from pseudopotential files):"/)')
     DO nwfc = 1, natomwfc
        WRITE(stdout,1000) &
             nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
             nlmchi(nwfc)%n, nlmchi(nwfc)%l, nlmchi(nwfc)%m
     ENDDO
1000 FORMAT (5x,"state #",i4,": atom ",i3," (",a3,"), wfc ",i2, &
                " (l=",i1," m=",i2,")")
     !
     ALLOCATE(idx(natomwfc), proj1 (natomwfc) )
     DO ik = 1, nkstot
        WRITE( stdout, '(/" k = ",3f14.10)') (xk (i, ik) , i = 1, 3)
        DO ibnd = 1, nbnd
           WRITE( stdout, '("==== e(",i4,") = ",f11.5," eV ==== ")') &
              ibnd, et (ibnd, ik) * rytoev
           !
           ! sort projections by magnitude, in decreasing order
           !
           DO nwfc = 1, natomwfc
              idx (nwfc) = 0
              proj1 (nwfc) = - proj (nwfc, ibnd, ik)
           ENDDO
           !
           ! projections differing by less than 1.d-4 are considered equal
           !
           CALL hpsort_eps (natomwfc, proj1, idx, eps4)
           !
           !  only projections that are larger than 0.001 are written
           !
           DO nwfc = 1, natomwfc
              proj1 (nwfc) = - proj1(nwfc)
              IF ( abs (proj1(nwfc)) < 0.001d0 ) GOTO 20
           ENDDO
           nwfc = natomwfc + 1
20         nwfc = nwfc -1
           !
           ! fancy (?!?) formatting
           !
           WRITE( stdout, '(5x,"psi = ",5(f5.3,"*[#",i4,"]+"))') &
                (proj1 (i), idx(i), i = 1, min(5,nwfc))
           DO j = 1, (nwfc-1)/5
              WRITE( stdout, '(10x,"+",5(f5.3,"*[#",i4,"]+"))') &
                   (proj1 (i), idx(i), i = 5*j+1, min(5*(j+1),nwfc))
           ENDDO
           psum = SUM ( proj(1:natomwfc, ibnd, ik) )
           WRITE( stdout, '(4x,"|psi|^2 = ",f5.3)') psum
           !
        ENDDO
     ENDDO
     DEALLOCATE (idx, proj1)
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE lprojwave
!
!-----------------------------------------------------------------------
FUNCTION compute_mj(j,l,m)
   !-----------------------------------------------------------------------
   USE kinds, ONLY: DP
   IMPLICIT NONE
   !
   REAL(DP) :: compute_mj, j
   INTEGER  :: l, m

   IF (abs(j-l-0.5d0)<1.d-4) THEN
       compute_mj=m+0.5d0
   ELSEIF (abs(j-l+0.5d0)<1.d-4) THEN
      compute_mj=m-0.5d0
   ELSE
      CALL errore('compute_mj','l and j not compatible',1)
   ENDIF

   RETURN
END FUNCTION compute_mj
!
!-----------------------------------------------------------------------
SUBROUTINE  write_proj (filename, lbinary, projs, lwrite_ovp, ovps )
  !-----------------------------------------------------------------------
  !
  USE kinds
  USE io_files,         ONLY : iun => iunsat, prefix, tmp_dir, postfix
  USE basis,            ONLY : natomwfc
  USE cell_base
  USE klist,            ONLY : wk, xk, nkstot, nelec
  USE noncollin_module, ONLY : noncolin
  USE lsda_mod,         ONLY : nspin, isk
  USE ener,             ONLY : ef
  USE wvfct,            ONLY : et, nbnd
  USE iotk_module
  IMPLICIT NONE

  CHARACTER(*),  INTENT(IN) :: filename
  LOGICAL,       INTENT(IN) :: lbinary
  COMPLEX(DP),   INTENT(IN) :: projs(natomwfc,nbnd,nkstot)
  LOGICAL,       INTENT(IN) :: lwrite_ovp
  COMPLEX(DP),   INTENT(IN) :: ovps(natomwfc,natomwfc,nkstot)
  !
  CHARACTER(256)          :: tmp
  CHARACTER(iotk_attlenx) :: attr
  INTEGER :: ik, ik_eff, isp, ia, ierr, num_k_points

!
! subroutine body
!

  tmp = trim( tmp_dir ) // trim( prefix ) // postfix //trim(filename)
  !
  IF ( lbinary ) THEN
      tmp = TRIM(tmp) // ".dat"
  ELSE
      tmp = TRIM(tmp) // ".xml"
  ENDIF
  !
  CALL iotk_open_write(iun, FILE=trim(tmp), ROOT="ATOMIC_PROJECTIONS", &
                       BINARY=lbinary, IERR=ierr )
  IF ( ierr /= 0 ) RETURN
  !
  !
  num_k_points = nkstot
  IF ( nspin == 2 ) num_k_points = nkstot / 2
  !
  CALL iotk_write_begin(iun, "HEADER")
  !
  CALL iotk_write_dat(iun, "NUMBER_OF_BANDS", nbnd)
  CALL iotk_write_dat(iun, "NUMBER_OF_K-POINTS", num_k_points )
  CALL iotk_write_dat(iun, "NUMBER_OF_SPIN_COMPONENTS", nspin)
  CALL iotk_write_dat(iun, "NON-COLINEAR_CALCULATION",noncolin)
  CALL iotk_write_dat(iun, "NUMBER_OF_ATOMIC_WFC", natomwfc)
  CALL iotk_write_dat(iun, "NUMBER_OF_ELECTRONS", nelec )
  CALL iotk_write_attr(attr, "UNITS", " 2 pi / a", FIRST=.true.  )
  CALL iotk_write_empty (iun,  "UNITS_FOR_K-POINTS", ATTR=attr)
  CALL iotk_write_attr(attr, "UNITS", "Rydberg", FIRST=.true.  )
  CALL iotk_write_empty (iun,  "UNITS_FOR_ENERGY", ATTR=attr)
  CALL iotk_write_dat(iun, "FERMI_ENERGY", ef )
  !
  CALL iotk_write_end(iun, "HEADER")
  !
  !
  CALL iotk_write_dat(iun, "K-POINTS", xk(:,1:num_k_points) , COLUMNS=3 )
  CALL iotk_write_dat(iun, "WEIGHT_OF_K-POINTS", wk(1:num_k_points), COLUMNS=8 )
  !
  CALL iotk_write_begin(iun, "EIGENVALUES")
  !
  DO ik=1,num_k_points
     !
     CALL iotk_write_begin( iun, "K-POINT"//trim(iotk_index(ik)) )
     !
     IF ( nspin == 2 ) THEN
        !
        ik_eff = ik + num_k_points
        !
        CALL iotk_write_dat( iun, "EIG.1", et(:,ik) )
        CALL iotk_write_dat( iun, "EIG.2", et(:,ik_eff) )
        !
     ELSE
        !
        CALL iotk_write_dat( iun, "EIG", et(:,ik) )
        !
     ENDIF
     !
     CALL iotk_write_end( iun, "K-POINT"//trim(iotk_index(ik)) )
     !
  ENDDO
  !
  CALL iotk_write_end(iun, "EIGENVALUES")

  !
  ! main loop atomic wfc
  !
  CALL iotk_write_begin(iun, "PROJECTIONS")
  !
  DO ik=1,num_k_points
     !
     CALL iotk_write_begin( iun, "K-POINT"//trim(iotk_index(ik)) )
     !
     IF ( nspin == 2 ) THEN
        !
        CALL iotk_write_begin ( iun, "SPIN.1" )
           !
           DO ia = 1, natomwfc
               CALL iotk_write_dat(iun, "ATMWFC"//trim(iotk_index(ia)), projs(ia,:,ik)  )
           ENDDO
           !
        CALL iotk_write_end ( iun, "SPIN.1" )
        !
        ik_eff = ik + num_k_points
        !
        CALL iotk_write_begin ( iun, "SPIN.2" )
           !
           DO ia = 1, natomwfc
               CALL iotk_write_dat(iun, "ATMWFC"//trim(iotk_index(ia)), projs(ia,:,ik_eff)  )
           ENDDO
           !
        CALL iotk_write_end ( iun, "SPIN.2" )
        !
     ELSE
        !
        DO ia = 1,natomwfc
            CALL iotk_write_dat(iun, "ATMWFC"//trim(iotk_index(ia)), projs(ia,:,ik)  )
        ENDDO
        !
     ENDIF
     !
     CALL iotk_write_end( iun, "K-POINT"//trim(iotk_index(ik)) )
     !
  ENDDO
  !
  CALL iotk_write_end(iun, "PROJECTIONS")

  !
  ! overlaps
  !
  IF ( lwrite_ovp ) THEN
      !
      CALL iotk_write_begin(iun, "OVERLAPS")
      !
      DO ik=1,num_k_points
          !
          CALL iotk_write_begin( iun, "K-POINT"//trim(iotk_index(ik)) )
          !
          DO isp = 1, nspin
              !
              ik_eff = ik + num_k_points * ( isp -1 )
              !
              CALL iotk_write_dat(iun, "OVERLAP"//trim(iotk_index(isp)), ovps(:,:,ik_eff)  )
              !
              !
          ENDDO
          !
          CALL iotk_write_end( iun, "K-POINT"//trim(iotk_index(ik)) )
          !
      ENDDO
      !
      CALL iotk_write_end(iun, "OVERLAPS")
      !
  ENDIF
  !
  ! closing the file
  !
  CALL iotk_close_write(iun)

END SUBROUTINE write_proj
!
    SUBROUTINE fill_nlmchi_dipole ( natom, atom_index, nl_tab, l_tab, natomwfc, nwfc, lmax_wfc )
      !
      USE kinds,      ONLY : DP
      USE projections
      USE ions_base, ONLY : ityp, nat
      USE uspp_param, ONLY: upf
      USE spin_orb, ONLY: lspinorb
      USE noncollin_module, ONLY: noncolin
      !
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: natom, atom_index(natom), nl_tab, l_tab(3)
      integer :: natomwfc
      INTEGER, INTENT (OUT) :: nwfc, lmax_wfc 
      !
      integer :: iatom
      INTEGER :: na, nt, n, n1, n2, l, m, ind
      REAL(dp) :: jj, fact(2)
      REAL(dp), EXTERNAL :: spinor
      CHARACTER(LEN=1) :: label
      CHARACTER(LEN=1) :: spdf(0:3) = ['S','P','D','F']
      INTEGER :: nn(0:3)
      !
      ALLOCATE (nlmchi(natomwfc))
      nwfc=0
      lmax_wfc = 0
      DO iatom = 1, natom
         na = atom_index(iatom)
         nt = ityp (na)
         n2 = 0
         nn = [1,2,3,4]
         DO n = 1, nl_tab
               l = l_tab(n)
               write(label,'(A1)') spdf(l)
               lmax_wfc = max (lmax_wfc, l )
               IF (lspinorb) THEN
                  IF (upf(nt)%has_so) THEN
                     jj = upf(nt)%jchi (n)
                     ind = 0
                     DO m = -l-1, l
                        fact(1) = spinor(l,jj,m,1)
                        fact(2) = spinor(l,jj,m,2)
                        IF (abs(fact(1)) > 1.d-8 .or. abs(fact(2)) > 1.d-8) THEN
                           nwfc = nwfc + 1
                           ind = ind + 1
                           nlmchi(nwfc)%na = na
                           nlmchi(nwfc)%n  =  n
                           nlmchi(nwfc)%l  =  l
                           nlmchi(nwfc)%m  =  m
                           nlmchi(nwfc)%ind = ind
                           nlmchi(nwfc)%jj  = jj
                           nlmchi(nwfc)%els = label
                        ENDIF
                     ENDDO
                  ELSE
                     DO n1 = l, l+1
                        jj= dble(n1) - 0.5d0
                        ind = 0
                        IF (jj>0.d0)  THEN
                           n2 = n2 + 1
                           DO m = -l-1, l
                              fact(1) = spinor(l,jj,m,1)
                              fact(2) = spinor(l,jj,m,2)
                              IF (abs(fact(1)) > 1.d-8 .or. abs(fact(2)) > 1.d-8) THEN
                                 nwfc = nwfc + 1
                                 ind = ind + 1
                                 nlmchi(nwfc)%na = na
                                 nlmchi(nwfc)%n  =  n2
                                 nlmchi(nwfc)%l  =  l
                                 nlmchi(nwfc)%m  =  m
                                 nlmchi(nwfc)%ind = ind
                                 nlmchi(nwfc)%jj  = jj
                                 nlmchi(nwfc)%els = label
                              ENDIF
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDIF
               ELSE
                  DO m = 1, 2 * l + 1
                     nwfc=nwfc+1
                     nlmchi(nwfc)%na = na
                     nlmchi(nwfc)%n  =  n
                     nlmchi(nwfc)%l  =  l
                     nlmchi(nwfc)%m  =  m
                     nlmchi(nwfc)%ind=  m
                     nlmchi(nwfc)%jj =  0.d0
                     nlmchi(nwfc)%els=  label
                  ENDDO
                  IF ( noncolin) THEN
                     DO m = 1, 2 * l + 1
                        nwfc=nwfc+1
                        nlmchi(nwfc)%na = na
                        nlmchi(nwfc)%n  =  n
                        nlmchi(nwfc)%l  =  l
                        nlmchi(nwfc)%m  =  m
                        nlmchi(nwfc)%ind=  m+2*l+1
                        nlmchi(nwfc)%jj =  0.d0
                        nlmchi(nwfc)%els = label
                     END DO
                  ENDIF
               ENDIF
         ENDDO
      ENDDO
      !
      IF (lmax_wfc > 3) CALL errore ('fill_nlmchi_dipole', 'l > 3 not yet implemented',1)
      IF (nwfc /= natomwfc) CALL errore ('fill_nlmchi_dipole','wrong # of atomic wfcs',1)
      
    END SUBROUTINE fill_nlmchi_dipole

