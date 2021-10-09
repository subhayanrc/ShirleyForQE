! ---------------------------------------------------------------------- 
  module hamq_shirley
! ---------------------------------------------------------------------- 

  implicit none

  private

  public :: nbasis
  public :: alloc_hamq_local, init_kin, &
            init_vnl_kgrid, get_vnl_kgrid_cart, &
            init_vnl_k, init_vnl_spline, &
            dump_hamq, init_stdout, init_qtrans, &
            build_hamq_local, build_hamq_projs, build_hamq_nlprojs, &
            build_dhamq_local, build_dhamq_projs, &
            build_hamq_kin, build_hamq_vloc, &
            read_hamq, write_hamq, bcast_hamdim, &
            trnlp2kin, trkin2nlp, &
            mom_hamq
  public :: build_hamq_atomprojs, resize_vnl
  public :: open_hamq, close_hamq, write_hamloc, write_hamprj, write_haminf, &
            read_hamloc, read_hamprj
  public :: matrix_list
  public :: natom, ntype, nproj, nproj_nl, &
            type_atom, nproj_type, nproj_type_nl, &
            index_betaq, index_nlproj_type
  public :: update_atomic_proj, write_nloper, read_nloper, ncpp
  public :: init_system, dump_system, init_ldaU
  public :: fhhaminf, fhhamloc, fhhamprj
  public :: lda_plus_u, natomproj, &
            Hubbard_lmax, Hubbard_l, Hubbard_U, Hubbard_alpha, &
            index_ldaUq, &
            update_atomic_proj_ldaU

  public :: kin0, kin1, vloc

  interface init_kin
    module procedure init_kin, init_kin_packed
  end interface init_kin
  interface dump_hamq
    module procedure dump_hamq, dump_hamq_packed
  end interface dump_hamq


  ! necessary parameters and variables
  integer,parameter :: dp = kind(1.d0)
  real(dp),parameter :: MByte = 2.d0**20
  complex(dp),parameter :: zero=cmplx(0.d0,0.d0)
  complex(dp),parameter :: one =cmplx(1.d0,0.d0)
  complex(dp),parameter :: two =cmplx(2.d0,0.d0)
  complex(dp),parameter :: half=cmplx(0.5d0,0.0d0)
  real(dp),parameter :: rytoev = 13.60569193d0

  ! i/o
  integer :: stdout = 6

  ! q-vector transformations between the kinetic energy Cartesian basis
  ! and the non-local potential crystal basis
  real(dp) :: trkin2nlp(3,3)
  real(dp) :: trnlp2kin(3,3)

  integer :: nbasis      ! no. of basis functions
  integer :: nproj       ! no. of projectors
  !
  integer :: nproj_nl    ! no. of strictly non-local projectors
  !
  ! details for non-local operators based on atom-centers
  !
  integer :: natom   ! no. of atoms
  integer :: ntype   ! no. of type of atoms
  integer :: nproj_type_max   ! max no. of projectors over all types of atoms
  !
  integer,allocatable :: type_atom(:)      ! the type of each atom
  integer,allocatable :: nproj_type(:)     ! no. of projectors of a given type
  integer,allocatable :: nproj_type_nl(:)  ! no. of strictly non-local projectors of a given type
  !
  ! typedefs for non-local operators
  type matrix_list
    real(dp),pointer :: matrix(:,:)
  end type
  
  ! flag for norm-conserving pseudos
  logical :: ncpp

  ! which projectors are strictly non-local by type
  integer,allocatable :: index_nlproj_type(:,:)
  ! mapping from atoms and projections to betaq
  integer,allocatable :: index_betaq(:,:)

  ! indices of those elements of beta that are strictly non-local
  integer,allocatable :: index_nlproj_betaq(:)

  ! system details - useful for normalization etc
  real(dp) :: nelec
  real(dp) :: alat, omega
  real(dp) :: at(3,3), bg(3,3)
  real(dp) :: tpiba
  integer :: nspin

  ! LDA+U
  logical :: lda_plus_u
  !
  integer :: natomproj   ! no. of atomic projectors - for LDA+U right now
  !
  integer ::  Hubbard_lmax
  integer,allocatable ::  Hubbard_l(:)
  real(dp),allocatable :: Hubbard_U(:)
  real(dp),allocatable :: Hubbard_alpha(:)
  integer,allocatable :: index_ldaUq(:,:)

  ! coefficients of q-dependent expansion
  complex(dp),allocatable :: kin0(:,:)
  complex(dp),allocatable :: kin1(:,:,:)
  complex(dp),allocatable :: vloc(:,:,:)

  ! non-local potential
  ! the grid in k-space of matrix elements of the non-local potential
  ! used to determine the B-spline interpolation
  complex(dp),allocatable :: vnl(:,:,:)
  ! grid in k-space
  integer :: nkx, nky, nkz, nktot
  integer,allocatable :: ikmap(:,:), ikmap_inv(:,:,:)
  real(dp),allocatable :: xkx(:), xky(:), xkz(:)
  ! splined terms
  integer :: kxord, kyord, kzord
  real(dp),allocatable :: xknot(:), yknot(:), zknot(:)
  integer :: nxcoef, nycoef, nzcoef
  real(dp),allocatable :: bscoefR(:,:,:), bscoefI(:,:,:)
  real(dp),allocatable :: ascoefR(:,:,:), ascoefI(:,:,:)

  ! file handles
  integer :: fhhamloc
  integer :: fhhamprj
  integer :: fhhaminf


  contains


! ---------------------------------------------------------------------- 
  subroutine init_stdout( stdout_ )
! ---------------------------------------------------------------------- 

  integer,intent(in) :: stdout_

  ! replace existing stdout unit with the input
  stdout = stdout_

  return
  end subroutine init_stdout


! ---------------------------------------------------------------------- 
  subroutine init_system( nelec_, alat_, omega_, at_, bg_, tpiba_, &
                          nspin_, lda_plus_u_ )
! ---------------------------------------------------------------------- 

  real(dp) :: nelec_, alat_, omega_, at_(3,3), bg_(3,3), tpiba_
  integer :: nspin_
  logical :: lda_plus_u_

  nelec = nelec_
  alat = alat_
  omega = omega_
  at = at_
  bg = bg_
  tpiba = tpiba_
  nspin = nspin_
  lda_plus_u = lda_plus_u_

  end subroutine init_system

! ---------------------------------------------------------------------- 
  subroutine dump_system( nelec_, alat_, omega_, at_, bg_, tpiba_, &
                          nspin_, lda_plus_u_ )
! ---------------------------------------------------------------------- 

  real(dp),intent(out) :: nelec_, alat_, omega_, at_(3,3), bg_(3,3), tpiba_
  integer :: nspin_
  logical :: lda_plus_u_

  nelec_ = nelec
  alat_ = alat
  omega_ = omega
  at_ = at
  bg_ = bg
  tpiba_ = tpiba
  nspin_ = nspin
  lda_plus_u_ = lda_plus_u

  end subroutine dump_system

! ---------------------------------------------------------------------- 
  subroutine init_ldaU( ntype_, Hubbard_l_, Hubbard_U_, Hubbard_alpha_ )
! ---------------------------------------------------------------------- 

  integer :: ntype_
  integer :: Hubbard_l_(ntype_)
  real(dp) :: Hubbard_U_(ntype_)
  real(dp) :: Hubbard_alpha_(ntype_)

  allocate( Hubbard_l(ntype_), Hubbard_U(ntype_), Hubbard_alpha(ntype_) )
  Hubbard_l = Hubbard_l_
  Hubbard_U = Hubbard_U_
  Hubbard_alpha = Hubbard_alpha_

  end subroutine init_ldaU


! ---------------------------------------------------------------------- 
  subroutine init_qtrans( trkin2nlp_, trnlp2kin_ )
! ---------------------------------------------------------------------- 

  real(dp),intent(in) ::  trkin2nlp_(3,3), trnlp2kin_(3,3)

  trkin2nlp = trkin2nlp_
  trnlp2kin = trnlp2kin_

  return
  end subroutine init_qtrans


! ---------------------------------------------------------------------- 
  subroutine alloc_hamq_local( nr, nc )
! ---------------------------------------------------------------------- 

  integer,intent(in) :: nr, nc
  integer :: ierr

  write(stdout,*) ' local ham: attempt to allocate ', &
    (dble(nr)*dble(nc)*5.d0)*16.d0/MByte, ' MBytes'

  if( allocated( vloc ) ) deallocate( vloc )
  allocate( vloc(nr,nc,nspin), stat=ierr )
  if( ierr/=0 ) call errore('alloc_hamq_local','unable to allocate vloc',1)

  if( allocated( kin0 ) ) deallocate( kin0 )
  allocate( kin0(nr,nc), stat=ierr )
  if( ierr/=0 ) call errore('alloc_hamq_local','unable to allocate kin0',1)

  if( allocated( kin1 ) ) deallocate( kin1 )
  allocate( kin1(nr,nc,3), stat=ierr )
  if( ierr/=0 ) call errore('alloc_hamq_local','unable to allocate kin1',1)

  return
  end subroutine alloc_hamq_local


! ---------------------------------------------------------------------- 
  subroutine alloc_knots( )
! ---------------------------------------------------------------------- 

  integer :: ierr

  write(stdout,*) ' non-local splines: attempt to allocate ', &
    ( dble(nxcoef+kxord) + dble(nycoef+kyord) + dble(nzcoef+kzord) &
    )*8.d0/MByte, ' MBytes'

  if( allocated(xknot) ) deallocate(xknot)
  if( allocated(yknot) ) deallocate(yknot)
  if( allocated(zknot) ) deallocate(zknot)

  allocate( xknot(nxcoef+kxord), &
            yknot(nycoef+kyord), &
            zknot(nzcoef+kzord), &
            stat=ierr )
  if( ierr /= 0 ) call errore('alloc_knots','unable to allocate knots',1)

  end subroutine alloc_knots


! ---------------------------------------------------------------------- 
  subroutine alloc_projs( n )
! ---------------------------------------------------------------------- 

  integer,intent(in) :: n
  integer :: ierr


  write(stdout,*) ' non-local splines: attempt to allocate ', &
    ( dble(nxcoef*nycoef*nzcoef)*dble(nproj)*dble(n)*2.d0 )*8.d0/MByte, ' MBytes'

  if( allocated(bscoefR) ) deallocate(bscoefR)
  if( allocated(bscoefI) ) deallocate(bscoefI)

  allocate( bscoefR(nxcoef*nycoef*nzcoef,nproj,n), &
            bscoefI(nxcoef*nycoef*nzcoef,nproj,n), &
            stat=ierr )
  if( ierr /= 0 ) call errore('alloc_projs','unable to allocate projector coefficients',1)

  write(stdout,*) 'alloc_projs: projector coefficients allocated'

  end subroutine alloc_projs


! ---------------------------------------------------------------------- 
  subroutine alloc_projs_ldaU( n )
! ---------------------------------------------------------------------- 

  integer,intent(in) :: n
  integer :: ierr


  write(stdout,*) ' non-local splines: attempt to allocate ', &
    ( dble(nxcoef*nycoef*nzcoef)*dble(natomproj)*dble(n)*2.d0 )*8.d0/MByte, ' MBytes'

  if( allocated(ascoefR) ) deallocate(ascoefR)
  if( allocated(ascoefI) ) deallocate(ascoefI)

  allocate( ascoefR(nxcoef*nycoef*nzcoef,natomproj,n), &
            ascoefI(nxcoef*nycoef*nzcoef,natomproj,n), &
            stat=ierr )
  if( ierr /= 0 ) call errore('alloc_projs_ldaU','unable to allocate ldaU projector coefficients',1)

  end subroutine alloc_projs_ldaU


! ---------------------------------------------------------------------- 
  subroutine alloc_atomic_proj( nproj_, nproj_nl_, &
                                natom_, ntype_, nproj_type_max_ )
! ---------------------------------------------------------------------- 

  integer,intent(in) :: nproj_
  integer,intent(inout) :: nproj_nl_
  integer,intent(in) :: natom_, ntype_, nproj_type_max_

  integer :: i, it, ierr

  if( nproj_ < 0 ) then
    call errore('alloc_atomic_proj','strange number of projectors',1)
  endif

  if( nproj_nl_ < 0 .or. nproj_nl_ > nproj_ ) nproj_nl_ = nproj_

  ! update module dimensions
  nproj = nproj_
  nproj_nl = nproj_nl_

  natom = natom_
  ntype = ntype_
  nproj_type_max = nproj_type_max_

  write(stdout,*) ' non-local projs: attempt to allocate ', &
    ( dble(natom_) + 2.d0*dble(ntype_) + dble(nproj_type_max_)*(dble(ntype_) &
    + dble(natom_)) + dble(nproj_nl) )*4.d0/MByte, ' MBytes'

  ! allocations
  if( allocated(type_atom) ) deallocate(type_atom)
  if( allocated(nproj_type) ) deallocate(nproj_type)
  if( allocated(nproj_type_nl) ) deallocate(nproj_type_nl)
 
  allocate( type_atom(natom_) )
  allocate( nproj_type(ntype_) )
  allocate( nproj_type_nl(ntype_) )

  if( allocated(index_nlproj_type) ) deallocate( index_nlproj_type )
  allocate( index_nlproj_type(nproj_type_max_,ntype_) )

  if( allocated(index_betaq) ) deallocate( index_betaq )
  allocate( index_betaq(nproj_type_max_,natom_) )

  if( allocated( index_nlproj_betaq ) ) deallocate( index_nlproj_betaq )
  allocate( index_nlproj_betaq(nproj_nl), stat=ierr )
  if( ierr/=0 ) call errore('alloc_atomic_proj','unable to allocate index_nlproj_betaq',1)

  end subroutine alloc_atomic_proj


! ---------------------------------------------------------------------- 
  subroutine alloc_atomic_proj_ldaU( natom_, Hubbard_lmax_ )
! ---------------------------------------------------------------------- 

  integer,intent(in) :: natom_, Hubbard_lmax_

  if( allocated(index_ldaUq) ) deallocate( index_ldaUq )
  allocate( index_ldaUq(2*Hubbard_lmax_+1,natom_) )

  end subroutine alloc_atomic_proj_ldaU


! ---------------------------------------------------------------------- 
  subroutine init_kin( nbasis_, lda, kin_, ixyz )
! ---------------------------------------------------------------------- 

  integer,intent(in) :: nbasis_, lda
  integer,intent(in),optional :: ixyz
  complex(dp),intent(in) :: kin_(lda,nbasis_)

  integer :: i, j

  if( nbasis_ /= nbasis ) call errore('init_kin','inconsistent size: nbasis',1)

  if( present(ixyz) ) then

    if( ixyz < 1 .or. ixyz > 3 ) call errore('init_kin','inconsistent size: ixyz',1)

    forall( i=1:nbasis, j=1:nbasis ) kin1(i,j,ixyz) = kin_(i,j)

  else

    forall( i=1:nbasis, j=1:nbasis ) kin0(i,j) = kin_(i,j)

  endif

  return
  end subroutine init_kin

! ---------------------------------------------------------------------- 
  subroutine init_kin_packed( nbasis_, kin_, ixyz )
! ---------------------------------------------------------------------- 

  integer,intent(in) :: nbasis_
  integer,intent(in),optional :: ixyz
  complex(dp),intent(in) :: kin_((nbasis_*(nbasis_+1))/2)

  integer :: i, j, ij

  if( nbasis_ /= nbasis ) call errore('init_kin','inconsistent size: nbasis',1)

  if( present(ixyz) ) then

    ij=0
    do j=1,nbasis
      do i=1,j
        ij=ij+1
        kin1(i,j,ixyz) = kin_(ij)
      enddo
    enddo

    do j=1,nbasis
      do i=j+1,nbasis
        kin1(i,j,ixyz) = conjg(kin1(j,i,ixyz))
      enddo
    enddo

  else

    ij=0
    do j=1,nbasis
      do i=1,j
        ij=ij+1
        kin0(i,j) = kin_(ij)
      enddo
    enddo

    do j=1,nbasis
      do i=j+1,nbasis
        kin0(i,j) = conjg(kin0(j,i))
      enddo
    enddo
  
  endif

  return
  end subroutine init_kin_packed


! ---------------------------------------------------------------------- 
  subroutine update_atomic_proj( nkb, nat, ntyp, ityp, nhm, nh, nhtol, local_channel )
! ---------------------------------------------------------------------- 

  ! store the index mapping from becp or vkb to atoms and projections 
  ! also store local info on atoms and projectors

  integer,intent(in) :: nkb
  integer,intent(in) :: nat, ntyp
  integer,intent(in) :: ityp(nat)
  integer,intent(in) :: nhm
  integer,intent(in) :: nh(ntyp), nhtol(nhm,ntyp)
  integer,intent(in) :: local_channel(ntyp)

  integer :: ijkb0, nt, na, ikb
  integer :: ip, ipnl, ikbnl

  ! assign local copies of these variables
  nproj = nkb

  ikb=0
  ikbnl=0
  do nt=1,ntyp
    do na=1,nat
      if( ityp(na) == nt ) then
        do ip=1,nh(nt)
          ikb=ikb+1
          if( nhtol(ip,nt) /= local_channel(nt) ) then
            ikbnl=ikbnl+1
          endif
        enddo
      endif
    enddo
  enddo

  ! check
  if( ikb /= nproj ) call errore('update_atomic_proj','number of projectors inconsistent with atoms and types',1)

  nproj_nl = ikbnl

  natom = nat
  ntype = ntyp
  nproj_type_max = nhm

  call alloc_atomic_proj( nproj, nproj_nl, natom, ntype, nproj_type_max )

  type_atom = ityp
  nproj_type = nh

  index_nlproj_type=0
  do nt=1,ntype
    ipnl=0
    nproj_type_nl(nt) = 0
    do ip=1,nproj_type(nt)
      if( nhtol(ip,nt) == local_channel(nt) ) cycle
      nproj_type_nl(nt) = nproj_type_nl(nt) + 1
      ipnl=ipnl+1
      index_nlproj_type(ipnl,nt) = ip
    enddo
  enddo

  index_betaq=0
  index_nlproj_betaq=0
  ijkb0=0
  ikbnl=0
  do nt=1,ntype
    do na=1,natom
      if( ityp(na) == nt ) then
        do ip=1,nproj_type(nt)
          ikb=ijkb0+ip
          index_betaq(ip,na) = ikb
          if( nhtol(ip,nt) /= local_channel(nt) ) then
            ikbnl=ikbnl+1
            index_nlproj_betaq(ikbnl) = ikb
          endif
        enddo
        ijkb0 = ijkb0 + nproj_type(nt)
      endif
    enddo
  enddo

  end subroutine update_atomic_proj


! ---------------------------------------------------------------------- 
  subroutine update_atomic_proj_ldaU( natomwfc, nat, ntyp, ityp, Hubbard_lmax_, Hubbard_l_, Hubbard_U_, Hubbard_alpha_ )
! ---------------------------------------------------------------------- 

  ! store the index mapping from becp or vkb to atoms and projections 
  ! also store local info on atoms and projectors

  integer,intent(in) :: natomwfc
  integer,intent(in) :: nat, ntyp
  integer,intent(in) :: ityp(nat)
  integer,intent(in) :: Hubbard_lmax_
  integer,intent(in) :: Hubbard_l_(ntyp)
  real(dp),intent(in) :: Hubbard_U_(ntyp), Hubbard_alpha_(ntyp)

  integer :: nt, na
  integer :: ikb, ijkb0, ip 

  ! copy Hubbard parameters
  Hubbard_lmax=Hubbard_lmax_
  allocate( Hubbard_l(ntyp), Hubbard_U(ntyp), Hubbard_alpha(ntyp) )
  Hubbard_l=Hubbard_l_
  Hubbard_U=Hubbard_U_
  Hubbard_alpha=Hubbard_alpha_

  ! how many LDA+U projectors out of the total number of atomic natomwfc
  natomproj=0
  do nt=1,ntyp
    do na=1,nat
      if( ityp(na) == nt ) then
        if( Hubbard_U(nt) /= 0.d0 .or. Hubbard_alpha(nt) /= 0.d0 ) then
          do ip=1,2*Hubbard_l(nt)+1
            natomproj=natomproj+1
          enddo
        endif
      endif
    enddo
  enddo

  ! check
  if( natomproj <= natomwfc ) then
    write(stdout,*) ' update_atomic_proj_ldaU: note that number of projectors for LDA+U is ', natomproj
    write(stdout,*) ' update_atomic_proj_ldaU: compare with the number of atomic projectors: ', natomwfc
  else
    call errore('update_atomic_proj_ldaU','number of LDA+U projectors inconsistent with atoms and types',1)
  endif

  call alloc_atomic_proj_ldaU( natom, Hubbard_lmax )

  index_ldaUq=0
  ijkb0=0
  do nt=1,ntype
    do na=1,natom
      if( ityp(na) == nt ) then
        if( Hubbard_U(nt) /= 0.d0 .or. Hubbard_alpha(nt) /= 0.d0 ) then
          do ip=1,2*Hubbard_l(nt)+1
            ikb=ijkb0+ip
            index_ldaUq(ip,na) = ikb
          enddo
          ijkb0 = ijkb0 + 2*Hubbard_l(nt)+1
        endif
      endif
    enddo
  enddo

  end subroutine update_atomic_proj_ldaU


! ---------------------------------------------------------------------- 
  subroutine init_vnl_kgrid( nkgrid, ikgrid, nkb, nbnd, nktot_out )
! ---------------------------------------------------------------------- 

  integer,intent(in) :: nkgrid(3), ikgrid(3)
  integer,intent(in) :: nkb, nbnd
  integer,intent(out) :: nktot_out

  integer :: i, j, k, ierr
  integer :: nkgrid_(3)

  if( any( nkgrid < 0 ) ) &
    call errore('hamq','nkgrid - number of grid points must be non-negative',1)

  if( any( ikgrid < 0 ) .or. any( ikgrid > 1 ) ) &
    call errore('hamq','ikgrid - grid shift values possible are 0 or 1',1)

  nkx = nkgrid(1) + ikgrid(1) + 1
  nky = nkgrid(2) + ikgrid(2) + 1
  nkz = nkgrid(3) + ikgrid(3) + 1
  nktot = nkx*nky*nkz

  allocate( xkx(nkx), xky(nky), xkz(nkz), &
            ikmap(3,nktot), ikmap_inv(nkx,nky,nkz) )

  nkgrid_=nkgrid
  where( nkgrid==0 ) nkgrid_=1 

  nktot=0
  do i=1,nkx
     do j=1,nky
        do k=1,nkz
           !  xk are the components of the complete grid in crystal axis
           xkx(i) = DBLE(i-1)/nkgrid_(1) - DBLE(ikgrid(1))/2/nkgrid_(1)
           xky(j) = DBLE(j-1)/nkgrid_(2) - DBLE(ikgrid(2))/2/nkgrid_(2)
           xkz(k) = DBLE(k-1)/nkgrid_(3) - DBLE(ikgrid(3))/2/nkgrid_(3)
           nktot=nktot+1
           ikmap(1:3,nktot) = (/ i, j, k /)
           ikmap_inv(i,j,k) = nktot
        end do
     end do
  end do
 
  allocate( vnl(nktot,nkb,nbnd), stat=ierr )
  if( ierr/=0 ) call errore('init_vnl_kgrid','unable to allocate vnl',1)
  write(stdout,*) '  nkb = ', nkb, '   nproj = ', nproj
  write(stdout,*) ' nbnd = ', nbnd, ' nbasis = ', nbasis

  nktot_out = nktot

  return
  end subroutine init_vnl_kgrid


! ---------------------------------------------------------------------- 
  subroutine resize_vnl( nk, nkb, nbnd )
! ---------------------------------------------------------------------- 

  integer,intent(in) :: nk, nkb, nbnd

  if( allocated(vnl) ) deallocate( vnl )
  allocate( vnl(nk,nkb,nbnd) )

  end subroutine resize_vnl


! ---------------------------------------------------------------------- 
  subroutine get_vnl_kgrid_cart( xk_cart )
! ---------------------------------------------------------------------- 

  real(dp),intent(out) :: xk_cart(3,nktot)

  real(dp) :: xk(3)
  integer :: ik

  write(stdout,'(/,a6,3a12,2x,3a12,/)') 'kpt', 'xk(1)','xk(2)','xk(3)', &
                                     'xk(x)','xk(y)','xk(z)'
  do ik=1,nktot
    xk = (/ xkx(ikmap(1,ik)), xky(ikmap(2,ik)), xkz(ikmap(3,ik)) /)
    xk_cart(:,ik) = matmul( trnlp2kin, xk )
    !
    write(stdout,'(i6,3f12.6,2x,3f12.6)') ik, xkx(ikmap(1,ik)), xky(ikmap(2,ik)), xkz(ikmap(3,ik)), xk_cart(1:3,ik)
  enddo
  write(stdout,*)

  return
  end subroutine get_vnl_kgrid_cart


! ---------------------------------------------------------------------- 
  subroutine init_vnl_k( nproj_, nbasis_, ik, vnl_ )
! ---------------------------------------------------------------------- 

  integer,intent(in) :: nproj_, nbasis_
  integer,intent(in) :: ik
  complex(dp),intent(in) :: vnl_(:,:)

  integer :: i, j

  if( nproj_ /= size(vnl,2) ) call errore('init_vnl_k','inconsistent size: nproj',1)
  if( nbasis_ /= size(vnl,3) ) call errore('init_vnl_k','inconsistent size: nbasis',1)

  forall( i=1:nproj_, j=1:nbasis_ ) vnl(ik,i,j) = vnl_(i,j)

  return
  end subroutine init_vnl_k

! ---------------------------------------------------------------------- 
  subroutine init_vnl_spline( kxord_, kyord_, kzord_ )
! ---------------------------------------------------------------------- 

  use bspline90_22, only : dbsnak, dbs3in

  integer,intent(in) :: kxord_, kyord_, kzord_

  integer :: ierr, i, j
  real(dp),allocatable :: fdata(:,:,:)
  real(dp),allocatable :: xkxp(:), xkyp(:), xkzp(:)
  complex(dp),allocatable :: vnlp(:,:,:)
  integer,allocatable :: ixp(:), iyp(:), izp(:)
  integer :: ic, jc, kc

  ! take input
  kxord = kxord_
  kyord = kyord_
  kzord = kzord_

  ! define the number of coefficients
  nxcoef = nkx
  nycoef = nky
  nzcoef = nkz

  ! defaults 
  if( kxord == 0 ) kxord = nxcoef
  if( kyord == 0 ) kyord = nycoef
  if( kzord == 0 ) kzord = nzcoef

  ! error check
  if( kxord < 0 ) call errore('init_vnl_spline','order of x-spline < 1',1)
  if( kyord < 0 ) call errore('init_vnl_spline','order of y-spline < 1',1)
  if( kzord < 0 ) call errore('init_vnl_spline','order of z-spline < 1',1)

  ! warnings
  if( kxord > nxcoef ) call errore('init_vnl_spline','order of x-spline > no. x-points. Reducing',-1)
  if( kyord < nycoef ) call errore('init_vnl_spline','order of y-spline < no. y-points. Reducing',-1)
  if( kzord < nzcoef ) call errore('init_vnl_spline','order of z-spline < no. z-points. Reducing',-1)

  kxord = min(kxord,nxcoef)
  kyord = min(kyord,nycoef)
  kzord = min(kzord,nzcoef)

  write(stdout,*) 'init_vnl_spline:', kxord, kyord, kzord

  ! important copying of basis dimension for projectors
  if( nproj /= size( vnl, 2 ) ) call errore('init_vnl_spline','nproj inconsistent with vnl',1)

  call alloc_knots()
  call alloc_projs( size(vnl,3) )

  allocate( xkxp(nxcoef) )
  allocate( xkyp(nycoef) )
  allocate( xkzp(nzcoef) )

  xkxp = xkx(1:nkx)
  xkyp = xky(1:nky)
  xkzp = xkz(1:nkz)
  call dbsnak( nxcoef, xkxp, kxord, xknot )
  call dbsnak( nycoef, xkyp, kyord, yknot )
  call dbsnak( nzcoef, xkzp, kzord, zknot )

  write(stdout,*) 'init_vnl_spline: knots made'

  allocate( fdata(nxcoef,nycoef,nzcoef), stat=ierr )
  if( ierr /= 0 ) call errore('init_vnl_spline','unable to allocate fdata',1)

  write(stdout,*) 'init_vnl_spline: fdata allocated'

  write(stdout,*) ' splining '
  do j=1,size(vnl,3)
  do i=1,nproj
    do kc=1,nzcoef
    do jc=1,nycoef
    do ic=1,nxcoef
      fdata(ic,jc,kc) =  real( vnl( ikmap_inv(ic,jc,kc), i, j ) )
    enddo
    enddo
    enddo
  
    call dbs3in( nxcoef, xkxp, nycoef, xkyp, nzcoef, xkzp, fdata, nxcoef, nycoef, &
                 kxord, kyord, kzord, xknot, yknot, zknot, bscoefR(1,i,j) )

    do kc=1,nzcoef
    do jc=1,nycoef
    do ic=1,nxcoef
      fdata(ic,jc,kc) = aimag( vnl( ikmap_inv(ic,jc,kc), i, j ) )
    enddo
    enddo
    enddo
  
    call dbs3in( nxcoef, xkxp, nycoef, xkyp, nzcoef, xkzp, fdata, nxcoef, nycoef, &
                 kxord, kyord, kzord, xknot, yknot, zknot, bscoefI(1,i,j) )
  enddo
  enddo

  deallocate( fdata )

  return
  end subroutine init_vnl_spline

! ---------------------------------------------------------------------- 
  subroutine build_hamq_kin( qvec, lkin, hamq_local )
! ---------------------------------------------------------------------- 

  use hamq_pool, only : local_cyclic_dims, ndiag_cyclic, diag_cyclic

  real(dp),intent(in) :: qvec(3)
  logical,intent(in) :: lkin  ! true (false) if qvec is in kinetic (nonlocal) basis
  complex(dp),intent(out) :: hamq_local(:,:)

  integer :: i, j, ierr
  integer :: nr, nc
  real(dp) :: tqvec(3), sqvec(3)
  complex(dp) :: zqvec(3), zqvec2

  if( lkin ) then
    tqvec = qvec
    sqvec = matmul( trkin2nlp, qvec )
  else
    tqvec = matmul( trnlp2kin, qvec )
    sqvec = qvec
  endif
  
  ! reduce sqvec to [0,1)
  sqvec = sqvec - floor( sqvec )
  ! redefine consistent cartesian coordinates
  tqvec = matmul( trnlp2kin, sqvec )

  zqvec = cmplx( tqvec, kind=dp )
  zqvec2 = cmplx( dot_product(tqvec,tqvec), kind=dp )
  
  ! construct the local Hamiltonian
  
  ! constant terms
  hamq_local = kin0
  
  ! linear terms
  call local_cyclic_dims( nr, nc )
  forall( i=1:nr, j=1:nc ) &
    hamq_local(i,j) = hamq_local(i,j) + two * sum( kin1(i,j,1:3) * zqvec(1:3) )
  
  ! quadratic terms
  if( ndiag_cyclic > 0 ) then
    forall( i=1:ndiag_cyclic ) &
      hamq_local(diag_cyclic(1,i),diag_cyclic(2,i)) &
        = hamq_local(diag_cyclic(1,i),diag_cyclic(2,i)) + zqvec2
  endif

  return
  end subroutine build_hamq_kin

! ---------------------------------------------------------------------- 
  subroutine build_hamq_vloc( qvec, lkin, ispin, hamq_local )
! ---------------------------------------------------------------------- 

  real(dp),intent(in) :: qvec(3)
  logical,intent(in) :: lkin  ! true (false) if qvec is in kinetic (nonlocal) basis
  integer,intent(in) :: ispin
  complex(dp),intent(out) :: hamq_local(:,:)

  ! constant terms
  hamq_local = vloc(:,:,ispin)
  
  return
  end subroutine build_hamq_vloc

! ---------------------------------------------------------------------- 
  subroutine build_hamq_local( qvec, lkin, ispin, hamq_local )
! ---------------------------------------------------------------------- 

  use hamq_pool, only : local_cyclic_dims, ndiag_cyclic, diag_cyclic

  real(dp),intent(in) :: qvec(3)
  logical,intent(in) :: lkin  ! true (false) if qvec is in kinetic (nonlocal) basis
  integer,intent(in) :: ispin
  complex(dp),intent(out) :: hamq_local(:,:)

  integer :: i, j, ierr
  integer :: nr, nc
  real(dp) :: tqvec(3), sqvec(3)
  complex(dp) :: zqvec(3), zqvec2

  if( lkin ) then
    tqvec = qvec
    sqvec = matmul( trkin2nlp, qvec )
  else
    tqvec = matmul( trnlp2kin, qvec )
    sqvec = qvec
  endif
  
  ! reduce sqvec to [0,1)
  sqvec = sqvec - floor( sqvec )
  ! redefine consistent cartesian coordinates
  tqvec = matmul( trnlp2kin, sqvec )

  zqvec = cmplx( tqvec, kind=dp )
  zqvec2 = cmplx( dot_product(tqvec,tqvec), kind=dp )
  
  ! construct the local Hamiltonian
  
  ! constant terms
  hamq_local = vloc(:,:,ispin) + kin0
  
  ! linear terms
  call local_cyclic_dims( nr, nc )
  forall( i=1:nr, j=1:nc ) &
    hamq_local(i,j) = hamq_local(i,j) + two * sum( kin1(i,j,1:3) * zqvec(1:3) )
  
  ! quadratic terms
  if( ndiag_cyclic > 0 ) then
    forall( i=1:ndiag_cyclic ) &
      hamq_local(diag_cyclic(1,i),diag_cyclic(2,i)) &
        = hamq_local(diag_cyclic(1,i),diag_cyclic(2,i)) + zqvec2
  endif

  return
  end subroutine build_hamq_local

! ---------------------------------------------------------------------- 
  subroutine build_dhamq_local( ideriv, qvec, lkin, dhamq_local )
! ---------------------------------------------------------------------- 

  use hamq_pool, only : local_cyclic_dims, ndiag_cyclic, diag_cyclic

  ! generate the derivative in the ideriv dimensionless reciprocal lattice direction

  integer,intent(in) :: ideriv
  real(dp),intent(in) :: qvec(3)
  logical,intent(in) :: lkin  ! true (false) if qvec is in kinetic (nonlocal) basis
  complex(dp),intent(out) :: dhamq_local(:,:)

  integer :: i, j, k, ierr
  integer :: nr, nc
  real(dp) :: tqvec(3), sqvec(3)
  complex(dp) :: zqvec(3)
  complex(dp) :: ztrnlp2kin(3,3)

  if( ideriv > 3 .or. ideriv < 1 ) &
    call errore( 'build_dhamq_local','ideriv must be 1,2, or 3',abs(ideriv))
 
  if( lkin ) then
    tqvec = qvec
    sqvec = matmul( trkin2nlp, qvec )
  else
    tqvec = matmul( trnlp2kin, qvec )
    sqvec = qvec
  endif
  
  ! reduce sqvec to [0,1)
  sqvec = sqvec - floor( sqvec )
  ! redefine consistent cartesian coordinates
  tqvec = matmul( trnlp2kin, sqvec )

  ! complex versions
  zqvec = cmplx( tqvec, kind=dp ) 
  ztrnlp2kin = cmplx( trnlp2kin, kind=dp )

  ! construct the derivative of the local Hamiltonian
  ! from linear terms and quadratic

  ! reciprocal lattice coords
  call local_cyclic_dims( nr, nc )
  forall(i=1:nr, j=1:nc) &
    dhamq_local(i,j) = two * sum( kin1(i,j,1:3)*ztrnlp2kin(1:3,ideriv) )
  if( ndiag_cyclic > 0 ) then
    forall( i=1:ndiag_cyclic ) &
      dhamq_local(diag_cyclic(1,i),diag_cyclic(2,i)) &
        = dhamq_local(diag_cyclic(1,i),diag_cyclic(2,i)) &
        + two * sum( zqvec(1:3)*ztrnlp2kin(1:3,ideriv) )
  endif
  
  return
  end subroutine build_dhamq_local

! ---------------------------------------------------------------------- 
  subroutine build_hamq_projs( qvec, lkin, betaq )
! ---------------------------------------------------------------------- 

  use mp_world, only : mpime

  real(dp),intent(in) :: qvec(3)
  logical,intent(in) :: lkin  ! true (false) if qvec is in kinetic (nonlocal) basis
  complex(dp),intent(out) :: betaq(:,:)

  integer :: i, j, ierr
  real(dp) :: sqvec(3)
  real(dp) :: fR, fI

  if( nproj > 0 ) then

  if( lkin ) then
    sqvec = matmul( trkin2nlp, qvec )
  else
    sqvec = qvec
  endif
  
  ! reduce sqvec to [0,1)
  sqvec = sqvec - floor( sqvec )
  
  ! fit nonlocal potential terms next
  ! uses qvec in crystal coordinates
  do j=1,nbasis
  do i=1,nproj
    fR = spline_evaluate( sqvec, bscoefR(1:size(bscoefR,1),i,j) )
    fI = spline_evaluate( sqvec, bscoefI(1:size(bscoefI,1),i,j) )
    betaq(i,j) = cmplx( fR, fI )
  enddo
  enddo

  endif ! nproj > 0

  return
  end subroutine build_hamq_projs

! ---------------------------------------------------------------------- 
  subroutine build_hamq_atomprojs( qvec, lkin, atomq )
! ---------------------------------------------------------------------- 

  use mp_world, only : mpime

  real(dp),intent(in) :: qvec(3)
  logical,intent(in) :: lkin  ! true (false) if qvec is in kinetic (nonlocal) basis
  complex(dp),intent(out) :: atomq(:,:)

  integer :: i, j, ierr
  real(dp) :: sqvec(3)
  real(dp) :: fR, fI

  if( natomproj > 0 ) then

  if( lkin ) then
    sqvec = matmul( trkin2nlp, qvec )
  else
    sqvec = qvec
  endif
  
  ! reduce sqvec to [0,1)
  sqvec = sqvec - floor( sqvec )
  
  ! fit nonlocal potential terms next
  ! uses qvec in crystal coordinates
  do j=1,nbasis
  do i=1,natomproj
    fR = spline_evaluate( sqvec, ascoefR(1:size(ascoefR,1),i,j) )
    fI = spline_evaluate( sqvec, ascoefI(1:size(ascoefI,1),i,j) )
    atomq(i,j) = cmplx( fR, fI )
  enddo
  enddo

  endif ! natomproj > 0

  return
  end subroutine build_hamq_atomprojs

! ---------------------------------------------------------------------- 
  subroutine build_dhamq_projs( ideriv, qvec, lkin, dbetaq )
! ---------------------------------------------------------------------- 

  ! generate the derivative in the ideriv dimensionless reciprocal lattice direction

  integer,intent(in) :: ideriv
  real(dp),intent(in) :: qvec(3)
  logical,intent(in) :: lkin  ! true (false) if qvec is in kinetic (nonlocal) basis
  complex(dp),intent(out) :: dbetaq(nproj,nbasis)

  integer :: i, j, ierr
  real(dp) :: sqvec(3)
  real(dp) :: fR, fI

  if( nproj > 0 ) then

  if( lkin ) then
    sqvec = matmul( trkin2nlp, qvec )
  else
    sqvec = qvec
  endif
  
  ! reduce sqvec to [0,1)
  sqvec = sqvec - floor( sqvec )
  
  ! fit nonlocal potential terms next
  ! uses qvec in crystal coordinates
  do j=1,nbasis
  do i=1,nproj
    fR = spline_derivative( ideriv, sqvec, bscoefR(1:size(bscoefR,1),i,j) )
    fI = spline_derivative( ideriv, sqvec, bscoefI(1:size(bscoefI,1),i,j) )
    dbetaq(i,j) = cmplx( fR, fI )
  enddo
  enddo

  endif ! nproj > 0

  return
  end subroutine build_dhamq_projs

! ---------------------------------------------------------------------- 
  subroutine build_hamq_nlprojs( qvec, lkin, betaq )
! ---------------------------------------------------------------------- 

  real(dp),intent(in) :: qvec(3)
  logical,intent(in) :: lkin  ! true (false) if qvec is in kinetic (nonlocal) basis
  complex(dp),intent(out) :: betaq(nproj_nl,nbasis)

  integer :: i, j, ierr
  real(dp) :: sqvec(3)
  real(dp) :: fR, fI

  if( nproj_nl > 0 ) then

  if( lkin ) then
    sqvec = matmul( trkin2nlp, qvec )
  else
    sqvec = qvec
  endif
  
  ! reduce sqvec to [0,1)
  sqvec = sqvec - floor( sqvec )
  
  ! fit nonlocal potential terms next
  ! uses qvec in crystal coordinates
  do j=1,nbasis
  do i=1,nproj_nl
    fR = spline_evaluate( sqvec, bscoefR(1:size(bscoefR,1),index_nlproj_betaq(i),j) )
    fI = spline_evaluate( sqvec, bscoefI(1:size(bscoefI,1),index_nlproj_betaq(i),j) )
    betaq(i,j) = cmplx( fR, fI )
  enddo
  enddo

  endif ! nproj_nl > 0

  return
  end subroutine build_hamq_nlprojs

! ---------------------------------------------------------------------- 
  subroutine mom_hamq( qvec, lkin, mom )
! ---------------------------------------------------------------------- 

  use hamq_pool, only : ndiag_cyclic, diag_cyclic

  real(dp),intent(in) :: qvec(3)
  logical,intent(in) :: lkin  ! true (false) if qvec is in kinetic (nonlocal) basis
  complex(dp),intent(out) :: mom(:,:,:)

  integer :: i, j, ixyz
  real(dp) :: tqvec(3), sqvec(3)
  complex(dp) :: zqvec(3)


  if( lkin ) then
    tqvec = qvec
    sqvec = matmul( trkin2nlp, qvec )
  else
    sqvec = qvec
  endif
  
  ! reduce sqvec to [0,1)
  sqvec = sqvec - floor( sqvec )
  ! redefine consistent cartesian coordinates
  tqvec = matmul( trnlp2kin, sqvec )

  zqvec = cmplx( tqvec, kind=dp )
  
  ! linear term in kinetic energy plus the q-vector
  mom = kin1 
  if( ndiag_cyclic > 0 ) then
    forall( i=1:ndiag_cyclic, ixyz=1:3 ) &
      mom(diag_cyclic(1,i),diag_cyclic(2,i),ixyz) &
        = mom(diag_cyclic(1,i),diag_cyclic(2,i),ixyz) + zqvec(ixyz)
  endif
  
  return
  end subroutine mom_hamq


! ---------------------------------------------------------------------- 
  subroutine vel_hamq( qvec, lkin, momcr, vel )
! ---------------------------------------------------------------------- 

  real(dp),intent(in) :: qvec(3)
  logical,intent(in) :: lkin  ! true (false) if qvec is in kinetic (nonlocal) basis
  complex(dp),intent(in) :: momcr(:,:,:)
  complex(dp),intent(out) :: vel(:,:,:)

  integer :: i, j, ixyz, ierr
  real(dp) :: tqvec(3), sqvec(3)
  complex(dp) :: zqvec(3)
  real(dp) :: fR, fI
  complex(dp),allocatable :: betaq(:,:)
  complex(dp),allocatable :: ztmp(:,:)


  ! momentum first
  call mom_hamq( qvec, lkin, vel )

  if( lkin ) then
    tqvec = qvec
    sqvec = matmul( trkin2nlp, qvec )
  else
    sqvec = qvec
  endif
  
  ! reduce sqvec to [0,1)
  sqvec = sqvec - floor( sqvec )
  ! redefine consistent cartesian coordinates
  tqvec = matmul( trnlp2kin, sqvec )

  zqvec = cmplx( tqvec, kind=dp )
  
  ! now the projector contribution
  if( nproj > 0 ) then

    ! fit nonlocal potential terms next
    ! uses qvec in crystal coordinates
    ! requires that the momentum core-repair matrix elements
    ! momcr are input

    allocate( betaq(nproj,nbasis), ztmp(nproj,nbasis), stat=ierr )
    if( ierr /= 0 ) call errore('vel_hamq','unable to allocate tmp arrays',1)
     
    do j=1,nbasis
    do i=1,nproj
      fR = spline_evaluate( sqvec, bscoefR(1:size(bscoefR,1),i,j) )
      fI = spline_evaluate( sqvec, bscoefI(1:size(bscoefI,1),i,j) )
      betaq(i,j) = cmplx( fR, fI )
    enddo
    enddo

    ! fold in the core-repair terms
    do ixyz=1,3
      CALL ZGEMM( 'N', 'N', nproj, nbasis, nproj, one, &
                  momcr, nproj, betaq, nproj, zero, ztmp, nproj )
      CALL ZGEMM( 'C', 'N', nbasis, nbasis, nproj, one, &
                  betaq, nproj, ztmp, nproj, one, vel(1,1,ixyz), nbasis )
    enddo
     
    deallocate( betaq, ztmp )
  
  endif ! projectors
  
  return
  end subroutine vel_hamq


! ---------------------------------------------------------------------- 
  subroutine open_hamq( prefix )
! ---------------------------------------------------------------------- 

  use mp_world, only : mpime, root, world_comm
  use mp, only : mp_bcast
  use mpio

  character(*),intent(in) :: prefix
  character(len(prefix)+7) :: hamloc_file, hamprj_file, haminf_file
  integer,external :: freeunit

  ! MPI-IO
  if( mpime==root ) then
    hamloc_file=trim(prefix)//'.hamloc'
    hamprj_file=trim(prefix)//'.hamprj'
  endif

  call mp_bcast( hamloc_file, root, world_comm )
  call mp_bcast( hamprj_file, root, world_comm )

  call mp_file_open_dp( hamloc_file, fhhamloc, root, world_comm )
  call mp_file_open_dp( hamprj_file, fhhamprj, root, world_comm )

  end subroutine open_hamq


! ---------------------------------------------------------------------- 
  subroutine close_hamq( )
! ---------------------------------------------------------------------- 

  use mp_world, only : mpime, root

  integer :: ierr

  call mpi_file_close( fhhamloc, ierr )
  call mpi_file_close( fhhamprj, ierr )

  end subroutine close_hamq


! ---------------------------------------------------------------------- 
  subroutine write_hamprj( posn )
! ---------------------------------------------------------------------- 

  use mp_world, only : mpime, root, world_comm
  use mpio, only : mp_file_scatter_write
  use parallel_include

  integer(kind=MPI_OFFSET_KIND),optional :: posn
  integer(kind=MPI_OFFSET_KIND),save :: fposn=0

  if( present(posn) ) then
    fposn=posn
  endif
  call mp_file_scatter_write( fhhamprj, fposn, bscoefR, root, world_comm )
  call mp_file_scatter_write( fhhamprj, fposn, bscoefI, root, world_comm )

  end subroutine write_hamprj


! ---------------------------------------------------------------------- 
  subroutine write_hamloc( fposn, n, ham )
! ---------------------------------------------------------------------- 

  use mp_world, only : root, world_comm
  use mp_scatt, only : mp_scatter_size, mp_scatter
  use mpio, only : mp_file_scatter_write
  use parallel_include

  integer(kind=MPI_OFFSET_KIND),intent(inout) :: fposn
  integer,intent(in) :: n
  complex(dp),intent(in) :: ham(:)

  call mp_file_scatter_write( fhhamloc, fposn, ham, root, world_comm )

  end subroutine write_hamloc


! ---------------------------------------------------------------------- 
  subroutine read_hamprj( fposn, proj )
! ---------------------------------------------------------------------- 

  use kinds, only : dp
  use mp_world, only : mpime, root
  use mp, only : mp_bcast, mp_gather
  use mp_scatt, only : mp_scatter_size, mp_scatter_displ
  use mpio, only : mp_file_scatter_read
  use hamq_pool, only : nproc_per_pool, &
                        nproj_global, &
                        natom_global, type_atom_global, &
                        mypoolid, mypoolroot, intra_pool_comm, &
                        mypool, rootpool, cross_pool_comm, &
                        index_betaq_global
  use parallel_include

  integer(kind=MPI_OFFSET_KIND),intent(inout) :: fposn
  ! proj dimensions are (coefs, projs, basis)
  real(dp),intent(out) :: proj(:,:,:)

  real(dp),allocatable :: projsct(:,:,:), projtmp(:,:,:)
  integer :: ierr, i, j, k
  integer :: iatom, iatom_global
  integer :: it, ip, ip_global, np, np_global, nproj_type_max
  integer :: length, displ, nk, npg
  integer :: id
  integer,allocatable :: length_id(:), displ_id(:)

  if( mpime==root ) write(*,*) ' read_hamprj '

  call mp_scatter_size( size(proj,3), length, mypoolroot, intra_pool_comm )
  call mp_scatter_displ( size(proj,3), displ, mypoolroot, intra_pool_comm )

  allocate( length_id(nproc_per_pool), displ_id(nproc_per_pool) )

  call mp_gather( length, length_id, mypoolroot, intra_pool_comm )
  call mp_gather( displ,  displ_id,  mypoolroot, intra_pool_comm )
  call mp_bcast( length_id, mypoolroot, intra_pool_comm )
  call mp_bcast( displ_id, mypoolroot, intra_pool_comm )

  nk=size(proj,1)
  npg = nproj_global
  allocate( projsct(nk,npg,length) )
#ifdef DEBUG
  write(stdout,*) ' before scattering: proj ', size(proj,1), size(proj,2), size(proj,3)
  write(stdout,*) '  after scattering: proj ', nk, npg, length, ' mpime ', mpime
#endif

  ! only do this on first pool
  if( mypool == rootpool ) then
    call mp_file_scatter_read( fhhamprj, fposn, projsct, &
                               mypoolroot, intra_pool_comm )
  endif
  ! send across pools - within slices
  call mp_bcast( projsct, rootpool, cross_pool_comm )

  ! redistribute within pool
  proj=zero
  do id=0,nproc_per_pool-1
    allocate( projtmp(nk,npg,length_id(id+1)) )
    if( mypoolid==id ) projtmp=projsct
    call mp_bcast( projtmp, id, intra_pool_comm )

    ip=0
    ip_global=0
    iatom=0
    do iatom_global=1,natom_global
      it=type_atom_global(iatom_global)
      np_global = nproj_type(it)

      if( mod(iatom_global-1,nproc_per_pool) == mypoolid ) then
        iatom=iatom+1
        it=type_atom(iatom)
        np = nproj_type(it)
        ! check
        if( np /= np_global ) call errore('read_hamprj','projector mismatch',1)
        
#ifdef DEBUG
        write(stdout,*) ' read_hamprj:'
        write(stdout,*) ' iatom_global, natom_global : ', iatom_global, natom_global
        write(stdout,*) ' iatom, it, np : ', iatom, it, np
        write(stdout,*) ' nk, np, length_id : ', nk, np, length_id(id+1)
        write(stdout,*) ' proj    : ', size(proj,1), size(proj,2), size(proj,3)
        write(stdout,*) ' projtmp : ', size(projtmp,1), size(projtmp,2), size(projtmp,3)
        write(stdout,*) index_betaq(1:np,iatom)
        write(stdout,*) index_betaq_global(1:np,iatom_global)
#endif

        forall( i=1:nk, j=1:np, k=1:length_id(id+1) ) &
          proj(i,index_betaq(j,iatom),displ_id(id+1)+k) &
          = projtmp(i,index_betaq_global(j,iatom_global),k)
        ip = ip + np
      endif

      ip_global = ip_global + np_global
    enddo
    if( iatom /= natom ) call errore('read_hamprj','problem with atoms',1)
    deallocate( projtmp )
  enddo

  deallocate( projsct )

  end subroutine read_hamprj



! ---------------------------------------------------------------------- 
  subroutine read_hamprj_ldaU( fposn, proj )
! ---------------------------------------------------------------------- 

  use kinds, only : dp
  use mp_world, only : mpime, root
  use mp, only : mp_bcast, mp_gather
  use mp_scatt, only : mp_scatter_size, mp_scatter_displ
  use mpio, only : mp_file_scatter_read
  use hamq_pool, only : nproc_per_pool, &
                        natomproj_global, &
                        natom_global, type_atom_global, &
                        mypoolid, mypoolroot, intra_pool_comm, &
                        mypool, rootpool, cross_pool_comm, &
                        index_betaq_global, index_ldaUq_global
  use parallel_include

  integer(kind=MPI_OFFSET_KIND),intent(inout) :: fposn
  ! proj dimensions are (coefs, projs, basis)
  real(dp),intent(out) :: proj(:,:,:)

  real(dp),allocatable :: projsct(:,:,:), projtmp(:,:,:)
  integer :: ierr, i, j, k
  integer :: iatom, iatom_global
  integer :: it, ip, ip_global, np, np_global, nproj_type_max
  integer :: length, displ, nk, npg
  integer :: id
  integer,allocatable :: length_id(:), displ_id(:)

  if( mpime==root ) write(*,*) ' read_hamprj_ldaU '

  call mp_scatter_size( size(proj,3), length, mypoolroot, intra_pool_comm )
  call mp_scatter_displ( size(proj,3), displ, mypoolroot, intra_pool_comm )

  allocate( length_id(nproc_per_pool), displ_id(nproc_per_pool) )

  call mp_gather( length, length_id, mypoolroot, intra_pool_comm )
  call mp_gather( displ,  displ_id,  mypoolroot, intra_pool_comm )
  call mp_bcast( length_id, mypoolroot, intra_pool_comm )
  call mp_bcast( displ_id, mypoolroot, intra_pool_comm )

  nk=size(proj,1)
  npg=natomproj_global
  allocate( projsct(nk,npg,length) )

  ! only do this on first pool
  if( mypool == rootpool ) then
    call mp_file_scatter_read( fhhamprj, fposn, projsct, &
                               mypoolroot, intra_pool_comm )
  endif
  ! send across pools - within slices
  call mp_bcast( projsct, rootpool, cross_pool_comm )

  ! redistribute within pool
  proj=zero
  do id=0,nproc_per_pool-1
    allocate( projtmp(nk,npg,length_id(id+1)) )
    if( mypoolid==id ) projtmp=projsct
    call mp_bcast( projtmp, id, intra_pool_comm )

    ip=0
    ip_global=0
    iatom=0
    do iatom_global=1,natom_global
      it=type_atom_global(iatom_global)
      if( Hubbard_U(it) /= 0.d0 .or. Hubbard_alpha(it) /= 0.d0 ) then
        np_global = 2*Hubbard_l(it)+1
      else
        np_global = 0
      endif

      if( mod(iatom_global-1,nproc_per_pool) == mypoolid ) then
        iatom=iatom+1
        it=type_atom(iatom)
        if( Hubbard_U(it) /= 0.d0 .or. Hubbard_alpha(it) /= 0.d0 ) then
          np = 2*Hubbard_l(it)+1
        else
          np = 0
        endif
        ! check
        if( np /= np_global ) call errore('read_hamprj','projector mismatch',1)
        
        forall( i=1:nk, j=1:np, k=1:length_id(id+1) ) &
          proj(i,index_ldaUq(j,iatom),displ_id(id+1)+k) &
          = projtmp(i,index_ldaUq_global(j,iatom_global),k)
        ip = ip + np
      endif

      ip_global = ip_global + np_global
    enddo
    if( iatom /= natom ) call errore('read_hamprj','problem with atoms',1)
    deallocate( projtmp )
  enddo

  deallocate( projsct )

  end subroutine read_hamprj_ldaU


! ---------------------------------------------------------------------- 
  subroutine read_hamloc( fposn, n, hmat )
! ---------------------------------------------------------------------- 

  use mp_world, only : root, nproc
  use mp_scatt, only : mp_scatter_size, mp_scatter_displ, mp_scatter
  use mp, only : mp_max, mp_sum, mp_bcast, mp_gather
  use mpio, only : mp_file_scatter_read
  use hamq_pool, only : nproc_per_pool, &
                        cyclic_localindex, &
                        mypool, rootpool, cross_pool_comm, &
                        mypoolid, mypoolroot, intra_pool_comm
  use parallel_include

  integer(kind=MPI_OFFSET_KIND),intent(inout) :: fposn
  integer,intent(in) :: n 
  complex(dp),intent(out) :: hmat(:,:)

  integer :: ij, i, j, il, jl, ijl
  logical :: islocal
  integer :: length, displ, id
  integer,allocatable :: length_id(:), displ_id(:)
  integer :: ierr
  complex(dp),allocatable :: ham(:), hamtmp(:)


  write(stdout,*) ' read_hamloc '
#ifdef DEBUG
  write(stdout,*) ' read_hamloc: ', fposn, n, size(hmat,1), size(hmat,2)
#endif

  call mp_scatter_size( (n*(n+1))/2, length, mypoolroot, intra_pool_comm )
  call mp_scatter_displ( (n*(n+1))/2, displ, mypoolroot, intra_pool_comm )

  allocate( length_id(nproc_per_pool), displ_id(nproc_per_pool) )

  call mp_gather( length, length_id, mypoolroot, intra_pool_comm )
  call mp_gather( displ,  displ_id,  mypoolroot, intra_pool_comm )
  call mp_bcast( length_id, mypoolroot, intra_pool_comm )
  call mp_bcast( displ_id, mypoolroot, intra_pool_comm )

#ifdef DEBUG
  write(stdout,*) ' read_hamloc: ', length, displ
#endif

  allocate( ham(length), stat=ierr )
  if( ierr/=0 ) call errore('read_hamloc','problem allocating space',abs(ierr))

  ! only do this on first pool
  if( mypool == rootpool ) then
    call mp_file_scatter_read( fhhamloc, fposn, ham, &
                               mypoolroot, intra_pool_comm )
  endif
  ! send across pools - within slices
  call mp_bcast( ham, rootpool, cross_pool_comm )

  ! need to unpack the ham to hmat within pool
  hmat=zero
  do id=0,nproc_per_pool-1
    allocate( hamtmp(length_id(id+1)) )
    if( mypoolid==id ) hamtmp=ham
    call mp_bcast( hamtmp, id, intra_pool_comm )
    ij=0
    do j=1,n
      do i=1,j
        ij=ij+1
        if( ij <= displ_id(id+1) .or. ij > displ_id(id+1)+length_id(id+1) ) cycle

        ijl=ij-displ_id(id+1)
        call cyclic_localindex( i, j, il, jl, islocal )
        if( islocal ) hmat(il,jl) = hamtmp(ijl)
      enddo
    enddo
    deallocate( hamtmp )
  enddo

  deallocate( ham )

  end subroutine read_hamloc


! ---------------------------------------------------------------------- 
  subroutine write_haminf()
! ---------------------------------------------------------------------- 

  namelist /system/ nelec, alat, omega, at, bg, tpiba, ncpp, &
                    trkin2nlp, trnlp2kin
  namelist /hamdim/ nbasis, nproj, nproj_nl
  namelist /nonloc/ natom, ntype, nproj_type_max

  
  ! system details
  write(fhhaminf,nml=system) 

  ! hamiltonian dimensions
  write(fhhaminf,nml=hamdim)

  if( nproj > 0 ) then

    ! atomic projectors
    write(fhhaminf,nml=nonloc) 

    ! atomic projector details
    write(fhhaminf,*) type_atom
    write(fhhaminf,*) nproj_type
    write(fhhaminf,*) nproj_type_nl

    write(fhhaminf,*) index_nlproj_type
    write(fhhaminf,*) index_betaq

    write(fhhaminf,*) index_nlproj_betaq

    ! spline info
    write(fhhaminf,*) kxord, kyord, kzord
    write(fhhaminf,*) nxcoef, nycoef, nzcoef
    write(fhhaminf,*) xknot, yknot, zknot

  endif

  return
  end subroutine write_haminf


! ---------------------------------------------------------------------- 
  subroutine write_hamq( iunit )
! ---------------------------------------------------------------------- 

  integer,intent(in) :: iunit

  integer :: i
  
  ! system details
  write(iunit) nelec, alat, omega, at, bg, tpiba, ncpp, nspin, lda_plus_u
  write(stdout,*) ' write_hamq: nspin = ', nspin

  ! transformation matrices
  write(iunit) trkin2nlp
  write(iunit) trnlp2kin

  ! size of basis
  write(iunit) nbasis

  ! non-local Hamiltonian
  write(iunit) nproj, nproj_nl

  if( nproj > 0 ) then

    ! atomic projectors
    write(iunit) natom, ntype, nproj_type_max

    ! read atomic projector details
    write(iunit) type_atom
    write(iunit) nproj_type
    write(iunit) nproj_type_nl

    write(iunit) index_nlproj_type
    write(iunit) index_betaq

    write(iunit) index_nlproj_betaq

    ! spline info
    write(iunit) kxord, kyord, kzord
    write(iunit) nxcoef, nycoef, nzcoef
    write(iunit) xknot, yknot, zknot

    ! LDA+U
    if( lda_plus_u ) then
      write(iunit) natomproj
      write(iunit) Hubbard_lmax
      write(iunit) Hubbard_l
      write(iunit) Hubbard_U
      write(iunit) Hubbard_alpha
      write(iunit) index_ldaUq
    endif
  endif

  return
  end subroutine write_hamq


! ---------------------------------------------------------------------- 
  subroutine read_hamq( iunit )
! ---------------------------------------------------------------------- 

  use mp_world, only : mpime
  use hamq_pool, only : nproc_per_pool, natom_global, &
                        nbasis_stripe, nbasis_stripe_max, nb_stripe, &
                        create_cyclic_desc, local_cyclic_dims, &
                        distrib_projs_byatom, distrib_ldaU_byatom, &
                        create_striped_desc, &
                        local_striped_dim
  use parallel_include


  integer,intent(in) :: iunit

  integer(kind=MPI_OFFSET_KIND) :: fposn
  integer :: ixyz, nr, nc, nstripe, ispin
  complex(dp),allocatable :: tmp(:,:)

  call read_hamdim( iunit )

  ! distribute hamiltonian information within pool
  call create_cyclic_desc( nbasis )
  call local_cyclic_dims( nr, nc )

  ! allocate space for local Hamiltonian
  call alloc_hamq_local( nr, nc )

  ! tmp space
  allocate( tmp(nr,nc) )

  ! read local Hamiltonian
  write(stdout,*) ' read_hamloc '
  fposn=0
  call read_hamloc( fposn, nbasis, kin0 )
  do ixyz=1,3
    ! just using vloc as tmp space
    call read_hamloc( fposn, nbasis, tmp )
    kin1(:,:,ixyz)=tmp
  enddo
  do ispin=1,nspin
    call read_hamloc( fposn, nbasis, tmp )
    vloc(:,:,ispin) = tmp
  enddo
  write(stdout,*) ' done read_hamloc '

  ! distribute projector information within pool
  call distrib_projs_byatom( natom, ntype, type_atom, nproj_type, &
                             nproj_type_nl, nproj_type_max, nproj, &
                             nproj_nl, index_betaq, index_nlproj_betaq )

  if( lda_plus_u ) then
    call distrib_ldaU_byatom( ntype, &
                              natomproj, Hubbard_lmax, Hubbard_l, &
                              Hubbard_U, Hubbard_alpha, index_ldaUq )
  endif

  ! allocate space for projectors
  call alloc_projs( nbasis )

  ! read projectors
  fposn=0
  call read_hamprj( fposn, bscoefR )
  call read_hamprj( fposn, bscoefI )

  if( lda_plus_u ) then
    call alloc_projs_ldaU( nbasis )
    call read_hamprj_ldaU( fposn, ascoefR )
    call read_hamprj_ldaU( fposn, ascoefI )
  endif

  ! set up striped distribution for calculating non-local potential
  call create_striped_desc( nbasis )

  end subroutine read_hamq

! ---------------------------------------------------------------------- 
  subroutine read_hamdim( iunit )
! ---------------------------------------------------------------------- 

  use mp_world, only : mpime, root

  integer,intent(in) :: iunit

  integer :: ierr
  integer :: i

  if( mpime==root ) then

  ! system details
  read(iunit) nelec, alat, omega, at, bg, tpiba, ncpp, nspin, lda_plus_u
  write(stdout,*) ' nspin = ', nspin

  ! transformation matrices
  read(iunit) trkin2nlp
  read(iunit) trnlp2kin

  ! size of basis
  read(iunit) nbasis

  ! report
  write(stdout,*) ' basis set size = ', nbasis
            
  ! non-local Hamiltonian
  read(iunit) nproj, nproj_nl

  write(stdout,*) ' non-local proj = ', nproj
  write(stdout,*) '    nl ham proj = ', nproj_nl

  if( nproj > 0 ) then

    ! atomic projectors
    read(iunit) natom, ntype, nproj_type_max

    ! allocate space
    call alloc_atomic_proj( nproj, nproj_nl, natom, ntype, nproj_type_max )

    ! read atomic projector details
    read(iunit) type_atom
    read(iunit) nproj_type
    read(iunit) nproj_type_nl

    read(iunit) index_nlproj_type
    read(iunit) index_betaq

    read(iunit) index_nlproj_betaq

    ! nonlocal spline info for projectors
    read(iunit) kxord, kyord, kzord
    read(iunit) nxcoef, nycoef, nzcoef

    ! report
    write(stdout,*) '    spline orders = ', kxord, kyord, kzord
    write(stdout,*) ' number of coeffs = ', nxcoef, nycoef, nzcoef
            
    ! allocate space
    call alloc_knots()

    ! knots
    read(iunit) xknot, yknot, zknot

    ! LDA+U
    if( lda_plus_u ) then
      read(iunit) natomproj
      read(iunit) Hubbard_lmax
      allocate( Hubbard_l(ntype), Hubbard_U(ntype), Hubbard_alpha(ntype), index_ldaUq(2*Hubbard_lmax+1,natom) )
      read(iunit) Hubbard_l
      read(iunit) Hubbard_U
      read(iunit) Hubbard_alpha
      read(iunit) index_ldaUq

      write(stdout,*) ' natomproj = ', natomproj
      write(stdout,*) ' Hubbard parameters: l         U     alpha'
      do i=1,ntype
        write(stdout,'(a,i4,3x,i4,2f10.6)') '  atom type ', i, Hubbard_l(i), &
                     Hubbard_U(i)*rytoev, Hubbard_alpha(i)*rytoev
      enddo
      do i=1,natom
        write(stdout,*) ' atom ', i, ' index: ', index_ldaUq(:,i)
      enddo
      write(stdout,*)
    endif
  endif ! nproj > 0

  endif ! mpime==root

  call bcast_hamdim()

  return
  end subroutine read_hamdim


! ---------------------------------------------------------------------- 
  subroutine bcast_hamdim()
! ---------------------------------------------------------------------- 

  ! this uses Cavazzoni's mp routines
  use mp, only : mp_bcast
  use mp_world, only : mpime, root, world_comm

  integer :: nbasis_, nproj_, nproj_nl_
  integer :: i

  ! system details
  call mp_bcast( nelec, root, world_comm )
  call mp_bcast( alat, root, world_comm )
  call mp_bcast( omega, root, world_comm )
  call mp_bcast( at, root, world_comm )
  call mp_bcast( bg, root, world_comm )
  call mp_bcast( tpiba, root, world_comm )
  call mp_bcast( ncpp, root, world_comm )
  call mp_bcast( nspin, root, world_comm )
  call mp_bcast( lda_plus_u, root, world_comm )

  ! transformation matrices
  call mp_bcast( trkin2nlp, root, world_comm )
  call mp_bcast( trnlp2kin, root, world_comm )

  ! size of basis
  call mp_bcast( nbasis, root, world_comm )

  call mp_bcast( nproj, root, world_comm )
  call mp_bcast( nproj_nl, root, world_comm )

  if( nproj > 0 ) then

    ! atomic projectors
    call mp_bcast(natom, root, world_comm )
    call mp_bcast(ntype, root, world_comm )
    call mp_bcast(nproj_type_max, root, world_comm )

    if( mpime/=root ) then
      ! allocate space
      nproj_ = nproj
      nproj_nl_ = nproj_nl
      call alloc_atomic_proj( nproj_, nproj_nl_, natom, ntype, nproj_type_max )
    endif

    call mp_bcast(type_atom, root, world_comm )
    call mp_bcast(nproj_type, root, world_comm )
    call mp_bcast(nproj_type_nl, root, world_comm )

    call mp_bcast(index_nlproj_type, root, world_comm )
    call mp_bcast(index_betaq, root, world_comm )
    call mp_bcast( index_nlproj_betaq, root, world_comm )

    ! nonlocal spline info
    call mp_bcast( kxord, root, world_comm )
    call mp_bcast( kyord, root, world_comm )
    call mp_bcast( kzord, root, world_comm )
    call mp_bcast( nxcoef, root, world_comm )
    call mp_bcast( nycoef, root, world_comm )
    call mp_bcast( nzcoef, root, world_comm )

    ! allocate space
    if( mpime /= root ) call alloc_knots()

    ! knots
    call mp_bcast( xknot, root, world_comm )
    call mp_bcast( yknot, root, world_comm )
    call mp_bcast( zknot, root, world_comm )

    ! LDA+U
    if( lda_plus_u ) then
      call mp_bcast( natomproj, root, world_comm )
      call mp_bcast( Hubbard_lmax, root, world_comm )
      if( mpime/=root ) allocate( Hubbard_l(ntype), Hubbard_U(ntype), Hubbard_alpha(ntype), index_ldaUq(2*Hubbard_lmax+1,natom) )
      call mp_bcast( Hubbard_l, root, world_comm )
      call mp_bcast( Hubbard_U, root, world_comm )
      call mp_bcast( Hubbard_alpha, root, world_comm )
      call mp_bcast( index_ldaUq, root, world_comm )
    endif
  endif

  return
  end subroutine bcast_hamdim


! ---------------------------------------------------------------------- 
  subroutine dump_hamq( iunit, n, m, lda, a )
! ---------------------------------------------------------------------- 

  integer,intent(in) :: iunit, n, m, lda
  complex(dp),intent(in) :: a(lda,m)

  integer :: i, j

  do j=1,m
    do i=1,n
      write(iunit,'(2i6,2e12.5)') i, j, a(i,j)
    enddo
    write(iunit,*)
  enddo

  return
  end subroutine dump_hamq

! ---------------------------------------------------------------------- 
  subroutine dump_hamq_packed( iunit, n, a )
! ---------------------------------------------------------------------- 

  integer,intent(in) :: iunit, n
  complex(dp),intent(in) :: a((n*(n+1))/2)

  integer :: ntot, ij, ijt, i, j

  ij=0
  do j=1,n
    do i=1,j
      ij = ij + 1
      write(iunit,'(2i6,2e12.5)') i, j, a(ij)
    enddo
    ijt=ij
    do i=j+1,n
      ijt = ijt + (i-1)
      write(iunit,'(2i6,2e12.5)') i, j, conjg(a(ijt))
    enddo
    write(iunit,*)
  enddo

  return
  end subroutine dump_hamq_packed


! ---------------------------------------------------------------------- 
  subroutine write_nloper( iunit, nloper )
! ---------------------------------------------------------------------- 

  integer,intent(in) :: iunit
  type(matrix_list) :: nloper(:)

  integer :: iatom

  write(iunit) natom
  do iatom=1,natom
    write(iunit) size(nloper(iatom)%matrix,1)
    write(iunit) nloper(iatom)%matrix
  enddo
  
  end subroutine write_nloper

! ---------------------------------------------------------------------- 
  subroutine read_nloper( iunit, nloper )
! ---------------------------------------------------------------------- 

  use hamq_pool, only : nproc_per_pool

  integer,intent(in) :: iunit
  type(matrix_list) :: nloper(natom)

  if( nproc_per_pool > 1 ) then
    call read_nloper_pool( iunit, nloper )
  else
    call read_nloper_serial( iunit, nloper )
  endif

  end subroutine read_nloper

! ---------------------------------------------------------------------- 
  subroutine read_nloper_serial( iunit, nloper )
! ---------------------------------------------------------------------- 

  use mp_world, only : mpime, root

  integer,intent(in) :: iunit
  type(matrix_list) :: nloper(natom)

  integer :: iatom
  integer :: natom_, nproj_

  if( mpime==root ) then

  read(iunit) natom_
  if( natom_ /= natom ) &
    call errore('read_nloper','inconsistent number of atoms',abs(natom_))
  do iatom=1,natom
    read(iunit) nproj_
    if( nproj_ /= size(nloper(iatom)%matrix,1) ) &
      call errore('read_nloper','inconsistent number of projectors',abs(nproj_))
    read(iunit) nloper(iatom)%matrix
  enddo
  
  endif

  call bcast_nloper( nloper )

  end subroutine read_nloper_serial


! ---------------------------------------------------------------------- 
  subroutine read_nloper_pool( iunit, nloper )
! ---------------------------------------------------------------------- 

  use hamq_pool, only : natom_global, type_atom_global, &
                        mypoolid, mypoolroot, intra_pool_comm, &
                        mypool, rootpool, cross_pool_comm, &
                        nproc_per_pool
  use parallel_include
  use mp_world, only : mpime, root
  use mp, only : mp_bcast

  integer,intent(in) :: iunit
  type(matrix_list) :: nloper(natom)

  integer :: iatom_global
  integer :: iatom
  integer :: natom_, nproj_
  real(dp),allocatable :: matrix(:,:)
  integer :: dest, ierr
  integer :: istat(MPI_STATUS_SIZE)

  if( mpime==root ) then
    read(iunit) natom_
    if( natom_ /= natom_global ) &
      call errore('read_nloper','inconsistent number of atoms',abs(natom_))
  endif

  iatom=0
  do iatom_global=1,natom_global
    if( mpime==root ) then
      read(iunit) nproj_
!      if( nproj_ /= size(nloper(iatom)%matrix,1) ) &
!        call errore('read_nloper','inconsistent number of projectors',abs(nproj_))
    endif
    if( mypoolid==mypoolroot ) call mp_bcast( nproj_, rootpool, cross_pool_comm )
    if( mypoolid==mypoolroot ) allocate( matrix(nproj_,nproj_) )
    if( mpime==root ) then
      read(iunit) matrix
    endif
    if( mypoolid==mypoolroot ) call mp_bcast( matrix, rootpool, cross_pool_comm )

    ! send nproj_ to the correct process in the pool
    dest=mod(iatom_global-1,nproc_per_pool)
    if( mypoolid==mypoolroot ) then
      if( dest/=mypoolid ) then
        call mpi_send( nproj_, 1, MPI_INTEGER, dest, &
                       natom_global+iatom_global-1, &
                       intra_pool_comm, ierr )
      endif
    endif
    if( dest==mypoolid ) then
      if( mypoolid/=mypoolroot ) then
        call mpi_recv( nproj_, 1, MPI_INTEGER, mypoolroot, &
                       natom_global+iatom_global-1, &
                       intra_pool_comm, istat, ierr )
      endif
    endif

    ! now send the matrix
    if( mypoolid==mypoolroot ) then
      if( dest/=mypoolid ) then
        call mpi_send( matrix, nproj_*nproj_, MPI_DOUBLE_PRECISION, dest, &
                       iatom_global-1, intra_pool_comm, ierr )
      endif
    endif
    if( dest==mypoolid ) then
      iatom=iatom+1
      if( nproj_ /= size(nloper(iatom)%matrix,1) ) &
        call errore('read_nloper','inconsistent number of projectors',1)
      if( mypoolid/=mypoolroot ) then
        allocate( matrix(nproj_,nproj_) )
        call mpi_recv( matrix, nproj_*nproj_, MPI_DOUBLE_PRECISION, mypoolroot, &
                       iatom_global-1, intra_pool_comm, istat, ierr )
        nloper(iatom)%matrix = matrix
      else
        nloper(iatom)%matrix = matrix
      endif
    endif
    if( allocated(matrix) ) deallocate( matrix )
  enddo
  
  end subroutine read_nloper_pool


! ---------------------------------------------------------------------- 
  subroutine bcast_nloper( nloper )
! ---------------------------------------------------------------------- 

  ! this uses Cavazzoni's mp routines
  use mp, only : mp_bcast
  use mp_world, only : mpime, root, world_comm

  type(matrix_list) :: nloper(natom)

  integer :: iatom

  do iatom=1,natom
    call mp_bcast( nloper(iatom)%matrix, root, world_comm )
  enddo
  
  end subroutine bcast_nloper

! ---------------------------------------------------------------------- 
  subroutine spline_make_knots( xkxp, xkyp, xkzp )
! ---------------------------------------------------------------------- 

  use bspline90_22, only : dbsnak

  real(dp),intent(in) :: xkxp(nxcoef)
  real(dp),intent(in) :: xkyp(nycoef)
  real(dp),intent(in) :: xkzp(nzcoef)

  call dbsnak( nxcoef, xkxp, kxord, xknot )
  call dbsnak( nycoef, xkyp, kyord, yknot )
  call dbsnak( nzcoef, xkzp, kzord, zknot )

  return
  end subroutine spline_make_knots


! ---------------------------------------------------------------------- 
  subroutine spline_interpolate( xkxp, xkyp, xkzp, fin, bscoef )
! ---------------------------------------------------------------------- 

  use bspline90_22, only : dbsint, dbs2in, dbs3in

  real(dp),intent(in) :: xkxp(nxcoef)
  real(dp),intent(in) :: xkyp(nycoef)
  real(dp),intent(in) :: xkzp(nzcoef)
  real(dp),intent(in) :: fin(nxcoef*nycoef*nzcoef)
  real(dp),intent(out) :: bscoef(nxcoef*nycoef*nzcoef)

  ! 0-D
  if( nxcoef==1 .and. nycoef==1 .and. nzcoef==1 ) then
    bscoef(1) = fin(1)
  ! 1-D
  else if( nycoef==1 .and. nzcoef==1 ) then
    call dbsint( nxcoef, xkxp, &
                 fin, &
                 kxord, xknot, bscoef(1) )
  else if( nxcoef==1 .and. nzcoef==1 ) then
    call dbsint( nycoef, xkyp, &
                 fin, &
                 kyord, yknot, bscoef(1) )
  else if( nxcoef==1 .and. nycoef==1 ) then
    call dbsint( nzcoef, xkzp, &
                 fin, &
                 kzord, zknot, bscoef(1) )
  ! 2-D
  else if( nzcoef==1 ) then
    call dbs2in( nxcoef, xkxp, nycoef, xkyp, &
                 fin, nxcoef, &
                 kxord, kyord, xknot, yknot, bscoef(1) )
  else if( nycoef==1 ) then
    call dbs2in( nxcoef, xkxp, nzcoef, xkzp, &
                 fin, nxcoef, &
                 kxord, kzord, xknot, zknot, bscoef(1) )
  else if( nxcoef==1 ) then
    call dbs2in( nycoef, xkyp, nzcoef, xkzp, &
                 fin, nycoef, &
                 kyord, kzord, yknot, zknot, bscoef(1) )
  ! 3-D
  else
    call dbs3in( nxcoef, xkxp, nycoef, xkyp, nzcoef, xkzp, &
                 fin, nxcoef, nycoef, &
                 kxord, kyord, kzord, xknot, yknot, zknot, bscoef(1) )
  endif

  
  end subroutine spline_interpolate


! ---------------------------------------------------------------------- 
  function spline_evaluate( qvec, bscoef )

  use mp_world, only : mpime
  use bspline90_22, only : dbsval, dbs2vl, dbs3vl

  real(dp),intent(in) :: qvec(3)
  real(dp),intent(in) :: bscoef(nxcoef*nycoef*nzcoef)
  real(dp) :: spline_evaluate
  real(dp) :: fout

  ! 0-D
  if( nxcoef==1 .and. nycoef==1 .and. nzcoef==1 ) then
    fout = bscoef(1)
  ! 1-D
  else if( nycoef==1 .and. nzcoef==1 ) then
    fout = dbsval(qvec(1),kxord,xknot,nxcoef,bscoef(1))
  else if( nxcoef==1 .and. nzcoef==1 ) then
    fout = dbsval(qvec(2),kyord,yknot,nycoef,bscoef(1))
  else if( nxcoef==1 .and. nycoef==1 ) then
    fout = dbsval(qvec(3),kzord,zknot,nzcoef,bscoef(1))
  ! 2-D
  else if( nzcoef==1 ) then
    fout = dbs2vl(qvec(1),qvec(2),kxord,kyord,xknot,yknot,nxcoef,nycoef,bscoef(1))
  else if( nycoef==1 ) then
    fout = dbs2vl(qvec(1),qvec(3),kxord,kzord,xknot,zknot,nxcoef,nzcoef,bscoef(1))
  else if( nxcoef==1 ) then
    fout = dbs2vl(qvec(2),qvec(3),kyord,kzord,yknot,zknot,nycoef,nzcoef,bscoef(1))
  ! 3-D
  else
    fout = dbs3vl( qvec(1), qvec(2), qvec(3), &
                   kxord, kyord, kzord, xknot, yknot, zknot, &
                   nxcoef, nycoef, nzcoef, bscoef(1) )
  endif
  spline_evaluate = fout

  return
  end function spline_evaluate

! ---------------------------------------------------------------------- 
  function spline_derivative( ideriv, qvec, bscoef )

  use bspline90_22, only : dbsder, dbs2dr, dbs3dr

  integer,intent(in) :: ideriv
  real(dp),intent(in) :: qvec(3)
  real(dp),intent(in) :: bscoef(nxcoef*nycoef*nzcoef)
  real(dp) :: spline_derivative
  real(dp) :: fout

  fout = 0.d0  ! default
  ! 0-D
  if( nxcoef==1 .and. nycoef==1 .and. nzcoef==1 ) then
    fout = 0.d0
  ! 1-D
  else if( nycoef==1 .and. nzcoef==1 .and. ideriv==1 ) then
    fout = dbsder(1,qvec(1),kxord,xknot,nxcoef,bscoef(1))
  else if( nxcoef==1 .and. nzcoef==1 .and. ideriv==2 ) then
    fout = dbsder(1,qvec(2),kyord,yknot,nycoef,bscoef(1))
  else if( nxcoef==1 .and. nycoef==1 .and. ideriv==3 ) then
    fout = dbsder(1,qvec(3),kzord,zknot,nzcoef,bscoef(1))
  ! 2-D
  else if( nzcoef==1 ) then
    if( ideriv==1 ) then
      fout = dbs2dr(1,0,qvec(1),qvec(2),kxord,kyord,xknot,yknot,nxcoef,nycoef,bscoef(1))
    else if ( ideriv==2 ) then
      fout = dbs2dr(0,1,qvec(1),qvec(2),kxord,kyord,xknot,yknot,nxcoef,nycoef,bscoef(1))
    endif
  else if( nycoef==1 ) then
    if( ideriv==1 ) then
      fout = dbs2dr(1,0,qvec(1),qvec(3),kxord,kzord,xknot,zknot,nxcoef,nzcoef,bscoef(1))
    else if( ideriv==3 ) then
      fout = dbs2dr(0,1,qvec(1),qvec(3),kxord,kzord,xknot,zknot,nxcoef,nzcoef,bscoef(1))
    endif
  else if( nxcoef==1 ) then
    if( ideriv==2 ) then
      fout = dbs2dr(1,0,qvec(2),qvec(3),kyord,kzord,yknot,zknot,nycoef,nzcoef,bscoef(1))
    else if( ideriv==3 ) then
      fout = dbs2dr(0,1,qvec(2),qvec(3),kyord,kzord,yknot,zknot,nycoef,nzcoef,bscoef(1))
    endif
  ! 3-D
  else
    if( ideriv==1 ) then
      fout = dbs3dr( 1,0,0,qvec(1), qvec(2), qvec(3), &
                     kxord, kyord, kzord, xknot, yknot, zknot, &
                     nxcoef, nycoef, nzcoef, bscoef(1) )
    else if( ideriv==2 ) then
      fout = dbs3dr( 0,1,0,qvec(1), qvec(2), qvec(3), &
                     kxord, kyord, kzord, xknot, yknot, zknot, &
                     nxcoef, nycoef, nzcoef, bscoef(1) )
    else if( ideriv==3 ) then
      fout = dbs3dr( 0,0,1,qvec(1), qvec(2), qvec(3), &
                     kxord, kyord, kzord, xknot, yknot, zknot, &
                     nxcoef, nycoef, nzcoef, bscoef(1) )
    endif
  endif
  spline_derivative = fout

  return
  end function spline_derivative


! ---------------------------------------------------------------------- 
  end module hamq_shirley
! ---------------------------------------------------------------------- 
