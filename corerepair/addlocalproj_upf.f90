!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
module upf_wrapper
  use upf

  contains

  subroutine dump_pseudo_upf( upfin )
  use funct, only: set_dft_from_name, get_iexch, get_icorr, get_igcx, get_igcc
  use pseudo_types, only: pseudo_upf
  type(pseudo_upf),intent(in) :: upfin
  character(len=26) :: dft_

  ! pp_info
  rel=0
  select case (upfin%rel)
    case ('no')
      rel=0
    case ('scalar')
      rel=1
    case ('full')
      rel=2
  end select
  rcloc=upfin%rcloc
  nwfs=upfin%nwfc
  allocate(oc(nwfs))
  allocate(rcut(nwfs))
  allocate(rcutus(nwfs))
  allocate(epseu(nwfs))
  allocate(els(nwfs))
  els=''
  allocate(lchi(nwfs))
  allocate(nns(nwfs))
  if(associated(upfin%oc) .and. .not. allocated(oc)) then
    oc=upfin%oc
  endif
  if(associated(upfin%rcut) .and. .not. allocated(rcut)) then
    rcut=upfin%rcut
  endif
  if(associated(upfin%rcutus) .and. .not. allocated(rcutus)) then
    rcutus=upfin%rcutus
  endif
  if(associated(upfin%epseu) .and. .not. allocated(epseu)) then
    epseu=upfin%epseu
  endif
  if(associated(upfin%els_beta) .and. .not. allocated(els)) then
    els=upfin%els_beta
  endif
  if(associated(upfin%lchi) .and. .not. allocated(lchi)) then
    lchi=upfin%lchi
  endif
  if(associated(upfin%nn) .and. .not. allocated(nns)) then
    nns=upfin%nn
  endif

  ! pp_header
  generated=upfin%generated
  date_author=upfin%date // ' ' // upfin%author
  comment=upfin%comment
  psd=upfin%psd
  pseudotype=upfin%typ
  !nv=upfin%nv
  call set_dft_from_name( upfin%dft )
  iexch=get_iexch()
  icorr=get_icorr()
  igcx=get_igcx()
  igcc=get_igcc()
  lmax=upfin%lmax
  mesh=upfin%mesh
  nbeta=upfin%nbeta
  ntwfc=upfin%nwfc
  nlcc=upfin%nlcc
  zp=upfin%zp
  ecutrho=upfin%ecutrho
  ecutwfc=upfin%ecutwfc
  etotps=upfin%etotps
  if(associated(upfin%oc) .and. .not. allocated(ocw)) then
    allocate(ocw(size(upfin%oc)))
    ocw=upfin%oc
  endif
  if(associated(upfin%els) .and. .not. allocated(elsw)) then
    allocate(elsw(size(upfin%els)))
    elsw=upfin%els
  endif
  if(associated(upfin%lchi) .and. .not. allocated(lchiw)) then
    allocate(lchiw(size(upfin%lchi)))
    lchiw=upfin%lchi
  endif

  ! pp_mesh
  if(associated(upfin%r) .and. .not. allocated(r)) then
    allocate(r(size(upfin%r)))
    r=upfin%r
  endif
  if(associated(upfin%rab) .and. .not. allocated(rab)) then
    allocate(rab(size(upfin%rab)))
    rab=upfin%rab
  endif
  
  ! pp_nlcc
  if(associated(upfin%rho_atc) .and. .not. allocated(rho_atc)) then
    allocate(rho_atc(size(upfin%rho_atc)))
    rho_atc=upfin%rho_atc
  endif

  ! pp_local
  if(associated(upfin%vloc) .and. .not. allocated(vloc0)) then
    allocate(vloc0(size(upfin%vloc)))
    vloc0=upfin%vloc
  endif

  ! pp_nonlocal
  ! pp_beta
  if(associated(upfin%beta) .and. .not. allocated(betar)) then
    allocate(betar(size(upfin%beta,1),size(upfin%beta,2)))
    betar=upfin%beta
  endif
  if(associated(upfin%lll) .and. .not. allocated(lll)) then
    allocate(lll(size(upfin%lll)))
    lll=upfin%lll
  endif
  if(associated(upfin%kbeta) .and. .not. allocated(ikk2)) then
    allocate(ikk2(size(upfin%kbeta)))
    ikk2=upfin%kbeta
  endif
  ! pp_dij
  if(associated(upfin%dion) .and. .not. allocated(dion)) then
    allocate(dion(size(upfin%dion,1),size(upfin%dion,2)))
    dion=upfin%dion
  endif
  ! pp_qij
  nqf=upfin%nqf
  nqlc=upfin%nqlc
  if(associated(upfin%rinner) .and. .not. allocated(rinner)) then
    allocate(rinner(size(upfin%rinner)))
    rinner=upfin%rinner
  endif
  if(associated(upfin%qqq) .and. .not. allocated(qqq)) then
    allocate(qqq(size(upfin%qqq,1),size(upfin%qqq,2)))
    qqq=upfin%qqq
  endif
  if(associated(upfin%qfuncl) .and. .not. allocated(qfunc)) then
    allocate(qfunc(size(upfin%qfuncl,1),size(upfin%qfuncl,2),size(upfin%qfuncl,3)))
    qfunc=upfin%qfuncl
  endif
  ! pp_qfcoef
  if(associated(upfin%qfcoef) .and. .not. allocated(qfcoef)) then
    allocate(qfcoef(size(upfin%qfcoef,1),size(upfin%qfcoef,2),size(upfin%qfcoef,3),size(upfin%qfcoef,4)))
    qfcoef=upfin%qfcoef
  endif
  
  ! pp_pswfc
  if(associated(upfin%chi) .and. .not. allocated(chi)) then
    allocate(chi(size(upfin%chi,1),size(upfin%chi,2)))
    chi=upfin%chi
  endif
  
  ! pp_rhoatom
  if(associated(upfin%rho_at) .and. .not. allocated(rho_at)) then
    allocate(rho_at(size(upfin%rho_at)))
    rho_at=upfin%rho_at
  endif

  end subroutine dump_pseudo_upf

end module upf_wrapper

!---------------------------------------------------------------------
program addlocalproj_upf
  !---------------------------------------------------------------------
  !
  !  Read pseudopotentials in the Unified Pseudopotential Format (UPF)
  !  Add a local projector to the UPF file
  !  Dump to a new file in UPF again
  !
  use upf_module
  use pseudo_types, only: pseudo_upf, nullify_pseudo_upf
  use upf_wrapper, only: dump_pseudo_upf
  USE radial_grids, ONLY: radial_grid_type, nullify_radial_grid
  USE environment, ONLY: environment_start, environment_end
  USE mp_global, ONLY: mp_startup, mp_global_end
  USE io_global, ONLY: ionode, stdout

  implicit none
  !
  real(8),parameter :: eps12=1.d-12
  real(8),parameter :: lambda=6.d0 ! exponent for local projector
  !
  integer :: ios
  integer :: iunps = 4
  character (len=256) :: filein
  type(pseudo_upf) :: upf
  !! derived type where is possible to store data on the radial mesh
  TYPE(radial_grid_type),TARGET :: grid
  integer :: nb, lloc, ichi_loc, i, ikk_loc
  integer,allocatable :: map(:)
  real(8) :: rcut_loc, epsloc, dion_loc
  !
  ! pp_beta
  real(8), allocatable :: betar_old(:,:)
  integer, allocatable:: lll_old(:), ikk2_old(:)
  ! pp_dij
  real(8), allocatable :: dion_old(:,:)
  ! pp_pswfc
  real(8), allocatable :: chi_loc(:)
  real(8), allocatable :: betar_loc(:)
  !
  real(8),external :: intradialprod
  !
  logical, external :: matches
  !
#if defined(__MPI)
   CALL mp_startup()
#endif
   CALL environment_start('addlocalproj_upf')
   IF (ionode) THEN

  !---------------------------------------------------------------------- 
  ! read
  print '(''  Input PP file in UPF format : '',$)'
  read (*, '(a)', end = 20, err = 20) filein
  open(unit=iunps,file=filein,status='old',form='formatted',iostat=ios)
  if (ios.ne.0) stop
  CALL nullify_pseudo_upf ( upf )
  CALL nullify_radial_grid ( grid )
  upf%grid => grid
  call read_upf(upf, unit=iunps, grid=grid, ierr=ios)
  close (unit=iunps)

  !---------------------------------------------------------------------- 
  ! do something to pseudo
  ! for some reason there is not a logical to note that we are norm-conserving
  if( .not. (matches ('NC', upf%typ) ) ) then
    write(*,*) 'pseudo type = ', upf%typ
    write(*,*) ' error: this code expects a norm-conserving pseudo'
    goto 20
  endif

  write(*,*) ' number of projectors = ', upf%nbeta
  do nb=1,upf%nbeta
    write(*,*) ' proj ', nb, ' l=', upf%lll(nb)
  enddo
  write(*,*) ' number of atomic wave functions = ', upf%nwfc
  do nb=1,upf%nwfc
    write(*,*) '  chi ', nb, ' l=', upf%lchi(nb)
  enddo
  write(*,*) ' number of projectors and wave functions differ by ', upf%nwfc-upf%nbeta
  if( upf%nwfc-upf%nbeta /= 1 ) then
    write(*,*) ' error: this code expects only one missing projector'
    goto 20
  endif
  !print '(''  Angular momentum of local channel : '',$)'
  !read (*,*,end=20,err=20) lloc
  allocate( map(upf%nwfc) )
  map = -1
  do i=1,upf%nwfc
    do nb=1,upf%nbeta
      if( upf%lll(nb) == upf%lchi(i) ) then
        write(*,*) ' beta ', nb, ' matches chi ', i
        map(i)=nb
      endif
    enddo
  enddo
  forall( i=1:upf%nwfc, map(i)==-1 ) ichi_loc=i
  lloc=upf%lchi(ichi_loc)
  write(*,*) ' Angular momentum of local channel chosen to be ', lloc

  ! check on existing wave functions
  ! this might be impossible
  !if( upf%nwfc == 0 ) then
  !  write(*,*) ' There are no atomic wave functions for generating this projector'
  !  if( lloc > 0 ) then
  !    write(*,*) ' hydrogenic wave implemented only for s-waves'
  !    goto 20
  !  endif
  !  write(*,*) ' Using a hydrogenic wave instead'
  !  ichi_loc=0
  !else
  !  do nb=1,upf%nwfc
  !    if( lloc==upf%lchi(nb) ) write(*,*) ' chi ', nb, ' has the correct ang mom'
  !  enddo
  !  print '(''  Pick which chi to use in constructing the local projector > '',$)'
  !  read(*,*,end=20,err=20) ichi_loc
  !  if( ichi_loc < 1 .or. ichi_loc > upf%nwfc ) then
  !    write(*,*) ' error: chosen chi is out of bounds' ; goto 20
  !  else if( upf%lchi(ichi_loc)/=lloc ) then
  !    write(*,*) ' error: chosen chi has the wrong ang momentum'; goto 20
  !  endif
  !endif

  ! cut-off
  ! recompute upf%kbeta
  do nb=1,upf%nbeta
    i=size(upf%beta,1)
    do while( i > 1 .and. abs(upf%beta(i,nb)) < eps12 )
      i=i-1
    enddo
    upf%kbeta(nb) = i
  enddo
  if( upf%nbeta > 0 ) then
    ikk_loc = maxval(upf%kbeta(1:upf%nbeta)) 
    rcut_loc = upf%r(ikk_loc)
  else
    rcut_loc = 0.5d0
    do i=1,size(upf%r)
      if( upf%r(i) > rcut_loc ) then
        ikk_loc = i-1
        exit
      endif
    enddo
  endif
  write(*,*) ' chosen cut-off radius = ', rcut_loc

  ! copy old projectors and resize
  if( upf%nbeta > 0 ) then
    ! pp_beta
    allocate( betar_old(size(upf%beta,1),size(upf%beta,2)), &
              lll_old(size(upf%lll,1)), &
              ikk2_old(size(upf%kbeta,1)) )
    betar_old = upf%beta
    lll_old = upf%lll
    ikk2_old = upf%kbeta
    deallocate( upf%beta ) ; allocate( upf%beta(size(betar_old,1),size(betar_old,2)+1) )
    deallocate( upf%lll ) ; allocate( upf%lll(size(lll_old,1)+1) )
    deallocate( upf%kbeta ) ; allocate( upf%kbeta(size(ikk2_old,1)+1) )

    ! pp_dij
    allocate( dion_old(size(upf%dion,1),size(upf%dion,2)) )
    dion_old = upf%dion
    deallocate( upf%dion ) ; allocate( upf%dion(size(dion_old,1)+1,size(dion_old,2)+1) )
  else
    ! if no projectors then make space for one
    allocate(upf%beta(size(betar_loc),1))
    allocate(upf%lll(1))
    allocate(upf%kbeta(1))
    allocate(upf%dion(1,1))
  endif

  ! generate chi_loc and betar_loc
  allocate( chi_loc(size(upf%r)), betar_loc(size(upf%r)) )
  if( ichi_loc<1 ) then
    ! this assumes that if we have no projector at all, which seems likely for the H atom only
    ! maybe this should be rewritten for a general Z value, not assuming Z=1

    ! compute chi_loc
    epsloc = (-log(eps12))
    ! compute ikk_loc
    i=size(upf%r)
    do while( i > 1 .and. upf%r(i) > epsloc )
      i=i-1
    enddo
    ikk_loc = i
    chi_loc = 0.d0
    do i=1,ikk_loc
      ! hydrogenic s-wave
      chi_loc(i) = exp( -upf%r(i) )
    enddo
    ! why the factor of 2 here - surely it's irrelevant
    chi_loc = 2.d0 * upf%r * chi_loc  ! dont forget that we store r*chi_loc

    ! compute betar_loc
    epsloc = ((-log(eps12))**(1.d0/lambda))*rcut_loc
    ! compute ikk_loc
    i=size(upf%r)
    do while( i > 1 .and. upf%r(i) > epsloc )
      i=i-1
    enddo
    ikk_loc = i
    betar_loc = 0.d0
    do i=1,ikk_loc
      ! hydrogenic s-wave
      betar_loc(i) = exp( -(upf%r(i)/rcut_loc)**lambda )*chi_loc(i)
    enddo

  else
    ! now, assuming that there was a chi_loc that matches the lloc value

    ! modify rcut_loc to give beta_loc=0 outside the cut-off radius
    rcut_loc = rcut_loc / ((-log(eps12) )**(1.d0/lambda))
    chi_loc = upf%chi(:,ichi_loc)
    betar_loc = 0.d0
    do i=1,ikk_loc
      betar_loc(i) = exp(-(upf%r(i)/rcut_loc)**lambda)*chi_loc(i)
    enddo
    
  endif

  ! upf%dion = 1 / integral beta * chi
  dion_loc = 1.d0 / intradialprod( ikk_loc, upf%r, betar_loc, chi_loc )

  ! update pseudo-variables
  ! to match the ordering of the chi wave functions
  upf%nbeta = upf%nbeta+1
  upf%dion=0.d0
  do nb=1,upf%nbeta
    if( map(nb) /= -1) then
      upf%beta(:,nb) = betar_old(:,map(nb))
      upf%lll(nb) = lll_old(map(nb))
      upf%kbeta(nb) = ikk2_old(map(nb))
      ! assumption for NCPP's that dion is diagonal
      upf%dion(nb,nb) = dion_old(map(nb),map(nb))
    else
      upf%beta(:,nb) = betar_loc
      upf%lll(nb) = lloc
      upf%kbeta(nb) = ikk_loc
      upf%dion(nb,nb) = dion_loc
    endif
  enddo
    
  write(*,*) ' new local projector has been made'

  !---------------------------------------------------------------------- 
  ! write
  print '(''  Output PP file in UPF format > '',$)'
  read (*, '(a)', end = 20, err = 20) filein
  open(unit=iunps,file=filein,form='formatted',iostat=ios)
  if (ios.ne.0) stop
  call dump_pseudo_upf(upf)
  call write_upf_v1(iunps)
  close (unit=iunps)

   END IF
   CALL environment_end('addlocalproj_upf')
#if defined(__MPI) 
   CALL mp_global_end()
#endif 

20 stop
end program addlocalproj_upf
!
