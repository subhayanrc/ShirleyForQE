  program valence_corerepair

! Read the pseudo-waves
! Ask for corresponding ae wave
! Generate a momentum matrix for each pair of channels
!   for ae and ps
! Subtract to give corerepair term
!   <ae_l|p|ae_m> - <ps_l|p|ps_m>

  use kinds
  use fileio
  use atomic_waves

  implicit none

  real(dp),parameter :: eps12=1.d-12

  character(255) :: fmtstr
  
  integer :: nwfc
  type(atomic),allocatable :: ps_wfc(:)
  type(atomic),allocatable :: ae_wfc(:)
  real(dp),allocatable :: rcut(:)
  real(dp) :: rc
  logical :: tnormcnsv
  real(dp) :: alpha
  real(dp),allocatable :: ps_mommat(:,:,:)
  real(dp),allocatable :: ae_mommat(:,:,:)
  complex(dp),allocatable :: corerep(:,:,:)

  integer :: ngrid, l, igrid
  integer :: iwfc, jwfc
  integer :: ixyz, i0, j0
  real(dp) :: rdum, rcut_eff, dx
  integer :: immax, jmmax, im, jm
  integer :: nproj, iproj, jproj
  integer :: nnonzero
  real(dp) :: crm

  ! header to explain type of core-repair
  write(*,'(2x,a)') 'momentum'

  ! read number of wave functions
  read(*,*) nwfc
  write(*,'(2x,2i4,t24,a)') nwfc, nwfc, '! nwfc1, nwfc2'

  ! read the wave functions
  allocate( ps_wfc(nwfc), ae_wfc(nwfc), rcut(nwfc) )
  do iwfc=1,nwfc
    read(*,*) ngrid, l, rc
    rcut(iwfc) = rc

    call init_atomic_wave( ngrid, l, ps_wfc(iwfc) )
    do igrid=1,ngrid
      read(*,*) ps_wfc(iwfc)%r(igrid), ps_wfc(iwfc)%f(igrid)
    enddo

    read(*,*) ngrid, l, rc
    ! check that the l-value matches
    if( l /= ps_wfc(iwfc)%l ) then
      write(*,*) ' error: angular momentum of ae_wfc does not match ps_wfc'
      stop
    endif
    ! check that this is the ae_wfc
    if( rc /= 0.d0 ) then
      write(*,*) ' error: this should be an all-electron wave function'
      write(*,*) '        the cut-off radius should be zero'
      stop
    endif

    call init_atomic_wave( ngrid, l, ae_wfc(iwfc) )
    do igrid=1,ngrid
      read(*,*) ae_wfc(iwfc)%r(igrid), ae_wfc(iwfc)%f(igrid)
    enddo
  enddo

  ! angular momenta
  write(fmtstr,*) '(2x,', nwfc, 'i2,t24,a)'
  write(*,fmtstr) ae_wfc(:)%l, '! lwfc1(1:nwfc1)'
  write(*,fmtstr) ae_wfc(:)%l, '! lwfc2(1:nwfc2)'
  ! Compute momentum matrices
  nproj=0
  do iwfc=1,nwfc
    nproj=nproj+2*ae_wfc(iwfc)%l+1
  enddo
  allocate( corerep(3,nproj,nproj) )
  corerep = 0.d0

  i0=0
  do iwfc=1,nwfc
    j0=0
    do jwfc=1,nwfc
      rcut_eff = max(rcut(iwfc),rcut(jwfc))
      immax = 2 * ps_wfc(iwfc)%l + 1
      jmmax = 2 * ps_wfc(jwfc)%l + 1
      allocate( ps_mommat(3,jmmax,immax), &
                ae_mommat(3,jmmax,immax) )

      call momentum_matrix( &
        ps_wfc(jwfc)%ngrid, ps_wfc(jwfc)%l, ps_wfc(jwfc)%r, ps_wfc(jwfc)%f, &
        ps_wfc(iwfc)%ngrid, ps_wfc(iwfc)%l, ps_wfc(iwfc)%r, ps_wfc(iwfc)%f, &
        rcut_eff, ps_mommat &
                            )

      call momentum_matrix( &
        ae_wfc(jwfc)%ngrid, ae_wfc(jwfc)%l, ae_wfc(jwfc)%r, ae_wfc(jwfc)%f, &
        ae_wfc(iwfc)%ngrid, ae_wfc(iwfc)%l, ae_wfc(iwfc)%r, ae_wfc(iwfc)%f, &
        rcut_eff, ae_mommat &
                          )
      ! Subtract ae-ps and dump
      do im=1,immax
      do jm=1,jmmax
        do ixyz=1,3
          corerep(ixyz,j0+jm,i0+im) = (ae_mommat(ixyz,jm,im)-ps_mommat(ixyz,jm,im))
        enddo
      enddo
      enddo

      deallocate( ps_mommat )
      deallocate( ae_mommat )
      j0=j0+jmmax
    enddo
    i0=i0+immax
  enddo
    
  ! multiply by factor of -i
  corerep = corerep * cmplx(0.d0,-1.d0)

  nnonzero=0
  do iproj=1,nproj
    do jproj=1,nproj
      do ixyz=1,3
        crm = conjg(corerep(ixyz,jproj,iproj))*corerep(ixyz,jproj,iproj)
        if( crm > eps12 ) then
          nnonzero=nnonzero+1
        endif
      enddo
    enddo
  enddo

  write(*,'(2x,i6,t24,a)') nnonzero, '! nonzero elements (i,j,ixyz,cR,cI)'

  do iproj=1,nproj
    do jproj=1,nproj
      do ixyz=1,3
        crm = conjg(corerep(ixyz,jproj,iproj))*corerep(ixyz,jproj,iproj)
        if( crm > eps12 ) then
          write(*,'(2x,3i3,x,2e20.12)') jproj, iproj, ixyz, corerep(ixyz,jproj,iproj)
        endif
      enddo
    enddo
  enddo

  stop

  end program valence_corerepair
