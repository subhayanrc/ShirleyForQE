  program valence_overlap

  ! David Prendergast, LBNL, July 2018

    use kinds
    use fileio
    use atomic_waves
    use splines
  
    implicit none
  
#ifdef __PGI
    integer,external :: iargc
#else
    integer,intrinsic :: iargc
#endif
  
    real(dp), parameter :: zero   =  0.d0
    real(dp), parameter :: eps12  =  1.d-12
  
    integer :: narg,  iunc,  iunv
    character(255) :: filename,  fmtstr
    
    type(spline_struct) :: spl

    type(atomic) :: core_wave,  tmp_wave
    type(atomic), allocatable :: ae_g(:), ps_g(:)
    type(atomic), allocatable :: ae_x(:), ps_x(:)
    real(dp) :: r, phi
    real(dp) :: rmax
    real(dp), allocatable :: rm(:)
    integer :: imax, imax1(1)
    real(dp), allocatable :: rgrid(:), fgrid(:), qgrid(:)
    real(dp), allocatable :: aeg(:,:), psg(:,:)
    real(dp), allocatable :: aex(:,:), psx(:,:)

    real(dp), allocatable :: qqq(:,:)
  
    integer, parameter :: ideg = 5
    type(spline_struct), allocatable :: spl_ae_g(:)
  
    integer :: ngrid,  l,  igrid
    integer :: ierr,  nnonzero
  
    integer :: nvalence_g, nvalence_x
    integer :: ivalence,  jvalence
  
    character(3) :: orbit(1 : 9)  =  [ character(3) :: 'x', 'y', 'z', 'x^2', 'y^2', 'z^2', 'xy', 'yz', 'zx' ]
  
    ! Get user input
    narg  =  iargc()
    if( narg /=  2 ) then
      write(0,  * ) 'usage: valence_overlap valence_g.dat valence_x.dat'
      stop
    endif
  
    call getarg( 1,  filename )
    ! open file
    iunv  =  freeunit()
    open(iunv, file = trim(filename), form = 'formatted', err = 911)
  
    ! Read the GS valence wave
#ifdef __DEBUG
    write( * ,  * ) 'reading GS valence wave'
#endif
    read(iunv,  * ) nvalence_g
  
    allocate( ae_g(nvalence_g), ps_g(nvalence_g) )

    do ivalence = 1, nvalence_g
    
      call read_wave( iunv, ps_g(ivalence) )
      call read_wave( iunv, ae_g(ivalence) )

    enddo

    close( iunv )


    ! do the same for the ps_x and ae_x waves
    ! I need the pseudo cut-offs for integration and then I need to integrate only up to those cut-offs or know that the functions will go to zero beyond them
    ! then spline on the best grid (pick rmax and ngrid max)
    ! then integrate 4-way products ae_g ae_x - ps_g ps_x = Q
    ! dump and leave


    call getarg( 2,  filename )
    ! open file
    iunv  =  freeunit()
    open(iunv, file = trim(filename), form = 'formatted', err = 911)
  
    ! Read the XS valence wave
#ifdef __DEBUG
    write( * ,  * ) 'reading XS valence wave'
#endif
    read(iunv,  * ) nvalence_x

    allocate( ae_x(nvalence_x), ps_x(nvalence_x) )

    do ivalence = 1, nvalence_x
    
      call read_wave( iunv, ps_x(ivalence) )
      call read_wave( iunv, ae_x(ivalence) )

    enddo

    close( iunv )

    allocate( rm(nvalence_g+nvalence_x) )
    do ivalence = 1, nvalence_g
      rm(ivalence) = ps_g(ivalence)%r(ps_g(ivalence)%ngrid)
    enddo
    do ivalence = 1, nvalence_x
      rm(ivalence+nvalence_g) = ps_x(ivalence)%r(ps_x(ivalence)%ngrid)
    enddo
    rmax = maxval( rm(:) )
    imax1 = maxloc( rm(:) )
    imax = imax1(1)
    if( imax > nvalence_g ) then
      ngrid = ps_x(imax-nvalence_g)%ngrid
      allocate( rgrid(ngrid) )
      rgrid = ps_x(imax-nvalence_g)%r
    else
      ngrid = ps_g(imax)%ngrid
      allocate( rgrid(ngrid) )
      rgrid = ps_g(imax)%r
    endif
    allocate( fgrid(size(rgrid)) )
    allocate( qgrid(size(rgrid)) )
 
    !write(*,*) " rmax = ", rmax, rgrid(1), rgrid(ngrid)

    allocate( psg(ngrid,nvalence_g), aeg(ngrid,nvalence_g), &
              psx(ngrid,nvalence_x), aex(ngrid,nvalence_x) )

    do ivalence = 1, nvalence_g
    
      call spline_fit( ideg, ps_g(ivalence)%ngrid, ps_g(ivalence)%r, ps_g(ivalence)%f, &
                       spl, ierr, zero, rmax )
      call spline_eval( spl, ngrid, rgrid, psg(:,ivalence), ierr )

      call spline_fit( ideg, ae_g(ivalence)%ngrid, ae_g(ivalence)%r(:), ae_g(ivalence)%f(:), &
                       spl, ierr, zero, rmax )
      call spline_eval( spl, ngrid, rgrid, aeg(:,ivalence), ierr )

    enddo

    do ivalence = 1, nvalence_x
    
      call spline_fit( ideg, ps_x(ivalence)%ngrid, ps_x(ivalence)%r, ps_x(ivalence)%f, &
                       spl, ierr, zero, rmax )
      call spline_eval( spl, ngrid, rgrid, psx(:,ivalence), ierr )

      call spline_fit( ideg, ae_x(ivalence)%ngrid, ae_x(ivalence)%r, ae_x(ivalence)%f, &
                       spl, ierr, zero, rmax )
      call spline_eval( spl, ngrid, rgrid, aex(:,ivalence), ierr )

    enddo

    ! allocate qqq
    allocate( qqq(nvalence_g,nvalence_x) )
    do ivalence=1,nvalence_g
    do jvalence=1,nvalence_x

      if( ae_g(ivalence)%l == ae_x(jvalence)%l ) then

        qgrid = aeg(:,ivalence) * aex(:,jvalence) - psg(:,ivalence) * psx(:,jvalence)
        call spline_fit( ideg, ngrid, rgrid, qgrid, spl, ierr, zero, rmax )
        qqq(ivalence,jvalence) = spline_integral( spl, zero, rmax )

      else

        qqq(ivalence,jvalence) = zero

      endif

    enddo
    write(fmtstr,'(a1,i6,a7)') '(', nvalence_x, 'e20.12)'
    write(*,fmtstr) qqq(ivalence,1:nvalence_x)  
    enddo

    stop

911 write(0,*) 'corevalence_position : unable to open file ', trim(filename)
    stop

  contains

    subroutine read_wave( iun, phi )

    integer,intent(in) :: iun
    type(atomic),intent(out) :: phi

    integer :: ngrid, l, igrid, jgrid
    real(dp) :: r, f
      

      ! read the pseudo wave first
      read(iunv,  * ) ngrid,  l
      read(iunv,  * ) r, f
      if( r > zero ) ngrid=ngrid+1
      call init_atomic_wave( ngrid,  l,  phi )
      if( r > zero ) then
        phi%r(1)=zero; phi%f(1)=zero
        phi%r(2)=r;    phi%f(2)=f
        jgrid=3
      else
        phi%r(1)=r; phi%f(1)=f
        jgrid=2
      endif
      do igrid=jgrid,ngrid
        read(iunv,  * ) phi%r(igrid), phi%f(igrid)
      enddo

      end subroutine read_wave
  end program valence_overlap
