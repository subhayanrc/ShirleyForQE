  module splines_module

! Interface to bpspline90_22
! assumes periodic splines

  use kinds, only : dp
  use bspline90_22

  implicit none

  type spline_dataset_type
    integer :: nx, ny, nz
    real(dp),pointer :: x(:), y(:), z(:)
    type(grid3D_type),pointer :: datasetR(:,:)
    type(grid3D_type),pointer :: datasetI(:,:)
    real(dp) :: transcart(3,3), transgrid(3,3)
    integer,pointer :: gridindex(:,:)
  end type

  type spline_fit_type
    real(dp),pointer :: xknot(:), yknot(:), zknot(:)
    integer :: xorder, yorder, zorder
    integer :: nxcoef, nycoef, nzcoef
    type(spline_coef_type),pointer :: coefR(:,:)
    type(spline_coef_type),pointer :: coefI(:,:)
  end type

  type grid3D_type
    real(dp),pointer :: gxyz(:,:,:)
  end type
    
  type spline_coef_type
    real(dp),pointer :: bscoef(:)
  end type

  
  contains

! ----------------------------------------------------------------------
  subroutine create_spline_dataset( n, m, ngrid, igrid, &
                                    transcart, transgrid, ds )
! ----------------------------------------------------------------------
  integer,intent(in) :: n, m
  integer,intent(in) :: ngrid(3), igrid(3)
  real(dp),intent(in) :: transcart(3,3), transgrid(3,3)
  type(spline_dataset_type) :: ds

  integer :: i, j, k, ijk

  ds%nx = ngrid(1) + igrid(1) + 1
  ds%ny = ngrid(2) + igrid(2) + 1
  ds%nz = ngrid(3) + igrid(3) + 1
  allocate( ds%x( ds%nx ) )
  allocate( ds%y( ds%ny ) )
  allocate( ds%z( ds%nz ) )

  allocate( ds%datasetR(n,m), ds%datasetI(n,m) )
  do j=1,m
  do i=1,n
    allocate( ds%datasetR(i,j)%gxyz(ds%nx,ds%ny,ds%nz) )
    allocate( ds%datasetI(i,j)%gxyz(ds%nx,ds%ny,ds%nz) )
  enddo
  enddo

  forall( i=1:ds%nx ) &
    ds%x(i) = dble(i-1)/ngrid(1) - dble(igrid(1))/2/ngrid(1)
  forall( i=1:ds%ny ) &
    ds%y(i) = dble(i-1)/ngrid(2) - dble(igrid(2))/2/ngrid(2)
  forall( i=1:ds%nz ) &
    ds%z(i) = dble(i-1)/ngrid(3) - dble(igrid(3))/2/ngrid(3)

  ijk=(ds%nx-1)*(ds%ny-1)*(ds%nz-1)
  allocate( ds%gridindex(3,ijk) )
  ijk=0
  do k=1,ds%nz-1
    do j=1,ds%ny-1
      do i=1,ds%nx-1
        ijk=ijk+1
        ds%gridindex(1:3,ijk) = (/ i, j, k /)
      enddo
    enddo
  enddo

  ds%transcart = transcart
  ds%transgrid = transgrid

  end subroutine create_spline_dataset

! ----------------------------------------------------------------------
  function numgridcoord_spline_dataset( ds )
! ----------------------------------------------------------------------

  type(spline_dataset_type) :: ds
  integer :: numgridcoord_spline_dataset

  numgridcoord_spline_dataset = size(ds%gridindex,2)

  end function numgridcoord_spline_dataset

! ----------------------------------------------------------------------
  subroutine getgridcoord_spline_dataset( ixyz, r, ds )
! ----------------------------------------------------------------------

  integer,intent(in) :: ixyz
  real(dp) :: r(3)
  type(spline_dataset_type) :: ds

  integer :: i,j,k

  i = ds%gridindex(1,ixyz)
  j = ds%gridindex(2,ixyz)
  k = ds%gridindex(3,ixyz)

  r(1) = ds%x(i)
  r(2) = ds%y(j)
  r(3) = ds%z(k)

  end subroutine getgridcoord_spline_dataset

! ----------------------------------------------------------------------
  subroutine getcartcoord_spline_dataset( ixyz, r, ds )
! ----------------------------------------------------------------------

  integer,intent(in) :: ixyz
  real(dp) :: r(3)
  type(spline_dataset_type) :: ds

  call getgridcoord_spline_dataset( ixyz, r, ds )
  r = matmul( ds%transcart, r )

  end subroutine getcartcoord_spline_dataset

! ----------------------------------------------------------------------
  subroutine putelement_spline_dataset( i, j, ixyz, zdata, ds )
! ----------------------------------------------------------------------

  integer,intent(in) :: i,j,ixyz
  complex(dp),intent(in) :: zdata
  type(spline_dataset_type) :: ds
  integer :: ix,iy,iz

  ix = ds%gridindex(1,ixyz)
  iy = ds%gridindex(2,ixyz)
  iz = ds%gridindex(3,ixyz)

  ds%datasetR(i,j)%gxyz(ix,iy,iz) =  real(zdata)
  ds%datasetI(i,j)%gxyz(ix,iy,iz) = dimag(zdata)

  ! assuming periodicity 1==n
  if( ix==1 .and. iy==1 .and. iz==1 ) then
    ds%datasetR(i,j)%gxyz(   ix,   iy,ds%nz) =  real(zdata)
    ds%datasetI(i,j)%gxyz(   ix,   iy,ds%nz) = dimag(zdata)
    ds%datasetR(i,j)%gxyz(   ix,ds%ny,   iz) =  real(zdata)
    ds%datasetI(i,j)%gxyz(   ix,ds%ny,   iz) = dimag(zdata)
    ds%datasetR(i,j)%gxyz(   ix,ds%ny,ds%nz) =  real(zdata)
    ds%datasetI(i,j)%gxyz(   ix,ds%ny,ds%nz) = dimag(zdata)
    ds%datasetR(i,j)%gxyz(ds%nx,   iy,   iz) =  real(zdata)
    ds%datasetI(i,j)%gxyz(ds%nx,   iy,   iz) = dimag(zdata)
    ds%datasetR(i,j)%gxyz(ds%nx,   iy,ds%nz) =  real(zdata)
    ds%datasetI(i,j)%gxyz(ds%nx,   iy,ds%nz) = dimag(zdata)
    ds%datasetR(i,j)%gxyz(ds%nx,ds%ny,   iz) =  real(zdata)
    ds%datasetI(i,j)%gxyz(ds%nx,ds%ny,   iz) = dimag(zdata)
    ds%datasetR(i,j)%gxyz(ds%nx,ds%ny,ds%nz) =  real(zdata)
    ds%datasetI(i,j)%gxyz(ds%nx,ds%ny,ds%nz) = dimag(zdata)
  else if( iy==1 .and. iz==1 ) then
    ds%datasetR(i,j)%gxyz(   ix,   iy,ds%nz) =  real(zdata)
    ds%datasetI(i,j)%gxyz(   ix,   iy,ds%nz) = dimag(zdata)
    ds%datasetR(i,j)%gxyz(   ix,ds%ny,   iz) =  real(zdata)
    ds%datasetI(i,j)%gxyz(   ix,ds%ny,   iz) = dimag(zdata)
    ds%datasetR(i,j)%gxyz(   ix,ds%ny,ds%nz) =  real(zdata)
    ds%datasetI(i,j)%gxyz(   ix,ds%ny,ds%nz) = dimag(zdata)
  else if( ix==1 .and. iz==1 ) then
    ds%datasetR(i,j)%gxyz(   ix,   iy,ds%nz) =  real(zdata)
    ds%datasetI(i,j)%gxyz(   ix,   iy,ds%nz) = dimag(zdata)
    ds%datasetR(i,j)%gxyz(ds%nx,   iy,   iz) =  real(zdata)
    ds%datasetI(i,j)%gxyz(ds%nx,   iy,   iz) = dimag(zdata)
    ds%datasetR(i,j)%gxyz(ds%nx,   iy,ds%nz) =  real(zdata)
    ds%datasetI(i,j)%gxyz(ds%nx,   iy,ds%nz) = dimag(zdata)
  else if( ix==1 .and. iy==1 ) then
    ds%datasetR(i,j)%gxyz(   ix,ds%ny,   iz) =  real(zdata)
    ds%datasetI(i,j)%gxyz(   ix,ds%ny,   iz) = dimag(zdata)
    ds%datasetR(i,j)%gxyz(ds%nx,   iy,   iz) =  real(zdata)
    ds%datasetI(i,j)%gxyz(ds%nx,   iy,   iz) = dimag(zdata)
    ds%datasetR(i,j)%gxyz(ds%nx,ds%ny,   iz) =  real(zdata)
    ds%datasetI(i,j)%gxyz(ds%nx,ds%ny,   iz) = dimag(zdata)
  else if( iz==1 ) then
    ds%datasetR(i,j)%gxyz(   ix,   iy,ds%nz) =  real(zdata)
    ds%datasetI(i,j)%gxyz(   ix,   iy,ds%nz) = dimag(zdata)
  else if( iy==1 ) then
    ds%datasetR(i,j)%gxyz(   ix,ds%ny,   iz) =  real(zdata)
    ds%datasetI(i,j)%gxyz(   ix,ds%ny,   iz) = dimag(zdata)
  else if( ix==1 ) then
    ds%datasetR(i,j)%gxyz(ds%nx,   iy,   iz) =  real(zdata)
    ds%datasetI(i,j)%gxyz(ds%nx,   iy,   iz) = dimag(zdata)
  endif

  end subroutine putelement_spline_dataset

! ----------------------------------------------------------------------
  subroutine getelement_spline_dataset( i, j, ixyz, zdata, ds )
! ----------------------------------------------------------------------

  integer,intent(in) :: i,j,ixyz
  complex(dp),intent(out) :: zdata
  type(spline_dataset_type) :: ds
  integer :: ix,iy,iz

  ix = ds%gridindex(1,ixyz)
  iy = ds%gridindex(2,ixyz)
  iz = ds%gridindex(3,ixyz)

  zdata = cmplx( ds%datasetR(i,j)%gxyz(ix,iy,iz), &
                 ds%datasetI(i,j)%gxyz(ix,iy,iz) )

  end subroutine getelement_spline_dataset

! ----------------------------------------------------------------------
  subroutine dump_spline_dataset( ds )
! ----------------------------------------------------------------------
  type(spline_dataset_type) :: ds
  integer :: ix, iy, iz

  do iz=1,ds%nz
  do iy=1,ds%ny
  do ix=1,ds%nx
    write(600+iz,'(3f12.5,1(2f12.5))') ds%x(ix), ds%y(iy), ds%z(iz), &
      ds%datasetR(1,1)%gxyz(ix,iy,iz), ds%datasetI(1,1)%gxyz(ix,iy,iz)
  enddo
    write(600+iz,*)
  enddo
  enddo

  end subroutine dump_spline_dataset

! ----------------------------------------------------------------------
  subroutine destroy_spline_dataset( ds )
! ----------------------------------------------------------------------
  type(spline_dataset_type) :: ds
  integer :: i,j

  do j=1,size(ds%datasetR,2)
  do i=1,size(ds%datasetR,1)
    deallocate( ds%datasetR(i,j)%gxyz )
    deallocate( ds%datasetI(i,j)%gxyz )
  enddo
  enddo
  deallocate( ds%datasetR, ds%datasetI )

  deallocate( ds%gridindex )
  deallocate( ds%x )
  deallocate( ds%y )
  deallocate( ds%z )

  end subroutine destroy_spline_dataset

! ----------------------------------------------------------------------
  subroutine fit_spline_to_dataset( ds, fit )
! ----------------------------------------------------------------------
  type(spline_dataset_type) :: ds
  type(spline_fit_type) :: fit
  integer :: i, j

  do j=1,size(ds%datasetR,2)
  do i=1,size(ds%datasetR,1)

    call dbs3in( ds%nx, ds%x, ds%ny, ds%y, ds%nz, ds%z, &
                 ds%datasetR(i,j)%gxyz, size(ds%datasetR(i,j)%gxyz,1), &
                 size(ds%datasetR(i,j)%gxyz,2), &
                 fit%xorder, fit%yorder, fit%zorder, &
                 fit%xknot, fit%yknot, fit%zknot, fit%coefR(i,j)%bscoef )

  enddo
  enddo
  
  do j=1,size(ds%datasetI,2)
  do i=1,size(ds%datasetI,1)

    call dbs3in( ds%nx, ds%x, ds%ny, ds%y, ds%nz, ds%z, &
                 ds%datasetI(i,j)%gxyz, size(ds%datasetI(i,j)%gxyz,1), &
                 size(ds%datasetI(i,j)%gxyz,2), &
                 fit%xorder, fit%yorder, fit%zorder, &
                 fit%xknot, fit%yknot, fit%zknot, fit%coefI(i,j)%bscoef )

  enddo
  enddo

  end subroutine fit_spline_to_dataset

! ----------------------------------------------------------------------
  subroutine create_spline_fit( n, m, order, ds, fit )
! ----------------------------------------------------------------------
  integer,intent(in) :: n,m,order(3)
  type(spline_dataset_type) :: ds
  type(spline_fit_type) :: fit
  integer :: i,j

  ! number of coefficients is equal to number of data points
  fit%nxcoef = ds%nx
  fit%nycoef = ds%ny
  fit%nzcoef = ds%nz

  ! fix spline order
  fit%xorder = order(1)
  fit%yorder = order(2)
  fit%zorder = order(3)

  ! check for default order
  if( fit%xorder <= 0 .or. fit%xorder > ds%nx ) fit%xorder = ds%nx
  if( fit%yorder <= 0 .or. fit%yorder > ds%ny ) fit%yorder = ds%ny
  if( fit%zorder <= 0 .or. fit%zorder > ds%nz ) fit%zorder = ds%nz

  allocate( fit%xknot( ds%nx + fit%xorder ) )
  allocate( fit%yknot( ds%ny + fit%yorder ) )
  allocate( fit%zknot( ds%nz + fit%zorder ) )

  allocate( fit%coefR(n,m), fit%coefI(n,m) ) 
  do j=1,m
  do i=1,n
    allocate( fit%coefR(i,j)%bscoef(ds%nx*ds%ny*ds%nz) )
    allocate( fit%coefI(i,j)%bscoef(ds%nx*ds%ny*ds%nz) )
  enddo
  enddo

  ! make the knots
  call dbsnak( ds%nx, ds%x, fit%xorder, fit%xknot )
  call dbsnak( ds%ny, ds%y, fit%yorder, fit%yknot )
  call dbsnak( ds%nz, ds%z, fit%zorder, fit%zknot )

  end subroutine create_spline_fit

! ----------------------------------------------------------------------
  subroutine destroy_spline_fit( fit )
! ----------------------------------------------------------------------
  type(spline_fit_type) :: fit
  integer :: i,j

  deallocate( fit%xknot )
  deallocate( fit%yknot )
  deallocate( fit%zknot )

  do j=1,size(fit%coefR,2)
  do i=1,size(fit%coefR,1)
    deallocate( fit%coefR(i,j)%bscoef )
    deallocate( fit%coefI(i,j)%bscoef )
  enddo
  enddo
  deallocate( fit%coefR, fit%coefI ) 

  end subroutine destroy_spline_fit

! ----------------------------------------------------------------------
  subroutine write_spline_fit( iun, fit )
! ----------------------------------------------------------------------
  integer,intent(in) :: iun
  type(spline_fit_type) :: fit
  integer :: i,j

  write(iun) fit%xorder, fit%yorder, fit%zorder
  write(iun) fit%nxcoef, fit%nycoef, fit%nzcoef
  write(iun) size(fit%coefR,1), size(fit%coefR,2)
  write(iun) fit%xknot, fit%yknot, fit%zknot
  do j=1,size(fit%coefR,2)
  do i=1,size(fit%coefR,1)
    write(iun) fit%coefR(i,j)%bscoef 
  enddo
  enddo
  do j=1,size(fit%coefI,2)
  do i=1,size(fit%coefI,1)
    write(iun) fit%coefI(i,j)%bscoef 
  enddo
  enddo

  end subroutine write_spline_fit

! ----------------------------------------------------------------------
  subroutine read_spline_fit( iun, fit )
! ----------------------------------------------------------------------
  integer,intent(in) :: iun
  type(spline_fit_type) :: fit
  integer :: n,m, i,j, nxyz

  read(iun) fit%xorder, fit%yorder, fit%zorder
  read(iun) fit%nxcoef, fit%nycoef, fit%nzcoef
  read(iun) n,m
  allocate( fit%xknot(fit%nxcoef+fit%xorder) )
  allocate( fit%yknot(fit%nycoef+fit%yorder) )
  allocate( fit%zknot(fit%nzcoef+fit%zorder) )
  read(iun) fit%xknot, fit%yknot, fit%zknot
  allocate( fit%coefR(n,m), fit%coefI(n,m) )
  nxyz=fit%nxcoef*fit%nycoef*fit%nzcoef
  do j=1,m
  do i=1,n
    allocate( fit%coefR(i,j)%bscoef(nxyz) )
    read(iun) fit%coefR(i,j)%bscoef
  enddo
  enddo
  do j=1,m
  do i=1,n
    allocate( fit%coefI(i,j)%bscoef(nxyz) )
    read(iun) fit%coefI(i,j)%bscoef
  enddo
  enddo
  
  end subroutine read_spline_fit
  
! ----------------------------------------------------------------------
  subroutine bcast_spline_fit( fit, mpime, root, comm )
! ----------------------------------------------------------------------
  use mp, only : mp_bcast, mp_barrier
  integer,intent(in) :: mpime, root, comm
  type(spline_fit_type) :: fit
  integer :: n,m, i,j, nxyz

  call mp_barrier( comm )
  call mp_bcast( fit%xorder, root, comm )
  call mp_bcast( fit%yorder, root, comm )
  call mp_bcast( fit%zorder, root, comm )
  call mp_bcast( fit%nxcoef, root, comm )
  call mp_bcast( fit%nycoef, root, comm )
  call mp_bcast( fit%nzcoef, root, comm )
  
  if( mpime/=root ) then
    allocate( fit%xknot(fit%nxcoef+fit%xorder) )
    allocate( fit%yknot(fit%nycoef+fit%yorder) )
    allocate( fit%zknot(fit%nzcoef+fit%zorder) )
  endif

  call mp_bcast( fit%xknot, root, comm )
  call mp_bcast( fit%yknot, root, comm )
  call mp_bcast( fit%zknot, root, comm )

  if( mpime==root ) then
    n=size(fit%coefR,1)
    m=size(fit%coefR,2)
  endif
  call mp_bcast( n, root, comm )
  call mp_bcast( m, root, comm )

  if( mpime/=root ) then
    allocate( fit%coefR(n,m), fit%coefI(n,m) )
    nxyz=fit%nxcoef*fit%nycoef*fit%nzcoef
    do j=1,m
    do i=1,n
      allocate( fit%coefR(i,j)%bscoef(nxyz) )
      allocate( fit%coefI(i,j)%bscoef(nxyz) )
    enddo
    enddo
  endif

  do j=1,m
  do i=1,n
    call mp_bcast( fit%coefR(i,j)%bscoef, root, comm )
    call mp_bcast( fit%coefI(i,j)%bscoef, root, comm )
  enddo
  enddo
  
  end subroutine bcast_spline_fit
  
! ----------------------------------------------------------------------
  subroutine getdims_spline_fit( n, m, fit )
! ----------------------------------------------------------------------
  integer,intent(out) :: n, m
  type(spline_fit_type) :: fit

  n = size(fit%coefR,1)
  m = size(fit%coefR,2)

  end subroutine getdims_spline_fit 

! ----------------------------------------------------------------------
  subroutine evaluate_spline_fit( r, fit, evfit )
! ----------------------------------------------------------------------
  real(dp),intent(in) :: r(3)
  type(spline_fit_type) :: fit
  complex(dp) :: evfit(size(fit%coefR,1),size(fit%coefR,2))
  real(dp) :: evR, evI
  integer :: i,j

  do j=1,size(fit%coefR,2)
  do i=1,size(fit%coefR,1)
    
    evR  = spline_eval( r, fit%nxcoef, fit%nycoef, fit%nzcoef, &
                        fit%xorder, fit%yorder, fit%zorder,    &
                        fit%xknot, fit%yknot, fit%zknot, &
                        fit%coefR(i,j)%bscoef )

    evI = spline_eval( r, fit%nxcoef, fit%nycoef, fit%nzcoef, &
                       fit%xorder, fit%yorder, fit%zorder,    &
                       fit%xknot, fit%yknot, fit%zknot, &
                       fit%coefI(i,j)%bscoef )

    evfit(i,j) = cmplx( evR, evI )

  enddo
  enddo

  end subroutine evaluate_spline_fit

! ----------------------------------------------------------------------
  function spline_eval( r, nxcoef, nycoef, nzcoef, xorder, yorder, zorder, &
                        xknot, yknot, zknot, bscoef )
! ----------------------------------------------------------------------
  real(dp) :: r(3)
  integer :: nxcoef, nycoef, nzcoef
  integer :: xorder, yorder, zorder
  real(dp) :: xknot(nxcoef+xorder), yknot(nycoef+yorder), zknot(nzcoef+zorder)
  real(dp) :: bscoef(nxcoef*nycoef*nzcoef)
  real(dp) :: sout, spline_eval

  ! 0-D
  if( nxcoef==1 .and. nycoef==1 .and. nzcoef==1 ) then
    sout = bscoef(1)
  ! 1-D
  else if( nycoef==1 .and. nzcoef==1 ) then
    sout = dbsval(r(1),xorder,xknot,nxcoef,bscoef)
  else if( nxcoef==1 .and. nzcoef==1 ) then
    sout = dbsval(r(2),yorder,yknot,nycoef,bscoef)
  else if( nxcoef==1 .and. nycoef==1 ) then
    sout = dbsval(r(3),zorder,zknot,nzcoef,bscoef)
  ! 2-D
  else if( nzcoef==1 ) then
    sout = dbs2vl(r(1),r(2),xorder,yorder,xknot,yknot,nxcoef,nycoef,bscoef)
  else if( nycoef==1 ) then
    sout = dbs2vl(r(1),r(3),xorder,zorder,xknot,zknot,nxcoef,nzcoef,bscoef)
  else if( nxcoef==1 ) then
    sout = dbs2vl(r(2),r(3),yorder,zorder,yknot,zknot,nycoef,nzcoef,bscoef)
  ! 3-D
  else
    sout = dbs3vl( r(1), r(2), r(3), &
                   xorder, yorder, zorder, xknot, yknot, zknot, &
                   nxcoef, nycoef, nzcoef, bscoef )
  endif
  spline_eval = sout

  return
  end function spline_eval


  end module splines_module

