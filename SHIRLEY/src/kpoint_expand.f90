! ---------------------------------------------------------------------- 
  subroutine kpoint_expand( ndim )
! ---------------------------------------------------------------------- 

  ! This is an important mapping of the available wave functions from
  ! their original k-points to k in [0,1]
  ! Note that we include the zone boundary at 1 aswell
  ! This means for example that the Gamma-point will be mapped
  ! to itself (000) and also (001), (010), (011), (100), (101), (110), (111)
  ! If we have a system of reduced dimensions (e.g. 1D-periodic) then the
  ! possible mappings is reduced, e.g. (000) and (001) for 1D along z

  ! I need to add inversion symmetry to this routine

  use kinds,     only : dp
  use klist,     only : xk, nks, nkstot
  use cell_base, only : at
  use symm_base,     only : s, nsym
  use io_global, only : stdout

  use wfc_shirley, only : xk_map

  implicit none

  real(dp),parameter :: eps=1.d-12

  integer,intent(in) :: ndim(3)

  integer :: ik, isym
  real(dp) :: xk_latt(3,nkstot)
  integer :: imap, ikm
  integer :: ig, jg, kg, izero(3)


  write(stdout,*)
  write(stdout,*) ' kpoint_expand '
  write(stdout,*) ' Expand original k-point set in the positive BZ: [0,1]^3'
  write(stdout,*)
  write(stdout,*) ' system dimensionality: ndim = ', ndim

  write(stdout,*) ' nks, nkstot, nsym = ', nks, nkstot, nsym
  write(stdout,*) ' k-points'
  write(stdout,'(6x,8x,a8,8x,2x,8x,a8,8x)') 'k (cart)', 'k (latt)'
  do ik=1,nkstot
    xk_latt(:,ik) = matmul( transpose(at), xk(:,ik) )
    write(stdout,'(i6,3f8.4,2x,3f8.4)') ik, xk(:,ik), xk_latt(:,ik)
  enddo

  if( nsym > 1 ) then
    call errore('kpoint_expand','Number of symmetries is > 1',-1)
    write(stdout,*) ' deliberately not using rotational or inversion symmetry to expand'
    write(stdout,*) ' If you have used symmetry you may have too few k-points'
    write(stdout,*) ' Set nosym=.true. in PWSCF input'
    write(stdout,*)
  endif

  ! now pad the BZ with edges
  ! and collect the final number of mapped kpoints
  if( allocated(xk_map) ) then
    do ik=1,size(xk_map)
      deallocate( xk_map(ik)%xmap )
    enddo
    deallocate(xk_map)
  endif
  allocate( xk_map(nkstot) )
  xk_map(:)%nmap = 0
  do ik=1,nkstot
    xk_map(ik)%nmap = xk_map(ik)%nmap + 1

    izero=0
    where( abs(xk_latt(:,ik)) < eps ) izero=1
    do kg=0,min(ndim(3),izero(3))
    do jg=0,min(ndim(2),izero(2))
    do ig=0,min(ndim(1),izero(1))
      if( ig==0 .and. jg==0 .and. kg==0 ) cycle
      xk_map(ik)%nmap = xk_map(ik)%nmap + 1
    enddo
    enddo
    enddo
  enddo

  ! allocate
  do ik=1,nkstot
    allocate( xk_map(ik)%xmap( 3, xk_map(ik)%nmap ) )
  enddo

  ! store displacements
  xk_map(:)%nmap = 0
  do ik=1,nkstot
    xk_map(ik)%nmap = xk_map(ik)%nmap + 1
    xk_map(ik)%xmap(:,xk_map(ik)%nmap) = 0.d0

    izero=0
    where( abs(xk_latt(:,ik)) < eps ) izero=1
    do kg=0,min(ndim(3),izero(3))
    do jg=0,min(ndim(2),izero(2))
    do ig=0,min(ndim(1),izero(1))
      if( ig==0 .and. jg==0 .and. kg==0 ) cycle
      xk_map(ik)%nmap = xk_map(ik)%nmap + 1
      xk_map(ik)%xmap(:,xk_map(ik)%nmap) = (/ ig, jg, kg /)
    enddo
    enddo
    enddo
  enddo

  ! report
  write(stdout,*) ' expanded set of k-points: '
  write(stdout,*) '   # symmetry operations = ', nsym
  write(stdout,*) ' this code excludes rotational and inversion symmetry'
  write(stdout,*) '   # original kpts = ', nkstot
  write(stdout,*) '   # expanded kpts = ', sum(xk_map(:)%nmap)
  write(stdout,'(2(2x,6x,12x,a12,12x))') '  expanded k', '  original k'
  ikm=0
  do ik=1,nkstot
    write(stdout,*)
    do imap=1,xk_map(ik)%nmap
      ikm=ikm+1
      write(stdout,'(2(2x,i6,3f12.6))') &
        ikm, xk_latt(:,ik)+xk_map(ik)%xmap(:,imap), &
        ik, xk_latt(:,ik)
    enddo
  enddo
  write(stdout,*)
 
  return
  
  end subroutine kpoint_expand

