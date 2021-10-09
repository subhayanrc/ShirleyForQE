  program xas_para

!#include "f_defs.h"
  use parallel_include
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_world, only : nproc, mpime, world_comm, mp_world_end
  use mp, only : mp_bcast, mp_barrier, mp_sum
  use mpio

  implicit none

  real(dp),parameter :: evtory=7.3498649d-2

  character(len=3) :: nodenumber

  integer :: narg

  integer :: nener
  real(dp) :: e1, e2, sigma, delta, efermi, dilation, evector(3)
  logical :: have_efermi, have_dilation, have_evector
  real(dp),allocatable :: ener(:), spec(:,:,:)
  complex(dp),allocatable :: xas(:,:,:)
  real(dp) :: xas_xyz(7)
  real(dp) :: xas_tensor(3,3)
  real(dp) :: spec_spin(14)

  integer :: i, j, k, n, m, ik, ispin
  real(dp) :: de
  real(dp),allocatable :: eigval(:)

  real(dp) :: ei, enorm

  logical :: cartesian
  character(255) :: grid_type
  real(dp),allocatable :: kvec(:,:), wk(:)
  real(dp) :: wktot, wkfac

  character(255) :: filename
  character(255) :: ce1
  character(255) :: ce2
  character(255) :: cnener
  character(255) :: csigma
  character(255) :: cdelta
  character(255) :: cefermi
  character(255) :: cdilation
  character(255) :: cevector(3)

  integer :: iunout
  character(255) :: fout
  integer :: iunstick
  character(255) :: stick_file, fmtstr
  integer :: iuninf, ierr
  integer :: fheigval, fhxmat
  character(255) :: eigval_file, xmat_file
  integer :: status(MPI_STATUS_SIZE)
  integer(kind=MPI_OFFSET_KIND) :: offset
  integer :: reqeigval, reqxmat

  integer,external :: freeunit
#ifdef __PGI
  integer,external :: iargc
#else
  integer,intrinsic :: iargc
#endif

  integer :: nk, nbnd, ncp, nbasis
  real(dp) :: nelec, alat, volume, &
              at(3,3), bg(3,3), tpiba, fermi_energy
  integer :: nspin
  logical :: lda_plus_u
  namelist /info/ nk, nbnd, ncp, nelec, nbasis, alat, volume, &
                  at, bg, tpiba, fermi_energy, nspin, lda_plus_u

  
  ! initialize mpi
  CALL start_shirley (nodenumber)

  if( ionode ) then
    narg = iargc()
  endif
  call mp_bcast( narg, ionode_id, world_comm )

  if( narg < 6 .or. narg > 11 ) then
    write(stdout,*) ' usage: xas_para e1 e2 nener sigma delta [efermi [dilation [Ex Ey Ez]]] filename'
    call mp_world_end
    stop
  endif

  if( ionode ) then

  have_efermi=.false.
  have_dilation=.false.
  have_evector=.false.

  call getarg( 1, ce1 )
  call getarg( 2, ce2 )
  call getarg( 3, cnener )
  call getarg( 4, csigma )
  call getarg( 5, cdelta )
  if( narg>=7 ) then
    have_efermi=.true.
    call getarg( 6, cefermi )
    if( narg>=8 ) then
      have_dilation=.true.
      call getarg( 7, cdilation )
      if( narg>=11 ) then
        have_evector=.true.
        call getarg( 8, cevector(1) )
        call getarg( 9, cevector(2) )
        call getarg( 10, cevector(3) )
      endif
    endif
  endif
  ! final argument should be the filename
  call getarg( narg, filename )

  ! default values
  dilation=1.d0
  evector=(/ 1.d0, 1.d0, 1.d0 /)

  read(ce1,*) e1
  read(ce2,*) e2
  read(cnener,*) nener
  read(csigma,*) sigma
  read(cdelta,*) delta
  if( have_efermi ) then
    read(cefermi,*) efermi
    if( have_dilation ) read(cdilation,*) dilation
    if( have_evector ) then
      read(cevector(1),*) evector(1)
      read(cevector(2),*) evector(2)
      read(cevector(3),*) evector(3)
    endif
  endif

  if( have_efermi ) then
    write(stdout,*) ' Excluding states with eigenenergies below efermi = ', efermi, ' eV'
    write(stdout,*) ' Dilation of energy scale = ', dilation
    if( have_evector ) then
      write(stdout,*) ' Electric field vector (will be normalized) = ( ', evector, ' )'
    else
      write(stdout,*) ' Assumed orientational averaging of all polarizations of sample'
    endif
  endif

  e1 = e1*evtory
  e2 = e2*evtory
  sigma = sigma*evtory
  delta = delta*evtory
  efermi = efermi*evtory

  ! normalize evector
  enorm = dot_product ( evector, evector )
  evector = evector / sqrt(enorm)

  endif

  call mp_bcast( e1, ionode_id, world_comm )
  call mp_bcast( e2, ionode_id, world_comm )
  call mp_bcast( nener, ionode_id, world_comm )
  call mp_bcast( sigma, ionode_id, world_comm )
  call mp_bcast( delta, ionode_id, world_comm )
  call mp_bcast( have_efermi, ionode_id, world_comm )
  call mp_bcast( efermi, ionode_id, world_comm )
  call mp_bcast( dilation, ionode_id, world_comm )
  call mp_bcast( evector, ionode_id, world_comm )
  call mp_bcast( filename, ionode_id, world_comm )

  if( ionode ) then
    iunout=freeunit()
    fout=trim(filename)//'.xas'
    open(iunout,file=trim(fout),form='formatted')
    write(stdout,*) '    output in '//trim(fout)

    iuninf=freeunit()
    open(iuninf,file=trim(filename)//'.info',form='formatted')
    read(iuninf,nml=info)
    allocate( wk(nk), kvec(1:3,nk) )
    read(iuninf,*) cartesian
    read(iuninf,*) grid_type
    read(iuninf,*) wk
    read(iuninf,*) kvec
    close(iuninf)
    write(stdout,*) ' reading info from ', trim(filename)
    write(stdout,*) ' fermi_energy (stored) = ', fermi_energy
    write(stdout,*) ' nbnd = ', nbnd
    write(stdout,*) ' ncp = ', ncp
  endif

  call mp_bcast( nk, ionode_id, world_comm )
  call mp_bcast( nbnd, ionode_id, world_comm )
  call mp_bcast( ncp, ionode_id, world_comm )
  call mp_bcast( nelec, ionode_id, world_comm )
  call mp_bcast( alat, ionode_id, world_comm )
  call mp_bcast( volume, ionode_id, world_comm )
  call mp_bcast( at, ionode_id, world_comm )
  call mp_bcast( bg, ionode_id, world_comm )
  call mp_bcast( tpiba, ionode_id, world_comm )
  call mp_bcast( fermi_energy, ionode_id, world_comm )
  call mp_bcast( nspin, ionode_id, world_comm )

  if( .not. ionode ) allocate( wk(nk), kvec(3,nk) )
  call mp_bcast( wk, ionode_id, world_comm )
  call mp_bcast( cartesian, ionode_id, world_comm )
  call mp_bcast( kvec, ionode_id, world_comm )


  de = (e2-e1)/dble(nener)
  allocate( ener(nener), spec(nener,7,nspin) )
  do i=1,nener
    ener(i) = e1 + dble(i-1)*de
  enddo
  spec(1:nener,1:7,1:nspin)=0.0d0


  wktot = sum(wk)
  wkfac=2.d0/wktot

  ! allocate
  allocate( eigval(nbnd) )
  allocate( xas(nbnd,ncp,3) )

  ! MPI-IO
  eigval_file=trim(filename)//'.eigval'
  xmat_file=trim(filename)//'.xmat'

  call mp_file_open_dp( eigval_file, fheigval, ionode_id, world_comm )
  call mp_file_open_dp( xmat_file, fhxmat, ionode_id, world_comm )


  iunstick=freeunit()
  stick_file=trim(filename)//'.stick'
  if( trim(nodenumber) /= '' ) then
    stick_file = trim(stick_file) // '.' // adjustl(nodenumber)
  endif
  open(iunstick,file=trim(stick_file),form='formatted')
  

  write(stdout,*) '    running...'

  do ispin=1,nspin
  do ik=1,nk
    if( mod(ik-1,nproc)/=mpime ) cycle

    write(stdout,*) ' reading ik= ', ik, ' of ', nk

    ! read eigenvalues
    offset = ((ispin-1)*nk + ik-1)*nbnd
    call mpi_file_read_at( fheigval, offset, &
                           eigval, nbnd, &
                           MPI_DOUBLE_PRECISION, status, ierr )

    ! read xas
    offset = ((ispin-1)*nk + ik-1)*nbnd*ncp*3*2
    call mpi_file_read_at( fhxmat, offset, &
                           xas, nbnd*ncp*3*2, &
                           MPI_DOUBLE_PRECISION, status, ierr )

    write(fmtstr,'(a,i6,a)') '(2i8,e20.12,e20.12,',6,'e20.12)'
    do i=1,nbnd
      if( have_efermi .and. eigval(i)<(efermi-sigma) ) cycle
      ! this line might need modification for L2/L3 branching ratio

      ! x*x, y*y, z*z
      forall(j=1:3) xas_xyz(j) = sum(conjg(xas(i,1:ncp,j))*xas(i,1:ncp,j))*wk(ik)
      ! x*y
      xas_xyz(4) = sum(conjg(xas(i,1:ncp,1))*xas(i,1:ncp,2))*wk(ik)
      ! x*z
      xas_xyz(5) = sum(conjg(xas(i,1:ncp,1))*xas(i,1:ncp,3))*wk(ik)
      ! y*z
      xas_xyz(6) = sum(conjg(xas(i,1:ncp,2))*xas(i,1:ncp,3))*wk(ik)

      ei = ((eigval(i)-efermi)*dilation+efermi)+delta

      forall(j=1:3) xas_tensor(j,j)=xas_xyz(j)
      xas_tensor(1,2) = xas_xyz(4)
      xas_tensor(2,1) = xas_xyz(4)
      xas_tensor(1,3) = xas_xyz(5)
      xas_tensor(3,1) = xas_xyz(5)
      xas_tensor(2,3) = xas_xyz(6)
      xas_tensor(3,2) = xas_xyz(6)

      if( have_evector ) then
        xas_xyz(7) = dot_product( evector, matmul( xas_tensor, evector ) )
      else
        ! make the orientational average
        xas_xyz(7) = (xas_tensor(1,1)+xas_tensor(2,2)+xas_tensor(3,3))/3.d0
      endif

      write(iunstick,trim(fmtstr)) i, ik, ei/evtory, xas_xyz(7), xas_xyz(1:6)

      ! I may want to split the spectrum over spin-channels
      ! polarized or averaged spectrum
      call add_gauss( nener, ener, ei, &
                      sigma, xas_xyz(7:7), spec(1:nener,1:1,ispin) )
      ! tensor components (xx, xy, xz, etc)
      call add_gauss( nener, ener, ei, &
                      sigma, xas_xyz(1:6), spec(1:nener,2:7,ispin) )

    enddo
  enddo
  enddo

  ! close 
  call mpi_file_close( fheigval, ierr )
  call mpi_file_close( fhxmat, ierr )

  close(iunstick)
  
  if( ionode ) write(stdout,*) '    waiting for other processors...'
  call mp_barrier( world_comm )

  call mp_sum( spec, world_comm )
  if( ionode ) then
    write(iunout,'(a,f12.5,a)') '# Applied a delta shift of ', delta/evtory, ' eV'
    if( nspin==1 ) then
      do i=1,nener
        write(iunout,'(8e14.5e3)') ener(i)/evtory,spec(i,1:size(spec,2),1)
      enddo
    else
      do i=1,nener
        do j=1,size(spec,2)
          spec_spin(2*j-1:2*j)=(/ spec(i,j,1), spec(i,j,2) /)
        enddo
        write(iunout,'(15e14.5e3)') ener(i)/evtory,spec_spin(1:2*size(spec,2))
      enddo
    endif
    close(iunout)
  endif

  call mp_barrier( world_comm ) 
  call mp_world_end
  stop

  contains

    subroutine add_gauss( nener, ener, e, sigma, weight, spec )

    integer :: nener
    real(dp) :: ener(nener), e, sigma, weight(:), spec(:,:)

    real(dp) :: its2, pre(size(weight))
    real(dp) :: arg(nener)
    integer :: i, j

    its2=2.d0*sigma*sigma
    pre=1.d0/sqrt(acos(-1.d0)*its2) * weight
    its2 = 1.d0/its2
    arg = ener - e
    arg = arg ** 2.d0
    arg = - arg * its2
    forall( i=1:nener, j=1:size(weight) ) &
      spec(i,j) = spec(i,j) + pre(j) * exp( arg(i) )
     
    end subroutine add_gauss

  end program xas_para

