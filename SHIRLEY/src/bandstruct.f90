  program bandstruct

!#include "f_defs.h"
  use parallel_include
  USE io_global,  ONLY : stdout, ionode, ionode_id
  use mp_world, only : nproc, mpime, world_comm, mp_world_end
  use mp, only : mp_barrier
  use mpio

  implicit none

  real(dp),parameter :: evtory=7.3498649d-2
  REAL(DP), PARAMETER :: rytoev=13.6058d0

  character(len=3) :: nodenumber

  integer :: narg
  character(255) :: ic

  integer :: nener
  real(dp) :: e1, e2, sigma, delta
  real(dp) :: efermi=0.d0
  logical :: have_efermi
  real(dp),allocatable :: ener(:), spec(:,:)
  complex(dp),allocatable :: xas(:,:,:)
  real(dp) :: xas_xyz(3)

  integer :: i, j, k, n, m, ik, ispin
  real(dp) :: de
  real(dp),allocatable :: eigval(:)

  logical :: cartesian
  character(255) :: grid_type
  real(dp),allocatable :: kvec(:,:), wk(:)
  real(dp) :: wktot, wkfac
  integer :: nksp
  real(dp),allocatable :: kpathlen(:)
  character(255),allocatable :: labelsp(:)
  real(dp),allocatable :: kpathlensp(:)
  real(dp) :: dk(3)

  character(255) :: filename
  character(255) :: ce1
  character(255) :: ce2
  character(255) :: cnener
  character(255) :: csigma
  character(255) :: cdelta
  character(255) :: cefermi

  integer :: iunout, iunplt
  character(255) :: fout, fplt
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

  integer :: nk, nbnd, ncp, nspin
  real(dp) :: nelec, alat, volume, &
              at(3,3), bg(3,3), tpiba, fermi_energy
  logical :: lda_plus_u
  namelist /info/ nk, nbnd, ncp, nelec, alat, volume, &
                  at, bg, tpiba, fermi_energy, nspin, lda_plus_u
  

  ! initialize mpi
  CALL start_shirley (nodenumber)
  
  if( nproc>1 ) then
    write(stdout,*) ' you should really only run this on one processor'
  endif

  if( mpime==0 ) then

  narg = iargc()
  if( narg /= 1 .and. narg /= 2 ) then
    write(stdout,*) ' usage: bandstruct [efermi] filename'
    stop
  endif

  if( narg==2 ) then
    have_efermi=.true.
    call getarg( 1, cefermi )
    call getarg( 2, filename )
  else
    have_efermi=.false.
    call getarg( 1, filename )
  endif

  if( have_efermi ) then
    read(cefermi,*) efermi
  endif

  if( have_efermi ) then
    write(stdout,*) ' Shifting energies to align the Fermi-level at zero: efermi = ', efermi, ' eV'
  endif

  efermi = efermi*evtory

  iunout=freeunit()
  fout=trim(filename)//'.bandstruct'
  open(iunout,file=trim(fout),form='formatted')
  write(stdout,*) '    output in '//trim(fout)

  iunplt=freeunit()
  fplt=trim(filename)//'.gnuplot'
  open(iunplt,file=trim(fplt),form='formatted')
  write(stdout,*) '    gnuplot in '//trim(fplt)

  iuninf=freeunit()
  open(iuninf,file=trim(filename)//'.info',form='formatted')
  read(iuninf,nml=info)
  allocate( wk(nk), kvec(1:3,nk) )
  read(iuninf,*) cartesian
  read(iuninf,*) grid_type
  read(iuninf,*) wk
  read(iuninf,*) kvec
  if( trim(grid_type) == 'bandstructure' ) then
    read(iuninf,*) nksp
    allocate( kpathlensp(nksp), labelsp(nksp) )
    do i=1,nksp
      read(iuninf,*) kpathlensp(i), &
                     labelsp(i)
    enddo
    allocate( kpathlen(nk) )
    read(iuninf,*) kpathlen
  else
    allocate( kpathlen(nk) )
    kpathlen(1)=0
    do ik=2,nk
      dk=kvec(:,ik)-kvec(:,ik-1)
      kpathlen(ik) =  vlength( dk, bg )
    enddo
  endif
  close(iuninf)
  write(stdout,*) ' reading info from ', trim(filename)
  write(stdout,*) ' fermi_energy (stored) = ', fermi_energy
  if( .not. have_efermi ) efermi = fermi_energy

  wktot = sum(wk)
  wkfac=2.d0/wktot

  ! allocate
  allocate( eigval(nbnd) )

  ! MPI-IO
  eigval_file=trim(filename)//'.eigval'

  call mp_file_open_dp( eigval_file, fheigval, ionode_id, world_comm )

  write(stdout,*) '    running...'

  do ispin=1,nspin
  do ik=1,nk
    write(stdout,*) ' reading ik= ', ik, ' of ', nk

    ! read eigenvalues
    offset = ((ispin-1)*nk + ik-1)*nbnd
    call mpi_file_read_at( fheigval, offset, &
                           eigval, nbnd, &
                           MPI_DOUBLE_PRECISION, status, ierr )

    write(fmtstr,'(a,i8,a)') '(i8,4e20.12,',nbnd,'e20.12)'
    write(iunout,fmtstr) ik, kvec(1:3,ik)/tpiba, &
                         kpathlen(ik), &
                         eigval(:)*rytoev
  enddo
    write(iunout,*)
  enddo

  if( trim(grid_type) == 'bandstructure' ) then
    write(iunplt,*) "set nokey"
    write(iunplt,*) "set style data l"
    write(iunplt,*) "bandcolor=-1"
    write(iunplt,*) "efermi=", efermi
    write(iunplt,*) "# Uncomment the lines below here for eps output"
    write(iunplt,*) "# set term postscript eps enhanced"
    write(iunplt,*) "# set output 'qdiag.eps'"
  
    write(iunplt,*) "set xtics ( \"
    do i=1,nksp-1
      write(iunplt,'(a,f20.12,a)') '"'//trim(labelsp(i))//'"', kpathlensp(i), ", \"
    enddo
    write(iunplt,'(a,f20.12,a)') '"'//trim(labelsp(nksp))//'"', kpathlensp(nksp), ")"
    do i=2,nksp-1
      write(iunplt,'(a,i8,a,f20.12,a,f20.12,a)') "set arrow ", i, &
        " from first ", kpathlensp(i), ", graph 0 to ", kpathlensp(i), &
        ", graph 1 lt bandcolor lw 1 nohead" 
    enddo

    write(iunplt,'(a, a, a)') "plot '", trim(fout), "' u 5:($6-efermi) lt bandcolor \"
    do i=2,nbnd-1
      write(ic,'(i8)') i+5 ; ic=adjustl(ic)
      write(iunplt,'(a, a, a)') ", '' u 5:($", trim(ic), "-efermi) lt bandcolor \"
    enddo
    write(ic,'(i8)') nbnd+5 ; ic=adjustl(ic)
    write(iunplt,'(a, a, a)') ", '' u 5:($", trim(ic), "-efermi) lt bandcolor"
  endif

  ! close 
  call mpi_file_close( fheigval, ierr )

  endif

  call mp_barrier( world_comm ) 
  call mp_world_end
  stop

  contains

  function vlength( r, a )
  real(dp) :: vlength
  real(dp),intent(in) :: r(3)
  real(dp),intent(in),optional :: a(3,3)
  real(dp) :: rp(3)

  if( present(a) ) then
    rp = matmul( r, a )
  else
    rp = r
  endif
  vlength = sum( rp(:)*rp(:) )
  end function vlength

  end program bandstruct

