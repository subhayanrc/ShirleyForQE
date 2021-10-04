  module mpio

  use parallel_include
  use kinds, only : dp

  implicit none

  public

  interface mp_file_scatter_write
    module procedure mp_file_scatter_write_r1
    module procedure mp_file_scatter_write_r3
    module procedure mp_file_scatter_write_c1
    module procedure mp_file_scatter_write_c2
  end interface mp_file_scatter_write

  interface mp_file_scatter_read
    module procedure mp_file_scatter_read_r3
    module procedure mp_file_scatter_read_c1
  end interface mp_file_scatter_read

  contains

  subroutine mp_file_open_dp( filename, fhandle, root, cm )

!  use mp, only : mp_bcast

  character(*) :: filename
  integer :: root
  ! specifying optional here leads to seg faults when I don't include cm
  !integer,optional :: cm
  integer :: cm
  integer :: fhandle

  integer :: comm, fmode, finfo, ierr
  integer(kind=MPI_OFFSET_KIND) :: disp

  integer :: errorcode, resultlen
  character(255) :: string

  comm=MPI_COMM_WORLD
  ! this caused problems
  !if( present(cm) ) comm=cm
  comm=cm

  fmode=IOR(MPI_MODE_RDWR,MPI_MODE_CREATE)

  call MPI_INFO_CREATE( finfo, ierr )
  if( ierr/=0 ) then
    errorcode=ierr
    call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
    call errore(string,errorcode)
  endif

  call mpi_barrier( comm, ierr )
  call MPI_FILE_OPEN( comm, filename, fmode, finfo, fhandle, ierr )
  if( ierr/=0 ) then
    errorcode=ierr
    call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
    call errore(string,errorcode)
  endif
  call mpi_barrier( comm, ierr )

  disp=0
  call MPI_FILE_SET_VIEW( fhandle, disp, MPI_DOUBLE_PRECISION, &
                          MPI_DOUBLE_PRECISION, 'native', finfo, ierr )
  if( ierr/=0 ) then
    errorcode=ierr
    call MPI_ERROR_STRING( errorcode, string, resultlen, ierr )
    call errore(string,errorcode)
  endif

  return
  end subroutine mp_file_open_dp


  subroutine mp_wait( request, status, ierr )

  integer,intent(in) :: request
  integer :: status(MPI_STATUS_SIZE)
  integer,intent(out) :: ierr
  
#ifdef __MPICH
  call MPIO_WAIT(request, status, ierr)
#else
  call MPI_WAIT(request, status, ierr)
#endif
  end subroutine mp_wait


  subroutine mp_file_scatter_write_r1( fhandle, fposn, msg, root, cm )

  use mp, only : mp_bcast, mp_stop

  integer,intent(in) :: fhandle
  integer(kind=MPI_OFFSET_KIND),intent(inout) :: fposn
  real(dp),intent(in) :: msg(:)
  integer,intent(in) :: root, cm

  integer :: status(MPI_STATUS_SIZE)
  integer :: msglen, group, numtask, taskid, ierr
  integer(kind=MPI_OFFSET_KIND) :: i
  integer,allocatable :: lengths(:)
  integer(kind=MPI_OFFSET_KIND),allocatable :: displs(:)
  integer(kind=MPI_OFFSET_KIND) :: offset

  ! message length
  msglen = SIZE(msg)
  ! group
  group = mpi_comm_world
  group = cm

  call mpi_comm_size(group,numtask,ierr)
  if( ierr/=0 ) call mp_stop(8400)
  call mpi_comm_rank(group,taskid,ierr)
  if( ierr/=0 ) call mp_stop(8401)

  if( taskid == root ) then
    allocate( lengths(numtask), displs(numtask) )
  endif
  call mpi_gather( size(msg), 1, MPI_INTEGER, &
                   lengths, 1, MPI_INTEGER, &
                   root, group, ierr )
  if( ierr/=0 ) call mp_stop(8402)
  if( taskid == root ) then
    displs(1)=fposn
    do i=2,numtask
      displs(i)=displs(i-1)+lengths(i-1)
    enddo
    ! update file position
    fposn=displs(numtask)+lengths(numtask)
  endif
  call mpi_scatter( displs, 1, MPI_INTEGER8, &
                    i, 1, MPI_INTEGER8, &
                    root, group, ierr )
  if( ierr/=0 ) call mp_stop(8403)
  ! here we convert between different kinds of integer
  offset=i
  call mp_bcast( fposn, root, group )
  if( taskid == root ) then
    deallocate( lengths, displs )
  endif

  call mpi_file_write_at( fhandle, offset, &
                          msg, msglen, &
                          MPI_DOUBLE_PRECISION, status, ierr )
  if( ierr/=0 ) call mp_stop(8404)

  return
  end subroutine mp_file_scatter_write_r1


  subroutine mp_file_scatter_write_r3( fhandle, fposn, msg, root, cm )

  use mp, only : mp_bcast, mp_stop

  integer,intent(in) :: fhandle
  integer(kind=MPI_OFFSET_KIND),intent(inout) :: fposn
  real(dp),intent(in) :: msg(:,:,:)
  integer,intent(in) :: root, cm

  integer :: status(MPI_STATUS_SIZE)
  integer :: msglen, group, numtask, taskid, ierr
  integer(kind=MPI_OFFSET_KIND) :: i
  integer,allocatable :: lengths(:)
  integer(kind=MPI_OFFSET_KIND),allocatable :: displs(:)
  integer(kind=MPI_OFFSET_KIND) :: offset

  ! message length
  msglen = SIZE(msg)
  ! group
  group = mpi_comm_world
  group = cm

  call mpi_comm_size(group,numtask,ierr)
  if( ierr/=0 ) call mp_stop(8400)
  call mpi_comm_rank(group,taskid,ierr)
  if( ierr/=0 ) call mp_stop(8401)

  if( taskid == root ) then
    allocate( lengths(numtask), displs(numtask) )
  endif
  call mpi_gather( size(msg), 1, MPI_INTEGER, &
                   lengths, 1, MPI_INTEGER, &
                   root, group, ierr )
  if( ierr/=0 ) call mp_stop(8402)
  if( taskid == root ) then
    displs(1)=fposn
    do i=2,numtask
      displs(i)=displs(i-1)+lengths(i-1)
    enddo
    ! update file position
    fposn=displs(numtask)+lengths(numtask)
  endif
  call mpi_scatter( displs, 1, MPI_INTEGER8, &
                    i, 1, MPI_INTEGER8, &
                    root, group, ierr )
  if( ierr/=0 ) call mp_stop(8403)
  ! here we convert between different kinds of integer
  offset=i
  call mp_bcast( fposn, root, group )
  if( taskid == root ) then
    deallocate( lengths, displs )
  endif

  call mpi_file_write_at( fhandle, offset, &
                          msg, msglen, &
                          MPI_DOUBLE_PRECISION, status, ierr )
  if( ierr/=0 ) call mp_stop(8404)

  return
  end subroutine mp_file_scatter_write_r3


  subroutine mp_file_scatter_write_c1( fhandle, fposn, msg, root, cm )

  use mp, only : mp_bcast, mp_stop

  integer,intent(in) :: fhandle
  integer(kind=MPI_OFFSET_KIND),intent(inout) :: fposn
  complex(dp),intent(in) :: msg(:)
  integer,intent(in) :: root, cm

  integer :: status(MPI_STATUS_SIZE)
  integer :: msglen, group, numtask, taskid, ierr
  integer(kind=MPI_OFFSET_KIND) :: i
  integer,allocatable :: lengths(:)
  integer(kind=MPI_OFFSET_KIND),allocatable :: displs(:)
  integer(kind=MPI_OFFSET_KIND) :: offset

  ! message length (in real units)
  msglen = SIZE(msg)*2
  ! group
  group = mpi_comm_world
  group = cm

  call mpi_comm_size(group,numtask,ierr)
  if( ierr/=0 ) call mp_stop(8400)
  call mpi_comm_rank(group,taskid,ierr)
  if( ierr/=0 ) call mp_stop(8401)

  if( taskid == root ) then
    allocate( lengths(numtask), displs(numtask) )
  endif
  call mpi_gather( msglen, 1, MPI_INTEGER, &
                   lengths, 1, MPI_INTEGER, &
                   root, group, ierr )
  if( ierr/=0 ) call mp_stop(8402)
  if( taskid == root ) then
    displs(1)=fposn
    do i=2,numtask
      displs(i)=displs(i-1)+lengths(i-1)
    enddo
    ! update file position
    fposn=displs(numtask)+lengths(numtask)
  endif
  call mpi_scatter( displs, 1, MPI_INTEGER8, &
                    i, 1, MPI_INTEGER8, &
                    root, group, ierr )
  if( ierr/=0 ) call mp_stop(8403)
  call mp_bcast( fposn, root, group )
  ! here we convert between different kinds of integer
  offset=i
  if( taskid == root ) then
    deallocate( lengths, displs )
  endif

  call mpi_file_write_at( fhandle, offset, &
                          msg, msglen, &
                          MPI_DOUBLE_PRECISION, status, ierr )
  if( ierr/=0 ) call mp_stop(8404)

  return
  end subroutine mp_file_scatter_write_c1


  subroutine mp_file_scatter_write_c2( fhandle, fposn, msg, root, cm )

  use mp, only : mp_bcast, mp_stop

  integer,intent(in) :: fhandle
  integer(kind=MPI_OFFSET_KIND),intent(inout) :: fposn
  complex(dp),intent(in) :: msg(:,:)
  integer,intent(in) :: root, cm

  integer :: status(MPI_STATUS_SIZE)
  integer :: msglen, group, numtask, taskid, ierr
  integer(kind=MPI_OFFSET_KIND) :: i
  integer,allocatable :: lengths(:)
  integer(kind=MPI_OFFSET_KIND),allocatable :: displs(:)
  integer(kind=MPI_OFFSET_KIND) :: offset

  ! message length (in real units)
  msglen = SIZE(msg)*2
  ! group
  group = mpi_comm_world
  group = cm

  call mpi_comm_size(group,numtask,ierr)
  if( ierr/=0 ) call mp_stop(8400)
  call mpi_comm_rank(group,taskid,ierr)
  if( ierr/=0 ) call mp_stop(8401)

  if( taskid == root ) then
    allocate( lengths(numtask), displs(numtask) )
  endif
  call mpi_gather( msglen, 1, MPI_INTEGER, &
                   lengths, 1, MPI_INTEGER, &
                   root, group, ierr )
  if( ierr/=0 ) call mp_stop(8402)
  if( taskid == root ) then
    displs(1)=fposn
    do i=2,numtask
      displs(i)=displs(i-1)+lengths(i-1)
    enddo
    ! update file position
    fposn=displs(numtask)+lengths(numtask)
  endif
  call mpi_scatter( displs, 1, MPI_INTEGER8, &
                    i, 1, MPI_INTEGER8, &
                    root, group, ierr )
  if( ierr/=0 ) call mp_stop(8403)
  call mp_bcast( fposn, root, group )
  ! here we convert between different kinds of integer
  offset=i
  if( taskid == root ) then
    deallocate( lengths, displs )
  endif

  call mpi_file_write_at( fhandle, offset, &
                          msg, msglen, &
                          MPI_DOUBLE_PRECISION, status, ierr )
  if( ierr/=0 ) call mp_stop(8404)

  return
  end subroutine mp_file_scatter_write_c2


  subroutine mp_file_scatter_read_r3( fhandle, fposn, msg, root, cm )

  use mp, only : mp_bcast, mp_stop

  integer,intent(in) :: fhandle
  integer(kind=MPI_OFFSET_KIND),intent(inout) :: fposn
  real(dp),intent(out) :: msg(:,:,:)
  integer,intent(in) :: root, cm

  integer :: status(MPI_STATUS_SIZE)
  integer :: msglen, group, numtask, taskid, ierr
  integer(kind=MPI_OFFSET_KIND) :: i
  integer,allocatable :: lengths(:)
  integer(kind=MPI_OFFSET_KIND),allocatable :: displs(:)
  integer(kind=MPI_OFFSET_KIND) :: offset

  ! message length (in real units)
  msglen = SIZE(msg)
  ! group
  group = mpi_comm_world
  group = cm

  call mpi_comm_size(group,numtask,ierr)
  if( ierr/=0 ) call mp_stop(8400)
  call mpi_comm_rank(group,taskid,ierr)
  if( ierr/=0 ) call mp_stop(8401)

  !if( taskid == root ) then
    allocate( lengths(numtask), displs(numtask) )
  !endif
  call mpi_gather( msglen, 1, MPI_INTEGER, &
                   lengths, 1, MPI_INTEGER, &
                   root, group, ierr )
  if( ierr/=0 ) call mp_stop(8402)
  if( taskid == root ) then
    displs(1)=fposn
    do i=2,numtask
      displs(i)=displs(i-1)+lengths(i-1)
    enddo
    ! update file position
    fposn=displs(numtask)+lengths(numtask)
  endif
  call mpi_scatter( displs, 1, MPI_INTEGER8, &
                    i, 1, MPI_INTEGER8, &
                    root, group, ierr )
  if( ierr/=0 ) call mp_stop(8403)
  call mp_bcast( fposn, root, group )
  ! here we convert between different kinds of integer
  offset=i
  !if( taskid == root ) then
    deallocate( lengths, displs )
  !endif

  call mpi_file_read_at( fhandle, offset, &
                         msg, msglen, &
                         MPI_DOUBLE_PRECISION, status, ierr )
  if( ierr/=0 ) call mp_stop(8404)

  return
  end subroutine mp_file_scatter_read_r3


  subroutine mp_file_scatter_read_c1( fhandle, fposn, msg, root, cm )

  use mp, only : mp_bcast, mp_stop

  integer,intent(in) :: fhandle
  integer(kind=MPI_OFFSET_KIND),intent(inout) :: fposn
  complex(dp),intent(out) :: msg(:)
  integer,intent(in) :: root, cm

  integer :: status(MPI_STATUS_SIZE)
  integer :: msglen, group, numtask, taskid, ierr
  integer(kind=MPI_OFFSET_KIND) :: i
  integer,allocatable :: lengths(:)
  integer(kind=MPI_OFFSET_KIND),allocatable :: displs(:)
  integer(kind=MPI_OFFSET_KIND) :: offset

  ! message length (in real units)
  msglen = SIZE(msg)*2
  ! group
  group = mpi_comm_world
  group = cm

  call mpi_comm_size(group,numtask,ierr)
  if( ierr/=0 ) call mp_stop(8400)
  call mpi_comm_rank(group,taskid,ierr)
  if( ierr/=0 ) call mp_stop(8401)

  !if( taskid == root ) then
    allocate( lengths(numtask), displs(numtask) )
  !endif
  call mpi_gather( msglen, 1, MPI_INTEGER, &
                   lengths, 1, MPI_INTEGER, &
                   root, group, ierr )
  if( ierr/=0 ) call mp_stop(8402)
  if( taskid == root ) then
    displs(1)=fposn
    do i=2,numtask
      displs(i)=displs(i-1)+lengths(i-1)
    enddo
    ! update file position
    fposn=displs(numtask)+lengths(numtask)
  endif
  call mpi_scatter( displs, 1, MPI_INTEGER8, &
                    i, 1, MPI_INTEGER8, &
                    root, group, ierr )
  if( ierr/=0 ) call mp_stop(8403)
  call mp_bcast( fposn, root, group )
  ! here we convert between different kinds of integer
  offset=i
  !if( taskid == root ) then
    deallocate( lengths, displs )
  !endif

  call mpi_file_read_at( fhandle, offset, &
                         msg, msglen, &
                         MPI_DOUBLE_PRECISION, status, ierr )
  if( ierr/=0 ) call mp_stop(8404)

  return
  end subroutine mp_file_scatter_read_c1


  end module mpio

