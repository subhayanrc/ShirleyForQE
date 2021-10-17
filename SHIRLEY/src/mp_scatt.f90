!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!------------------------------------------------------------------------------!
    MODULE mp_scatt
!------------------------------------------------------------------------------!
      USE kinds,     ONLY : DP
      USE parallel_include
      use mp, only : mp_stop
      !
      IMPLICIT NONE
!
      public :: mp_scatter, mp_scatter_size, mp_scatter_displ
!
      interface mp_scatter
        module procedure mp_scatter_i, mp_scatter_r, mp_scatter_c
        module procedure mp_scatter_i1, mp_scatter_r1, mp_scatter_c1
        module procedure mp_scatter_i2, mp_scatter_r2, mp_scatter_c2
        module procedure mp_scatter_iv, mp_scatter_rv, mp_scatter_cv
      end interface
!
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
! routines added by davegp to allow distributing of data across resources - scattering
!
!..mp_scatter_size
      subroutine mp_scatter_size( size_sour, size_dest, root, gid )
      ! for a given size_sour on root, provides a uniform distribution size_dest across procs in gid
      implicit none
      integer :: size_sour
      integer,intent(out) :: size_dest
      integer,intent(in) :: root
      integer,optional,intent(in) :: gid
      integer :: msglen, numtask, taskid, group, ierr, i
      integer,allocatable :: lengths(:)
#if defined (__MPI)
      ! message length
      msglen = 1
      ! group
      group = mpi_comm_world
      if( present(gid) ) group = gid
      call mpi_comm_size(group,numtask,ierr)
      if( ierr/=0 ) call mp_stop(8400)
      call mpi_comm_rank(group,taskid,ierr)
      if( ierr/=0 ) call mp_stop(8400)

      if( taskid==root ) then
        allocate( lengths(numtask) )      
        lengths = size_sour / numtask
        do i=1,size_sour - lengths(numtask)*numtask
          lengths(i)=lengths(i)+1
        enddo
      else
        allocate( lengths(1) )
      endif

      call mp_scatter( lengths, size_dest, root, group )

      deallocate( lengths )
#else
      msglen = 1
      size_dest = size_sour
#endif
      return
      end subroutine mp_scatter_size


!------------------------------------------------------------------------------!
!..mp_scatter_displ
      subroutine mp_scatter_displ( size_sour, displ_dest, root, gid )
      ! generates the displacement table for a uniform distribution of size_sour across procs in gid
      implicit none
      integer :: size_sour
      integer,intent(out) :: displ_dest
      integer,intent(in) :: root
      integer,optional,intent(in) :: gid
      integer :: msglen, numtask, taskid, group, ierr, i
      integer,allocatable :: lengths(:), displs(:)
#if defined (__MPI)
      ! message length
      msglen = 1
      ! group
      group = mpi_comm_world
      if( present(gid) ) group = gid
      call mpi_comm_size(group,numtask,ierr)
      if( ierr/=0 ) call mp_stop(8400)
      call mpi_comm_rank(group,taskid,ierr)
      if( ierr/=0 ) call mp_stop(8400)

      if( taskid==root ) then
        allocate( lengths(numtask), displs(numtask) )      
        lengths = size_sour / numtask
        do i=1,size_sour - lengths(numtask)*numtask
          lengths(i)=lengths(i)+1
        enddo
        displs(1)=0
        do i=2,numtask
          displs(i)=displs(i-1)+lengths(i-1)
        enddo
      else
        allocate( lengths(1), displs(1) )
      endif

      call mp_scatter( displs, displ_dest, root, group )

      deallocate( lengths, displs )
#else
      msglen = 1
      displ_dest = 0
#endif
      return
      end subroutine mp_scatter_displ


!------------------------------------------------------------------------------!
!..mp_scatter_i
      subroutine mp_scatter_i( msg_sour, msg_dest, root, gid )
      implicit none
      integer :: msg_sour(:)
      integer :: msg_dest
      integer,intent(in) :: root
      integer,optional,intent(in) :: gid

      integer :: msglen, group, numtask, ierr
#if defined (__MPI)
      ! message length
      msglen = 1
      ! group
      group = mpi_comm_world
      if( present(gid) ) group = gid
      ! check that size of group is consistent
      call mpi_comm_size(group,numtask,ierr)
      if( ierr/=0 ) call mp_stop(8400)

      call mpi_scatter( msg_sour, 1, MPI_INTEGER, &
                        msg_dest, 1, MPI_INTEGER, &
                        root, group, ierr )
      if( ierr/=0 ) call mp_stop(8403)
# else
      msglen = 1
      msg_dest = msg_sour(1)
# endif
      return
      end subroutine mp_scatter_i


!------------------------------------------------------------------------------!
!..mp_scatter_r
      subroutine mp_scatter_r( msg_sour, msg_dest, root, gid )
      implicit none
      real(dp) :: msg_sour(:)
      real(dp) :: msg_dest
      integer,intent(in) :: root
      integer,optional,intent(in) :: gid

      integer :: msglen, group, numtask, ierr
#if defined (__MPI)
      ! message length
      msglen = 1
      ! group
      group = mpi_comm_world
      if( present(gid) ) group = gid
      ! check that size of group is consistent
      call mpi_comm_size(group,numtask,ierr)
      if( ierr/=0 ) call mp_stop(8400)

      call mpi_scatter( msg_sour, 1, MPI_DOUBLE_PRECISION, &
                        msg_dest, 1, MPI_DOUBLE_PRECISION, &
                        root, group, ierr )
      if( ierr/=0 ) call mp_stop(8403)
# else
      msglen = 1
      msg_dest = msg_sour(1)
# endif
      return
      end subroutine mp_scatter_r

!------------------------------------------------------------------------------!
!..mp_scatter_c
      subroutine mp_scatter_c( msg_sour, msg_dest, root, gid )
      implicit none
      complex(dp) :: msg_sour(:)
      complex(dp) :: msg_dest
      integer,intent(in) :: root
      integer,optional,intent(in) :: gid

      integer :: msglen, group, numtask, ierr
#if defined (__MPI)
      ! message length
      msglen = 1
      ! group
      group = mpi_comm_world
      if( present(gid) ) group = gid
      ! check that size of group is consistent
      call mpi_comm_size(group,numtask,ierr)
      if( ierr/=0 ) call mp_stop(8400)

      call mpi_scatter( msg_sour, 1, MPI_DOUBLE_COMPLEX, &
                        msg_dest, 1, MPI_DOUBLE_COMPLEX, &
                        root, group, ierr )
      if( ierr/=0 ) call mp_stop(8403)
# else
      msglen = 1
      msg_dest = msg_sour(1)
# endif
      return
      end subroutine mp_scatter_c


!------------------------------------------------------------------------------!
!..mp_scatter_i1
      subroutine mp_scatter_i1( msg_sour, msg_dest, &
                                root, gid )
      implicit none
      integer :: msg_sour(:)
      integer :: msg_dest(:)
      integer,intent(in) :: root
      integer,optional,intent(in) :: gid

      integer :: msglen, group, numtask, taskid, ierr, i
      integer,allocatable :: lengths(:), displs(:)
#if defined (__MPI)
      ! message length
      msglen = SIZE(msg_dest)
      ! group
      group = mpi_comm_world
      if( present(gid) ) group = gid

      call mpi_comm_size(group,numtask,ierr)
      if( ierr/=0 ) call mp_stop(8400)
      call mpi_comm_rank(group,taskid,ierr)
      if( ierr/=0 ) call mp_stop(8400)

      if( taskid == root ) then
        allocate( lengths(numtask), displs(numtask) )
      else
        allocate( lengths(1), displs(1) )
      endif
      call mpi_gather( size(msg_dest), 1, MPI_INTEGER, &
                       lengths, 1, MPI_INTEGER, &
                       root, group, ierr )
      if( taskid == root ) then
        displs(1)=0
        do i=2,numtask
          displs(i)=displs(i-1)+lengths(i-1)
        enddo
      endif

      call mpi_scatterv( msg_sour, lengths, displs, MPI_INTEGER, &
                         msg_dest, size(msg_dest), MPI_INTEGER, &
                         root, group, ierr )
      if( ierr/=0 ) call mp_stop(8403)

      deallocate( lengths, displs )
# else
      msglen = SIZE(msg_dest)
      msg_dest = msg_sour
# endif
      return
      end subroutine mp_scatter_i1


!------------------------------------------------------------------------------!
!..mp_scatter_r1
      subroutine mp_scatter_r1( msg_sour, msg_dest, &
                                root, gid )
      implicit none
      real(dp) :: msg_sour(:)
      real(dp) :: msg_dest(:)
      integer,intent(in) :: root
      integer,optional,intent(in) :: gid

      integer :: msglen, group, numtask, taskid, ierr, i
      integer,allocatable :: lengths(:), displs(:)
#if defined (__MPI)
      ! message length
      msglen = SIZE(msg_dest)
      ! group
      group = mpi_comm_world
      if( present(gid) ) group = gid

      call mpi_comm_size(group,numtask,ierr)
      if( ierr/=0 ) call mp_stop(8400)
      call mpi_comm_rank(group,taskid,ierr)
      if( ierr/=0 ) call mp_stop(8400)

      if( taskid == root ) then
        allocate( lengths(numtask), displs(numtask) )
      else
        allocate( lengths(1), displs(1) )
      endif
      call mpi_gather( size(msg_dest), 1, MPI_INTEGER, &
                       lengths, 1, MPI_INTEGER, &
                       root, group, ierr )
      if( taskid == root ) then
        displs(1)=0
        do i=2,numtask
          displs(i)=displs(i-1)+lengths(i-1)
        enddo
      endif

      call mpi_scatterv( msg_sour, lengths, displs, MPI_DOUBLE_PRECISION, &
                         msg_dest, size(msg_dest), MPI_DOUBLE_PRECISION, &
                         root, group, ierr )
      if( ierr/=0 ) call mp_stop(8403)

      deallocate( lengths, displs )
# else
      msglen = SIZE(msg_dest)
      msg_dest = msg_sour
# endif
      return
      end subroutine mp_scatter_r1

!------------------------------------------------------------------------------!
!..mp_scatter_c1
      subroutine mp_scatter_c1( msg_sour, msg_dest, &
                                root, gid )
      implicit none
      complex(dp) :: msg_sour(:)
      complex(dp) :: msg_dest(:)
      integer,intent(in) :: root
      integer,optional,intent(in) :: gid

      integer :: msglen, group, numtask, taskid, ierr, i
      integer,allocatable :: lengths(:), displs(:)
#if defined (__MPI)
      ! message length
      msglen = SIZE(msg_dest)
      ! group
      group = mpi_comm_world
      if( present(gid) ) group = gid

      call mpi_comm_size(group,numtask,ierr)
      if( ierr/=0 ) call mp_stop(8400)
      call mpi_comm_rank(group,taskid,ierr)
      if( ierr/=0 ) call mp_stop(8400)

      if( taskid == root ) then
        allocate( lengths(numtask), displs(numtask) )
      else
        allocate( lengths(1), displs(1) )
      endif
      call mpi_gather( size(msg_dest), 1, MPI_INTEGER, &
                       lengths, 1, MPI_INTEGER, &
                       root, group, ierr )
      if( taskid == root ) then
        displs(1)=0
        do i=2,numtask
          displs(i)=displs(i-1)+lengths(i-1)
        enddo
      endif

      call mpi_scatterv( msg_sour, lengths, displs, MPI_DOUBLE_COMPLEX, &
                         msg_dest, size(msg_dest), MPI_DOUBLE_COMPLEX, &
                         root, group, ierr )
      if( ierr/=0 ) call mp_stop(8403)

      deallocate( lengths, displs )
# else
      msglen = SIZE(msg_dest)
      msg_dest = msg_sour
# endif
      return
      end subroutine mp_scatter_c1


!------------------------------------------------------------------------------!
!..mp_scatter_i2
      subroutine mp_scatter_i2( msg_sour, msg_dest, &
                                root, gid )
      implicit none
      integer :: msg_sour(:,:)
      integer :: msg_dest(:,:)
      integer,intent(in) :: root
      integer,optional,intent(in) :: gid

      integer :: msglen, group, numtask, taskid, ierr, i
      integer,allocatable :: lengths(:), displs(:)
#if defined (__MPI)
      ! message length
      msglen = SIZE(msg_dest)
      ! group
      group = mpi_comm_world
      if( present(gid) ) group = gid

      call mpi_comm_size(group,numtask,ierr)
      if( ierr/=0 ) call mp_stop(8400)
      call mpi_comm_rank(group,taskid,ierr)
      if( ierr/=0 ) call mp_stop(8400)

      if( taskid == root ) then
        allocate( lengths(numtask), displs(numtask) )
      else
        allocate( lengths(1), displs(1) )
      endif
      call mpi_gather( size(msg_dest), 1, MPI_INTEGER, &
                       lengths, 1, MPI_INTEGER, &
                       root, group, ierr )
      if( taskid == root ) then
        displs(1)=0
        do i=2,numtask
          displs(i)=displs(i-1)+lengths(i-1)
        enddo
      endif

      call mpi_scatterv( msg_sour, lengths, displs, MPI_INTEGER, &
                         msg_dest, size(msg_dest), MPI_INTEGER, &
                         root, group, ierr )
      if( ierr/=0 ) call mp_stop(8403)

      deallocate( lengths, displs )
# else
      msglen = SIZE(msg_dest)
      msg_dest = msg_sour
# endif
      return
      end subroutine mp_scatter_i2


!------------------------------------------------------------------------------!
!..mp_scatter_r2
      subroutine mp_scatter_r2( msg_sour, msg_dest, &
                                root, gid )
      implicit none
      real(dp) :: msg_sour(:,:)
      real(dp) :: msg_dest(:,:)
      integer,intent(in) :: root
      integer,optional,intent(in) :: gid

      integer :: msglen, group, numtask, taskid, ierr, i
      integer,allocatable :: lengths(:), displs(:)
#if defined (__MPI)
      ! message length
      msglen = SIZE(msg_dest)
      ! group
      group = mpi_comm_world
      if( present(gid) ) group = gid

      call mpi_comm_size(group,numtask,ierr)
      if( ierr/=0 ) call mp_stop(8400)
      call mpi_comm_rank(group,taskid,ierr)
      if( ierr/=0 ) call mp_stop(8400)

      if( taskid == root ) then
        allocate( lengths(numtask), displs(numtask) )
      else
        allocate( lengths(1), displs(1) )
      endif
      call mpi_gather( size(msg_dest), 1, MPI_INTEGER, &
                       lengths, 1, MPI_INTEGER, &
                       root, group, ierr )
      if( taskid == root ) then
        displs(1)=0
        do i=2,numtask
          displs(i)=displs(i-1)+lengths(i-1)
        enddo
      endif

      call mpi_scatterv( msg_sour, lengths, displs, MPI_DOUBLE_PRECISION, &
                         msg_dest, size(msg_dest), MPI_DOUBLE_PRECISION, &
                         root, group, ierr )
      if( ierr/=0 ) call mp_stop(8403)

      deallocate( lengths, displs )
# else
      msglen = SIZE(msg_dest)
      msg_dest = msg_sour
# endif
      return
      end subroutine mp_scatter_r2

!------------------------------------------------------------------------------!
!..mp_scatter_c2
      subroutine mp_scatter_c2( msg_sour, msg_dest, &
                                root, gid )
      implicit none
      complex(dp) :: msg_sour(:,:)
      complex(dp) :: msg_dest(:,:)
      integer,intent(in) :: root
      integer,optional,intent(in) :: gid

      integer :: msglen, group, numtask, taskid, ierr, i
      integer,allocatable :: lengths(:), displs(:)
#if defined (__MPI)
      ! message length
      msglen = SIZE(msg_dest)
      ! group
      group = mpi_comm_world
      if( present(gid) ) group = gid

      call mpi_comm_size(group,numtask,ierr)
      if( ierr/=0 ) call mp_stop(8400)
      call mpi_comm_rank(group,taskid,ierr)
      if( ierr/=0 ) call mp_stop(8400)

      if( taskid == root ) then
        allocate( lengths(numtask), displs(numtask) )
      else
        allocate( lengths(1), displs(1) )
      endif
 
      call mpi_gather( size(msg_dest), 1, MPI_INTEGER, &
                       lengths, 1, MPI_INTEGER, &
                       root, group, ierr )
      if( taskid == root ) then
        displs(1)=0
        do i=2,numtask
          displs(i)=displs(i-1)+lengths(i-1)
        enddo
        print *, "taskid = ", taskid
        print *, "mp_scatter_c2 displs = ", displs, " lengths = ", lengths
      endif

      call mpi_scatterv( msg_sour, lengths, displs, MPI_DOUBLE_COMPLEX, &
                         msg_dest, size(msg_dest), MPI_DOUBLE_COMPLEX, &
                         root, group, ierr )
      if( ierr/=0 ) call mp_stop(8403)

      deallocate( lengths, displs )
# else
      msglen = SIZE(msg_dest)
      msg_dest = msg_sour
# endif
      return
      end subroutine mp_scatter_c2

!------------------------------------------------------------------------------!
!..mp_scatter_iv
      subroutine mp_scatter_iv( msg_sour, lengths, displs, msg_dest, &
                                root, gid )
      implicit none
      integer :: msg_sour(:)
      integer :: lengths(:), displs(:)
      integer :: msg_dest(:)
      integer,intent(in) :: root
      integer,optional,intent(in) :: gid

      integer :: msglen, group, ierr
#if defined (__MPI)
      ! message length
      msglen = SIZE(msg_dest)
      ! group
      group = mpi_comm_world
      if( present(gid) ) group = gid
      ! check that size of group is consistent
      call mpi_scatterv( msg_sour, lengths, displs, MPI_INTEGER, &
                         msg_dest, size(msg_dest), MPI_INTEGER, &
                         root, group, ierr )
      if( ierr/=0 ) call mp_stop(8403)
# else
      msglen = SIZE(msg_dest)
      msg_dest = msg_sour(displs(1)+1:displs(1)+lengths(1))
# endif
      return
      end subroutine mp_scatter_iv


!------------------------------------------------------------------------------!
!..mp_scatter_rv
      subroutine mp_scatter_rv( msg_sour, lengths, displs, msg_dest, &
                                root, gid )
      implicit none
      real(dp) :: msg_sour(:)
      integer :: lengths(:), displs(:)
      real(dp) :: msg_dest(:)
      integer,intent(in) :: root
      integer,optional,intent(in) :: gid

      integer :: msglen, group, ierr
#if defined (__MPI)
      ! message length
      msglen = SIZE(msg_dest)
      ! group
      group = mpi_comm_world
      if( present(gid) ) group = gid
      ! check that size of group is consistent
      call mpi_scatterv( msg_sour, lengths, displs, MPI_DOUBLE_PRECISION, &
                         msg_dest, size(msg_dest), MPI_DOUBLE_PRECISION, &
                         root, group, ierr )
      if( ierr/=0 ) call mp_stop(8403)
# else
      msglen = SIZE(msg_dest)
      msg_dest = msg_sour(displs(1)+1:displs(1)+lengths(1))
# endif
      return
      end subroutine mp_scatter_rv

!------------------------------------------------------------------------------!
!..mp_scatter_cv
      subroutine mp_scatter_cv( msg_sour, lengths, displs, msg_dest, &
                                root, gid )
      implicit none
      complex(dp) :: msg_sour(:)
      integer :: lengths(:), displs(:)
      complex(dp) :: msg_dest(:)
      integer,intent(in) :: root
      integer,optional,intent(in) :: gid

      integer :: msglen, group, ierr
#if defined (__MPI)
      ! message length
      msglen = SIZE(msg_dest)
      ! group
      group = mpi_comm_world
      if( present(gid) ) group = gid
      ! check that size of group is consistent
      call mpi_scatterv( msg_sour, lengths, displs, MPI_DOUBLE_COMPLEX, &
                         msg_dest, size(msg_dest), MPI_DOUBLE_COMPLEX, &
                         root, group, ierr )
      if( ierr/=0 ) call mp_stop(8403)
# else
      msglen = SIZE(msg_dest)
      msg_dest = msg_sour(displs(1)+1:displs(1)+lengths(1))
# endif
      return
      end subroutine mp_scatter_cv

! end of scatter routines
! davegp

!------------------------------------------------------------------------------!
    END MODULE mp_scatt
!------------------------------------------------------------------------------!

