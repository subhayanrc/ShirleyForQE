!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
subroutine stop_shirley
  !--------------------------------------------------------------------
  !
  ! Synchronize processes before stopping.
  !
  use control_flags, only: twfcollect
  use io_files, only: iunwfc
  use mp_global, only: mp_global_end
  use environment, only : environment_end
  USE parallel_include

  character(len=9) :: code = 'SHIRLEY'
#ifdef __MPI
  integer :: info
  logical :: op

  inquire ( iunwfc, opened = op )

  if ( op ) then
     if (twfcollect) then
        close (unit = iunwfc, status = 'delete')
     else
        close (unit = iunwfc, status = 'keep')
     end if
  end if 

  call print_clock_shirley()

  !call mp_barrier()

  ! call mpi_finalize (info)
#endif
 
  CALL environment_end( code )
  !
  call mp_global_end()

#ifdef __T3E
  !
  ! set streambuffers off
  !

  call set_d_stream (0)
#endif

  stop
end subroutine stop_shirley
