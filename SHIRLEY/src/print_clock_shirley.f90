!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE print_clock_shirley()
   !---------------------------------------------------------------------------
   !
   ! ... this routine prints out the clocks at the end of the run
   ! ... it tries to construct the calling tree of the program.
   !
   USE io_global,     ONLY : stdout
   USE control_flags, ONLY : isolve
   USE force_mod,     ONLY : lforce, lstres
!   USE mp_global,     ONLY : mpime, root
   !
   IMPLICIT NONE
   !
   !
!   IF ( mpime /= root ) &
!      OPEN( UNIT = stdout, FILE = '/dev/null', STATUS = 'UNKNOWN' )
   !
   WRITE( stdout, * )
   !
   CALL print_clock( 'shirley' )
   !
   WRITE( stdout, '(5X,"General routines")' )
   !
   CALL print_clock( 'firstfft' )
   CALL print_clock( 'secondfft' )
!   CALL print_clock( 'ccalbec' )
   CALL print_clock( 'cft3' )
   CALL print_clock( 'cft3s' )
!   CALL print_clock( 'interpolate' )
   CALL print_clock( 'davcio' )
   !    
   WRITE( stdout, * )
   !
#if defined (__MPI)
   WRITE( stdout, '(5X,"Parallel routines")' )
   !
   CALL print_clock( 'reduce' )
   CALL print_clock( 'fft_scatter' )
   CALL print_clock( 'poolreduce' )
#endif
   !
#ifdef EXX
   WRITE( stdout, '(5X,"EXX routines")' )
   !
   CALL print_clock( 'exx_grid' )
   CALL print_clock( 'exxinit' )
   CALL print_clock( 'vexx' )
   CALL print_clock( 'exxenergy' )
   CALL print_clock( 'exxen2' )
   CALL print_clock ('cycleig')
#endif
   !
   RETURN
   !
END SUBROUTINE print_clock_shirley
