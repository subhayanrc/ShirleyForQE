  module fileio

  integer,parameter :: stdin  = 5
  integer,parameter :: stdout = 6

  contains

! ---------------------------------------------------------------------- 
  function freeunit(iunit_seed)
! ---------------------------------------------------------------------- 

  integer :: freeunit
  integer,optional,intent(in) :: iunit_seed
  logical :: opnd
  integer :: iunit

  ! Avoid the units used for input (5) and output (6)
  ! and any erroneous negative values
  if( present(iunit_seed) ) then
    iunit=iunit_seed
    if( iunit<7 ) iunit=7
  else
    iunit=7
  endif

  opnd=.true.
  do while(opnd)
    iunit=iunit+1
    inquire(unit=iunit,opened=opnd)
  enddo

  freeunit=iunit

  end function freeunit


  function strlen_int( i )

  integer :: strlen_int
  integer,intent(in) :: i

  integer :: l

  l = log10(real(abs(i)))+1
  if( i < 0 ) l=l+1

  strlen_int = l  

  end function strlen_int


  end module fileio
