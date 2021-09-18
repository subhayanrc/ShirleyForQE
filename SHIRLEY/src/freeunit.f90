! ---------------------------------------------------------------------- 
  function freeunit()
! ---------------------------------------------------------------------- 

  integer :: freeunit
  logical :: opnd
  integer :: iun

  iun=10
  opnd=.true.
  free_unit_loop: do while(opnd)
    iun=iun+1
    inquire(unit=iun,opened=opnd)
  enddo free_unit_loop
  freeunit=iun

  end function freeunit

