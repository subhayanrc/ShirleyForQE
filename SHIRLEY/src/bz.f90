  module bz

  contains

! ---------------------------------------------------------------------- 
  subroutine map2bz( x )
! ---------------------------------------------------------------------- 
  ! maps vectors to the positive octant [0,1)^3
  use kinds, only : dp

  implicit none

  real(dp),intent(inout)  :: x(3)

  x = x - floor( x )

  return

  contains

  ! ---------------------------------------------------------------------- 
    elemental function floor_fn( f )
  ! ---------------------------------------------------------------------- 
    ! floor function
    real(dp),intent(in) :: f
    real(dp) :: floor_fn
    floor_fn = int(f - 0.5d0 + sign(0.5d0,f))

    end function floor_fn


  end subroutine map2bz


  end module bz
