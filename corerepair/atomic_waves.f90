  module atomic_waves

  use kinds

  implicit none

  type atomic
    integer :: ngrid, l
    real(dp),pointer :: r(:) => null()
    real(dp),pointer :: f(:) => null()
  end type atomic

  contains

  subroutine init_atomic_wave( ngrid, l, wave )

  integer :: ngrid, l
  type(atomic) :: wave

  if( ngrid<1 ) return
  call destroy_atomic_wave( wave )
  wave%ngrid = ngrid
  wave%l = l
  allocate( wave%r(ngrid) )
  allocate( wave%f(ngrid) )

  return
  end subroutine init_atomic_wave

  subroutine destroy_atomic_wave( wave )

  type(atomic) :: wave

  if( associated(wave%r) ) deallocate(wave%r)
  if( associated(wave%f) ) deallocate(wave%f)

  return
  end subroutine destroy_atomic_wave

  end module atomic_waves
