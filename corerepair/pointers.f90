  module pointers

  ! various utilities for handling and managing pointers and their
  ! memory management

  ! David Prendergast, 2005DEC11

  use kinds

  implicit none

  private

  public :: pointer_destroy

  interface pointer_destroy
    module procedure pointer_i_1d_destroy, &
                     pointer_d_1d_destroy, &
                     pointer_c_1d_destroy
    module procedure pointer_i_2d_destroy, &
                     pointer_d_2d_destroy, &
                     pointer_c_2d_destroy
  end interface pointer_destroy
  
  contains


  subroutine pointer_i_1d_destroy( p )

  integer,dimension(:),pointer :: p

  if( associated(p) ) then
    deallocate(p)
  else
    nullify(p)
  endif

  end subroutine pointer_i_1d_destroy


  subroutine pointer_d_1d_destroy( p )

  real(dp),dimension(:),pointer :: p

  if( associated(p) ) then
    deallocate(p)
  else
    nullify(p)
  endif

  end subroutine pointer_d_1d_destroy


  subroutine pointer_c_1d_destroy( p )

  character,dimension(:),pointer :: p

  if( associated(p) ) then
    deallocate(p)
  else
    nullify(p)
  endif

  end subroutine pointer_c_1d_destroy


  subroutine pointer_i_2d_destroy( p )

  integer,dimension(:,:),pointer :: p

  if( associated(p) ) then
    deallocate(p)
  else
    nullify(p)
  endif

  end subroutine pointer_i_2d_destroy


  subroutine pointer_d_2d_destroy( p )

  real(dp),dimension(:,:),pointer :: p

  if( associated(p) ) then
    deallocate(p)
  else
    nullify(p)
  endif

  end subroutine pointer_d_2d_destroy


  subroutine pointer_c_2d_destroy( p )

  character,dimension(:,:),pointer :: p

  if( associated(p) ) then
    deallocate(p)
  else
    nullify(p)
  endif

  end subroutine pointer_c_2d_destroy

  end module pointers
