  module shirley_constants

  use kinds

  integer,parameter :: maxchar = 255

  real(dp),parameter :: sqrt2 = 1.414213562d0
  real(dp),parameter :: pi    = 3.141592654d0

  complex(dp),parameter :: ONE  = cmplx( 1.d0, 0.d0 )
  complex(dp),parameter :: ZERO = cmplx( 0.d0, 0.d0 )
  complex(dp),parameter :: IOTA = cmplx( 0.d0, 1.d0 )


  real(dp),parameter :: rytoev=13.6058d0
  real(dp),parameter :: kelvin2rydberg = 6.3336303d-6
  real(dp),parameter :: speed_of_light = 13.7035999074d0

  end module shirley_constants
