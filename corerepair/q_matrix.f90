  SUBROUTINE q_matrix( n_f_, l_f, r_f_, f_, &
                       n_g_, l_g, r_g_, g_, &
                       rmax, qqq )

  ! Compute the Q overlap term between two sets of atomic waves
   
  ! David Prendergast, LBNL July 2018 
  
    USE kinds
    USE splines

    IMPLICIT NONE
  
    INTEGER,  PARAMETER :: ideg = 5
    REAL(DP), PARAMETER :: zero = 0.d0
    REAL(DP), PARAMETER :: eps12 = 1.d-12
    INTEGER, PARAMETER :: max_l = 4

    ! 3 dimensions of Jmat: ixyz, m_f, m_g; l_f and l_g are fixed
    ! ixyz ranges from 1..9 = x, y, z, x^2, y^2, z^2, xy, yz, zx
    ! m ranges from -l..l
    REAL(DP), ALLOCATABLE :: Jmat(:, :, :) 
    REAL(DP) :: RY3_22, RY3_20, RY3_00
    INTEGER :: m_f, m_g, ixyz

    ! in corevalence_position, f for valence and g for core
    INTEGER :: n_f, l_f
    INTEGER :: n_f_
    REAL(DP) :: r_f_(n_f_), f_(n_f_)
    REAL(DP),allocatable :: r_f(:), f(:)
    INTEGER :: n_g, l_g
    INTEGER :: n_g_
    REAL(DP) :: r_g_(n_g_), g_(n_g_)
    REAL(DP),allocatable :: r_g(:), g(:)
    REAL(DP) :: rmax
    REAL(DP) :: posmat(9, 2 * l_f + 1, 2 * l_g + 1)
  
    INTEGER :: i
    INTEGER :: n_fmax
    INTEGER :: ierr
    type(spline_struct) :: spl
    ! vector
    REAL(DP),allocatable :: gbr(:)
    ! tensor
    REAL(DP),allocatable :: gbrr(:)
    REAL(DP) :: fgbr_int, fgbrr_int
    INTEGER :: inonzero
    INTEGER :: imap
    
    ALLOCATE(Jmat(9, -l_f : l_f, -l_g : l_g))
    Jmat = 0

    DO m_f = -l_f, l_f
    DO m_g = -l_g, l_g

      ! x
      Jmat(1, m_f, m_g) = SQRT(4 * PI / 3) * RY3(l_f, l_g, 1, m_f, m_g,  1)

      ! y
      Jmat(2, m_f, m_g) = SQRT(4 * PI / 3) * RY3(l_f, l_g, 1, m_f, m_g, -1)

      ! z
      Jmat(3, m_f, m_g) = SQRT(4 * PI / 3) * RY3(l_f, l_g, 1, m_f, m_g,  0)

      RY3_22 = RY3(l_f, l_g, 2, m_f, m_g, 2) ! lm=(2,2): dx2-y2 = SQRT(15/pi)/4(x^2-y^2)
      RY3_20 = RY3(l_f, l_g, 2, m_f, m_g, 0) ! lm=(2,0): dz2 = SQRT(5/pi)/4(-x^2-y^2+2z^2)
      RY3_00 = RY3(l_f, l_g, 0, m_f, m_g, 0) ! lm=(0,0): s = SQRT(1/pi)/2

      ! x^2 = 1 / 2 (x^2-y^2) - 1 / 6 (-x^2-y^2+2z^2) + 1 / 3 (x^2+y^2+z^2)
      Jmat(4, m_f, m_g) = 1 / 2.0_DP * SQRT(PI / 15) * 4 * RY3_22 & 
                        - 1 / 6.0_DP * SQRT(PI / 5)  * 4 * RY3_20 &
                        + 1 / 3.0_DP * SQRT(PI / 1)  * 2 * RY3_00

      ! y^2 =-1 / 2 (x^2-y^2) - 1 / 6 (-x^2-y^2+2z^2) + 1 / 3 (x^2+y^2+z^2)
      Jmat(5, m_f, m_g) =-1 / 2.0_DP * SQRT(PI / 15) * 4 * RY3_22 & 
                        - 1 / 6.0_DP * SQRT(PI / 5)  * 4 * RY3_20 &
                        + 1 / 3.0_DP * SQRT(PI / 1)  * 2 * RY3_00

      ! z^2 = 1 / 3 (-x^2-y^2+2z^2) + 1 / 3 (x^2+y^2+z^2)
      Jmat(6, m_f, m_g) = 1 / 3.0_DP * SQRT(PI / 5) * 4 * RY3_20 &
                        + 1 / 3.0_DP * SQRT(PI / 1) * 2 * RY3_00

      ! xy
      Jmat(7, m_f, m_g) = SQRT(4 * PI / 15) * RY3(l_f, l_g, 2, m_f, m_g, -2)

      ! yz
      Jmat(8, m_f, m_g) = SQRT(4 * PI / 15) * RY3(l_f, l_g, 2, m_f, m_g, -1)

      ! zx
      Jmat(9, m_f, m_g) = SQRT(4 * PI / 15) * RY3(l_f, l_g, 2, m_f, m_g, 1)

    END DO ! im_g
    END DO ! im_f

    ! quick check to see IF integrals are zero by angular symmetry
    IF ( all( abs(Jmat(:, :, :)) < eps12 ) ) THEN
      posmat = zero
      RETURN
    ENDIF

    ! add zeros
    allocate( r_f(n_f_ + 1), f(n_f_ + 1) )
    allocate( r_g(n_g_ + 1), g(n_g_ + 1) )

    IF( r_f_(1) > zero ) THEN
      r_f = (/ zero, r_f_ /)
      f = (/ zero, f_ /)
      n_f = n_f_ + 1
    ELSE
      r_f = r_f_
      f = f_
      n_f = n_f_ 
    ENDIF
  
    IF( r_g_(1) > zero ) THEN
      r_g = (/ zero, r_g_ /)
      g = (/ zero, g_ /)
      n_g = n_g_ + 1
    ELSE
      r_g = r_g_
      g = g_
      n_g = n_g_ 
    ENDIF
  
    ! check that rmax inside limits of radial grids
    IF( rmax > r_f(n_f) .or. rmax > r_g(n_g) ) THEN
      WRITE(0,*) ' position_matrix : rmax too large for these radial grids'
      WRITE(0,*) ' rmax = ', rmax
      WRITE(0,*) ' r_f(n_f) = ', r_f(n_f)
      WRITE(0,*) ' r_g(n_g) = ', r_g(n_g)
      STOP 1
    ENDIF
  
    DO i = 1, n_f
      IF( r_f(i) > rmax ) exit
    ENDDO
    IF( i>=n_f ) THEN
      n_fmax = n_f
    ELSE
      n_fmax = i+1
    ENDIF
  
    ! should I modify rmax
    !IF( r_g(n_g) < r_f(n_fmax) ) THEN
    !  n_fmax = n_fmax - 1
    !  rmax = r_f(n_fmax)
    !ENDIF
    ! WRITE(*,*) 'rmax = ', rmax
  
    allocate( gbr(n_fmax) )
    allocate( gbrr(n_fmax) )
  
    ! spline g
    call spline_fit( ideg, n_g, r_g, g, spl, ierr, zero, rmax )
    IF( ierr > 0 ) THEN
      WRITE(0,*) ' spline error 1', ierr
      STOP
    ENDIF
  
    ! evaluate g on f grid
    ! WRITE(*,*) 'check :', n_fmax, r_f(1), r_f(n_fmax), zero, rmax
    call spline_eval( spl, n_fmax, r_f, gbr, ierr )
    IF( ierr > 0 ) THEN
      WRITE(0,*) ' spline error 2', ierr
      STOP
    ENDIF
  
    ! compute f g/r
    do i=1,n_fmax
      if( r_f(i) > zero ) then
        gbr(i)  = f(i) * gbr(i) * r_f(i)
        gbrr(i) = f(i) * gbr(i) * r_f(i) * r_f(i)
      else
        gbr(i)  = zero
        gbrr(i) = zero
      endif
    enddo 

    ! position vector
    call spline_fit( ideg, n_fmax, r_f, gbr, spl, ierr , zero, r_f(n_fmax) )
    IF( ierr > 0 ) THEN
      WRITE(0,*) ' position vector: spline error 5', ierr
      STOP
    ENDIF
  
    fgbr_int = spline_integral( spl, zero, rmax )
     ! WRITE(*,*) 'fgbr_int = ', fgbr_int
    
    ! position tensor
    call spline_fit( ideg, n_fmax, r_f, gbrr, spl, ierr , zero, r_f(n_fmax) )
    IF( ierr > 0 ) THEN
      WRITE(0,*) ' position tensor: spline error 5', ierr
      STOP
    ENDIF
  
    fgbrr_int = spline_integral( spl, zero, rmax )
    ! WRITE(*,*) 'fgbrr_int = ', fgbrr_int
  
    ! WRITE(*,*) 'angular contributions'
    DO m_f = -l_f, l_f
    DO m_g = -l_g, l_g
      ! vector
      DO ixyz = 1, 3 
        posmat(ixyz, imap(l_f, m_f), imap(l_g, m_g)) = fgbr_int  * Jmat( ixyz, m_f, m_g )
      ENDDO
      ! tensor
      DO ixyz = 4, 9
        posmat(ixyz, imap(l_f, m_f), imap(l_g, m_g)) = fgbrr_int * Jmat( ixyz, m_f, m_g )
      ENDDO
      ! WRITE(*,'(2x,4i4,9e14.7)') l_f, m_f, l_g, m_g, (Jmat(ixyz,m_f,l_f,m_g,l_g), ixyz=1,9)
    ENDDO
    ENDDO
  
    DEALLOCATE(Jmat)

  END SUBROUTINE position_matrix

  ! from conventional (l, m) in Quantum Mechanics to positive-definite m
  FUNCTION imap(l, m)

    INTEGER :: l, m, imap

    imap = ABS(m) * 2
    IF ( m <= 0 ) THEN
      imap = imap + 1
    ENDIF

  END FUNCTION

  ! Backward
  FUNCTION mimap(l, i)
    
    INTEGER :: l, i, mimap
  
    SELECT CASE(MOD(i,2))
    CASE(1)
      mimap = -(i - 1) / 2
    CASE(0)
      mimap = i / 2
    END SELECT

  END FUNCTION
 
