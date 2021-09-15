  MODULE WIGNER3J_MODULE

  ! evaluate the Wigner 3j symbol using the Racah formula
  ! the most brute-force, compatible way, suitable for medium lm value

  ! Yufeng Liang, LBNL, Dec 2014

  USE kinds, ONLY : DP 
  USE constants, ONLY : PI, GSMALL, SQRT2

  IMPLICIT NONE

  ! Only WIGNER3J_3 is used. WIGNER3J and WIGNER3J_2 are NOT STABLE.
  PUBLIC :: WIGNER3J, WIGNER3J_2, WIGNER3J_3, Y3, RY3

  CONTAINS

    ! Gaunt coefficient, defined as the integral of the product of three Ylm
    ! Integrate over Yl1m1 Yl2m2 Yl3m3
    ! Ylm here are ALL COMPLEX, not directly applicable to conventional atomic
    ! orbits, which are representd by REAL Ylm
    ! See http://en.wikipedia.org/wiki/Table_of_spherical_harmonics

    FUNCTION Y3( l1, l2, l3, m1, m2, m3 )

      INTEGER, INTENT(IN) :: l1, l2, l3, m1, m2, m3
      REAL(DP) :: Y3      

      Y3 = SQRT((2 * l1 + 1) * (2 * l2 + 1) * (2 * l3 + 1) / 4.0 / PI) &
         * WIGNER3J_3( l1, l2, l3,  0,  0,  0 ) &
         * WIGNER3J_3( l1, l2, l3, m1, m2, m3 ) 

    END FUNCTION Y3

    ! Gaunt coefficient for REAL RYlm
    ! Integrate over REAL RYl1m1 RYl2m2 RYl3m3
    ! Idea: expand real RYlm as a linear combinations of complex Ylm

    FUNCTION RY3( l1, l2, l3, m1, m2, m3 )

      INTEGER, INTENT(IN) :: l1, l2, l3, m1, m2, m3
      INTEGER :: M(3), J, J1, J2, J3
      COMPLEX(DP) :: Coeff(3, 2) = 0, res
      COMPLEX(DP), PARAMETER :: i = CMPLX(0, 1)
      REAL(DP) :: RY3

      M = (/ m1, m2, m3 /)
      ! Coefficient: from complex to real
      ! See real form in http://en.wikipedia.org/wiki/Spherical_harmonics
      DO J = 1, 3
      IF     ( M(J) == 0 ) THEN
        Coeff(J, 1 : 2) = (/ 1, 0 /) 
      ELSEIF ( M(J) <  0 ) THEN
        Coeff(J, 1 : 2) = i / SQRT2 * (/ -(-1) ** M(J), 1 /)
      ELSE
        Coeff(J, 1 : 2) = 1 / SQRT2 * (/ 1, (-1) ** M(J) /)
      END IF
      END DO

      res = 0.0
      DO J1 = 1, 2
      DO J2 = 1, 2
      DO J3 = 1, 2
        IF ( ABS(  Coeff(1, J1) * Coeff(2, J2) * Coeff(3, J3) ) < GSMALL ) CYCLE
        res = res + Coeff(1, J1) * Coeff(2, J2) * Coeff(3, J3) &
            * Y3(l1, l2, l3, (2 * J1 - 3) * m1, (2 * J2 - 3) * m2, (2 * J3 - 3) * m3)
      ENDDO
      ENDDO
      ENDDO
      RY3 = REAL(res)

    END FUNCTION RY3

    ! Wigner 3J symbol: method 1
    ! reference from the Eq. (7) on the website
    ! http://mathworld.wolfram.com/Wigner3j-Symbol.html

    FUNCTION WIGNER3J( j1, j2, j3, m1, m2, m3 )
    
      INTEGER, INTENT(IN) :: j1, j2, j3, m1, m2, m3
      INTEGER :: mu, t
      REAL(DP) :: x, res
      REAL(DP) :: WIGNER3J

      mu = MIN(j1 + m1, j1 - m1, j2 + m2, j2 - m2, j3 + m3, j3 - m3, &
               j2 + j3 - j1, j3 + j1 - j2, j1 + j2 - j3)

      res = 0.0
      DO t = 0, 2 * (j1 + j2 + j3)
        x = FAC(t) * FAC(j3 - j2 + t + m1) * FAC(j3 - j1 + t - m2) &
          * FAC(j1 + j2 - j3 - t) * FAC(j1 - t - m1) * FAC(j2 - t + m2)
        IF (ABS(x) < gsmall) CYCLE
        res = res + (-1) ** t / x
      END DO

      WIGNER3J = (-1) ** (j1 - j2 - m3) * SQRT(DELTA(j1, j2, j3)) &
               * SQRT(FAC(j1 + m1) * FAC(j1 - m1)) &
               * SQRT(FAC(j2 + m2) * FAC(j2 - m2)) &
               * SQRT(FAC(j3 + m3) * FAC(j3 - m3)) &
               * res
   
    END FUNCTION WIGNER3J

    ! Wigner 3J symbol: method 2
    ! http://www.ams.org/journals/mcom/1996-65-216/S0025-5718-96-00774-0/S0025-5718-96-00774-0.pdf
    ! Something's wrong

    FUNCTION WIGNER3J_2( j1, j2, j3, m1, m2, m3 )

      INTEGER, INTENT(IN) :: j1, j2, j3, m1, m2, m3
      INTEGER :: k
      REAL(DP) :: x, y, res
      REAL(DP) :: WIGNER3J_2

      res = 0.0
      IF ( m1 + m2 + m3 /= 0 ) THEN
        WIGNER3J_2 = 0.0
        RETURN
      END IF

      DO k = 0, 2 * (j1 + j2 + j3)
        x = FAC(j2 + j3 - m1 - k) * FAC(j1 - m1 + k)
        y = FAC(k) * FAC(j3 - j1 + j2 - k) * FAC(j3 - m3 - k) * FAC(k + j1 - j2 + m3)
        IF (ABS(y) < gsmall) CYCLE
        res = res + (-1) ** (k + j2 + m2) * x / y
      END DO

      WIGNER3J_2 = (-1) ** (j1 - j2 - m3) * SQRT(DELTA(j1, j2, j3)) &
                 / SQRT(FAC(j1 + m1) * FAC(j1 - m1)) &
                 / SQRT(FAC(j2 + m2) * FAC(j2 - m2)) &
                 * SQRT(FAC(j3 + m3) * FAC(j3 - m3)) &
                 * res
      
    END FUNCTION WIGNER3J_2

    ! Wigner 3J symbol: method 3
    ! http://epubs.siam.org/doi/pdf/10.1137/S1064827503422932
    ! same as in Mathematica

    FUNCTION WIGNER3J_3( j1, j2, j3, m1, m2, m3)
    
      INTEGER, INTENT(IN) :: j1, j2, j3, m1, m2, m3
      INTEGER :: k, kmin, kmax
      REAL(DP) :: x, res
      REAL(DP) :: WIGNER3J_3
      
      IF ( m1 + m2 + m3 /= 0 ) THEN
        WIGNER3J_3 = 0.0
        RETURN
      ENDIF

      res = 0.0
      kmax = min(j1 + j2 - j3, j1 - m1, j2 + m2)
      kmin = max(-j3 + j2 - m1, -j3 + j1 + m2, 0)
      DO k = kmin, kmax
        x = FAC(k) * FAC(j1 + j2 - j3 - k) &
          * FAC(j1 - m1 - k) * FAC(j2 + m2 - k) &
          * FAC(j3 - j2 + m1 + k) * FAC(j3 - j1 - m2 + k)
        res = res + (-1) ** k / x
      END DO
      
      WIGNER3J_3 = (-1) ** (j1 - j2 - m3) * SQRT(DELTA(j1, j2, j3)) &
                 * SQRT(FAC(j1 + m1) * FAC(j1 - m1)) &
                 * SQRT(FAC(j2 + m2) * FAC(j2 - m2)) &
                 * SQRT(FAC(j3 + m3) * FAC(j3 - m3)) &
                 * res

    END FUNCTION WIGNER3J_3

    FUNCTION DELTA( A, B, C )
    
      INTEGER, INTENT(IN) :: A, B, C
      REAL(DP) :: DELTA 

      DELTA = FAC(A + B - C) * FAC(A - B + C) * FAC(-A + B + C) &
            / FAC(A + B + C + 1)

    END FUNCTION DELTA

    RECURSIVE FUNCTION FAC(N) RESULT(FAC_N)
      
      INTEGER, INTENT(IN) :: N
      REAL(DP) :: FAC_N
   
      SELECT CASE(N)
      CASE(0)
        FAC_N = 1.0
      CASE(1:)
        FAC_N = N * FAC(N - 1)
      CASE DEFAULT
        FAC_N = 0.0
      END SELECT

    END FUNCTION FAC

  END MODULE WIGNER3J_MODULE
