    !----------------------------------------------------------------------------
    SUBROUTINE gk_l2gmap( ngm, ig_l2g, ngk, igk, igk_l2g )
      !----------------------------------------------------------------------------
      !
      ! ... This subroutine maps local G+k index to the global G vector index
      ! ... the mapping is used to collect wavefunctions subsets distributed
      ! ... across processors.
      ! ... Written by Carlo Cavazzoni
      !
      IMPLICIT NONE
      !
      ! ... Here the dummy variables
      !
      INTEGER, INTENT(IN)  :: ngm, ngk, igk(ngk), ig_l2g(ngm)
      INTEGER, INTENT(OUT) :: igk_l2g(ngk)
      INTEGER              :: ig
      !
      ! ... input: mapping between local and global G vector index
      !
      DO ig = 1, ngk
         !
         igk_l2g(ig) = ig_l2g(igk(ig))
         !
      END DO
      !
      RETURN
      !
    END SUBROUTINE gk_l2gmap

