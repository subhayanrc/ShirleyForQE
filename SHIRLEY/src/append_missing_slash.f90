  subroutine append_missing_slash( path )

  character(*),intent(inout) :: path

  integer :: l
 
  ! ... verify if path ends with /, add one if needed
  !
  l = LEN_TRIM( path )
  !
  IF ( path(l:l) /= '/' ) THEN
    !
    IF ( l > 0 .AND. l < LEN( path ) ) THEN
      !
      path(l+1:l+1) = '/'
      !
    ELSE
      !
      CALL errore( 'path: ', path // ' truncated or empty', 1 )
      !
    END IF
    !
  END IF
  !
  end subroutine append_missing_slash
