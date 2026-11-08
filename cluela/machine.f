*-----------------------------------------------------------------------
* machine.f - front end to FORTRAN 95 machine constant routines
*-----------------------------------------------------------------------
      FUNCTION rmachine (in)
       IMPLICIT NONE
       INTEGER in
       REAL rmachine, r
       DATA r / 1.0 /          ! avoid compiler whimpering

       rmachine = 0.
       IF (in < 1 .OR. in > 3) RETURN
       IF (in .EQ. 1) THEN
         rmachine = EPSILON(r)  ! machine precision
       ELSE IF (in .EQ. 2) THEN
         rmachine = TINY(r)     ! smallest positive number
       ELSE IF (in .EQ. 3) THEN
         rmachine = HUGE(r)     ! largest positive number
       END IF
       RETURN
      END  ! of rmachine


*-----------------------------------------------------------------------
      FUNCTION dmachine (in)
       IMPLICIT NONE
       INTEGER in
       DOUBLE PRECISION dmachine, d
       DATA d / 1D0 /

       dmachine = 0D0
       IF (in < 1 .OR. in > 3) RETURN
       IF (in .EQ. 1) THEN
         dmachine = EPSILON(d)  ! machine precision
       ELSE IF (in .EQ. 2) THEN
         dmachine = TINY(d)     ! smallest positive number
       ELSE IF (in .EQ. 3) THEN
         dmachine = HUGE(d)     ! largest positive number
       END IF
       RETURN
      END  ! of dmachine


*-----------------------------------------------------------------------
      FUNCTION imachine (in)
       IMPLICIT NONE
       INTEGER in, imachine, i
       REAL r
       DOUBLE PRECISION d
       DATA i, r, d / 1, 1.0, 1D0 /

       imachine = 0
       IF (in < 1 .OR. in > 10) RETURN
       IF (in .EQ. 1) THEN
         imachine = RADIX(i)
       ELSE IF (in .EQ. 2) THEN
         imachine = BIT_SIZE(i)
       ELSE IF (in .EQ. 3) THEN
         imachine = HUGE(i)
       ELSE IF (in .EQ. 4) THEN
         imachine = RADIX(r)
       ELSE IF (in .EQ. 5) THEN
         imachine = DIGITS(r)
       ELSE IF (in .EQ. 6) THEN
         imachine = MINEXPONENT(r)
       ELSE IF (in .EQ. 7) THEN
         imachine = MAXEXPONENT(r)
       ELSE IF (in .EQ. 8) THEN
         imachine = DIGITS(d)
       ELSE IF (in .EQ. 9) THEN
         imachine = MINEXPONENT(d)
       ELSE IF (in .EQ. 10) THEN
         imachine = MAXEXPONENT(d)
       END IF
       RETURN
      END  ! of imachine
*++++++++++++++++++++++++++ End of file machine.f ++++++++++++++++++++++
