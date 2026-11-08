*-----------------------------------------------------------------------
* dump.f - Write CSV scatter plot, with clusters by column
* by Andy Allinger, 2017, released to the public domain
* This program may be used by any person for any purpose.
*
*___Name_____Type______I/O___Description________________________________
*   R(N,M)   Real      In    Data points projected into two dimensions
*   N        Integer   In    Number of objects
*   M        Integer   In    Projection dimension, 2, always
*   Z(N)     Integer   In    Cluster ID's
*   IOU      Integer   In    Unit number of an open I/O unit
*   IFAULT   Integer   Out   Error code:
*                               0 = no errors
*                              -1 = Z out of bounds
*                              -2 = M is not 2
*                               other as returned by WRITE statements
*-----------------------------------------------------------------------
      SUBROUTINE WCSVSP (R, N, M, Z, IOU, IFAULT)
       IMPLICIT NONE
       INTEGER N, M, Z, IOU, IFAULT
       REAL R
       DIMENSION R(N,M), Z(N)

       INTEGER I, L

*          Validate
       DO I = 1, N
         IF (Z(I) < 1) THEN
           IFAULT = -1
           RETURN
         END IF
       END DO

       IF (M .NE. 2) THEN
         IFAULT = -2
         RETURN
       END IF

*           One record for each object
       DO I = 1, N
  70     FORMAT (G9.3, $)
         WRITE (UNIT=IOU, FMT=70, IOSTAT=IFAULT) R(I,1)  ! First coordinate
         IF (IFAULT .NE. 0) THEN
           PRINT *, 'wcsvsp: trouble writing X-coord'
           RETURN
         END IF

         DO L = 1, Z(I)               ! One comma per cluster ID number
  80       FORMAT (',', $)
           WRITE (UNIT=IOU, FMT=80, IOSTAT=IFAULT)
           IF (IFAULT .NE. 0) THEN
             PRINT *, 'wcsvsp: trouble writing comma'
             RETURN
           END IF
         END DO

  90     FORMAT (G9.3)
         WRITE (UNIT=IOU, FMT=90, IOSTAT=IFAULT) R(I,2)  ! Second coordinate
         IF (IFAULT .NE. 0) THEN
           PRINT *, 'wcsvsp: trouble writing Y-coord'
           RETURN
         END IF
       END DO

       RETURN
      END  ! of WCSVSP
*++++++++++++++++++++++++++ End of file dump.f +++++++++++++++++++++++++
