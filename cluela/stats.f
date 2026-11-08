*-----------------------------------------------------------------------
* stats.f - averages and deviations with missing data support
* by Andy Allinger, 2017, released to the public domain
* This program may be used by any person for any purpose.
*-----------------------------------------------------------------------

************************************************************************
* AVEDEV - summary statistics:  mean and standard deviation
* Returns 0.0 when not enough data is available
************************************************************************
*
*___Name_________Type______In/Out____Description________________________
*   X(PLIM,N)    Real      In        Data matrix
*   NX(PLIM,N)   Integer   In        1 if data exists for this observation
*   PLIM         Integer   In        Row bound of X and NX
*   P            Integer   In        Number of variables
*   N            Integer   In        Number of objects
*   NA(P)        Integer   Neither   Count of data available
*   AVE(P)       Real      Out       Means
*   DEV(P)       Real      Out       Standard deviations
************************************************************************
      SUBROUTINE AVEDEV (X, NX, PLIM, P, N, NA, AVE, DEV)
       IMPLICIT NONE
       INTEGER NX, PLIM, P, N, NA
       REAL X, AVE, DEV
       DIMENSION X(PLIM,N), NX(PLIM,N), NA(P), AVE(P), DEV(P)

       INTEGER I, J

*            begin
       DO J = 1, P                     ! clear
         AVE(J) = 0.
         NA(J) = 0
       END DO
       DO I = 1, N                     ! add
         DO J = 1, P
           IF (NX(J,I) > 0) THEN
             AVE(J) = AVE(J) + X(J,I)
             NA(J) = NA(J) + 1
           END IF
         END DO
       END DO
       DO J = 1, P                     ! divide
         IF (NA(J) > 0) AVE(J) = AVE(J) / FLOAT(NA(J))
       END DO

       DO J = 1, P                     ! clear
         DEV(J) = 0.
       END DO
       DO I = 1, N                     ! add
         DO J = 1, P
           IF (NX(J,I) > 0) DEV(J) = DEV(J) + (X(J,I) - AVE(J)) **2
         END DO
       END DO
       DO J = 1, P                     ! divide
         IF (NA(J) > 1) THEN
           DEV(J) = SQRT(DEV(J) / FLOAT(NA(J) - 1))
         ELSE
           DEV(J) = 0.
         END IF
       END DO
       RETURN
      END  ! of avedev


************************************************************************
* MEDMAD - summary statistics:  median and median absolute deviation (MAD)
* Returns 0.0 when not enough data is available
************************************************************************
*
*___Name_________Type______In/Out____Description________________________
*   X(PLIM,N)    Real      In        Data matrix
*   NX(PLIM,N)   Integer   In        1 if data exists for this observation
*   PLIM         Integer   In        Row bound of X and NX
*   P            Integer   In        Number of variables
*   N            Integer   In        Number of objects
*   RWORK(N)     Real      Neither   Workspace
*   AVE(P)       Real      Out       Medians
*   DEV(P)       Real      Out       Median absolute deviations
************************************************************************
      SUBROUTINE MEDMAD (X, NX, PLIM, P, N, RWORK, AVE, DEV)
       IMPLICIT NONE
       INTEGER NX, PLIM, P, N
       REAL X, RWORK, AVE, DEV
       DIMENSION X(PLIM,N), NX(PLIM,N), RWORK(N), AVE(P), DEV(P)

       INTEGER I, J, K

*            begin
       DO J = 1, P                         ! copy
         K = 0
         DO I = 1, N
           IF (NX(J,I) > 0) THEN
             K = K + 1
             RWORK(K) = X(J,I)
           END IF
           IF (K > 0) THEN
             CALL MEDIAN (RWORK, K, AVE(J))
           ELSE
             AVE(J) = 0.
           END IF
         END DO
       END DO

       DO J = 1, P                        ! clear deviations
         K = 0
         DO I = 1, N
           IF (NX(J,I) > 0) THEN
             K = K + 1
             RWORK(K) = ABS(X(J,I) - AVE(J))
           END IF
         END DO
         IF (K > 0) THEN
           CALL MEDIAN (RWORK, K, DEV(J))
         ELSE
           DEV(J) = 0.
         END IF
       END DO
       RETURN
      END  ! of medmad


************************************************************************
* LEAGRE - summary statistics:  least and greatest
* Returns 0.0 when not enough data is available
************************************************************************
*
*___Name_________Type______In/Out____Description________________________
*   X(PLIM,N)    Real      In        Data matrix
*   NX(PLIM,N)   Integer   In        1 if data exists for this observation
*   PLIM         Integer   In        Row bound of X and NX
*   P            Integer   In        Number of variables
*   N            Integer   In        Number of objects
*   NA(P)        Integer   Neither   Count of data available
*   LEA(P)       Real      Out       Least values
*   GRE(P)       Real      Out       Greatest values
************************************************************************
      SUBROUTINE LEAGRE (X, NX, PLIM, P, N, NA, LEA, GRE)
       IMPLICIT NONE
       INTEGER PLIM, P, N
       INTEGER NX(PLIM,N), NA(P)
       REAL X(PLIM,N), LEA(P), GRE(P)

       REAL BIG
       PARAMETER (BIG = 1.E36)
       INTEGER I, J

*             begin
       DO J = 1, P
         LEA(J) = +BIG
         GRE(J) = -BIG
         NA(J) = 0
       END DO
       DO I = 1, N
         DO J = 1, P
           IF (NX(J,I) > 0) THEN
             IF (X(J,I) < LEA(J)) LEA(J) = X(J,I)
             IF (X(J,I) > GRE(J)) GRE(J) = X(J,I)
             NA(J) = NA(J) + 1
           END IF
         END DO
       END DO
       DO J = 1, P
         IF (NA(J) .EQ. 0) THEN
           LEA(J) = 0.
           GRE(J) = 0.
         END IF
       END DO
       RETURN
      END  ! of leagre
