************************************************************************
* t_round.f - a test program for CLUELA
* By Andy Allinger, 2017-2020, released to the public domain
* This program may be used by any person for any purpose.
************************************************************************
      PROGRAM ROUND
       IMPLICIT NONE
       INTEGER IOU, N, PLIM
       REAL TWOPI
       CHARACTER*32 OFILE
       PARAMETER (OFILE='round_output.csv',
     $            IOU = 1,
     $            N = 30 000,
     $            PLIM = 2,
     $            TWOPI = 6.2831853)
       INTEGER KMIN, KMAX, KLIM, NP
       PARAMETER (KMIN = 1,
     $            KMAX = 10,
     $            KLIM = 1000,
     $            NP = N * PLIM)

       INTEGER I, J, L, ITER, Z(N), IFAULT, P, K, POP(KLIM),
     $  IWORK(4*N+25*KLIM), NX(PLIM,N), NC(PLIM,KLIM)
       REAL X(PLIM,N), C(PLIM,KLIM), R, THETA, RNUM, W(PLIM),
     $  F(KMIN:KMAX), U, URAND, RWORK(2*N+9*KLIM), THRESH,
     $  XT(N,PLIM), CT(KLIM,PLIM), TARRAY(2), RESULT
       LOGICAL ROBUST, SOSDO
       DATA W / PLIM * 1. /
       DATA NX / NP * 1 /

*             Begin
       CALL RNDINI ()
       DO I = 1, N                      ! random points over unit circle
         RNUM = URAND()
         R = SQRT(RNUM)                      ! radius
         RNUM = URAND()
         THETA = RNUM * TWOPI                ! angle
         X(1,I) = R * COS(THETA)          ! Cartesian coordinates
         X(2,I) = R * SIN(THETA)
       END DO

       K = 6
       P = PLIM
       ROBUST = .FALSE.                  ! edit this line for robust mode
       SOSDO = .NOT. ROBUST

*             Call subroutine
       CALL CLUELA (X, NX, PLIM, P, N, W, C, NC, KMIN, K, KMAX,
     $  Z, POP, F, ROBUST, SOSDO, IWORK, RWORK, U, IFAULT)
       PRINT *, 'cluela returns with error code ', IFAULT,
     $  ' made ', K, ' clusters with residual ', U

*         Write CSV file with Y-coordinates in different columns by cluster
       OPEN (UNIT=1, FILE=OFILE, STATUS='NEW', IOSTAT=IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 't_round: trouble opening file'

  10   FORMAT (F8.4, $)
  20   FORMAT (',', $)
  30   FORMAT (F8.4)
       DO I = 1, N
         WRITE (UNIT=1, FMT=10, IOSTAT=IFAULT) X(1,I)
         IF (IFAULT .NE. 0) PRINT *, 't_round: trouble writing X-coord'
         DO L = 1, Z(I)                     ! one comma per cluster ID
           WRITE (UNIT=1, FMT=20, IOSTAT=IFAULT)
           IF (IFAULT .NE. 0) PRINT *, 't_round: trouble writing comma'
         END DO
         WRITE (UNIT=1, FMT=30, IOSTAT=IFAULT) X(2,I)
         IF (IFAULT .NE. 0) PRINT *, 't_round: trouble writing Y-coord'
       END DO

*           Write the centroids in the far column
       DO L = 1, K
         WRITE (UNIT=1, FMT=10, IOSTAT=IFAULT) C(1,L)
         IF (IFAULT .NE. 0) PRINT *, 't_round: trouble writing X-coord'
         DO I = 1, K+1
           WRITE (UNIT=1, FMT=20, IOSTAT=IFAULT)
           IF (IFAULT .NE. 0) PRINT *, 't_round: trouble writing comma'
         END DO
         WRITE (UNIT=1, FMT=30, IOSTAT=IFAULT) C(2,L)
         IF (IFAULT .NE. 0) PRINT *, 't_round: trouble writing Y-coord'
       END DO

*         Write the results for F
  40   FORMAT (I3, ', ', F8.4)
       DO L = KMIN, KMAX
         WRITE (UNIT=1, FMT=40, IOSTAT=IFAULT) L, F(L)
         IF (IFAULT .NE. 0) PRINT *, 't_round: trouble writing F'
       END DO
       CLOSE (UNIT=1)

*          Robust mode isn't comparable to KMNS!
       IF (ROBUST) STOP 'round: done'

*             Time for CLUSTE
       SOSDO = .FALSE.
       K = KLIM
       CALL DTIME (TARRAY, RESULT)
       CALL CLUSTE (X, NX, PLIM, P, N, W, C, NC, KMIN, K, KLIM,
     $  Z, POP, ROBUST, SOSDO, IWORK, RWORK, U, IFAULT)
       CALL DTIME (TARRAY, RESULT)
       PRINT *, 'cluste returns #', IFAULT, ' in ', RESULT, ' seconds.',
     $  ' Made ', K, ' clusters with residual ', U

*             Time for KMNS
       K = KLIM
       ITER = 5000
       CALL DINIT (X, NX, PLIM, P, N, W, C, NC, K, ROBUST,
     $  RWORK(1), RWORK(N+1), THRESH)
       DO J = 1, P
         DO I = 1, N
           XT(I,J) = X(J,I)
         END DO
         DO L = 1, K
           CT(L,J) = C(J,L)
         END DO
       END DO
       CALL DTIME (TARRAY, RESULT)
       CALL KMNS (XT, N, P, CT, K, Z, IWORK(1), POP, RWORK(N+1),
     $  RWORK(K+N+1), IWORK(N+1), RWORK(1), IWORK(N+K+1),
     $  IWORK(N+2*K+1), ITER, RWORK(N+2*K+1), IFAULT)
       U = 0.
       DO L = 1, K
         U = U + RWORK(N+2*K+L)
       END DO
       CALL DTIME (TARRAY, RESULT)
       PRINT *, 'kmns resturns #', IFAULT, ' in ', RESULT, ' seconds.',
     $  ' Made ', K, ' clusters with residual ', U

      END  ! of t_round.f **********************************************
