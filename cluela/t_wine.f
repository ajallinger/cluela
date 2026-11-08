*=======================================================================
* Cluster the Wine data
* by Andy Allinger, 2011-2017, public domain
*=======================================================================
      PROGRAM  WINE
       IMPLICIT NONE
       INTEGER G, IOU, KLIM, M, N, P
       PARAMETER (
     $            G = 3,          ! # correct classes
     $            IOU = 1,        ! input/output unit #
     $            KLIM = 20,      ! bound for max. clusters
     $            M = 2,          ! projection dimension
     $            N = 178,        ! # objects
     $            P = 13)         ! # variables
       INTEGER KMIN, KMAX
       LOGICAL SOSDO
       PARAMETER (KMIN=1,
     $            KMAX=KLIM,
     $            SOSDO = .FALSE.)
       INTEGER PN
       PARAMETER (PN = P * N)

       REAL X(P,N), W(P), WSPEC(P), C(P,KLIM), F(KMIN:KMAX), U,
     $       RWORK(2*N+9*KMAX), AVE(P), DEV(P)
       INTEGER K, POP(KLIM), Y(N), Z(N), IWORK(4*N+25*KMAX),
     $    NX(P,N), NC(P,KMAX), NA(P)
       INTEGER I, J, L, IFAULT
       REAL ACCARI, ACCVI, VALTAB(G*KLIM+G+KLIM)
       LOGICAL FEXIST, ROBUST
       CHARACTER*32 FNAME
       LOGICAL SAFDIV

*        workspace for principal component analysis
       INTEGER Q
       REAL E(P), V(P,P), R(N,M)

       DATA NX / PN * 1 /
       DATA W / P * 1. /
       DATA WSPEC / 2.76, 0.648, 1.42, 0.581, 0.695, 0.0516, 1.92,
     $              0.05, 0.05, 1.66, 1.90, 3.11, 4.32 /

*=======================================================================
*        begin
*=======================================================================
       CALL RNDINI ()
       IFAULT = 0

*        get input file
       FNAME = 'wine.data'
       INQUIRE (FILE=FNAME, EXIST=FEXIST, IOSTAT=IFAULT)
       IF (.NOT. FEXIST) STOP 'File not found'
       IF (IFAULT .NE. 0) STOP 'Error inquiring of input file'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='OLD',
     $        ACCESS='SEQUENTIAL', IOSTAT=IFAULT)
       IF (IFAULT .NE. 0) STOP 'trouble opening input file'

*       read in data
       DO I = 1, N
         READ (UNIT=1, FMT=*, END=30) Y(I), (X(J,I), J = 1, P)
       END DO
  30   CLOSE (UNIT=1)

*=======================================================================
*              Call the subroutine
*=======================================================================
       PRINT *, 'on the raw data: '
       ROBUST = .FALSE.
       K = 3
       CALL CLUELA (X, NX, P, P, N, W, C, NC, KMIN, K, KMAX, Z, POP,
     $  F, ROBUST, SOSDO, IWORK, RWORK, U, IFAULT)
       PRINT *, 'cluela returns #', IFAULT, ' Made ', K, ' clusters'

*        accuracy calculation
       CALL ARIND (Y, Z, N, G, K, VALTAB, ACCARI)
       CALL VINFO (Y, Z, N, G, K, VALTAB, ACCVI)
       PRINT *, 'residual: ', U, ' ARI: ', ACCARI, ' VI: ', ACCVI

*=======================================================================
* Standardize:  attempt to replicate Pham, Dimov & Nguyen
*=======================================================================
       CALL AVEDEV (X, NX, P, P, N, NA, AVE, DEV)
       DO J = 1, P
         IF (SAFDIV(1., DEV(J))) THEN
           DEV(J) = 1. / DEV(J)
         ELSE
           DEV(J) = 0.
         END IF
       END DO

       DO I = 1, N
         DO J = 1, P
           X(J,I) = (X(J,I) - AVE(J)) * DEV(J)
         END DO
       END DO

*=======================================================================
*              Call the subroutine
*=======================================================================
       K = 3
       PRINT *, 'attempting to replicate Pham...'
       CALL CLUELA (X, NX, P, P, N, W, C, NC, KMIN, K, KMAX, Z, POP,
     $  F, ROBUST, SOSDO, IWORK, RWORK, U, IFAULT)
       PRINT *, 'cluela returns #', IFAULT, ' Made ', K, ' clusters'

*        accuracy calculation
       CALL ARIND (Y, Z, N, G, K, VALTAB, ACCARI)
       CALL VINFO (Y, Z, N, G, K, VALTAB, ACCVI)
       PRINT *, 'residual: ', U, ' ARI: ', ACCARI, ' VI: ', ACCVI

*         write F values
       FNAME = 'wine_fvals_stand.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'wine: trouble opening ', FNAME
  80   FORMAT (I5, ',', G9.3)
       DO L = KMIN, KMAX
         WRITE (UNIT=IOU, FMT=80, IOSTAT=IFAULT) L, F(L)
         IF (IFAULT .NE. 0) PRINT *, 'wine: trouble writing F stats'
       END DO
       CLOSE (UNIT=IOU)

*=======================================================================
*           Now using specially obtained importance weights
*=======================================================================
       K = 3
       PRINT *, 'with special weights...'
       CALL CLUELA (X, NX, P, P, N, WSPEC, C, NC, KMIN, K, KMAX, Z, POP,
     $  F, ROBUST, SOSDO, IWORK, RWORK, U, IFAULT)
       PRINT *, 'cluela returns #', IFAULT, ' Made ', K, ' clusters'

*        accuracy calculation
       CALL ARIND (Y, Z, N, G, K, VALTAB, ACCARI)
       CALL VINFO (Y, Z, N, G, K, VALTAB, ACCVI)
       PRINT *, 'residual: ', U, ' ARI: ', ACCARI, ' VI: ', ACCVI

*         write F values
       FNAME = 'wine_fvals_weighted.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'wine: trouble opening ', FNAME
       DO L = KMIN, KMAX
         WRITE (UNIT=IOU, FMT=80, IOSTAT=IFAULT) L, F(L)
         IF (IFAULT .NE. 0) PRINT *, 'wine: trouble writing F stats'
       END DO

*-----------------------------------------------------------------------
*             project to 2-D
*-----------------------------------------------------------------------
       DO I = 1, N                            ! multiply in weights
         DO J = 1, P
           X(J,I) = X(J,I) * WSPEC(J)
         END DO
       END DO

       FNAME = 'wine_classes_weighted.csv'
       CALL PCA (X,P,P,N,ROBUST,R,N,M,Q,AVE,E,V,DEV,IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'pca: error #', IFAULT
       IF (Q < M) PRINT *, 'only ', Q, ' dimensions'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'wine: trouble opening ', FNAME
       CALL WCSVSP (R, N, M, Y, IOU, IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'wcsvsp: error #', IFAULT
       CLOSE (UNIT=IOU)

       FNAME = 'wine_classes_robust.csv'
       ROBUST = .TRUE.
       CALL PCA (X,P,P,N,ROBUST,R,N,M,Q,AVE,E,V,DEV,IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'pca: error #', IFAULT
       IF (Q < M) PRINT *, 'only ', Q, ' dimensions'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'wine: trouble opening ', FNAME
       CALL WCSVSP (R, N, M, Y, IOU, IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'wcsvsp: error #', IFAULT
       CLOSE (UNIT=IOU)

       STOP 'program complete'
      END  ! of WINE

