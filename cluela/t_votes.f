*-----------------------------------------------------------------------
* test clustering on the House votes '84 data set
* by Andy Allinger, 2011-2017, public domain
*-----------------------------------------------------------------------
      PROGRAM VOTES
       IMPLICIT NONE

*             fixed parameters
       INTEGER G, IOU, M, NLIM, PLIM
       REAL XBOGUS
       PARAMETER (
     $            G = 2,         ! # correct classes
     $            IOU = 1,       ! I/O unit number
     $            XBOGUS = -9.,  ! filler value for missing data
     $            M = 2,         ! projection dimension for visualization
     $            NLIM = 435,    ! Number of objects
     $            PLIM = 16)     ! # attibutes

*             adjustable parameters
       INTEGER KMIN, KMAX
       PARAMETER (KMIN = 1,            ! fewest clusters to make
     $            KMAX = 20)           ! most clusters to make

*             variables in test program
       INTEGER
     $         FERR,                     ! file i/o error
     $         I, J, L,                  ! counters
     $         LL,
     $         Y(NLIM)                   ! known classes

       REAL
     $      ACCARI, ACCVI,         ! accuracy vs. external standard
     $      R(NLIM,M)              ! projected coordinates

       LOGICAL FEXIST,            ! file exists
     $         ROBUST,            ! use L1 statistics
     $         SOSDO              ! minimize sum-of-squares error
       CHARACTER RAW(PLIM,NLIM)   ! input data
       CHARACTER*32 FNAME
       CHARACTER*256 TEXT

*             external function
       INTEGER LTRIM

*             arguments to clustering subroutine
       INTEGER K, N, P, Z(NLIM), POP(KMAX), IWORK(4*NLIM+25*KMAX),
     $         NX(PLIM,NLIM), NC(PLIM,KMAX), IFAULT
       REAL X(PLIM,NLIM), W(PLIM), C(PLIM,KMAX), RWORK(2*NLIM+9*KMAX)
       REAL F(KMIN:KMAX), U

*             workspace for validation
       REAL VALTAB(G*KMAX+G+KMAX)

*             arguments for principal component analysis
       INTEGER Q
       REAL A(PLIM), E(PLIM), V(PLIM,PLIM), T(PLIM)

       DATA W / PLIM * 1. /
       DATA FNAME / 'house-votes-84.data' /
*-----------------------------------------------------------------------
* 1984 United States Congressional Voting Records Database
* From the University of California at Irvine Machine Learning repository
*
* Source:
*      Congressional Quarterly Almanac, 98th Congress,
*      2nd session 1984, Volume XL: Congressional Quarterly Inc.
*      Washington, D.C., 1985.
*
* Reference:
*      J. C. Schlimmer
*      Concept acquisition through representational adjustment
*      Doctoral dissertation, Department of Information and Computer Science,
*      University of California, Irvine, CA., 1987.
*
* Classes are:
*     democrat, republican
*
* Attributes are:
*
*___Attribute________________________________Type______#___Values_______
*   handicapped-infants                      Boolean   2   y,n,?
*   water-project-cost-sharing               Boolean   2   y,n,?
*   adoption-of-the-budget-resolution        Boolean   2   y,n,?
*   physician-fee-freeze                     Boolean   2   y,n,?
*   el-salvador-aid                          Boolean   2   y,n,?
*   religious-groups-in-schools              Boolean   2   y,n,?
*   anti-satellite-test-ban                  Boolean   2   y,n,?
*   aid-to-nicaraguan-contras                Boolean   2   y,n,?
*   mx-missile                               Boolean   2   y,n,?
*   immigration                              Boolean   2   y,n,?
*   synfuels-corporation-cutback             Boolean   2   y,n,?
*   education-spending                       Boolean   2   y,n,?
*   superfund-right-to-sue                   Boolean   2   y,n,?
*   crime                                    Boolean   2   y,n,?
*   duty-free-exports                        Boolean   2   y,n,?
*   export-administration-act-south-africa   Boolean   2   y,n,?
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*          Begin
*-----------------------------------------------------------------------
       P = PLIM
       N = NLIM
       CALL RNDINI ()
       INQUIRE (FILE=FNAME, EXIST=FEXIST, IOSTAT=FERR)
       IF (FERR .NE. 0) STOP 'Error inquiring of input file'
       IF (.NOT. FEXIST) STOP 'File not found'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='OLD', ACCESS='SEQUENTIAL',
     $        IOSTAT=FERR)
       IF (FERR .NE. 0) STOP 'trouble opening input file'

*-----------------------------------------------------------------------
*        Read in data
*-----------------------------------------------------------------------
       DO I = 1, N
         READ (UNIT=IOU, FMT=10, END=30, IOSTAT=FERR) TEXT
         IF (FERR .NE. 0) PRINT *, 'trouble reading line ', I
  10     FORMAT (A)
         IF (TEXT(1:1) .EQ. 'd') Y(I) = 1       ! class labels
         IF (TEXT(1:1) .EQ. 'r') Y(I) = 2
         LL = LTRIM(TEXT)
         L = 0
         DO J = 1, LL-1
           IF (TEXT(J:J) .EQ. ',') THEN
             L = L + 1
             RAW(L,I) = TEXT(J+1:J+1)
           END IF
         END DO ! next j
       END DO
  30   CLOSE (UNIT = IOU)

*-----------------------------------------------------------------------
* Consider the problem as binary data with some values missing.
*-----------------------------------------------------------------------
       DO I = 1, N
         DO J = 1, P
           IF (RAW(J,I) .EQ. 'y') THEN
             X(J,I) = 1.
             NX(J,I) = 1
           ELSE IF (RAW(J,I) .EQ. 'n') THEN
             X(J,I) = 0.
             NX(J,I) = 1
           ELSE IF (RAW(J,I) .EQ. '?') THEN
             X(J,I) = XBOGUS
             NX(J,I) = 0
           ELSE
             STOP 'votes: trouble recoding input!'
           END IF
         END DO
       END DO

*-----------------------------------------------------------------------
*              Call the subroutine
*-----------------------------------------------------------------------
        K = 3
        ROBUST = .FALSE.
        SOSDO = .TRUE.
        CALL CLUELA (X, NX, PLIM, P, N, W, C, NC, KMIN, K, KMAX, Z, POP,
     $         F, ROBUST, SOSDO, IWORK, RWORK, U, IFAULT)
          PRINT *, 'cluela returns #', IFAULT, ' made ', K, ' clusters'

*-----------------------------------------------------------------------
* Reconsider the problem as scalar data, with no missing values.
*-----------------------------------------------------------------------
       DO I = 1, N
         DO J = 1, P
           IF (RAW(J,I) .EQ. 'y') X(J,I) = 1.0
           IF (RAW(J,I) .EQ. 'n') X(J,I) = 0.0
           IF (RAW(J,I) .EQ. '?') X(J,I) = 0.5
           NX(J,I) = 1
         END DO
       END DO

*-----------------------------------------------------------------------
*              Call the subroutine
*-----------------------------------------------------------------------
       K = 3
       CALL CLUELA (X, NX, PLIM, P, N, W, C, NC, KMIN, K, KMAX, Z, POP,
     $         F, ROBUST, SOSDO, IWORK, RWORK, U, IFAULT)
       PRINT *, 'cluela returns #', IFAULT, ' made ', K, ' clusters'

*            Find accuracy
       IF (IFAULT.EQ. 0) THEN
         CALL ARIND (Y, Z, N, G, K, VALTAB, ACCARI)  ! adjusted Rand index
         CALL VINFO (Y, Z, N, G, K, VALTAB, ACCVI)   ! variation of information
       END IF

*             project to 2-D
       CALL PCA (X,PLIM,P,N,ROBUST,R,NLIM,M,Q,A,E,V,T,IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'pca: error #', IFAULT
       IF (Q < M) PRINT *, 'pca: only ', Q, ' directions!'

*-----------------------------------------------------------------------
*             Write output
*-----------------------------------------------------------------------
       FNAME = 'votes_clusters_sosdo.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'votes: trouble opening ', FNAME
       CALL WCSVSP (R, N, M, Z, IOU, IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'wcsvsp: error #', IFAULT
       CLOSE (UNIT=IOU)

       FNAME = 'votes_classes.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'votes: trouble opening ', FNAME
       CALL WCSVSP (R, N, M, Y, IOU, IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'wcsvsp: error #', IFAULT
       CLOSE (UNIT=IOU)

*         Write F values
       FNAME = 'votes_fvals.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'votes: trouble opening ', FNAME
  80   FORMAT (I5, ',', G9.3)
       DO L = KMIN, KMAX
         WRITE (UNIT=IOU, FMT=80, IOSTAT=IFAULT) L, F(L)
         IF (IFAULT .NE. 0) PRINT *, 'votes: trouble writing F stats'
       END DO
       CLOSE (UNIT=IOU)

       PRINT *, 'Returned residual: ', U
       PRINT *, 'Adjusted Rand Index: ', ACCARI
       PRINT *, ' Variation of Information:', ACCVI

*-----------------------------------------------------------------------
*             Repeat in robust mode
*-----------------------------------------------------------------------
       ROBUST = .TRUE.
       SOSDO = .FALSE.
       K = 3
       CALL CLUELA (X, NX, PLIM, P, N, W, C, NC, KMIN, K, KMAX, Z, POP,
     $         F, ROBUST, SOSDO, IWORK, RWORK, U, IFAULT)
       PRINT *, 'cluela returns #', IFAULT, ' made ', K, ' clusters'

       IF (IFAULT.EQ. 0) THEN
         CALL ARIND (Y, Z, N, G, K, VALTAB, ACCARI)  ! adjusted Rand index
         CALL VINFO (Y, Z, N, G, K, VALTAB, ACCVI)   ! variation of information
       END IF

       CALL PCA (X,PLIM,P,N,ROBUST,R,NLIM,M,Q,A,E,V,T,IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'pca: error #', IFAULT
       IF (Q < M) PRINT *, 'pca: only ', Q, ' directions!'

       FNAME = 'votes_clusters_robust.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'votes: trouble opening ', FNAME
       CALL WCSVSP (R, N, M, Z, IOU, IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'wcsvsp: error #', IFAULT
       CLOSE (UNIT=IOU)

       PRINT *, 'Returned residual: ', U
       PRINT *, 'Adjusted Rand Index: ', ACCARI
       PRINT *, ' Variation of Information:', ACCVI
       STOP 'program complete'
      END  ! of votes
