*=======================================================================
* test clustering on the credit approval data
* by Andy Allinger, 2017-2020, public domain
*=======================================================================
      PROGRAM CRX
       IMPLICIT NONE

*             fixed parameters
       INTEGER CATMAX, G, IOU, M, N, NS, O, P
       REAL XBOGUS
       PARAMETER (CATMAX = 14,     ! Most categories any attribute
     $            G = 2,           ! # correct classes
     $            IOU = 1,         ! I/O unit number
     $            XBOGUS = -9.,    ! filler value for missing data
     $            M = 2,           ! projection dimension for visualization
     $            N = 690,         ! Number of objects
     $            NS = 12,         ! node size of binary tree
     $            O = 15,          ! # input attributes
     $            P = 38)          ! # attibutes after encoding
       INTEGER NO
       PARAMETER (NO = N * O)

*             adjustable parameters
       INTEGER KMIN, KDEF, KMAX
       PARAMETER (KMIN = 1,             ! fewest clusters to make
     $            KDEF = 2,             ! default clusters to make
     $            KMAX = 20)            ! most clusters to make

*             variables in test program
       INTEGER
     $         CAT(O),             ! categories each input attribute
     $         FERR,               ! file i/o error
     $         I, J, JBEG, L,      ! counters
     $         IPERM(N),           ! index array
     $         JPERM(N),           ! index array
     $         LL,                 ! line length
     $         NX0(N),             ! workspace for IREID
     $         Y(N)                ! known classes

       REAL
     $      ACCARI, ACCVI,         ! accuracy vs. external standard
     $      S(CATMAX-1,CATMAX)     ! regular simplex coordinates

       LOGICAL FEXIST,             ! file exists
     $         ROBUST,             ! use L1 metric?
     $         SOSDO               ! optimize sum of square distance?
       CHARACTER CLASS(N), CH1
       CHARACTER*2 SYMBOL(N,O), HWORK
       CHARACTER*32 FNAME, RAW
       CHARACTER*256 TEXT

*             external functions
       INTEGER LTRIM
       LOGICAL SAFDIV

*             arguments to clustering subroutine
       INTEGER K, Z(N), POP(KMAX), IFAULT, NX(P,N), NC(P,KMAX),
     $  LIVEC(KMAX), LIVEX(N), TREE(NS,2*KMAX)
       REAL X(P,N), W(P), C(P,KMAX), U, SHORTC(KMAX), SHORTX(N),
     $  CRIT(4,2*KMAX), THRESH, WORK(N)
       LOGICAL MISS

*             backup of input, not to be recoded
       INTEGER NXBAK(P,N)
       REAL WBAK(P), XBAK(P,N)

*             workspace for stats
       INTEGER NA(P)
       REAL LEA(P), GRE(P), R

*             workspace for validation
       REAL VALTAB(G*KMAX+G+KMAX)

*             workspace for multi-dimensional scaling
       INTEGER Q
       REAL E(N), T(N), V(N,N)

       DATA FNAME / 'crx.data' /
*=======================================================================
* Credit card application data
* From the University of California at Irvine Machine Learning repository
*
* Reference:
*      Ross Quinlan
*      "Simplifying decision trees"
*      Int J Man-Machine Studies v.27, pp.221-234
*      Dec 1987
*
* This data was provided by an anonymous donor.  The attributes have been
* changed to arbitrary values.
*
* Classes are:
*    +, -
*
* Attributes are:
*
*___Attribute___Type___________#___Values_______________________________
*    1          Categorical    2   b, a
*    2          Continuous         13.75 ... 80.25
*    3          Continuous         0 ... 28.0
*    4          Categorical    4   u, y, l, t
*    5          Categorical    3   g, p, gg
*    6          Categorical   14   c, d, cc, i, j, k, m, r, q, w, x, e, aa, ff
*    7          Categorical    9   v, h, bb, j, n, z, dd, ff, o
*    8          Continuous         0 ... 28.5
*    9          Boolean        2   t, f
*   10          Boolean        2   t, f
*   11          Continuous         0 ... 67
*   12          Boolean        2   t, f
*   13          Categorical    3   g, p, s
*   14          Continuous         0 ... 2000
*   15          Continuous         0 ... 100,000
*=======================================================================

*        info above is summarized thus:
       DATA CAT / 2, 0, 0, 4, 3,         !  1-5
     $           14, 9, 0, 2, 2,         !  6-10
     $            0, 2, 3, 0, 0 /        ! 11-15

*             equal weights
       DATA W / P * 1. /
       DATA SYMBOL / NO * '$$' /         ! bogus initial value

*=======================================================================
*             begin
*=======================================================================
       CALL RNDINI ()               ! initialize random number generator
       INQUIRE (FILE=FNAME, EXIST=FEXIST, IOSTAT=FERR)  ! get input file
       IF (.NOT. FEXIST) STOP 'File not found'
       IF (FERR .NE. 0) STOP 'Error seeking input file'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='OLD',
     &        ACCESS='SEQUENTIAL', IOSTAT=FERR)
       IF (FERR .NE. 0) PRINT *, 'crx: trouble opening input file'

*=======================================================================
*             read in data
*=======================================================================
       DO I = 1, N
         READ (UNIT=IOU, FMT=10, END=30, IOSTAT=FERR) TEXT
         IF (FERR .NE. 0) PRINT *, 'trouble reading from file, line', I
  10     FORMAT (A)
         LL = LTRIM(TEXT)
         JBEG = 1
         L = 0
         DO 20 J = 1, LL                       ! split on commas
           IF (TEXT(J:J) .EQ. ',') THEN
             L = L + 1
             RAW = TEXT(JBEG:J-1)
           ELSE IF (J .EQ. LL) THEN
             L = L + 1
             RAW = TEXT(JBEG:J)
           ELSE
             GO TO 20
           END IF

           IF (L > O+1) STOP '!trouble reading L > O+1'
           IF (L .EQ. 16) THEN                           ! answer key
             CLASS(I) = RAW(1:1)
           ELSE IF (INDEX(RAW, '?') .EQ. 0) THEN          ! non-missing
             IF (CAT(L) .EQ. 0) THEN                   ! numerical?
               READ (UNIT=RAW, FMT=*, IOSTAT=FERR) X(L,I)
               IF (FERR .NE. 0) PRINT*,'error reading X, i: ',I,' L: ',L
             ELSE
               SYMBOL(I,L) = RAW(1:2)
             END IF
             NX(L,I) = 1
           ELSE                                 ! missing
             X(L,I) = XBOGUS
             NX(L,I) = 0
           END IF
           JBEG = J+1
  20     CONTINUE  ! next J
       END DO  ! next I
  30   CLOSE (UNIT=IOU)

*=======================================================================
* Preprocess input.
*=======================================================================
*             convert text labels to integers
       DO L = 1, O
         IF (CAT(L) > 0) THEN
           DO I = 1, N
             NX0(I) = NX(L,I)
           END DO
           CALL IREID (SYMBOL(1,L),Z,NX0,IPERM,JPERM,N,HWORK,IFAULT)
           IF (IFAULT .NE. 0) PRINT *, 'ireid: error #', IFAULT
           DO I = 1, N
             X(L,I) = FLOAT(Z(I))
           END DO
         END IF
       END DO

*-----------------------------------------------------------------------
* Variable 15 has a standard deviation greater than its semi-interquartile
* range by a factor of 26.  Apply a logarithm.
*-----------------------------------------------------------------------
       DO I = 1, N
          X(15,I) = LOG(1. + X(15,I))
        END DO

*             rescale by range
       CALL LEAGRE (X, NX, P, P, N, NA, LEA, GRE)
       DO L = 1, O
         IF (CAT(L) .EQ. 0) THEN
           R = GRE(L) - LEA(L)
           IF (SAFDIV(1., R)) THEN
             R = 1. / R
           ELSE
             R = 0.
             PRINT *, 'crx: no variance column ', L
           END IF
           DO I = 1, N
             X(L,I) = (X(L,I) - LEA(L)) * R
           END DO
         END IF
       END DO

*             back up input
       DO I = 1, N
         DO L = 1, O
           XBAK(L,I) = X(L,I)
           NXBAK(L,I) = NX(L,I)
         END DO
       END DO
       DO L = 1, O
         WBAK(L) = W(L)
       END DO

*             recode categoricals with corresponding column of simplex
       ROBUST = .FALSE.
       CALL RSCODE (X,NX,P,P,N,W,CAT,O,S,CATMAX,ROBUST,IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'rscode: error #', IFAULT

*             relabel class ID's
       DO I = 1, N
         NX0(I) = 1
       END DO
       CALL IREID (CLASS, Y, NX0, IPERM, JPERM, N, CH1, IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'ireid on class, error #', IFAULT

*=======================================================================
*              call the subroutine
*=======================================================================
       K = KDEF
       SOSDO = .TRUE.
       MISS = .TRUE.
       CALL DINIT (X,NX,P,P,N,W,C,NC,K,ROBUST,SHORTX,WORK,THRESH)
       CALL KCLUST (X, NX, P, P, N, W, C, NC, KMIN, K, KMAX, Z, POP,
     $  MISS, ROBUST, SOSDO, THRESH, LIVEC, LIVEX, IPERM, TREE, CRIT,
     $  SHORTC, SHORTX, WORK, U, IFAULT)
       PRINT *, 'kclust returns #', IFAULT, ' made ', K, ' clusters'

*             try reducing the threshhold
       CALL DINIT (X,NX,P,P,N,W,C,NC,K,ROBUST,SHORTX,WORK,THRESH)
       THRESH = THRESH / 2.
       CALL KCLUST (X, NX, P, P, N, W, C, NC, KMIN, K, KMAX, Z, POP,
     $  MISS, ROBUST, SOSDO, THRESH, LIVEC, LIVEX, IPERM, TREE, CRIT,
     $  SHORTC, SHORTX, WORK, U, IFAULT)
       PRINT *, 'with half the threshhold distance, kclust returns #',
     $  IFAULT, ' made ', K, ' clusters'

*            find accuracy
       IF (IFAULT.EQ. 0) THEN
         CALL ARIND (Y, Z, N, G, K, VALTAB, ACCARI)  ! adjusted Rand index
         CALL VINFO (Y, Z, N, G, K, VALTAB, ACCVI)   ! variation of information
       END IF

*             project to 2-D
       CALL MDS (XBAK,NXBAK,P,O,N,WBAK,CAT,ROBUST,V,N,Q,E,T,IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'mds: error #', IFAULT
       IF (Q < M) PRINT *, 'mds: only ', Q, ' directions!'

*-----------------------------------------------------------------------
* Write CSV file with coordinates in different columns by cluster
*-----------------------------------------------------------------------
       FNAME = 'crx_clusters.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'crx: trouble opening ', FNAME
       CALL WCSVSP (V, N, M, Z, IOU, IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'wcsvsp: error #', IFAULT
       CLOSE (UNIT=IOU)

*-----------------------------------------------------------------------
* Write CSV file with coordinates in different columns by class
*-----------------------------------------------------------------------
       FNAME = 'crx_classes_euclid.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'crx: trouble opening ', FNAME
       CALL WCSVSP (V, N, M, Y, IOU, IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'wcsvsp: error #', IFAULT
       CLOSE (UNIT=IOU)

*-----------------------------------------------------------------------
* Robust projection
*-----------------------------------------------------------------------
       ROBUST = .TRUE.
       CALL MDS (XBAK,NXBAK,P,O,N,WBAK,CAT,ROBUST,V,N,Q,E,T,IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'mds: error #', IFAULT
       IF (Q < M) PRINT *, 'mds: only ', Q, ' directions!'

       FNAME = 'crx_classes_robust.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'crx: trouble opening ', FNAME
       CALL WCSVSP (V, N, M, Y, IOU, IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'wcsvsp: error #', IFAULT
       CLOSE (UNIT=IOU)

       PRINT *, 'Returned residual: ', U
       PRINT *, 'Adjusted Rand Index: ', ACCARI
       PRINT *, ' Variation of Information:', ACCVI
       STOP 'program complete'
      END  ! of t_crx.f ************************************************
