*=======================================================================
* test clustering on the horse colic data
* by Andy Allinger, 2013-2020, public domain
*=======================================================================
      PROGRAM HORSE
       IMPLICIT NONE

*             fixed parameters
       INTEGER CATMAX, G, IOU, MVIS, N, PIN, O, P, PLIM, Q, QN
       REAL WINIT, XBOGUS
       PARAMETER (CATMAX = 346,       ! Most categories any attribute
     $            G = 62,             ! number of distinct diagnoses
     $            IOU = 1,            ! I/O unit number
     $            MVIS = 2,           ! projection dimension for visualization
     $            N = 368,            ! Number of objects
     $            PIN = 25,           ! # input attributes
     $            O = 86,             ! # after recoding multi-attribute
     $            P = 440,            ! # attibutes after recoding categoricals
     $            Q = 3,              ! # concurrent diagnoses
     $            XBOGUS = -9.)       ! filler value for missing data
       PARAMETER (PLIM = P,           ! row bound of X
     $            QN = Q * N,         ! size of array GNOSIS
     $            WINIT = 1.0 / PIN)  ! weights to sum to one

*             adjustable parameters
       INTEGER ITER, KMIN, KMAX, KDEF
       PARAMETER (ITER = 1000,          ! iterations allowed in KMNS
     $            KMIN = 1,             ! fewest clusters to make
     $            KDEF = 2,             ! default clusters
     $            KMAX = 20)            ! most clusters to make

*             parameters for simulated annealing
       INTEGER NIT, PEX
       REAL SBEG, SEND, TBEG, TEND, SMIDGE
       PARAMETER (NIT = 10 000,          ! iterations
     $            PEX = 3,               ! penalty exponent
     $            SBEG = 0.1,            ! jitter amount, beginning
     $            SEND = 0.001,          ! jitter amount, end
     $            TBEG = 1.E-3,          ! temperature, beginning
     $            TEND = 1.E-6,          ! temperature, end
     $            SMIDGE = 1. / PIN**2)  ! minimal weight

*             variables in test program
       INTEGER
     $         CAT(O),                ! categories each before calling RSCODE
     $         FERR,                  ! file i/o error
     $         GNOSIS(Q,N),           ! diagnoses
     $         HOSPTL(N),             ! hospital ID number
     $         H, HBEG, I, J, L,      ! counters
     $         IPERM(QN),             ! index array
     $         JPERM(QN),             ! index array
     $         KBEST,                 !
     $         LL,                    ! line length
     $         NX0(QN)                ! workspace for JSCODE

       REAL
     $      BIG,                   ! arbitrary large number
     $      R,                     ! temporary scalar
     $      S(CATMAX-1,CATMAX),    ! regular simplex coordinates
     $      THRESH,                ! threshhold distance for kclust
     $      UBEST,                 ! lowest residual from KMNS
     $      X0(N)                  ! workspace for ORDNAL

       LOGICAL
     $         FEXIST,      ! file exists
     $         ORDREQ(PIN)  ! ordinal rescaling required

       CHARACTER*32 FNAME, RAW
       CHARACTER*256 TEXT

*             arguments to clustering subroutine
       INTEGER K, POP(KMAX), IWORK(4*N+25*KMAX), NX(PLIM,N), NXBAK(O,N),
     $          NC(PLIM,KMAX), Z(N), ZBEST(N), IFAULT
       REAL X(PLIM,N), XBAK(O,N), W(PLIM), WIN(PIN), C(PLIM,KMAX),
     $       RWORK(2*N+9*KMAX), F(KMIN:KMAX), U, WBAK(O)
       LOGICAL MISS, ROBUST, SOSDO

*             arguments for simulated annealing
       INTEGER IFLAG,
     $         LIS                    ! last improving step
       REAL PNLTY, PNLMIN, WALL, WKEEP(PIN), WBEST(PIN)

*             workspace for KMNS
       INTEGER NCP(KDEF), ITRAN(KDEF), LIVE(KDEF)
       REAL CT(KDEF,N), AN1(KDEF), AN2(KDEF), WSS(KDEF)

*             workspace for stats
       INTEGER NA(PIN)
       REAL AVE(PIN), DEV(PIN)

*             workspace for multi-dimensional scaling
       INTEGER M
       REAL E(N), V(N,N), T(N)

*             timing info
       REAL TARRAY(2), RESULT

*             external functions
       INTEGER LTRIM
       LOGICAL ALMEQ, SAFDIV

*             initial values
       DATA FNAME / 'horse-colic.both' /
*=======================================================================
* Horse colic data
* From the University of California at Irvine Machine Learning repository
*
* Creators: Mary McLeish & Matt Cecile
*           Department of Computer Science
*           University of Guelph
*           Guelph, Ontario, Canada N1G 2W1
*           mdmcleish@water.waterloo.edu
*
* Classes: As the analyst may choose.  Possibilities are outcome
*          (lived / died / was euthanized) or type of lesion (one of 62
*           diagnosis numbers).
*
* Attributes are:
*
*________Attribute__________Type____________#____Values__________
*    1   surgery            Boolean         2   1:  yes
*                                               2:  no
*
*    2   age                Categorical     2   1:  adult
*                                               9:  young     ! 9?
*
*    3   hospital number    Categorical   346   ID number assigned to horse
* If the horse lives, then this might be considered a categorical
* variable.  Since a horse only dies once, for horses in the groups
* Died and Was Euthanized, this is a contra-categorical variable, since
* the presence of a horse in one of these groups means that no other
* instances of this horse may be present.
*
*    4   rectal temp        Continuous          [degrees Celsius]
*
*    5   pulse              Continuous          [beats per minute]
*
*    6   respiratory rate   Continuous
*
*    7   extremities temp   Ordinal             1:  Normal
*                                               2:  Warm
*                                               3:  Cool
*                                               4:  Cold
*                                                                Wrong order!
*
*    8   peripheral pulse   Ordinal             1:  normal
*                                               2:  increased
*                                               3:  reduced
*                                               4:  absent
*                                                                Wrong order!
*
*    9   mucous membranes   Categorical     6   1:  normal pink
*                                               2:  bright pink
*                                               3:  pale pink
*                                               4:  pale cyanotic
*                                               5:  bright red / injected
*                                               6:  dark cyanotic
*
*   10   capillary refill   Categorical     3   1:  < 3 seconds
*                                               2:  >= 3 seconds
*                                               3:  undocumented!
*
*   11   pain               Continuous          1:  alert, no pain
*                                               2:  depressed
*                                               3:  intermittent mild pain
*                                               4:  intermittent severe pain
*                                               5:  continuous severe pain
* The file 'horse-colic.names' says: "should NOT be treated as a ordered or
* discrete variable!"  Obey this instruction and treat pain as continuous.
*
*   12   peristalsis        Ordinal             1:  hypermotile
*                                               2:  normal
*                                               3:  hypomotile
*                                               4:  absent
*
*   13   abdom distension   Ordinal             1:  none
*                                               2:  slight
*                                               3:  moderate
*                                               4:  severe
*
*   14   nasogastric tube   Ordinal             1:  none
*                                               2:  slight
*                                               3:  significant
*
*   15   naso. reflux       Ordinal             1:  none
*                                               2:  > 1 liter
*                                               3:  < 1 liter
*                                                                Wrong order!
*
*   16   reflux PH          Continuous          0 ... 14
*
*   17   rectal feces       Ordinal             1:  normal
*                                               2:  increased
*                                               3:  decreased
*                                               4:  absent
*                                                                Wrong order!
*
*   18   abdomen            Categorical     5   1:  normal
*                                               2:  other
*                                               3:  firm feces large intestine
*                                               4:  distended small intestine
*                                               5:  distended large intestine
*
*   19   packed cell vol.   Continuous          red cells by volume in blood
*
*   20   total protein      Continous           [gms/dL]
*
*   21   abdominocentesis   Categorical     3   1:  clear
*                                               2:  cloudy
*                                               3:  serosanguinous
*
*   22   abdom. protein     Continuous          [gms/dL]
*
*   23   outcome            Categorical     3   1:  lived
*                                               2:  died
*                                               3:  was euthanized
*
*   24   surgical lesion?   Categorical     2   1:  Yes
*                                               2:  No
*
*   25   type of lesion     Multi          62
*   26
*   27
* Three columns contain either 00000 or one of 62 distinct numeric codes
* for a diagnosis.
*
*  28    cp_data            Boolean             1:  Yes
*                                               2:  No
* This column is irrelevant and shall be omitted.
*
*=======================================================================

*             info above is summarized thus:
       DATA CAT / 2, 2, 346, 0, 0, 0, 0, 0, 6, 3,  !  1 ... 10
     $            0, 0, 0, 0, 0, 0, 0, 5, 0, 0,    ! 11 ... 20
     $            3, 0, 3, 2,                      ! 21 ... 24
     $  G * 2 /                                    ! 25 ... 86

       DATA ORDREQ / .FALSE., .FALSE., .FALSE., .FALSE., .FALSE.,
     $               .FALSE.,  .TRUE.,  .TRUE., .FALSE., .FALSE.,
     $               .FALSE.,  .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,
     $               .FALSE.,  .TRUE., .FALSE., .FALSE., .FALSE.,
     $               .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. /

       DATA WIN / PIN * WINIT /      ! equal weights
       DATA NX0 / QN * 1 /
       DATA BIG / 1.E27 /
       DATA ROBUST / .FALSE. /

*=======================================================================
*             begin
*=======================================================================
       CALL RNDINI ()             ! initialize random number generator
       INQUIRE (FILE=FNAME, EXIST=FEXIST, IOSTAT=FERR)  ! get input file
       IF (.NOT. FEXIST) STOP 'File not found'
       IF (FERR .NE. 0) STOP 'Error seeking input file'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='OLD',
     $        ACCESS='SEQUENTIAL', IOSTAT=FERR)
       IF (FERR .NE. 0) PRINT *, 'horse: trouble opening input file'

*=======================================================================
*             Read input data
*=======================================================================
       DO I = 1, N
  10     FORMAT (A)
         READ (UNIT=IOU, FMT=10, END=30, IOSTAT=FERR) TEXT
         IF (FERR .NE. 0) PRINT *, 'trouble reading from file, line', I
         LL = LTRIM(TEXT)
         HBEG = 1
         L = 0
         J = 0
         DO 20 H = 1, LL
           IF (TEXT(H:H) .EQ. ' ') THEN          ! split on spaces
             L = L + 1
             RAW = TEXT(HBEG:H-1)
           ELSE IF (H .EQ. LL) THEN
             L = L + 1
             RAW = TEXT(HBEG:H)
           ELSE
             GO TO 20
           END IF

           IF (L .EQ. 3) THEN  ! ID number
             J = J + 1
             NX(J,I) = 1       ! nothing missing
             READ (UNIT=RAW, FMT=*, IOSTAT=FERR) HOSPTL(I)
             IF (FERR .NE. 0) PRINT *, 'error reading ID line ', I
           ELSE IF (L .EQ. 25) THEN
             READ (UNIT=RAW, FMT=*, IOSTAT=FERR) GNOSIS(1,I)
             IF (FERR .NE. 0) PRINT *, 'error read attrib1 object ', I
           ELSE IF (L .EQ. 26) THEN
             READ (UNIT=RAW, FMT=*, IOSTAT=FERR) GNOSIS(2,I)
             IF (FERR .NE. 0) PRINT *, 'error read attrib2 object ', I
           ELSE IF (L .EQ. 27) THEN
             READ (UNIT=RAW, FMT=*, IOSTAT=FERR) GNOSIS(3,I)
             IF (FERR .NE. 0) PRINT *, 'error read attrib3 object ', I
           ELSE IF (L .EQ. 28) THEN            ! ignore cp_data
             CONTINUE
           ELSE IF (INDEX(RAW, '?') .EQ. 0) THEN      ! data present?
             J = J + 1
             READ (UNIT=RAW, FMT=*, IOSTAT=FERR) X(J,I)
             IF (FERR .NE. 0) PRINT *,'error reading X, i: ', I,' L:', L
             NX(J,I) = 1
           ELSE                                  ! missing
             J = J + 1
             X(J,I) = XBOGUS
             NX(J,I) = 0
           END IF
           HBEG = H + 1
  20     CONTINUE  ! next H
       END DO  ! next i
  30   CLOSE (UNIT=IOU)

*=======================================================================
*             preprocess data
*=======================================================================
*             relabel the young / adult variable
       DO I = 1, N
         IF (X(2,I) > 1.5) X(2,I) = 2.
       END DO

*             relabel case ID's
       CALL JSCODE (HOSPTL, NX0, IPERM, JPERM, N)
       DO I = 1, N
         X(3,I) = FLOAT(HOSPTL(I))
       END DO

*             fix wrong ordered ordinals
       DO I = 1, N
         IF (ALMEQ(X(7,I), 1.)) THEN   ! put subjective temperature in order
           X(7,I) = 2.
         ELSE IF (ALMEQ(X(7,I), 2.)) THEN
           X(7,I) = 1.
         END IF

         IF (ALMEQ(X(8,I), 1.)) THEN   ! put subjective pulse in order
           X(8,I) = 2.
         ELSE IF (ALMEQ(X(8,I), 2.)) THEN
           X(8,I) = 1.
         END IF

         IF (ALMEQ(X(15,I), 2.)) THEN            ! put reflux in order
           X(15,I) = 3.
         ELSE IF (ALMEQ(X(15,I), 3.)) THEN
           X(15,I) = 2.
         END IF

         IF (ALMEQ(X(17,I), 1.)) THEN            ! put feces in order
           X(17,I) = 2.
         ELSE IF (ALMEQ(X(17,I), 2.)) THEN
           X(17,I) = 1.
         END IF
       END DO

       IF (ROBUST) THEN
         CALL MEDMAD (X, NX, PLIM, PIN, N, X0, AVE, DEV)
       ELSE
         CALL AVEDEV (X, NX, PLIM, PIN, N, NA, AVE, DEV)
       END IF

       DO J = 1, PIN
         IF (ORDREQ(J)) THEN                     ! rescale ordinals
           DO I = 1, N
             X0(I) = X(J,I)
             NX0(I) = NX(J,I)
           END DO
           CALL ORDNAL (X0, NX0, IPERM, JPERM, N)
           DO I = 1, N
             X(J,I) = X0(I)
           END DO
         ELSE IF (CAT(J) .EQ. 0) THEN            ! standardize continuous
           IF (SAFDIV(1., DEV(J))) THEN
             R = 1. / DEV(J)
           ELSE
             PRINT *, 'variable ',J, ': no variance'
             R = 0.
           END IF
           DO I = 1, N
             IF (NX(J,I) > 0) X(J,I) = (X(J,I) - AVE(J)) * R
           END DO
         END IF
       END DO

*             convert diagnosis codes to small consecutive integers
       DO I = 1, QN
         NX0(I) = 1
       END DO
       CALL JSCODE (GNOSIS, NX0, IPERM, JPERM, QN)
       DO I = 1, N
         DO H = 1, Q
           GNOSIS(H,I) = GNOSIS(H,I) - 1         ! zero is null attribute!
         END DO
       END DO

*             Recode the multi-attribute
       DO I = 1, N
         DO J = 25, O
           X(J,I) = 2.                           ! no
           NX(J,I) = 1
         END DO
       END DO
       DO I = 1, N
         DO H = 1, Q
           IF (GNOSIS(H,I) .NE. 0) THEN
             J = 24 + GNOSIS(H,I)
             X(J,I) = 1.                         ! yes
           END IF
         END DO
       END DO

*             Back up input.
       DO I = 1, N
         DO J = 1, O
           XBAK(J,I) = X(J,I)
           NXBAK(J,I) = NX(J,I)
         END DO
       END DO
       DO J = 1, PIN-1
         W(J) = WIN(J)
       END DO
       DO J = PIN, O              ! spread weight across multi-attribute
         W(J) = WINIT / FLOAT(G)
       END DO
       DO J = 1, O
         WBAK(J) = W(J)
       END DO

*             recode categoricals with corresponding column of simplex
       CALL RSCODE (X,NX,PLIM,P,N,W,CAT,O,S,CATMAX,ROBUST,IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'rscode: error #', IFAULT

*=======================================================================
*             call CLUELA on the full problem - equal weights
*=======================================================================
       PRINT *, 'Original problem:'
       K = KDEF
       SOSDO = .NOT. ROBUST
       CALL DTIME (TARRAY, RESULT)
       CALL CLUELA (X, NX, PLIM, P, N, W, C, NC, KMIN, K, KMAX, Z,
     $  POP, F, ROBUST, SOSDO, IWORK, RWORK, U, IFAULT)
       CALL DTIME (TARRAY, RESULT)
       PRINT *, 'cluela returns #', IFAULT, ' in ', RESULT, ' seconds',
     $  ' made ', K, ' clusters with residual ', U

*         write F values
       FNAME = 'horse_fvals_equal.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=FERR)
       IF (FERR .NE. 0) PRINT *, 'horse: trouble opening ', FNAME
  40   FORMAT (I5, ',', G9.3)
       DO L = KMIN, KMAX
         WRITE (UNIT=IOU, FMT=40, IOSTAT=FERR) L, F(L)
         IF (FERR .NE. 0) PRINT *, 'horse: trouble writing F stats'
       END DO
       CLOSE (UNIT=IOU)

*             project to lower dimension
       CALL MDS (XBAK,NXBAK,O,O,N,WBAK,CAT,ROBUST,V,N,M,E,T,IFAULT)
       PRINT *, 'mds returns #', IFAULT, ' with ', M, ' directions'

*             Scatter plot of CLUELA result
       FNAME = 'horse_clusters_equal.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=FERR)
       IF (FERR .NE. 0) PRINT *, 'horse: trouble opening ', FNAME
       CALL WCSVSP (V, N, MVIS, Z, IOU, IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'wcsvsp: error #', IFAULT
       CLOSE (UNIT=IOU)

*=======================================================================
*             Call KCLUST directly.  Optimize weights
*=======================================================================
       MISS = .TRUE.
       IFLAG = 1
       PNLMIN = BIG
       PNLTY = 0.
       LIS = 0
       DO L = 0, NIT
         K = KDEF
         CALL OSIMAN (PIN, WIN, PNLTY,
     $                  SBEG, SEND, TBEG, TEND, NIT, WKEEP, IFLAG)
         DO I = 1, N
           DO J = 1, O
             X(J,I) = XBAK(J,I)
             NX(J,I) = NXBAK(J,I)
           END DO
         END DO
         WALL = 0.
         DO J = 1, PIN
           IF (WIN(J) < SMIDGE) WIN(J) = SMIDGE
           WALL = WALL + WIN(J)
         END DO
         DO J = 1, PIN
           WIN(J) = WIN(J) / WALL
           W(J) = WIN(J)
         END DO
         DO J = PIN, O
           W(J) = WIN(PIN) / FLOAT(G)
         END DO

         PNLTY = 0.
         DO J = 1, PIN              ! weight penalty
           PNLTY = PNLTY + WIN(J)**PEX
         END DO
         CALL RSCODE (X,NX,PLIM,P,N,W,CAT,O,S,CATMAX,ROBUST,IFAULT)
         IF (IFAULT .NE. 0) PRINT *, 'rscode: error #', IFAULT
         CALL DINIT (X, NX, PLIM, P, N, W, C, NC, K, ROBUST,
     $               RWORK, RWORK(N+1), THRESH)
         CALL KCLUST (X, NX, PLIM, P, N, W, C, NC, KMIN, K, KMAX, Z,
     $     POP, MISS, ROBUST, SOSDO, BIG, IWORK(1), IWORK(KMAX+1),
     $     IWORK(KMAX+N+1), IWORK(KMAX+2*N+1), RWORK(KMAX+2*N+1),
     $     RWORK(1), RWORK(KMAX+1), RWORK(KMAX+N+1), U, IFAULT)
         PNLTY = PNLTY * U
  50     FORMAT ('iter ', I7, ' kclust returns #', I3, '.  Made ', I3,
     $     ' clusters with residual ', E11.3, ' penalty: ', E11.3, $)
  60     FORMAT (' <---- New low penalty')
  70     FORMAT ()
         PRINT 50, L, IFAULT, K, U, PNLTY
         IF (IFAULT .EQ. 0 .AND. K .EQ. KDEF) THEN
           IF (PNLTY < PNLMIN) THEN
             PNLMIN = PNLTY
             UBEST = U
             DO I = 1, N
               ZBEST(I) = Z(I)
             END DO
             DO J = 1, PIN
               WBEST(J) = WIN(J)
             END DO
             PRINT 60
             LIS = L
           ELSE
             PRINT 70
           END IF
         ELSE
           PNLTY = BIG
           PRINT 70
         END IF
       END DO
  80   FORMAT (25F7.3)
       PRINT *, 'penalty exponent ', PEX, ' and robust=', ROBUST,
     $  ' yield optimized weights: '
       PRINT 80, WBEST
       PRINT *, 'best value iter ', LIS, ' ubest: ', UBEST

*=======================================================================
*             project to lower dimension
*=======================================================================
       DO J = 1, PIN-1
         W(J) = WBEST(J)
       END DO
       DO J = PIN, O
         W(J) = WBEST(PIN) / FLOAT(G)
       END DO
       CALL MDS (XBAK,NXBAK,O,O,N,W,CAT,ROBUST,V,N,M,E,T,IFAULT)
       PRINT *, 'mds returns #', IFAULT, ' with ', M, ' directions'

*             Scatter plot of KCLUST result
       FNAME = 'horse_clusters_kclust.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=FERR)
       IF (FERR .NE. 0) PRINT *, 'horse: trouble opening ', FNAME
       CALL WCSVSP (V, N, MVIS, ZBEST, IOU, IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'wcsvsp: error #', IFAULT
       CLOSE (UNIT=IOU)

*             Transpose V into X.  There is no longer any missing data.
*             Remember that dimension is now M.
       DO I = 1, N
         DO J = 1, M
           X(J,I) = V(I,J)
           NX(J,I) = 1
         END DO
       END DO
       DO J = 1, M                ! remember to change W
         W(J) = 1.
       END DO
       ROBUST = .FALSE.           ! not for embedded data

*=======================================================================
* Call CLUELA on the projected data.
*=======================================================================
       K = KDEF
       SOSDO = .TRUE.
       CALL CLUELA (X, NX, PLIM, M, N, W, C, NC, KMIN, K, KMAX, Z,
     $  POP, F, ROBUST, SOSDO, IWORK, RWORK, U, IFAULT)
       PRINT *, 'cluela returns #', IFAULT, '  Residual: ', U

       FNAME = 'horse_clusters_weighted.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=FERR)
       IF (FERR .NE. 0) PRINT *, 'horse: trouble opening ', FNAME
       CALL WCSVSP (V, N, MVIS, Z, IOU, IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'wcsvsp: error #', IFAULT
       CLOSE (UNIT=IOU)

*         write F values
       FNAME = 'horse_fvals_weighted.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=FERR)
       IF (FERR .NE. 0) PRINT *, 'horse: trouble opening ', FNAME
       DO L = KMIN, KMAX
         WRITE (UNIT=IOU, FMT=40, IOSTAT=FERR) L, F(L)
         IF (FERR .NE. 0) PRINT *, 'horse: trouble writing F stats'
       END DO
       CLOSE (UNIT=IOU)

*=======================================================================
* Call alternate interface
*=======================================================================
       CALL DTIME (TARRAY, RESULT)
       SOSDO = .FALSE.
       UBEST = BIG
       KBEST = 0
       DO H = 1, 10
         K = KDEF
         CALL CLUSTE (X, NX, PLIM, M, N, W, C, NC, KMIN, K, KMAX,
     $  Z, POP, ROBUST, SOSDO, IWORK, RWORK, U, IFAULT)
         IF (IFAULT .EQ. 0) THEN
           IF (K .EQ. KDEF) THEN
             KBEST = K
             IF (U < UBEST) THEN
               UBEST = U
               DO I = 1, N
                 ZBEST(I) = Z(I)
               END DO
             END IF
           END IF
         ELSE
           PRINT *, 'cluste: error #', IFAULT
         END IF
       END DO
       CALL DTIME (TARRAY, RESULT)
       IF (KBEST .NE. KDEF)
     $    PRINT *, 'cluste did not find requested number of clusters!'
       PRINT *, 'cluste returns in ', RESULT, ' seconds.  Residual: ',
     $            UBEST

       FNAME = 'horse_clusters_cluste.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=FERR)
       IF (FERR .NE. 0) PRINT *, 'horse: trouble opening ', FNAME
       CALL WCSVSP (V, N, MVIS, ZBEST, IOU, IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'wcsvsp: error #', IFAULT
       CLOSE (UNIT=IOU)

*=======================================================================
* Call KMNS on projected data.  Remember to initalize centers.
* Matrices in KMNS are the transpose of those in CLUELA!
*=======================================================================
       CALL DTIME (TARRAY, RESULT)
       K = KDEF
       UBEST = BIG
       DO H = 1, 10
         CALL DINIT (X, NX, PLIM, M, N, W, C, NC, K, ROBUST,
     $                RWORK, RWORK(N+1), THRESH)
         DO L = 1, K
           DO J = 1, M            ! transpose center
             CT(L,J) = C(J,L)
           END DO
         END DO                   ! next center to initialize
         CALL KMNS (V, N, M, CT, K, Z, IPERM, POP, AN1, AN2, NCP,
     $               RWORK, ITRAN, LIVE, ITER, WSS, IFAULT)

         U = 0.
         DO L = 1, K              ! combine per-cluster residuals
           U = U + WSS(L)
         END DO
         IF (IFAULT .EQ. 0) THEN
           IF (U < UBEST) THEN
             UBEST = U
             DO I = 1, N
               ZBEST(I) = Z(I)
             END DO
           END IF
         ELSE
           PRINT *, 'kmns: error #', IFAULT
         END IF
       END DO
       CALL DTIME (TARRAY, RESULT)
       PRINT *, 'kmns returns in ', RESULT, 'seconds. Residual: ', UBEST

       FNAME = 'horse_clusters_kmns.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=FERR)
       IF (FERR .NE. 0) PRINT *, 'horse: trouble opening ', FNAME
       CALL WCSVSP (V, N, MVIS, ZBEST, IOU, IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'wcsvsp: error #', IFAULT
       CLOSE (UNIT=IOU)

*=======================================================================
* Call the neural gas routine.
*=======================================================================
       CALL DTIME (TARRAY, RESULT)
       K = KDEF
       CALL SHUFFL (IWORK, N)             ! choose random centers
       DO L = 1, K
         I = IWORK(L)
         DO J = 1, M
           C(J,L) = X(J,I)
         END DO
       END DO
       CALL NEUGAS (X, PLIM, M, N, C, K, Z, IWORK, RWORK, U)
       CALL DTIME (TARRAY, RESULT)
       PRINT *, 'neugas returns in ', RESULT, ' seconds.  Residual: ', U

       FNAME = 'horse_clusters_ng.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=FERR)
       IF (FERR .NE. 0) PRINT *, 'horse: trouble opening ', FNAME
       CALL WCSVSP (V, N, MVIS, Z, IOU, IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'wcsvsp: error #', IFAULT
  90   FORMAT (G9.3, ',,,', G9.3)
       WRITE (UNIT=IOU, FMT=90) C(1,1), C(2,1)
       WRITE (UNIT=IOU, FMT=90) C(1,2), C(2,2)
       CLOSE (UNIT=IOU)

       STOP 'program complete'
      END  ! of horse
