*=======================================================================
* test clustering on the Soybean data
* by Andy Allinger, 2017, public domain
*=======================================================================
      PROGRAM SOY
       IMPLICIT NONE

*           fixed parameters
       INTEGER CATMAX, G, IOU, M, NLIM, O, PLIM
       REAL XBOGUS
       PARAMETER (CATMAX = 5,      ! Most categories any attribute
     $            G = 19,          ! # correct classes
     $            IOU = 1,         ! I/O unit number
     $            XBOGUS = -9.,    ! filler value for missing data
     $            M = 2,           ! projection dimension for visualization
     $            NLIM = 683,      ! Number of objects
     $            O = 35,          ! # input attributes
     $            PLIM = 54)       ! # attibutes after encoding
       INTEGER NN
       PARAMETER (NN = NLIM * NLIM)

*          adjustable parameters
       INTEGER KMIN, KDEF, KMAX
       LOGICAL ROBUST, SOSDO
       PARAMETER (KMIN = 1,             ! fewest clusters to make
     $            KDEF = G,             ! default number of clusters
     $            KMAX = 50,            ! most clusters to make
     $            ROBUST = .TRUE.,      ! use L1 metric?
     $            SOSDO = .FALSE.)      ! minimize sum-of-squares dissimilarity

*              variables in test program
       INTEGER
     $         CATC(O),            ! categories each input attribute
     $         CATD(O),            ! interpret some attributes as dichotomous
     $         FERR,               ! file i/o error
     $         I, J, JBEG, L,      ! counters
     $         IPERM(NLIM),        ! index array
     $         JPERM(NLIM),        ! index array
     $         LL,                 ! line length
     $         NX0(NLIM),          ! work array
     $         Y(NLIM)             ! known classes

       REAL
     $      ACCARI, ACCVI,         ! accuracy vs. external standard
     $      S(CATMAX-1,CATMAX),    ! regular simplex coordinates
     $      XLEA, XGRE             ! extreme values of column of X

       LOGICAL
     $         FEXIST,             ! file exists
     $         ISDIX(O),           ! may be treated as dichotomous
     $         ORDREQ(O)           ! ordinal scaling required
       CHARACTER*32 FNAME, RAW, CLASS(NLIM)
       CHARACTER*256 TEXT

*        external function
       INTEGER LTRIM

*        arguments to clustering subroutine
       INTEGER P, K, N, Z(NLIM), POP(KMAX), IWORK(4*NLIM+25*KMAX),
     $          IFAULT, NX(PLIM,NLIM), NC(PLIM,KMAX)
       REAL X(PLIM,NLIM), W(PLIM), C(PLIM,KMAX), RWORK(2*NLIM+9*KMAX),
     $       F(KMIN:KMAX), U

*        arrays for embedded data
       INTEGER CATR(PLIM), NV(NLIM,NLIM), NCV(NLIM,KMAX)
       REAL CV(NLIM,KMAX), WV(NLIM)

*        workspace for validation
       REAL VALTAB(G*KMAX+G+KMAX), DIF, WORST

*        workspace for multidimensional scaling
       INTEGER MC, MD, MR
       REAL E(NLIM), T(NLIM)
       REAL VC(NLIM,NLIM), VD(NLIM,NLIM), VR(NLIM,NLIM), VDT(NLIM,NLIM)

       DATA FNAME / 'soybean.both' /
*=======================================================================
* Large Soybean Database
* From the University of California at Irvine Machine Learning repository
*
* R.S. Michalski and R.L. Chilausky
* "Learning by Being Told and Learning from Examples:
* An Experimental Comparison of the Two Methods of Knowledge Acquisition
* in the Context of Developing an Expert System for Soybean Disease Diagnosis",
* International Journal of Policy Analysis and Information Systems,
* Vol. 4, No. 2, 1980.
*
* Classes are:
*
* diaporthe-stem-canker, charcoal-rot, rhizoctonia-root-rot, phytophthora-rot,
* brown-stem-rot, powdery-mildew, downy-mildew, brown-spot, bacterial-blight,
* bacterial-pustule, purple-seed-stain, anthracnose, phyllosticta-leaf-spot,
* alternarialeaf-spot, frog-eye-leaf-spot, diaporthe-pod-&-stem-blight,
* cyst-nematode, 2-4-d-injury, herbicide-injury
*
* Attributes are:
*
*___Attribute_________Type__________#___Values_______________________
*   date              Continuous        april,may,june,july,august,
*                                         september,october,?
*   plant-stand       Categorical   2   normal,lt-normal,?
*   precip            Ordinal           lt-norm,norm,gt-norm,?
*   temp              Ordinal           lt-norm,norm,gt-norm,?
*   hail              Categorical   2   yes,no,?
*                     (Dichotomous)
*   crop-hist         Ordinal           diff-lst-year,same-lst-yr,
*                                         same-lst-two-yrs,same-lst-sev-yrs,?
*   area-damaged      Categorical   4   scattered,low-areas,upper-areas,
*                                         whole-field,?
*   severity          Ordinal           minor,pot-severe,severe,?
*   seed-tmt          Categorical   3   none,fungicide,other,?
*   germination       Ordinal           90-100%,80-89%,lt-80%,?
*   plant-growth      Categorical   2   norm,abnorm,?
*   leaves            Categorical   2   norm,abnorm.
*   leafspots-halo    Categorical   3   absent,yellow-halos,no-yellow-halos,?
*   leafspots-marg    Categorical   3   w-s-marg,no-w-s-marg,dna,?
*   leafspot-size     Categorical   3   lt-1/8,gt-1/8,dna,?
*   leaf-shread       Categorical   2   absent,present,?
*                     (Dichotomous)
*   leaf-malf         Categorical   2   absent,present,?
*                     (Dichotomous)
*   leaf-mild         Categorical   3   absent,upper-surf,lower-surf,?
*   stem              Categorical   2   norm,abnorm,?
*   lodging           Categorical   2   yes,no,?
*                     (Dichotomous)
*   stem-cankers      Categorical   4   absent,below-soil,above-soil,
*                                         above-sec-nde,?
*   canker-lesion     Categorical   4   dna,brown,dk-brown-blk,tan,?
*   fruiting-bodies   Categorical   2   absent,present,?
*                     (Dichotomous)
*   external decay    Categorical   3   absent,firm-and-dry,watery,?
*   mycelium          Categorical   2   absent,present,?
*                     (Dichotomous)
*   int-discolor      Categorical   3   none,brown,black,?
*   sclerotia         Categorical   2   absent,present,?
*                     (Dichotomous)
*   fruit-pods        Categorical   4   norm,diseased,few-present,dna,?
*   fruit spots       Categorical   5   absent,colored,brown-w/blk-specks,
*                                         distort,dna,?
*   seed              Categorical   2   norm,abnorm,?
*                     (Dichotomous)
*   mold-growth       Categorical   2   absent,present,?
*                     (Dichotomous)
*   seed-discolor     Categorical   2   absent,present,?
*                     (Dichotomous)
*   seed-size         Categorical   2   norm,lt-norm,?
*                     (Dichotomous)
*   shriveling        Categorical   2   absent,present,?
*                     (Dichotomous)
*   roots             Categorical   3   norm,rotted,galls-cysts,?
*=======================================================================

*        info above is summarized thus:
       DATA CATC / 0, 2, 0, 0, 2, 0, 4, 0, 3, 0,           !  1-10
     $             2, 2, 3, 3, 3, 2, 2, 3, 2, 2,           ! 11-20
     $             4, 4, 2, 3, 2, 3, 2, 4, 5, 2,           ! 21-30
     $             2, 2, 2, 2, 3 /                         ! 31-35

       DATA CATD / 0, 2, 0, 0, 1, 0, 4, 0, 3, 0,           !  1-10
     $             2, 2, 3, 3, 3, 1, 1, 3, 2, 1,           ! 11-20
     $             4, 4, 1, 3, 1, 3, 1, 4, 5, 1,           ! 21-30
     $             1, 1, 1, 1, 3 /                         ! 31-35

       DATA CATR / PLIM * 0 /

       DATA ORDREQ / .FALSE., .FALSE., .TRUE. , .TRUE. , .FALSE.,
     $               .TRUE. , .FALSE., .TRUE. , .FALSE., .TRUE. ,
     $               25 * .FALSE. /

       DATA ISDIX / .FALSE., .FALSE., .FALSE., .FALSE., .TRUE. ,   !  1 -  5
     $              .FALSE., .FALSE., .FALSE., .FALSE., .FALSE.,   !  6 - 10
     $              .FALSE., .FALSE., .FALSE., .FALSE., .FALSE.,   ! 11 - 15
     $              .TRUE. , .TRUE. , .FALSE., .FALSE., .TRUE. ,   ! 16 - 20
     $              .FALSE., .FALSE., .TRUE. , .FALSE., .TRUE. ,   ! 21 - 25
     $              .FALSE., .TRUE. , .FALSE., .FALSE., .TRUE. ,   ! 26 - 30
     $              .TRUE. , .TRUE. , .TRUE. , .TRUE. , .FALSE. /  ! 31 - 35

*           weights (not to be taken seriously)
       DATA W / 2.742, 1.944, 6.565, 2.532, 2.200,
     $          2.227, 1.386, 2.997, 1.378, 2.186,
     $          3.578, 4.254, 3.271, 4.313, 4.405,
     $          2.795, 5.595, 12.57, 3.092, 3.847,
     $          2.163, 2.124, 3.362, 2.623, 12.17,
     $          12.57, 12.57, 3.162, 2.544, 4.050,
     $          4.626, 4.764, 4.301, 5.215, 3.464,
     $            19 * 0. /

       DATA WV / NLIM * 1. /
       DATA NV / NN * 1 /         ! no missing data after embedding

*=======================================================================
*          begin
*=======================================================================
       CALL RNDINI ()             ! initialize random number generator
       N = NLIM
       P = PLIM
       INQUIRE (FILE=FNAME, EXIST=FEXIST, IOSTAT=FERR)  ! get input file
       IF (.NOT. FEXIST) STOP 'File not found'
       IF (FERR .NE. 0) STOP 'Error seeking input file'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='OLD',
     &        ACCESS='SEQUENTIAL', IOSTAT=FERR)
       IF (FERR .NE. 0) PRINT *, 'soy: trouble opening input file'

*=======================================================================
*       read in data
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

           IF (L .EQ. 1) THEN                           ! answer key
             READ (UNIT=RAW, FMT=*, IOSTAT=FERR) CLASS(I)
             IF (FERR .NE. 0) THEN
               PRINT *, 'error reading class label'
               PRINT *, 'raw was ', RAW
             END IF
           ELSE IF (INDEX(RAW, '?') .EQ. 0) THEN          ! non-missing
             READ (UNIT=RAW, FMT=*, IOSTAT=FERR) X(L-1,I)
             IF (FERR .NE. 0) PRINT*,'error reading X, i: ',I,' L:',L
             X(L-1,I) = X(L-1,I) + 1.
             NX(L-1,I) = 1
           ELSE
             X(L-1,I) = XBOGUS
             NX(L-1,I) = 0
           END IF
           JBEG = J+1
           IF (L > O+1) PRINT *, '!trouble reading L>o '
  20     CONTINUE  !  next J
       END DO  ! next I
  30   CLOSE (UNIT=IOU)

*=======================================================================
*          preprocess input
*=======================================================================
*             change date attribute to unit range
       DO I = 1, N
         X(1,I) = X(1,I) / 6.
       END DO

*            rescale ordinal variables to replicate DAISY
       DO L = 1, O
         IF (ORDREQ(L)) THEN
           XLEA = 1.E33
           XGRE = 0.
           DO I = 1, N
             IF (NX(L,I) > 0) THEN
               XLEA = MIN(XLEA, X(L,I))
               XGRE = MAX(XGRE, X(L,I))
             END IF
           END DO
           DO I = 1, N
             IF (NX(L,I) > 0) X(L,I) = (X(L,I) - XLEA) / (XGRE - XLEA)
           END DO
         END IF
       END DO

*             embed
       CALL MDS (X,NX,PLIM,O,N,W,CATC,ROBUST,VC,NLIM,MC,E,T,IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'mds on categorical #', IFAULT
       IF (MC < M) PRINT *, 'mds on categorical: ', MC, ' directions!'

*             embed using dictomous variables
       DO I = 1, N
         DO L = 1, O
           IF (ISDIX(L)) X(L,I) = X(L,I) - 1.  ! must be {0,1} instead of {1,2}
         END DO
       END DO
       CALL MDS (X,NX,PLIM,O,N,W,CATD,ROBUST,VD,NLIM,MD,E,T,IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'mds with dichotomous #', IFAULT
       IF (MD < M) PRINT *, 'mds with dichotomous: ', MD, ' directions!'
       DO I = 1, N
         DO L = 1, O
           IF (ISDIX(L)) X(L,I) = X(L,I) + 1.
         END DO
       END DO

*             recode categoricals with corresponding column of simplex
       CALL RSCODE (X,NX,PLIM,P,N,W,CATC,O,S,CATMAX,ROBUST,IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'rscode: error #', IFAULT

*=======================================================================
*             call the subroutine
*=======================================================================
       K = KDEF
       CALL CLUELA (X, NX, PLIM, P, N, W, C, NC, KMIN, K, KMAX, Z, POP,
     $               F, ROBUST, SOSDO, IWORK, RWORK, U, IFAULT)
       PRINT *, 'original problem: cluela returns #', IFAULT,
     $              ' made ', K, ' clusters'

*            relabel class ID's
       DO I = 1, N
         NX0(I) = 1
       END DO
       CALL IREID (CLASS, Y, NX0, IPERM, JPERM, N, RAW, IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'ireid: error #', IFAULT

*            find accuracy
       IF (IFAULT.EQ. 0) THEN
         CALL ARIND (Y, Z, N, G, K, VALTAB, ACCARI)  ! adjusted Rand index
         CALL VINFO (Y, Z, N, G, K, VALTAB, ACCVI)   ! variation of information
       END IF
       PRINT *, 'Returned residual: ', U
       PRINT *, 'Adjusted Rand Index: ', ACCARI
       PRINT *, 'Variation of Information:', ACCVI

*=======================================================================
* Write CSV file with coordinates in different columns by cluster
*=======================================================================
       FNAME = 'soy_clusters.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'soy: trouble opening ', FNAME
       CALL WCSVSP (VC, N, M, Z, IOU, IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'wcsvsp: error #', IFAULT
       CLOSE (UNIT=IOU)

*=======================================================================
* Write CSV file with coordinates in different columns by class
*=======================================================================
       FNAME = 'soy_classes.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'soy: trouble opening ', FNAME
       CALL WCSVSP (VC, N, M, Y, IOU, IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'wcsvsp: error #', IFAULT
       CLOSE (UNIT=IOU)

*=======================================================================
*         write F values
*=======================================================================
       FNAME = 'soy_fvals.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'soy: trouble opening ', FNAME
  80   FORMAT (I5, ',', G9.3)
       DO L = KMIN, KMAX
         WRITE (UNIT=IOU, FMT=80, IOSTAT=IFAULT) L, F(L)
         IF (IFAULT .NE. 0) PRINT *, 'soy: trouble writing F stats'
       END DO
       CLOSE (UNIT=IOU)

*=======================================================================
* Repeat clustering interpreting 12 variables as dichotomous.
*=======================================================================
       K = KDEF
       DO J = 1, N
         DO I = 1, N
           VDT(I,J) = VD(J,I)
         END DO
       END DO
       CALL CLUELA (VDT, NV, NLIM, MD, N, WV, CV, NCV, KMIN, K, KMAX,
     $               Z, POP, F, .FALSE., SOSDO, IWORK, RWORK, U, IFAULT)
       PRINT *, 'dichotomous problem: cluela returns #', IFAULT,
     $          ' made ', K, ' clusters'

       IF (IFAULT.EQ. 0) THEN
         CALL ARIND (Y, Z, N, G, K, VALTAB, ACCARI)
         CALL VINFO (Y, Z, N, G, K, VALTAB, ACCVI)
       END IF
       PRINT *, 'Returned residual: ', U
       PRINT *, 'Adjusted Rand Index: ', ACCARI
       PRINT *, 'Variation of Information:', ACCVI

       FNAME = 'soy_dichot_clusters.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'soy: trouble opening ', FNAME
       CALL WCSVSP (VD, N, M, Z, IOU, IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'wcsvsp: error #', IFAULT
       CLOSE (UNIT=IOU)

       FNAME = 'soy_dichot_classes.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'soy: trouble opening ', FNAME
       CALL WCSVSP (VD, N, M, Y, IOU, IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'wcsvsp: error #', IFAULT
       CLOSE (UNIT=IOU)

       FNAME = 'soy_dichot_fvals.csv'
       OPEN (UNIT=IOU, FILE=FNAME, STATUS='NEW', IOSTAT=IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'soy: trouble opening ', FNAME
       DO L = KMIN, KMAX
         WRITE (UNIT=IOU, FMT=80, IOSTAT=IFAULT) L, F(L)
         IF (IFAULT .NE. 0) PRINT *, 'soy: trouble writing F stats'
       END DO
       CLOSE (UNIT=IOU)

*             embed the recoded data
       CALL MDS (X,NX,PLIM,P,N,W,CATR,ROBUST,VR,NLIM,MR,E,T,IFAULT)
       IF (IFAULT .NE. 0) PRINT *, 'mds on recoded #', IFAULT
       IF (MR < M) PRINT *, 'mds on recoded: ', MR, ' directions!'

       WORST = 0.
       DO J = 1, N
         DO I = 1, N
           DIF = MIN(ABS(VC(I,J) + VR(I,J)), ABS(VC(I,J) - VR(I,J)))
           WORST = MAX(WORST, DIF)
         END DO
       END DO
       PRINT *, 'greatest difference between original ',
     $   'and recoded embeddings: ', WORST
       STOP 'program complete'
      END  ! of t_soy.f
