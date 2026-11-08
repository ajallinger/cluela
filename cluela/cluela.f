*-----------------------------------------------------------------------
* cluela.f - Data clustering by the k-means method
* By Andy Allinger, 2011-2017, released to the public domain.
* This program may be used by any person for any purpose.
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        BIT_SIZE replacement for old compilers
*-----------------------------------------------------------------------
      FUNCTION BITSZ ()
       IMPLICIT NONE
       INTEGER BITSZ
       INTEGER I, J

       J = 0
       J = NOT(J)
       BITSZ = 0
       DO I = 1, 1024
         J = ISHFT(J, -1)
         IF (.NOT. BTEST(J, 0)) THEN
           BITSZ = I
           RETURN
         END IF
       END DO
       RETURN      ! error!
      END ! function BITSZ


*-----------------------------------------------------------------------
*         almost equality test
*-----------------------------------------------------------------------
      FUNCTION ALMEQ (R, S)
       IMPLICIT NONE
       REAL R, S
       LOGICAL ALMEQ

       REAL AR, AS, TWOEPS, LIT
       LOGICAL FIRST
       REAL RMACHINE                 ! external function
       SAVE TWOEPS, LIT, FIRST
       DATA FIRST / .TRUE. /

       IF (FIRST) THEN
         TWOEPS = 2. * RMACHINE(1)
         LIT = RMACHINE(2)
         FIRST = .FALSE.
       END IF
       AR = ABS(R)
       AS = ABS(S)
       IF (AR < LIT .AND. AS < LIT) THEN      ! all denormal numbers are equal
         ALMEQ = .TRUE.
       ELSE
         ALMEQ = ABS(R - S) < TWOEPS * MAX(AR, AS)
       END IF
       RETURN
      END ! of almeq


*-----------------------------------------------------------------------
* Random permutation of 1...N
* origin:  Durstenfeld, 1964
*-----------------------------------------------------------------------
      SUBROUTINE SHUFFL (I, N)
       IMPLICIT NONE
       INTEGER I, N
       DIMENSION I(N)
       INTEGER J, K, L          ! local variables
       INTEGER HAT              ! external function

       DO J = 1, N              ! begin
         I(J) = J
       END DO
       DO J = N, 2, -1
         K = HAT(J)
         L = I(K)
         I(K) = I(J)
         I(J) = L
       END DO
       RETURN
      END ! of shuffl


*-----------------------------------------------------------------------
* Bisection search for the first element greater than a value
*
*___Name_______Type_____In/Out_____Description__________________________
*   X(N)       REAL     In         Array (in increasing order)
*   N          INTEGER  In         Length of X
*   V          REAL     In         The value sought
*   I          INTEGER  Out        Index where X > V found
*-----------------------------------------------------------------------
      SUBROUTINE BSGTR (X, N, V, I)
       IMPLICIT NONE
       INTEGER N, I
       REAL X, V
       DIMENSION X(N)

       INTEGER IL, IR

*       check input
       IF (N < 1) RETURN

*        V is greater than X(N) or less than X(1)
       IF (V < X(1)) THEN
         I = 1
         RETURN
       ELSE IF (V > X(N)) THEN
         I = N
         RETURN
       END IF

*        bisection search
       IL = 1
       IR = N

  10   I = (IL + IR) / 2
       IF (I .EQ. IL) THEN
         IF (X(I) > V) THEN
           I = IL
         ELSE
           I = IR
         END IF
         RETURN
       END IF

       IF (X(I) .LE. V) THEN
         IL = I
       ELSE
         IR = I
       END IF
       GO TO 10
      END  ! of bsgtr


*-----------------------------------------------------------------------
* Passive sort of unsigned integers by an unstable radix algorithm
*
*___Name_____Type______In/Out___Description_____________________________
*   IX(N)    Integer   In       Input array
*   IND(N)   Integer   Out      Permutation vector
*   N        Integer   In       Length of array
*   M        Integer   In       Largest value in array
*-----------------------------------------------------------------------
      SUBROUTINE USORTU (IX, IND, N, M)
       IMPLICIT NONE
       INTEGER IX, IND, N, M
       DIMENSION IX(N), IND(N)

       INTEGER B, I, J, K, L, IL(0:31), IG(0:31), LG(0:31)  ! local variables
       INTEGER BITSZ                                        ! external function

*         Begin
       DO I = 1, N                ! initialize
         IND(I) = I
       END DO
       IF (N < 2) RETURN          ! validate

       J = BITSZ() - 1             ! greatest bit position
       B = J                       ! avoid compiler warning
       DO I = J, 0, -1
         B = I
         IF (BTEST(M, I)) GO TO 10
       END DO

  10   L = B                 ! most significant bit
       IL(L) = 1
       IG(L) = N
       LG(L) = -1

*            Sort
  20   I = IL(L)
       J = IG(L)
       IF (I .GE. J) GO TO 50    ! if empty or single, back up

*           Find a set bit at the front of the partition
  30   IF (I .LE. J) THEN
         IF (.NOT. BTEST(IX(IND(I)), L)) THEN
           I = I + 1
           GO TO 30
         END IF
       END IF

*           Find a clear bit at the rear of the partition
  40   IF (I .LE. J) THEN
         IF (BTEST(IX(IND(J)), L)) THEN
           J = J - 1
           GO TO 40
         END IF
       END IF

*             Swap these elements
       IF (I < J) THEN
         K = IND(I)
         IND(I) = IND(J)
         IND(J) = K
         I = I + 1
         J = J - 1
         GO TO 30
       END IF

*              Traverse tree
       IF (L .EQ. 0) GO TO 50  ! At the least significant bit, go back up
       IL(L-1) = IL(L)
       IG(L-1) = J
       LG(L) = -1
       L = L - 1                     ! Descend
       GO TO 20

  50   L = L + 1                     ! Back up
       IF (L > B) RETURN
       IF (LG(L) > 0) GO TO 50  ! Both subpartitions sorted; go back up
       IL(L-1) = IG(L-1) + 1         ! Do the greater sub-region
       IG(L-1) = IG(L)
       LG(L) = +1                    ! Mark as greater
       L = L - 1                     ! Descend
       GO TO 20
      END   ! of usortu


*-----------------------------------------------------------------------
* Dissimilarity function with missing data support and
* an optional correction to minimize within-cluster square distance.
*
*___Name_____Type______In/Out___Description________________________
*   X(P)     Real      In       Data vector of object
*   NX(P)    Integer   In       Data exists this observation
*   C(P)     Real      In       Data vector of center
*   NC(P)    Integer   In       Count of data present this center
*   P        Integer   In       Dimension of the data
*   W(P)     Real      In       Importance weights
*   ROBUST   Logical   In       Use robust statistics
*   ISIN     Integer   In       -1 delete point, 0 neutral, +1 insert point
*   DSM      Real      Out      Distance to point
*-----------------------------------------------------------------------
      FUNCTION DSM (X, NX, C, NC, P, W, ROBUST, ISIN)
       IMPLICIT NONE
       INTEGER NX, NC, P, ISIN
       REAL X, C, W, DSM
       LOGICAL ROBUST
       DIMENSION X(P), NX(P), C(P), NC(P), W(P)

*            local variables
       INTEGER J, NCJ
       REAL BIG, DIF, FACT, TOT, WSUM
       LOGICAL SAFDIV            ! external function
       PARAMETER (BIG = 1.E30)

*          begin
       TOT = 0.
       WSUM = 0.
       DO J = P, 1, -1
         IF (NX(J) > 0 .AND. NC(J) > 0) THEN
           DIF = X(J) - C(J)
           IF (ROBUST) THEN
             TOT = TOT + W(J) * ABS(DIF)          ! Manhattan distance
           ELSE
             NCJ = NC(J)
             IF (NCJ + ISIN .EQ. 0) THEN
               FACT = 1.
             ELSE
               FACT = FLOAT(NCJ) / FLOAT(NCJ + ISIN)
             END IF
             TOT = TOT + W(J) * FACT * DIF **2    ! Pythagorean distance
           END IF
           WSUM = WSUM + W(J)
         END IF  ! data available ?
       END DO  ! next J

*           adjust for missing values
       IF (SAFDIV(TOT, WSUM)) THEN
         DSM = TOT / WSUM
       ELSE
         DSM = BIG
       END IF
       RETURN
      END ! of dsm


*-----------------------------------------------------------------------
* Count cluster populations
*
*___Name__________Type______In/Out_____Description______________________
*   Z(N)          Integer   In         Membership of each point
*   N             Integer   In         Number of points
*   POP(K)        Integer   Out        Populations each cluster
*   K             Integer   In         Number of clusters
*-----------------------------------------------------------------------
      SUBROUTINE CNTPOP (Z, N, POP, K)
       IMPLICIT NONE
       INTEGER Z, N, POP, K
       DIMENSION Z(N), POP(K)
       INTEGER I, L

*           Begin
       DO L = 1, K         ! clear
         POP(L) = 0
       END DO
       DO I = 1, N         ! add
         L = Z(I)
         POP(L) = POP(L) + 1
       END DO
       RETURN
      END  ! of cntpop


*-----------------------------------------------------------------------
* Let each center be the average of its points
*
*___Name_________Type______In/Out____Description_________________________
*   X(PLIM,N)    Real      In        Data points
*   NX(PLIM,N)   Integer   In        Data exists
*   PLIM         Integer   In        Row dimension of X
*   P            Integer   In        Variables each datum
*   N            Integer   In        # of data points
*   C(PLIM,K)    Real      Out       Centers
*   NC(PLIM,K)   Integer   Out       Center data available
*   K            Integer   In        Number of clusters
*   Z(N)         Integer   In        Cluster membership, each datum
*   POP(K)       Integer   In        Population, each cluster
*   ROBUST       Logical   In        Median instead of mean
*   LIVEC(K)     Integer   In        .GE. 1 => Need to recompute this center
*   IPERM(N)     Integer   Neither   Workspace
*   WORK(N)      Real      Neither   Workspace
*-----------------------------------------------------------------------
      SUBROUTINE CENT (X, NX, PLIM, P, N, C, NC, K, Z,
     $  POP, ROBUST, LIVEC, IPERM, WORK)

       IMPLICIT NONE
       INTEGER NX, PLIM, P, N, NC, K, Z, POP, LIVEC, IPERM
       REAL X, C, WORK
       LOGICAL ROBUST
       DIMENSION X(PLIM,N), NX(PLIM,N), C(PLIM,K), NC(PLIM,K),
     $  Z(N), POP(K), LIVEC(K), IPERM(N), WORK(N)

       INTEGER H, I, I1, I5, I9, J, L

       IF (ROBUST) THEN
         CALL USORTU (Z, IPERM, N, K)            ! sort by cluster ID
         DO J = 1, P
           I9 = 0                 ! index of last object in cluster
           DO 50 L = 1, K
             I1 = I9 + 1          ! index of first object in cluster
             I9 = I9 + POP(L)
             IF (LIVEC(L) .LE. 0) GO TO 50       ! skip dead
             H = 0
             DO I5 = I1, I9
               I = IPERM(I5)
               IF (NX(J,I) > 0) THEN             ! fill work array
                 H = H + 1
                 WORK(H) = X(J,I)
               END IF
             END DO    ! next point
             NC(J,L) = H
             IF (H > 0) CALL MEDIAN (WORK, H, C(J,L))
  50       CONTINUE  ! next cluster
         END DO  ! next variable

       ELSE
         DO L = 1, K
           DO J = 1, P
             C(J,L) = 0.
             NC(J,L) = 0
           END DO
         END DO
         DO I = 1, N
           L = Z(I)
           DO J = 1, P
             IF (NX(J,I) > 0) THEN
               C(J,L) = C(J,L) + X(J,I)
               NC(J,L) = NC(J,L) + 1
             END IF
           END DO
         END DO
         DO L = 1, K
           DO J = 1, P
             IF (NC(J,L) > 0) C(J,L) = C(J,L) / FLOAT(NC(J,L))
           END DO
         END DO
       END IF

       RETURN
      END  ! of cent


*-----------------------------------------------------------------------
* Initialize cluster centers based on distance
*
* Refer to:
*    "kmeans++: the advantages of careful seeding"
*    David Arthur and Sergei Vassilvitskii
*    Proceedings of the eighteenth annual ACM-SIAM Symposium
*      on Discrete Algorithms, 2007
*
*___Name_________Type______I/O_______Description________________________
*   X(PLIM,N)    Real      In        Data points
*   NX(PLIM,N)   Integer   In        Data exists this observation
*   PLIM         Integer   In        Row dimension of X
*   P            Integer   In        Dimension of the data
*   N            Integer   In        Number of points
*   W(P)         Real      In        Importance weights
*   C(PLIM,K)    Real      Out       Center points of clusters
*   NC(PLIM,K)   Integer   Out       Count of center data available
*   K            Integer   In        Number of clusters to make
*   ROBUST       Logical   In        Use medians instead of means
*   SHORTX(N)    Real      Neither   Squared distance to each point
*   CUMD(N)      Real      Neither   Cumulative squared distance
*   THRESH       Real      Out       Threshhold distance for dp-means
*-----------------------------------------------------------------------
      SUBROUTINE DINIT (X, NX, PLIM, P, N, W, C, NC, K, ROBUST,
     $                    SHORTX, CUMD, THRESH)

       IMPLICIT NONE
       INTEGER NX, PLIM, P, N, NC, K
       REAL X, W, C, SHORTX, CUMD, THRESH
       LOGICAL ROBUST
       DIMENSION X(PLIM,N), NX(PLIM,N), W(P), C(PLIM,K), NC(PLIM,K),
     $  SHORTX(N), CUMD(N)

*          Local variables
       INTEGER I, J, L
       REAL BIG, D, PAD, TOT, U

*          External functions
       INTEGER HAT
       REAL DSM, URAND

       DATA BIG / 1.E30 /

* Begin.  Choose each center with probability proportional to its
* least squared distance from any existing center.
       DO I = 1, N
         SHORTX(I) = BIG
       END DO

       I = HAT(N)  ! choose first center at random
       DO J = 1, P                         ! copy
         C(J,1) = X(J,I)
         NC(J,1) = NX(J,I)
       END DO

       DO L = 2, K
         PAD = 0.
         DO I = 1, N             ! compare each point to previous center
           D = DSM (X(1,I), NX(1,I), C(1,L-1), NC(1,L-1), P,W,ROBUST, 0)
           IF (D < SHORTX(I)) SHORTX(I) = D
           IF (SHORTX(I) < PAD) PAD = SHORTX(I)
         END DO
         TOT = 0.
         DO I = 1, N
           TOT = TOT + SHORTX(I) - PAD
           CUMD(I) = TOT
         END DO                      ! next data point

         U = URAND() * TOT    ! uniform at random over cumulative distance
         CALL BSGTR (CUMD, N, U, I)
         DO J = 1, P         ! assign center
           C(J,L) = X(J,I)
           NC(J,L) = NX(J,I)
         END DO
       END DO               ! next center to initialize

*   Threshhold distance is longest distance of any point to its center
       THRESH = 0.
       DO I = 1, N
         IF (SHORTX(I) > THRESH) THRESH = SHORTX(I)
       END DO
       RETURN
      END  ! of dinit


*-----------------------------------------------------------------------
* Partition clustering for mixed and missing data
* origin:  Hugo Steinhaus, 1956
*
* Refer to:
*
*    "Revisiting K-means:  New Algorithms via Bayesian Nonparametrics"
*    Brian Kulis and Michael I. Jordan
*    Proceedings of the 29th International Conference on Machine Learning, 2012
*
*___Name_____________Type______I/O_______Description____________________
*   X(PLIM,N)         Real      In        Data points
*   NX(PLIM,N)        Integer   In        Data exists this observation
*   PLIM              Integer   In        Row dimension of X
*   P                 Integer   In        Dimension of the data
*   N                 Integer   In        Number of points
*   W(P)              Real      In        Importance weights
*   C(PLIM,KMAX)      Real      Both      Center points of clusters
*   NC(PLIM,KMAX)     Integer   Both      Count of center data available
*   KMIN              Integer   In        Fewest clusters to make
*   K                 Integer   Both      Default number of clusters to make
*   KMAX              Integer   In        Most clusters to make
*   Z(N)              Integer   Out       What cluster a point is in
*   POP(KMAX)         Integer   Out       # of points in each cluster
*   MISS              Logical   In        Some data is missing
*   ROBUST            Logical   In        Use medians instead of means
*   SOSDO             Logical   In        Sum of squared distance objective
*   THRESH            Real      In        Threshhold for forming new cluster
*   LIVEC(KMAX)       Integer   Neither   .GE. 1 => cluster needs updating
*   LIVEX(N)          Integer   Neither   If each point is in the live set
*   IPERM(N)          Integer   Neither   Permutation vector
*   TREE(NS,2*KMAX)   Integer   Neither   generalized hyperplane tree
*   CRIT(4,2*KMAX)    Real      Neither   critical distances in tree
*   SHORTC(KMAX)      Real      Neither   Shortest distance to another center
*   SHORTX(N)         Real      Neither   Shortest distance to current center
*   WORK(N)           Real      Neither   For averages
*   U                 Real      Out       Residual
*   IFAULT            Integer   Out       Error code 0-5 (see CLUELA)
*-----------------------------------------------------------------------
      SUBROUTINE KCLUST (X, NX, PLIM, P, N, W, C, NC, KMIN, K, KMAX,
     $ Z, POP, MISS, ROBUST, SOSDO, THRESH, LIVEC, LIVEX, IPERM, TREE,
     $ CRIT, SHORTC, SHORTX, WORK, U, IFAULT)

       IMPLICIT NONE
       INTEGER NX, PLIM, P, N, NC, KMIN, K, KMAX, Z, POP, LIVEX, LIVEC,
     $  IPERM, TREE, IFAULT
       REAL X, W, C, THRESH, CRIT, SHORTC, SHORTX, WORK, U
       LOGICAL MISS, ROBUST, SOSDO
       INTEGER NS
       PARAMETER (NS = 12)
       DIMENSION X(PLIM,N), NX(PLIM,N), W(P), C(PLIM,KMAX),
     $  NC(PLIM,KMAX), Z(N), POP(KMAX), LIVEC(KMAX), LIVEX(N),
     $  IPERM(N), TREE(NS,2*KMAX), CRIT(4,2*KMAX), SHORTC(KMAX),
     $  SHORTX(N), WORK(N)

*               constants
       INTEGER HS, ITER
       REAL BIG
       PARAMETER (BIG = 1.E30, ! arbitrary large number
     $            HS = 13,     ! length of history hash buffer
     $            ITER = 1000) ! maximum iterations

*                local variables
       INTEGER
     $         H, I, I1, J, L, M,  ! counters
     $         HL,                 ! length in bytes of Z()
     $         IFAR,               ! farthest outlying point
     $         ISIN,               ! -1 if point is in cluster compared to
     $         K0,                 ! temp. value of K
     $         LOLD, LNEW          ! present, new cluster ID

       INTEGER*4 HASH, HIST(HS)   ! hashes of prior Z()

       REAL BEST,                 ! shortest distance to a center
     $      D,                    ! a dissimilarity
     $      DFAR,                 ! farthest distance of any point
     $      WALL                  ! sum of weights

       LOGICAL CHANGE,      ! whether any points have been reassigned
     $         METRIC       ! use a generalized hyperplane tree

*                       external functions
       INTEGER BITSZ
       REAL DSM
       LOGICAL ALMEQ

*-----------------------------------------------------------------------
*           begin
*-----------------------------------------------------------------------
       IFAULT = 0                          ! validate, initialize
       IF (K < 1 .OR. K > N) THEN
         IFAULT = 3
         RETURN
       END IF
       HL = N * BITSZ() / 8
       METRIC = .NOT. (MISS .OR. SOSDO)
       DO J = 1, P
         IF (W(J) < 0.) METRIC = .FALSE.
       END DO

*          initialize
       DO L = 1, K
         LIVEC(L) = 2
         SHORTC(L) = BIG
       END DO
       DO I = 1, N
         LIVEX(I) = 1
         Z(I) = 1                 ! arbitrary valid ID
         SHORTX(I) = BIG
       END DO

*-----------------------------------------------------------------------
*                      main loop
*-----------------------------------------------------------------------
       DO H = 1, ITER       ! main loop
         CHANGE = .FALSE.
         CALL SHUFFL (IPERM, N)

*-----------------------------------------------------------------------
*          Put centers into generalized hyperplane tree
*-----------------------------------------------------------------------
         IF (METRIC) THEN
           CALL GHTINI (K, TREE, I)
           DO L = 1, K
             CALL GHTINS (PLIM, P, W, C, K, ROBUST, L, TREE, CRIT, I)
           END DO

*-----------------------------------------------------------------------
* If the distance from a point to its center is less than or equal to half
* the distance of the center to its nearest center, no transfer can occur.
*-----------------------------------------------------------------------
           DO L = 1, K          ! nearest center each center
             BEST = BIG
             CALL GHTSEL (C, PLIM, P, K, W, C, K, ROBUST, .TRUE., L,
     $                        LIVEC, 1, LNEW, BEST, TREE, CRIT)
             IF (ROBUST) THEN
               SHORTC(L) = 0.5 * BEST
             ELSE                            ! quarter of squared distance
               SHORTC(L) = 0.25 * BEST
             END IF
           END DO
         END IF

*-----------------------------------------------------------------------
* Find nearest center for each point.  Points are in the live set if the
* distance to the nearest center has increased.
* Clusters updated on this step are set LIVEC = 2 and will be recomputed.
* Clusters updated on the prior step have LIVEC = 1 and are still live.
* If neither the point nor the cluster are in their respective live sets,
* skip the comparison.
*-----------------------------------------------------------------------
         DFAR = -BIG
         IFAR = 1
         DO 40 I1 = 1, N
           I = IPERM(I1)
           LOLD = Z(I)                ! old assignment
           LNEW = LOLD
           BEST = SHORTX(I)

           IF (METRIC) THEN
             IF (BEST .LE. SHORTC(LOLD)) GO TO 40   ! don't change point
             CALL GHTSEL (X, PLIM, P, N, W, C, K, ROBUST, .FALSE., I,
     $                      LIVEC, LIVEX(I), LNEW, BEST, TREE, CRIT)

           ELSE                                    ! exhaustive search
             DO 20 L = 1, K
               IF (LIVEX(I) .LE. 0 .AND. LIVEC(L) .LE. 0) GO TO 20
               IF ( (.NOT. SOSDO) .OR. (H .EQ. 1) ) THEN
                 ISIN = 0
               ELSE IF (L .EQ. LOLD) THEN
                 ISIN = -1
               ELSE
                 ISIN = +1
               END IF
               D = DSM (X(1,I), NX(1,I), C(1,L), NC(1,L),
     $                                 P, W, ROBUST, ISIN)
               IF (D < BEST) THEN
                 BEST = D
                 LNEW = L
               END IF
  20         CONTINUE  ! next L
           END IF  ! use generalized hyperplane tree?
           SHORTX(I) = BEST
           IF (BEST > DFAR) THEN
             DFAR = BEST
             IFAR = I
           END IF

*               check for stray points
           IF (ALMEQ(BEST, BIG)) THEN
             IFAULT = 4
             RETURN
           END IF

*-----------------------------------------------------------------------
*              reassign point
*-----------------------------------------------------------------------
           IF (LOLD .NE. LNEW) THEN
             LIVEC(LOLD) = 2                    ! recenter old cluster
             LIVEC(LNEW) = 2                    ! recenter new cluster
             Z(I) = LNEW
             CHANGE = .TRUE.
             IF (SOSDO .AND. H > 1) THEN
               DO J = 1, P                        ! update centers
                 IF (NX(J,I) > 0) THEN
                   M = MAX(1, NC(J,LOLD) - 1)       ! don't divide by zero
                   C(J,LOLD) = C(J,LOLD) -
     $                    (X(J,I) - C(J,LOLD)) / FLOAT(M)
                   NC(J,LOLD) = NC(J,LOLD) - 1
                   C(J,LNEW) = C(J,LNEW) +
     $                    (X(J,I) - C(J,LNEW)) / FLOAT(NC(J,LNEW) + 1)
                   NC(J,LNEW) = NC(J,LNEW) + 1
                 END IF  ! data present?
               END DO  ! next J
             END IF  ! SOS mode?
           END IF   ! assignment changed?
  40     CONTINUE  ! next i

*-----------------------------------------------------------------------
*        if shortest distance exceeds a threshhold, make new cluster
*-----------------------------------------------------------------------
         IF (.NOT. CHANGE .AND. DFAR > THRESH) THEN
           IF (K .GE. KMAX) THEN
             IFAULT = 5
             RETURN
           END IF
           K = K + 1
           LNEW = K
           LOLD = Z(IFAR)
           LIVEC(LOLD) = 2
           LIVEC(LNEW) = 2
           Z(IFAR) = LNEW
           CHANGE = .TRUE.
         END IF

         IF (.NOT. CHANGE .AND. H > 1) GO TO 80   ! algorithm halts

*-----------------------------------------------------------------------
*               delete empty clusters
*-----------------------------------------------------------------------
         CALL CNTPOP (Z, N, POP, K)      ! count populations
         K0 = K
         DO L = K0, 1, -1
           IF (POP(L) .EQ. 0) THEN
             IF (K .LE. KMIN) THEN
               IFAULT = 1
               RETURN
             END IF

             DO J = 1, P
               C(J,L) = C(J,K)
               NC(J,L) = NC(J,K)
             END DO
             DO I = 1, N
               IF (Z(I) .EQ. K) Z(I) = L
             END DO
             POP(L) = POP(K)
             LIVEC(L) = LIVEC(K)
             SHORTC(L) = SHORTC(K)
             POP(K) = 0
             K = K - 1
           END IF
         END DO

*-----------------------------------------------------------------------
*           find cluster centers
*-----------------------------------------------------------------------
         DO L = 1, K                            ! decrement LIVEC
           LIVEC(L) = LIVEC(L) - 1
         END DO
         CALL CENT (X, NX, PLIM, P, N, C, NC, K, Z, POP, ROBUST,
     $                   LIVEC, IPERM, WORK)

*-----------------------------------------------------------------------
*            check if anything is changing
*-----------------------------------------------------------------------
         HASH = -2128831035           ! canonical arbitrary initialization
         CALL FNV32 (Z, HL, HASH)        ! hash of membership vector
         DO J = 1, MIN(H-1, HS)
           IF (HASH .EQ. HIST(J)) GO TO 80      ! exit on loop criterion
         END DO
         J = MOD(H, HS)
         IF (J .EQ. 0) J = HS
         HIST(J) = HASH                          ! remember

*           Determine live points
         DO I = 1, N         ! no live points
           LIVEX(I) = 0
         END DO
         IF (SOSDO) THEN
           ISIN = -1
         ELSE
           ISIN = 0
         END IF
         DO I = 1, N
           L = Z(I)
           IF (LIVEC(L) > 0) THEN
             D = DSM (X(1,I), NX(1,I), C(1,L), NC(1,L),
     $                  P, W, ROBUST, ISIN)
             IF (D > SHORTX(I)) LIVEX(I) = 1  ! points now farther are live
             SHORTX(I) = D
           END IF
         END DO

       END DO                 ! next iteration
       IFAULT = 2               ! too many iterations
       RETURN

*-----------------------------------------------------------------------
*              compute residual
*-----------------------------------------------------------------------
  80   U = 0.
       DO I = 1, N
         L = Z(I)
         D = DSM (X(1,I), NX(1,I), C(1,L), NC(1,L), P, W, ROBUST, 0)
         U = U + D
       END DO
       WALL = 0.
       DO J = 1, P
         WALL = WALL + W(J)
       END DO
       U = U * WALL
       RETURN
      END ! of kclust


*-----------------------------------------------------------------------
* Validate input, cluster several tries, and keep the best result.
* Choose the number of clusters.
*
*___Variable__________Type______I/O_______Description___________________
*   X(PLIM,N)         Real      In        Data points in R^P
*   NX(PLIM,N)        Integer   In        1 indicates data is present for
*                                            the corresponding X value.
*   PLIM              Integer   In        Row dimension of X
*   P                 Integer   In        Measurements per point
*   N                 Integer   In        Number of points
*   W(P)              Real      In        Importance weights of variables
*                                              if in doubt use /P * 1./
*   C(PLIM,KMAX)      Real      Out       Cluster centroids
*   NC(PLIM,KMAX)     Integer   Out       Count of centroid data available
*   KMIN              Integer   In        Fewest clusters to make
*   K                 Integer   Both      Default/actual clusters to make
*   KMAX              Integer   In        Most clusters to make
*   Z(N)              Integer   Out       What cluster each point is in
*   POP(KMAX)         Integer   Out       Number of points in each cluster
*   F(KMIN:KMAX)      Real      Out       Low value indicates good # clusters
*   ROBUST            Logical   In        Use L1 statistics
*   SOSDO             Logical   In        Sum of square distance objective
*   IPERM(N)          Integer   Neither   Orders points by cluster
*   ZLOF(N)           Integer   Neither   Saved memberships with low F value
*   ZTRY(N)           Integer   Neither   Membership array passed to KCLUST
*   LIVEX(N)          Integer   Neither   Whether each point is live
*   LIVEC(KMAX)       Integer   Neither   Whether each cluster is live
*   TREE(NS,2*KMAX)   Integer   Neither   generalized hyperplane tree
*   SHORTX(N)         Real      Neither   Shortest distance each point
*   CUMD(N)           Real      Neither   Cumulative distribution of distance
*   SHORTC(KMAX)      Real      Neither   Shortest distance, each center
*   CRIT(4,2*KMAX)    Real      Neither   critical distances in tree
*   U                 Real      Out       Residual error
*   IFAULT            Integer   Out       Error code
*-----------------------------------------------------------------------
      SUBROUTINE CLUE77 (X, NX, PLIM, P, N, W, C, NC,
     $  KMIN, K, KMAX, Z, POP, F, ROBUST, SOSDO, IPERM, ZLOF, ZTRY,
     $  LIVEX, LIVEC, TREE, SHORTX, CUMD, SHORTC, CRIT, U, IFAULT)

       IMPLICIT NONE
       INTEGER NX, PLIM, P, N, NC, KMIN, K, KMAX, Z, POP, IPERM,
     $  ZLOF, ZTRY, LIVEX, LIVEC, TREE, IFAULT
       REAL X, W, C, F, CRIT, SHORTX, CUMD, SHORTC, U
       LOGICAL ROBUST, SOSDO
       INTEGER NS
       PARAMETER (NS = 12)
       DIMENSION X(PLIM,N), NX(PLIM,N), W(P), C(PLIM,KMAX),
     $  NC(PLIM,KMAX), Z(N), POP(KMAX), F(KMIN:KMAX), IPERM(N),
     $  ZLOF(N), ZTRY(N), LIVEX(N), LIVEC(KMAX), TREE(NS,2*KMAX),
     $  CRIT(4,2*KMAX), SHORTX(N), CUMD(N), SHORTC(KMAX)

*-----------------------------------------------------------------------
*           constants
*-----------------------------------------------------------------------
       INTEGER C_TRIES
       REAL BIG, CUT
       PARAMETER (BIG = 1.E30,      ! arbitrary large number
     $            C_TRIES = 10,     ! # repeats of k-means
     $            CUT = 0.85)       ! limit for forming clusters

*-----------------------------------------------------------------------
*            local variables
*-----------------------------------------------------------------------
       INTEGER
     $        G, H, I,          ! count trials, count points
     $        J,                ! count dimensions
     $        K0,               ! one less than KMIN
     $        KBEST,            ! preferable value for K
     $        KTRY,             ! trial value for K
     $        L                 ! count clusters

       REAL
     $     A, FBEST,           ! measure for selecting k
     $     THRESH,             ! threshhold distance for new clusters
     $     UBEST, UPREV,       ! residual from clustering routines
     $     URET,
     $     WALL                ! total weights

       LOGICAL KNOWN,        ! number of clusters is given?
     $         MISS,         ! some data is missing
     $         OK            ! at least one successful result to return

*-----------------------------------------------------------------------
*         external functions
*-----------------------------------------------------------------------
       LOGICAL ALMEQ                 ! almost-equal function

*-----------------------------------------------------------------------
*       Begin executable statements.  Validate input.
*-----------------------------------------------------------------------
       IFAULT = 0                       ! clear error codes
       KBEST = K                        ! avoid compiler warning
       OK = .TRUE.                      ! avoid compiler warning
       IF (KMIN < 1 .OR. K < KMIN .OR. K > KMAX .OR. KMAX > N) THEN
         IFAULT = 6
         RETURN
       END IF

       IF (ROBUST .AND. SOSDO) THEN
         IFAULT = 7
         RETURN
       END IF

       IF (P < 1 .OR. N < 1) THEN       ! non-null data array?
         IFAULT = 8
         RETURN
       END IF

*        at least one element of W() must be nonzero
       DO J = 1, P
         IF (.NOT. ALMEQ(W(J), 0.)) GO TO 10
       END DO
       IFAULT = 9
       RETURN

*        check for data with no values present
  10   DO 20 I = 1, N
         DO J = 1, P
           IF (NX(J,I) > 0) GO TO 20
         END DO
         IFAULT = 10
         RETURN
  20   CONTINUE

*        look for missing data
       MISS = .FALSE.
       DO I = 1, N
         DO J = 1, P
           IF (NX(J,I) .LE. 0) THEN
             MISS = .TRUE.
             GO TO 30
           END IF
         END DO
       END DO

*        weights summing to zero prevent comparisons
  30   WALL = 0.
       DO J = 1, P
         WALL = WALL + W(J)
       END DO
       IF (ALMEQ(WALL, 0.) .AND. .NOT. MISS) THEN
         IFAULT = 11
         RETURN
       END IF

*-----------------------------------------------------------------------
*        case K = N
*-----------------------------------------------------------------------
       IF (KMIN .EQ. N) THEN
         DO I = 1, N
           DO J = 1, P
             C(J,I) = X(J,I)
             NC(J,I) = NX(J,I)
           END DO
         END DO
         DO I = 1, N
           Z(I) = I
         END DO
         DO L = 1, K
           POP(L) = 1
         END DO
         U = 0.
         RETURN
       END IF

*-----------------------------------------------------------------------
*         case K = 1
*-----------------------------------------------------------------------
       IF (KMAX .EQ. 1) THEN
         DO I = 1, N
           Z(I) = 1
         END DO
         GO TO 90                ! and compute center
       END IF

*-----------------------------------------------------------------------
*                initialize alpha and fk (see below)
*-----------------------------------------------------------------------
       KNOWN = KMIN .EQ. KMAX
       IF (KNOWN) THEN
         K0 = K
         F(K) = 1.
       ELSE
         K0 = MAX(1, KMIN-1)              ! bound of DO-loop
       END IF
       FBEST = BIG
       UBEST = 0.
       URET = 0.
       A  = 1. - 3. / (4. * P)
       DO G = 3, K0-1
         A = A + (1. - A) / 6.
       END DO

*-----------------------------------------------------------------------
*               Loop over possible values of k
*-----------------------------------------------------------------------
       DO 60 G = K0, KMAX
         UPREV = UBEST
         UBEST = BIG
         OK = .FALSE.
         DO H = 1, C_TRIES            ! seek least residual
           KTRY = G
           CALL DINIT (X, NX, PLIM, P, N, W, C, NC, KTRY, ROBUST,
     $                    SHORTX, CUMD, THRESH)
           THRESH = BIG
           CALL KCLUST (X,NX,PLIM,P,N,W,C,NC,KTRY,KTRY,KTRY,ZTRY,POP,
     $ MISS,ROBUST,SOSDO,THRESH,LIVEC,LIVEX,IPERM,TREE,CRIT,SHORTC,
     $ SHORTX,CUMD,U,IFAULT)
           IF (IFAULT .EQ. 0 .AND. U < UBEST) THEN
             UBEST = U
             OK = .TRUE.
             IF (KNOWN) KBEST = KTRY
             DO I = 1, N
               Z(I) = ZTRY(I)
             END DO
             IF (KTRY .EQ. 1) GO TO 50    ! save computation for 1 cluster
           END IF
         END DO  ! next H
         IF (.NOT. OK) THEN               ! failure
           K = G
           RETURN
         END IF

  50     IF (KNOWN) THEN        ! ignore f criterion
           K = KBEST
           U = UBEST
           GO TO 90             ! bail out of G loop
         END IF

*-----------------------------------------------------------------------
* Choose the number of clusters.  refer to:
*     "Selection of K in K-means clustering"
*      D T Pham, S S Dimov, and C D Nguyen
*      Journal of Mechanical Engineering Science, v.219, 2005 [2004]
*
* They propose the formula:
*             { 1,                            k = 1
*      f(k) = { S[k] / (a[k] @ S[k-1]),       S[k-1] .NE. 0, k > 1
*             { 1,                            S[k-1] = 0
*
* where
*             { 1 - 3 / (4 @ n),              k = 2, n > 1
*      a[k] = {
*             { a[k-1] + (1 - a[k]) / 6       k > 2, n > 1
*
* and k is the number of clusters, n is the number of dimensions, and
* S[k] is the residual from clustering with k clusters.
*-----------------------------------------------------------------------
         IF (G .EQ. 1) THEN
           A = 0.
         ELSE IF (G .EQ. 2) THEN
           A = 1. - 3. / (4. * P)
         ELSE
           A = A + (1. - A) / 6.
         END IF
         IF (G < KMIN) GO TO 60        ! don't select too few
         IF (G .EQ. 1 .OR. ALMEQ(UPREV, 0.)) THEN
           F(G) = 1.
         ELSE
           F(G) = UBEST / (A * UPREV)
         END IF

*-----------------------------------------------------------------------
* Low F indicates distinct clusters.  Accept the default number of clusters
* if the lowest value of F is not less than the CUT parameter.
*-----------------------------------------------------------------------
         IF ( (F(G) < CUT .AND. F(G) < FBEST)
     $          .OR. (FBEST .GE. CUT .AND. G .EQ. K) ) THEN
           FBEST = MIN(FBEST, F(G))
           KBEST = G
           URET = UBEST
           DO I = 1, N                 ! stash best Z for this K
             ZLOF(i) = Z(i)
           END DO
         END IF
  60   CONTINUE  ! next G

       K = KBEST                   ! set return values
       U = URET
       DO I = 1, N
         Z(i) = ZLOF(i)
       END DO

*-----------------------------------------------------------------------
*   Count cluster populations, return centers
*-----------------------------------------------------------------------
  90   IF (OK) IFAULT = 0               ! clear error codes
       CALL CNTPOP (Z, N, POP, K)       ! count populations
       DO L = 1, K                      ! all clusters live
         LIVEC(L) = 1
       END DO
       CALL CENT (X, NX, PLIM, P, N, C, NC, K, Z, POP, ROBUST,
     $                 LIVEC, IPERM, CUMD)
       RETURN
      END  ! of clue77


*-----------------------------------------------------------------------
* CLUELA - CLUstering with ELAborations
*
*___Variable________Type______I/O_______Description_____________________
*   X(PLIM,N)       Real      In        Data points in R^P
*   NX(PLIM,N)      Integer   In        1 indicates data is present for
*                                            the corresponding X value.
*                                       0 indicates data is missing for
*                                            the corresponding X value.
*   PLIM            Integer   In        Row dimension of X
*   P               Integer   In        Measurements per point
*   N               Integer   In        Number of points
*   W(P)            Real      In        Importance weights of variables
*                                              if in doubt use /P * 1./
*   C(PLIM,KMAX)    Real      Out       Cluster centroids
*   NC(PLIM,KMAX)   Integer   Out       Count of data available each center
*   KMIN            Integer   In        Fewest clusters to make
*   K               Integer   Both      Default/actual clusters to make
*   KMAX            Integer   In        Most clusters to make
*   Z(N)            Integer   Out       What cluster each point is in
*   POP(KMAX)       Integer   Out       Number of points in each cluster
*   F(KMIN:KMAX)    Real      Out       Low value indicates good # clusters
*   ROBUST          Logical   In        Use L1 statistics
*   SOSDO           Logical   In        Sum of squared distance objective
*   IWORK(4*N+25*KMAX)
*                   Integer   Neither   Workspace
*   RWORK(2*N+9*KMAX)
*                   Real      Neither   Workspace
*   U               Real      Out       Residual error
*   IFAULT          Integer   Out       Error code:
*                                          0  Normal return
*                                          1  Empty cluster
*                                          2  Max iterations exceeded
*                                          3  K out of bounds (in KCLUST)
*                                          4  Point cannot be compared
*                                          5  Formed too many clusters
*                                          6  K out of bounds (in CLUE77)
*                                          7  Conflicting options
*                                          8  Array limits out of bounds
*                                          9  All weights are zero
*                                         10  Datum with no values
*                                         11  Weights sum to zero
*-----------------------------------------------------------------------
      SUBROUTINE CLUELA (X, NX, PLIM, P, N, W, C, NC, KMIN, K, KMAX,
     $  Z, POP, F, ROBUST, SOSDO, IWORK, RWORK, U, IFAULT)

       IMPLICIT NONE
       INTEGER NX, PLIM, P, N, NC, KMIN, K, KMAX, Z, POP, IWORK, IFAULT
       REAL X, W, C, F, RWORK, U
       LOGICAL ROBUST, SOSDO
       INTEGER NS
       PARAMETER (NS = 12)
       DIMENSION X(PLIM,N), NX(PLIM,N), W(P), C(PLIM,KMAX),
     $  NC(PLIM,KMAX), Z(N), POP(KMAX), F(KMIN:KMAX),
     $  IWORK(4*N+(2*NS+1)*KMAX), RWORK(2*N+9*KMAX)

*          Partition workspace
       CALL CLUE77 (X, NX, PLIM, P, N, W, C, NC, KMIN, K, KMAX, Z, POP,
     $  F, ROBUST, SOSDO, IWORK(1), IWORK(N+1), IWORK(2*N+1),
     $  IWORK(3*N+1), IWORK(4*N+1), IWORK(4*N+KMAX+1), RWORK(1),
     $  RWORK(N+1), RWORK(2*N+1), RWORK(2*N+KMAX+1), U, IFAULT)

       RETURN
      END  ! of cluela
