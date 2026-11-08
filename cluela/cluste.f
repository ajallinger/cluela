*-----------------------------------------------------------------------
* cluste.f - CLUstering with Superior Theoretical Efficiency
* by Andy Allinger, 2017, released to the public domain
* This program may be used by any person for any purpose.
*-----------------------------------------------------------------------

*=======================================================================
* hclust - Hyperplane partitioning
* Divisive clustering based on partial sorting, a reduced form of
* Anthony Hoare's quick sort algorithm.
*
*___Name_________Type______I/O____Description___________________________
*   X(PLIM,N)    Real      In     Data array
*   NX(PLIM,N)   Integer   In     1 => data present for this measurement
*   PLIM         Integer   In     Row dimension of X
*   P            Integer   In     Dimension of data
*   N            Integer   In     Number of data points
*   W(P)         Real      In     Importance weights
*   K            Integer   In     Desired number of clusters
*   Z(N)         Integer   Out    What cluster each point is in
*   POP(K)       Integer   Out    Cluster populations
*   ROBUST       Logical   In     Use L1 distance?
*   IND(N)       Integer   None   Index array that sorts X into
*                                     output order from input order
*   JND(K)       Integer   None   Sorting workspace
*   BOUND(K)     Integer   None   Work array
*=======================================================================
      SUBROUTINE HCLUST (X, NX, PLIM, P, N, W, K, Z, POP, ROBUST,
     $                    IND, JND, BOUND)

       IMPLICIT NONE
       INTEGER NX, PLIM, P, N, K, Z, POP, IND, JND, BOUND
       REAL X, W
       LOGICAL ROBUST
       DIMENSION X(PLIM,N), NX(PLIM,N), W(P), Z(N), POP(K),
     $  IND(N), JND(K), BOUND(K)

*          Local variables
       INTEGER
     $          I, L,              ! counters
     $          MAXPOP,            ! population of largest partition
     $          LL,                ! index of largest population
     $          LO, HI,            ! indices at end of current partition
     $          IHI, ILO,          ! indices for sorting current partition
     $          KLO, KMID, KHI,    ! indices of randomly chosen pivots
     $          ISWAP              ! temporary scalar

       REAL DHL, DHM, DML          ! dissimilarities hi-lo, hi-mid, mid-lo

*        External functions
       REAL DSM
       INTEGER HAT

*         initialize
       DO I = 1, N
         IND(I) = I
       END DO
       DO L = 1, K
         BOUND(L) = N
       END DO
       POP(1) = N
       LL = 1

*=======================================================================
*       Partition a section of the array such that elements closer to the first
*       element are moved to the front, and elements closer to the last element
*       are moved to the rear.  Continue until enough clusters have been made.
*=======================================================================
       DO L = 2, K
         MAXPOP = 0              ! find partition with greatest population
         DO I = 1, L-1
           IF (POP(I) > MAXPOP) THEN
             MAXPOP = POP(I)
             LL = I
           END IF
         END DO

         HI = BOUND(LL)               ! find limits of partition
         LO = HI - MAXPOP + 1
         IHI = HI - 1
         ILO = LO + 1
         IF (MAXPOP .EQ. 2) GO TO 70

*          swap random pivot values to ends of partition (introduce chance)
         KHI = HAT(MAXPOP) - 1 + LO
  20     KLO = HAT(MAXPOP) - 1 + LO
         IF (KLO .EQ. KHI) GO TO 20

*           a third pivot value in case first two are near each other
  30     KMID = HAT(MAXPOP) - 1 + LO
         IF (KMID .EQ. KHI .OR. KMID .EQ. KLO) GO TO 30
         DHL = DSM (X(1,IND(KHI)), NX(1,IND(KHI)),
     $           X(1,IND(KLO)), NX(1,IND(KLO)), P, W, ROBUST, 0)
         DHM = DSM (X(1,IND(KHI)), NX(1,IND(KHI)),
     $           X(1,IND(KMID)), NX(1,IND(KMID)), P, W, ROBUST, 0)
         DML = DSM (X(1,IND(KMID)), NX(1,IND(KMID)),
     $           X(1,IND(KLO)), NX(1,IND(KLO)), P, W, ROBUST, 0)

         IF (DHM > DHL .OR. DML > DHL) THEN
           IF (DHM > DML) THEN
             KLO = KMID
           ELSE
             KHI = KMID
           END IF
         END IF

         ISWAP = IND(KHI)
         IND(KHI) = IND(HI)
         IND(HI) = ISWAP
         ISWAP = IND(KLO)
         IND(KLO) = IND(LO)
         IND(LO) = ISWAP

*=======================================================================
*            sort
*=======================================================================
  50     IF (ILO .LE. IHI) THEN  ! look for object that belongs in HI side
           IF (DSM(X(1,IND(ILO)), NX(1,IND(ILO)),
     $           X(1,IND(LO)), NX(1,IND(LO)),
     $           P, W, ROBUST, 0) .LE.
     $         DSM(X(1,IND(ILO)), NX(1,IND(ILO)),
     $           X(1,IND(HI)), NX(1,IND(HI)),
     $           P, W, ROBUST, 0)) THEN
             ILO = ILO + 1
             GO TO 50
           END IF
         END IF

  60     IF (IHI .GE. ILO) THEN  ! look for object that belongs in LO side
           IF (DSM(X(1,IND(IHI)), NX(1,IND(IHI)),
     $           X(1,IND(HI)), NX(1,IND(HI)),
     $           P, W, ROBUST, 0) .LE.
     $         DSM(X(1,IND(IHI)), NX(1,IND(IHI)),
     $           X(1,IND(LO)), NX(1,IND(LO)),
     $           P, W, ROBUST, 0)) THEN
             IHI = IHI - 1
             GO TO 60
           END IF
         END IF

         IF (ILO < IHI) THEN         ! swap these objects
           ISWAP = IND(ILO)
           IND(ILO) = IND(IHI)
           IND(IHI) = ISWAP
           ILO = ILO + 1
           IHI = IHI - 1
           GO TO 50
         END IF

*=======================================================================
  70   CONTINUE      ! Add this partition boundary to the bounds array
*=======================================================================
*             equal split if all distances were the same
         IF (ILO .EQ. HI) THEN
           IHI = (LO + HI) / 2
           ILO = IHI + 1
         END IF
         POP(LL) = HI - IHI           ! reduce population
         POP(L) = ILO - LO            ! population in new partition
         BOUND(L) = IHI               ! append
       END DO     ! next partition

*=======================================================================
*  Recover memberships and populations from the index array
*=======================================================================
       CALL QSORTI (BOUND, JND, K)
       CALL IORDER (BOUND, JND, K)
       L = 1
       DO I = 1, N
         IF (I > BOUND(L)) L = L + 1
         Z(IND(I)) = L
       END DO

       I = 0
       DO L = 1, K
         POP(L) = BOUND(L) - I
         I = BOUND(L)
       END DO
       RETURN
      END  ! of hclust


*-----------------------------------------------------------------------
* Speedy clustering with theoretical efficiency less than O(iNKP)
* if there is NO MISSING DATA, NO SUM-OF-SQUARES MINIMIZATION, and
* NO NEGATIVE WEIGHTS.  These conditions are not checked!
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
*   ROBUST            Logical   In        Use L1 statistics
*   SOSDO             Logical   In
*   IWORK(2*N+25*KMAX)
*                     Integer   Neither   Workspace
*   RWORK(2*N+9*KMAX)
*                     Real      Neither   Shortest distance each point
*   U                 Real      Out       Residual error
*   IFAULT            Integer   Out       Error code:
*-----------------------------------------------------------------------
      SUBROUTINE CLUSTE (X, NX, PLIM, P, N, W, C, NC, KMIN, K, KMAX, Z,
     $                    POP, ROBUST, SOSDO, IWORK, RWORK, U, IFAULT)

       IMPLICIT NONE
       INTEGER NX, PLIM, P, N, NC, KMIN, K, KMAX, Z, POP, IWORK, IFAULT
       REAL X, W, C, RWORK, U
       LOGICAL ROBUST, SOSDO
       INTEGER NS
       PARAMETER (NS = 12)
       DIMENSION X(PLIM,N), NX(PLIM,N), W(P), C(PLIM,KMAX),
     $  NC(PLIM,KMAX), Z(N), POP(KMAX),
     $  IWORK(2*N+(2*NS+1)*KMAX), RWORK(2*N+9*KMAX)

*            local variables
       INTEGER
     $        I,                ! count points
     $        J,                ! count dimensions
     $        L                 ! count clusters

       REAL
     $     D,                  ! a distance
     $     THRESH,             ! threshhold distance for new clusters
     $     WALL                ! sum of weights

       LOGICAL
     $         MISS            ! missing data?

*           External functions
       REAL DSM              ! dissimilarity
       LOGICAL ALMEQ         ! almost equal

*-----------------------------------------------------------------------
*       Begin executable statements.  Validate input.
*-----------------------------------------------------------------------
       IFAULT = 0                       ! clear error codes
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

*         look for missing data
       MISS = .FALSE.
       DO I = 1, N
         DO J = 1, P
           IF (NX(J,I) .LE. 0) THEN
             MISS = .TRUE.
             GO TO 30
           END IF
         END DO
       END DO

  30   WALL = 0.
       DO J = 1, P
         WALL = WALL + W(J)
       END DO
       IF (ALMEQ(WALL, 0.) .AND. .NOT. MISS) THEN
         IFAULT = 11
         RETURN
       END IF

*-----------------------------------------------------------------------
*           Choose initial centers
*-----------------------------------------------------------------------
       CALL HCLUST (X, NX, PLIM, P, N, W, K, Z, POP, ROBUST,
     $    IWORK(1), IWORK(N+1), IWORK(N+K+1))
       DO L = 1, K
         IWORK(L) = 1
       END DO
       CALL CENT (X, NX, PLIM, P, N, C, NC, K, Z, POP,
     $    ROBUST, IWORK(1), IWORK(K+1), RWORK(1))
       THRESH = 0.
       DO I = 1, N
         L = Z(I)
         D = DSM (X(1,I), NX(1,I), C(1,L), NC(1,L),
     $                P, W, ROBUST, 0)
         IF (D > THRESH) THRESH = D
       END DO
       CALL KCLUST (X, NX, PLIM, P, N, W, C, NC, KMIN, K, KMAX,
     $    Z, POP, MISS, ROBUST, SOSDO, THRESH, IWORK(1),
     $    IWORK(KMAX+1), IWORK(N+KMAX+1), IWORK(2*N+KMAX+1),
     $    RWORK(2*N+KMAX+1), RWORK(2*N+1), RWORK(1), RWORK(N+1), U,
     $    IFAULT)
       RETURN
      END  ! of cluste
