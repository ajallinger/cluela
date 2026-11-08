*-----------------------------------------------------------------------
* KMNS - data clustering by the K-Means algorithm
* by J.A. Hartigan & M.A. Wong
*
* Code in this file based on Applied Statistics algorithms
* (C) Royal Statistical Society 1979
*
* According to Statlib (http://lib.stat.cmu.edu/apstat/):
* "The Royal Statistical Society holds the copyright to these routines,
* but has given its permission for their distribution provided that
* no fee is charged."
*
* added comments and declarations:  AJA, 2011-2017
*
*-----------------------------------------------------------------------
C
C     ALGORITHM AS 136  APPL. STATIST. (1979) VOL.28, NO.1
C
C     Divide M points in N-dimensional space into K clusters so that
C     the within cluster sum of squares is minimized.
C
*                         TRANSFER LIST
*___Name_______Type______I/O_______Description__________________________
*   A(M,N)     Real      In        the data matrix
*   M          Integer   In        the number of points
*   N          Integer   In        the number of dimensions
*   C(K,N)     Real      input:    the matrix of initial cluster centers
*                        output:   the matrix of final cluster centers
*   K          Integer   In        the number of clusters
*   IC1(M)     Integer   Out       the cluster each point belongs to
*   IC2(M)     Integer   Neither   this array is used to remember the cluster
*                                    which each point is most likely to be
*                                    tranferred to at each step
*   NC(K)      Integer   Out       the number of points in each cluster
*   AN1(K)     Real      Neither   Ratio that corrects a square distance
*                                    to a center if a point were removed
*   AN2(K)     Real      Neither   Ratio that corrects a square distance
*                                    to a center if a point were added
*   NCP(K)     Integer   Neither   Step the cluster is no longer updated
*   D(M)       Real      Neither   Distance from point to its nearest center
*   ITRAN(K)   Integer   Neither   Relates to if there is convergence
*   LIVE(K)    Integer   Neither   Relates to if there is convergence
*   ITER       Integer   In        the maximum number of iterations allowed
*   WSS(K)     Real      Out       the within-cluster sum of squares
*                                    of each cluster
*   IFAULT     Integer   Out       Error code
*
*                        FAULT DIAGNOSTICS
*   IFAULT = 0   No fault
*   IFAULT = 1   At least one cluster is empty after the initial assignment.
*                (A better set of initial cluster centres is called for)
*   IFAULT = 2   The allowed maximum number of iterations is exceeded
*   IFAULT = 3   K is less than or equal to 1 or greater than or equal to M
*
*-----------------------------------------------------------------------
      SUBROUTINE KMNS(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,
     $     ITRAN, LIVE, ITER, WSS, IFAULT)
      IMPLICIT NONE
      INTEGER M,N,K,ITER,IFAULT
      INTEGER IC1(M), IC2(M), NC(K), NCP(K), ITRAN(K), LIVE(K)
      REAL    A(M,N), D(M), C(K,N), AN1(K), AN2(K), WSS(K), DT(2)
      INTEGER I,J,L,INDX,II,IJ,IL
      REAL AA, DA, DB, DC, TEMP
      REAL    BIG, ZERO, ONE

C     Define BIG to be a very large positive number

      DATA BIG /1.E30/, ZERO /0.0/, ONE /1.0/

*-----------------------------------------------------------------------
* Begin.  Validate K, return error code 3 for K < 2 or K > M-1
*-----------------------------------------------------------------------
      IFAULT = 3
      IF (K .LE. 1 .OR. K .GE. M) RETURN
      IFAULT = 0

*-----------------------------------------------------------------------
* Initial values.  C() must have initial values already assigned.
*-----------------------------------------------------------------------
C     For each point I, find its two closest centres, IC1(I) and
C     IC2(I).     Assign it to IC1(I).

      DO 50 I = 1, M
      IC1(I) = 1
      IC2(I) = 2
      DO 10 IL = 1, 2  ! square distance from centroids 1 and 2 to point
        DT(IL) = ZERO
        DO 10 J = 1, N
          DA = A(I,J) - C(IL,J)
          DT(IL) = DT(IL) + DA*DA
   10   CONTINUE
      IF (DT(1) .GT. DT(2)) THEN  ! swap the square distances if appropriate
        IC1(I) = 2
        IC2(I) = 1
        TEMP = DT(1)
        DT(1) = DT(2)
        DT(2) = TEMP
      END IF
      DO 50 L = 3, K          ! compare the point to the centroids 3...K
        DB = ZERO
        DO 30 J = 1, N
          DC = A(I,J) - C(L,J)
          DB = DB + DC*DC
          IF (DB .GE. DT(2)) GO TO 50
   30     CONTINUE
        IF (DB .LT. DT(1)) GO TO 40
        DT(2) = DB                           ! better point than IC2
        IC2(I) = L
        GO TO 50
   40     DT(2) = DT(1)                      ! better point than IC1
        IC2(I) = IC1(I)
        DT(1) = DB
        IC1(I) = L
   50 CONTINUE

C     Update cluster centres to be the average of points contained
C     within them.

      DO 70 L = 1, K                      ! set previous values to zero
      NC(L) = 0
      DO 60 J = 1, N
   60   C(L,J) = ZERO
   70 CONTINUE
      DO 90 I = 1, M
      L = IC1(I)                  ! index of cluster nearest point I
      NC(L) = NC(L) + 1           ! increment count of points in cluster
      DO 80 J = 1, N
   80   C(L,J) = C(L,J) + A(I,J)
   90 CONTINUE

C     Check to see if there is any empty cluster at this stage

      DO 120 L = 1, K
      IF (NC(L) .EQ. 0) THEN
        IFAULT = 1
        RETURN
      END IF
      AA = NC(L)
      DO 110 J = 1, N
  110   C(L,J) = C(L,J) / AA

C     Initialize AN1, AN2, ITRAN & NCP
C     AN1(L) = NC(L) / (NC(L) - 1)
C     AN2(L) = NC(L) / (NC(L) + 1)
C     ITRAN(L) = 1 if cluster L is updated in the quick-transfer stage,
C              = 0 otherwise
C     In the optimal-transfer stage, NCP(L) stores the step at which
C     cluster L is last updated.
C     In the quick-transfer stage, NCP(L) stores the step at which
C     cluster L is last updated plus M.

      AN2(L) = AA / (AA + ONE)
      AN1(L) = BIG
      IF (AA .GT. ONE) AN1(L) = AA / (AA - ONE)
      ITRAN(L) = 1
      NCP(L) = -1
  120 CONTINUE
      INDX = 0

*-----------------------------------------------------------------------
*         Main loop
*-----------------------------------------------------------------------
      DO 140 IJ = 1, ITER

C     In this stage, there is only one pass through the data.   Each
C     point is re-allocated, if necessary, to the cluster that will
C     induce the maximum reduction in within-cluster sum of squares.

      CALL OPTRA(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,
     $        ITRAN, LIVE, INDX)

C     Stop if no transfer took place in the last M optimal transfer
C     steps.

      IF (INDX .EQ. M) GO TO 150

C     Each point is tested in turn to see if it should be re-allocated
C     to the cluster to which it is most likely to be transferred,
C     IC2(I), from its present cluster, IC1(I).   Loop through the
C     data until no further change is to take place.

      CALL QTRAN(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,
     $       ITRAN, INDX)

C     If there are only two clusters, there is no need to re-enter the
C     optimal transfer stage.

      IF (K .EQ. 2) GO TO 150

C     NCP has to be set to 0 before entering OPTRA.

      DO 130 L = 1, K
  130   NCP(L) = 0
  140 CONTINUE  ! end of main loop

C     Since the specified number of iterations has been exceeded, set
C     IFAULT = 2.   This may indicate unforeseen looping.

      IFAULT = 2

*-----------------------------------------------------------------------
*     Compute final centers
*-----------------------------------------------------------------------
  150 DO 160 L = 1, K
      WSS(L) = ZERO
      DO 160 J = 1, N
        C(L,J) = ZERO
  160 CONTINUE
      DO 170 I = 1, M
      II = IC1(I)
      DO 170 J = 1, N
        C(II,J) = C(II,J) + A(I,J)
  170 CONTINUE
      DO 190 J = 1, N
      DO 180 L = 1, K
  180   C(L,J) = C(L,J) / FLOAT(NC(L))

C     Compute within-cluster sum of squares for each cluster.

      DO 190 I = 1, M
        II = IC1(I)
        DA = A(I,J) - C(II,J)
        WSS(II) = WSS(II) + DA*DA
  190 CONTINUE

      RETURN
      END  ! of KMNS


*-----------------------------------------------------------------------
C     ALGORITHM AS 136.1  APPL. STATIST. (1979) VOL.28, NO.1
C
C     This is the optimal transfer stage.
C
C     Each point is re-allocated, if necessary, to the cluster that
C     will induce a maximum reduction in the within-cluster sum of
C     squares.
*-----------------------------------------------------------------------
      SUBROUTINE OPTRA(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,
     $      ITRAN, LIVE, INDX)
      IMPLICIT NONE
      INTEGER M,N,K,INDX
      INTEGER IC1(M), IC2(M), NC(K), NCP(K), ITRAN(K), LIVE(K)
      REAL    A(M,N), D(M), C(K,N), AN1(K), AN2(K)

      INTEGER L,I,L1,L2,LL,J
      REAL DE,DF,DA,DB,R2,RR,DC,DD,AL1,ALW,AL2,ALT
      REAL BIG, ZERO, ONE

C     Define BIG to be a very large positive number.

      DATA BIG /1.0E30/, ZERO /0.0/, ONE/1.0/

C     If cluster L is updated in the last quick-transfer stage, it
C     belongs to the live set throughout this stage.   Otherwise, at
C     each step, it is not in the live set if it has not been updated
C     in the last M optimal transfer steps.

      DO 10 L = 1, K
      IF (ITRAN(L) .EQ. 1) LIVE(L) = M + 1
   10 CONTINUE
      DO 100 I = 1, M
      INDX = INDX + 1
      L1 = IC1(I)               ! index of nearest cluster
      L2 = IC2(I)               ! index of second-nearest cluster
      LL = L2

C     If point I is the only member of cluster L1, no transfer.

      IF (NC(L1) .EQ. 1) GO TO 90

C     If L1 has not yet been updated in this stage, no need to
C     re-compute D(I).

      IF (NCP(L1) .EQ. 0) GO TO 30
      DE = ZERO
      DO 20 J = 1, N
        DF = A(I,J) - C(L1,J)
        DE = DE + DF*DF
   20   CONTINUE
      D(I) = DE * AN1(L1)   ! weighted sum of squares from point I to L1

C     Find the cluster with minimum R2.

   30   DA = ZERO
      DO 40 J = 1, N
        DB = A(I,J) - C(L2,J)
        DA = DA + DB*DB
   40   CONTINUE
      R2 = DA * AN2(L2)     ! weighted sum of squares from point I to L2
      DO 60 L = 1, K

C     If I >= LIVE(L1), then L1 is not in the live set.   If this is
C     true, we only need to consider clusters that are in the live set
C     for possible transfer of point I.   Otherwise, we need to consider
C     all possible clusters.

        IF (I .GE. LIVE(L1) .AND. I .GE. LIVE(L) .OR. L .EQ. L1 .OR.
     $        L .EQ. LL) GO TO 60
        RR = R2 / AN2(L)
        DC = ZERO
        DO 50 J = 1, N
          DD = A(I,J) - C(L,J)
          DC = DC + DD*DD
          IF (DC .GE. RR) GO TO 60
   50     CONTINUE
        R2 = DC * AN2(L)           ! a new second-nearest cluster
        L2 = L
   60     CONTINUE
        IF (R2 .LT. D(I)) GO TO 70

C     If no transfer is necessary, L2 is the new IC2(I).

        IC2(I) = L2
        GO TO 90

C     Update cluster centres, LIVE, NCP, AN1 & AN2 for clusters L1 and
C     L2, and update IC1(I) & IC2(I).

   70     INDX = 0
        LIVE(L1) = M + I
        LIVE(L2) = M + I
        NCP(L1) = I
        NCP(L2) = I
        AL1 = NC(L1)
        ALW = AL1 - ONE
        AL2 = NC(L2)
        ALT = AL2 + ONE
        DO 80 J = 1, N
          C(L1,J) = (C(L1,J) * AL1 - A(I,J)) / ALW ! change center of L1
          C(L2,J) = (C(L2,J) * AL2 + A(I,J)) / ALT ! change center of L2
   80     CONTINUE
        NC(L1) = NC(L1) - 1       ! change number of points in L1 and L2
        NC(L2) = NC(L2) + 1
        AN2(L1) = ALW / AL1       ! change AN1 and AN2
        AN1(L1) = BIG
        IF (ALW .GT. ONE) AN1(L1) = ALW / (ALW - ONE)
        AN1(L2) = ALT / AL2
        AN2(L2) = ALT / (ALT + ONE)
        IC1(I) = L2             ! update the nearest clusters to point I
        IC2(I) = L1
   90   CONTINUE
      IF (INDX .EQ. M) RETURN
  100 CONTINUE
      DO 110 L = 1, K

C     ITRAN(L) = 0 before entering QTRAN.   Also, LIVE(L) has to be
C     decreased by M before re-entering OPTRA.

      ITRAN(L) = 0
      LIVE(L) = LIVE(L) - M
  110 CONTINUE

      RETURN
      END  ! of OPTRA


*-----------------------------------------------------------------------
C     ALGORITHM AS 136.2  APPL. STATIST. (1979) VOL.28, NO.1
C
C     This is the quick transfer stage.
C     IC1(I) is the cluster which point I belongs to.
C     IC2(I) is the cluster which point I is most likely to be
C         transferred to.
C     For each point I, IC1(I) & IC2(I) are switched, if necessary, to
C     reduce within-cluster sum of squares.  The cluster centres are
C     updated after each step.
*-----------------------------------------------------------------------
      SUBROUTINE QTRAN(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,
     $    ITRAN, INDX)
      IMPLICIT NONE
      INTEGER M,N,K,INDX
      INTEGER IC1(M), IC2(M), NC(K), NCP(K), ITRAN(K)
      REAL    A(M,N), D(M), C(K,N), AN1(K), AN2(K)

      INTEGER ICOUN,ISTEP,I,L1,L2,J
      REAL DA,DB,DD,AL1,ALW,AL2,ALT,R2,DE
      REAL BIG, ZERO, ONE

C     Define BIG to be a very large positive number

      DATA BIG /1.0E30/, ZERO /0.0/, ONE /1.0/

C     In the optimal transfer stage, NCP(L) indicates the step at which
C     cluster L is last updated.   In the quick transfer stage, NCP(L)
C     is equal to the step at which cluster L is last updated plus M.

      ICOUN = 0
      ISTEP = 0
   10 DO 70 I = 1, M
      ICOUN = ICOUN + 1
      ISTEP = ISTEP + 1
      L1 = IC1(I)
      L2 = IC2(I)

C     If point I is the only member of cluster L1, no transfer.

      IF (NC(L1) .EQ. 1) GO TO 60

C     If ISTEP > NCP(L1), no need to re-compute distance from point I to
C     cluster L1.   Note that if cluster L1 is last updated exactly M
C     steps ago, we still need to compute the distance from point I to
C     cluster L1.

      IF (ISTEP .GT. NCP(L1)) GO TO 30
      DA = ZERO
      DO 20 J = 1, N
        DB = A(I,J) - C(L1,J)
        DA = DA + DB*DB
   20   CONTINUE
      D(I) = DA * AN1(L1)

C     If ISTEP >= both NCP(L1) & NCP(L2) there will be no transfer of
C     point I at this step.

   30   IF (ISTEP .GE. NCP(L1) .AND. ISTEP .GE. NCP(L2)) GO TO 60
      R2 = D(I) / AN2(L2)
      DD = ZERO
      DO 40 J = 1, N
        DE = A(I,J) - C(L2,J)
        DD = DD + DE*DE
        IF (DD .GE. R2) GO TO 60
   40   CONTINUE

C     Update cluster centres, NCP, NC, ITRAN, AN1 & AN2 for clusters
C     L1 & L2.   Also update IC1(I) & IC2(I).   Note that if any
C     updating occurs in this stage, INDX is set back to 0.

      ICOUN = 0
      INDX = 0
      ITRAN(L1) = 1
      ITRAN(L2) = 1
      NCP(L1) = ISTEP + M
      NCP(L2) = ISTEP + M
      AL1 = NC(L1)
      ALW = AL1 - ONE
      AL2 = NC(L2)
      ALT = AL2 + ONE
      DO 50 J = 1, N                                !  new centers
        C(L1,J) = (C(L1,J) * AL1 - A(I,J)) / ALW
        C(L2,J) = (C(L2,J) * AL2 + A(I,J)) / ALT
   50   CONTINUE
      NC(L1) = NC(L1) - 1
      NC(L2) = NC(L2) + 1
      AN2(L1) = ALW / AL1                           ! update AN1 and AN2
      AN1(L1) = BIG
      IF (ALW .GT. ONE) AN1(L1) = ALW / (ALW - ONE)
      AN1(L2) = ALT / AL2
      AN2(L2) = ALT / (ALT + ONE)
      IC1(I) = L2
      IC2(I) = L1

C     If no re-allocation took place in the last M steps, return.

   60   IF (ICOUN .EQ. M) RETURN
   70 CONTINUE
      GO TO 10
      END  ! of QTRAN
*+++++++++++++++++++++ End of file kmns.f ++++++++++++++++++++++++++++++
