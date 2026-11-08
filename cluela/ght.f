*-----------------------------------------------------------------------
* ght.f - a generalized hyperplane tree for the k-means problem
* by Andy Allinger, 2013-2020, released to the public domain
* This program may be used by any person for any purpose.
*
*-----------------------------------------------------------------------
*   The generalized hyperplane tree organizes a data set by partitioning,
* like any binary tree.  Unlike the k-d tree, there is no requirement for
* the data to be represented in R^N Euclidean space.  A data set which has
* a dissimilarity between each datum may be split by choosing two points,
* and assigning the remaining points to a partition by whether they are
* nearer the first point or the second point.  If the data and dissimilarity
* function satisfy the requirements of a metric space, then the triangle
* inequality may be used to efficiently search the tree for nearest
* neighbors of a point.
*
* Early research was:
*     I. Kalantari and G. McDonald
*     "A data structure and an algorithm for the nearest point problem"
*     IEEE Transactions on Software Engineering, v.9 p.5, 1983
*
* cited in:
*     Edgar Chavez, Gonzallo Navarro, and Jose Luis Marroquin
*     "Searching in Metric Spaces"
*     Association for Computing Machinery, 1999
*
* It was discussed also in:
*     Jeffrey K. Uhlmann
*     "Satisfying General Proximity/Similarity Queries with Metric Trees"
*     Information Processing Letters, v.40 pp.175-179, November 1991
*
*-----------------------------------------------------------------------
*     Suppose the pivots are O and P, a query point is Q.  There may be
* a hypothetical point S near the boundary surface.  The distance QS must
* be bounded to compare it to r, the radius of the ball containing the known
* neighbors of Q.  Construct triangle PQS, then by the triangle inequality:
*       QS + PS .GE. PQ             =>      QS .GE. PQ - PS
* Construct triangle OQS, then
*       QS + OQ .GE. OS             =>      QS .GE. OS - OQ
* Add
*                          2 @ QS .GE. (PQ - OQ) + (OS - PS)
* OS = PS by definition, so
*                              QS .GE. (PQ - OQ) / 2
*                               r < (PQ - OQ) / 2
*
*-----------------------------------------------------------------------
*      Also keep the farthest distance from a pivot to any object beyond
* it in the tree.  Call this distance f.  The ball around P of radius f
* has no intersection with the ball around Q with radius r if
*      r + f < PQ                =>      r < PQ - f
* F may be arbitrarily close to zero or arbitrarily large.  In the condition
* considered previously, OQ may be arbitrarily close to zero or arbitrarily
* large.  It follows that neither condition implies the other, thus both
* conditions may be used to prune the search.
*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
* ghtdsm - a dissimilarity function
*
*___Name_____Type______In/Out___Description_____________________________
*   A(P)     Real      In       Data vector
*   B(P)     Real      In       Data vector
*   P        Integer   In       Dimension of the data
*   W(P)     Real      In       Importance weights
*   ROBUST   Logical   In       Use robust statistics
*   ODWA     Real      In       1 / SUM w[j]
*   GHTDSM   Real      Out      Manhattan distance or squared Euclidean
*-----------------------------------------------------------------------
      FUNCTION GHTDSM (A, B, P, W, ROBUST, ODWA)
       IMPLICIT NONE
       INTEGER P
       REAL A(P), B(P), W(P), ODWA, GHTDSM
       LOGICAL ROBUST
       INTEGER J

*          begin
       GHTDSM = 0.
       IF (ROBUST) THEN                               ! Manhattan
         DO J = 1, P
           GHTDSM = GHTDSM + W(J) * ABS(A(J) - B(J))
         END DO
       ELSE                                           ! squared Euclidean
         DO J = 1, P
           GHTDSM = GHTDSM + W(J) * (A(J) - B(J))**2
         END DO
       END IF  ! robust ?
       GHTDSM = GHTDSM * ODWA
       RETURN
      END  ! of ghtdsm


*-----------------------------------------------------------------------
* ghtdbo - a dissimilarity function with bailout
*
*___Name_____Type______In/Out___Description_____________________________
*   A(P)     Real      In       Data vector
*   B(P)     Real      In       Data vector
*   P        Integer   In       Dimension of the data
*   W(P)     Real      In       Importance weights
*   ROBUST   Logical   In       Use robust statistics
*   ODWA     Real      In       1 / SUM w[j]
*   BAIL     Real      In       Bail out if this value is exceeded
*   GHTDBO   Real      Out      Manhattan distance or squared Euclidean
*-----------------------------------------------------------------------
      FUNCTION GHTDBO (A, B, P, W, ROBUST, ODWA, BAIL)
       IMPLICIT NONE
       INTEGER P
       REAL A(P), B(P), W(P), ODWA, BAIL, GHTDBO
       LOGICAL ROBUST
       INTEGER J

*          begin
       GHTDBO = 0.
       IF (ROBUST) THEN                               ! Manhattan
         DO J = 1, P
           GHTDBO = GHTDBO + W(J) * ABS(A(J) - B(J))
           IF (GHTDBO > BAIL) GO TO 10
         END DO
       ELSE                                           ! squared Euclidean
         DO J = 1, P
           GHTDBO = GHTDBO + W(J) * (A(J) - B(J))**2
           IF (GHTDBO > BAIL) GO TO 10
         END DO
       END IF  ! robust ?
  10   GHTDBO = GHTDBO * ODWA
       RETURN
      END  ! of ghtdbo


*-----------------------------------------------------------------------
* ghtini - initialize a generalized hyperplane tree
*
*___Name___________Type______I/O_____Description________________________
*   K              Integer   In      Number of items to be indexed
*   TREE(NS,2*K)   Integer   Both    The generalized hyperplane tree
*   LAST           Integer   Out     Last filled index in tree
*-----------------------------------------------------------------------
      SUBROUTINE GHTINI (K, TREE, LAST)
       IMPLICIT NONE
       INTEGER K, TREE, LAST
       INTEGER NS
       PARAMETER (NS = 12)
       DIMENSION TREE(NS,2*K)

*           constants
       INTEGER STEM, LEAF, BUCK0, BUCK9
       PARAMETER (STEM = 1,
     $            LEAF = 2,
     $            BUCK0 = 3,
     $            BUCK9 = NS)

*         local variables
       INTEGER J

*          begin
       LAST = 1
       TREE(STEM,LAST) = -1       ! is root
       TREE(LEAF,LAST) = -1       ! is leaf
       DO J = BUCK0, BUCK9
         TREE(J,LAST) = -1        ! no data
       END DO
       RETURN
      END  ! of ghtini


*-----------------------------------------------------------------------
* ghtins - insert a cluster center into a generalized hyperplane tree
*
*___Name___________Type______I/O____Description_________________________
*   PLIM           Integer   In     Row dimension of matrix
*   P              Integer   In     # variables each datum
*   W(P)           Real      In     Importance weights
*   C(PLIM,K)      Real      In     Center data
*   K              Integer   In     # of clusters
*   ROBUST         Logical   In     Use Manhattan distance
*   Q              Integer   In     Index into C of a center to insert
*   TREE(NS,2*K)   Integer   Both   The generalized hyperplane tree
*   CRIT(4,2*K)    Real      Out    Critical distances each node
*   LAST           Integer   Both   Last filled index in tree
*-----------------------------------------------------------------------
      SUBROUTINE GHTINS (PLIM, P, W, C, K, ROBUST, Q, TREE, CRIT, LAST)
       IMPLICIT NONE
       INTEGER PLIM, P, K, Q, TREE, LAST
       REAL W(P), C(PLIM,K), CRIT(4,2*K)
       LOGICAL ROBUST
       INTEGER NS
       PARAMETER (NS = 12)             ! node size
       DIMENSION TREE(NS,2*K)

*             constants
       INTEGER POINTS, POINTL, POINTR, PIVOTL, PIVOTR, FARL, FARR,
     $   LEAF, BUCK0, BUCK9
       REAL BIG
       PARAMETER (POINTS = 1,     ! index of stem node
     $            POINTL = 2,     ! index of left branch node
     $            POINTR = 3,     ! index of right branch node
     $            PIVOTL = 4,     ! index of left pivot datum
     $            PIVOTR = 5,     ! index of right pivot datum
     $            LEAF = 2,       ! negative one indicates a leaf node
     $            BUCK0 = 3,      ! index of first datum in bucket
     $            BUCK9 = NS,     ! index of last datum in bucket
     $            FARL = 1,       ! farthest distance from left pivot
     $            FARR = 2,       ! farthest distance from right pivot
     $            BIG = 1.E30)    ! arbitrary large number

*             local variables
       INTEGER
     $   G, H, I, J,              ! counters
     $   L,                       ! index into tree of active node
     $   NL, NR,                  ! indices of newly made nodes
     $   PL, PR                   ! indices of left and right pivots

       REAL D, DFAR, DL, DR, FL, FR,             ! distances
     $      ODWA, WALL                           ! inverse, sum of weights

       REAL GHTDSM                               ! external function

*-----------------------------------------------------------------------
*             Begin.
*-----------------------------------------------------------------------
       WALL = 0.
       DO J = 1, P
         WALL = WALL + W(J)
       END DO
       ODWA = 1. / WALL
       L = 1                                     ! the root node
  10   IF (TREE(LEAF,L) .EQ. -1) THEN            ! leaf node?
         DO H = BUCK0, BUCK9
           IF (TREE(H,L) .EQ. -1) THEN           ! space remaining in node
             TREE(H,L) = Q                       ! insert the point
             RETURN                              ! success
           END IF
         END DO
         GO TO 50                                ! node is full, so split
       END IF

*         measure
       PL = TREE(PIVOTL,L)
       PR = TREE(PIVOTR,L)
       DL = GHTDSM(C(1,PL), C(1,Q), P, W, ROBUST, ODWA)
       DR = GHTDSM(C(1,PR), C(1,Q), P, W, ROBUST, ODWA)

       IF (DL < DR .OR. (DL .EQ. DR .AND. BTEST(Q,0))) THEN
         IF (DL > CRIT(FARL,L)) CRIT(FARL,L) = DL
         L = TREE(POINTL,L)       ! climb left branch
       ELSE
         IF (DR > CRIT(FARR,L)) CRIT(FARR,L) = DR
         L = TREE(POINTR,L)       ! climb right branch
       END IF
       GO TO 10

*-----------------------------------------------------------------------
*        Split existing node; create two new leaf nodes.
*-----------------------------------------------------------------------
  50   DFAR = -BIG             ! find longest distance between points in node
       PL = TREE(BUCK0,L)      ! avoid uninitialized warning
       PR = TREE(BUCK9,L)
       DO J = BUCK0, BUCK9-1
         NR = TREE(J,L)
         DO I = J+1, BUCK9
           NL = TREE(I,L)
           D = GHTDSM(C(1,NL), C(1,NR), P, W, ROBUST, ODWA)
           IF (D > DFAR) THEN
             DFAR = D
             PL = NL
             PR = NR
           END IF
         END DO
       END DO

       LAST = LAST + 1            ! make new nodes
       NL = LAST
       TREE(POINTS,NL) = L
       TREE(LEAF,NL) = -1
       DO H = BUCK0, BUCK9
         TREE(H,NL) = -1
       END DO

       LAST = LAST + 1
       NR = LAST
       TREE(POINTS,NR) = L
       TREE(LEAF,NR) = -1
       DO H = BUCK0, BUCK9
         TREE(H,NR) = -1
       END DO

*              sort
       I = BUCK0
       J = BUCK0
       FL = 0.
       FR = 0.
       DO H = BUCK0, BUCK9
         G = TREE(H,L)
         IF (G .EQ. PL) THEN                     ! force pivot into its node
           DL = -1.0
         ELSE
           DL = GHTDSM(C(1,G), C(1,PL), P, W, ROBUST, ODWA)
         END IF
         IF (G .EQ. PR) THEN
           DR = -1.0
         ELSE
           DR = GHTDSM(C(1,G), C(1,PR), P, W, ROBUST, ODWA)
         END IF

         IF (DL < DR .OR. (DL .EQ. DR .AND. BTEST(G,0))) THEN
           TREE(I,NL) = G                        ! copy datum to left branch
           IF (DL > FL) FL = DL
           I = I + 1
         ELSE
           TREE(J,NR) = G                        ! copy datum to right branch
           IF (DR > FR) FR = DR
           J = J + 1
         END IF
       END DO

*                current node becomes a stem node
       TREE(POINTL,L) = NL
       TREE(POINTR,L) = NR
       TREE(PIVOTL,L) = PL
       TREE(PIVOTR,L) = PR
       CRIT(FARL,L) = FL
       CRIT(FARR,L) = FR
       GO TO 10
      END  ! of ghtins


*-----------------------------------------------------------------------
* ghtsel - select nearest center from a generalized hyperplane tree
*
* DIST should be set to a very large value on input, unless a distance to
* a center is already known.  Centers farther away will not be considered.
*
*___Name___________Type______I/O____Description_________________________
*   X(PLIM,N)      Real      In     Data array
*   PLIM           Integer   In     Row dimension of arrays
*   P              Integer   In     # variables each datum
*   N              Integer   In     # points in X
*   W(P)           Real      In     Importance weights of variables
*   C(PLIM,K)      Real      In     Coordinates of centers
*   K              Integer   In     Number of centers
*   ROBUST         Logical   In     Use Manhattan distance
*   SELF           Logical   In     Don't permit to select NEAR = Q
*   Q              Integer   In     Index into X of a point to query
*   LIVEC(K)       Integer   In     > 0 if each cluster is in the live set
*   LIVEX          Integer   In     > 0 if query point is in the live set
*   NEAR           Integer   Out    Index into C of the nearest center
*   DIST           Real      Both   Distance to the nearest neighbor
*   TREE(NS,2*K)   Integer   Both   The generalized hyperplane tree
*   CRIT(4,2*K)    Real      Both   Critical distances each node
*-----------------------------------------------------------------------
      SUBROUTINE GHTSEL (X, PLIM, P, N, W, C, K, ROBUST, SELF, Q,
     $                    LIVEC, LIVEX, NEAR, DIST, TREE, CRIT)
       IMPLICIT NONE
       INTEGER PLIM, P, N, K
       INTEGER Q, LIVEC(K), LIVEX, NEAR, TREE
       REAL X(PLIM,N), W(P), C(PLIM,K), DIST, CRIT(4,2*K)
       LOGICAL ROBUST, SELF
       INTEGER NS
       PARAMETER (NS = 12)
       DIMENSION TREE(NS,2*K)

*             constants
       INTEGER POINTS, POINTL, POINTR, PIVOTL, PIVOTR, FARL, FARR, MARK,
     $   DIR, KEEPY, KEEPZ, LEAF, BUCK0, BUCK9
       REAL SMIDGE, SLTH, SLTO
       PARAMETER (POINTS = 1,         ! index of stem node
     $            POINTL = 2,         ! index of left branch node
     $            POINTR = 3,         ! index of right branch node
     $            PIVOTL = 4,         ! index of left pivot datum
     $            PIVOTR = 5,         ! index of right pivot datum
     $            MARK = 6,           ! which branches have been visited
     $            DIR = 7,            ! which branch is nearer the test point
     $            LEAF = 2,           ! -1 here indicates a leaf node
     $            BUCK0 = 3,          ! index of first datum in bucket
     $            BUCK9 = NS,         ! index of last datum in bucket
     $            FARL = 1,           ! farthest distance from left pivot
     $            FARR = 2,           ! farthest distance from right pivot
     $            KEEPY = 3,          ! PQ - f
     $            KEEPZ = 4,          ! (1/2) (PQ - OQ)
     $            SMIDGE = 1.E-4)
       PARAMETER (SLTH = 1. / (2. + SMIDGE),   ! smidge less than half
     $            SLTO = 1. / (1. + SMIDGE))   ! smidge less than one

*                local variables
       INTEGER
     $   I, J,                     ! counters
     $   L,                        ! index into tree of active node
     $   NB, FB,                   ! indices of near branch and far branch
     $   PL, PR                    ! indices of left and right pivots

       REAL
     $   BAIL,       ! sum that exceeds shortest distance
     $   D2,         ! possibly squared distance
     $   DF,         ! a far distance
     $   DL, DR,     ! Metric (not squared) distances
     $   ODWA,       ! inverse of sum of weights
     $   PQ,         ! distance from P to Q
     $   R,          ! radius around query point to search
     $   WALL,       ! sum of weights
     $   Y,          ! distance from query to point in far partition
     $   Z           ! distance from query to the generalized hyperplane

*           external functions
       REAL GHTDSM, GHTDBO

*-----------------------------------------------------------------------
*           Begin.
*-----------------------------------------------------------------------
       WALL = 0.
       DO J = 1, P
         WALL = WALL + W(J)
       END DO
       BAIL = DIST * WALL
       ODWA = 1. / WALL
       IF (ROBUST) THEN
         R = DIST
       ELSE
         R = SQRT(DIST)
       END IF
       L = 1              ! begin at root

*-----------------------------------------------------------------------
*         Come here when climbing.
*-----------------------------------------------------------------------
  10   IF (TREE(LEAF,L) .EQ. -1) THEN            ! inspect a leaf
         DO 20 J = BUCK0, BUCK9
           I = TREE(J,L)
           IF (I .EQ. -1) GO TO 30               ! reached bottom of bucket
           IF (SELF .AND. I .EQ. Q) GO TO 20     ! don't select zero distance
           IF (LIVEX .LE. 0 .AND. LIVEC(I) .LE. 0) GO TO 20  ! skip dead
           D2 = GHTDBO(X(1,Q), C(1,I), P, W, ROBUST, ODWA, BAIL)
           IF (D2 < DIST) THEN
             NEAR = I
             DIST = D2
             BAIL = DIST * WALL
             IF (ROBUST) THEN
               R = DIST
             ELSE
               R = SQRT(DIST)
             END IF
           END IF      ! shorter distance?
  20     CONTINUE      ! next datum in leaf

  30     L = TREE(POINTS,L)                      ! back to stem node

       ELSE
         TREE(MARK,L) = 0                        ! clear mark in stem node
         TREE(DIR,L) = 0
       END IF

*-----------------------------------------------------------------------
*         Come here when backing down.
*-----------------------------------------------------------------------
  50   IF (L .EQ. -1) RETURN                     ! is this the root node?

*                measure as required
       IF (TREE(DIR,L) .EQ. -1) THEN             ! left branch nearer
         NB = POINTL
         FB = POINTR
         Y = CRIT(KEEPY,L)
         Z = CRIT(KEEPZ,L)

       ELSE IF (TREE(DIR,L) .EQ. 1) THEN         ! right branch nearer
         NB = POINTR
         FB = POINTL
         Y = CRIT(KEEPY,L)
         Z = CRIT(KEEPZ,L)

       ELSE                                      ! unknown
         PL = TREE(PIVOTL,L)
         PR = TREE(PIVOTR,L)
         DL = GHTDSM(X(1,Q), C(1,PL), P, W, ROBUST, ODWA)
         DR = GHTDSM(X(1,Q), C(1,PR), P, W, ROBUST, ODWA)
         IF (.NOT. ROBUST) THEN
           DL = SQRT(DL)
           DR = SQRT(DR)
         END IF
         PQ = MAX(DL, DR)
         Z = ABS(DL - DR) * SLTH
         CRIT(KEEPZ,L) = Z     ! remember

         IF (DL < DR .OR. (DL .EQ. DR .AND. BTEST(Q,0))) THEN
           NB = POINTL
           FB = POINTR
           TREE(DIR,L) = -1         ! remember
           DF = CRIT(FARR,L)
         ELSE
           NB = POINTR
           FB = POINTL
           TREE(DIR,L) = +1
           DF = CRIT(FARL,L)
         END IF
         IF (.NOT. ROBUST) DF = SQRT(DF)
         Y = (PQ - DF) * SLTO
         CRIT(KEEPY,L) = Y
       END IF               ! need to measure ?

*             Explore the near branch
       IF (TREE(MARK,L) .EQ. 0) THEN             ! near branch unmarked?
         TREE(MARK,L) = 1                        ! mark
         L = TREE(NB,L)
         GO TO 10                                ! and inspect the node
       END IF

*             Consider far branch.
       IF (R > Y .AND. R > Z .AND. TREE(MARK,L) .EQ. 1) THEN
         TREE(MARK,L) = 2                        ! mark
         L = TREE(FB,L)
         GO TO 10                                ! and inspect the node
       END IF

       L = TREE(POINTS,L)                        ! backtrack
       GO TO 50
      END  ! of ghtsel
