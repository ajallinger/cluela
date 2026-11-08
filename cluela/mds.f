*-----------------------------------------------------------------------
* mds.f - multi-dimensional scaling and principal component analysis
* by Andy Allinger, 2013-2020, released to the public domain
* This program may be used by any person for any purpose.
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
* mds - Multi-Dimensional Scaling
* Given a data set, find an embedding with the same distances.
*
*___Name_________Type______In/Out____Description________________________
*   X(PLIM,N)    Real      In        Data of P dimensions by N objects
*   NX(PLIM,N)   Integer   In        1 indicates data is present for
*                                       the correspoding X value
*                                    0 indicates data is missing for
*                                       the corresponding X value
*   PLIM         Integer   In        Array bound of X and NX
*   P            Integer   In        Number of variables
*   N            Integer   In        Number of objects
*   W(P)         Real      In        Importance weights
*   CAT(P)       Integer   In        Number of categories each variable
*                                      0 indicates a continuous variable
*                                      1 indicates a dichotomous variable
*                                        values {0.0, 1.0} for absent, present
*                                     .GE. 2 indicates a categorical variable
*                                         values {1.0, 2.0, ...}
*                                     -1 indicates a cyclic variable,
*                                             period of one
*   ROBUST       Logical   In        Use Manhattan dissimilarity
*   V(NLIM,N)    Real      Out       Embedded data matrix, N by M
*                                       Transposed relative to X!
*   NLIM         Integer   In        Array bound of V
*   M            Integer   Out       Number of directions found.  M .LE. N
*   E(N)         Real      Out       Eigenvalues
*   T(N)         Real      Neither   Work array
*   IFAULT       Integer   Out       Error code
*
*___Return values_______________________________________________________
*    0   no errors
*   +n   eigenvectors incomplete
*   -1   error in eigensolver
*   -2   N or P less than 1
*   -3   fractional category value
*   -4   category value out of bounds, dichotomous
*   -5   category value out of bounds, nominal
*   -6   all zero weights
*   -7   bad value in CAT()
*   -8   no comparisons possible
*-----------------------------------------------------------------------
      SUBROUTINE MDS (X, NX, PLIM, P, N, W, CAT, ROBUST,
     $                 V, NLIM, M, E, T, IFAULT)
       IMPLICIT NONE
       INTEGER PLIM, P, N, NLIM
       INTEGER NX(PLIM,N), CAT(P), M, IFAULT
       REAL X(PLIM,N), W(P), V(NLIM,N), E(N), T(N)
       LOGICAL ROBUST

       INTEGER J, K            ! counters
       REAL
     $      S,                    ! a temp scalar
     $      TOT                   ! a total


*-----------------------------------------------------------------------
* Begin.
*-----------------------------------------------------------------------
       CALL DSMATX (X,NX,PLIM,P,N,W,CAT,ROBUST,V,NLIM,IFAULT)
       IF (IFAULT .NE. 0) RETURN

       DO K = 1, N                               ! column sums
         T(K) = 0.
         DO J = 1, N
           T(K) = T(K) + V(J,K)
         END DO
       END DO

       TOT = 0.                                  ! overall total
       S = 1. / FLOAT(N)
       DO J = 1, N
         T(J) = T(J) * S
         TOT = TOT + T(J)
       END DO
       TOT = TOT * S

*-----------------------------------------------------------------------
*           center the matrix
*-----------------------------------------------------------------------
*                       1                1                 1
* -2 C[i,j] = A[i,j] - --- SUM A[i,k] - --- SUM A[k,j] + ----- SUM SUM A[k,l]
*                       N   k            N   k            N^2   k   l
*-----------------------------------------------------------------------
       DO K = 1, N
         DO J = K, N
           V(J,K) = -0.5 * (V(J,K) - T(J) - T(K) + TOT)
         END DO
       END DO

       CALL SEIGV (V, NLIM, N, E, V, T, IFAULT)  ! factor
       IF (IFAULT .NE. 0) RETURN

       DO K = 1, N/2                             ! decreasing order
         S = E(K)
         E(K) = E(N-K+1)
         E(N-K+1) = S
         DO J = 1, N
           S = V(J,K)
           V(J,K) = V(J,N-K+1)
           V(J,N-K+1) = S
         END DO
       END DO

       DO K = 1, N                               ! embed
         IF (E(K) > 0.) M = K
         S = SQRT(ABS(E(K)))
         DO J = 1, N
           V(J,K) = V(J,K) * S
         END DO
       END DO

       RETURN
      END  ! of mds


*-----------------------------------------------------------------------
* dsmatx - Construct a dissimilarity matrix
*
*___Name_________Type______In/Out____Description________________________
*   X(PLIM,N)    Real      In        Data of P dimensions by N objects
*   NX(PLIM,N)   Integer   In        1 indicates data is present for
*                                       the correspoding X value
*                                    0 indicates data is missing for
*                                       the corresponding X value
*   PLIM         Integer   In        Array bound of X and NX
*   P            Integer   In        Number of variables
*   N            Integer   In        Number of objects
*   W(P)         Real      In        Importance weights
*   CAT(P)       Integer   In        Number of categories each variable
*                                      0 indicates a continuous variable
*                                      1 indicates a dichotomous variable
*                                        values {0.0, 1.0} for absent, present
*                                     .GE. 2 indicates a categorical variable
*                                         values {1.0, 2.0, ...}
*                                     -1 indicates a cyclic variable,
*                                             period of one
*   ROBUST       Logical   In        Use Manhattan dissimilarity
*   D(NLIM,N)    Real      Out       Dissimilarity matrix
*   NLIM         Integer   In        Array bound of D
*   IFAULT       Integer   Out       Error code
*-----------------------------------------------------------------------
      SUBROUTINE DSMATX (X,NX,PLIM,P,N,W,CAT,ROBUST,D,NLIM,IFAULT)
       IMPLICIT NONE
       INTEGER PLIM, P, N, NLIM
       INTEGER NX(PLIM,N), CAT(P), IFAULT
       REAL X(PLIM,N), W(P), D(NLIM,N)
       LOGICAL ROBUST

       REAL CTOL, SMALL
       PARAMETER (CTOL = 1.E-3,   ! tolerance on category labels
     $            SMALL = 1.E-36)

       INTEGER I, J, K,           ! counters
     $         ND                 ! data present in D
       REAL DAV,                  ! average distance
     $      DIF,                  ! a difference
     $      TOT,                  ! a total
     $      WSUM                  ! sum of weights
       LOGICAL AZ                 ! all zero weights

       LOGICAL SAFDIV            ! external function

*-----------------------------------------------------------------------
* Begin.
*-----------------------------------------------------------------------
       IFAULT = 0
       IF (P < 1 .OR. N < 1) THEN     ! validate
         IFAULT = -2
         RETURN
       END IF

*             check for valid category labels
       DO I = 1, N
         DO J = 1, P
           IF (NX(J,I) > 0) THEN
             K = 0                               ! avoid compiler warning
             IF (CAT(J) > 0) THEN
               K = NINT(X(J,I))
               IF (ABS(X(J,I) - FLOAT(K)) > CTOL) THEN
                 IFAULT = -3
                 RETURN
               END IF
             END IF
             IF (CAT(J) .EQ. 1) THEN
               IF (.NOT. (K .EQ. 0 .OR. K .EQ. 1)) THEN
                 IFAULT = -4
                 RETURN
               ENDIF
             END IF
             IF (CAT(J) > 1) THEN
               IF (K < 1 .OR. K > CAT(J)) THEN
                 IFAULT = -5
                 RETURN
               END IF
             END IF
           END IF
         END DO
       END DO

*             require a nonzero weight
       AZ = .TRUE.
       DO J = 1, P
         IF (ABS(W(J)) > SMALL) AZ = .FALSE.
       END DO
       IF (AZ) THEN
         IFAULT = -6
         RETURN
       END IF

*             dissimilarity matrix
       ND = 0
       DAV = 0.
       DO K = 1, N-1
         DO J = K+1, N
           TOT = 0.
           WSUM = 0.
           DO 10 I = 1, P
             IF (NX(I,J) .LE. 0 .OR. NX(I,K) .LE. 0) GO TO 10
             IF (CAT(I) .EQ. 0) THEN
               DIF = X(I,J) - X(I,K)
             ELSE IF (CAT(I) > 1) THEN           ! categorical
               IF (ABS(X(I,J) - X(I,K)) > 0.5) THEN
                 DIF = 1.
               ELSE
                 DIF = 0.
               END IF
             ELSE IF (CAT(I) .EQ. 1) THEN        ! dichotomous
               IF (X(I,J) < 0.5 .AND. X(I,K) < 0.5) GO TO 10
               IF (X(I,J) < 0.5 .OR. X(I,K) < 0.5) THEN
                 DIF = 1.
               ELSE
                 DIF = 0.
               END IF
             ELSE IF (CAT(I) .EQ. -1) THEN       ! cyclic
               DIF = MOD(ABS(X(I,J) - X(I,K)), 1.)
               DIF = 2. * MIN(DIF, 1. - DIF)
             ELSE
               IFAULT = -7
               RETURN
             END IF

             WSUM = WSUM + W(I)
             IF (ROBUST) THEN
               TOT = TOT + W(I) * ABS(DIF)       ! Manhattan
             ELSE
               TOT = TOT + W(I) * DIF**2         ! squared Euclidean
             END IF
  10       CONTINUE  ! next I

           IF (SAFDIV(TOT, WSUM)) THEN           ! adjust for missing values
             ND = ND + 1
             D(J,K) = TOT / WSUM
             DAV = DAV + D(J,K)
           ELSE
             D(J,K) = -1.                        ! sentinel to be replaced
           END IF
         END DO
       END DO

       IF (ND .EQ. 0) THEN
         IFAULT = -8
         RETURN
       END IF
       DAV = DAV / FLOAT(ND)

       DO K = 1, N                               ! zero on the diagonal
         D(K,K) = 0.
         DO J = K+1, N
           IF (D(J,K) < 0.) D(J,K) = DAV         ! replace null distances
           D(K,J) = D(J,K)
         END DO
       END DO

       RETURN
      END  ! of dsmatx


*-----------------------------------------------------------------------
* pca - Principal Component Analysis.  No support for missing data.
*
*___Name_________Type______In/Out____Description________________________
*   X(PLIM,N)    Real      In        Data matrix
*   PLIM         Integer   In        Row bound of X
*   P            Integer   In        Variables each object
*   N            Integer   In        Number of objects
*   ROBUST       Logical   In        Use comedian instead of covariance
*   R(NLIM,M)    Real      Out       Projected data
*   NLIM         Integer   In        Row bound of R
*   M            Integer   In        # of dimensions requested to project
*   Q            Integer   Out       # of dimensions found
*   A(P)         Real      Out       Averages each variable
*   E(P)         Real      Out       Eigenvalues of covariance matrix
*   V(PLIM,P)    Real      Out       Eigenvectors of covariance matrix
*   T(P)         Real      Neither   Workspace
*   IFAULT       Integer   Out       Error code
*-----------------------------------------------------------------------
      SUBROUTINE PCA (X,PLIM,P,N,ROBUST,R,NLIM,M,Q,A,E,V,T,IFAULT)
       IMPLICIT NONE
       INTEGER PLIM, P, N, NLIM, M, Q, IFAULT
       REAL X(PLIM,N), R(NLIM,M), A(P), E(P), V(PLIM,P), T(P)
       LOGICAL ROBUST

       INTEGER I, J, K, K1        ! counters
       REAL TOT,                  ! a total
     $      XJ, XK                ! weighted X values

*-----------------------------------------------------------------------
*           Begin.  Validate.
*-----------------------------------------------------------------------
       IF (M < 1 .OR. M > P) THEN
         IFAULT = -1
         RETURN
       END IF
       IF (P < 1 .OR. N < 2) THEN  ! need two for variance
         IFAULT = -2
         RETURN
       END IF
       IFAULT = 0

*-----------------------------------------------------------------------
*           Compute averages.
*-----------------------------------------------------------------------
       DO J = 1, P
         IF (ROBUST) THEN
           DO I = 1, N                           ! fill work array
             R(I,1) = X(J,I)
           END DO
           CALL MEDIAN (R, N, A(J))
         ELSE
           TOT = 0.
           DO I = 1, N
             TOT = TOT + X(J,I)
           END DO
           A(J) = TOT / FLOAT(N)
         END IF
       END DO

*-----------------------------------------------------------------------
* Covariance / Comedian matrix.
*
* The comedian matrix is a highly robust estimate of covariance
* that is easy to compute.  Refer to:
*      Michael Falk
*      "On MAD and Comedians"
*      Ann. Inst. Statist. Math., v.49 n.4 pp.615-644, 1997
*
* The covariance matrix generalizes the variance.
* The comedian matrix generalizes the median absolute deviation.
*-----------------------------------------------------------------------
       DO K = 1, P
         DO J = K, P
           TOT = 0.
           DO I = 1, N
             XJ = X(J,I) - A(J)
             XK = X(K,I) - A(K)
             R(I,1) = XJ * XK
           END DO
           IF (ROBUST) THEN
             CALL MEDIAN (R, N, V(J,K))
           ELSE
             TOT = 0.
             DO I = 1, N
               TOT = TOT + R(I,1)
             END DO
             V(J,K) =  TOT / FLOAT(N - 1)
           END IF
         END DO
       END DO

*-----------------------------------------------------------------------
*         Factor covariance matrix
*-----------------------------------------------------------------------
       CALL SEIGV (V, PLIM, P, E, V, T, IFAULT)
       IF (IFAULT .NE. 0) RETURN                 ! failure

*-----------------------------------------------------------------------
* Step through columns of V in reverse order so that the first column
* of R has the greatest variance.
*             R[i,k] = X[j,i] V[j,p-k+1]
*-----------------------------------------------------------------------
       Q = 0
       DO K = 1, M
         K1 = P - K + 1
         IF (E(K1) .LE. 0.) RETURN
         Q = K
         DO I = 1, N
           R(I,K) = 0.
           DO J = 1, P
             R(I,K) = R(I,K) + (X(J,I) - A(J)) * V(J,K1)
           END DO
         END DO
       END DO
       RETURN
      END  ! of pca
