*-----------------------------------------------------------------------
* eigen.f - Solve real symmetric eigenproblems.  An extract of EISPACK,
* as distributed in the National Bureau of Standards Core Mathematical
* Library.  This file is public domain.
*-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C   TRED2
C   DATE WRITTEN   760101   (YYMMDD)
C   REVISION DATE  830518   (YYMMDD)
C   CATEGORY NO.  D4C1B1
C   AUTHOR  Smith, B. T., et al.
C   PURPOSE  Reduce real symmetric matrix to symmetric tridiagonal
C            matrix using and accumulating orthogonal transformation
C   DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure TRED2,
C     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
C     Handbook for Auto. Comp., Vol.II-Linear Algebra, 212-226(1971).
C
C     This subroutine reduces a REAL SYMMETRIC matrix to a
C     symmetric tridiagonal matrix using and accumulating
C     orthogonal similarity transformations.
C
C     On Input
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrix.
C
C        A contains the real symmetric input matrix.  Only the
C          lower triangle of the matrix need be supplied.
C
C     On Output
C
C        D contains the diagonal elements of the tridiagonal matrix.
C
C        E contains the subdiagonal elements of the tridiagonal
C          matrix in its last N-1 positions.  E(1) is set to zero.
C
C        Z contains the orthogonal transformation matrix
C          produced in the reduction.
C
C        A and Z may coincide.  If distinct, A is unaltered.
C
C     Questions and comments should be directed to B. S. Garbow,
C     Applied Mathematics Division, Argonne National Laboratory
C     ------------------------------------------------------------------
C   REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema, C. B. Moler, *Matrix Eigen-
C                 system Routines - EISPACK Guide*, Springer-Verlag,
C                 1976.
C   ROUTINES CALLED  (NONE)
C
C-----------------------------------------------------------------------
      SUBROUTINE TRED2(NM,N,A,D,E,Z)
      IMPLICIT NONE
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL A(NM,N),D(N),E(N),Z(NM,N)
      REAL F,G,H,HH,SCALE

C   First executable statement
      DO 100 I = 1, N

         DO 100 J = 1, I
            Z(I,J) = A(I,J)
  100 CONTINUE

      IF (N .EQ. 1) GO TO 320
C     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0
         SCALE = 0.0
         IF (L .LT. 2) GO TO 130
C     .......... Scale row (ALGOL TOL then not needed) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(Z(I,K))

         IF (SCALE .NE. 0.0) GO TO 140
  130    E(I) = Z(I,L)
         GO TO 290

  140    DO 150 K = 1, L
            Z(I,K) = Z(I,K) / SCALE
            H = H + Z(I,K) * Z(I,K)
  150    CONTINUE

         F = Z(I,L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         Z(I,L) = F - G
         F = 0.0

         DO 240 J = 1, L
            Z(J,I) = Z(I,J) / H
            G = 0.0
C     .......... Form element of A*U ..........
            DO 180 K = 1, J
  180       G = G + Z(J,K) * Z(I,K)

            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220

            DO 200 K = JP1, L
  200       G = G + Z(K,J) * Z(I,K)
C     .......... Form element of P ..........
  220       E(J) = G / H
            F = F + E(J) * Z(I,J)
  240    CONTINUE

         HH = F / (H + H)
C     .......... Form reduced A ..........
         DO 260 J = 1, L
            F = Z(I,J)
            G = E(J) - HH * F
            E(J) = G

            DO 260 K = 1, J
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)
  260    CONTINUE

  290    D(I) = H
  300 CONTINUE

  320 D(1) = 0.0
      E(1) = 0.0
C     .......... Accumulation of transformation matrices ..........
      DO 500 I = 1, N
         L = I - 1
         IF (D(I) .EQ. 0.0) GO TO 380

         DO 360 J = 1, L
            G = 0.0

            DO 340 K = 1, L
  340       G = G + Z(I,K) * Z(K,J)

            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * Z(K,I)
  360    CONTINUE

  380    D(I) = Z(I,I)
         Z(I,I) = 1.0
         IF (L .LT. 1) GO TO 500

         DO 400 J = 1, L
            Z(I,J) = 0.0
            Z(J,I) = 0.0
  400    CONTINUE

  500 CONTINUE

      RETURN
      END  ! of TRED2


C-----------------------------------------------------------------------
C PYTHAG - Finds sqrt(A**2+B**2) without overflow or destructive underflow
C-----------------------------------------------------------------------
      REAL FUNCTION PYTHAG(A,B)
      IMPLICIT NONE
      REAL A,B

      REAL P,Q,R,S,T
C   First executable statement
      P = AMAX1(ABS(A),ABS(B))
      Q = AMIN1(ABS(A),ABS(B))
      IF (Q .EQ. 0.0) GO TO 20
   10 CONTINUE
         R = (Q/P)**2
         T = 4.0 + R
         IF (T .EQ. 4.0) GO TO 20
         S = R/T
         P = P + 2.0*P*S
         Q = Q*S
      GO TO 10
   20 PYTHAG = P
      RETURN
      END  ! of PYTHAG


C-----------------------------------------------------------------------
C   TQL2
C   DATE WRITTEN   760101   (YYMMDD)
C   REVISION DATE  830518   (YYMMDD)
C   CATEGORY NO.  D4A5,D4C2A
C   AUTHOR  Smith, B. T., et al.
C   PURPOSE  Compute eigenvalues and eigenvectors of symmetric
C            tridiagonal matrix.
C   DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure TQL2,
C     Num. Math. 11, 293-306(1968) by Bowdler, Martin, Reinsch, and
C     Wilkinson.
C     Handbook for Auto. Comp., Vol.II-Linear Algebra, 227-240(1971).
C
C     This subroutine finds the eigenvalues and eigenvectors
C     of a SYMMETRIC TRIDIAGONAL matrix by the QL method.
C     The eigenvectors of a FULL SYMMETRIC matrix can also
C     be found if  TRED2  has been used to reduce this
C     full matrix to tridiagonal form.
C
C     On Input
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrix.
C
C        D contains the diagonal elements of the input matrix.
C
C        E contains the subdiagonal elements of the input matrix
C          in its last N-1 positions.  E(1) is arbitrary.
C
C        Z contains the transformation matrix produced in the
C          reduction by  TRED2, if performed.  If the eigenvectors
C          of the tridiagonal matrix are desired, Z must contain
C          the identity matrix.
C
C      On Output
C
C        D contains the eigenvalues in ascending order.  If an
C          error exit is made, the eigenvalues are correct but
C          unordered for indices 1,2,...,IERR-1.
C
C        E has been destroyed.
C
C        Z contains orthonormal eigenvectors of the symmetric
C          tridiagonal (or full) matrix.  If an error exit is made,
C          Z contains the eigenvectors associated with the stored
C          eigenvalues.
C
C        IERR is set to
C          Zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after 30 iterations.
C
C     Calls PYTHAG(A,B) for SQRT(A**2 + B**2).
C
C     Questions and comments should be directed to B. S. Garbow,
C     Applied Mathematics Division, Argonne National Laboratory
C     ------------------------------------------------------------------
C   REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema, C. B. Moler, *Matrix Eigen-
C                 system Routines - EISPACK Guide*, Springer-Verlag,
C                 1976.
C   ROUTINES CALLED  PYTHAG
C
C-----------------------------------------------------------------------
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
      IMPLICIT NONE
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      REAL D(N),E(N),Z(NM,N)
      REAL B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2
      REAL PYTHAG

C   First executable statement
      IERR = 0
      IF (N .EQ. 1) GO TO 1001

      DO 100 I = 2, N
  100 E(I-1) = E(I)

      F = 0.0
      B = 0.0
      E(N) = 0.0
      S2 = 0.0
      C3 = 0.0

      DO 240 L = 1, N
         J = 0
         H = ABS(D(L)) + ABS(E(L))
         IF (B .LT. H) B = H
C     .......... Look for small sub-diagonal element ..........
         DO 110 M = L, N
            IF (B + ABS(E(M)) .EQ. B) GO TO 120
C     .......... E(N) is always zero, so there is no exit
C                through the bottom of the loop ..........
  110    CONTINUE

*            added 2011/06/13 AJA
         IF (M > N) THEN
           IERR = -1
           M = N
         END IF

  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... Form shift ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0 * E(L))
         R = PYTHAG(P,1.0)
         D(L) = E(L) / (P + SIGN(R,P))
         D(L1) = E(L) * (P + SIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145

         DO 140 I = L2, N
  140    D(I) = D(I) - H

  145    F = F + H
C     .......... QL transformation ..........
         P = D(M)
         C = 1.0
         C2 = C
         EL1 = E(L1)
         S = 0.0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150
            C = E(I) / P
            R = SQRT(C*C+1.0)
            E(I+1) = S * P * R
            S = C / R
            C = 1.0 / R
            GO TO 160
  150       C = P / E(I)
            R = SQRT(C*C+1.0)
            E(I+1) = S * E(I) * R
            S = 1.0 / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
C     .......... Form vector ..........
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE

  200    CONTINUE

         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         IF (B + ABS(E(L)) .GT. B) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
C     .......... Order eigenvalues and eigenvectors ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)

         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE

         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P

         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE

  300 CONTINUE

      GO TO 1001
C     .......... Set error -- No convergence to an
C                eigenvalue after 30 iterations ..........
 1000 IERR = L
 1001 RETURN
      END  ! of TQL2


C-----------------------------------------------------------------------
C        Eigenvalues and eigenvectors of symmetric real matrices
C-----------------------------------------------------------------------
*
*___Name______Type______I/O_______Description___________________________
*   A(KA,N)   Real      In        The matrix
*   KA        Integer   In        Row dimension of A
*   N         Integer   In        Order of A
*   W(N)      Real      Out       Eigenvalues
*   Z(KA,N)   Real      Out       Eigenvectors
*   T(N)      Real      Neither   Workspace
*   IERR      Integer   Out       Error code
*-----------------------------------------------------------------------
      SUBROUTINE SEIGV (A, KA, N, W, Z, T, IERR)
      IMPLICIT NONE
      INTEGER KA, N, IERR
      REAL A(KA,N), W(N), Z(KA,N), T(N)

      CALL TRED2 (KA, N, A, W, T, Z)
      CALL TQL2 (KA, N, W, T, Z, IERR)
      RETURN
      END  ! of SEIGV
*++++++++++++++++++ End of file eigen.f ++++++++++++++++++++++++++++++++
