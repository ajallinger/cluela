      PROGRAM EMBED  ! embed cyclic data
       IMPLICIT NONE
       INTEGER P, N, M, J
       PARAMETER (P = 1,
     $            N = 60)
       INTEGER NX(P,N), CAT(P), IFAULT
       REAL X(P,N), W(P), V(N,N), E(N), T(N)
       LOGICAL ROBUST

       W(1) = 1.0
       ROBUST = .FALSE.
       CAT(1) = -1
       DO J = 1, N
         X(P,J) = J / FLOAT(N)
         NX(P,J) = 1
       END DO

       CALL MDS (X, NX, P, P, N, W, CAT, ROBUST,
     $            V, N, M, E, T, IFAULT)
       PRINT *, 'mds returns #', ifault, ' with ', M, ' real directions'

       PRINT *, 'coordinates:'
       DO J = 1, N
         PRINT *, V(J,1), V(J,2)
       END DO
      END

