************************************************************************
*  valid.f - measures for clustering validation
*    by Andy Allinger, 2014, released to the public domain
*    This program may be used by any person for any purpose.
************************************************************************

************************************************************************
*  ARIND - adjusted Rand index
*   refer to:
*    "Comparing Partitions"
*    Lawrence Hubert and Phipps Arabie
*    Journal of Classification v.2 pp.193-218, 1985.
*
*                            I/O LIST
*___Name____________Type_________In/Out____Description________________
*   Y(N)            Integer      In        Correct class memberships
*   Z(N)            Integer      In        Computed cluster labels
*   N               Integer      In        Length of Y and Z
*   G               Integer      In        Number of correct classes
*   K               Integer      In        Number of discovered clusters
*   WORK(G*K+G+K)   Real         Neither   Workspace
*   ACC             Real         Out       Accuracy
************************************************************************
      SUBROUTINE ARIND (Y, Z, N, G, K, WORK, ACC)
       IMPLICIT NONE
       INTEGER Y, Z, N, G, K
       REAL WORK, ACC
       DIMENSION Y(N), Z(N), WORK(G*K+G+K)

*         partition workspace
       CALL ARI01 (Y, Z, N, G, K,
     $              WORK(1), WORK(G*K+1), WORK(G*K+G+1), ACC)

       RETURN
      END  ! of ARIND


************************************************************************
      SUBROUTINE ARI01 (Y, Z, N, G, K, M, A, B, ACC)
       IMPLICIT NONE
       INTEGER Y, Z, N, G, K
       REAL M, A, B, ACC
       DIMENSION Y(N), Z(N), M(G,K), A(G), B(K)

*            local variables
       INTEGER I, J, L              ! counters
       REAL RI, AS, BS, ADJ

*         validate
       ACC = 0.
       DO I = 1, N
         IF (Y(I) < 1 .OR. Y(I) > G .OR. Z(I) < 1 .OR. Z(I) > K) RETURN
       END DO

       DO J = 1, K                ! clear
         DO I = 1, G
           M(I,J) = 0.
         END DO
       END DO
       DO I = 1, N                ! count contingency table
         J = Y(I)
         L = Z(I)
         M(J,L) = M(J,L) + 1.
       END DO

*          row sums
       DO I = 1, G
         A(I) = 0.
         DO J = 1, K
           A(I) = A(I) + M(I,J)
         END DO
       END DO

*          column sums
       DO J = 1, K
         B(J) = 0.
         DO I = 1, G
           B(J) = B(J) + M(I,J)
         END DO
       END DO

*               main term
       RI = 0.
       DO J = 1, K
         DO I = 1, G
           L = NINT(M(I,J))
           RI = RI + L * (L-1) / 2
         END DO
       END DO

*               adjust term
       AS = 0.
       DO I = 1, G
         L = NINT(A(I))
         AS = AS + L * (L-1) / 2
       END DO

       BS = 0.
       DO J = 1, K
         L = NINT(B(J))
         BS = BS + L * (L-1) / 2
       END DO

       ADJ = AS * BS / ( N * (N-1) / 2)

*               finish
       ACC = (RI - ADJ) / (0.5 * (AS + BS) - ADJ)
       RETURN
      END  ! of ARI01


**********************************************************************
*  VINFO - variation of information
*   refer to:
*    "Comaparing Clusterings"
*    Marina Meila
*    University of Washington Technical report, 2002
*
*                            I/O LIST
*___Name____________Type_________In/Out____Description________________
*   Y(N)            Integer      In        Correct class memberships
*   Z(N)            Integer      In        Computed cluster labels
*   N               Integer      In        Length of Y and Z
*   G               Integer      In        Number of correct classes
*   K               Integer      In        Number of discovered clusters
*   WORK(G*K+G+K)   Real         Neither   Workspace
*   ACC             Real         Out       Accuracy
***********************************************************************
      SUBROUTINE VINFO (Y, Z, N, G, K, WORK, ACC)
       IMPLICIT NONE
       INTEGER Y, Z, N, G, K
       REAL WORK, ACC
       DIMENSION Y(N), Z(N), WORK(G*K+G+K)

*         partition workspace
       CALL VI01 (Y, Z, N, G, K,
     $             WORK(1), WORK(G*K+1), WORK(G*K+G+1), ACC)
       RETURN
      END  ! of VINFO


***********************************************************************
      SUBROUTINE VI01 (Y, Z, N, G, K, M, A, B, ACC)
       IMPLICIT NONE
       INTEGER Y, Z, N, G, K
       REAL M, A, B, ACC
       DIMENSION Y(N), Z(N), M(G,K), A(G), B(K)

*            local variables
       INTEGER I, J, L              ! counters
       REAL HA, HB, IMUT, SMALL  ! entropies, mutual information, small #
       SAVE SMALL
       DATA SMALL / 1E-9 /

*         validate
       ACC = 0.
       DO I = 1, N
         IF (Y(I) < 1 .OR. Y(I) > G .OR. Z(I) < 1 .OR. Z(I) > K) RETURN
       END DO

       DO J = 1, K                        ! clear
         DO I = 1, G
           M(I,J) = 0.
         END DO
       END DO
       DO I = 1, N                        ! count
         J = Y(I)
         L = Z(I)
         M(J,L) = M(J,L) + 1.
       END DO

*          row sums
       DO I = 1, G
         A(I) = 0.
         DO J = 1, K
           A(I) = A(I) + M(I,J)
         END DO
         A(I) = A(I) / FLOAT(N)    ! make a probability
       END DO

*          column sums
       DO J = 1, K
         B(J) = 0.
         DO I = 1, G
           B(J) = B(J) + M(I,J)
         END DO
         B(J) = B(J) / FLOAT(N)     ! make a probability
       END DO

*         change counts to a joint probability
       DO J = 1, K
         DO I = 1, G
           M(I,J) = M(I,J) / FLOAT(N)
         END DO
       END DO

*         entropy A
       HA = 0.
       DO 10 I = 1, G
         IF (A(I) < SMALL) GO TO 10
         HA = HA - A(I) * LOG(A(I))
  10   CONTINUE

*         entropy B
       HB = 0.
       DO 20 J = 1, K
         IF (B(J) < SMALL) GO TO 20
         HB = HB - B(J) * LOG(B(J))
  20   CONTINUE

*         mutual information
       IMUT = 0.
       DO J = 1, K
         DO 30 I = 1, G
           IF (M(I,J) < SMALL) GO TO 30
           IMUT = IMUT + M(I,J) * LOG( M(I,J) / (A(I) * B(J)) )
  30     CONTINUE
       END DO

*             variation of information
       ACC = HA + HB - 2. * IMUT
       ACC = ACC / LOG(2.)       ! convert to bits
       RETURN
      END  ! of VI01
*+++++++++++++++++++++++++ End of file valid.f +++++++++++++++++++++++++
