************************************************************************
* neugas.f - neural gas, a data clustering algorithm
* by Andy Allinger, 2013-2017, released to the public domain
* This program may be used by any person for any purpose.
*
* Refer to:
*      Thomas Martinetz and Klaus Schulten
*      A "Neural-Gas" Network Learns Topologies
*      in Artificial Neural Networks, Ed. T. Kohonen, K. Makisara,
*      O. Simula, and J. Kangas.  Holland: Elsevier, 1991.
*
* Also the description by Bernd Fritzke at:
*      www.demogng.de/JavaPaper/node16.html
*
* This routine follows Fritzke and omits the topology.
*
*___Name________Type______In/Out____Description_______________________
*   X(PLIM,N)   Real      In        Data Matrix
*   PLIM        Integer   In        Row dimension of X
*   P           Integer   In        Measurements per object
*   N           Integer   In        Number of objects
*   C(PLIM,K)   Real      Both      Cluster centers
*   K           Integer   In        Number of clusters
*   Z(N)        Integer   Out       The cluster each point belongs to
*   IND(K)      Integer   Neither   Index array
*   D(K)        Real      Neither   Distances to each center
*   U           Real      Out       Total squared distance of each point
*                                     to its nearest center
************************************************************************
      SUBROUTINE NEUGAS (X, PLIM, P, N, C, K, Z, IND, D, U)
       IMPLICIT NONE
       INTEGER PLIM, N, K
       INTEGER P, Z(N), IND(K)
       REAL X(PLIM,N), C(PLIM,K), D(K), U

*             constants
       INTEGER TAU
       REAL BIG, LAMI, LAMF, EPSI, EPSF

*             revised parameters
        PARAMETER (BIG = 1.E30,   ! arbitrary large number
     $             LAMI = 10.,    ! range of adjustment, initial
     $             LAMF = 0.01,   ! range of adjustment, final
     $             EPSI = 0.3,    ! adjustment, initial
     $             EPSF = 0.005,  ! adjustment, final
     $             TAU = 40 000)  ! iterations

*             local variables
       INTEGER H, I, J, L         ! counters
       INTEGER HAT                ! external function

       REAL
     $      BEST,                 ! a distance
     $      LFACT, EFACT,         ! update factors for EPS, LAMBDA
     $      EPS,                  ! adjustment factor
     $      LAMBDA,               ! range factor
     $      TOT,                  ! a total
     $      Y                     ! factor in update formula

************************************************************************
*             begin
************************************************************************
       LAMBDA = LAMI
       EPS = EPSI
       LFACT = (LAMF / LAMI) ** (1. / FLOAT(TAU - 1))
       EFACT = (EPSF / EPSI) ** (1. / FLOAT(TAU - 1))

************************************************************************
*             main loop
************************************************************************
       DO H = 1, TAU
         I = HAT(N)               ! choose a datum at random

*             find distance to each center
         DO L = 1, K
           TOT = 0.
           DO J = 1, P
             TOT = TOT + (X(J,I) - C(J,L)) **2
           END DO
           D(L) = TOT
         END DO
         CALL QSORTR (D, IND, K)            ! sort distances
         DO L = 1, K
           D(IND(L)) = FLOAT(L - 1)             ! rank
         END DO

*             adapt the cluster centers
         DO L = 1, K
           Y = EPS * EXP(-D(L) / LAMBDA)
           DO J = 1, P
             C(J,L) = C(J,L) + Y * (X(J,I) - C(J,L))
           END DO
         END DO

*             adjust parameters
         LAMBDA = LAMBDA * LFACT
         EPS = EPS * EFACT
       END DO  ! next iteration

*             assignments and residual
       U = 0.
       DO I = 1, N
         BEST = BIG
         DO L = 1, K
           TOT = 0.
           DO J = 1, P
             TOT = TOT + (X(J,I) - C(J,L))**2
           END DO
           IF (TOT < BEST) THEN
             BEST = TOT
             Z(I) = L
           END IF
         END DO
         U = U + BEST
       END DO

       RETURN
      END  ! of neugas
