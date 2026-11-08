*-----------------------------------------------------------------------
* prep.f - Routines for preprocessing data
* by Andy Allinger, 2011-2017, released to the public domain
* This program may be used by any person for any purpose.
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*       substitute for LEN_TRIM for old compilers
*-----------------------------------------------------------------------
      FUNCTION LTRIM (STRING)
       IMPLICIT NONE
       CHARACTER*(*) STRING
       INTEGER LTRIM
       INTEGER J, L

       L = LEN(STRING)
       DO 50 J = L, 1, -1
         IF (STRING(J:J) .EQ. ' ') GO TO 50
         LTRIM = J
         RETURN
  50   CONTINUE
       LTRIM = 0
       RETURN
      END ! of ltrim


*-----------------------------------------------------------------------
* Sort non-missing data to the front, by a single iteration of quick sort.
*___Name_____Type______In/Out___Description_____________________________
*   NX(N)    Integer   In       Input array, 0 = missing, 1 = present
*   IND(N)   Integer   Out      Permutation vector
*   N        Integer   In       Length of array
*-----------------------------------------------------------------------
      SUBROUTINE QSORTM (NX, IND, N)
       IMPLICIT NONE
       INTEGER NX, IND, N
       DIMENSION NX(N), IND(N)

*        Local variables
       INTEGER I, J, K

*         Begin.
       DO I = 1, N                ! initialize
         IND(I) = I
       END DO
       IF (N < 2) RETURN          ! validate
       I = 1
       J = N

*           Find a missing datum at the front of the array.
  50   IF (I < J) THEN
         IF (NX(IND(I)) > 0) THEN
           I = I + 1
           GO TO 50
         END IF
       END IF

*           Find a non-missing datum at the rear of the array
  60   IF (I < J) THEN
         IF (NX(IND(J)) .LE. 0) THEN
           J = J - 1
           GO TO 60
         END IF
       END IF

*          Swap these elements
       IF (I < J) THEN
         K = IND(I)
         IND(I) = IND(J)
         IND(J) = K
         I = I + 1
         J = J - 1
         GO TO 50
       END IF
       RETURN
      END   ! of qsortm


*-----------------------------------------------------------------------
* jscode - succint code:  change arbitrary integers to consecutive integers
*
*___Name______Type______In/Out____Description___________________________
*   JNAM(N)   Integer   Both      Data array
*   NX(N)     Integer   In        1 if data exists, each datum
*   IND(N)    Integer   Neither   Permutation vector
*   JND(N)    Integer   Neither   Permutation vector
*   N         Integer   In        Length of JNAM
*-----------------------------------------------------------------------
      SUBROUTINE JSCODE (JNAM, NX, IND, JND, N)
       IMPLICIT NONE
       INTEGER JNAM, NX, IND, JND, N
       DIMENSION JNAM(N), NX(N), IND(N), JND(N)

*           local variables
       INTEGER I, J, K, L, M

*            Begin. Move missing data to end.
       CALL QSORTM (NX, JND, N)
       M = 0                      ! count data available
       DO I = 1, N
         IF (NX(I) > 0) M = M + 1
       END DO
       IF (M .EQ. 0) RETURN
       IF (M .EQ. 1) THEN
         JNAM(JND(1)) = 1
         RETURN
       END IF
       CALL IORDER (JNAM, JND, N)

*              Sort implicitly.
       CALL QSORTI (JNAM, IND, M)

*             Compare adjacent entries.
       K = JNAM(IND(1))
       L = 1
       JNAM(IND(1)) = 1
       DO I = 2, M
         J = JNAM(IND(I))
         IF (J .NE. K) THEN
           K = J
           L = L + 1
         END IF
         JNAM(IND(I)) = L
       END DO

*            Restore original order.
       DO I = 1, N
         IND(JND(I)) = I
       END DO
       CALL IORDER (JNAM, IND, N)

       RETURN
      END ! of jscode


*-----------------------------------------------------------------------
* ireid - replace text names with integers 1, 2, 3... for each distinct item
*
*___Name_____Type____________I/O_______Description______________________
*   H(N)     Character*(*)   In        Items or names
*   ID(N)    Integer         Out       Successive small integers
*   NX(N)    Integer         In        1 if datum exists
*   IND(N)   Integer         Neither   Permutation vector
*   JND(N)   Integer         Neither   Permutation vector
*   N        Integer         In        Length of array
*   HWORK    Character*(*)   Neither   Used by HPSORT, same length as H
*   IERR     Integer         Out       Error code
*-----------------------------------------------------------------------
      SUBROUTINE IREID (H, ID, NX, IND, JND, N, HWORK, IERR)     ! re-id
       IMPLICIT NONE
       INTEGER ID, NX, IND, JND, N, IERR
       CHARACTER*(*) H, HWORK
       DIMENSION H(N), ID(N), NX(N), IND(N), JND(N)

*           local variables
       INTEGER I, J, L, M

*            Begin.
       IERR = 0
       CALL QSORTM (NX, JND, N)         ! move missing data to end
       M = 0                            ! count data available
       DO I = 1, N
         IF (NX(I) > 0) M = M + 1
       END DO
       IF (M .EQ. 0) RETURN
       IF (M .EQ. 1) THEN
         ID(IND(1)) = 1
         RETURN
       END IF
       CALL HPPERM (H, N, JND, HWORK, IERR)       ! apply permutation
       IF (IERR .NE. 0) RETURN

*            Sort to bring duplicate values adjacent.
       L = LEN(H(1))
       CALL HPSORT(H, M, 1, L, IND, 2, HWORK, IERR)
       IF (IERR .NE. 0) RETURN
       CALL IORDER (JND, IND, M)       ! apply permutation

*             compare adjacent strings
       J = 1
       L = 1
       ID(J) = L
       DO I = 2, M
         IF (H(I) .NE. H(J)) THEN
           J = I
           L = L + 1
         END IF
         ID(I) = L
       END DO

*               restore original order
       DO I = 1, N
         IND(JND(I)) = I
       END DO
       CALL HPPERM (H, N, IND, HWORK, IERR)
       CALL IORDER (ID, IND, N)
       IF (IERR .NE. 0) RETURN

* Change labels to be the order of first occurance in the input array
       DO I = 1, N
         IND(I) = 0
       END DO
       L = 0
       DO I = 1, N                ! form replacement table in IND
         IF (NX(I) > 0) THEN
           IF (IND(ID(I)) .EQ. 0) THEN
             L = L + 1
             IND(ID(I)) = L
           END IF
         END IF
       END DO
       DO I = 1, N                  ! replace
         IF (NX(I) > 0) ID(I) = IND(ID(I))
       END DO
       RETURN
      END  ! of ireid


*-----------------------------------------------------------------------
* ordnal - replace data with its ordinal rank, rescaled to 0...1
*
*___Name_____Type______I/O________Description___________________________
*   X(N)     Real      Both       Data vector
*   NX(N)    Integer   In         1 if data exists
*   IND(N)   Integer   Neither    Permutation vector
*   JND(N)   Integer   Neither    Permutation vector
*   N        Integer   In         Length of vector
*-----------------------------------------------------------------------
      SUBROUTINE ORDNAL (X, NX, IND, JND, N)
       IMPLICIT NONE
       INTEGER NX, IND, JND, N
       REAL X
       DIMENSION X(N), NX(N), IND(N), JND(N)

*           local variables
       INTEGER I, J, K,        ! indices
     $         M               ! count of data present

       REAL              V                   ! new value
       LOGICAL ALMEQ             ! almost-equal function

*             begin
       CALL QSORTM (NX, JND, N)         ! move missing data to end
       M = 0                            ! count data available
       DO I = 1, N
         IF (NX(I) > 0) M = M + 1
       END DO
       IF (M .EQ. 0) RETURN
       IF (M .EQ. 1) THEN
         X(IND(1)) = 0.
         RETURN
       END IF
       CALL RORDER (X, JND, N)       ! apply permutation

       CALL QSORTR (X, IND, M)         ! bring duplicate values adjacent
       CALL RORDER (X, IND, M)         ! apply permutation
       CALL IORDER (JND, IND, M)       ! apply permutation

*        compare adjacent values
       J = 1
       DO 10 I = 2, M
         IF (ALMEQ(X(I), X(J))) GO TO 10
         V = 0.5 * FLOAT(J + I-1) - 1.
         DO K = J, I-1          ! replace
           X(K) = V
         END DO
         J = I
  10   CONTINUE

*             last stretch
       V = 0.5 * FLOAT(J + M) - 1.
       DO K = J, M
         X(K) = V
       END DO

*               rescale 0...1
       V = 1. / FLOAT(M - 1)
       DO I = 1, M
         X(I) = X(I) * V
       END DO

*               restore original order
       DO I = 1, N
         IND(JND(I)) = I
       END DO
       CALL RORDER (X, IND, N)
       RETURN
      END ! of ordnal


*-----------------------------------------------------------------------
* Construct vertices of a regular simplex according to Manhattan metric
*
*___Name_______Type______I/O___Description______________________________
*   M          Integer   In    Number of vertices
*   S(M-1,M)   Real      Out   Coordinates
*-----------------------------------------------------------------------
      SUBROUTINE RSCL1 (M, S)     ! regular simplex coordinates, Lebesgue 1
       IMPLICIT NONE
       INTEGER M
       REAL S
       DIMENSION S(M-1,M)

*        local variables
       INTEGER I, J, K

*        begin
       IF (M < 2) RETURN
       S(1,1) = 0.
       S(1,2) = 1.

       DO K = 3, M            ! each vertex
         S(1,K) = 0.5          ! first entry always 1/2

         DO I = 2, K-2          ! zeros
           S(I,K) = 0.
         END DO

         DO J = 1, K-1        ! fill in zeros
           S(K-1,J) = 0.
         END DO

         S(K-1,K) = 0.5         ! last entry always 1/2
       END DO  ! next k
       RETURN
      END  ! of rscl1


*-----------------------------------------------------------------------
* Construct vertices of a regular simplex according to Euclidean metric
*
*___Name_______Type______I/O___Description______________________________
*   M          Integer   In    Number of vertices
*   S(M-1,M)   Real      Out   Coordinates
*-----------------------------------------------------------------------
      SUBROUTINE RSCL2 (M, S)     ! regular simplex coordinates, Lebesgue 2
       IMPLICIT NONE
       INTEGER M
       REAL S
       DIMENSION S(M-1,M)

*        local variables
       INTEGER I, J, K
       REAL T

*        begin
       IF (M < 2) RETURN
       S(1,1) = 0.
       S(1,2) = 1.

       DO K = 3, M            ! each vertex
         DO I = 1, K-2         ! average each row
           T = 0.
           DO J = 1, K-1
             T = T + S(I,J)
           END DO
           S(I,K) = T / FLOAT(K-1)
         END DO  ! next I

         DO J = 1, K-1        ! fill in zeros
           S(K-1,J) = 0.
         END DO

         T = 0.
         DO I = 1, K-2           ! make unit length
           T = T + S(I,K) **2
         END DO
         S(K-1,K) = SQRT(1. - T)
       END DO  ! next k
       RETURN
      END  ! of rscl2


*-----------------------------------------------------------------------
* rscode - Recode categorical variables as column of regular simplex
* Also scale down importance weights to preserve original importances.
*
*___Name_________Type______I/O_______Description________________________
*   X(PLIM,N)    Real      Both      Data array
*   NX(PLIM,N)   Integer   Both      1 if data is present for this measurement
*   PLIM         Integer   In        Row bound of X and EX
*   P            Integer   In        Number of coded variables
*   N            Integer   In        Number of objects
*   W(P)         Real      Both      Importance weights
*   CAT(O)       Integer   In        Number of categories each variable
*   O            Integer   In        Number of input variables
*   S(CATMAX-1,CATMAX)
*                Real      Neither   Simplex coordinates
*   CATMAX       Integer   In        Maximum number of categories
*   ROBUST       Logical   In        Use L1 statistics?
*   IFAULT       Integer   Out       Error code
*-----------------------------------------------------------------------
      SUBROUTINE RSCODE (X,NX,PLIM,P,N,W,CAT,O,S,CATMAX,ROBUST,IFAULT)
       IMPLICIT NONE
       INTEGER NX, PLIM, P, N, CAT, O, CATMAX, IFAULT
       REAL X, W, S
       LOGICAL ROBUST
       DIMENSION X(PLIM,N), NX(PLIM,N), W(P), CAT(O), S(CATMAX-1,CATMAX)

       REAL XBOGUS                   ! constant
       PARAMETER (XBOGUS = -999.)

       INTEGER CLM1, H, I, J, K, L       ! local variables
       REAL FCLM1

*-----------------------------------------------------------------------
*          Begin.  Validate
*-----------------------------------------------------------------------
       IFAULT = 0
       J = 0
       K = 0
       DO L = 1, O
         IF (CAT(L) .EQ. 1) THEN
           IFAULT = 1
           RETURN
         ELSE IF (CAT(L) < 0) THEN
           IFAULT = -1
           RETURN
         END IF
         J = J + 1 + DIM(CAT(L), 2)
         K = MAX(K, CAT(L))
       END DO
       IF (J .NE. P) THEN
         IFAULT = 2
         RETURN
       END IF
       IF (K > CATMAX) THEN
         IFAULT = 3
         RETURN
       END IF

*-----------------------------------------------------------------------
*        change weights
*-----------------------------------------------------------------------
       J = P + 1
       DO L = O, 1, -1
         IF (CAT(L) .EQ. 0) THEN
           J = J - 1
           W(J) = W(L)
         ELSE
           CLM1 = CAT(L) - 1
           FCLM1 = FLOAT(CLM1)
           DO H = CLM1, 1, -1
             J = J - 1
             W(J) = W(L) / FCLM1
           END DO
         END IF
       END DO

*-----------------------------------------------------------------------
*              Recode.
*-----------------------------------------------------------------------
       IF (ROBUST) THEN                  ! construct simplex
         CALL RSCL1 (CATMAX, S)
       ELSE
         CALL RSCL2 (CATMAX, S)
       END IF

       DO I = N, 1, -1
         J = P + 1
         DO L = O, 1, -1
           IF (CAT(L) .EQ. 0) THEN
             J = J - 1
             IF (NX(L,I) > 0) THEN
               X(J,I) = X(L,I)
               NX(J,I) = 1
             ELSE
               X(J,I) = XBOGUS
               NX(J,I) = 0
             END IF
           ELSE
             IF (NX(L,I) > 0) THEN
               K = NINT(X(L,I))
               IF (K < 1 .OR. K > CAT(L)) THEN
                 IFAULT = 4
                 RETURN
               END IF
             END IF
             CLM1 = CAT(L) - 1
             FCLM1 = FLOAT(CLM1)
             IF (.NOT. ROBUST) FCLM1 = SQRT(FCLM1)
             DO H = CLM1, 1, -1
               J = J - 1
               IF (NX(L,I) > 0) THEN
                 X(J,I) = S(H,K) * FCLM1
                 NX(J,I) = 1
               ELSE
                 X(J,I) = XBOGUS
                 NX(J,I) = 0
               END IF
             END DO
           END IF
         END DO
       END DO

       RETURN
      END  ! of rscode
