*-----------------------------------------------------------------------
* sort.f - Routines for sorting and ordering, from the SLATEC and
* Naval Surface Warfare Center libraries.  This file is public domain.
*-----------------------------------------------------------------------

C***********************************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C   This subroutine uses an order N*log(N) quick sort to sort an integer
C array X into increasing order.  The algorithm is as follows.  IND is
C initialized to the ordered sequence of indices 1,...,N, and all
C interchanges are applied to IND.  X is divided into two portions by
C picking a central element T.  The first and last elements are compared
C with T, and interchanges are applied as necessary so that the three
C values are in ascending order.  Interchanges are then applied so that
C all elements greater than T are in the upper portion of the array and
C all elements less than T are in the lower portion.  The upper and
C lower indices of one of the portions are saved in local arrays, and
C the process is repeated iteratively on the other portion.  When a
C portion is completely sorted, the process begins again by retrieving
C the indices bounding another unsorted portion.
C
C Input parameters -   N - Length of the array X.
C
C                      X - Vector of length N to be sorted.
C
C                    IND - Vector of length .GE. N.
C
C N and X are not altered by this routine.
C
C Output parameter - IND - Sequence of indices 1,...,N
C                          permuted in the same fashion as X
C                          would be.  Thus, the ordering on
C                          X is defined by Y(I) = X(IND(I)).
C
C Intrinsic functions called by QSORTI - IFIX, FLOAT
C
C***********************************************************************
C
C NOTE -- IU and IL must be dimensioned .GE. log(N) where log has base 2.
C
C***********************************************************************
      SUBROUTINE QSORTI (X, IND, N)
      INTEGER N, X(N), IND(N)
      INTEGER IU(32), IL(32)
      INTEGER M, I, J, K, L, IJ, IT, ITT, INDX, T
      REAL    R

C Local parameters -
C
C IU,IL =  Temporary storage for the upper and lower
C            indices of portions of the array X
C M =      Index for IU and IL
C I,J =    Lower and upper indices of a portion of X
C K,L =    Indices in the range I,...,J
C IJ =     Randomly chosen index between I and J
C IT,ITT = Temporary storage for interchanges in IND
C INDX =   Temporary index for X
C R =      Pseudo random number for generating IJ
C T =      Central element of X

      IF (N .LE. 0) RETURN

C Initialize IND, M, I, J, and R
      DO 1 I = 1,N
    1   IND(I) = I
      M = 1
      I = 1
      J = N
      R = .375

C Top of loop
    2 IF (I .GE. J) GO TO 10
      IF (R .GT. .5898437) GO TO 3
      R = R + .0390625
      GO TO 4
    3 R = R - .21875

C Initialize K
    4 K = I

C Select a central element of X and save it in T
      IJ = I + IFIX(R*FLOAT(J-I))
      IT = IND(IJ)
      T = X(IT)

C If the first element of the array is greater than T,
C   interchange it with T
      INDX = IND(I)
      IF (X(INDX) .LE. T) GO TO 5
      IND(IJ) = INDX
      IND(I) = IT
      IT = INDX
      T = X(IT)

C Initialize L
    5 L = J

C If the last element of the array is less than T,
C   interchange it with T
      INDX = IND(J)
      IF (X(INDX) .GE. T) GO TO 7
      IND(IJ) = INDX
      IND(J) = IT
      IT = INDX
      T = X(IT)

C If the first element of the array is greater than T,
C   interchange it with T
      INDX = IND(I)
      IF (X(INDX) .LE. T) GO TO 7
      IND(IJ) = INDX
      IND(I) = IT
      IT = INDX
      T = X(IT)
      GO TO 7

C Interchange elements K and L
    6 ITT = IND(L)
      IND(L) = IND(K)
      IND(K) = ITT

C Find an element in the upper part of the array which is
C   not larger than T
    7 L = L - 1
      INDX = IND(L)
      IF (X(INDX) .GT. T) GO TO 7

C Find an element in the lower part of the array which is
C   not smaller than T
    8 K = K + 1
      INDX = IND(K)
      IF (X(INDX) .LT. T) GO TO 8

C If K .LE. L, interchange elements K and L
      IF (K .LE. L) GO TO 6

C Save the upper and lower subscripts of the portion of the
C   array yet to be sorted
      IF (L-I .LE. J-K) GO TO 9
      IL(M) = I
      IU(M) = L
      I = K
      M = M + 1
      GO TO 11

    9 IL(M) = K
      IU(M) = J
      J = L
      M = M + 1
      GO TO 11

C Begin again on another unsorted portion of the array
   10 M = M - 1
      IF (M .EQ. 0) RETURN
      I = IL(M)
      J = IU(M)

   11 IF (J-I .GE. 11) GO TO 4
      IF (I .EQ. 1) GO TO 2
      I = I - 1

C Sort elements I+1,...,J.  Note that 1 .LE. I < J and J-I < 11.
   12 I = I + 1
      IF (I .EQ. J) GO TO 10
      INDX = IND(I+1)
      T = X(INDX)
      IT = INDX
      INDX = IND(I)
      IF (X(INDX) .LE. T) GO TO 12
      K = I

   13 IND(K+1) = IND(K)
      K = K - 1
      INDX = IND(K)
      IF (T .LT. X(INDX)) GO TO 13
      IND(K+1) = IT
      GO TO 12
      END  ! of QSORTI


C***********************************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C   This subroutine uses an order N*log(N) quick sort to sort an integer
C array X into increasing order.  The algorithm is as follows.  IND is
C initialized to the ordered sequence of indices 1,...,N, and all
C interchanges are applied to IND.  X is divided into two portions by
C picking a central element T.  The first and last elements are compared
C with T, and interchanges are applied as necessary so that the three
C values are in ascending order.  Interchanges are then applied so that
C all elements greater than T are in the upper portion of the array and
C all elements less than T are in the lower portion.  The upper and
C lower indices of one of the portions are saved in local arrays, and
C the process is repeated iteratively on the other portion.  When a
C portion is completely sorted, the process begins again by retrieving
C the indices bounding another unsorted portion.
C
C Input parameters -   N - length of the array X.
C
C                      X - Vector of length N to be sorted.
C
C                    IND - Vector of length .GE. N.
C
C N and X are not altered by this routine.
C
C Output parameter - IND - Sequence of indices 1,...,N
C                          permuted in the same fashion as X
C                          would be.  Thus, the ordering on
C                          X is defined by Y(I) = X(IND(I)).
C
C Intrinsic functions called by QSORTR - IFIX, FLOAT
C
C***********************************************************************
C
C NOTE -- IU and IL must be dimensioned .GE. log(N) where log has base 2.
C
C***********************************************************************
      SUBROUTINE QSORTR (X, IND, N)
      INTEGER N, IND(N)
      REAL    X(N)

      INTEGER IU(32), IL(32)
      INTEGER M, I, J, K, L, IJ, IT, ITT, INDX
      REAL    R, T

C Local parameters -
C
C IU,IL =  Temporary storage for the upper and lower
C            indices of portions of the array X
C M =      Index for IU and IL
C I,J =    Lower and upper indices of a portion of X
C K,L =    Indices in the range I,...,J
C IJ =     Randomly chosen index between I and J
C IT,ITT = Temporary storage for interchanges in IND
C INDX =   Temporary index for X
C R =      Pseudo random number for generating IJ
C T =      Central element of X

      IF (N .LE. 0) RETURN

C Initialize IND, M, I, J, and R
      DO 1 I = 1,N
    1   IND(I) = I
      M = 1
      I = 1
      J = N
      R = .375

C Top of loop
    2 IF (I .GE. J) GO TO 10
      IF (R .GT. .5898437) GO TO 3
      R = R + .0390625
      GO TO 4
    3 R = R - .21875

C Initialize K
    4 K = I

C Select a central element of X and save it in T
      IJ = I + IFIX(R*FLOAT(J-I))
      IT = IND(IJ)
      T = X(IT)

C If the first element of the array is greater than T,
C   interchange it with T
      INDX = IND(I)
      IF (X(INDX) .LE. T) GO TO 5
      IND(IJ) = INDX
      IND(I) = IT
      IT = INDX
      T = X(IT)

C Initialize L
    5 L = J

C If the last element of the array is less than T,
C   interchange it with T
      INDX = IND(J)
      IF (X(INDX) .GE. T) GO TO 7
      IND(IJ) = INDX
      IND(J) = IT
      IT = INDX
      T = X(IT)

C If the first element of the array is greater than T,
C   interchange it with T
      INDX = IND(I)
      IF (X(INDX) .LE. T) GO TO 7
      IND(IJ) = INDX
      IND(I) = IT
      IT = INDX
      T = X(IT)
      GO TO 7

C Interchange elements K and L
    6 ITT = IND(L)
      IND(L) = IND(K)
      IND(K) = ITT

C Find an element in the upper part of the array which is
C   not larger than T
    7 L = L - 1
      INDX = IND(L)
      IF (X(INDX) .GT. T) GO TO 7

C Find an element in the lower part of the array which is
C   not smaller than T
    8 K = K + 1
      INDX = IND(K)
      IF (X(INDX) .LT. T) GO TO 8

C If K .LE. L, interchange elements K and L
      IF (K .LE. L) GO TO 6

C Save the upper and lower subscripts of the portion of the
C   array yet to be sorted
      IF (L-I .LE. J-K) GO TO 9
      IL(M) = I
      IU(M) = L
      I = K
      M = M + 1
      GO TO 11
    9 IL(M) = K
      IU(M) = J
      J = L
      M = M + 1
      GO TO 11

C Begin again on another unsorted portion of the array
   10 M = M - 1
      IF (M .EQ. 0) RETURN
      I = IL(M)
      J = IU(M)

   11 IF (J-I .GE. 11) GO TO 4
      IF (I .EQ. 1) GO TO 2
      I = I - 1

C Sort elements I+1,...,J.  Note that 1 .LE. I < J and J-I < 11.
   12 I = I + 1
      IF (I .EQ. J) GO TO 10
      INDX = IND(I+1)
      T = X(INDX)
      IT = INDX
      INDX = IND(I)
      IF (X(INDX) .LE. T) GO TO 12
      K = I

   13 IND(K+1) = IND(K)
      K = K - 1
      INDX = IND(K)
      IF (T .LT. X(INDX)) GO TO 13
      IND(K+1) = IT
      GO TO 12
      END  ! of QSORTR


C***********************************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C   This routine applies a set of permutations to a vector.
C
C Input parameters -  N - Length of A and IP.
C
C                    IP - Vector containing the sequence of
C                         integers 1,...,N permuted in the
C                         same fashion that a is to be permuted.
C
C                     A - Vector to be permuted.
C
C N and IP are not altered by this routine.
C
C Output parameter -  A - Reordered vector reflecting the
C                         permutations defined by IP.
C
C***********************************************************************
      SUBROUTINE IORDER (A, IP, N)
      INTEGER N, A(N), IP(N)

      INTEGER NN, K, J, IPJ, TEMP

C Local parameters -
C
C NN =   Local copy of N
C K =    Index for IP and for the first element of A in a permutation
C J =    Index for IP and A, J .GE. K
C IPJ =  IP(J)
C TEMP = Temporary storage for A(K)

      NN = N
      IF (NN .LT. 2) RETURN
      K = 1

C Loop on permutations
    1 J = K
      TEMP = A(K)

C Apply permutation to A.  IP(J) is marked (made negative)
C   as being included in the permutation.
    2 IPJ = IP(J)
      IP(J) = -IPJ
      IF (IPJ .EQ. K) GO TO 3
      A(J) = A(IPJ)
      J = IPJ
      GO TO 2
    3 A(J) = TEMP

C Search for an unmarked element of IP
    4 K = K + 1
      IF (K .GT. NN) GO TO 5
      IF (IP(K) .GT. 0) GO TO 1
      GO TO 4

C All permutations have been applied.  Unmark IP.
    5 DO 6 K = 1,NN
    6   IP(K) = -IP(K)
      RETURN
      END  ! of IORDER


C***********************************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C   This routine applies a set of permutations to a vector.
C
C Input parameters -  N - Length of A and IP.
C
C                    IP - Vector containing the sequence of
C                         integers 1,...,N permuted in the
C                         same fashion that a is to be permuted.
C
C                     A - Vector to be permuted.
C
C N and IP are not altered by this routine.
C
C Output parameter -  A - Reordered vector reflecting the
C                         permutations defined by IP.
C
C***********************************************************************
      SUBROUTINE RORDER (A, IP, N)
      INTEGER N, IP(N)
      REAL    A(N)

      INTEGER NN, K, J, IPJ
      REAL    TEMP

C Local parameters -
C
C NN =   Local copy of N
C K =    Index for IP and for the first element of A in a permutation
C J =    INDEX for IP and A, J .GE. K
C IPJ =  IP(J)
C TEMP = Temporary storage for A(K)

      NN = N
      IF (NN .LT. 2) RETURN
      K = 1

C Loop on permutations
    1 J = K
      TEMP = A(K)

C Apply permutation to A.  IP(J) is marked (made negative)
C   as being included in the permutation.
    2 IPJ = IP(J)
      IP(J) = -IPJ
      IF (IPJ .EQ. K) GO TO 3
      A(J) = A(IPJ)
      J = IPJ
      GO TO 2
    3 A(J) = TEMP

C Search for an unmarked element of IP
    4 K = K + 1
      IF (K .GT. NN) GO TO 5
      IF (IP(K) .GT. 0) GO TO 1
      GO TO 4

C All permutations have been applied.  Unmark IP.
    5 DO 6 K = 1,NN
    6   IP(K) = -IP(K)
      RETURN
      END  ! of RORDER


C***********************************************************************
C HPSORT
C   PURPOSE  Return the permutation vector generated by sorting a
C            substring within a character array and, optionally,
C            rearrange the elements of the array.  The array may be
C            sorted in forward or reverse lexicographical order.  A
C            slightly modified quicksort algorithm is used.
C   LIBRARY   SLATEC
C   CATEGORY  N6A1C, N6A2C
C   KEYWORDS  PASSIVE SORTING, SINGLETON QUICKSORT, SORT, STRING SORTING
C   AUTHOR  Jones, R. E., (SNLA)
C           Rhoads, G. S., (NBS)
C           Sullivan, F. E., (NBS)
C           Wisniewski, J. A., (SNLA)
C
C   DESCRIPTION
C   HPSORT returns the permutation vector IPERM generated by sorting
C   the substrings beginning with the character STRBEG and ending with
C   the character STREND within the strings in array HX and, optionally,
C   rearranges the strings in HX.   HX may be sorted in increasing or
C   decreasing lexicographical order.  A slightly modified quicksort
C   algorithm is used.
C
C   IPERM is such that HX(IPERM(I)) is the Ith value in the
C   rearrangement of HX.  IPERM may be applied to another array by
C   calling IPPERM, SPPERM, DPPERM or HPPERM.
C
C   An active sort of numerical data is expected to execute somewhat
C   more quickly than a passive sort because there is no need to use
C   indirect references. But for the character data in HPSORT, integers
C   in the IPERM vector are manipulated rather than the strings in HX.
C   Moving integers may be enough faster than moving character strings
C   to more than offset the penalty of indirect referencing.
C
C   Description of Parameters
C      HX - input/output -- array of type character to be sorted.
C           For example, to sort a 80 element array of names,
C           each of length 6, declare HX as character HX(100)*6.
C           If ABS(KFLAG) = 2, then the values in HX will be
C           rearranged on output; otherwise, they are unchanged.
C      N  - input -- number of values in array HX to be sorted.
C      STRBEG - input -- the index of the initial character in
C               the string HX that is to be sorted.
C      STREND - input -- the index of the final character in
C               the string HX that is to be sorted.
C      IPERM - output -- permutation array such that IPERM(I) is the
C              index of the string in the original order of the
C              HX array that is in the Ith location in the sorted
C              order.
C      KFLAG - input -- control parameter:
C            =  2  means return the permutation vector resulting from
C                  sorting HX in lexicographical order and sort HX also.
C            =  1  means return the permutation vector resulting from
C                  sorting HX in lexicographical order and do not sort
C                  HX.
C            = -1  means return the permutation vector resulting from
C                  sorting HX in reverse lexicographical order and do
C                  not sort HX.
C            = -2  means return the permutation vector resulting from
C                  sorting HX in reverse lexicographical order and sort
C                  HX also.
C      WORK - character variable which must have a length specification
C             at least as great as that of HX.
C      IER - output -- error indicator:
C          =  0  if no error,
C          =  1  if N is zero or negative,
C          =  2  if KFLAG is not 2, 1, -1, or -2,
C          =  3  if work array is not long enough,
C          =  4  if string beginning is beyond its end,
C          =  5  if string beginning is out-of-range,
C          =  6  if string end is out-of-range.
C
C     E X A M P L E  O F  U S E
C
C      CHARACTER*2 HX, W
C      INTEGER STRBEG, STREND
C      DIMENSION HX(10), IPERM(10)
C      DATA (HX(I),I=1,10)/ '05','I ',' I','  ','Rs','9R','R9','89',
C     1     ',*','N"'/
C      DATA STRBEG, STREND / 1, 2 /
C      CALL HPSORT (HX,10,STRBEG,STREND,IPERM,1,W)
C      PRINT 100, (HX(IPERM(I)),I=1,10)
C 100 FORMAT (2X, A2)
C      STOP
C      END
C
C   REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
C                 for sorting with minimal storage, Communications of
C                 the ACM, 12, 3 (1969), pp. 185-187.
C   REVISION HISTORY  (YYMMDD)
C   761101  DATE WRITTEN
C   761118  Modified by John A. Wisniewski to use the Singleton
C           quicksort algorithm.
C   811001  Modified by Francis Sullivan for string data.
C   850326  Documentation slightly modified by D. Kahaner.
C   870423  Modified by Gregory S. Rhoads for passive sorting with the
C           option for the rearrangement of the original data.
C   890620  Algorithm for rearranging the data vector corrected by R.
C           Boisvert.
C   890622  Prologue upgraded to Version 4.0 style by D. Lozier.
C   920507  Modified by M. McClain to revise prologue text.
C   920818  Declarations section rebuilt and code restructured to use
C           IF-THEN-ELSE-ENDIF.  (SMR, WRB)
C   170523  Deleted calls to XERMSG, increased array bounds.  (AJA)
C***********************************************************************
      SUBROUTINE HPSORT (HX, N, STRBEG, STREND, IPERM, KFLAG, WORK, IER)
C     .. Scalar Arguments ..
      INTEGER IER, KFLAG, N, STRBEG, STREND
      CHARACTER * (*) WORK
C     .. Array Arguments ..
      INTEGER IPERM(*)
      CHARACTER * (*) HX(*)
C     .. Local Scalars ..
      REAL R
      INTEGER I, IJ, INDX, INDX0, IR, ISTRT, J, K, KK, L, LM, LMT, M,
     $        NN, NN2
C     .. Local Arrays ..
      INTEGER IL(32), IU(32)
C     .. Intrinsic Functions ..
      INTRINSIC ABS, INT, LEN
C   First executable statement
      IER = 0
      NN = N
      IF (NN .LT. 1) THEN
         IER = 1  ! The number of values to be sorted, N, is not positive.
         RETURN
      ENDIF
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
         IER = 2  ! The sort control parameter, KFLAG, is not 2, 1, -1, or -2.
         RETURN
      ENDIF
      IF(LEN(WORK) .LT. LEN(HX(1))) THEN
         IER = 3  ! The length of the work variable, WORK, is too short.
         RETURN
      ENDIF
      IF (STRBEG .GT. STREND) THEN
         IER = 4  ! The string beginning, STRBEG, is beyond its end, STREND.
         RETURN
      ENDIF
      IF (STRBEG .LT. 1 .OR. STRBEG .GT. LEN(HX(1))) THEN
         IER = 5  ! The string beginning, STRBEG, is out-of-range.
         RETURN
      ENDIF
      IF (STREND .LT. 1 .OR. STREND .GT. LEN(HX(1))) THEN
         IER = 6  ! The string end, STREND, is out-of-range.
         RETURN
      ENDIF

C     Initialize permutation vector
      DO 10 I=1,NN
         IPERM(I) = I
   10 CONTINUE

C     Return if only one value is to be sorted
      IF (NN .EQ. 1) RETURN

C     Sort HX only
      M = 1
      I = 1
      J = NN
      R = .375E0

   20 IF (I .EQ. J) GO TO 70
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF

   30 K = I

C     Select a central element of the array and save it in location L
      IJ = I + INT((J-I)*R)
      LM = IPERM(IJ)

C     If first element of array is greater than LM, interchange with LM
      IF (HX(IPERM(I))(STRBEG:STREND) .GT. HX(LM)(STRBEG:STREND)) THEN
         IPERM(IJ) = IPERM(I)
         IPERM(I) = LM
         LM = IPERM(IJ)
      ENDIF
      L = J

C     If last element of array is less than LM, interchange with LM
      IF (HX(IPERM(J))(STRBEG:STREND) .LT. HX(LM)(STRBEG:STREND)) THEN
         IPERM(IJ) = IPERM(J)
         IPERM(J) = LM
         LM = IPERM(IJ)

C If first element of array is greater than LM, interchange with LM
         IF (HX(IPERM(I))(STRBEG:STREND) .GT. HX(LM)(STRBEG:STREND))
     $      THEN
               IPERM(IJ) = IPERM(I)
               IPERM(I) = LM
               LM = IPERM(IJ)
         ENDIF
      ENDIF
      GO TO 50
   40 LMT = IPERM(L)
      IPERM(L) = IPERM(K)
      IPERM(K) = LMT

C Find an element in the second half of the array which is smaller
C     than LM
   50 L = L-1
      IF (HX(IPERM(L))(STRBEG:STREND) .GT. HX(LM)(STRBEG:STREND))
     $    GO TO 50

C     Find an element in the first half of the array which is greater
C     than LM
   60 K = K+1
      IF (HX(IPERM(K))(STRBEG:STREND) .LT. HX(LM)(STRBEG:STREND))
     $   GO TO 60

C     Interchange these elements
      IF (K .LE. L) GO TO 40

C     Save upper and lower subscripts of the array yet to be sorted
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 80

C     Begin again on another portion of the unsorted array
   70 M = M-1
      IF (M .EQ. 0) GO TO 110
      I = IL(M)
      J = IU(M)

   80 IF (J-I .GE. 1) GO TO 30
      IF (I .EQ. 1) GO TO 20
      I = I-1

   90 I = I+1
      IF (I .EQ. J) GO TO 70
      LM = IPERM(I+1)
      IF (HX(IPERM(I))(STRBEG:STREND) .LE. HX(LM)(STRBEG:STREND))
     $   GO TO 90
      K = I

  100 IPERM(K+1) = IPERM(K)
      K = K-1

      IF (HX(LM)(STRBEG:STREND) .LT. HX(IPERM(K))(STRBEG:STREND))
     $    GO TO 100
      IPERM(K+1) = LM
      GO TO 90

C     Clean up
  110 IF (KFLAG .LE. -1) THEN

C        Alter array to get reverse order, if necessary
         NN2 = NN/2
         DO 120 I=1,NN2
           IR = NN-I+1
           LM = IPERM(I)
           IPERM(I) = IPERM(IR)
           IPERM(IR) = LM
  120    CONTINUE
      ENDIF

C     Rearrange the values of HX if desired
      IF (KK .EQ. 2) THEN

C        Use the IPERM vector as a flag.
C        If IPERM(I) < 0, then the I-th value is in correct location

         DO 140 ISTRT=1,NN
            IF (IPERM(ISTRT) .GE. 0) THEN
               INDX = ISTRT
               INDX0 = INDX
               WORK = HX(ISTRT)
  130          IF (IPERM(INDX) .GT. 0) THEN
                  HX(INDX) = HX(IPERM(INDX))
                  INDX0 = INDX
                  IPERM(INDX) = -IPERM(INDX)
                  INDX = ABS(IPERM(INDX))
                  GO TO 130
               ENDIF
               HX(INDX0) = WORK
            ENDIF
  140    CONTINUE

C        Revert the signs of the IPERM values
         DO 150 I=1,NN
            IPERM(I) = -IPERM(I)
  150    CONTINUE

      ENDIF

      RETURN
      END  ! of HPSORT


C***********************************************************************
C HPPERM
C   PURPOSE  Rearrange a given array according to a prescribed
C            permutation vector.
C   LIBRARY   SLATEC
C   CATEGORY  N8
C   AUTHOR  McClain, M. A., (NIST)
C           Rhoads, G. S., (NBS)
C
C   DESCRIPTION
C         HPPERM rearranges the data vector HX according to the
C         permutation IPERM: HX(I) <--- HX(IPERM(I)).  IPERM could come
C         from one of the sorting routines IPSORT, SPSORT, DPSORT or
C         HPSORT.
C
C     Description of Parameters
C         HX - input/output -- character array of values to be rearranged.
C         N - input -- number of values in character array HX.
C         IPERM - input -- permutation vector.
C         WORK - character variable which must have a length
C                   specification at least as great as that of HX.
C         IER - output -- error indicator:
C             =  0  if no error,
C             =  1  if N is zero or negative,
C             =  2  if work array is not long enough,
C             =  3  if IPERM is not a valid permutation.
C
C   REVISION HISTORY  (YYMMDD)
C   901004  DATE WRITTEN
C   920507  Modified by M. McClain to revise prologue text and to add
C           check for length of work array.
C   170604  Deleted calls to XERMSG
C***********************************************************************
      SUBROUTINE HPPERM (HX, N, IPERM, WORK, IER)
      INTEGER N, IPERM(*), I, IER, INDX, INDX0, ISTRT
      CHARACTER*(*) HX(*), WORK

C   First executable statement
      IER=0
      IF (N .LT. 1) THEN
         IER = 1  ! The number of values to be rearranged, N, is not positive.
         RETURN
      ENDIF
      IF (LEN(WORK) .LT. LEN(HX(1))) THEN
         IER = 2  ! The length of the work variable, WORK, is too short.
         RETURN
      ENDIF

C     Check whether IPERM is a valid permutation
      DO 100 I = 1, N
         INDX = ABS(IPERM(I))
         IF ((INDX .GE. 1) .AND. (INDX .LE. N)) THEN
            IF (IPERM(INDX) .GT. 0) THEN
               IPERM(INDX) = -IPERM(INDX)
               GOTO 100
            ENDIF
         ENDIF
         IER = 3  ! The permutation vector, IPERM, is not valid.
         RETURN
  100 CONTINUE

C     Rearrange the values of HX
C     Use the IPERM vector as a flag.
C     If IPERM(I) > 0, then the I-th value is in correct location
      DO 330 ISTRT = 1, N
         IF (IPERM(ISTRT) .GT. 0) GOTO 330
         INDX = ISTRT
         INDX0 = INDX
         WORK = HX(ISTRT)
  320    CONTINUE
         IF (IPERM(INDX) .GE. 0) GOTO 325
            HX(INDX) = HX(-IPERM(INDX))
            INDX0 = INDX
            IPERM(INDX) = -IPERM(INDX)
            INDX = IPERM(INDX)
            GOTO 320
  325    CONTINUE
         HX(INDX0) = WORK
  330 CONTINUE

      RETURN
      END  ! of HPPERM
*+++++++++++++++++++++++++ End of file sort.f ++++++++++++++++++++++++++
