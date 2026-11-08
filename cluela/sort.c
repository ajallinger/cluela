/*----------------------------------------------------------------------
sort.c - Routines for sorting and ordering, translated from the SLATEC
and Naval Surface Warfare Center libraries.  This file is public domain.
----------------------------------------------------------------------*/
#include <stdlib.h>
#include <string.h>

/***********************************************************************

                                              Robert Renka
                                      Oak Ridge Natl. Lab.

  This subroutine uses an order N*log(N) quick sort to sort an integer
array X into increasing order.  The algorithm is as follows.  IND is
initialized to the ordered sequence of indices 1,...,N, and all
interchanges are applied to IND.  X is divided into two portions by
picking a central element T.  The first and last elements are compared
with T, and interchanges are applied as necessary so that the three
values are in ascending order.  Interchanges are then applied so that
all elements greater than T are in the upper portion of the array and
all elements less than T are in the lower portion.  The upper and
lower indices of one of the portions are saved in local arrays, and
the process is repeated iteratively on the other portion.  When a
portion is completely sorted, the process begins again by retrieving
the indices bounding another unsorted portion.

Input parameters -   N - Length of the array X.

                     X - Vector of length N to be sorted.

                   IND - Vector of length .GE. N.

N and X are not altered by this routine.

Output parameter - IND - Sequence of indices 1,...,N
                         permuted in the same fashion as X
                         would be.  Thus, the ordering on
                         X is defined by Y(I) = X(IND(I)).

Intrinsic functions called by QSORTI - IFIX, FLOAT

***********************************************************************

NOTE -- IU and IL must be dimensioned .GE. log(N) where log has base 2.

***********************************************************************/
void qsorti (int *x, int *ind, int n)
{
/* Local variables */
	int i, j, k, l, m;
	float r;
	int t, ij, il[32], it, iu[32], itt, indx;

/***********************************************************************
Local parameters -
IU,IL =  Temporary storage for the upper and lower
           indices of portions of the array X
M =      Index for IU and IL
I,J =    Lower and upper indices of a portion of X
K,L =    Indices in the range I,...,J
IJ =     Randomly chosen index between I and J
IT,ITT = Temporary storage for interchanges in IND
INDX =   Temporary index for X
R =      Pseudo random number for generating IJ
T =      Central element of X
***********************************************************************/

/* Parameter adjustments */
	--x;
	--ind;
/* Function Body */
	if (n <= 0) return;
/* Initialize IND, M, I, J, and R */
	for (i = 1; i <= n; ++i) ind[i] = i;
	m = 1;
	i = 1;
	j = n;
	r = .375f;
/* Top of loop */
L2:
	if (i >= j) goto L10;
	if (r > .5898437f) {
		r -= .21875f;
	} else {
		r += .0390625f;
	}
/* Initialize K */
L4:
	k = i;
/* Select a central element of X and save it in T */
	ij = i + (int) (r * (float) (j - i));
	it = ind[ij];
	t = x[it];
/* If the first element of the array is greater than T, interchange it with T */
	indx = ind[i];
	if (x[indx] <= t) goto L5;
	ind[ij] = indx;
	ind[i] = it;
	it = indx;
	t = x[it];
/* Initialize L */
L5:
	l = j;
/* If the last element of the array is less than T, interchange it with T */
	indx = ind[j];
	if (x[indx] >= t) goto L7;
	ind[ij] = indx;
	ind[j] = it;
	it = indx;
	t = x[it];
/* If the first element of the array is greater than T, interchange it with T */
	indx = ind[i];
	if (x[indx] <= t) goto L7;
	ind[ij] = indx;
	ind[i] = it;
	it = indx;
	t = x[it];
	goto L7;
/* Interchange elements K and L */
L6:
	itt = ind[l];
	ind[l] = ind[k];
	ind[k] = itt;
/* Find an element in the upper part of the array which is not larger than T */
L7:
	do {
		--l;
		indx = ind[l];
	} while (x[indx] > t);
/* Find an element in the lower part of the array which is not smaller than T */
	do {
		++k;
		indx = ind[k];
	} while (x[indx] < t);
/* If K .LE. L, interchange elements K and L */
	if (k <= l) goto L6;
/* Save the upper and lower subscripts of the portion of the
   array yet to be sorted */
	if (l - i <= j - k) {
		il[m-1] = k;
		iu[m-1] = j;
		j = l;
		++m;
	} else {
		il[m-1] = i;
		iu[m-1] = l;
		i = k;
		++m;
	}
	goto L11;
/* Begin again on another unsorted portion of the array */
L10:
	--m;
	if (0 == m) return;
	i = il[m-1];
	j = iu[m-1];
L11:
	if (j - i >= 11) goto L4;
	if (1 == i) goto L2;
	--i;
/* Sort elements I+1,...,J.  Note that 1 .LE. I < J and J-I < 11. */
L12:
	++i;
	if (i == j) goto L10;
	indx = ind[i+1];
	t = x[indx];
	it = indx;
	indx = ind[i];
	if (x[indx] <= t) goto L12;
	k = i;
	do {
		ind[k+1] = ind[k];
		--k;
		indx = ind[k];
	} while (t < x[indx]);
	ind[k+1] = it;
	goto L12;
} /* end of qsorti */

/***********************************************************************

                                              Robert Renka
                                      Oak Ridge Natl. Lab.

  This subroutine uses an order N*log(N) quick sort to sort an integer
array X into increasing order.  The algorithm is as follows.  IND is
initialized to the ordered sequence of indices 1,...,N, and all
interchanges are applied to IND.  X is divided into two portions by
picking a central element T.  The first and last elements are compared
with T, and interchanges are applied as necessary so that the three
values are in ascending order.  Interchanges are then applied so that
all elements greater than T are in the upper portion of the array and
all elements less than T are in the lower portion.  The upper and
lower indices of one of the portions are saved in local arrays, and
the process is repeated iteratively on the other portion.  When a
portion is completely sorted, the process begins again by retrieving
the indices bounding another unsorted portion.

Input parameters -   N - length of the array X.

                     X - Vector of length N to be sorted.

                   IND - Vector of length .GE. N.

N and X are not altered by this routine.

Output parameter - IND - Sequence of indices 1,...,N
                         permuted in the same fashion as X
                         would be.  Thus, the ordering on
                         X is defined by Y(I) = X(IND(I)).

Intrinsic functions called by QSORTR - IFIX, FLOAT

***********************************************************************

NOTE -- IU and IL must be dimensioned .GE. log(N) where log has base 2.

***********************************************************************/
void qsortr (float *x, int *ind, int n)
{
/* Local variables */
	int i, j, k, l, m;
	float r, t;
	int ij, il[32], it, iu[32], itt, indx;

/***********************************************************************
Local parameters -

IU,IL =  Temporary storage for the upper and lower
           indices of portions of the array X
M =      Index for IU and IL
I,J =    Lower and upper indices of a portion of X
K,L =    Indices in the range I,...,J
IJ =     Randomly chosen index between I and J
IT,ITT = Temporary storage for interchanges in IND
INDX =   Temporary index for X
R =      Pseudo random number for generating IJ
T =      Central element of X
***********************************************************************/

/* Parameter adjustments */
	--x;
	--ind;
/* Function Body */
	if (n <= 0) return;
/* Initialize IND, M, I, J, and R */
	for (i = 1; i <= n; ++i) ind[i] = i;
	m = 1;
	i = 1;
	j = n;
	r = .375f;
/* Top of loop */
L2:
	if (i >= j) goto L10;
	if (r > .5898437f) {
		r -= .21875f;
	} else {
		r += .0390625f;
	}
/* Initialize K */
L4:
	k = i;
/* Select a central element of X and save it in T */
	ij = i + (int) (r * (float) (j - i));
	it = ind[ij];
	t = x[it];
/* If the first element of the array is greater than T, interchange it with T */
	indx = ind[i];
	if (x[indx] <= t) goto L5;
	ind[ij] = indx;
	ind[i] = it;
	it = indx;
	t = x[it];
/* Initialize L */
L5:
	l = j;
/* If the last element of the array is less than T, interchange it with T */
	indx = ind[j];
	if (x[indx] >= t) goto L7;
	ind[ij] = indx;
	ind[j] = it;
	it = indx;
	t = x[it];
/* If the first element of the array is greater than T, interchange it with T */
	indx = ind[i];
	if (x[indx] <= t) goto L7;
	ind[ij] = indx;
	ind[i] = it;
	it = indx;
	t = x[it];
	goto L7;
/* Interchange elements K and L */
L6:
	itt = ind[l];
	ind[l] = ind[k];
	ind[k] = itt;
/* Find an element in the upper part of the array which is not larger than T */
L7:
	do {
		--l;
		indx = ind[l];
	} while (x[indx] > t);
/* Find an element in the lower part of the array which is not smaller than T */
	do {
		++k;
		indx = ind[k];
	} while (x[indx] < t);
/* If K .LE. L, interchange elements K and L */
	if (k <= l) goto L6;
/* Save the upper and lower subscripts of the portion of the
   array yet to be sorted */
	if (l - i <= j - k) {
		il[m-1] = k;
		iu[m-1] = j;
		j = l;
		++m;
	} else {
		il[m-1] = i;
		iu[m-1] = l;
		i = k;
		++m;
	}
	goto L11;
/* Begin again on another unsorted portion of the array */
L10:
	--m;
	if (0 == m) return;
	i = il[m-1];
	j = iu[m-1];
L11:
	if (j - i >= 11) goto L4;
	if (1 == i) goto L2;
	--i;
/* Sort elements I+1,...,J.  Note that 1 .LE. I < J and J-I < 11. */
L12:
	++i;
	if (i == j) goto L10;
	indx = ind[i+1];
	t = x[indx];
	it = indx;
	indx = ind[i];
	if (x[indx] <= t) goto L12;
	k = i;
	do {
		ind[k+1] = ind[k];
		--k;
		indx = ind[k];
	} while (t < x[indx]);
	ind[k+1] = it;
	goto L12;
} /* end of qsortr */

/***********************************************************************

                                              Robert Renka
                                      Oak Ridge Natl. Lab.

  This routine applies a set of permutations to a vector.

Input parameters -  N - Length of A and IP.

                   IP - Vector containing the sequence of
                        integers 1,...,N permuted in the
                        same fashion that a is to be permuted.

                    A - Vector to be permuted.

N and IP are not altered by this routine.

Output parameter -  A - Reordered vector reflecting the
                        permutations defined by IP.

***********************************************************************/
void iorder (int *a, int *ip, int n)
{
/* Local variables */
    int j, k, ipj, temp;

/***********************************************************************
Local parameters -

K =    Index for IP and for the first element of A in a permutation
J =    Index for IP and A, J .GE. K
IPJ =  IP(J)
TEMP = Temporary storage for A(K)
***********************************************************************/

/* Parameter adjustments */
	--a;
	--ip;

/* Function Body */
	if (n < 2) return;
	k = 1;
/* Loop on permutations */
L1:
	j = k;
	temp = a[k];
/* Apply permutation to A.  IP(J) is marked (made negative)
   as being included in the permutation. */
L2:
	ipj = ip[j];
	ip[j] = -ipj;
	if (ipj == k) goto L3;
	a[j] = a[ipj];
	j = ipj;
	goto L2;
L3:
	a[j] = temp;
/* Search for an unmarked element of IP */
L4:
	++k;
	if (k > n) goto L5;
	if (ip[k] > 0) goto L1;
	goto L4;
/* All permutations have been applied.  Unmark IP. */
L5:
	for (k = 1; k <= n; ++k) ip[k] = -ip[k];
	return;
} /* end of iorder */

/***********************************************************************

                                              Robert Renka
                                      Oak Ridge Natl. Lab.

  This routine applies a set of permutations to a vector.

Input parameters -  N - Length of A and IP.

                   IP - Vector containing the sequence of
                        integers 1,...,N permuted in the
                        same fashion that a is to be permuted.

                    A - Vector to be permuted.

N and IP are not altered by this routine.

Output parameter -  A - Reordered vector reflecting the
                        permutations defined by IP.

***********************************************************************/
void rorder (float *a, int *ip, int n)
{
/* Local variables */
	int j, k, ipj;
	float temp;

/***********************************************************************
Local parameters -

K =    Index for IP and for the first element of A in a permutation
J =    INDEX for IP and A, J .GE. K
IPJ =  IP(J)
TEMP = Temporary storage for A(K)
***********************************************************************/

/* Parameter adjustments */
	--a;
	--ip;
/* Function Body */
	if (n < 2) return;
	k = 1;
/* Loop on permutations */
L1:
	j = k;
	temp = a[k];
/* Apply permutation to A.  IP(J) is marked (made negative)
   as being included in the permutation. */
L2:
	ipj = ip[j];
	ip[j] = -ipj;
	if (ipj == k) goto L3;
	a[j] = a[ipj];
	j = ipj;
	goto L2;
L3:
	a[j] = temp;
/* Search for an unmarked element of IP */
L4:
	++k;
	if (k > n) goto L5;
	if (ip[k] > 0) goto L1;
	goto L4;
/* All permutations have been applied.  Unmark IP. */
L5:
	for (k = 1; k <= n; ++k) ip[k] = -ip[k];
	return;
} /* end of rorder */

/***********************************************************************
HPSORT
  PURPOSE  Return the permutation vector generated by sorting a
           substring within a character array and, optionally,
           rearrange the elements of the array.  The array may be
           sorted in forward or reverse lexicographical order.  A
           slightly modified quicksort algorithm is used.
  LIBRARY   SLATEC
  CATEGORY  N6A1C, N6A2C
  KEYWORDS  PASSIVE SORTING, SINGLETON QUICKSORT, SORT, STRING SORTING
  AUTHOR  Jones, R. E., (SNLA)
          Rhoads, G. S., (NBS)
          Sullivan, F. E., (NBS)
          Wisniewski, J. A., (SNLA)

  DESCRIPTION
  HPSORT returns the permutation vector IPERM generated by sorting
  the substrings beginning with the character STRBEG and ending with
  the character STREND within the strings in array HX and, optionally,
  rearranges the strings in HX.   HX may be sorted in increasing or
  decreasing lexicographical order.  A slightly modified quicksort
  algorithm is used.

  IPERM is such that HX(IPERM(I)) is the Ith value in the
  rearrangement of HX.  IPERM may be applied to another array by
  calling IPPERM, SPPERM, DPPERM or HPPERM.

  An active sort of numerical data is expected to execute somewhat
  more quickly than a passive sort because there is no need to use
  indirect references. But for the character data in HPSORT, integers
  in the IPERM vector are manipulated rather than the strings in HX.
  Moving integers may be enough faster than moving character strings
  to more than offset the penalty of indirect referencing.

  Description of Parameters
     HX - input/output -- array of type character to be sorted.
          For example, to sort a 80 element array of names,
          each of length 6, declare HX as character HX(100)*6.
          If ABS(KFLAG) = 2, then the values in HX will be
          rearranged on output; otherwise, they are unchanged.
     N  - input -- number of values in array HX to be sorted.
     HX_LEN - input -- length of each entry in HX
     STRBEG - input -- the index of the initial character in
              the string HX that is to be sorted.
     STREND - input -- the index of the final character in
              the string HX that is to be sorted.
     IPERM - output -- permutation array such that IPERM(I) is the
             index of the string in the original order of the
             HX array that is in the Ith location in the sorted
             order.
     KFLAG - input -- control parameter:
           =  2  means return the permutation vector resulting from
                 sorting HX in lexicographical order and sort HX also.
           =  1  means return the permutation vector resulting from
                 sorting HX in lexicographical order and do not sort
                 HX.
           = -1  means return the permutation vector resulting from
                 sorting HX in reverse lexicographical order and do
                 not sort HX.
           = -2  means return the permutation vector resulting from
                 sorting HX in reverse lexicographical order and sort
                 HX also.
     WORK - character variable which must have a length specification
            at least as great as that of HX.
     RETURN VALUE -
         =  0  if no error,
         =  1  if N is zero or negative,
         =  2  if KFLAG is not 2, 1, -1, or -2,
         =  3  if work array is not long enough,
         =  4  if string beginning is beyond its end,
         =  5  if string beginning is out-of-range,
         =  6  if string end is out-of-range.

  REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
                for sorting with minimal storage, Communications of
                the ACM, 12, 3 (1969), pp. 185-187.
  REVISION HISTORY  (YYMMDD)
  761101  DATE WRITTEN
  761118  Modified by John A. Wisniewski to use the Singleton
          quicksort algorithm.
  811001  Modified by Francis Sullivan for string data.
  850326  Documentation slightly modified by D. Kahaner.
  870423  Modified by Gregory S. Rhoads for passive sorting with the
          option for the rearrangement of the original data.
  890620  Algorithm for rearranging the data vector corrected by R.
          Boisvert.
  890622  Prologue upgraded to Version 4.0 style by D. Lozier.
  920507  Modified by M. McClain to revise prologue text.
  920818  Declarations section rebuilt and code restructured to use
          IF-THEN-ELSE-ENDIF.  (SMR, WRB)
  170523  Deleted calls to XERMSG, increased array bounds.  (AJA)
***********************************************************************/
int hpsort (char **hx, int n, int hx_len, int strbeg, int strend,
            int *iperm, int kflag, char *work)
{
/* Local variables */
	int i, j, k, l, m;
	float r;
	int ij, il[32], kk, lm, ir, iu[32], nn2, lmt, indx, indx0, istrt;
	int bm1, embp1;

/* Parameter adjustments */
	--iperm;
/* Function Body */
	if (n < 1) return 1;
	kk = abs(kflag);
	if (kk != 1 && kk != 2) return 2;
	if (hx_len < 1) return 3;
	if (strbeg > strend) return 4;
	if (strbeg < 1 || strbeg > hx_len) return 5;
	if (strend < 1 || strend > hx_len) return 6;
/* Initialize permutation vector */
	for (i = 1; i <= n; ++i) iperm[i] = i;
	bm1 = strbeg - 1; /* offset of substring to sort */
	embp1 = strend - strbeg + 1; /* length of substring to sort */
/* Return if only one value is to be sorted */
	if (1 == n) return 0;
/* Sort HX only */
	m = 1;
	i = 1;
	j = n;
	r = .375f;
L20:
	if (i == j) goto L70;
	if (r <= .5898437f) {
		r += .0390625f;
	} else {
		r -= .21875f;
	}
L30:
	k = i;
/* Select a central element of the array and save it in location L */
	ij = i + (int) ((j - i) * r);
	lm = iperm[ij];
/* If first element of array is greater than LM, interchange with LM */
	if (strncmp(hx[iperm[i]-1]+bm1, hx[lm-1]+bm1, embp1) > 0) {
		iperm[ij] = iperm[i];
		iperm[i] = lm;
		lm = iperm[ij];
	}
	l = j;
/* If last element of array is less than LM, interchange with LM */
	if (strncmp(hx[iperm[j]-1]+bm1, hx[lm-1]+bm1, embp1) < 0) {
		iperm[ij] = iperm[j];
		iperm[j] = lm;
		lm = iperm[ij];
/* If first element of array is greater than LM, interchange with LM */
		if (strncmp(hx[iperm[i]-1]+bm1, hx[lm-1]+bm1, embp1) > 0) {
			iperm[ij] = iperm[i];
			iperm[i] = lm;
			lm = iperm[ij];
		}
	}
	goto L50;
L40:
	lmt = iperm[l];
	iperm[l] = iperm[k];
	iperm[k] = lmt;
/* Find an element in the second half of the array which is smaller than LM */
L50:
	do {
		--l;
	} while (strncmp(hx[iperm[l]-1]+bm1, hx[lm-1]+bm1, embp1) > 0);
/* Find an element in the first half of the array which is greater than LM */
	do {
		++k;
	} while (strncmp(hx[iperm[k]-1]+bm1, hx[lm-1]+bm1, embp1) < 0);
/* Interchange these elements */
	if (k <= l) goto L40;
/* Save upper and lower subscripts of the array yet to be sorted */
	if (l - i > j - k) {
		il[m-1] = i;
		iu[m-1] = l;
		i = k;
		++m;
	} else {
		il[m-1] = k;
		iu[m-1] = j;
		j = l;
		++m;
	}
	goto L80;
/* Begin again on another portion of the unsorted array */
L70:
	--m;
	if (0 == m) goto L110;
	i = il[m-1];
	j = iu[m-1];
L80:
	if (j - i >= 1) goto L30;
	if (1 == i) goto L20;
	--i;
L90:
	++i;
	if (i == j) goto L70;
	lm = iperm[i+1];
	if (strncmp(hx[iperm[i]-1]+bm1, hx[lm-1]+bm1, embp1) <= 0) goto L90;
	k = i;
	do {
		iperm[k+1] = iperm[k];
		--k;
	} while (strncmp(hx[lm-1]+bm1, hx[iperm[k]-1]+bm1, embp1) < 0);
	iperm[k+1] = lm;
	goto L90;
/* Clean up */
L110:
	if (kflag <= -1) {
/*Alter array to get reverse order, if necessary */
		nn2 = n / 2;
		for (i = 1; i <= nn2; ++i) {
			ir = n - i + 1;
			lm = iperm[i];
			iperm[i] = iperm[ir];
			iperm[ir] = lm;
		}
	}
/* Rearrange the values of HX if desired */
	if (2 == kk) {
/*Use the IPERM vector as a flag.
If IPERM(I) < 0, then the I-th value is in correct location */
		for (istrt = 1; istrt <= n; ++istrt) {
			if (iperm[istrt] >= 0) {
				indx = istrt;
				indx0 = indx;
				strncpy(work, hx[istrt-1], hx_len);
				while (iperm[indx] > 0) {
					strncpy(hx[indx-1], hx[iperm[indx]-1], hx_len);
					indx0 = indx;
					iperm[indx] = -iperm[indx];
					indx = abs(iperm[indx]);
				}
				strncpy(hx[indx0-1], work, hx_len);
			}
		}
/*Revert the signs of the IPERM values */
		for (i = 1; i <= n; ++i) iperm[i] = -iperm[i];
	}
	return 0;
} /* end of hpsort */

/***********************************************************************
HPPERM
  PURPOSE  Rearrange a given array according to a prescribed
           permutation vector.
  LIBRARY   SLATEC
  CATEGORY  N8
  AUTHOR  McClain, M. A., (NIST)
          Rhoads, G. S., (NBS)

  DESCRIPTION
        HPPERM rearranges the data vector HX according to the
        permutation IPERM: HX(I) <--- HX(IPERM(I)).  IPERM could come
        from one of the sorting routines IPSORT, SPSORT, DPSORT or
        HPSORT.

    Description of Parameters
        HX - input/output -- character array of values to be rearranged.
        N - input -- number of values in character array HX.
        HX - input -- length of each entry in HX
        IPERM - input -- permutation vector.
        WORK - character variable which must have a length
                  specification at least as great as that of HX.
        RETURN VALUE -
            =  0  if no error,
            =  1  if N is zero or negative,
            =  2  if work array is not long enough,
            =  3  if IPERM is not a valid permutation.

  REVISION HISTORY  (YYMMDD)
  901004  DATE WRITTEN
  920507  Modified by M. McClain to revise prologue text and to add
          check for length of work array.
  170604  Deleted calls to XERMSG
***********************************************************************/
int hpperm (char **hx, int n, int hx_len, int *iperm, char *work)
{
/* Local variables */
	int i, indx, indx0, istrt;
/* Parameter adjustments */
	--iperm;
/* Function Body */
	if (n < 1) return 1;
	if (hx_len < 1) return 2;
/* Check whether IPERM is a valid permutation */
	for (i = 1; i <= n; ++i) {
		indx = abs(iperm[i]);
		if (indx >= 1 && indx <= n) {
			if (iperm[indx] > 0) {
				iperm[indx] = -iperm[indx];
				continue;
			}
		}
		return 3;
	}
/* Rearrange the values of HX.  Use the IPERM vector as a flag.
 If IPERM(I) > 0, then the I-th value is in correct location */
	for (istrt = 1; istrt <= n; ++istrt) {
		if (iperm[istrt] > 0) continue;
		indx = istrt;
		indx0 = indx;
		strncpy(work, hx[istrt-1], hx_len);
L320:
		if (iperm[indx] >= 0) {
			strncpy(hx[indx0-1], work, hx_len);
		} else {
			strncpy(hx[indx-1], hx[(-iperm[indx])-1], hx_len);
			indx0 = indx;
			iperm[indx] = -iperm[indx];
			indx = iperm[indx];
			goto L320;
		}
	}
	return 0;
} /* end of hpperm */
