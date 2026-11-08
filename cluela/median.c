/*----------------------------------------------------------------------
median.f90 -> median.c
by Alan Miller, placed in the public domain
translated to C, by AJA 23 April 2017
for the original, see http://www.jblevins.org/mirror/amiller/
------------------------------------------------------------------------

___Name___Type_____In/Out___Description_________________________________
   x[n]   float    In       Data array
   n      int      In       Length of array
   xmed   float*   Out      The median
----------------------------------------------------------------------*/
#include <math.h>
#include <stdbool.h>
void median (float *x, int n, float *xmed)
{
/* Local variables */
	int i, j, k, hi, lo;
	int mid, nby2, nby2p1;
	bool odd;
	float xhi, xlo, temp, xmin, xmax;

/* Find the median of x[1], ... , x[n], using as much of the quicksort
 algorithm as is needed to isolate it.
 N.B. On exit, the array x is partially ordered.
     Latest revision - 26 November 1996 */
	nby2 = n / 2;
	nby2p1 = nby2 + 1;
	odd = true;

/* hi & lo are position limits encompassing the median. */
	if (n == nby2 << 1) odd = false;
	lo = 1;
	hi = n;
	if (n < 3) {
		if (n < 1) {
			*xmed = 0.f;
			return;
		}
		*xmed = x[0];
		if (1 == n) return;
		*xmed = (*xmed + x[1]) * .5f;
		return;
	}

/* Find median of 1st, middle & last values. */
L10:
	mid = (lo + hi) / 2;
	*xmed = x[mid-1];
	xlo = x[lo-1];
	xhi = x[hi-1];

/* Swap xhi & xlo */
	if (xhi < xlo) {
		temp = xhi;
		xhi = xlo;
		xlo = temp;
	}
	if (*xmed > xhi) {
		*xmed = xhi;
	} else if (*xmed < xlo) {
		*xmed = xlo;
	}

/* The basic quicksort algorithm to move all values <= the sort key (xmed)
 to the left-hand end, and all higher values to the other end. */
	i = lo;
	j = hi;
L50:
	while (true) {
		if (x[i-1] >= *xmed) break;
		++i;
	}
	while (true) {
		if (x[j-1] <= *xmed) break;
		--j;
	}
	if (i < j) {
		temp = x[i-1];
		x[i-1] = x[j-1];
		x[j-1] = temp;
		++i;
		--j;

/* Decide which half the median is in. */
		if (i <= j) goto L50;
	}

	if (! odd) {
		if (j == nby2 && i == nby2p1) goto L130;
		if (j < nby2) lo = i;
		if (i > nby2p1) hi = j;
		if (i != j) goto L100;
		if (i == nby2) lo = nby2;
		if (j == nby2p1) hi = nby2p1;
	} else {
		if (j < nby2p1) lo = i;
		if (i > nby2p1) hi = j;
		if (i != j) goto L100;

/* Test whether median has been isolated. */
		if (i == nby2p1) return;
	}
L100:
	if (lo < hi - 1) goto L10;
	if (! odd) {
		*xmed = (x[nby2-1] + x[nby2p1-1]) * .5f;
		return;
	}
	temp = x[lo-1];
	if (temp > x[hi-1]) {
		x[lo-1] = x[hi-1];
		x[hi-1] = temp;
	}
	*xmed = x[nby2p1-1];
	return;

/* Special case, n even, j = n/2 & i = j + 1, so the median is
 between the two halves of the series.   Find max. of the first
 half & min. of the second half, then average. */
L130:
	xmax = x[0];
	for (k = lo; k <= j; ++k) xmax = fmaxf(xmax, x[k-1]);
	xmin = x[n-1];
	for (k = i; k <= hi; ++k) xmin = fminf(xmin, x[k-1]);
	*xmed = (xmin + xmax) * .5f;
	return;
} /* end of median */
