/*-----------------------------------------------------------------------
prep.c - Routines for preprocessing data
by Andy Allinger, 2011-2017, released to the public domain
This program may be used by any person for any purpose.
-----------------------------------------------------------------------*/
#include <stdbool.h>
#include <string.h>
#include <math.h>

void qsorti (int *, int *, int); /* defined in sort.c */
void qsortr (float *, int *, int);
int hpsort (char **, int, int, int, int, int *, int, char *);
void iorder (int *, int *, int);
void rorder (float *, int *, int);
int hpperm (char **, int, int, int *, char *);
bool almeq (float, float); /* defined in cluela.c */

/*----------------------------------------------------------------------
Sort non-missing data to the front, by a single iteration of quick sort.

___Name_____Type___In/Out____Description________________________________
   nx[n]    int    In*       Input array, 0 = missing, 1 = present
   ind[n]   int    Out*      Permutation vector
   n        int    In        Length of array
----------------------------------------------------------------------*/
void qsortm (int *nx, int *ind, int n)
{
/* Local variables */
	int i, j, k;
/* Parameter adjustments */
	--nx;
	--ind;
/* Function Body */
	for (i = 1; i <= n; ++i) ind[i] = i;
/* validate */
	if (n < 2) return;
	i = 1;
	j = n;
/* Find a missing datum at the front of the array */
L50:
	if (i < j) {
		if (nx[ind[i]] > 0) {
			++i;
			goto L50;
		}
	}
/* Find a non-missing datum at the rear of the array */
L60:
	if (i < j) {
		if (nx[ind[j]] <= 0) {
			--j;
			goto L60;
		}
	}
/* Swap these elements */
	if (i < j) {
		k = ind[i];
		ind[i] = ind[j];
		ind[j] = k;
		++i;
		--j;
		goto L50;
	}
	return;
} /* end of qsortm */

/*----------------------------------------------------------------------
jscode - succint code:  change arbitrary integers to consecutive integers

___Name______Type___In/Out____Description_______________________________
   jnam[n]   int*   Both      Data array
   nx[n]     int*   In        1 if data exists, each datum
   ind[n]    int*   Neither   Permutation vector
   jnd[n]    int*   Neither   Permutation vector
   n         int    In        Length of JNAM
----------------------------------------------------------------------*/
void jscode (int *jnam, int *nx, int *ind, int *jnd, int n)
{
/* Local variables */
	int i, j, k, l, m;
/* Parameter adjustments */
	--jnam;
	--nx;
	--ind;
	--jnd;
/* Function Body */
/* Begin. Move missing data to end. */
	qsortm (&nx[1], &jnd[1], n);
/* count data available */
	m = 0;
	for (i = 1; i <= n; ++i) if (nx[i] > 0) ++m;
	if (0 == m) return;
	if (1 == m) {
		jnam[jnd[1]] = 1;
		return;
	}
	iorder (&jnam[1], &jnd[1], n);
/* Sort implicitly. */
	qsorti (&jnam[1], &ind[1], m);
/* Compare adjacent entries. */
	k = jnam[ind[1]];
	l = 1;
	jnam[ind[1]] = 1;
	for (i = 2; i <= m; ++i) {
		j = jnam[ind[i]];
		if (j != k) {
			k = j;
			++l;
		}
		jnam[ind[i]] = l;
	}
/* Restore original order. */
	for (i = 1; i <= n; ++i) ind[jnd[i]] = i;
	iorder (&jnam[1], &ind[1], n);
	return;
} /* end of jscode */

/*----------------------------------------------------------------------
ireid - replace text names with integers 1, 2, 3... for each distinct item

___Name_____Type_____I/O_______Description______________________________
   h[n]     char**   In        Items or names
   id[n]    int*     Out       Successive small integers
   nx[n]    int*     In        1 if datum exists
   ind[n]   int*     Neither   Permutation vector
   jnd[n]   int*     Neither   Permutation vector
   n        int      In        Length of array
   h_len    int      In        Length of each entry in H
   hwork    char*    Neither   Used by HPSORT, same length as H
----------------------------------------------------------------------*/
int ireid (char **h, int *id, int *nx, int *ind, int *jnd, int n,
           int h_len, char *hwork)
{
/* Constants */
	static int ONE = 1;
	static int TWO = 2;
/* Local variables */
	int i, j, l, m, ierr;
/* Parameter adjustments */
	--id;
	--nx;
	--ind;
	--jnd;
/* Function Body */
	ierr = 0;
/* move missing data to end */
	qsortm (&nx[1], &jnd[1], n);
/* count data available */
	m = 0;
	for (i = 1; i <= n; ++i) if (nx[i] > 0) ++m;
	if (0 == m) return 0;
	if (1 == m) {
		id[ind[1]] = 1;
		return 0;
	}
/* apply permutation */
	ierr = hpperm (h, n, h_len, &jnd[1], hwork);
	if (ierr != 0) return ierr;
/* Sort to bring duplicate values adjacent. */
	ierr = hpsort (h, m, h_len, ONE, h_len, &ind[1], TWO, hwork);
	if (ierr != 0) return ierr;
/* apply permutation */
	iorder (&jnd[1], &ind[1], m);
/* compare adjacent strings */
	j = 1;
	l = 1;
	id[j] = l;
	for (i = 2; i <= m; ++i) {
		if (strncmp(h[i-1], h[j-1], h_len)) { /* strings differ */
			j = i;
			++l;
		}
		id[i] = l;
	}
/* restore original order */
	for (i = 1; i <= n; ++i) ind[jnd[i]] = i;
	ierr = hpperm (h, n, h_len, &ind[1], hwork);
	iorder (&id[1], &ind[1], n);
	if (ierr != 0) return ierr;
/* Change labels to be the order of first occurance in the input array */
	for (i = 1; i <= n; ++i) ind[i] = 0;
	l = 0;
/* form replacement table in IND */
	for (i = 1; i <= n; ++i) {
		if (nx[i] > 0) {
			if (0 == ind[id[i]]) {
				++l;
				ind[id[i]] = l;
			}
		}
	}
/* replace */
	for (i = 1; i <= n; ++i) if (nx[i] > 0) id[i] = ind[id[i]];
	return 0;
} /* end of ireid */

/*----------------------------------------------------------------------
ordnal - replace data with its ordinal rank, rescaled to 0...1

___Name_____Type_____I/O_______Description______________________________
   x[n]     float*   Both      Data vector
   nx[n]    int*     In        1 if data exists
   ind[n]   int*     Neither   Permutation vector
   jnd[n]   int*     Neither   Permutation vector
   n        int      In        Length of vector
----------------------------------------------------------------------*/
void ordnal (float *x, int *nx, int *ind, int *jnd, int n)
{
/* Local variables */
	int i, j, k, m;
	float v;
/* Parameter adjustments */
	--x;
	--nx;
	--ind;
	--jnd;
/* Function Body */
/* move missing data to end */
	qsortm (&nx[1], &jnd[1], n);
/* count data available */
	m = 0;
	for (i = 1; i <= n; ++i) if (nx[i] > 0) ++m;
	if (0 == m) return;
	if (1 == m) {
		x[ind[1]] = 0.f;
		return;
	}
/* apply permutation */
	rorder (&x[1], &jnd[1], n);
/* bring duplicate values adjacent */
	qsortr (&x[1], &ind[1], m);
/* apply permutation */
	rorder (&x[1], &ind[1], m);
/* apply permutation */
	iorder (&jnd[1], &ind[1], m);
/* compare adjacent values */
	j = 1;
	for (i = 2; i <= m; ++i) {
		if (almeq (x[i], x[j])) continue;
		v = (float) (j + i - 1) * .5f - 1.f;
/* replace */
		for (k = j; k <= i-1; ++k) x[k] = v;
		j = i;
	}
/* last stretch */
	v = (float) (j + m) * .5f - 1.f;
	for (k = j; k <= m; ++k) x[k] = v;
/* rescale 0...1 */
	v = 1.f / (float) (m - 1);
	for (i = 1; i <= m; ++i) x[i] *= v;
/* restore original order */
	for (i = 1; i <= n; ++i) ind[jnd[i]] = i;
	rorder (&x[1], &ind[1], n);
	return;
} /* end of ordnal */

/*----------------------------------------------------------------------
Construct vertices of a regular simplex according to Manhattan metric

___Name________Type______I/O___Description______________________________
   m           int       In    Number of vertices
   s[m][m-1]   float**   Out   Coordinates
----------------------------------------------------------------------*/
/* regular simplex coordinates, Lebesgue 1 */
void rscl1 (int m, float **s)
{
/* Local variables */
	int i, j, k;
/* Function Body */
	if (m < 2) return;
	s[0][0] = 0.f;
	s[1][0] = 1.f;
/* each vertex */
	for (k = 3; k <= m; ++k) {
/* first entry always 1/2 */
		s[k-1][0] = .5f;
/* zeros */
		for (i = 2; i <= k-2; ++i) s[k-1][i-1] = 0.f;
/* fill in zeros */
		for (j = 1; j <= k-1; ++j) s[j-1][k-2] = 0.f;
/* last entry always 1/2 */
		s[k-1][k-2] = .5f;
	}
	return;
} /* end of rscl1 */

/*----------------------------------------------------------------------
Construct vertices of a regular simplex according to Euclidean metric

___Name________Type______I/O___Description______________________________
   m           int       In    Number of vertices
   s[m][m-1]   float**   Out   Coordinates
----------------------------------------------------------------------*/
/* regular simplex coordinates, Lebesgue 2 */
void rscl2 (int m, float **s)
{
/* Local variables */
	int i, j, k;
	float r, t;
/* Function Body */
	if (m < 2) return;
	s[0][0] = 0.f;
	s[1][0] = 1.f;
/* each vertex */
	for (k = 3; k <= m; ++k) {
/* average each row */
		for (i = 1; i <= k-2; ++i) {
			t = 0.f;
			for (j = 1; j <= k-1; ++j) t += s[j-1][i-1];
			s[k-1][i-1] = t / (float) (k-1);
		}
/* fill in zeros */
		for (j = 1; j <= k-1; ++j) s[j-1][k-2] = 0.f;
/* make unit length */
		t = 0.f;
		for (i = 1; i <= k-2; ++i) {
			r = s[k-1][i-1];
			t += r * r;
		}
		s[k-1][k-2] = sqrtf(1.f - t);
	}
	return;
} /* end of rscl2 */

/*======================================================================
rscode - Recode categorical variables as column of regular simplex
Also scale down importance weights to preserve original importances.

___Name_______Type______I/O_______Description___________________________
   x[n][p]    float**   Both      Data array
   nx[n][p]   int**     Both      1 if data is present for this measurement
   n          int       In        Number of objects
   p          int       In        Number of coded variables
   w[p]       float*    Both      Importance weights
   cat[o]     int*      In        Number of categories each variable
   o          int       In        Number of input variables
   s[catmax][catmax-1]
              float**   Neither   Simplex coordinates
   catmax     int       In        Maximum number of categories
   robust     bool      In        Use L1 statistics?
=====================================================================*/
int rscode (float **x, int **nx, int n, int p, float *w,
            int *cat, int o, float **s, int catmax, bool robust)
{
/* Local variables */
	int h, i, j, l, clm1;
	long k;
	float fclm1;
/* Function Body */
/* Validate */
	j = 0;
	k = 0;
	for (l = 0; l < o; ++l) {
		if (1 == cat[l]) return 1;
		if (cat[l] < 0) return -1;
		j = j + 1 + (cat[l]-2 > 0 ? cat[l]-2 : 0);
		if (cat[l] > k) k = cat[l];
	}
	if (j != p) return 2;
	if (k > catmax) return 3;

/*======================================================================
        change weights
======================================================================*/
	j = p;
	for (l = o-1; l >= 0; --l) {
		if (0 == cat[l]) {
			--j;
			w[j] = w[l];
		} else {
			clm1 = cat[l] - 1;
			fclm1 = (float) clm1;
			for (h = clm1-1; h >= 0; --h) {
				--j;
				w[j] = w[l] / fclm1;
			}
		}
	}

/*======================================================================
              Recode.
======================================================================*/
/* construct simplex */
	if (robust) {
		rscl1 (catmax, s);
	} else {
		rscl2 (catmax, s);
	}
	for (i = n-1; i >= 0; --i) {
		j = p;
		for (l = o-1; l >= 0; --l) {
			if (0 == cat[l]) {
				--j;
				if (nx[i][l] > 0) {
					x[i][j] = x[i][l];
					nx[i][j] = 1;
				} else {
					x[i][j] = -999.f;
					nx[i][j] = 0;
				}
			} else {
				if (nx[i][l] > 0) {
					k = lrintf(x[i][l]);
					if (k < 1 || k > cat[l]) return 4;
				}
				clm1 = cat[l] - 1;
				fclm1 = (float) clm1;
				if (!robust) fclm1 = sqrtf(fclm1);
				for (h = clm1-1; h >= 0; --h) {
					--j;
					if (nx[i][l] > 0) {
						x[i][j] = s[k-1][h] * fclm1;
						nx[i][j] = 1;
					} else {
						x[i][j] = -999.f;
						nx[i][j] = 0;
					}
				}
			}
		}
	}
	return 0;
} /* end of rscode */
