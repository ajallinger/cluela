/***********************************************************************
 valid.c - measures for clustering validation
   by Andy Allinger, 2014-2017, released to the public domain
   This program may be used by any person for any purpose.
***********************************************************************/
#include <math.h>

/***********************************************************************
 ARIND - adjusted Rand index
  refer to:
   "Comparing Partitions"
   Lawrence Hubert and Phipps Arabie
   Journal of Classification v.2 pp.193-218, 1985.

___Name______Type______In/Out____Description____________________________
   y[n]      int*      In        Correct class memberships
   z[n]      int*      In        Computed cluster labels
   n         int       In        Length of Y and Z
   g         int       In        Number of correct classes
   k         int       In        Number of discovered clusters
   m[g][k]   float**   Neither   Confusion matrix
   a[g]      float*    Neither   Class marginals
   b[k]      float*    Neither   Cluster marginals
   arind     float     Out       Adjusted Rand index
***********************************************************************/
float arind (int *y, int *z, int n, int g, int k,
            float **m, float *a, float *b)
{
/* Local variables */
	int i, j;
	long l;
	float as, bs, ri, adj, acc;
/* Function Body */
/* validate */
	acc = 0.f;
	for (i = 0; i < n; ++i) {
		if (y[i] < 1 || y[i] > g || z[i] < 1 || z[i] > k) return 0.f;
	}
/* clear */
	for (j = 0; j < k; ++j) {
		for (i = 0; i < g; ++i) {
			m[i][j] = 0.f;
		}
	}
/* count contingency table */
	for (i = 0; i < n; ++i) {
		j = y[i] - 1;
		l = z[i] - 1;
		m[j][l] += 1.f;
	}
/* row sums */
	for (i = 0; i < g; ++i) {
		a[i] = 0.f;
		for (j = 0; j < k; ++j) a[i] += m[i][j];
	}
/* column sums */
	for (j = 0; j < k; ++j) {
		b[j] = 0.f;
		for (i = 0; i < g; ++i) b[j] += m[i][j];
	}
/* main term */
	ri = 0.f;
	for (j = 0; j < k; ++j) {
		for (i = 0; i < g; ++i) {
			l = lrintf(m[i][j]);
			ri += l * (l-1) / 2;
		}
	}
/* adjust term */
	as = 0.f;
	for (i = 0; i < g; ++i) {
		l = lrintf(a[i]);
		as += l * (l-1) / 2;
	}
	bs = 0.f;
	for (j = 0; j < k; ++j) {
		l = lrintf(b[j]);
		bs += l * (l-1) / 2;
	}
	adj = as * bs / (n * (n-1) / 2);
/* finish */
	acc = (ri - adj) / ((as + bs) * .5f - adj);
	return acc;
} /* end of arind */

/***********************************************************************
 VINFO - variation of information
  refer to:
   "Comaparing Clusterings"
   Marina Meila
   University of Washington Technical report, 2002

___Name______Type______In/Out____Description____________________________
   y[n]      int*      In        Correct class memberships
   z[n]      int*      In        Computed cluster labels
   n         int       In        Length of Y and Z
   g         int       In        Number of correct classes
   k         int       In        Number of discovered clusters
   m[g][k]   float**   Neither   Confusion matrix
   a[g]      float*    Neither   Class marginals
   b[k]      float*    Neither   Cluster marginals
   vinfo     float     Out       Variation of information
***********************************************************************/
float vinfo (int *y, int *z, int n, int g, int k,
            float **m, float *a, float *b)
{
/* Initialized data */
	static float small = 1e-9f;
/* Local variables */
	int i, j, l;
	float ha, hb, imut, acc;
/* Function Body */
/* validate */
	acc = 0.f;
	for (i = 0; i < n; ++i) {
		if (y[i] < 1 || y[i] > g || z[i] < 1 || z[i] > k) return 0.f;
	}
/* clear */
	for (j = 0; j < k; ++j) {
		for (i = 0; i < g; ++i) {
			m[i][j] = 0.f;
		}
	}
/* count */
	for (i = 0; i < n; ++i) {
		j = y[i] - 1;
		l = z[i] - 1;
		m[j][l] += 1.f;
	}
/* row sums */
	for (i = 0; i < g; ++i) {
		a[i] = 0.f;
		for (j = 0; j < k; ++j) a[i] += m[i][j];
/* make a probability */
		a[i] /= (float) (n);
	}
/* column sums */
	for (j = 0; j < k; ++j) {
		b[j] = 0.f;
		for (i = 0; i < g; ++i) b[j] += m[i][j];
/* make a probability */
		b[j] /= (float) (n);
	}
/* change counts to a joint probability */
	for (j = 0; j < k; ++j) {
		for (i = 0; i < g; ++i) {
			m[i][j] /= (float) (n);
		}
	}
/* entropy A */
	ha = 0.f;
	for (i = 0; i < g; ++i) {
		if (a[i] < small) continue;
		ha -= a[i] * logf(a[i]);
	}
/* entropy B */
	hb = 0.f;
	for (j = 0; j < k; ++j) {
		if (b[j] < small) continue;
		hb -= b[j] * logf(b[j]);
	}
/* mutual information */
	imut = 0.f;
	for (j = 0; j < k; ++j) {
		for (i = 0; i < g; ++i) {
			if (m[i][j] < small) continue;
			imut += m[i][j] * logf(m[i][j] / (a[i] * b[j]));
		}
	}
/* variation of information */
	acc = ha + hb - imut * 2.f;
/* convert to bits */
	acc /= (float)log(2.f);
	return acc;
} /* end of vinfo */
