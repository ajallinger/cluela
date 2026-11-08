/*----------------------------------------------------------------------
mds.c - multi-dimensional scaling and principal component analysis
by Andy Allinger, 2013-2020, released to the public domain
This program may be used by any person for any purpose.
----------------------------------------------------------------------*/
#include <stdbool.h>
#include <math.h>
#include <float.h>

int seigv (float **, int, float *, float **, float *); /* defined in eigen.c */
void median (float *, int, float *); /* defined in median.c */
bool safdiv (float, float); /* defined in safdiv.c */

/*----------------------------------------------------------------------
dsmatx - Construct a dissimilarity matrix
Given a data set, find an embedding with the same distances.

___Name_______Type______In/Out____Description___________________________
   x[n][p]    float**   In        Data of P dimensions by N objects
   nx[n][p]   int**     In        1 indicates data is present for
                                    the correspoding X value
                                  0 indicates data is missing for
                                    the corresponding X value
   n          int       In        Number of objects
   p          int       In        Number of variables
   w[p]       float*    In        Importance weights
   cat[p]     int*      In        Number of categories each variable
                                  0 indicates a continuous variable
                                  1 indicates a dichotomous variable
                                    values {0.0, 1.0} for absent, present
                                  .GE. 2 indicates a categorical variable
                                         values {1.0, 2.0, ...}
                                  -1 indicates a cyclic variable,
                                     period of one
   robust     bool      In        Use Manhattan dissimilarity
   d[n][n]    float**   Out       Dissimilarity matrix
----------------------------------------------------------------------*/
int dsmatx (float **x, int **nx, int n, int p, float *w, int *cat, bool robust,
             float **d)
{
/* Constants */
	static float ctol = 0.001f;
	static float small = 1.e-36f;
/* Local variables */
	int i, j, k, nd;
	float dif, tot, dav, wsum;
	bool az;
/* Function Body */
	if (p < 1 || n < 1) return -2;
/* check for valid category labels */
	for (i = 0; i < n; ++i) {
		for (j = 0; j < p; ++j) {
			if (nx[i][j] > 0) {
				k = 0; /* avoid compiler warning */
				if (cat[j] > 0) {
					k = roundf(x[i][j]);
					if (fabsf(x[i][j] - (float) k) > ctol) return -3;
				}
				if (1 == cat[j]) {
					if (! (0 == k || 1 == k)) return -4;
				}
				if (cat[j] > 1) {
					if (k < 1 || k > cat[j]) return -5;
				}
			}
		}
	}
/* require a nonzero weight */
	az = true;
	for (j = 0; j < p; ++j) if (fabsf(w[j]) > small) az = false;
	if (az) return -6;
/* dissimilarity matrix */
	nd = 0;
	dav = 0.f;
	for (k = 0; k < n-1; ++k) {
		for (j = k+1; j < n; ++j) {
			tot = 0.f;
			wsum = 0.f;
			for (i = 0; i < p; ++i) {
				if (nx[j][i] <= 0 || nx[k][i] <= 0) continue;
				if (0 == cat[i]) {
					dif = x[j][i] - x[k][i];
				} else if (cat[i] > 1) {
					dif = (fabsf(x[j][i] - x[k][i]) > 0.5f) ? 1.f : 0.f;
				} else if (1 == cat[i]) {
					if (x[j][i] < 0.5f && x[k][i] < 0.5f) continue;
					dif = (x[j][i] < 0.5f || x[k][i] < 0.5f) ? 1.f : 0.f;
				} else if (-1 == cat[i]) {
					dif =  fmodf(fabsf(x[j][i] - x[k][i]), 1.f);
					dif = 2.f * fminf(dif, 1.f - dif);
				} else {
					return -7;
				}
				wsum += w[i];
				tot += robust ? w[i] * fabsf(dif) : w[i] * dif * dif;
			}
			if (safdiv(tot, wsum)) { /* adjust for missing values*/
				++nd;
				d[j][k] = tot / wsum;
				dav += d[j][k];
			} else {
				d[j][k] = -1.f; /* sentinel to be replaced */
			}
		}
	}
	if (0 == nd) return -8;
	dav /= (float) nd;
	for (k = 0; k < n; ++k) {
		d[k][k] = 0.f; /* zero on the diagonal */
		for (j = k+1; j < n; ++j) {
			if (d[j][k] < 0.f) d[j][k] = dav;
			d[k][j] = d[j][k];
		}
	}
	return 0;
} /* end of dsmatx */

/*----------------------------------------------------------------------
mds - Multi-Dimensional Scaling
Given a data set, find an embedding with the same distances.

___Name_______Type______In/Out____Description___________________________
   x[n][p]    float**   In        Data of P dimensions by N objects
   nx[n][p]   int**     In        1 indicates data is present for
                                    the correspoding X value
                                  0 indicates data is missing for
                                    the corresponding X value
   n          int       In        Number of objects
   p          int       In        Number of variables
   w[p]       float*    In        Importance weights
   cat[p]     int*      In        Number of categories each variable
                                  0 indicates a continuous variable
                                  1 indicates a dichotomous variable
                                    values {0.0, 1.0} for absent, present
                                  .GE. 2 indicates a categorical variable
                                         values {1.0, 2.0, ...}
                                  -1 indicates a cyclic variable,
                                     period of one
   robust     bool      In        Use Manhattan dissimilarity
   v[n][n]    float**   Out       Embedded data matrix, N by M
   m          int*      Out       Number of directions found.  M .LE. N
   e[n]       float*    Out       Eigenvalues
   t[n]       float*    Neither   Work array

___Return values________________________________________________________
    0   no errors
   +n   eigenvectors incomplete
   -1   error in eigensolver
   -2   N or P less than 1
   -3   fractional category value
   -4   category value out of bounds, dichotomous
   -5   category value out of bounds, nominal
   -6   all zero weights
   -7   bad value in CAT()
   -8   no comparisons possible
----------------------------------------------------------------------*/
int mds (float **x, int **nx, int n, int p, float *w, int *cat, bool robust,
         float **v, int *m, float *e, float *t)
{
/* Local variables */
	int j, k, ifault;
	float s, tot;
/* Function Body */
	ifault = dsmatx(x, nx, n, p, w, cat, robust, v);
	if (ifault != 0) return ifault;
	for (k = 0; k < n; ++k) { /* column sums */
		t[k] = 0.f;
		for (j = 0; j < n; ++j) t[k] += v[j][k];
	}
	tot = 0.f; /* overall total */
	s = 1.f / (float) (n);
	for (j = 0; j < n ; ++j) {
		t[j] *= s;
		tot += t[j];
	}
	tot *= s;
/*----------------------------------------------------------------------
   center the matrix

                      1                1                 1
-2 C[i,j] = A[i,j] - --- SUM A[i,k] - --- SUM A[k,j] + ----- SUM SUM A[k,l]
                      N   k            N   k            N^2   k   l

----------------------------------------------------------------------*/
	for (k = 0; k < n ; ++k) {
		for (j = k; j < n; ++j) {
			v[j][k] = -0.5f * (v[j][k] - t[j] - t[k] + tot);
		}
	}
	ifault = seigv (v, n, e, v, t); /* factor */
	if (ifault != 0) return ifault;
	for (k = 0; k < n/2; ++k) { /* decreasing order */
		s = e[k];
		e[k] = e[n-k-1];
		e[n-k-1] = s;
		for (j = 0; j < n; ++j) {
			s = v[j][k];
			v[j][k] = v[j][n-k-1];
			v[j][n-k-1] = s;
		}
	}
	for (k = 0; k < n; ++k) { /* embed */
		if (e[k] > 0.f) *m = k + 1;
		s = sqrtf(fabsf(e[k]));
		for (j = 0; j < n; ++j) v[j][k] *= s;
	}
	return 0;
} /* end of mds */

/*----------------------------------------------------------------------
pca - Principal Component Analysis.  No support for missing data.

___Name______Type______In/Out____Description____________________________
   x[n][p]   float**   In        Data matrix
   n         int       In        Number of objects
   p         int       In        Variables each object
   robust    bool      In        Use comedian instead of covariance
   r[n][m]   float**   Out       Projected data
   m         int       In        # of dimensions requested to project
   q         int       Out       # of dimensions found
   a[p]      float*    Out       Averages each variable
   e[p]      float*    Out       Eigenvalues of covariance matrix
   v[p][p]   float**   Out       Eigenvectors of covariance matrix
   t[p]      float*    Neither   Workspace
   u[n]      float*    Neither   Workspace
----------------------------------------------------------------------*/
int pca (float **x, int n, int p, bool robust, float **r, int m, int *q,
     float *a, float *e, float **v, float *t, float *u)
{
/* Local variables */
	int i, j, k, k1, ifault;
	float tot, xj, xk;
/* Function Body */
	if (m < 1 || m > p) return -1;
	if (p < 1 || n < 2) return -2; /* need two for variance */
/* Compute averages. */
	for (j = 0; j < p; ++j) {
		if (robust) {
			for (i = 0; i < n; ++i) u[i] = x[i][j];
			median (u, n, &a[j]);
		} else {
			tot = 0.f;
			for (i = 0; i < n; ++i) tot += x[i][j];
			a[j] = tot / (float) n;
		}
	}
/*----------------------------------------------------------------------
Covariance / Comedian matrix.

The comedian matrix is a highly robust estimate of covariance
that is easy to compute.  Refer to:
     Michael Falk
     "On MAD and Comedians"
     Ann. Inst. Statist. Math., v.49 n.4 pp.615-644, 1997

The covariance matrix generalizes the variance.
The comedian matrix generalizes the median absolute deviation.
----------------------------------------------------------------------*/
	for (k = 0; k < p; ++k) {
		for (j = k; j < p; ++j) {
			for (i =0; i < n; ++i) {
				xj = x[i][j] - a[j];
				xk = x[i][k] - a[k];
				u[i] = xj * xk;
			}
			if (robust) {
				median (u, n, &v[j][k]);
			} else {
				tot = 0.f;
				for (i = 0; i < n; ++i) tot += u[i];
				v[j][k] = tot / (float) (n - 1);
			}
		}
	}
/* Factor covariance matrix */
	ifault = seigv (v, p, e, v, t);
	if (0 != ifault) return ifault;
/*----------------------------------------------------------------------
Step through columns of V in reverse order so that the first column
of R has the greatest variance.
        R[i,k] = X[j,i] V[j,p-k+1]
----------------------------------------------------------------------*/
	*q = 0;
	for (k = 0; k < m; ++k) {
		k1 = p - k - 1;
		if (e[k1] <= 0.f) return 0;
		*q = k + 1;
		for (i = 0; i < n; ++i) {
	    		r[i][k] = 0.f;
			for (j = 0; j < p; ++j) r[i][k] += (x[i][j] - a[j]) * v[j][k1];
	    }
	}
	return 0;
} /* end of pca */
