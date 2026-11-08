/*----------------------------------------------------------------------
   cluste.c - CLUstering with Superior Theoretical Efficiency
   by Andy Allinger, 2017, released to the public domain
   This program may be used by any person for any purpose.
----------------------------------------------------------------------*/
#include <stdbool.h>

void qsorti (int *, int *, int); /* defined in sort.c */
void iorder (int *, int *, int);
bool almeq (float, float); /* defined in cluela.c */
int hat (int);
float dsm (float *, int *, float *, int *, int, float *, bool, int);
void cent (float **, int **, int, int, float **, int **, int, int *, int *,
       bool, int *, int *, float *);
int kclust (float **, int **, int, int, float *, float **, int **,
      int, int *, int, int *, int *, bool, bool, bool, float, int *, int *,
      int *, int **, float **, float *, float *, float *, float *);

/*======================================================================
hclust - Hyperplane partitioning
Divisive clustering based on partial sorting, a reduced form of
Anthony Hoare's quick sort algorithm.

___Name_______Type______I/O____Description______________________________
   x[n][p]    float**   In     Data array
   nx[n][p]   int**     In     1 => data present for this measurement
   n          int       In     Number of data points
   p          int       In     Dimension of data
   w[p]       float*    In     Importance weights
   k          int       In     Desired number of clusters
   z[n]       int*      Out    What cluster each point is in
   pop[k]     int*      Out    Cluster populations
   robust     bool      In     Use L1 distance?
   ind[n]     int*      None   Index array that sorts X into
                                    output order from input order
   jnd[k]     int*      None   Sorting workspace
   bound[k]   int*      None   Work array
======================================================================*/
void hclust (float **x, int **nx, int n, int p, float *w, int k, int *z,
       int *pop, bool robust, int *ind, int *jnd, int *bound)
{
/* Constants */
	static int ZERO = 0;
/* Local variables */
	int i, l, hi, ll, lo, ihi, khi, ilo, klo, kmid, iswap, maxpop;
	float dhl, dhm, dml;
/* Function Body */
	for (i = 0; i < n; ++i) ind[i] = i;
	for (l = 0; l < k; ++l) bound[l] = n-1;
	pop[0] = n;
	ll = 0;

/*======================================================================
   Partition a section of the array such that elements closer to the first
   element are moved to the front, and elements closer to the last element
   are moved to the rear.  Continue until enough clusters have been made.
======================================================================*/
	for (l = 1; l < k; ++l) {
/* find partition with greatest population */
		maxpop = 0;
		for (i = 0; i < l; ++i) {
			if (pop[i] > maxpop) {
				maxpop = pop[i];
				ll = i;
			}
		}
/* find limits of partition */
		hi = bound[ll];
		lo = hi - maxpop + 1;
		ihi = hi - 1;
		ilo = lo + 1;
		if (2 == maxpop) goto L70;
/* swap random pivot values to ends of partition (introduce chance) */
		khi = hat(maxpop) - 1 + lo;
L20:
		klo = hat(maxpop) - 1 + lo;
		if (klo == khi) goto L20;
/* a third pivot value in case first two are near each other */
L30:
		kmid = hat(maxpop) - 1 + lo;
		if (kmid == khi || kmid == klo) goto L30;
		dhl = dsm (x[ind[khi]], nx[ind[khi]], x[ind[klo]], nx[ind[klo]], p, w, robust, ZERO);
		dhm = dsm (x[ind[khi]], nx[ind[khi]], x[ind[kmid]], nx[ind[kmid]], p, w, robust, ZERO);
		dml = dsm (x[ind[kmid]], nx[ind[kmid]], x[ind[klo]], nx[ind[klo]], p, w, robust, ZERO);
		if (dhm > dhl || dml > dhl) {
			if (dhm > dml) {
				klo = kmid;
			} else {
				khi = kmid;
			}
		}
		iswap = ind[khi];
		ind[khi] = ind[hi];
		ind[hi] = iswap;
		iswap = ind[klo];
		ind[klo] = ind[lo];
		ind[lo] = iswap;

/*======================================================================
           sort
======================================================================*/
L50: /* look for object that belongs in HI side */
		if (ilo <= ihi) {
			if (dsm(x[ind[ilo]], nx[ind[ilo]], x[ind[lo]], nx[ind[lo]], p, w, robust, ZERO)
			  <= dsm(x[ind[ilo]], nx[ind[ilo]], x[ind[hi]], nx[ind[hi]], p, w, robust, ZERO)) {
				++ilo;
				goto L50;
			}
		}
L60: /* look for object that belongs in LO side */
		if (ihi >= ilo) {
			if (dsm(x[ind[ihi]], nx[ind[ihi]], x[ind[hi]], nx[ind[hi]], p, w, robust, ZERO)
			  <= dsm(x[ind[ihi]], nx[ind[ihi]], x[ind[lo]], nx[ind[lo]], p, w, robust, ZERO)) {
				--ihi;
				goto L60;
			}
		}
/* swap these objects */
		if (ilo < ihi) {
			iswap = ind[ilo];
			ind[ilo] = ind[ihi];
			ind[ihi] = iswap;
			++ilo;
			--ihi;
			goto L50;
		}
L70: /* Add this partition boundary to the bounds array */
		if (ilo == hi) { /* equal split if all distances were the same */
			ihi = (lo + hi) / 2;
			ilo = ihi + 1;
		}
		pop[ll] = hi - ihi; /* reduce population */
		pop[l] = ilo - lo; /* population in new partition */
		bound[l] = ihi; /* append */
	}

/*======================================================================
   Recover memberships and populations from the index array
======================================================================*/
	qsorti (bound, jnd, k);
	iorder (bound, jnd, k);
	l = 0;
	for (i = 0; i < n; ++i) {
		if (i > bound[l]) ++l;
		z[ind[i]] = l + 1;
	}
	i = -1;
	for (l = 0; l < k; ++l) {
		pop[l] = bound[l] - i;
		i = bound[l];
	}
	return;
} /* end of hclust */

/*----------------------------------------------------------------------
Speedy clustering with theoretical efficiency faster than O(iNPK)
if there is NO MISSING DATA, NO SUM-OF-SQUARES MINIMIZATION, and
NO NEGATIVE WEIGHTS.  These conditions are not checked!  It is OK
to call the function anyway, but efficiency will be worse.

___Variable________Type______I/O_______Description______________________
   x[n][p]         float**   In        Data points in R^P
   nx[n][p]        int**     In        1 indicates data is present for
                                           the corresponding X value.
   n               int       In        Number of points
   p               int       In        Measurements per point
   w[p]            float*    In        Importance weights of variables
   c[kmax][p]      float**   Out       Cluster centroids
   nc[kmax][p]     int**     Out       Count of centroid data available
   kmin            int       In        Fewest clusters to make
   k               int*      Both      Default/actual clusters to make
   kmax            int       In        Most clusters to make
   z[n]            int*      Out       What cluster each point is in
   pop[kmax]       int*      Out       Number of points in each cluster
   robust          bool      In        Use L1 statistics
   sosdo           bool      In        Sum of squares distance objective
   iwork[2*n+kmax]
                   int*      Neither   Workspace
   rwork[2*n+kmax]
                   float*    Neither   Workspace
   tree[2*k][12]   int**     Neither   Workspace
   crit[2*k][4]    float**   Neither   Workspace
   u               float*    Out       Residual error
----------------------------------------------------------------------*/
int cluste (float **x, int **nx, int n, int p, float *w, float **c, int **nc,
      int kmin, int *k, int kmax, int *z, int *pop, bool robust,
      bool sosdo, int *iwork, float *rwork, int **tree, float **crit, float *u)
{
/* Constants */
	static int ZERO = 0;
	static int FZERO = 0.f;
/* Local variables */
	int i, j, k0, l, ifault;
	float d, thresh, wall;
	bool miss;
/* Function Body */
	k0 = *k;
	if (kmin < 1 || k0 < kmin || k0 > kmax || kmax > n) return 6;
	if (robust && sosdo) return 7;
	if (p < 1 || n < 1) return 8; /* non-null data array? */
/* at least one element of W() must be nonzero */
	for (j = 0; j < p; ++j) if (! almeq (w[j], FZERO)) goto L10;
	return 9;
L10: /* check for data with no values present */
	for (i = 0; i < n; ++i) {
		for (j = 0; j < p; ++j) {
			if (nx[i][j] > 0) goto L20;
		}
		return 10;
L20:
		;
	}
/* look for missing data */
	miss = false;
	for (i = 0; i < n; ++i) {
		for (j = 0; j < p; ++j) {
			if (nx[i][j] <= 0) {
				miss = true;
				goto L30;
			}
		}
	}
L30:
	wall = FZERO;
	for (j = 0; j < p; ++j) wall += w[j];
	if (almeq(wall, FZERO) && !miss) return 11;
/* Choose initial centers */
	hclust (x, nx, n, p, w, k0, z, pop, robust, iwork, &iwork[n], &iwork[n+k0]);
	for (l = 0; l < k0; ++l) iwork[l] = 1;
	cent (x, nx, n, p, c, nc, k0, z, pop, robust, iwork, &iwork[k0], rwork);
	thresh = 0.f;
	for (i = 0; i < n; ++i) {
		l = z[i] - 1;
		d = dsm (x[i], nx[i], c[l], nc[l], p, w, robust, ZERO);
		if (d > thresh) thresh = d;
	}
/* cluster */
	ifault = kclust (x, nx, n, p, w, c, nc, kmin, k, kmax, z, pop,
	  miss, robust, sosdo, thresh, iwork, &iwork[kmax], &iwork[n+kmax], tree, crit,
	  rwork, &rwork[kmax], &rwork[n+kmax], u);
	return ifault;
} /* end of cluste */
