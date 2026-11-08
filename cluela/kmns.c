/*----------------------------------------------------------------------
kmns.c - a translation of kmns.f
Use and distribution is subject to the same conditions.
KMNS - data clustering by the K-Means algorithm
by J.A. Hartigan & M.A. Wong

Code in this file based on Applied Statistics algorithms
(C) Royal Statistical Society 1979

According to Statlib (http://lib.stat.cmu.edu/apstat/):
"The Royal Statistical Society holds the copyright to these routines,
but has given its permission for their distribution provided that
no fee is charged."

added comments and declarations:  AJA, 2011-2017
----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
    ALGORITHM AS 136.1  APPL. STATIST. (1979) VOL.28, NO.1

    This is the optimal transfer stage.

    Each point is re-allocated, if necessary, to the cluster that
    will induce a maximum reduction in the within-cluster sum of
    squares.
----------------------------------------------------------------------*/
static void optra (float **a, int m, int n, float **c,
	int k, int *ic1, int *ic2, int *nc, float *an1, float *an2,
	int *ncp, float *d, int *itran, int *live, int *indx)
{
/* Initialized data */
/* Define BIG to be a very large positive number. */
	static float big = 1e30f;
	static float zero = 0.f;
	static float one = 1.f;

/* Local variables */
	int i, j, l, l1, l2;
	float r2, da, db, dc, dd, de, df;
	int ll;
	float rr, al1, al2, alt, alw;

/* Function Body */
/* If cluster L is updated in the last quick-transfer stage, it
   belongs to the live set throughout this stage.   Otherwise, at
   each step, it is not in the live set if it has not been updated
   in the last M optimal transfer steps. */
	for (l = 1; l <= k; ++l) {
		if (1 == itran[l]) live[l] = m + 1;
	}
	for (i = 1; i <= m; ++i) {
		++(*indx);
		l1 = ic1[i]; /* index of nearest cluster */
		l2 = ic2[i]; /* index of second-nearest cluster */
		ll = l2;
/* If point I is the only member of cluster L1, no transfer. */
		if (1 == nc[l1]) goto L90;
/* If L1 has not yet been updated in this stage, no need to re-compute D(I). */
		if (0 == ncp[l1]) goto L30;
		de = zero;
		for (j = 1; j <= n; ++j) {
			df = a[i-1][j-1] - c[l1-1][j-1];
			de += df * df;
		}
/* weighted sum of squares from point I to L1 */
		d[i] = de * an1[l1];
/* Find the cluster with minimum R2. */
L30:
		da = zero;
		for (j = 1; j <= n; ++j) {
			db = a[i-1][j-1] - c[l2-1][j-1];
			da += db * db;
		}
/* weighted sum of squares from point I to L2 */
		r2 = da * an2[l2];
		for (l = 1; l <= k; ++l) {
/* If I >= LIVE(L1), then L1 is not in the live set.   If this is
   true, we only need to consider clusters that are in the live set
   for possible transfer of point I.   Otherwise, we need to consider
   all possible clusters. */
			if ((i >= live[l1] && i >= live[l]) || l == l1 || l == ll) goto L60;
			rr = r2 / an2[l];
			dc = zero;
			for (j = 1; j <= n; ++j) {
				dd = a[i-1][j-1] - c[l-1][j-1];
				dc += dd * dd;
				if (dc >= rr) goto L60;
			}
			r2 = dc * an2[l]; /* a new second-nearest cluster */
			l2 = l;
L60:
			;
		}
		if (r2 < d[i]) goto L70;
/* If no transfer is necessary, L2 is the new IC2(I). */
		ic2[i] = l2;
		goto L90;
/* Update cluster centres, LIVE, NCP, AN1 & AN2 for clusters L1 and
    L2, and update IC1(I) & IC2(I). */
L70:
		*indx = 0;
		live[l1] = m + i;
		live[l2] = m + i;
		ncp[l1] = i;
		ncp[l2] = i;
		al1 = (float) nc[l1];
		alw = al1 - one;
		al2 = (float) nc[l2];
		alt = al2 + one;
		for (j = 1; j <= n; ++j) {
/* change center of L1 */
			c[l1-1][j-1] = (c[l1-1][j-1] * al1 - a[i-1][j-1]) / alw;
/* change center of L2 */
			c[l2-1][j-1] = (c[l2-1][j-1] * al2 + a[i-1][j-1]) / alt;
		}
/* change number of points in L1 and L2 */
		--nc[l1];
		++nc[l2];
/* change AN1 and AN2 */
		an2[l1] = alw / al1;
		an1[l1] = big;
		if (alw > one) an1[l1] = alw / (alw - one);
		an1[l2] = alt / al2;
		an2[l2] = alt / (alt + one);
/* update the nearest clusters to point I */
		ic1[i] = l2;
		ic2[i] = l1;
L90:
		if (*indx == m) return;
	}
	for (l = 1; l <= k; ++l) {
/* ITRAN(L) = 0 before entering QTRAN.   Also, LIVE(L) has to be
    decreased by M before re-entering OPTRA. */
		itran[l] = 0;
		live[l] -= m;
	}
	return;
} /* end of optra */

/*----------------------------------------------------------------------
    ALGORITHM AS 136.2  APPL. STATIST. (1979) VOL.28, NO.1

    This is the quick transfer stage.
    IC1(I) is the cluster which point I belongs to.
    IC2(I) is the cluster which point I is most likely to be
        transferred to.
    For each point I, IC1(I) & IC2(I) are switched, if necessary, to
    reduce within-cluster sum of squares.  The cluster centres are
    updated after each step.
----------------------------------------------------------------------*/
static void qtran (float **a, int m, int n, float **c,
	int k, int *ic1, int *ic2, int *nc, float *an1, float *an2,
	int *ncp, float *d, int *itran, int *indx)
{
/* Initialized data */
/* Define BIG to be a very large positive number */
	static float big = 1e30f;
	static float zero = 0.f;
	static float one = 1.f;

/* Local variables */
	int i, j, l1, l2;
	float r2, da, db, dd, de, al1, al2, alt, alw;
	int icoun, istep;

/* Function Body */
/* In the optimal transfer stage, NCP(L) indicates the step at which
   cluster L is last updated.   In the quick transfer stage, NCP(L)
   is equal to the step at which cluster L is last updated plus M. */
	icoun = 0;
	istep = 0;
L10:
	for (i = 1; i <= m; ++i) {
		++icoun;
		++istep;
		l1 = ic1[i];
		l2 = ic2[i];
/* If point I is the only member of cluster L1, no transfer. */
		if (1 == nc[l1]) goto L60;
/* If ISTEP > NCP(L1), no need to re-compute distance from point I to
   cluster L1.   Note that if cluster L1 is last updated exactly M
   steps ago, we still need to compute the distance from point I to
   cluster L1. */
		if (istep > ncp[l1]) goto L30;
		da = zero;
		for (j = 1; j <= n; ++j) {
			db = a[i-1][j-1] - c[l1-1][j-1];
			da += db * db;
		}
		d[i] = da * an1[l1];
/* If ISTEP >= both NCP(L1) & NCP(L2) there will be no transfer of
   point I at this step. */
L30:
		if (istep >= ncp[l1] && istep >= ncp[l2]) goto L60;
		r2 = d[i] / an2[l2];
		dd = zero;
		for (j = 1; j <= n; ++j) {
			de = a[i-1][j-1] - c[l2-1][j-1];
			dd += de * de;
			if (dd >= r2) goto L60;
		}
/* Update cluster centres, NCP, NC, ITRAN, AN1 & AN2 for clusters
   L1 & L2.   Also update IC1(I) & IC2(I).   Note that if any
   updating occurs in this stage, INDX is set back to 0. */
		icoun = 0;
		*indx = 0;
		itran[l1] = 1;
		itran[l2] = 1;
		ncp[l1] = istep + m;
		ncp[l2] = istep + m;
		al1 = (float) nc[l1];
		alw = al1 - one;
		al2 = (float) nc[l2];
		alt = al2 + one;
/* new centers */
		for (j = 1; j <= n; ++j) {
			c[l1-1][j-1] = (c[l1-1][j-1] * al1 - a[i-1][j-1]) / alw;
			c[l2-1][j-1] = (c[l2-1][j-1] * al2 + a[i-1][j-1]) / alt;
		}
		--nc[l1];
		++nc[l2];
/* update AN1 and AN2 */
		an2[l1] = alw / al1;
		an1[l1] = big;
		if (alw > one) an1[l1] = alw / (alw - one);
		an1[l2] = alt / al2;
		an2[l2] = alt / (alt + one);
		ic1[i] = l2;
		ic2[i] = l1;
/* If no re-allocation took place in the last M steps, return. */
L60:
		if (icoun == m) return;
	}
	goto L10;
} /* end of qtran */

/*----------------------------------------------------------------------

    ALGORITHM AS 136  APPL. STATIST. (1979) VOL.28, NO.1

    Divide M points in N-dimensional space into K clusters so that
    the within cluster sum of squares is minimized.

                          TRANSFER LIST
___Name_______Type______I/O_______Description___________________________
   a[m][n]    float**   In        the data matrix
   m          int       In        the number of points
   n          int       In        the number of dimensions
   c[k][n]    float**   input:    the matrix of initial cluster centers
                        output:   the matrix of final cluster centers
   k          int       In        the number of clusters
   ic1[m]     int*      Out       the cluster each point belongs to
   ic2[m]     int*      Neither   this array is used to remember the cluster
                                    which each point is most likely to be
                                    tranferred to at each step
   nc[k]      int*      Out       the number of points in each cluster
   an1[k]     float*    Neither   Ratio that corrects a square distance
                                    to a center if a point were removed
   an2[k]     float*    Neither   Ratio that corrects a square distance
                                    to a center if a point were added
   ncp[k]     int*      Neither   Step the cluster is no longer updated
   d[m]       float*    Neither   Distance from point to its nearest center
   itran[k]   int*      Neither   Relates to if there is convergence
   live[k]    int*      Neither   Relates to if there is convergence
   iter       int       In        the maximum number of iterations allowed
   wss[k]     float*    Out       the within-cluster sum of squares
                                    of each cluster

                        FAULT DIAGNOSTICS
   0   No fault
   1   At least one cluster is empty after the initial assignment.
         (A better set of initial cluster centres is called for)
   2   The allowed maximum number of iterations is exceeded
   3   K is less than or equal to 1 or greater than or equal to M

----------------------------------------------------------------------*/
int kmns (float **a, int m, int n, float **c, int k,
          int *ic1, int *ic2, int *nc, float *an1, float *an2,
	     int *ncp, float *d, int *itran, int *live, int iter, float *wss)
{
/* Define BIG to be a very large positive number */
	static float big = 1e30f;
	static float zero = 0.f;
	static float one = 1.f;

/* Local variables */
	int i, j, l;
	float aa, da, db, dc;
	int ii, ij, il;
	float dt[2];
	int indx;
	float temp;
	int ifault;

/* Parameter adjustments */
	--ic1;
	--ic2;
	--nc;
	--an1;
	--an2;
	--ncp;
	--d;
	--itran;
	--live;
	--wss;

/* Function Body */
/*----------------------------------------------------------------------
Begin.  Validate K, return error code 3 for K < 2 or K > M-1
----------------------------------------------------------------------*/
	if (k <= 1 || k >= m) return 3;
	ifault = 0;

/*----------------------------------------------------------------------
Initial values.  C must have initial values already assigned.
----------------------------------------------------------------------*/
/* For each point I, find its two closest centres, IC1(I) and
   IC2(I).  Assign it to IC1(I). */
	for (i = 1; i <= m; ++i) {
		ic1[i] = 1;
		ic2[i] = 2;
 /* square distance from centroids 1 and 2 to point */
		for (il = 1; il <= 2; ++il) {
			dt[il-1] = zero;
			for (j = 1; j <= n; ++j) {
				da = a[i-1][j-1] - c[il-1][j-1];
				dt[il-1] += da * da;
			}
		}
/* swap the square distances if appropriate */
		if (dt[0] > dt[1]) {
			ic1[i] = 2;
			ic2[i] = 1;
			temp = dt[0];
			dt[0] = dt[1];
			dt[1] = temp;
		}
/* compare the point to the centroids 3...K */
		for (l = 3; l <= k; ++l) {
			db = zero;
			for (j = 1; j <= n; ++j) {
				dc = a[i-1][j-1] - c[l-1][j-1];
				db += dc * dc;
				if (db >= dt[1]) goto L50;
			}
			if (db < dt[0]) goto L40;
/* better point than IC2 */
			dt[1] = db;
			ic2[i] = l;
			goto L50;
/* better point than IC1 */
L40:
			dt[1] = dt[0];
			ic2[i] = ic1[i];
			dt[0] = db;
			ic1[i] = l;
L50:
			;
		}
}
/* Update cluster centres to be the average of points contained
within them. */
	for (l = 1; l <= k; ++l) { /* set previous values to zero */
		nc[l] = 0;
		for (j = 1; j <= n; ++j)	c[l-1][j-1] = zero;
	}
	for (i = 1; i <= m; ++i) {
		l = ic1[i]; /* index of cluster nearest point I */
		++nc[l]; /* increment count of points in cluster */
		for (j = 1; j <= n; ++j) c[l-1][j-1] += a[i-1][j-1];
	}
/* Check to see if there is any empty cluster at this stage */
	for (l = 1; l <= k; ++l) {
		if (0 == nc[l]) return 1;
		aa = (float) nc[l];
		for (j = 1; j <= n; ++j) c[l-1][j-1] /= aa;
/* Initialize AN1, AN2, ITRAN & NCP
   AN1(L) = NC(L) / (NC(L) - 1)
   AN2(L) = NC(L) / (NC(L) + 1)
   ITRAN(L) = 1 if cluster L is updated in the quick-transfer stage,
            = 0 otherwise
   In the optimal-transfer stage, NCP(L) stores the step at which
   cluster L is last updated.
   In the quick-transfer stage, NCP(L) stores the step at which
   cluster L is last updated plus M. */
		an2[l] = aa / (aa + one);
		an1[l] = big;
		if (aa > one) an1[l] = aa / (aa - one);
		itran[l] = 1;
		ncp[l] = -1;
	}
	indx = 0;
/*----------------------------------------------------------------------
        Main loop
----------------------------------------------------------------------*/
	for (ij = 1; ij <= iter; ++ij) {
/*     In this stage, there is only one pass through the data.   Each
    point is re-allocated, if necessary, to the cluster that will
    induce the maximum reduction in within-cluster sum of squares. */
	optra (a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, itran, live, &indx);
/*     Stop if no transfer took place in the last M optimal transfer steps. */
		if (indx == m) goto L150;
/*     Each point is tested in turn to see if it should be re-allocated
    to the cluster to which it is most likely to be transferred,
    IC2(I), from its present cluster, IC1(I).   Loop through the
    data until no further change is to take place. */
		qtran (a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, itran, &indx);
/*    If there are only two clusters, there is no need to re-enter the
    optimal transfer stage. */
		if (2 == k) goto L150;
/*     NCP has to be set to 0 before entering OPTRA. */
		for (l = 1; l <= k; ++l) ncp[l] = 0;
	} /* end of main loop */
/*     Since the specified number of iterations has been exceeded, set
    IFAULT = 2.   This may indicate unforeseen looping. */
	ifault = 2;
/*----------------------------------------------------------------------
    Compute final centers
----------------------------------------------------------------------*/
L150:
	for (l = 1; l <= k; ++l) {
		wss[l] = zero;
		for (j = 1; j <= n; ++j) c[l-1][j-1] = zero;
	}
	for (i = 1; i <= m; ++i) {
		ii = ic1[i];
		for (j = 1; j <= n; ++j) c[ii-1][j-1] += a[i-1][j-1];
	}
	for (j = 1; j <= n; ++j) {
		for (l = 1; l <= k; ++l) c[l-1][j-1] /= (float) nc[l];
/* Compute within-cluster sum of squares for each cluster. */
		for (i = 1; i <= m; ++i) {
			ii = ic1[i];
			da = a[i-1][j-1] - c[ii-1][j-1];
			wss[ii] += da * da;
		}
	}
	return ifault;
} /* end of kmns */
