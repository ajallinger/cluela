/*----------------------------------------------------------------------
cluela.c - Data clustering by the k-means method
By Andy Allinger, 2011-2017, released to the public domain.
This program may be used by any person for any purpose.
----------------------------------------------------------------------*/
#include <stdbool.h>
#include <float.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

void median (float *, int, float *); /* defined in median.c */
void fnv32 (void *, int, int32_t *); /* defined in fnv32.c */
bool safdiv (float, float); /* defined in safdiv.c */
void ghtini (int **, int *); /* defined in ght.c */
void ghtins (int, float *, float **, int, bool, int, int **, float **, int *);
void ghtsel (float **, int, int, float *, float **, int, bool, bool, int,
       int *, int, int *, float *, int **, float **);

/*----------------------------------------------------------------------
        almost equality test
----------------------------------------------------------------------*/
bool almeq (float r, float s)
{
/* Local variables */
	float ar, as;
	bool retval;
/* Begin. */
	ar = fabsf(r);
	as = fabsf(s);
/* all denormal numbers are equal */
	if (ar < FLT_MIN && as < FLT_MIN) {
		retval = true;
	} else {
		retval = fabsf(r - s) < 2.f * FLT_EPSILON * fmaxf(ar, as);
	}
	return retval;
} /* end of almeq */

/*----------------------------------------------------------------------
Pull a number out of a hat -- Uniform distribution 1 to n.
----------------------------------------------------------------------*/
int hat (int n)
{
	int j;
	float u;
	u = ((float) rand()) / ((float) RAND_MAX);
	j = ((int)(u * n)) + 1;
	if (j > n) j = n;
	return j;
} /* end of hat */

/*----------------------------------------------------------------------
Random permutation of 1...N
origin:  Durstenfeld, 1964
----------------------------------------------------------------------*/
void shuffl (int *i, int n)
{
/* Local variables */
	int j, k, l;
/* Parameter adjustments */
	--i;
/* Function Body */
	for (j = 1; j <= n; ++j) i[j] = j;
	for (j = n; j >= 2; --j) {
		k = hat(j);
		l = i[k];
		i[k] = i[j];
		i[j] = l;
	}
	return;
} /* end of shuffl */

/*----------------------------------------------------------------------
Bisection search for the first element greater than a value

___Name_______Type_____In/Out___Description_____________________________
   x[n]       float*   In       Array (in increasing order)
   n          int      In       Length of X
   v          float    In       The value sought
   i          int*     Out      Index where X > V found
----------------------------------------------------------------------*/
void bsgtr (float *x, int n, float v, int *i)
{
	int il, ir;
/* Function Body */
	if (n < 1) {
		*i = 0;
		return;
	}
/*V is greater than X(N) or less than X(1) */
	if (v < x[0]) {
		*i = 0;
		return;
	} else if (v > x[n-1]) {
		*i = n - 1;
		return;
	}
/*bisection search */
	il = 0;
	ir = n - 1;
L10:
	*i = (il + ir) / 2;
	if (*i == il) {
		if (x[*i] <= v) *i = ir;
		return;
	}
	if (x[*i] <= v) {
		il = *i;
	} else {
		ir = *i;
	}
	goto L10;
} /* end of bsgtr */

/*----------------------------------------------------------------------
Passive sort of unsigned integers by an unstable radix algorithm

___Name_____Type___In/Out___Description_________________________________
   ix[n]    int*   In       Input array
   ind[n]   int*   Out      Permutation vector
   n        int    In       Length of array
   m        int    In       Largest value in array
----------------------------------------------------------------------*/
void usortu (int *ix, int *ind, int n, int m)
{
/* Local variables */
	int b, i, j, k, l, ig[32], lg[32], il[32];
/* Parameter adjustments */
	--ix;
	--ind;
/* Initialize */
	for (i = 1; i <= n; ++i) ind[i] = i;
	if (n < 2) return;
	j = 8 * sizeof(int) - 1; /* greatest bit position */
	b = j; /* avoid compiler warning */
	for (i = j; i >= 0; --i) {
		b = i;
		if (m & (1 << i)) break;
	}
	l = b; /* most significant bit */
	il[l] = 1;
	ig[l] = n;
	lg[l] = -1;
L20: /* Sort */
	i = il[l];
	j = ig[l];
	if (i >= j) goto L50; /* if empty or single, back up */
L30: /* Find a set bit at the front of the partition */
	if (i <= j) {
		if (!(ix[ind[i]] & (1 << l))) {
			++i;
			goto L30;
		}
	}
L40: /* Find a clear bit at the rear of the partition */
	if (i <= j) {
		if (ix[ind[j]] & (1 << l)) {
			--j;
			goto L40;
		}
	}
	if (i < j) { /* Swap these elements */
		k = ind[i];
		ind[i] = ind[j];
		ind[j] = k;
		++i;
		--j;
		goto L30;
	}
/* Traverse tree */
	if (0 == l) goto L50; /* At the least significant bit, go back up */
	il[l-1] = il[l];
	ig[l-1] = j;
	lg[l] = -1;
	--l; /* Descend */
	goto L20;
L50: /* Back up */
	++l;
	if (l > b) return;
	if (lg[l] > 0) goto L50; /* Both subpartitions sorted; go back up */
/* Do the greater sub-region */
	il[l-1] = ig[l-1] + 1;
	ig[l-1] = ig[l];
	lg[l] = +1; /* Mark as greater */
	--l; /* Descend */
	goto L20;
} /* end of usortu */

/*----------------------------------------------------------------------
Dissimilarity function with missing data support and
an optional correction to minimize within-cluster square distance.

___Name_____Type_____In/Out___Description_______________________________
   x[p]     float*   In       Data vector of object
   nx[p]    int*     In       Data exists this observation
   c[p]     float*   In       Data vector of center
   nc[p]    int*     In       Count of data present this center
   p        int      In       Dimension of the data
   w[p]     float*   In       Importance weights
   robust   bool     In       Use robust statistics
   isin     int      In       -1 delete point, 0 neutral, +1 insert point
----------------------------------------------------------------------*/
float dsm (float *x, int *nx, float *c, int *nc, int p,
           float *w, bool robust, int isin)
{
/* Constants */
	static float BIG = 1e30f;
/* Local variables */
	int j, ncj;
	float dif, tot, fact, wsum, retval;
/* Function Body */
	tot = 0.f;
	wsum = 0.f;
	for (j = p-1; j >= 0; --j) {
		if (nx[j] > 0 && nc[j] > 0) {
			dif = x[j] - c[j];
			if (robust) { /* Manhattan distance */
				tot += w[j] * fabsf(dif);
			} else { /* Pythagorean distance */
				ncj = nc[j];
				if (0 == ncj + isin) {
					fact = 1.f;
				} else {
					fact = (float) ncj / (float) (ncj + isin);
				}
				tot += w[j] * fact * (dif * dif);
			}
			wsum += w[j];
		}
	}
/* adjust for missing values */
	if (safdiv (tot, wsum)) {
		retval = tot / wsum;
	} else {
		retval = BIG;
	}
	return retval;
} /* end of dsm */

/*----------------------------------------------------------------------
Count cluster populations

___Name_____Type___In/Out___Description_________________________________
   z[n]     int*   In       Membership of each point
   n        int    In       Number of points
   pop[k]   int*   Out      Populations each cluster
   k        int    In       Number of clusters
----------------------------------------------------------------------*/
void cntpop (int *z, int n, int *pop, int k)
{
/* Local variables */
	int i, l;
/* Parameter adjustments */
	--z;
	--pop;
/* Function Body */
	for (l = 1; l <= k; ++l) pop[l] = 0;
	for (i = 1; i <= n; ++i) {
		l = z[i];
		++pop[l];
	}
	return;
} /* end of cntpop */

/*----------------------------------------------------------------------
Let each center be the average of its points

___Name_______Type______In/Out____Description___________________________
   x[n][p]    float**   In        Data points
   nx[n][p]   int**     In        Data exists
   n          int       In        # of data points
   p          int       In        Variables each datum
   c[k][p]    float**   Out       Centers
   nc[k][p]   int**     Out       Center data available
   k          int       In        Number of clusters
   z[n]       int*      In        Cluster membership, each datum
   pop[k]     int*      In        Population, each cluster
   robust     bool      In        Median instead of mean
   livec[k]   int*      In        .GE. 1 => Need to recompute this center
   iperm[n]   int*      Neither   Workspace
   work[n]    float*    Neither   Workspace
----------------------------------------------------------------------*/
void cent (float **x, int **nx, int n, int p, float **c, int **nc, int k,
           int *z, int *pop, bool robust, int *livec, int *iperm, float *work)
{
/* Local variables */
	int h, i, j, l, i1, i5, i9;
/* Function Body */
	if (robust) {
		usortu (z, iperm, n, k); /* Sort implicitly by cluster ID */
		for (j = 0; j < p; ++j) { /* averages */
			i9 = -1; /* index of last object in cluster */
			for (l = 0; l < k; ++l) {
				i1 = i9 + 1; /* index of first object in cluster */
				i9 += pop[l];
				if (livec[l] <= 0) continue; /* skip dead */
				h = 0;
				for (i5 = i1; i5 <= i9; ++i5) {
					i = iperm[i5] - 1;
					if (nx[i][j] > 0) {
						work[h] = x[i][j];
						++h;
					}
				}
				nc[l][j] = h;
				if (h > 0) median(work, h, &c[l][j]);
			}
		}
	}
	else {
		for (l = 0; l < k; ++l) {
			for (j = 0; j < p; ++j) {
				c[l][j] = 0.f;
				nc[l][j] = 0;
			}
		}
		for (i = 0; i < n; ++i) {
			l = z[i] - 1;
			for (j = 0; j < p; ++j) {
				if (nx[i][j] > 0) {
					c[l][j] += x[i][j];
					++nc[l][j];
				}
			}
		}
		for (l = 0; l < k; ++l) {
			for (j = 0; j < p; ++j) {
				if (nc[l][j] > 0) c[l][j] /= (float)nc[l][j];
			}
		}
	}
	return;
} /* end of cent */

/*----------------------------------------------------------------------
Initialize cluster centers based on distance

Refer to:
   "kmeans++: the advantages of careful seeding"
   David Arthur and Sergei Vassilvitskii
   Proceedings of the eighteenth annual ACM-SIAM Symposium
     on Discrete Algorithms, 2007

___Name________Type______I/O_______Description__________________________
   x[n][p]     float**   In        Data points
   nx[n][p]    int**     In        Data exists this observation
   n           int       In        Number of points
   p           int       In        Dimension of the data
   w[p]        float*    In        Importance weights
   c[k][p]     float**   Out       Center points of clusters
   nc[k][p]    int**     Out       Count of center data available
   k           int       In        Number of clusters to make
   robust      bool      In        Use medians instead of means
   shortx[n]   float*    Neither   Squared distance to each point
   cumd[n]     float*    Neither   Cumulative squared distance
   thresh      float*    Out       Threshhold distance for dp-means
----------------------------------------------------------------------*/
void dinit (float **x, int **nx, int n, int p, float *w, float **c, int **nc,
           int k, bool robust, float *shortx, float *cumd, float *thresh)
{
/* Initialized data */
	static float BIG = 1e30f;
/* Local variables */
	int i, j, l;
	float d, pad, tot, u;
/* Begin.  Choose each center with probability proportional to its
   least squared distance from any existing center. */
	for (i = 0; i < n; ++i) shortx[i] = BIG;
/* choose first center at random */
	i = hat(n) - 1;
	for (j = 0; j < p; ++j) { /* copy */
		c[0][j] = x[i][j];
		nc[0][j] = nx[i][j];
	}
	for (l = 1; l < k; ++l) {
		pad = 0.f;
		for (i = 0; i < n; ++i) { /* compare each point to previous center */
			d = dsm (x[i], nx[i], c[l-1], nc[l-1], p, w, robust, 0);
			if (d < shortx[i]) shortx[i] = d;
			if (shortx[i] < pad) pad = shortx[i];
		}
		tot = 0.f;
		for (i = 0; i < n; ++i) {
			tot += shortx[i] - pad;
			cumd[i] = tot;
		}
/* uniform at random over cumulative distance */
		u = (float) rand() / (float) RAND_MAX * tot;
		bsgtr (cumd, n, u, &i);
/* assign center */
		for (j = 0; j < p; ++j) {
			c[l][j] = x[i][j];
			nc[l][j] = nx[i][j];
		}
	}
/* Threshhold distance is longest distance of any point to its center */
	*thresh = 0.f;
	for (i = 0; i < n; ++i) if (shortx[i] > *thresh) *thresh = shortx[i];
	return;
} /* end of dinit */

/*----------------------------------------------------------------------
Partition clustering for mixed and missing data
origin:  Hugo Steinhaus, 1956

Refer to:

   "Revisiting K-means:  New Algorithms via Bayesian Nonparametrics"
   Brian Kulis and Michael I. Jordan
   Proceedings of the 29th International Conference on Machine Learning, 2012

___Name____________Type______I/O_______Description______________________
   x[n][p]         float**   In        Data points
   nx[n][p]        int**     In        Data exists this observation
   n               int       In        Number of points
   p               int       In        Dimension of the data
   w[p]            float*    In        Importance weights
   c[kmax][p]      float**   Both      Center points of clusters
   nc[kmax][p]     int**     Both      Count of center data available
   kmin            int       In        Fewest clusters to make
   k               int*      Both      Default number of clusters to make
   kmax            int       In        Most clusters to make
   z[n]            int*      Out       What cluster a point is in
   pop[kmax]       int*      Out       # of points in each cluster
   miss            bool      In        Some data is missing
   robust          bool      In        Use medians instead of means
   sosdo           bool      In        Sum of squared distance objective
   thresh          float     In        Threshhold for forming new cluster
   livec[kmax]     int*      Neither   .GE. 1 => cluster needs updating
   livex[n]        int*      Neither   If each point is in the live set
   iperm[n]        int*      Neither   Permutation vector
   tree[2*k][12]   int**     Neither   generalized hyperplane tree
   crit[2*k][4]    float**   Neither   critical distances in tree
   shortc[kmax]    float*    Neither   Shortest distance to another center
   shortx[n]       float*    Neither   Shortest distance to current center
   work[n]         float*    Neither   For averages
   u               float*    Out       Residual
----------------------------------------------------------------------*/
int kclust (float **x, int **nx, int n, int p, float *w, float **c, int **nc,
  int kmin, int *k, int kmax, int *z, int *pop, bool miss, bool robust, bool sosdo,
  float thresh, int *livec, int *livex, int *iperm, int **tree, float **crit,
  float *shortc, float *shortx, float *work, float *u)
{
/* Constants */
	static int ZERO = 0;
	static int ONE = 1;
	static int HS = 13;
	static int ITER = 1000;
	static float BIG = 1e30f;
/* Local variables */
	int h, i, j, l, m, i1, k0, hl;
	int hash, ifar, lold, lnew, isin, hist[HS];
	float d, dfar, best, wall;
	bool change, metric;
/* Function Body */
	if (*k < 1 || *k > n) return 3;
	hl = n * sizeof(int);
	metric = ! (miss || sosdo);
	for (j = 0; j < p; ++j) {
		if (w[j] < 0.f)  metric = false;
	}
/* initialize */
	k0 = *k;
	for (l = 0; l < k0; ++l) {
		livec[l] = 2;
		shortc[l] = BIG;
	}
	for (i = 0; i < n; ++i) {
		livex[i] = 1;
		z[i] = 1; /* arbitrary valid ID */
		shortx[i] = BIG;
	}

/*----------------------------------------------------------------------
                 main loop
----------------------------------------------------------------------*/
	for (h = 0; h < ITER; ++h) {
		change = false;
		shuffl (iperm, n);
		k0 = *k;

/*----------------------------------------------------------------------
     Put centers into generalized hyperplane tree
----------------------------------------------------------------------*/
		if (metric) {
			ghtini (tree, &i);
			for (l = 0; l < k0; ++l) ghtins (p, w, c, k0, robust, l, tree, crit, &i);

/*----------------------------------------------------------------------
If the distance from a point to its center is less than or equal to half
the distance of the center to its nearest center, no transfer can occur.
----------------------------------------------------------------------*/
			for (l = 0; l < k0; ++l) { /* nearest center each center */
				best = BIG;
				ghtsel (c, k0, p, w, c, k0, robust, true, l, livec, ONE, &lnew, &best, tree, crit);
				shortc[l] = robust ? best * .50f : best * .25f;
			}
		}

/*----------------------------------------------------------------------
Find nearest center for each point.  Points are in the live set if the
distance to the nearest center has increased.
Clusters updated on this step are set LIVEC = 2 and will be recomputed.
Clusters updated on the prior step have LIVEC = 1 and are still live.
If neither the point nor the cluster are in their respective live sets,
skip the comparison.
----------------------------------------------------------------------*/
		dfar = -BIG;
		ifar = 0;
		for (i1 = 0; i1 < n; ++i1) {
			i = iperm[i1] - 1;
			lold = z[i] - 1; /* old assignment */
			lnew = lold;
			best = shortx[i];
			if (metric) {
				if (best <= shortc[lold]) continue; /* don't change point */
				ghtsel (x, n, p, w, c, k0, robust, false, i, livec, livex[i], &lnew, &best, tree, crit);
			} else {
				for (l = 0; l < k0; ++l) { /* exhaustive search */
					if (livex[i] <= 0 && livec[l] <= 0) continue;
					if ((!sosdo) || h == 1) {
						isin = 0;
					} else if (l == lold) {
						isin = -1;
					} else {
						isin = +1;
					}
					d = dsm (x[i], nx[i], c[l], nc[l], p, w, robust, isin);
					if (d < best) {
						best = d;
						lnew = l;
					}
				}
			}
			shortx[i] = best;
			if (best > dfar) {
				dfar = best;
				ifar = i;
			}
/* check for stray points */
			if (almeq (best, BIG)) return 4;

/*----------------------------------------------------------------------
             reassign point
----------------------------------------------------------------------*/
			if (lold != lnew) {
				livec[lold] = 2; /* recenter old cluster */
				livec[lnew] = 2; /* recenter new cluster */
				z[i] = lnew + 1;
				change = true;
				if (sosdo && h > 0) { /* update centers */
					for (j = 0; j < p; ++j) {
						if (nx[i][j] > 0) {
							m = nc[lold][j] - 1;
							if (m < 1) m = 1; /* don't divide by zero */
							c[lold][j] -= (x[i][j] - c[lold][j]) / (float) m;
							--nc[lold][j];
							c[lnew][j] += (x[i][j] - c[lnew][j]) / (float) (nc[lnew][j] + 1);
							++nc[lnew][j];
						}
					}
				}
			}
		}

/*----------------------------------------------------------------------
       if shortest distance exceeds a threshhold, make new cluster
----------------------------------------------------------------------*/
		if ( (!change) && (dfar > thresh) ) {
			if (k0 >= kmax) return 5;
			++(*k);
			lnew = *k - 1;
			lold = z[ifar] - 1;
			livec[lold] = 2;
			livec[lnew] = 2;
			z[ifar] = lnew + 1;
			change = true;
		}

/* algorithm halts */
		if ( (!change) && (h > 0) ) goto L80;

/*----------------------------------------------------------------------
              delete empty clusters
----------------------------------------------------------------------*/
		k0 = *k; /* count populations */
		cntpop (z, n, pop, k0);
		for (l = k0-1; l >= 0; --l) {
			if (0 == pop[l]) {
				if (*k <= kmin) return 1;
				for (j = 0; j < p; ++j) {
					c[l][j] = c[*k-1][j];
					nc[l][j] = nc[*k-1][j];
				}
				for (i = 0; i < n; ++i) if (z[i] == *k) z[i] = l + 1;
				pop[l] = pop[*k-1];
				livec[l] = livec[*k-1];
				shortc[l] = shortc[*k-1];
				pop[*k-1] = 0;
				--(*k);
			}
		}

/*----------------------------------------------------------------------
          find cluster centers
----------------------------------------------------------------------*/
		k0 = *k; /* decrement LIVEC */
		for (l = 0; l < k0; ++l) --livec[l];
		cent (x, nx, n, p, c, nc, k0, z, pop, robust, livec, iperm, work);

/*----------------------------------------------------------------------
           check if anything is changing
----------------------------------------------------------------------*/
		hash = -2128831035; /* canonical arbitrary initialization */
		fnv32 (z, hl, &hash); /* hash of membership vector */
/* exit on loop criterion */
		i = (h < HS) ? h : HS;
		for (j = 0; j < i; ++j) if (hash == hist[j]) goto L80;
		j = h % HS;
		hist[j] = hash; /* remember */
/* Determine live points */
		for (i = 0; i < n; ++i) livex[i] = 0;
		isin = sosdo ? -1 : 0;
		for (i = 0; i < n; ++i) {
			l = z[i] - 1;
			if (livec[l] > 0) {
				d = dsm (x[i], nx[i], c[l], nc[l], p, w, robust, isin);
/* points now farther are live */
				if (d > shortx[i]) livex[i] = 1;
				shortx[i] = d;
			}
		}
	} /* next iteration */
	return 2;

/*----------------------------------------------------------------------
             compute residual
----------------------------------------------------------------------*/
L80:
	*u = 0.f;
	for (i = 0; i < n; ++i) {
		l = z[i] - 1;
		d = dsm (x[i], nx[i], c[l], nc[l], p, w, robust, ZERO);
		*u += d;
	}
	wall = 0.f;
	for (j = 0; j < p; ++j) wall += w[j];
	*u *= wall;
	return 0;
} /* end of kclust */

/*----------------------------------------------------------------------
Validate input, cluster several tries, and keep the best result.
Choose the number of clusters.

___Variable___________Type______I/O_______Description___________________
   x[n][p]            float**   In        Data points in R^P
   nx[n][p]           int**     In        1 indicates data is present for
                                           the corresponding X value.
   n                  int       In        Number of points
   p                  int       In        Measurements per point
   w[p]               float*    In        Importance weights of variables
   c[kmax][p]         float**   Out       Cluster centroids
   nc[kmax][p]        int**     Out       Count of centroid data available
   kmin               int       In        Fewest clusters to make
   k                  int*      Both      Default/actual clusters to make
   kmax               int       In        Most clusters to make
   z[n]               int*      Out       What cluster each point is in
   pop[kmax]          int*      Out       Number of points in each cluster
   f[kmax-kmin+1]     float*    Out       Low value indicates good # clusters
   robust             bool      In        Use L1 statistics
   sosdo              bool      In        Sum of square distance objective
   iperm[n]           int*      Neither   Orders points by cluster
   zlof[n]            int*      Neither   Saved memberships with low F value
   ztry[n]            int*      Neither   Membership array passed to KCLUST
   livex[n]           int*      Neither   Whether each point is live
   livec[kmax]        int*      Neither   Whether each cluster is live
   tree[2*kmax][12]   int**     Neither   generalized hyperplane tree
   crit[2*kmax][4]    float**   Neither   critical distances in tree
   shortx[n]          float*    Neither   Shortest distance each point
   cumd[n]            float*    Neither   Cumulative distribution of distance
   shortc[kmax]       float*    Neither   Shortest distance, each center
   u                  float*    Out       Residual error
----------------------------------------------------------------------*/
int clue99 (float **x, int **nx, int n, int p, float *w, float **c, int **nc,
      int kmin, int *k, int kmax, int *z, int *pop, float *f, bool robust,
      bool sosdo, int *iperm, int *zlof, int *ztry, int *livex, int *livec,
      int **tree, float **crit, float *shortx, float *cumd, float *shortc, float *u)
{
/* Constants */
	static int C_TRIES = 10;
	static float FZERO = 0.f;
	static float BIG = 1e30f;
	static float CUT = 0.85f;
/* Local variables */
	int g, h, i, j, l, k0, ifault, ktry, kbest;
	float a, fbest, ubest, uret, uprev, thresh, wall;
	bool ok, miss, known;
/* Function Body */
	kbest = *k;
	ok = true;
	if (kmin < 1 || *k < kmin || *k > kmax || kmax > n) return 6;
	if (robust && sosdo) return 7;
	if (p < 1 || n < 1) return 8;
	for (j = 0; j < p; ++j) if (! almeq(w[j], FZERO)) goto L10;
	return 9;
L10: /* check for data with no values present */
	for (i = 0; i < n; ++i) {
		for (j = 0; j < p; ++j) if (nx[i][j] > 0) goto L20;
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

/*----------------------------------------------------------------------
        case K = N
----------------------------------------------------------------------*/
	if (kmin == n) {
		for (i = 0; i < n; ++i) {
			for (j = 0; j < p; ++j) {
				c[i][j] = x[i][j];
				nc[i][j] = nx[i][j];
			}
		}
		for (i = 0; i < n; ++i) z[i] = i + 1;
		k0 = *k;
		for (l = 0; l < k0; ++l) pop[l] = 1;
		*u = 0.f;
		return 0;
	}

/*----------------------------------------------------------------------
        case K = 1
----------------------------------------------------------------------*/
	if (1 == kmax) {
		for (i = 0; i < n; ++i) z[i] = 1;
		goto L90; /* and compute center */
	}

/*----------------------------------------------------------------------
               initialize alpha and fk (see below)
----------------------------------------------------------------------*/
	known = (kmin == kmax);
	if (known) {
		k0 = *k;
		f[*k-kmin] = 1.f;
	} else {
		k0 = kmin - 1; /* bound of loop */
		if (k0 < 1) k0 = 1;
	}
	fbest = BIG;
	ubest = 0.f;
	uret = 0.f;
	a = 1.f - 3.f / (p * 4.f);
	for (g = 3; g < k0; ++g) a += (1.f - a) / 6.f;

/*----------------------------------------------------------------------
              Loop over possible values of k
----------------------------------------------------------------------*/
	for (g = k0; g <= kmax; ++g) {
		uprev = ubest;
		ubest = BIG;
		ok = false;
/* seek least residual */
		for (h = 0; h < C_TRIES; ++h) {
			ktry = g;
			dinit (x, nx, n, p, w, c, nc, ktry, robust, shortx, cumd,
			  &thresh);
			thresh = BIG;
			ifault = kclust (x, nx, n, p, w, c, nc, ktry, &ktry, ktry, ztry,
			  pop, miss, robust, sosdo, thresh, livec, livex, iperm, tree, crit,
			  shortc, shortx, cumd, u);
			if (0 == ifault && *u < ubest) {
				ubest = *u;
				ok = true;
				if (known) kbest = ktry;
				for (i = 0; i < n; ++i) z[i] = ztry[i];
				if (1 == ktry) break; /* save computation for 1 cluster */
			}
		}
/* failure */
		if (! ok) {
			*k = g;
			return ifault;
		}
		if (known) { /* ignore f criterion */
			*k = kbest;
			*u = ubest;
			goto L90; /* bail out of G loop */
		}

/*----------------------------------------------------------------------
Choose the number of clusters.  refer to:
    "Selection of K in K-means clustering"
     D T Pham, S S Dimov, and C D Nguyen
     Journal of Mechanical Engineering Science, v.219, 2005 [2004]

They propose the formula:
             1,                            k = 1
     f(k) =  S[k] / (a[k] @ S[k-1]),       S[k-1] .NE. 0, k > 1
             1,                            S[k-1] = 0

where
             1 - 3 / (4 @ n),              k = 2, n > 1
     a[k] =
             a[k-1] + (1 - a[k]) / 6       k > 2, n > 1

and k is the number of clusters, n is the number of dimensions, and
S[k] is the residual from clustering with k clusters.
----------------------------------------------------------------------*/
		if (1 == g) {
			a = 0.f;
		} else if (2 == g) {
			a = 1.f - 3.f / (p * 4.f);
		} else {
			a += (1.f - a) / 6.f;
		}
/* don't select too few */
		if (g < kmin) continue;
		if (1 == g || almeq (uprev, FZERO)) {
			f[g-kmin] = 1.f;
		} else {
			f[g-kmin] = ubest / (a * uprev);
		}

/*----------------------------------------------------------------------
Low F indicates distinct clusters.  Accept the default number of clusters
if the lowest value of F is not less than the CUT parameter.
----------------------------------------------------------------------*/
		if ( (f[g-kmin] < CUT && f[g-kmin] < fbest) || (fbest >= CUT && g == *k) ) {
			fbest = fminf(fbest, f[g-kmin]);
			kbest = g;
			uret = ubest;
/* stash best Z for this K */
			for (i = 0; i < n; ++i) zlof[i] = z[i];
		}
	}
/* set return values */
	*k = kbest;
	*u = uret;
	for (i = 0; i < n; ++i) z[i] = zlof[i];

/*----------------------------------------------------------------------
  Count cluster populations, return centers
----------------------------------------------------------------------*/
L90:
	if (ok) ifault = 0; /* clear error codes */
	k0 = *k; /* count populations */
	cntpop (z, n, pop, k0);
	for (l = 0; l < k0; ++l) livec[l] = 1; /* all clusters live */
	cent (x, nx, n, p, c, nc, k0, z, pop, robust, livec, iperm, cumd);
	return ifault;
} /* end of clue99 */

/*----------------------------------------------------------------------
CLUELA - CLUstering with ELAborations

___Variable______Type______I/O_______Description_____________________
   x[n][p]       float**   In        Data points in R^P
   nx[n][p]      int**     In        1 indicates data is present for
                                           the corresponding X value.
                                      0 indicates data is missing for
                                           the corresponding X value.
   n             int       In        Number of points
   p             int       In        Measurements per point
   w[p]          float*    In        Importance weights of variables
   c[kmax][p]    float*    Out       Cluster centroids
   nc[kmax][p]   int*      Out       Count of data available each center
   kmin          int       In        Fewest clusters to make
   k             int*      Both      Default/actual clusters to make
   kmax          int       In        Most clusters to make
   z[n]          int*      Out       What cluster each point is in
   pop[kmax]     int*      Out       Number of points in each cluster
   f[kmax-kmin+1]
                 float*    Out       Low value indicates good # clusters
   robust        bool      In        Use L1 statistics
   sosdo         bool      In        Sum of squared distance objective
   iwork[4*n+kmax]
                 int*      Neither   Workspace
   rwork[2*n+kmax]
                 float*    Neither   Workspace
   tree[2*kmax][12]
                 int**     Neither   Workspace
   crit[2*kmax][4]
                 float**   Neither   Workspace
   u             float*    Out       Residual error

RETURN VALUE:
                                         0  Normal return
                                         1  Empty cluster
                                         2  Max iterations exceeded
                                         3  K out of bounds (in kclust)
                                         4  Point cannot be compared
                                         5  Formed too many clusters
                                         6  K out of bounds (in clue99)
                                         7  Conflicting options
                                         8  Array limits out of bounds
                                         9  All weights are zero
                                        10  Datum with no values
                                        11  Weights sum to zero
----------------------------------------------------------------------*/
int cluela (float **x, int **nx, int n, int p, float *w, float **c, int **nc,
      int kmin, int *k, int kmax, int *z, int *pop, float *f, bool robust,
      bool sosdo, int *iwork, float *rwork, int **tree, float **crit, float *u)
{
/* Function Body */
	int ifault;
	ifault = clue99 (x, nx, n, p, w, c, nc, kmin, k, kmax, z, pop, f,
	  robust, sosdo, iwork, &iwork[n], &iwork[2*n], &iwork[3*n], &iwork[4*n],
	  tree, crit, rwork, &rwork[n], &rwork[2*n], u);
	return ifault;
} /* end of cluela */
