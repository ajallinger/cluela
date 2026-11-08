/*----------------------------------------------------------------------
ght.c - a generalized hyperplane tree for the k-means problem
by Andy Allinger, 2013-2020, released to the public domain
This program may be used by any person for any purpose.

------------------------------------------------------------------------
     The generalized hyperplane tree organizes a data set by partitioning,
like any binary tree.  Unlike the k-d tree, there is no requirement for
the data to be represented in R^N Euclidean space.  A data set which has
a dissimilarity between each datum may be split by choosing two points,
and assigning the remaining points to a partition by whether they are
nearer the first point or the second point.  If the data and dissimilarity
function satisfy the requirements of a metric space, then the triangle
inequality may be used to efficiently search the tree for nearest
neighbors of a point.

Early research was:
     I. Kalantari and G. McDonald
     "A data structure and an algorithm for the nearest point problem"
     IEEE Transactions on Software Engineering, v.9 p.5, 1983

cited in:
     Edgar Chavez, Gonzallo Navarro, and Jose Luis Marroquin
     "Searching in Metric Spaces"
     Association for Computing Machinery, 1999

It was discussed also in:
     Jeffrey K. Uhlmann
     "Satisfying General Proximity/Similarity Queries with Metric Trees"
     Information Processing Letters, v.40 pp.175-179, November 1991

------------------------------------------------------------------------
     Suppose the pivots are O and P, a query point is Q.  There may be
a hypothetical point S near the boundary surface.  The distance QS must
be bounded to compare it to r, the radius of the ball containing the known
neighbors of Q.  Construct triangle PQS, then by the triangle inequality:
     QS + PS .GE. PQ             =>      QS .GE. PQ - PS
Construct triangle OQS, then
     QS + OQ .GE. OS             =>      QS .GE. OS - OQ
Add
                     2 @ QS .GE. (PQ - OQ) + (OS - PS)
OS = PS by definition, so
                         QS .GE. (PQ - OQ) / 2
                          r < (PQ - OQ) / 2

------------------------------------------------------------------------
     Also keep the farthest distance from a pivot to any object beyond
it in the tree.  Call this distance f.  The ball around P of radius f
has no intersection with the ball around Q with radius r if
     r + f < PQ                =>      r < PQ - f
F may be arbitrarily close to zero or arbitrarily large.  In the condition
considered previously, OQ may be arbitrarily close to zero or arbitrarily
large.  It follows that neither condition implies the other, thus both
conditions may be used to prune the search.
----------------------------------------------------------------------*/
#include <stdbool.h>
#include <math.h>

/*----------------------------------------------------------------------
ghtdsm - a dissimilarity function

___Name_____Type_____In/Out___Description_______________________________
   a[p]     float*   In       Data vector
   b[p]     float*   In       Data vector
   p        int      In       Dimension of the data
   w[p]     float*   In       Importance weights
   robust   bool     In       Use robust statistics
   odwa     float    In       1 / w[j]
   ghtdsm   float    Out      Manhattan distance or squared Euclidean
----------------------------------------------------------------------*/
float ghtdsm (float *a, float *b, int p, float *w, bool robust, float odwa)
{
/* Local variables */
	int j;
	float r, retval;
/* Function Body */
	retval = 0.f;
	if (robust) { /* Manhattan */
		for (j = 0; j < p; ++j) retval += w[j] * fabsf(a[j] - b[j]);
	} else { /* squared Euclidean */
		for (j = 0; j < p; ++j) {
			r = a[j] - b[j];
			retval += w[j] * (r * r);
		}
	} /* robust ? */
	retval *= odwa;
	return retval;
} /* end of ghtdsm */

/*----------------------------------------------------------------------
ghtdbo - a dissimilarity function with bailout

___Name_____Type_____In/Out___Description_______________________________
   a[p]     float*   In       Data vector
   b[p]     float*   In       Data vector
   p        int      In       Dimension of the data
   w[p]     float*   In       Importance weights
   robust   bool     In       Use robust statistics
   odwa     float    In       1 / w[j]
   bail     float    In       Bail out if this value is exceeded
   ghtdbo   float    Out      Manhattan distance or squared Euclidean
----------------------------------------------------------------------*/
float ghtdbo (float *a, float *b, int p, float *w,
              bool robust, float odwa, float bail)
{
/* Local variables */
	int j;
	float r, retval;
/* Function Body */
	retval = 0.f;
	if (robust) { /* Manhattan */
		for (j = 0; j < p; ++j) {
			retval += w[j] * fabsf(a[j] - b[j]);
			if (retval > bail) goto L10;
		}
	} else { /* squared Euclidean */
		for (j = 0; j < p; ++j) {
			r = a[j] - b[j];
			retval += w[j] * (r * r);
			if (retval > bail) goto L10;
		}
	} /* robust ? */
L10:
	retval *= odwa;
	return retval;
} /* end of ghtdbo */

/*-----------------------------------------------------------------------
ghtini - initialize a generalized hyperplane tree

___Name____________Type_____I/O_____Description_________________________
   tree[2*k][12]   int**    Both    The generalized hyperplane tree
   last            int*     Out     Last filled index in tree
-----------------------------------------------------------------------*/
void ghtini (int **tree, int *last)
{
/* Constants */
	enum {STEM, LEAF, BUCK0, BUCK9 = 11};
/* Local variables */
	int j;
/* Function Body */
	*last = 0;
	tree[*last][STEM] = -1;  /* is root */
	tree[*last][LEAF] = -1;  /* is leaf */
	for (j = BUCK0; j <= BUCK9; ++j) tree[*last][j] = -1;  /* no data */
	return;
} /* end of ghtini */

/*-----------------------------------------------------------------------
ghtins - insert a cluster center into a generalized hyperplane tree

___Name____________Type______I/O____Description_________________________
   p               int       In     # variables each datum
   w[p]            float*    In     Importance weights
   c[k][p]         float**   In     Center data
   k               int       In     # of clusters
   robust          bool      In     Use Manhattan distance
   q               int       In     Index into c of a center to insert
   tree[2*k][12]   int**     Both   The generalized hyperplane tree
   crit[2*k][4]    float**   Both   Critical distances in tree
   last            int*      Both   Last filled index in tree
-----------------------------------------------------------------------*/
void ghtins (int p, float *w, float **c, int k, bool robust, int q,
             int **tree, float **crit, int *last)
{
/* Constants */
	enum {POINTS, POINTL, POINTR, PIVOTL, PIVOTR,
	      LEAF = 1, BUCK0 = 2, BUCK9 = 11};
	enum {FARL, FARR};
	static float BIG = 1.e30f;
/* Local variables */
	int g, h, i, j, l;
	int nl, pl, nr, pr;
	float d, dfar, dl, dr, fl, fr, odwa, wall;
/* Function Body */
	wall = 0.f;
	for (j = 0; j < p; ++j) wall += w[j];
	odwa = 1.f / wall;
	l = 0; /* the root node */
L10:
	if (-1 == tree[l][LEAF]) { /* leaf node? */
		for (h = BUCK0; h <= BUCK9; ++h) {
			if (-1 == tree[l][h]) { /* space remaining in node ? */
				tree[l][h] = q; /* insert the point */
				return; /* success */
			}
		}
		goto L50; /* node is full, so split */
	}
/* measure */
	pl = tree[l][PIVOTL];
	pr = tree[l][PIVOTR];
	dl = ghtdsm (c[pl], c[q], p, w, robust, odwa);
	dr = ghtdsm (c[pr], c[q], p, w, robust, odwa);
	if (dl < dr || ( (dl == dr) && (q % 2) )) {
		if (dl > crit[l][FARL]) crit[l][FARL] = dl;
		l = tree[l][POINTL]; /* climb left branch */
	} else {
		if (dr > crit[l][FARR]) crit[l][FARR] = dr;
		l = tree[l][POINTR]; /* climb right branch */
	}
	goto L10;
/* Split existing node; create two new leaf nodes. */
L50:
	dfar = -BIG; /* find longest distance between points in node */
	pl = tree[l][BUCK0]; /* avoid uninitialized warning */
	pr = tree[l][BUCK9];
	for (j = BUCK0; j <= BUCK9-1; ++j) {
		nr = tree[l][j];
		for (i = j+1; i <= BUCK9; ++i) {
			nl = tree[l][i];
			d = ghtdsm (c[nl], c[nr], p, w, robust, odwa);
			if (d > dfar) {
				dfar = d;
				pl = nl;
				pr = nr;
			}
		}
	}
/* make new nodes */
	++(*last);
	nl = *last;
	tree[nl][POINTS] = l;
	tree[nl][LEAF] = -1;
	for (h = BUCK0; h <= BUCK9; ++h) tree[nl][h] = -1;
	++(*last);
	nr = *last;
	tree[nr][POINTS] = l;
	tree[nr][LEAF] = -1;
	for (h = BUCK0; h <= BUCK9; ++h) tree[nr][h] = -1;
/* Sort */
	i = BUCK0;
	j = BUCK0;
	fl = 0.f;
	fr = 0.f;
	for (h = BUCK0; h <= BUCK9; ++h) {
		g = tree[l][h];
/* force pivot into its node */
		dl = (g == pl) ? -1.f : ghtdsm (c[g], c[pl], p, w, robust, odwa);
		dr = (g == pr) ? -1.f : ghtdsm (c[g], c[pr], p, w, robust, odwa);
		if (dl < dr || ( (dl == dr) && (g % 2) )) { /* copy datum to left branch */
			tree[nl][i] = g;
			if (dl > fl) fl = dl;
			++i;
		} else { /* copy datum to right branch */
			tree[nr][j] = g;
			if (dr > fr) fr = dr;
			++j;
		}
	}
/* current node becomes a stem node */
	tree[l][POINTL] = nl;
	tree[l][POINTR] = nr;
	tree[l][PIVOTL] = pl;
	tree[l][PIVOTR] = pr;
	crit[l][FARL] = fl;
	crit[l][FARR] = fr;
	goto L10;
} /* end of ghtins */

/*-----------------------------------------------------------------------
ghtsel - select nearest center from a generalized hyperplane tree

dist should be set to a very large value on input, unless a distance to
a center is already known.  Centers farther away will not be considered.

___Name____________Type______I/O____Description_________________________
   x[n][p]         float**   In     Data array
   n               int       In     # points in x
   p               int       In     # variables each datum
   w[p]            float*    In     Importance weights of variables
   c[k][p]         float**   In     Coordinates of centers
   k               int       In     Number of centers
   robust          bool      In     Use Manhattan distance
   self            bool      In     Don't permit to select NEAR = Q
   q               int       In     Index into x of a point to query
   livec[k]        int*      In     > 0 if each cluster is in the live set
   livex           int       In     > 0 if query point is in the live set
   near            int*      Out    Index into c of the nearest center
   dist            float*    Both   Distance to the nearest neighbor
   tree[2*k][12]   int**     Both   The generalized hyperplane tree
   crit[2*k][4]    float**   Both   Critical distances in tree
-----------------------------------------------------------------------*/
void ghtsel (float **x, int n, int p, float *w, float **c, int k,
             bool robust, bool self, int q, int *livec, int livex,
             int *near, float *dist, int **tree, float **crit)
{
/* Constants */
	enum {POINTS, POINTL, POINTR, PIVOTL, PIVOTR, MARK, DIR,
	      LEAF = 1, BUCK0 = 2, BUCK9 = 11};
	enum {FARL, FARR, KEEPY, KEEPZ};
	static float SLTH = 0.499975f, SLTO = 0.9999f;
/* Local variables */
	int i, j, l, nb, fb, pl, pr;
	float bail, d2, df, dl, dr, odwa, pq, r, wall, y, z;
/* Function Body */
	wall = 0.f;
	for (j = 0; j < p; ++j) wall += w[j];
	odwa = 1.f / wall;
	bail = *dist * wall;
	r = robust ? *dist : sqrtf(*dist);
	l = 0; /* begin at root */
L10: /* Come here when climbing. */
	if (-1 == tree[l][LEAF]) { /* inspect a leaf */
		for (j = BUCK0; j <= BUCK9; ++j) {
			i = tree[l][j];
			if (-1 == i) break; /* reached bottom of bucket */
			if (self && i == q) continue; /* don't select zero distance */
			if (livex <= 0 && livec[i] <= 0) continue; /* skip dead */
			d2 = ghtdbo (x[q], c[i], p, w, robust, odwa, bail);
			if (d2 < *dist) {
				*near = i;
				*dist = d2;
				bail = d2 * wall;
				r = robust ? d2 : sqrtf(d2);
			}
		}
		l = tree[l][POINTS];  /* back to stem node */
	} else {
		tree[l][MARK] = 0; /* clear mark in stem node */
		tree[l][DIR] = 0;
	}
L50: /* Come here when backing down. */
	if (-1 == l) return; /* is this the root node ? */
/* measure as required */
	if (-1 == tree[l][DIR]) { /* left branch nearer */
		nb = POINTL;
		fb = POINTR;
		y = crit[l][KEEPY];
		z = crit[l][KEEPZ];
	} else if (+1 == tree[l][DIR]) { /* right branch nearer */
		nb = POINTR;
		fb = POINTL;
		y = crit[l][KEEPY];
		z = crit[l][KEEPZ];
	} else { /* unknown */
		pl = tree[l][PIVOTL];
		pr = tree[l][PIVOTR];
		dl = ghtdsm (x[q], c[pl], p, w, robust, odwa);
		dr = ghtdsm (x[q], c[pr], p, w, robust, odwa);
		if (!robust) {
			dl = sqrtf(dl);
			dr = sqrtf(dr);
		}
		pq = fmaxf(dl, dr);
		z = fabsf(dl - dr) * SLTH;
		crit[l][KEEPZ] = z; /* remember */
		if (dl < dr || ( (dl == dr) && (q % 2) )) {
			nb = POINTL;
			fb = POINTR;
			tree[l][DIR] = -1; /* remember */
			df = crit[l][FARR];
		} else {
			nb = POINTR;
			fb = POINTL;
			tree[l][DIR] = +1;
			df = crit[l][FARL];
		}
		if (!robust) df = sqrtf(df);
		y = (pq - df) * SLTO;
		crit[l][KEEPY] = y;
	} /* need to measure ? */
/* Explore the near branch */
	if (0 == tree[l][MARK]) { /* near branch unmarked? */
		tree[l][MARK] = 1; /* mark */
		l = tree[l][nb];
		goto L10; /* and inspect the node */
	}
/* Consider far branch. */
	if (r > y && r > z && 1 == tree[l][MARK]) {
		tree[l][MARK] = 2; /* mark */
		l = tree[l][fb];
		goto L10; /* and inspect the node */
	}
/* backtrack */
	l = tree[l][POINTS];
	goto L50;
} /* end of ghtsel */
