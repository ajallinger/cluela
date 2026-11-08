/***********************************************************************
neugas.c - neural gas, a data clustering algorithm
by Andy Allinger, 2013-2020, released to the public domain
This program may be used by any person for any purpose.

Refer to:
     Thomas Martinetz and Klaus Schulten
     A "Neural-Gas" Network Learns Topologies
     in Artificial Neural Networks, Ed. T. Kohonen, K. Makisara,
     O. Simula, and J. Kangas.  Holland: Elsevier, 1991.

Also the description by Bernd Fritzke at:
     www.demogng.de/JavaPaper/node16.html

This routine follows Fritzke and omits the topology.

___Name______Type______In/Out____Description____________________________
   x[n][p]   float**   In        Data Matrix
   n         int       In        Number of objects
   p         int       In        Measurements per object
   c[k][p]   float**   Both      Cluster centers
   k         int       In        Number of clusters
   z[n]      int*      Out       The cluster each point belongs to
   ind[k]    int*      Neither   Index array
   d[k]      float*    Neither   Distances to each center
   u         float*    Out       Total squared distance of each point
                                    to its nearest center
***********************************************************************/
#include <stdlib.h>
#include <math.h>
void qsortr(float *x, int *ind, int n);

void neugas (float **x, int n, int p, float **c, int k, int *z,
             int *ind, float *d, float *u)
{
/* Revised parameters */
	const float BIG = 1.e30f,   /* arbitrary large number */
	            LAMI = 10.f,    /* range of adjustment, initial */
	            LAMF = 0.01f,   /* range of adjustment, final */
	            EPSI = 0.3f,    /* adjustment, initial */
	            EPSF = 0.005f;  /* adjustment, final */
	const int TAU = 40000;     /* iterations */
/* Local variables */
	int h, i, j, l;
	float best, eps, efact, lfact, lambda, r, tot, y;
/* Function Body */
	lambda = LAMI;
	eps = EPSI;
	lfact = powf( (LAMF / LAMI), (1.f / (float)(TAU - 1)) );
	efact = powf( (EPSF / EPSI), (1.f / (float)(TAU - 1)) );
/* main loop */
	for (h = 0; h < TAU; ++h) {
		r = (float)rand() / (float)RAND_MAX;  /* choose a datum at random */
		i = (int)(r * n);
		if (i > (n - 1)) i = n - 1;
/* find distance to each center */
		for (l = 0; l < k; ++l) {
			tot = 0.f;
			for (j = 0; j < p; ++j) {
				r = x[i][j] - c[l][j];
				tot += r * r;
			}
			d[l] = tot;
		}
		qsortr (d, ind, k);  /* sort distances */
		for (l = 0; l < k; ++l) d[ind[l]-1] = (float)l;  /* rank */
/* adapt the cluster centers */
		for (l = 0; l < k; ++l) {
			y = eps * expf(-d[l] / lambda);
			for (j = 0; j < p; ++j) c[l][j] += y * (x[i][j] - c[l][j]);
		}
/* adjust parameters */
		lambda *= lfact;
		eps *= efact;
	}
/* assignments and residual */
	*u = 0.f;
	for (i = 0; i < n; ++i) {
		best = BIG;
		for (l = 0; l < k; ++l) {
			tot = 0.f;
			for (j = 0; j < p; ++j) {
				r = x[i][j] - c[l][j];
				tot += r * r;
			}
			if (tot < best) {
				best = tot;
				z[i] = l + 1;
			}
		}
		*u += best;
	}
	return;
} /* end of neugas */
