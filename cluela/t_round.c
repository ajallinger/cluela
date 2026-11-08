/***********************************************************************
t_round.c - a test program for CLUELA
By Andy Allinger, 2017, released to the public domain
This program may be used by any person for any purpose.
***********************************************************************/
#include <stdbool.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int **new_int_mtx (int m, int n); /* defined in array.c */
float **new_flt_mtx (int m, int n);
int del_int_mtx (int **b, int m);
int del_flt_mtx (float **x, int m);
int cluela (float **x, int **nx, int n, int p, float *w, float **c, int **nc, /* defined in cluela.c */
      int kmin, int *k, int kmax, int *z, int *pop, float *f, bool robust,
      bool sosdo, int *iwork, float *rwork, int **tree, float **crit, float *u);
void dinit (float **x, int **nx, int n, int p, float *w, float **c, int **nc,
      int k, bool robust, float *shortx, float *cumd, float *thresh);
int cluste (float **x, int **nx, int n, int p, float *w, float **c, int **nc, /* defined in cluste.c */
      int kmin, int *k, int kmax, int *z, int *pop, bool robust,
      bool sosdo, int *iwork, float *rwork, int **tree, float **crit, float *u);
int kmns (float **a, int m, int n, float **c, int k, int *ic1, int *ic2,     /* defined in kmns.c */
      int *nc, float *an1, float *an2, int *ncp, float *d, int *itran,
      int *live, int iter, float *wss);

/* Main program */
int main (void)
{
/* Constants */
#define KMAX 10
#define KLIM 1000
#define N 30000
#define NS 12
#define P 2
	const int kmin = 1;
	static float w[P] = { 1.f,1.f };
	const float twopi = 6.2831853f;
 /* Local variables */
	int i, j, k, l, iter, z[N], iwork[4*N+KLIM], ifault, ferr, pop[KLIM];
	float r, u, f[KMAX], rnum, theta, rwork[2*N+KLIM], thresh;
	bool robust, sosdo;
	char ofile[32] = "round_output.csv";
	FILE *iou;
/* Workspace for kmns */
	int ic2[N], ncp[KLIM], itran[KLIM], live[KLIM];
	float an1[KLIM], an2[KLIM], d[N], wss[KLIM];
/* Timing variables */
	clock_t before, after;
	float result;
/* Allocate matrices */
	int **nc, **nx, **tree;
	float **c, **crit, **x;
	nc = new_int_mtx(KLIM, P);
	nx = new_int_mtx(N, P);
	tree = new_int_mtx(2*KLIM, NS);
	c = new_flt_mtx(KLIM, P);
	crit = new_flt_mtx(2*KLIM, 4);
	x = new_flt_mtx(N, P);
	if (!nc || !nx || !tree || !c || !crit || !x) goto L90;
/* initial values */
	for (i = 0; i < N; ++i) {
		for (j = 0; j < P; ++j) {
			nx[i][j] = 1;
		}
	}
	srand( (unsigned int)time(NULL));
	for (i = 0; i < N; ++i) { /* random points over unit circle */
		rnum = (float)rand() / (float)RAND_MAX;
		r = sqrtf(rnum); /* radius */
		rnum = (float)rand() / (float)RAND_MAX;
		theta = rnum * twopi; /* angle */
		x[i][0] = r * cosf(theta); /* Cartesian coordinates */
		x[i][1] =r * sinf(theta);
	}
	k = 6;
	robust = false;
	sosdo = ! robust;
/* Call subroutine */
	ifault = cluela (x, nx, N, P, w, c, nc, kmin, &k, KMAX, z, pop,
	  f, robust, sosdo, iwork, rwork, tree, crit, &u);
	printf("cluela returns #%d, made %d clusters with residual %f\n", ifault, k, u);
/* Write CSV file with Y-coordinates in different columns by cluster */
	iou = fopen(ofile, "w");
	if (NULL == iou) printf("t_round: trouble opening file %s\n", ofile);
	for (i = 0; i < N; ++i) {
		ferr = fprintf(iou, "%8.4f", x[i][0]);
		if (0 == ferr) printf("t_round: trouble writing X-coord\n");
/* one comma per cluster ID */
		for (l = 0; l < z[i]; ++l) {
			ferr = fprintf(iou, ",");
			if (0 == ferr) printf("t_round: trouble writing comma\n");
		}
		ferr = fprintf(iou, "%8.4f\n", x[i][1]);
		if (0 == ferr) printf("t_round: trouble writing Y-coord\n");
	}
/* Write the centroids in the far column */
	for (l = 0; l < k; ++l) {
		ferr = fprintf(iou, "%8.4f", c[l][0]);
		if (0 == ferr) printf("t_round: trouble writing X-coord\n");
		for (i = 0; i < k+1; ++i) {
			ferr = fprintf(iou, ",");
			if (0 == ferr) printf("t_round: trouble writing comma\n");
		}
		ferr = fprintf(iou, "%8.4f\n", c[l][1]);
		if (0 == ferr) printf("t_round: trouble writing Y-coord\n");
	}
/* Write the results for F */
	for (l = kmin; l <= KMAX; ++l) {
		ferr = fprintf(iou, "%d, %8.4f\n", l, f[l-kmin]);
		if (0 == ferr) printf("t_round: trouble writing F\n");
	}
	ferr = fclose(iou);
	if (0 != ferr) printf("t_round: trouble closing file.\n");
/* Robust mode is not comparable to KMNS! */
	if (robust) goto L90;
/* Try cluste */
	robust = false;
	sosdo = false;
	k = KLIM;
	before = clock();
	ifault = cluste (x, nx, N, P, w, c, nc, kmin, &k, KLIM, z, pop,
	  robust, sosdo, iwork, rwork, tree, crit, &u);
	after = clock();
	result = (float)(after - before) / (float)CLOCKS_PER_SEC;
	printf("cluste returns #%i, in %f seconds.  Made %i clusters with residual %f\n",
	  ifault, result, k, u);
/* Try kmns */
	k = KLIM;
	iter = 5000;
	dinit (x, nx, N, P, w, c, nc, k, robust, &rwork[0], &rwork[N], &thresh);
	before = clock();
	ifault = kmns (x, N, P, c, k, z, ic2, pop, an1, an2, ncp,
	  d, itran, live, iter, wss);
	u = 0;
	for (l = 0; l < k; ++l) u += wss[l];
	after = clock();
	result = (float)(after - before) / (float)CLOCKS_PER_SEC;
	printf("kmns returns #%i in %f seconds.  Made %i clusters with residual %f\n",
	  ifault, result, k, u);
L90: /* release memory */
	if (nc) del_int_mtx(nc, KLIM);
	if (nx) del_int_mtx(nx, N);
	if (tree) del_int_mtx(tree, 2*KLIM);
	if (c) del_flt_mtx(c, KLIM);
	if (crit) del_flt_mtx(crit, 2*KLIM);
	if (x) del_flt_mtx(x, N);
	return 0;
} /* end of main */

