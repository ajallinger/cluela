/*======================================================================
t_wine.c - Cluster the Wine data
by Andy Allinger, 2011-2017, public domain
This program may be used by any person for any purpose.
======================================================================*/
#include <stdbool.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int **new_int_mtx (int m, int n); /* defined in array.c */
float **new_flt_mtx (int m, int n);
char **new_ch_mtx (int m, int n);
int del_int_mtx (int **b, int m);
int del_flt_mtx (float **x, int m);
int del_ch_mtx (char **h, int m);
void avedev (float **x, int **nx, int n, int p, /* defined in stats.c */
             int *na, float *ave, float *dev);
bool safdiv (float numer, float denom); /* defined in cluela.c */
int cluela (float **x, int **nx, int n, int p, float *w, float **c, int **nc,
      int kmin, int *k, int kmax, int *z, int *pop, float *f, bool robust,
      bool sosdo, int *iwork, float *rwork, int **tree, float **crit, float *u);
int pca (float **x, int n, int p, bool robust, float **r, int m, int *q, /* defined in mds.c */
         float *a, float *e, float **v, float *t, float *u);
int wcsvsp (float **r, int n, int *z, FILE *iou); /* defined in dump.c */
float arind (int *y, int *z, int n, int g, int k, /* defined in valid.c */
            float **m, float *a, float *b);
float vinfo (int *y, int *z, int n, int g, int k,
            float **m, float *a, float *b);

/* Main program */
int main (void)
{
/* Constants */
#define FNLEN 32   /* length of file name */
#define G 3        /* # classes */
#define M 2        /* projection dimension for graphs */
#define N 178      /* # of objects */
#define NS 12      /* node size of binary tree */
#define P 13       /* # of attributes */
#define KMAX 20    /* most clusters to make */
#define TXTLEN 256 /* length of input buffer */
	const int kmin = 1;
	bool sosdo = false;
/* Local variables */
	int i, j, k, l, q, y[N], z[N], na[P], pop[KMAX], iwork[4*N+KMAX], ferr, ifault;
	float u, e[P], f[KMAX], ave[P], dev[P], rwork[2*N+KMAX], aval[G], bval[KMAX];
	float accari = -1.f,
	      accvi = -1.f;
	bool robust = false;
	char text[TXTLEN];
	char *sptr;
	FILE *iou;
/* Allocate matrices */
	int **nc, **nx, **tree;
	float **c, **crit, **mval, **r, **v, **x;
	nc = new_int_mtx(N, KMAX);
	nx = new_int_mtx(N, P);
	tree = new_int_mtx(2*KMAX, NS);
	c = new_flt_mtx(N, KMAX);
	crit = new_flt_mtx(2*KMAX, 4);
	mval = new_flt_mtx(G, KMAX);
	r = new_flt_mtx(N, M);
	v = new_flt_mtx(P, P);
	x = new_flt_mtx(N, P);
	if (!nc || !nx || !tree || !c || !crit || !mval || !r || !x) goto L90;

	static char fname[FNLEN] = "wine.data";
/*----------------------------------------------------------------------
Wine recognition data
From the University of California at Irvine Machine Learning repository

       Forina, M. et al, PARVUS - An Extendible Package for Data
       Exploration, Classification and Correlation. Institute of Pharmaceutical
       and Food Analysis and Technologies, Via Brigata Salerno,
       16147 Genoa, Italy.

Attributes:

 	1) Alcohol
 	2) Malic acid
 	3) Ash
	4) Alcalinity of ash
 	5) Magnesium
	6) Total phenols
 	7) Flavanoids
 	8) Nonflavanoid phenols
 	9) Proanthocyanins
	10)Color intensity
 	11)Hue
 	12)OD280/OD315 of diluted wines
 	13)Proline
----------------------------------------------------------------------*/

	static float w[P] = { 1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f };
	static float wspec[P] = { 2.76f,.648f,1.42f,.581f,.695f,.0516f,1.92f,.05f,
	  .05f,1.66f,1.9f,3.11f,4.32f };
	for (i = 0; i < N; ++i) {
		for (j = 0; j < P; ++j) {
			nx[i][j] = 1;
		}
	}
	srand( (unsigned int)time(NULL));
	ifault = 0;
/* open input file */
	iou = fopen(fname, "r");
	if (NULL == iou) {
		printf("trouble opening input file\n");
		return EXIT_FAILURE;
	}
/* read in data */
	for (i = 0; i < N; ++i) {
		sptr = fgets(text, TXTLEN, iou);
		if (NULL == sptr) printf("trouble reading input file, line %d\n", i);
		sptr = strtok(text, ",");
		if (NULL == sptr) printf("nothing to read, line %d\n", i);
		ferr = sscanf(sptr, " %d", &y[i]);
		if (0 == ferr) printf("trouble scanning class label, line %d\n", i);
		for (j = 0; j < P; ++j) {
			sptr = strtok(NULL, ",");
			if (NULL == sptr) printf("nothing to read, line %d, position %d\n", i, j);
			ferr = sscanf(sptr, " %f", &x[i][j]);
			if (0 == ferr) printf("trouble scanning data,line %d, position %d\n", i, j);
		}
	}
	ferr = fclose(iou);
	if (0 != ferr) printf("trouble closing input file\n");
/* Call the subroutine */
	printf("on the raw data:\n");
	k = 3;
	ifault = cluela (x, nx, N, P, w, c, nc, kmin, &k, KMAX,
	z, pop, f, robust, sosdo, iwork, rwork, tree, crit, &u);
	printf("cluela returns #%d, Made %d clusters\n", ifault, k);
/* accuracy calculation */
	accari = arind (y, z, N, G, k, mval, aval, bval);
	accvi = vinfo(y, z, N, G, k, mval, aval, bval);
	printf("residual: %f, ARI: %f, VI: %f\n", u, accari, accvi);
/* Standardize:  attempt to replicate Pham, Dimov & Nguyen */
	avedev (x, nx, N, P, na, ave, dev);
	for (j = 0; j < P; ++j) {
		if (safdiv(1.0f, dev[j])) {
			dev[j] = 1.f / dev[j];
		} else {
			dev[j] = 0.f;
		}
	}
	for (i = 0; i < N; ++i) {
		for (j = 0; j < P; ++j) {
			x[i][j] = (x[i][j] - ave[j]) * dev[j];
		}
	}
/* Call the subroutine */
	k = 3;
	printf("on standardized data:\n");
	ifault = cluela (x, nx, N, P, w, c, nc, kmin, &k, KMAX,
	  z, pop, f, robust, sosdo, iwork, rwork, tree, crit, &u);
	printf("cluela returns #%d, Made %d clusters\n", ifault, k);
/* accuracy calculation */
	accari = arind (y, z, N, G, k, mval, aval, bval);
	accvi = vinfo (y, z, N, G, k, mval, aval, bval);
	printf("residual: %f, ARI: %f, VI: %f\n", u, accari, accvi);
/* write F values */
	strncpy(fname, "wine_fvals_stand.csv", FNLEN);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("wine: trouble opening %s\n", fname);
	for (l = kmin; l <= KMAX; ++l) {
		ferr = fprintf(iou, "%i, %9.3g\n", l, f[l-kmin]);
		if (0 == ferr) printf("wine: trouble writing F stats\n");
	}
	ferr = fclose(iou);
	if (0 != ferr) printf("wine: trouble closing %s\n", fname);
/* Now using specially obtained importance weights */
	k = 3;
	printf("with special weights...\n");
	ifault = cluela (x, nx, N, P, wspec, c, nc, kmin, &k, KMAX,
	  z, pop, f, robust, sosdo, iwork, rwork, tree, crit, &u);
	printf("cluela returns #%d, Made %d clusters\n", ifault, k);
/* accuracy calculation */
	accari = arind (y, z, N, G, k, mval, aval, bval);
	accvi = vinfo (y, z, N, G, k, mval, aval, bval);
	printf("residual: %f, ARI: %f, VI: %f\n", u, accari, accvi);
/* write F values */
	strncpy(fname, "wine_fvals_weighted.csv", FNLEN);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("wine: trouble opening %s\n", fname);
	for (l = kmin; l <= KMAX; ++l) {
		ferr = fprintf(iou, "%i, %9.3g\n", l, f[l-kmin]);
		if (0 == ferr) printf("wine: trouble writing F stats\n");
	}
	ferr = fclose(iou);
	if (0 != ferr) printf("trouble closing file %s\n", fname);
/* project to 2-D */
	for (i = 0; i < N; ++i) {
		for (j = 0; j < P; ++j) {
			x[i][j] *= wspec[j];
		}
	}
	strncpy(fname, "wine_classes_weighted.csv", FNLEN);
	ifault = pca (x, N, P, robust, r, M, &q, ave, e, v, dev, rwork);
	if (ifault != 0) printf("pca: error #%d\n", ifault);
	if (q < M) printf("only %i dimensions!", q);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("wine: trouble opening %s\n", fname);
	ifault = wcsvsp (r, N, y, iou);
	if (ifault != 0) printf ("wcsvsp: error #%d\n", ifault);
	ferr = fclose(iou);
	if (0 != ferr) printf("wine: trouble closing file %s\n", fname);
	strncpy(fname, "wine_classes_robust.csv", FNLEN);
	robust = true;
	ifault = pca (x, N, P, robust, r, M, &q, ave, e, v, dev, rwork);
	if (ifault != 0) printf("pca: error #%d\n", ifault);
	if (q < M) printf("only %i dimensions!", q);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("wine: trouble opening %s\n", fname);
	ifault = wcsvsp (r, N, y, iou);
	if (ifault != 0) printf ("wcsvsp: error #%d\n", ifault);
	ferr = fclose(iou);
	if (0 != ferr) printf("wine: trouble closing file %s\n", fname);
L90: /* release memory */
	if (nc) del_int_mtx(nc, N);
	if (nx) del_int_mtx(nx, N);
	if (tree) del_int_mtx(tree, 2*KMAX);
	if (c) del_flt_mtx(c, N);
	if (crit) del_flt_mtx(crit, 2*KMAX);
	if (mval) del_flt_mtx(mval, G);
	if (r) del_flt_mtx(r, N);
	if (v) del_flt_mtx(v, P);
	if (x) del_flt_mtx(x, N);
	printf("program complete\n");
	return 0;
} /* end of main */
