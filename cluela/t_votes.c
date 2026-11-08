/*----------------------------------------------------------------------
test clustering on the House votes '84 data set
by Andy Allinger, 2011-2017, public domain
This program may be used by any person for any purpose.
----------------------------------------------------------------------*/
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
int cluela (float **x, int **nx, int n, int p, float *w, float **c, int **nc,  /* defined in cluela.c */
      int kmin, int *k, int kmax, int *z, int *pop, float *f, bool robust,
      bool sosdo, int *iwork, float *rwork, int **tree, float **crit, float *u);
int pca (float **x, int n, int p, bool robust, float **r, int m, int *q,  /* defined in mds.c */
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
#define G 2        /* # political parties */
#define M 2        /* projection dimension for graphs */
#define N 435      /* # of Representatives */
#define NS 12      /* node size of binary tree */
#define P 16       /* # of bills */
#define KMIN 1     /* fewest clusters to make */
#define KMAX 20    /* most clusters to make */
#define TXTLEN 256 /* length of input buffer */
	const float xbogus = -9.f;
/* Local variables */
	int i, j, k, l, q, y[N], z[N], ll, iwork[4*N+KMAX], pop[KMAX], ferr, ifault;
	float a[P],  e[P], f[KMAX-KMIN+1], t[P], u, rwork[2*N+KMAX], aval[G], bval[KMAX];
	float accari = -1.f,
	  accvi = -1.f;
	bool robust, sosdo;
	char text[TXTLEN];
	char *sptr;
	FILE *iou;
/* Allocate matrices */
	int **nc, **nx, **tree;
	float **c, **crit, **mval, **r, **v, **x;
	char **raw;
	nc = new_int_mtx(N, KMAX);
	nx = new_int_mtx(N, P);
	tree = new_int_mtx(2*KMAX, NS);
	c = new_flt_mtx(N, KMAX);
	crit = new_flt_mtx(2*KMAX, 4);
	mval = new_flt_mtx(G, KMAX);
	r = new_flt_mtx(N, M);
	v = new_flt_mtx(P, P);
	x = new_flt_mtx(N, P);
	raw = new_ch_mtx(N, P);
	if (!nc || !nx || !tree || !c || !crit || !mval || !r || !x || !raw) goto L90;

	static char fname[FNLEN] = "house-votes-84.data";
/*----------------------------------------------------------------------
1984 United States Congressional Voting Records Database
From the University of California at Irvine Machine Learning repository

Source:
     Congressional Quarterly Almanac, 98th Congress,
     2nd session 1984, Volume XL: Congressional Quarterly Inc.
     Washington, D.C., 1985.

Reference:
     J. C. Schlimmer
     Concept acquisition through representational adjustment
     Doctoral dissertation, Department of Information and Computer Science,
     University of California, Irvine, CA., 1987.

Classes are:
    democrat, republican

Attributes are:

___Attribute________________________________Type______#___Values_______
   handicapped-infants                      Boolean   2   y,n,?
   water-project-cost-sharing               Boolean   2   y,n,?
   adoption-of-the-budget-resolution        Boolean   2   y,n,?
   physician-fee-freeze                     Boolean   2   y,n,?
   el-salvador-aid                          Boolean   2   y,n,?
   religious-groups-in-schools              Boolean   2   y,n,?
   anti-satellite-test-ban                  Boolean   2   y,n,?
   aid-to-nicaraguan-contras                Boolean   2   y,n,?
   mx-missile                               Boolean   2   y,n,?
   immigration                              Boolean   2   y,n,?
   synfuels-corporation-cutback             Boolean   2   y,n,?
   education-spending                       Boolean   2   y,n,?
   superfund-right-to-sue                   Boolean   2   y,n,?
   crime                                    Boolean   2   y,n,?
   duty-free-exports                        Boolean   2   y,n,?
   export-administration-act-south-africa   Boolean   2   y,n,?
----------------------------------------------------------------------*/

	static float w[P] = { 1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,
	  1.f,1.f,1.f };
	srand((unsigned int)time(NULL));
/* open file */
	iou = fopen(fname, "r");
	if (NULL == iou) {
		printf("trouble opening input file\n");
		return EXIT_FAILURE;
	}
/* Read in data */
	for (i = 0; i < N; ++i) {
		sptr = fgets(text, TXTLEN, iou);
		if (NULL == sptr) printf("trouble reading line %d\n", i);
		if ('d' == text[0]) y[i] = 1; /* class labels */
		if ('r' == text[0]) y[i] = 2;
		ll= (int) strlen(text);
		l = 0;
		for (j = 0; j < ll-1; ++j) {
			if (',' == text[j]) {
				raw[i][l] = text[j+1];
				++l;
			}
		}
	}
	ferr = fclose(iou);
	if (0 != ferr) printf("trouble closing file %s", fname);

/*----------------------------------------------------------------------
Consider the problem as binary data with some values missing.
----------------------------------------------------------------------*/
	for (i = 0; i < N; ++i) {
		for (j = 0; j < P; ++j) {
			if ('y' == raw[i][j]) {
				x[i][j] = 1.f;
				nx[i][j] = 1;
			} else if ('n' == raw[i][j]) {
				x[i][j] = 0.f;
				nx[i][j] = 1;
			} else if ('?' == raw[i][j]) {
				x[i][j] = xbogus;
				nx[i][j] = 0;
			} else {
				printf("tvotes: trouble recoding input\n");
				return EXIT_FAILURE;
			}
		}
	}
	k = 3;
	robust = false;
	sosdo = true;
	ifault = cluela (x, nx, N, P, w, c, nc, KMIN, &k, KMAX, z, pop, f,
	  robust, sosdo, iwork, rwork, tree, crit, &u);
	printf("cluela returns #%d, made %d clusters\n", ifault, k);

/*----------------------------------------------------------------------
Reconsider the problem as scalar data, with no missing values
----------------------------------------------------------------------*/
	for (i = 0; i < N; ++i) {
		for (j = 0; j < P; ++j) {
			if ('y' == raw[i][j]) x[i][j] = 1.f;
			if ('n' == raw[i][j]) x[i][j] = 0.f;
			if ('?' == raw[i][j]) x[i][j] = 0.5f;
			nx[i][j] = 1;
		}
	}
	k = 3;
	ifault = cluela (x, nx, N, P, w, c, nc, KMIN, &k, KMAX, z, pop, f,
	  robust, sosdo, iwork, rwork, tree, crit, &u);
	printf("cluela returns #%d, made %d clusters\n", ifault, k);
/* Find accuracy */
	if (0 == ifault) {
		accari = arind (y, z, N, G, k, mval, aval, bval); /* adjusted Rand index */
		accvi = vinfo (y, z, N, G, k, mval, aval, bval); /* variation of information */
	}
/* project to 2-D */
	ifault = pca (x, N, P, robust, r, M, &q, a, e, v, t, rwork);
	if (ifault != 0) printf("pca: error#%d\n", ifault);
	if (q < M) printf("pca: only %d directions!\n", q);
/* Write output */
	strncpy(fname, "votes_clusters_sosdo.csv", FNLEN);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("tvotes: trouble opening %s\n", fname);
	ifault = wcsvsp (r, N, z, iou);
	if (ifault != 0) printf("wcsvsp: error #%d\n", ifault);
	ferr = fclose(iou);
	if (0 != ferr) printf("trouble closing file %s\n", fname);
	strncpy(fname, "votes_classes.csv", FNLEN);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("tvotes: trouble opening %s\n", fname);
	ifault = wcsvsp (r, N, y, iou);
	if (ifault != 0) printf("wcsvsp: error #%d\n", ifault);
	ferr = fclose(iou);
	if (0 != ferr) printf("trouble closing file %s\n", fname);
/* Write F values */
	strncpy(fname, "votes_fvals.csv", FNLEN);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("tvotes: trouble opening %s\n", fname);
	for (l = KMIN; l <= KMAX; ++l) {
		ferr = fprintf(iou, "%i, %9.3g\n", l, f[l-KMIN]);
		if (0 == ferr) printf("tvotes: trouble writing F stats\n");
	}
	ferr = fclose(iou);
	if (0 != ferr) printf("trouble closing file %s\n", fname);
	printf("Returned residual: %f\n", u);
	printf("Adjusted Rand Index: %f\n", accari);
	printf("Variation of Information: %f\n", accvi);

/*----------------------------------------------------------------------
              Repeat in robust mode
----------------------------------------------------------------------*/
	k = 3;
	robust = true;
	sosdo = false;
	ifault = cluela (x, nx, N, P, w, c, nc, KMIN, &k, KMAX, z, pop, f,
	  robust, sosdo, iwork, rwork, tree, crit, &u);
	printf("cluela returns #%d, made %d clusters\n", ifault, k);
/* Find accuracy */
	if (0 == ifault) {
		accari = arind (y, z, N, G, k, mval, aval, bval); /* adjusted Rand index */
		accvi = vinfo (y, z, N, G, k, mval, aval, bval); /* variation of information */
	}
/* project to 2-D */
	ifault = pca (x, N, P, robust, r, M, &q, a, e, v, t, rwork);
	if (ifault != 0) printf("pca: error#%d\n", ifault);
	if (q < M) printf("pca: only %d directions!\n", q);
/* Write output */
	strncpy(fname, "votes_clusters_robust.csv", FNLEN);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("tvotes: trouble opening %s\n", fname);
	ifault = wcsvsp (r, N, z, iou);
	if (ifault != 0) printf("wcsvsp: error #%d\n", ifault);
	ferr = fclose(iou);
	if (0 != ferr) printf("trouble closing file %s\n", fname);
	printf("Returned residual: %f\n", u);
	printf("Adjusted Rand Index: %f\n", accari);
	printf("Variation of Information: %f\n", accvi);

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
	if (raw) del_ch_mtx(raw, N);
	printf("program complete\n");
	return EXIT_SUCCESS;
} /* end of main */

