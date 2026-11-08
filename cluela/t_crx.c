/*======================================================================
test clustering on the credit approval data
by Andy Allinger, 2017, public domain
This program may be used by any person for any purpose.
======================================================================*/
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

int **new_int_mtx (int, int); /* defined in array.c */
float **new_flt_mtx (int, int);
char **new_ch_mtx (int, int);
int del_int_mtx (int **, int);
int del_flt_mtx (float **, int);
int del_ch_mtx (char **, int);
bool safdiv (float numer, float denom);  /* defined in safdiv.c */
void leagre (float **x, int **nx, int n, int p,  /* defined in stats.c */
             int *na, float *lea, float *gre);
int ireid (char **h, int *id, int *nx, int *ind, int *jnd, int n, int h_len, char *hwork);  /* defined in prep.c */
int rscode (float **x, int **nx, int n, int p, float *w,
            int *cat, int o, float **s, int catmax, bool robust);
int mds (float **x, int **nx, int n, int p, float *w, int *cat, bool robust,  /* defined in mds.c */
         float **v, int *m, float *e, float *t);
void dinit (float **x, int **nx, int n, int p, float *w, float **c, int **nc,  /* defined in cluela.c */
           int k, bool robust, float *shortx, float *cumd, float *thresh);
int kclust (float **x, int **nx, int n, int p, float *w, float **c, int **nc,
  int kmin, int *k, int kmax, int *z, int *pop, bool miss, bool robust, bool sosdo,
  float thresh, int *livec, int *livex, int *iperm, int **tree, float **crit,
  float *shortc, float *shortx, float *work, float *u);
float arind (int *y, int *z, int n, int g, int k, float **m, float *a, float *b);  /* defined in valid.c */
float vinfo (int *y, int *z, int n, int g, int k, float **m, float *a, float *b);
int wcsvsp (float **, int, int *, FILE *); /* defined in dump.c */

/* Main program */
int main (void)
{
/* Constants */
#define CATMAX 14  /* most categories any attribute */
#define FNLEN 32   /* length of a file name */
#define G 2        /* # correct classes */
#define HLEN 32    /* length of a string */
#define KMIN 1     /* fewest clusters to make */
#define KDEF 2     /* default clusters to make */
#define KMAX 20    /* most clusters to make */
#define LBLLEN 1   /* length of the class label */
#define M 2        /* dimensions to embed into */
#define N 690      /* number of objects */
#define NO 10350   /* N @ O */
#define NS 12      /* node size of binary tree */
#define O 15       /* input attributes */
#define P 38       /* attributes after encoding */
#define SYMLEN 2   /* length of symbol strings */
#define TXTLEN 256 /* length of line of input file */

	const float fone = 1.f;
	const float xbogus = -9.f; /* arbitrary value for missing data */

/* Local variables */
	int i, j, k, l, q, y[N], z[N], ll, na[P], nx0[N], pop[KMAX], jbeg, ferr,
	  cplen, iperm[N], jperm[N], ifault, livex[N], livec[KMAX];
	float aval[G], bval[KMAX], e[N], lea[P], gre[P], r,
	  shortx[N], shortc[KMAX], thresh, u, wbak[P], work[N];
	float accvi = -1.f,
	      accari = -1.f;
	bool miss, robust, sosdo;
	char ch1[LBLLEN], hwork[SYMLEN], empty[HLEN], raw[HLEN], text[TXTLEN];
	char *sptr;
	FILE *iou;

/* Allocate matrices */
	int **nc, **nx, **nxbak, **tree;
	float **c, **crit, **mval, **s, **v, **x, **xbak;
	char **clslbl, **symbol;
	nc = new_int_mtx (KMAX, P);
	nx = new_int_mtx (N, P);
	nxbak = new_int_mtx (N, P);
	tree = new_int_mtx (2*KMAX, NS);
	c = new_flt_mtx (KMAX, P);
	crit = new_flt_mtx (2*KMAX, 4);
	mval = new_flt_mtx (G, KMAX);
	s = new_flt_mtx (CATMAX, CATMAX-1);
	v = new_flt_mtx (N, N);
	x = new_flt_mtx (N, P);
	xbak = new_flt_mtx (N, P);
	clslbl = new_ch_mtx (N, LBLLEN);
	symbol = new_ch_mtx (N*O, SYMLEN);
	if (!nc || !nx || !nxbak || !tree || !c || !crit || !mval || !s ||
	    !v || !x || !xbak || !clslbl || !symbol) goto L90;

	char fname[FNLEN] = "crx.data";
/*======================================================================
Credit card application data
From the University of California at Irvine Machine Learning repository

Reference:
     Ross Quinlan
     "Simplifying decision trees"
     Int J Man-Machine Studies v.27, pp.221-234
     Dec 1987

This data was provided by an anonymous donor.  The attributes have been
changed to arbitrary values.

Classes are:
   +, -

Attributes are:

___Attribute___Type___________#___Values_______________________________
   1          Categorical    2   b, a
   2          Continuous
   3          Continuous
   4          Categorical    4   u, y, l, t
   5          Categorical    3   g, p, gg
   6          Categorical   14   c, d, cc, i, j, k, m, r, q, w, x, e, aa, ff
   7          Categorical    9   v, h, bb, j, n, z, dd, ff, o
   8          Continuous
   9          Boolean        2   t, f
  10          Boolean        2   t, f
  11          Continuous
  12          Boolean        2   t, f
  13          Categorical    3   g, p, s
  14          Continuous
  15          Continuous
======================================================================*/

/* Info above is summarized thus: */
	static int cat[O] = { 2,0,0,4,3,14,9,0,2,2,0,2,3,0,0 };
/* Equal weights */
	float w[P] = { 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f,
	               1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f,
	               1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f,
	               1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f };
/* bogus initial value */
	for (i = 0; i < NO; ++i) {
		for (j = 0; j < SYMLEN; ++j) {
			symbol[i][j] = '$';
		}
	}
	for (i = 0; i < HLEN; ++i) empty[i] = '\0';
/* initialize random number generator */
	srand( (unsigned int)time(NULL));
/* open input file */
	iou = fopen(fname, "r");
	if (NULL == iou) printf("crx: trouble opening input file\n");
/* read in data */
	for (i = 0; i < N; ++i) {
		sptr = fgets(text, TXTLEN, iou);
		if (NULL == sptr) printf("trouble reading from file, line %d\n", i);
		ll = (int) strlen(text);
		jbeg = 0;
		l = 0;
/* split on commas */
		for (j = 0; j < ll; ++j) {
			if (',' == text[j]) {
				++l;
				strncpy(raw, empty, HLEN);
				cplen = j - jbeg;
				if (cplen > HLEN) cplen = HLEN;
				strncpy(raw, &text[jbeg], cplen);
			} else if (j == ll-1) {
				++l;
				strncpy(raw, empty, HLEN);
				cplen = j - jbeg + 1;
				if (cplen > HLEN) cplen = HLEN;
				strncpy(raw, &text[jbeg], cplen);
			} else {
				continue;
			}
			if (l > O+1) printf("trouble reading l > O+1\n");
			if (16 == l) { /* answer key */
				clslbl[i][0] = raw[0];
			} else if (NULL == strstr(raw, "?")) { /* non-missing */
				if (0 == cat[l-1]) { /* numerical? */
					ferr = sscanf(raw, " %f", &x[i][l-1]);
					if (0 == ferr) printf("error scanning X, i: %d, l: %d, string: %s \n", i, l, raw);
				} else {
					symbol[(l-1)*N+i][0] = raw[0];
					symbol[(l-1)*N+i][1] = raw[1];
				}
				nx[i][l-1] = 1;
			} else { /* missing */
				x[i][l-1] = xbogus;
				nx[i][l-1] = 0;
			}
			jbeg = j + 1;
		}
	}
	ifault = fclose(iou);
	if (0 != ifault) printf("trouble closing file %s\n", fname);
/* preprocess input */
	for (l = 0; l < O; ++l) { /* convert text labels to integers */
		if (cat[l] > 0) {
			for (i = 0; i < N; ++i) nx0[i] = nx[i][l];
			ifault = ireid (&symbol[l*N], z, nx0, iperm, jperm, N, SYMLEN, hwork);
			if (ifault != 0) printf("ireid: error#%d\n", ifault);
			for (i = 0; i < N; ++i) x[i][l] = (float) z[i];
		}
	}
/* Variable 15 has a standard deviation greater than its semi-interquartile range
by a factor of 26.  Apply a logarithm. */
	for (i = 0; i < N; ++i) x[i][14] = logf(1.f + x[i][14]);
/* rescale by range */
	leagre (x, nx, N, P, na, lea, gre);
	for (l = 0; l < O; ++l) {
		if (0 == cat[l]) {
			r = gre[l] - lea[l];
			if (safdiv(fone, r)) {
				r = 1. / r;
			} else {
				r = 0.f;
				printf("crx: no variance column%i\n", l);
			}
			for (i = 0; i < N; ++i) x[i][l] = (x[i][l] - lea[l]) * r;
		}
	}
/* back up input */
	for (i = 0; i < N; ++i) {
		for (l = 0; l < O; ++l) {
			xbak[i][l] = x[i][l];
			nxbak[i][l] = nx[i][l];
		}
	}
	for (l = 0; l < O; ++l) wbak[l] = w[l];
/* recode categoricals with corresponding column of simplex */
	robust = false;
	ifault = rscode (x, nx, N, P, w, cat, O, s, CATMAX, robust);
	if (ifault != 0) printf("rscode: error #%d\n", ifault);
/* relabel class ID's */
	for (i = 0; i < N; ++i) nx0[i] = 1;
	ifault = ireid (clslbl, y, nx0, iperm, jperm, N, LBLLEN, ch1);
	if (ifault != 0) printf("ireid on class, error #%d\n", ifault);
/* call the subroutine */
	k = KDEF;
	sosdo = true;
	miss = true;
	dinit (x, nx, N, P, w, c, nc, k, robust, shortx, work, &thresh);
	ifault = kclust (x, nx, N, P, w, c, nc, KMIN, &k, KMAX,
	  z, pop, miss, robust, sosdo, thresh,
	  livec, livex, iperm, tree, crit, shortc, shortx, work, &u);
	printf("kclust returns #%d, made %d clusters\n", ifault, k);
/* repeat with half the threshhold */
	dinit (x, nx, N, P, w, c, nc, k, robust, shortx, work, &thresh);
	thresh /= 2.f;
	ifault = kclust (x, nx, N, P, w, c, nc, KMIN, &k, KMAX,
	  z, pop, miss, robust, sosdo, thresh,
	  livec, livex, iperm, tree, crit, shortc, shortx, work, &u);
	printf("with half the threshhold distance, kclust returns #%d, made %d clusters\n", ifault, k);
/* find accuracy */
	if (ifault == 0) {
		accari = arind (y, z, N, G, k, mval, aval, bval); /* adjusted Rand index */
		accvi = vinfo (y, z, N, G, k, mval, aval, bval); /* variation of information */
	}
/* project to 2-D */
	ifault = mds (xbak, nxbak, N, P, wbak, cat, robust, v, &q, e, work);
	if (ifault != 0) printf("mds: error#%d\n", ifault);
	if (q < M) printf("mds: only %d directions!\n", q);
/* Write CSV file with coordinates in different columns by cluster */
	strncpy(fname, "crx_clusters.csv", FNLEN);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("crx: trouble opening %s\n", fname);
	ifault = wcsvsp (v, N, z, iou);
	if (ifault != 0) printf("wcsvsp: error#%d\n", ifault);
	ferr = fclose(iou);
	if (0 != ferr) printf("trouble closing file %s\n", fname);
/* Write CSV file with coordinates in different columns by class */
	strncpy(fname, "crx_classes_euclid.csv", FNLEN);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("crx: trouble opening %s\n", fname);
	ifault = wcsvsp (v, N, y, iou);
	if (0 != ifault) printf("wcsvsp: error#%d\n", ifault);
	ferr = fclose(iou);
	if (0 != ferr) printf("trouble closing file %s\n", fname);
/* Robust projection */
	robust = true;
	ifault = mds (xbak, nxbak, N, P, wbak, cat, robust, v, &q, e, work);
	if (ifault != 0) printf("mds: error#%d\n", ifault);
	if (q < M) printf("mds: only %d directions!\n", q);
	strncpy(fname, "crx_classes_robust.csv", FNLEN);
	iou = fopen(fname, "w");
	if (ifault != 0) printf("crx: trouble opening %s\n", fname);
	ifault = wcsvsp (v, N, y, iou);
	if (ifault != 0) printf("wcsvsp: error#%d\n", ifault);
	ferr = fclose(iou);
	if (0 != ferr) printf("trouble closing file %s\n", fname);
/* Print results */
	printf("Returned residual:  %f\n", u);
	printf("Adjusted Rand Index:  %f\n", accari);
	printf("Variation of Information:  %f\n", accvi);
L90: /* release memory */
	if (nc) del_int_mtx (nc, KMAX);
	if (nx) del_int_mtx (nx, N);
	if (nxbak) del_int_mtx (nxbak, N);
	if (tree) del_int_mtx (tree, 2*KMAX);
	if (c) del_flt_mtx (c, KMAX);
	if (crit) del_flt_mtx (crit, 2*KMAX);
	if (mval) del_flt_mtx (mval, G);
	if (s) del_flt_mtx (s, CATMAX);
	if (v) del_flt_mtx (v, P);
	if (x) del_flt_mtx (x, N);
	if (xbak) del_flt_mtx (xbak, N);
	if (clslbl) del_ch_mtx (clslbl, N);
	if (symbol) del_ch_mtx (symbol, N*O);
	printf("program complete\n");
	return 0;
} /* end of main */
