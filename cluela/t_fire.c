/*----------------------------------------------------------------------
t_fire.c - Cluster the forest fire data
by Andy Allinger, 2016-2020, public domain
----------------------------------------------------------------------*/
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>

int **new_int_mtx (int m, int n);  /* defined in array.c */
float **new_flt_mtx (int m, int n);
int del_int_mtx (int **b, int m);
int del_flt_mtx (float **x, int m);
void leagre (float **x, int **nx, int n, int p,  /* defined in stats.c */
             int *na, float *lea, float *gre);
int mds (float **x, int **nx, int n, int p, float *w, int *cat, bool robust,  /* defined in mds.c */
         float **v, int *m, float *e, float *t);
int cluela (float **x, int **nx, int n, int p, float *w, float **c, int **nc,  /* defined in cluela.c */
      int kmin, int *k, int kmax, int *z, int *pop, float *f, bool robust,
      bool sosdo, int *iwork, float *rwork, int **tree, float **crit, float *u);
float arind (int *y, int *z, int n, int g, int k,  /* defined in valid.c */
            float **m, float *a, float *b);
float vinfo (int *y, int *z, int n, int g, int k,
            float **m, float *a, float *b);
int wcsvsp (float **r, int n, int *z, FILE *iou);  /* defined in dump.c */

/* fixed parameters */
#define HLEN 32     /* length of filename */
#define P 13        /* # attributes */
#define N 517       /* # objects */
#define NS 12       /* node size of binary tree */
#define PLIM N      /* so embedding will fit in same arrays */
#define TXTLEN 256  /* length of a line of input */

/* adjustable parameters */
#define KMIN 1   /* fewest clusters to make */
#define KDEF 3   /* default clusters to make */
#define KMAX 10  /* most clusters to make */

int main(void)
{
/* Initialized data */
	bool robust = false;
/* Local variables */
	int i, j, l, m, ifault, iwork[4*N+KMAX], kembed, korig, kpse,
	    na[P], pop[KMAX], ze[N], zo[N], zs[N];
	float a[KMAX], b[KMAX], accari, accvi, e[N], f[KMAX-KMIN+1],
	      lea[P], gre[P], rwork[2*N+KMAX], t[N], u, w[PLIM];
	bool ok, sosdo;
	char day[4], month[4];  /* three characters plus a null byte */
	char *sptr, text[TXTLEN];
	FILE *iou;
/* allocate matrices */
	int **nc, **nx, **tree;
	float **c, **confus, **crit, **v, **x;
	nc = new_int_mtx(KMAX, PLIM);
	nx = new_int_mtx(N, PLIM);
	tree = new_int_mtx(2*KMAX, NS);
	c = new_flt_mtx(KMAX, PLIM);
	confus = new_flt_mtx(KMAX, KMAX);
	crit = new_flt_mtx(2*KMAX, 4);
	v = new_flt_mtx(N, N);
	x = new_flt_mtx(N, PLIM);
	if (!nc || !nx || !tree || !c || !confus || !crit || !v || !x) goto L90;
/* the data */
	char fname[HLEN] = "forestfires.csv";
/*----------------------------------------------------------------------
The Forest Fire data is described in:

     Paulo Cortez and Anibal Morais
     "A Data Mining Approach to Predict Forest Fires using Meteorological Data"
     In J. Neves, M. F. Santos and J. Machado Eds.,
     New Trends in Artificial Intelligence,
     In Proceedings of the 13th EPIA 2007
     Portuguese Conference on Artificial Intelligence,
     Guimaraes, Portugal, pp. 512-523, APPIA, ISBN-13 978-989-95618-0-9.
     December, 2007. (http://www.dsi.uminho.pt/~pcortez/fires.pdf)

Number of Instances: 517
Number of Attributes: 12 + output attribute

____#___Attribute_________________Type_________Values__________________
    0   spacial coordinate X      Continuous   1 ... 9
    1   spacial coordinate Y      Continuous   2 ... 9
    2   month                     Cyclic       "jan" ... "dec"
    3   day                       Cyclic       "mon" ... "sun"
    4   Fine Fuel Moisture Code   Continuous   18.7 ... 96.20
    5   Duff Moisture Code DMC    Continuous   1.1 ... 291.3
    6   Drought Code              Continuous   7.9 ... 860.6
    7   Initial Spread Index      Continuous   0.0 ... 56.10
    8   temperature, Celsius      Continuous   2.2 ... 33.30
    9   relative humidity         Continuous   15.0 ... 100
   10   wind speed, km/h          Continuous   0.40 ... 9.40
   11   recent rain, mm           Continuous   0.0 ... 6.4
   12   burned area, ha           Continuous   0.0 ... 1090.84
----------------------------------------------------------------------*/
	int cat[P] = { 0,0,-1,-1,0,0,0,0,0,0,0,0,0 };
/* Begin. */
	srand(time(NULL));
/*get input file */
	iou = fopen(fname, "r");
	if (NULL == iou) printf("trouble opening input file: %s\n", fname);
/* Read in data. */
	sptr = fgets(text, TXTLEN, iou);  /* burn header line */
	if (NULL == sptr) printf("trouble reading header line\n");
	for (i = 0; i < N; ++i) {
		sptr = fgets(text, TXTLEN, iou);
		if (NULL == sptr) printf("trouble reading line %i\n", i+1);
		ifault = sscanf(text, "%f,%f,%3s,%3s,%f,%f,%f,%f,%f,%f,%f,%f,%f",
		    &x[i][0], &x[i][1], month, day, &x[i][4], &x[i][5], &x[i][6],
		    &x[i][7], &x[i][8], &x[i][9], &x[i][10], &x[i][11], &x[i][12]);
		if (ifault < 0) printf("trouble scanning line %i\n", i+1);
		if (0 == strncmp(month, "jan", 3)) {
			x[i][2] = 1.f;
		} else if (0 == strncmp(month, "feb", 3)) {
			x[i][2] = 2.f;
		} else if (0 == strncmp(month, "mar", 3)) {
			x[i][2] = 3.f;
		} else if (0 == strncmp(month, "apr", 3)) {
			x[i][2] = 4.f;
		} else if (0 == strncmp(month, "may", 3)) {
			x[i][2] = 5.f;
		} else if (0 == strncmp(month, "jun", 3)) {
			x[i][2] = 6.f;
		} else if (0 == strncmp(month, "jul", 3)) {
			x[i][2] = 7.f;
		} else if (0 == strncmp(month, "aug", 3)) {
			x[i][2] = 8.f;
		} else if (0 == strncmp(month, "sep", 3)) {
			x[i][2] = 9.f;
		} else if (0 == strncmp(month, "oct", 3)) {
			x[i][2] = 10.f;
		} else if (0 == strncmp(month, "nov", 3)) {
			x[i][2] = 11.f;
		} else if (0 == strncmp(month, "dec", 3)) {
			x[i][2] = 12.f;
		} else {
			printf("fire: bad month %s\n", month);
			exit(EXIT_FAILURE);
		}
		if (0 == strncmp(day, "sun", 3)) {
			x[i][3] = 1.f;
		} else if (0 == strncmp(day, "mon", 3)) {
			x[i][3] = 2.f;
		} else if (0 == strncmp(day, "tue", 3)) {
			x[i][3] = 3.f;
		} else if (0 == strncmp(day, "wed", 3)) {
			x[i][3] = 4.f;
		} else if (0 == strncmp(day, "thu", 3)) {
			x[i][3] = 5.f;
		} else if (0 == strncmp(day, "fri", 3)) {
			x[i][3] = 6.f;
		} else if (0 == strncmp(day, "sat", 3)) {
			x[i][3] = 7.f;
		} else {
			printf("fire: bad day %s\n", day);
			exit(EXIT_FAILURE);
		}
	}
	fclose(iou);
	for (i = 0; i < N; ++i) for (j = 0; j < PLIM; ++j) nx[i][j] = 1;  /* nothing missing */
/*----------------------------------------------------------------------
Take the logarithm of skewed variable "area burned"
Scale month and day to unit period.  Scale the other variables to unit range.
Then embed.
----------------------------------------------------------------------*/
	for (i = 0; i < N; ++i) {
		x[i][12] = log(x[i][12] + 1.f);  /* recommended by the authors */
	}
	leagre (x, nx, N, P, na, lea, gre);
	for (i = 0; i < N; ++i) {
		for (j = 0; j < 2; ++j) x[i][j] = (x[i][j] - lea[j]) / (gre[j] - lea[j]);
		x[i][2] /= 12.f;
		x[i][3] /= 7.f;
		for (j = 4; j < P; ++j) x[i][j] = (x[i][j] - lea[j]) / (gre[j] - lea[j]);
	}
	for (j = 0; j < PLIM; ++j) w[j] = 1.f;  /* equal weights */
	ifault = mds (x, nx, N, P, w, cat, robust, v, &m, e, t);
	printf("mds returns #%d with %i directions.\n", ifault, m);
/* Cluster the original data in which X is interpreted as continuous. */
	ok = true;
	korig = KDEF;
	sosdo = ! robust;
	ifault = cluela (x, nx, N, P, w, c, nc, KMIN, &korig, KMAX, zo, pop, f,
	                robust, sosdo, iwork, rwork, tree, crit, &u);
	printf("original problem:  cluela returns #%i made %i clusters\n", ifault, korig);
	ok = ok && (0 == ifault);
/* Cluster the embedded data, real directions only */
	for (i = 0; i < N; ++i) for (j = 0; j < N; ++j) x[i][j] = v[i][j];  /* copy V into X */
	kembed = KDEF;
	robust = false ;  /* don't use robust mode on embedded data */
	sosdo = true;
	ifault = cluela (x, nx, N, m, w, c, nc, KMIN, &kembed, KMAX, ze, pop, f,
	                 robust, sosdo, iwork, rwork, tree, crit, &u);
	printf("embedded problem:  cluela returns #%i made %i clusters\n", ifault, kembed);
	ok = ok && (0 == ifault);
/* Cluster the embedded data, pseudo-Euclidean space */
	for (j = m; j < N; ++j) w[j] = -1.f;  /* imaginary! */
	kpse = KDEF;
	ifault = cluela (x, nx, N, N, w, c, nc, KMIN, &kpse, KMAX, zs, pop, f,
	                 robust, sosdo, iwork, rwork, tree, crit, &u);
	printf("pseudo-Euclidean problem: cluela returns #%i made %i clusters\n", ifault, kpse);
	ok = ok && (0 == ifault);
/* compare */
	if (ok) {
		accari = arind (zo, zs, N, korig, kpse, confus, a, b);
		accvi = vinfo (zo, zs, N, korig, kpse, confus, a, b);
		printf("Comparing the linear against pseudo-Euclidean:  "
		        "Adjusted Rand Index: %f Variation of Information: %f\n", accari, accvi);
		accari = arind (ze, zs, N, kembed, kpse, confus, a, b);
		accvi = vinfo (ze, zs, N, kembed, kpse, confus, a, b);
		printf("Comparing the real embedding against pseudo-Euclidean:  "
		        "Adjusted Rand Index: %f Variation of Information: %f\n", accari, accvi);
	}
/* Write CSV file */
	strncpy(fname, "fire_clusters.csv", HLEN);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("fire: trouble opening %s\n", fname);
	ifault = wcsvsp(v, N, zs, iou);
	if (0 != ifault) printf("wcsvsp: error #%i\n", ifault);
	fclose(iou);
/* write F values */
	strncpy(fname, "fire_fvals.csv", HLEN);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("fire: trouble opening %s\n", fname);
	for (l = 1; l <= KMAX; ++l) {
		ifault = fprintf(iou, "%5i,%9.3g\n", l, (double)f[l-KMIN]);
		if (ifault < 0) printf("fire: trouble writing F stats\n");
	}
	fclose(iou);
     printf("program complete\n");
L90:  /* release memory */
	if (nc) del_int_mtx(nc, KMAX);
	if (nx) del_int_mtx(nx, N);
	if (tree) del_int_mtx(tree, 2*KMAX);
	if (c) del_flt_mtx(c, KMAX);
	if (confus) del_flt_mtx(confus, KMAX);
	if (crit) del_flt_mtx(crit, 2*KMAX);
	if (x) del_flt_mtx(x, N);
	if (v) del_flt_mtx(v, N);
	return 0;
} /* end of main */
