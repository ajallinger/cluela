/*======================================================================
test clustering on the Soybean data
by Andy Allinger, 2017, public domain
This program may be used by any person for any purpose.
======================================================================*/
#include <stdbool.h> /* standard C99 */
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

int **new_int_mtx (int, int); /* defined in array.c */
float **new_flt_mtx (int, int);
char **new_ch_mtx (int, int);
int del_int_mtx (int **, int);
int del_flt_mtx (float **, int);
int del_ch_mtx (char **, int);
int mds (float **, int **, int, int, float *, int *, bool, /* defined in mds.c */
         float **, int *, float *, float *);
int ireid (char **, int *, int *, int *, int *, int, int, char *); /* defined in prep.c */
int rscode (float **, int **, int, int, float *, int *, int, float **, int, bool);
int cluela (float **, int **, int, int, float *, float **, int **, int, int *, int, /* defined in cluela.c */
            int *, int *, float *, bool, bool, int *, float *, int **, float **, float *);
int wcsvsp (float **, int, int *, FILE *); /* defined in dump.c */
float arind (int *, int *, int, int, int, float **, float *, float *); /* defined in valid.c */
float vinfo (int *, int *, int, int, int, float **, float *, float *);

int main (void)
{
/* Constants */
#define CATMAX 5   /* most categories any variable */
#define G 19       /* # of correct classes */
#define HLEN 32    /* length of string for label */
#define KMIN 1     /* fewest clusters to make */
#define KDEF G     /* a default number of clusters */
#define KMAX 50    /* most clusters to make */
#define M 2        /* projection dimension for visualization */
#define N 683      /* # of objects */
#define NS 12      /* node size of binary tree */
#define O 35       /* # input variables */
#define P 54       /* # variables after recoding */
#define TXTLEN 256 /* length of string for input text */

	const float xbogus = -9.f; /* arbitrary value for missing data */
	const bool robust = true,  /* use Manhattan distance and median averages */
	           sosdo = false;  /* do not seek sum of squares distance objective */
/* Local variables */
	int i, j, k, l, ll, mc, md, mr, cplen, y[N], z[N], nx0[N], iperm[N], jperm[N],
	    pop[KMAX], jbeg, ferr, iwork[4*N+KMAX], ifault;
	float rwork[2*N+KMAX], accari = 0.f, accvi = 999.f, aval[G], bval[KMAX],
	      t[N], u, e[N], f[KMAX], wv[N], xlea, xgre, dif, worst;
	char raw[HLEN], empty[HLEN], text[TXTLEN], *sptr;
	FILE *iou;

/* allocate matrices */
	int **nx, **nc, **nv, **ncv, **tree;
	float **x, **c, **crit, **cv, **s, **mval, **vc, **vd, **vr;
	char **clslbl;
	nx = new_int_mtx(N, P);
	nc = new_int_mtx(KMAX, P);
	nv = new_int_mtx(N, N);
	ncv = new_int_mtx(KMAX, N);
	tree = new_int_mtx(2*KMAX, NS);
	x = new_flt_mtx(N, P);
	c = new_flt_mtx(KMAX, P);
	crit = new_flt_mtx(2*KMAX, 4);
	cv = new_flt_mtx(KMAX, N);
	s = new_flt_mtx(CATMAX, CATMAX-1);
	mval = new_flt_mtx(G, KMAX);
	vc = new_flt_mtx(N, N);
	vd = new_flt_mtx(N, N);
	vr = new_flt_mtx(N, N);
	clslbl = new_ch_mtx(N, HLEN);
	if (!nx || !nc || !nv || !ncv || !tree || !nx || !c || !crit || !cv || !s ||
	     !mval || !vc || !vd || !vr || !clslbl) goto L90;

	char fname[HLEN] = "soybean.both";
/*======================================================================
Large Soybean Database
From the University of California at Irvine Machine Learning repository

R.S. Michalski and R.L. Chilausky
"Learning by Being Told and Learning from Examples:
An Experimental Comparison of the Two Methods of Knowledge Acquisition
in the Context of Developing an Expert System for Soybean Disease Diagnosis",
International Journal of Policy Analysis and Information Systems,
Vol. 4, No. 2, 1980.

Classes are:

diaporthe-stem-canker, charcoal-rot, rhizoctonia-root-rot, phytophthora-rot,
brown-stem-rot, powdery-mildew, downy-mildew, brown-spot, bacterial-blight,
bacterial-pustule, purple-seed-stain, anthracnose, phyllosticta-leaf-spot,
alternarialeaf-spot, frog-eye-leaf-spot, diaporthe-pod-&-stem-blight,
cyst-nematode, 2-4-d-injury, herbicide-injury

Attributes are:

___Attribute_________Type__________#___Values_______________________
   date              Continuous        april,may,june,july,august,
                                         september,october,?
   plant-stand       Categorical   2   normal,lt-normal,?
   precip            Ordinal           lt-norm,norm,gt-norm,?
   temp              Ordinal           lt-norm,norm,gt-norm,?
   hail              Categorical   2   yes,no,?
                     (Dichotomous)
   crop-hist         Ordinal           diff-lst-year,same-lst-yr,
                                         same-lst-two-yrs,same-lst-sev-yrs,?
   area-damaged      Categorical   4   scattered,low-areas,upper-areas,
                                         whole-field,?
   severity          Ordinal           minor,pot-severe,severe,?
   seed-tmt          Categorical   3   none,fungicide,other,?
   germination       Ordinal           90-100%,80-89%,lt-80%,?
   plant-growth      Categorical   2   norm,abnorm,?
   leaves            Categorical   2   norm,abnorm.
   leafspots-halo    Categorical   3   absent,yellow-halos,no-yellow-halos,?
   leafspots-marg    Categorical   3   w-s-marg,no-w-s-marg,dna,?
   leafspot-size     Categorical   3   lt-1/8,gt-1/8,dna,?
   leaf-shread       Categorical   2   absent,present,?
                     (Dichotomous)
   leaf-malf         Categorical   2   absent,present,?
                     (Dichotomous)
   leaf-mild         Categorical   3   absent,upper-surf,lower-surf,?
   stem              Categorical   2   norm,abnorm,?
   lodging           Categorical   2   yes,no,?
                     (Dichotomous)
   stem-cankers      Categorical   4   absent,below-soil,above-soil,
                                         above-sec-nde,?
   canker-lesion     Categorical   4   dna,brown,dk-brown-blk,tan,?
   fruiting-bodies   Categorical   2   absent,present,?
                     (Dichotomous)
   external decay    Categorical   3   absent,firm-and-dry,watery,?
   mycelium          Categorical   2   absent,present,?
                     (Dichotomous)
   int-discolor      Categorical   3   none,brown,black,?
   sclerotia         Categorical   2   absent,present,?
                     (Dichotomous)
   fruit-pods        Categorical   4   norm,diseased,few-present,dna,?
   fruit spots       Categorical   5   absent,colored,brown-w/blk-specks,
                                         distort,dna,?
   seed              Categorical   2   norm,abnorm,?
                     (Dichotomous)
   mold-growth       Categorical   2   absent,present,?
                     (Dichotomous)
   seed-discolor     Categorical   2   absent,present,?
                     (Dichotomous)
   seed-size         Categorical   2   norm,lt-norm,?
                     (Dichotomous)
   shriveling        Categorical   2   absent,present,?
                     (Dichotomous)
   roots             Categorical   3   norm,rotted,galls-cysts,?
======================================================================*/

/* info above is summarized thus: */
	static int catc[O] = { 0,2,0,0,2,0,4,0,3,0,2,2,3,3,3,2,2,3,2,2,4,4,2,
	  3,2,3,2,4,5,2,2,2,2,2,3 };
	static int catd[O] = { 0,2,0,0,1,0,4,0,3,0,2,2,3,3,3,1,1,3,2,1,4,4,1,
	  3,1,3,1,4,5,1,1,1,1,1,3 };
	static int catr[P] = {
	  0,0,0,0,0,0,0,0,0,0,
	  0,0,0,0,0,0,0,0,0,0,
	  0,0,0,0,0,0,0,0,0,0,
	  0,0,0,0,0,0,0,0,0,0,
	  0,0,0,0,0,0,0,0,0,0,
	  0,0,0,0 };
	static bool ordreq[O] = { false,false,true,true,false,true,
	  false,true,false,true,false,false,false,false,false,
	  false,false,false,false,false,false,false,false,false,
	  false,false,false,false,false,false,false,false,false,
	  false,false };
	static bool isdix[O] = { false,false,false,false,true,false,
	  false,false,false,false,false,false,false,false,false,true,
	  true,false,false,true,false,false,true,false,true,false,true,
	  false,false,true,true,true,true,true,false };
/* weights (not to be taken seriously) */
	static float w[P] = { 2.742f,1.944f,6.565f,2.532f,2.2f,2.227f,1.386f,
	  2.997f,1.378f,2.186f,3.578f,4.254f,3.271f,4.313f,4.405f,2.795f,
	  5.595f,12.57f,3.092f,3.847f,2.163f,2.124f,3.362f,2.623f,12.17f,
	  12.57f,12.57f,3.162f,2.544f,4.05f,4.626f,4.764f,4.301f,5.215f,
	  3.464f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,
	  0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f };
/* begin */
	srand( (unsigned int) time(NULL));
	for (i = 0; i < HLEN; ++i) empty[i] = '\0';
/* open the file */
	iou = fopen(fname, "r");
	if (NULL == iou) {
		printf("trouble opening file: %s \n", fname);
		return -1;
	}
/* read in data */
	for (i = 0; i < N; ++i) {
		sptr = fgets(text, TXTLEN, iou);
		if (NULL == sptr) printf("trouble reading from file, line %i \n", i+1);
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
/* scan each token */
			if (1 == l) { /* answer key */
				strncpy(clslbl[i], raw, HLEN);
			} else if (NULL == strstr(raw, "?")) { /* non-missing */
				ferr = sscanf(raw, " %f", &x[i][l-2]);
				if (0 == ferr) printf("trouble reading X, i: %i, l: %i \n", i, l);
				x[i][l-2] = x[i][l-2] + 1.f;
				nx[i][l-2] = 1;
			} else {
				x[i][l-2] = xbogus;
				nx[i][l-2] = 0;
			}
			jbeg = j + 1;
			if (l > O+1) printf("!trouble reading, l>o! \n");
		}
	}
	ferr = fclose(iou);
	if (0 != ferr) printf("trouble closing file %s \n", fname);

/*======================================================================
  preprocess input
======================================================================*/
/* scale date to unit range */
	for (i = 0; i < N; ++i) x[i][0] /= 6.f;
/* rescale ordinal variables according to the method of DAISY */
	for (l = 0; l < O; ++l) {
		if (ordreq[l]) {
			xlea = 1.e33f;
			xgre = 0.f;
			for (i = 0; i < N; ++i) {
				if (nx[i][l] > 0) {
					xlea = fminf(xlea, x[i][l]);
					xgre = fmaxf(xgre, x[i][l]);
				}
			}
			for (i = 0; i < N; ++i) {
				if (nx[i][l] > 0) {
					x[i][l] = (x[i][l] - xlea) / (xgre - xlea);
				}
			}
		}
	}
/* embed */
	ifault = mds (x, nx, N, O, w, catc, robust, vc, &mc, e, t);
	if (0 != ifault) printf("mds on categorical error# %i \n", ifault);
	if (mc < M) printf("mds on categorical %i dimensions!\n", mc);
/* embed using dichotomous variables */
	for (i = 0; i < N; ++i) {
		for (l = 0; l < O; ++l) {
			if (isdix[l]) x[i][l] -= 1.f;
		}
	}
	ifault = mds (x, nx, N, O, w, catd, robust, vd, &md, e, t);
	if (0 != ifault) printf("mds with dichotomous error# %i \n", ifault);
	if (md < M) printf("mds with dichotomous %i dimensions!\n", md);
	for (i = 0; i < N; ++i) {
		for (l = 0; l < O; ++l) {
			if (isdix[l]) x[i][l] += 1.f;
		}
	}
/* recode categoricals with corresponding column of simplex */
	ifault = rscode (x, nx, N, P, w, catc, O, s, CATMAX, robust);
	if (0 != ifault) printf("rscode: error# %i \n", ifault);
/* relabel class ID's */
	for (i = 0; i < N; ++i) nx0[i] = 1;
	ifault = ireid (clslbl, y, nx0, iperm, jperm, N, HLEN, raw);
	if (0 != ifault) printf("ireid: error #%i \n", ifault);
/* call the subroutine */
	k = KDEF;
	ifault = cluela (x, nx, N, P, w, c, nc, KMIN, &k, KMAX, z, pop, f,
	                 robust, sosdo, iwork, rwork, tree, crit, &u);
/* find accuracy */
	printf("cluela returns #%i, made %i clusters\n", ifault, k);
	if (0 == ifault) {
		accari = arind (y, z, N, G, k, mval, aval, bval); /* adjusted Rand index */
		accvi = vinfo (y, z, N, G, k, mval, aval, bval); /* variation of information */
	}
	printf("Returned residual: %12.4g\n", u);
	printf("Adjusted Rand Index: %5.3f\n", accari);
	printf("Variation of Information: %7.3f\n", accvi);
/* Write CSV file with coordinates in different columns by cluster */
	strncpy (fname, "soy_clusters.csv", HLEN);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("soy: trouble opening %s\n", fname);
	ifault = wcsvsp (vc, N, z, iou);
	if (0 != ifault) printf("wcsvsp: error #%i\n", ifault);
	ferr = fclose(iou);
	if (0 != ferr) printf("trouble closing file %s\n", fname);
/* Write CSV file with coordinates in different columns by class */
	strncpy(fname, "soy_classes.csv", HLEN);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("soy: trouble opening %s\n", fname);
	ifault = wcsvsp (vc, N, y, iou);
	if (0 != ifault) printf("wcsvsp: error #%i\n", ifault);
	ferr = fclose(iou);
	if (0 != ferr) printf("trouble closing file %s\n", fname);
/* write F values */
	strncpy (fname, "soy_fvals.csv", HLEN);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("soy: trouble opening %s\n", fname);
	for (l = KMIN; l <= KMAX; ++l) {
		ifault = fprintf(iou, "%i, %9.3g\n", l, f[l-KMIN]);
		if (0 == ifault) printf("soy: trouble writing f stats\n");
	}
	ferr = fclose(iou);
	if (0 != ferr) printf("soy: trouble closing file %s\n", fname);
/* repeat for dichotomous variables */
	k = KDEF;
	for (i = 0; i < N; ++i) {
		wv[i] = 1.f;
		for (j = 0; j < N; ++j) {
			nv[i][j] = 1;
		}
	}
	ifault = cluela (vd, nv, N, md, wv, cv, ncv, KMIN, &k, KMAX, z, pop, f,
	                 false, sosdo, iwork, rwork, tree, crit, &u);
	printf("cluela returns #%i, made %i clusters\n", ifault, k);
	if (0 == ifault) {
		accari = arind (y, z, N, G, k, mval, aval, bval); /* adjusted Rand index */
		accvi = vinfo (y, z, N, G, k, mval, aval, bval); /* variation of information */
	}
	printf("Returned residual: %12.4g\n", u);
	printf("Adjusted Rand Index: %5.3f\n", accari);
	printf("Variation of Information: %7.3f\n", accvi);
	strncpy (fname, "soy_dichot_clusters.csv", HLEN);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("soy: trouble opening %s\n", fname);
	ifault = wcsvsp (vd, N, z, iou);
	if (0 != ifault) printf("wcsvsp: error #%i\n", ifault);
	ferr = fclose(iou);
	if (0 != ferr) printf("trouble closing file %s\n", fname);
	strncpy(fname, "soy_dichot_classes.csv", HLEN);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("soy: trouble opening %s\n", fname);
	ifault = wcsvsp (vd, N, y, iou);
	if (0 != ifault) printf("wcsvsp: error #%i\n", ifault);
	ferr = fclose(iou);
	if (0 != ferr) printf("trouble closing file %s\n", fname);
	strncpy (fname, "soy_dichot_fvals.csv", HLEN);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("soy: trouble opening %s\n", fname);
	for (l = KMIN; l <= KMAX; ++l) {
		ifault = fprintf(iou, "%i, %9.3g\n", l, f[l-KMIN]);
		if (0 == ifault) printf("soy: trouble writing f stats\n");
	}
	ferr = fclose(iou);
	if (0 != ferr) printf("soy: trouble closing file %s\n", fname);
/* embed the recoded data */
	ifault = mds (x, nx, N, P, w, catr, robust, vr, &mr, e, t);
	if (0 != ifault) printf("mds on recoded data error# %i \n", ifault);
	if (mr < M) printf("mds on recoded data %i dimensions!\n", mr);
/* compare embeddings */
	worst = 0.f;
	for (i = 0; i < N; ++i) {
		for (j = 0; j < N; ++j) {
			dif = fminf(fabsf(vc[i][j] + vr[i][j]),
			            fabsf(vc[i][j] - vr[i][j]));
			worst = fmaxf(worst, dif);
		}
	}
	printf("Greatest difference between original and recoded embeddings: %9.2e \n", worst);
L90: /* Release memory */
	if (nx) del_int_mtx(nx, N);
	if (nc) del_int_mtx(nc, KMAX);
	if (nv) del_int_mtx(nv, N);
	if (ncv) del_int_mtx(ncv, KMAX);
	if (tree) del_int_mtx(tree, 2*KMAX);
	if (x) del_flt_mtx(x, N);
	if (c) del_flt_mtx(c, KMAX);
	if (crit) del_flt_mtx(crit, 2*KMAX);
	if (cv) del_flt_mtx(cv, KMAX);
	if (s) del_flt_mtx(s, CATMAX);
	if (mval) del_flt_mtx(mval, G);
	if (vc) del_flt_mtx(vc, N);
	if (vd) del_flt_mtx(vd, N);
	if (vr) del_flt_mtx(vr, N);
	if (clslbl) del_ch_mtx(clslbl, N);
	printf("program complete\n");
	return 0;
} /* end of main */
