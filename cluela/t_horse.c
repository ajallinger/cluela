/*======================================================================
test clustering on the horse colic data
by Andy Allinger, 2013-2020, public domain
======================================================================*/
#include <stdbool.h>  /* standard C99 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

int kmns (float **a, int m, int n, float **c, int k,  /* defined in kmns.c */
          int *ic1, int *ic2, int *nc, float *an1, float *an2,
	     int *ncp, float *d, int *itran, int *live, int iter, float *wss);
bool safdiv (float numer, float denom);  /* defined in safdiv.c */
int **new_int_mtx (int m, int n);  /* defined in array.c */
float **new_flt_mtx (int m, int n);
int del_int_mtx (int **b, int m);
int del_flt_mtx (float **x, int m);
void avedev (float **x, int **nx, int n, int p,  /* defined in stats.c */
             int *na, float *ave, float *dev);
void medmad (float **x, int **nx, int n, int p,
             float *rwork, float *ave, float *dev);
void jscode (int *jnam, int *nx, int *ind, int *jnd, int n);  /* defined in prep.c */
void ordnal (float *x, int *nx, int *ind, int *jnd, int n);
int rscode (float **x, int **nx, int n, int p, float *w,
            int *cat, int o, float **s, int catmax, bool robust);
int mds (float **x, int **nx, int n, int p, float *w, int *cat, bool robust,  /* defined in mds.c */
         float **v, int *m, float *e, float *t);
bool almeq (float r, float s);  /* defined in cluela.c */
void shuffl (int *i, int n);
void dinit (float **x, int **nx, int n, int p, float *w, float **c, int **nc,
            int k, bool robust, float *shortx, float *cumd, float *thresh);
int kclust (float **x, int **nx, int n, int p, float *w, float **c, int **nc,
            int kmin, int *k, int kmax, int *z, int *pop, bool miss,
            bool robust, bool sosdo, float thresh, int *livec, int *livex,
            int *iperm, int **tree, float **crit, float *shortc, float *shortx,
            float *work, float *u);
int cluela (float **x, int **nx, int n, int p, float *w, float **c, int **nc,
            int kmin, int *k, int kmax, int *z, int *pop, float *f,
            bool robust, bool sosdo, int *iwork, float *rwork, int **tree,
            float **crit, float *u);
int cluste (float **x, int **nx, int n, int p, float *w, float **c, int **nc,  /* defined in cluste.c */
            int kmin, int *k, int kmax, int *z, int *pop, bool robust,
            bool sosdo, int *iwork, float *rwork, int **tree, float **crit, float *u);
void osiman (int n, float *x, float e, float sbeg, float send,  /* defined in osiman.c */
             float tbeg, float tend, int nit, float *xkeep, int *iflag);
void neugas (float **x, int n, int p, float **c, int k, int *z, int *ind,
             float *d, float *u);  /* defined in neugas.c */
int wcsvsp (float **r, int n, int *z, FILE *iou);  /* defined in dump.c */

/* fixed parameters */
#define CATMAX 346  /* most categories any attribute */
#define HLEN 32     /* length of filename strings */
#define G 62        /* # distinct diagnoses */
#define PIN 25      /* # input attributes */
#define N 368       /* # objects */
#define NS 12       /* node size of binary tree */
#define O 86        /* # attributes after recoding multi-attribute */
#define P 440       /* # attributes after recoding categoricals */
#define Q 3         /* max concurrent diagnoses */
#define QN 1104     /* Q*N */
#define TXTLEN 256  /* length of a line of input file */

/* adjustable parameters */
#define KMIN 1      /* fewest clusters to make */
#define KDEF 2      /* default # clusters to make */
#define KMAX 20     /* most clusters to make */

/* Main program */
int main (void)
{
/* Constants */
	const int ITER = 1000,
	          NIT = 10000;
	const float BIG = 1e27f,
	            FONE = 1.f,
	            FTWO = 2.f,
	            FTHREE = 3.f,
	            PEX = 3.f,
	            SBEG = 0.1f,
	            SEND = 0.001f,
	            TBEG = 1.e-3f,
	            TEND = 1.e-6f,
	            XBOGUS = -9.f;
/* Local variables */
	int h, i, j, k, l, m, cplen, ferr, gnosis[QN], hbeg, hosptl[N], iperm[QN], ifault, iflag,
	    itran[KDEF], iwork[4*N+KMAX], jperm[QN], kbest, lis, live[KDEF], ll,
	    na[PIN], ncp[KDEF], nx0[QN], pop[KMAX], z[N], zbest[N];
	float an1[KDEF], an2[KDEF], ave[PIN], dev[PIN], e[N], f[KMAX-KMIN+1], pnlty,
	      pnlmin, r, result, rwork[2*N+KMAX], smidge, t[N], thresh, u, ubest, w[P], wall,
	      wbak[O], wbest[PIN], win[PIN], winit, wkeep[PIN], wss[KDEF], x0[N];
	bool miss, robust, sosdo;
	char empty[HLEN], raw[HLEN], text[TXTLEN], *sptr;
	FILE *iou;
	clock_t after, before;
/* allocate matrices */
	int **nc, **nx, **nxbak, **tree;
	float **crit, **x, **c, **xbak, **s, **v;
	nc = new_int_mtx (KMAX, P);
	nx = new_int_mtx (N, P);
	nxbak = new_int_mtx (N, O);
	tree = new_int_mtx (2*KMAX, NS);
	crit = new_flt_mtx (2*KMAX, 4);
	x = new_flt_mtx (N, P);
	c = new_flt_mtx (KMAX, P);
	xbak = new_flt_mtx (N, O);
	s = new_flt_mtx (CATMAX, CATMAX-1);
	v = new_flt_mtx (N, N);
	if (!nc || !nx || !nxbak || !tree || !crit ||
	     !x || !c || !xbak || !s || !v) {
		printf("trouble allocating arrays\n");
		goto L90;
	}

	char fname[HLEN] = "horse-colic.both";
/*=======================================================================
Horse colic data
From the University of California at Irvine Machine Learning repository

Creators: Mary McLeish & Matt Cecile
          Department of Computer Science
          University of Guelph
          Guelph, Ontario, Canada N1G 2W1
          mdmcleish@water.waterloo.edu

Classes: As the analyst may choose.  Possibilities are outcome
         (lived / died / was euthanized) or type of lesion (one of 62
          diagnosis numbers).

Attributes are:

________Attribute__________Type____________#___Values___________________
    0   surgery            Boolean         2   1:  yes
                                               2:  no

    1   age                Categorical     2   1:  adult
                                               9:  young     ! 9?

    2   hospital number    Categorical   346   ID number assigned to horse
If the horse lives, then this might be considered a categorical
variable.  Since a horse only dies once, for horses in the groups
Died and Was Euthanized, this is a contra-categorical variable, since
the presence of a horse in one of these groups means that no other
instances of this horse may be present.

    3   rectal temp        Continuous          [degrees Celsius]

    4   pulse              Continuous          [beats per minute]

    5   respiratory rate   Continuous

    6   extremities temp   Ordinal             1:  Normal
                                               2:  Warm
                                               3:  Cool
                                               4:  Cold
                                                               Wrong order!

    7   peripheral pulse   Ordinal             1:  normal
                                               2:  increased
                                               3:  reduced
                                               4:  absent
                                                               Wrong order!

    8   mucous membranes   Categorical     6   1:  normal pink
                                               2:  bright pink
                                               3:  pale pink
                                               4:  pale cyanotic
                                               5:  bright red / injected
                                               6:  dark cyanotic

    9   capillary refill   Categorical     3   1:  < 3 seconds
                                               2:  >= 3 seconds
                                               3:  undocumented!

   10   pain               Continuous          1:  alert, no pain
                                               2:  depressed
                                               3:  intermittent mild pain
                                               4:  intermittent severe pain
                                               5:  continuous severe pain
The file 'horse-colic.names' says: "should NOT be treated as a ordered or
discrete variable!"  Obey this instruction and treat pain as continuous.

   11   peristalsis        Ordinal             1:  hypermotile
                                               2:  normal
                                               3:  hypomotile
                                               4:  absent

   12   abdom distension   Ordinal             1:  none
                                               2:  slight
                                               3:  moderate
                                               4:  severe

   13   nasogastric tube   Ordinal             1:  none
                                               2:  slight
                                               3:  significant

   14   naso. reflux       Ordinal             1:  none
                                               2:  > 1 liter
                                               3:  < 1 liter
                                                               Wrong order!

   15   reflux PH          Continuous          0 ... 14

   16   rectal feces       Ordinal             1:  normal
                                               2:  increased
                                               3:  decreased
                                               4:  absent
                                                               Wrong order!

   17   abdomen            Categorical     5   1:  normal
                                               2:  other
                                               3:  firm feces large intestine
                                               4:  distended small intestine
                                               5:  distended large intestine

   18   packed cell vol.   Continuous          red cells by volume in blood

   19   total protein      Continous           [gms/dL]

   20   abdominocentesis   Categorical     3   1:  clear
                                               2:  cloudy
                                               3:  serosanguinous

   21   abdom. protein     Continuous          [gms/dL]

   22   outcome            Categorical     3   1:  lived
                                               2:  died
                                               3:  was euthanized

   23   surgical lesion?   Categorical     2   1:  Yes
                                               2:  No

   24   type of lesion     Multi          62
   25
   26
Three columns contain either 00000 or one of 62 distinct numeric codes
for a diagnosis.

        cp_data            Boolean             1:  Yes
                                               2:  No
This column is irrelevant and shall be omitted.

=====================================================================*/
/* info above is summarized thus: */
	int cat[O] = { 2,2,346,0,0,0,0,0,6,3,0,0,0,0,0,0,0,5,0,0,3,0,3,2,
	   2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
	   2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2 };
	const bool ordreq[PIN] = { false, false, false, false, false,
	                           false,  true,  true, false, false,
	                           false,  true,  true,  true,  true,
	                           false,  true, false, false, false,
	                           false, false, false, false, false };
/* set initial values */
	robust = false;
	winit = 1.f / (float)PIN;
	for (j = 0; j < PIN; ++j) win[j] = winit;
	for (i = 0; i < HLEN; ++i) empty[i] = '\0';
	for (i = 0; i < QN; ++i) nx0[i] = 1;
	smidge = 1.f / (float)(PIN * PIN);
	srand((unsigned int)time(NULL)); /* initialize random number generator */
/* read data file */
	iou = fopen(fname, "r");
	if (NULL == iou) {
		printf("trouble opening input file\n");
		return -1;
	}
	for (i = 0; i < N; ++i) {
		sptr = fgets(text, TXTLEN, iou);
		if (NULL == sptr) printf("trouble reading from file line %d \n", i+1);
		ll = (int)strlen(text);
		hbeg = 0;
		l = 0;
		j = 0;
		for (h = 0; h < ll; ++h) {
			if (' ' == text[h]) {  /* split on spaces */
				++l;
				strncpy(raw, empty, HLEN);
				cplen = h - hbeg;
				if (cplen > HLEN) cplen = HLEN;
				strncpy(raw, &text[hbeg], cplen);
			} else if (h == ll-1) {
				++l;
				strncpy(raw, empty, HLEN);
				cplen = h - hbeg + 1;
				if (cplen > HLEN) cplen = HLEN;
				strncpy(raw, &text[hbeg], cplen);
			} else {
				continue;
			}
			if (3 == l) {  /* ID number */
				nx[i][j] = 1;  /* nothing missing */
				++j;
				ferr = sscanf(raw, " %d", &hosptl[i]);
				if (0 == ferr) printf("error reading ID line %d\n", i);
			} else if (25 == l) {
				ferr = sscanf(raw, " %d", &gnosis[Q*i]);
				if (0 == ferr) printf("error reading attrib1 object %d\n", i);
			} else if (26 == l) {
				ferr = sscanf(raw, " %d", &gnosis[Q*i+1]);
				if (0 == ferr) printf("error reading attrib2 object %d\n", i);
			} else if (27 == l) {
				ferr = sscanf(raw, " %d", &gnosis[Q*i+2]);
				if (0 == ferr) printf("error reading attrib3 object %d\n", i);
			} else if (28 == l) {
				; /* ignore cp_data */
			} else if (NULL == strstr(raw, "?")) {
 				ferr = sscanf(raw, " %f", &x[i][j]);
 				if (0 == ferr) printf("error reading X, i: %d l: %d\n", i, l);
 				nx[i][j] = 1;
				++j;
			} else {  /* missing */
				x[i][j] = XBOGUS;
				nx[i][j] = 0;
				++j;
			}
			hbeg = h + 1;
		}
	}
	ferr = fclose(iou);
	if (0 != ferr) printf("trouble closing file %s\n", fname);
/* preprocess data */
	for (i = 0; i < N; ++i) if (x[i][1] > 1.5f) x[i][1] = 2.f;  /* relabel the young / adult variable */
	jscode (hosptl, nx0, iperm, jperm, N);  /* relabel case ID's */
	for (i = 0; i < N; ++i) x[i][2] = (float)hosptl[i];
/* fix wrong ordered ordinals */
	for (i = 0; i < N; ++i) {
		if (almeq(x[i][6], FONE)) {  /* put subjective temperature in order */
			x[i][6] = FTWO;
		} else if (almeq(x[i][6], FTWO)) {
			x[i][6] = FONE;
		}
		if (almeq(x[i][7], FONE)) {  /* put subjective pulse in order */
			x[i][7] = FTWO;
		} else if (almeq(x[i][7], FTWO)) {
			x[i][7] = FONE;
		}
		if (almeq(x[i][14], FTWO)) {  /* put reflux in order */
			x[i][14] = FTHREE;
		} else if (almeq(x[i][14], FTHREE)) {
			x[i][14] = FTWO;
		}
		if (almeq(x[i][16], FONE)) {  /* put feces in order */
			x[i][16] = FTWO;
		} else if (almeq(x[i][16], FTWO)) {
			x[i][16] = FONE;
		}
	}
	if (robust) {  /* averages */
		medmad (x, nx, N, PIN, rwork, ave, dev);
	} else {
		avedev (x, nx, N, PIN, na, ave, dev);
	}
	for (j = 0; j < PIN; ++j) {
		if (ordreq[j]) {
			for (i = 0; i < N; ++i) {
				x0[i] = x[i][j];
				nx0[i] = nx[i][j];
			}
			ordnal (x0, nx0, iperm, jperm, N);
			for (i = 0; i < N; ++i) x[i][j] = x0[i];
		} else if (0 == cat[j]) {
			if (safdiv(FONE, dev[j])) {
				r = FONE / dev[j];
			} else {
				printf("variable %d: no variance\n", j);
				r = 0.f;
			}
			for (i = 0; i < N; ++i) if (nx[i][j] > 0) x[i][j] = (x[i][j] - ave[j]) * r;
		}
	}
/* convert diagnosis codes to small consecutive integers */
	for (i = 0; i < QN; ++i) nx0[i] = 1;
	jscode (gnosis, nx0, iperm, jperm, QN);
	for (i = 0; i < QN; ++i) --gnosis[i];  /* zero is null attribute */
/* Recode the multi-attribute */
	for (i = 0; i < N; ++i) {
		for (j = 24; j < O; ++j) {
			x[i][j] = 2.f;  /* no */
			nx[i][j] = 1;
		}
	}
	for (i = 0; i < N; ++i) {
		for (h = 0; h < Q; ++h) {
			if (gnosis[Q*i+h] > 0) {
				j = 23 + gnosis[Q*i+h];
				x[i][j] = 1.f;  /* yes */
			}
		}
	}
/* Back up input. */
	for (i = 0; i < N; ++i) {
		for (j = 0; j < O; ++j) {
			xbak[i][j] = x[i][j];
			nxbak[i][j] = nx[i][j];
		}
	}
	for (j = 0; j < PIN; ++j) w[j] = win[j];
	for (j = PIN-1; j < O; ++j) w[j] = winit / (float)G;  /* spread weight across multi-attribute >*/
	for (j = 0; j < O; ++j) wbak[j] = w[j];
/* recode categoricals with corresponding column of simplex */
	ifault = rscode (x, nx, N, P, w, cat, O, s, CATMAX, robust);
	if (0 != ifault) printf("rscode: error #%d\n", ifault);
/* call CLUELA on the full problem - equal weights */
	printf("Original problem:\n");
	k = KDEF;
	sosdo = !robust;
	before = clock();
	ifault = cluela (x, nx, N, P, w, c, nc, KMIN, &k, KMAX, z, pop, f,
	                  robust, sosdo, iwork, rwork, tree, crit, &u);
	after = clock();
	result = (float)(after - before) / (float)CLOCKS_PER_SEC;
	printf("cluela returns #%d in %f seconds.  Made %d clusters with residual %f\n",
	       ifault, (double)result, k, (double)u);
/* write F values */
	strncpy(fname, "horse_fvals_equal.csv", HLEN);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("horse: trouble opening %s \n", fname);
	for (l = KMIN; l <= KMAX; ++l) {
		ifault = fprintf(iou, "%5i, %9.3g\n", l, (double)f[l-KMIN]);
		if (ifault <= 0) printf("horse: trouble writing F stats\n");
	}
	fclose(iou);
/* project to lower dimension */
	ifault = mds (xbak, nxbak, N, O, wbak, cat, robust, v, &m, e, t);
	printf("mds returns #%d with %d directions\n", ifault, m);
/* Scatter plot of CLUELA result */
	strncpy(fname, "horse_clusters_equal.csv", HLEN);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("horse: trouble opening %s\n", fname);
	ifault = wcsvsp (v, N, z, iou);
	if (0 != ifault) printf("wcsvsp: error #%d\n", ifault);
	fclose(iou);
/* Call KCLUST directly.  Optimize weights */
	miss = true;
	iflag = 1;
	ubest = 0.f;
	pnlmin = BIG;
	pnlty = 0.f;
	lis = 0;
	for (l = 0; l <= NIT; ++l) {
		k = KDEF;
		osiman (PIN, win, pnlty, SBEG, SEND, TBEG, TEND, NIT, wkeep, &iflag);
		for (i = 0; i < N; ++i) {
			for (j = 0; j < O; ++j) {
				x[i][j] = xbak[i][j];
				nx[i][j] = nxbak[i][j];
			}
		}
		wall = 0.f;
		for (j = 0; j < PIN; ++j) {
			if (win[j] < smidge) win[j] = smidge;
			wall += win[j];
		}
		for (j = 0; j < PIN; ++j) {
			win[j] = win[j] / wall;
			w[j] = win[j];
		}
		for (j = PIN-1; j < O; ++j) w[j] = win[PIN-1] / (float)G;
		pnlty = 0.f;
		for (j = 0; j < PIN; ++j) pnlty += powf(win[j], PEX);  /* weight penalty */
		ifault = rscode(x, nx, N, P, w, cat, O, s, CATMAX, robust);
		if (0 != ifault) printf("rscode: error #%d\n", ifault);
		dinit (x, nx, N, P, w, c, nc, k, robust, &rwork[0], &rwork[N], &thresh);
		ifault = kclust (x, nx, N, P, w, c, nc, KMIN, &k, KMAX, z, pop,
		  miss, robust, sosdo, BIG, &iwork[0], &iwork[KMAX], &iwork[KMAX+N],
		  tree, crit, &rwork[0], &rwork[KMAX], &rwork[KMAX+N], &u);
		pnlty *= u;
		printf("iter %d kclust returns #%d.  Made %d clusters with residual %11.3e penalty %11.3e",
		  l, ifault, k, (double)u, (double)pnlty);
		if (0 == ifault && KDEF == k) {
			if (pnlty < pnlmin) {
				pnlmin = pnlty;
				ubest = u;
				for (i = 0; i < N; ++i) zbest[i] = z[i];
				for (j = 0; j < PIN; ++j) wbest[j] = win[j];
				printf(" <---- New low penalty\n");
				lis = l;
			} else {
				printf("\n");
			}
		} else {
			pnlty = BIG;
			printf("\n");
		}
	}
	printf("penalty exponent %f and robust=%s yield optimized weights: ",
	       (double)PEX, robust?"T":"F");
	for (j = 0; j < PIN; ++j) printf("%7.3f", wbest[j]);
	printf("\nbest value iter %d ubest: %f\n", lis, (double)ubest);
/* project to lower dimension */
	for (j = 0; j < PIN; ++j) w[j] = wbest[j];
	for (j = PIN-1; j < O; ++j) w[j] = wbest[PIN-1] / (float)G;
	ifault = mds (xbak, nxbak, N, O, w, cat, robust, v, &m, e, t);
	printf ("mds returns #%d with %d directions\n", ifault, m);
/* Scatter plot of KCLUST result */
	strncpy(fname, "horse_clusters_kclust.csv", HLEN);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("horse: trouble opening %s\n", fname);
	ifault = wcsvsp (v, N, zbest, iou);
	if (0 != ifault) printf("wcsvsp: error #%d\n", ifault);
	fclose(iou);
/* There is no longer any missing data.  Remember that dimension is now M */
	for (i = 0; i < N; ++i) for (j = 0; j < m; ++j) nx[i][j] = 1;
	for (j = 0; j < m; ++j) w[j] = 1.f;  /* remember to change W */
	robust = false;
/* Call CLUELA on the projected data. */
	k = KDEF;
	sosdo = true;
	cluela (v, nx, N, m, w, c, nc, KMIN, &k, KMAX, z,	pop, f,
	         robust, sosdo, iwork, rwork, tree, crit, &u);
	printf("cluela returns #%d  Residual: %f\n", ifault, (double)u);
	strncpy(fname, "horse_clusters_weighted.csv", HLEN);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("horse: trouble opening %s\n", fname);
	ifault = wcsvsp (v, N, z, iou);
	if (0 != ifault) printf("wcsvsp: error #%d\n", ifault);
	fclose(iou);
	strncpy(fname, "horse_fvals_weighted.csv", HLEN);  /* write F values */
	iou = fopen(fname, "w");
	if (NULL == iou) printf("horse: trouble opening %s\n", fname);
	for (l = KMIN; l <= KMAX; ++l) {
		ifault = fprintf(iou, "%5i, %9.3g\n", l, (double)f[l-KMIN]);
		if (ifault <= 0) printf("horse: trouble writing F stats\n");
	}
	fclose(iou);
/* Call alternate interface */
	before = clock();
	sosdo = false;
	ubest = BIG;
	kbest = 0;
	for (h = 0; h < 10; ++h) {
		k = KDEF;
		ifault = cluste (v, nx, N, m, w, c, nc, KMIN, &k, KMAX, z, pop,
	                      robust, sosdo, iwork, rwork, tree, crit, &u);
		if (0 == ifault) {
			if (KDEF == k) {
				kbest = k;
				if (u < ubest) {
					ubest = u;
					for (i = 0; i < N; ++i) zbest[i] = z[i];
				}
			}
		} else {
			printf("cluste: error #%d\n", ifault);
		}
	}
	after = clock();
	if (KDEF != kbest) printf("cluste did not find requested number of clusters\n");
	result = (float)(after - before) / (float)CLOCKS_PER_SEC;
	printf("cluste returns in %f seconds.  Residual: %f\n", (double)result, (double)ubest);
	strncpy(fname, "horse_clusters_cluste.csv", HLEN);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("horse: trouble opening %s\n", fname);
	ifault = wcsvsp (v, N, zbest, iou);
	if (0 != ifault) printf("wcsvsp: error #%d\n", ifault);
	fclose(iou);
/* Call KMNS on projected data.  Remember to initalize centers. */
	before = clock();
	k = KDEF;
	ubest = BIG;
	for (h = 0; h < 10; ++h) {
		dinit (v, nx, N, m, w, c, nc, k, robust, &rwork[0], &rwork[N], &thresh);
		ifault = kmns (v, N, m, c, k, z, iperm, pop,
		                an1, an2, ncp, rwork, itran, live, ITER, wss);
		u = 0.f;
		for (l = 0; l < k; ++l) u += wss[l];  /* combine per-cluster residuals */
		if (0 == ifault) {
			if (u < ubest) {
				ubest = u;
				for (i = 0; i < N; ++i) zbest[i] = z[i];
			}
		} else {
			printf("kmns: error #%d\n", ifault);
		}
	}
	after = clock();
	result = (float)(after - before) / (float)CLOCKS_PER_SEC;
	printf("kmns returns in %f seconds. Residual: %f\n", (double)result, (double)ubest);
	strncpy(fname, "horse_clusters_kmns.csv", HLEN);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("horse: trouble opening %s\n", fname);
	ifault = wcsvsp (v, N, zbest, iou);
	if (0 != ifault) printf("wcsvsp: error#%d\n", ifault);
	fclose(iou);
/* Call the neural gas routine. */
	before = clock();
	k = KDEF;
	shuffl (iwork, N);  /* choose random centers */
	for (l = 0; l < k; ++l) {
		i = iwork[l] - 1;
		for (j = 0; j < m; ++j) c[l][j] = x[i][j];
	}
	neugas(v, N, m, c, k, z, iwork, rwork, &u);
	after = clock();
	result = (float)(after - before) / (float)CLOCKS_PER_SEC;
	printf("neugas returns in %f seconds.  Residual: %f\n", (double)result, (double)u);
	strncpy(fname, "horse_clusters_ng.csv", HLEN);
	iou = fopen(fname, "w");
	if (NULL == iou) printf("horse: trouble opening %s\n", fname);
	ifault = wcsvsp (v, N, z, iou);
	if (0 != ifault) printf("wcsvsp: error #%d\n", ifault);
	fprintf(iou, "%9.3g,,,%9.3g\n", (double)c[0][0], (double)c[0][1]);
	fprintf(iou, "%9.3g,,,%9.3g\n", (double)c[1][0], (double)c[1][1]);
	fclose(iou);
	printf("program complete\n");
L90:
	if (nx) del_int_mtx (nx, N);
	if (nc) del_int_mtx (nc, KMAX);
	if (nxbak) del_int_mtx (nxbak, N);
	if (tree) del_int_mtx (tree, 2*KMAX);
	if (crit) del_flt_mtx (crit, 2*KMAX);
	if (x) del_flt_mtx (x, N);
	if (c) del_flt_mtx (c, KMAX);
	if (xbak) del_flt_mtx (xbak, N);
	if (s) del_flt_mtx (s, CATMAX);
	if (v) del_flt_mtx (v, N);
	return 0;
} /* end of horse */
