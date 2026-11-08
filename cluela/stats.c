/*----------------------------------------------------------------------
stats.c - averages and deviations with missing data support
by Andy Allinger, 2017, released to the public domain
This program may be used by any person for any purpose.
----------------------------------------------------------------------*/
#include <math.h>
int median (float *, int, float *);

/***********************************************************************
summary statistics:  mean and standard deviation
Returns 0.0 when not enough data is available

___Name_______Type______In/Out____Description___________________________
   x[n][p]    float**   In        Data matrix
   nx[n][p]   int**     In        1 if data exists for this observation
   n          int       In        Number of objects
   p          int       In        Number of variables
   na[p]      int*      Neither   Count of data available
   ave[p]     float*    Out       Means
   dev[p]     float*    Out       Standard deviations
***********************************************************************/
void avedev (float **x, int **nx, int n, int p,
             int *na, float *ave, float *dev)
{
/* Local variables */
	int i, j;
	float dif;
/* Function Body */
	for (j = 0; j < p; ++j) { /* clear */
		ave[j] = 0.f;
		na[j] = 0;
	}
	for (i = 0; i < n; ++i) { /* add */
		for (j = 0; j < p; ++j) {
			if (nx[i][j] > 0) {
				ave[j] += x[i][j];
				++na[j];
			}
		}
	}
	for (j = 0; j < p; ++j) { /* divide */
		if (na[j] > 0) ave[j] /= (float) na[j];
	}

	for (j = 0; j < p; ++j) dev[j] = 0.f; /* clear */
	for (i = 0; i < n; ++i) { /* add */
		for (j = 0; j < p; ++j) {
			if (nx[i][j] > 0) {
				dif = x[i][j] - ave[j];
				dev[j] += dif * dif;
			}
		}
	}
	for (j = 0; j < p; ++j) { /* divide */
		if (na[j] > 1) {
			dev[j] = sqrtf(dev[j] / (float) (na[j] - 1));
		} else {
			dev[j] = 0.f;
		}
	}
	return;
} /* end of avedev */

/***********************************************************************
summary statistics:  median and median absolute deviation (MAD)
Returns 0.0 when not enough data is available

___Name_______Type______In/Out____Description___________________________
   x[n][p]    float**   In        Data matrix
   nx[n][p]   int**     In        1 if data exists for this observation
   n          int       In        Number of objects
   p          int       In        Number of variables
   rwork[n]   float*    Neither   Workspace
   ave[p]     float*    Out       Medians
   dev[p]     float*    Out       Median absolute deviations
***********************************************************************/
void medmad (float **x, int **nx, int n, int p,
             float *rwork, float *ave, float *dev)
{
/* Local variables */
    int i, j, k;
/* Function Body */
	for (j = 0; j < p; ++j) {
		k = 0;
		for (i = 0; i < n; ++i) { /* fill work array */
			if (nx[i][j] > 0) {
				rwork[k] = x[i][j];
				++k;
			}
		}
		if (k > 0) {
			(void) median (rwork, k, &ave[j]);
		} else {
			ave[j] = 0.f;
		}
	}
	for (j = 0; j < p; ++j) {
		k = 0;
		for (i = 0; i < n; ++i) {
			if (nx[i][j] > 0) { /* absolute differences */
				rwork[k] = fabsf(x[i][j] - ave[j]);
				++k;
			}
		}
		if (k > 0) {
		    (void) median (rwork, k, &dev[j]);
		} else {
		    dev[j] = 0.f;
		}
	}
    return;
} /* end of medmad */

/***********************************************************************
summary statistics:  least and greatest
Returns 0.0 when not enough data is available

___Name_______Type______In/Out____Description__________________________
   x[n][p]    float**   In        Data matrix
   nx[n][p]   int**     In        1 if data exists for this observation
   n          int       In        Number of objects
   p          int       In        Number of variables
   na[p]      int*      Neither   Count of data available
   lea[p]     float*    Out       Least values
   gre[p]     float*    Out       Greatest values
***********************************************************************/
void leagre (float **x, int **nx, int n, int p,
             int *na, float *lea, float *gre)
{
/* Constant */
	static float BIG = 1.e36f;
/* Local variables */
	int i, j;
/* Function Body */
	for (j = 0; j < p; ++j) {
		lea[j] = +BIG;
		gre[j] = -BIG;
		na[j] = 0;
	}
	for (i = 0; i < n; ++i) {
		for (j = 0; j < p; ++j) {
			if (nx[i][j] > 0) {
				if (x[i][j] < lea[j]) lea[j] = x[i][j];
				if (x[i][j] > gre[j]) gre[j] = x[i][j];
				++na[j];
			}
		}
	}
	for (j = 0; j < p; ++j) {
		if (0 == na[j]) {
			lea[j] = 0.f;
			gre[j] = 0.f;
		}
	}
	return;
} /* end of leagre */
