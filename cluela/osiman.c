/*----------------------------------------------------------------------
osiman.c - optimize an array of positive real numbers by simulated annealing
by Andy Allinger, 2016-2020, released to the public domain
This program may be used by any person for any purpose.
----------------------------------------------------------------------*/
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
bool safdiv (float, float);

/*----------------------------------------------------------------------
ranorm - Random Approximate Normal Generator
origin:  G.E.P. Box and M.E. Muller, 1950.

refer to:
 Hand-Book on Statistical Distributions for experimentalists
 Christian Walck, Particle Physics Group, University of Stockholm
 10 September 2007

The random number generator is assumed to already be initialized.
----------------------------------------------------------------------*/
float ranorm (void)
{
/* Constants */
	static float TWOPI = 6.2831853;
/* Initialized data */
	static bool odd = false;
/* Local variables */
	float r, s;
	static float z1, z2;
/* begin */
	odd = ! odd;
	if (odd) {
L10:
		r = (float) rand() / (float) RAND_MAX;
		if (r <= 0.f || r > 1.f) goto L10;
		r = sqrtf(logf(r) * -2.f);
		if (r != r) goto L10;  /* just in case */
		s = TWOPI * (float) rand() / (float) RAND_MAX;
		z1 = r * sinf(s);
		z2 = r * cosf(s);
		return z1;
	} else {
		return z2;
	}
} /* end of ranorm */

/*----------------------------------------------------------------------
osiman - optimize an array of positive real numbers by simulated annealing

This is a reverse-communication interface.  The calling program must
call osiman to generate the next trial value of x[], then evaluate the
function e.

___Name_______Type_____In/Out____Description____________________________
   n          int      In        Size of problem
   x[n]       float*   Both      Parameters to optimize
   e          float    In        Evaluation of x
   sbeg       float    In        Initial scale (nondimensional)
   send       float    In        Final scale
   tbeg       float    In        Initial temperature (same units as e)
   tend       float    In        Final temperature
   nit        int      In        Number of iterations
   xkeep[n]   float*   Neither   Copy of x
   iflag      int*     Both      1 on input:  begin new problem
                                 0 on input:  continue
                                 0 on output:  successful step
----------------------------------------------------------------------*/
void osiman (int n, float *x, float e, float sbeg, float send,
    float tbeg, float tend, int nit, float *xkeep, int *iflag)
{
/* Constants */
static float BIG = 1e36f;
/* Local variables */
int j;
float h, prob;
static float blog, ekeep, r1, r2, s, t, samt, tfact;
/* Function Body */
	if (*iflag < 0) return;
	if (*iflag > 0) {
		samt = (send - sbeg) / (float)(nit - 1);  /* shrinkage schedule */
		r1 = tend / tbeg;
		r2 = 1.f / (float) (nit - 1);
		tfact = powf(r1, r2);  /* cooling schedule */
		blog = logf(BIG);
		s = sbeg;
		t = tbeg;
		ekeep = BIG;
		*iflag = 0;
		return;
	}
	h = ekeep - e;
	if (safdiv(h, t)) {
		h /= t;
		if (h < blog) {
			prob = expf(h);
		} else {
			prob = 1.f;
		}
	} else {
		if (h < 0.f) {
			prob = 0.f;
		} else {
			prob = 1.f;
		}
	}
	if ((float) rand() / (float) RAND_MAX < prob) {
		ekeep = e;
		for (j = 0; j < n; ++j) xkeep[j] = x[j];
	}
/* new trial point */
	for (j = 0; j < n; ++j) x[j] = xkeep[j] * expf(ranorm() * s);
	s += samt;
	t *= tfact;
	return;
} /* end of osiman */
