/*----------------------------------------------------------------------
eigen.c - Solve real symmetric eigenproblems.  Translated from EISPACK,
as distributed in the National Bureau of Standards Core Mathematical
Library.  This file is public domain.
----------------------------------------------------------------------*/
#include <math.h>

/*----------------------------------------------------------------------
  TRED2
  DATE WRITTEN   760101   (YYMMDD)
  REVISION DATE  830518   (YYMMDD)
  CATEGORY NO.  D4C1B1
  AUTHOR  Smith, B. T., et al.
  PURPOSE  Reduce real symmetric matrix to symmetric tridiagonal
           matrix using and accumulating orthogonal transformation
  DESCRIPTION

    This subroutine is a translation of the ALGOL procedure TRED2,
    NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
    Handbook for Auto. Comp., Vol.II-Linear Algebra, 212-226(1971).

    This subroutine reduces a REAL SYMMETRIC matrix to a
    symmetric tridiagonal matrix using and accumulating
    orthogonal similarity transformations.

    On Input

       N is the order of the matrix.

       A contains the real symmetric input matrix.  Only the
         lower triangle of the matrix need be supplied.

    On Output

       D contains the diagonal elements of the tridiagonal matrix.

       E contains the subdiagonal elements of the tridiagonal
         matrix in its last N-1 positions.  E(0) is set to zero.

       Z contains the orthogonal transformation matrix
         produced in the reduction.

       A and Z may coincide.  If distinct, A is unaltered.

    Questions and comments should be directed to B. S. Garbow,
    Applied Mathematics Division, Argonne National Laboratory
    ------------------------------------------------------------------
  REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
                Y. Ikebe, V. C. Klema, C. B. Moler, *Matrix Eigen-
                system Routines - EISPACK Guide*, Springer-Verlag,
                1976.
  ROUTINES CALLED  (NONE)
----------------------------------------------------------------------- */
void tred2 (int n, float **a, float *d, float *e, float **z)
{
/* Local variables */
	float f, g, h;
	int i, j, k, l;
	float hh;
	int ii, jp1;
	float scale;
	float tmp;
/* Function Body */
	for (i = 0; i < n; ++i) {
		for (j = 0; j <= i; ++j) {
			z[i][j] = a[i][j];
		}
	}
	if (1 == n) goto L320;
/* .......... FOR I=N STEP -1 UNTIL 2 DO -- .......... */
	for (ii = 2; ii <= n; ++ii) {
		i = n + 2 - ii;
		l = i - 1;
		h = 0.f;
		scale = 0.f;
		if (l < 2) goto L130;
/* .......... Scale row (ALGOL TOL then not needed) .......... */
		for (k = 1; k <= l; ++k) scale += fabsf(z[i-1][k-1]);
		if (scale != 0.f)  goto L140;
L130:
		e[i-1] = z[i-1][l-1];
		goto L290;
L140:
		for (k = 1; k <= l; ++k) {
			z[i-1][k-1] /= scale;
			tmp = z[i-1][k-1];
			h += tmp * tmp;
		}
		f = z[i-1][l-1];
		g = -copysignf(sqrtf(h), f);
		e[i-1] = scale * g;
		h -= f * g;
		z[i-1][l-1] = f - g;
		f = 0.f;
		for (j = 1; j <= l; ++j) {
			z[j-1][i-1] = z[i-1][j-1] / h;
			g = 0.f;
/* .......... Form element of A*U .......... */
			for (k = 1; k <= j; ++k) g += z[j-1][k-1] * z[i-1][k-1];
			jp1 = j + 1;
			if (l < jp1) goto L220;
			for (k = jp1; k <= l; ++k) g += z[k-1][j-1] * z[i-1][k-1];
/* .......... Form element of P .......... */
L220:
			e[j-1] = g / h;
			f += e[j-1] * z[i-1][j-1];
		}
		hh = f / (h + h);
/* .......... Form reduced A .......... */
		for (j = 1; j <= l; ++j) {
			f = z[i-1][j-1];
			g = e[j-1] - hh * f;
			e[j-1] = g;
			for (k = 1; k <= j; ++k)
			    z[j-1][k-1] = z[j-1][k-1] - f * e[k-1] - g * z[i-1][k-1];
		}
L290:
		d[i-1] = h;
	}
L320:
	d[0] = 0.f;
	e[0] = 0.f;
/* .......... Accumulation of transformation matrices .......... */
	for (i = 1; i <= n; ++i) {
		l = i - 1;
		if (0.f == d[i-1]) goto L380;
		for (j = 1; j <= l; ++j) {
			g = 0.f;
			for (k = 1; k <= l; ++k) g += z[i-1][k-1] * z[k-1][j-1];
			for (k = 1; k <= l; ++k) z[k-1][j-1] -= g * z[k-1][i-1];
		}
L380:
		d[i-1] = z[i-1][i-1];
		z[i-1][i-1] = 1.f;
		if (l < 1) continue;
		for (j = 1; j <= l; ++j) {
			z[i-1][j-1] = 0.f;
			z[j-1][i-1] = 0.f;
		}
	}
	return;
} /* end of tred2 */


/*----------------------------------------------------------------------
  TQL2
  DATE WRITTEN   760101   (YYMMDD)
  REVISION DATE  830518   (YYMMDD)
  CATEGORY NO.  D4A5,D4C2A
  AUTHOR  Smith, B. T., et al.
  PURPOSE  Compute eigenvalues and eigenvectors of symmetric
           tridiagonal matrix.
  DESCRIPTION

    This subroutine is a translation of the ALGOL procedure TQL2,
    Num. Math. 11, 293-306(1968) by Bowdler, Martin, Reinsch, and
    Wilkinson.
    Handbook for Auto. Comp., Vol.II-Linear Algebra, 227-240(1971).

    This subroutine finds the eigenvalues and eigenvectors
    of a SYMMETRIC TRIDIAGONAL matrix by the QL method.
    The eigenvectors of a FULL SYMMETRIC matrix can also
    be found if  TRED2  has been used to reduce this
    full matrix to tridiagonal form.

    On Input

       N is the order of the matrix.

       D contains the diagonal elements of the input matrix.

       E contains the subdiagonal elements of the input matrix
         in its last N-1 positions.  E(0) is arbitrary.

       Z contains the transformation matrix produced in the
         reduction by  TRED2, if performed.  If the eigenvectors
         of the tridiagonal matrix are desired, Z must contain
         the identity matrix.

     On Output

       D contains the eigenvalues in ascending order.  If an
         error exit is made, the eigenvalues are correct but
         unordered for indices 1,2,...,IERR-1.

       E has been destroyed.

       Z contains orthonormal eigenvectors of the symmetric
         tridiagonal (or full) matrix.  If an error exit is made,
         Z contains the eigenvectors associated with the stored
         eigenvalues.

       Return value is set to
         Zero       for normal return,
         J          if the J-th eigenvalue has not been
                    determined after 30 iterations.

    Questions and comments should be directed to B. S. Garbow,
    Applied Mathematics Division, Argonne National Laboratory
    ------------------------------------------------------------------
  REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
                Y. Ikebe, V. C. Klema, C. B. Moler, *Matrix Eigen-
                system Routines - EISPACK Guide*, Springer-Verlag,
                1976.
----------------------------------------------------------------------*/
int tql2 (int n, float *d, float *e, float **z)
{
/* Local variables */
	float b, c, f, g, h;
	int i, j, k, l, m;
	float p, r, s, c2, c3;
	int l1, l2;
	float s2;
	int ii;
	float dl1, el1;
	int mml;
	float one;
/* Function Body */
	one = 1.f;
	if (1 == n) return 0;
	for (i = 2; i <= n; ++i) e[i-2] = e[i-1];
	f = 0.f;
	b = 0.f;
	e[n-1] = 0.f;
	s2 = 0.f;
	c3 = 0.f;
	for (l = 1; l <= n; ++l) {
		j = 0;
		h = fabsf(d[l-1]) + fabsf(e[l-1]);
		if (b < h) b = h;
/* .......... Look for small sub-diagonal element .......... */
		for (m = l; m <= n; ++m) {
			if (b + fabsf(e[m-1]) == b) break;
/* .......... E(N) is always zero, so there is no exit */
/*            through the bottom of the loop .......... */
/* added 2011/06/13 AJA */
			if (m > n) return -1;
		}
		if (m == l)  goto L220;
L130:
/* .......... Set error -- No convergence to an */
/*            eigenvalue after 30 iterations .......... */
		if (30 == j) return l;
		++j;
/* .......... Form shift .......... */
		l1 = l + 1;
		l2 = l1 + 1;
		g = d[l-1];
		p = (d[l1-1] - g) / (e[l-1] * 2.f);
		r = hypotf(p, one);
		d[l-1] = e[l-1] / (p + copysignf(r, p));
		d[l1-1] = e[l-1] * (p + copysignf(r, p));
		dl1 = d[l1-1];
		h = g - d[l-1];
		if (l2 > n) goto L145;
		for (i = l2; i <= n; ++i) d[i-1] -= h;
L145:
		f += h;
/* .......... QL transformation .......... */
		p = d[m-1];
		c = 1.f;
		c2 = c;
		el1 = e[l1-1];
		s = 0.f;
		mml = m - l;
/* .......... FOR I=M-1 STEP -1 UNTIL L DO -- .......... */
		for (ii = 1; ii <= mml; ++ii) {
			c3 = c2;
			c2 = c;
			s2 = s;
			i = m - ii;
			g = c * e[i-1];
			h = c * p;
			if (fabsf(p) < fabsf(e[i-1])) goto L150;
			c = e[i-1] / p;
			r = sqrtf(c * c + 1.f);
			e[i] = s * p * r;
			s = c / r;
			c = 1.f / r;
			goto L160;
L150:
			c = p / e[i-1];
			r = sqrtf(c * c + 1.f);
			e[i] = s * e[i-1] * r;
			s = 1.f / r;
			c *= s;
L160:
			p = c * d[i-1] - s * g;
			d[i] = h + s * (c * g + s * d[i-1]);
/* .......... Form vector .......... */
			for (k = 1; k <= n; ++k) {
				h = z[k-1][i];
				z[k-1][i] = s * z[k-1][i-1] + c * h;
				z[k-1][i-1] = c * z[k-1][i-1] - s * h;
			}
		}
		p = -s * s2 * c3 * el1 * e[l-1] / dl1;
		e[l-1] = s * p;
		d[l-1] = c * p;
		if (b + fabsf(e[l-1]) > b) goto L130;
L220:
		d[l-1] += f;
	}
/* .......... Order eigenvalues and eigenvectors .......... */
	for (ii = 2; ii <= n; ++ii) {
		i = ii - 1;
		k = i;
		p = d[i-1];
		for (j = ii; j <= n; ++j) {
			if (d[j-1] >= p) continue;
			k = j;
			p = d[j-1];
		}
		if (k == i) continue;
		d[k-1] = d[i-1];
		d[i-1] = p;
		for (j = 1; j <= n; ++j) {
			p = z[j-1][i-1];
			z[j-1][i-1] = z[j-1][k-1];
			z[j-1][k-1] = p;
		}
	}
	return 0;
} /* end of tql2 */

/*----------------------------------------------------------------------
       Eigenvalues and eigenvectors of symmetric real matrices

___Name______Type______I/O_______Description____________________________
   a[n][n]   float**   In        The matrix
   n         int       In        Order of a
   w[n]      float*    Out       Eigenvalues
   z[n][n]   float**   Out       Eigenvectors
   t[n]      float*    Neither   Workspace
----------------------------------------------------------------------*/
int seigv (float **a, int n, float *w, float **z, float *t)
{
/* Local variables */
	int ierr;
/* Function Body */
	tred2 (n, a, w, t, z);
	ierr = tql2(n, w, t, z);
	return ierr;
} /* end of seigv */
