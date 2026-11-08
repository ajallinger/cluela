/* t_embed.c - embed cyclic data */
#include <stdbool.h>
#include <stdio.h>

int **new_int_mtx (int m, int n);  /* defined in array.c */
float **new_flt_mtx (int m, int n);
int del_int_mtx (int **b, int m);
int del_flt_mtx (float **x, int m);
int mds (float **x, int **nx, int n, int p, float *w, int *cat, bool robust,  /* defined in mds.c */
         float **v, int *m, float *e, float *t);

int main (void)
{
/* Constants */
#define P 1
#define N 60

/* Local variables */
	int cat[P], j, m, ifault;
	float e[N], t[N], w[P];
	bool robust;

/* Allcate matrices */
	int **nx;
	float **x, **v;
	nx = new_int_mtx(N, P);
	x = new_flt_mtx(N, P);
	v = new_flt_mtx(N, N);
	if (!nx || !x || !v) {
		printf("trouble allocating arrays!\n");
		goto L90;
	}
	w[0] = 1.f;
	robust = false;
	cat[0] = -1;
	for (j = 0; j < N; ++j) {
		x[j][0] = j / (float)N;
		nx[j][0] = 1;
	}
	ifault = mds (x, nx, N, P, w, cat, robust, v, &m, e, t);
	printf("mds returns #%i with %i real directions\n", ifault, m);
	printf("coordinates:\n");
	for (j = 0; j < N; ++j) printf("%f, %f\n", v[j][0], v[j][1]);

/* Release memory */
L90:
	if (nx) del_int_mtx(nx, N);
	if (x) del_flt_mtx(x, N);
	if (v) del_flt_mtx(v, N);
	return 0;
} /* end of t_embed */
