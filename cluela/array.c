/* array.c - allocate and deallocate arrays
by Andy Allinger, 2017, released to the public domain
This program may be used by any person for any purpose. */

#include <stdlib.h>
int **new_int_mtx (int m, int n)
{
	int i;
	int **b;
	b = (int **) malloc(m * sizeof(int *));
	for (i = 0; i < m; ++i) {
		b[i] = malloc(n * sizeof(int));
		if (!b[i]) break;
	}
	if (i < m) {
		for (i = i; i >= 0; --i) free (b[i]);
		free (b);
		b = NULL;
	}
	return b;
}

float **new_flt_mtx (int m, int n)
{
	int i;
	float **x;
	x = (float **) malloc(m * sizeof(float *));
	for (i = 0; i < m; ++i) {
		x[i] = malloc(n * sizeof(float));
		if (!x[i]) break;
	}
	if (i < m) {
		for (i = i; i >= 0; --i) free (x[i]);
		free (x);
		x = NULL;
	}
	return x;
}

char **new_ch_mtx (int m, int n)
{
	int i;
	char **h;
	h = (char **) malloc(m * sizeof(char *));
	for (i = 0; i < m; ++i) {
		h[i] = malloc(n * sizeof(char));
		if (!h[i]) break;
	}
	if (i < m) {
		for (i = i; i >= 0; --i) free (h[i]);
		free (h);
		h = NULL;
	}
	return h;
}

int del_int_mtx (int **b, int m)
{
	int i;
	for (i = 0; i < m; ++i) free (b[i]);
	free (b);
	return 0;
}

int del_flt_mtx (float **x, int m)
{
	int i;
	for (i = 0; i < m; ++i) free (x[i]);
	free (x);
	return 0;
}

int del_ch_mtx (char **h, int m)
{
	int i;
	for (i = 0; i < m; ++i) free (h[i]);
	free (h);
	return 0;
}
