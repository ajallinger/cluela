/*----------------------------------------------------------------------
dump.c - Write CSV scatter plot, with clusters by column
by Andy Allinger, 2017, released to the public domain
This program may be used by any person for any purpose.

___Name______Type______I/O___Description________________________________
   r[n][m]   float**   In    Data points projected into two dimensions
   n         int       In    Number of objects
   z[n]      int*      In    Cluster ID's
   iou       FILE*     In    Unit number of an open I/O unit

              Return value:
                              0 = no errors
                              1 = Z out of bounds
                              other as returned by WRITE statements
----------------------------------------------------------------------*/
#include <stdio.h>

int wcsvsp (float **r, int n, int *z, FILE *iou)
{
/* Local variables */
	int i, l, ifault;
/* Function Body */
	for (i = 0; i < n; ++i) if (z[i] < 1) return 1;
/*   One record for each object */
	for (i = 0; i < n; ++i) {
/* First coordinate */
		ifault = fprintf (iou, "%9.3g", (double)r[i][0]);
		if (ifault < 0) {
			printf ("wcsvsp: trouble writing X-coordinate!\n");
			return ifault;
		}
/* One comma per cluster ID number */
		for (l = 0; l < z[i]; ++l) {
			ifault = fprintf (iou, ",");
			if (ifault < 0) {
				printf ("wcsvsp: trouble writing comma\n");
				return ifault;
			}
	}
/* Second coordinate */
		ifault = fprintf (iou, "%9.3g\n", (double)r[i][1]);
		if (ifault < 0) {
			printf ("wcsvsp: trouble writing Y-coordinate!\n");
			return ifault;
		}
	}
	return 0;
} /* end of wcsvsp */
