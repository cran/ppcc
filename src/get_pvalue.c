#include <R.h>

/* get the p-value Pr(r < R), from nmc * r values */
double get_pvalue(double *r, double *ppcc, int *nmc)
{
  int i, count = 0, mc = nmc[0];

  for (i = 0; i < mc; i++)
      {
	if (r[i] < ppcc[0] )
	  {
	    count += 1 ;
	  }
      }
    return (double) count / (double) mc;
}
