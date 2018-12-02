/* Copyright (C) 2017 Thorsten Pohlert GPL-3 */
#include <R.h>
#include <Rmath.h>

void pmcor(double *x, double *y, int *n, double *res);

double get_pvalue(double *r, double *ppcc, int *nmc);

/* Probability Plot Correlation Coefficient Test
   for the GEV distribution */ 
void ppcctest_gev(double *p, double *ppcc, double *shape, 
		      int *n, int *nmc, double *pval)
{
  double *res; 
  int i , j;
  int mc = nmc[0];
  int nn = n[0];
  double *M, *qq, *r;
  
 
  /* allocate pointers */
  M = (double *) R_alloc(nn,  sizeof(double));
  qq = (double *) R_alloc(nn,  sizeof(double));
  r = (double *) R_alloc(mc,  sizeof(double));
  res = (double *) R_alloc(1,  sizeof(double));

  /* qgev */
  /* initialize the M */
  for (i = 0; i < nn ; i++)
    {
      M[i] = 1.0 / shape[0] * 
	(1.0 - R_pow(-log(p[i]), shape[0]));
    }
 
  /* Start Monte Carlo loop */
  GetRNGstate() ;
  for ( i = 0; i < mc; ++i)
    { 
      /* get random gev variates */
      for ( j = 0; j < nn; j++)
	{
	  qq[j] = 1.0 / shape[0] * 
	    (1.0 - R_pow(-log(unif_rand()), shape[0]));
	}
      /* sort the array */
      R_rsort(qq, nn);
      /* calculate correlation */
      pmcor(qq, M, n, res) ;
      r[i] = res[0];
    }
  PutRNGstate();
  
  /* estimate p value */
  pval[0] = get_pvalue(r, ppcc, nmc);
}
