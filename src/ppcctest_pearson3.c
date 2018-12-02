/* Copyright (C) 2017 Thorsten Pohlert GPL-3 
 * Modified 2018-09-15
 * -- fix abs() to fabs()
 */
#include <R.h>
#include <Rmath.h>

void pmcor(double *x, double *y, int *n, double *res);

double get_pvalue(double *r, double *ppcc, int *nmc);

/* Probability Plot Correlation Coefficient Test
   for the Pearson 3 distribution */ 
void ppcctest_pearson3(double *p, double *ppcc, double *shape, 
		      int *n, int *nmc, double *pval)
{
  double *res; 
  int i , j;
  int mc = nmc[0];
  int nn = n[0];
  double *M, *qq, *r;
  double alpha = 4.0 / R_pow_di(shape[0], 2);
  double beta = 1.0 / 2.0 * fabs(shape[0]);
 
  /* allocate pointers */
  M = (double *) R_alloc(nn,  sizeof(double));
  qq = (double *) R_alloc(nn,  sizeof(double));
  r = (double *) R_alloc(mc,  sizeof(double));
  res = (double *) R_alloc(1,  sizeof(double));

  /* qpearson3 */
  /* initialize the M */
  if (shape[0] > 0.0) {
    for (i = 0; i < nn ; i++)
      {
	      M[i] = -alpha * beta + qgamma(p[i], alpha, beta, 1, 0);
      }
  } else {
    for (i = 0; i < nn ; i++)
      {
	      M[i] =  alpha * beta - qgamma(1.0 - p[i], alpha, beta, 1, 0);
      }
  }

  /* test once */
  if (shape[0] > 0.0) {

    /* Start Monte Carlo loop */
    GetRNGstate() ;
    for ( i = 0; i < mc; ++i)
      { 
	    /* get random pearson3 variates */
	     for ( j = 0; j < nn; j++)
	        {
	            qq[j] = -alpha * beta + qgamma(unif_rand(), 
					             alpha, beta, 1, 0);
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

  } else {
    /* Start Monte Carlo loop */
    GetRNGstate() ;
    for ( i = 0; i < mc; ++i)
      { 
	    /* get random pearson3 variates */
	   for ( j = 0; j < nn; j++)
	        {
	          qq[j] = alpha * beta - qgamma(1.0 - unif_rand(), 
					          alpha, beta, 1, 0);
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
}
