#include <R.h>
#include <Rmath.h>

/*
double ppccsum(double *x, int *n){
  double xsum;
  xsum = 0.0;
  for (int i = 0; i < n[0]; ++i)
    {
      xsum += x[i];
    }
  return xsum ;
}
*/
/* mean calculation */
/*
double ppccsum(double *x, int *n);
double ppccmean(double *x, int *n)
{
  return ppccsum(x, n) / (double)n[0];
} 
*/
/*variance calculation */
/*
double ppccmean(double *x, int *n);
double ppccvar(double *x, int *n)
{
  double sum;
  double xmean = ppccmean(x, n);
  
  sum = 0.0;
  for (int i = 0; i < n[0] ; ++i)
    {
      sum += R_pow_di(x[i] - xmean, 2);
    }
  return sum;
}
*/

/*covariance calculation */
/*
double ppccmean(double *x, int *n);
double ppcccov(double *x, double *y, int *n)
{

  double sum;
  double xmean;
  double ymean;

  xmean = ppccmean(x, n);
  ymean = ppccmean(y, n);

  sum = 0.0;
  for (int i = 0; i < n[0]; ++i)
    {
      sum += (x[i] - xmean) * (y[i] - ymean);
    }

  return sum ;
}
*/

/* Pearson's corelation coefficient */ 
double ppcccov(double *x, double *y, int *n);
double ppccvar(double *x, int *n);
void pmcor(double *x, double *y, int *n, double *res)
{
  double xysum = 0.0, xsum = 0.0, ysum = 0.0;
  double x2sum = 0.0 , y2sum = 0.0 ;
  double deno, num;
  int i ;
  int nn = n[0];

  for ( i = 0; i < nn; ++i)
    {
      xysum += x[i] * y[i];
      xsum += x[i] ;
      ysum += y[i] ;
      x2sum += R_pow_di(x[i], 2) ;
      y2sum += R_pow_di(y[i], 2) ;
    }

  num = (double)nn * xysum - xsum * ysum;
  deno = ((double)nn * x2sum - R_pow_di(xsum, 2)) * 
    ((double)nn * y2sum - R_pow_di(ysum, 2));

  res[0] = num / sqrt (deno) ;
  

  /*res[0] = ppcccov(x, y, n) / sqrt(ppccvar(x, n) * ppccvar(y, n));*/
}


