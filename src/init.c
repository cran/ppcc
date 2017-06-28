#include "ppcc.h"
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

/* Copyright (C) 2017 Thorsten Pohlert, GPL-3 */

/* define argument types */
static R_NativePrimitiveArgType pmcor_t[] = {
  REALSXP, REALSXP, INTSXP, REALSXP};

static R_NativePrimitiveArgType locscale_t[] = {
  REALSXP, REALSXP, INTSXP,  INTSXP, REALSXP};

static R_NativePrimitiveArgType locscaleshp_t[] = {
  REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP};

/* define Centry points, their names, nr of arguments and their types */

static const R_CMethodDef CEntries[]  = {
  {"pmcor", (DL_FUNC) &pmcor, 4, pmcor_t},
  {"ppcctest_norm", (DL_FUNC) &ppcctest_norm, 5, locscale_t},
  {"ppcctest_lnorm", (DL_FUNC) &ppcctest_lnorm, 5, locscale_t},
  {"ppcctest_unif", (DL_FUNC) &ppcctest_unif, 5, locscale_t},
  {"ppcctest_rayleigh", (DL_FUNC) &ppcctest_rayleigh, 5, locscale_t},
  {"ppcctest_exp", (DL_FUNC) &ppcctest_exp, 5, locscale_t},
  {"ppcctest_logis", (DL_FUNC) &ppcctest_logis, 5, locscale_t},
  {"ppcctest_cauchy", (DL_FUNC) &ppcctest_cauchy, 5, locscale_t},
  {"ppcctest_gumbel", (DL_FUNC) &ppcctest_gumbel, 5, locscale_t},
  {"ppcctest_pearson3", (DL_FUNC) &ppcctest_pearson3, 6, locscaleshp_t},
  {"ppcctest_weibull", (DL_FUNC) &ppcctest_weibull, 6, locscaleshp_t},
  {"ppcctest_gev", (DL_FUNC) &ppcctest_gev, 6, locscaleshp_t},
  {"ppcctest_glogis", (DL_FUNC) &ppcctest_glogis, 6, locscaleshp_t},
  {"ppcctest_kappa2", (DL_FUNC) &ppcctest_kappa2, 6, locscaleshp_t},
  {NULL, 0}
};


void attribute_visible R_init_ppcc(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
