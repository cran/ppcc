#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "ppcc.h"

static const R_CMethodDef CEntries[] = {
    {"pmcor",             (DL_FUNC) &pmcor,             4},
    {"ppcctest_cauchy",   (DL_FUNC) &ppcctest_cauchy,   5},
    {"ppcctest_exp",      (DL_FUNC) &ppcctest_exp,      5},
    {"ppcctest_gev",      (DL_FUNC) &ppcctest_gev,      6},
    {"ppcctest_glogis",   (DL_FUNC) &ppcctest_glogis,   6},
    {"ppcctest_gumbel",   (DL_FUNC) &ppcctest_gumbel,   5},
    {"ppcctest_kappa2",   (DL_FUNC) &ppcctest_kappa2,   6},
    {"ppcctest_lnorm",    (DL_FUNC) &ppcctest_lnorm,    5},
    {"ppcctest_logis",    (DL_FUNC) &ppcctest_logis,    5},
    {"ppcctest_norm",     (DL_FUNC) &ppcctest_norm,     5},
    {"ppcctest_pearson3", (DL_FUNC) &ppcctest_pearson3, 6},
    {"ppcctest_rayleigh", (DL_FUNC) &ppcctest_rayleigh, 5},
    {"ppcctest_unif",     (DL_FUNC) &ppcctest_unif,     5},
    {"ppcctest_weibull",  (DL_FUNC) &ppcctest_weibull,  6},
    {NULL, NULL, 0}
};

void R_init_ppcc(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
