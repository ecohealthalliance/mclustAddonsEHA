#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*
Routines registration obtained with 

tools::package_native_routine_registration_skeleton(".", character_only = FALSE)
 
FIXME: Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _mclustAddonsEHA_powerTransform(SEXP, SEXP, SEXP);
extern SEXP _mclustAddonsEHA_rangepowerTransformDeriv_lb(SEXP, SEXP, SEXP, SEXP);
extern SEXP _mclustAddonsEHA_rangepowerTransformDeriv_lub(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _mclustAddonsEHA_rangepowerTransformDeriv_unb(SEXP, SEXP);
extern SEXP _mclustAddonsEHA_rangeTransform(SEXP, SEXP, SEXP);
extern SEXP _mclustAddonsEHA_rowMax(SEXP);
extern SEXP _mclustAddonsEHA_rowSum(SEXP);
extern SEXP _mclustAddonsEHA_colMax(SEXP);
extern SEXP _mclustAddonsEHA_colSum(SEXP);
extern SEXP _mclustAddonsEHA_logsumexp_Rcpp(SEXP, SEXP);
extern SEXP _mclustAddonsEHA_softmax_Rcpp(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_mclustAddonsEHA_powerTransform",               (DL_FUNC) &_mclustAddonsEHA_powerTransform,               3},
    {"_mclustAddonsEHA_rangepowerTransformDeriv_lb",  (DL_FUNC) &_mclustAddonsEHA_rangepowerTransformDeriv_lb,  4},
    {"_mclustAddonsEHA_rangepowerTransformDeriv_lub", (DL_FUNC) &_mclustAddonsEHA_rangepowerTransformDeriv_lub, 6},
    {"_mclustAddonsEHA_rangepowerTransformDeriv_unb", (DL_FUNC) &_mclustAddonsEHA_rangepowerTransformDeriv_unb, 2},
    {"_mclustAddonsEHA_rangeTransform",               (DL_FUNC) &_mclustAddonsEHA_rangeTransform,               3},
    {"_mclustAddonsEHA_rowMax",                       (DL_FUNC) &_mclustAddonsEHA_rowMax,                       1},
    {"_mclustAddonsEHA_rowSum",                       (DL_FUNC) &_mclustAddonsEHA_rowSum,                       1},
    {"_mclustAddonsEHA_colMax",                       (DL_FUNC) &_mclustAddonsEHA_colMax,                       1},
    {"_mclustAddonsEHA_colSum",                       (DL_FUNC) &_mclustAddonsEHA_colSum,                       1},
    {"_mclustAddonsEHA_logsumexp_Rcpp",               (DL_FUNC) &_mclustAddonsEHA_logsumexp_Rcpp,               2},
    {"_mclustAddonsEHA_softmax_Rcpp",                 (DL_FUNC) &_mclustAddonsEHA_softmax_Rcpp,                 2},
    {NULL, NULL, 0}
};

void R_init_mclustAddonsEHA(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
