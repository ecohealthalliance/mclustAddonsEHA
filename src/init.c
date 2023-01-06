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
extern SEXP _mclustAddons_powerTransform(SEXP, SEXP, SEXP);
extern SEXP _mclustAddons_rangepowerTransformDeriv_lb(SEXP, SEXP, SEXP, SEXP);
extern SEXP _mclustAddons_rangepowerTransformDeriv_lub(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _mclustAddons_rangepowerTransformDeriv_unb(SEXP, SEXP);
extern SEXP _mclustAddons_rangeTransform(SEXP, SEXP, SEXP);
extern SEXP _mclustAddons_rowMax(SEXP);
extern SEXP _mclustAddons_rowSum(SEXP);
extern SEXP _mclustAddons_colMax(SEXP);
extern SEXP _mclustAddons_colSum(SEXP);
extern SEXP _mclustAddons_logsumexp_Rcpp(SEXP, SEXP);
extern SEXP _mclustAddons_softmax_Rcpp(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_mclustAddons_powerTransform",               (DL_FUNC) &_mclustAddons_powerTransform,               3},
    {"_mclustAddons_rangepowerTransformDeriv_lb",  (DL_FUNC) &_mclustAddons_rangepowerTransformDeriv_lb,  4},
    {"_mclustAddons_rangepowerTransformDeriv_lub", (DL_FUNC) &_mclustAddons_rangepowerTransformDeriv_lub, 6},
    {"_mclustAddons_rangepowerTransformDeriv_unb", (DL_FUNC) &_mclustAddons_rangepowerTransformDeriv_unb, 2},
    {"_mclustAddons_rangeTransform",               (DL_FUNC) &_mclustAddons_rangeTransform,               3},
    {"_mclustAddons_rowMax",                       (DL_FUNC) &_mclustAddons_rowMax,                       1},
    {"_mclustAddons_rowSum",                       (DL_FUNC) &_mclustAddons_rowSum,                       1},
    {"_mclustAddons_colMax",                       (DL_FUNC) &_mclustAddons_colMax,                       1},
    {"_mclustAddons_colSum",                       (DL_FUNC) &_mclustAddons_colSum,                       1},
    {"_mclustAddons_logsumexp_Rcpp",               (DL_FUNC) &_mclustAddons_logsumexp_Rcpp,               2},
    {"_mclustAddons_softmax_Rcpp",                 (DL_FUNC) &_mclustAddons_softmax_Rcpp,                 2},
    {NULL, NULL, 0}
};

void R_init_mclustAddons(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
