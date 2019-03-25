
#include "QTLRelR.h"
#include <R_ext/Rdynload.h> //R_CMethodDef
//#include <Rdefines.h>

//ibsFn.c
extern void ibsPrc(double*, int*, int*, double*);
extern void ibsFnc(int*, int*, int*, double*);
extern void deltaFnc(int*, int*, int*, double*);
//idcoef.c
extern void llints(int*);
extern void getsize(int*);
extern void kinshipc(int*, int*, int*, double*);
extern void phicw(int*, int*, int*, int*, int*, int*, char**, char**);
extern void phicr(int*, int*, int*, int*, int*, int*, char**, double*, int*);
extern void gen_Matrix(double*, int*, int*, int*, double*, double*, double*, double*, double*);
//imFn.c
extern void conGenoPrc(int*, int*, double*, double*, int*, int*, int*, int*, double*, int*);
//ks.c
extern void pkolmogorov2x(double *x, Sint *n);
extern void psmirnov2x(double *x, Sint *m, Sint *n);
//qqplot.c
extern void Fnc(double*, int*, double*, int*);
extern void kolm(double*, int*);
extern void qFnc(double*, int*, double*, int*);
//rgdat.c
extern void rgdata(int*, int*, int*, int*, int*, double*);
extern void rgdata2(int*, int*, int*, int*, int*, double*, int*);

static const R_CMethodDef cMethods[] = {
    {"ibsPrc",          (DL_FUNC) &ibsPrc,             4},
    {"ibsFnc",          (DL_FUNC) &ibsFnc,             4},
    {"deltaFnc",        (DL_FUNC) &deltaFnc,           4},
    {"llints",          (DL_FUNC) &llints,             1},
    {"getsize",         (DL_FUNC) &getsize,            1},
    {"kinshipc",        (DL_FUNC) &kinshipc,           4},
    {"phicw",           (DL_FUNC) &phicw,              8},
    {"phicr",           (DL_FUNC) &phicr,              9},
    {"gen_Matrix",      (DL_FUNC) &gen_Matrix,         9},
    {"conGenoPrc",      (DL_FUNC) &conGenoPrc,        10},
    {"pkolmogorov2x",   (DL_FUNC) &pkolmogorov2x,      2},
    {"psmirnov2x",      (DL_FUNC) &psmirnov2x,         3},
    {"Fnc",             (DL_FUNC) &Fnc,                4},
    {"kolm",            (DL_FUNC) &kolm,               2},
    {"qFnc",            (DL_FUNC) &qFnc,               4},
    {"rgdata",          (DL_FUNC) &rgdata,             6},
    {"rgdata2",         (DL_FUNC) &rgdata2,            7},
    {NULL, NULL, 0}
};

static const
R_FortranMethodDef fMethods[] = {
   {"sc10", (DL_FUNC) &F77_NAME(sc10), 23},
   {"sc11", (DL_FUNC) &F77_NAME(sc11), 25},
   {"sc20", (DL_FUNC) &F77_NAME(sc20), 24},
   {"sc21", (DL_FUNC) &F77_NAME(sc21), 26},
   {"dsyev", (DL_FUNC) &F77_NAME(dsyev), 9},
   {"dsyevr", (DL_FUNC) &F77_NAME(dsyevr), 21},
   {NULL, NULL, 0}
};

static const
R_CallMethodDef callMethods[] = {
   //{"Cdqrls", (DL_FUNC) &Cdqrls, 4},
   NULL
};

void R_init_QTLRel(DllInfo *dll)
{
    R_registerRoutines(dll, cMethods, callMethods, fMethods, NULL);
    R_useDynamicSymbols(dll, FALSE);
//    R_forceSymbols(dll, TRUE);
}

