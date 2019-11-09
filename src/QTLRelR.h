
/**************************
  compile with option:
       g++  -D_FILE_OFFSET_BITS=64 compareFiles.cc
       R CMD SHLIB -D_FILE_OFFSET_BITS=64 idcoef.cc
  or add on top "#define _FILE_OFFSET_BITS 64"
  --which seems crucial for large file writing:
       g++ compareFiles.cc
       R CMD SHLIB idcoef.cc
  long long type can hold an integer as large as
  ULLONG_MAX 18446744073709551615ULL
**************************/

#define _FILE_OFFSET_BITS 64 //must on top

#include <R_ext/Error.h>
#include <R_ext/Memory.h>
#ifdef ENABLE_NLS
   #include <libintl.h>
   #define _(String) dgettext ("stats", String)
#else
   #define _(String) (String)
#endif
#include <R.h>
#include <Rmath.h>
#include <R_ext/Print.h> // for Rprintf
#include <R_ext/Utils.h> //R_CheckUserInterrupt(void)
#include <R_ext/RS.h>
//#include <R_ext/BLAS.h>
//#include <R_ext/Lapack.h>

#include "stdlib.h"
#include "time.h"
#include "math.h"
#include "limits.h"
#include "stdio.h"
#include "string.h"
#include "signal.h"
//#include "iostream"
//#include "fstream"
//#include "iomanip"
//#include "string"
//#include "cmath"
//using namespace std;

typedef long long LONGLONG;

void F77_NAME(sc10)(double* y, int* n, double* x, int* p,
   double* xg, int* ng,
   double* coef, double* tt, double* pval, double* v, int* opt,
   int* pvt, double* b, double* r0, double* r1,
   double* xx, double* qty, double* qraux, double* work);

void F77_NAME(sc11)(double* y, int* n, double* x, int* p, int* nc,
   double* xg, int* ng, double* gcv,
   double* coef, double* tt, double* pval, double* v, int* opt,
   int* pvt, double* b, double* r0, double* r1,
   double* xx, double* qty, double* qraux, double* work);

void F77_NAME(sc20)(double* y, int* n, double* x, int* p,
   double* xg, int* ng, int* nl,
   double* coef, double* tt, double* pval, double* v, int* opt,
   int* pvt, double* b, double* r0, double* r1,
   double* xx, double* qty, double* qraux, double* work);

void F77_NAME(sc21)(double* y, int* n, double* x, int* p, int* nc,
   double* xg, int* ng, int* nl, double* gcv,
   double* coef, double* tt, double* pval, double* v, int* opt,
   int* pvt, double* b, double* r0, double* r1,
   double* xx, double* qty, double* qraux, double* work);

void F77_NAME(dsyev)(char* JOBZ, char* UPLO, int* N, double* A, int* LDA,
   double* W, double* WORK, int* LWORK, int* INFO);

void F77_NAME(dsyevr)(char* JOBZ, char* RANGE, char* UPLO, int* N, double* A, int* LDA,
   double* VL, double* VU, int* IL, int* IU, double* ABSTOL, int* M, double* W,
   double* Z, int* LDZ, int* ISUPPZ, double* WORK, int* LWORK, int* IWORK, int* LIWORK,
   int* INFO);

