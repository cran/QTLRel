
/***************************************************
  compile with option:
       g++  -D_FILE_OFFSET_BITS=64 compareFiles.cc
       R CMD SHLIB -D_FILE_OFFSET_BITS=64 idcoef.cc
  or add on top "#define _FILE_OFFSET_BITS 64"
  --which seems crucial for large file writing:
       g++ compareFiles.cc
       R CMD SHLIB idcoef.cc
  long long type can hold an integer as large as
  ULLONG_MAX 18446744073709551615ULL
****************************************************/

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

