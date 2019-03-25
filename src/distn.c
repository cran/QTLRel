#include <R.h>
#include "nmath.h"
#include "dpq.h"
#include <limits.h>

/*
  Note: it doesn't work to pass back a value of the function
 */
void F77_SUB(pchif)(double* p, double* q, double* df, int* lower, int* log){
   *p = pchisq(*q, *df, *lower, *log);
   return;
}

void F77_SUB(pff)(double* p, double* q, double* df1, double* df2, int* lower, int* log){
   *p = pf(*q, *df1, *df2, *lower, *log);
   return;
}

double F77_SUB(normrnd)(void){//not working!!!
   return norm_rand();
}

