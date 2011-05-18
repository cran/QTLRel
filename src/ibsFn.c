/* 
   estimate densed identity coefficients from genetype data
*/

/***************************************************
  compile with option:
     g++  -D_FILE_OFFSET_BITS=64 ibsFn.cc
     R CMD SHLIB -D_FILE_OFFSET_BITS=64 ibsFn.cc
****************************************************/

#include "xxx.h"

void deltaFn();

/*
   gdat: nr by nc array, genetoype data, 1-AA,2-AB,3-BB,0-missing
   ibs: nr(nr+1)/2 by 9 array, condensed identity coefficeints
   delta: nr(nr+1)/2 by 5 array (ksp,dlt1,dlt2,dlt35,dlt7)
*/
//extern "C"{
   void deltaFnc(int* gdat,int* nr,int* nc,double* delta){
         int nn;
         nn = (*nr)*((*nr)+1)/2;
         int* gdatp[(*nr)]; for(int i=0;i<(*nr);i++) gdatp[i]=gdat+i*(*nc);
         double* deltap[nn]; for(int i=0;i<nn;i++) deltap[i]=delta+i*5;
         deltaFn(gdatp,(*nr),(*nc),deltap);
   }
//}

//***********************************************
double pr(int *x,int k,int o){
   double s=1.0;
   if(o==1) for(int i=0;i<k;i++){
      s *= (3.0-x[i])/2; //both A alleles
   }else if(o==2) for(int i=0;i<k;i++){
      s *= (x[i]-1.0)/2; //both B alleles
   }else{
      Rprintf("o in pr: 1 or 2 only.\n");
      exit(0);
   }

   return s;
}

double _phi_2_(int a, int b, int** gdat,int j){
   int x[2];
   double s=0.0;
   x[0] = gdat[a][j];
   x[1] = gdat[b][j];
   if(x[0]*x[1] != 0){
      s += pr(x,2,1) + pr(x,2,2);
   }

   return s;
}

double dlt1(int a, int b, int** gdat,int j){
   double s=0.0;
   if((gdat[a][j]==1 && gdat[b][j]==1) || (gdat[a][j]==3 && gdat[b][j]==3)) s++;

   return s;
}

double dlt2(int a, int b, int** gdat,int j){
   double s=0.0;
   if((gdat[a][j]==1 && gdat[b][j]==3) || (gdat[a][j]==3 && gdat[b][j]==1)) s++;

   return s;
}

double dlt35(int a, int b, int** gdat,int j){
   double s=0.0;
   if((gdat[a][j]==2 || gdat[b][j]==2) && (gdat[a][j]!=gdat[b][j])) s++;

   return s;
}

double dlt7(int a, int b, int** gdat,int j){
   double s=0.0;
   if(gdat[a][j]==2 && gdat[b][j]==2) s++;

   return s;
}

void deltaFn(int** gdat,int nr,int nc,double** delta){
   int ii=0;
   int npairs = nr*(nr+1)/2;
   double s[5];

//   Rprintf("There are %d pairs:\n", npairs);
   for(int a=0;a<nr;a++){
      for(int b=a;b<nr;b++){
         for(int t=0; t<5; t++) s[t] = 0.0;
         for(int j=0; j<nc; j++){
            s[0] += _phi_2_(a,b,gdat,j);
            s[1] += dlt1(a,b,gdat,j);
            s[2] += dlt2(a,b,gdat,j);
            s[3] += dlt35(a,b,gdat,j);
            s[4] += dlt7(a,b,gdat,j);
         }
         for(int t=0; t<5; t++) delta[ii][t] = s[t]/nc;

         ii++;
//         Rprintf("%d/%d\r", ii, npairs);
      }
   }
}

