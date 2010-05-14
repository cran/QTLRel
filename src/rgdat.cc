
/***************************************************
  compile with option:
     g++  -D_FILE_OFFSET_BITS=64 rgdat.cc
     R CMD SHLIB -D_FILE_OFFSET_BITS=64 rgdat.cc
****************************************************/

#include <R_ext/Error.h>
#include <R_ext/Memory.h>
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("stats", String)
#else
#define _(String) (String)
#endif

#include "iostream"
#include "fstream"
#include "iomanip"
#include "string"
#include "stdlib.h"
#include "time.h"
#include "math.h"
#include "limits.h"
#include "cmath"
#include <stdio.h>
using namespace std;
#define _FILE_OFFSET_BITS 64

void rgeno1(int* data,int nr,int nc,int ninit,int* pedd,double* rr,unsigned long seed=0);
void rgeno2(int* data,int nr,int nc,int ninit,int* pedd,double* rr,bool xchr=true,unsigned long seed=0);
extern "C"{
   void rgdata(int* gdata,int &nr,int &nc,int& ninit,int* pedigree,double* recomb,unsigned long &seed){
      rgeno1(gdata,nr,nc,ninit,pedigree,recomb,seed);
   }
   void rgdata2(int* gdata,int &nr,int &nc,int& ninit,int* pedigree,double* recomb,bool &xchr,unsigned long &seed){
      rgeno2(gdata,nr,nc,ninit,pedigree,recomb,xchr,seed);
   }
}

/*---------------------------------------------------
 generate genotype data using pedigree
 gdata: nr by 2*nc arrary
    gdata[i,] -- (L1A1,L1A2,...,LjA1,LjA2,...) for individual i
    gdata[1,] and gdata[2,] represent F1 parents with known genotypes
 pedd: nr by 4 array, pedd[i,] -- (id,sire,dam,sex)
    id=1,2,...,nr
 recomb: vector of length nc, recomb[j] is the recombination rate
    between SNP[j-1] and SNP[j]
 seed: take system time as seed if 0
 ---------------------------------------------------*/
void rgeno1(int* data,int nr,int nc,int ninit,int* pedd,double* rr,unsigned long seed){
   double u;
   int i,j;
   int ii,jj;
   if(seed == 0) srand ( time(NULL) ); else srand(seed); //printf("%d\n",seed);
   if(nr<2){
//      printf("pedigree: at least 2 rows.\n"); exit(0);
      error(_("pedigree: at least 2 rows.\n"));
   }else if(nr>INT_MAX){
//      printf("pedigree: too many rows.\n"); exit(0);
      error(_("pedigree: too many rows.\n"));
   }
   if(nc<1){
//      printf("recombinaton rate: at least 1 SNP.\n"); exit(0);
      error(_("recombinaton rate: at least 1 SNP.\n"));
   }else if(nc>(INT_MAX-1)/2){
//      printf("recombinaton rate: too many SNPs.\n"); exit(0);
      error(_("recombinaton rate: too many SNPs.\n"));
   }
   for(i=ninit;i<nr;i++){
      //paternal gamete
      ii = pedd[i*4+1] - 1;
      if(ii>=0){
         u =  rand()/((double)RAND_MAX + 1);
         if(u<0.5) jj = 0; else jj = 1;
         data[i*nc*2] = data[ii*nc*2+jj];
         for(j=1;j<nc;j++){
            u =  rand()/((double)RAND_MAX + 1);
            if(u<rr[j]) jj += 1; jj %= 2;
            data[i*nc*2+j*2] = data[ii*nc*2+j*2+jj];
         }
      }
      //maternal gamete
      ii = pedd[i*4+2] - 1;
      if(ii>=0){
         u =  rand()/((double)RAND_MAX + 1);
         if(u<0.5) jj = 0; else jj = 1;
         data[i*nc*2+1] = data[ii*nc*2+jj];
         for(j=1;j<nc;j++){
            u =  rand()/((double)RAND_MAX + 1);
            if(u<rr[j]) jj += 1; jj %= 2;
            data[i*nc*2+j*2+1] = data[ii*nc*2+j*2+jj];
         }
      }
   }
}


/*---------------------------------------------------
 generate genotype data for ONE chromosome using pedigree
 gdata: nr by 2*nc arrary
    gdata[i,] -- (L1A1,L1A2,...,LjA1,LjA2,...) for individual i
    gdata[1,] and gdata[2,] represent F1 parents with known genotypes
 pedd: nr by 4 array, pedd[i,] -- (id,sire,dam,sex)
    id=1,2,...,nr
 recomb: vector of length nc, recomb[j] is the recombination rate
    between SNP[j-1] and SNP[j]
 xchr: true if it is x-chromosome
 seed: take system time as seed if 0
 ---------------------------------------------------*/

void rgeno2(int* data,int nr,int nc,int ninit,int* pedd,double* rr,bool xchr,unsigned long seed){
   double u;
   int i,j;
   int ii,jj;
   if(seed == 0) srand ( time(NULL) ); else srand(seed); //printf("%d\n",seed);
   if(nr<2){
//      printf("pedigree: at least 2 rows.\n"); exit(0);
      error(_("pedigree: at least 2 rows.\n"));
   }else if(nr>INT_MAX){
//      printf("pedigree: too many rows.\n"); exit(0);
      error(_("pedigree: too many rows.\n"));
   }
   if(nc<1){
//      printf("recombinaton rate: at least 1 SNP.\n"); exit(0);
      error(_("recombinaton rate: at least 1 SNP.\n"));
   }else if(nc>(INT_MAX-1)/2){
//      printf("recombinaton rate: too many SNPs.\n"); exit(0);
      error(_("recombinaton rate: too many SNPs.\n"));
   }
   if(xchr){
      for(i=ninit;i<nr;i++){
         //paternal gamete
         ii = pedd[i*4+1] - 1;
         if(ii>=0){
            if(pedd[i*4+3]==0) jj = 1; else jj = 0;
            data[i*nc*2] = data[ii*nc*2+jj];
            for(j=1;j<nc;j++){
               data[i*nc*2+j*2] = data[ii*nc*2+j*2+jj];
            }
         }
         //maternal gamete
         ii = pedd[i*4+2] - 1;
         if(ii>=0){
            u =  rand()/((double)RAND_MAX + 1);
            if(u<0.5) jj = 0; else jj = 1;
            data[i*nc*2+1] = data[ii*nc*2+jj];
            for(j=1;j<nc;j++){
               u =  rand()/((double)RAND_MAX + 1);
               if(u<rr[j]) jj += 1; jj %= 2;
               data[i*nc*2+j*2+1] = data[ii*nc*2+j*2+jj];
            }
         }
      }
   }else{
      for(i=ninit;i<nr;i++){
         //paternal gamete
         ii = pedd[i*4+1] - 1;
         if(ii>=0){
            u =  rand()/((double)RAND_MAX + 1);
            if(u<0.5) jj = 0; else jj = 1;
            data[i*nc*2] = data[ii*nc*2+jj];
            for(j=1;j<nc;j++){
               u =  rand()/((double)RAND_MAX + 1);
               if(u<rr[j]) jj += 1; jj %= 2;
               data[i*nc*2+j*2] = data[ii*nc*2+j*2+jj];
            }
         }
         //maternal gamete
         ii = pedd[i*4+2] - 1;
         if(ii>=0){
            u =  rand()/((double)RAND_MAX + 1);
            if(u<0.5) jj = 0; else jj = 1;
            data[i*nc*2+1] = data[ii*nc*2+jj];
            for(j=1;j<nc;j++){
               u =  rand()/((double)RAND_MAX + 1);
               if(u<rr[j]) jj += 1; jj %= 2;
               data[i*nc*2+j*2+1] = data[ii*nc*2+j*2+jj];
            }
         }
      }
   }
}


