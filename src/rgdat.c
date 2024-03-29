
#include "QTLRelR.h"
#include <R_ext/Utils.h> //R_CheckUserInterrupt(void)

void rgeno1(int* data,int nr,int nc,int ninit,int* pedd,double* rr);
void rgeno2(int* data,int nr,int nc,int ninit,int* pedd,double* rr,int xchr);
//extern "C"{
   void rgdata(int* gdata,int *nr,int *nc,int* ninit,int* pedigree,double* recomb){
      rgeno1(gdata,*nr,*nc,*ninit,pedigree,recomb);
   }
   void rgdata2(int* gdata,int *nr,int *nc,int* ninit,int* pedigree,double* recomb,int *xchr){
      rgeno2(gdata,*nr,*nc,*ninit,pedigree,recomb,*xchr);
   }
//}

/*---------------------------------------------------
 generate genotype data using pedigree
 gdata: nr by 2*nc arrary
    gdata[i,] -- (L1A1,L1A2,...,LjA1,LjA2,...) for individual i
    gdata[1,] and gdata[2,] represent F1 parents with known genotypes
 pedd: nr by 4 array, pedd[i,] -- (id,sire,dam,sex)
    id=1,2,...,nr
 recomb: vector of length nc, recomb[j] is the recombination rate
    between SNP[j-1] and SNP[j]
 ---------------------------------------------------*/
void rgeno1(int* data,int nr,int nc,int ninit,int* pedd,double* rr){
   double u;
   int i,j;
   int ii,jj;
   if(nr<2){
      error(_("pedigree: at least 2 rows.\n"));
   }else if(nr>INT_MAX){//useless, R won't pass int larger than INT_MAX
      error(_("pedigree: too many rows.\n"));
   }
   if(nc<1){
      error(_("recombinaton rate: at least 1 SNP.\n"));
   }else if(nc>(INT_MAX-1)/2){
      error(_("recombinaton rate: too many SNPs.\n"));
   }
   for(i=ninit;i<nr;i++){
      R_CheckUserInterrupt();
      //paternal gamete
      ii = pedd[i*4+1] - 1;
      if(ii>=0){
         GetRNGstate();
         u =  unif_rand();
         PutRNGstate();
         if(u<0.5) jj = 0; else jj = 1;
         data[i*nc*2] = data[ii*nc*2+jj];
         for(j=1;j<nc;j++){
            GetRNGstate();
            u =  unif_rand();
            PutRNGstate();
            if(u<rr[j]){jj += 1; jj %= 2;}
            data[i*nc*2+j*2] = data[ii*nc*2+j*2+jj];
         }
      }
      //maternal gamete
      ii = pedd[i*4+2] - 1;
      if(ii>=0){
         GetRNGstate();
         u =  unif_rand();
         PutRNGstate();
         if(u<0.5) jj = 0; else jj = 1;
         data[i*nc*2+1] = data[ii*nc*2+jj];
         for(j=1;j<nc;j++){
            GetRNGstate();
            u =  unif_rand();
            PutRNGstate();
            if(u<rr[j]){jj += 1;jj %= 2;}
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
 xchr: 1 if it is x-chromosome
 ---------------------------------------------------*/

void rgeno2(int* data,int nr,int nc,int ninit,int* pedd,double* rr,int xchr){
   double u;
   int i,j;
   int ii,jj;
   if(nr<2){
      error(_("pedigree: at least 2 rows.\n"));
   }else if(nr>INT_MAX){
      error(_("pedigree: too many rows.\n"));
   }
   if(nc<1){
      error(_("recombinaton rate: at least 1 SNP.\n"));
   }else if(nc>(INT_MAX-1)/2){
      error(_("recombinaton rate: too many SNPs.\n"));
   }
   if(xchr){
      for(i=ninit;i<nr;i++){
         R_CheckUserInterrupt();
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
            GetRNGstate();
            u =  unif_rand();
            PutRNGstate();
            if(u<0.5) jj = 0; else jj = 1;
            data[i*nc*2+1] = data[ii*nc*2+jj];
            for(j=1;j<nc;j++){
               GetRNGstate();
               u =  unif_rand();
               PutRNGstate();
               if(u<rr[j]){jj += 1; jj %= 2;}
               data[i*nc*2+j*2+1] = data[ii*nc*2+j*2+jj];
            }
         }
      }
   }else{
      for(i=ninit;i<nr;i++){
         R_CheckUserInterrupt();
         //paternal gamete
         ii = pedd[i*4+1] - 1;
         if(ii>=0){
            GetRNGstate();
            u =  unif_rand();
            PutRNGstate();
            if(u<0.5) jj = 0; else jj = 1;
            data[i*nc*2] = data[ii*nc*2+jj];
            for(j=1;j<nc;j++){
               GetRNGstate();
               u =  unif_rand();
               PutRNGstate();
               if(u<rr[j]){jj += 1; jj %= 2;}
               data[i*nc*2+j*2] = data[ii*nc*2+j*2+jj];
            }
         }
         //maternal gamete
         ii = pedd[i*4+2] - 1;
         if(ii>=0){
            GetRNGstate();
            u =  unif_rand();
            PutRNGstate();
            if(u<0.5) jj = 0; else jj = 1;
            data[i*nc*2+1] = data[ii*nc*2+jj];
            for(j=1;j<nc;j++){
               GetRNGstate();
               u =  unif_rand();
               PutRNGstate();
               if(u<rr[j]){jj += 1; jj %= 2;}
               data[i*nc*2+j*2+1] = data[ii*nc*2+j*2+jj];
            }
         }
      }
   }
}


