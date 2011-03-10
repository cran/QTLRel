/* 
   calculate generalized identity coefficients:
   Ph(ij), Phi(ijk), Phi(ijkl),Phi(ij,kl)
   using Karigl algorithm
*/

/***************************************************
  compile with option:
     g++  -D_FILE_OFFSET_BITS=64 idcoef.cc
     R CMD SHLIB -D_FILE_OFFSET_BITS=64 idcoef.cc
  the following seem fine with "#define _FILE_OFFSET_BITS 64",
  --which seems crucial for large file writing:
       g++ idcoef.cc
       R CMD SHLIB idcoef.cc
  long long type can hold an integer as large as
  ULLONG_MAX 18446744073709551615ULL
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
#include "string.h"
#include "stdlib.h"
using namespace std;
#define _FILE_OFFSET_BITS 64

typedef int INT;

int fseekerr;
size_t frwsize;
int counter;
double buff;
long long jj;
long long o0[4];
long long o[4];

inline void checkages(INT &a, INT &b);
void kinship(INT** ped, int nr, double** kc);
double phi(INT a, INT b, INT** ped, int* top, FILE** ifs);
double phi(INT a, INT b, INT c, INT** ped, int* top, FILE** ifs);
double phi(INT a, INT b, INT c, INT d, INT** ped, int* top, FILE** ifs);
double phi22(INT a, INT b, INT c, INT d, INT** ped, int* top, FILE** ifs);
void idcoef(INT** ped,INT nr,INT* id,INT nid, int* top, FILE** ifs,FILE** ofs);
void idcoef(INT** ped,INT nr,INT* id,INT nid, int* top, FILE** ifs,double* idcf,int verbose);
void gen_Matrix(double** idcf, int nn, double** ksp, double** DD, double** AD, double** HH, double** MH);

/*************************************************************************
   pedigree: nr by nc array with (id,father,mother,...)
   id: vector of ids for which the identity coefficients are calculated
   outfs[4]: specify out file names
**************************************************************************/

extern "C"{
   void getsize(int& n){
      n=sizeof(long long);
   }

   void kinship(INT* pedigree, INT &nr, INT& nc, double* ksp){
      INT** ped=new INT*[nr]; for(INT i=0;i<nr;i++) ped[i]=pedigree+i*nc;
      double** kc=new double*[nr]; for(INT i=0;i<nr;i++) kc[i]=ksp+i*nr;
      kinship(ped,nr,kc);
   }

   //write to file outfs[4]: phi(a,b), phi(a,b,c), phi(a,b,c,d) and phi22(a,b,c,d)
   void phicw(INT* pedigree,INT& nr,int& nc,INT* id,INT& nid, int* top, char** infs, char** outfs){
         FILE* ifs[4];
         if(top[0]!=-999) for(int i=0; i<4; i++){
            ifs[i] = fopen64(infs[i],"rb+");
            if(!ifs[i]){
               error(_("In_file failed to open.\n"));
            }
         }
         FILE* ofs[4];
         for(int i=0; i<4; i++){
            ofs[i] = fopen64(outfs[i],"wb");
            if(!ofs[i]){
               error(_("Out_file failed to open.\n"));
            }
         }

         INT** ped=new INT*[nr]; for(INT i=0;i<nr;i++) ped[i]=pedigree+i*nc;
         idcoef(ped,nr,id,nid,top,ifs,ofs);

         delete[] ped;
         for(int i=0; i<4; i++)
            fclose(ofs[i]);
         if(top[0]!=-999) for(int i=0; i<4; i++){
            fclose(ifs[i]);
            remove(infs[i]);
         }
   }

   //store in idcf[,9]
   void phicr(INT* pedigree,INT& nr,int& nc,INT* id,INT& nid, int* top, char** infs, double* idcf,int& verbose){
         FILE* ifs[4];
         if(top[0]!=-999) for(int i=0; i<4; i++){
            ifs[i] = fopen64(infs[i],"rb+");
            if(!ifs[i]){
               error(_("In_file failed to open.\n"));
            }
         }

         INT** ped=new INT*[nr]; for(INT i=0;i<nr;i++) ped[i]=pedigree+i*nc;
         idcoef(ped,nr,id,nid,top,ifs,idcf,verbose);

         delete[] ped;
         if(top[0]!=-999) for(int i=0; i<4; i++){
            fclose(ifs[i]);
            remove(infs[i]);
         }
   }

   void gen_Matrix(double* idcf, int& nr, int& nc, int& nn,
      double* ksp, double* DD, double* AD, double* HH, double* MH){
      double** idc=new double*[nr]; for(int i=0;i<nr;i++) idc[i]=idcf+i*nc;
      double** ks=new double*[nn]; for(int i=0;i<nn;i++) ks[i]=ksp+i*nn;
      double** dd=new double*[nn]; for(int i=0;i<nn;i++) dd[i]=DD+i*nn;
      double** ad=new double*[nn]; for(int i=0;i<nn;i++) ad[i]=AD+i*nn;
      double** hh=new double*[nn]; for(int i=0;i<nn;i++) hh[i]=HH+i*nn;
      double** mh=new double*[nn]; for(int i=0;i<nn;i++) mh[i]=MH+i*nn;

      gen_Matrix(idc, nn, ks, dd, ad, hh, mh);

      delete[] idc; delete[] ks; delete[] dd; delete[] ad; delete[] hh; delete[] mh;
   }
}

// write to ofs[k]
void idcoef(INT** ped,INT nr,INT* id,INT nid, int* top, FILE** ifs, FILE** ofs){
   for(INT i=0;i<nid;i++){
      for(INT j=0;j<=i;j++){
         buff = phi(id[i],id[j],ped,top,ifs);
         frwsize = fwrite(&buff,sizeof(double),1,ofs[0]);
         if(frwsize!=1){
            error(_("Data writing errors.\n"));
         }
      }
   }

   for(INT i=0;i<nid;i++){
      for(INT j=0;j<=i;j++){
         for(INT k=0;k<=j;k++){
            buff = phi(id[i],id[j],id[k],ped,top,ifs);
            frwsize = fwrite(&buff,sizeof(double),1,ofs[1]);
            if(frwsize!=1){
               error(_("Data writing errors.\n"));
            }
         }
      }
   }

   for(INT i=0;i<nid;i++){
      for(INT j=0;j<=i;j++){
         for(INT k=0;k<=j;k++){
            for(INT l=0;l<=k;l++){
               buff = phi(id[i],id[j],id[k],id[l],ped,top,ifs);
               frwsize = fwrite(&buff,sizeof(double),1,ofs[2]);
               if(frwsize!=1){
                  error(_("Data writing errors.\n"));
               }
            }
         }
      }
   }

   for(INT i=0;i<nid;i++){
      for(INT j=0;j<=i;j++){
         for(INT k=0;k<=i;k++){
            for(INT l=0;l<=k;l++){
               buff = phi22(id[i],id[j],id[k],id[l],ped,top,ifs);
               frwsize = fwrite(&buff,sizeof(double),1,ofs[3]);
               if(frwsize!=1){
                  error(_("Data writing errors.\n"));
               }
            }
         }
      }
   }
}

// store in idcf[i][j]
void idcoef(INT** ped,INT nr,INT* id,INT nid, int* top, FILE** ifs, double* idcf,int verbose){
   long long ii;
   double aa,bb,ab,aab,abb,aabb,aaxbb,abxab;

   ii = 0;
   for(INT i=0;i<nid;i++){
      if(verbose) cout<<endl<<i+1<<flush;
      for(INT j=0;j<=i;j++){
         if(verbose) cout<<"."<<flush;

         aa = 2.0*phi(id[i],id[i],ped,top,ifs);
         bb = 2.0*phi(id[j],id[j],ped,top,ifs);
         ab = 4.0*phi(id[i],id[j],ped,top,ifs);
         aab = 8.0*phi(id[i],id[i],id[j],ped,top,ifs);
         abb = 8.0*phi(id[i],id[j],id[j],ped,top,ifs);
         aabb = 16.0*phi(id[i],id[i],id[j],id[j],ped,top,ifs); 
         aaxbb = 4.0*phi22(id[i],id[i],id[j],id[j],ped,top,ifs);
         abxab = 16.0*phi22(id[i],id[j],id[i],id[j],ped,top,ifs);

//from eqn (7) in Karigl 81
         idcf[ii] =   0.0 + 0.0*aa + 0.0*bb + 0.25*ab - 0.25*aab - 0.25*abb + 0.25*aabb + 0.0*aaxbb + 0.0*abxab; ii++;
         idcf[ii] =   1.0 - 1.0*aa - 1.0*bb - 0.25*ab + 0.25*aab + 0.25*abb - 0.25*aabb + 1.0*aaxbb + 0.0*abxab; ii++;
         idcf[ii] =   0.0 + 0.0*aa + 0.0*bb - 1.00*ab + 1.00*aab + 0.50*abb - 0.50*aabb + 0.0*aaxbb + 0.0*abxab; ii++;
         idcf[ii] =  -2.0 + 2.0*aa + 1.0*bb + 1.00*ab - 1.00*aab - 0.50*abb + 0.50*aabb - 1.0*aaxbb + 0.0*abxab; ii++;
         idcf[ii] =   0.0 + 0.0*aa + 0.0*bb - 1.00*ab + 0.50*aab + 1.00*abb - 0.50*aabb + 0.0*aaxbb + 0.0*abxab; ii++;
         idcf[ii] =  -2.0 + 1.0*aa + 2.0*bb + 1.00*ab - 0.50*aab - 1.00*abb + 0.50*aabb - 1.0*aaxbb + 0.0*abxab; ii++;
         idcf[ii] =   0.0 + 0.0*aa + 0.0*bb + 0.00*ab + 0.00*aab + 0.00*abb - 0.50*aabb + 0.0*aaxbb + 0.5*abxab; ii++;
         idcf[ii] =   0.0 + 0.0*aa + 0.0*bb + 4.00*ab - 2.00*aab - 2.00*abb + 2.00*aabb + 0.0*aaxbb - 1.0*abxab; ii++;
         idcf[ii] =   4.0 - 2.0*aa - 2.0*bb - 4.00*ab + 2.00*aab + 2.00*abb - 1.50*aabb + 1.0*aaxbb + 0.5*abxab; ii++;
      }
   }
   if(verbose) cout<<endl;
}


void gen_Matrix(double** idcf, int nn, double** ksp, double** DD, double** AD, double** HH, double** MH){
   int ii=0;
   for(int i=0;i<nn;i++){
      for(int j=0;j<=i;j++){
         ksp[i][j] = idcf[ii][0] + (idcf[ii][2] + idcf[ii][4] + idcf[ii][6])/2.0 + idcf[ii][7]/4.0;
            ksp[j][i] = ksp[i][j];
         DD[i][j] = idcf[ii][6];
            DD[j][i] = DD[i][j];
         AD[i][j] = 4*idcf[ii][0] + idcf[ii][2] + idcf[ii][4];
            AD[j][i] = AD[i][j];
         HH[i][j] = idcf[ii][0];
            HH[j][i] = HH[i][j];
         MH[i][j] = idcf[ii][0] + idcf[ii][1];
            MH[j][i] = MH[i][j];
         ii++;
      }
   }
   for(int i=0;i<nn;i++){
      for(int j=0;j<=i;j++){
         MH[i][j] -= (2.0*ksp[i][i] -1.0)*(2.0*ksp[j][j] -1.0);
            MH[j][i] = MH[i][j];
      }
   }
}

/*----------------------------------------
 sort: sort an array x of length n
 returned by arr of length n
 default: increasing=true
 ----------------------------------------*/
template <class T>
void sort(T* x,INT n,T* arr,bool increasing){
   T tmp;

   for(INT i=0;i<n;i++){
      arr[i]=x[i];
   }
   if(increasing==true){
      for(INT i=0;i<n-1;i++){
         for(INT j=i+1;j<n;j++){
            if(arr[j]<arr[i]){
               tmp=arr[i];
               arr[i]=arr[j];
               arr[j]=tmp;
            }
         }
      }
   }else if(increasing==false){
      for(INT i=0;i<n-1;i++){
         for(INT j=i+1;j<n;j++){
            if(arr[j]>arr[i]){
               tmp=arr[i];
               arr[i]=arr[j];
               arr[j]=tmp;
            }
         }
      } 
   }
}

/*----------------------------------------
 sort: sort two pairs x=(ab,cd) => (AB,CD)
 such that AB from ab or cd and A>=B,
 A>=C,C>=D
 ----------------------------------------*/
template <class T>
void sort22(T* x,INT n,T* arr){
   T tmp;

   if(n != 4){
      error(_("n should be 4.\n"));
   }

   for(INT i=0;i<n;i++){
      arr[i]=x[i];
   }
   if(arr[0]<arr[1]){
      tmp=arr[0];
      arr[0]=arr[1];
      arr[1]=tmp;
   }
   if(arr[2]<arr[3]){
      tmp=arr[2];
      arr[2]=arr[3];
      arr[3]=tmp;
   }
   if(arr[0]<arr[2]){
      tmp=arr[0];
      arr[0]=arr[2];
      arr[2]=tmp;
      tmp=arr[1];
      arr[1]=arr[3];
      arr[3]=tmp;
  }
}

/*----------------
 choose(n+k-1,k)
 accuracy won't come without trick!!!
 ----------------*/
long long fn2(long long n){
   long long s;
   s = n*(n+1)/2;
   return s;
}
long long fn3(long long n){
   long long s;
   s = fn2(n)*(n+2)/3;
   return s;
}
long long fn4(long long n){//fn4(1000) = 41917125250, which is 312.3069 Gb
   long long s;
   s = fn3(n)*(n+3)/4;
   return s;
}
// special one
long long fn4_2(long long n){//fn4_2(1000) = 125417041750, which is 934.4298 Gb.
   long long s;
   s = fn3(n)*(3*n+1)/4;
   return s;
}

long long s2(long long* x){
   long long s;
   s = fn2(x[0]-1) + x[1] - 1;
   return s;
}
long long s3(long long* x){
   long long s;
   s = fn3(x[0]-1) + fn2(x[1]-1) + x[2] - 1;
   return s;
}
long long s4(long long* x){
   long long s;
   s = fn4(x[0]-1) + fn3(x[1]-1) + fn2(x[2]-1) + x[3] - 1;
   return s;
}
long long s22(long long* x){
   long long s;
   s = fn4_2(x[0]-1) + (x[1]-1)*fn2(x[0]) + fn2(x[2]-1) + x[3] - 1;
   return s;
}

inline void checkages(INT &a, INT &b)
{
   if (a < b){
      INT tmp = a;
      a = b;
      b = tmp;
   }
}

double phi(INT a, INT b, INT** ped, double** kc)
{
   if( a == 0 || b == 0){
      return 0;
   }
   if(a == b){
//      buff = 0.5 + 0.5*kc[ped[a-1][1]-1][ped[a-1][2]-1];
      buff = 0.5 + 0.5*phi(ped[a-1][1], ped[a-1][2], ped, kc);
   }else if(a < b){
      if(ped[b-1][1]==0 && ped[b-1][2]==0) buff = 0;
      else if(ped[b-1][1]==0) buff = (kc[a-1][ped[b-1][2]-1])/2.0;
      else if(ped[b-1][2]==0) buff = (kc[a-1][ped[b-1][1]-1])/2.0;
      else buff = (kc[a-1][ped[b-1][1]-1] + kc[a-1][ped[b-1][2]-1])/2.0;
   }else{
      if(ped[a-1][1]==0 && ped[a-1][2]==0) buff = 0;
      else if(ped[a-1][1]==0) buff = (kc[ped[a-1][2]-1][b-1])/2.0;
      else if(ped[a-1][2]==0) buff = (kc[ped[a-1][1]-1][b-1])/2.0;
      else buff = (kc[ped[a-1][1]-1][b-1] + kc[ped[a-1][2]-1][b-1])/2.0;
   }

   return buff;
}

void kinship(INT** ped, int nr, double** kc)
{
   for(int i=0; i<nr; i++){
      for(int j=0; j<=i; j++){
         kc[i][j] = phi(i+1, j+1, ped, kc);
            kc[j][i] = kc[i][j];
      }
   }
}

double phi(INT a, INT b, INT** ped, int* top, FILE** ifs)
{
   if( a == 0 || b == 0)
      return 0;
   if(top[0]!=-999) if(top[a-1]==1 && top[b-1]==1){
      o0[0]=a; o0[1]=b;
      sort(o0,2,o,false);
      jj = s2(o);
      fseekerr = fseeko64(ifs[0],jj*sizeof(double),SEEK_SET);
/*
      counter = 3;
      while(fseekerr && counter>0){
         fseekerr = fseeko64(ifs[0],jj*sizeof(double),SEEK_SET);
         counter--;
      }
      if(fseekerr){
         error(_("Seeking errors (1) occurred repeatedly.\n"));
      }
*/
      frwsize = fread(&buff,sizeof(double),1,ifs[0]);
/*
      counter = 3;
      while(frwsize!=1 && counter>0){
         frwsize = fread(&buff,sizeof(double),1,ifs[0]);
         counter--;
      }
      if(frwsize!=1){
         error(_("Reading errors (1) occurred repeatedly.\n"));
      }
*/

      return(buff);
   }
   if(a == b)
      return (0.5 + 0.5*phi(ped[a-1][1], ped[a-1][2], ped, top, ifs));
   else{
      if(a < b){
         INT tmp = a;
         a = b;
         b = tmp;
      }
      return((phi(ped[a-1][1],b, ped, top, ifs) + phi(ped[a-1][2],b, ped, top, ifs))/2.0);
   }
}

double phi(INT a, INT b, INT c, INT** ped, int* top, FILE** ifs)
{
   if (a == 0 || b == 0 || c == 0) // case 0
      return 0.0;
   if(top[0]!=-999) if(top[a-1]==1 && top[b-1]==1 && top[c-1]==1){
      o0[0]=a; o0[1]=b; o0[2]=c;
      sort(o0,3,o,false);
      jj = s3(o);

      fseekerr = fseeko64(ifs[1],jj*sizeof(double),SEEK_SET);
/*
      counter = 3;
      while(fseekerr && counter>0){
         cout<<"."<<flush;
         fseekerr = fseeko64(ifs[1],jj*sizeof(double),SEEK_SET);
         counter--;
      }
      if(fseekerr){
         error(_("Seeking errors (2) occurred repeatedly.\n"));
      }
*/
      frwsize = fread(&buff,sizeof(double),1,ifs[1]);
/*
      counter = 3;
      while(frwsize!=1 && counter>0){
         cout<<"."<<flush;
         frwsize = fread(&buff,sizeof(double),1,ifs[1]);
         counter--;
      }
      if(frwsize!=1){
         error(_("Reading errors (2) occurred repeatedly.\n"));
      }
*/

      return(buff);
   }

   if (a == b && a == c) // case 1
      return ((1.0 + 3.0 * phi(ped[a-1][1], ped[a-1][2], ped, top, ifs)) / 4.0);

   checkages(a,c);
   checkages(b,c);
   if(a == b) // case 2
      return ((phi(a, c, ped, top, ifs) + phi(ped[a-1][1], ped[a-1][2], c, ped, top, ifs)) / 2.0);

   checkages(a,b);
   return ((phi(ped[a-1][1], b, c, ped, top, ifs) + phi(ped[a-1][2], b, c, ped, top, ifs)) / 2.0);
}


double phi(INT a, INT b, INT c, INT d, INT** ped, int* top, FILE** ifs)
{
   if (a == 0 || b == 0 || c == 0 || d == 0)
      return 0.0;
   if(top[0]!=-999) if(top[a-1]==1 && top[b-1]==1 && top[c-1]==1 && top[d-1]==1){
      o0[0]=a; o0[1]=b; o0[2]=c; o0[3]=d;
      sort(o0,4,o,false);
      jj = s4(o);

      fseekerr = fseeko64(ifs[2],jj*sizeof(double),SEEK_SET);
/*
      counter = 3;
      while(fseekerr && counter>0){
         fseekerr = fseeko64(ifs[2],jj*sizeof(double),SEEK_SET);
         counter--;
      }
      if(fseekerr){
         error(_("Seeking errors (3) occurred repeatedly.\n"));
      }
*/
      frwsize = fread(&buff,sizeof(double),1,ifs[2]);
/*
      counter = 3;
      while(frwsize!=1 && counter>0){
         frwsize = fread(&buff,sizeof(double),1,ifs[2]);
         counter--;
      }
      if(frwsize!=1){
         error(_("Reading errors (3) occurred repeatedly.\n"));
      }
*/

      return(buff);
   }

   if (a == b && a == c && a == d)   // case 1
      return ((1.0 + 7.0 * phi(ped[a-1][1], ped[a-1][2], ped, top, ifs)) / 8.0);

   checkages(a,d);
   checkages(b,d);
   checkages(c,d);
   if(a == b && b == c)   // case 2
      return ((phi(a, d, ped, top, ifs) + 3.0 * phi(ped[a-1][1], ped[a-1][2], d, ped, top, ifs)) / 4.0);

   checkages(a,c);
   checkages(b,c);
   if(a == b)   // case 3
      return ((phi(a, c, d, ped, top, ifs) + phi(ped[a-1][1], ped[a-1][2], c, d, ped, top, ifs)) / 2.0);

   checkages(a,b);
   return ((phi(ped[a-1][1], b, c, d, ped, top, ifs) + phi(ped[a-1][2], b, c, d, ped, top, ifs)) / 2.0);
}

double phi22(INT a, INT b, INT c, INT d, INT** ped, int* top, FILE** ifs)
{
   if( a == 0 || b == 0 || c == 0 || d == 0 )
      return 0.0;
   if(top[0]!=-999) if(top[a-1]==1 && top[b-1]==1 && top[c-1]==1 && top[d-1]==1){
      o0[0]=a; o0[1]=b; o0[2]=c; o0[3]=d;
      sort22(o0,4,o);
      jj = s22(o);

      fseekerr = fseeko64(ifs[3],jj*sizeof(double),SEEK_SET);
/*
      counter = 3;
      while(fseekerr && counter>0){
         fseekerr = fseeko64(ifs[3],jj*sizeof(double),SEEK_SET);
         counter--;
      }
      if(fseekerr){
         error(_("Seeking errors (4) occurred repeatedly.\n"));
      }
*/
      frwsize = fread(&buff,sizeof(double),1,ifs[3]);
/*
      counter = 3;
      while(frwsize!=1 && counter>0){
         frwsize = fread(&buff,sizeof(double),1,ifs[3]);
         counter--;
      }
      if(frwsize!=1){
         error(_("Reading errors (4) occurred repeatedly.\n"));
      }
*/

      return(buff);
   }

   if( a == b && a == c && a == d ) // case 1 
      return ((1.0 + 3.0 * phi(ped[a-1][1], ped[a-1][2], ped, top, ifs)) / 4.0);

   checkages(a,b);
   checkages(c,d);
   if( a == c)
      checkages(b,d);
   if( a == b && a == c ) // case 2
      return ((phi(a, d, ped, top, ifs) + phi(ped[a-1][1], ped[a-1][2], d, ped, top, ifs)) / 2.0);

   if( a < c ){
      INT tmp = a; a = c; c = tmp;
      tmp = b; b = d; d = tmp;
   }
   if( a == b ) // case 3
      return ((phi(c, d, ped, top, ifs) + phi22(ped[a-1][1], ped[a-1][2], c, d, ped, top, ifs)) / 2.0);

   if( a == c ) // case 4
      return ((2.0 * phi(a, b, d, ped, top, ifs) + phi22(ped[a-1][1], b, ped[a-1][2], d, ped, top, ifs) +
            phi22(ped[a-1][2], b, ped[a-1][1], d, ped, top, ifs)) / 4.0);

   return ((phi22(ped[a-1][1],b,c,d,ped, top, ifs) + phi22(ped[a-1][2],b,c,d,ped, top, ifs)) / 2.0);
}


