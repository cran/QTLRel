
#include "xxx.h"

/*------------------------
 Kolmogorov distribution
 -------------------------*/
extern "C"{
   void kolm(double* x, int&n){
      double xx, s, nmax2;
      unsigned long nmax;
      long int jj;

      for(int j=0; j < n; j++){
         if(x[j] <= 0) s = 0.0;
         else{
            xx = -2*x[j]*x[j];
            nmax = static_cast<unsigned long> (10.0/x[j] + 1);
            nmax2 = pow(2.0, static_cast<double> (sizeof(unsigned long)*8-1));
            if(nmax > nmax2) nmax = static_cast<unsigned long>(nmax2);

//nmax = 500;
            s = 1.0;
            jj = 2;
            for(unsigned long i=1; i <= nmax; i++, jj *= -1){
               s -= exp(xx*i*i)*jj;
            }
         }
         if(s < -1e-8) cout<<"Kolmogorov: negative..."<<endl;
         else if(s<0) s = 0.0;
         x[j] = s;
      }
   }

   void Fn(double* t, int& nt, double* x, int& nx){
      double s;
      for(int i=0; i<nt; i++){
         s = 0.0;
         for(int j=0; j<nx; j++){
            if(x[j] <= t[i]) s++;
         }
         t[i] = s/nx;
      }
   }

   void qFn(double* t, int& nt, double* x, int& nx){
      int jj;
      for(int i=0; i<nt; i++){
         if(t[i] <= 0.0) t[i] = -1e+300;
         else if(t[i] >= 1.0) t[i] = 1e+300;
         else{
            for(int j=0; j<nx; j++){
               if((j+1.0)/nx >= t[i]){jj = j; break;}
            }
            t[i] = x[jj];
         }
      }
   }
}


