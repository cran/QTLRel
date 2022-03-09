c
c  nl = 1, nc = 0
c
c  y: response (n)
c  x: covariate matrix (n x p)
c  nc: number of interactive covariates
c  xg: genotype data matrix (n x ng), numeric
c  nl: number of levels to be tested
c  gcv: sqrt(inv(GRM))
c  coef: coefficients (ng x ?)
c  tt: test stat
c  pval: test stats or p-values
c  v: variation explained by putative QTL
c  opt: 1 LRT, 2 p-value of F-test, 3 p-value of LRT (chisq)
c
      subroutine sc10(y, n, x, p, xg, ng, coef, tt, pval, v, opt,
     .   pvt, b, r0, r1, xx, qty, qraux, work)
      integer n, p, opt
      integer ng
      double precision y(n), x(n,p), xg(n,ng),
     .   coef(ng,p+1), tt(ng), pval(ng), v(ng)
c
      integer pvt(p + 1)
      double precision b(p+1), r0(n), r1(n), xx(n,p+1),
     .   qty(n), qraux(p+1), work(p+1,2)
c
c      double precision pchisq, pf, normrnd
c      external pchisq, pf, normrnd
      double precision PI
      parameter (PI = atan(1.d0) * 4)
      double precision llk0, llk1, llk, lrt, Ftest, df1, df2,
     .   tss, rss, tol, tmp
      integer k0, k1, np
c      character*10 tm
c      call date_and_time(time=tm)
c      write(*,*) tm
c
      tol = 1e-8
      llk = -0.5 * n * (log(2.0 * PI) + 1.0 - log(real(n)))
c      WRITE(*,*) PI, log(real(n))
      k0 = 0
      k1 = 0
      np = p + 1
c
c     TSS: total sum of squares
      tmp = 0.0
      do 3 i=1, n
         tmp = tmp + y(i)
   3  continue
      tmp = tmp/n
      tss = 0.0
      do 5 i=1, n
         tss = tss + (y(i)-tmp) ** 2.0
   5  continue
c
c     xx = gcv*x
      do 30 j=1, p
         do 20 i=1, n
            xx(i,j) = x(i,j)
   20    continue
   30 continue
c
c     null model
      do 85 k=1,np
         pvt(k) = k
   85 continue
      call mydqrls(xx,n,p,y,tol,b,r0,qty,k0,pvt,qraux,work)
c
c     scan the genome xg
      do 390 j=1, ng
         do 130 i=1, n
            do 120 l=1, p
               xx(i,l) = x(i,l)
  120       continue
            xx(i,p+1) = xg(i,j)
  130    continue
c
         do 320 k=1, np
            pvt(k) = k
  320    continue
         call mydqrls(xx,n,np,y,tol,b,r1,qty,k1,pvt,qraux,work)
c
         do 330 k=1, k1
            coef(j,pvt(k)) = b(k)
  330    continue
c
         Ftest = 0.0
         rss = 0.0
         do 360 i=1, n
            Ftest = Ftest + r0(i) ** 2.0
            rss = rss + r1(i) ** 2.0
  360    continue
         v(j) = (Ftest - rss)/tss
c
         if (opt .EQ. 1 .OR. opt .EQ. 3) then
            llk0 = 0.0
            llk1 = 0.0
            do 370 i=0, n
               llk0 = llk0 + r0(i) ** 2.0
               llk1 = llk1 + r1(i) ** 2.0
  370       continue
            llk0 = -0.5 * n * log(llk0) + llk
            llk1 = -0.5 * n * log(llk1) + llk
            lrt = 2.0 * (llk1 - llk0)
            tt(j) = lrt
            if (opt .EQ. 3) then
               df1 = k1 - k0
               call pchif(pval(j),lrt,df1,0,0)
            endif
c            WRITE(*, *) llk, llk0, llk1, pval(j)
         else if (opt .EQ. 2) then
            if(k1 .EQ. k0) then
               tt(j) = 0.0
               pval(j) = 1.0
            else
               Ftest = (Ftest-rss)/(k1-k0)
               rss = rss/(n-k1)
                  Ftest = Ftest/rss
               df1 = k1 - k0
               df2 = n - k1
               tt(j) = Ftest
               call pff(pval(j),Ftest,df1,df2,0,0)
            endif
c            WRITE(*, *) "10", Ftest, k1-k0, n-k1, pval(j), normrnd()
         endif
c
  390 continue
c      call date_and_time(time=tm)
c      write(*,*) tm
      return
      end

