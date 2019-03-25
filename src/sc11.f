c
c  nl = 1, nc > 0
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
      subroutine sc11(y, n, x, p, nc, xg, ng, gcv,
     .   coef, tt, pval, v, opt,
     .   pvt, xtx, qr, b, r0, r1, x1, x2, xty, qty,
     .   qraux, work)
      integer n, p, nc, opt
      integer ng
      double precision y(n), x(n,p), xg(n,ng), gcv(n,n),
     .   coef(ng,p+nc+1), tt(ng), pval(ng), v(ng)
c
      integer pvt(p + nc + 1)
      double precision xtx(p,p), qr(p+nc+1, p+nc+1),
     .   b(p+nc+1), r0(n), r1(n), x1(n,p), x2(n,nc+1),
     .   xty(p), qty(p+nc+1), qraux(p+nc+1), work(p+nc+1,2)
c
c      double precision pchisq, pf, normrnd
c      external pchisq, pf, normrnd
      parameter (PI = atan(1.d0) * 4)
      double precision llk0, llk1, llk, lrt, Ftest, df1, df2,
     .   tss, rss, tol, tmp
      integer k0, k1, ny, np
c      character*10 tm
c      call date_and_time(time=tm)
c      write(*,*) tm
c
      tol = 1e-8
      llk = -0.5 * n * (log(2.0 * PI) + 1.0 - log(real(n)))
c      WRITE(*,*) PI, log(real(n))
      k0 = 0
      k1 = 0
      ny = 1
      np = p + nc + 1
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
c     x1 = gcv*x
      do 30 j=1, p
         do 20 i=1, n
            tmp = 0.0
            do 10 m=1, n
               tmp = tmp + gcv(m,i)*x(m,j)
   10       continue
            x1(i,j) = tmp
   20    continue
   30 continue
c
c     xtx = x1'*x1
      do 60 i=1, p
         do 50 j=1, p
            tmp = 0.0
            do 40 m=1, n
               tmp = tmp + x1(m,i)*x1(m,j)
   40       continue
            xtx(i,j) = tmp
            qr(i,j) = tmp
   50    continue
   60 continue
c
c     xty = x1'y
      do 80 i=1, p
         tmp = 0.0
         do 70 m=1, n
            tmp = tmp + x1(m,i)*y(m)
   70    continue
         xty(i) = tmp
         qty(i) = tmp
   80 continue
c
c     null model
      do 85 k=1,np
         pvt(k) = k
   85 continue
      call dqrdc2(qr(1:p,1:p), p, p, p, tol, k0,
     .            qraux(1:p), pvt(1:p), work(1:p,1:2))
      call dqrcf(qr(1:p,1:p), p, k0, qraux(1:p),
     .            qty(1:p), ny, b(1:p), info)
c
      do 100 i=1, n
         tmp = 0.0
         do 90 j=1, k0
            tmp = tmp + x1(i,pvt(j))*b(j)
   90    continue
         r0(i) = y(i) - tmp
  100 continue
c
c     scan the genome xg
      do 390 j=1, ng
         do 120 i=1, p
            do 110 k=1, p
               qr(i,k) = xtx(i,k)
  110       continue
            qty(i) = xty(i)
  120    continue
         do 170 k=1, nc
            do 150 i=1, n
               tmp = 0.0
               do 140 m=1, n
                  tmp = tmp + gcv(m,i)*x(m,p-nc+k)*xg(m,j)
  140          continue
               x2(i,k) = tmp
  150       continue
  170    continue
         do 190 i=1, n
            tmp = 0.0
            do 180 m=1, n
               tmp = tmp + gcv(m,i)*xg(m,j)
  180       continue
            x2(i,nc+1) = tmp
  190    continue
         do 230 k=1, p
            do 220 l=p+1, np
               tmp = 0.0
               do 210 i=1, n
                  tmp = tmp + x1(i,k)*x2(i,l-p)
  210          continue
               qr(k,l) = tmp
               qr(l,k) = tmp
  220       continue
  230    continue
         do 290 k=p+1, np
            do 280 l=k, np
               tmp = 0.0
               do 270 i=1, n
                  tmp = tmp + x2(i,k-p)*x2(i,l-p)
  270          continue
               qr(k,l) = tmp
               qr(l,k) = tmp
  280       continue
  290    continue
         do 310 k=p+1, np
            tmp = 0.0
            do 300 i=1, n
               tmp = tmp + x2(i,k-p)*y(i)
  300       continue
            qty(k) = tmp
  310    continue
c
         do 320 k=1, np
            pvt(k) = k
  320    continue
         call dqrdc2(qr, np, np, np, tol, k1, qraux, pvt, work)
         call  dqrcf(qr, np, k1, qraux, qty, ny, b, info)
c
         do 330 k=1, k1
            coef(j,pvt(k)) = b(k)
  330    continue
         do 350 i=1, n
            tmp = 0.0
            do 340 k=1, k1
               if(pvt(k) > p) then
                  tmp = tmp + x2(i,pvt(k)-p)*b(k)
               else
                  tmp = tmp + x1(i,pvt(k))*b(k)
               endif
  340       continue
            r1(i) = y(i) - tmp
  350    continue
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
               pval(j) = 0.0
            else
               Ftest = (Ftest-rss)/(k1-k0)
               rss = rss/(n-k1)
                  Ftest = Ftest/rss
               df1 = k1 - k0
               df2 = n - k1
               tt(j) = Ftest
               call pff(pval(j),Ftest,df1,df2,0,0)
            endif
c            WRITE(*, *) Ftest, k1-k0, n-k1, pval(j), normrnd()
         endif
c
  390 continue
c      call date_and_time(time=tm)
c      write(*,*) tm
      return
      end

