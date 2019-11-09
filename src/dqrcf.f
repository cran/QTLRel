c     given qr of x, solve x*b = y
c
c     on entry
c
c        x      double precision(ldx,p).
c               x contains the output of dqrdc.
c
c        n      integer.
c               n is the number of rows of the matrix xk.  it must
c               have the same value as n in dqrdc.
c
c        k      integer.
c               k is the number of columns of the matrix xk.  k
c               must nnot be greater than min(n,p), where p is the
c               same as in the calling sequence to dqrdc.
c
c        qraux  double precision(p).
c               qraux contains the auxiliary output from dqrdc.
c
c        y      double precision(n)
c               y contains an n-vector that is to be manipulated
c               by dqrsl.
c
c        ny     integer.
c               n is the number of columns of the matrix y
c
c     on return
c
c        b      double precision(k)
c               b contains the solution of the least squares problem
c
c                    minimize norm2(y - xk*b),
c
c        info   integer.
c               info is zero unless the computation of b has
c               been requested and r is exactly singular.  in
c               this case, info is the index of the first zero
c               diagonal element of r and b is left unaltered.
c
      subroutine dqrcf(x, n, k, qraux, y, ny, b, info)

      integer n, k, ny, info
      double precision x(n,k), qraux(k), y(n,ny), b(k,ny)
      integer j
      double precision dummy(1)
      do 10 j = 1,ny
          call dqrsl(x, n, n, k, qraux, y(1,j), dummy,
     .               y(1,j), b(1,j), dummy, dummy, 100, info)
   10 continue
      return
      end

