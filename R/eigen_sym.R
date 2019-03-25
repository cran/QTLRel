##
## eigen decomposion of real symmetric matrix
##
.eigen.sym<- function(x){
   jobz<- "V"
   ul<- "L"
   n<- dim(x)
   if(length(n) != 2 || n[1] != n[2]){
      stop("x: not a square matrix?", call.=FALSE)
   }else n<- n[1]
   w<- rep(-99999.0, n)
   work<- 0.0
   lwork<- -1
   info<- 0
   ## inquire optimal size of work arrays
   ot<- .Fortran("dsyev",
      JOBZ = as.character(jobz),
      UPLO = as.character(ul),
      N = as.integer(n),
      A = x,
      LDA = as.integer(n),
      W = as.double(w),
      WORK = as.double(work),
      LWORK = as.integer(lwork),
      INFO = as.integer(info),
      PACKAGE = "QTLRel"
   )
   if(ot$INFO != 0)
      stop(paste("Error from Lapack 'dsyev': ", ot$INFO, sep=""))
   lwork = as.integer(ot$WORK);
   ot<- .Fortran("dsyev",
      JOBZ = as.character(jobz),
      UPLO = as.character(ul),
      N = as.integer(n),
      A = x,
      LDA = as.integer(n),
      W = as.double(w),
      WORK = mat.or.vec(lwork,1),
      LWORK = as.integer(lwork),
      INFO = as.integer(info),
      PACKAGE = "QTLRel"
   )
   if(ot$INFO != 0)
      stop(paste("Error from Lapack 'dsyev': ", ot$INFO, sep=""))
   ot<- ot[c("W","Z")]
   names(ot)<- c("values","vectors")
   od<- order(ot$values,decreasing=TRUE)
   ot$values<- ot$values[od]
   ot$vectors<- ot$vectors[,od,drop=FALSE]

   rm(list=setdiff(ls(),"ot")); gc()
   ot
}

##
## adapted from R "Lapack.c"; faster if eigenvectors are desired
##   No faster than R
##
eigen.sym<- function(x){
   jobz<- "V"
   rng<- "A"
   ul<- "L"
   n<- dim(x)
   if(length(n) != 2 || n[1] != n[2]){
      stop("x: not a square matrix?", call.=FALSE)
   }else n<- n[1]
   vl<- 0.0
   vu<- 0.0
   il<- 0
   iu<- 0
   abstol<- 0.0
   w<- rep(-99999.0, n)
   z<- matrix(-99999.0,nrow=n,ncol=n)
   isuppz<- rep(-9,n*2)
   work<- 0.0
   lwork<- -1
   iwork<- 1
   liwork<- -1
   info<- 0
   ## inquire optimal size of work arrays
   ot<- .Fortran("dsyevr",
      JOBZ = as.character(jobz),
      RANGE = as.character(rng),
      UPLO = as.character(ul),
      N = as.integer(n),
      A = x,
      LDA = as.integer(n),
      VL = as.double(vl),
      VU = as.double(vu),
      IL = as.integer(il),
      IU = as.integer(iu),
      ABSTOL = as.double(abstol),
      M = as.integer(n),
      W = as.double(w),
      Z = z,
      LDZ = as.integer(n),
      ISUPPZ = as.integer(isuppz),
      WORK = as.double(work),
      LWORK = as.integer(lwork),
      IWORK = as.integer(iwork),
      LIWORK = as.integer(liwork),
      INFO = as.integer(info),
      PACKAGE = "QTLRel"
   )
   if(ot$INFO != 0)
      stop(paste("Error from Lapack 'dsyevr': ", ot$INFO, sep=""))
   lwork = as.integer(ot$WORK);
   liwork = ot$IWORK
   ot<- .Fortran("dsyevr",
      JOBZ = as.character(jobz),
      RANGE = as.character(rng),
      UPLO = as.character(ul),
      N = as.integer(n),
      A = x,
      LDA = as.integer(n),
      VL = as.double(vl),
      VU = as.double(vu),
      IL = as.integer(il),
      IU = as.integer(iu),
      ABSTOL = as.double(abstol),
      M = as.integer(n),
      W = as.double(w),
      Z = z,
      LDZ = as.integer(n),
      ISUPPZ = as.integer(isuppz),
      WORK = mat.or.vec(lwork,1),
      LWORK = as.integer(lwork),
      IWORK = as.integer(mat.or.vec(liwork,1)),
      LIWORK = as.integer(liwork),
      INFO = as.integer(info),
      PACKAGE = "QTLRel"
   )
   if(ot$INFO != 0)
      stop(paste("Error from Lapack 'dsyevr': ", ot$INFO, sep=""))
   ot<- ot[c("W","Z")]
   names(ot)<- c("values","vectors")
   od<- order(ot$values,decreasing=TRUE)
   ot$values<- ot$values[od]
   ot$vectors<- ot$vectors[,od,drop=FALSE]

   rm(list=setdiff(ls(),"ot")); gc()
   ot
}

