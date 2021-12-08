
W.inv<- function(W, symmetric=TRUE,inverse=TRUE){
   eW <- eigen(W, symmetric=symmetric)
   d <- eW$values
   if (min(d) <0  && abs(min(d))>sqrt(.Machine$double.eps))
       stop("'W' is not positive definite.", call.=FALSE)
   else d[d<=0]<- ifelse(inverse, Inf, 0)
   A <- sweep(t(eW$vector), 1, d^ifelse(inverse, -0.5, 0.5), "*")
   A # t(A)%*%A = W^{-1}
}

# adapted from lm.gls in MASS
lmGls<- function (formula, data, A, contrasts = NULL, ...) {
    mf<- match.call(expand.dots = FALSE)
        m<- match(c("formula", "data"), names(mf), 0L)
        mf<- mf[c(1L, m)]
        mf$drop.unused.levels<- TRUE
        mf[[1L]]<- quote(stats::model.frame)
        mf<- eval(mf, parent.frame())
    Terms <- attr(mf, "terms")
    yy <- model.response(mf, "numeric")
       y<- A%*%yy
    xx <- model.matrix(Terms, mf, contrasts)
       x<- A%*%xx
    dtf<- data.frame(y=y,x)

    fit<- lm(y~.-1, data=dtf, ...)
    names(fit$coefficients)[1]<- "(Intercept)"

    fit
}

# extract model matrix
mdlMtr<- function (formula, data, contrasts = NULL, ...) {
    mf<- match.call(expand.dots = FALSE)
        m<- match(c("formula", "data"), names(mf), 0L)
        mf<- mf[c(1L, m)]
        mf$drop.unused.levels<- TRUE
        mf[[1L]]<- quote(stats::model.frame)
        mf<- eval(mf, parent.frame())
    Terms <- attr(mf, "terms")
    xx <- model.matrix(Terms, mf, contrasts)

    xx
}

# generalized least squares test
scanOne.0 <-
   function(y,
            x,
            prdat,
            hinvCov,
            intc = NULL,
            test = c("None","F","LRT"))
{
# prdat$pr: n by ? by ? matrix, allele probabilities
   test<- match.arg(test)

   snp<- prdat$snp
   chr<- prdat$chr
   dist<- prdat$dist
   opt<- switch(test,
       "None" = 1,
          "F" = 2,
      "LRT" = 3
   )

   if(is.null(intc)){
      if(!missing(x)){
         oTmp<- data.frame(y=y,x)
      }else{
         oTmp<- data.frame(y=y)
      }
      xx<- mdlMtr(y~., data=oTmp)
      xg<- as.matrix(prdat$pr[,1,] - prdat$pr[,3,])

      cns2<- "g"
      str<- colnames(xx)
         str<- c(str,cns2)

      rm(prdat, oTmp, cns2); gc()

      y<- as.vector(y)
      n<- nrow(xx)
      p<- ncol(xx)
      ng<- dim(xg)[2]
      ot<- .Fortran("sc10",
         y = hinvCov%*%y,
         n = as.integer(n),
         x = hinvCov%*%xx,
         p = as.integer(p),
         xg = hinvCov%*%xg,
         ng = as.integer(ng),
         coefficients = mat.or.vec(ng,p+1),
         tt = double(ng),
         pval = double(ng),
         v = double(ng),
         opt = as.integer(opt),
         pvt = 1L:(p+1),
         b = double(p+1),
         r0 = double(n),
         r1 = double(n),
         xx = mat.or.vec(n,p+1),
         qty = double(n),
         qraux = double(p+1),
         work = mat.or.vec(p+1,2),
         PACKAGE = "QTLRel"
      )[c("tt","pval","v", "coefficients")]
   }else{
      if(!missing(x)){
         oTmp<- data.frame(y=y,x,intc)
      }else{
         oTmp<- data.frame(y=y,intc)
      }
      xx<- mdlMtr(y~., data=oTmp)
      xg<- as.matrix(prdat$pr[,1,] - prdat$pr[,3,])

      oTmp<- data.frame(y=y,intc)
      cns1<- colnames(mdlMtr(y~., data=oTmp))
      cns2<- "g"
      str<- colnames(xx)
      for(s1 in cns1[-1]){
         str<- c(str, paste(s1,cns2,sep=":"))
      }
      str<- c(str,cns2)
      nc<- length(cns1)-1

      rm(x, prdat, oTmp, cns1, cns2, s1); gc()

      y<- as.vector(y)
      n<- nrow(xx)
      p<- ncol(xx)
      ng<- dim(xg)[2]
      ot<- .Fortran("sc11",
         y = hinvCov%*%y,
         n = as.integer(n),
         x = xx,
         p = as.integer(p),
         nc = as.integer(nc),
         xg = xg,
         ng = as.integer(ng),
         gcv = t(hinvCov),
         coefficients = mat.or.vec(ng,p+nc+1),
         tt = double(ng),
         pval = double(ng),
         v = double(ng),
         opt = as.integer(opt),
         pvt = 1L:(p+nc+1),
         b = double(p+nc+1),
         r0 = double(n),
         r1 = double(n),
         xx = mat.or.vec(n,p+nc+1),
         qty = double(n),
         qraux = double(p+nc+1),
         work = mat.or.vec(p+nc+1,2),
         PACKAGE = "QTLRel"
      )[c("tt","pval","v", "coefficients")]
   }
   names(ot$tt)<- 
      names(ot$pval)<- 
      names(ot$v)<- 
      rownames(ot$coefficients)<- snp
   nms<- c("LRT","pval","v","parameters")
   if(test == "F") nms[1]<- "F"
   names(ot)<- nms
   colnames(ot$parameters)<- str
   if(test == "None") ot[2]<- NULL

   ot
}

scanOne.1 <-
   function(y,
            x,
            prdat,
            hinvCov,
            intc = NULL,
            test = c("None","F","LRT"))
{
# prdat$pr: n by 3 by ? matrix, conditional probabilities
   test<- match.arg(test)

   snp<- prdat$snp
   chr<- prdat$chr
   dist<- prdat$dist
   opt<- switch(test,
       "None" = 1,
          "F" = 2,
        "LRT" = 3
   )

   if(is.null(intc)){
      if(!missing(x)){
         oTmp<- data.frame(y=y,x)
      }else{
         oTmp<- data.frame(y=y)
      }
      xx<- mdlMtr(y~., data=oTmp)
      xg<- prdat$pr
         xg<- aperm(xg, perm=c(1,3,2))
         xg[,,1]<- prdat$pr[,1,]-prdat$pr[,3,]
         xg[,,2]<- prdat$pr[,2,]
         xg<- array(xg[,,1:2],dim=c(dim(prdat$pr)[c(1,3)],2))
      for(k in 1:dim(xg)[3])
         xg[,,k]<- hinvCov%*%xg[,,k]

      cns2<- c("a", "d")
      str<- colnames(xx)
         str<- c(str,cns2)

      rm(k, prdat, oTmp, cns2); gc()

      y<- as.vector(y)
      n<- nrow(xx)
      p<- ncol(xx)
      ng<- dim(xg)[2]
      nl<- dim(xg)[3]
      ot<- .Fortran("sc20",
         y = hinvCov%*%y,
         n = as.integer(n),
         x = hinvCov%*%xx,
         p = as.integer(p),
         xg = xg,
         ng = as.integer(ng),
         nl = as.integer(nl),
         coefficients = mat.or.vec(ng,p+nl),
         tt = double(ng),
         pval = double(ng),
         v = double(ng),
         opt = as.integer(opt),
         pvt = 1L:(p+nl),
         b = double(p+nl),
         r0 = double(n),
         r1 = double(n),
         xx = mat.or.vec(n,p+nl),
         qty = double(n),
         qraux = double(p+nl),
         work = mat.or.vec(p+nl,2),
         PACKAGE = "QTLRel"
      )[c("tt","pval","v", "coefficients")]
   }else{
      if(!missing(x)){
         oTmp<- data.frame(y=y,x,intc)
      }else{
         oTmp<- data.frame(y=y,intc)
      }
      xx<- mdlMtr(y~., data=oTmp)
      xg<- prdat$pr
         xg<- aperm(xg, perm=c(1,3,2))
         xg[,,1]<- prdat$pr[,1,]-prdat$pr[,3,]
         xg[,,2]<- prdat$pr[,2,]
         xg<- array(xg[,,1:2],dim=c(dim(prdat$pr)[c(1,3)],2))

      oTmp<- data.frame(y=y,intc)
      cns1<- colnames(mdlMtr(y~., data=oTmp))
      cns2<- c("a", "d")
      str<- colnames(xx)
      for(s1 in cns1[-1]){
         for(s2 in cns2){
            str<- c(str, paste(s1,s2,sep=":"))
         }
      }
      str<- c(str,cns2)
      nc<- length(cns1)-1

      rm(prdat, oTmp, cns1, cns2, s1, s2); gc()

      y<- as.vector(y)
      n<- nrow(xx)
      p<- ncol(xx)
      ng<- dim(xg)[2]
      nl<- dim(xg)[3]
      ot<- .Fortran("sc21",
         y = hinvCov%*%y,
         n = as.integer(n),
         x = xx,
         p = as.integer(p),
         nc = as.integer(nc),
         xg = xg,
         ng = as.integer(ng),
         nl = as.integer(nl),
         gcv = t(hinvCov),
         coefficients = mat.or.vec(ng,p+nc*nl+nl),
         tt = double(ng),
         pval = double(ng),
         v = double(ng),
         opt = as.integer(opt),
         pvt = 1L:(p+nc*nl+nl),
         b = double(p+nc*nl+nl),
         r0 = double(n),
         r1 = double(n),
         xx = mat.or.vec(n,p+nc*nl+nl),
         qty = double(n),
         qraux = double(p+nc*nl+nl),
         work = mat.or.vec(p+nc*nl+nl,2),
         PACKAGE = "QTLRel"
      )[c("tt","pval","v", "coefficients")]
   }
   names(ot$tt)<- 
      names(ot$pval)<- 
      names(ot$v)<- 
      rownames(ot$coefficients)<- snp
   nms<- c("LRT","pval","v","parameters")
   if(test == "F") nms[1]<- "F"
   names(ot)<- nms
   colnames(ot$parameters)<- str
   if(test == "None") ot[2]<- NULL
   ot$chr<- chr
   ot$dist<- dist

   ot
}

scanOne.2 <-
   function(y,
            x,
            gdat,
            hinvCov,
            intc = NULL,
            numGeno = FALSE,
            test = c("None","F","LRT"))
{
# gdat: n by ? matrix, marker data. Markers in columes!!!
# intc: covariates that interact with QTL
   test<- match.arg(test)
   opt<- switch(test,
       "None" = 1,
          "F" = 2,
        "LRT" = 3
   )
 
  if(is.null(intc)){
      if(!missing(x)){
         oTmp<- data.frame(y=y,x)
      }else{
         oTmp<- data.frame(y=y)
      }
      xx<- mdlMtr(y~., data=oTmp)
      if(numGeno){
         xg<- as.matrix(gdat)
            xg<- hinvCov%*%xg
         nl<- 1
      }else{
         fct<- function(x){
            nlevels(as.factor(x))
         }
         nl<- sapply(as.data.frame(gdat),fct)
         tbl<- table(nl)
         if(length(tbl) != 1){
            print(tbl)
            stop("Levels across scanning loci NOT identical!",call.=FALSE)
         }else
            nl<- max(nl)
         xg<- array(NA, dim=c(nrow(gdat),ncol(gdat),nl))
         for(j in 1:ncol(gdat)){
            xt<- as.factor(gdat[,j])
            idx<- !is.na(xt)
            xt<- mdlMtr(~xt)
            xg[idx,j,1:ncol(xt)]<- xt
         }
         for(k in 1:nl)
            xg[,,k]<- hinvCov%*%xg[,,k]
         if(any(dim(xg)<2)) xg<- array(xg[,,-1], dim=c(dim(gdat),nl-1)) else
            xg<- xg[,,-1]
         nl<- nl-1

         rm(fct,j,k,xt,idx)
      }

      if(nl == 1) cns2<- "g" else
         cns2<- paste("g",1:nl+1, sep="")
      str<- colnames(xx)
         str<- c(str,cns2)

      cns<- colnames(gdat)

      rm(gdat, oTmp, cns2); gc()

      y<- as.vector(y)
      n<- nrow(xx)
      p<- ncol(xx)
      ng<- dim(xg)[2]
      #save(y,n,p,ng,xx,xg,hinvCov,file="tt.RData"); q("no")
      if(nl < 1){
         stop("'gdat': something wrong?")
      }else if(nl == 1) ot<- .Fortran("sc10",
         y = hinvCov%*%y,
         n = as.integer(n),
         x = hinvCov%*%xx,
         p = as.integer(p),
         xg = xg,
         ng = as.integer(ng),
         coefficients = mat.or.vec(ng,p+1),
         tt = double(ng),
         pval = double(ng),
         v = double(ng),
         opt = as.integer(opt),
         pvt = 1L:(p+1),
         b = double(p+1),
         r0 = double(n),
         r1 = double(n),
         xx = mat.or.vec(n,p+1),
         qty = double(n),
         qraux = double(p+1),
         work = mat.or.vec(p+1,2),
         PACKAGE = "QTLRel"
      )[c("tt","pval","v", "coefficients")] else ot<- .Fortran("sc20",
         y = hinvCov%*%y,
         n = as.integer(n),
         x = hinvCov%*%xx,
         p = as.integer(p),
         xg = xg,
         ng = as.integer(ng),
         nl = as.integer(nl),
         coefficients = mat.or.vec(ng,p+nl),
         tt = double(ng),
         pval = double(ng),
         v = double(ng),
         opt = as.integer(opt),
         pvt = 1L:(p+nl),
         b = double(p+nl),
         r0 = double(n),
         r1 = double(n),
         xx =  mat.or.vec(n,p+nl),
         qty = double(n),
         qraux = double(p+nl),
         work = mat.or.vec(p+nl,2),
         PACKAGE = "QTLRel"
      )[c("tt","pval","v", "coefficients")]
   }else{
      if(!missing(x)){
         oTmp<- data.frame(y=y,x,intc)
      }else{
         oTmp<- data.frame(y=y,intc)
      }
      xx<- mdlMtr(y~., data=oTmp)
      if(numGeno){
         xg<- as.matrix(gdat)
         nl<- 1
      }else{
         fct<- function(x){
            nlevels(as.factor(x))
         }
         nl<- sapply(as.data.frame(gdat),fct)
            nl<- max(nl)
         xg<- array(NA, dim=c(nrow(gdat),ncol(gdat),nl))
         for(j in 1:ncol(gdat)){
            xt<- as.factor(gdat[,j])
            idx<- !is.na(xt)
            xt<- mdlMtr(~xt)
            xg[idx,j,1:ncol(xt)]<- xt
         }
         if(any(dim(xg)<2)) xg<- array(xg[,,-1], dim=c(dim(gdat),nl-1)) else
            xg<- xg[,,-1]
         nl<- nl-1

         rm(fct, j, xt, idx)
      }

      oTmp<- data.frame(y=y,intc)
      cns1<- colnames(mdlMtr(y~., data=oTmp))
      if(nl == 1) cns2<- "g" else
         cns2<- paste("g",1:nl+1, sep="")
      str<- colnames(xx)
      for(s1 in cns1[-1]){
         for(s2 in cns2){
            str<- c(str, paste(s1,s2,sep=":"))
         }
      }
      str<- c(str,cns2)
      nc<- length(cns1)-1

      cns<- colnames(gdat)

      rm(gdat, oTmp, cns1, cns2, s1, s2); gc()

      y<- as.vector(y)
      n<- nrow(xx)
      p<- ncol(xx)
      ng<- dim(xg)[2]
      if(nl < 1){
         stop("'gdat': something wrong?")
      }else if(nl ==1) ot<- .Fortran("sc11",
         y = hinvCov%*%y,
         n = as.integer(n),
         x = xx,
         p = as.integer(p),
         nc = as.integer(nc),
         xg = xg,
         ng = as.integer(ng),
         gcv = t(hinvCov),
         coefficients = mat.or.vec(ng,p+nc+1),
         tt = double(ng),
         pval = double(ng),
         v = double(ng),
         opt = as.integer(opt),
         pvt = 1L:(p+nc+1),
         b = double(p+nc+1),
         r0 = double(n),
         r1 = double(n),
         xx = mat.or.vec(n,p+nc+1),
         qty = double(n),
         qraux = double(p+nc+1),
         work = mat.or.vec(p+nc+1,2),
         PACKAGE = "QTLRel"
      )[c("tt","pval","v", "coefficients")] else ot<- .Fortran("sc21",
         y = hinvCov%*%y,
         n = as.integer(n),
         x = xx,
         p = as.integer(p),
         nc = as.integer(nc),
         xg = xg,
         ng = as.integer(ng),
         nl = as.integer(nl),
         gcv = t(hinvCov),
         coefficients = mat.or.vec(ng,p+nc*nl+nl),
         tt = double(ng),
         pval = double(ng),
         v = double(ng),
         opt = as.integer(opt),
         pvt = 1L:(p+nc*nl+nl),
         b = double(p+nc*nl+nl),
         r0 = double(n),
         r1 = double(n),
         xx =  mat.or.vec(n,p+nc*nl+nl),
         qty = double(n),
         qraux = double(p+nc*nl+nl),
         work = mat.or.vec(p+nc*nl+nl,2),
         PACKAGE = "QTLRel"
      )[c("tt","pval","v", "coefficients")]
   }
   names(ot$tt)<- 
      names(ot$pval)<- 
      names(ot$v)<- 
      rownames(ot$coefficients)<- cns
   nms<- c("LRT","pval","v","parameters")
   if(test == "F") nms[1]<- "F"
   names(ot)<- nms
   colnames(ot$parameters)<- str
   if(test == "None") ot[2]<- NULL

   ot
}

scanOne<- 
   function(y,
            x,
            gdat,
            prdat = NULL,
            vc = NULL,
            intc = NULL,
            numGeno = FALSE,
            test = c("None","F","LRT"),
            minorGenoFreq = 0,
            rmv = TRUE)

{
   if(!all(is.finite(y)))
      stop("y: non-numeric or infinite data points not allowed.", call.=FALSE)
   if(!missing(x))
      if(any(sapply(x,is.infinite) | sapply(x,is.na)))
         stop("x: missing or infinite data points not allowed.", call.=FALSE)
   UseMethod("scanOne")
}

scanOne.default<- 
   function(y,
            x,
            gdat,
            prdat = NULL,
            vc = NULL,
            intc = NULL,
            numGeno = FALSE,
            test = c("None","F","LRT"),
            minorGenoFreq = 0,
            rmv = TRUE)
{
   if(!is.null(vc)){
      if(is.element("VC",attr(vc,"class"))){
         if(is.null(vc$hinvCov)){
            nv<- length(vc$v)
            nb<- length(vc$par) - nv
            nr<- nrow(vc$y)
            cov<- matrix(0,nrow=nr,ncol=nr)
            for(i in 1:nv)
               cov<- cov + vc$v[[i]]*vc$par[nb+i]
            hinvCov<- W.inv(cov)
         }else hinvCov<- vc$hinvCov
      }else{
         if(is.data.frame(vc)) vc<- as.matrix(vc)
         if(!is.matrix(vc)) stop("'vc' should be a matrix.", call.=FALSE)
         if(!is.numeric(vc)) stop("'vc' should be a numeric matrix.", call.=FALSE)
         hinvCov<- W.inv(vc)
      }
   }else hinvCov<- diag(nrow(as.matrix(y)))
   if(!is.null(prdat)){
      if(is.element("addEff",class(prdat))){
         pv<- scanOne.0(y=y,x=x,prdat=prdat,hinvCov=hinvCov,intc=intc,test=test)
      }else{
         pv<- scanOne.1(y=y,x=x,prdat=prdat,hinvCov=hinvCov,intc=intc,test=test)
      }
   }else{
      gdat<- as.matrix(gdat)
      if(any(is.na(gdat)))
         stop("There are missing genotypes...", call.=FALSE)
      if(numGeno){
         idx<- FALSE
      }else{
         tb<- sort(union(gdat,NULL))
         tbf<- NULL
         for(ii in tb) tbf<- rbind(tbf,colSums(gdat==ii))
            if(sum(tbf)!=nrow(gdat)*ncol(gdat)) stop("Error occurred.", call.=FALSE)
         tbf<- apply(tbf,2,min)
         idx<- (tbf < nrow(gdat)*minorGenoFreq)
      }
      if(sum(idx)>0){
         if(rmv){
            tb<- sort(union(gdat,NULL))
            tbf<- NULL
            for(ii in tb) tbf<- rbind(tbf,colSums(gdat==ii))
               if(sum(tbf)!=nrow(gdat)*ncol(gdat)) stop("Error occurred.", call.=FALSE)
            tbf<- apply(tbf,2,min)
            idx<- (tbf < nrow(gdat)*minorGenoFreq)
            gdat<- as.matrix(gdat[,!idx])
            rm(tb,tbf,ii,idx)
         }else{
            cat("   Minor genotype frequency is too small at one or more SNPs.\a\n")
            return(NULL)
        }
      }

      if(!numGeno)
         gdat<- as.data.frame(gdat)
      pv<- scanOne.2(y=y,x=x,gdat=gdat,hinvCov=hinvCov,intc=intc,numGeno=numGeno,test=test)
   }

   class(pv)<- c("scanOne",match.arg(test))
   pv
}

print.scanOne <-
   function(x,...)
{
   tt<- x; class(tt)<- NULL
      tt$parameters<- NULL
      tt<- as.data.frame(tt)
   if(length(tt$LRT)>0){
      idx<- 1:min(5,length(tt$LRT))

      cat("Test statistic (LRT):\n")
      print(tt$LRT[idx])
      cat("... ...\n\n")

      cat("Variance explained:\n")
      print(tt$v[idx])
      cat("... ...\n\n")

      cat("Coefficients:\n")
      print(x$par[idx,])
      cat("... ...\n\n")
   }else if(length(tt$pval)>0){
      idx<- 1:min(5,length(tt$pval))

      cat("Test statistic (P-value):\n")
      print(tt$pval[idx])
      cat("... ...\n\n")

      cat("Variance explained:\n")
      print(tt$v[idx])
      cat("... ...\n\n")

      cat("Coefficients:\n")
      print(x$par[idx,])
      cat("... ...\n\n")
   }else{
      cat("Test statistic:\n")
      print(tt)
      cat("\n")

      cat("Coefficients:\n")
      print(x$par)
   }
}

scanTwo.1 <-
   function(y,
            x,
            prdat,
            cov)
{
   diag.cov<- diag(as.matrix(cov))
   if( max( abs( cov-diag(diag.cov) ) ) < min(1e-5,1e-5*max(diag.cov)) ){
      if( max(diag.cov-min(diag.cov)) < min(1e-5,1e-5*max(diag.cov)) ){
         weights<- NULL
      }else weights<- 1/diag.cov
   }else weights<- NA
   gcv<- W.inv(cov)
   nsnp<- dim(prdat$pr)[3]
   P<- matrix(NA,nrow=nsnp,ncol=nsnp)
      rownames(P)<- colnames(P)<- prdat$snp
   if(nsnp<1) return(NULL)

   if(!missing(x)){
      oTmp.xy<- data.frame(y=y,x)
   }else{
      oTmp.xy<- data.frame(y=y)
   }
   for(i in 1:(nsnp-1)){
      for(k in (i+1):nsnp){
         xTmp<- data.frame(a1=prdat$pr[,1,i] - prdat$pr[,3,i],
                           d1=prdat$pr[,2,i],
                           a2=prdat$pr[,1,k] - prdat$pr[,3,k],
                           d2=prdat$pr[,2,k])
         oTmp<- cbind(oTmp.xy, oTmp)

         if( !is.null(weights[1]) && is.na(weights[1]) ){
            g0<- lmGls(y~.,data=oTmp,A=gcv)
         }else{
            g0<- lm(y~.,data=oTmp,weights=weights)
         }
         if( !is.null(weights[1]) && is.na(weights[1]) ){
            g<- lmGls(y~(a1+d1)*(a2+d2) + .,data=oTmp,A=gcv)
         }else{
            g<- lm(y~(a1+d1)*(a2+d2) + .,data=oTmp,weights=weights)
         }
         P[i,k]<- 2*(logLik(g)-logLik(g0))
      }
   }

   P
}

scanTwo.2 <-
   function(y,
            x,
            gdat,
            cov,
            numGeno)
{
   if(numGeno){
      num.geno<- I
   }else num.geno<- as.factor
   diag.cov<- diag(as.matrix(cov))
   if( max( abs( cov-diag(diag.cov) ) ) < min(1e-5,1e-5*max(diag.cov)) ){
      if( max(diag.cov-min(diag.cov)) < min(1e-5,1e-5*max(diag.cov)) ){
         weights<- NULL
      }else weights<- 1/diag.cov
   }else weights<- NA
   gcv<- W.inv(cov)

   nsnp<- dim(gdat)[2]
   P<- matrix(NA,nrow=nsnp,ncol=nsnp)
      rownames(P)<- colnames(P)<- colnames(gdat)
   if(nsnp<1) return(NULL)

   if(!missing(x)){
      oTmp.xy<- data.frame(y=y,x)
   }else{
      oTmp.xy<- data.frame(y=y)
   }
   for(i in 1:(nsnp-1)){
      for(j in (i+1):nsnp){
         oTmp<- data.frame(snp1=num.geno(gdat[,i]),
                           snp2=num.geno(gdat[,j]))
         oTmp<- cbind(oTmp.xy, oTmp)

         if( !is.null(weights[1]) && is.na(weights[1]) ){
            g0<- lmGls(y~.,data=oTmp,A=gcv)
         }else{
            g0<- lm(y~.,data=oTmp,weights=weights)
         }
         if( !is.null(weights[1]) && is.na(weights[1]) ){
            g<- lmGls(y~snp1*snp2 + .,data=oTmp,A=gcv)
         }else{
            g<- lm(y~snp1*snp2 + .,data=oTmp,weights=weights)
         }
         P[i,j]<- 2*(logLik(g)-logLik(g0))
      }
   }

   P
}

scanTwo<- 
   function(y,
            x,
            gdat,
            prdat = NULL,
            vc = NULL,
            numGeno = FALSE,
            minorGenoFreq = 0,
            rmv = TRUE)

{
   if(!all(is.finite(y)))
      stop("y: non-numeric or infinite data points not allowed.", call.=FALSE)
   if(!missing(x))
      if(any(sapply(x,is.infinite) | sapply(x,is.na)))
         stop("x: missing or infinite data points not allowed.", call.=FALSE)
   UseMethod("scanTwo")
}

scanTwo.default<- 
   function(y,
            x,
            gdat,
            prdat = NULL,
            vc = NULL,
            numGeno = FALSE,
            minorGenoFreq = 0,
            rmv = TRUE)
{
   if(!is.null(vc)){
      if(is.element("VC",attr(vc,"class"))){
         nv<- length(vc$v)
         nb<- length(vc$par) - nv
         nr<- nrow(vc$y)
         cov<- matrix(0,nrow=nr,ncol=nr)
         for(i in 1:nv)
            cov<- cov + vc$v[[i]]*vc$par[nb+i]
      }else{
         if(is.data.frame(vc)) vc<- as.matrix(vc)
         if(!is.matrix(vc)) stop("'vc' should be a matrix.", call.=FALSE)
         if(!is.numeric(vc)) stop("'vc' should be a numeric matrix.", call.=FALSE)
         cov<- vc
      }
   }else cov<- diag(nrow(as.matrix(y)))
   if(!is.null(prdat)){
      pv<- scanTwo.1(y=y,x=x,prdat=prdat,cov=cov)
   }else{
      gdat<- as.matrix(gdat)
      if(any(is.na(gdat)))
         stop("There are missing genotypes...", call.=FALSE)
      if(numGeno){
         idx<- FALSE
      }else{
         tb<- sort(union(gdat,NULL))
         tbf<- NULL
         for(ii in tb) tbf<- rbind(tbf,colSums(gdat==ii))
            if(sum(tbf)!=nrow(gdat)*ncol(gdat)) stop("Error occurred.", call.=FALSE)
         tbf<- apply(tbf,2,min)
         idx<- (tbf < nrow(gdat)*minorGenoFreq)
      }
      if(sum(idx)>0){
         if(rmv){
            tb<- sort(union(gdat,NULL))
            tbf<- NULL
            for(ii in tb) tbf<- rbind(tbf,colSums(gdat==ii))
               if(sum(tbf)!=nrow(gdat)*ncol(gdat)) stop("Error occurred.", call.=FALSE)
            tbf<- apply(tbf,2,min)
            idx<- (tbf < nrow(gdat)*minorGenoFreq)
            gdat<- as.matrix(gdat[,!idx])
            rm(tb,tbf,ii,idx)
         }else{
            cat("   Minor genotype frequency is too small at one or more SNPs.\a\n")
            return(NULL)
        }
      }

      if(!numGeno)
         gdat<- as.data.frame(gdat)
      pv<- scanTwo.2(y=y,x=x,gdat=gdat,cov=cov,numGeno=numGeno)
   }

   class(pv)<- "scanTwo"
   pv
}

# generalized least squares estimates
gls<- function(formula,data=NULL,vc=NULL,test=c("none","F")){
   cl <- match.call()
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
   mf <- mf[c(1L, m)]
      mf$drop.unused.levels <- TRUE
      mf[[1L]] <- quote(stats::model.frame)
      mf <- eval(mf, parent.frame())
   mt <- attr(mf, "terms")
   yy <- model.response(mf, "numeric")
   xx <- model.matrix(mt, mf)

   nr<- nrow(xx)
   if(!is.null(vc)){
      if(is.element("VC",attr(vc,"class"))){
         nv<- length(vc$v)
         nb<- length(vc$par) - nv
         nr<- nrow(vc$y)
         cov<- matrix(0,nrow=nr,ncol=nr)
         for(i in 1:nv)
            cov<- cov + vc$v[[i]]*vc$par[nb+i]
      }else{
         if(is.data.frame(vc)) vc<- as.matrix(vc)
         if(!is.matrix(vc)) stop("'vc' should be a matrix.", call.=FALSE)
         if(!is.numeric(vc)) stop("'vc' should be a numeric matrix.", call.=FALSE)
         cov<- vc
      }
   }else cov<- diag(nrow(as.matrix(yy)))
   A<- W.inv(cov)

   y<- A%*%yy; attributes(y)<- attributes(yy)
   x<- A%*%xx; attributes(x)<- attributes(xx)
   dtf<- data.frame(y=y,x)

   mdl<- lm.fit(x, y, singular.ok = TRUE)
   class(mdl) <- "lm"
   #mdl$data<- dtf
   mdl$xlevels <- .getXlevels(mt, mf)
   mdl$call <- match.call()
   mdl$terms <- mt
   mdl$model<- mf

#   print(logLik(mdl))
   test<- match.arg(test)
   if(test=="none") summary(mdl)$coeff else
      anova(mdl, test=test)
}

