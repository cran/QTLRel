
machineEps<- (.Machine$double.eps)^(2/4)
inf<- max(1e+38,sqrt(.Machine$double.xmax))

# one of the main functions, Nelder-Mead method
estVC <-
   function(y,
            x,
            v = list(E=diag(length(y))),
            initpar,
            nit = 25,
            control = list(),
            hessian = FALSE)
{
   UseMethod("estVC")
}

estVC.default <-
   function(y,
            x,
            v,
            initpar,
            nit = 25,
            control = list(),
            hessian = FALSE)
{
   if(!all(is.finite(y)))
      stop("y: non-numeric or infinite data points not allowed.", call.=FALSE)
   if(!missing(x))
      if(any(sapply(x,is.infinite) | sapply(x,is.na)))
         stop("x: missing or infinite data points not allowed.", call.=FALSE)

   if(!is.null(dim(y))){
      if(length(dim(y))>2) stop("y: may be wrong.", call.=FALSE)
      if(dim(y)[2]>1)
         cat("   Warning: y: only the fisrt column will be analyzed.\a\n")
      y<- y[,1]
   }
   xyf<- function(formula, data, contrasts = NULL){
      mf<- match.call(expand.dots = FALSE)
         m<- match(c("formula", "data"), names(mf), 0L)
         mf<- mf[c(1L, m)]
         mf$drop.unused.levels<- TRUE
         mf[[1L]]<- quote(stats::model.frame)
         mf<- eval(mf, parent.frame())
      mt<- attr(mf, "terms")
      x<- model.matrix(mt, mf, contrasts)
      y<- model.response(mf, "numeric")

      list(y=y, x=x)
   }
   if(!missing(x)){
      oTmp<- data.frame(y=y,x)
   }else{
      oTmp<- data.frame(y=y)
   }
   xy<- xyf(y~., data=oTmp)

   est.VC(y = xy$y,
          x = xy$x,
          v = v,
          initpar = initpar,
          nit = nit,
          control = control,
          hessian = hessian)
}

est.VC <-
   function(y,
            x,
            v,
            initpar,
            nit,
            control,
            hessian)
{
# estimate all background genetic variance (VC)
# y: vector, response
# x: covariates
# v: list of variance components (optionally with residual E)
# initpar: initial parameters, will be initilized automatically if missing
# nit: number of iterations to call optim()
# control: A list of control parameters
# hessian: logical. should a numerically differentiated Hessian matrix be returned?
   if(is.matrix(v)) v<- list(A=v)
   if(!is.null(v[["EE"]])){
      cat("   We now use 'E' (not 'EE') for residual variance matrix; see documentation.\a\n")
   }
   control$fnscale<- 1

   ginv<- function(W, symmetric=TRUE){
      eW <- eigen(W, symmetric=symmetric)
      d <- eW$values
      if(min(d) < .Machine$double.eps){
         pd<- FALSE
         cat("'W' not positive definite.\n")

         d[d < .Machine$double.eps]<- Inf
      }
      A <- sweep(eW$vector,2,d,"/") %*% t(eW$vector)
      A
   }

   ppfct<- function(y,x,v){# see Kang et al 2008 Genetics
      ny<- length(y)

      n<- match("E",names(v))
      eH<- eigen(v[[-n]], symmetric=TRUE)

      S<- x%*%ginv(t(x)%*%x)%*%t(x)
         S<- diag(1,ny) - S
      eS<- eigen(S%*%v[[-n]]%*%S, symmetric=TRUE)
         idx<- abs(eS$values) > machineEps
         eS$values<- eS$values[idx]
         eS$vectors<- eS$vectors[,idx]

      U<- t(eS$vector)
      eta<- U%*%y

      list(ny=ny, eta=eta, eH=eH, eS=eS)
   }

   optfct.2<- function(par,pp){# assume E = diag(1,n)
      dd<- pp$eS$value*exp(par)+1
      tmp<- -pp$ny*log(pp$ny/2/pi) + pp$ny
         tmp<- tmp + pp$ny*log(sum(pp$eta^2/dd))
         tmp<- tmp + sum(log(pp$eH$value*exp(par)+1))
      tmp/2
   }

   optfct<- function(par,y,x,v){
      ny<- length(y)
      nv<- length(v)
      nb<- length(par) - nv
      b<- par[1:nb]

      S<- matrix(0,nrow=ny,ncol=ny)
      for(i in 1:nv){
         S<- S + v[[i]]*exp(par[nb+i])
      }
      if(!all(is.finite(S)))
         return(inf)

      u<- x%*%b
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr))
         return(inf)
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp))
         return(inf)

      tmp
   }
   optfct.b<- function(vp,y,x,v){
      ny<- length(y)
      nv<- length(v)

      S<- matrix(0,nrow=ny,ncol=ny)
      for(i in 1:nv){
         S<- S + v[[i]]*exp(vp[i])
      }
      if(!all(is.finite(S)))
         return(inf)

      # it is slower to use lmGls
      dd<- eigen(S,symmetric=T)
         uu<- dd$vec
         dd<- abs(dd$val)
         dd[dd<machineEps]<- machineEps
      yy<- t(uu)%*%y
         yy<- as.matrix(yy)
      xx<- t(uu)%*%x
         xx<- as.matrix(xx)
      b<- lm.wfit(x=xx,y=yy,w=1/dd)$coef
         b[!is.finite(b)]<- 0
      b
   }
   optfct.v<- function(par,bp,y,x,v){
      ny<- length(y)
      nv<- length(v)

      S<- matrix(0,nrow=ny,ncol=ny)
      for(i in 1:nv){
         S<- S + v[[i]]*exp(par[i])
      }
      if(!all(is.finite(S)))
         return(inf)

      u<- x%*%bp
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr))
         return(inf)
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp))
         return(inf)

      tmp
   }

   ISDIAG<- FALSE
   EDIAG<- 1
   ny<- length(y)
   nb<- ncol(x)
   if(is.null(v[["E"]])){# add residual variance
      v$E<- diag(ny)
      ISDIAG<- TRUE
   }
   if(!ISDIAG){
      EDIAG<- diag(v$E)
      if(any(EDIAG < machineEps))
         stop("E: not positive!", call.=FALSE)
      EDIAG<- mean(EDIAG)
      if(max(abs(v$E - diag(EDIAG,ny))) < machineEps){
         ISDIAG<- TRUE
         if(abs(EDIAG -1) > machineEps)
            v$E<- diag(ny)
      }
   }
   olm<- lm(y~.-1, data=data.frame(y=y,x))
   vr<- var(olm$res)*(ny-1)/ny
   nv<- length(v)
      idx<- rep(TRUE,nv)
      for(n in 1:nv)
         if(is.null(v[[n]])) idx[n]<- FALSE
      v<- v[idx]
      nv<- length(v)
   rm(idx)

   if(missing(initpar)){
      initpar<- c(olm$coef, rep(vr,nv))
      if(nv > 2) for(n in 1:nv){# skip nv = 2 due to optfct.2()
         if(match("E",names(v))==n)
            next
         initpar[nb+n]<- initpar[nb+n]/mean(diag(v[[n]]))
      }
   }else{
      if(length(initpar) < nb+length(v))
         initpar<- c(initpar,vr)
   }

   if(nv != 2){
      if(ny > 1500){
         warning("The sample size is large so it may take time to finish...
         You might consider using other software.")
      }else if(ny > 900){
         warning("This may take a while...Please be patient.")
      }
      for(i in 1:nv){
         initpar[nb+i]<- log(initpar[nb+i])
      }
      rm(i)

      lmt<- initpar[-(1:nb)]
         names(lmt)<- names(v)
      pparTmp<- initpar
         oo<- nlminb(pparTmp[-c(1:nb)], optfct.v, bp=pparTmp[1:nb], y=y, x=x, v=v,
            lower=-19, upper=pmin(12,lmt+7))
         pparTmp[-c(1:nb)]<- oo$par
      lmt<- pparTmp[-(1:nb)]
         names(lmt)<- names(v)
      val<- inf
      if(missing(nit)) cnt<- 25 else
         cnt<- nit
      while(cnt>0){
         pparTmp[1:nb]<- optfct.b(pparTmp[-(1:nb)], y, x, v)
         oo<- nlminb(pparTmp[-c(1:nb)], optfct.v, bp=pparTmp[1:nb], y=y, x=x, v=v,
            lower=-19, upper=pmin(12,lmt+7))
            oo$value<- oo$objective
         if(abs(val-oo$val) < machineEps && max(abs(pparTmp[-c(1:nb)]-oo$par)) < machineEps^0.75)
            cnt<- 0
         pparTmp[-c(1:nb)]<- oo$par
         val<- oo$value
         cnt<- cnt-1
      }
      if(cnt == 0){
         cat("   Warning: optimization possibly failed. A larger 'nit' should be helpful.\a\n")
      }

      control$maxit<- 100
      oo<- optim(pparTmp, optfct, y=y, x=x, v=v, method="Nelder-Mead",
            control=control, hessian=FALSE)
         for(i in 1:nv){
            oo$par[nb+i]<- exp(oo$par[nb+i])
         }
   }else{
      pp<- ppfct(y=y,x=x,v=v)
      oo<- optimize(optfct.2, interval=c(-25,12), pp, maximum=FALSE, tol=machineEps)

      n<- match("E",names(v))
      dlt<- exp(oo$minimum) # note: this can't be too large
      H<- v[[-n]]*dlt + v[[n]]
      # rH<- solve(H)
      dd<- pp$eH$value*dlt + 1
      rH<- sweep(pp$eH$vectors, 2, dd, "/") %*% t(pp$eH$vectors)
      pXrH<- t(x)%*%rH
      b<- solve(pXrH%*%x,pXrH%*%y)
      R<- y - x%*%b
         R<- t(R)%*%rH%*%R
      ptmp<- initpar[-(1:nb)]
         ptmp[n]<- R/ny
         ptmp[-n]<- dlt*ptmp[n]
      initpar[1:nb]<- as.vector(b)
      initpar[-(1:nb)]<- ptmp

      oo<- list(value=oo$objective, par=initpar, eH=pp$eH)

      rm(dlt,n,H,dd,rH,pXrH,b,R,ptmp)
   }

   if(ISDIAG){
      n<- match("E",names(v))
      oo$par[nb+n]<- oo$par[nb+n]/EDIAG
      rm(n)
   }

   if(hessian){
      hessf<- function(par,y,x,v){
         ny<- length(y)
         nv<- length(v)
         nb<- length(par) - nv
         b<- par[1:nb]

         S<- matrix(0,nrow=ny,ncol=ny)
         for(i in 1:nv){
            S<- S + v[[i]]*par[nb+i]
         }
         #if(!all(is.finite(S))) return(inf)

         u<- x%*%b
         tmp<- qr(S)
         #if(tmp$rank<ncol(tmp$qr)) return(inf)
         ddtmp<- abs(diag(tmp$qr))
         tmp<- log(2*pi)*(ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
         #if(!is.finite(tmp)) return(inf)

         tmp
      }
      ppar<- oo$par
      for(i in 1:nv)
         ppar[nb+i]<- pmax(ppar[nb+i], 1e-5)
      control$ndeps<- rep(5e-6,length(ppar))
      oo$hessian<- -optimHess(ppar, hessf, y=y, x=x, v=v, control=control)
   }

   if(nv == 2){
      n<- match("E",names(v))
      pTmp<- oo$par[-c(1:nb)]
      dd<- pp$eH$value*pTmp[-n] + pTmp[n]
      oo$hinvCov<- sweep(t(oo$eH$vectors), 1, dd^(-0.5), "*")
      rm(n,pTmp,dd)
   }

   oo$value<- -oo$value
   attributes(oo$value)<- NULL
   names(oo$par)[1:nb]<- colnames(x)
   names(oo$par)[nb+1:nv]<- names(v)
   oo$y<- as.matrix(y)
   oo$x<- as.matrix(x)
   oo$v<- v
   class(oo)<- "VC"
   oo
}

blup<- function(object)
{
   UseMethod("blup")
}

blup.VC<- function(object){
# best linear unbiased prediction (BLUP) for all effects
# y: ny by 1 matrix
# x: design matrix including overall mean !!!
# object: object from estBVG
# v: list of variance components corresponding object$par[-c(1:nb)]
   ny<- nrow(object$y)
   nb<- ncol(object$x)
   nv<- length(object$v)
   idx<- NULL
   for(i in 1:nv)
      if(!is.null(object$v[[i]])) idx<- c(idx,i)
   object$v<- object$v[idx]
   nv<- length(object$v)

   S<- matrix(0,nrow=ny,ncol=ny)
   for(i in 1:nv)
      S<- S + object$v[[i]]*object$par[nb+i]

   out<- vector("list",nv)
   rr<- object$y-object$x%*%object$par[1:nb]
   rd<- solve(S)%*%rr
   out[[1]]<- sweep(object$x,2,object$par[1:nb],"*")
      names(out)[1]<- "fixed"
   for(i in 1:nv){
      out[[1+i]]<- object$v[[i]]%*%rd*object$par[nb+i]
      names(out)[1+i]<- names(object$v)[i]
   }

   out
}

print.VC<- function(x,...){
   cat("value:\n"); print(x$value)
   cat("parameters:\n"); print(x$par)
}

