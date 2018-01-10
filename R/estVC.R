
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
   if(is.matrix(v)) v<- list(v)
   if(!is.null(v[["EE"]])){
      cat("   We now use 'E' (not 'EE') for residual variance matrix; see documentation.\a\n")
   }
   control$fnscale<- 1

   .gls<- function(formula,data=NULL,vc=NULL){
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
         if(is.list(vc)){
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

      mdl
   }

   ginv<- function(W, symmetric=TRUE){
      eW <- eigen(W, symmetric=symmetric)
      d <- eW$values
      if(min(d) < .Machine$double.eps){
         pd<- FALSE
         cat("'W' not positive definite.\n")

         d[d < .Machine$double.eps]<- Inf
      }
      A <- eW$vector%*%diag(1/d) %*% t(eW$vector)
      A
   }

   ppfct<- function(y,x,v){# see Kang et al 2008 Genetics
      ny<- length(y)

      n<- match("E",names(v))
      Q<- chol(v[[n]],pivot=FALSE)
      rQ<- solve(Q)
      eH<- eigen(t(rQ)%*%v[[-n]]%*%rQ, symmetric=TRUE)

      S<- x%*%ginv(t(x)%*%x)%*%t(x)
         S<- diag(1,ny) - S
      eQ<- eigen(S, symmetric=TRUE) # note: S*S = S
      rQ<- eQ$vector
      eS<- eigen(t(rQ)%*%(S%*%v[[-n]]%*%S)%*%rQ, symmetric=TRUE)

      U<- t(eS$vector)%*%t(rQ)
      eta<- U%*%y

      list(ny=ny, eta=eta, eH=eH, eS=eS, eQ=eQ)
   }

   optfct.2<- function(par,pp){
      dd<- pp$eS$value*exp(par)+pp$eQ$value
      idx<- dd > machineEps
      tmp<- -pp$ny*log(pp$ny/2/pi) + pp$ny
         tmp<- tmp + pp$ny*log(sum(pp$eta[idx]^2/dd[idx]))
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
      nb<- length(par) - nv

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

   ny<- length(y)
   nb<- ncol(x)
   if(is.null(v[["E"]])){# add residual variance
      v$E<- diag(ny)
      olm<- lm(y~.-1, data=data.frame(y=y,x))
   }else{
      olm<- .gls(y~.-1, data=data.frame(y=y,x), vc=v$E)
   }
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
         #oo<- nlminb(c(-9,log(vr)), optfct.v, bp=olm$coef, y=y, x=x, v=c(v[n],v["E"]),
         #   lower=log(vr)-c(12,0), upper=log(vr)+c(3,0))
         #if(oo$objective > inf-1)
         #   cat("   Optimization possibly failed.\a\n")
         #initpar[nb+n]<- exp(oo$par[1])
         #rm(oo,n)
         ### the following seems to save a little bit of time
         oo<- .gls(y~.-1, data=data.frame(y=y,x), vc=v[[n]]+v[["E"]])
         vrTmp<- var(oo$res)*(ny-1)/ny
            vrTmp<- ifelse(vrTmp + machineEps > vr, vr*machineEps,
               (vr-vrTmp)*mean(diag(v[["E"]]))/mean(diag(v[[n]]))
            )
         initpar[nb+n]<- min(max(vrTmp,1e-5), 1e+5)
         rm(oo,vrTmp,n)
      }
   }else{
      if(length(initpar) < nb+length(v))
         initpar<- c(initpar,vr)
   }

   olm<- .gls(y~.-1, data=data.frame(y=y,x), vc=list(par=initpar,y=as.matrix(y),v=v))
   initpar<- c(olm$coef, initpar[-(1:nb)]*var(olm$res)*(ny-1)/ny)

   rm(olm,vr)

   for(i in 1:nv){
      initpar[nb+i]<- log(initpar[nb+i])
   }
   rm(i)

if(nv > 2){
   lmt<- initpar[-(1:nb)]
      names(lmt)<- names(v)
   pparTmp<- initpar
      oo<- nlminb(pparTmp[-c(1:nb)], optfct.v, bp=pparTmp[1:nb], y=y, x=x, v=v,
         lower=-19, upper=pmax(7,lmt+5))
      pparTmp[-c(1:nb)]<- oo$par
   lmt<- pparTmp[-(1:nb)]
      names(lmt)<- names(v)
   val<- inf
   if(missing(nit)) cnt<- 25 else
      cnt<- nit
   while(cnt>0){
      pparTmp[1:nb]<- optfct.b(pparTmp[-(1:nb)], y, x, v)
      oo<- nlminb(pparTmp[-c(1:nb)], optfct.v, bp=pparTmp[1:nb], y=y, x=x, v=v,
         lower=-19, upper=pmax(7,lmt+5))
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
}else if(nv < 2){
   oo<- optim(initpar, optfct, y=y, x=x, v=v, method="Nelder-Mead",
         control=list(fnscale=1,maxit=3), hessian=FALSE)
      for(i in 1:nv){
         oo$par[nb+i]<- exp(oo$par[nb+i])
      }
}else{
   oo<- optimize(optfct.2, interval=c(-25,25), ppfct(y=y,x=x,v=v),
      maximum=FALSE, tol=machineEps)

   n<- match("E",names(v))
   dlt<- exp(oo$minimum)
   H<- v[[-n]]*dlt + v[[n]]
   rH<- solve(H)
   pXrH<- t(x)%*%rH
   b<- solve(pXrH%*%x,pXrH%*%y)
   R<- y - x%*%b
      R<- t(R)%*%rH%*%R
   ptmp<- initpar[-(1:nb)]
      ptmp[n]<- R/ny
      ptmp[-n]<- dlt*ptmp[n]
   initpar[1:nb]<- as.vector(b)
   initpar[-(1:nb)]<- ptmp

   oo<- list(value=oo$objective, par=initpar)

   rm(dlt,n,H,rH,pXrH,b,R,ptmp)
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

est.VC.1 <-
   function(y,
            x,
            v,
            initpar,
            nit,
            control,
            hessian)
{### Fast; Trouble if E parameter is (nearly) 0 !!!
# estimate all background genetic variance (VC)
# y: vector, response
# x: covariates
# v: list of variance components
# initpar: initial parameters, will be initilized automatically if missing
# nit: number of iterations to call optim()
# control: A list of control parameters
# hessian: logical. should a numerically differentiated Hessian matrix be returned?
   if(!is.null(v[["EE"]])){
      cat("   We now use 'E' (not 'EE') for residual variance matrix; see documentation.\a\n")
   }
   control$fnscale<- 1

   ny<- length(y)
   nb<- ncol(x)
   if(is.null(v[["E"]])){ # add residual variance
      v$E<- diag(ny)
      if(!missing(initpar))
         initpar<- c(initpar,var(y))
   }
   nv<- length(v)
      idx<- rep(TRUE,nv)
      for(n in 1:nv)
         if(is.null(v[[n]])) idx[n]<- FALSE
      v<- v[idx]
      nv<- length(v)
   if(missing(initpar)){
      initpar<- c(rep(mean(y),nb),rep(var(y),nv))
      for(n in 1:nv){
         initpar[nb+n]<- initpar[nb+n]/mean(diag(v[[n]]))
      }
   }
   if(nv>1){
      n<- match("E",names(v))
      pTmp<- initpar[-(1:nb)]
         pTmp[-n]<- pTmp[-n]/pTmp[n]
      initpar[-(1:nb)]<- pTmp
      rm(pTmp)
   }else stop("There is no random effect; use other method instead.", call.=FALSE)
   rm(idx,n)

   for(i in 1:nv){
      initpar[nb+i]<- log(initpar[nb+i])
   }
   rm(i)

   Winv<- function(W, symmetric=TRUE, inverse=TRUE){
      eW <- eigen(W, symmetric=symmetric)
      d <- eW$values
      if (min(d) <0  && abs(min(d))>machineEps)
          stop("'W' is not positive definite.", call.=FALSE)
      else d[d<=0]<- ifelse(inverse, Inf, 0)
      A <- diag(d^ifelse(inverse, -0.5, 0.5)) %*% t(eW$vector)
      list(A=A, d=d) # t(A)%*%A = W^{-1}
   }
   Sf<- function(par,v){
   # par: not include E
      S<- v$E
      j<- match("E",names(v))
      i<- 0
      for(t in 1:length(v)){
         if(t != j){
            i<- i+1
            S<- S + v[[i]]*exp(par[i])
         }
      }
      S
   }
   optf<- function(par,y,x,v,opt=FALSE){
   # par: not include E
      S<- Sf(par,v)
      WD<- Winv(S)
      dtf<- data.frame(y=y,x)
      o<- lmGls(y~.-1,data=dtf,A=WD$A)

      s2<- sum(o$res^2)/ny
      n<- match("E",names(v))
      coef<- rep(NA,length(v))
         coef[-n]<- exp(par)*s2
         coef[n]<- s2
      names(coef)<- names(v)
      coef<- c(o$coef, coef)

      ny<- length(y)
      lik<- log(2*pi)*(ny/2) + ny/2*(log(s2)+1) + sum(log(WD$d))/2

      if(opt){
         lik
      }else{
         list(value=-lik, parameters=coef)
      }
   }

   pTmp<- initpar[-(1:nb)]
      idx<- match("E",names(v))
      pTmp<- pTmp[-idx]
   val<- inf
   if(missing(nit)) cnt<- 25 else
      cnt<- nit
   while(cnt > 0){
      if(nv > 2){
         oo<- optim(pTmp, optf, y=y, x=x, v=v, opt=TRUE, method="Nelder-Mead",
            control=control, hessian=FALSE)
      }else if(nv == 2){
         oo<- nlminb(pTmp, optf, y=y, x=x, v=v, opt=TRUE, lower=-19, upper=19)
         oo$value<- oo$objective
      }
      if(abs(val-oo$val) < machineEps && max(abs(pTmp-oo$par)) < machineEps^0.75)
         cnt<- 0
      pTmp<- oo$par
      val<- oo$value
      cnt<- cnt-1
   }

   oo<- optf(pTmp,y=y,x=x,v=v,opt=FALSE)

if(hessian){
   hessf<- function(par,y,x,v){
      ny<- length(y)
      nb<- ncol(x)
      nv<- length(v)

      b<- par[1:nb]
      S<- matrix(0,nrow=ny,ncol=ny)
      for(i in 1:nv){
         S<- S + v[[i]]*par[nb+i]
      }
      if(!all(is.finite(S))) return(inf)

      u<- x%*%b
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr)) return(inf)
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp)) return(inf)

      tmp
   }
   ppar<- oo$par
   for(i in 1:nv)
      ppar[nb+i]<- pmax(ppar[nb+i], 1e-5)
   control$ndeps<- rep(5e-6,length(ppar))
   oo$hessian<- -optimHess(ppar, hessf, y=y, x=x, v=v, control=control)
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

est.VC.2 <-
   function(y,
            x,
            v,
            initpar,
            nit,
            control,
            hessian)
{### slow !!!
# estimate all background genetic variance (VC)
# y: vector, response
# x: covariates
# v: list of variance components
# initpar: initial parameters, will be initilized automatically if missing
# nit: number of iterations to call optim()
# control: A list of control parameters
# hessian: logical. should a numerically differentiated Hessian matrix be returned?
   if(!is.null(v[["EE"]])){
      cat("   We now use 'E' (not 'EE') for residual variance matrix; see documentation.\a\n")
   }
   control$fnscale<- 1

   ny<- length(y)
   nb<- ncol(x)
   if(is.null(v[["E"]])){ # add residual variance
      v$E<- diag(ny)
      if(!missing(initpar))
         initpar<- c(initpar,var(y))
   }
   nv<- length(v)
      idx<- rep(TRUE,nv)
      for(n in 1:nv)
         if(is.null(v[[n]])) idx[n]<- FALSE
      v<- v[idx]
      nv<- length(v)
   if(missing(initpar)){
      initpar<- c(rep(mean(y),nb),rep(var(y),nv))
      for(n in 1:nv){
         initpar[nb+n]<- initpar[nb+n]/mean(diag(v[[n]]))
      }
   }
   if(nv>1){
      n<- match("E",names(v))
      pTmp<- initpar[-(1:nb)]
         pTmp[-n]<- pTmp[-n]/pTmp[n]
      initpar[-(1:nb)]<- pTmp
      rm(pTmp)
   }else stop("There is no random effect; use other method instead.", call.=FALSE)
   rm(idx,n)

   for(i in 1:nv){
      initpar[nb+i]<- log(initpar[nb+i])
   }
   rm(i)

   Winv<- function(W, symmetric=TRUE, inverse=TRUE){
      eW <- eigen(W, symmetric=symmetric)
      d <- eW$values
      if (min(d) <0  && abs(min(d))>machineEps)
          stop("'W' is not positive definite.", call.=FALSE)
      else d[d<=0]<- ifelse(inverse, Inf, 0)
      A <- diag(d^ifelse(inverse, -0.5, 0.5)) %*% t(eW$vector)
      list(A=A, d=d) # t(A)%*%A = W^{-1}
   }
   optf<- function(par,y,x,v,opt=FALSE){
   # par: not include E
      ny<- length(y)
      nv<- length(v)
      cv<- log(var(y)) - 10
      if(RPAR < cv){
         S<- matrix(0, nrow=ny, ncol=ny)
      }else{
         S<- v$E
      }
      j<- match("E",names(v))
      i<- 0
      for(t in 1:length(v)){
         if(t != j){
            i<- i+1
            S<- S + v[[i]]*exp(par[i])
         }
      }
      rm(i,j,t)
      WD<- Winv(S)
      dtf<- data.frame(y=y,x)
      o<- lmGls(y~.-1,data=dtf,A=WD$A)

      s2<- sum(o$res^2)/ny
      n<- match("E",names(v))
      coef<- rep(NA,length(v))
         coef[-n]<- exp(par)
         if(RPAR >= cv){
            coef[-n]<- coef[-n]*s2
            coef[n]<- RPAR<<- s2
         }else coef[n]<- RPAR<<- cv
      names(coef)<- names(v)
      coef<- c(o$coef, coef)

      if(RPAR < cv){
         u<- x%*%o$coef
         tmp<- qr(S)
         if(tmp$rank<ncol(tmp$qr)) return(inf)
         ddtmp<- abs(diag(tmp$qr))
         lik<- log(2*pi)*(ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
         if(!is.finite(tmp)) return(inf)
      }else{
         lik<- log(2*pi)*(ny/2) + ny/2*(log(s2)+1) + sum(log(WD$d))/2
      }
      if(opt){
         lik
      }else{
         list(value=-lik, parameters=coef)
      }
   }

   pTmp<- initpar[-(1:nb)]
      idx<- match("E",names(v))
      RPAR<- pTmp[idx]
      pTmp<- pTmp[-idx]
   val<- inf
   if(missing(nit)) cnt<- 25 else
      cnt<- nit
   while(cnt > 0){
      if(nv > 2){
         oo<- optim(pTmp, optf, y=y, x=x, v=v, opt=TRUE, method="Nelder-Mead",
            control=control, hessian=FALSE)
      }else if(nv == 2){
         oo<- nlminb(pTmp, optf, y=y, x=x, v=v, opt=TRUE, lower=-19, upper=19)
         oo$value<- oo$objective
      }
      if(abs(val-oo$val) < machineEps && max(abs(pTmp-oo$par)) < machineEps^0.75)
         cnt<- 0
      pTmp<- oo$par
      val<- oo$value
      cnt<- cnt-1
   }

   oo<- optf(pTmp,y=y,x=x,v=v,opt=FALSE)

if(hessian){
   hessf<- function(par,y,x,v){
      ny<- length(y)
      nb<- ncol(x)
      nv<- length(v)

      b<- par[1:nb]
      S<- matrix(0,nrow=ny,ncol=ny)
      for(i in 1:nv){
         S<- S + v[[i]]*par[nb+i]
      }
      if(!all(is.finite(S))) return(inf)

      u<- x%*%b
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr)) return(inf)
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp)) return(inf)

      tmp
   }
   ppar<- oo$par
   for(i in 1:nv)
      ppar[nb+i]<- pmax(ppar[nb+i], 1e-5)
   control$ndeps<- rep(5e-6,length(ppar))
   oo$hessian<- -optimHess(ppar, hessf, y=y, x=x, v=v, control=control)
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

est.VC.3 <-
   function(y,
            x,
            v,
            initpar,
            nit,
            control,
            hessian)
{### Not very accurate !!!
# estimate all background genetic variance (VC)
# y: vector, response
# x: covariates
# v: list of variance components
# initpar: initial parameters, will be initilized automatically if missing
# nit: number of iterations to call optim()
# control: A list of control parameters
# hessian: logical. should a numerically differentiated Hessian matrix be returned?
   if(!is.null(v[["EE"]])){
      cat("   We now use 'E' (not 'EE') for residual variance matrix; see documentation.\a\n")
   }
   control$fnscale<- 1

   ny<- length(y)
   nb<- ncol(x)
   if(is.null(v[["E"]])){ # add residual variance
      v$E<- diag(ny)
      if(!missing(initpar))
         initpar<- c(initpar,var(y))
   }
   nv<- length(v)
      idx<- rep(TRUE,nv)
      for(n in 1:nv)
         if(is.null(v[[n]])) idx[n]<- FALSE
      v<- v[idx]
      nv<- length(v)
   if(missing(initpar)){
      initpar<- c(rep(mean(y),nb),rep(var(y),nv))
      for(n in 1:nv){
         initpar[nb+n]<- initpar[nb+n]/mean(diag(v[[n]]))
      }
   }
   if(nv>1){
      n<- match("E",names(v))
      pTmp<- initpar[-(1:nb)]
         pTmp[-n]<- pTmp[-n]/pTmp[n]
      initpar[-(1:nb)]<- pTmp
      rm(pTmp)
   }else stop("There is no random effect; use other method instead.", call.=FALSE)
   rm(idx,n)

   for(i in 1:nv){
      initpar[nb+i]<- log(initpar[nb+i])
   }
   rm(i)

   Winv<- function(W, symmetric=TRUE, inverse=TRUE){
      eW <- eigen(W, symmetric=symmetric)
      d <- eW$values
      if (min(d) <0  && abs(min(d))>machineEps)
          stop("'W' is not positive definite.", call.=FALSE)
      else d[d<=0]<- ifelse(inverse, Inf, 0)
      A <- diag(d^ifelse(inverse, -0.5, 0.5)) %*% t(eW$vector)
      list(A=A, d=d) # t(A)%*%A = W^{-1}
   }
   optf<- function(par,y,x,v,opt=FALSE){
      ny<- length(y)
      nb<- ncol(x)
      nv<- length(v)

      S<- matrix(0,nrow=ny,ncol=ny)
      for(i in 1:nv){
         S<- S + v[[i]]*exp(par[i])
      }
      if(!all(is.finite(S))) return(inf)
      WD<- Winv(S)
      dtf<- data.frame(y=y,x)
      o<- lmGls(y~.-1,data=dtf,A=WD$A)

      s2<- sum(o$res^2)/ny
      ny<- length(y)
      lik<- log(2*pi)*(ny/2) + ny/2*(log(s2)+1) + sum(log(WD$d))/2

      if(opt){
         lik
      }else{
         coef<- exp(par)
         names(coef)<- names(v)
         coef<- c(o$coef, coef)

         list(value=-lik, parameters=coef)
      }
   }

   pTmp<- initpar[-(1:nb)]
   val<- inf
   if(missing(nit)) cnt<- 25 else
      cnt<- nit
   while(cnt > 0){
      if(nv > 1){
         oo<- optim(pTmp, optf, y=y, x=x, v=v, opt=TRUE, method="Nelder-Mead",
            control=control, hessian=FALSE)
      }else{
         oo<- nlminb(pTmp, optf, y=y, x=x, v=v, opt=TRUE, lower=-19, upper=19)
         oo$value<- oo$objective
      }
      if(abs(val-oo$val) < machineEps && max(abs(pTmp-oo$par)) < machineEps^0.75)
         cnt<- 0
      pTmp<- oo$par
      val<- oo$value
      cnt<- cnt-1
   }

   oo<- optf(pTmp,y=y,x=x,v=v,opt=FALSE)

if(hessian){
   hessf<- function(par,y,x,v){
      ny<- length(y)
      nb<- ncol(x)
      nv<- length(v)

      b<- par[1:nb]
      S<- matrix(0,nrow=ny,ncol=ny)
      for(i in 1:nv){
         S<- S + v[[i]]*par[nb+i]
      }
      if(!all(is.finite(S))) return(inf)

      u<- x%*%b
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr)) return(inf)
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp)) return(inf)

      tmp
   }
   ppar<- oo$par
   for(i in 1:nv)
      ppar[nb+i]<- pmax(ppar[nb+i], 1e-5)
   control$ndeps<- rep(5e-6,length(ppar))
   oo$hessian<- -optimHess(ppar, hessf, y=y, x=x, v=v, control=control)
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

