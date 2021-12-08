
machineEps<- (.Machine$double.eps)^(2/3)
inf<- min(1e+38,sqrt(.Machine$double.xmax))

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
         cat("Warning: matrix not in full rank!\n")

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

   fct.2<- function(y, x, v, pp, intv=c(-25,12), tol=machineEps){
      oo<- optimize(optfct.2, interval=intv, pp, maximum=FALSE, tol=machineEps)

      ny<- pp$ny
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
      ptmp<- rep(NA, 2)
         ptmp[n]<- R/ny
         ptmp[-n]<- dlt*ptmp[n]
      ptmp<- c(b, ptmp)

      oo<- list(value=oo$objective, par=ptmp, eH=pp$eH)

      rm(dlt,n,H,dd,rH,pXrH,b,R,ptmp)

      oo
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
      if(FALSE){
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
      }else{
         invS<- solve(S)
         tXinvS<- t(x)%*%invS
         b<- solve(tXinvS%*%x, tXinvS%*%y)
      }

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

   initparf<- function(y, x, v, nit=nit, opt=1){
      olm<- lm(y~.-1, data=data.frame(y=y,x))
      vr<- var(olm$res)*(ny-1)/ny
      initpar<- c(olm$coef, rep(vr,nv))
      if(opt==1){
         for(i in 1:nv){
            if(match("E",names(v))==i)
               next
            initpar[nb+i]<- initpar[nb+i]/mean(diag(v[[i]]))
         }
         rm(i)
      }else{
         j<- match("E", names(v))
         pt<- 0
         for(i in 1:nv){
            if(i == j) next
            vt<- list(Tmp=v[[i]], E=v[[j]])
            ppt<- ppfct(y=y, x=x, v=vt)
            ot<- fct.2(y=y, x=x, v=vt, pp=ppt, intv=c(-25,12), tol=machineEps)
            initpar[nb+i]<- ot$par[nb+1]
            pt<- c(pt, ot$par[nb+2])
            rm(vt,ppt,ot)
         }
         initpar[nb+j]<- mean(pt, na.rm = FALSE)
         rm(j,pt)
      }
      rm(olm,vr)
      ppar<- initpar
      for(i in 1:nv){
         ppar[nb+i]<- log(ppar[nb+i]+machineEps)
         rm(i)
      }

      ppar
   }

   oof<- function(ppar, y, x, v, nb=nb, nv=nv, nit=nit){
      val<- inf
      if(missing(nit)) cnt<- 25 else
         cnt<- nit
      lmt<- ppar[-(1:nb)]
         names(lmt)<- names(v)
      while(cnt>nit*1/2){
         # the following is much better than optim("Nelder-Mead")
         ppar[1:nb]<- optfct.b(ppar[-(1:nb)], y, x, v)
         oo<- nlminb(ppar[-c(1:nb)], optfct.v, bp=ppar[1:nb], y=y, x=x, v=v,
            lower=-19, upper=pmin(12,lmt+7))
         if(length(grep("Error",oo))>0){
            break
         }
         oo$value<- oo$objective
         cnt<- cnt-1
         if(abs(val-oo$val) < machineEps^0.75 && max(abs(ppar[-c(1:nb)]-oo$par)) < machineEps^0.5)
            break
         val<- oo$value
         ppar[-c(1:nb)]<- oo$par
      }

      cntr<- control
         cntr$maxit<- 500
         cntr$ndeps<- rep(1e-5, nb+nv)
         cntr$factr<- 1e5
      mtd<- "L-BFGS-B"
      while(cnt>0){
         lmt<- ppar[-(1:nb)]
         lw<- pmax(lmt - 5, -19)
            lw<- c(ppar[(1:nb)] - 25, lw)
         up<- pmin(lmt + 5, 12)
            up<- c(ppar[(1:nb)] + 25, up)
         oo<- try(
            optim(ppar, optfct, y=y, x=x, v=v, method=mtd, control=cntr, hessian=FALSE,
               lower=lw, upper=up), silent=TRUE
         )
         if(length(grep("Error",oo))>0){
            mtd<- "Nelder-Mead"
            cntr$maxit<- 100
            cntr$reltol<- machineEps
            next
         }
         cnt<- cnt-1
         if(abs(val-oo$val) < machineEps^0.75 && max(abs(ppar-oo$par)) < machineEps^0.5)
            break
         val<- oo$value
         ppar<- oo$par
      }
      if(cnt == 0){
         cat("   Optimization might have failed. Another run with 'initpar' being the resulting parameter values should be helpful.\a\n")
      }

      for(i in 1:nv){
         oo$par[nb+i]<- exp(oo$par[nb+i])
      }

      oo
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
   nv<- length(v)
      idx<- rep(TRUE,nv)
      for(n in 1:nv)
         if(is.null(v[[n]])) idx[n]<- FALSE
      v<- v[idx]
      nv<- length(v)
   rm(idx)

   if(nv == 2){
      if(ny > 5000){
         cat("   May take a long time due to large sample size...\n")
      }else if(ny > 3000){
         cat("   May take a while due to large sample size. Please be patient...\n")
      }
      pp<- ppfct(y=y,x=x,v=v)
      oo<- fct.2(y=y, x=x, v=v, pp=pp, intv=c(-25,12), tol=machineEps)
   }else{
      if(ny > 1500){
         cat("   May take a long time due to large sample size...\n")
      }else if(ny > 900){
         cat("   May take a while due to large sample size. Please be patient...\n")
      }
      if(missing(initpar)){
         pparTmp<- initparf(y, x, v, nit=nit, opt=1)
         oo<- oof(pparTmp, y, x, v, nb=nb, nv=nv, nit=nit)
      }else{
         if(length(initpar) != nb+nv)
            stop(paste("initpar: ", nb+nv, " parameters only!", sep=""), call.=FALSE)
         for(i in 1:nv) if(initpar[nb+i] < 0)
            stop("initpar: negative variance components?!", call.=FALSE)
         rm(i)

         pparTmp<- initpar
         for(i in 1:nv){
            pparTmp[nb+i]<- log(pparTmp[nb+i]+machineEps)
            rm(i)
         }
         oo<- oof(pparTmp, y, x, v, nb=nb, nv=nv, nit=nit)
      }
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
      cntr<- control
      if(is.null(cntr$ndeps)) cntr$ndeps<- rep(5e-6,length(ppar)) else
         cntr$ndeps<- pmin(cntr$ndeps, 5e-6)
      oo$hessian<- -optimHess(ppar, hessf, y=y, x=x, v=v, control=cntr)
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
# object: object from estVC
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
   rd<- solve(S,rr)
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

