
machineEps<- (.Machine$double.eps)^(2/4)
inf<- max(1e+38,sqrt(.Machine$double.xmax))

# extract info from specified variance components
fv <- function(vv){
# vv: list of variance components -- list(AA,DD,HH,AD,MH,EE,...) (Abney 200 pp635)
   if(any(!is.element(c("AA","DD","HH","AD","MH","EE"),names(vv)))){
      cat("Assume components are in correct order AA,DD,HH,AD,MH,EE...\n")
   }else {
      vv[1:6]<- list(AA=vv$AA,DD=vv$DD,HH=vv$HH,AD=vv$AD,MH=vv$MH,EE=vv$EE)
      names(vv)[1:6]<- c("AA","DD","HH","AD","MH","EE")
   }
   if(is.null(vv[[2]])){
      if(!is.null(vv[[3]])){
         vv[3]<- list(NULL)
         cat("v[[3]] is set to null because v[[2]] is null.\n")
      }
      if(!is.null(vv[[5]])){
         vv[5]<- list(NULL)
         cat("v[[5]] is set to null because v[[2]] is null.\n")
      }
   }
   if(is.null(vv[[1]])){
      if(!is.null(vv[[4]])){
        vv[4]<- list(NULL)
        cat("v[[4]] is set to null because v[[1]] is null.\n")
      }
   }
   if(is.null(vv[[3]])){
      if(!is.null(vv[[4]])){
         vv[4]<- list(NULL)
         cat("v[[4]] is set to null because v[[3]] is null.\n")
      }
   }
   nv<- length(vv)
   if(nv>0){
      nnl<- NULL
      for(i in 1:nv) nnl<- c(nnl,!is.null(vv[[i]]))
      nn<- cumsum(nnl)
   }else stop("At least the environmental variance component should be included.\n")

   list(v=vv,nv=nv,nnl=nnl,nn=nn)
}

# one of the main functions, Nelder-Mead method
estVC <-
   function(y,
            x,
            v,
            initpar,
            nit,
            method,
            control,
            hessian)
{
   UseMethod("estVC")
}

estVC.default <-
   function(y,
            x,
            v = vector("list",6),
            initpar,
            nit = 25,
            method = c("Nelder-Mead", "BFGS", "CG", "SANN"),
            control = list(),
            hessian = FALSE)
{
   if(any(is.infinite(y)))
      stop("y: missing or infinite data points not allowed.")
   if(!missing(x))
      if(any(is.infinite(x)))
         stop("x: missing or infinite data points not allowed.")

   if(!missing(x)){
      oTmp<- data.frame(y=y,x)
   }else oTmp<- data.frame(y=y)
   oTmp<- model.frame(y~.,oTmp)
   y<- model.response(oTmp)
   x<- model.matrix(y~.,oTmp)
   method<-  match.arg(method)

   estVC.2(y = y,
            x = x,
            v = v,
            initpar = initpar,
            nit = nit,
            method = method,
            control = control,
            hessian = hessian)
}

estVC.1 <-
   function(y,
            x,
            v,
            initpar,
            nit,
            method,
            control,
            hessian)
{
# estimate all background genetic variance (bgv)
# y: ny by 1 matrix, response
# x: covariates
# v: list of variance components -- list(AA,DD,HH,AD,MH,EE,...) (Abney 200 pp635)
# initpar: initial parameters, will be initilized automatically if missing
# nit: number of iterations to call optim()
# method: the method to be used
# control: A list of control parameters
# hessian: logical. should a numerically differentiated Hessian matrix be returned?
   fs<- control$fnscale
      if(!is.null(fs) && fs<0) fs<- -1 else fs<- 1

   if(!is.null(dim(y))){
      if(length(dim(y))>2) stop("y: may be wrong.\n")
      ny<- nrow(y)
   }else ny<- length(y)
   nb<- ncol(x)
   ov<- fv(v)
   if(missing(initpar)) initpar<- c(rep(mean(y),ncol(x)),rep(var(y),sum(ov$nnl)))

   ppar<- initpar
   vval<- inf
   optfct<- function(par,a=list(nb=nb,ny=ny,ov=ov)){
      b<- par[1:a$nb]
      if(a$ov$nnl[4]){
         tmp<- sqrt(exp(par[a$nb+a$ov$nn[1]]+par[a$nb+a$ov$nn[3]])/2)
         if(abs(par[a$nb+a$ov$nn[4]])>tmp) par[a$nb+a$ov$nn[4]]<- sign(par[a$nb+a$ov$nn[4]])*tmp
      }

      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
         if(a$ov$nnl[i]){
            if(i!=4){
               S<- S + a$ov$v[[i]]*exp(par[a$nb+a$ov$nn[i]])
            }else S<- S + a$ov$v[[i]]*par[a$nb+a$ov$nn[i]]
         }
      }
      if(any(is.infinite(S))) return(fs*(vval + machineEps))

      u<- x%*%b
#      tmp<- determinant(S,logarithm = TRUE)
#      if(!(tmp$sign>0) || abs(tmp$mod)>inf) return(fs*(vval + machineEps))
#      tmp.mod<- tmp$mod
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr)) return(fs*(vval + machineEps))
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(a$ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp)) return(fs*(vval + machineEps))
      if(abs(tmp)>inf) return(fs*(vval + machineEps))

      # necessary because of constraints
      if(tmp < vval){
         ppar<<- par; vval<<- tmp
      }
      fs*tmp
   }

   if(method != "Nelder-Mead"){
      contrTmp<- control
         contrTmp$abstol<- min(contrTmp$abstol,1e-6)
      oo<- optim(ppar,optfct,gr=NULL,a=list(nb=nb,ny=ny,ov=ov),method=method,control=contrTmp)
         oo$par<- ppar; oo$value<- vval
   }
   val1<- val2<- inf
   while(nit>0){
      oo<- optim(ppar,optfct,gr=NULL,a=list(nb=nb,ny=ny,ov=ov),method="Nelder-Mead",control=control,hessian=hessian)
         oo$par<- ppar; oo$value<- vval

      nit<- nit-1

      val1<- val2
      val2<- oo$value
      if(abs(val2-val1)<machineEps && abs(val2-val1)<abs(val1)*machineEps)
         break
   }

   oo$value<- as.numeric(oo$value)

   oo$value<- -oo$value #-fs*oo$value
   attributes(oo$value)<- NULL
   for(i in 1:ov$nv){
      if(ov$nnl[i]){
         if(i!=4) oo$par[nb+ov$nn[i]]<- exp(oo$par[nb+ov$nn[i]])
      }
   }
   names(oo$par)[1:nb]<- colnames(x)
   names(oo$par)[nb+1:sum(ov$nnl)]<- names(ov$v[ov$nnl])
   oo$y<- as.matrix(y)
   oo$x<- as.matrix(x)
   oo$v<- v
   oo$nv<- ov$nv
   oo$nnl<- ov$nnl
   oo$nn<- ov$nn
   class(oo)<- "bgv"
   oo
}

# NOTES: as accurate as but about 2.5 times as fast as estVC.1
estVC.2 <-
   function(y,
            x,
            v,
            initpar,
            nit,
            method,
            control,
            hessian)
{
# estimate all background genetic variance (bgv)
# y: ny by 1 matrix, response
# x: covariates
# v: list of variance components -- list(AA,DD,HH,AD,MH,EE,...) (Abney 200 pp635)
# initpar: initial parameters, will be initilized automatically if missing
# nit: number of iterations to call optim()
# method: the method to be used
# control: A list of control parameters
# hessian: logical. should a numerically differentiated Hessian matrix be returned?
   fs<- control$fnscale
      if(!is.null(fs) && fs<0) fs<- -1 else fs<- 1

   if(!is.null(dim(y))){
      if(length(dim(y))>2) stop("y: may be wrong.\n")
      ny<- nrow(y)
   }else ny<- length(y)
   nb<- ncol(x)
   ov<- fv(v)
   if(missing(initpar)) initpar<- c(rep(mean(y),ncol(x)),rep(var(y),sum(ov$nnl)))

   ppar<- initpar
   vval<- inf
   optfct<- function(par,a=list(nb=nb,ny=ny,ov=ov)){
      b<- par[1:a$nb]
      if(a$ov$nnl[4]){
         tmp<- sqrt(exp(par[a$nb+a$ov$nn[1]]+par[a$nb+a$ov$nn[3]])/2)
         if(abs(par[a$nb+a$ov$nn[4]])>tmp) par[a$nb+a$ov$nn[4]]<- sign(par[a$nb+a$ov$nn[4]])*tmp
      }

      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
         if(a$ov$nnl[i]){
            if(i!=4){
               S<- S + a$ov$v[[i]]*exp(par[a$nb+a$ov$nn[i]])
            }else S<- S + a$ov$v[[i]]*par[a$nb+a$ov$nn[i]]
         }
      }
      if(any(is.infinite(S))) return(fs*(vval + machineEps))

      u<- x%*%b
#      tmp<- determinant(S,logarithm = TRUE)
#      if(!(tmp$sign>0) || abs(tmp$mod)>inf) return(fs*(vval + machineEps))
#      tmp.mod<- tmp$mod
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr)) return(fs*(vval + machineEps))
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(a$ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp)) return(fs*(vval + machineEps))
      if(abs(tmp)>inf) return(fs*(vval + machineEps))

      # necessary because of constraints
      if(tmp < vval){
         ppar<<- par; vval<<- tmp
      }
      fs*tmp
   }
   optfct.b<- function(a=list(par=ppar,nb=nb,ny=ny,ov=ov)){
      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
         if(a$ov$nnl[i]){
            if(i!=4){
               S<- S + a$ov$v[[i]]*exp(a$par[a$nb+a$ov$nn[i]])
            }else S<- S + a$ov$v[[i]]*a$par[a$nb+a$ov$nn[i]]
         }
      }
      if(any(is.infinite(S))) return(a$par[1:a$nb])

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
   optfct.v<- function(par,a=list(par=ppar,nb=nb,ny=ny,ov=ov)){
      b<- a$par[1:a$nb]
      if(a$ov$nnl[4]){
         tmp<- sqrt(exp(par[a$ov$nn[1]]+par[a$ov$nn[3]])/2)
         if(abs(par[a$ov$nn[4]])>tmp) par[a$ov$nn[4]]<- sign(par[a$ov$nn[4]])*tmp
      }

      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
         if(a$ov$nnl[i]){
            if(i!=4){
               S<- S + a$ov$v[[i]]*exp(par[a$ov$nn[i]])
            }else S<- S + a$ov$v[[i]]*par[a$ov$nn[i]]
         }
      }
      if(any(is.infinite(S))) return(fs*(vval + machineEps))

      u<- x%*%b
#      tmp<- determinant(S,logarithm = TRUE)
#      if(!(tmp$sign>0) || abs(tmp$mod)>inf) return(fs*(vval + machineEps))
#      tmp.mod<- tmp$mod
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr)) return(fs*(vval + machineEps))
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(a$ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp)) return(fs*(vval + machineEps))
      if(abs(tmp)>inf) return(fs*(vval + machineEps))

      # necessary because of constraints
      if(tmp < vval){
         ppar<- a$par; ppar[-c(1:a$nb)]<<- par; vval<<- tmp
      }
      fs*tmp
   }

   val1<- val2<- inf
   pparTmp<- ppar
   if(length(ppar) < nb+2){
      while(nit>0){
         oo<- optimize(optfct.v,interval=c(-1000,1000),a=list(par=pparTmp,nb=nb,ny=ny,ov=ov),maximum=ifelse(fs<0,TRUE,FALSE))
            oo<- list(par=ppar[-c(1:nb)], value=vval)
         pparTmp<- ppar; pparTmp[1:nb]<- optfct.b(a=list(par=ppar,nb=nb,ny=ny,ov=ov))

         val1<- val2
         val2<- oo$value
         if(abs(val2-val1)<machineEps && abs(val2-val1)<abs(val1)*machineEps)
            break
         nit<- nit-1
      }
   }else{
      if(method != "Nelder-Mead"){
         while(nit>0){
            contrTmp<- control
               contrTmp$abstol<- min(contrTmp$abstol,1e-6)
            oo<- optim(pparTmp[-c(1:nb)],optfct.v,gr=NULL,a=list(par=pparTmp,nb=nb,ny=ny,ov=ov),method=method,control=contrTmp)
               oo$par<- ppar[-c(1:nb)]; oo$value<- vval
            pparTmp<- ppar; pparTmp[1:nb]<- optfct.b(a=list(par=ppar,nb=nb,ny=ny,ov=ov))

            val1<- val2
            val2<- oo$value
            if(abs(val2-val1)<machineEps)
               break
            nit<- nit-1
         }
      }
      while(nit>0){
         oo<- optim(pparTmp[-c(1:nb)],optfct.v,gr=NULL,a=list(par=pparTmp,nb=nb,ny=ny,ov=ov),method="Nelder-Mead",control=control)
            oo$par<- ppar[-c(1:nb)]; oo$value<- vval
         pparTmp<- ppar; pparTmp[1:nb]<- optfct.b(a=list(par=ppar,nb=nb,ny=ny,ov=ov))

         val1<- val2
         val2<- oo$value
         if(abs(val2-val1)<machineEps && abs(val2-val1)<abs(val1)*machineEps)
            break
         nit<- nit-1
      }
   }
   oo<- optim(ppar,optfct,gr=NULL,a=list(nb=nb,ny=ny,ov=ov),method="Nelder-Mead",control=control,hessian=hessian)
      oo$par<- ppar; oo$value<- vval
   oo$value<- as.numeric(oo$value)

   oo$value<- -oo$value
   attributes(oo$value)<- NULL
   for(i in 1:ov$nv){
      if(ov$nnl[i]){
         if(i!=4) oo$par[nb+ov$nn[i]]<- exp(oo$par[nb+ov$nn[i]])
      }
   }
   names(oo$par)[1:nb]<- colnames(x)
   names(oo$par)[nb+1:sum(ov$nnl)]<- names(ov$v[ov$nnl])
   oo$y<- as.matrix(y)
   oo$x<- as.matrix(x)
   oo$v<- ov$v
   oo$nv<- ov$nv
   oo$nnl<- ov$nnl
   oo$nn<- ov$nn
   class(oo)<- "bgv"
   oo
}

# NOTES: as accurate as but about 5 times as fast as estVC.1
# as accurate as but about 2 times as fast as estVC.2; may Not stable!!!
estVC.3 <-
   function(y,
            x,
            v,
            initpar,
            nit,
            method,
            control,
            hessian)
{
# estimate all background genetic variance (bgv)
# y: ny by 1 matrix, response
# x: desig matrix including overall mean !!!
# v: list of variance components -- list(AA,DD,HH,AD,MH,EE,...) (Abney 200 pp635)
# initpar: initial parameters, will be initilized automatically if missing
# nit: number of iterations
# ...: options passed to nlminb()
   if(!is.null(dim(y))){
      if(length(dim(y))>2) stop("y: may be wrong.\n")
      ny<- nrow(y)
   }else ny<- length(y)
   nb<- ncol(x)
   ov<- fv(v)
   if(missing(initpar)) initpar<- c(rep(mean(y),ncol(x)),rep(var(y),sum(ov$nnl)))

   ppar<- initpar
   vval<- inf
   optfct<- function(par,a=list(nb=nb,ny=ny,ov=ov)){
      b<- par[1:a$nb]
      if(a$ov$nnl[4]){
         tmp<- sqrt(exp(par[a$nb+a$ov$nn[1]]+par[a$nb+a$ov$nn[3]])/2)
         if(abs(par[a$nb+a$ov$nn[4]])>tmp) par[a$nb+a$ov$nn[4]]<- sign(par[a$nb+a$ov$nn[4]])*tmp
      }

      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
         if(a$ov$nnl[i]){
            if(i!=4){
               S<- S + a$ov$v[[i]]*exp(par[a$nb+a$ov$nn[i]])
            }else S<- S + a$ov$v[[i]]*par[a$nb+a$ov$nn[i]]
         }
      }
      if(any(is.infinite(S))) return(vval + machineEps)

      u<- x%*%b
#      tmp<- determinant(S,logarithm = TRUE)
#      if(!(tmp$sign>0) || abs(tmp$mod)>inf) return(vval + machineEps)
#      tmp.mod<- tmp$mod
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr)) return(vval + machineEps)
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(a$ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp)) return(vval + machineEps)
      if(abs(tmp)>inf) return(vval + machineEps)

      # necessary because of constraints
      if(tmp < vval){
         ppar<<- par; vval<<- tmp
      }
      tmp
   }
   optfct.b<- function(a=list(par=par,nb=nb,ny=ny,ov=ov)){
      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
         if(a$ov$nnl[i]){
            if(i!=4){
               S<- S + a$ov$v[[i]]*exp(a$par[a$nb+a$ov$nn[i]])
            }else S<- S + a$ov$v[[i]]*a$par[a$nb+a$ov$nn[i]]
         }
      }
      if(any(is.infinite(S))) return(a$par[1:a$nb])

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
   optfct.v<- function(par,a=list(par=par,nb=nb,ny=ny,ov=ov)){
      b<- a$par[1:a$nb]
      if(a$ov$nnl[4]){
         tmp<- sqrt(exp(par[a$ov$nn[1]]+par[a$ov$nn[3]])/2)
         if(abs(par[a$ov$nn[4]])>tmp) par[a$ov$nn[4]]<- sign(par[a$ov$nn[4]])*tmp
      }

      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
         if(a$ov$nnl[i]){
            if(i!=4){
               S<- S + a$ov$v[[i]]*exp(par[a$ov$nn[i]])
            }else S<- S + a$ov$v[[i]]*par[a$ov$nn[i]]
         }
      }
      if(any(is.infinite(S))) return(vval + machineEps)

      u<- x%*%b
#      tmp<- determinant(S,logarithm = TRUE)
#      if(!(tmp$sign>0) || abs(tmp$mod)>inf) return(vval + machineEps)
#      tmp.mod<- tmp$mod
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr)) return(vval + machineEps)
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(a$ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp)) return(vval + machineEps)
      if(abs(tmp)>inf) return(vval + machineEps)

      # necessary because of constraints
      if(tmp < vval){
         ppar<- a$par; ppar[-c(1:a$nb)]<<- par; vval<<- tmp
      }
      tmp
   }

   val1<- val2<- inf
   pparTmp<- ppar
   while(nit>0){
      oo<- nlminb(pparTmp[-c(1:nb)],optfct.v,a=list(par=pparTmp,nb=nb,ny=ny,ov=ov))
         oo$par<- ppar[-c(1:nb)]; oo$objective<- vval
      pparTmp<- ppar; pparTmp[1:nb]<- optfct.b(a=list(par=ppar,nb=nb,ny=ny,ov=ov))
      val1<- val2
      val2<- oo$objective
      if(abs(val2-val1)<machineEps && abs(val2-val1)<abs(val1)*machineEps)
         break
      nit<- nit-1
   }
   oo<- nlminb(ppar,optfct,a=list(nb=nb,ny=ny,ov=ov))
      oo$par<- ppar
      oo$objective<- vval
      names(oo)[match(c("objective"),names(oo))]<- c("value")
   oo<- optim(ppar,optfct,gr=NULL,a=list(nb=nb,ny=ny,ov=ov),method="Nelder-Mead",control=control,hessian=hessian)
      oo$par<- ppar; oo$value<- vval
   oo$value<- as.numeric(oo$value)

   oo$value<- -oo$value
   attributes(oo$value)<- NULL
   for(i in 1:ov$nv){
      if(ov$nnl[i]){
         if(i!=4) oo$par[nb+ov$nn[i]]<- exp(oo$par[nb+ov$nn[i]])
      }
   }
   names(oo$par)[1:nb]<- colnames(x)
   names(oo$par)[nb+1:sum(ov$nnl)]<- names(ov$v[ov$nnl])
   oo$y<- as.matrix(y)
   oo$x<- as.matrix(x)
   oo$v<- ov$v
   oo$nv<- ov$nv
   oo$nnl<- ov$nnl
   oo$nn<- ov$nn
   class(oo)<- "bgv"
   oo
}

blup<- function(object)
{
   UseMethod("blup")
}

blup.bgv<- function(object){
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

print.bgv<- function(object){
   cat("\nvalue:\n"); print(object$value)
   cat("\nparameters:\n"); print(object$par)
}

