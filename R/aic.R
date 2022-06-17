
machineEps<- (.Machine$double.eps)^(2/4)
inf<- max(1e+38,sqrt(.Machine$double.xmax))

estRR <-
   function(yy, xx, dd)
{
   g<- lm(yy~xx-1, weights=1/dd)
   b<- g$coef
   s2<- sum(g$res^2/dd)/length(g$res)
   loglik<- logLik(g)
   if(loglik > 1e+12) loglik<- -1e+12

   lst<- list(value=loglik,parameters=c(b,s2))
   lst
}

print.aic <- function(x,...)
{
   cat("AIC: ",x$aic,"\n", sep="")   
#   cat("maximum likelihood:",x$model$val,"\n")
   cat("Model parameters:\n"); print(x$model$par)
}

aicVC <- 
   function(y,
            x,
            v = list(E=diag(length(y))),
            initpar,
            k = 2,
            init = 1,
            keep = 1,
            direction = c("forward","backward"),
            nit = 25,
            msg = FALSE,
            control = list(),
            hessian = FALSE)
{
   if(!all(is.finite(y)))
      stop("y: non-numeric or infinite data points not allowed.", call.=FALSE)
   if(!missing(x))
      if(any(sapply(x,is.infinite) | sapply(x,is.na)))
         stop("x: missing or infinite data points not allowed.", call.=FALSE)
   UseMethod("aicVC")
}

aicVC.default <- 
   function(y,
            x,
            v = list(E=diag(length(y))),
            initpar,
            k = 2,
            init = 1,
            keep = 1,
            direction = c("forward","backward"),
            nit = 25,
            msg = FALSE,
            control = list(),
            hessian = FALSE)
{
   direction<- match.arg(direction)
   if(direction=="forward"){
      out<- aicVCf(y = y,
            x = x,
            v = v,
            initpar = initpar,
            k = k,
            init = init,
            nit = nit,
            msg = msg,
            control = control,
            hessian = hessian)
   }else if(direction=="backward"){
      out<- aicVCb(y = y,
            x = x,
            v = v,
            initpar = initpar,
            k = k,
            keep = keep,
            nit = nit,
            msg = msg,
            control = control,
            hessian = hessian)
   }else stop("Only forward or backward is considered.", call.=FALSE)
   out
}

# forward
aicVCf <-
   function(y,
            x,
            v = list(E=diag(length(y))),
            initpar,
            k = 2,
            init = 1,
            nit = 25,
            msg = FALSE,
            control = list(),
            hessian = FALSE)
{
# estimate all background genetic variance (VC)
# y: ny by 1 matrix, response
# x: covariates
# v: list of variance components
#    serve as the largest possible model
# initpar: initial parameters, will be initilized automatically if missing
# k: penalty on each parameter
# init: indicator of the initial model. only E by default
# nit: number of iterations to call optim()
# control: A list of control parameters
# hessian: logical. should a numerically differentiated Hessian matrix be returned?
   fv <- function(vv){
   # vv: list of variance components
      if(is.matrix(vv)){
         vv<- list(A=vv,E=diag(nrow(vv)))
      }

      nv<- length(vv)
      if(nv>0){
         nnl<- NULL
         for(i in 1:nv){
            nnl<- c(nnl,!is.null(vv[[i]]))
            if(!is.null(vv[[i]])){
               if(!all(is.finite(vv[[i]])))
                  stop("Only finite numerical elements are allowed in variance matrices!", call.=FALSE)
            }
         }
         nn<- cumsum(nnl)
      }

      list(v=vv,nv=nv,nnl=nnl,nn=nn)
   }

   if(!is.element("E",names(v))){
      v$E<- diag(length(y))
      init<- c(init,length(v))

      n<- match("E",names(v))
      init<- union(init,n)
      rm(n)
   }

   v0<- vector("list",length(v))
      names(v0)<- names(v)
   for(i in 1:length(init))
      v0[init[i]]<- list(v[[init[i]]])
   ov<- fv(v0)

   ob<- list(aic = NULL,
             model = NULL,
             lik = NULL,
             trace = NULL)

   if(msg)
      cat("\nModel with variance components: ",names(ov$v)[ov$nnl],"\n", sep="")
   oo<- estVC(y=y,x=x,v=ov$v,initpar=initpar,nit=nit,control=control,hessian=hessian)
   ob$aic<- -2*oo$value+length(oo$par)*k
   ob$model<- oo
   ob$lik<- c(lik=oo$value,ob$lik)
   ob$trace<- rbind(ov$nnl,ob$trace)
   if(msg){
      cat("*** initial model:",ov$nnl,"\n")
      cat("*** aic:",ob$aic,"\n\n")
   }

   go<- FALSE
   if(sum(!ov$nnl)>0) go<- TRUE
   while(go){
      go<- FALSE
      oo0<- oo
      ns<- c(1:length(v))[!ov$nnl]
      for(i in 1:sum(!ov$nnl)){
         if(ns[i]==6) next
         v0<- ov$v
         v0[ns[i]]<- list(v[[ns[i]]])
         ov0<- fv(v0)
         if(all(ov$nnl == ov0$nnl)) next
         if(msg)
            cat("Model with variance components: ",names(ov0$v)[ov0$nnl],"\n", sep="")
         o0<- estVC(y=y,x=x,v=ov0$v,nit=nit,control=control,hessian=hessian)
         ob$lik<- c(lik=o0$value,ob$lik)
         ob$trace<- rbind(ov0$nnl,ob$trace)
         if(-2*o0$value+length(o0$par)*k < -2*oo0$value+length(oo0$par)*k - 1e-5){
            oo0<- o0
            vv0<- ov0
         }
         if(msg) cat("AIC:",-2*o0$value+length(o0$par)*k,"\n", sep="")
      }
      if(-2*oo0$value+length(oo0$par)*k < ob$aic - 1e-5){
         oo<- oo0
         ov<- vv0
         ob$aic<- -2*oo$value+length(oo$par)*k
         ob$model<- oo
#         ob$lik<- c(lik=oo$value,ob$lik)
#         ob$trace<- rbind(ov$nnl,ob$trace)
         go<- TRUE
         if(sum(!ov$nnl)==0) go<- FALSE

         if(msg){
            cat("*** selected model:",ov$nnl,"\n")
            cat("*** aic:",ob$aic,"\n\n")
         }
      }
   }

   names(ob$lik)<- paste("model",length(ob$lik):1,sep="")
   rownames(ob$trace)<- paste("model",nrow(ob$trace):1,sep="")
   colnames(ob$trace)<- names(ov$v)
   class(ob)<- "aic"
   ob
}

# backward
aicVCb <-
   function(y,
            x,
            v = list(E=diag(length(y))),
            initpar,
            k = 2,
            keep = 1,
            nit = 25,
            msg = FALSE,
            control = list(),
            hessian = FALSE)
{
# estimate all background genetic variance (VC)
# y: ny by 1 matrix, response
# x: covariates
# v: list of variance components
#    serve as the large possible model
# initpar: initial parameters, will be initilized automatically if missing
# k: penalty on each parameter
# keep: terms to be kept. E is kept (won't drop) by default
# nit: number of iterations to call optim()
# control: A list of control parameters
# hessian: logical. should a numerically differentiated Hessian matrix be returned?
   fv <- function(vv){
   # vv: list of variance components
      if(is.matrix(vv)){
         vv<- list(A=vv,E=diag(nrow(vv)))
      }

      nv<- length(vv)
      if(nv>0){
         nnl<- NULL
         for(i in 1:nv){
            nnl<- c(nnl,!is.null(vv[[i]]))
            if(!is.null(vv[[i]])){
               if(!all(is.finite(vv[[i]])))
                  stop("Only finite numerical elements are allowed in variance matrices!", call.=FALSE)
            }
         }
         nn<- cumsum(nnl)
      }

      list(v=vv,nv=nv,nnl=nnl,nn=nn)
   }

   if(!is.element("E",names(v))){
      v$E<- diag(length(y))
      keep<- c(keep,length(v))

      n<- match("E",names(v))
      keep<- union(keep,n)
      rm(n)
   }

   ov<- fv(v)
   ob<- list(aic = NULL,
             model = NULL,
             lik = NULL,
             trace = NULL)

   if(msg)
      cat("\nModel with variance components: ",names(ov$v)[ov$nnl],"\n", sep="")
   oo<- estVC(y=y,x=x,v=ov$v,initpar=initpar,nit=nit,control=control,hessian=hessian)
   ob$aic<- -2*oo$value+length(oo$par)*k
   ob$model<- oo
   ob$lik<- c(lik=oo$value,ob$lik)
   ob$trace<- rbind(ov$nnl,ob$trace)
   if(msg){
      cat("*** initial model:",ov$nnl,"\n")
      cat("*** aic:",ob$aic,"\n\n")
   }

   kpt<- rep(TRUE,length(v)); kpt[keep]<- FALSE
   go<- FALSE
   if(sum(ov$nnl & kpt)>1) go<- TRUE
   while(go){
      go<- FALSE
      oo0<- oo
      ns<- c(1:length(v))[ov$nnl & kpt]
      for(i in 1:sum(ov$nnl & kpt)){
         v0<- ov$v
         v0[ns[i]]<- list(NULL)
         ov0<- fv(v0)
         if(msg)
            cat("Model with variance components: ",names(ov0$v)[ov0$nnl],"\n", sep="")
         o0<- estVC(y=y,x=x,v=ov0$v,nit=nit,control=control,hessian=hessian)
         ob$lik<- c(lik=o0$value,ob$lik)
         ob$trace<- rbind(ov0$nnl,ob$trace)
         if(-2*o0$value+length(o0$par)*k < -2*oo0$value+length(oo0$par)*k - 1e-5){
            oo0<- o0
            vv0<- ov0
         }
         if(msg) cat("AIC: ",-2*o0$value+length(o0$par)*k,"\n",sep="")
      }
      if(-2*oo0$value+length(oo0$par)*k < ob$aic - 1e-5){
         oo<- oo0
         ov<- vv0
         ob$aic<- -2*oo$value+length(oo$par)*k
         ob$model<- oo
         go<- TRUE
         if(sum(ov$nnl & kpt)<1) go<- FALSE

         if(msg){
            cat("*** selected model: ",ov$nnl,"\n",sep="")
            cat("*** AIC: ",ob$aic,"\n\n",sep="")
         }
      }
   }

   if(!is.null(ob$trace))
      if(sum(ob$trace[1,])>1){
         v0<- ov$v
         v0[1:(length(v0)-1)]<- vector("list",length(v0)-1)
         ov0<- fv(v0)
         o0<- estVC(y=y,x=x,v=ov0$v,nit=nit,control=control,hessian=hessian)
         ob$lik<- c(lik=o0$value,ob$lik)
         ob$trace<- rbind(ov0$nnl,ob$trace)
      }

   names(ob$lik)<- paste("model",length(ob$lik):1,sep="")
   rownames(ob$trace)<- paste("model",nrow(ob$trace):1,sep="")
   colnames(ob$trace)<- names(ov$v)
   class(ob)<- "aic"
   ob
}

blup.aic<- function(object){
# best linear unbiased prediction (BLUP) for all effects
# object: object from aicVC
   blup(object$model)
}

# backward elemination
AICb <-
   function(aic,
            xin,
            yy,
            xx,
            gdat,
            uu,
            dd,
            kk,
            msg = FALSE)
{
   xx<- t(uu)%*%xx
   oo<- NULL
   go<- TRUE
   while(go){
      go<- FALSE
      if(sum(xin)>0){
         nx<- length(xin)
         ii<- 1:nx
         ii<- ii[xin]
         oo1<- oo
         aic1<- aic
         xin1<- xin
         if(msg) cat("b")
         for(i in 1:length(ii)){
            xin0<- xin
            xin0[ii[i]]<- FALSE
            if(sum(xin0)>0){
               oTmp<- data.frame(gdat[,xin0])
               oTmp<- model.frame(~.,oTmp)
               xxx<- model.matrix(~.,oTmp)[,-1]
               xxx<- t(uu)%*%xxx
            }else xxx<- NULL
            xxx<- cbind(xx,xxx)
            oo0<- estRR(yy,xxx,dd)
            aic0<- -2*oo0$value + kk*length(oo0$par)
            if(aic0 < aic1 - 1e-5){
               oo1<- oo0
               aic1<- aic0
               xin1<- xin0
               go<- TRUE
               if(msg) cat(".")
            }
         }
         if(go){
            if(sum(xin1)!=sum(xin)-1) stop("Drop stage in AICb: something wrong.", call.=FALSE)
            oo<- oo1
            aic<- aic1
            xin<- xin1

            if(msg) cat("-")
         }
      }
   }

   list(oo=oo,aic=aic,xin=xin)
}

# forward selection
AICf <-
   function(aic,
            xin,
            yy,
            xx,
            gdat,
            uu,
            dd,
            kk,
            scope,
            msg = FALSE)
{
   xx<- t(uu)%*%xx
   if(missing(scope)) scope<- rep(TRUE,length(xin))
   oo<- NULL
   go<- FALSE
   if(sum(!xin)>0){
      nx<- length(xin)
      ii<- 1:nx
      ii<- ii[!xin]
      oo2<- oo
      aic2<- aic
      xin2<- xin
      if(msg) cat("f")
      for(i in 1:length(ii)){
         if(!scope[ii[i]]) next
         xin0<- xin
         xin0[ii[i]]<- TRUE
         oTmp<- data.frame(gdat[,xin0])
         oTmp<- model.frame(~.,oTmp)
         xxx<- model.matrix(~.,oTmp)[,-1]
            xxx<- t(uu)%*%xxx
            xxx<- cbind(xx,xxx)
         oo0<- estRR(yy,xxx,dd)
         aic0<- -2*oo0$value + kk*length(oo0$par)
         if(aic0 < aic2 - 1e-5){
            oo2<- oo0
            aic2<- aic0
            xin2<- xin0
            go<- TRUE
            if(msg) cat(".")
         }
      }
      if(go){
         if(sum(xin2)!=sum(xin)+1) stop("Move stage in ICf: something wrong.", call.=FALSE)
         oo<- oo2
         aic<- aic2
         xin<- xin2
         if(msg) cat("+")
      }
   }

   list(oo=oo,aic=aic,xin=xin)
}

# move QTL to the best locations
AICmove <-
   function(aic,
            xin,
            yy,
            xx,
            gdat,
            chrIdx,
            uu,
            dd,
            kk,
            msg=FALSE)
{
   xx<- t(uu)%*%xx
   go<- FALSE
   if(sum(xin)>0) go<- TRUE
   oo<- NULL
   while(go){
      go<- FALSE
      nx<- length(xin)
      for(i in 1:sum(xin)){
         yes<- FALSE
         ii<- 1:nx
            ii<- ii[xin] # QTL location index
         mn<- ifelse(i>1,ii[i-1]+1,1)
         mx<- ifelse(i<sum(xin),ii[i+1]-1,nx)
         if(!missing(chrIdx)) scope<- chrIdx==chrIdx[ii[i]] else scope<- rep(TRUE,length(xin))

         xin1<- xin
         oo1<- oo
         aic1<- aic
         if(msg) cat("<")
         for(j in mn:mx){
            if(i==ii[i] || !scope[j]) next
            xin0<- xin
            xin0[ii[i]]<- FALSE
            xin0[j]<- TRUE
            oTmp<- data.frame(gdat[,xin0])
            oTmp<- model.frame(~.,oTmp)
            xxx<- model.matrix(~.,oTmp)[,-1]
               xxx<- t(uu)%*%xxx
               xxx<- cbind(xx,xxx)
            oo0<- estRR(yy,xxx,dd)
            aic0<- -2*oo0$value + kk*length(oo0$par)
            if(aic0 < aic1 - 1e-5){
               oo1<- oo0
               aic1<- aic0
               xin1<- xin0
               go<- TRUE
               if(msg) cat(".")
            }
         }

         if(go){
            oo<- oo1
            aic<- aic1
            xin<- xin1
         }
      }
   }

   list(oo=oo,aic=aic,xin=xin)
}

# multiple QTL model selection
mAIC.default <-
   function(y,
            x,
            gdat,
            vc = NULL,
            chrIdx,
            xin,
            k = 2,
            direction = c("both", "backward", "forward"),
            ext = FALSE,
            msg = FALSE)
{
# y: single trait
# x: covariates
# gdat: n by ? matrix, marker data. Markers in columes!!!
# chrIdx: chromsome index of markers in columns of gdat
# xin: indicator of gdat columns at start
# k: panelty,e.g., 0.05 threshold/number of parameters
   if(!is.null(vc)){
      if(is.element("VC",attr(vc,"class"))){
         nv<- length(vc$v)
         nb<- length(vc$par) - nv
         nr<- nrow(vc$y)
         vcov<- matrix(0,nrow=nr,ncol=nr)
         for(i in 1:nv)
            vcov<- vcov + vc$v[[i]]*vc$par[nb+i]
      }else{
         if(is.data.frame(vc)) vc<- as.matrix(vc)
         if(!is.matrix(vc)) stop("'vc' should be a matrix.", call.=FALSE)
         if(!is.numeric(vc)) stop("'vc' should be a numeric matrix.", call.=FALSE)
         vcov<- vc
      }
   }else vcov<- diag(0,nrow(as.matrix(y)))
#   dd<- svd(vcov); uu<- dd$u; dd<- dd$d
   dd<- eigen(vcov,symmetric=T); uu<- dd$vec; dd<- dd$val
      if(min(dd)<0 && abs(min(dd))>sqrt(.Machine$double.eps)) stop("Variance-covariance may not be positive definite.", call.=FALSE)
   dd<- abs(dd)

   if(is.matrix(y) || is.data.frame(y)){
      if(dim(y)[2]>1)
         cat("   Warning: y: only the first column will be analyzed.\a\n")
      y<- y[,1]
   }
   if(!missing(x)){
      oTmp<- data.frame(y=y,x)
   }else oTmp<- data.frame(y=y)
   oTmp<- model.frame(y~.,oTmp)
   y<- model.response(oTmp)
   x<- model.matrix(y~.,oTmp)
   yy<- t(uu)%*%y
      yy<- as.matrix(yy)

   if(any(is.na(gdat)))
      stop("There are missing genotypes...", call.=FALSE)
   if(is.numeric(gdat)){
      gdat<- matrix(as.character(gdat),nrow=nrow(gdat))
   }
   if(missing(xin)) xin<- rep(F,ncol(gdat))
   if(any(xin)){
      oTmp<- data.frame(gdat[,xin])
      oTmp<- model.frame(~.,oTmp)
      xxx<- model.matrix(~.,oTmp)[,-1]
         xxx<- cbind(x,xxx)
   }else xxx<- x
      xxx<- t(uu)%*%xxx
   oo<- estRR(yy,xxx,dd)
   aic<- -2*oo$value + k*length(oo$par)

   loci.idx<- 1:length(xin)
if(ext){
   if(msg)
      cat("   Remove extraneous QTL...\n")
   tmp<- AICb(aic=aic,xin=xin,yy=yy,xx=x,gdat=gdat,uu=uu,dd=dd,kk=k,msg=msg)
   if(msg) cat("Done.\n")
   if(!is.null(tmp$oo)){
      oo<- tmp$oo
      aic<- tmp$aic
      xin=tmp$xin
   }
   if(msg){
      cat("   There are ",sum(xin)," candidate QTL: ",loci.idx[xin],"\n",sep="")
      cat("   Log-likelihood: ",oo$val,"\n",sep="")
      cat(date(),"\n\n")
   }

   if(msg)
      cat("   Move QTL to the best locations...\n")
   tmp<- AICmove(aic=aic,xin=xin,yy=yy,xx=x,gdat=gdat,chrIdx=chrIdx,uu=uu,dd=dd,kk=k,msg=msg)
   if(msg) cat("Done.\n")
   if(!is.null(tmp$oo)){
      oo<- tmp$oo
      aic<- tmp$aic
      xin=tmp$xin
   }
   if(msg){
      cat("   There are ",sum(xin)," candidate QTL: ",loci.idx[xin],"\n",sep="")
      cat("   Log-likelihood: ",oo$val,"\n",sep="")
      cat(date(),"\n\n")
   }
}

   direction<- match.arg(direction)
   forward<- direction=="both" || direction=="forward"
   backward<- direction=="both" || direction=="backward"
   go<- TRUE
   if(msg){
      if(forward && backward) cat("   Stepwise selection...\n")
      else if(forward) cat("   Forward selection...\n")
      else cat("   Backward selection...\n")
   }
   while(go){
      go<- FALSE
      if(backward){
         tmp<- AICb(aic=aic,xin=xin,yy=yy,xx=x,gdat=gdat,uu=uu,dd=dd,kk=k,msg=msg)
         if(!is.null(tmp$oo)){
            oo<- tmp$oo
            aic<- tmp$aic
            xin=tmp$xin

            if(msg){
               cat("   There are ",sum(xin)," candidate QTL: ",loci.idx[xin],"\n",sep="")
               cat("   Log-likelihood: ",oo$val,"\n",sep="")
               cat("   AIC: ",aic,"\n",sep="")
               cat(date(),"\n\n")
            }
            go<- TRUE
            next
         }
      }

      if(forward){# scope will be missing
         tmp<- AICf(aic=aic,xin=xin,yy=yy,xx=x,gdat=gdat,uu=uu,dd=dd,kk=k,msg=msg)
         if(!is.null(tmp$oo)){
            oo<- tmp$oo
            aic<- tmp$aic
            xin=tmp$xin

            if(msg){
               cat("   There are ",sum(xin)," candidate QTL: ",loci.idx[xin],"\n",sep="")
               cat("   Log-likelihood: ",oo$val,"\n",sep="")
               cat("   AIC: ",aic,"\n",sep="")
               cat(date(),"\n\n")
            }
            go<- TRUE
         }
      }
   }
   if(msg){
      cat("   Done.\n")
      cat("   There are ",sum(xin)," putative QTL.\n",sep="")
      cat("   Log-likelihood: ",oo$val,"\n\n",sep="")
   }

   list(model=oo,aic=aic,snp=colnames(gdat)[xin],xin=xin)
}

# ------------------------------
# Haley-Knott method
#
fxx.hk<- function(prdat,xin){
   xx<- NULL
   if(sum(xin)>0){
      ii<- 1:length(xin)
         ii<- ii[xin]
      for(k in ii){
         xxx<- cbind(prdat$pr[,1,k]-prdat$pr[,3,k],prdat$pr[,2,k])
         xx<- cbind(xx,xxx)
      }
   }
   xx
}

# backward elemination
AICb.HK <-
   function(aic,
            xin,
            yy,
            xx,
            prdat,
            uu,
            dd,
            kk,
            msg = FALSE)
{
   xx<- t(uu)%*%xx
   oo<- NULL
   go<- TRUE
   while(go){
      go<- FALSE
      if(sum(xin)>0){
         nx<- length(xin)
         ii<- 1:nx
         ii<- ii[xin]
         oo1<- oo
         aic1<- aic
         xin1<- xin
         if(msg) cat("b")
         for(i in 1:length(ii)){
            xin0<- xin
            xin0[ii[i]]<- FALSE
            if(sum(xin0)>0){
               xxx<- fxx.hk(prdat,xin0)
               xxx<- t(uu)%*%xxx
            }else xxx<- NULL
            xxx<- cbind(xx,xxx)
            oo0<- estRR(yy,xxx,dd)
            aic0<- -2*oo0$value + kk*length(oo0$par)
            if(aic0 < aic1 - 1e-5){
               oo1<- oo0
               aic1<- aic0
               xin1<- xin0
               go<- TRUE
               if(msg) cat(".")
            }
         }
         if(go){
            if(sum(xin1)!=sum(xin)-1)
               stop("Drop stage in AICb: something wrong.", call.=FALSE)
            oo<- oo1
            aic<- aic1
            xin<- xin1

            if(msg) cat("-")
         }
      }
   }

   list(oo=oo,aic=aic,xin=xin)
}

# forward selection
AICf.HK <-
   function(aic,
            xin,
            yy,
            xx,
            prdat,
            uu,
            dd,
            kk,
            scope,
            msg = FALSE)
{
   xx<- t(uu)%*%xx
   if(missing(scope)) scope<- rep(TRUE,length(xin))
   oo<- NULL
   go<- FALSE
   if(sum(!xin)>0){
      nx<- length(xin)
      ii<- 1:nx
      ii<- ii[!xin]
      oo2<- oo
      aic2<- aic
      xin2<- xin
      if(msg) cat("f")
      for(i in 1:length(ii)){
         if(!scope[ii[i]]) next
         xin0<- xin
         xin0[ii[i]]<- TRUE
         xxx<- fxx.hk(prdat,xin0)
            xxx<- t(uu)%*%xxx
            xxx<- cbind(xx,xxx)
         oo0<- estRR(yy,xxx,dd)
         aic0<- -2*oo0$value + kk*length(oo0$par)
         if(aic0 < aic2 - 1e-5){
            oo2<- oo0
            aic2<- aic0
            xin2<- xin0
            go<- TRUE
            if(msg) cat(".")
         }
      }
      if(go){
         if(sum(xin2)!=sum(xin)+1) stop("Move stage in ICf: something wrong.", call.=FALSE)
         oo<- oo2
         aic<- aic2
         xin<- xin2
         if(msg) cat("+")
      }
   }

   list(oo=oo,aic=aic,xin=xin)
}

# move QTL to the best locations
AICmove.HK <-
   function(aic,
            xin,
            yy,
            xx,
            prdat,
            uu,
            dd,
            kk,
            msg = FALSE)
{
   xx<- t(uu)%*%xx
   go<- FALSE
   if(sum(xin)>0) go<- TRUE
   oo<- NULL
   while(go){
      go<- FALSE
      nx<- length(xin)
      for(i in 1:sum(xin)){
         yes<- FALSE
         ii<- 1:nx
            ii<- ii[xin] # QTL location index
         mn<- ifelse(i>1,ii[i-1]+1,1)
         mx<- ifelse(i<sum(xin),ii[i+1]-1,nx)
         scope<- prdat$chr==prdat$chr[ii[i]]

         xin1<- xin
         oo1<- oo
         aic1<- aic
         if(msg) cat("<")
         for(j in mn:mx){
            if(i==ii[i] || !scope[j]) next
            xin0<- xin
            xin0[ii[i]]<- FALSE
            xin0[j]<- TRUE
            xxx<- fxx.hk(prdat,xin0)
               xxx<- t(uu)%*%xxx
               xxx<- cbind(xx,xxx)
            oo0<- estRR(yy,xxx,dd)
            aic0<- -2*oo0$value + kk*length(oo0$par)
            if(aic0 < aic1 - 1e-5){
               oo1<- oo0
               aic1<- aic0
               xin1<- xin0
               go<- TRUE
               if(msg) cat(".")
            }
         }

         if(go){
            oo<- oo1
            aic<- aic1
            xin<- xin1
         }
      }
   }

   list(oo=oo,aic=aic,xin=xin)
}

# multiple QTL model selection
mAIC.HK <-
   function(y,
            x,
            prdat,
            vc = NULL,
            xin,
            k = 2,
            direction = c("both", "backward", "forward"),
            ext = FALSE,
            msg = FALSE)
{
# y: single trait
# x: covariates
# prdat: prdat$pr: n by 3 by ? matrix, conditional probabilities
# xin: indicator of prdat$pr columns at start
# k: panelty,e.g., 0.05 threshold/number of parameters
   if(!is.null(vc)){
      if(is.element("VC",attr(vc,"class"))){
         nv<- length(vc$v)
         nb<- length(vc$par) - nv
         nr<- nrow(vc$y)
         vcov<- matrix(0,nrow=nr,ncol=nr)
         for(i in 1:nv)
            vcov<- vcov + vc$v[[i]]*vc$par[nb+i]
      }else{
         if(is.data.frame(vc)) vc<- as.matrix(vc)
         if(!is.matrix(vc)) stop("'vc' should be a matrix.", call.=FALSE)
         if(!is.numeric(vc)) stop("'vc' should be a numeric matrix.", call.=FALSE)
         vcov<- vc
      }
   }else vcov<- diag(0,nrow(as.matrix(y)))
#   dd<- svd(vcov); uu<- dd$u; dd<- dd$d
   dd<- eigen(vcov,symmetric=T); uu<- dd$vec; dd<- dd$val
      if(min(dd)<0 && abs(min(dd))>sqrt(.Machine$double.eps)) stop("Variance-covariance: may not be positive definite.", call.=FALSE)
   dd<- abs(dd)

   if(is.matrix(y) || is.data.frame(y)){
      if(dim(y)[2]>1)
         cat("   Warning: y: only the first column will be analyzed.\a\n")
      y<- y[,1]
   }
   if(!missing(x)){
      oTmp<- data.frame(y=y,x)
   }else oTmp<- data.frame(y=y)
   oTmp<- model.frame(y~.,oTmp)
   y<- model.response(oTmp)
   x<- model.matrix(y~.,oTmp)
   yy<- t(uu)%*%y
      yy<- as.matrix(yy)

   if(missing(xin)) xin<- rep(F,dim(prdat$pr)[3])
   if(any(xin)){
      xxx<- fxx.hk(prdat,xin)
         xxx<- cbind(x,xxx)
   }else xxx<- x
   xxx<- t(uu)%*%xxx
   oo<- estRR(yy,xxx,dd)
   aic<- -2*oo$value + k*length(oo$par)

   loci.idx<- 1:length(xin)
if(ext){
   if(msg)
      cat("   Remove extraneous QTL...\n")
   tmp<- AICb.HK(aic=aic,xin=xin,yy=yy,xx=x,prdat=prdat,uu=uu,dd=dd,kk=k,msg=msg)
   if(msg) cat("   Done.\n")
   if(!is.null(tmp$oo)){
      oo<- tmp$oo
      aic<- tmp$aic
      xin=tmp$xin
   }
   if(msg){
      cat("   There are ",sum(xin)," candidate QTL: ",loci.idx[xin],"\n",sep="")
      cat("   Log-likelihood: ",oo$val,"\n",sep="")
      cat(date(),"\n\n")
   }

   if(msg)
      cat("   Move QTL to the best locations...\n")
   tmp<- AICmove.HK(aic=aic,xin=xin,yy=yy,xx=x,prdat=prdat,uu=uu,dd=dd,kk=k,msg=msg)
   if(msg) cat("Done.\n")
   if(!is.null(tmp$oo)){
      oo<- tmp$oo
      aic<- tmp$aic
      xin=tmp$xin
   }
   if(msg){
      cat("   There are ",sum(xin)," candidate QTL: ",loci.idx[xin],"\n",sep="")
      cat("   Log-likelihood: ",oo$val,"\n")
      cat(date(),"\n\n")
   }
}

   direction<- match.arg(direction)
   forward<- direction=="both" || direction=="forward"
   backward<- direction=="both" || direction=="backward"
   go<- TRUE
   if(msg){
      if(forward && backward) cat("   Stepwise selection...\n")
      else if(forward) cat("   Forward selection...\n")
      else cat("   Backward selection...\n")
   }
   while(go){
      go<- FALSE
      if(backward){
         tmp<- AICb.HK(aic=aic,xin=xin,yy=yy,xx=x,prdat=prdat,uu=uu,dd=dd,kk=k,msg=msg)
         if(!is.null(tmp$oo)){
            oo<- tmp$oo
            aic<- tmp$aic
            xin=tmp$xin

            if(msg){
               cat("   There are ",sum(xin)," candidate QTL: ",loci.idx[xin],"\n",sep="")
               cat("   Log-likelihood: ",oo$val,"\n",sep="")
               cat("   AIC: ",aic,"\n",sep="")
               cat(date(),"\n\n")
            }
            go<- TRUE
            next
         }
      }

      if(forward){# scope will be missing
         tmp<- AICf.HK(aic=aic,xin=xin,yy=yy,xx=x,prdat=prdat,uu=uu,dd=dd,kk=k,msg=msg)
         if(!is.null(tmp$oo)){
            oo<- tmp$oo
            aic<- tmp$aic
            xin=tmp$xin

            if(msg){
               cat("   There are ",sum(xin)," candidate QTL: ",loci.idx[xin],"\n",sep="")
               cat("   Log-likelihood: ",oo$val,"\n",sep="")
               cat("   AIC: ",aic,"\n",sep="")
               cat(date(),"\n\n")
            }
            go<- TRUE
         }
      }
   }
   if(msg){
      cat("   Done.\n")
      cat("   There are ",sum(xin)," putative QTL.\n",sep="")
      cat("   Log-likelihood: ",oo$val,"\n\n",sep="")
   }

   list(model=oo,aic=aic,snp=prdat$snp[xin],xin=xin)
}

mAIC<-
   function(y,
            x,
            gdat,
            prdat = NULL,
            vc = NULL,
            chrIdx,
            xin,
            k = 2,
            direction = c("both", "backward", "forward"),
            ext = FALSE,
            msg = FALSE)
{
   if(!all(is.finite(y)))
      stop("y: non-numeric or infinite data points not allowed.", call.=FALSE)
   if(!missing(x))
      if(any(sapply(x,is.infinite) | sapply(x,is.na)))
         stop("x: missing or infinite data points not allowed.", call.=FALSE)
   if(is.null(prdat) || !is.element("Pr",class(prdat))){
      mAIC.default(y = y,
                   x = x,
                   gdat = gdat,
                   vc = vc,
                   chrIdx = chrIdx,
                   xin = xin,
                   k = k,
                   direction = direction,
                   ext = ext,
                   msg = msg)
   }else{
      mAIC.HK(y = y,
              x = x,
              prdat = prdat,
              vc = vc,
              xin = xin,
              k = k,
              direction = direction,
              ext = ext,
              msg = msg)
   }
}

