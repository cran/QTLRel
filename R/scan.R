
machineEps<- (.Machine$double.eps)^(2/4)
inf<- max(1e+38,sqrt(.Machine$double.xmax))

scanOne<- 
   function(y,
            x,
            gdat,
            prdat = NULL,
            vc = NULL,
            intcovar = NULL,
            minorGenoFreq = 0,
            rmv = TRUE,
            nit = 25)
{
   UseMethod("scanOne")
}

scanOne.default<- 
   function(y,
            x,
            gdat,
            prdat = NULL,
            vc = NULL,
            intcovar = NULL,
            minorGenoFreq = 0,
            rmv = TRUE,
            nit = 25)
{
   if(!is.null(vc)){
      if(any(is.infinite(y)))
         stop("y: missing or infinite data points not allowed.")
      if(!missing(x))
          if(any(is.infinite(x)))
             stop("x: missing or infinite data points not allowed.")

      if(is.element("bgv",attr(vc,"class"))){
         nb<- length(vc$par) - sum(vc$nnl)
         nr<- nrow(vc$y)
         gcov<- matrix(0,nrow=nr,ncol=nr)
         for(i in 1:vc$nv)
            if(vc$nnl[i] && names(vc$v[i])!="EE") gcov<- gcov + vc$v[[i]]*vc$par[nb+vc$nn[i]]
      }else{
         if(is.data.frame(vc)) vc<- as.matrix(vc)
         if(!is.matrix(vc)) stop("vc should be a matrix.")
         if(!is.numeric(vc)) stop("vc should be a numeric matrix.")
         gcov<- vc
      }

      if(!is.null(prdat)){
         scanOne.llrHK(y = y,
            x = x,
            gcov = gcov,
            prdat = prdat,
            intcovar = intcovar,
            nit = nit)
      }else{
         scanOne.llr(y = y,
            x = x,
            gcov = gcov,
            gdat = gdat,
            intcovar = intcovar,
            minorGenoFreq = minorGenoFreq,
            rmv = rmv,
            nit = nit)
      }
   }else{
      if(!missing(x)){
         dat<- data.frame(y=y,x)
      }else dat<- data.frame(y=y)
      llrReg(dat=dat,
             gdat=gdat,
             intcovar=intcovar,
             minorGenoFreq=minorGenoFreq,
             rmv=rmv)
   }
}

print.llr <-
   function(x,...)
{
   tt<- x
   tt<- data.frame(snp=tt$snp,lr=tt$lr)
   if(nrow(tt)>5){
      tt<- tt[1:5,]
      print(tt)
      cat("   ...   ...\n")
   }else print(tt)
}

print.llrHK <-
   function(x,...)
{
   tt<- x
   tt<- data.frame(snp=tt$snp,chr=tt$chr,dist=tt$dist,lr=tt$lr)
   if(nrow(tt)>5){
      tt<- tt[1:5,]
      print(tt)
      cat("   ...   ...\n")
   }else print(tt)
}

# maximum likelihood method
scanOne.llr <-
   function(y,
            x,
            gcov,
            gdat,
            intcovar = NULL,
            minorGenoFreq = 0,
            rmv = TRUE,
            nit = 25)
{
# x: covariates
# gcov: n by n genetic variance-covariance matrix
# gdat: n by ? matrix, marker data. Markers in columes!!!
# intcover: covariates that interact with QTL
# minorGenoFreq: minor genotype frequency at a SNP
# rmv: TRUE -- remove SNPs where minorGenoFreq is smaller, FALSE -- return(NULL)
   tb<- sort(union(as.matrix(gdat),NULL))
   tbf<- NULL
   for(ii in tb) tbf<- rbind(tbf,colSums(gdat==ii))
      if(sum(tbf)!=nrow(gdat)*ncol(gdat)) stop("Error occurred.\n")
   tbf<- apply(tbf,2,min)
   idx<- (tbf < nrow(gdat)*minorGenoFreq)
   if(sum(idx)>0){
      if(rmv){
         gdat<- as.matrix(gdat)
         tb<- sort(union(as.matrix(gdat),NULL))
         tbf<- NULL
         for(ii in tb) tbf<- rbind(tbf,colSums(gdat==ii))
            if(sum(tbf)!=nrow(gdat)*ncol(gdat)) stop("Error occurred.\n")
         tbf<- apply(tbf,2,min)
         idx<- (tbf < nrow(gdat)*minorGenoFreq)
         gdat<- gdat[,!idx]
         rm(tb,tbf,ii,idx)
      }else{
         cat("minor genotype frequency is too small at one or more SNPs.\n")
         return(NULL)
     }
   }

#   dd<- svd(gcov); uu<- dd$u; dd<- dd$d
   dd<- eigen(gcov,symmetric=T); uu<- dd$vec; dd<- dd$val
      if(any(dd < -machineEps)) stop("gcov: may not be positive definite.")
   dd[dd<machineEps]<- machineEps

   if(!missing(x)){
      if(!is.null(intcovar)){
         oTmp<- data.frame(x,intcovar,y=y)
      }else oTmp<- data.frame(x,y=y)
   }else{
      if(!is.null(intcovar)){
         oTmp<- data.frame(intcovar,y=y)
      }else oTmp<- data.frame(y=y)
   }
   oTmp<- model.frame(y~.,oTmp)
   yy<- t(uu)%*%model.response(oTmp)
      yy<- as.matrix(yy)
   xx<- t(uu)%*%model.matrix(y~.,oTmp)
      xx<- as.matrix(xx)
   nb<- ncol(xx)
   oo<- estRR(yy,xx,dd,nit=nit)
   llk0<- oo$value

   initpar<- oo$par
   nsnp<- ncol(gdat)
   llk<- rep(Inf,nsnp)
      names(llk)<- colnames(gdat)
   model.par<- vector("list",nsnp)
      names(model.par)<- colnames(gdat)
   if(!is.null(intcovar)) nint<- ncol(as.matrix(intcovar))
   for(j in 1:nsnp){
#      cat(j,"/",nsnp,"\r")
      if(!missing(x)){
         if(!is.null(intcovar)){
            oTmp<- data.frame(x,intcovar,snp=as.factor(gdat[,j]))
         }else oTmp<- data.frame(x,snp=as.factor(gdat[,j]))
      }else{
         if(!is.null(intcovar)){
            oTmp<- data.frame(intcovar,snp=as.factor(gdat[,j]))
         }else oTmp<- data.frame(snp=as.factor(gdat[,j]))
      }
      oTmp<- model.frame(~.,oTmp)

      if(!is.null(intcovar)){
         nc<- ncol(oTmp)
         str<- paste(colnames(oTmp)[nc-(nint:1)],colnames(oTmp)[nc],collapse="+",sep=":")
         str<- paste("~.+",str,sep="")
      }else{
         str<- paste("~.",sep="")
      }
      xxx<- model.matrix(formula(str),oTmp)
      vns<- c(colnames(xxx),"g","e")
         xxx<- t(uu)%*%xxx

      initp<- c(initpar[1:ncol(xx)],rep(0,ncol(xxx)-ncol(xx)),initpar[1:2+ncol(xx)])
      oo<- estRR(yy,xxx,dd,initpar=initp,nit=nit)

      llk[j]<- 2*(oo$value - llk0)
      model.par[[j]]<- oo$par; names(model.par[[j]])<- vns
   }

   out<- list(snp=colnames(gdat),lr=llk,parameters=model.par)
   class(out)<- c("scanOne","llr")
   out
}


# maximum likelihood method: Haley-Knott method
scanOne.llrHK <-
   function(y,
            x,
            gcov,
            prdat,
            intcovar = NULL,
            nit = 25)
{
# x: covariates
# gcov: n by n genetic variance-covariance matrix
# prdat$pr: n by 3 by ? matrix, conditional probabilities
   dd<- eigen(gcov,symmetric=T); uu<- dd$vec; dd<- dd$val
      if(any(dd < -machineEps)) stop("gcov: may not be positive definite.")
   dd[dd<machineEps]<- machineEps

   if(!missing(x)){
      if(!is.null(intcovar)){
         oTmp<- data.frame(x,intcovar,y=y)
      }else oTmp<- data.frame(x,y=y)
   }else{
      if(!is.null(intcovar)){
         oTmp<- data.frame(intcovar,y=y)
      }else oTmp<- data.frame(y=y)
   }
   oTmp<- model.frame(y~.,oTmp)
   nb<- ncol(model.matrix(y~.,oTmp))
   yy<- t(uu)%*%model.response(oTmp)
      yy<- as.matrix(yy)
   xx<- t(uu)%*%model.matrix(y~.,oTmp)
      xx<- as.matrix(xx)

   oo<- estRR(yy,xx,dd,nit=nit)
   llk0<- oo$value

   initpar<- oo$par
   nsnp<- dim(prdat$pr)[3]
   llk<- rep(Inf,nsnp)
   model.par<- vector("list",nsnp)
      names(model.par)<- prdat$snp
   if(!is.null(intcovar)) nint<- ncol(as.matrix(intcovar))
   for(k in 1:nsnp){
      if(!missing(x)){
         if(!is.null(intcovar)){
            oTmp<- data.frame(x,intcovar,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
         }else oTmp<- data.frame(x,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
      }else{
         if(!is.null(intcovar)){
            oTmp<- data.frame(intcovar,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
         }else oTmp<- data.frame(a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
      }
      oTmp<- model.frame(~.,oTmp)

      if(!is.null(intcovar)){
         nc<- ncol(oTmp)
         str<- paste(paste("(",paste(colnames(oTmp)[nc-1-(nint:1)],collapse="+"),")",sep=""),
                     paste("(",paste(colnames(oTmp)[(nc-1):nc],collapse="+"),")",sep=""),
                     sep=":")
         str<- paste("~.+",str,sep="")
      }else{
         str<- paste("~.",sep="")
      }
      xxx<- model.matrix(formula(str),oTmp)
      vns<- c(colnames(xxx),"g","e")
         xxx<- t(uu)%*%xxx

      initp<- c(initpar[1:ncol(xx)],rep(0,ncol(xxx)-ncol(xx)),initpar[1:2+ncol(xx)])
      oo<- estRR(yy,xxx,dd,initpar=initp,nit=nit)

      llk[k]<- 2*(oo$value - llk0)
      model.par[[k]]<- oo$par; names(model.par[[k]])<- vns
#      cat(k,"/",nsnp,"\r")
   }

   out<- list(snp=prdat$snp,chr=prdat$chr,dist=prdat$dist,lr=llk,parameters=model.par)
   class(out)<- c("scanOne","llrHK")
   out
}

# regression method
llrReg<- 
   function(dat,
            gdat,
            intcovar=NULL,
            minorGenoFreq=0,
            rmv=TRUE){
# gdat: n by ? matrix, marker data. Markers in columes!!!
   tb<- sort(union(as.matrix(gdat),NULL))
   tbf<- NULL
   for(ii in tb) tbf<- rbind(tbf,colSums(gdat==ii))
      if(sum(tbf)!=nrow(gdat)*ncol(gdat)) stop("Error occurred.\n")
   tbf<- apply(tbf,2,min)
   idx<- (tbf < nrow(gdat)*minorGenoFreq)
   if(sum(idx)>0){
      if(rmv){
         gdat<- as.matrix(gdat)
         tb<- sort(union(as.matrix(gdat),NULL))
         tbf<- NULL
         for(ii in tb) tbf<- rbind(tbf,colSums(gdat==ii))
            if(sum(tbf)!=nrow(gdat)*ncol(gdat)) stop("Error occurred.\n")
         tbf<- apply(tbf,2,min)
         idx<- (tbf < nrow(gdat)*minorGenoFreq)
         gdat<- gdat[,!idx]
         rm(tb,tbf,ii,idx)
      }else{
         cat("minor genotype frequency is too small at one or more SNPs.\n")
         return(NULL)
     }
   }

   if(!is.data.frame(dat) || !is.element("y",colnames(dat)))
      stop("dat: should be a data frame include variable 'y'.")
   nsnp<- ncol(gdat)
   llk<- rep(-Inf,nsnp)
      names(llk)<- colnames(gdat)
   model.par<- vector("list",nsnp)
      names(model.par)<- colnames(gdat)
   if(!is.null(intcovar)){
      ncol.dat<- ncol(dat)
      datTmp<- data.frame(dat,intcovar)
         str<- paste("snp",colnames(datTmp)[-c(1:ncol.dat)],sep=":",collapse="+")
      for(j in 1:nsnp){
         tmpdat<- data.frame(datTmp,snp=as.factor(gdat[,j]))
         tmp<- lm(paste("y~.",str,sep="+"),data=tmpdat)
         llk[j]<- logLik(tmp)
         model.par[[j]]<- tmp$coef
      }
   }else{
      for(j in 1:nsnp){
         tmpdat<- data.frame(dat,snp=as.factor(gdat[,j]))
         tmp<- lm(y~.,data=tmpdat)
         llk[j]<- logLik(tmp)
         model.par[[j]]<- tmp$coef
      }
   }

   tmp<- lm(y~.,data=dat)
   llk<- 2*(llk-logLik(tmp))

   out<- list(snp=colnames(gdat),lr= llk,parameters=model.par)
   class(out)<- c("scanOne","llr")
   out
}


# function that is used by llr
estRR <-
   function(yy,
            xx,
            dd,
            initpar,
            nit = 25)
{
   UseMethod("estRR")
}

estRR.default <-
   function(yy,
            xx,
            dd,
            initpar,
            nit = 25)
{
   estRR.3(yy,
           xx,
           dd,
           initpar,
           nit)
}

estRR.1 <-
   function(yy,
            xx,
            dd,
            initpar,
            nit = 25)
{
# yy: ny by 1 matrix, response
# xx: desig matrix including overall mean !!!
# dd: diagnal matrix representing a random variance component
# initpar: initial parameters, will be initilized automatically if missing
# nit: number of iterations to call optim()
   ny<- nrow(yy)
   nb<- ncol(xx)
   if(missing(initpar)) initpar<- c(rep(mean(yy),ncol(xx)),rep(var(yy)*10,2))

   optfct<- function(par){
      b<- par[1:nb]
      sR<- par[nb+1]
      sE<- par[nb+2]

      rr<- yy - xx%*%b
      ds<- dd*exp(sR) + exp(sE)
      if(any(!is.finite(ds))) return(inf)
      if(any(ds<=0)) return(inf)
      tmp<- log(2*pi)*(ny/2) + sum(log(ds))/2 + t(rr)%*%((1/ds)*rr)/2
      # sum is better than prod!!!
      if(!is.finite(tmp)) return(inf)

      tmp
   }
   optfct.b<- function(par){
      sR<- par[nb+1]
      sE<- par[nb+2]
      ds<- dd*exp(sR) + exp(sE)
      b<- lm.wfit(x=xx,y=yy,w=1/ds)$coef
         b[!is.finite(b)]<- 0
      b
   }
   optfct.v<- function(par,a=list(par=initpar,nb=nb)){
      b<- a$par[1:a$nb]
      sR<- par[1]
      sE<- par[2]

      rr<- yy - xx%*%b
      ds<- dd*exp(sR) + exp(sE)
      if(any(!is.finite(ds))) return(inf)
      if(any(ds<=0)) return(inf)
      tmp<- log(2*pi)*(ny/2) + sum(log(ds))/2 + t(rr)%*%((1/ds)*rr)/2
      # sum is better than prod!!!
      if(!is.finite(tmp)) return(inf)

      tmp
   }

   val1<- val2<- inf
   pparTmp<- initpar
   while(nit>0){
      pparTmp[1:nb]<- optfct.b(par=pparTmp)
      oo<- optim(pparTmp[-c(1:nb)],optfct.v,a=list(par=pparTmp,nb=nb),method="Nelder-Mead")
      pparTmp[-c(1:nb)]<- oo$par

      val1<- val2
      val2<- oo$val
      if(abs(val2-val1)<machineEps && abs(val2-val1)<abs(val1)*machineEps) break
      nit<- nit-1
   }

   oo<- optim(pparTmp,optfct,method="Nelder-Mead")

   oo$par[nb+1:2]<- exp(oo$par[nb+1:2])
   oo$value<- -oo$value

   oo
}

# NOTE: nlm can lead to warnings "NA/Inf replaced by maximum positive value"
# from the source src/main/optimize.c
estRR.2 <-
   function(yy,
            xx,
            dd,
            initpar,
            nit=25)
{
# yy: ny by 1 matrix, response
# xx: desig matrix including overall mean !!!
# dd: diagnal matrix representing a random variance component
# initpar: initial parameters, will be initilized automatically if missing
# nit: number of iterations
   ny<- nrow(yy)
   nb<- ncol(xx)
   if(missing(initpar)) initpar<- c(rep(mean(yy),ncol(xx)),rep(var(yy)*10,2))

   optfct<- function(par){
      b<- par[1:nb]
      sR<- par[nb+1]
      sE<- par[nb+2]

      rr<- yy - xx%*%b
      ds<- dd*exp(sR) + exp(sE)
      if(any(!is.finite(ds))) return(inf)
      if(any(ds<=0)) return(inf)
      tmp<- log(2*pi)*(ny/2) + sum(log(ds))/2 + t(rr)%*%((1/ds)*rr)/2
      # sum is better than prod!!! or warnings
      if(!is.finite(tmp)) return(inf)

      tmp
   }
   optfct.b<- function(par){
      sR<- par[nb+1]
      sE<- par[nb+2]
      ds<- dd*exp(sR) + exp(sE)
      b<- lm.wfit(x=xx,y=yy,w=1/ds)$coef
         b[!is.finite(b)]<- 0
      b
   }
   optfct.v<- function(par,a=list(par=initpar,nb=nb)){
      b<- a$par[1:a$nb]
      sR<- par[1]
      sE<- par[2]

      rr<- yy - xx%*%b
      ds<- dd*exp(sR) + exp(sE)
      if(any(!is.finite(ds))) return(inf)
      if(any(ds<=0)) return(inf)
      tmp<- log(2*pi)*(ny/2) + sum(log(ds))/2 + t(rr)%*%((1/ds)*rr)/2
      # sum is better than prod!!! or warnings
      if(!is.finite(tmp)) return(inf)

      tmp
   }

   val1<- val2<- inf
   pparTmp<- initpar
   while(nit>0){
      pparTmp[1:nb]<- optfct.b(par=pparTmp)
      oo<- nlm(optfct.v,pparTmp[-c(1:nb)],a=list(par=pparTmp,nb=nb))
         oo$par<- oo$estimate; oo$value<- oo$minimum
      pparTmp[-c(1:nb)]<- oo$par

      val1<- val2
      val2<- oo$value
      if(abs(val2-val1)<machineEps && abs(val2-val1)<abs(val1)*machineEps)
         break
      nit<- nit-1
   }

   oo<- nlm(optfct,pparTmp)
      oo$par<- oo$estimate; oo$value<- oo$minimum
      oo$par[nb + 1:2]<- exp(oo$par[nb + 1:2])

   oo$value<- -oo$value

   oo
}

estRR.3 <-
   function(yy,
            xx,
            dd,
            initpar,
            nit=25)
{
# yy: ny by 1 matrix, response
# xx: desig matrix including overall mean !!!
# dd: diagnal matrix representing a random variance component
# initpar: initial parameters, will be initilized automatically if missing
# nit: number of iterations
# ...: other options passed to nlm()
   ny<- nrow(yy)
   nb<- ncol(xx)
   if(missing(initpar)) initpar<- c(rep(mean(yy),ncol(xx)),rep(var(yy)*10,2))

   optfct<- function(par){
      b<- par[1:nb]
      sR<- par[nb+1]
      sE<- par[nb+2]

      rr<- yy - xx%*%b
      ds<- dd*exp(sR) + exp(sE)
      if(any(!is.finite(ds))) return(inf)
      if(any(ds<=0)) return(inf)
      tmp<- log(2*pi)*(ny/2) + sum(log(ds))/2 + t(rr)%*%((1/ds)*rr)/2
      # sum is better than prod!!! or warnings
      if(!is.finite(tmp)) return(inf)

      tmp
   }
   optfct.b<- function(par){
      sR<- par[nb+1]
      sE<- par[nb+2]
      ds<- dd*exp(sR) + exp(sE)
      b<- lm.wfit(x=xx,y=yy,w=1/ds)$coef
         b[!is.finite(b)]<- 0
      b
   }
   optfct.v<- function(par,a=list(par=initpar,nb=nb)){
      b<- a$par[1:a$nb]
      sR<- par[1]
      sE<- par[2]

      rr<- yy - xx%*%b
      ds<- dd*exp(sR) + exp(sE)
      if(any(!is.finite(ds))) return(inf)
      if(any(ds<=0)) return(inf)
      tmp<- log(2*pi)*(ny/2) + sum(log(ds))/2 + t(rr)%*%((1/ds)*rr)/2
      # sum is better than prod!!! or warnings
      if(!is.finite(tmp)) return(inf)

      tmp
   }

   upper<- rep(25,2)
      upper<- c(rep(Inf,nb),upper)
   val1<- val2<- inf
   pparTmp<- initpar
   while(nit>0){
      pparTmp[1:nb]<- optfct.b(par=pparTmp)
      oo<- nlminb(pparTmp[-c(1:nb)],optfct.v,a=list(par=pparTmp,nb=nb),upper=upper[-c(1:nb)])
         oo$value<- oo$objective
      pparTmp[-c(1:nb)]<- oo$par

      val1<- val2
      val2<- oo$value
      if(abs(val2-val1)<machineEps && abs(val2-val1)<abs(val1)*machineEps)
         break
      nit<- nit-1
   }

   oo<- nlminb(pparTmp,optfct,upper=upper)
      oo$value<- oo$objective
      oo$par[nb + 1:2]<- exp(oo$par[nb + 1:2])

   oo$value<- -oo$value

   oo
}


W.inv<- function(W, symmetric=TRUE,inverse=TRUE){
   eW <- eigen(W, symmetric=symmetric)
   d <- eW$values
   if (min(d) <0  && abs(min(d))>sqrt(.Machine$double.eps))
       stop("'W' is not positive definite")
   else d[d<=0]<- ifelse(inverse, Inf, 0)
   A <- diag(d^ifelse(inverse, -0.5, 0.5)) %*% t(eW$vector)
   A # t(A)%*%A = W^{-1}
}

# adapted from lm.gls in MASS
lmGls<- function (formula, data, A, subset, na.action, method = "qr", ...) {
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    m$A <- m$method <- m$... <- NULL
    m[[1L]] <- as.name("model.frame")
    m <- eval.parent(m)
    if (method == "model.frame") 
        return(m)
    Terms <- attr(m, "terms")
    Y <- model.response(m)
    X <- model.matrix(Terms, m, contrasts)
    fit <- lm.fit(A %*% X, A %*% Y, method = method, ...)
    class(fit)<- "lm"
    fit
}

# generalized least squares test
scanOne4p.1 <-
   function(y,
            x,
            prdat,
            cov,
            intcovar = NULL,
            test = c("F","Chisq","None"))
{
# prdat$pr: n by 3 by ? matrix, conditional probabilities
# vc: object from estVC or aicVC
# test: “Chisq”, “F” or “Cp”
   gcv<- W.inv(cov)
   test<- match.arg(test)

   nsnp<- dim(prdat$pr)[3]
   if(!is.null(intcovar)) nint<- ncol(as.matrix(intcovar))
   model.par<- vector("list",nsnp)
      names(model.par)<- prdat$snp
   P<- rep(Inf,nsnp)
      names(P)<- prdat$snp
   if(is.null(intcovar)){
      if(!missing(x)){
         oTmp<- data.frame(y=y,x)
      }else{
         oTmp<- data.frame(y=y)
      }
      g0<- lmGls(y~.,data=oTmp,A=gcv)
      if(test=="None"){
         P0<- logLik(g0)
         for(k in 1:nsnp){
            if(!missing(x)){
               oTmp<- data.frame(y=y,x,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
            }else{
               oTmp<- data.frame(y=y,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
            }

            g<- lmGls(y~.,data=oTmp,A=gcv)
            model.par[[k]]<- g$coef
            P[k]<- logLik(g)
         }
         P<- 2*(P-P0)
      }else{
         for(k in 1:nsnp){
            if(!missing(x)){
               oTmp<- data.frame(y=y,x,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
            }else{
               oTmp<- data.frame(y=y,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
            }

            g<- lmGls(y~.,data=oTmp,A=gcv)
            model.par[[k]]<- g$coef
            P[k]<- anova(g0,g,test=test)$P[2]
         }
      }
   }else{
      if(!missing(x)){
         oTmp<- data.frame(y=y,x,intcovar)
      }else{
         oTmp<- data.frame(y=y,intcovar)
      }
      g0<- lmGls(y~.,data=oTmp,A=gcv)
      if(test=="None"){
         P0<- logLik(g0)
         for(k in 1:nsnp){
            if(!missing(x)){
               oTmp<- data.frame(y=y,x,intcovar,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
            }else{
               oTmp<- data.frame(y=y,intcovar,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
            }
            nc<- ncol(oTmp)
            str<- paste(paste("(",paste(colnames(oTmp)[nc-1-(nint:1)],collapse="+"),")",sep=""),
                        paste("(",paste(colnames(oTmp)[(nc-1):nc],collapse="+"),")",sep=""),
                        sep=":")
            str<- paste("y~.+",str,sep="")

            g<- lmGls(formula(str),data=oTmp,A=gcv)
            model.par[[k]]<- g$coef
            P[k]<- logLik(g)
         }
         P<- 2*(P-P0)
      }else{
         for(k in 1:nsnp){
            if(!missing(x)){
               oTmp<- data.frame(y=y,x,intcovar,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
            }else{
               oTmp<- data.frame(y=y,intcovar,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
            }
            nc<- ncol(oTmp)
            str<- paste(paste("(",paste(colnames(oTmp)[nc-1-(nint:1)],collapse="+"),")",sep=""),
                        paste("(",paste(colnames(oTmp)[(nc-1):nc],collapse="+"),")",sep=""),
                        sep=":")
            str<- paste("y~.+",str,sep="")

            g<- lmGls(formula(str),data=oTmp,A=gcv)
            model.par[[k]]<- g$coef
            P[k]<- anova(g0,g,test=test)$P[2]
         }
      }
   }

   list(P=P, parameters=model.par,chr=prdat$chr,dist=prdat$dist)
}

scanOne4p.2 <-
   function(y,
            x,
            gdat,
            cov,
            intcovar = NULL,
            test = c("F","Chisq","None"))
{
# gdat: n by ? matrix, marker data. Markers in columes!!!
# vc: object from estVC or aicVC
# intcover: covariates that interact with QTL
# test: “Chisq”, “F” or “Cp”
   gcv<- W.inv(cov)
   test<- match.arg(test)

   nsnp<- dim(gdat)[2]
   if(!is.null(intcovar)) nint<- ncol(as.matrix(intcovar))
   model.par<- vector("list",nsnp)
      names(model.par)<- colnames(gdat)
   P<- rep(Inf,nsnp)
      names(P)<- colnames(gdat)
   if(is.null(intcovar)){
      if(!missing(x)){
         oTmp<- data.frame(y=y,x)
      }else{
         oTmp<- data.frame(y=y)
      }
      g0<- lmGls(y~.,data=oTmp,A=gcv)
      if(test=="None"){
         P0<- logLik(g0)
         for(j in 1:nsnp){
            if(!missing(x)){
               oTmp<- data.frame(y=y,x,snp=as.factor(gdat[,j]))
            }else{
               oTmp<- data.frame(y=y,snp=as.factor(gdat[,j]))
            }

            g<- lmGls(y~.,data=oTmp,A=gcv)
            model.par[[j]]<- g$coef
            P[j]<- logLik(g)
         }
         P<- 2*(P - P0)
      }else{
         for(j in 1:nsnp){
            if(!missing(x)){
               oTmp<- data.frame(y=y,x,snp=as.factor(gdat[,j]))
            }else{
               oTmp<- data.frame(y=y,snp=as.factor(gdat[,j]))
            }

            g<- lmGls(y~.,data=oTmp,A=gcv)
            model.par[[j]]<- g$coef
            P[j]<- anova(g0,g,test=test)$P[2]
         }
      }
   }else{
      if(!missing(x)){
         oTmp<- data.frame(y=y,x,intcovar)
      }else{
         oTmp<- data.frame(y=y,intcovar)
      }
      g0<- lmGls(y~.,data=oTmp,A=gcv)
      if(test=="None"){
         P0<- logLik(g0)
         for(j in 1:nsnp){
            if(!missing(x)){
               oTmp<- data.frame(y=y,x,intcovar,snp=as.factor(gdat[,j]))
            }else{
               oTmp<- data.frame(y=y,intcovar,snp=as.factor(gdat[,j]))
            }

            nc<- ncol(oTmp)
            str<- paste(colnames(oTmp)[nc-(nint:1)],colnames(oTmp)[nc],collapse="+",sep=":")
            str<- paste("y~.+",str,sep="")

            g<- lmGls(formula(str),data=oTmp,A=gcv)
            model.par[[j]]<- g$coef
            P[j]<- logLik(g)
         }
         P<- 2*(P-P0)
      }else{
         for(j in 1:nsnp){
            if(!missing(x)){
               oTmp<- data.frame(y=y,x,intcovar,snp=as.factor(gdat[,j]))
            }else{
               oTmp<- data.frame(y=y,intcovar,snp=as.factor(gdat[,j]))
            }

            nc<- ncol(oTmp)
            str<- paste(colnames(oTmp)[nc-(nint:1)],colnames(oTmp)[nc],collapse="+",sep=":")
            str<- paste("y~.+",str,sep="")

            g<- lmGls(formula(str),data=oTmp,A=gcv)
            model.par[[j]]<- g$coef
            P[j]<- anova(g0,g,test=test)$P[2]
         }
      }
   }

   list(P=P, parameters=model.par)
}

scanOne4p<- 
   function(y,
            x,
            gdat,
            prdat = NULL,
            vc = NULL,
            intcovar = NULL,
            test = c("F","Chisq","None"),
            minorGenoFreq = 0,
            rmv = TRUE)

{
   UseMethod("scanOne4p")
}

scanOne4p.default<- 
   function(y,
            x,
            gdat,
            prdat = NULL,
            vc = NULL,
            intcovar = NULL,
            test = c("F","Chisq","None"),
            minorGenoFreq = 0,
            rmv = TRUE)
{
   if(!is.null(vc)){
      if(is.element("bgv",attr(vc,"class"))){
         nb<- length(vc$par) - sum(vc$nnl)
         nr<- nrow(vc$y)
         cov<- matrix(0,nrow=nr,ncol=nr)
         for(i in 1:vc$nv)
            if(vc$nnl[i]) cov<- cov + vc$v[[i]]*vc$par[nb+vc$nn[i]]
      }else{
         if(is.data.frame(vc)) vc<- as.matrix(vc)
         if(!is.matrix(vc)) stop("vc should be a matrix.")
         if(!is.numeric(vc)) stop("vc should be a numeric matrix.")
         cov<- vc
      }
   }else cov<- diag(nrow(as.matrix(y)))
   if(!is.null(prdat)){
      pv<- scanOne4p.1(y=y,x=x,prdat=prdat,cov=cov,intcovar=intcovar,test=test)
   }else{
      tb<- sort(union(as.matrix(gdat),NULL))
      tbf<- NULL
      for(ii in tb) tbf<- rbind(tbf,colSums(gdat==ii))
         if(sum(tbf)!=nrow(gdat)*ncol(gdat)) stop("Error occurred.\n")
      tbf<- apply(tbf,2,min)
      idx<- (tbf < nrow(gdat)*minorGenoFreq)
      if(sum(idx)>0){
         if(rmv){
            gdat<- as.matrix(gdat)
            tb<- sort(union(as.matrix(gdat),NULL))
            tbf<- NULL
            for(ii in tb) tbf<- rbind(tbf,colSums(gdat==ii))
               if(sum(tbf)!=nrow(gdat)*ncol(gdat)) stop("Error occurred.\n")
            tbf<- apply(tbf,2,min)
            idx<- (tbf < nrow(gdat)*minorGenoFreq)
            gdat<- gdat[,!idx]
            rm(tb,tbf,ii,idx)
         }else{
            cat("minor genotype frequency is too small at one or more SNPs.\n")
            return(NULL)
        }
      }

      gdat<- as.data.frame(gdat)
      pv<- scanOne4p.2(y=y,x=x,gdat=gdat,cov=cov,intcovar=intcovar,test=test)
   }

   class(pv)<- c("scanOne4p",test)
   pv
}

print.scanOne4p <-
   function(x,...)
{
   tt<- x
   if(length(tt$P)>5){
      cat("P values:\n")
      print(tt$P[1:5])
      cat("... ...\n\n")

      cat("Coefficients:\n")
      print(tt$par[1:5])
      cat("... ...\n\n")
   }else print(tt)
}

