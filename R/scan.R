
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
scanOne.1 <-
   function(y,
            x,
            prdat,
            cov,
            intcovar = NULL,
            test = c("None","F","Chisq"))
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

   list(chr=prdat$chr,
        dist=prdat$dist,
        p=P,
        parameters=model.par)
}

scanOne.2 <-
   function(y,
            x,
            gdat,
            cov,
            intcovar = NULL,
            test = c("None","F","Chisq"))
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

   list(p=P,
        parameters=model.par)
}

scanOne<- 
   function(y,
            x,
            gdat,
            prdat = NULL,
            vc = NULL,
            intcovar = NULL,
            test = c("None","F","Chisq"),
            minorGenoFreq = 0,
            rmv = TRUE)

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
            test = c("None","F","Chisq"),
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
      pv<- scanOne.1(y=y,x=x,prdat=prdat,cov=cov,intcovar=intcovar,test=test)
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
      pv<- scanOne.2(y=y,x=x,gdat=gdat,cov=cov,intcovar=intcovar,test=test)
   }

   class(pv)<- c("scanOne")
   pv
}

print.scanOne <-
   function(x,...)
{
   tt<- x; class(tt)<- NULL
      tt$parameters<- NULL
      tt<- as.data.frame(tt)
   if(length(tt$p)>5){
      cat("Test statistic:\n")
      print(tt[1:5,])
      cat("... ...\n\n")

      cat("Coefficients:\n")
      print(x$par[1:5])
      cat("... ...\n\n")
   }else{
      cat("Test statistic:\n")
      print(tt)
      cat("\n")

      cat("Coefficients:\n")
      print(x$par)
   }
}

