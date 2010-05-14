
smpl<- function(pr,x){
   n<- nrow(pr)
   y<- rep(NA,n)
   for(i in 1:n) y[i]<- sample(x,size=1,prob=pr[i,])
   y
}

# estimate variances of (genetic) components
gvar<- function (par, ov){
# par: parameter estimates
# ov: ov$v is a list of genetic matrices
   str<- c("AA", "DD", "HH", "AD", "MH", "EE")
   if(!is.null(par)) par<- par[str]
   ov<- ov$v[str]
   vv <- rep(NA, 6)
      names(vv)<- str
   for (i in 1:6){
      v0<- ov[[i]]*par[i]
         v0<- diag(v0)
      vv[i]<- mean(v0)
   }
   vv
}

# estimate QTL variances
qtlVar<- function(lrt,prdat,simulation=FALSE,nsim=25){
# lrt: data frame (a,d,g,e,...)
# probs: prDat$pr
# estimated genetic variance-covariance matrix
   vv<- rep(NA,nrow(lrt))
   for(ii in 1:nrow(lrt)){
      tmp<- lrt[ii,]
      prd<- prdat[,,ii]

      if(!simulation){
         tmp1<- sweep(prd,2,c(tmp$a,tmp$d,-tmp$a),"*")
            tmp1<- rowSums(tmp1) # mean
         tmp2<- cbind(tmp$a-tmp1,tmp$d-tmp1,-tmp$a-tmp1)
            tmp2<- tmp2^2
            tmp2<- tmp2*prd
            tmp2<- rowSums(tmp2) # variance

         vv[ii]<- var(tmp1) + mean(tmp2)
      }else{# simulation method -- takes time
         vr<- rep(NA,100)
         for(i in 1:100){
            vr.<- rep(NA,nrow(prd))
            pr.<- runif(nrow(prd), min=0, max=1)
            idx1<- pr. <= prd[,1]
            idx2<- (!idx1) & (pr. <= prd[,1]+prd[,2])
            idx3<- !idx1 & !idx2
            if(any(idx1)) vr.[idx1]<- tmp$a
            if(any(idx2)) vr.[idx2]<- tmp$d
            if(any(idx3)) vr.[idx3]<- -tmp$a
            vr[i]<- var(vr.)
         }
         vv[ii]<- mean(vr)
      }
   }
   vv
}

###################
# no need for now #
###################

qtlvarPr<- function(y,lrt,prdat,lodci,gcov=NULL){
# one QTL only
   if(missing(lodci)){
      idx<- lrt$lr==max(lrt$lr)
      lodci<- data.frame(chr=prdat$chr[idx][1],
                         lower=0,
                         upper=Inf)
      rm(idx)
   }
   out<- NULL
   if(nrow(lodci)>0) for(n in 1:nrow(lodci)){
      ii<- prdat$chr==lodci$chr[n]
         ii<- ii & lrt$lr==max(lrt$lr[ii])
         ii<- ii & prdat$dist >= lodci$lower[n]
         ii<- ii & prdat$dist <= lodci$upper[n]
         jj<- c(1:length(ii))[ii]
         ii<- jj[1]
      tmp<- lrt[ii,]
      prd<- prdat$pr[,,ii]

      if(!is.null(gcov)){
         gcv<- gcov*tmp$g
         diag(gcv)<- diag(gcv) + tmp$e
         iv<- chol(gcv,pivot=TRUE)
            iv<- t(iv[,order(attr(iv, "pivot"))])
            iv<- diag(sqrt(diag(gcv)))%*%solve(iv)
         yy<- iv%*%as.matrix(y)
      }
      tmp1<- sweep(prd,2,c(tmp$a,tmp$d,-tmp$a),"*")
         tmp1.0<- rowSums(tmp1)
         if(!is.null(gcov)) tmp1<- iv%*%tmp1
         tmp1<- rowSums(tmp1)
      v<- var(tmp1)
#      tmp2<- cbind(tmp$a-tmp1.0,tmp$d-tmp1.0,-tmp$a-tmp1.0)
#         if(!is.null(gcov)) tmp2<- iv%*%tmp2
#         tmp2<- tmp2^2
#         tmp2<- tmp2*prd
#         tmp2<- rowSums(tmp2)
#         v<- v + mean(tmp2) # no need due to llk procedure
      out<- c(out,v/var(yy))
   }
   out*100
}


