
nullSim<- function(y,x,gdat,ped,gmap,ids,method=c("permutation","gene dropping"),
   vc=NULL,intcovar=NULL,minorGenoFreq=0.05,rmv=TRUE,nit=25,ntimes=10){
   matr<- rep(NA,ntimes)
   method<- match.arg(method)
   if(method=="gene dropping"){
      pedR<- pedRecode(ped)
      for(n in 1:ntimes){
         gdatTmp<- genoSim(pedR,gmap,ids=ids,method="Haldane",
            recode.pedigree=FALSE)

         llkTmp<- scanOne(y=y,x=x,gdat=gdatTmp,vc=vc,intcovar=intcovar,
            minorGenoFreq=minorGenoFreq,rmv=rmv,nit=nit)
         if(is.null(llkTmp)) next

         matr[n]<- max(llkTmp$lr)
      }
   }else if(method=="permutation"){
      if(missing(gdat)) stop("gdat must be supplied.")
      idx<- sample(1:nrow(gdat),replace=FALSE)
      llkTmp<- scanOne(y=y,x=x,gdat=gdat[idx,],vc=vc,intcovar=intcovar,
         minorGenoFreq=minorGenoFreq,rmv=rmv,nit=nit)
      if(is.null(llkTmp)) next

      matr[n]<- max(llkTmp$lr)
   }else stop("permutation or gene dropping.")

   matr
}


