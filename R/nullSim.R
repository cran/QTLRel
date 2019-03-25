
nullSim<- function(y, x, gdat, prdat, ped, gmap, hap,
   method=c("permutation","gene dropping"), vc=NULL, intc=NULL,
   test = c("None","F","Chisq"), minorGenoFreq=0.05, rmv=TRUE,
   gr=2, ntimes=10){
   matr<- NULL
   method<- match.arg(method)
   test<- match.arg(test)
   if(method=="gene dropping"){
      pedR<- pedRecode(ped,all=TRUE,msg=FALSE)
      if(missing(prdat)){
         ids<- rownames(gdat)
         if(any(!is.element(ids, pedR$old)))
            stop("Not all sample IDs in both 'prdat' and 'ped'?", call.=FALSE)
         pos<- gmap
         for(n in 1:ntimes){
            gdatTmp<- genoSim(pedR, gmap=gmap, ids=ids, hap=hap, method="Haldane")
            llkTmp<- scanOne(y=y, x=x, gdat=gdatTmp, vc=vc, intc=intc,
               test=test, minorGenoFreq=minorGenoFreq, rmv=rmv)

            if(minorGenoFreq <= 0 && rmv) matr<- rbind(matr, llkTmp$p) else
            if(test != "None") stop("Test should be 'None'.", call.=FALSE) else
            matr<- rbind(matr, max(llkTmp$p))
         }
      }else{
         ids<- dimnames(prdat$pr)[[1]]
         if(any(!is.element(ids, pedR$old)))
            stop("Not all sample IDs in both 'prdat' and 'ped'?", call.=FALSE)
         pos<- data.frame(snp=prdat$snp, chr=prdat$chr, dist=prdat$dist)
         for(n in 1:ntimes){
            gdatTmp<- genoSim(pedR, gmap=gmap, ids=ids, hap=hap, method="Haldane")
            prd<- genoProb(gdatTmp, gmap, gr=gr, pos=pos, method="Haldane", msg = FALSE)
            llkTmp<- scanOne(y=y, x=x, prdat=prd, vc=vc, intc=intc,
               test=test, minorGenoFreq=minorGenoFreq, rmv=rmv)

            if(minorGenoFreq <= 0 && rmv) matr<- rbind(matr, llkTmp$p) else
            if(test != "None") stop("Test should be 'None'.", call.=FALSE)
            else matr<- rbind(matr, max(llkTmp$p))
         }
      }
   }else if(method=="permutation"){
      if(missing(prdat)){
         if(missing(gdat)) stop("Either 'gdat' or 'prdat' should be provided.", call.=FALSE)
         for(n in 1:ntimes){
            idx<- sample(1:nrow(gdat),replace=FALSE)
            llkTmp<- scanOne(y=y, x=x, gdat=gdat[idx,], vc=vc, intc=intc,
               test=test, minorGenoFreq=minorGenoFreq, rmv=rmv)

            matr<- rbind(matr, llkTmp$p)
         }
      }else{
         prd<- prdat
         for(n in 1:ntimes){
            idx<- sample(1:dim(prdat$pr)[1], replace=FALSE)
            prd$pr<- prdat$pr[idx,,]
            llkTmp<- scanOne(y=y, x=x, prdat=prd, vc=vc, intc=intc,
               test=test, minorGenoFreq=minorGenoFreq, rmv=rmv)

            matr<- rbind(matr, llkTmp$p)
         }
      }
   }else stop("Permutation or gene dropping.", call.=FALSE)

   matr
}


