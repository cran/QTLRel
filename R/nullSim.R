
nullSim<- function(y, x, gdat, prdat, ped, gmap, hap,
   method = c("permutation","gene dropping"), vc = NULL, intc = NULL,
   numGeno = FALSE, test = c("None","F","LRT"), minorGenoFreq = 0.05,
   rmv = TRUE, gr = 2, ntimes = 1000){
   matr<- NULL
   method<- match.arg(method)
   test<- match.arg(test)
   if(method=="gene dropping"){
      if(missing(prdat)){
         ids<- rownames(gdat)
         if(any(!is.element(ids, ped$id)))
            stop("Not all sample IDs in both 'prdat' and 'ped'?", call.=FALSE)
         pos<- gmap
         for(n in 1:ntimes){
            gdatTmp<- genoSim(ped, gmap=gmap, ids=ids, hap=hap, method="Haldane")
            llkTmp<- scanOne(y=y, x=x, gdat=gdatTmp, vc=vc, intc=intc,
               numGeno=numGeno, test=test, minorGenoFreq=minorGenoFreq, rmv=rmv)

            if(test == "None"){
               matr<- c(matr, max(llkTmp$LRT))
            }else{
               matr<- c(matr, max(-log10(llkTmp$pval)))
            }
         }
      }else{
         ids<- dimnames(prdat$pr)[[1]]
         if(any(!is.element(ids, ped$id)))
            stop("Not all sample IDs in both 'prdat' and 'ped'?", call.=FALSE)
         pos<- data.frame(snp=prdat$snp, chr=prdat$chr, dist=prdat$dist)
         for(n in 1:ntimes){
            gdatTmp<- genoSim(ped, gmap=gmap, ids=ids, hap=hap, method="Haldane")
            prd<- genoProb(gdatTmp, gmap, gr=gr, pos=pos, method="Haldane", msg = FALSE)
            llkTmp<- scanOne(y=y, x=x, prdat=prd, vc=vc, intc=intc,
               numGeno=numGeno, test=test, minorGenoFreq=minorGenoFreq, rmv=rmv)

            if(test == "None"){
               matr<- c(matr, max(llkTmp$LRT))
            }else{
               matr<- c(matr, max(-log10(llkTmp$pval)))
            }
         }
      }
   }else if(method=="permutation"){
      if(missing(prdat)){
         if(missing(gdat)) stop("Either 'gdat' or 'prdat' should be provided.", call.=FALSE)
         for(n in 1:ntimes){
            idx<- sample(1:nrow(gdat),replace=FALSE)
            llkTmp<- scanOne(y=y, x=x, gdat=gdat[idx,], vc=vc, intc=intc,
               numGeno=numGeno, test=test, minorGenoFreq=minorGenoFreq, rmv=rmv)

            if(test == "None"){
               matr<- c(matr, max(llkTmp$LRT))
            }else{
               matr<- c(matr, max(-log10(llkTmp$pval)))
            }
         }
      }else{
         prd<- prdat
         for(n in 1:ntimes){
            idx<- sample(1:dim(prdat$pr)[1], replace=FALSE)
            prd$pr<- prdat$pr[idx,,]
            llkTmp<- scanOne(y=y, x=x, prdat=prd, vc=vc, intc=intc,
               numGeno=numGeno, test=test, minorGenoFreq=minorGenoFreq, rmv=rmv)

            if(test == "None"){
               matr<- c(matr, max(llkTmp$LRT))
            }else{
               matr<- c(matr, max(-log10(llkTmp$pval)))
            }
         }
      }
   }else stop("Permutation or gene dropping.", call.=FALSE)

   matr
}


