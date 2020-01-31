
###########################################
# generate genotype data by gene dropping #
###########################################

#######################################
# generate haplotype or genotype data #
#######################################
# the follow code generate genotype data chrosomome by chromosome

.hapSim<- function(ped,gmap,ids,hap,method=c("Haldane","Kosambi"),genotype=FALSE){
   ped<- pedRecode(ped,ids=ids,all=TRUE,msg=FALSE)
   idx<- ped$father <= 0 | ped$mother <= 0 # founders
   if(any(idx)){
      pedTmp<- ped[idx,]
      idTmp<- trim(sapply(pedTmp$old.id,as.character))
   }else{
      stop("No founders? Unbelievable.", call.=FALSE)
   }
   rm(idx)
   if(!missing(hap)){
      if(any(hap<0))
         stop("hap: Use non-negative integers for alleles.", call.=FALSE)
      if(any(hap>2^14))
         stop("hap: Too large integers.", call.=FALSE)
      idx<- match(idTmp,trim(rownames(hap)))
      if(any(is.na(idx))){
         print(idTmp[is.na(idx)])
         stop("Above founders' haplotype data are not supplied. Make sure row names of 'hap' are founders' IDs... You may run 'pedRecode' to find out which are founders.", call.=FALSE)
      }
      hap<- hap[idx,]
      if(sum(idx) < 2)
         hap<- matrix(hap, nrow=1)
      rm(idx)
   }else{
      nr<- nrow(pedTmp)
      nc<- nrow(gmap)
      hap<- matrix(NA, nrow=nr, ncol=nc*2)
         rownames(hap)<- trim(pedTmp$old.id)
      for(n in 1:nr){
         hap[n,]<- rep(n,nc*2)
         if(is.element(substr(tolower(trim(pedTmp$sex[n])),1,1), c("1","m"))){
            idx<- (1:nc)[tolower(gmap$chr)=="x"]*2-1
            hap[n,idx]<- 0
         }
      }
      rm(nr,nc,n,idx)
   }
   rm(pedTmp,idTmp)
   if(!is.data.frame(gmap) || any(!is.element(c("snp","chr","dist"),colnames(gmap)))){
      stop("Genetic map should be a data frame (snp,chr,dist,...).", call.=FALSE)
   }
   if(!missing(ids)){
      ids<- trim(ids)
      if(length(ids)==0)
         stop("ids not correctly specified.", call.=FALSE)
      ii<- match(ids,ped$old.id)
      if(any(is.na(ii)))
         stop("Check ids for error.", call.=FALSE)
   }else ii<- 1:length(ped$old.id)

   gmap$chr<- reorder(factor(gmap$chr))
      ord<- order(gmap$chr,gmap$dist)
      gmap<- gmap[ord,]
   if(is.null(gmap$recRate)){
      if(!is.numeric(gmap$dist)){
         stop("gmap$dist should be numeric.", call.=FALSE)
      }
      method<- match.arg(method)
      dstTmp<- diff(gmap$dist)
         dstTmp[dstTmp<0]<- 0
         dstTmp<- c(0,dstTmp)/100
      gmap$recRate<- mappingFuncInv(dstTmp,method=method)
   }
   chr<- unique(gmap$chr)
      chr<- reorder(chr)
      chr<- sort(chr)
   nchr<- length(chr)
   out<- NULL
   for(n in 1:nchr){
      tmp<- .hapSim0(ped,gmap,chr[n],hap,genotype)
      out<- cbind(out,tmp)
      rm(tmp)
      gc(reset=TRUE)
   }
   out<- out[ii,];
      out<- as.matrix(out)
   rownames(out)<- ped$old.id[ii]
   idx<- order(ord)
   if(genotype){
      as.matrix(out[,idx])
   }else{
      idx<- rbind(idx*2-1,idx*2)
      as.matrix(out[,idx])
   }
}

.hapSim0<- function(pedd,gmap,chr,hap,genotype){
   pedd<- pedd[,c("id","father","mother","sex")]
   if(is.numeric(pedd[,"sex"])){
      pedd[,"sex"]<- pedd[,"sex"]==1
   }else{
      pedd[,"sex"]<- pedd[,"sex"]=="M" | pedd[,"sex"]=="Male" |
                     pedd[,"sex"]=="m" | pedd[,"sex"]=="male"
   }
   idx<- gmap$chr==chr
   rr<- gmap$recRate[idx]
   nr<- nrow(pedd)
   nc<- sum(idx)
   xchr<- FALSE
      if(chr=="x" || chr=="X")
         xchr<- TRUE
   gdat<- matrix(-99,nrow=nr,ncol=2*nc)
   if(!missing(hap)){
      ninit<- nrow(hap)
      gdat[1:ninit,]<- hap[,rep(idx,rep(2,length(idx)))]
   }
   if(genotype && nrow(hap) > 2){
      stop("You need to get genotypes manually unless there are more than two founders.", call.=FALSE)
   }
   out<- .C("rgdata2",
            gdata = as.integer(t(gdat)),
               nr = as.integer(nr),
               nc = as.integer(nc),
            ninit = as.integer(ninit),
         pedigree = as.integer(t(pedd)),
           recomb = as.double(rr),
             xchr = as.logical(xchr),
          PACKAGE = "QTLRel")$gdata
   out[out==-99]<- NA
   out<- matrix(out,nrow=nr,byrow=TRUE)
      storage.mode(out)<- "integer"

   if(genotype){
      oo<- out[,2*(1:nc)-1] + out[,2*(1:nc)]; oo<- as.matrix(oo)
      if(xchr){
         ii<- as.logical(pedd[,"sex"])
         if(any(out[ii,2*(1:nc)-1]!=0))
            stop("Paternal wrong! check pedigree for sex errors.", call.=FALSE)
         if(any(out[!ii,2*(1:nc)-1]==0))
            stop("Maternal wrong! check pedigree for sex errors.", call.=FALSE)
         oo[ii,]<- 2*oo[ii,]
      }
      oo<- oo-1
      storage.mode(oo)<- "integer"
      colnames(oo)<- gmap$snp[gmap$chr==chr]
      oo
   }else{
      out
   }
}

hapSim<- function(ped,gmap,ids,hap,method=c("Haldane","Kosambi")){
# ped: recoded pedigree
# gmap: genetic map (snp,chr,recom,dist,...)
# hap: founders' haplotypes
# ids: only output data for individuals with ID ids
   if(missing(gmap)){
      gmap<- data.frame(snp="N",chr="N",dist=0)
      if(!missing(hap)){
         hap<- as.matrix(hap)
         if(ncol(hap) < 2)
            stop("'hap' should have at least 2 columns.", call.=FALSE)
         hap<- hap[,1:2]
      }
   }
   out<- .hapSim(ped = ped,
                gmap = gmap,
                 ids = ids,
                 hap = hap,
              method = method,
            genotype = FALSE)
   out
}

genoSim<- function(ped,gmap,ids,hap,method=c("Haldane","Kosambi")){
# ped: recoded pedigree
# gmap: genetic map (snp,chr,recom,dist,...)
# output: individuals by SNPs
#   for each SNP 1--AA, 2--AB,3--BB
# ids: only output data for individuals with ID ids
   if(missing(gmap)){
      gmap<- data.frame(snp="N",chr="N",dist=0)
      if(!missing(hap)){
         hap<- as.matrix(hap)
         if(ncol(hap) < 2)
            stop("'hap' should have at least 2 columns.", call.=FALSE)
         hap<- hap[,1:2]
      }
   }
   out<- .hapSim(ped = ped,
                gmap = gmap,
                 ids = ids,
                 hap = hap,
              method = method,
            genotype = TRUE)

   out
}

