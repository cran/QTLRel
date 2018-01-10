
#######################
# recode the pedigree #
#######################

pedRecode <- function(ped,ids,all=TRUE,msg=TRUE){
# ped: pedigree (id, sire, dam,...) or  (id, generation, sire, dam,...)
#    missing values in sire or dam represented by 0 or NA
# ids: IDs of interest if not missing
# output: id=1,2,...
   strf<- function(ped, str){
      if(is.numeric(ped[,str])){
         idx<- sapply(ped[,str],Im) != 0
      }else{
         sr<- sapply(ped[,str],as.character)
            sr<- trim(sr)
         nch<- sapply(sr,nchar)
         idx<- tolower(substr(sr,nch,nch)) == "i"
         rm(sr,nch)
      }
      usr<- sapply(unique(ped[idx,str]), as.character)
         usr<- usr[!is.na(usr)]
      ii<- match(sapply(ped[,str],as.character),usr)
      idx<- !is.na(ii)

      list(idx=idx, values=-ii[idx])
   }

   ped<- as.data.frame(ped)
   if(is.null(ped$id)){
      stop("'id' missing...", call.=FALSE)
   }else{
      if(is.numeric(ped$id) | is.complex(ped$id)){
         idTmp<- sapply(Re(ped$id),as.character)
         idx<- Im(ped$id) != 0
         if(any(idx)){
            idTmp[idx]<- sapply(ped$id[idx],as.character)
         }
         ped$id<- idTmp
      }
      ped$id<- trim(ped$id)
   }
   if(is.null(ped$sire)){
      stop("'sire' missing...", call.=FALSE)
   }else{
      if(is.numeric(ped$sire) | is.complex(ped$sire)){
         idTmp<- sapply(Re(ped$sire),as.character)
         idx<- Im(ped$sire) != 0
         if(any(idx)){
            idTmp[idx]<- sapply(ped$sire[idx],as.character)
         }
         ped$sire<- idTmp
      }
      ped$sire<- trim(ped$sire)
   }
   if(is.null(ped$dam)){
      stop("'dam' missing...", call.=FALSE)
   }else{
      if(is.numeric(ped$dam) | is.complex(ped$dam)){
         idTmp<- sapply(Re(ped$dam),as.character)
         idx<- Im(ped$dam) != 0
         if(any(idx)){
            idTmp[idx]<- sapply(ped$dam[idx],as.character)
         }
         ped$dam<- idTmp
      }
      ped$dam<- trim(ped$dam)
   }

   if(all){
      idEx<- union(sapply(ped$sire,as.character),sapply(ped$dam,as.character))
         idEx<- setdiff(idEx,sapply(ped$id,as.character))
         idEx<- idEx[!is.na(idEx) & idEx != 0 & idEx != "0"]
         nch<- sapply(idEx,nchar)
         idx<- tolower( substr(idEx,nch,nch) ) != "i"
         idEx<- idEx[idx]
      if(length(idEx) > 0){
         cns<- colnames(ped)
         mtr<- matrix(NA,nrow=length(idEx),ncol=length(cns))
         if(!is.null(ped$generation)){
            idx<- match("generation",cns)
            mtr[,idx]<- rep("-999", length(idEx))
         }
         if(!is.null(ped$sex)){
            sx<- rep("M", length(idEx))
               idx<- is.element(idEx,sapply(ped$dam,as.character))
               sx[idx]<- "F"
            idx<- match("sex",cns)
            mtr[,idx]<- sx
            rm(sx)
         }
         idx<- match("id",cns)
         mtr[,idx]<- idEx
         idx<- match(c("sire","dam"),cns)
         mtr[,idx]<- 0
         ped<- rbind(mtr,as.matrix(ped))
         rownames(ped)<- NULL
         rm(cns,mtr,idx)
      }
      rm(idEx)
      ped<- as.data.frame(ped)
   }
   pedSave<- ped

   idx<- is.na(ped$id) | ped$id==0 | ped$id=="0"
   if(any(idx) && msg){
      print(ped[idx,])
      cat("   Above individuals with N/A IDs were removed.\a\n")
      ped<- ped[!idx,]
   }
   rm(idx)

   if(is.numeric(ped$id)){
      im<- sapply(ped$id,Im)
      idx<- im != 0
      rm(im)
   }else{
      idTmp<- sapply(ped$id,as.character)
         idTmp<- trim(idTmp)
         nch<- sapply(idTmp,nchar)
         idTmp<- tolower( substr(idTmp,nch,nch) )
      idx<- idTmp == "i"
      rm(idTmp,nch)
   }
   if(any(idx) && msg){
      print(ped[idx,])
      cat("   Above founder lines were removed.\a\n")
      ped<- ped[!idx,]
   }
   rm(idx)

   if(!missing(ids)){# discard irrelevant ones
      ids<- trim(ids)
      if(length(ids)==0) stop("IDs not correctly specified.", call.=FALSE)
      idx<- !is.element(ids,ped$id)
      if(any(idx) && msg){
         print(ids[idx])
         cat("   Above IDs are out of range and ignored!\a\n")
      }
      ids<- ids[!idx]

      idTmp<- ids
      idx<- rep(FALSE,nrow(ped))
      while(TRUE){
         idxTmp<- is.element(ped$id,idTmp)
         idx<- idx | idxTmp
         idx<- idx | is.element(ped$id,ped$sire[idxTmp]) |
                     is.element(ped$id,ped$dam[idxTmp])
         idTmp<- ped$id[idx]

         if(sum(idxTmp)==sum(idx)) break
      }
      ped<- ped[idx,]
      rm(idTmp,idx,idxTmp)
   }

   if(is.null(ped$generation))
      ped<- pedRecode.0(ped,msg)
   ped$generation<- trim(ped$generation)
   ids<- paste(ped$generation,ped$id,sep="~")
   uids<- unique(ids)
      idx<- match(uids,ids)
   if(length(uids)<length(ids) && msg){
      cat("   The following samples are repeated:\a\n")
      print(ped[-idx,c("id","generation","sire","dam")])
      cat("   Repeated IDs were excluded!\a\n")
   }
   ped<- ped[idx,]
   rm(idx,ids)

   # new code
   ped$generation<- reorder(factor(ped$generation))
   ped$id<- reorder(factor(ped$id))
   ped$sire<- reorder(factor(ped$sire))
   ped$dam<- reorder(factor(ped$dam))
   ped<- ped[order(ped$generation,ped$id,ped$sire,ped$dam),]
   idd<- data.frame(index=c(0,0),id=c(NA,0))
      idd<- rbind(idd,data.frame(index=1:nrow(ped),id=ped$id))
   ### recode here

   # recode IDs
   if(is.null(ped$sex)){
      idx<- match(sapply(ped$sire,as.character),idd$id)
         idx[is.na(idx)]<- 1
      sire<- idd$index[idx]
         str<- strf(ped,"sire")
         if(any(str$idx))
            sire[str$idx]<- str$val-10000
      idx<- match(sapply(ped$dam,as.character),idd$id)
         idx[is.na(idx)]<- 1
      dam<- idd$index[idx]
         str<- strf(ped,"dam")
         if(any(str$idx))
            dam[str$idx]<- str$val-20000
      ped<- data.frame(id=idd$index[-c(1:2)],
                       sire=sire,
                       dam=dam,
                       generation=ped$generation,
                       old.id=ped$id)
   }else{
      idx<- match(sapply(ped$sire,as.character),idd$id)
         idx[is.na(idx)]<- 1
      sire<- idd$index[idx]
         str<- strf(ped,"sire")
         if(any(str$idx))
            sire[str$idx]<- str$val-10000
      idx<- match(sapply(ped$dam,as.character),idd$id)
         idx[is.na(idx)]<- 1
      dam<- idd$index[idx]
         str<- strf(ped,"dam")
         if(any(str$idx))
            dam[str$idx]<- str$val-20000
      ped<- data.frame(id=idd$index[-c(1:2)],
                       sire=sire,
                       dam=dam,
                       sex=trim(ped$sex),
                       generation=ped$generation,
                       old.id=ped$id)
      ii<- match(sapply(ped$sire,as.character),ped$id)
         ii<- ii[!is.na(ii)]
      if(length(ii)>0){
         idx<- !is.element(ped$sex[ii], c("0", "1", "M", "Male"))
            idx<- unique(ii[idx])
         if(any(idx) && msg){
            cat("   Suppose 1, M or Male stands for male...\n")
            cat("   --------------------------------------\n")
            print(ped[idx,])
            cat("   --------------------------------------\n")
            cat("   Above should be sire(s)...\a\n\n")
         }
      }
      jj<- match(sapply(ped$dam,as.character),ped$id)
         jj<- jj[!is.na(jj)]
      if(length(jj)>0){
         idx<- is.element(ped$sex[jj], c("1", "M", "Male"))
            # idx<- idx | is.na(ped$sex[jj])
            idx<- unique(jj[idx])
         if(any(idx) && msg){
            cat("   Suppose !0, !1, !M and !Male stands for female...\n")
            cat("   --------------------------------------\n")
            print(ped[idx,])
            cat("   --------------------------------------\n")
            cat("   Above should be dam(s)...\a\n\n")
         }
      }
   }
   idx<- (ped$sire > ped$id) | (ped$dam > ped$id)
   if(any(idx)){
      idx<- match(sapply(ped$old[idx],as.character),pedSave$id)
      print(pedSave[idx,][1:min(sum(idx),3),])
      cat("... ...\n")
      stop("Check the above for errors...", call.=FALSE)
   }
   idx.s<- ped$generation[match(sapply(ped$sire,as.character), ped$id)] == ped$generation
      idx.s<- (1:length(idx.s))[idx.s]
      idx.s<- idx.s[!is.na(idx.s)]
   idx.d<- ped$generation[match(sapply(ped$dam,as.character), ped$id)] == ped$generation
      idx.d<- (1:length(idx.d))[idx.d]
      idx.d<- idx.d[!is.na(idx.d)]
   if(length(idx.s)>0 || length(idx.d)>0){
      idx<- c(idx.s, idx.d)
         idx<- sort(unique(idx))
      idx<- c(sort(unique(c(ped$sire[idx.s], ped$dam[idx.d]))),idx)
      pedTmp<- ped[idx,]
         pedTmp$id<- ped$old[match(sapply(pedTmp$id,as.character),ped$id)]
         pedTmp$sire<- ped$old[match(sapply(pedTmp$sire,as.character),ped$id)]
         pedTmp$dam<- ped$old[match(sapply(pedTmp$dam,as.character),ped$id)]
      pedTmp$old.id<- NULL
      print(pedTmp)
      stop("Check the above for errors regarding generations...", call.=FALSE)
   }
   rownames(ped)<- 1:nrow(ped)

   ped
}

# create "generation"
pedRecode.0<- function(ped,msg){
# ped: data frame (id,sire,dam,...)
   ped<- as.data.frame(ped)
      ped$generation<- NULL
   if(is.null(ped$id)){
      stop("'id' missing...", call.=FALSE)
   }else ped$id<- trim(ped$id)
   if(is.null(ped$sire)){
      stop("'sire' missing...", call.=FALSE)
   }else ped$sire<- trim(ped$sire)
   if(is.null(ped$dam)){
      stop("'dam' missing...", call.=FALSE)
   }else ped$dam<- trim(ped$dam)

   idx<- is.na(ped$id) | ped$id==0 | ped$id=="0"
   if(any(idx) && msg){
      print(ped[idx,])
      cat("   Above individuals with N/A IDs were removed.\n")
      ped<- ped[!idx,]
   }

   nr<- nrow(ped)
   ii<- matrix(TRUE,nrow=1,ncol=nr)
   ii0<- ii
   while(1){
      ii0<- is.element(ped$id,ped$sire[ii0]) |
            is.element(ped$id,ped$dam[ii0])
      if(any(ii0)){
         ii<- rbind(ii0,ii)
      }else break
   }
# no offspring or parents
#   idx<- is.element(ped$id,ped$sire) |
#         is.element(ped$id,ped$dam)  |
#         is.element(ped$sire,ped$id) |
#         is.element(ped$dam,ped$id)
#   ii[1,]<- ii[1,] | !idx

   idx<- ii[1,]
   idx0<- idx
   jj<- 0
   out<- cbind(generation=jj,ped[idx,])
   if(nrow(ii)>1){
      for(n in 2:nrow(ii)){
         idx0<- idx0 | ii[n-1,]
         idx<- ii[n,] & !idx0
         if(any(idx)){
            jj<- jj+1
            out<- rbind(out,cbind(generation=jj,ped[idx,]))
         }
      }
   }
   jj<- match("id",colnames(out))
   out<- cbind(id=out[,jj],generation=out$generation,out[-c(1,jj)])
   rownames(out)<- 1:nrow(out)

   out$generation<- reorder(factor(out$generation))
   out$id<- reorder(factor(out$id))
   out$sire<- reorder(factor(out$sire))
   out$dam<- reorder(factor(out$dam))

   out<- out[order(out$generation,out$id,out$sire,out$dam),]
   out
}

