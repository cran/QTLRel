
#######################
# recode the pedigree #
#######################

chechPedColNames<- function(ped){
   if(is.element("sire",colnames(ped))){
      if(is.element("father",colnames(ped)))
         stop("Only either 'father' or 'sire' represents paternal parent", call.=FALSE)
      idx<- match("sire",tolower(colnames(ped)))
      colnames(ped)[idx]<- "father"
   }else if(!is.element("father",colnames(ped))){
      stop("'father/sire' missing", call.=FALSE)
   }
   if(is.element("dam",colnames(ped))){
      if(is.element("mother",colnames(ped)))
         stop("Only either 'mother' or 'dam' represents maternal parent", call.=FALSE)
      idx<- match("dam",tolower(colnames(ped)))
      colnames(ped)[idx]<- "mother"
   }else if(!is.element("mother",colnames(ped))){
      stop("'mother/dam' missing", call.=FALSE)
   }

   if(!is.element("id",colnames(ped))){
      stop("'id' missing...", call.=FALSE)
   }

   as.data.frame(ped)
}

pedRecode <- function(ped,ids,all=TRUE,msg=TRUE){
# ped: pedigree (id, father, mother,...) or  (id, generation, father, mother,...)
#    missing values in father or mother represented by 0 or NA
# ids: IDs of interest if not missing
# output: id=1,2,...
   strf<- function(ped, str){
      if(is.numeric(ped[,str])){
         idx<- sapply(ped[,str],Im) != 0
      }else{
         sr<- sapply(ped[,str],as.character)
            sr<- trimws(sr)
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

   ped<- chechPedColNames(ped)
   if(is.numeric(ped$id) | is.complex(ped$id)){
      idTmp<- sapply(Re(ped$id),as.character)
      idx<- Im(ped$id) != 0
      if(any(idx)){
         idTmp[idx]<- sapply(ped$id[idx],as.character)
      }
      ped$id<- idTmp
      rm(idTmp, idx)
   }
   ped$id<- trimws(ped$id)
   if(is.null(ped$father)){
      if(is.null(ped$father)){
         stop("'father/sire' missing...", call.=FALSE)
      }else{
         nTmp<- match("father",colnames(ped))
         colnames(ped)[nTmp]<- "father"
         rm(nTmp)
      }
   }else{
      if(is.numeric(ped$father) | is.complex(ped$father)){
         idTmp<- sapply(Re(ped$father),as.character)
         idx<- Im(ped$father) != 0
         if(any(idx)){
            idTmp[idx]<- sapply(ped$father[idx],as.character)
         }
         ped$father<- idTmp
         rm(idTmp, idx)
      }
      ped$father<- trimws(ped$father)
   }
   if(is.null(ped$mother)){
      if(is.null(ped$mother)){
         stop("'mother/dam' missing...", call.=FALSE)
      }else{
         nTmp<- match("mother",colnames(ped))
         colnames(ped)[nTmp]<- "mother"
         rm(nTmp)
      }
   }else{
      if(is.numeric(ped$mother) | is.complex(ped$mother)){
         idTmp<- sapply(Re(ped$mother),as.character)
         idx<- Im(ped$mother) != 0
         if(any(idx)){
            idTmp[idx]<- sapply(ped$mother[idx],as.character)
         }
         ped$mother<- idTmp
         rm(idTmp, idx)
      }
      ped$mother<- trimws(ped$mother)
   }

   if(all){
      idEx<- union(sapply(ped$father,as.character),sapply(ped$mother,as.character))
         idEx<- setdiff(idEx,sapply(ped$id,as.character))
         idEx<- idEx[!is.na(idEx) & idEx != 0 & idEx != "0"]
         nch<- sapply(idEx,nchar)
         idx<- tolower( substr(idEx,nch,nch) ) != "i"
         idEx<- idEx[idx]
         rm(nch, idx)
      if(length(idEx) > 0){
         cns<- colnames(ped)
         mtr<- matrix(NA,nrow=length(idEx),ncol=length(cns))
         if(!is.null(ped$generation)){
            idx<- match("generation",cns)
            mtr[,idx]<- rep("-999", length(idEx))
            rm(idx)
         }
         if(!is.null(ped$sex)){
            sx<- rep("M", length(idEx))
               idx<- is.element(idEx,sapply(ped$mother,as.character))
               sx[idx]<- "F"
            idx<- match("sex",cns)
            mtr[,idx]<- sx
            rm(sx, idx)
         }
         idx<- match("id",cns)
         mtr[,idx]<- idEx
         idx<- match(c("father","mother"),cns)
         mtr[,idx]<- 0
         ped<- rbind(mtr,as.matrix(ped))
         rownames(ped)<- NULL
         rm(cns, mtr, idx)
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
         idTmp<- trimws(idTmp)
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
      ids<- trimws(ids)
      if(length(ids)==0) stop("IDs not correctly specified.", call.=FALSE)
      idx<- !is.element(ids,ped$id)
      if(any(idx) && msg){
         print(ids[idx])
         cat("   Above IDs are out of range and ignored!\a\n")
      }
      ids<- ids[!idx]
      rm(idx)

      idTmp<- ids
      idx<- rep(FALSE,nrow(ped))
      while(TRUE){
         idxTmp<- is.element(ped$id,idTmp)
         idx<- idx | idxTmp
         idx<- idx | is.element(ped$id,ped$father[idxTmp]) |
                     is.element(ped$id,ped$mother[idxTmp])
         idTmp<- ped$id[idx]

         if(sum(idxTmp)==sum(idx)) break
      }
      ped<- ped[idx,]
      rm(idTmp, idx, idxTmp)
   }

   if(is.null(ped$generation))
      ped<- pedRecode.0(ped,msg)
   ped$generation<- trimws(ped$generation)
   ids<- paste(ped$generation,ped$id,sep="~")
   uids<- unique(ids)
      idx<- match(uids,ids)
   if(length(uids)<length(ids) && msg){
      cat("   The following samples are repeated:\a\n")
      print(ped[-idx,c("id","generation","father","mother")])
      cat("   Repeated IDs were excluded!\a\n")
   }
   ped<- ped[idx,]
   rm(idx, ids)

   # new code
   ped$generation<- reorder(factor(ped$generation))
   ped$id<- reorder(factor(ped$id))
   ped$father<- reorder(factor(ped$father))
   ped$mother<- reorder(factor(ped$mother))
   ped<- ped[order(ped$generation,ped$id,ped$father,ped$mother),]
   idd<- data.frame(index=c(0,0),id=c(NA,0))
      idd<- rbind(idd,data.frame(index=1:nrow(ped),id=ped$id))
   ### recode here

   # recode IDs
   if(is.null(ped$sex)){
      idx<- match(sapply(ped$father,as.character),idd$id)
         idx[is.na(idx)]<- 1
      father<- idd$index[idx]
         str<- strf(ped,"father")
         if(any(str$idx))
            father[str$idx]<- str$val-10000
      idx<- match(sapply(ped$mother,as.character),idd$id)
         idx[is.na(idx)]<- 1
      mother<- idd$index[idx]
         str<- strf(ped,"mother")
         if(any(str$idx))
            mother[str$idx]<- str$val-20000
      ped<- data.frame(id=idd$index[-c(1:2)],
                       father=father,
                       mother=mother,
                       generation=ped$generation,
                       old.id=ped$id)
   }else{
      idx<- match(sapply(ped$father,as.character),idd$id)
         idx[is.na(idx)]<- 1
      father<- idd$index[idx]
         str<- strf(ped,"father")
         if(any(str$idx))
            father[str$idx]<- str$val-10000
      idx<- match(sapply(ped$mother,as.character),idd$id)
         idx[is.na(idx)]<- 1
      mother<- idd$index[idx]
         str<- strf(ped,"mother")
         if(any(str$idx))
            mother[str$idx]<- str$val-20000
      ped<- data.frame(id=idd$index[-c(1:2)],
                       father=father,
                       mother=mother,
                       sex=trimws(ped$sex),
                       generation=ped$generation,
                       old.id=ped$id)
      rm(idx, father, mother)
      ii<- match(sapply(ped$father,as.character),ped$id)
         ii<- ii[!is.na(ii)]
      if(length(ii)>0){
         idx<- !is.element(ped$sex[ii], c("0", "1", "M", "Male"))
            idx<- unique(ii[idx])
         if(any(idx) && msg){
            cat("   Suppose 1, M or Male stands for male...\n")
            cat("   --------------------------------------\n")
            print(ped[idx,])
            cat("   --------------------------------------\n")
            cat("   Above should be male(s)...\a\n\n")
         }
         rm(idx)
      }
      rm(ii)
      jj<- match(sapply(ped$mother,as.character),ped$id)
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
            cat("   Above should be female(s)...\a\n\n")
         }
         rm(idx)
      }
      rm(jj)
      fms<- intersect(sapply(ped$father, as.character), sapply(ped$mother, as.character))
         fms<- setdiff(fms, c(0,NA))
      if(length(fms) > 0){
         print(fms)
         cat("   Above are both father and mother?\a\n\n")
      }
   }
   idx<- (ped$father > ped$id) | (ped$mother > ped$id)
   if(any(idx)){
      idx<- match(sapply(ped$old[idx],as.character),pedSave$id)
      print(pedSave[idx,][1:min(sum(idx),3),])
      cat("... ...\n")
      stop("Check the above for errors...", call.=FALSE)
   }
   idx.s<- ped$generation[match(sapply(ped$father,as.character), ped$id)] == ped$generation
      idx.s<- (1:length(idx.s))[idx.s]
      idx.s<- idx.s[!is.na(idx.s)]
   idx.d<- ped$generation[match(sapply(ped$mother,as.character), ped$id)] == ped$generation
      idx.d<- (1:length(idx.d))[idx.d]
      idx.d<- idx.d[!is.na(idx.d)]
   if(length(idx.s)>0 || length(idx.d)>0){
      idx<- c(idx.s, idx.d)
         idx<- sort(unique(idx))
      idx<- c(sort(unique(c(ped$father[idx.s], ped$mother[idx.d]))),idx)
      pedTmp<- ped[idx,]
         pedTmp$id<- ped$old[match(sapply(pedTmp$id,as.character),ped$id)]
         pedTmp$father<- ped$old[match(sapply(pedTmp$father,as.character),ped$id)]
         pedTmp$mother<- ped$old[match(sapply(pedTmp$mother,as.character),ped$id)]
      pedTmp$old.id<- NULL
      print(pedTmp)
      stop("Check the above for errors regarding generations...", call.=FALSE)
   }
   rownames(ped)<- 1:nrow(ped)

   ped
}

# create "generation"
pedRecode.0<- function(ped,msg){
# ped: data frame (id,father,mother,...)
   ped<- as.data.frame(ped)
      ped$generation<- NULL
   if(is.null(ped$id)){
      stop("'id' missing...", call.=FALSE)
   }else ped$id<- trimws(ped$id)
   if(is.null(ped$father)){
      stop("'father/sire' missing...", call.=FALSE)
   }else ped$father<- trimws(ped$father)
   if(is.null(ped$mother)){
      stop("'mother/dam' missing...", call.=FALSE)
   }else ped$mother<- trimws(ped$mother)

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
      ii0<- is.element(ped$id,ped$father[ii0]) |
            is.element(ped$id,ped$mother[ii0])
      if(any(ii0)){
         ii<- rbind(ii0,ii)
      }else break
   }
# no offspring or parents
#   idx<- is.element(ped$id,ped$father) |
#         is.element(ped$id,ped$mother)  |
#         is.element(ped$father,ped$id) |
#         is.element(ped$mother,ped$id)
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
   out<- cbind(id=out[,jj],generation=out$generation,out[,-c(1,jj)])
      idx<- (is.na(out$father) | out$father == 0) & 
            (is.na(out$mother) | out$mother == 0)
      if(any(idx)) out$generation[idx]<- 0
      rm(idx)

   out$generation<- reorder(factor(out$generation))
   out$id<- reorder(factor(out$id))
   out$father<- reorder(factor(out$father))
   out$mother<- reorder(factor(out$mother))

   out<- out[order(out$generation,out$id,out$father,out$mother),]
   rownames(out)<- NULL
   out
}

