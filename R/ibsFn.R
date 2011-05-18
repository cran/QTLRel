
########################################
# genetic matrices from genotypic data #
#           identity by state          #
########################################

# ksp, delta1, delta2,delta3+delta5,delta7
genMatrix.default<- function(x){
# x: genotype data ("AA","AB","BB", or, 1,2,3); no missing data
# rows represent individuals, columns represent SNPs
   gdata<- as.matrix(x)
   if(!is.numeric(gdata))
      gdata<- (gdata=="AA")*1 +
              (gdata=="AB")*2 +
              (gdata=="BB")*3
   if(any(!is.element(unique(c(gdata)),c(1,2,3))))
      stop("x: wrong input...")
   nr<- nrow(gdata)
   nc<- ncol(gdata)
   npairs<- nr*(nr+1)/2

   deltac<- matrix(-1,nrow=npairs,ncol=5)
   out<- .C("deltaFnc",
            x=as.integer(t(gdata)),
            nr=as.integer(nr),
            nc=as.integer(nc),
            deltac=as.double(t(deltac)))

   delta<- matrix(out$deltac,ncol=5,byrow=TRUE)
      colnames(delta)<- c("ksp","delta1","delta2","delta35","delta7")

   ids<- rownames(x)
   if(is.null(ids)){
      ids<- 1:nrow(x)
   }else{
      ids<- trim(as.character(ids))
   }
   ksp<- matrix(-99,nrow=nr,ncol=nr)
      rownames(ksp)<- colnames(ksp)<- ids
   dlt1<- dlt2<- dlt35<- dlt7<- ksp
   ii<- 0
   for(i in 1:nr){
      for(j in i:nr){
         ii<- ii+1
         ksp[i,j]<- ksp[j,i]<- delta[ii,1]
         dlt1[i,j]<- dlt1[j,i]<- delta[ii,2]
         dlt2[i,j]<- dlt2[j,i]<- delta[ii,3]
         dlt35[i,j]<- dlt35[j,i]<- delta[ii,4]
         dlt7[i,j]<- dlt7[j,i]<- delta[ii,5]
      }
   }

   ib<- 2*diag(ksp)-1
      names(ib)<-  ids
   AA<- 2*ksp
      rownames(AA)<- colnames(AA)<- ids
   DD<- dlt7
      rownames(DD)<- colnames(DD)<- ids
   AD<- 4*dlt1 + dlt35
      rownames(AD)<- colnames(AD)<- ids
   HH<- dlt1
      rownames(HH)<- colnames(HH)<- ids
   MH<- dlt1 + dlt2 - ib%o%ib
      rownames(MH)<- colnames(MH)<- ids

   list(ib=ib,
        AA=AA,
        DD=DD,
        AD=AD,
        HH=HH,
        MH=MH)
}

genMatrix<- function(x){
   UseMethod("genMatrix")
}


