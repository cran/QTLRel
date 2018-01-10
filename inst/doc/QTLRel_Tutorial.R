

#
# NOTE: working folder: "QTLRel Tutorial/"
#

# create inbred line by selfing
slf<- function(ped, id, generations=99){
   pd<- cbind(id=character(0), sex=character(0), generation=character(0),sire=character(0), dam=character(0))
   idTmp1<- idTmp2<- id
   cnt<- 1
   while(cnt < generations){
      idTmp2<- paste(id,"EXT",cnt,sep="")
      pd<- rbind(
         cbind(id=idTmp1, sex=0, generation=-cnt, sire=idTmp2, dam=idTmp2),
         pd
      )
      idTmp1<- idTmp2
      cnt<- cnt+1
   }
   pd<- rbind(
      cbind(id=idTmp1, sex=0, generation=-cnt, sire=0, dam=0),
      pd
   )
   pd
}

# create inbred line by brother-sister mating
bsf<- function(ped, id, generations=99){
   sx<- sort(unique(sapply(ped$sex,as.character)), decreasing=TRUE)
   pd<- cbind(id=character(0), sex=character(0), generation=character(0),sire=character(0), dam=character(0))
   for(id in id){
      cnt<- 1
      ft1<- ft2<- paste(id,"EXTf",0,sep="")
      mt1<- mt2<- paste(id,"EXTm",0,sep="")
      pd<- rbind(
         cbind(id=id, sex=ifelse(is.element(id,ped$dam),sx[1],sx[2]), generation=-cnt, sire=ft2, dam=mt2),
         pd
      )
      cnt<- cnt+1
      while(cnt < generations){
         ft2<- unique(paste(id,"EXTf",cnt,sep=""))
         mt2<- unique(paste(id,"EXTm",cnt,sep=""))
         pd<- rbind(
            cbind(id=ft1, sex=sx[1], generation=-cnt, sire=ft2, dam=mt2),
            cbind(id=mt1, sex=sx[2], generation=-cnt, sire=ft2, dam=mt2),
            pd
         )
         ft1<- ft2
         mt1<- mt2
         cnt<- cnt+1
      }
      pd<- rbind(
         cbind(id=ft1, sex=sx[1], generation=-cnt, sire=0, dam=0),
         cbind(id=mt1, sex=sx[2], generation=-cnt, sire=0, dam=0),
         pd
      )
   }
   pd<- as.data.frame(pd)
      pd$sex<- reorder(pd$sex)
      pd$generation<- reorder(pd$generation)
   pd<- pd[order(pd$generation,pd$sex),]
   pd
}

pedF8<- read.table("data/pedigree.txt", header=TRUE, check.names=FALSE)
   pedF8<- pedF8[-c(1:2),]
   pedF8$sire[1:2]<- "1i"
   pedF8$dam[1:2]<- "2i"
rownames(pedF8)<- NULL
pedF8.1<- pedF8
   pedF8.1$sire[pedF8.1$gen=="F2"]<- 1
   pedF8.1$dam[pedF8.1$gen=="F2"]<- 1
   pedF8.1<- rbind(c(1,0,"F1",0,0,"F1"),as.matrix(pedF8.1[-(1:2),]))
   pedF8.1<- as.data.frame(pedF8.1)
rownames(pedF8.1)<- NULL
pedF8.2<- pedF8
   pedF8.2$sire[1:2]<- "fd.001"
   pedF8.2$dam[1:2]<- "fd.002"
   pedF8.2<- rbind(
         cbind(bsf(pedF8.2,c("fd.001", "fd.002"),150),family="NA"),
         as.matrix(pedF8.2)
   )
   pedF8.2<- as.data.frame(pedF8.2)
rownames(pedF8.2)<- NULL
gmapF8<- read.csv("data/genetic map F8.csv", header=TRUE, check.names=FALSE)
   colnames(gmapF8)<- c("snp","chr","dist","phyPos37")
gdatF8<- read.csv("data/genotype F8.csv", header=TRUE, check.names=FALSE)
   rownames(gdatF8)<- sapply(gdatF8$id,as.character)
   gmapF8<- gmapF8[is.element(gmapF8$snp,intersect(colnames(gdatF8),gmapF8$snp)),]
   idx<- is.element(colnames(gdatF8),intersect(colnames(gdatF8),gmapF8$snp))
   gdatF8<- gdatF8[,idx]
pdatF8<- read.table("data/phenotype F8.txt", sep="\t", header=TRUE, check.names=FALSE)
sum(!is.element(pdatF8$id,pedF8$id))
sum(!is.element(gdatF8$id,pedF8$id))
ids<- sapply(pdatF8$id[is.element(pdatF8$id,pedF8$id)], as.character)
id<- sample(ids, 500) # randomly choose 500 individuals

gdatF8<- gdatF8[match(id,rownames(gdatF8)),]
pdatF8<- pdatF8[match(id,pdatF8$id),]
   pdatF8<- pdatF8[,-1]
   rownames(pdatF8)<- id
pdatF8.<- pdatF8[,c("Sex","FCAge", "WeightD1", "Cage ID")]
colnames(pdatF8.)<- c("sex","age","bwt","cage")
pdatF8.[1:5,]
pdatF8<- pdatF8.

save(gdatF8, pdatF8, gmapF8, pedF8, pedF8.1, pedF8.2, file="QTLRelEx.RData")

#########################
# identity coefficients #
#########################

library(QTLRel)

load("QTLRelEx.RData")

# identical kinship
ksp<- kinship(pedF8, ids=rownames(pdatF8), all=TRUE, msg=FALSE)
ksp.1<- kinship(pedF8.1, ids=rownames(pdatF8), all=TRUE, msg=FALSE)
ksp.2<- kinship(pedF8.2, ids=rownames(pdatF8), all=TRUE, msg=FALSE)

dim(ksp) #500 500

sum(abs(ksp-ksp.1)) # 0
sum(abs(ksp-ksp.2)) # 6.869755e-10
sum(abs(ksp.1-ksp.2)) # 6.869755e-10

# only interested individuals in F8
id<- rownames(pdatF8)
id[1:5]
# check if phenotype and genotype data are from the same sample
sum(!is.element(id,pedF8$id))
sum(!is.element(id,rownames(gdatF8)))
sum(!is.element(rownames(gdatF8),id))

# running
idcf<- cic(pedF8,ids=id,df=3,ask=TRUE,msg=TRUE)
# extract genetic matrices
gmF8<- genMatrix(idcf); names(gmF8)
dim(gmF8$AA)
gmF8$AA[1:3,1:5]

max(abs(ksp - gmF8$AA/2))

save(gdatF8, pdatF8, gmapF8, pedF8, pedF8.1, pedF8.2, gmF8, file="QTLRelEx.RData")

#######################
# variance components #
#######################

idx<- !is.na(pdatF8[,"bwt"])
pdatTmp<- pdatF8[idx,] # remove missing values
gdatTmp<- gdatF8[match(rownames(pdatTmp),rownames(gdatF8)),] #matching

ii<- match(rownames(pdatTmp),rownames(gmF8$AA))
vc<- estVC(y = pdatTmp[,"bwt"], x = pdatTmp[,c("sex","age")],
   v=list(AA=gmF8$AA[ii,ii],DD=gmF8$DD[ii,ii]))
vc$value # likelihood
vc$par # parameter estimates

ranMtr<- rem(~cage, data=pdatTmp) # matrice for random effect (cage)
names(ranMtr)
dim(ranMtr$cage)
vc.cage<- estVC(y = pdatTmp[,"bwt"], x = pdatTmp[,c("sex","age")],
   v=list(AA=gmF8$AA[ii,ii],DD=gmF8$DD[ii,ii], cage=ranMtr$cage))

2*(vc.cage$value - vc$value) # likelihood ratio test for cage effect

###############
# genome scan #
###############

sum(is.na(gdatTmp)) # number of missing genotypes
gdatTmpImputed<- genoImpute(gdatTmp,gmapF8,gr=8,na.str=NA)

# genome scan
lrt<- scanOne(y=pdatTmp[,"bwt"], x=pdatTmp[,c("sex","age")],
    gdat=gdatTmpImputed, vc=vc)
pdf("bwt1.pdf")
   plot(lrt,gmap=gmapF8,main="Body Weight") # plotting
dev.off()

# genome scan: interactive sex
lrt.sex<- scanOne(y=pdatTmp[,"bwt"], x=pdatTmp[,c("age")],
    gdat=gdatTmpImputed, intcovar=pdatTmp[,c("sex")],vc=vc)
pdf("bwt2.pdf")
  plot(lrt.sex,gmap=gmapF8,main="Body Weight") # plotting
dev.off()

# QTL by sex interaction
lrtTmp<- lrt
   lrtTmp$p<- lrt.sex$p - lrt$p
pdf("bwt3.pdf")
   plot(lrtTmp,gmap=gmapF8,main="Body Weight: QTL by Sex Effect") # plotting
dev.off()

# Haley-Knott interval mapping
gdTmp<- gdatTmp
   gdTmp[is.na(gdTmp)]<- 0
   unique(c(as.matrix(gdTmp)))
prDat<- genoProb(gdat=gdTmp,gmap=gmapF8,step=3,gr=8)

# genome scan: Haley-Knott method
lrtHK<- scanOne(y=pdatTmp[,"bwt"], x=pdatTmp[,c("sex","age")],
    prd=prDat, vc=vc)
pdf("bwt4.pdf")
   plot(lrtHK,main="Body Weight (HK Method)") # plotting
dev.off()

# genome scan: Haley-Knott method, interactive sex
lrtHK.sex<- scanOne(y=pdatTmp[,"bwt"], x=pdatTmp[,c("age")],
    prd=prDat, intcovar=pdatTmp[,c("sex")],vc=vc)
pdf("bwt5.pdf")
   plot(lrtHK.sex,main="Body Weight (HK Method)") # plotting
dev.off()

# Haley-Knott method, QTL by sex interaction
lrtHKTmp<- lrtHK
   lrtHKTmp$p<- lrtHK.sex$p - lrtHK$p
pdf("bwt6.pdf")
   plot(lrtHKTmp,main="Body Weight (HK Method): QTL by Sex Effect") # plotting
dev.off()

##############
# thresholds #
##############

nn<- nrow(gdatTmpImputed) # sample size
ntimes<- 1000 # number of simulations
cvMtr<- NULL # matrix to save results of permuted data
for(n in 1:ntimes){
   idx<- sample(1:nn,replace=FALSE) # permutation
   tmp<- scanOne(y=pdatTmp[,"bwt"],x=pdatTmp[,c("sex","age")],
       gdat=gdatTmpImputed[idx,], vc=vc)
   cvMtr<- rbind(cvMtr,tmp$p)
   cat(n,"/",ntimes,"\r") # track process
}

# 0.05 thresholds
quantile(apply(cvMtr,1,max),0.95)
quantile(apply(cvMtr,1,max),0.95)/(2*log(10))

# for HK method...
cvMtrHK<- NULL
for(n in 1:ntimes){
   idx<- sample(1:nn,replace=FALSE)
   prdTmp<- prDat
   prdTmp$pr<- prdTmp$pr[idx,,]
   tmp<- scanOne(y=pdatTmp[,"bwt"], x=pdatTmp[,c("sex","age")],
       prd=prdTmp, vc=vc)
   cvMtrHK<- rbind(cvMtrHK,tmp$p)
   cat(n,"/",ntimes,"\r")
}

# gene dropping
pedR<- pedRecode(pedF8) # recode the pedigree
ids<- rownames(gdatTmpImputed) # relevant individual IDs
cvMtrGD<- NULL
for(n in 1:ntimes){
   gdatTmp<- genoSim(pedR, gmapF8, ids=ids)
   tmp<- scanOne(y=pdatTmp[,"bwt"],x=pdatTmp[,c("sex","age")],
       gdat=gdatTmp, vc=vc)
   cvMtrGD<- rbind(cvMtrGD,tmp$p)
   cat(n,"/",ntimes,"\r")
}

############
# plotting #
############

pdf("bwt7.pdf")
   plot(lrt,cv=3.2,gmap=gmapF8,main="Body Weight",cex=1)
dev.off()

idx<- match(colnames(gdatTmpImputed),gmapF8$snp)
Tmp<- data.frame(chr=gmapF8$chr[idx],
                 dist=gmapF8$dist[idx],
                 y=lrt$p)
   Tmp<- Tmp[order(Tmp$chr,Tmp$dist),] # order by chromosome and distance
pdf("bwt8.pdf")
plotit(Tmp, cv=15, main="Mapping Plot of Body Weight", xlab="Chromosome",
   ylab="LRT", col=as.integer(Tmp$ch)%%2+2,type="p",lty=2)
dev.off()

library(lattice) # dependency package
pdf("bwt9.pdf")
tmp<- plotit(Tmp,cv=12,type="p",lty=2,col=4,cex=0.65,
      xlab="Genetic Position (cM)",ylab="LRT",
      main="Body Weight",bychr=TRUE)
print(tmp)
dev.off()

# save.image("ttmp.RData")

#################
# miscellaneous #
#################

# load("ttmp.RData")

# lod ci
Tmp<- data.frame(chr=lrtHK$chr,
                 dist=lrtHK$dist,
                 y=lrtHK$p/(2*log(10))) # convert to LOD
   Tmp$chr<- reorder(Tmp$chr)
   Tmp<- Tmp[order(Tmp$chr,Tmp$dist),] # order by chromosome and distance
lc<- lodci(Tmp,cv=3.2,lod=1.5,drop=1.5)
lc

# multiple QTL model estimates
dtfTmp<- data.frame(
   y=pdatTmp[,"bwt"],
   age=pdatTmp[,"age"],
   sex=pdatTmp[,"sex"],
   a1=prDat$pr[,1,lc$index[1]]-prDat$pr[,3,lc$index[1]],
   d1=prDat$pr[,2,lc$index[1]],
   a2=prDat$pr[,1,lc$index[2]]-prDat$pr[,3,lc$index[2]],
   d2=prDat$pr[,2,lc$index[2]]
)
est1<- gls(y~age+sex+a1+d1+a2+d2,data=dtfTmp,vc=vc)
est2<- gls(y~age+sex*(a1+d1+a2+d2),data=dtfTmp,vc=vc)
est1 # with age and sex being additive covariates
est2 # with sex being an interactive covariate

# herit
ii<- match(rownames(pdatTmp),rownames(gmF8$AA))
vc0<- estVC(y = pdatTmp[,"bwt"], v=list(AA=gmF8$AA[ii,ii], DD=gmF8$DD[ii,ii]))
nb<- length(vc0$par) - length(vc0$v)
nr<- nrow(vc0$y)
cov<- matrix(0,nrow=nr,ncol=nr)
for(i in 1:length(vc0$v))
   cov<- cov + vc0$v[[i]]*vc0$par[nb+i]
tv<- mean(diag(cov)) # total variation
eff<- NULL # QTL effects
for(n in 1:length(lrtHK$par)){
   eff<- rbind(eff,lrtHK$par[[n]][c("a","d")])
}
eff<- data.frame(eff) # data frame!
qv<- qtlVar(eff,prDat$pr) # per QTL variation
qv[lc$index]/tv*100 # per QTL heritability

# QTLxQTL interaction
date()
qqInt.1<- scanTwo(y=pdatTmp[,"bwt"], x=pdatTmp[,c("sex","age")],
    gdat=gdatTmpImputed, vc=vc)
date()
qqInt.2<- scanTwo(y=pdatTmp[,"bwt"], x=pdatTmp[,c("sex","age")],
    prdat=prDat, vc=vc)
date()
# save.image("ttmp.RData")

#plotting 
png("bwt_qqint.png",width=640)
par(mfrow=c(1,2))
plot(qqInt.1/(2*log(10)),
   gmap=gmapF8[match(rownames(qqInt.1),gmapF8$snp),],
   main="Body Weight: Epistasis (A)\n\n",xlab="",ylab="")
plot(qqInt.2/(2*log(10)),
   gmap=data.frame(snp=prDat$snp,chr=prDat$chr,dist=prDat$dist),
   main="Body Weight: Epistasis (B)\n\n",xlab="",ylab="")
dev.off()


###################################################
# the end #
###########

