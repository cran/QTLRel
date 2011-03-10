
plot.scanOne<- function(x,...){
   xTmp<- list(...)
   if(is.null(xTmp$cex)) cex<- 0.3
   if(is.null(xTmp$main)) main<- ""
   if(is.null(xTmp$xlab)) xlab<- "Chromosome"
   cv<- xTmp$cv

   lrt<- data.frame(y=x$p)
   if(is.null(x$chr)){
      if(is.null(xTmp$gmap)) stop("need: gmap...")
      gmap<- xTmp$gmap
      idx<- match(names(x$p),gmap$snp)
      lrt$chr<- gmap$chr[idx]
      lrt$dist<- gmap$dist[idx]
   }else{
      lrt$chr=x$chr
      lrt$dist=x$dist
   }
   if(is.element("None",class(x))){
      lrt$y<- x$p/(2*log(10))
      plot.lrt(lrt,cv,cex=cex,main=main,xlab=xlab,ylab="LOD")
   }else{
      lrt$y<- -log10(x$p)
      plot.lrt(lrt,cv,cex=cex,main=main,xlab=xlab,ylab=expression(paste(-log[10],"(p-value)")))
   }
}

plot.lrt<- function(lrt,cv,...){
# lrt: data.frame(y,chr,dist,...)
   lrt$chr<- as.factor(lrt$chr)
      lrt$chr<- reorder(lrt$chr)
   lrt<- lrt[order(lrt$chr,lrt$dist),]

   chr<- unique(lrt$chr)
      nchr<- length(chr)
   for(i in 1:nchr){
      idx<- lrt$chr==chr[i]
      lrt$dist[idx]<- diff(c(0,lrt$dist[idx]))
   }
   lrt$dist<- cumsum(lrt$dist)

   plot(range(lrt$dist),range(lrt$y),type="n",xaxt="n",...)

   for(i in 1:nchr){
      idx<- lrt$chr==chr[i]
      lines(lrt$dist[idx],lrt$y[idx],type="p",col=i%%2 + 3,...)
   }
   idx<- c(TRUE,lrt$chr[-length(lrt$chr)]!=lrt$chr[-1])
   min.p<- lrt$dist[idx]
   min.p<- c(min.p,max(lrt$dist))

   if(!missing(cv)) abline(h=cv,col=2,lty=4,...)
   abline(v=min.p,lty=3,lwd=0.1)
   pos<- (min.p[-1]+min.p[-c(nchr+1)])/2
   axis(1,at=pos,labels=chr,tick=T,las=2)
}


plotit<- function(lrt,cv,bychr=FALSE,chr.labels=TRUE,
   type="p",lty=NULL,col=NULL,pch=NULL,cex=NULL,...){
# lrt: data.frame(y,chr,dist,group,...)
   hh<- NULL
   if(!missing(cv) && !is.null(cv)) hh<- cv

   if(bychr){
      groups<- NULL
      if(!is.null(lrt$group)){
         groups=lrt$group
      }

      if(length(unique(lrt$chr))>1){
         xyplot(y~dist|chr,data=lrt,
            groups=groups,
            panel=function(x,y,...){
               panel.xyplot(x,y,...)
               if(!is.null(hh)) panel.abline(h=hh,...)
               if(!is.null(lrt$group)) panel.superpose(x,y,...)
            },
            type=type,
            lty=lty,
            col=col,
            pch=if(!is.null(pch)) pch else 1,
            cex=cex,
            ...
         )
      }else{
         xyplot(y~dist,data=lrt,
            groups=groups,
            panel=function(x,y,...){
               panel.xyplot(x,y,...)
               if(!is.null(hh)) panel.abline(h=hh,...)
               if(!is.null(lrt$group)) panel.superpose(x,y,...)
            },
            type=type,
            lty=lty,
            col=col,
            pch=if(!is.null(pch)) pch else 1,
            cex=cex,
            ...
         )
      }
   }else{
      if(is.null(lrt$group)) lrt$group<- 1
      lrt<- as.data.frame(lrt)
      lrt$chr<- as.factor(lrt$chr)
         lrt$chr<- reorder(lrt$chr)
      lrt$group<- as.factor(lrt$group)
         lrt$group<- reorder(lrt$group)
      lrt<- lrt[order(lrt$group,lrt$chr,lrt$dist),]

      groups<- lrt$group
         groups<- sort(unique(groups))
      ngr<- length(groups)
      if(ngr>1) cat("  Groups:",as.character(groups),"\n")
      chr<- unique(lrt$chr)
         nchr<- length(chr)
      for(i in 1:ngr){
         idx<- lrt$group==groups[i]
         lrt0<- lrt[idx,]
         chr0<- unique(lrt0$chr)
            nchr0<- length(chr0)
         for(i in 1:nchr0){
            idx0<- lrt0$chr==chr0[i]
            lrt0$dist[idx0]<- diff(c(0,lrt0$dist[idx0]))
         }
         lrt0$dist<- cumsum(lrt0$dist)
         lrt$dist[idx]<- lrt0$dist
      }
      if(chr.labels){
         plot(range(lrt$dist),range(lrt$y),type="n",xaxt="n",...)
      }else{
         plot(range(lrt$dist),range(lrt$y),type="n",...)
      }

      if(ngr>1){
         if(is.null(lty)) lty<- 1:5; lty<- rep(lty,ngr)[1:ngr]
         if(is.null(col)) col<- 1:6; col<- rep(col,ngr)[1:ngr]
         pch<- rep(pch,ngr); pch<- pch[1:ngr]
         cex<- rep(cex,ngr); cex<- cex[1:ngr]
      }else if(is.null(col)) col<- rep(1,length(lrt$dist))
      colTmp<- matrix(col,ncol=ngr,byrow=FALSE)

      min.p<- matrix(NA,nrow=ngr,ncol=nchr)
      for(g in 1:ngr){
         idx0<- lrt$group==groups[g]
         lrt0<- lrt[idx0,]
         for(i in 1:nchr){
            idx<- lrt0$chr==chr[i]
            lines(lrt0$dist[idx],lrt0$y[idx],type=type,lty=lty[g],col=colTmp[idx,g],pch=pch[g],cex=cex[g],...)
         }
         idx<- c(TRUE,lrt0$chr[-length(lrt0$chr)]!=lrt0$chr[-1])
         min.p[g,]<- lrt0$dist[idx]
      }
      min.p<- apply(min.p,2,min)
      min.p<- c(min.p,max(lrt$dist))

      if(!missing(cv)) abline(h=cv,lty=lty,col=col)
      abline(v=min.p,lty=3,lwd=0.1)
      if(chr.labels){
         pos<- (min.p[-1]+min.p[-c(nchr+1)])/2
         axis(1,at=pos,labels=chr,tick=F)
      }
   }
}


