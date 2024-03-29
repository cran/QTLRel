
plot.scanOne<- function(x, ...){
   plot.lrt<- function(lrt, cv, ...){
   # lrt: data.frame(y,chr,dist,...)
      lrt$chr<- reorder(factor(lrt$chr))
      lrt<- lrt[order(lrt$chr,lrt$dist),]

      chr<- unique(lrt$chr)
         nchr<- length(chr)
      for(i in 1:nchr){
         idx<- lrt$chr==chr[i]
         lrt$dist[idx]<- diff(c(0,lrt$dist[idx]))
      }
      lrt$dist<- cumsum(lrt$dist)

      plot(range(lrt$dist),range(lrt$y),type="n",xaxt="n",...)

      points(lrt$dist,lrt$y,type="p",...)
      idx<- c(TRUE,lrt$chr[-length(lrt$chr)]!=lrt$chr[-1])
      min.p<- lrt$dist[idx]
      min.p<- c(min.p,max(lrt$dist))

      if(!is.null(cv)) abline(h=cv,col=2,lty=4)
      abline(v=min.p,lty=3,lwd=0.1)
      pos<- (min.p[-1]+min.p[-c(nchr+1)])/2
      axis(1,at=pos,labels=chr,tick=T,las=2)
   }

   xTmp<- list(...)
   if(is.null(xTmp$cex)) xTmp$cex<- 0.5
   if(is.null(xTmp$main)) xTmp$main<- ""
   if(is.null(xTmp$xlab)) xTmp$xlab<- "Chromosome"
   if(!is.null(xTmp$ylab)){
      cat("   Note: ylab supplied but ignored...\n")
   }

   if(is.null(x$LRT)) lrt<- data.frame(y=x$pval) else
      lrt<- data.frame(y=x$LRT)
   if(is.null(x$chr)){
      if(is.null(xTmp$gmap)) stop("Need: gmap...", call.=FALSE)
      gmap<- xTmp$gmap
      idx<- match(rownames(lrt),gmap$snp)
      lrt$chr<- gmap$chr[idx]
      lrt$dist<- gmap$dist[idx]
   }else{
      lrt$chr=x$chr
      lrt$dist=x$dist
   }
   if(is.null(xTmp$col)){
      col<- NULL
      chr<- unique(lrt$chr)
         nchr<- length(chr)
      for(i in 1:nchr){
         idx<- lrt$chr==chr[i]
         col<- c(col,rep(i%%2 + 3,sum(idx)))
      }
   }else{
      col<- xTmp$col
   }
   if(is.element("None",class(x))){
      lrt$y<- x$LRT/(2*log(10))
      plot.lrt(lrt, cv=xTmp$cv, cex=xTmp$cex, col=col, ylim=xTmp$ylim, main=xTmp$main,
         xlab=xTmp$xlab, ylab="LOD")
   }else{
      lrt$y<- -log10(x$pval)
      plot.lrt(lrt, cv=xTmp$cv, cex=xTmp$cex, col=col, ylim=xTmp$ylim, main=xTmp$main,
         xlab=xTmp$xlab, ylab=expression(paste(-log[10],"(p-value)")))
   }
}

plotit<- function(lrt, cv, bychr=FALSE, chr.labels=TRUE,
   type="p", lty=NULL, col=NULL, pch=NULL, cex=NULL, ...){
# lrt: data.frame(y,chr,dist,group,...)
   hh<- NULL
   if(!missing(cv) && !is.null(cv)) hh<- cv

   if(bychr){
      groups<- NULL
      if(!is.null(lrt$group)){
         groups=lrt$group
      }

      if(length(unique(lrt$chr))>1){
         lrt$chr<- reorder(factor(lrt$chr))
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
      lrt$chr<- reorder(factor(lrt$chr))
      lrt$group<- reorder(factor(lrt$group))
      lrt<- lrt[order(lrt$group,lrt$chr,lrt$dist),]

      groups<- lrt$group
         groups<- sort(unique(groups))
      ngr<- length(groups)
      if(ngr>1) cat("  Groups: ",as.character(groups),"\n", sep="")
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
         if(is.null(lty)){
            lty<- 1:5
            lty<- rep(lty,ngr)[1:ngr]
         }
         if(is.null(col)){
            col<- 1:ngr
         }
         pch<- rep(pch,ngr); pch<- pch[1:ngr]
         cex<- rep(cex,ngr); cex<- cex[1:ngr]
      }else{
         if(is.null(col)) col<- 1
      }
      if(length(col) < length(lrt$dist))
         colTmp<- col[match(lrt$group, groups)]
      else colTmp<- col

      min.p<- matrix(NA,nrow=ngr,ncol=nchr)
      for(g in 1:ngr){
         idx0<- lrt$group==groups[g]
         lrt0<- lrt[idx0,]
         col0<- colTmp[idx0]
         for(i in 1:nchr){
            idx<- lrt0$chr==chr[i]
            lines(lrt0$dist[idx],lrt0$y[idx],type=type,lty=lty[g],col=col0[idx],pch=pch[g],cex=cex[g],...)
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

plot.scanTwo<- function(x, ...){
# x: object of scanTwo
# a genetic map 'gmap' is needed
   lst<- list(...)
   if(is.null(lst$gmap))
      stop("Need a genetic map 'gmap'.", call.=FALSE)
   qqint<- x
   gmap<- lst$gmap
   v<- as.matrix(qqint)
   rv<- range(v,na.rm=TRUE)
      rv<- seq(floor(rv[1]),ceiling(rv[2]),by=0.5)
   dst<- gmap$dist
      dst<- diff(dst)
      dst[dst<0]<- 0
      dst[dst==0]<- 1e-5
      dst<- c(0,dst)
      dst<- cumsum(dst)
   chrs<- unique(gmap$chr)
   xat<- NULL
   for(chr in chrs){
      xat<- c(xat,mean(range(dst[gmap$chr==chr])))
   }

   if(is.null(lst$col)){
      col<- terrain.colors(12)
   }else{
      col<- lst$col
   }
   scale<- max(dst)/(max(rv)-min(rv))*0.75
   image(x=dst,y=dst,z=v,axes=F,xlab=lst$xlab,xlim=c(-2,max(dst)),
      ylim=c(0,max(dst)+2),ylab=lst$ylab,
      main=lst$main,col=col)
   image(x=max(dst)*c(19/20,1),y=(rv-min(rv))*scale,z=matrix(rv,nrow=1),
      col=col,add=TRUE)
   axis(4,at=(rv-min(rv))*scale,labels=rv,pos=par("usr")[2]*0.985,tick=FALSE)
   axis(2,at=xat,labels=chrs,tick=F,line=-0.75)
   axis(3,at=xat,labels=chrs,tick=F,line=-0.75)
   mtext("Chromosomes",2,line=1.5)
   mtext("Chromosomes",3,line=1.25,at=0)
   col<- 0
   for(chr in chrs){
      col<- col+1
      lines(c(-3,-3),range(dst[gmap$chr==chr]),col=col,lwd=5)
      lines(range(dst[gmap$chr==chr]),max(dst)-c(-3,-3),
      col=col,lwd=5)
   }
}

