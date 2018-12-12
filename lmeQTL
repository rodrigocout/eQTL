################################################
### simple eQTL analysis
################################################
# for gxdata and gtypedata, see getGxGtypesCohort1 and getGxGtypesCohort2 functions!
# geneSnpList is a named list of SNPs (names are genes to be tested)
lmEQTL<-function(geneSnpList, gxdata, gtypedata, gxids=NULL, cutoff=NULL, dosort=TRUE, 
                 global=FALSE, setFormat=FALSE, verbose=TRUE){
  geneSnpList<-geneSnpList[lapply(geneSnpList, length)>0]
  markers<-unlist(geneSnpList,use.names = FALSE)
  if(ncol(gxdata)!=ncol(gtypedata))stop("NOTE: please provide a sample matched 'gxdata' and 'gtypedata'!")
  if(!all(colnames(gxdata)==colnames(gtypedata)))stop("NOTE: please provide a sample matched 'gxdata' and 'gtypedata'!")
  if(!all(markers%in%rownames(gtypedata))){
    stop("NOTE: all markers in 'geneSnpList' should be available in 'gtypedata'!")
  }
  genes<-names(geneSnpList)
  if(!all(genes%in%rownames(gxdata))){
    stop("NOTE: all genes in 'geneSnpList' should be available in 'gxdata'!")
  }
  if(!is.null(gxids)){
    if(!is.data.frame(gxids))stop("NOTE: 'gxids' should be a data frame object!")
    if(!all(genes%in%rownames(gxids)))stop("NOTE: all ENTREZ IDs in 'geneSnpList' should be listed in 'gxids' rownames! ")
  }
  names(geneSnpList)<-paste("ID",names(geneSnpList),sep="")
  rownames(gxdata)<-paste("ID",rownames(gxdata),sep="")
  names(genes)<-names(geneSnpList)
  #---
  gxgtype<-as.data.frame(cbind(t(gxdata),t(gtypedata)))
  if(verbose) pb <- txtProgressBar(style=3)
  if(global){
    eqtls<-NULL
    for(nm in names(genes)){
      if(verbose) setTxtProgressBar(pb, which(names(genes)==nm)/length(genes))
      fm<-formula(paste( nm, " ~ ", paste(geneSnpList[[nm]], collapse = "+") , sep="" ) )
      fit <- lm(fm, data=gxgtype)
      res<-summary(fit)
      cf <- coef(res)
      if(nrow(cf)>=2){
        res<-cf[-1,,drop=FALSE]
      } else {
        res<-c(NA,NA,NA,1) #.. when MAF is 0 for one marker's case!
      }
      res<-data.frame(Gene=as.character(genes[nm]), Marker=rownames(res), res, 
                      check.names = FALSE, stringsAsFactors = FALSE)
      eqtls<-rbind(eqtls,res)
    }
  } else {
    eqtls<-NULL
    for(nm in names(genes)){
      if(verbose) setTxtProgressBar(pb, which(names(genes)==nm)/length(genes))
      mks<-geneSnpList[[nm]]
      res<-t(sapply(mks, function(mk){
        fm <- formula(paste( nm, " ~ ", mk , sep="" ) )
        fit <- lm(fm, data=gxgtype)
        res <- summary(fit)
        cf <- coef(res)
        if(nrow(cf)==2){
          res<-cf[-1,]
        } else {
          res<-c(NA,NA,NA,1) #.. when MAF is 0!
        }
        res
      }))
      res<-data.frame(Gene=as.character(genes[nm]), Marker=rownames(res), res, check.names = FALSE, stringsAsFactors = FALSE)
      eqtls<-rbind(eqtls,res)
    }
  }
  if(verbose) close(pb)
  colnames(eqtls)<-c("Gene","Marker","Coef","stdError","tvalue","pvalue")
  rownames(eqtls)<-paste(eqtls$Gene," ~ ",eqtls$Marker, sep = "")
  if(dosort) eqtls<-eqtls[sort.list(eqtls[,"pvalue"]),]
  eqtls$AdjPvalue<-p.adjust(eqtls$pvalue, method = "BH")
  if(!is.null(cutoff)){
    eqtls<-eqtls[eqtls[,"AdjPvalue"]<cutoff,]
  }
  if(setFormat){
    eqtls$Coef<-round(eqtls$Coef, digits = 3)
    eqtls$stdError<-round(eqtls$stdError, digits = 3)
    eqtls$tvalue<-round(eqtls$tvalue, digits = 2)
    eqtls$pvalue<-format(eqtls$pvalue, scientific=TRUE, digits = 2)
    eqtls$AdjPvalue<-format(eqtls$AdjPvalue, scientific=TRUE, digits = 2)
  }
  if(!is.null(gxids)){
    res<-cbind(gxids[eqtls$Gene,],eqtls)
    rownames(res)<-rownames(eqtls)
    eqtls<-res
  }
  eqtls
}

################################################
### return a list with genes mapped to snps
################################################
gene2snpCandidates<-function(annotgene,annotsnps,maxgap=25e4, onlyCandicates=TRUE){
  #---keep candicates
  if(onlyCandicates){
    annotgene<-overlapGRangesModel(annotgene, annotsnps, maxgap=maxgap)
    annotsnps<-overlapGRangesModel(annotsnps, annotgene, maxgap=maxgap)
  }
  #---
  geneSnpList<-list()
  for(i in 1:length(annotgene)){
    gene_id<-annotgene@elementMetadata$gene_id[i]
    tp<-overlapGRangesModel(annotsnps, annotgene[i,], maxgap=maxgap)
    geneSnpList[[gene_id]]<-tp@elementMetadata$rsid
  }
  geneSnpList<-geneSnpList[lapply(geneSnpList, length)>0]
  geneSnpList
}

################################################
### simple overlaping analysis with GRanges objects
################################################
countRangesOverlap<-function(snpRanges, geneRanges, maxgap=500000){
  nov<-overlapGRangesModel(snpRanges, geneRanges, maxgap=maxgap)
  length(nov)/length(snpRanges)
}
overlapGRangesModel<-function(GR, GRmodel, maxgap=0){
  idx<-overlapsAny(GR,GRmodel,maxgap=maxgap)
  GR[idx,]
}
getTxDb<-function(EntrezIDs,getTSS=FALSE,removeDuplicates=FALSE){
  tr<-get("TxDb.Hsapiens.UCSC.hg19.knownGene")
  tr <- transcripts(tr, columns = c("gene_id","tx_name"))
  tr <- checkChrName(tr, addChr = TRUE)
  tr <- tr[elementNROWS(values(tr)$gene_id) > 0]
  names(tr) <- values(tr)$tx_name
  tss <- checkChrName(tr, addChr = TRUE)
  if(getTSS)tss <- flank(tr, width = -1, start = TRUE)
  tss$gene_id<-as.character(tss$gene_id)
  if(!missing(EntrezIDs))tss<-tss[as.logical(tss$gene_id%in%EntrezIDs)]
  if(removeDuplicates){
    tss<-as.data.frame(tss[,"gene_id"])
    start<-tss[,c("start","gene_id")]
    end<-tss[,c("end","gene_id")]
    start<-aggregate(. ~ gene_id, data=start, mean, na.rm=TRUE)
    end<-aggregate(. ~ gene_id, data=end, mean, na.rm=TRUE)
    tp<-data.frame(gene_id=start$gene_id,start=start$start,end=end$end,stringsAsFactors=FALSE)
    #--
    tss<-tss[!duplicated(tss$gene_id),]
    idx<-match(tss$gene_id,tp$gene_id)
    tss$start<-tp$start[idx]
    tss$end<-tp$end[idx]
    #--
    tss <- with(tss, GRanges(seqnames, IRanges(start,end,names=gene_id), strand=tss$strand))
    tss$gene_id<-names(tss)
    seqlengths(tss)<-seqlengths(tr)
  }
  tss
}

## add chr to the chromosome names
checkChrName <- function(grange, addChr=TRUE, keepUnsual=FALSE, keepSex=FALSE) {
  if (is(grange, 'GRanges')) {
    chrName <- seqlevels(grange)
  } else if (is(grange, 'character')) {
    chrName <- grange
  } else if (is(grange, 'GenoSet')) {
    chrName <- chrNames(grange)
  } else if (is(grange, 'RangedData')) {
    chrName <- levels(space(grange))
  } else {
    chrName <- names(grange)
    if (is.null(chrName)) 
      stop('Un-supported data types!')
  }
  if (any(grepl('^chr[0-9XY][0-9]?', chrName))) {
    if (!addChr) {
      chrName <- sub('^chr', '', chrName)
    }
  } else if (addChr) {
    ind <- grep('^[0-9XY][0-9]?', chrName)
    chrName[ind] <- paste('chr', chrName[ind], sep='')
  }
  if (is(grange, 'GRanges')) {
    seqlevels(grange) <- chrName
  } else if (is(grange, 'character')) {
    grange <- chrName
  } else if (is(grange, 'GenoSet')) {
    names(locData(grange)) <- chrName
  } else {
    names(grange) <- chrName
  }
  if(!keepUnsual){
    sx<-c("chrX","chrY")
    chr<-c(paste("chr",1:22,sep=""),if(keepSex) sx)
    idx<-seqnames(grange)%in%chr
    grange<-grange[idx,]
    seqlevels(grange)<-chr
  }
  return(grange)
}



########################################
### plotqt
########################################
ploteqtl<-function(geneid, snpid, gxdata, gtypedata, eqtls, dirname=".",
                   labpheno="", ylab="expression", ylim=NULL, plotpdf=FALSE){
  cid<-paste(geneid,snpid,sep=" ~ ")
  if(!cid%in%rownames(eqtls))stop("NOTE: please provide a valid #geneid!")
  if(ncol(gxdata)!=ncol(gtypedata))stop("NOTE: please provide a sample matched 'gxdata' and 'gtypedata'!")
  if(!all(colnames(gxdata)==colnames(gtypedata)))stop("NOTE: please provide a sample matched 'gxdata' and 'gtypedata'!")
  gexp<-gxdata[geneid, ]
  gtype<-gtypedata[snpid, ]
  #---
  if(is.null(ylim))ylim<-range(gexp, na.rm=TRUE)
  #---
  pval<-eqtls[cid,"pvalue"]
  if("SYMBOL"%in%colnames(eqtls)){
    genesymbol<-eqtls[cid,"SYMBOL"]
  } else {
    genesymbol<-geneid
  }
  dt<-data.frame(gtype=gtype,gexp=gexp)
  dt<-dt[complete.cases(dt),]
  if(plotpdf)pdf(file=paste(dirname,"/",genesymbol,"_",snpid,".pdf",sep=""), width = 4, height = 5)
  boxplot(gexp~gtype,data=dt, las=2,ylim=ylim, ylab=paste(genesymbol," (",ylab,")",sep=""),
          xlab=snpid, pars=list(tcl=-0.1), axes=FALSE )
  axis(1,tcl=-0.2, labels = c("1","2","3"), at=c(1,2,3))
  axis(2,tcl=-0.2, las=2)
  x<-sinaplot(dt$gexp,dt$gtype, plot=FALSE)
  points(x$x,x$y, col=adjustcolor("blue", alpha.f = 0.3), bg=adjustcolor("cyan", alpha.f = 0.3), pch=21 )
  legend("topleft", legend = paste("P-value = ",format(pval,digits = 2, scientific = TRUE),sep=""),bty = "n")
  mtext(labpheno,side=3, adj = 1)
  if(plotpdf)dev.off()
}

########################################
### sinaplot code corrected for our analysis!
########################################
sinaplot <- function(x,
                     groups,
                     method = c("density", "neighbourhood"),
                     groupwiseScale = TRUE,
                     yFraction = 0.02,
                     neighbLimit = 1,
                     adjust = 3/4,
                     xSpread = 0.1,
                     labels = NULL,
                     plot = TRUE,
                     #Plot parameters
                     main = "",
                     ylab = "",
                     bw = FALSE,
                     shape = 16,
                     size = 2,
                     color = NULL
)
{
  
  ###Check input arguments
  if (length(x) != length(groups))
    stop("x and groups must be of the same length.")
  
  if (neighbLimit < 1 | !is.numeric(neighbLimit)){
    warning("Invalid neighbLimit value. Neighblimit was set to 1.")
    neighbLimit = 1
  }
  
  if (!is.null(labels) & length(labels) != length(unique(groups)))
    stop("labels and unique 'groups' values must be of the same length.")
  
  if (yFraction <= 0)
    stop("yFraction must be > 0")
  
  if (xSpread <= 0)
    stop("xSpread must be > 0")
  
  method <- match.arg(method)
  ###end
  
  #remove redundant labels
  groups <- factor(groups)
  
  yBins <- .binY(x, yFraction)
  
  #calculate new x coordinates
  x <- .getXcoord(x, groups, yBins, xSpread, groupwiseScale, neighbLimit,
                  adjust, method)
  
  if (is.null(labels))
    labels <- levels(groups)
  
  newGroups <- factor(rep(levels(groups), unlist(lapply(x, nrow))))
  
  #unlist
  x <- do.call(rbind, x)
  
  x$groups <- newGroups
  
  #order the data frame based on the input order
  x <- x[order(x$idx), ]
  x <- x[, c("x", "y", "groups")]
  x
}

.binY <- function(data, yFraction) {
  #get y value range
  ymin <- min(data)
  ymax <- max(data)
  
  #window width
  window_size <- (ymax - ymin) * yFraction
  
  yBins <- c()
  for (i in 0:ceiling(1/yFraction)) {
    yBins <- c(yBins, ymin + i * window_size)
  }
  
  yBins
}

.getXcoord <- function(data, groups, yBins, xSpread, groupwiseScale,
                       neighbLimit, adjust, method){
  
  #number of groups
  ngroups <- length(unique(groups))
  
  xyArray <- c()
  
  ###find the densiest area in the plot
  neighbours <- list()
  
  #set maxNeighbours
  maxNeighbours <- 0
  
  #method == "density"
  if (method == "density"){
    maxDensity <- 0
    densities <- c()
  }
  
  #keep an index of the original data order
  idx <- 1:length(data)
  
  for (j in 1:ngroups){
    
    #extract samples per group and store them in a data.frame
    keep <- groups == levels(groups)[j]
    cur_xyArray <- as.data.frame(cbind(rep(j, sum(keep)),
                                       as.numeric(data[keep]), idx[keep]))
    colnames(cur_xyArray) <- c("x", "y", "idx")
    
    #find the densiest neighbourhood in the current group and compare it
    #with the global max
    cur_neighbours <- table(findInterval(cur_xyArray$y, yBins))
    cur_maxNeighbours <- max(cur_neighbours)
    
    if (cur_maxNeighbours > maxNeighbours)
      maxNeighbours <- cur_maxNeighbours
    
    if (method == "density"){
      #find the highest density value
      cur_density <- density(cur_xyArray$y, adjust = adjust)
      cur_maxDensity <- max(cur_density$y)
      
      if (cur_maxDensity > maxDensity)
        maxDensity <- cur_maxDensity
      
      densities <- c(densities, list(cur_density))
    }
    
    #store neighbour and sample per group data frames in lists
    neighbours <- c(neighbours, list(cur_neighbours))
    xyArray <- c(xyArray, list(cur_xyArray))
  }
  
  relScalingFactor <- 1
  
  for (j in 1:ngroups){
    
    #confine the sample
    if (method == "density"){
      if (max(densities[[j]]$y) > 0.48)
        scalingFactor <- 0.48 / max(densities[[j]]$y)
      else
        scalingFactor <- 1
    }else {
      #if the space required to spread the samples in a neighbourhood exceeds
      #1, create a scaling factor to compress the points
      if (max(neighbours[[j]]) > 1/xSpread){
        scalingFactor <- (1/xSpread) / max(neighbours[[j]])
      }else
        scalingFactor <- 1
    }
    
    #scale all neighbourhoods based on their density relative to the
    #densiest neighbourhood
    if (groupwiseScale == TRUE)
      relScalingFactor <- max(neighbours[[j]]) / maxNeighbours
    
    for (i in names(neighbours[[j]])){
      #examine neighbourhoods with more than 'neighbLimit' samples
      if (neighbours[[j]][i] > neighbLimit && i<length(yBins)){
        cur_bin <- yBins[ as.integer(i) : (as.integer(i) + 1)]
        
        #find samples in the current bin and translate their X coord.
        points <- findInterval(xyArray[[j]]$y, cur_bin) == 1
        
        #compute the border margin for the current bin
        if (method == "density")
          xMax <- mean(densities[[j]]$y[findInterval(densities[[j]]$x,
                                                     cur_bin) == 1])
        else
          xMax <- xSpread*neighbours[[j]][i] / 2
        
        #assign the samples uniformely within the specified range
        xTranslate <- runif(neighbours[[j]][i], -xMax, xMax )
        
        #store new x coordinates
        xyArray[[j]]$x[points] <- xyArray[[j]]$x[points] +
          (xTranslate * scalingFactor * relScalingFactor)
      }
    }
  }
  xyArray
}
