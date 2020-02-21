
require(VariantAnnotation)
require(BSgenome)

get.ranges.from.file <-function(file, chrom=NA, extra.colnames=c()) {
  t <- read.table(file)
  colnames(t) = c(c("chr","start","end"), extra.colnames)
  if (!is.na(chrom)) {
    t <- t[t$chr==chrom,]
  }
  gr <- GRanges(seqnames = t$chr, IRanges(start=t$start, end=t$end))
  mcols(gr) <- t[,extra.colnames]
  return(gr)
}

get.roh.ranges <- function(file, chrom=NA, size.thr=100000) {
  gr = get.ranges.from.file(file, chrom, c("LEN","CNT","SCORE"))
  return(gr[which(gr$LEN > size.thr)])
}

get.chromosome.ranges <- function(chrom=NA, genome="hg38") {
  if (genome!="hg38") {
    Exception("Not implemented")
  } 
  
  g  <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")
  gr <- GRanges(seqnames=seqnames(g),IRanges(start=0, end=seqlengths(g)))
  
  if (!is.na(chrom)) {
    return(gr[which(seqnames(gr)==chrom)])
  }
  return(gr)
}

running.mean <- function(a, bin.width=1000) {
  b=bin.width/2
  return(sapply(c(1:length(a)), function(x) mean(a[max(0,(x-b)):min(length(a),(x+b))])))
}

plot.chrom <- function(chrom, vcf.path, image.path, roh.bed=NA, genome="hg38") {
  chr.gr   <- get.chromosome.ranges(chrom, genome)
  centr.gr <- get.ranges.from.file("centromers_hg38.bed", chrom=chrom)
  v <- readVcf(TabixFile(vcf.path), genome, param=ScanVcfParam(which=chr.gr, fixed="FILTER", geno=c("DP","AD")))
  # filter PASS variants
  v <- v[fixed(v)$FILTER == "PASS"]
  
  dp.track <- running.mean(geno(v)$DP, 1000)
  rd = sapply(geno(v)$AD, function(x) x[1])
  ad = sapply(geno(v)$AD, function(x) x[2])
  
  png(image.path, width=1200, height = 800, pointsize=16)
  par(mar=c(c(5.1,4.1,4.1,4.1)))
  plot(x=start(v), y=ad/(ad+rd), cex=0.2, pch=18, col="blue",
       xlim=c(0,max(start(v))), ylim=c(0,1), ylab='AD/DP', xlab='position')
  
  #add DP
  max.dp = max(dp.track)
  par(new=T)
  plot(x=start(v), y=dp.track/100, col=2, type='l', lwd=3, 
       axes=F, xlab=NA, ylab=NA, ylim=c(0,1), xlim=c(0, max(start(v))))
  axis(side=4, at = c(0:10)/10, labels = c(0:10)*10, 
       las=2, col=2, col.ticks = 2, col.axis=2)
  mtext(side = 4, line = 3, 'DP', col=2)
  
  #add centromeres
  abline(v=start(centr.gr)[1], col="orange", lwd=3, lty=3)
  abline(v=end(centr.gr)[2], col="orange", lwd=3, lty=3)
  
  
  if (!is.na(roh.bed)) {
    roh.gr   <- get.roh.ranges(roh.bed, chrom)
    arrows(x0=start(roh.gr), y0=1.01, x1=end(roh.gr), y1=1.01, 
           length = 0, lwd = 4, col = "darkgreen")
  }
  title(chrom)
  
  dev.off()
}

setwd('~/Work/projects/MNM/repos/variant_vis/')
output.prefix='~/Work/projects/MNM/repos/variant_vis/plots/'
chroms <- paste0("chr", c(1:22,"X","Y"))
for (chrom in chroms) {
  plot.chrom(chrom, 'test.vcf.gz', paste0(output.prefix,chrom,".png"), 
             roh.bed = 'test.norm.roh.bed')  
}


