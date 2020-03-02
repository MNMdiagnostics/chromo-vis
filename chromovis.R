#!/usr/bin/env Rscript

###
#
# Parse CLI options
#

suppressMessages(require(optparse))
option_list = list(
  make_option(c("-v", "--vcf"), type="character", default=NULL, 
              help="VCF file (required)", metavar="vcf"),
  make_option(c("-r", "--roh_bed"), type="character", default=NULL, 
              help="BED file with regions of homozygozity", metavar="bed"),
  make_option(c("-o", "--output_prefix"), type="character", default=NULL, 
              help="Output prefix name [default=%input_vcf]", metavar="prefix"),
  make_option(c("-c", "--chromosomes"), type="character", 
              default=paste0("chr",c(as.character(c(1:22)), "X","Y"),collapse=","), 
              help="Comma separated list of chromosomes[default=\"%default\"]", metavar="chr")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# test if there is at least the VCF argument: if not, return an error
if (is.null(opt$vcf)) {
  stop("Missing required VCF argument. Try --help for usage instructions.", call.=FALSE)
} 

# get dirname of the script
cmd = unlist(strsplit(opt_parser@usage, " "))[2]
script.dirname = dirname(cmd)
if (!startsWith(cmd, "/")) {
  script.dirname = paste(getwd(), script.dirname, sep="/")
}

if (is.null(opt$output_prefix)) {
  opt$output_prefix <- opt$vcf 
}

#######
#
# Script logic
#
#

suppressMessages(require(VariantAnnotation))
suppressMessages(require(BSgenome))

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

plot.chrom <- function(chrom, centr.gr, vcf.path, image.path, roh.bed=NA, genome="hg38") {
  chr.gr   <- get.chromosome.ranges(chrom, genome)
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


##################
#
# Plotting
#
########


for (chrom in unlist(strsplit(opt$chromosomes, ","))) {
  png = paste0(opt$output_prefix,".",chrom,".png")
  cat(paste0("Plotting chromosome ", chrom, " to file ", png, "\n"))
  centr.gr <- get.ranges.from.file(paste0(script.dirname,"/centromers_hg38.bed"), chrom=chrom)
  plot.chrom(chrom, centr.gr, opt$vcf, png, roh.bed = opt$roh_bed)  
}
