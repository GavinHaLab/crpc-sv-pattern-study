## This version of the script is designed to work for input generated using snakemake: combineSvabaTitan.snakefile.

## requires R-3.3 #
library(VariantAnnotation)
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(ggrepel)


args <- commandArgs(TRUE)

options(stringsAsFactors=F, width=175, scipen=999, bitmapType="cairo")

svDir <- args[1]  # path to combineSVABAandTITAN
minSPAN <- args[2] #1e+07
sampleFile <- args[3]
scoresFile <- args[4] # 0 if no score to plot 
keyFile <- args[5] # will not be used if scoresFile is 0
geneFile <- args[6] # will not be used if scoresFile is 0
outPrefix <- args[7]

outImage <- paste0(outPrefix, ".RData")

use.SV.for.NNI <- TRUE ### USES SV FOR COMPUTING NNI INSTEAD OF SEGS

cnsegs.suffix <- ".svabaTitan.cn.txt"
sv.suffix <- ".svabaTitan.sv.txt"


chrs <- paste0("chr", c(1:22, "X"))

segFiles <- list.files(svDir, pattern = cnsegs.suffix, recursive=T, full.names=T)
ids <- gsub(cnsegs.suffix, "", basename(segFiles))
names(segFiles) <- ids


svFiles <- list.files(svDir, pattern = sv.suffix, recursive=T, full.names=T)
ids <- gsub(sv.suffix, "", basename(svFiles))
names(svFiles) <- ids


samples <- fread(sampleFile)
svFiles <- svFiles[samples[[1]]]
segFiles <- segFiles[samples[[1]]]
numSamples <- length(svFiles)

save.image(file = outImage)

segsAll <- NULL
svAll <- NULL

for (i in 1:numSamples){ 
  id <- names(svFiles)[i]
  message("Loading data for ", id)

  ## load copy number data from TITAN ##
  segs <- fread(segFiles[i])
  segs <- segs[Chromosome %in% chrs]
  segs[, Length.snp. := Length.snp.]
  segs[, length := End - Start + 1]
  segs[, length.SV := SV_start_2 - SV_start_1 + 1]
  segs[, Segment_Mean := Median_logR]
  segs[, ID := Sample]
  
  segs <- segs[, .(ID, SEG.id, Chromosome, Start, End, Length.snp., length, length.SV,
                  Segment_Mean, Median_Ratio, 
                  Copy_Number, MinorCN, Clonal_Cluster, Cellular_Prevalence,
                  Corrected_Copy_Number, SV_overlap_type, SV.overlap.id, SV_start_1, SV_start_2)]
  
  ## load sv data from SvABA ##
  svSampleFile <- grep(id, svFiles, value=T)

  if (length(svSampleFile) == 1){
  	sv <- fread(svSampleFile)
    if (nrow(sv) > 0){
			setnames(sv, "Sample", "ID")
			sv <- sv %>% select(ID, everything())
			sv <- sv[SPAN >= minSPAN]
		}
    svAll <- rbind(svAll, sv)
  }else{
  	stop("Can't find SV file for ", id)
  }
  
  # store all filtered segments into one data.table
  segsAll <- rbind(segsAll, segs)
  
} # end for each sample

segs <- copy(segsAll)
segs[, SEG.id := 1:nrow(segs)]
setkey(segs, SEG.id)
sv <- copy(svAll)
sv[, SV.all.id := 1:nrow(sv)]


## USE SV EVENTS TO COMPUTE NNI 
## DIFFERENT THAN IN VISWANATHAN, HA ET AL. CELL, 2018. WHICH USES SEGS
## get SV that may be tandem dups ##
## renames columns in sv and then assigned to segs so that it can use the exact same code below
# pre-condition: only intra-chr events
# assume start_1 < start_2
if (use.SV.for.NNI){
	setnames(sv, c("chromosome_1","start_1","start_2","SPAN","SV.all.id","CN_overlap_type"), c("Chromosome","Start","End", "length.SV","SEG.id","SV_overlap_type"))
	segs <- copy(sv[, TD := SV_overlap_type == "TandemDup" ])
}

## analyze dispersion ##
# assumes sorted by coordinate within chromosomes
# get minimum distance to left and right event (passing size threshold before)
## assign distance to beginning or end of chromosomes ##
segs[, minDistToChrEnd := {chrStart = Start[1]
            chrEnd = End[.N]
            pmin(Start - chrStart, chrEnd - End)
            }, by=c("Chromosome", "ID")]
## filter segs.gain ##
segs.gain <- segs[SV_overlap_type=="TandemDup"]
segs.gain[, interDupDist := {
            next.start = c(Start[-1], NA) - End          
            }, by=c("Chromosome", "ID")]
# assign first and last event of each chr with dist to end of chr
segs.gain[is.na(interDupDist), interDupDist := minDistToChrEnd]
segs[segs.gain$SEG.id, interDupDist := segs.gain$interDupDist]


## nearest neighbor index (NNI) ##
# need to compute from original seg dt because it has full chromosome coordinates
segs[, NNI := {
          nni <- NA
          if (sum(!is.na(interDupDist)) > 1){
            nni <- mean(interDupDist, na.rm=T) / (0.5*((End[.N]-Start[1])/sum(!is.na(interDupDist))))
          }
          as.numeric(nni)
          }, by=c("Chromosome","ID")]
segs[is.na(NNI), NNI := 0]
# reassign NNI to segs.gain
segs.gain[, NNI := segs[segs.gain$SEG.id, NNI]]


#setkey(ulp.gain.all, ID, chrom)
## summary by sample ##
segs[SV_overlap_type!="TandemDup"|is.na(SV_overlap_type), length := NA]
segs[SV_overlap_type!="TandemDup"|is.na(SV_overlap_type), length.SV := NA]
segs.gain.summary <- segs[, .(Num.Gains=sum(SV_overlap_type=="TandemDup",na.rm=T), Median.Length=mean(length, na.rm=T), Median.Length.SV=mean(length.SV,na.rm=T), Median.InterDupDist=as.double(median(interDupDist,na.rm=T))), by=ID]

## summary by sample ##

## summarize by chromosome and sample
# group by ID and chrom
segs.gain.per.chr <- segs[, .(Num.Gains.per.chr=sum(SV_overlap_type=="TandemDup",na.rm=T), Median.Length.per.chr=as.numeric(median(length,na.rm=T)), Median.Length.SV.per.chr=as.numeric(median(length.SV,na.rm=T)), Median.InterDupDist.per.chr=as.numeric(median(interDupDist,na.rm=T)), NNI=as.numeric(mean(NNI,na.rm=T))), by=c("Chromosome","ID")]
# per chr column names
perChrCols <- c("Num.Gains.per.chr", "Median.Length.per.chr", "Median.InterDupDist.per.chr", "NNI")

segs.gain.per.chr.summary <- segs.gain.per.chr[, .(Num.Gains.per.chr=sum(Num.Gains.per.chr) / length(chrs), 
    Median.Length.per.chr=median(Median.Length.per.chr,na.rm=T), 
    Median.Length.SV.per.chr=median(Median.Length.SV.per.chr,na.rm=T),
    Median.InterDupDist.per.chr=median(Median.InterDupDist.per.chr,na.rm=T),
    NNI=sum(NNI)/length(chrs)), by=ID]
segs.gain.summary[, (perChrCols) := segs.gain.per.chr.summary[, perChrCols, with=F]]

# add back to original summary which is grouped by ID
#segs.gain.summary[order(NNI, decreasing=T)][Median.Length <= 3.5e6 & Num.Gains.per.chr >= 2 & Median.InterDupDist < 20e6,]

outFile <- paste0(outPrefix, "_summary.txt")
write.table(segs.gain.summary, file=outFile, col.names=T, row.names=F, quote=F, sep="\t")

outFile <- paste0(outPrefix, "_TDsegs.txt")
write.table(segs, file=outFile, col.names=T, row.names=F, quote=F, sep="\t")

save.image(file = outImage)

mat <- copy(segs.gain.summary)
if (use.SV.for.NNI){
	mat[, Median.Length := Median.Length.SV]
}
mat[, ID.label := ID]
mat[is.na(Num.Gains), Num.Gains := 0]; mat[is.na(NNI), NNI := 0]; 
mat[is.na(Num.Gains.per.chr), Num.Gains.per.chr := 0]
mat[is.na(Median.Length.SV), Median.Length.SV := 0]
gp <- ggplot(mat, aes(x=Median.Length.SV, y=NNI, size=Num.Gains.per.chr)) +
      geom_point(alpha = 0.5) +  #scale_y_continuous(trans="log10") +
      geom_text_repel(aes(label=ID.label), size=3) +
      xlab("Median TD Segment Length (Mb)") +
      ylab("TDP Dispersion Score (NNI)") +
      theme_bw() + 
      theme(axis.text.y=element_text(size=20), axis.text.x=element_text(size=16),
          axis.title=element_text(size=20), 
          plot.title=element_text(size=24, face="bold")) 
outPlot <- paste0(outPrefix, "_nni-vs-length.pdf")
ggsave(gp, file=outPlot, width=8, height=7)

mat <- copy(segs.gain.summary)
if (use.SV.for.NNI){
	mat[, Median.Length := Median.Length.SV]
}
mat[, ID.label := ID]
mat[is.na(Num.Gains), Num.Gains := 0]; mat[is.na(NNI), NNI := 0]; 
mat[is.na(Num.Gains.per.chr), Num.Gains.per.chr := 0]
mat[is.na(Median.Length), Median.Length := 0]

gp <- ggplot(mat, aes(x=Median.Length, y=Num.Gains, size=Num.Gains.per.chr)) +
      geom_point(alpha = 0.5) +  #scale_y_continuous(trans="log10") +
      geom_text_repel(aes(label=ID.label), size=3) +
      ylab("Number of TD Segments") +
      xlab("Median TD Segment Length (Mb)") +
      scale_x_continuous(trans="log10", limits=c(1e4,10e6),
         breaks=c(1e4,1e5,1e6,2e6,3e6,4e6,5e6,10e6), labels=c(0.01,0.1,1,2,3,4,5,10)) +
      scale_y_continuous(trans="log10") +
      theme_bw() + 
      theme(axis.text.y=element_text(size=20), axis.text.x=element_text(size=16),
          axis.title=element_text(size=20), 
          plot.title=element_text(size=24, face="bold")) 
outPlot <- paste0(outPrefix, "_numTD-vs-length.pdf")
ggsave(gp, file=outPlot, width=9, height=7)


