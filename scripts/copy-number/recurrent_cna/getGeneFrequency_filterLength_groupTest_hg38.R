
library(data.table)
library(GenomicRanges)
library(plyr)
library(ggplot2)
library(plyranges)
library(dplyr)
library(ggrepel)


options(stringsAsFactors=FALSE)

args <- commandArgs(TRUE)

matFile <- args[1]
sampleList <- args[2]
annotFile <- args[3] 
groupType <- args[4]
headerType <- args[5] # Gene or chrPosn
lengthThreshold <- as.numeric(args[6])
cnaType <- args[7] #{'extreme','overall'}
p.adjustMethod <- args[8]
gisticFile <- args[9]
geneFile <- args[10]
outRoot <- args[11]
outImage <- paste0(outRoot, ".RData")
#save.image(outImage)

genomeStyle <- "UCSC"
genomeBuild <- "hg38"
seqinfo <- Seqinfo(genome=genomeBuild)
numCores <- 10
chrs <- c(1:22, "X")
clonalThres <- 0.85
signifLevel <- 0.05
ylim <- c(0,6)

load(matFile)# loads all gene matrices
geneCN <- geneCNmat
geneLen <- geneLenmat

#print(head(geneLen))

genes <- fread(annotFile)
if (headerType == "Gene"){
	setnames(genes, c("cdsStart", "cdsEnd"), c("Start", "End"))
	genes <- genes[, .(Chr = Chr[1], Start = min(Start), End = max(End), Karyotype_band=Karyotype_band[1], strand=strand[1]), by=Gene]
	genes[, Length := End - Start]
	genes <- genes[Length > 1]
	genes[, chrPosn := paste0(Chr,":", Start, "-", End)]
}else{ #if (!headerType %in% colnames(genes)){
   colnames(genes)[1:3] <- c("Chr", "Start", "End")
}
seqlevelsStyle(genes$Chr) <- genomeStyle
seqlevelsStyle(chrs) <- genomeStyle
seqinfo <- seqinfo[chrs]
genes <- genes[Chr %in% chrs]
genes$Chr <- factor(genes$Chr, levels = chrs)
genes <- genes[order(Chr, Start)]
numGenes <- nrow(genes)

# order the genes in the matrices
geneCN <- geneCN[, genes[[headerType]]]
geneLen <- geneLen[, genes[[headerType]]]

numGenes <- ncol(geneCN)
numSamples <- nrow(geneCN)

if (nrow(genes) != numGenes){
	stop("Number of genes in matrices don't match gene list.")
}

## load sample list if given ##
if (sampleList != "0"){
	groups <- fread(sampleList)
}else{
	groups <- data.table(Sample = rownames(geneCN), 1)
	colnames(groups)[2] <- groupType
}
groups <- groups[, c(1,which(names(groups)==groupType)), with=F]
types <- unique(groups[, get(groupType)])
numGroups <- length(types)
samples <- groups[[1]]

geneCN <- geneCN[rownames(geneCN) %in% samples, ]


#save.image(outImage)
####################################################
######## CORRECT COPY NUMBER BY MEDIAN CN ##########
####################################################
auto.chr.genes <- genes[!grep("X", Chr), get(headerType)]
cn.autoChrGenes.median <- t(apply(geneCN[, auto.chr.genes], 1, function(x){ x - median(x, na.rm=T) }))
chrX.genes <- genes[grep("X", Chr), get(headerType)]
cn.chrXGenes.median <- t(apply(geneCN[, chrX.genes], 1, function(x){ x - median(x, na.rm=T) }))
cn.median.norm <- cbind(cn.autoChrGenes.median, cn.chrXGenes.median)
geneCN <- cn.median.norm

####################################################
###### FUNCTION TO COMPUTE FREQUENCY ###############
####################################################
getFrequency <- function(mat, neut = NULL, loss=1, gain=3, na.is.neut = FALSE, numSamples = NULL, ind = TRUE){
	if (is.null(numSamples)){
	  numSamples <- nrow(mat)
	}

	if (is.null(neut)){
		indGain <- (mat >= gain) & ind
		indLoss <- (mat <= loss) & ind
		neut <- (gain - loss)
	}else{
		indGain <- (mat >= neut) & ind
		indLoss <- (mat <= neut) & ind
	}
	if (na.is.neut){
	  #mat[is.na(mat)] <- neut
	  indGain[is.na(indGain)] <- FALSE
	  indLoss[is.na(indLoss)] <- FALSE
	}

	gainFreq <- colSums(indGain, na.rm=TRUE) / numSamples	
	lossFreq <- colSums(indLoss, na.rm=TRUE) / numSamples
	return(list(gainFreq=gainFreq, lossFreq=lossFreq))
}

#gain <- 3
#loss <- 1
# relative copy number to neutral = 0
# the above correction subtracts the median CN from all bins/genes
gain <- 1 
loss <- -1
neut <- NULL
if (cnaType == "extreme"){
	gain <- 3
	loss <- -2	
}

## Overall (clonal + subclonal) copy number frequency ##
freqAll <- getFrequency(geneCN, neut = neut, loss=loss, gain=gain, numSamples = nrow(groups), na.is.neut=T)

## copy number filtered by length threshold ##
freqLen <- getFrequency(geneCN, neut = neut, loss=loss, gain=gain, numSamples = nrow(groups), na.is.neut=T, 
    ind = geneLen[samples, ] <= lengthThreshold)



outMat <- data.table(genes[, 1:3], All_loss=freqAll$lossFreq, All_gain=freqAll$gainFreq)
outFile <- paste(outRoot, "_geneFreq.txt", sep="")
## output to text file ##
write.table(format(outMat,digits=4), file=outFile, col.names=T, row.names=F, quote=F, sep="\t")

outMat$pos <- rowMeans(cbind(outMat$Start, outMat$End))

genes.df <- as.data.frame(genes)
genes.df$Chr <- as.character(genes.df$Chr)
genes.gr <- makeGRangesFromDataFrame(genes.df, keep.extra.columns=TRUE)


#################################################
################ PLOT GENE FREQUENCIES ##########

##load GISTIC data to overlay hits over frequency plot
gistic <- read.table(gisticFile, sep="\t", head=T, stringsAsFactors=F)
gistic_peaks_only <- gistic[!grepl("CN",gistic[,1]),]
cn_type <- unlist(lapply(strsplit(gistic_peaks_only$Unique.Name,' '),'[[',1)) #Amplification or Deletion

Gain_idx <- which(cn_type == "Amplification")
Loss_idx <- which(cn_type == "Deletion")
Peak_label <- gsub(" ","",paste(cn_type, gistic_peaks_only$Descriptor,sep="_"))

## load CRPC gene list to annotate to GISTIC peak
CRPC_gene_list <- read.table(geneFile, sep="\t", head=T)
regions <- makeGRangesFromDataFrame(CRPC_gene_list, TRUE)

peak_type <- "Wide.Peak.Limits"


Peak <- gsub("\\(probes.*$","",gsub(" ", "", gistic_peaks_only[,peak_type]))
Peak_chr <- unlist(lapply(strsplit(Peak,":"),'[[',1))
Peak_start <- as.numeric(unlist(lapply(strsplit(sub(".*:", "", Peak),"-"),'[[',1)))
Peak_end <- as.numeric(unlist(lapply(strsplit(sub(".*:", "", Peak),"-"),'[[',2)))

Gain_regions.gr <- GRanges(seqnames=Peak_chr[Gain_idx], ranges=IRanges(start = Peak_start[Gain_idx], end =Peak_end[Gain_idx]),type=Peak_label[Gain_idx])
Loss_regions.gr <- GRanges(seqnames=Peak_chr[Loss_idx], ranges=IRanges(start = Peak_start[Loss_idx], end =Peak_end[Loss_idx]),type=Peak_label[Loss_idx])

Gain_intersectRegion <- join_overlap_inner(Gain_regions.gr, regions)
Loss_intersectRegion <- join_overlap_inner(Loss_regions.gr, regions)
gene_annotation <- rbind(as.data.frame(Gain_intersectRegion), as.data.frame(Loss_intersectRegion))

## output to text file ##
#write.table(gene_annotation,sprintf("GISTIC_%s_with_gene_annotation.txt", peak_type), sep='\t', col.names=T, row.names=F, quote=F)

## save image ##
#save.image(file=outImage)

############################
#find overlap between amp peak & bins

gistic_intersectGainRegion <- as.data.frame(subsetByOverlaps(genes.gr, Gain_regions.gr))
colnames(gistic_intersectGainRegion)[1:3] <- c("Chr","Start","End")

#find overlap between amp peak related genes and bins
gistic.gain.gr <- makeGRangesFromDataFrame(gistic_intersectGainRegion, TRUE)

Gain.overlap <- findOverlaps(gistic.gain.gr, Gain_intersectRegion)
Gain.matched <- gistic.gain.gr[queryHits(Gain.overlap)]
# Add the metadata 
mcols(Gain.matched) <- cbind.data.frame(
    mcols(Gain.matched),
    mcols(Gain_intersectRegion[subjectHits(Gain.overlap)]))

annot_intersectGainRegion = as.data.frame(Gain.matched)
colnames(annot_intersectGainRegion)[1:3] = c("Chr","Start","End")

#combine gene name for same peak regions
max_end_annot_gain <- as.data.frame( annot_intersectGainRegion %>% 
  group_by(Chr,Start,End,chrPosn,symbol) %>% 
  ungroup() %>%
  group_by(Chr,Start,End,chrPosn) %>%
  summarize(Series = paste(symbol, collapse = ',')) %>%
  group_by(Series) %>% top_n(1, End) )


############################
# draw for Gain
############################

data_plot <- outMat %>%
  group_by(Chr) %>%
  summarize(chr_len=max(End)) %>%
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  left_join(outMat, ., by="Chr") %>%
  arrange(Chr, pos) %>%
  mutate(pos_cs=pos+tot) %>%
  left_join(.,gistic_intersectGainRegion, by=c("Chr","Start","End")) %>%
  left_join(.,max_end_annot_gain, by=c("Chr","Start","End")) 

axis_x <- data_plot %>%
  group_by(Chr) %>%
  summarize(center=.5*(max(pos_cs)+min(pos_cs)), end=max(pos_cs))


mht_plot <- ggplot(data_plot, aes(x=pos_cs, y=All_gain)) +
  geom_vline(xintercept=c(0, axis_x$end), linetype="dashed", color="#88888888") +
  #geom_bar(color="grey60") +
  #geom_point(data=subset(data_plot, score<=fdr_thres), color="red") +
  geom_bar(data=data_plot, color="red", stat='identity') +
  geom_bar(data=data_plot[!is.na(data_plot$chrPosn.x),], stat='identity', color="black") +
 # geom_bar(data=data_plot[!is.na(data_plot$chrPosn),],stat='identity', color="black") +
  scale_x_continuous(name="Chromosome", label=c(1:22, "X"),
                     breaks=axis_x$center, expand=c(0, 0)) +
  scale_y_continuous(name="Frequency",limits = c(0, 1)) +
  geom_text_repel(data=subset(data_plot, !is.na(Series)), aes(label=Series,
                  size=4), box.padding=1.00, point.padding=0.20,
                  nudge_y=0.3, nudge_x=0.5, min.segment.length=0, segment.alpha=0.5) +
  theme_bw() +
  theme(legend.position="none", panel.border = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

pdf(sprintf("CRPC10X_WCDT_100kb_geneFreq_with_GISTIC_%s_Gain.pdf",peak_type), height=4, width=12)
print(mht_plot)
dev.off()


############################
#find overlap between del peak & bins

gistic_intersectLossRegion = as.data.frame(subsetByOverlaps(genes.gr, Loss_regions.gr))
colnames(gistic_intersectLossRegion)[1:3] = c("Chr","Start","End")

#find overlap between del peak related genes and bins
gistic.loss.gr <- makeGRangesFromDataFrame(gistic_intersectLossRegion,TRUE)

Loss.overlap <- findOverlaps(gistic.loss.gr, Loss_intersectRegion)
Loss.matched <- gistic.loss.gr[queryHits(Loss.overlap)]
# # Add the metadata 
mcols(Loss.matched) <- cbind.data.frame(
    mcols(Loss.matched),
    mcols(Loss_intersectRegion[subjectHits(Loss.overlap)]))

annot_intersectLossRegion <- as.data.frame(Loss.matched)
colnames(annot_intersectLossRegion)[1:3] = c("Chr","Start","End")

#combine gene name for same peak regions
max_end_annot_loss <- as.data.frame( annot_intersectLossRegion %>% 
	group_by(Chr,Start,End,chrPosn,symbol) %>% 
    ungroup() %>%
	group_by(Chr,Start,End,chrPosn) %>%
	summarize(Series = paste(symbol, collapse = ',')) %>%
	group_by(Series) %>% top_n(1, End) )

max_end_annot_loss <- max_end_annot_loss[max_end_annot_loss$Series!="CDH1,ZFHX3,FANCA,FANCA",]

############################
# draw for loss
############################

data_plot <- outMat %>%
  group_by(Chr) %>%
  summarize(chr_len=max(End)) %>%
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  left_join(outMat, ., by="Chr") %>%
  arrange(Chr, pos) %>%
  mutate(pos_cs=pos+tot) %>%
  left_join(.,gistic_intersectLossRegion, by=c("Chr","Start","End")) %>%
  left_join(.,max_end_annot_loss, by=c("Chr","Start","End")) 

axis_x <- data_plot %>%
  group_by(Chr) %>%
  summarize(center=.5*(max(pos_cs)+min(pos_cs)), end=max(pos_cs))


mht_plot <- ggplot(data_plot, aes(x=pos_cs, y=-1*All_loss)) +
  geom_vline(xintercept=c(0, axis_x$end), linetype="dashed", color="#88888888") +
  geom_bar(data=data_plot, color="blue",stat='identity') +
  geom_bar(data=data_plot[!is.na(data_plot$chrPosn.x),],stat='identity', color="black") +
  scale_x_continuous(name="Chromosome", label=c(1:22, "X"),
                     breaks=axis_x$center, expand=c(0, 0)) +
  scale_y_continuous(name="Frequency",limits = c(-1,0)) +
  geom_text_repel(data=subset(data_plot, !is.na(Series)), aes(label=Series,
                  size=4), box.padding=1.00, point.padding=0.20,
                  nudge_x=0.2, nudge_y=-0.4, min.segment.length=1, segment.alpha=0.5) +
  theme_bw() +
  theme(legend.position="none", panel.border = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

pdf(sprintf("CRPC10X_WCDT_100kb_geneFreq_with_GISTIC_%s_Loss.pdf",peak_type), height=4, width=12)
print(mht_plot)
dev.off()


