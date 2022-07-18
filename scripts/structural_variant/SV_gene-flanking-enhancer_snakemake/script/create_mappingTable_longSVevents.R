
# Treat bp1 and bp2 separately for this analysis
# mapping gene info (flanking region--1mb nearby gene body) of bp1 and then find parter genes (within 1mb) of bp2

library(dplyr)
library(tidyr)
library(data.table)
library(GenomicRanges)
library(RCircos)

args <- commandArgs(TRUE)


inBedFile <- args[1] # "long_sv_minus_gene_body_transection_full.txt"
regionFile <- args[2] 
geneAnnotFile <- args[3]
chrEnd_rds <- args[4]
outPrefix <- args[5] 


transectType <- "single"

genomeStyle <- "UCSC"
genomeBuild <- "hg38"
seqinfo <- Seqinfo(genome=genomeBuild)
chrs <- c(1:22,"X")
seqlevelsStyle(chrs) <- genomeStyle


annotTypes <- c("EnhancerUp", "EnhancerDown")
outImageAnalysis <- paste0("./results/longSV_analysis/",outPrefix,"_flank_with_longSV_",transectType,".RData")


inBed <- fread(inBedFile)
samples <- unique(inBed$Sample)

#gene annotation for flanking region
region <- fread(regionFile)

#consider chromosomes 1:22 and chrX only
region <- region[chromosome %in% paste0("chr",c(1:22,"X")),]

chrEnds <- readRDS(chrEnd_rds)

#region <- region[symbol != "AREnhancer",]

buffer <- 1e6 #1mb definition

#####################
## This gene annotation is going to be used for the partners
#####################

#geneAnnotFile <- "gencode.v33.annotation.genes.tsv"
genes <- fread(geneAnnotFile)
protein_coding_genes <- genes [gene_type=="protein_coding",c("gene_name","seqid","start","end","strand")]
protein_coding_genes[ ,Start := min(start), by=gene_name]
protein_coding_genes[ ,End :=max(end), by=gene_name]
genes_filtered <- unique(protein_coding_genes[,c("seqid","Start","End","gene_name")])

##create for up/dn stream
genes_filtered_up = copy(genes_filtered)
genes_filtered_up[, End := min(Start)-1, by=gene_name]
genes_filtered_up[, Start :=  as.integer(End - buffer)]
genes_filtered_up[, gene_name := paste0(gene_name,"_up")]
genes_filtered_up[Start < 0, Start := 1]

genes_filtered_dn = copy(genes_filtered)
genes_filtered_dn[, Start := max(End)+1, by=gene_name]
genes_filtered_dn[, End :=  as.integer(Start + buffer)]
genes_filtered_dn[, gene_name := paste0(gene_name,"_dn")]
genes_filtered_dn[End < 0, End := 1]

genes_filtered_all = rbind(genes_filtered, genes_filtered_up, genes_filtered_dn)

genes.gr <- as(rbind(genes_filtered, genes_filtered_up, genes_filtered_dn), "GRanges")


bed.all <- NULL

for (i in samples){

	## load sample bedpe files ##
	message("Analyzing sample ", i)
	sample_bed <- inBed[Sample==i,]
	
	bed.1 <- copy(sample_bed) 
	setnames(bed.1, c("chrom1","start1","end1"), c("chr","start","end"))
	bed.2 <- copy(sample_bed)
	setnames(bed.2, c("chrom2","start2","end2"), c("chr","start","end"))
	# bed.gr.1 <- as(bed.1[, -c("chrom2","start2","end2")], "GRanges")
	# bed.gr.2 <- as(bed.2[, -c("chrom1","start1","end1")], "GRanges")
	bed.gr.1 <- as(bed.1, "GRanges")
	bed.gr.2 <- as(bed.2, "GRanges")

	
	#Enhancer region only
	for (annotType in annotTypes){
		print(annotType)
		#expand gene region by expanding -1MB for EnhancerUp
		if (annotType == "EnhancerUp"){ #flankingLeft
			region[, new_end := start-1]
			region[, new_start :=  as.integer(new_end - buffer)]
			region[new_start < 0, new_start := 1]	

		#expand gene region by expanding +1MB for EnhancerDown
		} else if (annotType == "EnhancerDown"){ #flankingRight
			region[, new_start := end+1]
			region[, new_end :=  as.integer(new_start + buffer)]
			region[new_end < 0, new_end := 1]
		}

		region_new <- region[,c("chromosome","new_start","new_end","symbol")]
		region_new <- region_new %>% 
					merge(chrEnds, by.x="chromosome", by.y="Chromosome")
		region_new[V1 < new_end, new_end := V1]

		region.gr <- as(region_new[,c("chromosome","new_start", "new_end")], "GRanges")

		#colName <- paste0(transectType)


		### findoverlap between breakpoint and gene
		hits.1 <- findOverlaps(subject=region.gr, query = bed.gr.1)
		hits.2 <- findOverlaps(subject=region.gr, query = bed.gr.2)
		
		hits_list = list()

		#there could be multiple genes involved in SV breakpoints since we are looking at 1MB region 
		#separate the rows for different genes (easier to create gene matrix later on)

		#find partner gene
		if( length(hits.1) >0 ){
			
			bed_hits.1 <- as.data.table(bed.gr.1[queryHits(hits.1)])
			bed_hits.1[, eval("gene") := region_new[subjectHits(hits.1)]$symbol] #assign gene for pair1
#			bed_hits.1.pair <- bed_hits.1[,c("chromosome_2","start_2","start_2")]
			bed_hits.1.pair <- bed_hits.1[,c("chrom2","start2","end2")]
			colnames(bed_hits.1.pair) <- c("chr","start","end")
			bed_hits.1.pair.gr <- as(bed_hits.1.pair, "GRanges")
			
			##find pair2 from bed_hits.1
			bed_hits.1.pair_gene <- as.data.table(findOverlaps(subject=bed_hits.1.pair.gr, query = genes.gr))
			
			#mapping partner to genes
			bed_hits.1.pair_gene <- bed_hits.1.pair_gene[, paste0(unique(genes_filtered_all[queryHits, gene_name]), collapse=","), by=subjectHits]
			bed_hits.1[bed_hits.1.pair_gene$subjectHits, eval("partner") := bed_hits.1.pair_gene$V1]
			
			#trim the data
			#if it matches to start_1 info
			bed_hits.1.sub = unique(bed_hits.1[  ,c("seqnames","start","orient_1","chrom2","start2","orient_2","CN_overlap_type","svID","SPAN","gene","partner")])
			setnames(bed_hits.1.sub, c("seqnames","orient_1","chrom2","start2","orient_2","gene","partner"),
				c("chrom","orient","chrom_partner","start_partner","orient_partner","gene_flanking","gene_partner"))

			hits_list[["hits.1"]] = bed_hits.1.sub


		}
		if( length(hits.2) >0 ){
			
			bed_hits.2 <- as.data.table(bed.gr.2[queryHits(hits.2)])
			bed_hits.2[, eval("gene") := region_new[subjectHits(hits.2)]$symbol]#assign gene for pair1
			bed_hits.2.pair <- bed_hits.2[,c("chrom1","start1","end1")]
			colnames(bed_hits.2.pair) <- c("chr","start","end")
			bed_hits.2.pair.gr <- as(bed_hits.2.pair, "GRanges")
			
			##find pair2 from bed_hits.2
			bed_hits.2.pair_gene <- as.data.table(findOverlaps(subject=bed_hits.2.pair.gr, query = genes.gr))
			
			#mapping partner to genes
			bed_hits.2.pair_gene <- bed_hits.2.pair_gene[, paste0(unique(genes_filtered_all[queryHits, gene_name]), collapse=","), by=subjectHits]
			bed_hits.2[bed_hits.2.pair_gene$subjectHits, eval("partner") := bed_hits.2.pair_gene$V1]
			
			#trim the data
			#if it matches to start_2 info
			bed_hits.2.sub = unique(bed_hits.2[ ,c("seqnames","start","orient_2","chrom1","start1","orient_1","CN_overlap_type","svID","SPAN","gene","partner")])
			setnames(bed_hits.2.sub, c("seqnames","orient_2","chrom1","start1","orient_1","gene","partner"),
				c("chrom","orient","chrom_partner","start_partner","orient_partner","gene_flanking","gene_partner"))

			hits_list[["hits.2"]] = bed_hits.2.sub
		}

		#keep records for all SV flanking gene and SV breakpoint partner gene information

		hits.all <- do.call(rbind, hits_list)
		bed.all <- rbind(bed.all, cbind(i, hits.all, annotType))
	}
}
## output files ##
outFile <- paste0("./results/longSV_analysis/",outPrefix, "_longSV_summary_gene_info_", transectType, ".txt")
colnames(bed.all)[1] = "Sample"
fwrite(bed.all, file=outFile, col.names=T, row.names=F, quote=F, sep="\t")

save.image(outImageAnalysis)
	
