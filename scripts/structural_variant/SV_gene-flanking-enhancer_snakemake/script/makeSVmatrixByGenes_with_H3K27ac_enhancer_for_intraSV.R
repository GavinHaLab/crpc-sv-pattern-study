# Use this for short SV events

# Overlapping definition: consider SVs fully contains H3k27ac peaks. 
# Procedure : 
# 	1. Filter H3k27ac peaks within gene region : findoverlap between flanking region (1MB up/down stream) & H3k27ac peaks (within) 
# 	2. Filter SVs within gene region : findoverlap between SV breakpoint and flanking region ("any" for single and "within" for double) 
# 	3. Findoverlap between SV arcs and filtered H3K27ac peaks (from 1) ; peaks must be fully contained within the SV interval (within)
# 	4. Get common query index (bed file) between (2) and (3) for each single and double case


library(GenomicRanges)
library(data.table)
library(stringr) 
library(dplyr)

args <- commandArgs(TRUE)
options(stringsAsFactors=F, width=160)

inBedFile <- args[1] #gene filtered bedfile "short_sv_minus_gene_body_transection_full.txt"
regionFile <- args[2] 
chrEnd_rds <- args[3]
H3K27ac_peaks <- args[4]
outPrefix <- args[5] 

transectType <- "single" 


outImageAnalysis <- paste0("./results/shortSV_enh_analysis/",outPrefix,"_",transectType,".RData")

genomeStyle <- "UCSC"
genomeBuild <- "hg38"

chrEnds <- readRDS(chrEnd_rds)

seqinfo <- Seqinfo(genome=genomeBuild)
chrs <- c(1:22,"X")
seqlevelsStyle(chrs) <- genomeStyle

annotTypes <- c("EnhancerUp", "EnhancerDown")#args[3] # genes or tss
sv.type.enhancer <- c("TandemDup"=1,"Deletion"=2,"Inversion-Unbal"=3, "Inversion-Bal"=4, "Inversion-FoldBack"=5, "Trans-Unbal"=6,"Trans-Bal"=7,"Unbalanced"=8,"Balanced"=9,"LR-shortDup"=10)

met_H3K27ac_peaks = read.table(H3K27ac_peaks, sep="\t")[,1:3]
colnames(met_H3K27ac_peaks) = c("chr","start","end")
met_H3K27ac_peaks.gr = makeGRangesFromDataFrame(met_H3K27ac_peaks)

buffer <- 1e6

inBed <- fread(inBedFile)
samples <- unique(inBed$Sample)
numSamples <- length(samples)

#gene annotation for flanking region (focusing on driver genes only)
regions <- fread(regionFile)
regions <- regions[chromosome %in% paste0("chr",c(1:22,"X")),]

regions <- regions[symbol != "AREnhancer",]
regions.uniq <- unique(regions$symbol)
numGenes <- length(regions.uniq)

enhancerUpMat <- matrix(NaN, nrow=numSamples, ncol=numGenes, dimnames=list(samples, regions.uniq))
enhancerDownMat <- matrix(NaN, nrow=numSamples, ncol=numGenes, dimnames=list(samples, regions.uniq))


save.image(outImageAnalysis)



bed.all <- NULL


for (annotType in annotTypes){
	message("annotType : ", annotType, ", transectType : ",transectType)

	#expand gene region by expanding -1MB for EnhancerUp
	if (annotType == "EnhancerUp"){ #flankingLeft
		regions[, new_end := start-1]
		regions[, new_start :=  as.integer(new_end - buffer)]
		regions[new_start < 0, new_start := 1]	

	#expand gene region by expanding +1MB for EnhancerDown
	} else if (annotType == "EnhancerDown"){ #flankingRight
		regions[, new_start := end+1]
		regions[, new_end :=  as.integer(new_start + buffer)]
		regions[new_end < 0, new_end := 1]
	}

	region_new = regions[,c("chromosome","new_start","new_end","symbol")]
	region_coord = region_new %>% 
				merge(chrEnds, by.x="chromosome", by.y="Chromosome")
	region_coord[V1 < new_end, new_end := V1]


	for (i in samples){

		## load sample bedpe files ##
		message("Analyzing sample ", i)
		bed <- inBed[Sample==i,]
		setnames(bed, c("chrom1","start1","start2"), c("chr","start","end"))
		print(bed[end<start,])
	
		#if we don't consider Translocation like events
		bed.combined <- bed

		bed.combined.gr <- as(bed.combined, "GRanges")

		for(r in regions.uniq){

			print(paste0("process : ",r))
			genes <- region_coord[symbol == r]
			genes.gr <- as(genes, "GRanges")

			peaks_filter <- findOverlaps(subject = genes.gr, query = met_H3K27ac_peaks.gr, type="within")
			peaks_filter_idx <- unique(queryHits(peaks_filter)) 
			
			##any contains even when entire bed region covers gene regions
			genes_to_bed_lookup_table <- findOverlaps(subject = genes.gr, query = bed.combined.gr, type="any")
			beds_containing_entire_gene_region <- findOverlaps(subject=bed.combined.gr, query=genes.gr, type="within")
			
			#swap subject and query order from beds_containing_entire_gene_region
			beds_containing_entire_gene_region_swap  <- cbind(subjectHits(beds_containing_entire_gene_region),queryHits(beds_containing_entire_gene_region))
			colnames(beds_containing_entire_gene_region_swap) = c("queryHits","subjectHits")
			#set difference between genes_to_bed_lookup_table and beds_containing_entire_gene_region_swap
			setdiff_genes_to_bed_lookup_table <- fsetdiff(as.data.table(genes_to_bed_lookup_table), as.data.table(beds_containing_entire_gene_region_swap))
			
			#get SV arcs within gene region
			#the query interval must be wholly contained within the subject interval
			bed_within_table <- findOverlaps(subject = genes.gr, query = bed.combined.gr, type="within")
			bed_within <- unique(queryHits(bed_within_table)) #when both arcs are within genes
			bed_without <- unique(queryHits(genes_to_bed_lookup_table)) 
			bed_contains <- unique(subjectHits(findOverlaps(subject = bed.combined.gr, query = met_H3K27ac_peaks.gr[peaks_filter_idx], type="within")))
			#}

			##########################################
			## get genes that overlap a breakpoint - transecting ##
			## Genes that are broken by break1 or break2 ##
			##########################################
			if (transectType == "single"){
					common_hit	<- intersect(bed_without, bed_contains)
					hits <- setdiff_genes_to_bed_lookup_table[queryHits %in% common_hit]

			##########################################
			## get genes that overlap both breakpoint - double transecting ##
			## Genes that are broken by break1 AND break2 ##
			##########################################
			} else if (transectType == "double"){
					common_hit	<- intersect(bed_within, bed_contains)
					hits <- as.data.table(bed_within_table[queryHits(bed_within_table) %in% common_hit,])

			######################################################
			## get genes that are fully covered by intra-chr SV ##
			## Genes that are contained between break1 and break2 (intrachr) ##
			######################################################
			}else if (transectType == "span"){ ## not working yet - bed.intra
				bed.intra <- bed[chromosome_1==chromosome_2]
				setnames(bed.intra, c("chrom1","start1","start2"), c("chr","start","end"))
				bed.intra.gr <- as(bed.intra, "GRanges")
				hits <- findOverlaps(query = genes.gr, subject = bed.intra.gr, type="within")
				hits <- as.data.table(hits)
			}

			colName <- paste0(annotType, ".", transectType)
			bed <- copy(bed.combined)
			
			if (nrow(hits) > 0){
				message("Assign ", annotType, " overlaps")
				hits <- hits[, paste0(unique(genes[subjectHits, symbol]), collapse=","), by=queryHits]

				bed[hits$queryHits, eval(colName) := hits$V1] #original
				bed_hit = cbind(na.omit(bed, cols=colName))
				print(bed_hit)		
				bed.all <- rbind(bed.all, bed_hit,fill=TRUE)
				#bed[bed.intra[hits$subjectHits, ID], eval(colName) := hits$V1]	
			}else{
				bed[, eval(colName) := NA]
			}
			
			######################################################
			## build matrix for transecting events
			######################################################
			colName <- paste0(annotType, ".", transectType)
			genes.overlap <- bed[!is.na(get(colName)), get(colName)]
			if (length(genes.overlap) > 0){
				svType <- bed[!is.na(get(colName)), CN_overlap_type]
				print(svType)
				span <- bed[!is.na(get(colName)), SPAN]
				for (g in 1:length(genes.overlap)){ ## each SV that overlaps a gene
						svNum <- sv.type.enhancer[svType[g]]

					if (!is.na(svNum)){
						gene.set <- strsplit(genes.overlap[g], ",")[[1]]
						switch(annotType,

							EnhancerUp={ enhancerUpMat[i, gene.set] <-  paste(enhancerUpMat[i, gene.set], svNum, sep="_") },
							EnhancerDown={ enhancerDownMat[i, gene.set] <- paste(enhancerDownMat[i, gene.set],svNum,sep="_") }

						)	

					}
				}
				}
			}	
		}	
	}

## output files ##

save.image(outImageAnalysis)

outFile <- paste0("./results/shortSV_enh_analysis/",outPrefix, "_summary_H3K27ac_", transectType, ".bedpe")

if(!is.null(bed.all)){
	fwrite(bed.all, file=outFile, col.names=T, row.names=F, quote=F, sep="\t")
}


# output the matrix to file
#Input: Matrix to output; output file name
writeMatrixToFile <- function(mat,outfile){
	outMat <- cbind(rownames(mat),mat)
	if (!is.null(colnames(outMat))){
		colnames(outMat)[1] <- "Sample"
	}
	write.table(outMat,file=outfile,row.names=T,col.names=T,quote=F,na="NaN",sep="\t")
}


outFile <- paste0("./results/shortSV_enh_analysis/",outPrefix, "_enhUpMatrix_H3K27ac_", transectType, ".txt")
writeMatrixToFile(enhancerUpMat, outFile)

outFile <- paste0("./results/shortSV_enh_analysis/",outPrefix, "_enhDownMatrix_H3K27ac_", transectType, ".txt")
writeMatrixToFile(enhancerDownMat, outFile)


outImageMats <- paste0("./results/shortSV_enh_analysis/",outPrefix, "_mats_H3K27ac_", transectType, ".RData")

save(enhancerUpMat, enhancerDownMat, file=outImageMats)
