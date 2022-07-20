###SV breakpoint analysis for ETS family 

library(data.table)
library(GenomicRanges)


buffer <- 0 #no buffer


#use genecode annotation--filter gene by protein coding only
geneFile <- "../structural_variant/SV_gene-flanking-enhancer_snakemake/data/gencode.v33.annotation.genes.tsv"
genes <- fread(geneFile)
protein_coding_genes <- genes [gene_type=="protein_coding",c("gene_name","seqid","start","end","strand")]
protein_coding_genes[ ,Start := min(start)-buffer, by=gene_name]
protein_coding_genes[ ,End :=max(end)+buffer, by=gene_name]
genes_filtered <- unique(protein_coding_genes[,c("seqid","Start","End","gene_name")])
genes.gr <- as(genes_filtered, "GRanges")

#load bed file
bedFile <- "../../metadata/CRPC10X42-PoNToolFilter_manual_edit_WCDT101_svaba-titan.min1kb.bedpe" 
bed <- fread(bedFile)


#prepare bed file
bed.1 <- copy(bed) 
setnames(bed.1, c("chrom1","start1","end1"), c("chr","start","end"))
bed.2 <- copy(bed)
setnames(bed.2, c("chrom2","start2","end2"), c("chr","start","end"))
bed.gr.1 <- as(bed.1[, -c("chrom2","start2","end2")], "GRanges")
bed.gr.2 <- as(bed.2[, -c("chrom1","start1","end1")], "GRanges")

#load ETS gene coordinates (28 genes in total, downloaded ETS gene family from hgnc website)
ets_gene_full <- fread("../../metadata/ETS_gene_coordinates_hg38.txt")

ets_gene_full[,Start := min(txStart)-buffer, by=name2] #get minimum txStart for same gene symbol
ets_gene_full[,End := max(txEnd)+buffer, by=name2] #get maximum txEnd for same gene symbol

ets_gene <- unique(ets_gene_full[,c("chrom","Start","End","name2")])
ets_gene.gr <- as(ets_gene, "GRanges")


### findoverlap between breakpoint and ETS gene
hits.1 <- findOverlaps(subject=ets_gene.gr, query = bed.gr.1)
hits.2 <- findOverlaps(subject=ets_gene.gr, query = bed.gr.2)
ets_gene.gr[unique(subjectHits(hits.1))]
ets_gene.gr[unique(subjectHits(hits.2))]

#find partner of ets gene hit
###bed_hits.1 and bed_hits.2 are subsetted bed based on ETS gene family
bed_hits.1 <- as.data.table(bed.gr.1[queryHits(hits.1)])
bed_hits.1[, eval("gene1") := ets_gene[subjectHits(hits.1)]$name2] #assign gene for pair1
bed_hits.1.pair <- bed_hits.1[,c("chromosome_2","start_2","start_2")]
colnames(bed_hits.1.pair) <- c("chr","start","end")
bed_hits.1.pair.gr <- as(bed_hits.1.pair, "GRanges")

bed_hits.2 <- as.data.table(bed.gr.2[queryHits(hits.2)])
bed_hits.2[, eval("gene1") := ets_gene[subjectHits(hits.2)]$name2]#assign gene for pair1
bed_hits.2.pair <- bed_hits.2[,c("chromosome_1","start_1","start_1")]
colnames(bed_hits.2.pair) <- c("chr","start","end")
bed_hits.2.pair.gr <- as(bed_hits.2.pair, "GRanges")

##find pair2 from bed_hits.1 and bed_hits.2
bed_hits.1.pair_gene <- as.data.table(findOverlaps(subject=bed_hits.1.pair.gr, query = genes.gr))
bed_hits.2.pair_gene <- as.data.table(findOverlaps(subject=bed_hits.2.pair.gr, query = genes.gr))

#mapping to gene info
bed_hits.1.pair_gene <- bed_hits.1.pair_gene[, paste0(unique(genes_filtered[queryHits, gene_name]), collapse=","), by=subjectHits]
bed_hits.2.pair_gene <- bed_hits.2.pair_gene[, paste0(unique(genes_filtered[queryHits, gene_name]), collapse=","), by=subjectHits]

bed_hits.1[bed_hits.1.pair_gene$subjectHits, eval("gene2") := bed_hits.1.pair_gene$V1]
bed_hits.2[bed_hits.2.pair_gene$subjectHits, eval("gene2") := bed_hits.2.pair_gene$V1]

fwrite(bed_hits.1, file=sprintf("bed_hits.1.ETS_family_buffer_%d.txt",buffer), col.names=T, row.names=F, quote=F, sep="\t")
fwrite(bed_hits.2, file=sprintf("bed_hits.2.ETS_family_buffer_%d.txt",buffer), col.names=T, row.names=F, quote=F, sep="\t")




